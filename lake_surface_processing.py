import argparse
import glob
import os
from datetime import datetime

import matplotlib.path as mpth
import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset

# --------- NOTE ---------- #
# --- WRF-native + calculated lake surface variables --- #
# --- Only needs T_LAKE3D, LAKE_ICEFRAC3D, LAKEMASK, XLAT, XLONG from wrfout --- #

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SHAPEFILE_DIR = os.path.join(SCRIPT_DIR, "lake_shapefiles")

# LBRM basin outline csvs to use for the 5 named Great Lakes basins.
# (ER/GEO/STC csvs in lake_shapefiles/ are a duplicate Erie subset, the Georgian
# Bay subbasin, and the St. Clair connecting-channel subbasin respectively --
# they are not part of the 5 basins requested and are intentionally unused here.)
BASIN_CSVS = {
    "SUPERIOR": "SUP_lbrm_outline.csv",
    "HURON": "HU_lbrm_outline.csv",
    "MICHIGAN": "MIC_lbrm_outline.csv",
    "ERIE": "ERI_lbrm_outline.csv",
    "ONTARIO": "ON_lbrm_outline.csv",
}
# priority order used to resolve the small number of gridpoints whose basin
# polygons overlap (adjacent-lake boundary cells)
BASIN_PRIORITY = ["SUPERIOR", "HURON", "MICHIGAN", "ERIE", "ONTARIO"]
NONE_LABEL = "NONE"

# order/key for mean_ice_cover_percent's lake_group axis
LAKE_GROUPS = ["FULL_BASIN", "SUPERIOR", "HURON", "MICHIGAN", "ERIE", "ONTARIO", "MICHIGAN-HURON"]

# CLM lake-model fill sentinel used in *_LAKE3D variables at non-lake gridpoints
# (not exposed as a proper _FillValue attribute in the wrfout files)
LAKE3D_FILL_VALUE = -999.0

# name of the lake-level dimension in T_LAKE3D / LAKE_ICEFRAC3D
LAKE_LEVEL_DIM = "soil_levels_or_lake_levels_stag"
# top (surface) layer index of the lake-level dimension
# (confirmed via Z_LAKE3D: level 0 has the shallowest depth, i.e. the surface layer)
SURFACE_LAKE_LEVEL = 0

# areal lake-ice fraction at/above which a gridpoint-day is counted as "ice covered"
# for the categorical ice-day / first-ice-day / last-ice-day statistics
ICE_COVER_THRESHOLD = 0.5

MONTHS = list(range(1, 13))
N_DOY = 365  # leap days are dropped so every year aligns to a 365-day calendar

# wrfout files are 6-hourly (4 timesteps/day); a chunk of 200 timesteps is
# ~50 days per dask chunk. Tune this (and the SLURM --mem for this script) if
# you see it thrashing or under-using available memory.
TIME_CHUNK = 200


def check_year(value):
    ivalue = int(value)
    if ivalue < 2000 or ivalue > 2099:
        raise argparse.ArgumentTypeError(f"{value} is an invalid year. Must be between 2000 and 2099.")
    return ivalue


def parse_file_datetime(fpath):
    timestamp = os.path.basename(fpath).split("_d01_")[1]
    return datetime.strptime(timestamp, "%Y-%m-%d_%H:%M:%S")


def build_clim_doy(time_index):
    """
    1-indexed DOY on a 365-day calendar for every entry in a DatetimeIndex.
    Leap-year days after Feb 28 are shifted back by 1 so all years align.
    Feb 29 entries are marked NaN (leap day ignored per spec, matching the
    original normalized_doy()/per-day skip behavior).
    """
    times = pd.DatetimeIndex(time_index)
    doy = times.dayofyear.to_numpy().astype(float)
    shift = (times.is_leap_year & (times.month > 2)).astype(float)
    doy = doy - shift
    doy[(times.month == 2) & (times.day == 29)] = np.nan
    return doy


def load_basin_masks(lon, lat):
    """
    Return dict[basin_name] -> boolean mask over the (south_north, west_east) grid.

    Each csv packs several dozen disjoint "subbasin" loops (see the `subbasin`
    column) back to back in one long point list. Building a single
    matplotlib.path.Path from the whole file (as lake_basin_segmenting.ipynb
    does) draws a spurious connecting edge between the end of one loop and the
    start of the next, which silently excludes a large fraction of the true
    lake surface from the polygon test. Grouping by `subbasin` and unioning
    each closed loop's own contains_points result avoids that.
    """
    pts = np.column_stack((lon.ravel(), lat.ravel()))
    masks = {}
    for name, fname in BASIN_CSVS.items():
        df = pd.read_csv(os.path.join(SHAPEFILE_DIR, fname))
        inside_any = np.zeros(pts.shape[0], dtype=bool)
        for _, grp in df.groupby("subbasin"):
            poly = grp[["long", "lat"]].to_numpy(dtype=float)
            if len(poly) < 3:
                continue
            path = mpth.Path(poly)
            inside_any |= path.contains_points(pts)
        masks[name] = inside_any.reshape(lon.shape)
    return masks


def build_lake_basin_vars(lon, lat):
    masks = load_basin_masks(lon, lat)
    lake_basin = np.full(lon.shape, NONE_LABEL, dtype="<U8")
    for name in BASIN_PRIORITY:
        unassigned = lake_basin == NONE_LABEL
        lake_basin[masks[name] & unassigned] = name
    lake_basin_flag = (lake_basin != NONE_LABEL).astype(np.int8)
    return lake_basin, lake_basin_flag


def basin_group_masks(lake_basin, lake_basin_flag):
    """dict[group_name] -> boolean mask, in LAKE_GROUPS order."""
    masks = {"FULL_BASIN": lake_basin_flag == 1}
    for name in BASIN_PRIORITY:
        masks[name] = lake_basin == name
    masks["MICHIGAN-HURON"] = (lake_basin == "MICHIGAN") | (lake_basin == "HURON")
    return masks


def _preprocess_lake(ds):
    """Attach a real datetime Time coordinate; keep only the surface lake layer (float32)."""
    times = pd.to_datetime(ds["Times"].astype(str).values, format="%Y-%m-%d_%H:%M:%S")
    ds = ds.assign_coords(Time=("Time", times))
    icefrac = ds["LAKE_ICEFRAC3D"].isel({LAKE_LEVEL_DIM: SURFACE_LAKE_LEVEL})
    tlake = ds["T_LAKE3D"].isel({LAKE_LEVEL_DIM: SURFACE_LAKE_LEVEL})
    out = xr.Dataset({"ICEFRAC_sfc": icefrac, "TLAKE_sfc": tlake})
    return out.astype("float32")


def main():
    parser = argparse.ArgumentParser(description="Process lake surface variables from wrfout files.")
    parser.add_argument(
        "-p",
        "--path",
        help="directory of wrfout files",
        default="../WRF/run",
    )
    parser.add_argument(
        "-y",
        "--year_start",
        type=check_year,
        help="The starting year (integer between 2000-2099)",
        default=2005,
    )
    parser.add_argument(
        "-s",
        "--scenario",
        help="simulation scenario (default: historical)",
        default="historical",
    )
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()
    home_dir = args.path
    start_year = args.year_start
    end_year = start_year + 10

    if os.path.isdir(home_dir):
        print(f"Accessing directory: {home_dir}")
    else:
        print(f"Error: {home_dir} is not a valid directory.")
        return

    all_files = sorted(glob.glob(os.path.join(home_dir, "wrfout_d01_20*")))
    files = [f for f in all_files if start_year <= parse_file_datetime(f).year < end_year]
    if not files:
        print(f"No wrfout files found for years [{start_year}, {end_year}) in {home_dir}")
        return
    print(f"Found {len(files)} wrfout files for years {start_year}-{end_year - 1}")

    sample = Dataset(files[0])
    lat = np.asarray(sample.variables["XLAT"][0, :, :], dtype=float)
    lon = np.asarray(sample.variables["XLONG"][0, :, :], dtype=float)
    water_mask_np = np.asarray(sample.variables["LAKEMASK"][0, :, :]) == 1
    ny, nx = lat.shape
    sample.close()

    print("Building lake basin classification from LBRM basin outlines...")
    lake_basin, lake_basin_flag = build_lake_basin_vars(lon, lat)
    group_masks = basin_group_masks(lake_basin, lake_basin_flag)

    years_range = list(range(start_year, end_year))

    # ── Lazily load surface ICEFRAC/TLAKE across the whole period ──────────
    ds = xr.open_mfdataset(
        files,
        engine="netcdf4",
        combine="nested",
        concat_dim="Time",
        preprocess=_preprocess_lake,
        chunks={"Time": TIME_CHUNK},
        # netCDF4/HDF5 isn't reliably thread-safe for concurrent file opens;
        # parallel=True has caused spurious "Unknown file format" errors.
        # This only affects the (cheap) metadata-open phase, not the
        # dask-parallel computation that follows.
        parallel=False,
    )

    water_mask = xr.DataArray(water_mask_np, dims=("south_north", "west_east"))

    icefrac = ds["ICEFRAC_sfc"].where(ds["ICEFRAC_sfc"] > LAKE3D_FILL_VALUE + 1.0)
    tlake = ds["TLAKE_sfc"].where(ds["TLAKE_sfc"] > LAKE3D_FILL_VALUE + 1.0)
    icefrac = icefrac.where(water_mask)
    tlake = tlake.where(water_mask)

    # ── Monthly climatology (all raw timesteps, Feb 29 included) ───────────
    print("Computing monthly ICEFRAC/TLAKE climatology...")
    icefrac_mean = (
        icefrac.groupby("Time.month").mean(dim="Time", skipna=True)
        .reindex(month=MONTHS).transpose("month", "south_north", "west_east")
    )
    tlake_mean = (
        tlake.groupby("Time.month").mean(dim="Time", skipna=True)
        .reindex(month=MONTHS).transpose("month", "south_north", "west_east")
    )
    tlake_p95 = (
        tlake.groupby("Time.month").quantile(0.95, dim="Time", skipna=True)
        .reindex(month=MONTHS).transpose("month", "south_north", "west_east")
    )
    tlake_p99 = (
        tlake.groupby("Time.month").quantile(0.99, dim="Time", skipna=True)
        .reindex(month=MONTHS).transpose("month", "south_north", "west_east")
    )

    # ── Day-aggregated (Feb 29 excluded) series for ice-day / basin-pct stats ──
    print("Computing daily-mean ice fraction and ice-day statistics...")
    day_mean_icefrac = icefrac.resample(Time="1D").mean(skipna=True)

    doy_all = build_clim_doy(day_mean_icefrac["Time"].values)
    valid_idx = np.nonzero(~np.isnan(doy_all))[0]
    doy_valid = doy_all[valid_idx].astype(int)

    day_mean_icefrac_valid = day_mean_icefrac.isel(Time=valid_idx)
    times_valid = pd.DatetimeIndex(day_mean_icefrac["Time"].values)[valid_idx]

    doy_coord = xr.DataArray(doy_valid, dims="Time", coords={"Time": day_mean_icefrac_valid["Time"]})
    year_coord = xr.DataArray(times_valid.year, dims="Time", coords={"Time": day_mean_icefrac_valid["Time"]})
    spring_half = xr.DataArray(times_valid.month < 7, dims="Time", coords={"Time": day_mean_icefrac_valid["Time"]})

    ice_flag = (day_mean_icefrac_valid >= ICE_COVER_THRESHOLD) & water_mask

    spring_doy_where = xr.where(ice_flag & spring_half, doy_coord, np.nan).assign_coords(year=("Time", year_coord.values))
    fall_doy_where = xr.where(ice_flag & ~spring_half, doy_coord, np.nan).assign_coords(year=("Time", year_coord.values))
    ice_flag_yr = ice_flag.astype("float32").assign_coords(year=("Time", year_coord.values))

    # "last" ice day each year = latest Jan-Jun occurrence (no default; NaN if none)
    last_ice_day_peryear = spring_doy_where.groupby("year").max(dim="Time", skipna=True).reindex(year=years_range)
    # earliest Jan-Jun occurrence each year -- used only to backfill the *previous*
    # year's first_ice_day when that year's own Jul-Dec search found nothing
    # (a winter ice season that started after Jan 1)
    spring_min_doy_peryear = spring_doy_where.groupby("year").min(dim="Time", skipna=True).reindex(year=years_range)
    # "first" ice day each year = earliest Jul-Dec occurrence (no default; NaN if none)
    first_ice_day_peryear = fall_doy_where.groupby("year").min(dim="Time", skipna=True).reindex(year=years_range)

    # cross-year backfill: shift(year=-1) brings year Y+1's spring-min value to
    # year Y's slot, filling any year whose own fall search was empty
    spring_min_next_year = spring_min_doy_peryear.shift(year=-1)
    first_ice_day_peryear = xr.where(
        np.isnan(first_ice_day_peryear) & ~np.isnan(spring_min_next_year),
        spring_min_next_year,
        first_ice_day_peryear,
    )

    ice_days_peryear = ice_flag_yr.groupby("year").sum(dim="Time", skipna=True).reindex(year=years_range)

    lake_ice_days = ice_days_peryear.mean(dim="year", skipna=True)
    first_ice_day = first_ice_day_peryear.mean(dim="year", skipna=True)
    last_ice_day = last_ice_day_peryear.mean(dim="year", skipna=True)

    lake_ice_days = lake_ice_days.where(water_mask)
    first_ice_day = first_ice_day.where(water_mask)
    last_ice_day = last_ice_day.where(water_mask)

    # ── Basin ice-cover percent climatology by day of year ─────────────────
    print("Computing basin ice-cover percent climatology...")
    group_pct_by_doy = []
    for name in LAKE_GROUPS:
        mask = xr.DataArray(group_masks[name], dims=("south_north", "west_east"))
        group_ts = day_mean_icefrac_valid.where(mask).mean(dim=["south_north", "west_east"], skipna=True) * 100.0
        group_ts = group_ts.assign_coords(clim_doy=("Time", doy_valid))
        group_clim = group_ts.groupby("clim_doy").mean(dim="Time", skipna=True).reindex(clim_doy=np.arange(1, N_DOY + 1))
        group_pct_by_doy.append(group_clim)
    mean_ice_cover_percent = xr.concat(group_pct_by_doy, dim=pd.Index(LAKE_GROUPS, name="lake_group"))
    mean_ice_cover_percent = mean_ice_cover_percent.transpose("clim_doy", "lake_group").rename({"clim_doy": "doy"})

    # --- write output ---
    print("Writing output...")
    ds_out = xr.Dataset(
        data_vars={
            "ICEFRAC_mean": (
                ("month", "south_north", "west_east"),
                icefrac_mean.data,
                {
                    "units": "-",
                    "long_name": "mean lake surface (top layer) ice fraction by month",
                    "_FillValue": np.nan,
                    "coordinates": "XLAT XLONG",
                },
            ),
            "TLAKE_sfc_mean": (
                ("month", "south_north", "west_east"),
                tlake_mean.data,
                {
                    "units": "K",
                    "long_name": "mean top-layer T_LAKE3D by month",
                    "_FillValue": np.nan,
                    "coordinates": "XLAT XLONG",
                },
            ),
            "TLAKE_sfc_p95": (
                ("month", "south_north", "west_east"),
                tlake_p95.data,
                {
                    "units": "K",
                    "long_name": "95th percentile of top-layer T_LAKE3D by month",
                    "_FillValue": np.nan,
                    "coordinates": "XLAT XLONG",
                },
            ),
            "TLAKE_sfc_p99": (
                ("month", "south_north", "west_east"),
                tlake_p99.data,
                {
                    "units": "K",
                    "long_name": "99th percentile of top-layer T_LAKE3D by month",
                    "_FillValue": np.nan,
                    "coordinates": "XLAT XLONG",
                },
            ),
            "lake_ice_days": (
                ("stat", "south_north", "west_east"),
                lake_ice_days.data[np.newaxis, :, :],
                {
                    "units": "days",
                    "long_name": f"mean number of days per year with ice fraction >= {ICE_COVER_THRESHOLD}",
                    "_FillValue": np.nan,
                    "coordinates": "XLAT XLONG",
                },
            ),
            "first_ice_day": (
                ("stat", "south_north", "west_east"),
                first_ice_day.data[np.newaxis, :, :],
                {
                    "units": "days",
                    "long_name": "mean first day-of-year (365-day calendar) of ice cover in the cold season",
                    "_FillValue": np.nan,
                    "coordinates": "XLAT XLONG",
                },
            ),
            "last_ice_day": (
                ("stat", "south_north", "west_east"),
                last_ice_day.data[np.newaxis, :, :],
                {
                    "units": "days",
                    "long_name": "mean last day-of-year (365-day calendar) of ice cover in the cold season",
                    "_FillValue": np.nan,
                    "coordinates": "XLAT XLONG",
                },
            ),
            "lake_basin": (
                ("stat", "south_north", "west_east"),
                lake_basin[np.newaxis, :, :],
                {
                    "long_name": "Great Lakes basin each gridpoint drains to",
                    "valid_values": ", ".join(BASIN_PRIORITY + [NONE_LABEL]),
                    "coordinates": "XLAT XLONG",
                },
            ),
            "lake_basin_flag": (
                ("stat", "south_north", "west_east"),
                lake_basin_flag[np.newaxis, :, :],
                {
                    "long_name": "1 if gridpoint falls within any of the 5 Great Lakes basins, else 0",
                    "coordinates": "XLAT XLONG",
                },
            ),
            "mean_ice_cover_percent": (
                ("doy", "lake_group"),
                mean_ice_cover_percent.data,
                {
                    "units": "%",
                    "long_name": "mean lake surface ice cover by day of year",
                    "order": ", ".join(LAKE_GROUPS),
                    "_FillValue": np.nan,
                },
            ),
        },
        coords={
            "south_north": np.arange(ny),
            "west_east": np.arange(nx),
            "month": MONTHS,
            "stat": [0],
            "doy": np.arange(1, N_DOY + 1),
            "lake_group": LAKE_GROUPS,
            "XLAT": (("south_north", "west_east"), lat, {"units": "degree_north"}),
            "XLONG": (("south_north", "west_east"), lon, {"units": "degree_east"}),
        },
        attrs={
            "description": "Lake surface postprocessed statistics from WRF wrfout files",
            "scenario": args.scenario,
            "start_year": start_year,
            "end_year": end_year - 1,
            "ice_cover_threshold": ICE_COVER_THRESHOLD,
            "source_n_files": len(files),
        },
    )

    outfile = f"output_new/stats_lake_surface_{args.scenario}_{start_year}-{end_year}.nc"
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    print(f"Writing outfile (triggers the lazy computation): {outfile}")
    ds_out.to_netcdf(outfile)
    print("Done.")


if __name__ == "__main__":
    main()
