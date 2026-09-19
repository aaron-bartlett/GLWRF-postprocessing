import argparse
import glob
import os
import numpy as np
import xarray as xr
import pandas as pd


# wrfxtrm files are already daily-resolution (one T2MAX/T2MIN record per
# calendar day), so a chunk of 200 timesteps is ~200 days per dask chunk.
# Tune this (and the SLURM --mem for this script) if you see it thrashing or
# under-using available memory.
TIME_CHUNK = 200

FREEZE_F = 32.0
MID_SUMMER_MONTH = 7

# Physically-plausible T2MAX/T2MIN range (K); wrfxtrm occasionally has
# fill/garbage values outside this and they're masked out, same as before.
T2_VALID_MIN_K = 220.0
T2_VALID_MAX_K = 330.0


def k_to_c(k):
    """Convert Kelvin to Celsius."""
    return k - 273.15
def c_to_f(c):
    """Convert Celsius to Fahrenheit."""
    return c * 9.0 / 5.0 + 32.0


def build_clim_doy(time_index):
    """
    1-indexed DOY on a 365-day calendar for every entry in a DatetimeIndex.
    Leap-year days after Feb 28 are shifted back by 1 so all years align.
    Feb 29 entries are marked NaN (excluded from climatology), matching the
    to_clim_doy() convention used elsewhere in this pipeline.
    """
    times = pd.DatetimeIndex(time_index)
    doy = times.dayofyear.to_numpy().astype(float)
    shift = (times.is_leap_year & (times.month > 2)).astype(float)
    doy = doy - shift
    doy[(times.month == 2) & (times.day == 29)] = np.nan
    return doy


def get_grid_coords(directory):
    """Extract south_north/west_east coords from the first available wrfout file."""
    files = sorted(glob.glob(os.path.join(directory, "wrfout_d01_20*")))
    ds = xr.open_dataset(files[0])
    s_n = ds.coords["south_north"]
    w_e = ds.coords["west_east"]
    ds.close()
    return s_n, w_e


def get_grid_latlon(directory):
    """Extract 2D lat/lon arrays from the first available wrfout file."""
    files = sorted(glob.glob(os.path.join(directory, "wrfout_d01_20*")))
    ds = xr.open_dataset(files[0])
    XLAT = ds["XLAT"].isel(Time=0)
    XLONG = ds["XLONG"].isel(Time=0)
    ds.close()
    return XLAT, XLONG


def _preprocess_xtrm(ds):
    """Attach a real datetime Time coordinate and keep only T2MAX/T2MIN (as float32)."""
    times = pd.to_datetime(ds["Times"].astype(str).values, format="%Y-%m-%d_%H:%M:%S")
    ds = ds.assign_coords(Time=("Time", times))
    return ds[["T2MAX", "T2MIN"]].astype("float32")


def load_daily_extrema(directory, start_year, end_year):
    """
    Lazily load T2MAX/T2MIN from all wrfxtrm_d01 files in [start_year, end_year)
    and build daily max/min/mean as dask-backed arrays. Unlike the old
    implementation, no file's data is ever materialized as a plain Python list
    across the whole decade -- xarray/dask process it chunk by chunk.

    Returns:
        tmax_c, tmin_c, tmean_c : dask-backed xr.DataArray
            dims (Time=day, south_north, west_east), Celsius
    """
    print("Discovering wrfxtrm files...")
    all_files = sorted(glob.glob(os.path.join(directory, "wrfxtrm_d01_20*")))

    def _year_of(fpath):
        return int(os.path.basename(fpath).split("_d01_")[1][:4])

    files = [f for f in all_files if start_year <= _year_of(f) < end_year]
    if not files:
        raise FileNotFoundError(
            f"No wrfxtrm files found for years [{start_year}, {end_year}) in {directory}"
        )
    print(f"Found {len(files)} wrfxtrm files for years {start_year}-{end_year - 1}")

    ds = xr.open_mfdataset(
        files,
        engine="netcdf4",
        combine="nested",
        concat_dim="Time",
        preprocess=_preprocess_xtrm,
        chunks={"Time": TIME_CHUNK},
        # netCDF4/HDF5 isn't reliably thread-safe for concurrent file opens;
        # parallel=True has caused spurious "Unknown file format" errors.
        # This only affects the (cheap) metadata-open phase, not the
        # dask-parallel computation that follows.
        parallel=False,
    )

    t2max_valid = ds["T2MAX"].where((ds["T2MAX"] >= T2_VALID_MIN_K) & (ds["T2MAX"] <= T2_VALID_MAX_K))
    t2min_valid = ds["T2MIN"].where((ds["T2MIN"] >= T2_VALID_MIN_K) & (ds["T2MIN"] <= T2_VALID_MAX_K))

    tmax_c = k_to_c(t2max_valid)
    tmin_c = k_to_c(t2min_valid)
    tmean_c = xr.where(tmax_c.notnull() & tmin_c.notnull(), (tmax_c + tmin_c) / 2, np.nan)

    print(f"Built daily series lazily: {tmax_c.sizes['Time']} days, "
          f"shape {tmax_c.shape}, chunks {tmax_c.chunks}")
    return tmax_c, tmin_c, tmean_c


def write_extreme_stats(directory, tmax_c, tmin_c, tmean_c, start_year, end_year, scenario):
    """Compute all extreme temperature statistics and write extreme_temperature_stats.nc.

    Everything below stays lazy (dask-backed) until the single ds.to_netcdf()
    call at the end, so the underlying T2MAX/T2MIN read is only ever
    materialized once, chunk by chunk, instead of holding the whole decade in
    memory up front.
    """
    print("Computing extreme temperature statistics...")
    years_range = list(range(start_year, end_year))

    # ── Feb-29-excluded view, used only for DOY climatology + freeze dates ──
    doy_all = build_clim_doy(tmax_c["Time"].values)
    valid_idx = np.nonzero(~np.isnan(doy_all))[0]
    doy_valid = doy_all[valid_idx].astype(int)

    tmax_doy = tmax_c.isel(Time=valid_idx).assign_coords(clim_doy=("Time", doy_valid))
    tmin_doy = tmin_c.isel(Time=valid_idx).assign_coords(clim_doy=("Time", doy_valid))

    # ── DOY climatology (365-day calendar) ─────────────────────────────────
    TMAX_day_mean = (
        tmax_doy.groupby("clim_doy").mean(dim="Time")
        .reindex(clim_doy=np.arange(1, 366))
        .transpose("south_north", "west_east", "clim_doy")
        .rename({"clim_doy": "doy"})
    )
    TMIN_day_mean = (
        tmin_doy.groupby("clim_doy").mean(dim="Time")
        .reindex(clim_doy=np.arange(1, 366))
        .transpose("south_north", "west_east", "clim_doy")
        .rename({"clim_doy": "doy"})
    )
    diurnal_range = (
        (tmax_doy - tmin_doy).groupby("clim_doy").mean(dim="Time")
        .reindex(clim_doy=np.arange(1, 366))
        .transpose("south_north", "west_east", "clim_doy")
        .rename({"clim_doy": "doy"})
    )

    # ── Monthly climatology (full series, Feb 29 included) ─────────────────
    TMAX_mean = (
        tmax_c.groupby("Time.month").mean(dim="Time")
        .reindex(month=np.arange(1, 13))
        .transpose("south_north", "west_east", "month")
    )
    TMIN_mean = (
        tmin_c.groupby("Time.month").mean(dim="Time")
        .reindex(month=np.arange(1, 13))
        .transpose("south_north", "west_east", "month")
    )

    # ── Absolute extremes over full period ────────────────────────────────
    TMAX_abs = tmax_c.max(dim="Time")
    TMIN_abs = tmin_c.min(dim="Time")

    # ── Summer (JJA) mean T2 ──────────────────────────────────────────────
    jja = tmean_c.sel(Time=tmean_c["Time"].dt.month.isin([6, 7, 8]))
    T2_summer_mean = jja.mean(dim="Time")

    # ── Growing season / freeze dates ─────────────────────────────────────
    # Spring freeze: last day in Jan-Jun with T2MIN < 32F (default DOY 1 if none)
    # Fall freeze  : first day in Jul-Dec with T2MIN < 32F (default DOY 365 if none)
    # Feb 29 is excluded here, matching the DOY climatology above.
    tmin_f_valid = c_to_f(tmin_c.isel(Time=valid_idx))
    times_valid = pd.DatetimeIndex(tmax_c["Time"].values)[valid_idx]

    doy_coord = xr.DataArray(doy_valid, dims="Time", coords={"Time": tmin_f_valid["Time"]})
    year_coord = xr.DataArray(times_valid.year, dims="Time", coords={"Time": tmin_f_valid["Time"]})
    spring_half = xr.DataArray(times_valid.month < MID_SUMMER_MONTH, dims="Time",
                                coords={"Time": tmin_f_valid["Time"]})

    freeze_mask = tmin_f_valid < FREEZE_F
    spring_doy_where = xr.where(freeze_mask & spring_half, doy_coord, np.nan)
    fall_doy_where = xr.where(freeze_mask & ~spring_half, doy_coord, np.nan)
    spring_doy_where = spring_doy_where.assign_coords(year=("Time", year_coord.values))
    fall_doy_where = fall_doy_where.assign_coords(year=("Time", year_coord.values))

    # fillna BEFORE reindex: a year with data-but-no-freeze gets the 1/365
    # default; a year wholly absent from the input stays NaN (and is excluded
    # by the skipna mean below), matching the original per-year loop exactly.
    spring_freeze_peryear = (
        spring_doy_where.groupby("year").max(dim="Time", skipna=True)
        .fillna(1.0).reindex(year=years_range)
    )
    fall_freeze_peryear = (
        fall_doy_where.groupby("year").min(dim="Time", skipna=True)
        .fillna(365.0).reindex(year=years_range)
    )

    growing_season = (fall_freeze_peryear - spring_freeze_peryear).clip(min=0).mean(dim="year", skipna=True)
    last_freeze_date = spring_freeze_peryear.mean(dim="year", skipna=True)
    first_freeze_date = fall_freeze_peryear.mean(dim="year", skipna=True)

    # ── Degree days / threshold counts (full series, Feb 29 included) ──────
    tmax_f = c_to_f(tmax_c)
    tmin_f = c_to_f(tmin_c)
    tmean_f = c_to_f(tmean_c)

    def mean_annual_degdays(daily_vals):
        """Average annual sum of per-day degree-day values."""
        yearly = daily_vals.groupby("Time.year").sum(dim="Time", skipna=True)
        return yearly.reindex(year=years_range).mean(dim="year", skipna=True)

    def mean_annual_count(cond):
        """Mean annual number of days satisfying a boolean condition."""
        yearly = cond.astype("float32").groupby("Time.year").sum(dim="Time", skipna=True)
        return yearly.reindex(year=years_range).mean(dim="year", skipna=True)

    growing_degdays = mean_annual_degdays((tmean_f - 50.0).clip(min=0))
    freezing_degdays = mean_annual_degdays((32.0 - tmean_f).clip(min=0))
    cooling_degdays = mean_annual_degdays((tmean_f - 65.0).clip(min=0))
    heating_degdays = mean_annual_degdays((65.0 - tmean_f).clip(min=0))

    cold_nights_0F = mean_annual_count(tmin_f < 0.0)
    cold_nights_28F = mean_annual_count(tmin_f < 28.0)
    cold_nights_32F = mean_annual_count(tmin_f < 32.0)
    cold_days_32F = mean_annual_count(tmax_f < 32.0)
    cold_days_20F = mean_annual_count(tmax_f < 20.0)

    hot_nights_60F = mean_annual_count(tmin_f > 60.0)
    hot_nights_70F = mean_annual_count(tmin_f > 70.0)
    hot_nights_75F = mean_annual_count(tmin_f > 75.0)
    hot_nights_80F = mean_annual_count(tmin_f > 80.0)
    hot_nights_85F = mean_annual_count(tmin_f > 85.0)
    hot_nights_90F = mean_annual_count(tmin_f > 90.0)

    hot_days_85F = mean_annual_count(tmax_f > 85.0)
    hot_days_86F = mean_annual_count(tmax_f > 86.0)
    hot_days_90F = mean_annual_count(tmax_f > 90.0)
    hot_days_95F = mean_annual_count(tmax_f > 95.0)
    hot_days_100F = mean_annual_count(tmax_f > 100.0)
    hot_days_105F = mean_annual_count(tmax_f > 105.0)
    hot_days_110F = mean_annual_count(tmax_f > 110.0)
    hot_days_115F = mean_annual_count(tmax_f > 115.0)

    # ── Assemble output (still lazy) ────────────────────────────────────────
    s_n = tmax_c["south_north"]
    w_e = tmax_c["west_east"]
    doys = list(range(1, 366))
    months = list(range(1, 13))

    def var(da, description, units):
        return (da.dims, da.data, {"description": description, "units": units})

    ds = xr.Dataset(
        coords={"south_north": s_n, "west_east": w_e, "doy": doys, "month": months},
        data_vars={
            "TMAX_day_mean": var(TMAX_day_mean, "Mean daily maximum 2m temperature by day of year (365-day calendar)", "C"),
            "TMIN_day_mean": var(TMIN_day_mean, "Mean daily minimum 2m temperature by day of year (365-day calendar)", "C"),
            "diurnal_range": var(diurnal_range, "Mean diurnal temperature range (T2MAX - T2MIN) by day of year", "C"),
            "TMAX_mean": var(TMAX_mean, "Mean of daily maximum 2m temperature by calendar month", "C"),
            "TMIN_mean": var(TMIN_mean, "Mean of daily minimum 2m temperature by calendar month", "C"),
            "TMAX": var(TMAX_abs, "Absolute maximum 2m temperature over full period", "C"),
            "TMIN": var(TMIN_abs, "Absolute minimum 2m temperature over full period", "C"),
            "T2_summer_mean": var(T2_summer_mean, "Mean daily mean 2m temperature during JJA", "C"),
            "growing_season": var(growing_season,
                "Mean growing season length: days between last spring freeze "
                "(last day Jan-Jun with T2MIN < 32F) and first fall freeze "
                "(first day Jul-Dec with T2MIN < 32F), averaged across years", "days"),
            "first_freeze_date": var(first_freeze_date,
                "Mean first fall freeze DOY: first day after Jul 1 each year "
                "with T2MIN < 32F, averaged across years (365-day calendar)", "DOY"),
            "last_freeze_date": var(last_freeze_date,
                "Mean last spring freeze DOY: last day Jan 1-Jun 30 each year "
                "with T2MIN < 32F, averaged across years (365-day calendar)", "DOY"),
            "growing_degdays": var(growing_degdays, "Mean annual growing degree days: sum of max(T2MEAN_F - 50, 0) per day", "degF-days"),
            "freezing_degdays": var(freezing_degdays, "Mean annual freezing degree days: sum of max(32 - T2MEAN_F, 0) per day", "degF-days"),
            "cooling_degdays": var(cooling_degdays, "Mean annual cooling degree days: sum of max(T2MEAN_F - 65, 0) per day", "degF-days"),
            "heating_degdays": var(heating_degdays, "Mean annual heating degree days: sum of max(65 - T2MEAN_F, 0) per day", "degF-days"),
            "cold_nights_0F": var(cold_nights_0F, "Mean annual number of nights with T2MIN < 0F", "days/year"),
            "cold_nights_28F": var(cold_nights_28F, "Mean annual number of nights with T2MIN < 28F", "days/year"),
            "cold_nights_32F": var(cold_nights_32F, "Mean annual number of nights with T2MIN < 32F", "days/year"),
            "cold_days_32F": var(cold_days_32F, "Mean annual number of days with T2MAX < 32F", "days/year"),
            "cold_days_20F": var(cold_days_20F, "Mean annual number of days with T2MAX < 20F", "days/year"),
            "hot_nights_60F": var(hot_nights_60F, "Mean annual number of nights with T2MIN > 60F", "days/year"),
            "hot_nights_70F": var(hot_nights_70F, "Mean annual number of nights with T2MIN > 70F", "days/year"),
            "hot_nights_75F": var(hot_nights_75F, "Mean annual number of nights with T2MIN > 75F", "days/year"),
            "hot_nights_80F": var(hot_nights_80F, "Mean annual number of nights with T2MIN > 80F", "days/year"),
            "hot_nights_85F": var(hot_nights_85F, "Mean annual number of nights with T2MIN > 85F", "days/year"),
            "hot_nights_90F": var(hot_nights_90F, "Mean annual number of nights with T2MIN > 90F", "days/year"),
            "hot_days_85F": var(hot_days_85F, "Mean annual number of days with T2MAX > 85F", "days/year"),
            "hot_days_86F": var(hot_days_86F, "Mean annual number of days with T2MAX > 86F", "days/year"),
            "hot_days_90F": var(hot_days_90F, "Mean annual number of days with T2MAX > 90F", "days/year"),
            "hot_days_95F": var(hot_days_95F, "Mean annual number of days with T2MAX > 95F", "days/year"),
            "hot_days_100F": var(hot_days_100F, "Mean annual number of days with T2MAX > 100F", "days/year"),
            "hot_days_105F": var(hot_days_105F, "Mean annual number of days with T2MAX > 105F", "days/year"),
            "hot_days_110F": var(hot_days_110F, "Mean annual number of days with T2MAX > 110F", "days/year"),
            "hot_days_115F": var(hot_days_115F, "Mean annual number of days with T2MAX > 115F", "days/year"),
        },
        attrs={
            "description": "Extreme temperature postprocessed statistics from WRF wrfxtrm files",
            "scenario": scenario,
            "start_year": start_year,
            "end_year": end_year - 1,
        },
    )

    sample_XLAT, sample_XLONG = get_grid_latlon(directory)
    ds = ds.assign(XLAT=sample_XLAT, XLONG=sample_XLONG)

    outfile = f"output_new/stats_extreme_temperature_{scenario}_{start_year}-{end_year}.nc"
    print(f"Writing (triggers the lazy computation): {outfile}")
    ds.to_netcdf(outfile)
    print(f"Written: {outfile}")


def check_year(value):
    ivalue = int(value)
    if ivalue < 2000 or ivalue > 2099:
        raise argparse.ArgumentTypeError(f"{value} is an invalid year. Must be between 2000 and 2099.")
    return ivalue

def main(start_year=2005, end_year=2015):

    parser = argparse.ArgumentParser(description="Process some files in a home directory.")

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
        exit(1)

    print("Begin Extreme Temperature Processing")

    s_n, w_e = get_grid_coords(home_dir)
    print(f"Grid shape: {s_n.shape} by {w_e.shape}")

    tmax_c, tmin_c, tmean_c = load_daily_extrema(home_dir, start_year, end_year)

    write_extreme_stats(home_dir, tmax_c, tmin_c, tmean_c, start_year, end_year, args.scenario)

    print("Done.")


if __name__ == "__main__":
    main()
