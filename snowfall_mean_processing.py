import argparse
import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta

MM_PER_INCH = 25.4


def get_grid_coords(directory):
    """Extract 2D lat/lon arrays from the first available WRF file."""
    files = sorted(glob.glob(os.path.join(directory, "wrfout_d01_20*")))
    ds = xr.open_dataset(files[0])
    XLAT  = ds.coords["XLAT"][0, :, :]   # (ny, nx)
    XLONG = ds.coords["XLONG"][0, :, :]  # (ny, nx)
    ds.close()
    return XLAT, XLONG


def get_accum_snowfall(directory, yr, mo, dy):
    """Return SNOWNC at Time=0 (0th hour) for the given date."""
    path = f"{directory}/wrfout_d01_{yr}-{mo:02}-{dy:02}_00:00:00"
    ds = xr.open_dataset(path)
    accum = ds["SNOWNC"].isel(Time=0).values  # (ny, nx)
    ds.close()
    return accum


def build_daily_snowfall(directory, start_year, end_year, verbose):
    """
    Compute daily total snowfall for each day in [start_year-01-01, end_year-01-01).

    Reads accumulated SNOWNC at the 0th hour of each day (including the
    first day of end_year to close the final interval), then differences
    consecutive days to get daily totals.

    Returns:
        dates : list of datetime, length n_days
        daily : ndarray shape (n_days, ny, nx), mm
    """
    print("Building daily snowfall array...")

    all_dates = []
    d = datetime(start_year-1, 12, 31)
    #d = datetime(start_year, 1, 1)
    #stop = datetime(start_year, 1, 31)
    stop = datetime(end_year, 1, 1)
    while d < stop:
        all_dates.append(d)
        d += timedelta(days=1)

    accums = []
    for d in all_dates:
        if verbose:
            print(f"Reading {d.strftime('%Y-%m-%d')}")
        snowfall = get_accum_snowfall(directory, d.year, d.month, d.day)
        accums.append(snowfall)

    dates = np.array(all_dates[1:])
    winter_mask = np.isin([dt.month for dt in dates], [11, 12, 1, 2, 3])
    
    daily = np.stack(
        [accums[i + 1] - accums[i] for i in range(len(accums) - 1)],
        axis=0,
    )  # (n_days, ny, nx)

    neg_mask = daily < 0
    if neg_mask.any():
        print(f"WARNING: {neg_mask.sum()} negative daily-snowfall values clipped to 0")
        print(daily[neg_mask])
        print(dates[neg_mask.any(axis=(1, 2))])
        daily = np.maximum(daily, 0.0)
    winter = daily[winter_mask]
    return dates, daily, winter


def write_annual_stats(daily, winter):
    """

    Variables (dims: lat x lon):
        yr_daily_snowfall      - mean daily snowfall across all days
        yr_daily_snowfall_p90  - 90th-pct daily snowfall across all days
        yr_daily_snowfall_p95  - 95th-pct daily snowfall across all days
        yr_daily_snowfall_p99  - 99th-pct daily snowfall across all days
        yr_tot_snowfall        - sum of mo_tot_snowfall across all 12 months
    """

    yr_daily_snowfall     = np.nanmean(daily, axis=0)
    winter_daily_snowfall = np.nanmean(winter, axis=0)    
    winter_daily_snowfall_p90 = np.nanpercentile(winter, 90, axis=0)
    winter_daily_snowfall_p95 = np.nanpercentile(winter, 95, axis=0)
    winter_daily_snowfall_p99 = np.nanpercentile(winter, 99, axis=0)
    yr_tot_snowfall       = np.nansum(daily, axis=0) / 10
    days_peryear_snowfall_1cm = np.sum(daily >= 10.0, axis=0).astype(float)

    print("Complete.")
    
    return {
        "yr_daily_snowfall": (
            ("lat", "lon"), yr_daily_snowfall,
            {"description": "Mean daily snowfall across all days",
             "units": "mm/day", "_FillValue": np.nan},
        ),
        "winter_daily_snowfall": (
            ("lat", "lon"), winter_daily_snowfall,
            {"description": "Mean daily snowfall in winter (November-April)",
             "units": "mm/day", "_FillValue": np.nan},
        ),
        "winter_daily_snowfall_p90": (
            ("lat", "lon"), winter_daily_snowfall_p90,
            {"description": "90th percentile of daily snowfall in winter (November-April)",
             "units": "mm/day", "_FillValue": np.nan},
        ),
        "winter_daily_snowfall_p95": (
            ("lat", "lon"), winter_daily_snowfall_p95,
            {"description": "95th percentile of daily snowfall in winter (November-April)", 
             "units": "mm/day", "_FillValue": np.nan},
        ),
        "winter_daily_snowfall_p99": (
            ("lat", "lon"), winter_daily_snowfall_p99,
            {"description": "99th percentile of daily snowfall in winter (November-April)",
             "units": "mm/day", "_FillValue": np.nan},
        ),
        "yr_tot_snowfall": (
            ("lat", "lon"), yr_tot_snowfall,
            {"description": "Mean annual total snowfall (sum of mo_tot_snowfall across 12 months)",
             "units": "mm", "_FillValue": np.nan},
        ),
        "yr_days_snowfall_1cm": (
            ("lat", "lon"), days_peryear_snowfall_1cm,
            {"description": "Mean annual number of days with snowfall exceeding 1cm", 
             "units": "days/year", "_FillValue": np.nan},
        ),
    }


def write_daily_stats(dates, daily):
    """

    Variables:
        daily_snowfall     - (lat, lon, day)  daily total snowfall
        max_snowfall_1day  - (lat, lon)       max single-day snowfall over full period
        max_snowfall_5day  - (lat, lon)       max 5-day total snowfall over full period
        max_snowfall_10day - (lat, lon)       max 10-day total snowfall over full period
    """
    n_days, ny, nx = daily.shape

    #daily_mean = daily.reshape(365, int(n_days/365), ny, nx).mean(axis=1)
    max_1day = np.nanmax(daily, axis=0)

    # Rolling sums via cumsum for efficiency
    cumsum    = np.concatenate([np.zeros((1, ny, nx)), np.cumsum(daily, axis=0)], axis=0)
    rolling_5  = cumsum[5:]  - cumsum[:-5]   # (n_days-4, ny, nx)
    rolling_10 = cumsum[10:] - cumsum[:-10]  # (n_days-9, ny, nx)
    max_5day   = rolling_5.max(axis=0)
    max_10day  = rolling_10.max(axis=0)

    print("Complete.")
    
    return {
        "max_snowfall_1day": (
            ("lat", "lon"), max_1day,
            {"description": "Maximum single-day total snowfall over full period",
             "units": "mm", "_FillValue": np.nan},
        ),
        "max_snowfall_5day": (
            ("lat", "lon"), max_5day,
            {"description": "Maximum 5-day total snowfall over full period",
             "units": "mm", "_FillValue": np.nan},
        ),
        "max_snowfall_10day": (
            ("lat", "lon"), max_10day,
            {"description": "Maximum 10-day total snowfall over full period",
             "units": "mm", "_FillValue": np.nan},
        ),
    }
    '''
        "daily_snowfall": (
            ("lat", "lon", "doy"),
            np.moveaxis(daily_mean, 0, -1),  # (ny, nx, n_days)
            {"description": "Daily total snowfall (SNOWNC, differenced at 0h of each day)",
             "units": "mm", "_FillValue": np.nan},
        ),
    '''


# wrfout files are 6-hourly (4 timesteps/day); a chunk of 200 timesteps is
# ~50 days per dask chunk. Tune this (and the SLURM --mem for this script) if
# you see it thrashing or under-using available memory.
SNOWDEPTH_TIME_CHUNK = 200


def _build_clim_doy_simple(time_index):
    """
    1-indexed DOY on a 365-day calendar, using the same simplified leap-year
    rule (year % 4 == 0) as the original per-file loop this replaces. Feb 29
    is marked NaN (excluded), matching the original's "Leap Days not
    included" doy=366 bucket that was always skipped downstream.
    """
    times = pd.DatetimeIndex(time_index)
    doy = times.dayofyear.to_numpy().astype(float)
    shift = ((times.year % 4 == 0) & (times.month > 2)).astype(float)
    doy = doy - shift
    doy[(times.month == 2) & (times.day == 29)] = np.nan
    return doy


def _preprocess_snowdepth(ds):
    """Attach a real datetime Time coordinate; keep only SNOW/SNOWH/SNOWC (float32)."""
    times = pd.to_datetime(ds["Times"].astype(str).values, format="%Y-%m-%d_%H:%M:%S")
    ds = ds.assign_coords(Time=("Time", times))
    return ds[["SNOW", "SNOWH", "SNOWC"]].astype("float32")


def write_snowdepth_stats(directory, start_year):
    """
    Dask-backed rewrite: the original loaded every hourly SNOW/SNOWH/SNOWC
    grid for the entire multi-year period into Python lists before reducing
    (the single largest memory user in this pipeline). Here the whole SNOW/
    SNOWH/SNOWC series stays lazy (dask-backed) and is only ever materialized
    chunk by chunk, once, when ds.to_netcdf() runs in main().
    """
    all_files = sorted(glob.glob(os.path.join(directory, 'wrfout_d01_20*')))

    def _year_of(fpath):
        return int(os.path.basename(fpath).split("_d01_")[1][:4])

    # NOTE: matches the original -- no upper year bound is applied here.
    files = [f for f in all_files if _year_of(f) >= start_year]
    if not files:
        raise FileNotFoundError(f"No wrfout files found for year >= {start_year} in {directory}")

    ds = xr.open_mfdataset(
        files,
        engine="netcdf4",
        combine="nested",
        concat_dim="Time",
        preprocess=_preprocess_snowdepth,
        chunks={"Time": SNOWDEPTH_TIME_CHUNK},
        # netCDF4/HDF5 isn't reliably thread-safe for concurrent file opens;
        # parallel=True has caused spurious "Unknown file format" errors.
        # This only affects the (cheap) metadata-open phase, not the
        # dask-parallel computation that follows.
        parallel=False,
    )

    doy_all = _build_clim_doy_simple(ds["Time"].values)
    valid_idx = np.nonzero(~np.isnan(doy_all))[0]
    doy_valid = doy_all[valid_idx].astype(int)

    def doy_climatology(da):
        da_valid = da.isel(Time=valid_idx).assign_coords(clim_doy=("Time", doy_valid))
        return (
            da_valid.groupby("clim_doy").mean(dim="Time", skipna=True)
            .reindex(clim_doy=np.arange(1, 366))
            .transpose("south_north", "west_east", "clim_doy")
            .rename({"clim_doy": "doy"})
        )

    mean_snow_365 = doy_climatology(ds["SNOW"])
    mean_snowh_365 = doy_climatology(ds["SNOWH"])

    # winter = doy 1-90 (Jan-Mar) plus doy 305-365 (Nov-Dec), same slicing as original
    mean_snow_winter = xr.concat(
        [mean_snow_365.isel(doy=slice(0, 90)), mean_snow_365.isel(doy=slice(304, 365))], dim="doy"
    ).mean(dim="doy", skipna=True)
    mean_snowh_winter = xr.concat(
        [mean_snowh_365.isel(doy=slice(0, 90)), mean_snowh_365.isel(doy=slice(304, 365))], dim="doy"
    ).mean(dim="doy", skipna=True)

    # snow-cover day count uses every timestep (Feb 29 included), matching original
    snowc_days = ds["SNOWC"].sum(dim="Time", skipna=True)
    mean_snowc_days_per_year = snowc_days / 40  # unchanged from original (not start_year/end_year-derived)

    print("Complete.")

    return {
        "snow_depth": (
            ("lat", "lon", "doy"), mean_snow_365.data,
            {"description": "Mean liquid-equivalent snow depth by gridpoint for each date",
            "units": "kg*m^-2", "_FillValue": np.nan},
        ),
        "snow_depth_height": (
            ("lat", "lon", "doy"), mean_snowh_365.data,
            {"description": "Mean snow depth height by gridpoint for each date",
            "units": "m", "_FillValue": np.nan},
        ),
        "winter_snow_depth": (
            ("lat", "lon"), mean_snow_winter.data,
            {"description": "Mean winter (November-April) liquid-equivalent snow depth by gridpoint",
             "units": "kg*m^-2", "_FillValue": np.nan},
        ),
        "winter_snow_depth_height": (
            ("lat", "lon"), mean_snowh_winter.data,
            {"description": "Mean winter (November-April) snow depth height by gridpoint",
                "units": "m", "_FillValue": np.nan},
        ),
        "snow_cover_days": (
            ("lat", "lon"), mean_snowc_days_per_year.data,
            {"description": "Mean days per year with snow cover by gridpoint",
             "units": "days/year", "_FillValue": np.nan},
        ),
    }


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

    print("Begin Precipitation Mean Processing")

    XLAT, XLONG = get_grid_coords(home_dir)
    print(f"Grid shape: {XLAT.shape}")

    dates, daily, winter = build_daily_snowfall(home_dir, start_year, end_year, args.verbose)
    print(f"Daily snowfall: {len(dates)} days {dates[0].date()} - {dates[-1].date()}, shape: {daily.shape}")
     
    full_data_vars = {}
    
    print("Computing Monthly/Season/Annual Statistics...", end="")
    
    annual_data_vars = write_annual_stats(daily, winter)
    full_data_vars.update(annual_data_vars)

    print("Computing Daily Statistics...", end="")

    daily_data_vars = write_daily_stats(dates, daily)
    full_data_vars.update(daily_data_vars)
    
    print("Computing Snow Depth Statistics...", end="")
    
    depth_data_vars = write_snowdepth_stats(home_dir, start_year)
    full_data_vars.update(depth_data_vars)
    
    ds = xr.Dataset(
        coords={
            "lat": (("lat", "lon"), XLAT.values),
            "lon": (("lat", "lon"), XLONG.values),
            "doy": range(1,366),
        },
        data_vars=full_data_vars,
        attrs={
            "description": "Snowfall/snow depth postprocessed statistics from WRF wrfout files",
            "scenario": args.scenario,
            "start_year": start_year,
            "end_year": end_year - 1,
        },
    )

    outfile = f"output_new/stats_snowfall_{args.scenario}_{start_year}-{end_year}.nc"
    ds.to_netcdf(outfile)

    print(f"Snowfall NetCDF Written to {outfile}")


if __name__ == "__main__":
    main()
