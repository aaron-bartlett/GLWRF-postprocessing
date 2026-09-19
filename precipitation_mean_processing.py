import argparse
import os
import glob
import numpy as np
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


def get_accum_precip(directory, yr, mo, dy):
    """Return RAINC+RAINNC at Time=0 (0th hour) for the given date."""
    path = f"{directory}/wrfout_d01_{yr}-{mo:02}-{dy:02}_00:00:00"
    ds = xr.open_dataset(path)
    bucket_mm = ds.attrs.get("BUCKET_MM", 100.0)
    rainc = ds["RAINC"].isel(Time=0).values
    rainnc = ds["RAINNC"].isel(Time=0).values
    i_rainc = ds["I_RAINC"].isel(Time=0).values
    i_rainnc = ds["I_RAINNC"].isel(Time=0).values

    total_precip = (rainc + i_rainc * bucket_mm) + (
        rainnc + i_rainnc * bucket_mm
    )

    ds.close()
    return total_precip


def build_daily_precip(directory, start_year, end_year, verbose):
    """
    Compute daily total precipitation for each day in [start_year-01-01, end_year-01-01).

    Reads accumulated RAINC+RAINNC at the 0th hour of each day (including the
    first day of end_year to close the final interval), then differences
    consecutive days to get daily totals.

    Returns:
        dates : list of datetime, length n_days
        daily : ndarray shape (n_days, ny, nx), mm
    """
    print("Building daily precipitation array...")


    all_dates = []
    d = datetime(start_year-1, 12, 31)
    #stop = datetime(start_year, 1, 31)
    stop = datetime(end_year, 1, 1)
    while d < stop:
        all_dates.append(d)
        d += timedelta(days=1)

    accums = []
    for d in all_dates:
        if verbose:
            print(f"Reading {d.strftime('%Y-%m-%d')}")
        accums.append(get_accum_precip(directory, d.year, d.month, d.day))

    dates = np.array(all_dates[1:])
    daily = np.stack(
        [accums[i + 1] - accums[i] for i in range(len(accums) - 1)],
        axis=0,
    )  # (n_days, ny, nx)

    neg_mask = daily < 0
    if neg_mask.any():
        print(f"WARNING: {neg_mask.sum()} negative daily-precip values clipped to 0")
        print(daily[neg_mask])
        print(dates[neg_mask.any(axis=(1, 2))])
        daily = np.maximum(daily, 0.0)

    print(f"Daily precip built: {len(dates)} days, shape {daily.shape}")
    return dates, daily


def write_monthly_stats(dates, daily, XLAT, XLONG):
    """

    Variables (dims: lat x lon x month):
        mo_daily_precip      - mean daily precip by calendar month
        mo_daily_precip_p90  - 90th-pct daily precip by calendar month
        mo_daily_precip_p95  - 95th-pct daily precip by calendar month
        mo_daily_precip_p99  - 99th-pct daily precip by calendar month
        mo_tot_precip        - mean monthly total precip by calendar month

    Returns mo_tot_precip array (ny, nx, 12) for reuse in write_annual_stats.
    """
    print("Computing monthly stats...")
    ny, nx = XLAT.shape
    months = list(range(1, 13))
    seasons = list(range(4))
    n_mo   = len(months)

    mo_idx = {m: [] for m in months}
    for i, d in enumerate(dates):
        mo_idx[d.month].append(i)

    shape = (ny, nx, n_mo)
    mo_daily_precip     = np.full(shape, np.nan)
    mo_daily_precip_p90 = np.full(shape, np.nan)
    mo_daily_precip_p95 = np.full(shape, np.nan)
    mo_daily_precip_p99 = np.full(shape, np.nan)
    mo_tot_precip       = np.full(shape, np.nan)

    for mi, m in enumerate(months):
        idxs = mo_idx[m]
        if not idxs:
            continue

        vals = daily[idxs]  # (n_days_for_month_across_all_years, ny, nx)

        mo_daily_precip[:, :, mi]     = np.nanmean(vals, axis=0)
        mo_daily_precip_p90[:, :, mi] = np.nanpercentile(vals, 90, axis=0)
        mo_daily_precip_p95[:, :, mi] = np.nanpercentile(vals, 95, axis=0)
        mo_daily_precip_p99[:, :, mi] = np.nanpercentile(vals, 99, axis=0)

        # Sum daily values within each calendar year for this month, then
        # average those per-year monthly totals across all years.
        yr_totals = {}
        for i in idxs:
            yr = dates[i].year
            if yr not in yr_totals:
                yr_totals[yr] = np.zeros((ny, nx))
            yr_totals[yr] += daily[i]

        mo_tot_precip[:, :, mi] = np.nanmean(
            np.stack(list(yr_totals.values()), axis=0), axis=0
        )

    sn_tot_precip = np.full((daily.shape[1], daily.shape[2], 4), np.nan)
    sn_daily_precip = np.full((daily.shape[1], daily.shape[2], 4), np.nan)
    for si, season_months in enumerate([[12, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]):
        sn_tot_precip[:,:,si] = np.nansum(mo_tot_precip[:, :, [m - 1 for m in season_months]], axis=2)
        sn_daily_precip[:,:,si] = np.nanmean(mo_daily_precip[:, :, [m - 1 for m in season_months]], axis=2)

    print("Complete.")
    return mo_tot_precip, {
            "mo_daily_precip": (
                ("lat", "lon", "month"), mo_daily_precip,
                {"description": "Mean daily precipitation by calendar month",
                 "units": "mm/day", "_FillValue": np.nan},
            ),
            "mo_daily_precip_p90": (
                ("lat", "lon", "month"), mo_daily_precip_p90,
                {"description": "90th percentile of daily precipitation by calendar month",
                 "units": "mm/day", "_FillValue": np.nan},
            ),
            "mo_daily_precip_p95": (
                ("lat", "lon", "month"), mo_daily_precip_p95,
                {"description": "95th percentile of daily precipitation by calendar month",
                 "units": "mm/day", "_FillValue": np.nan},
            ),
            "mo_daily_precip_p99": (
                ("lat", "lon", "month"), mo_daily_precip_p99,
                {"description": "99th percentile of daily precipitation by calendar month",
                 "units": "mm/day", "_FillValue": np.nan},
            ),
            "mo_tot_precip": (
                ("lat", "lon", "month"), mo_tot_precip,
                {"description": "Mean monthly total precipitation by calendar month",
                 "units": "mm", "_FillValue": np.nan},
            ),
            "sn_daily_precip": (
                ("lat", "lon", "season"), sn_daily_precip,
                {"description": "Mean daily precipitation by season {0: DJF, 1: MAM, 2: JJA, 3: SON}",
                 "units": "mm/day", "_FillValue": np.nan},
            ),
            "sn_tot_precip": (
                ("lat", "lon", "season"), sn_tot_precip,
                {"description": "Mean total seasonal precipitation by season {0: DJF, 1: MAM, 2: JJA, 3: SON}",
                 "units": "mm/day", "_FillValue": np.nan},
            ),
        }


def write_annual_stats(daily, mo_tot_precip, XLAT, XLONG):
    """

    Variables (dims: lat x lon):
        yr_daily_precip      - mean daily precip across all days
        yr_daily_precip_p90  - 90th-pct daily precip across all days
        yr_daily_precip_p95  - 95th-pct daily precip across all days
        yr_daily_precip_p99  - 99th-pct daily precip across all days
        yr_tot_precip        - sum of mo_tot_precip across all 12 months
    """
    print("Computing annual stats...")
    yr_daily_precip     = np.nanmean(daily, axis=0)
    yr_daily_precip_p90 = np.nanpercentile(daily, 90, axis=0)
    yr_daily_precip_p95 = np.nanpercentile(daily, 95, axis=0)
    yr_daily_precip_p99 = np.nanpercentile(daily, 99, axis=0)
    yr_tot_precip       = np.nansum(mo_tot_precip, axis=2)  # sum over 12 months -> (ny, nx)

    print("Complete.")
    
    return {
        "yr_daily_precip": (
            ("lat", "lon"), yr_daily_precip,
            {"description": "Mean daily precipitation across all days",
             "units": "mm/day", "_FillValue": np.nan},
        ),
        "yr_daily_precip_p90": (
            ("lat", "lon"), yr_daily_precip_p90,
            {"description": "90th percentile of daily precipitation across all days",
             "units": "mm/day", "_FillValue": np.nan},
        ),
        "yr_daily_precip_p95": (
            ("lat", "lon"), yr_daily_precip_p95,
            {"description": "95th percentile of daily precipitation across all days",
             "units": "mm/day", "_FillValue": np.nan},
        ),
        "yr_daily_precip_p99": (
            ("lat", "lon"), yr_daily_precip_p99,
            {"description": "99th percentile of daily precipitation across all days",
             "units": "mm/day", "_FillValue": np.nan},
        ),
        "yr_tot_precip": (
            ("lat", "lon"), yr_tot_precip,
            {"description": "Mean annual total precipitation (sum of mo_tot_precip across 12 months)",
             "units": "mm", "_FillValue": np.nan},
        ),
    }


def write_daily_stats(dates, daily, XLAT, XLONG):
    """

    Variables:
        daily_precip     - (lat, lon, day)  daily total precip
        max_precip_1day  - (lat, lon)       max single-day precip over full period
        max_precip_5day  - (lat, lon)       max 5-day total precip over full period
        max_precip_10day - (lat, lon)       max 10-day total precip over full period
    """
    print("Computing daily stats...")
    n_days, ny, nx = daily.shape

    max_1day = np.nanmax(daily, axis=0)

    # Rolling sums via cumsum for efficiency
    cumsum    = np.concatenate([np.zeros((1, ny, nx)), np.cumsum(daily, axis=0)], axis=0)
    rolling_5  = cumsum[5:]  - cumsum[:-5]   # (n_days-4, ny, nx)
    rolling_10 = cumsum[10:] - cumsum[:-10]  # (n_days-9, ny, nx)
    max_5day   = rolling_5.max(axis=0)
    max_10day  = rolling_10.max(axis=0)

    print("Complete.")
    return {
        "daily_precip": (
            ("lat", "lon", "day"),
            np.moveaxis(daily, 0, -1),  # (ny, nx, n_days)
            {"description": "Daily total precipitation (RAINC+RAINNC, differenced at 0h of each day)",
             "units": "mm", "_FillValue": np.nan},
        ),
        "max_precip_1day": (
            ("lat", "lon"), max_1day,
            {"description": "Maximum single-day total precipitation over full period",
             "units": "mm", "_FillValue": np.nan},
        ),
        "max_precip_5day": (
            ("lat", "lon"), max_5day,
            {"description": "Maximum 5-day total precipitation over full period",
             "units": "mm", "_FillValue": np.nan},
        ),
        "max_precip_10day": (
            ("lat", "lon"), max_10day,
            {"description": "Maximum 10-day total precipitation over full period",
             "units": "mm", "_FillValue": np.nan},
        ),
    }


def write_event_stats(dates, daily, XLAT, XLONG, start_year, end_year):
    """

    Per-decade variables (lat x lon) - total count over full period:
        precip_days_perdecade_1in / _2in / _3in / _4in

    Per-year mean variables (lat x lon) - mean annual day count:
        precip_days_1mm
        precip_days_1in / _2in / _3in / _4in
    """
    print("Computing event stats...")
    ny, nx = XLAT.shape

    t1mm = 1.0
    t1in = 1.0 * MM_PER_INCH
    t2in = 2.0 * MM_PER_INCH
    t3in = 3.0 * MM_PER_INCH
    t4in = 4.0 * MM_PER_INCH

    # Totals over the entire period (per decade)
    precip_days_perdecade_1in = np.sum(daily >= t1in, axis=0).astype(float)
    precip_days_perdecade_2in = np.sum(daily >= t2in, axis=0).astype(float)
    precip_days_perdecade_3in = np.sum(daily >= t3in, axis=0).astype(float)
    precip_days_perdecade_4in = np.sum(daily >= t4in, axis=0).astype(float)

    # Mean annual counts
    years  = list(range(start_year, end_year))
    yr_idx = {yr: [] for yr in years}
    for i, d in enumerate(dates):
        if d.year in yr_idx:
            yr_idx[d.year].append(i)

    def mean_annual_count(threshold):
        counts = []
        for yr in years:
            idxs = yr_idx[yr]
            if idxs:
                counts.append(np.sum(daily[idxs] >= threshold, axis=0).astype(float))
            else:
                counts.append(np.zeros((ny, nx)))
        return np.nanmean(np.stack(counts, axis=0), axis=0)

    precip_days_1mm = mean_annual_count(t1mm)
    precip_days_1in = mean_annual_count(t1in)
    precip_days_2in = mean_annual_count(t2in)
    precip_days_3in = mean_annual_count(t3in)
    precip_days_4in = mean_annual_count(t4in)

    print("Complete.")

    return {
        "precip_days_perdecade_1in": (
            ("lat", "lon"), precip_days_perdecade_1in,
            {"description": "Total days with precip >= 1 inch over full period",
             "units": "days", "_FillValue": np.nan},
        ),
        "precip_days_perdecade_2in": (
            ("lat", "lon"), precip_days_perdecade_2in,
            {"description": "Total days with precip >= 2 inches over full period",
             "units": "days", "_FillValue": np.nan},
        ),
        "precip_days_perdecade_3in": (
            ("lat", "lon"), precip_days_perdecade_3in,
            {"description": "Total days with precip >= 3 inches over full period",
             "units": "days", "_FillValue": np.nan},
        ),
        "precip_days_perdecade_4in": (
            ("lat", "lon"), precip_days_perdecade_4in,
            {"description": "Total days with precip >= 4 inches over full period",
             "units": "days", "_FillValue": np.nan},
        ),
        "precip_days_1mm": (
            ("lat", "lon"), precip_days_1mm,
            {"description": "Mean annual number of days with precip >= 1 mm",
             "units": "days/year", "_FillValue": np.nan},
        ),
        "precip_days_1in": (
            ("lat", "lon"), precip_days_1in,
            {"description": "Mean annual number of days with precip >= 1 inch",
             "units": "days/year", "_FillValue": np.nan},
        ),
        "precip_days_2in": (
            ("lat", "lon"), precip_days_2in,
            {"description": "Mean annual number of days with precip >= 2 inches",
             "units": "days/year", "_FillValue": np.nan},
        ),
        "precip_days_3in": (
            ("lat", "lon"), precip_days_3in,
            {"description": "Mean annual number of days with precip >= 3 inches",
             "units": "days/year", "_FillValue": np.nan},
        ),
        "precip_days_4in": (
            ("lat", "lon"), precip_days_4in,
            {"description": "Mean annual number of days with precip >= 4 inches",
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

    months = list(range(1, 13))
    seasons = list(range(4)) 
    XLAT, XLONG = get_grid_coords(home_dir)
    print(f"Grid shape: {XLAT.shape}")

    dates, daily = build_daily_precip(home_dir, start_year, end_year, args.verbose)
    print(f"Daily precip: {len(dates)} days {dates[0].date()} - {dates[-1].date()}, shape: {daily.shape}")
   
    full_data_vars = {}

    print("Computing Monthly/Season Precipitation Statistics...", end="")
    
    mo_tot_precip, monthly_data_vars = write_monthly_stats(dates, daily, XLAT, XLONG)
    full_data_vars.update(monthly_data_vars)
   
    print("Computing Annual Precipitation Statistics...", end="")
    
    annual_data_vars = write_annual_stats(daily, mo_tot_precip, XLAT, XLONG)
    full_data_vars.update(annual_data_vars)

    print("Computing Daily Precipitation Statistics...", end="")
    
    daily_data_vars = write_daily_stats(dates, daily, XLAT, XLONG)
    full_data_vars.update(daily_data_vars)
    
    print("Computing Precipitation Event Statistics...", end="")
    
    event_data_vars = write_event_stats(dates, daily, XLAT, XLONG, start_year, end_year)
    full_data_vars.update(event_data_vars)
    
    #events_outfile = f"events_precipitation_{args.scenario}_{start_year}-{end_year}.nc"
    #events_ds.to_netcdf(events_outfile)
    
    ds = xr.Dataset(
        coords={
            "lat": (("lat", "lon"), XLAT.values),
            "lon": (("lat", "lon"), XLONG.values),
            "month": months,
            "season": seasons,
        },
        data_vars=full_data_vars,
        attrs={
            "description": "Precipitation postprocessed statistics from WRF wrfout files",
            "scenario": args.scenario,
            "start_year": start_year,
            "end_year": end_year - 1,
        },
    )
    
    outfile = f"output_new/stats_precipitation_{args.scenario}_{start_year}-{end_year}.nc"
    ds.to_netcdf(outfile)

    print(f"Precipitation NetCDF Written to {outfile}")

if __name__ == "__main__":
    main()
