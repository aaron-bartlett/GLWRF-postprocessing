import argparse
import os
import sys
import glob
import numpy as np
import xarray as xr
from datetime import datetime
from wrf import getvar
from netCDF4 import Dataset

# --------- NOTE ---------- #
# --- Only works for vars with ONLY --- #
# --- north_south, west_east, time dimensions --- #
#varnames = ['T2','D2','Q2','SLP','U10','V10','CLDFRA','PW','GRDFLX','HFX','LH']
varnames = ['T2','D2','Q2','SLP','U10','V10','PW','GRDFLX','HFX','LH']
#varnames = ['T2','U10','V10']

def get_wrfvar(ds, varname):

    if varname in ['D2', 'SLP','PW']:
        if varname == "D2":
            arr = getvar(ds, "td2").values  # dewpoint
        elif varname == "SLP":
            arr = getvar(ds, "slp").values  # sea-level pressure
        elif varname == "PW":
            arr = getvar(ds, "pw").values  # total column precipitable water
    else:
        if varname not in ds.variables:
            raise KeyError(f"Variable {varname} not found in WRF file")

        arr = np.asarray(ds[varname][:]).astype(float)

        if np.ma.isMaskedArray(arr):
            arr = arr.filled(np.nan)

        #arr = ds[varname].values

    if arr.ndim == 2:
        arr = arr[np.newaxis, :, :]

    return arr

def check_year(value):
    ivalue = int(value)
    if ivalue < 2000 or ivalue > 2099:
        raise argparse.ArgumentTypeError(f"{value} is an invalid year. Must be between 2000 and 2099.")
    return ivalue

def main():

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

    args = parser.parse_args()
    home_dir = args.path
    start_year = args.year_start
    end_year = start_year + 10

    if os.path.isdir(home_dir):
        print(f"Accessing directory: {home_dir}")
    else:
        print(f"Error: {home_dir} is not a valid directory.")



    months = range(1,13)

    files = sorted(glob.glob(os.path.join(home_dir, "wrfout_d01_20*")))

    sample = xr.open_dataset(files[0])
    s_n = sample.coords["south_north"]
    w_e = sample.coords["west_east"]
    
    data_vars = {}

    for v in varnames:
        data_vars[f"mean_{v}"] = (
            ("south_north", "west_east", "month"),
            np.full((s_n.shape[0], w_e.shape[0], len(months)), np.nan, dtype=float),
            {"_FillValue": np.nan}
        )
        data_vars[f"p95_{v}"] = (
            ("south_north", "west_east", "month"),
            np.full((s_n.shape[0], w_e.shape[0], len(months)), np.nan, dtype=float),
            {"_FillValue": np.nan}
        )
        data_vars[f"p99_{v}"] = (
            ("south_north", "west_east", "month"),
            np.full((s_n.shape[0], w_e.shape[0], len(months)), np.nan, dtype=float),
            {"_FillValue": np.nan}
        )

    ds = xr.Dataset(
        coords={
            "south_north": s_n,
            "west_east": w_e,
            "month": months,
        },
        data_vars=data_vars,
        attrs={
            "description": "Monthly mean/p95/p99 surface meteorological statistics from WRF wrfout files",
            "scenario": args.scenario,
            "start_year": start_year,
            "end_year": end_year - 1,
        },
    )

    ds = ds.assign(
        XLAT=sample["XLAT"].isel(Time=0),
        XLONG=sample["XLONG"].isel(Time=0),
    )
    sample.close()

    
    # Partition files by calendar month up front (filename parsing only, no
    # file opens). This replaces the old "rescan the whole file list once per
    # month" loop -- each file is now opened exactly once in total instead of
    # up to 12 times -- while still processing (and discarding) one month's
    # data at a time, so peak memory stays bounded to a single month's worth
    # of grids rather than accumulating across all 12 months.
    files_by_month = {m: [] for m in months}
    for f in files:
        timestamp = f.split("_d01_")[1]
        file_dt = datetime.strptime(timestamp, "%Y-%m-%d_%H:%M:%S")
        if file_dt.year >= start_year:
            files_by_month[file_dt.month].append(f)

    for i, mo in enumerate(months):

        print(f"Processing month {mo:02d}")

        monthly_values = {v: [] for v in varnames}

        for f in files_by_month[mo]:
            dsfile = Dataset(f)
            print(os.path.basename(f))

            for v in varnames:
                data = get_wrfvar(dsfile, v)
                monthly_values[v].append(data)

            dsfile.close()

        for v in varnames:
            # Skip if no data
            if len(monthly_values[v]) == 0:
                print(f"No data for {v} for month {mo:02d}")
                continue

            # Stack all values for the month
            clean_list = [a.filled(np.nan)
                    if np.ma.isMaskedArray(a) else a
                    for a in monthly_values[v]]
            allvals = np.concatenate(clean_list, axis=0)

            # Compute statistics
            mean_val = np.nanmean(allvals, axis=0)
            p95_val = np.nanpercentile(allvals, 95, axis=0)
            p99_val = np.nanpercentile(allvals, 99, axis=0)

            # Store in dataset
            ds[f"mean_{v}"].isel(month=i)[:, :] = mean_val
            ds[f"p95_{v}"].isel(month=i)[:, :] = p95_val
            ds[f"p99_{v}"].isel(month=i)[:, :] = p99_val

        # Free this month's raw grids before moving on to the next month.
        del monthly_values

    outfile = f"output_new/stats_monthly_mean_{args.scenario}_{start_year}-{start_year+10}.nc"
    print(f"Writing outfile: {outfile}")
    ds.to_netcdf(outfile)


if __name__ == "__main__":
    main()
