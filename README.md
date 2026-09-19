# GLWRF Postprocessing Scripts

Postprocessing pipeline for the Great Lakes WRF (GLWRF) regional climate
simulations. Each script reads raw WRF output (`wrfout_d01_*` / `wrfxtrm_d01_*`
files) for a ~10-year window of a given scenario and writes a single
compact NetCDF file of derived statistics (climatologies, extremes, degree
days, etc.) to `output_new/`, so downstream analysis never has to re-touch
the raw decade-scale model output.

For the full list of output variables — description, units, WRF source
variable, and data shape — see **[VARIABLES.md](VARIABLES.md)**.

For notes on memory/runtime characteristics and the SLURM parallelization
strategy, see **[PARALLELIZATION.md](PARALLELIZATION.md)**.

## Pipeline overview

```
wrfout_d01_YYYY-MM-DD_HH:MM:SS      (3-hourly/6-hourly raw WRF output)
wrfxtrm_d01_YYYY-MM-DD_HH:MM:SS     (daily T2MAX/T2MIN extrema output)
        │
        ▼
  processing script (this repo)
        │
        ▼
output_new/stats_<name>_<scenario>_<start>-<end>.nc
```

All six statistics scripts share the same CLI shape and decade-window
convention:

| Flag | Meaning | Default |
|---|---|---|
| `-p`, `--path` | Directory of `wrfout`/`wrfxtrm` files | `../WRF/run` |
| `-y`, `--year_start` | First year of the processing window (2000–2099) | `2005` |
| `-s`, `--scenario` | Scenario label written into output attrs/filename | `historical` |
| `-v`, `--verbose` | Print per-file progress | off |

The processing window is `[year_start, year_start + 10)` for every script
except `fvcom_forcing_processing.py`, which uses `[year_start, year_start + 12)`.

Grid dimensions are the native WRF domain — some scripts label them
`south_north`/`west_east`, others `lat`/`lon` (both are the same projected
model grid, not a geographic lat/lon vector); `VARIABLES.md` uses `X,Y` for
both. Sample domain size in this repo's test outputs is 145 × 174.

Five scenarios are processed in production (`HH`=historical,
`M2`/`M5`=mid-century SSP2-4.5/SSP5-8.5, `L2`/`L5`=late-century
SSP2-4.5/SSP5-8.5), each with its own start year and WRF run directory —
see `postproc_common.sh`.

## Scripts

### `extreme_processing.py`
Computes extreme/derived 2 m temperature statistics from **daily**
`wrfxtrm_d01_*` files (`T2MAX`/`T2MIN`). Builds a dask-backed (lazy)
daily max/min/mean series over the full decade, then derives:
day-of-year and monthly climatologies, absolute max/min, JJA mean,
growing-season length and freeze dates, heating/cooling/growing/freezing
degree days, and mean-annual counts of hot/cold day and night thresholds.
T2MAX/T2MIN values outside a physically plausible 220–330 K range are
masked before use. Output: `output_new/stats_extreme_temperature_<scenario>_<start>-<end>.nc`.

### `extreme_processing_sixhour.py`
Identical statistics set to `extreme_processing.py`, but derived from
**6-hourly** `wrfout_d01_*` files by resampling instantaneous `T2` to
daily max/min/mean first. Used when native `wrfxtrm` daily-extrema files
aren't available for a run. Output:
`output_new/stats_6hextreme_temperature_<scenario>_<start>-<end>.nc`.

### `lake_surface_processing.py`
Computes lake-surface ice and temperature statistics from the CLM lake
model's surface layer (`LAKE_ICEFRAC3D`, `T_LAKE3D`, top level only),
masked to water gridpoints via `LAKEMASK`. Classifies each lake gridpoint
into one of the five Great Lakes basins (Superior, Huron, Michigan, Erie,
Ontario) using LBRM basin-outline polygons in `lake_shapefiles/`, resolving
boundary overlaps by basin priority order. Produces monthly ice-fraction
and lake-temperature (mean/p95/p99) climatologies, mean annual ice-cover
duration and first/last ice dates (ice defined as areal fraction ≥ 0.5),
the basin/flag classification grid, and a day-of-year ice-cover-percent
climatology per basin group (including a combined Michigan-Huron group,
since they are hydraulically one lake). Output:
`output_new/stats_lake_surface_<scenario>_<start>-<end>.nc`.

### `monthly_mean_processing.py`
Computes monthly mean/p95/p99 statistics for a fixed set of surface
meteorological variables (`T2`, `D2`, `Q2`, `SLP`, `U10`, `V10`, `PW`,
`GRDFLX`, `HFX`, `LH`) from 6-hourly `wrfout_d01_*` files. `D2` (dewpoint),
`SLP` (sea-level pressure), and `PW` (precipitable water) are derived
diagnostics computed via `wrf-python`'s `getvar`; the rest are read
directly. Files are partitioned by calendar month up front (filename
parsing only) so each file is opened exactly once rather than rescanned
per month. Output:
`output_new/stats_monthly_mean_<scenario>_<start>-<end>.nc`.

### `precipitation_mean_processing.py`
Builds a daily total-precipitation array (mm) from the accumulated
bucket variables `RAINC`/`RAINNC` (plus `I_RAINC`/`I_RAINNC` bucket-reset
counters, scaled by the file's `BUCKET_MM` attribute), read once per day
at the 0th hour and differenced. From that daily series it derives:
monthly and seasonal mean/percentile/total precipitation, annual
mean/percentile/total precipitation, daily precipitation itself plus
maximum 1/5/10-day totals, and both per-decade totals and mean-annual
counts of days exceeding 1 mm / 1–4 inch thresholds. Output:
`output_new/stats_precipitation_<scenario>_<start>-<end>.nc`.

### `snowfall_mean_processing.py`
Builds a daily total-snowfall array (mm) from accumulated `SNOWNC`,
differenced the same way as precipitation. Derives annual/winter
(Nov–Apr) mean and percentile daily snowfall, annual total snowfall,
mean-annual count of days with snowfall over threshold, and maximum
1/5/10-day totals. Separately, a dask-backed pass over hourly
`SNOW`/`SNOWH`/`SNOWC` computes a day-of-year snow-depth climatology
(liquid-water-equivalent and physical height), a winter-mean snow depth,
and mean annual snow-cover days. Output:
`output_new/stats_snowfall_<scenario>_<start>-<end>.nc`.

### `fvcom_forcing_processing.py`
Not a statistics script — builds a single timeseries NetCDF forcing file
for the FVCOM hydrodynamic model directly from raw `wrfout_d01_*` files,
in the "WRF grid (structured)" surface-forcing format FVCOM expects
(matched via the `source` global attribute). Streams one file at a time
(never holds more than one timestep's grids in memory) and writes wind,
2 m air temperature, relative humidity (derived from `T2`/`Q2`/`PSFC` via
the Bolton 1980 formula), surface pressure, longwave/shortwave radiation,
precipitation rate (differenced `RAINC`+`RAINNC`), and an evaporation
rate (derived from latent heat flux `LH`) at every raw timestep. Output:
`output/fvcom_forcing_<scenario>_<start>-<end>.nc`.

## Running the pipeline

`submit_all.sh` submits six independent SLURM job arrays (one per
script), each with `--array=0-4` indexing the five scenarios defined in
`postproc_common.sh`:

```
./submit_all.sh ALL <path> <scenario>
```

To rerun or retune a single script/scenario combination, submit its
`.sbatch` file directly, e.g. `sbatch submit_snowfall.sbatch`. See
[PARALLELIZATION.md](PARALLELIZATION.md) for memory-tuning knobs
(`TIME_CHUNK`, `SNOWDEPTH_TIME_CHUNK`) and the rationale behind the
per-script job split.

## Dependencies

`numpy`, `pandas`, `xarray` (+ `dask` for the `open_mfdataset`-based
scripts), `netCDF4`, `wrf-python` (`monthly_mean_processing.py` only),
and `matplotlib` (`lake_surface_processing.py`, for basin-polygon
point-in-polygon tests).
