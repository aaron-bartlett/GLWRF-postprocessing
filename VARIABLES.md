# Postprocessed Variable Reference

Full listing of every variable written by each postprocessing script. `X,Y`
denotes the native WRF grid (`south_north` × `west_east`, sometimes labeled
`lat`/`lon`); other dimension names (`month`, `doy`, `season`, `day`, `time`,
`stat`, `lake_group`) are 1-D index/coordinate axes. See
[README.md](README.md) for what each script does and how to run it.

## `extreme_processing.py` / `extreme_processing_sixhour.py`

Both scripts write the identical set of variables to
`stats_extreme_temperature_*.nc` and `stats_6hextreme_temperature_*.nc`
respectively. The only difference is the WRF source: `extreme_processing.py`
reads daily `T2MAX`/`T2MIN` from `wrfxtrm` files directly, while
`extreme_processing_sixhour.py` derives the same daily max/min by
resampling instantaneous `T2` from 6-hourly `wrfout` files.

| Variable | Description | Units | WRF Source Variable | Shape |
|---|---|---|---|---|
| `TMAX_day_mean` | Mean daily maximum 2 m temperature by day of year (365-day calendar) | C | `T2MAX` (or `T2` resampled) | X,Y,doy |
| `TMIN_day_mean` | Mean daily minimum 2 m temperature by day of year (365-day calendar) | C | `T2MIN` (or `T2` resampled) | X,Y,doy |
| `diurnal_range` | Mean diurnal temperature range (TMAX − TMIN) by day of year | C | `T2MAX`, `T2MIN` | X,Y,doy |
| `TMAX_mean` | Mean of daily maximum 2 m temperature by calendar month | C | `T2MAX` | X,Y,month |
| `TMIN_mean` | Mean of daily minimum 2 m temperature by calendar month | C | `T2MIN` | X,Y,month |
| `TMAX` | Absolute maximum 2 m temperature over full period | C | `T2MAX` | X,Y |
| `TMIN` | Absolute minimum 2 m temperature over full period | C | `T2MIN` | X,Y |
| `T2_summer_mean` | Mean daily-mean 2 m temperature during JJA | C | `T2MAX`, `T2MIN` (averaged) | X,Y |
| `growing_season` | Mean growing season length: days between last spring freeze (last day Jan–Jun with TMIN < 32°F) and first fall freeze (first day Jul–Dec with TMIN < 32°F), averaged across years | days | `T2MIN` | X,Y |
| `first_freeze_date` | Mean first fall freeze day-of-year, averaged across years (365-day calendar) | DOY | `T2MIN` | X,Y |
| `last_freeze_date` | Mean last spring freeze day-of-year, averaged across years (365-day calendar) | DOY | `T2MIN` | X,Y |
| `growing_degdays` | Mean annual growing degree days: Σ max(T2MEAN_F − 50, 0) per day | degF-days | `T2MAX`, `T2MIN` | X,Y |
| `freezing_degdays` | Mean annual freezing degree days: Σ max(32 − T2MEAN_F, 0) per day | degF-days | `T2MAX`, `T2MIN` | X,Y |
| `cooling_degdays` | Mean annual cooling degree days: Σ max(T2MEAN_F − 65, 0) per day | degF-days | `T2MAX`, `T2MIN` | X,Y |
| `heating_degdays` | Mean annual heating degree days: Σ max(65 − T2MEAN_F, 0) per day | degF-days | `T2MAX`, `T2MIN` | X,Y |
| `cold_nights_0F` | Mean annual number of nights with TMIN < 0°F | days/year | `T2MIN` | X,Y |
| `cold_nights_28F` | Mean annual number of nights with TMIN < 28°F | days/year | `T2MIN` | X,Y |
| `cold_nights_32F` | Mean annual number of nights with TMIN < 32°F | days/year | `T2MIN` | X,Y |
| `cold_days_32F` | Mean annual number of days with TMAX < 32°F | days/year | `T2MAX` | X,Y |
| `cold_days_20F` | Mean annual number of days with TMAX < 20°F | days/year | `T2MAX` | X,Y |
| `hot_nights_60F` | Mean annual number of nights with TMIN > 60°F | days/year | `T2MIN` | X,Y |
| `hot_nights_70F` | Mean annual number of nights with TMIN > 70°F | days/year | `T2MIN` | X,Y |
| `hot_nights_75F` | Mean annual number of nights with TMIN > 75°F | days/year | `T2MIN` | X,Y |
| `hot_nights_80F` | Mean annual number of nights with TMIN > 80°F | days/year | `T2MIN` | X,Y |
| `hot_nights_85F` | Mean annual number of nights with TMIN > 85°F | days/year | `T2MIN` | X,Y |
| `hot_nights_90F` | Mean annual number of nights with TMIN > 90°F | days/year | `T2MIN` | X,Y |
| `hot_days_85F` | Mean annual number of days with TMAX > 85°F | days/year | `T2MAX` | X,Y |
| `hot_days_86F` | Mean annual number of days with TMAX > 86°F | days/year | `T2MAX` | X,Y |
| `hot_days_90F` | Mean annual number of days with TMAX > 90°F | days/year | `T2MAX` | X,Y |
| `hot_days_95F` | Mean annual number of days with TMAX > 95°F | days/year | `T2MAX` | X,Y |
| `hot_days_100F` | Mean annual number of days with TMAX > 100°F | days/year | `T2MAX` | X,Y |
| `hot_days_105F` | Mean annual number of days with TMAX > 105°F | days/year | `T2MAX` | X,Y |
| `hot_days_110F` | Mean annual number of days with TMAX > 110°F | days/year | `T2MAX` | X,Y |
| `hot_days_115F` | Mean annual number of days with TMAX > 115°F | days/year | `T2MAX` | X,Y |
| `XLAT` | Gridpoint latitude | degree_north | `XLAT` | X,Y |
| `XLONG` | Gridpoint longitude (west negative) | degree_east | `XLONG` | X,Y |

## `lake_surface_processing.py`

Output: `stats_lake_surface_*.nc`.

| Variable | Description | Units | WRF Source Variable | Shape |
|---|---|---|---|---|
| `ICEFRAC_mean` | Mean lake surface (top layer) ice fraction by month | - (fraction) | `LAKE_ICEFRAC3D` (surface level), masked by `LAKEMASK` | month,X,Y |
| `TLAKE_sfc_mean` | Mean top-layer lake temperature by month | K | `T_LAKE3D` (surface level) | month,X,Y |
| `TLAKE_sfc_p95` | 95th percentile of top-layer lake temperature by month | K | `T_LAKE3D` (surface level) | month,X,Y |
| `TLAKE_sfc_p99` | 99th percentile of top-layer lake temperature by month | K | `T_LAKE3D` (surface level) | month,X,Y |
| `lake_ice_days` | Mean number of days per year with ice fraction ≥ 0.5 | days | `LAKE_ICEFRAC3D` (derived, daily-resampled) | stat,X,Y |
| `first_ice_day` | Mean first day-of-year of ice cover in the cold season (365-day calendar) | days (DOY) | `LAKE_ICEFRAC3D` (derived) | stat,X,Y |
| `last_ice_day` | Mean last day-of-year of ice cover in the cold season (365-day calendar) | days (DOY) | `LAKE_ICEFRAC3D` (derived) | stat,X,Y |
| `lake_basin` | Great Lakes basin each gridpoint drains to (`SUPERIOR`/`HURON`/`MICHIGAN`/`ERIE`/`ONTARIO`/`NONE`) | categorical | derived from `XLAT`/`XLONG` + LBRM basin outline shapefiles | stat,X,Y |
| `lake_basin_flag` | 1 if gridpoint falls within any of the 5 Great Lakes basins, else 0 | flag (0/1) | derived from `XLAT`/`XLONG` + LBRM basin outline shapefiles | stat,X,Y |
| `mean_ice_cover_percent` | Mean lake surface ice cover by day of year, per basin group (`FULL_BASIN`, 5 named lakes, and combined `MICHIGAN-HURON`) | % | `LAKE_ICEFRAC3D` (derived, basin-averaged) | doy,lake_group |
| `XLAT` | Gridpoint latitude | degree_north | `XLAT` | X,Y |
| `XLONG` | Gridpoint longitude | degree_east | `XLONG` | X,Y |

## `monthly_mean_processing.py`

Output: `stats_monthly_mean_*.nc`. For each variable `v` in
`{T2, D2, Q2, SLP, U10, V10, PW, GRDFLX, HFX, LH}`, three statistics are
written: `mean_v`, `p95_v`, `p99_v` (all dims X,Y,month).

| Variable prefix | Description | Units | WRF Source Variable | Shape |
|---|---|---|---|---|
| `mean_T2` / `p95_T2` / `p99_T2` | Monthly mean / 95th / 99th percentile of 2 m temperature | K | `T2` | X,Y,month |
| `mean_D2` / `p95_D2` / `p99_D2` | Monthly mean / 95th / 99th percentile of 2 m dewpoint temperature | deg C | `td2` diagnostic (`wrf-python getvar`, from `T2`/`Q2`/`PSFC`) | X,Y,month |
| `mean_Q2` / `p95_Q2` / `p99_Q2` | Monthly mean / 95th / 99th percentile of 2 m water vapor mixing ratio | kg/kg | `Q2` | X,Y,month |
| `mean_SLP` / `p95_SLP` / `p99_SLP` | Monthly mean / 95th / 99th percentile of sea-level pressure | hPa | `slp` diagnostic (`wrf-python getvar`, from pressure/temperature/geopotential columns) | X,Y,month |
| `mean_U10` / `p95_U10` / `p99_U10` | Monthly mean / 95th / 99th percentile of 10 m U-wind | m/s | `U10` | X,Y,month |
| `mean_V10` / `p95_V10` / `p99_V10` | Monthly mean / 95th / 99th percentile of 10 m V-wind | m/s | `V10` | X,Y,month |
| `mean_PW` / `p95_PW` / `p99_PW` | Monthly mean / 95th / 99th percentile of total column precipitable water | kg/m² | `pw` diagnostic (`wrf-python getvar`, from full-column moisture profile) | X,Y,month |
| `mean_GRDFLX` / `p95_GRDFLX` / `p99_GRDFLX` | Monthly mean / 95th / 99th percentile of ground heat flux | W/m² | `GRDFLX` | X,Y,month |
| `mean_HFX` / `p95_HFX` / `p99_HFX` | Monthly mean / 95th / 99th percentile of upward surface sensible heat flux | W/m² | `HFX` | X,Y,month |
| `mean_LH` / `p95_LH` / `p99_LH` | Monthly mean / 95th / 99th percentile of surface latent heat flux | W/m² | `LH` | X,Y,month |
| `XLAT` | Gridpoint latitude | degree_north | `XLAT` | X,Y |
| `XLONG` | Gridpoint longitude | degree_east | `XLONG` | X,Y |

## `precipitation_mean_processing.py`

Output: `stats_precipitation_*.nc`. Daily total precipitation is built by
differencing accumulated `RAINC`+`RAINNC` (each corrected for bucket resets
via `I_RAINC`/`I_RAINNC` × the file's `BUCKET_MM` attribute) at the 0th
hour of consecutive days.

| Variable | Description | Units | WRF Source Variable | Shape |
|---|---|---|---|---|
| `mo_daily_precip` | Mean daily precipitation by calendar month | mm/day | `RAINC`, `RAINNC`, `I_RAINC`, `I_RAINNC` (derived) | X,Y,month |
| `mo_daily_precip_p90` | 90th percentile of daily precipitation by calendar month | mm/day | derived (as above) | X,Y,month |
| `mo_daily_precip_p95` | 95th percentile of daily precipitation by calendar month | mm/day | derived (as above) | X,Y,month |
| `mo_daily_precip_p99` | 99th percentile of daily precipitation by calendar month | mm/day | derived (as above) | X,Y,month |
| `mo_tot_precip` | Mean monthly total precipitation by calendar month | mm | derived (as above) | X,Y,month |
| `sn_daily_precip` | Mean daily precipitation by season (0=DJF, 1=MAM, 2=JJA, 3=SON) | mm/day | derived (as above) | X,Y,season |
| `sn_tot_precip` | Mean total seasonal precipitation by season (0=DJF, 1=MAM, 2=JJA, 3=SON) | mm/day *(as labeled in file; a total, despite the per-day unit string)* | derived (as above) | X,Y,season |
| `yr_daily_precip` | Mean daily precipitation across all days | mm/day | derived (as above) | X,Y |
| `yr_daily_precip_p90` | 90th percentile of daily precipitation across all days | mm/day | derived (as above) | X,Y |
| `yr_daily_precip_p95` | 95th percentile of daily precipitation across all days | mm/day | derived (as above) | X,Y |
| `yr_daily_precip_p99` | 99th percentile of daily precipitation across all days | mm/day | derived (as above) | X,Y |
| `yr_tot_precip` | Mean annual total precipitation (sum of `mo_tot_precip` across 12 months) | mm | derived (as above) | X,Y |
| `daily_precip` | Daily total precipitation, one value per calendar day in the period | mm | derived (as above) | X,Y,day |
| `max_precip_1day` | Maximum single-day total precipitation over full period | mm | derived (as above) | X,Y |
| `max_precip_5day` | Maximum 5-day total precipitation over full period | mm | derived (as above) | X,Y |
| `max_precip_10day` | Maximum 10-day total precipitation over full period | mm | derived (as above) | X,Y |
| `precip_days_perdecade_1in` | Total days with precip ≥ 1 inch over full period | days | derived (as above) | X,Y |
| `precip_days_perdecade_2in` | Total days with precip ≥ 2 inches over full period | days | derived (as above) | X,Y |
| `precip_days_perdecade_3in` | Total days with precip ≥ 3 inches over full period | days | derived (as above) | X,Y |
| `precip_days_perdecade_4in` | Total days with precip ≥ 4 inches over full period | days | derived (as above) | X,Y |
| `precip_days_1mm` | Mean annual number of days with precip ≥ 1 mm | days/year | derived (as above) | X,Y |
| `precip_days_1in` | Mean annual number of days with precip ≥ 1 inch | days/year | derived (as above) | X,Y |
| `precip_days_2in` | Mean annual number of days with precip ≥ 2 inches | days/year | derived (as above) | X,Y |
| `precip_days_3in` | Mean annual number of days with precip ≥ 3 inches | days/year | derived (as above) | X,Y |
| `precip_days_4in` | Mean annual number of days with precip ≥ 4 inches | days/year | derived (as above) | X,Y |
| `lat` | Gridpoint latitude | (unitless coord) | `XLAT` | X,Y |
| `lon` | Gridpoint longitude | (unitless coord) | `XLONG` | X,Y |

## `snowfall_mean_processing.py`

Output: `stats_snowfall_*.nc`. Daily total snowfall is built by
differencing accumulated `SNOWNC` at the 0th hour of consecutive days.

| Variable | Description | Units | WRF Source Variable | Shape |
|---|---|---|---|---|
| `yr_daily_snowfall` | Mean daily snowfall across all days | mm/day | `SNOWNC` (derived) | X,Y |
| `winter_daily_snowfall` | Mean daily snowfall in winter (Nov–Mar) | mm/day | `SNOWNC` (derived) | X,Y |
| `winter_daily_snowfall_p90` | 90th percentile of daily snowfall in winter (Nov–Mar) | mm/day | `SNOWNC` (derived) | X,Y |
| `winter_daily_snowfall_p95` | 95th percentile of daily snowfall in winter (Nov–Mar) | mm/day | `SNOWNC` (derived) | X,Y |
| `winter_daily_snowfall_p99` | 99th percentile of daily snowfall in winter (Nov–Mar) | mm/day | `SNOWNC` (derived) | X,Y |
| `yr_tot_snowfall` | Mean annual total snowfall (sum of daily snowfall ÷ 10 years) | mm | `SNOWNC` (derived) | X,Y |
| `yr_days_snowfall_1cm` | Mean annual number of days with snowfall ≥ 10 mm (liquid-equivalent) | days/year | `SNOWNC` (derived) | X,Y |
| `max_snowfall_1day` | Maximum single-day total snowfall over full period | mm | `SNOWNC` (derived) | X,Y |
| `max_snowfall_5day` | Maximum 5-day total snowfall over full period | mm | `SNOWNC` (derived) | X,Y |
| `max_snowfall_10day` | Maximum 10-day total snowfall over full period | mm | `SNOWNC` (derived) | X,Y |
| `snow_depth` | Mean liquid-water-equivalent snow depth by gridpoint, by day of year (365-day calendar) | kg/m² | `SNOW` | X,Y,doy |
| `snow_depth_height` | Mean physical snow depth by gridpoint, by day of year (365-day calendar) | m | `SNOWH` | X,Y,doy |
| `winter_snow_depth` | Mean winter (Nov–Apr) liquid-water-equivalent snow depth | kg/m² | `SNOW` (derived) | X,Y |
| `winter_snow_depth_height` | Mean winter (Nov–Apr) physical snow depth | m | `SNOWH` (derived) | X,Y |
| `snow_cover_days` | Mean days per year with snow cover | days/year | `SNOWC` (derived) | X,Y |
| `lat` | Gridpoint latitude | (unitless coord) | `XLAT` | X,Y |
| `lon` | Gridpoint longitude | (unitless coord) | `XLONG` | X,Y |

## `fvcom_forcing_processing.py`

Not a statistics file — writes one record per raw WRF timestep, in the
NetCDF layout FVCOM's structured-grid surface-forcing reader expects.
Output: `fvcom_forcing_*.nc`.

| Variable | Description | Units | WRF Source Variable | Shape |
|---|---|---|---|---|
| `Times` | Timestamp string for each record | `YYYY-MM-DD_HH:MM:SS` | `Times` | time |
| `time` | Time coordinate | days since 1858-11-17 (Modified Julian Day) | derived from `Times` | time |
| `XLAT` | Gridpoint latitude | degree_north | `XLAT` | X,Y |
| `XLONG` | Gridpoint longitude | degree_east | `XLONG` | X,Y |
| `LANDMASK` | Land mask (1 = land, 0 = water) | (unitless) | `LANDMASK` | X,Y |
| `uwind_speed` | 10 m U-component wind speed | m/s | `U10` | time,X,Y |
| `vwind_speed` | 10 m V-component wind speed | m/s | `V10` | time,X,Y |
| `air_temperature` | 2 m air temperature | deg C | `T2` (converted from K) | time,X,Y |
| `relative_humidity` | 2 m relative humidity | % | derived from `T2`, `Q2`, `PSFC` (Bolton 1980 saturation vapor pressure) | time,X,Y |
| `air_pressure` | Surface air pressure | Pa | `PSFC` | time,X,Y |
| `long_wave` | Downward longwave radiation flux | W/m² | `GLW` | time,X,Y |
| `short_wave` | Downward shortwave radiation flux | W/m² | `SWDOWN` | time,X,Y |
| `precip` | Precipitation rate | m/s | `RAINC`, `RAINNC` (differenced, derived) | time,X,Y |
| `evap` | Evaporation rate (negative = water loss to the atmosphere) | m/s | `LH` (derived, Rho_water/latent-heat-of-vaporization scaling) | time,X,Y |
