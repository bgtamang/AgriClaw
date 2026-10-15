# Demo data

- `daily_weather_champaign_2025.csv`: produced by `claw-weather --demo` (1 April – 31 October 2025, near Champaign, Illinois) from an Open-Meteo Historical Weather API response saved on 2026-10-15.
  - `precip_mm` and `et0_mm` are Open-Meteo's daily precipitation and FAO-56 Penman-Monteith reference ET.
  - `rh_min_pct` is daily minimum relative humidity.
  - `wind_2m_m_s` is the daily mean 10 m wind converted to 2 m with FAO-56 eq. 47.
  - Licence: Open-Meteo data are CC BY 4.0 (https://open-meteo.com/en/licence). Weather data by Open-Meteo.com (https://doi.org/10.5281/zenodo.7970649); generated using Copernicus Climate Change Service information (ERA5).
- `soil_drummer/tables/aws_profile.csv` and `soil_drummer/tables/mapunit_summary.csv`: produced by `claw-soil --demo` from USDA-NRCS Soil Data Access responses for map unit 242963, Drummer silty clay loam, 0 to 2 percent slopes (Champaign County, IL019), saved on 2026-10-15. SSURGO is a US Government work in the public domain.

Demo settings: field maize (`maize_grain`), planted 2025-04-25, root depth 1.2 m (a demo choice within FAO-56 Table 22's 1.0–1.7 m), FAO-56 Table 11 stage lengths, rainfed, root zone at field capacity at planting. Drummer is poorly drained with a seasonal water table near 15 cm and is normally tile-drained in Illinois. The report therefore warns that ignoring capillary rise probably overestimates stress, which is the behaviour to expect on such soils.
