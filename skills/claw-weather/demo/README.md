# Demo data

`openmeteo_champaign_2025.json` is an unmodified response from the Open-Meteo Historical Weather API (archive endpoint, default "best match" model), saved on 2026-10-15 for this request:

```
https://archive-api.open-meteo.com/v1/archive?latitude=40.08&longitude=-88.24&timezone=auto&wind_speed_unit=ms
  &start_date=2025-01-01&end_date=2025-12-31
  &daily=temperature_2m_max,temperature_2m_min,precipitation_sum,et0_fao_evapotranspiration,
         relative_humidity_2m_min,wind_speed_10m_mean,shortwave_radiation_sum
```

The demo summarises 1 April – 31 October 2025 near Champaign, Illinois, USA. The data are gridded model output for the grid cell centred at 40.105°N, 88.235°W, not station observations.

The files in `tests/fixtures/` are responses saved the same day: the same 2025 year, an ERA5-only request for 1–14 October 2026 whose last 6 days were not yet published, a forecast-API `past_days=14` response, and the HTTP 400 body returned for an end date in the future.

Licence and attribution: Open-Meteo API data are offered under the Creative Commons Attribution 4.0 International licence (CC BY 4.0), https://open-meteo.com/en/licence. Weather data by Open-Meteo.com (Zippenfenig 2023, https://doi.org/10.5281/zenodo.7970649); generated using Copernicus Climate Change Service information (ERA5 and ERA5-Land; Hersbach et al. 2020; Muñoz-Sabater et al. 2021).
