---
name: claw-weather
description: Summarise growing-season weather for a field location or station - daily and monthly temperature, rainfall/precipitation, reference evapotranspiration (ET0), a simple climatic water balance (P - ET0), frost and hard-freeze days, heat days, heavy rain, and dry spells - from free Open-Meteo historical data (ERA5 reanalysis grid) or the user's weather-station CSV. Missing days are counted and reported, never treated as zero, and in-season requests stop at the last complete day. Use when the user asks what the weather was like at a trial site, how much it rained, how many hot or frost days there were, how long the longest dry spell was, or wants weather context for a field season or a site comparison.
license: MIT
compatibility: Python >= 3.9 with matplotlib; standard library otherwise. Internet access to open-meteo.com for --lat/--lon runs; --demo and --input run offline.
metadata:
  version: "1.0.0"
  domain: field
  language: Python
  status: tested
  summary: "Growing-season weather: temperature, rain, reference ET, heat and dry spells"
  repository: https://github.com/bgtamang/AgriClaw
---

# Weather season summary

## What it does

Builds a daily weather table for a period and summarises it: season and monthly temperature, precipitation and FAO-56 reference ET0; a simple climatic water balance (P − ET0); counts of frost, hard-freeze, hot and heavy-rain days with sourced thresholds; and dry and wet spells. Data come from the Open-Meteo archive (gridded reanalysis/model data) or from the user's station CSV. Missing values are counted per variable and never filled with zeros.

## When to use / when not to use

- Use for weather context of a trial or field season, comparing sites or years, and checking for heat, frost, heavy rain or drought periods.
- Use **claw-gdd** for growing degree days, corn stage dates and freeze dates relative to planting.
- Use **claw-water-balance** for crop water use (ETc), soil-water depletion and water stress. P − ET0 here is a climate summary, not a drought index or crop water use.
- Don't use it for forecasts. It never requests days after yesterday.

## Inputs

Ask the user for anything missing. Never guess the period, the location, or the units of a station file.

| Input | Flag | Notes |
|---|---|---|
| Period | `--start-date`, `--end-date` (YYYY-MM-DD) | End defaults to yesterday (Open-Meteo) or the last CSV date. Future end dates are clamped. |
| Location | `--lat`, `--lon` (decimal degrees) | For Open-Meteo. Also recorded in the report for CSV runs. |
| Station data (optional) | `--input FILE.csv` | `date` plus any of `tmax`, `tmin` (needs `--temp-units C\|F`), `precip`, `et0` (needs `--precip-units mm\|in`), `rh_min_pct`, `wind_2m_m_s` or `wind_10m_m_s`, `srad_mj_m2`. Suffixed names (`tmax_c`, `precip_in`, ...) must agree with the flags. |
| Data model (optional) | `--model best_match\|era5\|era5_seamless\|era5_land\|ecmwf_ifs` | `era5` is temporally consistent for multi-year comparisons. `era5_land` has no precipitation or ET0. |

Ask:

- For a station file, which units the columns are in, and the missing-value code. Blank cells must replace codes such as `-9999`; the script stops on out-of-range values.
- Whether the user wants the whole calendar season or the crop season (planting to harvest). The indices count days in the chosen period only.
- For multi-year or climate-trend questions, suggest `--model era5`. Open-Meteo advises ERA5/ERA5-Land, not the IFS-based best match, for consistency over decades.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
python "${CLAUDE_SKILL_DIR}/scripts/weather.py" --demo
python "${CLAUDE_SKILL_DIR}/scripts/weather.py" --lat 40.08 --lon -88.24 --start-date 2025-04-01 --end-date 2025-10-31 --out results/weather
python "${CLAUDE_SKILL_DIR}/scripts/weather.py" --lat 40.08 --lon -88.24 --start-date 2026-05-01   # in-season: stops at yesterday
python "${CLAUDE_SKILL_DIR}/scripts/weather.py" --input station.csv --temp-units F --precip-units in --start-date 2025-04-15 --end-date 2025-10-15
```

The demo is Champaign, Illinois, 1 April – 31 October 2025, from a bundled Open-Meteo response.

## Outputs

- `report.md`: warnings first, then the data-source caveat, season summary with missing-day counts, monthly table with P − ET0, threshold indices with rules and sources, dry spells, figures and references.
- `tables/daily_weather.csv`: `date, tmax_c, tmin_c, precip_mm, et0_mm, rh_min_pct, wind_2m_m_s, srad_mj_m2, source`. Missing values are blank. This file is the `--weather` input of **claw-water-balance**.
- `tables/season_summary.csv` and `tables/monthly_summary.csv`: totals and means over days with data, `n_missing_*` columns, and `complete` (yes/no).
- `tables/threshold_indices.csv`: index, rule, days, missing days, source. `tables/threshold_days.csv` lists the dates.
- `tables/dry_spells.csv`: start, end, length and whether the spell touches missing data or the period edge.
- `tables/missing_days.csv`.
- Figures: `figures/temp_precip.png`, `water_balance.png` and `monthly_summary.png`.
- `raw/openmeteo_responses.json` (API runs) and `reproducibility/`.

## Expert decisions

- **Missing data.** Totals and means use days with data only. Each variable's missing days are counted and shown next to the total, and months with gaps are marked `complete = no`. A total with missing days is a lower bound.
  - A missing day ends a dry or wet spell. Spells that border missing data are flagged, because their true length is unknown.
  - Absent dates in a station file count as missing.
- **Threshold indices.** Counts cover the analysed period only (ETCCDI defines them per year):
  - FD frost days, TN < 0 °C; SU summer days, TX > 25 °C; TR tropical nights, TN > 20 °C.
  - R10mm / R20mm, RR ≥ 10 / 20 mm; Rx1day, the largest 1-day total.
  - CDD and CWD: the longest run of dry (RR < 1 mm) and wet (RR ≥ 1 mm) days.
  - These follow ETCCDI exactly (Zhang et al. 2011; ETCCDI list of 27 indices).
  - Hard-freeze days use TN ≤ −2.22 °C (28 °F), the NOAA NCEI threshold (Arguez 2012).
  - Hot days use TX > 32.2 °C (90 °F) and very hot days TX ≥ 35 °C (95 °F). Neild & Newman (1987, NCH-40) state that rain-fed corn starts to stress above 90 °F and loses yield for each day at or above 95 °F during the critical period.
  - The heat thresholds are crop-specific; adjust them with `--hot-c` and `--very-hot-c` for other crops. `--heavy-rain-mm` adds an ETCCDI-style Rnnmm count.
  - The dry-spell listing length (`--min-dry-spell`, default 7 days) is a rule of thumb for the table only; CDD is always reported.
- **P − ET0.** Monthly precipitation minus FAO-56 grass reference ET0 is presented only as a simple climatic water balance. It is not a drought index (no soil storage, runoff or crop coefficient). The old "mild/moderate/severe deficit" labels were removed because their thresholds had no source.
- **Data source.** The Open-Meteo archive "best match" combines ECMWF IFS (9 km, from 2017) with the ERA5 (0.25°, ~25 km) and ERA5-Land (0.1°, ~11 km) reanalyses. It is a grid cell, not a station.
  - Open-Meteo computes ET0 with FAO-56 Penman–Monteith.
  - Daily mean 10 m wind is converted to 2 m with FAO-56 eq. 47 (u2 = u10 × 0.748). The tests check this against FAO-56 Example 14.
- **In-season dates.** The end date is clamped to yesterday (Open-Meteo returns HTTP 400 for later dates). ERA5 is published about 5 days late. Trailing days with no archive data are filled from the forecast API (`&past_days=`, at most 92 days), labelled `forecast:past_days` and flagged as provisional.

## Checks before reporting

- Read the Warnings section of `report.md` and pass every warning on.
- Stop and tell the user if:
  - any core variable (tmax, tmin, precip, ET0) has many missing days, or a month is incomplete in the period they care about;
  - the source is a grid cell but the question is about field-specific events: one storm, frost in a low spot, irrigation timing;
  - `forecast:past_days` rows are present: those values are provisional;
  - station-file units were not stated by the user.
- Heat-day thresholds are corn-specific. For other crops, ask which thresholds apply before interpreting heat days.
- Never present P − ET0 as crop water stress. Run **claw-water-balance** instead.

## References

- Zhang, X., Alexander, L., Hegerl, G.C., Jones, P., Klein Tank, A., Peterson, T.C., Trewin, B. & Zwiers, F.W. (2011). Indices for monitoring changes in extremes based on daily temperature and precipitation data. WIREs Climate Change 2(6), 851–870. https://doi.org/10.1002/wcc.147 (definitions: https://etccdi.pacificclimate.org/list_27_indices.shtml)
- Neild, R.E. & Newman, J.E. (1987; Purdue reprint RR 4/90). Growing Season Characteristics and Requirements in the Corn Belt. National Corn Handbook NCH-40. https://www.extension.purdue.edu/extmedia/nch/nch-40.html
- Arguez, A. (2012). NOAA's 1981–2010 Climate Normals: Supplemental Normals. NOAA NCDC technical report, in NOAA NCEI Normals Calculation Methodology 2020. https://www.ncei.noaa.gov/data/normals-annualseasonal/1991-2020/doc/Normals_Calculation_Methodology_2020.pdf
- Allen, R.G., Pereira, L.S., Raes, D. & Smith, M. (1998). Crop evapotranspiration: guidelines for computing crop water requirements. FAO Irrigation and Drainage Paper 56. FAO, Rome. https://www.fao.org/4/x0490e/x0490e00.htm
- Zippenfenig, P. (2023). Open-Meteo.com Weather API. Zenodo. https://doi.org/10.5281/zenodo.7970649 (data sources and update delays: https://open-meteo.com/en/docs/historical-weather-api)
- Hersbach, H. et al. (2020). The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society 146(730), 1999–2049. https://doi.org/10.1002/qj.3803
- Muñoz-Sabater, J. et al. (2021). ERA5-Land: a state-of-the-art global reanalysis dataset for land applications. Earth System Science Data 13(9), 4349–4383. https://doi.org/10.5194/essd-13-4349-2021
