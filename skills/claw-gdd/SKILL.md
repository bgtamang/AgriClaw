---
name: claw-gdd
description: Compute growing degree days (GDD, heat units, GDUs, thermal time) from planting with the average or modified (86/50 °F) method, in both °C-days and °F-days, from Open-Meteo weather or a station CSV. Estimates corn (maize) growth stage dates (V6, VT/V16, silking R1, dent, black layer R6) from the cited NCH-40 table in °F-days, optionally scaled to the hybrid's GDD-to-black-layer rating. Reports daylength for soybean instead of guessing R stages, and last spring / first fall freeze and hard-freeze dates. Use when the user asks how many GDD have accumulated since planting, when corn will silk or reach black layer, whether an in-season crop is ahead or behind, or for frost dates at a field site.
license: MIT
compatibility: Python >= 3.9 with matplotlib; standard library otherwise. Internet access to open-meteo.com for --lat/--lon runs; --demo and --input run offline.
metadata:
  version: "1.0.0"
  domain: field
  language: Python
  status: tested
  repository: https://github.com/bgtamang/AgriClaw
---

# Growing degree days, corn stages and freeze dates

## What it does

Computes daily and cumulative growing degree days (GDD) from a planting date, using daily maximum and minimum air temperature from Open-Meteo or the user's station CSV. For corn it dates the stages in the Neild & Newman (1987) NCH-40 table, comparing °F-day thresholds with °F-day accumulation. For soybean it reports GDD and daylength but does not predict stages. It also finds freeze and hard-freeze dates before and after planting.

## When to use / when not to use

- Use for thermal time since planting, corn stage timing, comparing sites or planting dates in heat units, and freeze dates.
- Use **claw-weather** for a whole-season weather summary (rain, ET0, heat and dry spells), and **claw-water-balance** for crop water use and stress.
- Don't use it to predict soybean V or R stages. They depend on photoperiod and maturity group, so point the user to a photothermal model (SoyStage, CROPGRO-Soybean) and give them the daylength table.
- Don't use it for forecasts: it never goes past the last day with observations.
- Canadian Crop Heat Units (CHU) are a different formula; this skill doesn't compute them.

## Inputs

Ask the user for anything missing. Never guess a planting date, crop, temperature unit or base temperature.

| Input | Flag | Notes |
|---|---|---|
| Planting date | `--plant-date YYYY-MM-DD` | Required. Counted as day 1. |
| Crop | `--crop corn\|soybean\|other` | `other` needs `--base-c` (and `--upper-c` for the modified method). |
| Location | `--lat`, `--lon` | Decimal degrees. `--lat` is always required (daylength, frost season). |
| Station data (optional) | `--input FILE.csv --temp-units C\|F` | Columns `date` (YYYY-MM-DD), `tmax`, `tmin` (or `tmax_c`/`tmax_f`, `tmin_c`/`tmin_f`). Units are **never** inferred from the data. |
| End date (optional) | `--end-date` | Default: planting + 210 days, then clamped to the last day with data. |
| Corn hybrid rating (optional) | `--hybrid-gdd 2600` | Seed-company GDD to black layer in **°F-days** (86/50 method). |

Questions to ask:

- For corn, "What is the hybrid's GDD (GDU) rating to black layer, and is it in °F GDD, °C GDD, or Canadian CHU?" Only °F GDD can be used directly. Multiply a °C rating by 1.8. CHU can't be used at all.
- For `--crop other`, "Which base temperature (and upper threshold) does your reference or phenology model use?" For example, wheat models commonly use 0 °C, while NCH-40 lists 40 °F (≈ 4.4 °C) for cool-season crops.
- For a station CSV, ask the temperature unit (°C or °F) and the missing-value code. Blank or `NA` cells are fine. Ask the user to blank out codes such as `-9999`; the script stops on them.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
python "${CLAUDE_SKILL_DIR}/scripts/gdd.py" --demo
python "${CLAUDE_SKILL_DIR}/scripts/gdd.py" --lat 40.08 --lon -88.24 --plant-date 2025-04-25 --crop corn --hybrid-gdd 2650 --out results/gdd
python "${CLAUDE_SKILL_DIR}/scripts/gdd.py" --lat 40.08 --lon -88.24 --plant-date 2025-05-20 --crop soybean
python "${CLAUDE_SKILL_DIR}/scripts/gdd.py" --input station.csv --temp-units F --lat 40.08 --plant-date 2025-04-25 --crop corn
python "${CLAUDE_SKILL_DIR}/scripts/gdd.py" --lat 52.1 --lon 5.2 --plant-date 2025-03-20 --crop other --base-c 0 --method average
```

The demo is corn planted 2025-04-25 near Champaign, Illinois, using the bundled Open-Meteo response.

## Outputs

- `report.md`: warnings first, then inputs, season totals, a corn stage table (or the soybean note), freeze dates, figures and references.
- `tables/daily_gdd.csv`: one row per day from planting. Columns: `gdd_c`, `gdd_f`, `cum_gdd_c`, `cum_gdd_f`, `tmax_c`/`tmin_c` (and °F), `daylength_h` (sunrise to sunset), `daylength_civil_twilight_h`, `filled` (interpolated variables) and `source` (archive or forecast).
- `tables/corn_stages.csv`: for each NCH-40 stage, the table threshold, the threshold used (scaled), the date reached and days after planting. Stages not reached show `not reached as of <last date>`.
- `tables/freeze_events.csv`: every day with Tmin ≤ the freeze threshold, flagged `hard_freeze`, with its period (`before_planting`, `after_planting_before_fall`, `fall`).
- `tables/filled_days.csv`: every interpolated value.
- `figures/gdd_accumulation.png` and `figures/temperature_daylength.png`.
- `raw/openmeteo_responses.json` (API runs only), plus `reproducibility/`.

## Expert decisions

- **Methods** (McMaster & Wilhelm 1997):
  - `average`: GDD = max(0, (Tmax + Tmin)/2 − Tbase), with no upper threshold.
  - `modified`: Tmax and Tmin are each limited to [Tbase, Tupper] before averaging. This is the US corn "86/50" cut-off method described in NCH-40.
  - Corn and soybean default to `modified`, 10/30 °C (50/86 °F).
  - This is not a sine (single- or double-triangle) method.
- **Units.** GDD are computed in °C and reported in both units. °F-days = 1.8 × °C-days exactly when the thresholds are equivalent, and the tests check this.
- **Corn stages.** Thresholds are NCH-40 Table 3 values for a 2700 °F-day hybrid (Neild & Newman 1987): V2 200, V6 475, V12 870, V16 (tip of tassel emerging) 1135, R1 silking 1400, R2 blister 1660, R4 dough 1925, R5 denting 2190, kernels dented 2450, R6 black layer 2700 °F-days.
  - They are compared with cumulative °F-days, counted from the planting day.
  - Stages are reported only for the 86/50 modified method. With any other method or threshold, the script skips them and says why.
  - Lee (2011, AGR-202) reproduces this table and labels the 1135 point "VT". Codes for the Hanway stages are approximate leaf-collar equivalents; the descriptions are NCH-40's own words.
- **Hybrid scaling** (`--hybrid-gdd`): thresholds × rating/2700. This is an adjustable rule of thumb, not a calibrated model. Ratings outside 1800–3600 °F-days stop the run, because they usually mean °C-days or CHU.
- **Soybean.** Development depends on photoperiod and maturity group as well as temperature (Setiyono et al. 2007; Salmerón & Purcell 2016), so no stage dates are produced. The 10/30 °C thresholds are only a thermal-time convention borrowed from corn (adjustable).
  - Reported instead: daylength from the CBM model of Forsythe et al. (1995), both sunrise–sunset (p = 0.8333°) and with civil twilight (p = 6°). It agrees with Open-Meteo's daylight duration within 5 minutes (tested).
  - V stages are omitted too: the only GDD phyllochron estimates we found (58–83 °C-days) are conference abstracts and vary with temperature and sowing date.
  - Refer users to SoyStage (dos Santos et al. 2019; maturity groups 3–6, US Midsouth) or CROPGRO-Soybean in DSSAT (Jones et al. 2003).
- **Freeze thresholds** (NOAA NCEI, Arguez 2012): freeze = Tmin ≤ 0 °C (32 °F); hard freeze = Tmin ≤ −2.22 °C (28 °F). Comparisons are "at or below", as in NCEI. Both are adjustable (`--freeze-c`, `--hard-freeze-c`); killing temperatures vary by crop and stage.
  - The last freeze before planting is searched from January 1 to the planting date.
  - Freezes after planting are listed separately.
  - The first fall freeze is the first one on or after August 1. NCEI computes freeze dates over an August 1 – July 31 year. The script mirrors this to February 1 in the southern hemisphere.
  - The spring/fall summary assumes a summer crop. For winter cereals, read the event table instead.
- **In-season dates.** The end date is clamped to yesterday, so no future day is ever requested; Open-Meteo returns HTTP 400 for end dates after today. ERA5 is published about 5 days late (Open-Meteo documentation). Trailing days with no data are trimmed and filled from the Open-Meteo forecast API (`&past_days=`, at most 92 days). Those days are labelled `forecast:past_days` and flagged as provisional. Nothing is forecast or extrapolated: unreached stages are "not reached as of" the last date.
- **Gaps.** Runs of 3 or fewer missing days with data on both sides are linearly interpolated and listed (an adjustable rule of thumb, `MAX_GAP_DAYS`). Longer or edge gaps inside the season stop the run. Gaps before planting only produce a warning about the freeze search.
- **Data source.** The Open-Meteo archive "best match" blends ECMWF IFS (9 km, from 2017), ERA5 (0.25°) and ERA5-Land (0.1°) model data for a grid cell. It is not a station. `--model era5` gives a temporally consistent reanalysis for multi-decade comparisons.

## Checks before reporting

- Read the Warnings section of `report.md` and pass every warning on.
- Stop and tell the user if:
  - the planting date or `--temp-units` had to be assumed;
  - a hybrid rating's unit is unclear (°F GDD vs °C GDD vs CHU);
  - stages show "not reached" but the user expected them (in-season, or a short `--end-date`);
  - many days were interpolated, or `forecast:past_days` rows fall in the period that matters;
  - they ask for soybean R-stage dates.
- Corn planted in late April or May in the central US Corn Belt usually silks in July. A silking date in August from normal weather suggests a unit or threshold problem.
- Freeze dates from gridded data can miss cold-air pooling in low fields, so state that they come from a grid cell.

## References

- McMaster, G.S. & Wilhelm, W.W. (1997). Growing degree-days: one equation, two interpretations. Agricultural and Forest Meteorology 87(4), 291–300. https://doi.org/10.1016/S0168-1923(97)00027-0
- Neild, R.E. & Newman, J.E. (1987; Purdue reprint RR 4/90). Growing Season Characteristics and Requirements in the Corn Belt. National Corn Handbook NCH-40. https://www.extension.purdue.edu/extmedia/nch/nch-40.html
- Lee, C. (2011). Corn Growth Stages and Growing Degree Days: A Quick Reference Guide. AGR-202, University of Kentucky Cooperative Extension. https://publications.ca.uky.edu/sites/publications.ca.uky.edu/files/agr202.pdf
- Forsythe, W.C., Rykiel, E.J., Stahl, R.S., Wu, H. & Schoolfield, R.M. (1995). A model comparison for daylength as a function of latitude and day of year. Ecological Modelling 80(1), 87–95. https://doi.org/10.1016/0304-3800(94)00034-F
- Arguez, A. (2012). NOAA's 1981–2010 Climate Normals: Supplemental Normals. NOAA NCDC technical report, reproduced in NOAA NCEI Normals Calculation Methodology 2020. https://www.ncei.noaa.gov/data/normals-annualseasonal/1991-2020/doc/Normals_Calculation_Methodology_2020.pdf
- Setiyono, T.D., Weiss, A., Specht, J., Bastidas, A.M., Cassman, K.G. & Dobermann, A. (2007). Understanding and modeling the effect of temperature and daylength on soybean phenology under high-yield conditions. Field Crops Research 100(2–3), 257–271. https://doi.org/10.1016/j.fcr.2006.07.011
- Salmerón, M. & Purcell, L.C. (2016). Simplifying the prediction of phenology with the DSSAT-CROPGRO-soybean model based on relative maturity group and determinacy. Agricultural Systems 148, 178–187. https://doi.org/10.1016/j.agsy.2016.07.016
- dos Santos, C., Salmerón, M. & Purcell, L.C. (2019). Soybean phenology prediction tool for the US Midsouth. Agricultural & Environmental Letters 4(1). https://doi.org/10.2134/ael2019.09.0036 (tool: https://soystage.uark.edu/)
- Jones, J.W., Hoogenboom, G., Porter, C.H., Boote, K.J., Batchelor, W.D., Hunt, L.A., Wilkens, P.W., Singh, U., Gijsman, A.J. & Ritchie, J.T. (2003). The DSSAT cropping system model. European Journal of Agronomy 18(3–4), 235–265. https://doi.org/10.1016/S1161-0301(02)00107-7
- Zippenfenig, P. (2023). Open-Meteo.com Weather API. Zenodo. https://doi.org/10.5281/zenodo.7970649 (historical API documentation: https://open-meteo.com/en/docs/historical-weather-api)
- Hersbach, H. et al. (2020). The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society 146(730), 1999–2049. https://doi.org/10.1002/qj.3803
