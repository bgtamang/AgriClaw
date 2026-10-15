---
name: claw-water-balance
description: Daily crop soil-water balance with the FAO-56 single crop coefficient method - Kc curve, crop evapotranspiration (ETc), total and readily available water (TAW, RAW), water-stress coefficient Ks, actual ET, root-zone depletion, deep percolation, water-stress days, and optional irrigation scheduling at a management-allowed depletion (MAD). Takes daily rain and reference ET0 (from claw-weather or a CSV) and soil water holding capacity (from claw-soil SSURGO data, field capacity and wilting point, or TAW). Use when the user asks how much water a crop used, when or whether a field was drought stressed, how many stress days occurred in a growth stage, how to schedule irrigation, or wants crop water use rather than just rainfall minus ET0.
license: MIT
compatibility: Python >= 3.9 with matplotlib; standard library otherwise. No internet needed (inputs come from claw-weather / claw-soil outputs or user files).
metadata:
  version: "1.0.0"
  domain: field
  language: Python
  status: tested
  summary: "FAO-56 daily soil-water balance: crop water use, water stress, irrigation timing"
  repository: https://github.com/bgtamang/AgriClaw
---

# FAO-56 root-zone water balance

## What it does

Runs the FAO-56 (Allen et al. 1998) single-crop-coefficient soil water balance day by day from planting to the end of the season. It builds the Kc curve from stage lengths and Kc_ini/Kc_mid/Kc_end, adjusting Kc_mid and Kc_end for wind and humidity. It computes ETc, the stress coefficient Ks from root-zone depletion, and ETc adj. Outputs are stress days by stage, deep percolation and, optionally, irrigation dates and net depths.

## When to use / when not to use

- Use for crop water use, water-stress timing (e.g. stress around silking or flowering), comparing rainfed seasons, or a first irrigation schedule.
- Run **claw-weather** first for daily precipitation and ET0, and **claw-soil** for SSURGO available water storage (US sites).
- Don't use it for paddy rice (FAO-56 handles flooded rice differently), or for fields with a shallow water table or tile drainage unless the user can supply capillary rise or accepts the limitation.
- Don't use it as a calibrated yield model. Pair it with yield-response tools only with local calibration.

## Inputs

| Input | Flag | Notes |
|---|---|---|
| Weather | `--weather FILE.csv` | Columns `date`, `precip_mm`, `et0_mm` (FAO-56 grass reference ET). Optional: `rh_min_pct`, `wind_2m_m_s` (for the Kc climate adjustment), `runoff_mm`, `irrigation_mm` (net, applied), `capillary_rise_mm`. claw-weather's `tables/daily_weather.csv` works directly. No season day may be missing. |
| Crop | `--crop` | `maize_grain`, `soybean`, `spring_wheat`, `winter_wheat`, `barley`, `sorghum_grain`, `cotton`, `dry_bean`, `sunflower`, `potato`, or `custom` (then give `--kc-ini --kc-mid --kc-end --stage-days --root-depth-m --p`). |
| Planting date | `--plant-date` | Day 1 of the Kc curve (green-up date for winter crops after dormancy if you model from spring). |
| Soil water | exactly one of `--soil-dir`, `--theta-fc` + `--theta-wp`, `--taw-mm` | `--soil-dir` is a claw-soil output folder (uses `tables/aws_profile.csv`); θ values are volumetric (m³/m³). |
| Overrides | `--stage-days 30,40,50,30`, `--root-depth-m`, `--p`, `--kc-*`, `--crop-height-m` | Local values beat FAO-56 defaults. |
| Start condition | `--initial-depletion-mm` | Default 0 = field capacity at planting. |
| Irrigation | `--irrigate [--mad 0.5]` | Schedules a refill to field capacity when depletion reaches MAD × TAW (default MAD = p, i.e. RAW). |

Ask the user; never guess:

- planting date, and observed stage dates if known (stage lengths drive seasonal ET);
- soil water holding: their soil data, or permission to use SSURGO for the field location;
- whether the field is irrigated, and whether it is tile-drained or has a shallow water table;
- whether the soil was at field capacity at planting (after a wet spring, usually yes; after a dry winter, no);
- root depth if a restrictive layer, shallow soil or young crop limits rooting.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
python "${CLAUDE_SKILL_DIR}/scripts/water_balance.py" --demo
python "${CLAUDE_SKILL_DIR}/scripts/water_balance.py" --weather results/weather/tables/daily_weather.csv --soil-dir results/soil --crop maize_grain --plant-date 2025-04-25 --root-depth-m 1.2 --stage-days 25,45,50,35 --out results/wb
python "${CLAUDE_SKILL_DIR}/scripts/water_balance.py" --weather my_weather.csv --crop soybean --plant-date 2025-05-20 --theta-fc 0.30 --theta-wp 0.14 --irrigate --mad 0.5
python "${CLAUDE_SKILL_DIR}/scripts/water_balance.py" --weather my_weather.csv --crop custom --kc-ini 0.4 --kc-mid 1.15 --kc-end 0.5 --stage-days 20,30,60,25 --root-depth-m 0.9 --p 0.5 --taw-mm 140 --plant-date 2025-05-20
```

The demo is rainfed field maize planted 2025-04-25 near Champaign, Illinois. It uses claw-weather's demo weather and claw-soil's Drummer silty clay loam AWS profile, with Zr = 1.2 m.

## Outputs

- `report.md`: warnings first, then set-up, crop-coefficient table (values and sources), season water balance (P, runoff, irrigation, ETc, ETc adj, deep percolation, start/end depletion, mass-balance residual), stress days by stage, irrigation events, figures and assumptions.
- `tables/daily_water_balance.csv`: per day `stage`, `kc`, `et0_mm`, `etc_mm`, `p`, `raw_mm`, `taw_mm`, `precip_mm`, `runoff_mm`, `eff_precip_mm`, `irrigation_mm`, `deep_percolation_mm`, `dr_start_mm`, `dr_after_wetting_mm`, `ks`, `etc_adj_mm`, `dr_end_mm`, `stressed`, `et_limited_by_storage`.
- `tables/season_summary.csv`: totals, `stress_days`, `first_stress_date`, `longest_stress_run_days`, `irrigation_events`, `mass_balance_residual_mm`.
- `tables/crop_parameters.csv`: every parameter with its source, including the u2/RHmin means used in eq. 62/65.
- `tables/irrigation_events.csv`: dates and net depths (divide by application efficiency for gross depths).
- `figures/depletion.png` (Dr with RAW and TAW, stress days, rain/irrigation) and `figures/kc_et.png`, plus `reproducibility/`.

## Expert decisions

- **Kc curve** (FAO-56 eq. 66). Kc_ini holds in the initial stage, rises linearly to Kc_mid during development, holds through mid-season, then falls linearly to Kc_end over the late season. Day 1 is the planting date.
- **Crop defaults.**
  - Kc and maximum crop height come from FAO-56 Table 12. Kc_ini is the crop group value, a general value that is higher with frequent wetting (Table 12 footnote 1).
  - Stage lengths are one Table 11 row per crop (listed in `CROPS` and in the report), for example maize 30/40/50/30 days (April planting, Spain/California) and soybean 20/35/60/25 (Central USA). They are generic: always ask for local values.
  - Where Table 12 gives an unexplained range, the midpoint is used (sorghum Kc_mid 1.05, cotton Kc_mid 1.175 and Kc_end 0.60, sunflower Kc_mid 1.075, soybean h 0.75 m).
  - Maize Kc_end is 0.35 for grain dried in the field; use 0.60 for high-moisture harvest (footnote 11).
- **Climate adjustment** (eq. 62 for Kc_mid, eq. 65 for Kc_end). It uses the mean u2 and RHmin over the stage's days, with FAO's limits (1 ≤ u2 ≤ 6 m/s, 20 ≤ RHmin ≤ 80%, 0.1 ≤ h ≤ 10 m). Kc_end is adjusted only when its table value is ≥ 0.45. Verified on FAO-56 Example 27 (1.07 and 1.30) and Box 15 (1.19).
- **Soil water.** TAW = 1000 (θFC − θWP) Zr (eq. 82), or SSURGO available water storage interpolated at Zr. RAW = p × TAW (eq. 83). p comes from Table 22 and is adjusted daily as p + 0.04 (5 − ETc), limited to 0.1–0.8 (Table 22 footnote 2); use `--no-p-adjust` to keep it constant, as FAO-56 Example 37 does.
- **Root depth.** The default is the Table 22 range end: the larger value for rainfed or stress modeling, the smaller for irrigation scheduling (Table 22 footnote 1). It is constant all season; root growth is not modeled.
- **Daily update** (eqs. 84–86): Ks = (TAW − Dr)/(TAW − RAW) when Dr > RAW, otherwise 1. ETc adj = Ks × Kc × ET0 (eq. 81), and 0 ≤ Dr ≤ TAW.
  - Rain, irrigation and capillary rise enter early in the day; water above field capacity is lost to deep percolation the same day; then Ks and ET apply. This is the convention of FAO-56 Example 38.
  - FAO-56 eq. 88 subtracts the day's ET before computing drainage, which gives less drainage on wetting days.
  - ETc adj is capped at the water left in the root zone, so the balance closes exactly. The tests check closure to 1e-9 mm.
- **Simplifications.** Runoff and capillary rise are 0 unless supplied as columns (FAO-56: capillary rise ≈ 0 when the water table is more than ~1 m below the root zone). Optionally, daily rain below a fraction of ET0 can be ignored (`--ignore-rain-below-et0-frac 0.2`, as FAO-56 suggests for small events); the default is 0.
- **Initial depletion.** Defaults to 0 (FAO-56: near field capacity after heavy rain or irrigation). This is flagged as an assumption in every report.
- **Irrigation scheduling.** Irrigate on a day that starts with Dr ≥ MAD × TAW (default MAD = p), with a net depth equal to Dr, refilling to field capacity (FAO-56: irrigate at or before RAW is depleted, net depth ≤ Dr). With MAD = p no stress day occurs (tested).
- **Validation.** Reproduces FAO-56 Example 37 (daily Ks 1.00 → 0.62 and Dr 61.0 → 104.5 mm, within the table's rounding), Example 36 (TAW/RAW), Example 28 (Kc at days 20/40/70/95), and the analytic stress onset and geometric decline of TAW − Dr.

## Checks before reporting

- Read the Warnings section of `report.md` and pass every warning on, especially generic stage lengths, assumed initial depletion, and a shallow water table.
- Stop and tell the user if:
  - soil water holding had to be assumed (the script refuses unless one source is given; never invent θ values);
  - SSURGO shows a shallow water table or poor drainage and the field may not be tile-drained. Stress is then probably overestimated, because capillary rise is ignored;
  - the root depth exceeds the described soil profile, or a restrictive layer is likely;
  - weather data stop before the season ends (results are partial);
  - many days hit `et_limited_by_storage` (very small TAW: check root depth and units).
- Stress days are relative to FAO-56's p threshold. Present them as model estimates, and compare them with field observations (leaf rolling, soil moisture) where possible.

## References

- Allen, R.G., Pereira, L.S., Raes, D. & Smith, M. (1998). Crop evapotranspiration: guidelines for computing crop water requirements. FAO Irrigation and Drainage Paper 56. FAO, Rome. ISBN 92-5-104219-5. https://www.fao.org/4/x0490e/x0490e00.htm (Chapter 6: Tables 11–12, eqs. 62, 65, 66, Examples 27–28, Box 15; Chapter 8: eqs. 81–88, Table 22, Examples 36–38)
- Soil Survey Staff, Natural Resources Conservation Service, USDA. Soil Survey Geographic (SSURGO) Database (available water storage), via claw-soil. https://sdmdataaccess.sc.egov.usda.gov/
