---
name: claw-vpd-transpiration
description: Transpiration response to vapour pressure deficit (VPD) from gravimetric phenotyping (pot weights on balances or load cells). Computes VPD from air temperature and relative humidity, whole-plant transpiration rate (TR) from pot-weight loss with irrigation, spike and gap screening, and fits linear versus two-segment (breakpoint) TR-VPD regressions chosen by AICc and the Davies test. Reports the VPD breakpoint (threshold) with its 95% CI, the slopes before and after, the % slope reduction, and night-time TR. Use when the user has pot weights or TR and VPD data from growth chambers or greenhouses and asks about limited transpiration (the LT trait), water-saving or slow-wilting traits, TR response to evaporative demand, a VPD breakpoint or threshold, nocturnal transpiration, or segmented regression of transpiration on VPD.
license: MIT
compatibility: R >= 4.2 with the segmented package (>= 2.0); openssl optional for SHA-256 checksums. No internet needed.
metadata:
  version: "1.0.0"
  domain: physiology
  language: R
  status: tested
  repository: https://github.com/bgtamang/AgriClaw
---

# Transpiration response to VPD (limited-transpiration trait)

## What it does

Turns pot-weight and air temperature/RH logs (or a table of TR and VPD) into a TR-VPD response per genotype or per pot. It computes VPD, estimates TR for each VPD step after screening out irrigations, spikes and logger gaps, then compares a straight line with a two-segment line. When the two-segment model wins, the breakpoint is the VPD above which TR rises more slowly: the "limited transpiration under high VPD" water-saving trait (Sinclair et al. 2017). The workflow follows the published gravimetric platform of Tamang & Sadok (2018), with an open-source model fit in place of GraphPad Prism.

## When to use / when not to use

- Use for growth-chamber or greenhouse experiments where VPD is stepped (usually by changing RH at constant temperature) and pots sit on balances; for screening genotypes for a VPD breakpoint; for night-time TR from dark periods.
- Use with `--input` when TR and VPD per pot and step are already computed.
- Not for leaf-level gas exchange (stomatal conductance from an LI-6800): use `claw-licor-6800`, then `claw-stomatal-models` for gs-VPD models.
- Not for transpiration response to soil drying (FTSW dry-down curves) or field sap flow: not covered.
- Not for unsealed pots: weight loss then includes soil evaporation and is not transpiration.

## Inputs

**Raw logger data** (all timestamps ISO 8601, e.g. `2026-07-01 13:05:00`, read as clock time):

| File | Required columns | Optional |
|---|---|---|
| `--weights` | `plant_id`, `timestamp`, `weight_g` (one row per reading) | |
| `--plants` | `plant_id`, `genotype` (or `--group`), `leaf_area_cm2` | `chamber` |
| `--env` | `timestamp`, `air_temp_c`, `rh_pct` | `chamber`, `sensor` |
| `--steps` | `step`, `start`, `end` | `phase` (`day`/`night`), `chamber` |

Without `--steps`, readings are binned into `--bin-min` windows and all are treated as daytime (a warning says so). If `env` has a `chamber` column, `plants` must too; VPD for a pot comes only from its own chamber's sensors, averaged over sensors.

**Precomputed data** (`--input`): one row per pot and step with `genotype` (or `--group`), `vpd_kpa` (or `air_temp_c` + `rh_pct`) and a TR column (`--tr-col`, default `tr_mg_m2_s`).

**Ask the user before running** (never guess):

1. Were pots sealed against soil evaporation, and were plants well watered? If not, stop: the LT trait is defined for well-watered plants with weight loss = transpiration.
2. Logging interval and balance resolution, so `--irrigation-g` exceeds the largest normal loss per reading (at least 5x; the script warns if not).
3. Do weight and climate loggers share one clock? Daylight-saving or logger-setup offsets of exactly 1 h are common; correct with `--env-offset-min`.
4. The step schedule: start and end of each VPD step, and how many minutes plants need to settle after each change (set the step windows or `--discard-first-min`). Which periods were dark (`phase = night`)?
5. Was VPD changed through RH at constant temperature? If temperature also changed, the response mixes temperature and VPD effects; say so in the report.
6. Leaf area: units (cm2), one-sided, and when it was measured relative to the TR measurements.
7. For `--input`: the TR units (put them in the column name).

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/vpd_transpiration.R" --demo --out demo_vpd
Rscript "${CLAUDE_SKILL_DIR}/scripts/vpd_transpiration.R" --weights weights.csv --plants plants.csv \
  --env env.csv --steps steps.csv --out results/
Rscript "${CLAUDE_SKILL_DIR}/scripts/vpd_transpiration.R" --input tr_vpd.csv --tr-col tr_mg_m2_s --out results/
Rscript "${CLAUDE_SKILL_DIR}/scripts/vpd_transpiration.R" --weights w.csv --plants p.csv --env e.csv \
  --steps s.csv --level plant --selection aicc      # per pot, AICc only (published workflow)
```

`--help` lists every option. Key ones: `--svp` (`buck1981_moist`, `buck1981`, `fao56`), `--selection` (`aicc+davies` or `aicc`), `--ci-method` (`score` or `delta`), `--irrigation-g`, `--settle-min`, `--outlier-k`, `--outlier-window`, `--discard-first-min`, `--env-offset-min`, `--min-points`, `--min-side`, `--seed`.

## Outputs

- `report.md`: warnings first, methods with the QC counts, the fit table, night-time TR, limitations, references.
- `tables/tr_vpd_fits.csv`: one row per genotype (or pot). `model` is `linear`, `segmented` or `segmented (slope increases)`; `breakpoint_kpa` with `breakpoint_lo`/`breakpoint_hi` (95% CI); `slope1`, `slope2` (TR units per kPa) with SE and CI; `slope_reduction_pct` = 100 (slope1 - slope2) / slope1; `lin_slope`, `lin_r2`; `aicc_linear`, `aicc_segmented`, `delta_aicc` (linear minus segmented; positive favours the breakpoint), `weight_segmented` (Akaike weight), `davies_p`; `n_below_bp`, `n_above_bp`; `note`. The two-segment columns are filled whenever that fit converged, **even when `model` is `linear`**: always read `model` first.
- `tables/step_tr_vpd.csv` (raw mode): per pot and step, `tr_mg_m2_s` (with its least-squares SE, which ignores autocorrelated logger noise and so is a lower bound), `tr_mmol_m2_s`, `tr_g_h` (whole plant), `vpd_kpa`, `vpd_sd_kpa` (within-step stability), `sensor_sd_kpa`, mean `air_temp_c` and `rh_pct`, `n_readings`, `n_excluded`, `note`.
- `tables/weight_flags.csv`: every excluded reading and its reason; `tables/night_tr.csv` when the steps file has night windows.
- `figures/tr_vpd_fits_<n>.png` (data, linear and two-segment fits, breakpoint CI band) and `figures/weight_qc.png` (pot weights with flagged readings).
- `reproducibility/`: `commands.sh`, `sessionInfo.txt`, checksums.

## Expert decisions

- **VPD** = SVP(T) x (1 - RH/100), computed for each climate reading and then averaged over the step, because SVP is non-linear in T. Default SVP (`buck1981_moist`): 0.61375 exp(17.502 T / (240.97 + T)) kPa, the coefficient of the published gravimetric workflow (Tamang & Sadok 2018 repository). It equals Buck's (1981) 0.61121 kPa times an enhancement factor of about 1.0042 for moist air near 1000 hPa (Buck tabulates 1.0040-1.0043 at 20-30 degC). `buck1981` drops the enhancement factor; `fao56` is 0.6108 exp(17.27 T / (T + 237.3)) (Allen et al. 1998). The three differ by < 0.7 % from 5 to 45 degC, far less than a +/-2 % RH sensor error (about +/-0.08 kPa at 30 degC). RH of 100-103 % is set to 100 % and counted.
- **TR** per pot and step = minus the least-squares slope of weight on time over the retained readings, with one intercept per segment so an irrigation inside a step moves the level but not the slope. For clean data this equals -delta weight / delta t; the regression uses every reading, not just the two ends. TR (mg H2O m-2 s-1) = loss (g/min) x 1000 / 60 / (leaf area, m2), the unit of Tamang & Sadok (2018).
- **Quality control** (adjustable rules of thumb; counts go in the report):
  - gap: an interval > `--gap-factor` (3) x the median logging interval starts a new segment;
  - jump: the median of the next h readings minus the median of the previous h (h = half of `--outlier-window`), corrected for the local trend, exceeds `--irrigation-g` (5 g): a rise is an irrigation, a drop a disturbance; a new segment starts and the next `--settle-min` (30) minutes are excluded (drainage, wet-surface evaporation, balance settling);
  - spike: within each segment, after removing the trend, readings more than `--outlier-k` (3) robust SDs (1.4826 x MAD) from a running median of `--outlier-window` (15) readings (Hampel 1974 identifier; threshold 2.5-3 per Leys et al. 2013). On clean data about 0.3-1 % of readings are flagged by chance; that removes noise, not bias.
  - a pot x step needs `--min-readings` (10) retained readings.
- **Model**: line TR = a + b VPD versus TR = a + b1 VPD + (b2 - b1)(VPD - psi)+, fitted with `segmented` (Muggeo 2003) using bootstrap restarts against local optima (Wood 2001) with a fixed `--seed`. `segmented` resets the random-number generator for its restarts; the script restores the caller's stream afterwards, so repeated fits stay independent.
- **AICc** = n ln(RSS/n) + 2k + 2k(k+1)/(n - k - 1) (Hurvich & Tsai 1989), k = 3 for the line and 5 for the two-segment model (each counts the residual variance). This reproduces GraphPad Prism's comparison in the published workflow exactly (test 4). Akaike weights follow Burnham & Anderson (2002).
- **Selection** (default `aicc+davies`): the breakpoint model needs lower AICc **and** Davies test p < `--alpha` (0.05; Davies 2002, the linear-model version in `segmented::davies.test`). The breakpoint is not defined under the linear null, so AICc alone over-selects breakpoints: on 200 simulated linear data sets (design below) AICc alone chose a breakpoint in 18.5 %, the default in 2.5 %. The price is power for weak responses: with slopes 30 then 18 the default detected the breakpoint in 82 % of 200 data sets, AICc alone in 92.5 %. `--selection aicc` reproduces the published Prism/AICc workflow; report which rule was used.
- A breakpoint fit is only tried with >= `--min-points` (10) points and >= 4 distinct VPD levels, and kept only with >= `--min-side` (3) points on each side (rules of thumb). With 7 points per pot, AICc's small-sample penalty makes a breakpoint nearly unselectable, so pooled genotype fits (the published practice) are the default.
- **Breakpoint CI**: score-based interval (Muggeo 2017; `confint.segmented(method = "score")`), which is asymmetric and wider. In 120 simulated data sets per case it covered the true breakpoint in 94 % (slopes 30 then 8) and 92 % (30 then 18) of runs, against 84 % and 74 % for the delta method (`--ci-method delta`, the `segmented` default).
- `slope2 > slope1` is reported as `segmented (slope increases)`: a breakpoint, but not limited transpiration.

## Checks before reporting

Stop and tell the user when:

- any warning appears at the top of `report.md` (pass every one on);
- `figures/weight_qc.png` shows steps or drifts that were not flagged, or > 20 % of readings were excluded;
- TR does not increase with VPD at low VPD (`lin_slope` <= 0), or daytime TR is <= 0: wrong units, unsealed pots, water stress or condensation;
- `vpd_sd_kpa` > 0.2 kPa within steps (VPD not stable) or chamber sensors disagree (`sensor_sd_kpa` > 0.3);
- the breakpoint CI is wider than about 1.5 kPa (three 0.5-kPa steps) or reaches the edge of the VPD range, or a group spans < 1.5 kPa of VPD: the threshold is not determined, whatever the model column says;
- `note` reports "AICc prefers a breakpoint but Davies p = ...": evidence is weak; report both criteria;
- night TR is a large fraction of daytime TR (> ~30 %): check the night window and lights.

Also report: genotype-level fits pool repeated steps of the same pots as independent points (as in the published workflow), so CIs are somewhat narrow; a breakpoint summarises a gradual stomatal response and depends on the VPD steps imposed.

## Validation

`tests/test_vpd_transpiration.R` (55 checks) covers: VPD against hand calculations and the FAO-56 worked examples (3 decimals); per-chamber VPD assignment checked by hand; TR unit conversion (0.12 g/min over 500 cm2 = 40 mg m-2 s-1); an irrigation (+100 g) and a spike found at the right minute and TR recovered within 2 %; AICc difference and Akaike weight equal to a Prism output from the published repository (3.786 vs 3.788; 0.869 vs 0.8692); and 30 independent simulated data sets per case (7 VPD steps x 6 plants, residual SD 4 mg m-2 s-1, breakpoint 2.2 kPa, slopes 30 then 8): breakpoint within +/-0.2 kPa in 25/30 and within +/-0.3 kPa in 29/30, all classified segmented, median slopes 29.7 and 7.7, score CI coverage 29/30; linear truth classified linear in 30/30 by default and 28/30 by AICc alone. In 200-data-set runs of the same design: within +/-0.2 kPa 86 %, within +/-0.3 kPa 96 %, median error 0.08 kPa. The demo recovers breakpoints of 2.00 and 2.61 kPa (true 2.0 and 2.6). Limits: a weak slope change (30 to 18, SD 5) was located within +/-0.2 kPa in only 42 % of 200 data sets (median error 0.24 kPa); per-pot fits and data with real stomatal lags were not validated.

## References

- Allen RG, Pereira LS, Raes D, Smith M (1998) Crop evapotranspiration: guidelines for computing crop water requirements. FAO Irrigation and Drainage Paper 56. FAO, Rome. https://www.fao.org/4/x0490e/x0490e00.htm
- Buck AL (1981) New equations for computing vapor pressure and enhancement factor. Journal of Applied Meteorology 20:1527-1532. https://doi.org/10.1175/1520-0450(1981)020<1527:NEFCVP>2.0.CO;2
- Burnham KP, Anderson DR (2002) Model Selection and Multimodel Inference, 2nd edn. Springer, New York. https://doi.org/10.1007/b97636
- Davies RB (1987) Hypothesis testing when a nuisance parameter is present only under the alternative. Biometrika 74:33-43. https://doi.org/10.1093/biomet/74.1.33
- Davies RB (2002) Hypothesis testing when a nuisance parameter is present only under the alternative: linear model case. Biometrika 89:484-489. https://doi.org/10.1093/biomet/89.2.484
- Fletcher AL, Sinclair TR, Allen LH (2007) Transpiration responses to vapor pressure deficit in well watered 'slow-wilting' and commercial soybean. Environmental and Experimental Botany 61:145-151. https://doi.org/10.1016/j.envexpbot.2007.05.004
- Hampel FR (1974) The influence curve and its role in robust estimation. Journal of the American Statistical Association 69:383-393. https://doi.org/10.1080/01621459.1974.10482962
- Hurvich CM, Tsai C-L (1989) Regression and time series model selection in small samples. Biometrika 76:297-307. https://doi.org/10.1093/biomet/76.2.297
- Leys C, Ley C, Klein O, Bernard P, Licata L (2013) Detecting outliers: do not use standard deviation around the mean, use absolute deviation around the median. Journal of Experimental Social Psychology 49:764-766. https://doi.org/10.1016/j.jesp.2013.03.013
- Muggeo VMR (2003) Estimating regression models with unknown break-points. Statistics in Medicine 22:3055-3071. https://doi.org/10.1002/sim.1545
- Muggeo VMR (2017) Interval estimation for the breakpoint in segmented regression: a smoothed score-based approach. Australian & New Zealand Journal of Statistics 59:311-322. https://doi.org/10.1111/anzs.12200
- Sadok W, Sinclair TR (2009) Genetic variability of transpiration response to vapor pressure deficit among soybean cultivars. Crop Science 49:955-960. https://doi.org/10.2135/cropsci2008.09.0560
- Schoppach R, Sadok W (2012) Differential sensitivities of transpiration to evaporative demand and soil water deficit among wheat elite cultivars indicate different strategies for drought tolerance. Environmental and Experimental Botany 84:1-10. https://doi.org/10.1016/j.envexpbot.2012.04.016
- Sinclair TR, Devi J, Shekoofa A, Choudhary S, Sadok W, Vadez V, Riar M, Rufty T (2017) Limited-transpiration response to high vapor pressure deficit in crop species. Plant Science 260:109-118. https://doi.org/10.1016/j.plantsci.2017.04.007
- Tamang BG, Sadok W (2018) Nightly business: links between daytime canopy conductance, nocturnal transpiration and its circadian control illuminate physiological trade-offs in maize. Environmental and Experimental Botany 148:192-202. https://doi.org/10.1016/j.envexpbot.2017.11.016 (scripts: https://github.com/bgtamang/Gravimetric-phenotyping-of-water-saving-traits)
- Tamang BG, Schoppach R, Monnens D, Steffenson BJ, Anderson JA, Sadok W (2019) Variability in temperature-independent transpiration responses to evaporative demand correlate with nighttime water use and its circadian control across diverse wheat populations. Planta 250:115-127. https://doi.org/10.1007/s00425-019-03151-0
- Wood SN (2001) Minimizing model fitting objectives that contain spurious local minima by bootstrap restarting. Biometrics 57:240-244. https://doi.org/10.1111/j.0006-341X.2001.00240.x
