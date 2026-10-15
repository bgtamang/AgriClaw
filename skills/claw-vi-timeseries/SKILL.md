---
name: claw-vi-timeseries
description: Smooth per-plot vegetation-index time series (NDVI, NDRE, EVI, GNDVI) from drone/UAV flights or satellites (Sentinel-2, Landsat, PlanetScope) and extract crop phenology metrics - peak NDVI and peak date, green-up, onset of senescence, senescence rate, stay-green, area under the curve (AUC) and season length - on a days-after-planting or growing-degree-day axis. Fits a double logistic (Beck or Elmore), LOESS or Savitzky-Golay curve, down-weights cloud and shadow drops with upper-envelope weighting, refuses curves for sparse series, and never extrapolates past the last flight. Use when the user has a long table of plot x date x VI (e.g. from claw-vi-extract) and asks for VI curves, phenology dates, senescence timing, canopy duration, NDVI AUC, or time-series features for yield prediction.
license: MIT
compatibility: R >= 4.2 with minpack.lm and signal (phenofit optional, used only by the tests to cross-check the curve formulas). Base R graphics. No internet needed.
metadata:
  version: "1.0.0"
  domain: remote-sensing
  language: R
  status: tested
  summary: "Smooth vegetation-index time series and extract phenology: green-up, peak, senescence, AUC"
  repository: https://github.com/bgtamang/AgriClaw
---

# VI time series and phenology metrics

## What it does

Fits one smooth curve per plot to a vegetation-index (VI) time series and reads phenology metrics off the curve, inside each plot's observed window only. The axis is days after planting (DAP), growing degree days (GDD) or, as a last resort, day of year. Plots with too few dates get model-free summaries instead of a curve, and suspect observations and suspect flight dates are flagged.

## When to use / when not to use

- Use for plot-level VI series from repeated UAV flights or satellite overpasses within one growing season.
- Use it before `claw-yield-prediction` to turn VI series into per-plot features (`metrics.csv` plus `feature_time.csv`).
- Extract plot VIs from imagery first with `claw-vi-extract`. Compute GDD with `claw-gdd`.
- Don't use it for multi-year pixel series that need season detection, or for QA-flag handling across many cycles. Use the R package phenofit (Kong et al. 2022) directly.
- Don't use it for single-date canopy cover (`claw-canopy-cover`) or for statistical comparison of metrics between genotypes (`claw-trial-analysis` on `metrics.csv`).

## Inputs

A long CSV with one row per plot x date:

| Column | Meaning |
|---|---|
| `plot_id` (`--plot-col`) | plot identifier |
| `date` (`--date-col`) | observation date, `YYYY-MM-DD` |
| `ndvi` (`--vi`) | the VI, unscaled (NDVI in -1..1, not x 10000) |
| `planting_date` (optional) | per-plot planting date; otherwise use `--planting-date` |
| QA weight (optional, `--weight-col`) | prior weight in [0, 1], e.g. 0 for cloud-masked pixels |

Optional: a daily GDD CSV (`--gdd`). Use `claw-gdd`'s `tables/daily_gdd.csv`, which has `date` and the daily column `gdd_c` (the `--gdd-col` default; use `gdd_f` for deg F days). For several sites, stack their files and add an environment column (`--env-col`, matching the VI table). Every day from planting to the last flight must be present.

**Ask the user, never guess:**
- **Planting date of each plot or trial.** If it is truly unknown, `--axis doy` is the honest fallback, but plots planted on different dates are then not comparable.
- **Data source.** UAV or satellite decides the weighting (see Expert decisions). Ask whether cloud or QA masks exist; if so, pass them with `--weight-col`.
- **Which VI, and whether it is scaled.** Some platforms store NDVI x 10000.
- **GDD base temperature and units** (deg C or deg F days) when using `--axis gdd`.
- **Early-season use.** If the metrics feed an early-season prediction, ask at what DAP the prediction must be made and run with `--cutoff`.
- **Repeated rows.** Ask how to combine several rows for the same plot and date (two flights in one day, several sensors). The script refuses them.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/vi_timeseries.R" --demo
Rscript "${CLAUDE_SKILL_DIR}/scripts/vi_timeseries.R" --input vi_long.csv --vi ndvi --planting-date 2025-05-12 --keep-cols env,genotype,rep --out results/vi_ts
# satellite series with clouds: upper-envelope weighting, plus QA weights if available
Rscript "${CLAUDE_SKILL_DIR}/scripts/vi_timeseries.R" --input s2_ndvi.csv --planting-date 2025-05-12 --weighting upper --weight-col qa_weight
# thermal-time axis (daily GDD from claw-gdd, one series per environment)
Rscript "${CLAUDE_SKILL_DIR}/scripts/vi_timeseries.R" --input vi_long.csv --axis gdd --gdd daily_gdd.csv --env-col env
# early-season features for prediction at DAP 60: only flights up to DAP 60 are used
Rscript "${CLAUDE_SKILL_DIR}/scripts/vi_timeseries.R" --input vi_long.csv --cutoff 60 --features-at 40,50,60 --out results/vi_dap60
```

Other options: `--method dl-elmore|loess|sg`, `--weighting robust`, `--threshold`, `--min-limb-frac`, `--min-dates`, `--span`, `--sg-order`, `--sg-window`. See `--help`.

## Outputs

- `report.md`: warnings first, then plot counts by status, the noise SD, metric medians and quartiles, metric definitions and limitations.
- `tables/metrics.csv`: one row per plot. Time columns carry the axis suffix (`_dap`, `_gdd` or `_doy`):
  - `status`: `fitted`, `sparse`, `limb_missing` or `fit_failed`;
  - `peak_value` and `peak_<axis>`, with `peak_date` and the near-peak window `near_peak_start/end_<axis>`;
  - `greenup_<axis>` and `greenup_date`, `senescence_<axis>` and `senescence_date`;
  - `season_length_<days|gdd>`;
  - `greenup_rate_per_<day|gdd>` and `senescence_rate_per_<day|gdd>` (VI units per day or per GDD);
  - `auc_vi_<days|gdd>` (fitted curve) and `auc_obs_vi_<days|gdd>` (raw observations);
  - observed window `t_first/t_last_<axis>`, `max_obs_value`, `rmse`, outlier counts, `notes`, and `<vi>_at_<axis><t>` for each `--features-at` value.
  - Sparse plots have only the model-free columns filled.
- `tables/observations.csv`: each observation with its fitted value, residual, the weight used, the residual from the robust reference fit, and `outlier` (`low` / `high` / blank).
- `tables/fitted_curves.csv`: daily fitted values inside each plot's window. `tables/dl_parameters.csv`: double-logistic parameters (sos/eos are the inflection points, rsp/rau the rates per axis unit).
- `tables/date_diagnostics.csv`: median residual per flight date across plots.
- `tables/feature_time.csv`: each metric column with `available_t`, the last axis value it depends on. Pass it to `claw-yield-prediction --feature-time`.
- `figures/curves_examples.png`, `figures/curves_all.png`, `figures/metric_distributions.png`. Reproducibility bundle in `reproducibility/`.

## Expert decisions

**Axis.**
- DAP is the default because plots planted on different dates are only comparable relative to planting.
- GDD (summed from the day after planting, so GDD = 0 at planting) further removes temperature differences between sites and years.
- Day of year is allowed only as an explicit fallback, with a warning.

**Curve models.**
- `dl-beck` (default) is the 6-parameter double logistic of Beck et al. (2006).
- `dl-elmore` adds a linear summer decline (7 parameters; Elmore et al. 2012), for canopies that fade slowly before senescence.
- The tests check both formulas against `phenofit::doubleLog.Beck` and `doubleLog.Elmore`.
- Fitting uses `minpack.lm::nlsLM` (Levenberg-Marquardt; Moré 1978) on a rescaled 0-1 time axis. It tries 4 starting values. Bounds keep the curve a single hump (amplitude >= 0, senescence after green-up).
- Logistic rates are capped so a 10-90% transition is no shorter than the median sampling interval (rule of thumb: faster changes cannot be resolved). A capped rate is reported in `notes` as "at a bound".
- phenofit's own fitter was not used: on simulated crop seasons its default parameter bounds, designed for annual satellite cycles, excluded the true green-up rate.

**Minimum data.**
- A double logistic needs at least p + 1 dates (7 Beck, 8 Elmore). This is a mathematical floor, not a target: 10 or more dates spread over both sides of the peak is a better rule of thumb.
- It also needs at least 2 dates on each side of the highest value (rule of thumb). Otherwise the status is `limb_missing`.
- LOESS and SG need 6 dates (rule of thumb). Below these numbers only model-free metrics are reported (`sparse`).

**LOESS** (Cleveland & Devlin 1988): local quadratic, span 0.4, with at least 5 points per neighbourhood (rule of thumb). In simulations, span 0.5-0.75 flattened peaks and shifted senescence.

**Savitzky-Golay** (Savitzky & Golay 1964) through `signal::sgolayfilt`.
- The filter needs equal spacing, so observations are first linearly interpolated onto a regular grid (step = median sampling interval / 5).
- The window is 3 x the median sampling interval (rule of thumb, best of 3-6 x in simulation). The order is 3, the package default.
- Weights act by replacing each observation with w·y + (1-w)·ŷ before the next pass (the replacement step of Chen et al. 2004).

**Weighting.**
- `none` (default) suits clean UAV series.
- `upper`: upper envelope, the TIMESAT idea (Jönsson & Eklundh 2002). Points below the curve get bisquare weights at a scale of 6 x median |residual| (Cleveland 1979); points on or above keep weight 1. Use it for satellite data, where clouds, haze and shadows only push a VI down.
- `robust`: the same bisquare applied to both sides. Use it when errors can go either way, e.g. a mis-registered plot polygon.
- **The upper envelope trusts high values by design**: an upward error is not down-weighted and pulls the curve up. It is still flagged (see below).

**Metrics.** All are computed from the fitted curve on a fine grid within the observed window.
- Green-up and senescence are where the curve crosses base + 0.5 x amplitude, with base the minimum on that side of the peak. This is the midpoint idea of White et al. (1997), as in the TIMESAT season start/end (Jönsson & Eklundh 2004). Change the fraction with `--threshold`.
- Green-up is the first upward crossing after the pre-peak minimum. Senescence is the first downward crossing after the peak.
- A date is reported only if that side rises or falls by at least half the larger side (`--min-limb-frac 0.5`, rule of thumb). A series that ends mid-senescence gets NA, not an invented date.
- Rates are the steepest slope. They are NA when the steepest point is at the window edge, because the true maximum lies outside the window.
- AUC uses the trapezoid rule over the observed window.
- The near-peak window is where the curve is within 5% of the amplitude of its peak (rule of thumb).

**Outlier flags.**
- Each fitted plot is refitted with bisquare weights, using a noise SD pooled over all plots as the fixed scale. The pooled noise SD is the robust MAD of residuals, rescaled by n / residual degrees of freedom.
- Observations more than 3.5 noise SD from that reference curve are flagged (rule of thumb).
- The reference is the same whatever `--weighting` you chose, so an unweighted fit dragged down by a cloud does not flag its neighbors. An upward error absorbed by an upper-envelope fit is still flagged.
- A flight date is flagged when the median residual over at least 3 plots exceeds 2 noise SD (rule of thumb). This catches whole-flight calibration or illumination problems.

### Validation (tests/test_vi_timeseries.R, 76 checks)

- **Formulas and metric definitions.** The Beck and Elmore functions equal phenofit's to 1e-12. On an analytic Gaussian curve, the metric definitions reproduce the closed-form crossings, slopes, AUC and near-peak window.
- **Recovery on 20 simulated seasons** (12-16 irregular dates, noise SD 0.02):

  | Fit | Peak date | Peak value | Green-up | Senescence | AUC |
  |---|---|---|---|---|---|
  | `dl-beck` | <= 4 d | <= 0.02 | <= 2.5 d | <= 2.5 d | <= 2.5% |
  | `dl-elmore` (Elmore truth) | | | <= 2.5 d | <= 2.5 d | <= 2.5% |
  | LOESS and SG | | | <= 3.5 d | <= 3.5 d | <= 3% |

  The observed maxima for `dl-beck` were 3.1 d, 1.7 d, 1.6 d and 1.8%. Peak dates from LOESS/SG were off by up to 15 d on flat-topped curves.
- **Clouds.** Satellite-like series had 25 dates and 6 cloud drops of 0.1-0.4. Unweighted fits underestimated AUC (mean error below -5%). Upper-envelope fits had mean |AUC error| < 2% (max < 6%), with >= 80% of cloud points weighted below 0.2.
- **Upward plateau error.** Robust weighting rejects it (weight < 0.05, peak within 0.03). The upper envelope keeps it (weight > 0.99, peak biased up), as documented.
- **Flagging.** On simulated sets with 4% gross errors, flags caught 84-96% of them with 0.2-1.4% false flags.
- **Pipeline checks:**
  - a whole-flight shift is flagged on the right date;
  - DAP metrics are identical for one season planted 21 days apart;
  - GDD metrics are exactly 10 x DAP metrics at 10 GDD/day;
  - a GDD gap, a missing planting date, repeated plot x date rows and NDVI x 10000 are each rejected;
  - sparse and one-sided series are refused for the double logistic;
  - curves never leave the observed window.
- **Not validated:**
  - real-data agreement with field-scored phenology (e.g. R7 dates);
  - multi-peak seasons (regrowth, double cropping);
  - the GDD axis beyond its linear-scaling check.

## Checks before reporting

- Read **Warnings** in `report.md` first and pass each one on.
- **Stop and ask the user** if any of these appear:
  - many plots are `limb_missing` or `fit_failed`, or the peak sits at the first or last date. The planting date, axis or flight window is probably wrong, or the season was not covered.
  - the noise SD exceeds about 0.05 NDVI units, or a flight date is flagged. Suspect calibration, misregistration or wrong plot polygons before interpreting biology.
  - VI values fall outside -1..1 or look scaled.
  - `--axis doy` was used while planting dates differ.
- Flagged observations are only flagged, never removed. Show the user the list. With satellite data and many low flags, rerun with `--weighting upper` (and QA weights). With UAV data, a few isolated flags in both directions suggest `--weighting robust`.
- With a wide near-peak window (more than about 2 weeks), don't interpret the peak date. Use green-up, senescence or the window.
- NA senescence or green-up means "not reached in the observed window". Never fill it in.
- If observed windows differ between plots (warning), compare AUC and threshold dates only among plots with similar windows. Otherwise rerun with a common `--cutoff`.
- For early-season prediction, the metrics must come from a run with `--cutoff` at the prediction time. A whole-season fit leaks later data into every metric, including "green-up".
- Parameters at a bound mean the curve shape is not determined by the data. Treat those plots with caution.

## References

- Beck PSA, Atzberger C, Høgda KA, Johansen B, Skidmore AK (2006). Improved monitoring of vegetation dynamics at very high latitudes: a new method using MODIS NDVI. Remote Sensing of Environment 100(3):321-334. https://doi.org/10.1016/j.rse.2005.10.021
- Elmore AJ, Guinn SM, Minsley BJ, Richardson AD (2012). Landscape controls on the timing of spring, autumn, and growing season length in mid-Atlantic forests. Global Change Biology 18(2):656-674. https://doi.org/10.1111/j.1365-2486.2011.02521.x
- Cleveland WS (1979). Robust locally weighted regression and smoothing scatterplots. Journal of the American Statistical Association 74(368):829-836. https://doi.org/10.1080/01621459.1979.10481038
- Cleveland WS, Devlin SJ (1988). Locally weighted regression: an approach to regression analysis by local fitting. Journal of the American Statistical Association 83(403):596-610. https://doi.org/10.1080/01621459.1988.10478639
- Savitzky A, Golay MJE (1964). Smoothing and differentiation of data by simplified least squares procedures. Analytical Chemistry 36(8):1627-1639. https://doi.org/10.1021/ac60214a047
- Chen J, Jönsson P, Tamura M, Gu Z, Matsushita B, Eklundh L (2004). A simple method for reconstructing a high-quality NDVI time-series data set based on the Savitzky-Golay filter. Remote Sensing of Environment 91(3-4):332-344. https://doi.org/10.1016/j.rse.2004.03.014
- Jönsson P, Eklundh L (2002). Seasonality extraction by function fitting to time-series of satellite sensor data. IEEE Transactions on Geoscience and Remote Sensing 40(8):1824-1832. https://doi.org/10.1109/TGRS.2002.802519
- Jönsson P, Eklundh L (2004). TIMESAT - a program for analyzing time-series of satellite sensor data. Computers & Geosciences 30(8):833-845. https://doi.org/10.1016/j.cageo.2004.05.006
- White MA, Thornton PE, Running SW (1997). A continental phenology model for monitoring vegetation responses to interannual climatic variability. Global Biogeochemical Cycles 11(2):217-234. https://doi.org/10.1029/97GB00330
- Kong D, McVicar TR, Xiao M, Zhang Y, Peña-Arancibia JL, Filippa G, Xie Y, Gu X (2022). phenofit: an R package for extracting vegetation phenology from time series remote sensing. Methods in Ecology and Evolution 13(7):1508-1527. https://doi.org/10.1111/2041-210X.13870
- Moré JJ (1978). The Levenberg-Marquardt algorithm: implementation and theory. In: Numerical Analysis, Lecture Notes in Mathematics 630:105-116. https://doi.org/10.1007/BFb0067700
- Elzhov TV, Mullen KM, Spiess A-N, Bolker B. minpack.lm: R interface to the Levenberg-Marquardt nonlinear least-squares algorithm. R package. https://doi.org/10.32614/CRAN.package.minpack.lm
- signal developers. signal: signal processing. R package. https://doi.org/10.32614/CRAN.package.signal
