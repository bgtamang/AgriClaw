---
name: claw-yield-prediction
description: Predict crop yield from remote-sensing features - NDVI/VI time-series metrics (from claw-vi-timeseries), per-date vegetation indices, canopy cover, weather summaries - with elastic net (glmnet) and random forest (ranger), and measure prediction accuracy with nested, leakage-safe cross-validation. Schemes are leave-one-environment-out, leave-genotype-out, leave-year-out and new genotypes in a new environment; random k-fold is shown only for comparison. Reports r, RMSE, relative RMSE, bias, both R2 definitions (squared correlation and 1 - SSE/SST), within-environment accuracy and permutation importance. Flags inflated accuracy, temporal leakage (late-season features in an early forecast) and models worse than the mean. Use when the user asks to predict or forecast yield from drone, UAV or satellite data, build a phenomic or machine-learning yield model, cross-validate one, or check whether a reported prediction accuracy is inflated by leakage.
license: MIT
compatibility: R >= 4.2 with glmnet and ranger (openssl optional, for SHA-256 checksums). Base R graphics. No internet needed.
metadata:
  version: "1.0.0"
  domain: remote-sensing
  language: R
  status: tested
  repository: https://github.com/bgtamang/AgriClaw
---

# Yield prediction with leakage-safe cross-validation

## What it does

Fits elastic-net and random-forest models that predict yield from per-plot features. It estimates how well they predict new data with nested cross-validation (CV): hyperparameters are tuned inside each training set, and every scheme holds out whole groups (environments, genotypes, years, or both). This mimics the real prediction task. Random k-fold CV is reported beside the grouped schemes so the user can see how much it overstates accuracy.

## When to use / when not to use

- Use to predict yield, or another plot trait, from remote-sensing features: phenology metrics, VIs at fixed days after planting, canopy cover, weather totals.
- Use to audit a reported accuracy. Shared genotypes or environments between training and test sets, or late-season features in an "early" model, inflate it.
- Get features with `claw-vi-extract` and `claw-vi-timeseries` first; that skill's `tables/feature_time.csv` plugs into `--feature-time`. Convert units with `claw-yield-units`.
- Don't use it for marker-based genomic prediction (`claw-genomic-prediction`), for estimating genotype BLUEs from one trial (`claw-trial-analysis`), or for G x E and stability analysis (`claw-met-gxe`).
- Don't treat predicted yields as phenotypes for GWAS or selection. They are shrunk towards the mean and inherit the features' genetics, not yield's.

## Inputs

A CSV with one row per plot, or per genotype x environment mean or BLUE:

| Column | Meaning |
|---|---|
| `yield_kg_ha` (`--target`) | the target, **units in the name** (kg/ha, bu/ac, ...) |
| `plot_id` (`--id-col`) | optional row identifier (must be unique) |
| `env` (`--env-col`) | environment = site x year (or trial) |
| `genotype` (`--genotype-col`) | entry / line / hybrid |
| `year` (`--year-col`) | season, if several |
| features | numeric columns (`--features a,b,c`, or the features listed in `--feature-time`) |

`--feature-time` takes a CSV with `feature` and `available_t`: the time (e.g. DAP) from which each feature is known.

**Ask the user (don't guess):**
- **What is the prediction for?** This picks the headline scheme (`--primary`):

  | Goal | Scheme |
  |---|---|
  | new lines in environments already in the data | `leave-genotype-out` |
  | known lines in a new environment | `leave-env-out` |
  | next season | `leave-year-out` |
  | new lines in a new environment | `leave-genotype-env-out` |

- **Is it an early-season forecast?** If so, at what DAP, and when is each feature known? Then use `--feature-time` with `--cutoff`. Features for an early forecast must come from data up to that date (`claw-vi-timeseries --cutoff`), not from a whole-season fit.
- **Yield units and type.** Plot yield or BLUE? Moisture-adjusted? For leave-env-out, BLUEs must be per environment: BLUEs pooled across environments already carry information from the held-out one.
- **What defines an environment**, and is there a year column?
- **Feature provenance.** Did any feature come from harvest, from yield itself, from field position (row/range), or from processing that used all plots at once (e.g. normalising over the whole data set)? Any of these leaks.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/yield_prediction.R" --demo
Rscript "${CLAUDE_SKILL_DIR}/scripts/yield_prediction.R" --input plots.csv --target yield_kg_ha --env-col env --genotype-col genotype --out results/yield_cv
# features from claw-vi-timeseries, early forecast at DAP 60
Rscript "${CLAUDE_SKILL_DIR}/scripts/yield_prediction.R" --input plots.csv --feature-time vi_dap60/tables/feature_time.csv --cutoff 60 --primary leave-env-out
# faster run: fewer trees, two schemes
Rscript "${CLAUDE_SKILL_DIR}/scripts/yield_prediction.R" --input plots.csv --features ndvi_peak,auc_vi_days,senescence_dap --cv random,leave-genotype-out --num-trees 200
```

Other options: `--models`, `--folds`, `--inner-folds`, `--seed`, `--alpha-grid`, `--enet-lambda`, `--mtry-grid`, `--min-node-grid`, `--perm-repeats`, `--drop-na-rows`. See `--help`. A few hundred rows with all 5 schemes takes a few minutes. The leave-genotype-env-out scheme costs the most (environments x genotype folds).

## Outputs

- `report.md`: warnings first, then the features used and excluded, the accuracy table, top permutation importance and limitations.
- `tables/cv_metrics.csv`: one row per scheme x model, plus the `train_mean` baseline. Columns:
  - `r`, `rmse`, `rrmse_pct`, `bias` (mean predicted − observed);
  - `r2_cor` (squared correlation) and `r2_sse` (1 − SSE/SST);
  - `r_within_env` (mean of within-environment correlations);
  - `flag` ("worse than predicting the mean").
- `tables/cv_predictions.csv`: every out-of-fold prediction with ids, scheme, model and fold. `tables/fold_metrics.csv`: the same metrics per held-out fold, e.g. per environment.
- `tables/tuning.csv`: the chosen alpha/lambda or mtry/min.node.size per fold, and how the inner folds were formed.
- `tables/permutation_importance.csv`: headline scheme only; mean held-out RMSE increase, SD over folds, rank. `tables/features.csv`: each candidate feature, used or not, and why.
- `figures/accuracy_by_scheme.png`, `figures/obs_vs_pred.png` (headline scheme, environment centroids labelled), `figures/importance.png`. Reproducibility bundle in `reproducibility/`.

## Expert decisions

**Elastic net** (Zou & Hastie 2005) through `glmnet::cv.glmnet` (Friedman et al. 2010).
- Alpha is searched over 0, 0.25, 0.5, 0.75 and 1 (rule of thumb: ridge to lasso).
- Lambda is chosen by inner CV with the one-standard-error rule. That is glmnet's default for `predict()`, and it favours the sparser model when n is small. Change it with `--enet-lambda min`.

**Random forest** (Breiman 2001) through `ranger` (Wright & Ziegler 2017).
- 500 trees, ranger's default.
- `mtry` is searched over sqrt(p), p/3 and 2p/3 (rule of thumb spanning ranger's default sqrt(p) and the classic regression default p/3).
- `min.node.size` is 5, ranger's regression default.
- A seeded ranger gives identical forests at any thread count, so fits use all cores.

**Nested CV** (Varma & Simon 2006). Tuning error estimated on the same folds as the reported error is optimistic, so every hyperparameter is tuned by inner CV inside the outer training rows only.
- Inner folds are grouped like the outer scheme (environments, years or genotypes) when the training rows hold at least 3 groups. Otherwise they are random, and `tuning.csv` and the report say so.

**Grouped CV schemes** (Roberts et al. 2017).
- Random k-fold puts rows of the same genotype and environment on both sides of the split. A model can then score well by recognising them, through replicated genotypes or environment-specific calibration offsets, without having learnt anything that transfers.
- `leave-genotype-out` is the "untested genotypes" problem (CV1 of Burgueño et al. 2012).
- `leave-genotype-env-out` holds out one environment and one genotype fold together. The training set excludes both, matching the untested-sites-and-lines scenarios of Jarquín et al. (2017).
- Folds: `--folds` (5, rule of thumb) for random and genotype folds; one environment or year per fold otherwise.

**No preprocessing on the full data.**
- glmnet standardises features inside every fit from that fit's own training rows (`standardize = TRUE`). Forests need no scaling.
- Nothing is imputed: rows with missing features stop the run unless `--drop-na-rows`, which drops and counts them.

**Metrics.**
- Pooled over all out-of-fold predictions. R² is reported both ways (Kvålseth 1985):
  - `r2_cor` ignores bias and scale and can be high for badly biased predictions;
  - `r2_sse` = 1 − SSE/SST is the honest test-set measure (Alexander et al. 2015), and is negative when the model is worse than the mean of the observations.
- rRMSE = 100 × RMSE / mean observed (NA if the mean is not positive).
- `train_mean` is the baseline every model must beat.
- `r_within_env` averages correlations inside each environment. Pooled r across environments mixes in how well environment means are predicted.

**Permutation importance** (Breiman 2001).
- Each feature is shuffled within each held-out fold of the headline scheme (`--perm-repeats` 5), and the increase in test RMSE is recorded.
- Correlated features share or hide importance (Strobl et al. 2008), and shuffling creates unrealistic rows (Hooker et al. 2021).
- A feature that is constant within the held-out group, e.g. season weather under leave-env-out, always scores 0. Its effect appears only as the accuracy gap between schemes.

**Temporal leakage guard.** With `--cutoff`, every candidate feature needs an `available_t`, and features known later are dropped and listed in `features.csv`.

**Warnings** (all rules of thumb):
- fewer than 50 rows;
- more features than 0.2 × rows;
- a grouped scheme with fewer than 4 held-out environments or years;
- a feature with |r| > 0.95 with the target;
- features constant within environments;
- feature names like row/range/rep/block;
- target name without units;
- random-vs-headline gap in r above 0.1;
- pooled r more than 0.3 above within-environment r;
- negative R² (1 − SSE/SST).

### Validation (tests/test_yield_prediction.R, 75 checks)

- **Metric formulas:** checked against hand calculations.
- **Folds:** every row is tested exactly once, and no environment, genotype or year is shared between training and test.
- **No information from held-out rows:** changing the held-out rows' features and yield leaves that fold's elastic-net coefficients and random-forest predictions identical.
- **Environment leakage** (simulated fingerprints, no transferable signal): the random forest scores r 0.90 (R² 0.80) under random CV but r 0.20 (R² −0.07) under leave-env-out. The test requires gaps of >= 0.5 in both r and R².
- **Genotype leakage** (replicated genotypes, identifying traits): random r 0.67 versus leave-genotype-out r −0.13 (R² −0.36).
- **Real linear signal** (oracle r 0.848):
  - elastic net r 0.843 under leave-env-out and 0.829 under leave-genotype-env-out (tolerance: r within 0.05 of the oracle, R² within 0.1 of oracle r²);
  - random forest R² 0.57 against an oracle of 0.72 (tolerance 0.25);
  - permutation importance ranks the two informative features first, for both models.
- **Pure-noise target:** every R² <= 0 (max −0.02), and both models are flagged.
- **Other checks:** same seed gives identical predictions; the `--cutoff` guard and input errors behave as documented; the demo runs end to end.
- **Not validated:**
  - accuracy on real field data (no public benchmark is bundled);
  - hyperparameter grids other than the defaults;
  - categorical or classification targets (not supported);
  - repeated random CV (one fold assignment per seed).

## Checks before reporting

- Read **Warnings** first and pass each one on.
- Report the scheme that matches the user's goal. **Never headline random k-fold.** Show it only to quantify how much leakage inflates accuracy.
- **Stop and tell the user** in any of these cases:
  - R² (1 − SSE/SST) is negative for the headline scheme. Say plainly that the model is worse than the mean, and that its predictions must not be used as phenotypes (GWAS, selection, covariates).
  - Random CV is far above the grouped schemes. Explain that the high number comes from shared genotypes or environments, not from transferable signal.
  - A feature looks like yield, a harvest measurement, or a field position.
  - An "early" model uses features whose `available_t` is after the forecast date, or features from a whole-season curve fit.
- Compare pooled r with `r_within_env`:
  - pooled much higher means the model mostly predicts environment means;
  - pooled negative with positive within-environment r means it ranks genotypes but cannot place a new environment.
- With fewer than about 4 environments or years, call grouped estimates anecdotal; one fold more or less can swing them.
- Read `fold_metrics.csv`: one bad environment can dominate the pooled number.
- Importance means "the model relies on this feature", not "this feature causes yield".

## References

- Zou H, Hastie T (2005). Regularization and variable selection via the elastic net. Journal of the Royal Statistical Society B 67(2):301-320. https://doi.org/10.1111/j.1467-9868.2005.00503.x
- Friedman J, Hastie T, Tibshirani R (2010). Regularization paths for generalized linear models via coordinate descent. Journal of Statistical Software 33(1):1-22. https://doi.org/10.18637/jss.v033.i01
- Breiman L (2001). Random forests. Machine Learning 45(1):5-32. https://doi.org/10.1023/A:1010933404324
- Wright MN, Ziegler A (2017). ranger: a fast implementation of random forests for high dimensional data in C++ and R. Journal of Statistical Software 77(1):1-17. https://doi.org/10.18637/jss.v077.i01
- Varma S, Simon R (2006). Bias in error estimation when using cross-validation for model selection. BMC Bioinformatics 7:91. https://doi.org/10.1186/1471-2105-7-91
- Roberts DR, Bahn V, Ciuti S, et al. (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. Ecography 40(8):913-929. https://doi.org/10.1111/ecog.02881
- Kaufman S, Rosset S, Perlich C, Stitelman O (2012). Leakage in data mining: formulation, detection, and avoidance. ACM Transactions on Knowledge Discovery from Data 6(4):15. https://doi.org/10.1145/2382577.2382579
- Burgueño J, de los Campos G, Weigel K, Crossa J (2012). Genomic prediction of breeding values when modeling genotype × environment interaction using pedigree and dense molecular markers. Crop Science 52(2):707-719. https://doi.org/10.2135/cropsci2011.06.0299
- Jarquín D, Lemes da Silva C, Gaynor RC, et al. (2017). Increasing genomic-enabled prediction accuracy by modeling genotype × environment interactions in Kansas wheat. The Plant Genome 10(2). https://doi.org/10.3835/plantgenome2016.12.0130
- Kvålseth TO (1985). Cautionary note about R². The American Statistician 39(4):279-285. https://doi.org/10.1080/00031305.1985.10479448
- Alexander DLJ, Tropsha A, Winkler DA (2015). Beware of R²: simple, unambiguous assessment of the prediction accuracy of QSAR and QSPR models. Journal of Chemical Information and Modeling 55(7):1316-1322. https://doi.org/10.1021/acs.jcim.5b00206
- Strobl C, Boulesteix A-L, Kneib T, Augustin T, Zeileis A (2008). Conditional variable importance for random forests. BMC Bioinformatics 9:307. https://doi.org/10.1186/1471-2105-9-307
- Hooker G, Mentch L, Zhou S (2021). Unrestricted permutation forces extrapolation: variable importance requires at least one more model, or there is no free variable importance. Statistics and Computing 31:82. https://doi.org/10.1007/s11222-021-10057-z
