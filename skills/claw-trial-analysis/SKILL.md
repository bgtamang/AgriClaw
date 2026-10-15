---
name: claw-trial-analysis
description: Analyse a single field trial with mixed models in statgenSTA - spatial P-spline models (SpATS) and design-based models (lme4) for RCBD, alpha-lattice / incomplete-block, row-column, augmented, p-rep and unreplicated designs. Gives genotype BLUEs and BLUPs, generalised heritability (Cullis or Oakey H2), CV, variance components, outliers, spatial-trend maps and a spatial vs non-spatial model comparison, and writes BLUEs with weights for multi-environment analysis or GWAS. Use when the user asks to analyse a yield trial, compute genotype means, BLUEs or BLUPs, heritability or repeatability, correct for field spatial trends, or prepare stage-1 means for a MET or GWAS.
license: MIT
compatibility: R >= 4.2 with statgenSTA, SpATS, lme4 and emmeans (openssl optional, for SHA-256 checksums). No internet needed.
metadata:
  version: "1.0.0"
  domain: trials
  language: R
  status: tested
  summary: "Single-trial mixed models (SpATS spatial or lme4): BLUEs, BLUPs, heritability, outliers"
  repository: https://github.com/bgtamang/AgriClaw
---

# Single-trial analysis (statgenSTA)

## What it does

For each trial in the input, builds a statgenSTA trial object (`createTD`) and fits the design model with `fitTD`: a spatial model (SpATS: the design terms plus a 2-D P-spline surface) and, by default, the classical design model (lme4). Genotype is fitted as fixed (BLUEs and their standard errors) and as random (BLUPs, genetic variance, heritability). The two models are compared, outliers are flagged with `outlierSTA`, and the BLUEs are written with weights (1/se^2) ready for `claw-met-gxe` or GWAS.

## When to use / when not to use

- Use on one or several trials with a single genotype (entry) factor: RCBD, alpha / resolvable incomplete blocks, row-column, Latin square, augmented, p-rep, diagonal-check and other unreplicated trials with field coordinates.
- Run `claw-fieldbook` first to merge and QC the data; its `clean_data.csv` is the expected input.
- Use `claw-met-gxe` for genotype x environment analysis across trials, with this skill's `met_input.csv` files.
- Don't use it for factorial treatment experiments (split-plot, strip-plot, factorial): they need an ANOVA with the right error strata. The script refuses designs labelled split, strip or factorial.

## Inputs

One row per plot (CSV). Default column names match `claw-trial-design` and `claw-fieldbook`:

| Column | Option | Meaning |
|---|---|---|
| `entry` (or `genotype`) | `--genotype` | genotype / entry name |
| traits | `--traits` | numeric, units in the name (`yield_kg_ha`) |
| `trial` | `--trial` | optional; each trial is analysed separately |
| `rep`, `block` | `--rep`, `--block` | replicate and incomplete block (block labels may restart in each replicate) |
| `row`, `col` | `--row`, `--col` | field coordinates, needed for spatial and row-column models |
| `is_check` | `--check` | optional TRUE/FALSE; used for augmented-type trials |
| `design` | | optional; written by `claw-trial-design`, used to choose the model |
| `is_filler` | | optional; FILLER plots are dropped |

**Ask the user, never guess:**

- the trait units, and the column holding the genotype if it is not `entry` or `genotype`;
- the design, when there is no `design` column and the columns are ambiguous (for example a row-column design looks like an RCBD unless `--design res.rowcol` is given);
- whether `row` and `col` are true field positions in plot units, with no gaps from renumbering;
- which entries are checks, in augmented or other unreplicated trials;
- whether the BLUEs will feed a MET or GWAS (export BLUEs) or a selection decision within this trial (use BLUPs).

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/trial_analysis.R" --demo
Rscript "${CLAUDE_SKILL_DIR}/scripts/trial_analysis.R" --input clean_data.csv --traits yield_kg_ha,plant_height_cm --out sta/
Rscript "${CLAUDE_SKILL_DIR}/scripts/trial_analysis.R" --input trial.csv --traits yield_kg_ha --genotype variety --design res.rowcol --engine both
Rscript "${CLAUDE_SKILL_DIR}/scripts/trial_analysis.R" --input trial.csv --traits yield_kg_ha --engine lme4 --design rcbd
```

The demo analyses `demo/trial_clean.csv`, a synthetic 100-entry alpha trial with a known spatial trend and a planted outlier (see `demo/README.md`).

## Outputs

- `tables/blues.csv` and `tables/blups.csv`: trial, genotype, trait, model (`SpATS` or `lme4`), estimate and standard error.
- `tables/model_summary.csv`: per trial x trait x model, the n, mean, `cv_pct`, `var_genetic`, `var_residual`, `h2` (definition below), `h2_statgenSTA`, `mean_se_blue`, residual df, outlier limit and `n_outliers`.
- `tables/model_comparison.csv`: spatial vs non-spatial H2, mean BLUE standard error, BLUE correlation, outlier counts, and which model is more precise.
- `tables/outliers.csv`: plot, genotype, position, value and standardised residual.
- `tables/met_input.csv`: one row per trial x genotype, with `<trait>` (BLUE), `se_<trait>` and `wt_<trait>` = 1/se^2. Row-bind these files across trials for `claw-met-gxe`.
- `figures/spatial_<trial>_<trait>.png`: raw data, fitted spatial trend and residuals on the field grid, with outliers ringed. `figures/blues_spatial_vs_nonspatial_*.png` compares the two models' BLUEs.
- `report.md`: warnings first, then models, results, comparison and outliers.

## Expert decisions

- **Models.** These are statgenSTA's `fitTD` models. Fixed terms: replicate (for `rcbd`, `res.ibd`, `res.rowcol`). Random terms: blocks within replicates (`res.ibd`), blocks (`ibd`), rows and columns within replicates (`res.rowcol`), or rows and columns (`rowcol`). SpATS adds `PSANOVA(col, row)` with segments = half the number of columns and rows (statgenSTA default; `--nseg` overrides). The design comes from `--design`, the `design` column (alpha -> `res.ibd`, rowcol -> `res.rowcol`, augmented -> `ibd`, p-rep / diagonal / CRD -> `rowcol`), or the columns present.
- **Spatial model by default.** SpATS models a smooth field trend without having to choose a covariance structure (Rodriguez-Alvarez et al. 2018). Spatial analysis of field trials is standard practice (Gilmour et al. 1997). The SpATS model contains the lme4 design terms, so when there is no trend it reproduces the lme4 results. In the tests (6 simulated fields, 100 genotypes x 2 replicates) H2 was 0.626 vs 0.630 without a trend; with a trend, the BLUE correlation with the true genetic values was 0.78 vs 0.68. `met_input.csv` therefore uses the SpATS BLUEs when coordinates exist (`--blues-from lme4` overrides).
- **BLUEs vs BLUPs.** Pass BLUEs (genotype fixed) with their weights to MET and GWAS. BLUPs are already shrunk towards the mean, and using them as responses shrinks twice (Smith et al. 2005; Piepho et al. 2012; Garrick et al. 2009 on de-regression). Use BLUPs to rank or select genotypes within the trial.
- **Heritability (`h2`).** SpATS: Oakey et al. (2006) generalised heritability, the effective dimension of the genotype term over its nominal dimension (`SpATS::getHeritability`, reported unrounded). lme4: Cullis et al. (2006) generalised heritability, 1 - mean PEV of BLUP differences / (2 sigma2_g), computed from lme4's full mixed-model equations; for balanced RCBD it equals sigma2_g / (sigma2_g + sigma2_e / r) exactly (tested). `h2_statgenSTA` for lme4 is statgenSTA's own sigma2_g / (sigma2_g + sigma2_e / r), or the plot-level ratio when there is no replicate column. It ignores block adjustments and imbalance, so prefer `h2` (see also Piepho & Mohring 2007).
- **Comparing models.** The comparison uses precision (mean standard error of BLUEs) and generalised H2, plus the correlation of the two sets of BLUEs. AIC and BIC are not compared, because SpATS and lme4 likelihoods are not on the same footing.
- **Outliers.** `outlierSTA` flags |standardised residual| > `qnorm(1 - 0.5 / rdf)`, kept between 2 and 4 (statgenSTA default; `--outlier-limit` overrides). Flags are reported, never removed. With 2 replicates, one bad plot drags its genotype's BLUE, so its twin plot is flagged with the opposite sign; check the raw values to see which plot is wrong.
- **CV** = 100 x sqrt(residual variance) / mean (statgenSTA). The SpATS residual excludes the trend, so its CV is lower than the lme4 or ANOVA CV; don't compare CVs across engines.
- **Unreplicated trials** need replicated checks or a spatial model: without either, the script stops. When checks are more replicated than test entries, the SpATS BLUP model fits them as fixed `checkId` effects (statgenSTA convention), so BLUPs and H2 refer to the test lines. statgenSTA's lme4 engine drops check plots under `checkId`, so with lme4 the checks stay random genotypes; the report warns. BLUEs are unaffected.
- **Coordinates** must be complete and unique for spatial and row-column models (the script stops otherwise). For block designs fitted only with `--engine lme4`, incomplete coordinates are ignored with a warning.
- **Rules of thumb** that trigger warnings: H2 < 0.2 (genotypes poorly separated), CV > 20% for yield, and zero genetic variance.

## Checks before reporting

- Always report CV, H2 and the number of outliers for every trait and trial. Pass on every warning.
- Look at the residual panel of the spatial figure: a pattern left there (stripes along rows or columns) means structure the model missed. Consider `--nseg` or a row-column design model.
- Zero genetic variance, or H2 near 0, means the trial cannot separate genotypes. Stop and tell the user before passing its BLUEs on.
- A large H2 gain from the spatial model means a strong field trend. Mention it, because the non-spatial results would be misleading.
- Outliers: ask the user to check field notes; refit only if they decide a value is wrong and remove it in `claw-fieldbook`.
- Confirm that BLUEs, not BLUPs, go to `claw-met-gxe` and GWAS.

## References

- van Rossum B (2025). statgenSTA: Single Trial Analysis (STA) of Field Trials. R package version 1.0.15. https://doi.org/10.32614/CRAN.package.statgenSTA
- Rodriguez-Alvarez MX, Boer MP, van Eeuwijk FA, Eilers PHC (2018). Correcting for spatial heterogeneity in plant breeding experiments with P-splines. Spatial Statistics 23:52-71. https://doi.org/10.1016/j.spasta.2017.10.003
- Bates D, Maechler M, Bolker B, Walker S (2015). Fitting linear mixed-effects models using lme4. Journal of Statistical Software 67(1):1-48. https://doi.org/10.18637/jss.v067.i01
- Cullis BR, Smith AB, Coombes NE (2006). On the design of early generation variety trials with correlated data. Journal of Agricultural, Biological, and Environmental Statistics 11(4):381-393. https://doi.org/10.1198/108571106X154443
- Oakey H, Verbyla A, Pitchford W, Cullis B, Kuchel H (2006). Joint modeling of additive and non-additive genetic line effects in single field trials. Theoretical and Applied Genetics 113(5):809-819. https://doi.org/10.1007/s00122-006-0333-z
- Piepho H-P, Mohring J (2007). Computing heritability and selection response from unbalanced plant breeding trials. Genetics 177(3):1881-1888. https://doi.org/10.1534/genetics.107.074229
- Gilmour AR, Cullis BR, Verbyla AP (1997). Accounting for natural and extraneous variation in the analysis of field experiments. Journal of Agricultural, Biological, and Environmental Statistics 2(3):269-293. https://doi.org/10.2307/1400446
- Smith AB, Cullis BR, Thompson R (2005). The analysis of crop cultivar breeding and evaluation trials: an overview of current mixed model approaches. Journal of Agricultural Science 143(6):449-462. https://doi.org/10.1017/S0021859605005587
- Piepho H-P, Mohring J, Schulz-Streeck T, Ogutu JO (2012). A stage-wise approach for the analysis of multi-environment trials. Biometrical Journal 54(6):844-860. https://doi.org/10.1002/bimj.201100219
- Garrick DJ, Taylor JF, Fernando RL (2009). Deregressing estimated breeding values and weighting information for genomic regression analyses. Genetics Selection Evolution 41:55. https://doi.org/10.1186/1297-9686-41-55
