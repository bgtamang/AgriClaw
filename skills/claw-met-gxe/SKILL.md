---
name: claw-met-gxe
description: Analyse multi-environment trials (MET) from per-environment genotype BLUEs with statgenGxE - genotype-by-environment (GxE) variance components and heritability across environments, MET genotype BLUPs (weighted by stage-1 standard errors when given), AMMI and GGE biplots (which-won-where), Finlay-Wilkinson regression, stability (Shukla stability variance, Wricke ecovalence, Lin-Binns superiority, static stability) and mega-environments. Use when the user has the same genotypes in several trials, locations or years and asks about GxE, adaptation, yield stability, the best genotype across or within environments, AMMI or GGE biplots, mega-environments, or heritability over environments.
license: MIT
compatibility: R >= 4.2 with statgenGxE, statgenSTA, lme4, emmeans and nlme (agricolae only for the tests; openssl optional). No internet needed.
metadata:
  version: "1.0.0"
  domain: trials
  language: R
  status: tested
  repository: https://github.com/bgtamang/AgriClaw
---

# Multi-environment trials: GxE, stability and biplots

## What it does

Stage 2 of a two-stage MET analysis. It takes one BLUE per genotype per environment from stage 1 (`claw-trial-analysis`) and estimates genotype and GxE variance and heritability across environments, predicts genotype performance over the environments (BLUPs), and describes GxE. The GxE descriptions are AMMI and GGE biplots, Finlay-Wilkinson sensitivity, stability measures and mega-environments. When the BLUEs carry standard errors, a stage-2 model with known stage-1 variances weights each BLUE by its precision.

## When to use / when not to use

- Use when the same genotypes were tested in two or more trials (locations, years, managements). AMMI, GGE, stability and mega-environments need at least 3.
- Run `claw-trial-analysis` on each trial first and stack its `tables/met_input.csv` files. Give this skill BLUEs, never BLUPs: shrinking in both stages ("double shrinkage") distorts variances and rankings.
- Don't use it on plot-level data: the input must have one value per genotype per environment, and the script refuses duplicates.
- Don't use it for a single trial (`claw-trial-analysis`), for GWAS (`claw-gwas-crop`) or for genomic prediction across environments (`claw-genomic-prediction`).

## Inputs

One or more CSVs (`--input a.csv,b.csv` are stacked), one row per genotype per environment:

| Column | Meaning |
|---|---|
| `trial` (`--env`) | environment ID (trial, location-year) |
| `genotype` (`--genotype`) | genotype ID, spelled identically in every trial |
| `<trait>` (`--trait`) | per-environment BLUE, units in the name (`yield_kg_ha`) |
| `wt_<trait>` or `se_<trait>` | optional: weight 1/SE^2 or the SE of each BLUE (`--weight`, `--no-weights`) |
| location, year, scenario | optional: `--location` + `--year` (each trial one location-year) or `--nesting` (trials nested in a region or scenario) |

**Ask the user, never guess:**

- the trait units, and whether higher or lower is better (`--best max|min`; use `min` for days to maturity or disease score);
- how environments are structured: locations x years, or trials nested in regions or management scenarios. This decides `--location/--year` or `--nesting`; leave both out if unsure;
- whether genotype names are harmonised across trials (the same line spelled two ways becomes two genotypes);
- whether the trials represent the target population of environments, because heritability and BLUPs refer only to the environments sampled;
- whether stage 1 used the same model in every trial, so that the SEs are comparable.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/met_gxe.R" --demo
Rscript "${CLAUDE_SKILL_DIR}/scripts/met_gxe.R" --input ames/tables/met_input.csv,lewis/tables/met_input.csv,nashua/tables/met_input.csv --trait yield_kg_ha --out met_yield/
Rscript "${CLAUDE_SKILL_DIR}/scripts/met_gxe.R" --input met_all.csv --trait days_to_maturity --best min --location location --year year
Rscript "${CLAUDE_SKILL_DIR}/scripts/met_gxe.R" --input met_all.csv --trait yield_kg_ha --nesting region --complete-only --npc 2
```

The demo analyses `demo/met_input.csv`: 40 synthetic genotypes in 8 trials with SEs and 3 missing cells; see `demo/README.md`. Expect about 20-30 s, mostly package loading. `gxeFw` slows down with many genotypes (about 1 min for 250 genotypes x 10 trials).

## Outputs

- `tables/genotype_blups.csv`: genotype, number of environments, raw mean, `blup` (statgenGxE, unweighted), `blup_weighted` (when SEs are given), rank. Ranking uses the weighted BLUP when available.
- `tables/variance_components.csv`, `tables/heritability.csv`, `tables/environment_summary.csv` (mean, SD and mean SE per environment).
- `tables/ammi_anova.csv`, `tables/pc_importance.csv` (% variance per PC, AMMI and GGE), `tables/ammi_scores.csv`, `tables/gge_scores.csv` (symmetric scaling).
- `tables/finlay_wilkinson.csv` (sensitivity, SE, mean, deviation MS), `tables/fw_anova.csv`, `tables/fw_environments.csv`.
- `tables/stability.csv`: mean, superiority, Wricke's ecovalence, Shukla's stability variance, static stability, and ranks.
- `tables/mega_environments.csv`: trials, mega-environment, winning genotype.
- `figures/ammi2_biplot.png`, `figures/gge_biplot.png` (which-won-where polygon), `figures/finlay_wilkinson.png`, `figures/stability_vs_mean.png`. Biplot axes show the % of variance explained.
- `report.md`: warnings first. Then the data and missing cells, variance components and H2, the weighted model, top genotypes, AMMI, GGE, FW, stability and mega-environments.

## Expert decisions

- **Two-stage analysis.** Stage 1 estimates genotype BLUEs within each trial; stage 2 combines them (Smith et al. 2005; Piepho et al. 2012). BLUEs go in, never BLUPs.
- **Variance components and BLUPs** (`statgenGxE::gxeVarComp`, lme4; Bates et al. 2015). Model: trial (fixed) + genotype (random) + residual. With one BLUE per genotype x trial, the GxE term is confounded with the residual, so statgenGxE drops it and the **residual = GxE + stage-1 error**. `--nesting` adds genotype x nesting group, and `--location/--year` adds genotype x location and genotype x year (statgenGxE's documented models). The location x year model needs one trial per location-year (the factorial structure statgenGxE documents); with several trials per location-year (for example two managements), use `--nesting`. Genotype BLUPs come from `predict(predictLevel = "genotype")`, the mean trial effect plus the genotype BLUP. emmeans degrees of freedom are set to asymptotic, which gives identical values faster.
- **Heritability** (`statgenGxE::herit`): H2 = sigma2_G / (sigma2_G + sum(sigma2_term / n_term) + sigma2_residual / (r n)). This is the entry-mean heritability of a genotype's mean over the environments here (Atlin et al. 2000); for unbalanced data the numbers of levels are harmonic means (Holland et al. 2003). Tested: within 0.03 of the true 0.914.
- **Weights.** statgenGxE's `useWt` passes weights to lme4, where they are relative: residual variance sigma2 / w, with sigma2 estimated. The residual then lives on the weight scale, and `herit()` reads it as a variance. In a simulation with true H2 = 0.89, it reported 0.54-0.69, while BLUP accuracy barely changed. So statgenGxE runs **unweighted**. When SEs are available, a second stage-2 model fixes each BLUE's residual variance at its SE^2 (`nlme::lme` with `varFixed` and `sigma = 1`, as in Piepho et al. 2012) and estimates the GxE variance separately. Tested: sigma2_GE recovered within 0.1 of 0.5, and weighted BLUPs were at least as accurate as unweighted ones when trial precision differs. Only the trial model is weighted. AMMI, GGE, FW and stability assume equal precision and stay unweighted.
- **AMMI and GGE** (`gxeAmmi`, `gxeGGE`; Gauch 1988, 2006; Yan et al. 2000). `--npc 2` by default, capped below the number of environments. % variance is of the GxE sum of squares (AMMI) or the G + GxE sum of squares (GGE), shown on the biplot axes. **Missing cells**: statgenGxE drops all-missing genotypes, stops if more than 30% of genotype x environment cells are empty (reported as skipped), and otherwise fills empty cells by iterative regression (`multMissing`) before the PCA. The report gives the count and %. `--complete-only` drops incomplete genotypes instead. Biplots use symmetric scaling. The GGE polygon and its perpendicular rays give the which-won-where view (Yan & Tinker 2006). Within-environment genotype rankings from inner products do not depend on the scaling, but the 2-D view is only as good as the PC1 + PC2 share.
- **Finlay-Wilkinson** (`gxeFw`; Finlay & Wilkinson 1963; deviations as in Eberhart & Russell 1966). Each genotype is regressed on the estimated environmental index. Slopes average 1; slope > 1 means a stronger response to good environments. Tested: correlation with the true slopes > 0.9 (0.975 to 0.978).
- **Stability.** It is computed only for genotypes observed in **every** environment, because `gxeStability` still divides by the full number of trials when a cell is absent. Measures:
  - Lin & Binns (1988) superiority: distance to the best genotype in each environment, combining level and consistency; low = good.
  - Wricke's ecovalence (statgenGxE `wricke`).
  - Shukla's (1972) stability variance: sigma2_i = p W_i / ((p-2)(q-1)) - sum(W) / ((p-1)(p-2)(q-1)), computed from the ecovalence W. It equals `agricolae::stability.par` "Sigma-square" divided by its replicate count, and its mean equals the GxE mean square (both tested).
  - statgenGxE's `static` is each genotype's **variance across environments**, i.e. static stability (Becker & Leon 1988). statgenGxE's help page calls it Shukla's stability variance, but it is not (source and test).
- **Stability is not performance.** Static stability favours unresponsive, often low-yielding genotypes, so always read stability next to the mean (the figure plots both).
- **Mega-environments** (`gxeMegaEnv`): trials are grouped by their AMMI2 winning genotype. A grouping is a mega-environment only if it repeats across years (Yan et al. 2000); the script warns when there are fewer than 8 trials. Tested: two simulated groups and their winners recovered.
- **Number of environments.** Rule of thumb: with fewer than about 5 environments, GxE variance, stability and biplots are weakly estimated; the report warns. Fewer than 3 skips AMMI, GGE, stability and mega-environments.

## Checks before reporting

- Read the Warnings in `report.md` first and pass every one on.
- Zero or near-zero variance components (singular fit) make H2 unreliable; say so.
- Many imputed cells (above about 10%, a rule of thumb) make AMMI and GGE partly model-made. Re-run with `--complete-only` and compare.
- If PC1 + PC2 explain little of GGE (rule of thumb: under about 60%), don't read which-won-where from the biplot; report the share.
- A low correlation between weighted and unweighted BLUPs (below about 0.9) means very unequal trial precision or odd SEs. Check the `mean_se` column of `environment_summary.csv`.
- Genotypes tested in one environment have BLUPs, sensitivities and scores resting on one trial; flag them.
- Never call a genotype "stable" without its mean, and never call a one-year grouping a mega-environment.
- Stop and ask if genotype names do not match across trials, or if the input looks like plot data. The script refuses duplicate genotype x environment rows.

## References

- van Rossum B (2025). statgenGxE: Genotype by Environment (GxE) Analysis. R package version 1.0.11. https://doi.org/10.32614/CRAN.package.statgenGxE
- Finlay KW, Wilkinson GN (1963). The analysis of adaptation in a plant-breeding programme. Australian Journal of Agricultural Research 14(6):742-754. https://doi.org/10.1071/AR9630742
- Eberhart SA, Russell WA (1966). Stability parameters for comparing varieties. Crop Science 6(1):36-40. https://doi.org/10.2135/cropsci1966.0011183X000600010011x
- Shukla GK (1972). Some statistical aspects of partitioning genotype-environmental components of variability. Heredity 29(2):237-245. https://doi.org/10.1038/hdy.1972.87
- Lin CS, Binns MR (1988). A superiority measure of cultivar performance for cultivar x location data. Canadian Journal of Plant Science 68(1):193-198. https://doi.org/10.4141/cjps88-018
- Becker HC, Leon J (1988). Stability analysis in plant breeding. Plant Breeding 101(1):1-23. https://doi.org/10.1111/j.1439-0523.1988.tb00261.x
- Gauch HG (1988). Model selection and validation for yield trials with interaction. Biometrics 44(3):705-715. https://doi.org/10.2307/2531585
- Gauch HG (2006). Statistical analysis of yield trials by AMMI and GGE. Crop Science 46(4):1488-1500. https://doi.org/10.2135/cropsci2005.07-0193
- Yan W, Hunt LA, Sheng Q, Szlavnics Z (2000). Cultivar evaluation and mega-environment investigation based on the GGE biplot. Crop Science 40(3):597-605. https://doi.org/10.2135/cropsci2000.403597x
- Yan W, Tinker NA (2006). Biplot analysis of multi-environment trial data: Principles and applications. Canadian Journal of Plant Science 86(3):623-645. https://doi.org/10.4141/P05-169
- Atlin GN, Baker RJ, McRae KB, Lu X (2000). Selection response in subdivided target regions. Crop Science 40(1):7-13. https://doi.org/10.2135/cropsci2000.4017
- Holland JB, Nyquist WE, Cervantes-Martinez CT (2003). Estimating and interpreting heritability for plant breeding: an update. Plant Breeding Reviews 22:9-112. https://doi.org/10.1002/9780470650202.ch2
- Smith AB, Cullis BR, Thompson R (2005). The analysis of crop cultivar breeding and evaluation trials: an overview of current mixed model approaches. Journal of Agricultural Science 143(6):449-462. https://doi.org/10.1017/S0021859605005587
- Piepho H-P, Mohring J, Schulz-Streeck T, Ogutu JO (2012). A stage-wise approach for the analysis of multi-environment trials. Biometrical Journal 54(6):844-860. https://doi.org/10.1002/bimj.201100219
- Bates D, Machler M, Bolker B, Walker S (2015). Fitting linear mixed-effects models using lme4. Journal of Statistical Software 67(1):1-48. https://doi.org/10.18637/jss.v067.i01
- Pinheiro JC, Bates DM (2000). Mixed-Effects Models in S and S-PLUS. Springer, New York (R package nlme). https://doi.org/10.1007/b98882
