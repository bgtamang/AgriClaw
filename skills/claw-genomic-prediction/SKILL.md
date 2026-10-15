---
name: claw-genomic-prediction
description: Genomic prediction (genomic selection) of breeding values from SNP markers with GBLUP and RR-BLUP (rrBLUP) and optional Bayesian models BayesB, BayesC and Bayesian ridge regression (BGLR), evaluated by leakage-aware cross-validation - repeated k-fold, family-stratified folds, leave-one-family-out, and CV1/CV2 for multi-environment trials. Reports predictive ability r(y, y_hat) per fold (mean and SD), accuracy r/sqrt(h2), prediction bias (slope of observed on predicted), relatedness between test and training lines, and GEBVs for unphenotyped selection candidates. Use when the user asks for genomic prediction, genomic selection, GEBVs or GEBV ranking, GBLUP, rrBLUP, ridge regression BLUP, BayesB, BayesC, Bayesian alphabet, prediction accuracy, cross-validation of genomic prediction, training population design, predicting new families or untested lines, or sparse testing across environments.
license: MIT
compatibility: R >= 4.2 with rrBLUP; BGLR for Bayesian models; coda (optional, MCMC diagnostics), data.table (optional, fast CSV reading), openssl (optional, SHA-256 checksums). No internet needed.
metadata:
  version: "1.0.0"
  domain: genomics
  language: R
  status: tested
  summary: "Genomic prediction (GBLUP, RR-BLUP, Bayesian) with leakage-safe cross-validation"
  repository: https://github.com/bgtamang/AgriClaw
---

# Genomic prediction with leakage-aware cross-validation

## What it does

Fits genomic prediction models to phenotypes (BLUEs) and marker genotypes, estimates how well they predict lines that were not used for training by cross-validation, and predicts every genotyped line, including unphenotyped selection candidates. It reports predictive ability, accuracy, bias and how related test lines are to training lines, so inflated cross-validation is visible.

## When to use / when not to use

- Use for genomic selection in inbred or partially inbred crop populations (RIL families, breeding lines, diversity panels) with one BLUE per line, or BLUEs per line and environment.
- Genotypes should first pass `claw-marker-qc` (call rate, MAF, duplicates). Phenotype BLUEs come from `claw-trial-analysis` (one trial) or `claw-met-gxe` (several environments); the entry-mean heritability for `--h2` comes from the same analysis.
- For finding loci use `claw-gwas-crop` (panels) or `claw-qtl-mapping` (biparental and multiparent crosses). To turn accuracy into genetic gain per year use `claw-selection-index`.
- Not for G x E modelling with environmental covariates or reaction norms, multi-trait models, or non-additive kernels; say so and suggest a dedicated model (for example BGLR multi-kernel, Jarquin et al. 2014).

## Inputs

- `--geno`: CSV, first column line ID, one column per marker. Coding `0/1/2` (allele dosage), `-1/0/1` (rrBLUP) or `0/1` (two homozygous classes); `auto` detects 0/1/2 and -1/0/1 and refuses an ambiguous 0/1 file. Missing calls are mean-imputed by `rrBLUP::A.mat` and counted in the report. Markers in rows: add `--geno-transposed`.
- `--pheno`: CSV with `--id-col` (default `id`) and the `--trait` column, units in the name (`yield_kg_ha`). One row per line, or one row per line x environment with `--env-col`. Optional `--family-col` (cross, family or population) for family-aware CV.
- `--h2`: entry-mean heritability of exactly these phenotypes, for accuracy = r/sqrt(h2).

**Ask the user** (never guess):

- Are the phenotypes **BLUEs** (genotype fixed) or BLUPs? BLUPs are already shrunk toward the mean, and training on them shrinks twice; use BLUEs, or deregressed values (Garrick et al. 2009). If only BLUPs exist, say the bias (slope) estimates are unreliable.
- What **h2** applies to these values (entry-mean heritability over the reps and environments averaged)? Without it, report predictive ability only.
- What is the **prediction target**: new lines from the same families (k-fold or family-stratified), new families or crosses (leave-one-family-out), untested lines in tested environments (CV1) or sparse testing (CV2)? Choose the CV scheme that matches; report that one.
- Which column identifies **families/crosses**, and are any lines duplicates or re-sampled seed lots?
- Was any marker or line selected using the phenotypes of all lines (for example GWAS hits)? That leaks test information into training: stop and redo selection inside each fold.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/genomic_prediction.R" --demo
Rscript "${CLAUDE_SKILL_DIR}/scripts/genomic_prediction.R" --demo-multienv
Rscript "${CLAUDE_SKILL_DIR}/scripts/genomic_prediction.R" --geno geno_012.csv --pheno blues.csv \
    --trait yield_kg_ha --family-col family --cv kfold,lofo --h2 0.62 --out gp_results/
Rscript "${CLAUDE_SKILL_DIR}/scripts/genomic_prediction.R" --geno geno.csv --pheno blues.csv \
    --trait yield_kg_ha --models gblup,bayesb --n-iter 12000 --burn-in 2000
Rscript "${CLAUDE_SKILL_DIR}/scripts/genomic_prediction.R" --geno geno.csv --pheno blues_by_env.csv \
    --trait yield_kg_ha --env-col env --cv cv1,cv2 --h2 0.45
```

Run time grows with lines x markers; Bayesian models refit the sampler in every fold (25 fits for 5 x 5 CV). Try `--reps 1` first on large data.

## Outputs

- `report.md`: warnings first, data summary, CV table by scheme and model, notes on GBLUP/RR-BLUP agreement, variance components and BGLR convergence.
- `tables/cv_summary.csv`: per scheme and model `predictive_ability_mean`, `predictive_ability_sd`, `accuracy`, `slope_mean`, `slope_sd`, `median_max_relatedness`, mean training and test sizes.
- `tables/cv_folds.csv`: one row per fold (per environment for CV1/CV2): `n_train`, `n_test`, `r`, `slope`, `rmse`.
- `tables/cv_predictions.csv`: every out-of-fold prediction with observed value; `tables/cv_by_environment.csv` (multi-environment).
- `tables/predictions_all_lines.csv`: models trained on all phenotyped lines, predicting every genotyped line (`phenotyped` FALSE = selection candidates). For phenotyped lines these are fitted values, not an accuracy.
- `tables/variance_components.csv`: REML Vg, Ve and genomic heritability (Legarra 2016 scaling); `tables/bglr_convergence.csv` and `figures/bglr_trace_<model>.png` for Bayesian models.
- `figures/cv_predictive_ability.png` (fold-wise r by model and scheme), `figures/cv_observed_vs_predicted.png`.
- `reproducibility/`: `commands.sh`, `sessionInfo.txt`, checksums. The seed is in `commands.sh` and the report.

## Expert decisions

- **Markers**: recoded to rrBLUP's -1/0/1; `A.mat` removes monomorphic markers and mean-imputes missing calls; A = WW'/c with W the centred markers and c = 2 sum p(1-p) (VanRaden 2008; Endelman & Jannink 2012). Imputation is reported, never silent.
- **GBLUP** (`rrBLUP::kin.blup`, K = A) and **RR-BLUP** (`rrBLUP::mixed.solve`, Z = W) use REML variance components and are equivalent models (Endelman 2011); the report shows their agreement as a check.
- **Bayesian models** (`BGLR`; Perez & de los Campos 2014): BayesB, BayesC (Habier et al. 2011) and BRR with BGLR default priors (R2 = 0.5, df0 = 5, probIn = 0.5, counts = 10). Default 12000 iterations, 2000 burn-in, thin 5: an adjustable rule of thumb, longer than BGLR's own 1500/500 default. Convergence of the residual variance is checked with the Geweke z (flag |z| > 2) and effective sample size (flag < 100) from coda: rules of thumb.
- **Cross-validation**: 5 folds x 5 repetitions by default (adjustable rule of thumb). Family-stratified folds put about 1/k of each family in every fold. Leave-one-family-out predicts each family from the others. CV1 holds out whole lines; CV2 spreads each line's records over different folds (Burgueno et al. 2012). The same folds are used for every model.
- **Metrics**: predictive ability r(y, y_hat) per fold, mean and SD over folds. Accuracy = r / sqrt(h2), the usual indirect estimate (Estaghvirou et al. 2013): r(y, y_hat) = r(g, y_hat) x h when test-set errors are independent of the predictions, so h2 must be the heritability of the values in y (entry means). It can be biased and is sensitive to h2 (Estaghvirou et al. 2013). Bias = slope of observed on predicted, 1 = unbiased; report accuracy, bias and test-training relatedness together (Daetwyler et al. 2013). Multi-environment r is computed within each environment and averaged.
- **Leakage diagnostics**: for each fold, the median over test lines of the largest genomic correlation with any training line. Random folds keep sibs in training and inflate r relative to new families (Habier et al. 2007; Legarra et al. 2008; Riedelsheimer et al. 2013). Near-identical lines (genomic correlation > 0.95, rule of thumb) are flagged.
- **Expected accuracy for planning**: r is about sqrt(N h2 / (N h2 + Me)), with N training lines and Me independent chromosome segments (Daetwyler et al. 2008); accuracy rises with training size and relatedness.
- **Warning thresholds** (adjustable rules of thumb): fewer than 100 phenotyped lines; families under 10 lines in leave-one-family-out; |mean slope - 1| > 0.2; k-fold r more than 0.1 above leave-one-family-out; values beyond 4 SD; more than 10% missing calls or a marker more than 50% missing.

**What the tests validate**: GBLUP = RR-BLUP (identical predictions, Vg = c Vu); A = WW'/c; fold designs (balanced, family-stratified, LOFO, CV1 never shows a test line, CV2 always does); seed reproducibility; GBLUP accuracy on simulated unlinked loci between Daetwyler's first- and second-order expectations and increasing with N (0.41/0.55/0.65 for N = 120/240/480, M = 400, h2 = 0.5); r/sqrt(h2) unbiased for unrelated test lines (mean error -0.01); BRR agrees with RR-BLUP (r = 0.999); BayesB and BayesC beat GBLUP on a 5-QTL trait (0.94 and 0.93 vs 0.84); CV2 > CV1 with correlated environments; demo true accuracy higher for k-fold (0.71) than LOFO (0.40). Not validated: BGLR results with other priors, and the multi-environment model beyond main effects.

## Checks before reporting

- Read the Warnings section of `report.md` first and pass every warning on to the user.
- Report the CV scheme that matches the user's prediction target, with mean and SD of r. Never report the correlation between fitted values and phenotypes of training lines as accuracy.
- Compare random k-fold with leave-one-family-out (or CV1 with CV2): a large gap means relatives or shared environments drive the k-fold value.
- **Stop and tell the user** if: phenotypes are BLUPs; accuracy exceeds 1 (h2 too low or on the wrong basis); markers or lines were chosen using all phenotypes; duplicates or near-identical lines are spread across folds; BGLR convergence flags appear; fewer than about 100 lines (estimates are very noisy).
- A slope far from 1 means GEBVs are over- or under-dispersed: rank candidates, but do not read predicted values as expected performance.
- `median_max_relatedness` near 1 for test lines means clones or duplicates are being predicted from themselves.

## References

- Endelman JB (2011) Ridge regression and other kernels for genomic selection with R package rrBLUP. Plant Genome 4:250-255. https://doi.org/10.3835/plantgenome2011.08.0024
- Endelman JB, Jannink J-L (2012) Shrinkage estimation of the realized relationship matrix. G3 2:1405-1413. https://doi.org/10.1534/g3.112.004259
- VanRaden PM (2008) Efficient methods to compute genomic predictions. J Dairy Sci 91:4414-4423. https://doi.org/10.3168/jds.2007-0980
- Meuwissen THE, Hayes BJ, Goddard ME (2001) Prediction of total genetic value using genome-wide dense marker maps. Genetics 157:1819-1829. https://doi.org/10.1093/genetics/157.4.1819
- Perez P, de los Campos G (2014) Genome-wide regression and prediction with the BGLR statistical package. Genetics 198:483-495. https://doi.org/10.1534/genetics.114.164442
- Habier D, Fernando RL, Kizilkaya K, Garrick DJ (2011) Extension of the Bayesian alphabet for genomic selection. BMC Bioinformatics 12:186. https://doi.org/10.1186/1471-2105-12-186
- Burgueno J, de los Campos G, Weigel K, Crossa J (2012) Genomic prediction of breeding values when modeling genotype x environment interaction using pedigree and dense molecular markers. Crop Sci 52:707-719. https://doi.org/10.2135/cropsci2011.06.0299
- Daetwyler HD, Villanueva B, Woolliams JA (2008) Accuracy of predicting the genetic risk of disease using a genome-wide approach. PLoS ONE 3:e3395. https://doi.org/10.1371/journal.pone.0003395
- Daetwyler HD, Pong-Wong R, Villanueva B, Woolliams JA (2010) The impact of genetic architecture on genome-wide evaluation methods. Genetics 185:1021-1031. https://doi.org/10.1534/genetics.110.116855
- Daetwyler HD, Calus MPL, Pong-Wong R, de los Campos G, Hickey JM (2013) Genomic prediction in animals and plants: simulation of data, validation, reporting, and benchmarking. Genetics 193:347-365. https://doi.org/10.1534/genetics.112.147983
- Habier D, Fernando RL, Dekkers JCM (2007) The impact of genetic relationship information on genome-assisted breeding values. Genetics 177:2389-2397. https://doi.org/10.1534/genetics.107.081190
- Riedelsheimer C, Endelman JB, Stange M, Sorrells ME, Jannink J-L, Melchinger AE (2013) Genomic predictability of interconnected biparental maize populations. Genetics 194:493-503. https://doi.org/10.1534/genetics.113.150227
- Legarra A, Robert-Granie C, Manfredi E, Elsen J-M (2008) Performance of genomic selection in mice. Genetics 180:611-618. https://doi.org/10.1534/genetics.108.088575
- Ould Estaghvirou SB, Ogutu JO, Schulz-Streeck T, Knaak C, Ouzunova M, Gordillo A, Piepho H-P (2013) Evaluation of approaches for estimating the accuracy of genomic prediction in plant breeding. BMC Genomics 14:860. https://doi.org/10.1186/1471-2164-14-860
- Garrick DJ, Taylor JF, Fernando RL (2009) Deregressing estimated breeding values and weighting information for genomic regression analyses. Genet Sel Evol 41:55. https://doi.org/10.1186/1297-9686-41-55
- Legarra A (2016) Comparing estimates of genetic variance across different relationship models. Theor Popul Biol 107:26-30. https://doi.org/10.1016/j.tpb.2015.08.005
- Jarquin D et al. (2014) A reaction norm model for genomic selection using high-dimensional genomic and environmental data. Theor Appl Genet 127:595-607. https://doi.org/10.1007/s00122-013-2243-1
- Plummer M, Best N, Cowles K, Vines K (2006) coda: convergence diagnosis and output analysis for MCMC. R package. https://CRAN.R-project.org/package=coda
