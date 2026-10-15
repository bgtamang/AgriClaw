---
name: claw-qtl-mapping
description: Map quantitative trait loci (QTL) in biparental and multiparent crop populations with R/qtl2 - recombinant inbred lines (RIL), F2, backcross, doubled haploids, and MAGIC or multi-parent RILs from a qtl2 control file. Computes genotype probabilities with a genotyping-error hidden Markov model, runs Haley-Knott regression or a LOCO linear mixed model genome scan, sets genome-wide LOD thresholds by permutation (5% and 10%, separate X thresholds), finds peaks with 1.5-LOD support or Bayes credible intervals (also in Mb when a physical map is given), estimates additive and dominance or founder allele effects with standard errors and percent variance explained, and reports marker QC (missing data, segregation distortion, duplicate lines). Use when the user asks for QTL mapping, interval mapping, a LOD scan or LOD profile, QTL support intervals, permutation thresholds, linkage mapping in a RIL, DH or F2 population, MAGIC or multiparent QTL analysis, or finding loci for a trait in a biparental cross.
license: MIT
compatibility: R >= 4.2 with qtl2; openssl optional (SHA-256 checksums). No internet needed; the qtl2 example data used by the tests ship with the package.
metadata:
  version: "1.0.0"
  domain: genomics
  language: R
  status: tested
  summary: "QTL mapping in biparental and multiparent populations with permutation thresholds (qtl2)"
  repository: https://github.com/bgtamang/AgriClaw
---

# QTL mapping with R/qtl2

## What it does

Runs single-QTL genome scans on an experimental population: genotype probabilities along the genome, a LOD score at every marker and pseudomarker, genome-wide significance thresholds by permutation, QTL peaks with intervals, effects and variance explained, plus marker and line QC. Biparental crosses can be given as plain CSV files; multiparent populations need an R/qtl2 control file.

## When to use / when not to use

- Use for populations derived from known parents: RIL (selfing or sib-mating), DH, F2, backcross, and qtl2 multiparent types (riself4/8/16, risib4/8, magic19, dh6, genril, DO, HS).
- Diversity panels and breeding lines without a cross design: use `claw-gwas-crop`.
- Clean markers first with `claw-marker-qc`; compute BLUEs with `claw-trial-analysis` or `claw-met-gxe`. For genes inside a QTL interval use `claw-gene-lookup` with the Mb interval; for polygenic prediction use `claw-genomic-prediction`.

## Inputs

**CSV route** (biparental): `--geno`, `--map`, `--pheno`, `--cross-type`.

- `--geno`: first column individual ID, one column per marker, codes such as `A`/`B` (parental homozygotes) and `H` (heterozygote, F2/backcross only). Missing: empty, `NA`, `-`, `.`, `N`. Codes are set with `--geno-codes`; defaults `A=1,B=2` (riself, risib, dh, haploid), `A=1,H=2,B=3` (f2), `A=1,H=2` (bc). `A` is the reference parent: a positive additive effect means the `B` allele increases the trait.
- `--map`: `marker,chr,pos` with positions in cM (genetic map). `--pmap` (optional): `marker,chr,pos` in Mb, used only to report intervals in Mb.
- `--pheno`: `--id-col` (default `id`) plus numeric traits, units in the names (`days_to_flowering`, `yield_kg_ha`). BLUEs are fine.
- `--covar` + `--addcovar`: additive covariates (for example a batch or environment column).

**Control-file route**: `--control file.yaml|.json|.zip` in the R/qtl2 format (needed for multiparent populations: founder genotypes and cross information). See the qtl2 input-file vignette.

**Ask the user** (never guess):

- The **cross type and generation**: RIL by selfing or sib-mating, DH, F2, backcross, or a multiparent design (and its funnel/cross information). A wrong cross type gives wrong probabilities.
- Whether the map is **genetic (cM) or physical (bp/Mb)**. QTL mapping needs a genetic map of this population or a published linkage map. Only if none exists, `--map-units bp|Mb --cm-per-mb X` converts with one genome-average rate; say clearly that this is a rough approximation.
- Which parent is `A` and which is `B`, and how missing and heterozygous calls are coded (RIL residual heterozygotes must be coded missing).
- Whether phenotypes are BLUEs/means and their units; whether any covariate (batch, environment, sex for animal data) should be included.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/qtl_mapping.R" --demo
Rscript "${CLAUDE_SKILL_DIR}/scripts/qtl_mapping.R" --geno geno.csv --map gmap.csv --pheno blues.csv \
    --cross-type riself --pmap pmap_mb.csv --out qtl_results/
Rscript "${CLAUDE_SKILL_DIR}/scripts/qtl_mapping.R" --geno f2_geno.csv --map gmap.csv --pheno pheno.csv \
    --cross-type f2 --geno-codes "AA=1,AB=2,BB=3" --interval bayes
Rscript "${CLAUDE_SKILL_DIR}/scripts/qtl_mapping.R" --control magic.yaml --kinship loco --traits yield_kg_ha
```

Permutations dominate run time (about 1000 scans per trait); with an X chromosome, X-specific permutations add many more. Use `--cores` to parallelize.

## Outputs

- `report.md`: warnings first, data and QC summary, method, thresholds, QTL table (position, interval in cM and Mb, LOD, significance, PVE, additive/dominance effects with SE), limitations.
- `tables/peaks.csv`: one row per peak and trait: `pos_cM`, `ci_lo_cM`/`ci_hi_cM` (and `_Mb`), `lod`, `significance` (genome-wide 5% or suggestive 10%), thresholds, `n`, `pve_pct`, nearest and interval-end markers, `a_effect`/`a_se` (= (BB - AA)/2), `d_effect`/`d_se` (F2), `AB_minus_AA` (backcross).
- `tables/peak_effects.csv`: genotype means or founder allele effects (centered) with SE at each peak.
- `tables/thresholds.csv`, `tables/lod_scores.csv` (every position x trait), `tables/marker_qc.csv`, `tables/individual_qc.csv`, `tables/phenotype_summary.csv`, `tables/duplicate_pairs.csv` (if any).
- `figures/lod_<trait>.png` (LOD curves with 5% and 10% thresholds), `figures/effects_<trait>_chr<c>_<pos>.png` (effects along the chromosome above the LOD curve).
- `qtl2_input/` (CSV route): the qtl2 files and `cross.yaml` built from the input, readable with `qtl2::read_cross2()`.
- `reproducibility/`: `commands.sh`, `sessionInfo.txt`, checksums.

## Expert decisions

- **Genotype probabilities**: `qtl2::calc_genoprob` (hidden Markov model; Broman et al. 2019) with genotyping error probability 0.002, the value used in every worked example of the qtl2 user guide and help pages (the function default is 1e-4). Raise it (for example to 0.01, rule of thumb) for noisy GBS calls. Haldane map function (qtl2 default; no crossover interference, as in the HMM). Pseudomarkers every 1 cM (user-guide example); `--step 0` scans markers only, adequate for dense SNP maps.
- **Scan**: Haley-Knott regression (Haley & Knott 1992), the qtl2 default, a fast approximation to interval mapping (Lander & Botstein 1989). `--kinship loco` fits a linear mixed model with leave-one-chromosome-out kinship (Broman et al. 2019) for populations with residual relatedness (multiparent designs, families); permutations are then much slower. Multiparent crosses are scanned with founder allele probabilities (additive allele model, as in the qtl2 user guide).
- **Thresholds**: 1000 permutations (`--n-perm`; qtl2 user-guide value; permutation thresholds of Churchill & Doerge 1994), genome-wide alpha = 0.05 (significant) and 0.10 (suggestive). With an X chromosome, separate autosome and X thresholds (Broman et al. 2006). Each trait is permuted separately with the same seed.
- **Intervals**: 1.5-LOD support intervals (R/qtl `lodint` default and qtl2 user-guide example) or 95% Bayes credible intervals (`--interval bayes`), both expanded to flanking markers. Manichaikul et al. (2006) showed both behave appropriately, unlike bootstrap intervals. Peaks above the 10% threshold; one per chromosome and trait unless `--peakdrop` is set (the user guide uses 1.8 "with caution").
- **Effects**: `qtl2::fit1` at the peak (proper standard errors; `scan1coef` SEs are not exact with zero-sum constraints). Additive a = (BB - AA)/2, dominance d = AB - (AA + BB)/2 (F2), AB - AA (backcross), centered founder allele effects (multiparent).
- **PVE** = 1 - 10^(-2 LOD / n): exact for Haley-Knott, since LOD = (n/2) log10(RSS0/RSS1); approximate for the mixed model. Per QTL, not additive across QTL.
- **QC rules of thumb** (adjustable): segregation distortion flagged at chi-square p < 0.05 / number of markers (Bonferroni); markers or individuals with > 20% missing calls; pairs of individuals with > 95% identical genotypes (`qtl2::compare_geno`); trait skewness > 1; values beyond 4 SD; fewer than 100 individuals (low power, Beavis effect); intervals wider than 20 cM.

**What the tests validate**: the forward RIL simulator matches Haldane & Waddington (1931) map expansion; the bundled iron data reproduce the qtl2 user guide exactly (peaks on chr 2, 7, 16, 8, 9 at 56.8, 50.1, 28.6, 13.6, 56.6 cM with LOD 4.9576 to 12.0632, all 1.5-LOD intervals, the two chr-7 peaks with peakdrop 1.8, and the chr-9 Bayes interval 53.6-61.2), and liver permutation thresholds within Monte Carlo error of the guide (100 permutations); fit1 reproduces the scan LOD; PVE equals the regression R^2; in 20 simulated RIL populations (150 lines, PVE 0.20) the 1.5-LOD interval contained the true QTL 20/20 times and the 95% Bayes interval 19/20, mean position error 2.0 cM, mean additive estimate 1.02 (true 1); a simulated 4-way RIL (riself4, control-file route) recovers the QTL and the founder allele effects, also with LOCO; a factor covariate reproduces qtl2's own covariate scan; QC flags duplicates and distorted markers; unit and code errors fail with advice. Not validated: multiparent cross types other than riself4, LOCO permutation thresholds (only the LOCO scan is checked), and X-chromosome-specific thresholds (too slow to test; they come directly from qtl2::scan1perm).

## Checks before reporting

- Read the Warnings section of `report.md` first and pass every warning on to the user.
- Look at `marker_qc.csv` and the LOD figures: a LOD spike at a single marker, or peaks only at distorted or high-missing markers, usually means genotyping errors.
- **Stop and tell the user** if: the map was converted from physical positions (results are provisional); many markers are distorted genome-wide (wrong cross type or mislabelled parents); duplicate lines are present; the phenotype is strongly skewed or has outliers; heterozygous calls appear in RIL/DH data.
- Report QTL above the 5% threshold as significant and those between the 10% and 5% thresholds only as suggestive.
- Several peaks for one trait on one chromosome, or wide intervals, need a multiple-QTL model before they are called separate QTL. This skill fits single-QTL scans only.
- With fewer than a few hundred individuals, warn that effects and PVE of detected QTL are overestimated (Beavis effect; Xu 2003).

## References

- Broman KW, Gatti DM, Simecek P, Furlotte NA, Prins P, Sen S, Yandell BS, Churchill GA (2019) R/qtl2: software for mapping quantitative trait loci with high-dimensional data and multiparent populations. Genetics 211:495-502. https://doi.org/10.1534/genetics.118.301595
- Broman KW. R/qtl2 user guide. https://kbroman.org/qtl2/assets/vignettes/user_guide.html
- Haley CS, Knott SA (1992) A simple regression method for mapping quantitative trait loci in line crosses using flanking markers. Heredity 69:315-324. https://doi.org/10.1038/hdy.1992.131
- Lander ES, Botstein D (1989) Mapping Mendelian factors underlying quantitative traits using RFLP linkage maps. Genetics 121:185-199. https://doi.org/10.1093/genetics/121.1.185
- Churchill GA, Doerge RW (1994) Empirical threshold values for quantitative trait mapping. Genetics 138:963-971. https://doi.org/10.1093/genetics/138.3.963
- Broman KW, Sen S, Owens SE, Manichaikul A, Southard-Smith EM, Churchill GA (2006) The X chromosome in quantitative trait locus mapping. Genetics 174:2151-2158. https://doi.org/10.1534/genetics.106.061176
- Manichaikul A, Dupuis J, Sen S, Broman KW (2006) Poor performance of bootstrap confidence intervals for the location of a quantitative trait locus. Genetics 174:481-489. https://doi.org/10.1534/genetics.106.061549
- Broman KW, Sen S (2009) A Guide to QTL Mapping with R/qtl. Springer, New York. https://doi.org/10.1007/978-0-387-92125-9
- Xu S (2003) Theoretical basis of the Beavis effect. Genetics 165:2259-2268. https://doi.org/10.1093/genetics/165.4.2259
- Haldane JBS, Waddington CH (1931) Inbreeding and linkage. Genetics 16:357-374. https://doi.org/10.1093/genetics/16.4.357 (RIL map expansion R = 2r/(1+2r), used by the tests to check the population simulator)
