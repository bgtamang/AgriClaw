---
name: claw-gwas-crop
description: Genome-wide association study (GWAS) for crop diversity panels and breeding lines with a mixed linear model (statgenGWAS). Uses a marker kinship matrix (Astle-Balding, VanRaden or IBS), optional leave-one-chromosome-out (LOCO) kinship, principal components (Q+K model) and covariates, and a MAF filter. Applies Bonferroni or Benjamini-Hochberg FDR thresholds and reports the genomic-control inflation factor lambda, with QQ and Manhattan plots. Groups significant SNPs into QTL regions using an LD distance. Use when the user has SNP genotypes (0/1/2 table or VCF) and one phenotype value (BLUE) per line, and asks for GWAS, association mapping, marker-trait associations, MLM, EMMA, Q+K or kinship models, significant SNPs, QTL regions or candidate windows, Manhattan or QQ plots, or why GWAS p-values look inflated.
license: MIT
compatibility: R >= 4.2 with statgenGWAS and vcfR (openssl optional, for SHA-256 checksums). No internet needed.
metadata:
  version: "1.0.0"
  domain: genomics
  language: R
  status: tested
  repository: https://github.com/bgtamang/AgriClaw
---

# Mixed-model GWAS for crop panels

## What it does

Tests every marker for association with each trait in a mixed linear model: y = covariates + marker + g + e, where the polygenic effect g ~ N(0, Vg K) uses a marker-based kinship matrix K. The model is fitted by statgenGWAS: EMMA REML variance components, then generalised least squares for each marker (Kang et al. 2008). The skill reports p-values, effects per ALT allele, genomic-control λ, QQ and Manhattan plots, a PCA of population structure, and QTL regions built from an LD distance.

## When to use / when not to use

- Use for association panels: diversity collections, breeding lines, cultivars and landraces with one value per line.
- Before running:
  - genotypes not yet QC'd, or in HapMap or samples-in-rows layout → `claw-marker-qc` first;
  - plot-level or multi-environment phenotypes → `claw-trial-analysis` or `claw-met-gxe` first, to get one BLUE per line;
  - LD window for QTL regions → `claw-ld-decay`.
- After running: genes in the QTL regions → `claw-gene-lookup`.
- Don't use for biparental populations (RIL, DH, F2): use `claw-qtl-mapping`, which models the cross. Don't use for predicting breeding values: use `claw-genomic-prediction`.

## Inputs

- `--geno`: the 0/1/2 marker table from `claw-marker-qc` (`marker, chr, pos`, optional `ref, alt`, then one column per line; values are copies of the ALT allele; `.csv`, `.tsv` or `.gz`), or a VCF with biallelic diploid GT calls. Positions in base pairs.
- `--pheno`: CSV with one row per line, an ID column (`--id-col`, default `genotype`) matching the genotype column names, numeric trait columns with units in the name (`yield_kg_ha`), and optional covariate columns (`--covar`). Repeated IDs stop the run.

**Ask the user (never guess):**

- Are the phenotypes BLUEs or means per line? Raw plot data must go through `claw-trial-analysis` first. Prefer BLUEs over BLUPs: BLUPs are already shrunk towards the mean (Piepho et al. 2008), which shrinks marker effects again.
- The units of each trait, and which traits to analyse.
- Is there known structure (subpopulations, breeding programmes, maturity groups)? Should flowering time or another covariate be fitted, for example for yield?
- Which genome assembly the marker positions use. `claw-gene-lookup` needs it.
- The LD distance for candidate windows. Take it from `claw-ld-decay` (`tables/ld_window.csv`) rather than a remembered value.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/gwas_crop.R" --demo
Rscript "${CLAUDE_SKILL_DIR}/scripts/gwas_crop.R" --geno qc/tables/genotypes_filtered_012.csv \
        --pheno blues.csv --traits yield_kg_ha,height_cm --ld-table ld/tables/ld_window.csv --out gwas/
# strongly structured panel, or lambda > 1.1: add principal components (Q+K)
Rscript "${CLAUDE_SKILL_DIR}/scripts/gwas_crop.R" --geno geno.csv --pheno blues.csv --traits yield_kg_ha \
        --n-pcs 3 --loco --compare-naive --out gwas_pcs/
```

The demo is a synthetic panel of 250 inbred lines from three subpopulations, with 1,347 markers. `height_cm` has three QTNs; `protein_pct` has no marker effects, only a subpopulation effect. The demo finds the three QTNs, reports nothing for `protein_pct`, and shows λ ≈ 0.95 with kinship against 3.2 to 5.8 without.

Other options:
- `--kinship astle|vanRaden|IBS`
- `--maf` (default 0.05)
- `--max-missing` (default 0.2)
- `--threshold bonferroni|fdr`, `--alpha` (default 0.05)
- `--ld-distance BP`
- `--n-cores`

## Outputs

- `report.md`: warnings first (inflation, deflation, low power, imputation, outlier lines), then the model, a per-trait table (λ, λ without kinship, hits, regions, detectable QTN size) and the QTL regions.
- `tables/gwas_<trait>.csv`: one row per marker with these columns:
  - `alt_allele_freq`, `maf`, `p_value`, `minus_log10_p`;
  - `effect_per_alt_allele` (trait units per ALT-allele copy; half the homozygote difference in inbred lines) and `effect_se`;
  - `r2_lr` (likelihood-ratio R², Sun et al. 2010), `q_value_bh` and `significant`.

  Markers below the MAF filter have `p_value` NA.
- `tables/qtl_regions.csv`: `trait, region_id, chr, start_bp, end_bp, lead_snp, lead_pos_bp, lead_p, n_significant_snps, ld_distance_bp, lead_effect_per_alt_allele, lead_alt_allele_freq`. This is the input for `claw-gene-lookup`.
- `tables/significant_snps.csv`, `tables/principal_components.csv` and `tables/diagnostics.csv`. The diagnostics table has `lambda_gc`, `lambda_naive`, `vg`, `ve`, `bonferroni_p` and `min_r2_80pct_power`.
- `figures/manhattan_<trait>.png`, `figures/qq_<trait>.png` (with a 95% band under the null) and `figures/pca_structure.png` (PC1 vs PC2 and a scree plot).

## Expert decisions

- **Kinship is essential in structured panels.** Without it, structure and relatedness create false associations. Mixed models with kinship (Yu et al. 2006; Kang et al. 2008) control this.
  - On the simulated structured panels in the tests (300 lines, 3 subpopulations, 12 seeds), a trait with no marker effects gave λ 7.7 to 9.6 without kinship and 0.94 to 1.08 with kinship.
  - `--compare-naive` reports the no-kinship λ, so the user can see this on their own data.
- **Kinship method:** `astle` (Astle & Balding 2009) by default, the statgenGWAS default; `vanRaden` (VanRaden 2008) and `IBS` are alternatives. Markers below the MAF filter are left out of the kinship matrix.
- **LOCO** (`--loco`): each chromosome is tested with a kinship matrix built from the other chromosomes. This avoids the power loss from including the tested marker in K (Rincent et al. 2014).
  - In our simulations, LOCO alone gave λ 1.10 to 1.21 on traits without marker effects. Adding 2 PCs brought it to 0.93 to 1.12.
  - So with `--loco` in a structured panel, also add PCs.
- **PCs as covariates** (`--n-pcs`, the Q+K model; Price et al. 2006; Yu et al. 2006). The default is 0, because kinship alone was calibrated in our tests. Add 2 to 5 PCs, choosing from the scree plot, when λ > 1.1 or the PCA shows discrete groups. The number of PCs is a judgement call: report it.
- **Missing genotypes:** statgenGWAS needs complete markers. Markers with more than 20% missing calls are dropped (`--max-missing`, rule of thumb), and the rest are mean-imputed (2 × ALT frequency). The report states the imputed fraction and recommends Beagle-type imputation above 10% (rule of thumb).
- **MAF ≥ 0.05** (adjustable rule of thumb): markers carried by few lines give unstable tests.
- **Thresholds:**
  - **Bonferroni** (default): family-wise α = 0.05 over the markers actually tested. It is conservative because markers in LD are not independent.
  - **FDR** (`--threshold fdr`): Benjamini-Hochberg q ≤ α over the tested markers (Benjamini & Hochberg 1995). FDR control holds under positive dependence (Benjamini & Yekutieli 2001); it accepts some false positives in exchange for power.
  - **statgenGWAS `thrType = "fdr"` is not used.** In version 1.0.14 its code compares the k-th cluster representative with α/k after an unadjusted p < α gate, instead of the documented k·α/m (Brzyski et al. 2017). On 12 simulated traits with no marker effects at all, it called 8 to 13 SNPs significant per trait (p up to 0.006); Benjamini-Hochberg called none.
- **Genomic-control λ** (Devlin & Roeder 1999), as computed by statgenGWAS from the median test statistic. The warning limits, > 1.1 for inflation and < 0.9 for deflation, are adjustable rules of thumb.
  - A highly polygenic trait in a large panel can raise λ without any confounding (Yang et al. 2011).
  - So investigate before calling it confounding, but never interpret hits from an inflated scan without checking.
- **Power:** below about 200 lines (rule of thumb) power is low. The report gives the smallest QTN, as a fraction of phenotypic variance, detectable with 80% power at the threshold used. It uses a 1-df test with non-centrality n·R²/(1−R²) (Sham & Purcell 2014), ignoring the extra loss when a QTN is confounded with structure. At Bonferroni α = 0.05:

  | lines | 5,000 markers | 50,000 markers |
  |---|---|---|
  | 100 | 22% | 25% |
  | 200 | 12% | 14% |
  | 500 | 5% | 6% |

- **QTL regions:** significant SNPs on a chromosome closer than the LD distance are joined, and the region is extended by that distance on both sides, starting at bp 1 at the chromosome end. The lead SNP has the smallest p. Without an LD distance, regions are not extended, and the report says so.

## Checks before reporting

- Read the Warnings in `report.md` and relay each one.
- **QQ plot:** the points should follow the diagonal and leave the band only at the tail. An early, whole-curve departure means inflation: stop, add PCs, check duplicates and outliers, and re-run before interpreting any hit.
- **Stop and tell the user** when:
  - λ > 1.1 or < 0.9;
  - fewer than 200 lines (state the detectable QTN size);
  - the phenotype file had repeated IDs or obvious outliers;
  - many hits for a trait with low heritability;
  - a hit rests on a SNP with MAF near the filter;
  - IDs failed to match between the two files.
- Compare λ with and without kinship (`--compare-naive`). A large gap confirms the panel is structured and that the kinship model matters.
- Report lead SNPs as markers in LD with causal variants, not as causal genes. Give region coordinates together with the genome assembly.
- Don't claim a hit replicates unless it was tested in an independent panel or environment.

## References

- Astle W, Balding DJ (2009) Population structure and cryptic relatedness in genetic association studies. Statistical Science 24:451-471. https://doi.org/10.1214/09-STS307
- Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society B 57:289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Benjamini Y, Yekutieli D (2001) The control of the false discovery rate in multiple testing under dependency. Annals of Statistics 29:1165-1188. https://doi.org/10.1214/aos/1013699998
- Brzyski D, Peterson CB, Sobczyk P, Candès EJ, Bogdan M, Sabatti C (2017) Controlling the rate of GWAS false discoveries. Genetics 205:61-75. https://doi.org/10.1534/genetics.116.193987
- Devlin B, Roeder K (1999) Genomic control for association studies. Biometrics 55:997-1004. https://doi.org/10.1111/j.0006-341X.1999.00997.x
- Kang HM, Zaitlen NA, Wade CM, Kirby A, Heckerman D, Daly MJ, Eskin E (2008) Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723. https://doi.org/10.1534/genetics.107.080101
- Piepho HP, Möhring J, Melchinger AE, Büchse A (2008; online 2007) BLUP for phenotypic selection in plant breeding and variety testing. Euphytica 161:209-228. https://doi.org/10.1007/s10681-007-9449-8
- Price AL, Patterson NJ, Plenge RM, Weinblatt ME, Shadick NA, Reich D (2006) Principal components analysis corrects for stratification in genome-wide association studies. Nature Genetics 38:904-909. https://doi.org/10.1038/ng1847
- Rincent R et al. (2014) Recovering power in association mapping panels with variable levels of linkage disequilibrium. Genetics 197:375-387. https://doi.org/10.1534/genetics.113.159731
- Sham PC, Purcell SM (2014) Statistical power and significance testing in large-scale genetic studies. Nature Reviews Genetics 15:335-346. https://doi.org/10.1038/nrg3706
- Sun G, Zhu C, Kramer MH, Yang SS, Song W, Piepho HP, Yu J (2010) Variation explained in mixed-model association mapping. Heredity 105:333-340. https://doi.org/10.1038/hdy.2010.11
- VanRaden PM (2008) Efficient methods to compute genomic predictions. Journal of Dairy Science 91:4414-4423. https://doi.org/10.3168/jds.2007-0980
- Yang J, Weedon MN, Purcell S, et al., the GIANT Consortium (2011) Genomic inflation factors under polygenic inheritance. European Journal of Human Genetics 19:807-812. https://doi.org/10.1038/ejhg.2011.39
- Yu J et al. (2006; online 2005) A unified mixed-model method for association mapping that accounts for multiple levels of relatedness. Nature Genetics 38:203-208. https://doi.org/10.1038/ng1702
- statgenGWAS R package (van Rossum BJ, Kruijer W). https://CRAN.R-project.org/package=statgenGWAS
