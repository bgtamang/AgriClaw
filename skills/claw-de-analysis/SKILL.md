---
name: claw-de-analysis
description: Differential expression analysis of RNA-seq counts with DESeq2 - design formula with batch or block, explicit control/reference levels, pairwise contrasts (Wald test) or likelihood-ratio tests for time courses and multi-level factors, independent filtering, Benjamini-Hochberg FDR, apeglm log2 fold-change shrinkage for ranking, MA plot, volcano plot and PCA. Refuses TPM/FPKM input, unreplicated groups and batch confounded with treatment. Works for any crop (soybean, maize, wheat, rice) from tximport (Salmon/kallisto) or featureCounts counts. Use when the user asks for DEGs, differentially expressed genes, DESeq2, fold changes, up- and down-regulated genes between treatments or genotypes, drought, heat or stress versus control comparisons, time-course RNA-seq, or volcano and MA plots.
license: MIT
compatibility: R >= 4.3 with Bioconductor DESeq2, apeglm, SummarizedExperiment and S4Vectors. Tested with R 4.6.1, DESeq2 1.52 and apeglm 1.34. No internet needed.
metadata:
  version: "1.0.0"
  domain: transcriptomics
  language: R
  status: tested
  summary: "Differential expression with DESeq2: contrasts, likelihood-ratio tests, shrunken fold changes"
  repository: https://github.com/bgtamang/AgriClaw
---

# Differential expression with DESeq2

## What it does

Fits DESeq2's negative-binomial model to raw RNA-seq counts (with tximport length offsets when available), tests the requested comparisons with Wald tests or a likelihood-ratio test (LRT), controls the false discovery rate with Benjamini-Hochberg and independent filtering, and shrinks log2 fold changes with apeglm for ranking and plotting. It stops on inputs and designs that cannot give valid results.

## When to use / when not to use

- Use for treatment vs control, genotype comparisons, multi-level factors and time courses with biological replicates.
- Input comes from **claw-rnaseq-import** (`txi.rds`) or a raw count matrix (featureCounts, HTSeq).
- Don't use without biological replicates: no valid DE test exists; describe fold changes only and say so.
- Don't use on TPM, FPKM, CPM or log values: go back to the counts (claw-rnaseq-import).
- Not for microarrays (use limma), single-cell data or differential isoform usage.
- Next steps: **claw-enrichment** (GO/KEGG enrichment; use `padj`-tested genes as background), **claw-coexpression** (WGCNA on `tables/vst.csv` with 15+ samples).

## Inputs

- Counts, one of:
  - `--txi txi.rds` from claw-rnaseq-import (estimated counts plus average transcript lengths, used as offsets);
  - `--counts counts.csv`: first column `gene_id`, then one column of raw integer counts per sample; optional `--lengths gene_lengths.csv` (same shape) to use tximport offsets.
- `--samples samples.csv`: a `sample` column matching the count columns, plus every design variable.
- `--design`: R formula with the variable of interest last, e.g. `~ batch + condition`.

**Ask the user, never guess:**
1. Which column is the treatment, and which level is the control (`--reference condition=control`)? Alphabetical order would make `mutant` the reference over `WT` and flip every fold change.
2. What batches exist: sequencing run or lane, RNA extraction day, growth chamber, planting date, field block or harvest time? Put each one that is not confounded with the treatment into the design.
3. Is there pairing or blocking (same plant or plot sampled twice, field blocks)? Add it as a term (e.g. `~ block + treatment`).
4. Which comparisons answer the question? Give them with `--contrast`.
5. For time courses: is the question "does the response over time differ between genotypes" (LRT: full `~ genotype + time + genotype:time`, reduced `~ genotype + time`) or "which genes change over time at all" (LRT: reduced `~ genotype`)?
6. Were any samples flagged by claw-rnaseq-import QC, and should they be excluded (`--exclude`, with a reason)?
7. Is the library 3' tag-seq (QuantSeq)? Then add `--no-length-offset`.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/de_analysis.R" --demo --out demo_de/
Rscript "${CLAUDE_SKILL_DIR}/scripts/de_analysis.R" --txi rnaseq_import/txi.rds --samples rnaseq_import/tables/samples.csv \
    --design "~ batch + condition" --reference condition=control --out de/
Rscript "${CLAUDE_SKILL_DIR}/scripts/de_analysis.R" --counts counts.csv --samples samples.csv --design "~ block + treatment" \
    --reference treatment=control --contrast "treatment,drought,control;treatment,heat,control;treatment,heat,drought"
Rscript "${CLAUDE_SKILL_DIR}/scripts/de_analysis.R" --counts counts.csv --samples samples.csv --test LRT --as-factor time \
    --design "~ genotype + time + genotype:time" --reduced "~ genotype + time" --reference genotype=WT,time=0
```

| Option | Default | Meaning |
|---|---|---|
| `--reference` | none | reference levels, `factor=level,...`; required for the variable of interest unless `--contrast` names both levels |
| `--contrast` | every level vs reference | `factor,numerator,denominator` or a coefficient name; `;` separates several |
| `--test` / `--reduced` | Wald / none | `LRT` with a reduced formula for multi-level factors, interactions and time courses |
| `--alpha` | 0.05 | FDR for calling genes; also the target of independent filtering |
| `--lfc-threshold` | 0 | test \|log2FC\| > threshold instead of > 0 |
| `--shrink` | apeglm | `none` to skip shrinkage |
| `--min-count` / `--min-samples` | 10 / smallest group | pre-filter |
| `--min-reps` | 3 | warn below this many replicates |
| `--exclude` | none | samples to drop (reported) |
| `--as-factor` | none | numeric columns to treat as groups (time points, blocks) |
| `--no-length-offset` | off | ignore tximport lengths (3' tag-seq) |
| `--allow-noninteger` | off | round tximport scaledTPM/lengthScaledTPM counts |

## Outputs

- `report.md`: warnings, design and group sizes, results per comparison with the top genes, methods text.
- `tables/de_<comparison>.csv`, every input gene: `baseMean`, `log2FoldChange` (MLE), `lfcSE`, `stat`, `pvalue`, `padj`, `log2FoldChange_shrunk`, `lfcSE_shrunk` (apeglm), `significant` (padj < alpha), `direction` (up/down/ns), `status` (tested; filtered_low_mean; count_outlier; prefiltered). Genes with a `padj` are the tested universe for claw-enrichment.
- `tables/de_summary.csv`: tested genes, up, down, filter threshold, Cook's outliers per comparison.
- `tables/vst.csv` (VST, `blind = FALSE`; for clustering and claw-coexpression), `tables/normalized_counts.csv`, `tables/size_factors.csv`, `tables/samples_used.csv`.
- `figures/`: `ma_*.png`, `volcano_*.png`, `pvalue_hist_*.png`, `pca.png`, `dispersion.png`.

## Expert decisions

- **Counts, not abundances.** tximport counts go in with their average transcript lengths as offsets (`DESeqDataSetFromTximport`; Soneson et al. 2015). Inputs whose columns sum to about 10^6 (TPM/CPM) or that hold non-integers without offsets (FPKM, log values) are refused (Zhao et al. 2020 explain why such normalised values mislead between-sample comparisons).
- **Pre-filter** ≥ 10 counts in at least the smallest group size: the DESeq2 vignette's suggestion, for speed and cleaner plots. The statistical filtering is DESeq2's independent filtering on the mean count, optimised for `alpha` (Bourgon et al. 2010); filtered genes get `padj = NA`.
- **Wald test on the unshrunken estimate**, BH-adjusted, alpha 0.05 (DESeq2 defaults except alpha, which `results()` must receive to tune the filter). Cook's-distance outliers follow DESeq2's defaults: counts are replaced when a group has 7 or more replicates; otherwise a gene with an outlier in a group of 3 or more replicates gets p = NA (`status = count_outlier`).
- **apeglm shrinkage is for ranking and plots only** (Zhu et al. 2019). The p-values come from the Wald test on the MLE; the shrunken column has no p-value. apeglm needs a model coefficient, so for a denominator that is not the reference level the factor is re-levelled and `nbinomWaldTest` re-run (DESeq2 vignette).
- **LRT** for 3+ level factors, interactions and time courses: it tests all dropped terms together; its `log2FoldChange` column describes only one coefficient (named in the report).
- **`--lfc-threshold`** tests against a fold-change threshold. Filtering significant genes by fold change afterwards does not control the FDR for that threshold.
- **Refusals.** A compared group with one sample (no biological replicate; an unreplicated level outside the comparison only gives a warning), a design matrix that is not full rank (batch completely confounded with treatment, nesting, empty interaction cells; DESeq2 vignette, "Model matrix not full rank"), no residual degrees of freedom, and a missing reference level all stop the run with an explanation. Partly confounded batches give a warning.
- **Replicates:** at least 2 per tested group are required; fewer than 3 gives a warning. Schurch et al. (2016) found that with 3 replicates most tools detect only 20-40% of the DE genes found with 42 replicates, and recommend at least 6 (12 to detect small fold changes).

### Validation (tests/test_de_analysis.R, tests/test_de_lrt_cli.R)

- Negative-binomial simulations (1200 genes, 4 vs 4 in two batches, 10% DE at 2- to 4-fold): mean empirical FDR 0.06 over 5 seeds at nominal 0.05 (0.065 over 15 seeds in a separate calibration), power 0.88. The test allows FDR ≤ 0.08.
- apeglm: shrunken fold changes of strongly expressed DE genes lie within 0.2 log2 units of the truth (median); noisy null estimates shrink by more than half; Wald p-values are unchanged.
- Reversed contrasts negate fold changes exactly; the re-levelled apeglm path mirrors the forward contrast (r > 0.999).
- tximport offsets: with simulated isoform switching (same abundance, doubled length), ≤ 10% of switching genes are called DE with offsets, versus ≥ 50% without.
- **LRT limitation:** with 3 replicates per genotype × time cell, DESeq2's LRT p-values were anti-conservative in AgriClaw simulations: among null genes, 2.0% had p < 0.01 and 0.34% had p < 0.001. Empirical FDR was about 0.10 at nominal 0.05 (mean 0.098 over 10 seeds) and about 0.06 at nominal 0.025 (0.058 over 8 seeds). The report warns in this situation. Lund et al. (2012) proposed quasi-likelihood F-tests with shrunken dispersions, which compared favourably for FDR estimation in their simulations (edgeR's glmQLFit/glmQLFTest), but they are not part of this skill.
- Refusals are tested: TPM, FPKM, log values, NA counts, unreplicated groups, confounded and nested designs, and a missing reference level.

## Checks before reporting

- Read the Warnings section of `report.md` first and pass every warning on.
- **p-value histogram** of tested genes: expect a flat histogram with a peak near 0. A U shape, a hump in the middle or a peak near 1 suggests a missing batch term, outlier samples or a wrong model. Stop and tell the user.
- **Dispersion plot:** gene estimates should scatter around the fitted trend. A cloud far above it points to outliers or mixed sample types.
- **PCA:** if samples group by an unmodelled factor (lane, date, block), stop and ask whether to add it to the design.
- More than about 1% of genes flagged as count outliers suggests a bad sample; go back to the claw-rnaseq-import QC.
- Never change `--alpha`, drop samples or switch the design after seeing results in order to obtain more DE genes. Report zero or few DE genes honestly, with group sizes.
- Report `padj` from the test. Use shrunken fold changes for ranking and effect sizes, and say which column was used.
- For LRT time courses with 3 or fewer replicates, quote the calibration warning and consider `--alpha 0.025`.

## References

- Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15:550. https://doi.org/10.1186/s13059-014-0550-8
- Zhu A, Ibrahim JG, Love MI (2019). Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics 35:2084-2092. https://doi.org/10.1093/bioinformatics/bty895
- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society B 57:289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Bourgon R, Gentleman R, Huber W (2010). Independent filtering increases detection power for high-throughput experiments. PNAS 107:9546-9551. https://doi.org/10.1073/pnas.0914005107
- Schurch NJ, Schofield P, Gierliński M, et al. (2016). How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? RNA 22:839-851. https://doi.org/10.1261/rna.053959.115
- Soneson C, Love MI, Robinson MD (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research 4:1521. https://doi.org/10.12688/f1000research.7563.1
- Zhao S, Ye Z, Stanton R (2020). Misuse of RPKM or TPM normalization when comparing across samples and sequencing protocols. RNA 26:903-909. https://doi.org/10.1261/rna.074922.120
- Lund SP, Nettleton D, McCarthy DJ, Smyth GK (2012). Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates. Statistical Applications in Genetics and Molecular Biology 11(5). https://doi.org/10.1515/1544-6115.1826
- DESeq2 vignette (Bioconductor): https://bioconductor.org/packages/DESeq2/
