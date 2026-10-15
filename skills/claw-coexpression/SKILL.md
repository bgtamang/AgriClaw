---
name: claw-coexpression
description: Weighted gene co-expression network analysis (WGCNA) for crop RNA-seq - filters low-expression and low-variance genes, flags outlier samples, optional ComBat batch removal, picks the soft-thresholding power for a signed network, finds co-expression modules (blockwiseModules), module eigengenes, module-trait correlations with Benjamini-Hochberg correction, hub genes by module membership (kME), and a Cytoscape edge list. Use when the user asks for WGCNA, co-expression modules, gene co-expression networks, hub genes, module-trait relationships, eigengenes, or which gene modules relate to traits such as yield, drought tolerance, stomatal conductance or nutrient content across many samples (at least 15, preferably 20 or more) of soybean, maize, wheat, rice or other plants.
license: MIT
compatibility: R >= 4.3 with WGCNA (CRAN), plus Bioconductor DESeq2 (for --counts) and sva (for --batch). Tested with R 4.6.1, WGCNA 1.74, DESeq2 1.52, sva 3.60. No internet needed.
metadata:
  version: "1.0.0"
  domain: transcriptomics
  language: R
  status: tested
  summary: "WGCNA co-expression modules and module-trait relationships"
  repository: https://github.com/bgtamang/AgriClaw
---

# Co-expression modules with WGCNA

## What it does

Builds a signed weighted co-expression network from variance-stabilized expression data, cuts it into modules of co-expressed genes, summarizes each module by its eigengene (first principal component), correlates eigengenes with sample traits, ranks hub genes by module membership (kME), and exports within-module edges for Cytoscape. Sample outliers, sample size and scale-free fit are checked along the way.

## When to use / when not to use

- Use with 15 or more samples (20+ preferred) that vary in genotype, environment, time or trait values, to find groups of co-regulated genes and relate them to traits.
- Input: `tables/vst.csv` from **claw-de-analysis** (VST, low counts already pre-filtered), a log-CPM matrix, or raw counts (`--counts`, filtered and transformed here).
- Don't use for two-group comparisons with few replicates (use **claw-de-analysis**), or on a list of DE genes only (see below).
- Module gene lists can go to **claw-enrichment**, with the genes in the network as the universe.

## Inputs

- `--expr`: genes × samples CSV, first column `gene_id`, on a log scale (VST or log-CPM). Raw or normalized counts are refused (values > 100).
- `--counts`: raw integer counts; genes with < 10 counts in more than 90% of samples are removed, then DESeq2 VST is applied.
- `--traits`: CSV with `sample` plus trait columns, units in names (`yield_g_plant`, `stomatal_conductance_mol_m2_s`). Text columns with 2 levels become 0/1 (the coding is reported); with more levels, one 0/1 column per level.
- `--batch` (and `--protect`): a traits column to remove with ComBat, protecting listed biological columns.

**Ask the user, never guess:**
1. How many samples, and what drives their differences (genotypes, environments, time points, tissues)? Mixed tissues or species should usually be analyzed separately (consensus analysis is not covered here).
2. Are there batches (sequencing run, harvest date, field block)? Batch effects should be removed before the network is built (WGCNA FAQ), but only when the batch is not confounded with the biology.
3. Which traits matter, measured on the same samples, and in what units?
4. Has the expression matrix been filtered by differential expression? It must not be: run on all expressed genes.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/coexpression.R" --demo --out demo_wgcna/
Rscript "${CLAUDE_SKILL_DIR}/scripts/coexpression.R" --expr de/tables/vst.csv --traits traits.csv --out wgcna/
Rscript "${CLAUDE_SKILL_DIR}/scripts/coexpression.R" --counts counts.csv --traits traits.csv --batch harvest_date --protect genotype
```

| Option | Default | Meaning |
|---|---|---|
| `--n-top-var` | 5000 | most variable genes kept (0 = all) |
| `--min-count` / `--min-fraction` | 10 / 0.1 | `--counts` filter: ≥ 10 counts in ≥ 10% of samples |
| `--min-mean` | none | `--expr` mean-expression filter |
| `--power` | auto | soft-threshold power |
| `--r2-cut` | 0.8 | signed scale-free fit target |
| `--cor-type` / `--max-p-outliers` | bicor / 0.05 | correlation; bicor outlier cap |
| `--min-module-size` | 30 | smallest module |
| `--merge-cut-height` | 0.25 | merge modules with eigengene correlation > 0.75 |
| `--deep-split` | 2 | tree-cut sensitivity |
| `--outlier-z` / `--remove-outliers` | -2.5 / off | sample connectivity flag; removal only on request |
| `--allow-small-n` | off | run with fewer than 15 samples (exploratory only) |
| `--hub-n` / `--hub-kme` | 10 / 0.8 | hub genes per module; kME needed to be called a hub |
| `--edge-threshold` / `--max-edges` | 0.1 / 2000 | Cytoscape export |

## Outputs

- `report.md`: warnings, data and filtering, chosen power and fit, modules with top hub, module-trait table, methods.
- `tables/gene_modules.csv` (`module` number, 0 = unassigned, which WGCNA labels `grey`; `module_color`; `kME_own`), `tables/kME_all_modules.csv`, `tables/module_sizes.csv`.
- `tables/module_eigengenes.csv`, `tables/module_trait_correlations.csv` (`r`, `n`, `pvalue`, `padj`, `significant`).
- `tables/hub_genes.csv`, `tables/cytoscape_edges.csv` (fromNode, toNode, weight = topological overlap, module) and `tables/cytoscape_nodes.csv`.
- `tables/soft_threshold.csv`, `tables/sample_connectivity.csv` (Z.k per sample).
- `figures/`: `sample_outliers.png`, `soft_threshold.png`, `gene_dendrogram.png`, `module_trait.png`, `kme_vs_gene_significance.png`.

## Expert decisions

- **Sample size:** the WGCNA FAQ advises against fewer than 15 samples and recommends at least 20; below 15 the script stops unless `--allow-small-n` (rule of thumb from the FAQ).
- **Filtering:** remove consistently low counts (< 10 in more than 90% of samples) and apply a variance-stabilizing transformation (WGCNA FAQ); then keep the 5000 most variable genes (adjustable rule of thumb that keeps one network block). **Never filter by differential expression:** it builds artificially few, highly correlated modules and invalidates the scale-free criterion (WGCNA FAQ).
- **Signed network, bicor with maxPOutliers 0.05:** the FAQ prefers signed networks and recommends the biweight midcorrelation with `maxPOutliers` 0.05-0.10 (Langfelder & Horvath 2012). Signed TOM (`TOMType = "signed"`).
- **Soft-threshold power:** the lowest power with signed scale-free fit R² ≥ 0.8 (Zhang & Horvath 2005 approximate scale-free topology criterion; 0.8 is the FAQ's working target, an adjustable rule of thumb). If no power up to 30 reaches it, the FAQ table for signed networks is used (18 for < 20 samples, 16 for 20-30, 14 for 30-40, 12 for > 40) and the report warns that a strong driver (batch, outliers, a dominant treatment) is likely.
- **Modules:** blockwiseModules with dynamic tree cut (Langfelder et al. 2008), minimum size 30, merge height 0.25, deepSplit 2, `pamRespectsDendro = FALSE` and `reassignThreshold = 0` (the WGCNA tutorial settings; adjustable).
- **Sample outliers:** standardized connectivity Z.k from a signed sample network ((1 + r)/2)^2 (Oldham et al. 2012, who used Z.k = -2 as an example cut-off; -2.5 here, an adjustable rule of thumb). A sample with an extreme but genuine profile can also fall below the cut-off: in the test simulation without any planted outlier one sample still reached Z.k -3.1, and in a similar simulation the lowest sample was the one with the most extreme module eigengene values. Flags are therefore not removals; `--remove-outliers` must be requested.
- **Batch:** ComBat (sva; Johnson et al. 2007), as the FAQ suggests for categorical batches, protecting `--protect` covariates. The script stops if the batch is completely confounded with a protected factor.
- **Module-trait statistics:** Pearson correlation of eigengenes with traits, Student p-values (`corPvalueStudent`), BH over all module × trait tests, excluding the unassigned module. The report states the smallest |r| that can reach p < 0.05 at the given n (0.40 at n = 24).
- **Hub genes:** ranked by kME (signed eigengene correlation, `signedKME`); `hub = TRUE` when kME ≥ 0.8 (adjustable rule of thumb).
- **Edges:** within-module signed topological overlap from the TOM that blockwiseModules computed; edges ≥ 0.1, the strongest 2000 per module (adjustable).

### Validation (tests/test_coexpression.R)

- WGCNA::simulateDatExpr data (800 genes, 5 planted modules, 40 samples): all 5 modules found, adjusted Rand index (Hubert & Arabie 1985) 0.95 among planted-module genes and 0.78 over all genes. The lower overall value comes from genes near the boundary between a module and the unassigned genes (true r about 0.4). Every gene with true kME ≥ 0.7 was assigned to its module.
- A trait planted on one eigengene (r about 0.8) is linked to the matching module (padj < 0.05); a null trait is linked to none.
- Soft-threshold rules, the FAQ fallback table, a degraded sample (Z.k -5.7, next lowest -1.1), ComBat removal of a +2 shift, and exported edge weights equal to `TOMsimilarityFromExpr` are all checked.
- The demo (24 samples) recovers its 4 planted modules and the planted stomatal-conductance and leaf-N links, with no module for the null trait.

## Checks before reporting

- Read the Warnings section first and pass every warning on.
- `figures/sample_outliers.png`: a sample far from the rest or a clear split into groups (batch, tissue) must be explained before modules are interpreted. **Stop and ask** if the clusters match a batch that is not in `--batch`.
- `figures/soft_threshold.png`: if the FAQ fallback power was used, say so and explain the likely driver.
- If most genes are unassigned or one module holds most genes, report it; the network has little structure (or one dominant driver).
- **Don't over-interpret module-trait correlations at small n:** with 20-30 samples, r estimates are imprecise and only |r| ≥ 0.36-0.44 reaches nominal p < 0.05. Report r, n and padj, and call the links hypotheses, not mechanisms.
- A module correlated with treatment largely restates the DE result. Traits measured on other plants or other dates than the RNA samples cannot be correlated.
- Hub genes are candidates for follow-up, not validated regulators.

## References

- Langfelder P, Horvath S (2008). WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9:559. https://doi.org/10.1186/1471-2105-9-559
- Zhang B, Horvath S (2005). A general framework for weighted gene co-expression network analysis. Statistical Applications in Genetics and Molecular Biology 4:Article17. https://doi.org/10.2202/1544-6115.1128
- Langfelder P, Zhang B, Horvath S (2008). Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R. Bioinformatics 24:719-720. https://doi.org/10.1093/bioinformatics/btm563
- Langfelder P, Horvath S (2012). Fast R functions for robust correlations and hierarchical clustering. Journal of Statistical Software 46(11). https://doi.org/10.18637/jss.v046.i11
- Oldham MC, Langfelder P, Horvath S (2012). Network methods for describing sample relationships in genomic datasets: application to Huntington's disease. BMC Systems Biology 6:63. https://doi.org/10.1186/1752-0509-6-63
- Johnson WE, Li C, Rabinovic A (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics 8:118-127. https://doi.org/10.1093/biostatistics/kxj037
- Hubert L, Arabie P (1985). Comparing partitions. Journal of Classification 2:193-218. https://doi.org/10.1007/BF01908075
- Langfelder P, Horvath S. WGCNA package: Frequently Asked Questions (last updated December 24, 2017). https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html (the live page returned HTTP 503 in October 2026; archived copy: http://web.archive.org/web/20230227152143/https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html)
