---
name: claw-enrichment
description: Gene set over-representation analysis (ORA) for crop RNA-seq results with clusterProfiler - GO (BP, MF, CC), KEGG, MapMan, PANTHER or any gene-to-term table from crop annotation files (Phytozome annotation_info, Ensembl Plants BioMart, GAF), no OrgDb package needed. Uses the tested (expressed) genes as the background universe, optional GO propagation to parent terms with GO.db, gene-set size limits, Benjamini-Hochberg correction, and an ID check that catches mismatched gene IDs (e.g. soybean Glyma.01G000100 vs Glyma01g00010, or transcript vs gene IDs). Use when the user asks for GO enrichment, KEGG pathway enrichment, functional enrichment, over-represented terms or pathways among DE genes, up- or down-regulated genes, or the genes of a co-expression module, in soybean, maize, wheat, rice or other plants.
license: MIT
compatibility: R >= 4.3 with Bioconductor clusterProfiler, GO.db and AnnotationDbi. Tested with R 4.6.1, clusterProfiler 4.20 and GO.db 3.23 (GO release 2026-01-23). No internet needed.
metadata:
  version: "1.0.0"
  domain: transcriptomics
  language: R
  status: tested
  repository: https://github.com/bgtamang/AgriClaw
---

# Gene set enrichment (over-representation)

## What it does

Tests whether GO terms, pathways or any user-defined gene sets contain more differentially expressed genes than expected, using the hypergeometric test in `clusterProfiler::enricher` with a gene-to-term table from the crop's own annotation. The background is the set of genes that were actually tested. GO annotations can be propagated to ancestor terms, and gene IDs are checked against the annotation before anything is tested.

## When to use / when not to use

- Use after **claw-de-analysis** (pass its full `tables/de_*.csv`) or on a module gene list from **claw-coexpression**.
- Works with any gene-to-term mapping: GO from Phytozome/Ensembl Plants/PLAZA/GAF files, KEGG pathways, MapMan bins, PANTHER families, custom lists.
- Don't use it to compare whole ranked lists (that is gene set enrichment analysis, GSEA, not covered here), or when almost all tested genes are DE.
- It doesn't correct for gene-length or expression-level bias in DE detection (GOseq does; Young et al. 2010). Say so when long or highly expressed gene families dominate the results.

## Inputs

- **Gene sets**, one of:
  - `--de-table`: the full DE results with `gene_id`, `padj` and optionally `log2FoldChange`. DE genes = `padj < --alpha`; the universe = every gene with a non-missing `padj`. Up, down and all DE genes are tested separately.
  - `--genes` + `--universe`: plain ID lists (one per line, or a CSV with `gene_id`). The universe is required.
- **Annotation** (`--annotation`):
  - long table: one gene-term pair per row;
  - wide table: one row per gene, terms separated by `,` `;` `|` or spaces (e.g. the `GO` column of Phytozome `annotation_info.txt`, with `--gene-col locusName --term-col GO`);
  - GO GAF file (`!gaf-version` header): NOT-qualified annotations are excluded.
  - `--term-names` (term_id, name) for non-GO sets; GO names come from GO.db.

**Ask the user, never guess:**
1. Which genome assembly and annotation release do the DE gene IDs come from, and is the annotation file from the same release? (Soybean Wm82.a1 `Glyma01g00010` and Wm82.a2+ `Glyma.01G000100` are different gene models; maize v3/v4/v5, rice MSU vs RAP-DB and wheat IWGSC versions differ too.)
2. Is the annotation per gene or per transcript/protein? Per transcript needs `--strip-isoform`.
3. Does the GO file list only the most specific terms (most crop files do)? If so, use `--propagate`.
4. Which gene list is the universe, if not a full DE table? It must be the genes that were tested, never the genome.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/enrichment.R" --demo --out demo_enrich/
Rscript "${CLAUDE_SKILL_DIR}/scripts/enrichment.R" --de-table de/tables/de_condition_drought_vs_control.csv \
    --annotation Gmax_508_Wm82.a4.v1.annotation_info.txt --gene-col locusName --term-col GO --propagate --out enrich/
Rscript "${CLAUDE_SKILL_DIR}/scripts/enrichment.R" --genes module_turquoise.txt --universe network_genes.txt \
    --annotation kegg_gene2pathway.csv --term-names kegg_pathway_names.csv --out enrich_kegg/
```

| Option | Default | Meaning |
|---|---|---|
| `--ontology` | BP,MF,CC | GO namespaces, each tested separately |
| `--propagate` | off | add is_a/part_of ancestors (GO only) |
| `--direction` | all,up,down | gene sets made from `--de-table` |
| `--alpha` | 0.05 | DE cut-off for `--de-table` and FDR for terms |
| `--lfc-min` | 0 | also require \|log2FoldChange\| ≥ this for DE genes |
| `--min-size` / `--max-size` | 10 / 500 | gene-set size limits within the annotated universe |
| `--min-overlap` | 0.2 | stop if a smaller share of universe or DE genes is in the annotation |
| `--strip-isoform` | off | merge `.1`, `.1.p`, `_T001`, `_P001` IDs to genes |
| `--format`, `--gene-col`, `--term-col` | auto | annotation layout |

## Outputs

- `report.md`: warnings (ID problems, obsolete GO IDs, unpropagated annotation), coverage, results per gene set and ontology.
- `tables/enrichment_<set>_<ontology>.csv`: `term_id`, `term_name`, `ontology`, `n_de_in_term` (k), `n_de_annotated` (n), `term_size` (M), `n_universe_annotated` (N), `gene_ratio`, `bg_ratio`, `fold_enrichment` = (k/n)/(M/N), `pvalue`, `padj`, `significant`, `genes`. Only terms with at least one DE gene are listed.
- `tables/annotation_coverage.csv`: share of universe and DE genes present in the annotation.
- `tables/run_summary.csv`: genes, terms tested and significant terms per analysis.
- `figures/dotplot_<set>_<ontology>.png`: top terms by p-value (filled = significant).

## Expert decisions

- **Universe = tested genes.** Genes that were never expressed, or were removed before testing, could never be called DE. Including them in the background makes any category of well-expressed genes look enriched (Timmons et al. 2015: the background should contain only genes that were measured and called expressed). Wijesooriya et al. (2022) found that 95% of published over-representation analyses used no appropriate background or did not report it. The tests show the effect: a well-expressed term unrelated to DE gets p < 1e-4 against the genome and p > 0.05 against the tested genes.
- **Hypergeometric test** (one-sided, over-representation) as implemented in clusterProfiler (Wu et al. 2021). Only annotated genes enter the test (N and n count annotated genes), which is clusterProfiler's behaviour.
- **BH over every size-eligible term.** clusterProfiler adjusts only over terms that contain at least one DE gene; this skill also counts eligible terms without DE genes (p = 1), which is the complete family of tests and slightly more conservative.
- **Size limits 10-500** (clusterProfiler's defaults, adjustable): tiny sets are unstable and huge sets (e.g. "biological process") are uninformative.
- **GO propagation over is_a and part_of only.** The GO annotation rules say positive annotations propagate up is_a and part_of, not regulates ("involved in regulation of glycolysis" does not mean "involved in glycolysis"; GO documentation). GO.db's `GOBPANCESTOR` map also follows regulates edges (e.g. it lists "gene expression" as an ancestor of "regulation of gene expression"), so this skill walks `GOBPPARENTS`/`GOMFPARENTS`/`GOCCPARENTS` restricted to is_a/part_of instead. Without propagation, broad terms are under-counted when the file lists only the most specific terms; the report measures this.
- **Secondary GO IDs** are mapped to primary IDs (GOSYNONYM); obsolete IDs, or IDs newer than the installed GO.db release, are dropped and counted.
- **ID check before testing:** below 20% overlap (adjustable rule of thumb) the run stops with a diagnosis: isoform suffixes, annotation-version suffixes (`.Wm82.a4.v1`), letter case, `gene:` prefixes, or different ID systems (soybean a1 vs a2+, maize v3/v4/v5, rice MSU vs RAP-DB, wheat 01G/02G/03G codes). Converting between releases needs an official correspondence table (SoyBase, MaizeGDB, Gramene, Phytozome); string edits are not a valid conversion.

### Validation (tests/test_enrichment.R)

- p-values equal `phyper(k - 1, M, N - M, n, lower.tail = FALSE)` and BH values computed by hand; size limits exclude the right sets.
- Planted enrichment (3 sets, 40% of members DE) is recovered at padj < 0.001; at most 1 of 47 null sets passes padj < 0.05 (0 in the test data).
- The universe effect described above, with p-values matching hand calculations.
- Propagation on a hand-built DAG adds is_a and part_of ancestors but not regulates targets; on GO.db, "regulation of gene expression" reaches the BP root but not "gene expression".
- The demo recovers exactly the planted drought and light-harvesting terms, and no null term.

## Checks before reporting

- Read the Warnings section first. **Stop and tell the user** if the run stopped on ID overlap, if fewer than about half of the tested genes have any annotation, or if obsolete GO IDs are numerous (old annotation file).
- Check `annotation_coverage.csv`: DE genes much less annotated than the universe (e.g. many novel or lineage-specific genes) biases results toward well-studied families.
- With fewer than about 10 annotated DE genes, say that the analysis has almost no power.
- Report fold enrichment and the number of DE genes in each term, not only p-values. Terms that share most genes (parents and children after propagation) are one finding, not several.
- Enrichment is evidence about gene sets, not about individual genes or causality; the background and annotation release must be stated in any write-up.
- Don't test GO terms for a species against another species' annotation unless the user explicitly chose ortholog-based annotation, and then say so.

## References

- Wu T, Hu E, Xu S, et al. (2021). clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. The Innovation 2:100141. https://doi.org/10.1016/j.xinn.2021.100141
- Ashburner M, Ball CA, Blake JA, et al. (2000). Gene Ontology: tool for the unification of biology. Nature Genetics 25:25-29. https://doi.org/10.1038/75556
- The Gene Ontology Consortium (2023). The Gene Ontology knowledgebase in 2023. Genetics 224:iyad031. https://doi.org/10.1093/genetics/iyad031
- Gene Ontology documentation: annotation propagation rules, https://geneontology.org/docs/go-annotations/ and relations, https://geneontology.org/docs/ontology-relations/
- Timmons JA, Szkop KJ, Gallagher IJ (2015). Multiple sources of bias confound functional enrichment analysis of global -omics data. Genome Biology 16:186. https://doi.org/10.1186/s13059-015-0761-7
- Wijesooriya K, Jadaan SA, Perera KL, Kaur T, Ziemann M (2022). Urgent need for consistent standards in functional enrichment analysis. PLOS Computational Biology 18:e1009935. https://doi.org/10.1371/journal.pcbi.1009935
- Young MD, Wakefield MJ, Smyth GK, Oshlack A (2010). Gene ontology analysis for RNA-seq: accounting for selection bias. Genome Biology 11:R14. https://doi.org/10.1186/gb-2010-11-2-r14
- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society B 57:289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
