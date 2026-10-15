---
name: claw-rnaseq-import
description: Import Salmon (quant.sf), kallisto (abundance.tsv) or featureCounts RNA-seq quantifications into gene-level counts with tximport, build the transcript-to-gene table (tx2gene) from a GFF3/GTF or a table, catch transcript ID mismatches (version suffixes, Phytozome vs Ensembl ID styles, wrong annotation release), and run sample QC - mapping rate, library size, detected genes, sample correlations, PCA on VST and outlier flags. Built for crops and non-model plants (soybean, maize, wheat, rice) whose annotation comes from GFF3 files rather than OrgDb packages. Use when the user has Salmon or kallisto output folders, quant.sf files or a featureCounts matrix and asks to import them, summarise transcripts to genes, make a gene count matrix, check RNA-seq sample quality, or find outlier samples before differential expression.
license: MIT
compatibility: R >= 4.3 with Bioconductor tximport, DESeq2 and SummarizedExperiment, plus jsonlite (CRAN); readr optional for faster reading. Tested with R 4.6.1, tximport 1.40, DESeq2 1.52. No internet needed.
metadata:
  version: "1.0.0"
  domain: transcriptomics
  language: R
  status: tested
  summary: "Import Salmon, kallisto or featureCounts output to gene counts, with sample QC"
  repository: https://github.com/bgtamang/AgriClaw
---

# RNA-seq import and sample QC

## What it does

Reads per-sample transcript quantifications from Salmon or kallisto (or a featureCounts gene matrix), maps transcripts to genes with the same GFF3/GTF annotation that built the index, and summarises them to genes with tximport. It keeps the estimated counts plus average transcript lengths, which claw-de-analysis passes to DESeq2 as offsets. It stops when transcript IDs don't match the annotation, then reports per-sample QC (mapping rate, library size, detected genes, correlations, PCA) and flags outlier samples.

## When to use / when not to use

- Use after Salmon, kallisto or featureCounts and before any differential expression or co-expression analysis.
- Use to check that a GFF3/GTF and a transcriptome index match (same annotation release).
- Don't use for read alignment or quantification from FASTQ: run Salmon/kallisto or an established pipeline such as nf-core/rnaseq first.
- Differential expression: **claw-de-analysis** (takes this skill's `txi.rds` and `tables/samples.csv`).
- Co-expression modules: **claw-coexpression**. GO/KEGG enrichment: **claw-enrichment**.
- Not for single-cell data or differential transcript usage (needs transcript-level methods).

## Inputs

**Sample sheet** (CSV, `--samples`), one row per library:

| Column | Meaning |
|---|---|
| `sample` | unique sample name (used as the column name everywhere downstream) |
| `path` | Salmon `quant.sf` or kallisto `abundance.tsv`, or the folder that holds it; relative paths are read relative to the sample sheet's folder. Not needed for featureCounts. |
| `column` | featureCounts only, optional: the matrix header for this sample (otherwise matched by BAM file name) |
| other columns | design variables, e.g. `condition`, `genotype`, `batch`, `time`; copied to `tables/samples.csv` |

Salmon writes `aux_info/meta_info.json` and kallisto writes `run_info.json` next to the quant file; keep the whole output folder so mapping rates can be read.

**Annotation** (one of):
- `--gff`: GFF3 or GTF (plain or `.gz`). For GFF3 the script finds features whose `Parent` is a gene (mRNA, transcript, lnc_RNA, ncRNA, ...) and tries the `ID` (with `gene:`/`transcript:` prefixes removed), `Name`, `transcript_id` and raw `ID` attributes, keeping the one that matches the quantified IDs. For example, Phytozome GFF3 transcripts carry an `ID` with an annotation-version suffix (`Glyma.01G000100.1.Wm82.a4.v1`) and a plain `Name` (`Glyma.01G000100.1`); which one matches depends on how the transcript FASTA for the index was made, so both are tried and the match rate of each is reported.
- `--tx2gene`: CSV/TSV with columns `transcript_id`, `gene_id` (or the first two columns), with a header row.

**Ask the user, never guess:**
1. Which genome assembly and annotation release built the Salmon/kallisto index (e.g. soybean Wm82.a4.v1 vs Wm82.a2.v1, maize B73 v5 vs v4, wheat IWGSC RefSeq v1.1 vs v2.1)? Use the GFF3/GTF from that exact release.
2. Which sample-sheet column is the treatment, which are batch or block, and how replicates are defined.
3. Is the library full-length (standard) or 3' tag-seq (e.g. QuantSeq)? 3' tag-seq needs `--three-prime`.
4. Were organelle genomes, spike-ins (ERCC) or extra sequences added to the index? They are expected to be unmatched.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/rnaseq_import.R" --demo --out demo_import/
Rscript "${CLAUDE_SKILL_DIR}/scripts/rnaseq_import.R" --samples samples.csv --gff Gmax_508_Wm82.a4.v1.gene.gff3.gz --out rnaseq_import/
Rscript "${CLAUDE_SKILL_DIR}/scripts/rnaseq_import.R" --samples samples.csv --tx2gene tx2gene.csv --out rnaseq_import/
Rscript "${CLAUDE_SKILL_DIR}/scripts/rnaseq_import.R" --samples samples.csv --type featurecounts --counts featureCounts.txt --fc-summary featureCounts.txt.summary
```

| Option | Default | Meaning |
|---|---|---|
| `--type` | auto | `salmon`, `kallisto` or `featurecounts`; auto-detected from file names |
| `--id-attr` | auto | force the GFF3 attribute for transcript IDs: `ID`, `Name`, `transcript_id` |
| `--strip-tx-version` | off | remove a trailing `.N` version on both sides; refused when `.N` is an isoform number |
| `--min-id-match` | 0.95 | stop if fewer quantified transcripts than this fraction are in tx2gene |
| `--min-mapping-rate` | 50 | flag samples below this mapping rate (%) |
| `--outlier-mad` | 3 | flag samples this many MADs below the median |
| `--detect-min-count` | 10 | reads for a gene to count as detected |
| `--transform` | vst | `vst` or `rlog` for QC plots |
| `--three-prime` | off | 3' tag-seq: no transcript-length offsets |
| `--color-by` | first design column | PCA colouring |

The demo (6 synthetic Salmon samples, Phytozome-style GFF3) shows a planted low-mapping sample (S6) and two chloroplast transcripts missing from the nuclear GFF3.

## Outputs

- `report.md`: warnings first, then the ID-matching summary, the QC table and the next step.
- `txi.rds`: the gene-level tximport list (`counts`, `abundance`, `length`, `countsFromAbundance = "no"`); input for claw-de-analysis (`--txi`). Not written for featureCounts.
- `tables/gene_counts.csv`: estimated counts (gene × sample; non-integer for Salmon/kallisto, which is expected).
- `tables/gene_lengths.csv`: average transcript length per gene and sample (the offsets).
- `tables/gene_tpm.csv`: TPM, for plotting expression only, **never for differential expression**.
- `tables/samples.csv`: the sample sheet plus `qc_flagged` and `qc_flags`.
- `tables/sample_qc.csv`: `mapping_rate_pct`, `assigned_reads`, `detected_genes`, `median_cor`, `PC1`, `PC2`, `flags`.
- `tables/tx2gene.csv` (mapping used), `tables/unmatched_transcripts.csv` (if any), `tables/sample_correlation.csv`.
- `figures/qc_overview.png`, `figures/pca.png`, `figures/sample_correlation.png`.
- `reproducibility/`: command, sessionInfo and checksums.

## Expert decisions

- **Counts plus offsets, not TPM.** tximport is run with `countsFromAbundance = "no"`: gene counts are sums of transcript estimated counts, and gene length is the abundance-weighted mean of transcript effective lengths. DESeq2 uses these lengths as offsets that correct for changes in isoform composition (Soneson et al. 2015; DESeq2 vignette). TPM and FPKM are not counts and must not go into DESeq2 (Zhao et al. 2020).
- **Zero-abundance genes** get the geometric mean of their lengths in the other samples (tximport's rule), so offsets stay defined.
- **ID matching is exact.** tximport silently drops transcripts missing from tx2gene; this skill reports how many and what share of reads they held, and stops below `--min-id-match` 0.95 (adjustable rule of thumb). Expected unmatched IDs are organelle or spike-in sequences.
- **Never use tximport's `ignoreTxVersion` on plant IDs that contain dots.** It deletes everything after the first dot, so `Glyma.01G000100.1` becomes `Glyma` and all soybean (and sorghum, poplar: `Sobic.`, `Potri.`) transcripts collapse into one gene. `--strip-tx-version` removes only a trailing `.N` and refuses when that merges isoforms. In soybean, Arabidopsis, rice or wheat IDs, `.1`/`.2` is an isoform number, not a version.
- **Nomenclature check.** When IDs don't match, the error names the ID systems detected (e.g. soybean `Glyma01g00010` is Wm82.a1, `Glyma.01G000100` is Wm82.a2 or later; maize `GRMZM2G...` v3, `Zm00001d` v4, `Zm00001eb` v5; rice MSU `LOC_Os` vs RAP-DB `Os01g0100100`). Different releases need an official cross-reference table; string edits are not valid conversions.
- **QC transform:** DESeq2 `vst(blind = TRUE)`, as the DESeq2 vignette recommends for sample QA; with fewer than 1000 genes of mean count > 5, `varianceStabilizingTransformation` is used instead, because `vst()` subsamples 1000 genes. PCA uses the 500 most variable genes, as `DESeq2::plotPCA` does.
- **Detected gene:** ≥ 10 estimated reads, the count the DESeq2 vignette uses for pre-filtering (adjustable rule of thumb).
- **Mapping-rate floor 50%:** adjustable rule of thumb. Conesa et al. (2016) expect 70-90% of reads to align to a genome and somewhat less to a transcriptome; plant transcriptomes with incomplete annotation often map 60-85%.
- **Outliers:** more than 3 MADs below the median for mapping rate, log10 library size, detected genes or median between-sample correlation. The MAD is used instead of the SD because it is not distorted by the outlier itself (Leys et al. 2013); the multiplier 3 is an adjustable rule of thumb. With fewer than 4 samples the rule has little power.
- **Flags are not removals.** No sample is dropped by this skill.
- **3' tag-seq** (`--three-prime`): counts have no length bias, so no length offsets are used (tximport vignette).

### Validation (tests/test_rnaseq_import.R)

- Synthetic Salmon and kallisto files (6 samples, about 240 transcripts from 120 soybean-style genes): gene counts and TPM equal hand sums (to 1e-8), gene lengths equal TPM-weighted mean effective lengths, and a zero-abundance gene gets tximport's geometric-mean length.
- A planted sample with 31% mapping is flagged by both the absolute and the MAD rule; no other sample gets a mapping flag.
- Version-suffix mismatches are caught in both directions and fixed by `--strip-tx-version`; stripping is refused when it would merge isoforms (`Gm.01G000100.1`/`.2`); soybean Wm82.a1 vs a2+ IDs are named in the diagnosis.
- Phytozome-style GFF3 (Name chosen automatically), gzipped Ensembl-style GFF3 (prefixes removed, ncRNA genes kept), GTF, URL-escaped attributes and featureCounts matrices are parsed correctly.
- The demo flags only S6, reports the two organelle transcripts, and its gene counts equal hand sums of the quant.sf files.

## Checks before reporting

- Read the Warnings section of `report.md` first and pass every warning to the user.
- **Stop and tell the user** when:
  - the script stops on ID matching: show the diagnosis and ask which annotation release built the index; don't edit IDs by hand to force a match;
  - samples were quantified with different indices (the script stops), or Salmon reports different library types across samples;
  - most samples are flagged, or every sample maps below about 50%: suspect the wrong reference species or release, rRNA or other contamination, or a strandedness problem;
  - unmatched transcripts hold more than a few percent of reads and are not organelle or spike-in sequences.
- For each flagged sample, look at `figures/pca.png` and `figures/sample_correlation.png`: does it cluster with its own group? Ask the user for lab notes (RNA quality, sequencing lane) before anyone removes it. If removed, re-run claw-de-analysis with `--exclude SAMPLE` and report the removal.
- If PC1 or PC2 separates sequencing batches or lanes rather than treatments, tell the user and include the batch column in the DE design.
- Assigned reads well below about 5-10 million per sample for full-length libraries reduce power (rule of thumb); say so.

## References

- Soneson C, Love MI, Robinson MD (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research 4:1521. https://doi.org/10.12688/f1000research.7563.1
- Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods 14:417-419. https://doi.org/10.1038/nmeth.4197
- Bray NL, Pimentel H, Melsted P, Pachter L (2016). Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology 34:525-527. https://doi.org/10.1038/nbt.3519
- Liao Y, Smyth GK, Shi W (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics 30:923-930. https://doi.org/10.1093/bioinformatics/btt656
- Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15:550. https://doi.org/10.1186/s13059-014-0550-8
- Conesa A, Madrigal P, Tarazona S, et al. (2016). A survey of best practices for RNA-seq data analysis. Genome Biology 17:13. https://doi.org/10.1186/s13059-016-0881-8
- Zhao S, Ye Z, Stanton R (2020). Misuse of RPKM or TPM normalization when comparing across samples and sequencing protocols. RNA 26:903-909. https://doi.org/10.1261/rna.074922.120
- Leys C, Ley C, Klein O, Bernard P, Licata L (2013). Detecting outliers: do not use standard deviation around the mean, use absolute deviation around the median. Journal of Experimental Social Psychology 49:764-766. https://doi.org/10.1016/j.jesp.2013.03.013
- tximport and DESeq2 vignettes (Bioconductor): https://bioconductor.org/packages/tximport/ and https://bioconductor.org/packages/DESeq2/
