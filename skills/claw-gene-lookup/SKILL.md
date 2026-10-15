---
name: claw-gene-lookup
description: Find the genes inside GWAS or QTL intervals and list their IDs, coordinates, strands, descriptions and distance to the peak SNP. Genes come from the Ensembl REST API, which serves Ensembl Plants genomes (soybean glycine_max, maize zea_mays, wheat triticum_aestivum, rice oryza_sativa, sorghum sorghum_bicolor, Arabidopsis arabidopsis_thaliana and others), or from a local GFF3 annotation for any assembly. The skill checks that the coordinates and the annotation use the same genome assembly, matches chromosome names (Gm01/Chr01/1), and converts 0-based BED coordinates. Use when the user asks which genes are near a significant SNP or in a QTL region, for candidate genes, gene annotation of GWAS hits, genes within an LD window, or to annotate qtl_regions.csv from claw-gwas-crop.
license: MIT
compatibility: R >= 4.2 with httr2 and jsonlite (Ensembl mode needs internet access to rest.ensembl.org; GFF3 mode runs offline). openssl optional for SHA-256 checksums.
metadata:
  version: "1.0.0"
  domain: genomics
  language: R
  status: tested
  summary: "Genes inside GWAS or QTL intervals from Ensembl Plants or a local GFF3, with assembly checks"
  repository: https://github.com/bgtamang/AgriClaw
---

# Genes in GWAS and QTL intervals

## What it does

Takes intervals (explicit start and end, or peak SNPs ± a window) and returns every annotated gene that overlaps them. Overlap means gene start ≤ interval end and gene end ≥ interval start, in 1-based inclusive coordinates. For each gene it reports the ID, name, coordinates, strand, biotype, description and the distance to the peak SNP. Genes come from Ensembl REST or a local GFF3 file. The genome assembly of the user's coordinates is checked against the gene source before anything is looked up.

## When to use / when not to use

- Use after `claw-gwas-crop` (pass `tables/qtl_regions.csv` directly) or `claw-qtl-mapping` (support intervals), or for a list of significant SNPs.
- Take the window around a SNP from `claw-ld-decay` (`--ld-table tables/ld_window.csv`), not from a remembered rule.
- Don't use it to judge which gene is causal. It lists positional candidates only. Functional evidence (expression, mutants, orthologues) is a separate step; `claw-enrichment` tests gene-set enrichment.
- Don't use the Ensembl mode for assemblies Ensembl does not serve. For example, Ensembl Plants serves soybean Glycine_max_v2.1 only, not Wm82.a4 or a6. Use `--gff` with the matching annotation (for example from Phytozome, SoyBase or MaizeGDB) instead.

## Inputs

- Regions, one of:
  - `--intervals FILE`: CSV with `chr`, `start`, `end` (or `start_bp`, `end_bp`), plus optional `region_id`/`name`, `lead_snp`, `peak_pos`/`lead_pos_bp` and `trait`. CSV coordinates are 1-based and inclusive; `--coords bed` declares 0-based starts. A `.bed` file is always read as 0-based, half-open, and its start is shifted by +1.
  - `--snps FILE`: CSV with `snp`, `chr`, `pos`, plus `--window-bp N` or `--ld-table FILE`.
- Genes, one of:
  - `--gff FILE`: a GFF3 file (1-based, inclusive). `--feature-type gene` by default; use `gene,pseudogene` to include pseudogenes.
  - `--ensembl-species NAME`: an Ensembl name such as `glycine_max`, `zea_mays`, `triticum_aestivum`, `oryza_sativa`, `sorghum_bicolor` or `arabidopsis_thaliana`.
- `--assembly NAME` (required): the assembly the user's positions come from.

**Ask the user (never guess):**

- **Which genome assembly** the SNP or interval positions use: ask which reference the reads or markers were mapped to. Examples:
  - soybean: Wm82.a2, a4 or a6;
  - maize: B73 RefGen_v4 or Zm-B73-REFERENCE-NAM-5.0;
  - wheat: IWGSC RefSeq v1.0 or v2.1;
  - rice: IRGSP-1.0.
- **The window**, if SNPs rather than intervals are given: use the LD distance from `claw-ld-decay`.
- **Whether a BED-style 0-based file** is involved, if the intervals did not come from an AgriClaw skill.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/gene_lookup.R" --demo
# GWAS regions on soybean Wm82.a2 coordinates, genes from Ensembl Plants
Rscript "${CLAUDE_SKILL_DIR}/scripts/gene_lookup.R" --intervals gwas/tables/qtl_regions.csv \
        --ensembl-species glycine_max --assembly Glycine_max_v2.1 --out genes/
# peak SNPs with per-chromosome LD windows and a local annotation (any assembly)
Rscript "${CLAUDE_SKILL_DIR}/scripts/gene_lookup.R" --snps hits.csv --ld-table ld/tables/ld_window.csv \
        --gff Gmax_Wm82.a4.v1.gene.gff3.gz --assembly Wm82.a4.v1 --out genes_a4/
```

The demo uses a synthetic 3-chromosome GFF3 and three QTL regions written in the `claw-gwas-crop` format, with chromosomes named `1` to `3` against the GFF3's `Chr01` to `Chr03`. It runs offline.

## Outputs

- `report.md`: warnings, then the assembly and gene source, notes on coordinate conversion and chromosome renaming, a per-interval gene count and the nearest gene to each peak.
- `tables/genes_in_intervals.csv`: one row per gene and interval. Columns:
  - `interval_id`, `trait`, `interval_chr`, `interval_start_bp`, `interval_end_bp`, `lead_snp`, `peak_pos_bp`;
  - `gene_id`, `gene_name`, `chr` (the source's chromosome name), `gene_start_bp`, `gene_end_bp`, `strand`, `biotype`, `description`;
  - `distance_to_peak_bp` (0 when the peak lies inside the gene, otherwise the gap to the nearer gene end) and `peak_in_gene`.

  Rows are sorted by distance to the peak.
- `tables/intervals.csv`: the intervals actually searched, after conversion and clipping, with gene counts.
- `tables/chromosome_map.csv`: the user's chromosome names mapped to the source's names, and how each was matched.
- `figures/gene_tracks.png`: genes by strand for up to 12 intervals, with the peak marked.
- `raw_responses/` (Ensembl mode): every JSON response, for reproducibility.

## Expert decisions

- **Assembly check: refuse to mix.**
  - **Ensembl mode:** `--assembly` must equal Ensembl's `assembly_name` or GenBank accession for the species (from `/info/assembly`), ignoring case and punctuation. Otherwise the run stops and names the served assembly.
    - On 2026-10-15 (release 116) Ensembl served `glycine_max` Glycine_max_v2.1 (GCA_000004515.4), annotated with the JGI Wm82.a2.v1 gene set. Its gene IDs are ENA locus tags such as `GLYMA_01G010400`, where Phytozome and SoyBase write `Glyma.01G010400`.
    - The other species: `zea_mays` Zm-B73-REFERENCE-NAM-5.0, `triticum_aestivum` IWGSC, `oryza_sativa` IRGSP-1.0, `sorghum_bicolor` Sorghum_bicolor_NCBIv3, `arabidopsis_thaliana` TAIR10. Releases change, so trust the run's report, not this list.
  - **GFF3 mode:** a `#!genome-build`, `##genome-build` or `##annot-version` header must match `--assembly`. Without such a header the run continues with a warning, and the agent must confirm the file's assembly with the user.
  - **Coordinate checks:** a position past the end of a chromosome stops the run, because it means another assembly. A window that only runs past the end is clipped.
- **Same gene IDs do not mean same coordinates.** SoyBase carries Wm82.a2 gene-model names forward to Wm82.a4 and later assemblies where applicable, so an identical soybean ID tells you nothing about which coordinates are meant.
- **Maize gene IDs encode the assembly.** `Zm00001d…` IDs are B73 RefGen_v4 (the Ensembl archive lists them under B73_RefGen_v4); `Zm00001eb…` IDs are Zm-B73-REFERENCE-NAM-5.0 (MaizeGDB). A clear conflict with `--assembly` stops the run.
- **Chromosome names.** Exact names are used first. Otherwise names are matched after removing `chr`/`chromosome`/`Gm`/`Zm`/`Os`/`Sb`/`At`/`Ta`/`Sl` prefixes and leading zeros. So Gm01, Chr01, chr1 and 1 match one another, while wheat 1A, 1B and 1D stay distinct. Ambiguous or unmatched names stop the run (use `--chr-map`); every renaming is listed in the report.
- **Coordinates.** GFF3 and VCF are 1-based and inclusive; BED is 0-based and half-open (UCSC). BED starts get +1, BED ends are unchanged.
- **Ensembl limits.** The overlap endpoint accepts at most 5 Mb per request, so longer intervals are split. Requests are throttled to 10 per second (Ensembl allows 55,000 per hour, about 15 per second) and retried on 429, 500, 502, 503 and 504, all of which occurred during development. The Ensembl source tag is removed from descriptions.
- **More than 300 genes in one interval** (rule of thumb) triggers a warning: the window is too wide to suggest candidates.

## Checks before reporting

- Confirm that the assembly line of `report.md` is what the user's positions came from. If the GFF3 had no build header, say so explicitly.
- Read `chromosome_map.csv`. A mapping such as `1 → Chr01` is expected; anything surprising (scaffolds, `Un`, organelles) needs the user's confirmation.
- **Stop and tell the user** when:
  - an interval has no genes (wrong chromosome naming, a gene desert, or a window that is too small);
  - an interval has hundreds of genes;
  - a mismatch error appears: never re-run with a different `--assembly` just to make it pass;
  - descriptions are mostly empty (common for soybean in Ensembl): offer the GFF3 or another annotation source.
- Present genes as positional candidates ranked by distance only as a convenience. LD may span the whole interval, and the nearest gene is often not the causal one.

## References

- Yates A, Beal K, Keenan S, McLaren W, Pignatelli M, Ritchie GRS, Ruffier M, Taylor K, Vullo A, Flicek P (2015; online 2014) The Ensembl REST API: Ensembl data for any language. Bioinformatics 31:143-145. https://doi.org/10.1093/bioinformatics/btu613
- Yates AD et al. (2022; online 2021) Ensembl Genomes 2022: an expanding genome resource for non-vertebrates. Nucleic Acids Research 50:D996-D1003. https://doi.org/10.1093/nar/gkab1007
- Schmutz J et al. (2010) Genome sequence of the palaeopolyploid soybean. Nature 463:178-183. https://doi.org/10.1038/nature08670
- Valliyodan B et al. (2019) Construction and comparison of three reference-quality genome assemblies for soybean. The Plant Journal 100:1066-1082. https://doi.org/10.1111/tpj.14500
- Jiao Y et al. (2017) Improved maize reference genome with single-molecule technologies. Nature 546:524-527. https://doi.org/10.1038/nature22971
- Hufford MB et al. (2021) De novo assembly, annotation, and comparative analysis of 26 diverse maize genomes. Science 373:655-662. https://doi.org/10.1126/science.abg5289
- Ensembl REST API rate limits. https://github.com/Ensembl/ensembl-rest/wiki/Rate-Limits
- Sequence Ontology. GFF3 specification. https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
- UCSC Genome Browser. BED format (0-based starts). https://genome.ucsc.edu/FAQ/FAQformat.html#format1
- Ensembl Plants soybean annotation page (Glycine_max_v2.1, Wm82.a2.v1 gene set). https://plants.ensembl.org/Glycine_max/Info/Annotation/
- SoyBase genome nomenclature ("Gene model names from Wm82.a2 assembly will be carried forward to the Wm82.a4 and later assemblies where applicable"). https://legacy.soybase.org/nomenclature/aboutgenomenomenclature.php
- MaizeGDB, Zm-B73-REFERENCE-NAM-5.0 README (official gene set Zm00001eb.1). https://download.maizegdb.org/Genomes/NAM_Founders/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.README
