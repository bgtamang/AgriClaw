---
name: claw-marker-qc
description: Quality control and format conversion for crop SNP genotype data from VCF, HapMap (TASSEL) or a 0/1/2 marker matrix. Reports marker call rate, minor allele frequency (MAF), heterozygosity and monomorphic markers; sample call rate, heterozygosity outliers (DNA contamination, seed mix or outcrossing in inbred lines) and duplicate or near-identical samples by identity-by-state (IBS); segregation distortion (chi-square; RIL/DH 1:1, F2 1:2:1, backcross 1:1) and parent consistency for biparental mapping populations. Writes a filtered 0/1/2 ALT-allele-count matrix (optionally VCF) ready for GWAS, LD decay, QTL mapping or genomic prediction. Use when the user has SNP-chip, GBS or sequencing genotypes and asks to clean, filter or QC markers, check missing data, find duplicated or mislabelled lines, check heterozygosity or off-types, test segregation ratios, or convert VCF/HapMap to a 0/1/2 matrix.
license: MIT
compatibility: R >= 4.2 with the vcfR package (openssl optional, for SHA-256 checksums). No internet needed.
metadata:
  version: "1.0.0"
  domain: genomics
  language: R
  status: tested
  repository: https://github.com/bgtamang/AgriClaw
---

# Marker QC for crop genotype data

## What it does

Reads a VCF, HapMap or 0/1/2 genotype matrix and computes per-marker statistics (call rate, ALT and minor allele frequency, observed heterozygosity, monomorphic markers) and per-sample statistics (call rate, heterozygosity, near-identical pairs by identity-by-state). For biparental populations it also tests segregation ratios and flags progeny that carry alleles neither parent has. It writes a filtered matrix coded as the number of ALT-allele copies (0/1/2), plus QC tables, plots and a report.

## When to use / when not to use

- Use first, before `claw-gwas-crop`, `claw-ld-decay`, `claw-qtl-mapping` or `claw-genomic-prediction`. Those skills read this skill's `genotypes_filtered_012.csv`.
- Use to find duplicated, mislabelled, contaminated or outcrossed lines, and markers with too much missing data.
- Don't use for read alignment or variant calling: start from a VCF made by an established pipeline (for example nf-core).
- Don't use for polyploid dosage data (for example autotetraploid potato with 0 to 4 dosages), multi-allelic markers or imputation. Multi-allelic sites are excluded and counted; missing calls stay missing.

## Inputs

| Format | Detected by | Notes |
|---|---|---|
| VCF (`.vcf`, `.vcf.gz`) | `##fileformat=VCF` | GT field; diploid biallelic calls; phased `0\|1` read as unphased |
| HapMap (`.hmp.txt`) | header `rs#  alleles ...` | 11 annotation columns; two-letter (`AG`) or IUPAC (`R`) calls; counts copies of the **second** allele in `alleles` |
| 0/1/2 matrix (`.csv`, `.tsv`, `.gz`) | anything else | codes 0/1/2 or -1/0/1; empty, `NA`, `N` = missing; hard calls only |

A matrix may have markers in rows or samples in rows. The script detects the layout from `chr`/`pos` columns, from IDs listed in `--map`, or from marker-like IDs (such as `S1_12345`, `Gm01_123456`, `ss715578788`). If it can't tell, it stops and asks for `--orientation`.

**Ask the user before running:**

- **Germplasm (required, never guess):** inbred lines (RILs, DH, pure-line cultivars of soybean, wheat, rice and similar) or outbred material (open-pollinated maize, hybrids, clones, F2 individuals)? The heterozygosity checks depend on it.
- For a mapping population: which type (RIL, DH, F2 or BC1), and which sample IDs are the parents? For BC1, which parent is the recurrent one?
- Whether flagged samples should be dropped or only reported. The default is to report them and keep them.
- Which genome assembly the VCF positions use. Record it in the report, because downstream skills (`claw-gene-lookup`) must not mix assemblies.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/marker_qc.R" --demo
Rscript "${CLAUDE_SKILL_DIR}/scripts/marker_qc.R" --input panel.vcf.gz --germplasm inbred --out qc/
Rscript "${CLAUDE_SKILL_DIR}/scripts/marker_qc.R" --input geno.csv --map map.csv --germplasm inbred \
        --orientation samples-in-rows --min-marker-call-rate 0.8 --out qc/
Rscript "${CLAUDE_SKILL_DIR}/scripts/marker_qc.R" --input rils.hmp.txt --germplasm inbred \
        --population ril --parents Williams82,PI437654 --out qc/
```

The demo is a synthetic inbred panel with 120 lines and 600 markers. Planted problems: 20 high-missing markers, 10 monomorphic, 10 rare, 4 low-call-rate samples, 3 contaminated samples and 3 duplicate pairs. All are found exactly.

Useful options: `--drop-duplicates` (keep the best-called sample per duplicate group), `--drop-het-outliers`, `--drop-distorted`, `--drop-nonparental`, `--write-vcf`, `--seed`. Run `--help` for the full list.

## Outputs

- `report.md`: warnings first (failed and flagged samples, duplicate pairs, distorted markers), then thresholds, marker and sample summaries.
- `tables/genotypes_filtered_012.csv`: markers in rows with `marker, chr, pos, ref, alt`, then one column per sample. Values are copies of the ALT allele; empty means missing. This is the input format for `claw-gwas-crop` and `claw-ld-decay`.
- `tables/marker_qc.csv`: `call_rate`, `alt_freq`, `maf`, `het_obs`, `monomorphic`, `fail_call_rate`, `fail_maf`, `pass`, `fail_reason`. With `--population` it adds `seg_n_A/H/B`, `seg_chisq`, `seg_p` and `distorted`; with `--parents`, `parent1`, `parent2`, `parents_polymorphic` and `nonparental_call_rate`.
- `tables/sample_qc.csv`: `call_rate`, `het_obs`, `fail_call_rate`, `het_flag`, `duplicate_group`, `keep`, `fail_reason`. With `--parents` it adds `nonparental_rate` and `nonparental_flag`.
- `tables/duplicate_pairs.csv`: `sample_1`, `sample_2`, `ibs`, `n_shared_markers`.
- `tables/qc_summary.csv`, `figures/marker_qc.png`, `figures/sample_qc.png`, `figures/ibs_pairs.png` and `figures/segregation.png` (with `--population`). `genotypes_filtered.vcf` is written with `--write-vcf`.

## Expert decisions

Order of operations, with sample QC first as in Anderson et al. (2010):
1. sample call rate on all markers;
2. marker statistics on the samples that pass;
3. sample heterozygosity and duplicates on the markers that pass.

| Setting | Default | Why | Source |
|---|---|---|---|
| marker call rate | ≥ 0.90 | High-missing markers are usually failed assays. Lower it (0.5 to 0.8) for GBS data that will be imputed. | adjustable rule of thumb |
| sample call rate | ≥ 0.80 | Low-call samples are usually poor DNA; they also bias heterozygosity. | adjustable rule of thumb |
| MAF | ≥ 0.05 | Rare alleles give unstable estimates in small panels: in 200 lines, MAF 0.05 is about 10 carriers. | adjustable rule of thumb |
| inbred heterozygosity | flag > 0.10 | Residual heterozygosity of an F_t-derived line is (1/2)^(t-1): F5 6.3%, F6 3.1%, F8 0.8%. A DNA mix of two lines, or an F1 outcross, is heterozygous at most markers where the lines differ. Raise it for early-generation material. | Mendelian expectation; threshold is a rule of thumb |
| outbred heterozygosity | flag outside mean ± 3 SD | Excess suggests contamination; deficit suggests inbreeding or selfing. | Anderson et al. 2010 |
| duplicates | IBS ≥ 0.98 on ≤ 5000 random passing markers (seeded) | Allows for genotyping error. Sister lines, sub-lines and near-isogenic lines can legitimately exceed 0.98, so pairs are reported for review. | adjustable rule of thumb |
| segregation distortion | chi-square, Bonferroni family-wise α = 0.05 | RIL and DH 1:1 (heterozygotes excluded and counted); F2 1:2:1; BC1 1:1 (heterozygote : recurrent parent). Untested: markers monomorphic between parents, or in the progeny when no parents are given. | Broman & Sen 2009 (data checking) |
| non-parental calls | flag progeny > 0.05 | Uses markers where both parents are homozygous for the same allele: any other call is non-parental. High rates mean seed mix, outcross or mislabel. | adjustable rule of thumb |
| Hardy-Weinberg test | not applied | Inbred lines and structured crop panels are not in Hardy-Weinberg proportions, so the test would flag good markers. | — |

- Distorted markers are **flagged, not dropped**, unless the user asks (`--drop-distorted`). Distortion can be real biology, such as gametophytic selection, and QTL mapping can even benefit from it (Xu 2008).
- Coding: ALT allele count for VCF; second listed allele for HapMap; -1/0/1 matrices are shifted to 0/1/2 with a note in the report. Matrices of imputed (non-integer) dosages are refused, because heterozygosity and call-rate QC need hard calls.

## Checks before reporting

- Read the Warnings section of `report.md` first and pass every item to the user.
- **Stop and ask** when:
  - more than about 10% of samples fail call rate: suspect a plate or batch failure; check `sample_qc.csv` for a pattern by plate;
  - many samples in an "inbred" panel exceed the heterozygosity limit: the material may not be inbred (early generation, hybrids), so confirm `--germplasm`;
  - duplicate pairs appear: ask whether they are the same cultivar under two names, sister lines or a sample swap; never merge or drop them silently;
  - a parent is heterozygous or missing at many markers, or most markers look distorted: the wrong parents or population type may have been given;
  - fewer than a few hundred markers remain, or the map columns are empty (GWAS and LD need chromosome and position).
- Report how many samples and markers went in and came out, and list flagged sample IDs. Don't just say "QC passed".

## References

- Anderson CA, Pettersson FH, Clarke GM, Cardon LR, Morris AP, Zondervan KT (2010) Data quality control in genetic case-control association studies. Nature Protocols 5:1564-1573. https://doi.org/10.1038/nprot.2010.116
- Broman KW, Sen Ś (2009) A Guide to QTL Mapping with R/qtl. Springer, New York. https://doi.org/10.1007/978-0-387-92125-9
- Bradbury PJ, Zhang Z, Kroon DE, Casstevens TM, Ramdoss Y, Buckler ES (2007) TASSEL: software for association mapping of complex traits in diverse samples. Bioinformatics 23:2633-2635. https://doi.org/10.1093/bioinformatics/btm308 (HapMap genotype format)
- Danecek P et al. (2011) The variant call format and VCFtools. Bioinformatics 27:2156-2158. https://doi.org/10.1093/bioinformatics/btr330
- Knaus BJ, Grünwald NJ (2017; online 2016) vcfR: a package to manipulate and visualize variant call format data in R. Molecular Ecology Resources 17:44-53. https://doi.org/10.1111/1755-0998.12549
- Xu S (2008) Quantitative trait locus mapping can benefit from segregation distortion. Genetics 180:2201-2208. https://doi.org/10.1534/genetics.108.090688
