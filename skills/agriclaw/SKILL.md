---
name: agriclaw
description: Catalog and router for the AgriClaw crop-research skills. Use when a user asks a broad agriculture, plant breeding, field trial, phenotyping, genotyping, GWAS/QTL, RNA-seq, drone or satellite imagery, or plant physiology question and you need to pick which AgriClaw skill(s) to run and in what order, for example trial design, then field book, then trial analysis, then GWAS, then gene lookup.
license: MIT
compatibility: No code of its own. The skills it routes to need R >= 4.2 or Python >= 3.9; see each skill's compatibility field.
metadata:
  version: "1.0.0"
  domain: meta
  language: none
  status: guide
  summary: "Which AgriClaw skill to use for a request, and how the skills chain together"
  repository: https://github.com/bgtamang/AgriClaw
---

# AgriClaw: which skill to use

## What it does

AgriClaw is a set of skills for crop research. Each one wraps established R or Python tools, encodes the decisions experts make, and checks its own results. This page maps common requests to skills and shows how the skills chain together. It runs no code itself.

## When to use

- The request spans several steps (e.g. "analyse my trial and find the genes behind the yield QTL").
- You aren't sure which AgriClaw skill fits.

When the request clearly matches one skill, use that skill directly.

## Ground rules for every AgriClaw skill

1. **Ask; don't guess.** Ask the user for:
   - units and the moisture basis;
   - the experimental design and coordinates;
   - genome assembly and chromosome naming;
   - planting date, crop and maturity group;
   - whether leaves were dark-adapted;
   - the band order of imagery.
2. **Read the warnings first.** Every skill writes `report.md` with warnings at the top. Pass them on; don't bury them.
3. **Use BLUEs between stages.** Pass per-environment BLUEs, not BLUPs, from trial analysis to multi-environment, GWAS, QTL and genomic-prediction skills, to avoid double shrinkage.
4. **Keep the reproducibility bundle.** Each run writes `reproducibility/commands.sh`, session info and checksums. Point the user to them.
5. **Stop on red flags.** Each skill's "Checks before reporting" section lists conditions that mean "tell the user before interpreting".

## Catalog

### Field and environment

| Need | Skill |
|---|---|
| Growing degree days, frost dates, corn growth stages | `claw-gdd` |
| Season weather summary, dry spells, heat or heavy-rain days | `claw-weather` |
| Soil texture, drainage, pH, organic matter, available water (USDA SSURGO) | `claw-soil` |
| Crop water use and water stress (FAO-56 water balance), irrigation timing | `claw-water-balance` |
| bu/ac ↔ kg/ha, moisture adjustment, plot weight to yield | `claw-yield-units` |

### Field trials

| Need | Skill |
|---|---|
| Randomise a trial: RCBD, alpha-lattice, augmented, p-rep, split-plot, strip-plot, Latin square, row-column | `claw-trial-design` |
| Field Book app export/import, merge collected data, data QC | `claw-fieldbook` |
| Analyse one trial: spatial model, BLUEs/BLUPs, heritability, outliers | `claw-trial-analysis` |
| Analyse many environments: GxE, stability, AMMI/GGE, Finlay–Wilkinson | `claw-met-gxe` |

### Genotypes and genetics

| Need | Skill |
|---|---|
| Genotype QC: call rate, MAF, heterozygosity, duplicates, VCF conversion | `claw-marker-qc` |
| Genome-wide association (mixed model with kinship) | `claw-gwas-crop` |
| LD decay, candidate-window size | `claw-ld-decay` |
| Genes under a GWAS/QTL peak (Ensembl Plants or a local GFF3) | `claw-gene-lookup` |
| Biparental or multiparent QTL mapping with permutation thresholds | `claw-qtl-mapping` |
| Genomic prediction (GBLUP, Bayesian) with honest cross-validation | `claw-genomic-prediction` |
| Breeder's equation, selection intensity, Smith–Hazel and desired-gains indices | `claw-selection-index` |

### Gene expression

| Need | Skill |
|---|---|
| Import Salmon/kallisto/featureCounts output; sample QC | `claw-rnaseq-import` |
| Differential expression (DESeq2) | `claw-de-analysis` |
| GO/KEGG/custom enrichment with the right background | `claw-enrichment` |
| Co-expression modules (WGCNA) | `claw-coexpression` |

### Drone and satellite imagery

| Need | Skill |
|---|---|
| Check an orthomosaic or DSM before analysis; processing guide | `claw-ortho-qc` |
| Make plot polygons from field corners and layout | `claw-plot-grid` |
| Per-plot vegetation indices (NDVI, NDRE, …) | `claw-vi-extract` |
| Canopy cover from RGB images | `claw-canopy-cover` |
| Smooth VI time series; peak, green-up, senescence, AUC | `claw-vi-timeseries` |
| Predict yield from imagery with leakage-safe validation | `claw-yield-prediction` |

### Plant physiology

| Need | Skill |
|---|---|
| Read and QC LI-COR LI-6800 files | `claw-licor-6800` |
| A-Ci curves: Vcmax, J, TPU (C3 and C4) | `claw-aci-fit` |
| Light-response curves: Asat, quantum yield, light compensation point | `claw-light-response` |
| Ball–Berry and Medlyn stomatal models, intrinsic WUE | `claw-stomatal-models` |
| Chlorophyll fluorescence: Fv/Fm, ΦPSII, NPQ, qL, ETR | `claw-fluorescence` |
| Transpiration response to VPD (breakpoints) from pot weights | `claw-vpd-transpiration` |
| Growth analysis: RGR, NAR, LAR, SLA | `claw-growth-analysis` |

## Common pipelines

**Breeding trial to candidate genes:**
1. `claw-trial-design` → `claw-fieldbook` (export).
2. Field season.
3. `claw-fieldbook` (import and QC) → `claw-trial-analysis`, run per environment, giving BLUEs.
4. `claw-met-gxe` across environments.
5. `claw-marker-qc` → `claw-gwas-crop` (or `claw-qtl-mapping` for biparental populations).
6. `claw-ld-decay` for window size → `claw-gene-lookup`.

**Selection decisions:** `claw-trial-analysis` / `claw-met-gxe` (variance components) → `claw-selection-index`; or `claw-marker-qc` → `claw-genomic-prediction`.

**Imagery to yield:**
1. `claw-ortho-qc` → `claw-plot-grid` → `claw-vi-extract` (and `claw-canopy-cover`), repeated per flight.
2. `claw-vi-timeseries` → `claw-yield-prediction`.

Predicted traits enter genetic analyses only if their out-of-sample R² is clearly positive.

**Environment context for a trial:** `claw-weather` + `claw-gdd` + `claw-soil` → `claw-water-balance`. Use `claw-yield-units` whenever yields cross unit systems.

**Expression study:** `claw-rnaseq-import` → `claw-de-analysis` → `claw-enrichment`. With about 15 or more samples, also run `claw-coexpression`.

**Leaf physiology:** `claw-licor-6800` (read and QC) → `claw-aci-fit`, `claw-light-response`, `claw-stomatal-models` and `claw-fluorescence`.

## Out of scope

- **Read alignment and variant calling:** use an established pipeline (e.g. nf-core). AgriClaw starts from a VCF or a genotype matrix.
- **Orthomosaic generation:** use OpenDroneMap or commercial software. `claw-ortho-qc` checks the result and includes a processing guide.

## References

- Agent Skills specification: https://agentskills.io/specification
- AgriClaw repository and conventions: https://github.com/bgtamang/AgriClaw
