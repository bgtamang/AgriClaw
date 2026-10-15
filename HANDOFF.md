# AgriClaw v1 rebuild: handoff

Working branch: `v1-rebuild`. `main` stays unchanged until the owner reviews the rebuild.
Started 2026-10-15.

## Goal

Turn the March 2026 scaffold into a tested, portable skill library for crop research:

- rebuild in the open Agent Skills format, one folder per skill;
- fix the original skills (the review found wrong texture classes, °F/°C mix-ups in corn growth stages, and invalid trial designs);
- add skills across field trials, genomics, transcriptomics, remote sensing and plant physiology.

**Every skill must pass tests that check its results against known answers.**

## How to resume (any person or agent session)

1. `git clone https://github.com/bgtamang/AgriClaw && cd AgriClaw && git checkout v1-rebuild`
2. Read `CONVENTIONS.md` (the rules), then the status table below.
3. Check the current state:
   ```
   python -m pip install -r requirements-dev.txt
   Rscript tools/install_r_deps.R
   python tools/validate_skills.py
   python -m pytest
   Rscript tests/run_r_tests.R
   ```
4. Continue with any skill that isn't `done`. Update this table in the same commit as the skill.

## Status

Status values: `planned`, `building`, `done` (validator clean and tests pass), `experimental` (runs, but validation is incomplete).

| Skill | Domain | Lang | Status | Notes |
|---|---|---|---|---|
| claw-gdd | field | Python | building | fix °F/°C stage thresholds, in-season dates, CSV planting filter, frost logic, citations |
| claw-weather | field | Python | building | fix missing-data handling, in-season dates, data-source caveats |
| claw-soil | field | Python | building | use SSURGO texture names; fix fallback classifier; add root-zone water storage |
| claw-water-balance | field | Python | building | FAO-56 daily soil-water balance with crop coefficients |
| claw-yield-units | field | Python | building | bu/ac ↔ kg/ha, moisture adjustment, plot-to-area scaling |
| claw-trial-design | trials | R | building | rebuild on agricolae / FielDHub; designs checked for validity |
| claw-fieldbook | trials | R | building | Field Book app export/import and data QC |
| claw-trial-analysis | trials | R | building | single-trial spatial analysis, BLUEs/BLUPs, heritability (statgenSTA/SpATS) |
| claw-met-gxe | trials | R | building | multi-environment analysis: GxE, stability, AMMI/GGE (statgenGxE) |
| claw-marker-qc | genomics | R | building | genotype QC: call rate, MAF, heterozygosity, duplicates |
| claw-gwas-crop | genomics | R | building | mixed-model GWAS with kinship (statgenGWAS) |
| claw-ld-decay | genomics | R | building | LD decay and critical distance for candidate windows |
| claw-gene-lookup | genomics | R | building | genes in QTL/GWAS intervals (Ensembl Plants REST or local GFF3), assembly checks |
| claw-qtl-mapping | genomics | R | building | biparental QTL mapping with permutation thresholds (qtl2) |
| claw-genomic-prediction | genomics | R | building | GBLUP / Bayesian prediction with leakage-safe cross-validation |
| claw-selection-index | genomics | R | building | breeder's equation, Smith-Hazel and desired-gains indices |
| claw-rnaseq-import | transcriptomics | R | building | Salmon/kallisto import (tximport) and sample QC |
| claw-de-analysis | transcriptomics | R | building | DESeq2 differential expression with contrasts, shrinkage, LRT |
| claw-enrichment | transcriptomics | R | building | GO/KEGG/custom enrichment with the correct background (clusterProfiler) |
| claw-coexpression | transcriptomics | R | building | WGCNA co-expression modules and module-trait links |
| claw-plot-grid | remote-sensing | R | building | plot polygons from field corners and layout |
| claw-vi-extract | remote-sensing | R | building | per-plot vegetation indices from multispectral imagery |
| claw-canopy-cover | remote-sensing | R | building | canopy cover from RGB imagery |
| claw-ortho-qc | remote-sensing | R | building | orthomosaic quality checks and processing guide |
| claw-vi-timeseries | remote-sensing | R | building | VI time-series smoothing and phenology metrics |
| claw-yield-prediction | remote-sensing | R | building | yield prediction with leakage-safe (grouped) cross-validation |
| claw-licor-6800 | physiology | R | building | read and QC LI-COR LI-6800 files |
| claw-aci-fit | physiology | R | building | C3/C4 A-Ci fitting (PhotoGEA) |
| claw-light-response | physiology | R | building | light-response curves (non-rectangular hyperbola) |
| claw-stomatal-models | physiology | R | building | Ball-Berry and Medlyn stomatal model fits |
| claw-vpd-transpiration | physiology | R | building | transpiration response to VPD with breakpoint models |
| claw-growth-analysis | physiology | R | building | RGR, NAR, LAR, SLA from sequential harvests |
| claw-fluorescence | physiology | R | building | PAM chlorophyll-fluorescence parameters |
| agriclaw | meta | none | done | catalog and pipeline router across all skills |

## Repo-level work

| Item | Status |
|---|---|
| Flat `skills/<name>/` layout, CONVENTIONS.md, templates, validator, R/Python test runners | done |
| `tools/install_r_deps.R` (all R packages) | done |
| Claude Code plugin + marketplace manifest (`.claude-plugin/`) | done (passes `claude plugin validate`) |
| GitHub Actions CI (validator + Python + R tests) | written (`.github/workflows/tests.yml`), not pushed until skills pass locally |
| README rewrite (honest status table, install for Claude Code / OpenClaw / others) | done; table filled by `tools/build_catalog.py` |
| CONTRIBUTING.md update, CHANGELOG.md, CITATION.cff | done |

## How the build is organised

Nine build agents worked in parallel on 2026-10-15, one per area: field (Python), trials, markers/GWAS, breeding genetics, RNA-seq, imagery, imagery time series, gas exchange, and plant water/growth/fluorescence. Each could edit only its own skill folders. The orchestrator reviews each area, runs the validator and tests, and commits area by area.

## Decisions log

- 2026-10-15: Follow the Agent Skills spec strictly. Only name, description, license, compatibility, metadata and allowed-tools appear at the top level; `metadata` values are strings. The OpenClaw-specific `metadata.openclaw` block was dropped, because the spec requires string values. OpenClaw still loads the skills.
- 2026-10-15: R for statistics, wrapping established packages (agricolae, FielDHub, statgen*, qtl2, rrBLUP, DESeq2, WGCNA, PhotoGEA, terra). Python only for the field/weather data plumbing.
- 2026-10-15: Removed the earlier fabricated or garbled citations (Federer 1961 HPR, "Cullis 2006 J Agric Sci", "Gao 1992"). Every reference is now checked against Crossref.
- 2026-10-15: Skills are flat folders (`skills/<name>/`) so they install directly in Claude Code, OpenClaw and other tools. Domains are recorded in `metadata.domain`.
- 2026-10-15: SKILL.md commands call scripts through `${CLAUDE_SKILL_DIR}`, because an installed skill runs from the user's project folder. Claude Code substitutes the variable; other tools are told to use the skill folder.
- Out of scope: read alignment and variant calling. Use established pipelines such as nf-core; AgriClaw starts from a VCF or genotype matrix.
