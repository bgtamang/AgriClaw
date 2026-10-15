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
| claw-gdd | field | Python | done | NCH-40 stages in °F-days (R1 Jul 9 for the Apr 25 demo); in-season end-date clamp; CSV planting filter; NCEI freeze logic; no soybean stages (photoperiod); 27 tests |
| claw-weather | field | Python | done | missing values counted, never zeroed; in-season clamp with forecast past_days; station CSV input; reanalysis caveats; 12 tests |
| claw-soil | field | Python | done | SSURGO texture class first, corrected SSM 2017 fallback (all 12 classes, 21 boundaries); AWS rebuilt from horizons matches muaggatt; 73 tests |
| claw-water-balance | field | Python | done | FAO-56 single-Kc daily balance; reproduces Examples 36-38 and Box 15; mass balance closes; 20 tests |
| claw-yield-units | field | Python | done | bu/ac and kg/ha with AH-697 bushel weights; moisture shrink checked against UF AG442 and FSA-1078 examples; 18 tests |
| claw-trial-design | trials | R | building | rebuild on agricolae / FielDHub; designs checked for validity |
| claw-fieldbook | trials | R | building | Field Book app export/import and data QC |
| claw-trial-analysis | trials | R | building | single-trial spatial analysis, BLUEs/BLUPs, heritability (statgenSTA/SpATS) |
| claw-met-gxe | trials | R | building | multi-environment analysis: GxE, stability, AMMI/GGE (statgenGxE) |
| claw-marker-qc | genomics | R | done | call rate, MAF, heterozygosity, duplicates, segregation distortion; planted problems found exactly; VCF/HapMap I/O; 64 checks |
| claw-gwas-crop | genomics | R | done | statgenGWAS mixed model with kinship, LOCO and PCs; 35/36 simulated QTNs found, null lambda 0.94-1.08; BH FDR (statgenGWAS fdr option is faulty); 50 checks |
| claw-ld-decay | genomics | R | done | r2 decay, Hill-Weir fit and critical distance; refuses spurious windows in low-LD panels; 50 checks |
| claw-gene-lookup | genomics | R | done | genes in intervals from Ensembl Plants REST or a local GFF3; assembly and gene-ID checks; 50 offline checks + 1 live |
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
| claw-vi-timeseries | remote-sensing | R | done | double logistic (Beck/Elmore), LOESS or Savitzky-Golay on days-after-planting or GDD; upper-envelope weights for cloudy dates; peak/green-up/senescence recovered within 4 d on simulated seasons; 76 checks |
| claw-yield-prediction | remote-sensing | R | done | elastic net and random forest with nested grouped CV (leave-environment/genotype/year-out); leakage shown (r 0.90 random vs 0.20 environment-out); 75 checks |
| claw-licor-6800 | physiology | R | building | read and QC LI-COR LI-6800 files |
| claw-aci-fit | physiology | R | building | C3/C4 A-Ci fitting (PhotoGEA) |
| claw-light-response | physiology | R | building | light-response curves (non-rectangular hyperbola) |
| claw-stomatal-models | physiology | R | building | Ball-Berry and Medlyn stomatal model fits |
| claw-vpd-transpiration | physiology | R | done | VPD (FAO-56 examples), irrigation/spike QC, segmented breakpoint with AICc and Davies test; breakpoint within 0.2 kPa in 86% of simulations; 55 checks |
| claw-growth-analysis | physiology | R | done | RGR, NAR, SLA, LAR, LMF (classical and spline); Hoffmann-Poorter bias reproduced; 40 checks |
| claw-fluorescence | physiology | R | done | Fv/Fm, PhiPSII, qP, qL, NPQ, ETR; match the LI-6800's own values on 424 rows; 124 checks |
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
