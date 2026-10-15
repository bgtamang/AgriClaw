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
| claw-trial-design | trials | R | done | agricolae and FielDHub designs laid out as compact field blocks; design-property tests (alpha concurrence and efficiency factors, checks in every augmented block, intact whole plots and strips, no dropped entries); 75 checks |
| claw-fieldbook | trials | R | done | Field Book app export/import, wide files and two-tablet merges; 15 planted error types caught; 69 checks |
| claw-trial-analysis | trials | R | done | statgenSTA (SpATS or lme4) BLUEs/BLUPs and Cullis H2; spatial model more accurate than lme4 on trended fields; 44 checks |
| claw-met-gxe | trials | R | done | statgenGxE variance components, Finlay-Wilkinson, AMMI/GGE, SE-weighted stage 2, stability (Shukla checked against agricolae); 58 checks |
| claw-marker-qc | genomics | R | done | call rate, MAF, heterozygosity, duplicates, segregation distortion; planted problems found exactly; VCF/HapMap I/O; 64 checks |
| claw-gwas-crop | genomics | R | done | statgenGWAS mixed model with kinship, LOCO and PCs; 35/36 simulated QTNs found, null lambda 0.94-1.08; BH FDR (statgenGWAS fdr option is faulty); 50 checks |
| claw-ld-decay | genomics | R | done | r2 decay, Hill-Weir fit and critical distance; refuses spurious windows in low-LD panels; 50 checks |
| claw-gene-lookup | genomics | R | done | genes in intervals from Ensembl Plants REST or a local GFF3; assembly and gene-ID checks; 50 offline checks + 1 live |
| claw-qtl-mapping | genomics | R | done | qtl2 scans with permutation thresholds and LOD/Bayes intervals; reproduces the qtl2 user guide (iron data); true QTL inside the 1.5-LOD interval in 20/20 simulations; 59 checks |
| claw-genomic-prediction | genomics | R | done | GBLUP/RR-BLUP and BGLR Bayesian models with family- or environment-grouped CV; accuracy tracks Daetwyler's expectations; random vs leave-family-out inflation shown (0.71 vs 0.40); 53 checks |
| claw-selection-index | genomics | R | done | selection intensities, Smith-Hazel, restricted and desired-gains indices; reproduces the Ceron-Rojas & Crossa (2018) worked examples; 65 checks |
| claw-rnaseq-import | transcriptomics | R | done | Salmon, kallisto or featureCounts import; tx2gene from GFF3/GTF; transcript-version and mapping-rate QC; 54 checks |
| claw-de-analysis | transcriptomics | R | done | DESeq2 Wald and LRT with apeglm shrinkage; empirical FDR 0.06 at nominal 0.05; refuses TPM, unreplicated or confounded designs; LRT anti-conservatism measured and reported; 56 checks |
| claw-enrichment | transcriptomics | R | done | over-representation with the tested genes as the universe; GO propagation via is_a/part_of; planted terms found, 0 of 47 null terms called; 44 checks |
| claw-coexpression | transcriptomics | R | done | WGCNA modules and module-trait links; planted modules recovered (ARI 0.95); 47 checks |
| claw-plot-grid | remote-sensing | R | done | plot polygons from RTK corners or origin and bearing, in UTM with the scale-factor correction; serpentine or row-wise numbering; field-book join with a position cross-check (works on claw-trial-design output); 80 checks |
| claw-vi-extract | remote-sensing | R | done | 7 indices checked against the primary formulas; coverage-weighted plot means; soil masking; reflectance-scale and band-order checks; 70 checks |
| claw-canopy-cover | remote-sensing | R | done | ExG/ExGR with one Otsu threshold; within 0.83 percentage points on synthetic plots of known cover; glare and shadow fractions; 45 checks |
| claw-ortho-qc | remote-sensing | R | done | 18 orthomosaic/DSM checks (CRS, rotation, DN vs reflectance, NoData, saturation, DSM spikes); every injected problem detected; processing guide; 52 checks |
| claw-vi-timeseries | remote-sensing | R | done | double logistic (Beck/Elmore), LOESS or Savitzky-Golay on days-after-planting or GDD; upper-envelope weights for cloudy dates; peak/green-up/senescence recovered within 4 d on simulated seasons; 76 checks |
| claw-yield-prediction | remote-sensing | R | done | elastic net and random forest with nested grouped CV (leave-environment/genotype/year-out); leakage shown (r 0.90 random vs 0.20 environment-out); 75 checks |
| claw-licor-6800 | physiology | R | done | reads LI-6800 Excel and plaintext logs (values match PhotoGEA's reader); stops on unevaluated Excel formulas; duplicate-observation and instrument-confounding checks; 55 checks |
| claw-aci-fit | physiology | R | done | PhotoGEA C3/C4 A-Ci fits; reproduces PhotoGEA's published Vcmax, J, TPU and Vpmax tables within 0.5%; simulated-curve recovery; 62 checks in 3 files |
| claw-light-response | physiology | R | done | non-rectangular hyperbola with compensation and saturation points; simulated recovery within 3.5 SE; 34 checks |
| claw-stomatal-models | physiology | R | done | Ball-Berry and Medlyn via PhotoGEA; reproduces PhotoGEA's Ball-Berry table; catches unit mistakes; 38 checks |
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
| GitHub Actions CI (validator + Python + R tests) | done (`.github/workflows/tests.yml`) |
| README rewrite (honest status table, install for Claude Code / OpenClaw / others) | done; table filled by `tools/build_catalog.py` from each skill's `metadata.summary` |
| CONTRIBUTING.md update, CHANGELOG.md, CITATION.cff | done |
| `tools/check_dois.py` (looks up every cited DOI in Crossref/DataCite) | done; needs internet, so not in CI |
| Plugin install check | done: installs from a local marketplace in a clean Claude Code config, all 34 skills listed; a `--plugin-dir` session ran `claw-gdd` through `${CLAUDE_SKILL_DIR}` |

## How the build is organised

Nine build agents worked in parallel on 2026-10-15, one per area: field (Python), trials, markers/GWAS, breeding genetics, RNA-seq, imagery, imagery time series, gas exchange, and plant water/growth/fluorescence. Each could edit only its own skill folders. The orchestrator reviewed each area (validator, tests, DOI check, link check, spot checks of the science) and committed area by area.

## Ideas for later versions

Each skill's SKILL.md lists its own limitations. The larger gaps:

- No real-data benchmarks yet for imagery (canopy cover, VI time series) or yield prediction; validation is on simulated or synthetic data.
- `claw-genomic-prediction`: no GxE kernel in the multi-environment model.
- `claw-qtl-mapping`: single-QTL scans only.
- `claw-water-balance`: constant root depth; no runoff or capillary rise unless supplied.
- `claw-enrichment`: over-representation only (no GSEA).
- Not covered: crop-model calibration (DSSAT/APSIM), hyperspectral imagery, root phenotyping, seed composition (NIR).
- Upstream package problems found during testing are listed in CHANGELOG.md; consider reporting them to the package maintainers.

## Decisions log

- 2026-10-15: Follow the Agent Skills spec strictly. Only name, description, license, compatibility, metadata and allowed-tools appear at the top level; `metadata` values are strings. The OpenClaw-specific `metadata.openclaw` block was dropped, because the spec requires string values. OpenClaw still loads the skills.
- 2026-10-15: R for statistics, wrapping established packages (agricolae, FielDHub, statgen*, qtl2, rrBLUP, DESeq2, WGCNA, PhotoGEA, terra). Python only for the field/weather data plumbing.
- 2026-10-15: Removed the earlier fabricated or garbled citations (Federer 1961 HPR, "Cullis 2006 J Agric Sci", "Gao 1992"). Every reference is now checked against Crossref.
- 2026-10-15: Skills are flat folders (`skills/<name>/`) so they install directly in Claude Code, OpenClaw and other tools. Domains are recorded in `metadata.domain`.
- 2026-10-15: SKILL.md commands call scripts through `${CLAUDE_SKILL_DIR}`, because an installed skill runs from the user's project folder. Claude Code substitutes the variable; other tools are told to use the skill folder.
- Out of scope: read alignment and variant calling. Use established pipelines such as nf-core; AgriClaw starts from a VCF or genotype matrix.
