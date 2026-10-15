<h3 align="center">AgriClaw</h3>

<p align="center">
  <strong>Tested AI-agent skills for crop research.</strong><br>
  Field trials · genotypes and genetics · gene expression · drone and satellite imagery · plant physiology
</p>

<p align="center">
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="MIT License"></a>
  <a href="https://github.com/bgtamang/AgriClaw/actions/workflows/tests.yml"><img src="https://github.com/bgtamang/AgriClaw/actions/workflows/tests.yml/badge.svg" alt="Tests"></a>
  <a href="https://agentskills.io"><img src="https://img.shields.io/badge/format-Agent%20Skills-blue" alt="Agent Skills format"></a>
</p>

---

AgriClaw is a library of [Agent Skills](https://agentskills.io): folders of instructions and code that teach an AI agent how to run a crop-research analysis the way an experienced researcher would. They work in Claude Code, OpenAI Codex, Gemini CLI, GitHub Copilot, Cursor, OpenClaw and other tools that read the Agent Skills format.

Each skill:

- **wraps established tools** instead of re-implementing them: agricolae, FielDHub, statgenSTA, statgenGxE, statgenGWAS, qtl2, rrBLUP, DESeq2, WGCNA, PhotoGEA, terra and others;
- **encodes the expert decisions**: which design, which model, which threshold, and why, with sources;
- **asks before guessing** units, experimental design, genome assembly, planting date or band order;
- **checks itself**, with tests that compare results against known answers (simulations with known truth, published worked examples, design-property checks);
- **leaves a record**: every run writes a report (warnings first), tables, figures, and a reproducibility bundle (`commands.sh`, session info, SHA-256 checksums).

## Why

You run a field trial: 200 genotypes, 3 locations, drone flights every two weeks. Then come the growing degree days, the vegetation indices, the BLUPs across environments, the GWAS and the gene lists. Each step has pitfalls that a new student spends weeks learning and that a postdoc knows but keeps in undocumented scripts.

An AI agent can write each analysis script in seconds. What it lacks is the field-specific judgement that decides whether the result is right. Some examples:

- A split-plot layout whose whole plots straddle two field rows.
- Corn growth-stage thresholds in °F-days compared against °C-day accumulation.
- BLUPs passed into a GWAS, which shrinks them twice.
- LI-6800 Excel files read before their formulas were ever calculated, so A and Ci come back as zeros.
- Random cross-validation that leaks genotypes between folds and inflates prediction accuracy.
- GO enrichment run against the whole genome instead of the expressed genes.

AgriClaw freezes that expert knowledge into skills, and tests it.

## Install

### Claude Code

```
/plugin install agriclaw --marketplace bgtamang/AgriClaw
```

On Claude Code versions before 2.1.275, add the marketplace first:

```
/plugin marketplace add bgtamang/AgriClaw
/plugin install agriclaw@agriclaw
```

### OpenClaw

```bash
git clone https://github.com/bgtamang/AgriClaw.git
openclaw skills install ./AgriClaw/skills/claw-trial-analysis --as claw-trial-analysis   # repeat for each skill you want
```

### Other Agent Skills tools

Copy the folders you want from `skills/` into your tool's skills directory. Each folder is self-contained.

### Software the skills run

- **R ≥ 4.2.** Install every package the R skills use with `Rscript tools/install_r_deps.R`.
- **Python ≥ 3.9.** Run `python -m pip install -r requirements.txt`.

On Windows the R installer doesn't put `Rscript` on the PATH. Add R's `bin` folder (for example `C:\Program Files\R\R-4.6.1\bin`) to your PATH so agents can run the R skills.

Each skill's `compatibility` field lists exactly what that skill needs.

## Try it

Once installed, ask in plain language and the agent picks the skills. For example:

| Ask | Skills used |
|---|---|
| "Design an alpha-lattice for 120 soybean entries in 3 reps, 12 plots per field row." | `claw-trial-design` |
| "My corn was planted April 25 near Champaign. How many GDD so far, and when will it silk?" | `claw-gdd` |
| "Fit these LI-6800 A-Ci curves and give me Vcmax and Jmax per plant." | `claw-licor-6800` → `claw-aci-fit` |
| "Run a GWAS on these BLUEs with this VCF and list the genes under the peaks." | `claw-marker-qc` → `claw-gwas-crop` → `claw-ld-decay` → `claw-gene-lookup` |
| "Make plot polygons from these four RTK corners and extract NDVI per plot." | `claw-plot-grid` → `claw-vi-extract` |

The agent asks for anything it can't know (units, design, genome assembly) instead of guessing. Each run writes `report.md` with the warnings at the top.

Every skill also runs on its own with bundled demo data:

```bash
python skills/claw-gdd/scripts/gdd.py --demo
Rscript skills/claw-trial-design/scripts/trial_design.R --demo
```

## Skills

Start with [`agriclaw`](skills/agriclaw/), the router skill. It maps requests to skills and shows how they chain together, for example trial design → field book → trial analysis → multi-environment analysis → GWAS → gene lookup.

<!-- skills:start -->
### Start here

| Skill | What it does | Language | Status |
|---|---|---|---|
| [`agriclaw`](skills/agriclaw/) | Which AgriClaw skill to use for a request, and how the skills chain together | none | guide |

### Field and environment

| Skill | What it does | Language | Status |
|---|---|---|---|
| [`claw-gdd`](skills/claw-gdd/) | Growing degree days, corn growth stages (NCH-40) and freeze dates from planting | Python | tested |
| [`claw-soil`](skills/claw-soil/) | USDA SSURGO soil at a field: texture, drainage, pH, organic matter, water storage | Python | tested |
| [`claw-water-balance`](skills/claw-water-balance/) | FAO-56 daily soil-water balance: crop water use, water stress, irrigation timing | Python | tested |
| [`claw-weather`](skills/claw-weather/) | Growing-season weather: temperature, rain, reference ET, heat and dry spells | Python | tested |
| [`claw-yield-units`](skills/claw-yield-units/) | Yield conversions (bu/ac, kg/ha, t/ha), moisture adjustment, plot-to-area scaling | Python | tested |

### Field trials

| Skill | What it does | Language | Status |
|---|---|---|---|
| [`claw-fieldbook`](skills/claw-fieldbook/) | Field Book app export and import, merging tablets, plot-level data QC | R | tested |
| [`claw-met-gxe`](skills/claw-met-gxe/) | Multi-environment trials: GxE variance, heritability, stability, Finlay-Wilkinson, AMMI/GGE | R | tested |
| [`claw-trial-analysis`](skills/claw-trial-analysis/) | Single-trial mixed models (SpATS spatial or lme4): BLUEs, BLUPs, heritability, outliers | R | tested |
| [`claw-trial-design`](skills/claw-trial-design/) | Randomise and lay out trials: RCBD, alpha-lattice, row-column, augmented, p-rep, split, strip | R | tested |

### Genotypes and genetics

| Skill | What it does | Language | Status |
|---|---|---|---|
| [`claw-gene-lookup`](skills/claw-gene-lookup/) | Genes inside GWAS or QTL intervals from Ensembl Plants or a local GFF3, with assembly checks | R | tested |
| [`claw-genomic-prediction`](skills/claw-genomic-prediction/) | Genomic prediction (GBLUP, RR-BLUP, Bayesian) with leakage-safe cross-validation | R | tested |
| [`claw-gwas-crop`](skills/claw-gwas-crop/) | Mixed-model GWAS with kinship, LOCO and principal components (statgenGWAS) | R | tested |
| [`claw-ld-decay`](skills/claw-ld-decay/) | LD decay curves and the critical distance for candidate-gene windows | R | tested |
| [`claw-marker-qc`](skills/claw-marker-qc/) | Genotype QC: call rate, MAF, heterozygosity, duplicates, distortion; VCF and HapMap | R | tested |
| [`claw-qtl-mapping`](skills/claw-qtl-mapping/) | QTL mapping in biparental and multiparent populations with permutation thresholds (qtl2) | R | tested |
| [`claw-selection-index`](skills/claw-selection-index/) | Response to selection, selection intensity, Smith-Hazel and desired-gains indices | R | tested |

### Gene expression

| Skill | What it does | Language | Status |
|---|---|---|---|
| [`claw-coexpression`](skills/claw-coexpression/) | WGCNA co-expression modules and module-trait relationships | R | tested |
| [`claw-de-analysis`](skills/claw-de-analysis/) | Differential expression with DESeq2: contrasts, likelihood-ratio tests, shrunken fold changes | R | tested |
| [`claw-enrichment`](skills/claw-enrichment/) | GO, KEGG, MapMan or custom term enrichment against the correct background | R | tested |
| [`claw-rnaseq-import`](skills/claw-rnaseq-import/) | Import Salmon, kallisto or featureCounts output to gene counts, with sample QC | R | tested |

### Drone and satellite imagery

| Skill | What it does | Language | Status |
|---|---|---|---|
| [`claw-canopy-cover`](skills/claw-canopy-cover/) | Per-plot green canopy cover from RGB drone orthomosaics | R | tested |
| [`claw-ortho-qc`](skills/claw-ortho-qc/) | Check a drone orthomosaic and DSM before plot extraction; processing guide | R | tested |
| [`claw-plot-grid`](skills/claw-plot-grid/) | Plot polygons from field corners and layout, joined to the field book | R | tested |
| [`claw-vi-extract`](skills/claw-vi-extract/) | Per-plot vegetation indices (NDVI, NDRE, GNDVI, OSAVI, EVI2, MCARI2) from multispectral imagery | R | tested |
| [`claw-vi-timeseries`](skills/claw-vi-timeseries/) | Smooth vegetation-index time series and extract phenology: green-up, peak, senescence, AUC | R | tested |
| [`claw-yield-prediction`](skills/claw-yield-prediction/) | Predict yield from imagery features with grouped (leakage-safe) cross-validation | R | tested |

### Plant physiology

| Skill | What it does | Language | Status |
|---|---|---|---|
| [`claw-aci-fit`](skills/claw-aci-fit/) | Fit C3 and C4 A-Ci curves: Vcmax, J, TPU, Vpmax (PhotoGEA) | R | tested |
| [`claw-fluorescence`](skills/claw-fluorescence/) | PAM chlorophyll fluorescence: Fv/Fm, PhiPSII, qP, qL, NPQ, ETR | R | tested |
| [`claw-growth-analysis`](skills/claw-growth-analysis/) | Growth analysis from sequential harvests: RGR, NAR, LAR, SLA, LMF | R | tested |
| [`claw-licor-6800`](skills/claw-licor-6800/) | Read and QC LI-COR LI-6800 log files (Excel and plain text) | R | tested |
| [`claw-light-response`](skills/claw-light-response/) | Light-response curves: Asat, quantum yield, compensation and saturation points | R | tested |
| [`claw-stomatal-models`](skills/claw-stomatal-models/) | Ball-Berry and Medlyn stomatal model fits and intrinsic water-use efficiency | R | tested |
| [`claw-vpd-transpiration`](skills/claw-vpd-transpiration/) | Transpiration response to VPD from pot weights, with breakpoint models | R | tested |

34 skills: 33 tested, 0 experimental, 1 guide.
<!-- skills:end -->

**Status:**

- `tested`: the tests check results against known answers, and they pass.
- `experimental`: the skill runs, but some parts aren't validated yet. Its SKILL.md says which.
- `guide`: instructions only.

## Quality checks

The 33 skills carry about 1,800 checks against known answers: 150 Python tests and 1,644 R checks. The answers come from published worked examples (FAO-56, NCH-40, the qtl2 user guide, PhotoGEA's article tables, Cerón-Rojas & Crossa 2018), simulations with known truth, and design-property checks. Every cited DOI was looked up in Crossref or DataCite (`tools/check_dois.py`).

Run these from the repository root:

```bash
python -m pip install -r requirements-dev.txt
python tools/validate_skills.py        # every SKILL.md against the Agent Skills spec and AgriClaw rules
python -m pytest                       # Python skill tests
Rscript tests/run_r_tests.R            # R skill tests, each in a fresh R session
```

GitHub Actions runs all three on every push.

## Contributing

New skills and fixes are welcome. Read [CONVENTIONS.md](CONVENTIONS.md) (the rules every skill follows) and [CONTRIBUTING.md](CONTRIBUTING.md), and start from `templates/r-skill` or `templates/python-skill`.

## Citation

If AgriClaw helps your research, please cite the software. [CITATION.cff](CITATION.cff) has the details, and GitHub shows a "Cite this repository" button. Please also cite the packages and methods each skill names in its report.

## License

MIT. See [LICENSE](LICENSE).
