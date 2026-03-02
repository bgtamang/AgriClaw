<h3 align="center">AgriClaw</h3>

<p align="center">
  <strong>Agriculture-native AI agent skill library.</strong><br>
  Built on <a href="https://github.com/openclaw/openclaw">OpenClaw</a>. Local-first. Reproducible.
</p>

<p align="center">
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="MIT License"></a>
  <a href="#skills"><img src="https://img.shields.io/badge/skills-8_planned-orange" alt="Skills"></a>
  <a href="https://github.com/bgtamang/AgriClaw/issues"><img src="https://img.shields.io/github/issues/bgtamang/AgriClaw" alt="Open Issues"></a>
</p>

---

## The Problem

You run a field trial. 200 genotypes, 3 locations, UAV flights every two weeks. Now you need to:

1. Compute growing degree days. Open Excel, find the weather station, download CSV, clean it, calculate GDD manually.
2. Extract vegetation indices from drone imagery. Which band combination? NDVI or NDRE? What about soil background?
3. Fit BLUPs across environments. Was it `lme4` or `sommer`? Do you need spatial correction? Which variance structure?
4. Run GWAS. GAPIT? TASSEL? Mixed model or GLM? What's the right significance threshold?
5. Annotate hits. Which genome assembly? Liftover needed? Where's the gene list?

Each step is a rabbit hole. A new student spends **weeks** learning the pitfalls. A postdoc has the knowledge but it lives in undocumented scripts on their laptop.

**AgriClaw freezes that expert knowledge into skills that an AI agent executes correctly every time.**

---

## What Is AgriClaw?

A **skill** is a domain expert's pipeline — frozen into code — that an AI agent runs reproducibly.

```
LLM alone        = guesses at agriculture analysis, hallucinates thresholds
AgriClaw skill    = a breeder's proven pipeline that the AI executes correctly
```

- **Local-first**: Your field data, genotypes, and imagery stay on your machine.
- **Reproducible**: Every analysis exports `commands.sh`, `environment.yml`, and SHA-256 checksums.
- **Modular**: Each skill is a self-contained directory (`SKILL.md` + scripts) organised by domain.
- **Expert-driven**: Skills encode the hard-won decisions — the QC thresholds, the model choices, the pitfalls to avoid.

### Why not just ask an LLM?

Ask an LLM to "run GWAS on my soybean data." It will:

- Use GLM instead of MLM, ignoring population structure
- Skip kinship correction entirely
- Set an arbitrary p-value threshold instead of FDR
- Produce results with genomic inflation factor > 10
- No reproducibility bundle, no audit trail

AgriClaw encodes the correct decisions so the agent gets it right first time.

---

## Skills

Skills are organised by domain:

### Phenomics

| Skill | Status | Description |
|-------|--------|-------------|
| [claw-gdd](skills/phenomics/claw-gdd/) | **MVP** | Growing degree days, frost dates, growth stage estimation from free weather APIs |
| [claw-vi-extract](skills/phenomics/claw-vi-extract/) | Planned | Vegetation index extraction from UAV/satellite multispectral imagery |

### Genomics

| Skill | Status | Description |
|-------|--------|-------------|
| [claw-gwas-crop](skills/genomics/claw-gwas-crop/) | **MVP** | Mixed-model GWAS with kinship, FDR control, and genomic inflation checks |
| [claw-blup](skills/genomics/claw-blup/) | Planned | Multi-environment BLUP computation with R² pre-filtering |
| [claw-qtl-annotator](skills/genomics/claw-qtl-annotator/) | Planned | SNP annotation with genome liftover and gene lookup (Phytozome/Ensembl Plants) |

### Field

| Skill | Status | Description |
|-------|--------|-------------|
| [claw-trial-design](skills/field/claw-trial-design/) | **MVP** | Field trial randomisation (RCBD, alpha-lattice, augmented designs) |
| [claw-soil](skills/field/claw-soil/) | Planned | Soil properties from USDA SSURGO — texture, OM, drainage, pH |
| [claw-weather](skills/field/claw-weather/) | Planned | Weather summaries, drought indices, precipitation from PRISM/Open-Meteo |

### Orchestrator

| Skill | Status | Description |
|-------|--------|-------------|
| [agri-orchestrator](skills/orchestrator/agri-orchestrator/) | Planned | Routes agriculture requests to the right specialist skill |

---

## Quick Start

### Prerequisites

- [OpenClaw](https://github.com/openclaw/openclaw) installed
- Python 3.9+

### Install and run

```bash
# Clone
git clone https://github.com/bgtamang/AgriClaw.git
cd AgriClaw

# Run a skill with demo data
python skills/phenomics/claw-gdd/gdd.py --demo

# Or use via OpenClaw
openclaw install skills/phenomics/claw-gdd
openclaw "Calculate GDD for soybeans planted May 15 in Champaign IL"
```

Every skill includes **demo data** — try it immediately without your own files.

---

## Architecture

```
User: "How many GDD have accumulated for my soybean trial?"
         |
  +------v--------+
  |  Agri          |  <- routes by data type + keywords
  |  Orchestrator  |
  +------+---------+
         |
  +------v------------------------------------------+
  |                                                  |
  GDD       GWAS      Trial      Soil      VI
  Calc      Crop      Design     Lookup    Extract  ...
  |                                                  |
  +------+------------------------------------------+
         |
  +------v--------+
  |  Markdown      |  <- report + figures + checksums
  |  Report        |     + reproducibility bundle
  +---------------+
```

Each skill works standalone — the orchestrator routes, but you can run any skill directly.

---

## Reproducibility

Every AgriClaw analysis produces a **reproducibility bundle**:

```
output/
  report.md              # Full analysis with figures and tables
  figures/               # Publication-quality plots
  tables/                # CSV data tables
  commands.sh            # Exact commands to reproduce
  environment.yml        # Conda/pip environment snapshot
  checksums.sha256       # SHA-256 of every input and output
```

A collaborator can reproduce your results without emailing you. A reviewer can verify your analysis in 30 seconds.

---

## Community Wanted Skills

We want skills from the agricultural research community. If you work with crops, livestock, soil, or remote sensing — **wrap your pipeline as a skill**.

| Skill | Domain | Your expertise |
|-------|--------|----------------|
| **claw-gxe** | Genomics | G x E interaction (AMMI, GGE biplot) |
| **claw-gs** | Genomics | Genomic selection / prediction |
| **claw-ld** | Genomics | LD analysis, haplotype blocks, pruning |
| **claw-canopy** | Phenomics | Canopy cover estimation from RGB imagery |
| **claw-yield-predict** | Phenomics | ML yield prediction from time-series imagery |
| **claw-spatial** | Field | Spatial field correction (SpATS, AR1xAR1) |
| **claw-irrigation** | Field | Irrigation scheduling from ET + soil data |

See [CONTRIBUTING.md](CONTRIBUTING.md) for the submission process and [templates/SKILL-TEMPLATE.md](templates/SKILL-TEMPLATE.md) for the skill template.

---

## Contributing

We welcome contributions! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

The short version:
1. Fork the repo
2. Create your skill using [templates/SKILL-TEMPLATE.md](templates/SKILL-TEMPLATE.md)
3. Include demo data and a `--demo` flag
4. Submit a PR

---

## Citation

If you use AgriClaw in your research, please cite:

```bibtex
@software{agriclaw_2026,
  author = {Tamang, Bishal},
  title = {AgriClaw: An Agriculture-Native AI Agent Skill Library for Reproducible Crop Research},
  year = {2026},
  url = {https://github.com/bgtamang/AgriClaw}
}
```

## Links

- [OpenClaw](https://github.com/openclaw/openclaw) — The agent platform
- [ClawBio](https://github.com/ClawBio/ClawBio) — Bioinformatics skill library (inspiration)

## License

MIT — clone it, run it, build a skill, submit a PR.
