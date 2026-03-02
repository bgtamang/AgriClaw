# Contributing to AgriClaw

Thanks for your interest in contributing! AgriClaw is community-driven — every skill comes from a researcher who solved a real problem and wants others to benefit.

## What Makes a Good Skill?

A skill is **not** a thin wrapper around an existing tool. A good skill encodes **expert decisions** — the things you learned the hard way:

- The QC thresholds that actually matter
- The model choice that avoids common pitfalls
- The parameter defaults that work for real agricultural data
- The sanity checks that catch silent failures

**Ask yourself**: "What do I know that a first-year grad student gets wrong?" That's your skill.

## How to Contribute a Skill

### 1. Pick a domain folder

```
skills/
├── field/              # Trial design, GDD, soil, weather, spatial analysis
├── genomics/           # GWAS, QTL, WGS, genomic selection, marker analysis
├── transcriptomics/    # RNA-seq, differential expression, pathway enrichment
├── remote-sensing/     # Orthomosaics, vegetation indices, canopy, UAV/satellite
└── orchestrator/       # Routing and integration
```

### 2. Create your skill directory

```bash
mkdir skills/your-domain/claw-your-skill
```

### 3. Write your SKILL.md

Use [templates/SKILL-TEMPLATE.md](templates/SKILL-TEMPLATE.md) as a starting point. The `SKILL.md` file is the most important part — it tells the AI agent exactly how to run your pipeline.

### 4. Include demo data

Every skill **must** include demo data and a `--demo` flag. If someone can't try your skill in 30 seconds without their own data, they won't try it at all.

```
skills/your-domain/claw-your-skill/
├── SKILL.md
├── your_script.py
└── demo_data/
    └── example_input.csv
```

### 5. Add a reproducibility bundle

Your skill's output should include:
- `commands.sh` — exact commands to reproduce
- `environment.yml` — dependencies
- `checksums.sha256` — SHA-256 of all inputs and outputs

### 6. Submit a PR

- Fork the repo
- Create a branch: `git checkout -b add-claw-your-skill`
- Add your skill
- Test: `python your_script.py --demo` should work cleanly
- Submit a PR with a description of what expert knowledge your skill encodes

## Code Standards

- Python 3.9+ compatible
- Minimal dependencies (prefer stdlib where possible)
- No hardcoded paths
- No cloud uploads without explicit user consent
- Include type hints for public functions
- Include a brief docstring for the main script

## Naming Convention

Skills follow the pattern `claw-{name}`:
- `claw-gdd` (growing degree days)
- `claw-gwas-crop` (crop GWAS)
- `claw-trial-design` (field trial design)

Keep names short, descriptive, and lowercase with hyphens.

## Questions?

Open an [issue](https://github.com/bgtamang/AgriClaw/issues) or start a [discussion](https://github.com/bgtamang/AgriClaw/discussions).
