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

An AI agent can write an analysis script in seconds. What it lacks is the field-specific judgement that decides whether the result is right. Some examples:

- A split-plot layout whose whole plots straddle two field rows.
- Corn growth-stage thresholds in °F-days compared against °C-day accumulation.
- BLUPs passed into a GWAS, which shrinks them twice.
- LI-6800 Excel files read before their formulas were ever calculated, so A and Ci come back as zeros.
- Random cross-validation that leaks genotypes between folds and inflates prediction accuracy.
- GO enrichment run against the whole genome instead of the expressed genes.

AgriClaw skills encode these decisions and test them.

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

Each skill's `compatibility` field lists exactly what that skill needs.

## Skills

Start with [`agriclaw`](skills/agriclaw/), the router skill. It maps requests to skills and shows how they chain together, for example trial design → field book → trial analysis → multi-environment analysis → GWAS → gene lookup.

<!-- skills:start -->
<!-- skills:end -->

**Status:**

- `tested`: the tests check results against known answers, and they pass.
- `experimental`: the skill runs, but some parts aren't validated yet. Its SKILL.md says which.
- `guide`: instructions only.

## Quality checks

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
