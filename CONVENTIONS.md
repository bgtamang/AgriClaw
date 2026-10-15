# AgriClaw skill conventions

Every AgriClaw skill follows these rules. They keep skills portable across AI agents, scientifically correct, and reproducible. Read this file in full before adding or changing a skill.

## 1. What a skill is

A skill is one folder, `skills/<name>/`:

```
skills/claw-example/
├── SKILL.md        instructions for the agent (required)
├── scripts/        the code the agent runs (R or Python)
├── demo/           small demo inputs, plus README.md saying where they came from
├── tests/          tests that check results against known answers
└── references/     optional longer method notes the agent reads only when needed
```

Skills get installed one folder at a time. **A skill must be self-contained**: it can't import code from another skill or from shared folders in this repo. Small helpers (argument parsing, the reproducibility bundle) are copied into each script from the templates.

Working templates: `templates/r-skill/` and `templates/python-skill/`. Copy one, rename it, and replace the example logic.

## 2. SKILL.md frontmatter (Agent Skills specification)

AgriClaw follows the open Agent Skills specification (<https://agentskills.io/specification>), so the same folder loads in Claude Code, OpenAI Codex, Gemini CLI, GitHub Copilot, Cursor, OpenClaw and other tools. Only these top-level keys are allowed:

| Key | Rule |
|---|---|
| `name` | Required. Lowercase letters, digits and single hyphens; at most 64 characters; **must equal the folder name**. |
| `description` | Required. At most 1024 characters. Says what the skill does **and** when to use it, using the words a user would actually type. Agents choose skills from this text alone. |
| `license` | `MIT`. |
| `compatibility` | At most 500 characters: language, required packages, internet needs. |
| `metadata` | A map of string keys to **string** values. AgriClaw requires `version`, `domain`, `language` and `status`. |
| `allowed-tools` | Optional and experimental in the spec. Don't use it unless needed. |

A top-level `version`, or any other key, is invalid: put extra information under `metadata`.

`metadata.status` must be one of:

- `tested`: has scripts and tests that check results against known answers, and the tests pass.
- `experimental`: runs, but validation is incomplete. SKILL.md must say exactly what is unvalidated.
- `guide`: instructions only, no scripts.

`metadata.domain` must be one of `field`, `trials`, `genomics`, `transcriptomics`, `remote-sensing`, `physiology` or `meta`.

```yaml
---
name: claw-aci-fit
description: Fit C3 and C4 CO2-response (A-Ci) curves from LI-COR LI-6800 gas-exchange data to estimate Vcmax, J, TPU and respiration with the Farquhar-von Caemmerer-Berry model. Use when the user has A-Ci or CO2-response curves and asks for Vcmax, Jmax, photosynthetic capacity or biochemical limitation.
license: MIT
compatibility: R >= 4.2 with the PhotoGEA package. No internet needed.
metadata:
  version: "1.0.0"
  domain: physiology
  language: R
  status: tested
  repository: https://github.com/bgtamang/AgriClaw
---
```

Always quote `version`. `python tools/validate_skills.py` checks all of these rules.

## 3. SKILL.md body

Write for an AI agent that knows general science but not the pitfalls of this method. Keep SKILL.md under about 300 lines and move long derivations to `references/`. Use these sections, in this order:

1. **What it does**: two or three sentences.
2. **When to use / when not to use**, naming the AgriClaw skill to use instead where relevant.
3. **Inputs**: file formats, required columns with units, and **what to ask the user** when information is missing. Never guess a planting date, maturity group, genome assembly, unit, or experimental design.
4. **How to run**: exact commands, with the demo first and then real-data examples. An installed skill runs in the user's project folder, not the skill folder, so call scripts through `${CLAUDE_SKILL_DIR}`. Start the section with this sentence: "Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md." Then write commands like `Rscript "${CLAUDE_SKILL_DIR}/scripts/<file>.R" --demo`. Paths to the user's own files stay as ordinary paths.
5. **Outputs**: the files, key columns, and how to read them.
6. **Expert decisions**: each default and threshold, why it was chosen, and its source.
7. **Checks before reporting**: diagnostics the agent must look at, and red flags that mean "stop and tell the user".
8. **References**: full citations with DOIs.

## 4. Scripts

- One main script per skill in `scripts/`, named after the skill with underscores (`scripts/aci_fit.R`). Helper scripts are allowed.
- Use R for statistics and anything with an established R package; use Python for data plumbing and API access. R scripts use base R style: base graphics, no tidyverse.
- **Wrap established, peer-reviewed packages** (agricolae, FielDHub, statgenSTA, statgenGxE, statgenGWAS, qtl2, rrBLUP, DESeq2, WGCNA, PhotoGEA, terra, …). Don't re-implement a published method when a maintained package exists. Writing small formulas directly (unit conversions, index equations) is fine, as long as a test covers them.
- Command line: `--demo` runs on the bundled demo data; `--input`, `--out` and skill-specific options use `--long-names`; `--help` prints usage. The default output folder is `agriclaw_output/<skill-name>/` in the current directory.
- Scripts must be loadable without running, so tests can call their functions. In R, end the script with `if (sys.nframe() == 0L) main()`. In Python, use `if __name__ == "__main__":`.
- Validate inputs and **fail loudly** with a message that says how to fix the problem (missing column, wrong units, too few replicates). Never drop or fill data silently. If the method needs imputation, say so in the report.
- Put units in column names: `yield_kg_ha`, `tmax_c`, `ppfd_umol_m2_s`.
- Set random seeds and record them in the output.
- Don't access the network unless the skill exists to fetch public data. Such skills write nothing outside `--out` and must have offline tests that use saved fixtures.

## 5. Outputs

Every run writes:

```
<out>/
├── report.md             what was done, key results, diagnostics, warnings, references
├── tables/*.csv
├── figures/*.png         base R graphics or matplotlib, 150 dpi
└── reproducibility/
    ├── commands.sh       exact command to re-run
    ├── sessionInfo.txt   R session info (Python skills write environment.txt)
    └── checksums.sha256  SHA-256 of inputs and outputs (R falls back to checksums.md5 without openssl)
```

Put warnings at the top of `report.md`, for example "λ = 1.4: p-values are inflated; see diagnostics". Don't bury them.

## 6. Tests

A skill needs tests to reach status `tested`.

- **R:** `tests/test_<name>.R`, plain base R, run with `Rscript tests/test_<name>.R`. Use the `check()` helper from the template, finish by printing `PASS: <n> checks`, and `stop()` on any failure so the exit code is non-zero.
- **Python:** `tests/test_<name>.py`, run with pytest.
- **Tests check results against known answers**, not just that the code runs. Use at least one of these:
  - simulate data with known true parameters and check the method recovers them within a stated tolerance;
  - reproduce a published worked example or a package's documented output;
  - check design properties, such as every entry appearing exactly once per replicate.
- Include one end-to-end test that runs the command line with `--demo` into a temporary folder and checks the expected files exist.
- Tests run offline, and each test file finishes in under about 2 minutes.

To run everything:

```
python tools/validate_skills.py
python -m pytest
Rscript tests/run_r_tests.R
```

## 7. Scientific honesty

- Every default, threshold and constant has a cited source, or is labelled in SKILL.md as an adjustable rule of thumb.
- **Every reference must be real.** Check each DOI (for example `https://api.crossref.org/works/<doi>`) before citing it. Leave a citation out rather than invent one.
- Don't claim validation that wasn't done. If the tests cover only part of the method, say which part is untested.
- State each method's limitations in SKILL.md and in `report.md`.

## 8. Demo data

Keep demo data small, ideally under 500 KB per skill. Use either synthetic data from a committed generator script with a fixed seed (`demo/make_demo.R`), or public data whose licence allows redistribution, with the source named in `demo/README.md`. Example data that ships inside an R package (PhotoGEA, qtl2, statgenSTA, …) can be loaded with `system.file()` instead of being copied.

## 9. Naming

Skill names are `claw-<short-name>`: lowercase with hyphens. Script and test file names use underscores.
