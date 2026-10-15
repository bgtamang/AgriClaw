---
name: claw-example
description: Summarise a numeric trait by group (n, mean, SD, CV) from a CSV and flag groups with too few observations. Template skill that shows the AgriClaw layout; use as a starting point for new skills, not for real analyses.
license: MIT
compatibility: Python >= 3.9 with matplotlib. No internet needed.
metadata:
  version: "1.0.0"
  domain: meta
  language: Python
  status: tested
  repository: https://github.com/bgtamang/AgriClaw
---

# Example: trait summary by group

## What it does

Reads a CSV with one row per observation and reports, for each group, the number of observations, mean, standard deviation and coefficient of variation, flagging groups with fewer than `--min-n` observations.

## When to use / when not to use

- Use as the template for a new AgriClaw Python skill.
- Don't use it for real trial analysis: use `claw-trial-analysis`, which accounts for the experimental design.

## Inputs

A CSV file with at least:

| Column | Meaning |
|---|---|
| `genotype` (or `--group`) | grouping variable |
| `yield_kg_ha` (or `--trait`) | numeric trait, units in the column name |

If the user's trait column has no units in its name, ask what the units are before reporting results.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
python "${CLAUDE_SKILL_DIR}/scripts/example.py" --demo
python "${CLAUDE_SKILL_DIR}/scripts/example.py" --input trial.csv --trait yield_kg_ha --group genotype --out results/
```

## Outputs

- `report.md`: warnings first, then the summary table.
- `tables/summary.csv`: one row per group with `n`, `mean`, `sd`, `cv_pct`, `low_n`.
- `figures/means.png`: dot chart of group means.
- `reproducibility/`: `commands.sh`, `environment.txt` and checksums.

## Expert decisions

- Missing values are excluded and **counted in the report**, never dropped silently.
- `sd` uses the n − 1 denominator (`statistics.stdev`).
- Groups with n < 3 (adjustable rule of thumb, `--min-n`) are flagged, because their SD and CV are unreliable.

## Checks before reporting

- Read the Warnings section of `report.md` first and pass every warning on to the user.
- A CV above about 30% for yield in a uniform trial usually means data-entry or unit problems. Ask the user before interpreting it.

## References

- Python Software Foundation. The Python Standard Library: statistics. https://docs.python.org/3/library/statistics.html
