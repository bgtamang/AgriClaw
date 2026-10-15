# Contributing to AgriClaw

Thanks for helping. AgriClaw is only as good as the expertise in its skills, so contributions from people who run these analyses for real are especially welcome.

## What makes a good skill

A skill isn't a thin wrapper around a function. It captures the decisions an experienced researcher makes and a newcomer, or an AI agent, gets wrong:

- which design or model to use, and when not to use it;
- the thresholds and defaults that matter, with their sources;
- the questions to ask before analysing (units, design, assembly, …);
- the diagnostics to check before believing a result.

Ask yourself: "What do I know about this analysis that a first-year graduate student would get wrong?"

## Adding a skill

1. Read [CONVENTIONS.md](CONVENTIONS.md). It covers the frontmatter, layout, outputs, tests and citation rules.
2. Copy `templates/r-skill` or `templates/python-skill` to `skills/claw-<your-skill>` and rename it.
3. Write the script so it wraps established packages, validates its inputs and fails loudly.
4. Add small demo data and a `--demo` option.
5. Write tests that check results against known answers (simulation with known truth, a published worked example, or property checks), plus one end-to-end `--demo` run.
6. Run the checks:
   ```bash
   python tools/validate_skills.py claw-<your-skill>
   python -m pytest skills/claw-<your-skill>            # Python skills
   Rscript tests/run_r_tests.R claw-<your-skill>        # R skills
   python tools/build_catalog.py                        # refresh catalog.json and the README table
   python tools/check_dois.py claw-<your-skill>         # look up every cited DOI (needs internet)
   ```
7. Open a pull request saying what expert knowledge the skill encodes and how the tests validate it.

## Fixing a skill

If a skill gives a wrong answer, please open an issue with a minimal input that shows it. Better still, send a fix with a test that would have caught it.

## Ground rules

- Every reference must be real and verified (check the DOI). No citation is better than an invented one.
- Don't claim validation you didn't do. Use `status: experimental` and say what's untested.
- Keep demo data small and freely redistributable.
- Never commit outputs; `agriclaw_output/` is git-ignored.
