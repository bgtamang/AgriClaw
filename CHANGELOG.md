# Changelog

## 1.0.0 (unreleased, branch `v1-rebuild`)

A rebuild of the March 2026 scaffold.

### Changed
- Skills follow the open Agent Skills specification, one self-contained folder each under `skills/`, so they load in Claude Code, Codex, Gemini CLI, Copilot, Cursor, OpenClaw and others.
- AgriClaw installs as a Claude Code plugin (`.claude-plugin/`).
- Statistics skills are written in R and wrap established packages. Field/weather data access stays in Python.

### Fixed
- `claw-soil`: the texture classifier misclassified common textures (silt loam reported as silty clay loam; sandy loam as loamy sand). It now uses SSURGO's own texture class, with a corrected fallback classifier.
- `claw-gdd`:
  - corn growth-stage thresholds in °F-days were compared against °C-days;
  - soybean reproductive stages were predicted from GDD alone;
  - in-season queries failed;
  - missing days were filled indefinitely;
  - frost dates were searched in the wrong window.
- `claw-trial-design`: the earlier Python generator produced invalid designs (broken strip-plots, pseudo alpha-lattices that dropped entries, augmented blocks without checks, split whole plots). It was replaced by an R skill built on agricolae and FielDHub, with design-property tests.
- Removed fabricated or garbled citations.

### Added
- New skills across field trials, genetics, gene expression, imagery and physiology. See the README table.
- `CONVENTIONS.md`, working skill templates, `tools/validate_skills.py`, `tools/build_catalog.py`, `tests/run_r_tests.R`, `tools/install_r_deps.R`, and GitHub Actions CI.

## 0.1.0 (2026-03-02)

Initial scaffold: README, skill template, and four Python skills (GDD, soil, weather, trial design).
