# Changelog

## 1.0.0 (2026-10-15)

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
- 33 tested skills plus the `agriclaw` router: field and environment (5), field trials (4), genotypes and genetics (7), gene expression (4), drone and satellite imagery (6) and plant physiology (7). See the README table.
- Tests that compare results with known answers: published worked examples (FAO-56, NCH-40, the qtl2 user guide, PhotoGEA's article tables, Cerón-Rojas & Crossa 2018), simulations with known truth, and design-property checks.
- `CONVENTIONS.md`, working skill templates, `tools/validate_skills.py`, `tools/build_catalog.py`, `tools/check_dois.py`, `tests/run_r_tests.R`, `tools/install_r_deps.R`, and GitHub Actions CI.

### Problems found in upstream packages (worked around in the skills)
- statgenGWAS 1.0.14: `thrType = "fdr"` compares the k-th cluster with α/k after an unadjusted p < α gate, so traits with no marker effects get significant SNPs. `claw-gwas-crop` uses Benjamini–Hochberg instead.
- statgenGxE 1.0.11: `gxeStability(method = "static")` is documented as Shukla's stability variance but returns each genotype's variance across environments. `claw-met-gxe` computes Shukla's measure itself.
- FielDHub 1.5.0: `row_column()` fails for some seeds ("best_efficiencies" not found). `claw-trial-design` retries with the next seed and records the seed used.
- tximport: `ignoreTxVersion` cuts IDs at the first dot, which merges every `Glyma.`/`Sobic.`/`Potri.` transcript. `claw-rnaseq-import` strips only a trailing version number.
- segmented: `seg.control(seed = )` resets the caller's random-number stream. `claw-vpd-transpiration` restores it.
- terra: 8- and 16-bit rasters are written with NoData 255 / 65535 by default, so saturated pixels become missing. `claw-canopy-cover` and `claw-ortho-qc` warn about it.

## 0.1.0 (2026-03-02)

Initial scaffold: README, skill template, and four Python skills (GDD, soil, weather, trial design).
