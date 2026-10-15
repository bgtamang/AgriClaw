---
name: claw-fieldbook
description: Move field-trial data between the design and the analysis. Exports a trial field book to a Field Book app (PhenoApps) field import CSV, and imports collected data (Field Book database or table exports, or any spreadsheet, from one or several tablets) back onto the design with QC - unknown or duplicated plot IDs, plots or traits with missing data, out-of-range values, robust (MAD) outliers, decimal commas and date problems - then writes a clean, analysis-ready table. Use when the user wants to load a trial into Field Book, merge phenotype data from tablets, check or clean field data before analysis, or find data-entry errors.
license: MIT
compatibility: R >= 4.2, base R only (openssl optional, for SHA-256 checksums). No internet needed.
metadata:
  version: "1.0.0"
  domain: trials
  language: R
  status: tested
  summary: "Field Book app export and import, merging tablets, plot-level data QC"
  repository: https://github.com/bgtamang/AgriClaw
---

# Field Book export, data import and QC

## What it does

`--mode export` turns a field book (one row per plot, from `claw-trial-design`) into a Field Book app import file: a unique plot ID plus primary and secondary order columns, with headers and IDs that Field Book accepts. `--mode import` reads collected data, merges it onto the design by plot ID, runs QC, and writes `clean_data.csv` (one row per plot, one column per trait) for `claw-trial-analysis`, plus a table of every flag and the action taken.

## When to use / when not to use

- Use before a season (export) and after data collection (import, QC), including when data come from several tablets or from a spreadsheet instead of Field Book.
- Use `claw-trial-design` to create the field book first, and `claw-trial-analysis` after QC for model-based analysis. Its outlier test, on residuals after removing genotype and spatial effects, is more sensitive than the screen here.
- Don't use it to fix the design itself (wrong randomization, missing plots in the design): regenerate the design instead.

## Inputs

- `--design FILE.csv`: one row per plot with a unique ID column (`plot_id` by default, `--design-id`), and `row` and `col` (needed for export and for the QC maps). `is_filler` is used if present.
- `--data FILES`: comma-separated. The format is detected per file:
  - Field Book **database** export: one row per observation; columns are the exported field columns, then `trait, value, timestamp, person, location, number` (newer versions add `attached_*` and `device_name`).
  - Field Book **table** export, or any **wide** CSV: one row per plot, one column per trait. Design columns repeated in the file are ignored.
  - The plot-ID column is `--id-col`, or else the design ID column name, or else the column whose values match the design IDs.
- `--limits`: plausible range per trait, as a CSV (`trait,min,max`) or inline (`'yield_kg_ha=0:12000,height_cm=20:200'`).
- `--planting-date`, `--season-start`, `--season-end` (YYYY-MM-DD); `--today` fixes the reference date for future-date checks.

**Ask the user, never guess:**

- the units of every trait, if the column names don't contain them (`yield_kg_ha`, not `yield`);
- plausible min/max per trait, for `--limits`; without them only the robust-outlier screen runs;
- the planting date and season window, for the date checks;
- which record to keep when a plot x trait was entered twice with different values (`--duplicates`), and whether repeated measures are intended;
- the date format, if dates are not YYYY-MM-DD (for example `--date-format %m/%d/%Y`).

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/fieldbook.R" --demo
Rscript "${CLAUDE_SKILL_DIR}/scripts/fieldbook.R" --mode export --design field_book.csv --out fb_export/
Rscript "${CLAUDE_SKILL_DIR}/scripts/fieldbook.R" --mode import --design field_book.csv --data tablet1.csv,tablet2.csv --limits "yield_kg_ha=0:12000,height_cm=20:200" --planting-date 2026-05-12 --out fb_qc/
```

The demo exports `demo/design_field_book.csv` and imports `demo/fieldbook_database_export.csv`, a simulated Field Book export with 12 planted errors (listed in `demo/README.md`); all of them appear in `qc_flags.csv`.

## Outputs

- `tables/fieldbook_import.csv` (export): `plot_id`, `row`, `col` first, then all other design columns. At import in Field Book choose `plot_id` as the unique identifier, `row` as primary order and `col` as secondary order.
- `tables/clean_data.csv` (import): one row per design plot (FILLER plots dropped). It has the design columns and one column per trait: numbers for numeric traits, ISO dates for date traits (plus `<trait>_dap`, days after planting, when a planting date is given), and text otherwise. `qc_issues` lists the flags on each plot.
- `tables/qc_flags.csv`: `plot_id, trait, issue, value, detail, action` for every problem. Issues: `unknown_plot_id`, `data_on_filler`, `plot_not_observed`, `missing_value`, `trait_not_observed`, `duplicate_identical`, `duplicate_conflict`, `non_numeric_value`, `unparseable_date`, `out_of_range`, `mad_outlier`, `future_date`, `before_planting`, `outside_season`, `bad_timestamp`.
- `tables/qc_summary.csv`: per trait, observed and missing counts, flags, min, median and max.
- `figures/qc_field_maps.png`: each numeric trait on the field grid (light = low, dark = high), with missing plots crossed and flagged plots boxed in orange. Spatial patterns (a dead row, a scale changed mid-field) show up here.
- `report.md`: warnings first, then the summary and flags.

## Expert decisions

- **Field Book format.** Field files need a unique identifier plus primary and secondary order columns, chosen at import (Field Book documentation). The import code (PhenoApps Field-Book, `FieldFileObject` and `DataHelper`) rejects duplicate or empty IDs and IDs containing `/` or `\`, and strips `[ ] ` " '` from headers. The export enforces the same rules and also replaces spaces and other symbols in headers with `_`. IDs must be unique across every field on the device, which is why `plot_id` includes the trial name.
- **Database vs table export.** Field Book's table export keeps only the first value of repeated measures and drops timestamps (Field Book documentation), so prefer the database export. With it, repeats are visible and the date checks can run.
- **Missing codes.** Blank, `NA`, `N/A`, `NaN`, `NULL`, `.` and `-` count as missing.
- **Trait type.** A trait is numeric if more than half of its values parse as numbers, and a date if more than half parse as dates. The remaining values are flagged (`non_numeric_value`, `unparseable_date`) and set to NA; a decimal comma is named in the detail. Everything else is text and passes through.
- **Duplicates (`--duplicates flag`).** Identical repeats collapse to one value. With conflicting repeats there is no safe automatic choice, so the value is set to NA and flagged. `last` or `first` (by timestamp, then file order) and `mean` are explicit opt-ins. Use `mean` only for intended repeated measures.
- **Range limits.** Values outside `--limits` are physically implausible, so they are set to NA and flagged; `--keep-out-of-range` keeps them, still flagged.
- **Robust outliers.** Robust z = (x - median) / (1.4826 x MAD), across all plots of a trait, after removing out-of-range values. |z| > 3.5 is flagged, following Iglewicz & Hoaglin's modified z-score cut-off; values are **kept**, because real genotype differences produce extremes too. This is a screen for data-entry errors and a rule of thumb: Leys et al. (2013) argue for MAD over SD and suggest 2.5 as a stricter default. The screen needs at least 10 values, also a rule of thumb, and is skipped when MAD = 0 (for example score traits with one dominant value).
- **Dates.** For **date values** (for example flowering date), a date before planting or after `--today` is impossible, so it is set to NA (kept with `--keep-out-of-range`). For **observation timestamps** the same problems are only flagged, because they usually mean a wrong tablet clock rather than a wrong value. A date before planting is reported as `before_planting` only, not also as `outside_season`.
- **FILLER plots** are dropped from `clean_data.csv`; any value recorded on them is flagged.
- **Nothing is dropped silently.** Every exclusion or change to NA is a row in `qc_flags.csv` and is counted in the report.

## Checks before reporting

- Read the Warnings in `report.md` and pass every one on.
- `unknown_plot_id`: more than a few means the wrong design file or a different trial. Stop and ask before analyzing.
- `plot_not_observed` or `trait_not_observed`: ask whether the plots were lost (record why) or the data are in another file.
- `duplicate_conflict`: resolve with the data collector; don't pick a value yourself unless the user chooses a policy.
- `mad_outlier`: check against field notes (lodging, damage, entry errors). Don't delete them; `claw-trial-analysis` re-tests outliers on model residuals.
- Look at `qc_field_maps.png` for spatial runs of missing or extreme values, such as a planter skip or a scale recalibrated mid-field.
- Confirm the trait units before passing `clean_data.csv` to `claw-trial-analysis`.

## References

- Rife TW, Poland JA (2014). Field Book: An open-source application for field data collection on Android. Crop Science 54(4):1624-1627. https://doi.org/10.2135/cropsci2013.08.0579
- Field Book documentation (PhenoApps): fields and export formats. https://field-book.readthedocs.io/ ; source code: https://github.com/PhenoApps/Field-Book
- Leys C, Ley C, Klein O, Bernard P, Licata L (2013). Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median. Journal of Experimental Social Psychology 49(4):764-766. https://doi.org/10.1016/j.jesp.2013.03.013
- Iglewicz B, Hoaglin DC (1993). How to Detect and Handle Outliers. ASQC Basic References in Quality Control: Statistical Techniques, vol. 16. Cut-off as summarized in the NIST/SEMATECH e-Handbook of Statistical Methods, section 1.3.5.17: https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
