---
name: claw-licor-6800
description: Read, combine and quality-check LI-COR LI-6800 photosynthesis (gas exchange) log files, plaintext or Excel .xlsx, with PhotoGEA, and write one tidy CSV for A-Ci, light-response and stomatal-model fitting. Stops on Excel logs whose formulas were never calculated (A, gsw, Ci read as 0 or blank), records file and instrument serial numbers, flags Ci < 0, Ci > Ca, negative gsw, implausible leaf temperature, unstable points and tiny CO2 differentials, and warns when genotypes or treatments were measured on different instruments (machine confounding). Use when the user has Licor 6800 / LI6800 files or gas-exchange data from several machines and wants to load, merge, clean or QC them, or before fitting Vcmax, light curves or Ball-Berry/Medlyn models.
license: MIT
compatibility: R >= 4.2 with PhotoGEA (tested with 1.4.0; it installs openxlsx). No internet needed.
metadata:
  version: "1.0.0"
  domain: physiology
  language: R
  status: tested
  summary: "Read and QC LI-COR LI-6800 log files (Excel and plain text)"
  repository: https://github.com/bgtamang/AgriClaw
---

# LI-COR LI-6800: read, combine and QC

## What it does

Reads LI-6800 log files with PhotoGEA's `read_gasex_file`, stops with an explanation if an Excel file's formulas were never evaluated, combines files with a file name and instrument ID on every row, adds QC flags, checks whether design factors are confounded with instrument, and writes `tables/licor_tidy.csv`, the input of `claw-aci-fit`, `claw-light-response` and `claw-stomatal-models`. Rows are never dropped: flags are columns.

## When to use / when not to use

- Use first, on raw LI-6800 logs, before any fitting skill.
- LI-6400, CIRAS or already-tabulated data: PhotoGEA reads only LI-6800 logs. Rename the columns to the tidy names below (units in the names) and pass the CSV straight to the fitting skill.
- Fitting A-Ci, light-response or stomatal curves: `claw-aci-fit`, `claw-light-response`, `claw-stomatal-models`.
- Chlorophyll fluorescence parameters: `claw-fluorescence` (planned), not this skill.

## Inputs

LI-6800 logs as the instrument wrote them: plaintext files (no extension, always written) or `.xlsx` (optional). Prefer plaintext: it stores computed values, while the `.xlsx` stores formulas that only Excel evaluates.

Ask the user, never guess:

- **Which user constants identify one curve or leaf** (`--curve-cols`, e.g. `genotype,plot`). Without them per-curve checks are skipped. If the files have no user constants, ask how curves were logged (one file per curve? remarks?) before going further.
- **Which columns hold the genotypes or treatments** (`--design-cols`) for the confounding check.
- **How instruments are identified** if the files lack serial numbers (the Remarks sheet of an `.xlsx` can be deleted): a user constant such as `instrument` (`--instrument-col`).
- **Whether every leaf filled the chamber aperture.** A, E, gsw and Ci are computed with the leaf area logged in the file (`leaf_area_cm2`); if a leaf did not fill the chamber, values must be recomputed with the true area before analysis. This skill does not recompute them.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/licor_6800.R" --demo
Rscript "${CLAUDE_SKILL_DIR}/scripts/licor_6800.R" --input-dir licor_logs/ --curve-cols genotype,plot --design-cols genotype --out results/licor
Rscript "${CLAUDE_SKILL_DIR}/scripts/licor_6800.R" --input "logs/2026-07-01-0930_logdata,logs/2026-07-02-1015_logdata" --curve-cols plot --instrument-col machine
```

The demo reads two A-Ci log files shipped with PhotoGEA (`system.file("extdata", package = "PhotoGEA")`), measured on two instruments; one has no serial numbers, so the demo uses the `instrument` user constant. Options: `--tleaf-min`/`--tleaf-max` (5/45 C), `--min-delta-co2` (1 umol/mol), `--tleaf-drift` (2 C). `--help` lists all.

## Outputs

- `tables/licor_tidy.csv`: one row per logged point. Identifiers first (`source_file`, `instrument_id`, `curve_id`, user constants), then measurements with units in the names, then QC flags. Key columns:

| Column | LI-6800 | Units |
|---|---|---|
| `a_umol_m2_s`, `ci_umol_mol`, `ca_umol_mol` | A, Ci, Ca | umol m-2 s-1, umol/mol |
| `gsw_mol_m2_s`, `gbw_mol_m2_s`, `e_mol_m2_s` | gsw, gbw, E | mol m-2 s-1 |
| `tleaf_c` | TleafCnd (Tleaf if absent) | C |
| `qin_umol_m2_s`, `qabs_umol_m2_s`, `leaf_absorptance_assumed` | Qin, Qabs, alpha | umol m-2 s-1, fraction |
| `co2_r_umol_mol`, `co2_s_umol_mol`, `co2_r_sp_umol_mol`, `delta_co2_umol_mol` | CO2_r, CO2_s, CO2_r_sp, CO2_r - CO2_s | umol/mol |
| `h2o_s_mmol_mol`, `rh_cham_pct`, `vpd_leaf_kpa` | H2O_s, RHcham, VPDleaf (leaf to chamber air) | mmol/mol, %, kPa |
| `pa_kpa`, `delta_pcham_kpa`, `oxygen_pct`, `leaf_area_cm2` | Pa, DeltaPcham, Oxygen, S | kPa, %, cm2 |
| `obs`, `elapsed_s`, `time_unix_s`, `date_local` | obs, elapsed, time, date | |
| `stable_n`, `stable_total`, `console_serial`, `head_serial` | Stable, Total, file header | |

  Full list: `references/tidy_columns.md`. QC columns: `qc_missing_values`, `qc_ci_negative`, `qc_ci_gt_ca`, `qc_gsw_negative`, `qc_tleaf_out_of_range`, `qc_duplicate` (these six make `qc_hard_fail`), `qc_unstable`, `qc_small_delta_co2` (cautions), and `qc_flags` (text).
- `tables/file_summary.csv` (rows, serials, file opened, software version per file), `tables/qc_flag_counts.csv`, `tables/curve_summary.csv` (points, files, instruments, leaf-temperature range and setpoint sequence per curve), `tables/confounding_summary.csv` and `tables/confounding_counts.csv`.
- `figures/qc_overview.png`: A vs Ci coloured by instrument with hard failures marked; leaf temperature by row.
- `report.md`: warnings first, then the tables above.

## Expert decisions

- **Reader.** PhotoGEA handles Excel and plaintext logs, files closed and reopened, oxygen changed mid-log, and Unicode column names. Units of every mapped column are checked against the contract; a mismatch stops the run.
- **Unevaluated Excel formulas.** The LI-6800 writes A, gsw, Ci and other computed columns as formulas evaluated only when Excel opens and saves the file; read programmatically they are 0 or blank (documented in `?PhotoGEA::read_licor_6800_Excel`). If A, gsw and Ci are all zero or empty while raw CO2/H2O readings are present, the run stops and recommends the plaintext log (or Excel > Formulas > Calculate Now, then save). PhotoGEA's own check looks only for zeros; this one also catches blanks.
- **Instrument ID** = `<Console s/n>/<Head s/n>` from the file header (the head holds the IRGAs), else `--instrument-col`, else `unknown:<file>` with a warning.
- **Twins and duplicates.** The LI-6800 can save each log twice (plaintext `X` and `X.xlsx`); when both are given (or found by `--input-dir`) the plaintext is read and the `.xlsx` skipped, with a note. Rows that repeat an observation already read from another file (same instrument and timestamp) are flagged `qc_duplicate`, a hard failure, so fits use the first copy. The check needs a known instrument ID, because two instruments can log in the same second.
- **Hard failures** are physically impossible or broken readings: missing A/Ci/gsw; Ci < 0; Ci > Ca while A > 0 (a leaf taking up CO2 must have Ci below Ca; von Caemmerer & Farquhar 1981); gsw < 0; leaf temperature outside 5-45 C (rule of thumb for a detached or faulty thermocouple). Ci > Ca with A < 0 (below the compensation point) is legitimate and not flagged.
- **Cautions** (rules of thumb): `qc_unstable` when fewer stability criteria were met than defined (`Stable < Total`; PhotoGEA articles drop such points for stomatal models); `qc_small_delta_co2` when |CO2_r - CO2_s| < 1 umol/mol, where A rests on a difference near analyser noise and chamber leaks weigh more (Long & Bernacchi 2003). It is expected near the CO2 compensation point of an A-Ci curve. Leaf-temperature drift within a curve > 2 C is reported.
- **Confounding** is counted in curves (rows if no `--curve-cols`): *completely confounded* when no instrument measured more than one level; *partially confounded* when some level was never measured on some instrument; *crossed* otherwise, with a warning when Cramer's V >= 0.5 (rule of thumb). The fix is design: balance instruments across genotypes/treatments, or block on instrument (fixed block, random effect or covariate). Complete confounding has no statistical fix.
- **Oxygen.** Rows with oxygen other than 21% are reported, because A-Ci fits use the logged value.

## Checks before reporting

- Read the Warnings section of `report.md` first and pass every warning on.
- **Stop and tell the user** when: the run stops on unevaluated formulas (point them to the plaintext logs); rows are flagged `qc_duplicate` (the same log was read twice, or files were merged by hand); a design factor is completely or partially confounded with instrument; a curve ID spans several files (the curve columns do not identify single curves); instrument IDs are `unknown` while several machines were used; oxygen is not 21% unexpectedly; many points fail hard checks in one file (often a thermocouple, leak or leaf-area problem).
- Look at `figures/qc_overview.png`: one instrument's points offset from the others for the same material suggests a calibration difference; ask whether the instruments were matched or cross-checked on the same leaves.
- Hard-failure points are excluded by the fitting skills by default; say so.

## References

- Lochocki EB, Salesse-Smith CE, McGrath JM (2025) PhotoGEA: an R package for closer fitting of photosynthetic gas exchange data with non-Gaussian confidence interval estimation. Plant, Cell & Environment 48: 5104-5119. https://doi.org/10.1111/pce.15501
- Long SP, Bernacchi CJ (2003) Gas exchange measurements, what can they tell us about the underlying limitations to photosynthesis? Procedures and sources of error. Journal of Experimental Botany 54: 2393-2401. https://doi.org/10.1093/jxb/erg262
- von Caemmerer S, Farquhar GD (1981) Some relationships between the biochemistry of photosynthesis and the gas exchange of leaves. Planta 153: 376-387. https://doi.org/10.1007/BF00384257
