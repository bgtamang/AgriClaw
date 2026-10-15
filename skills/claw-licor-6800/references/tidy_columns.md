# Tidy CSV columns (claw-licor-6800 output)

`tables/licor_tidy.csv` has one row per logged LI-6800 point. It is the input
contract for `claw-aci-fit`, `claw-light-response` and `claw-stomatal-models`.
Data from other instruments can use the same column names and units.

Columns appear only when the source column exists in the log. User constants
(for example `genotype`, `plot`, `instrument`) keep their names; a user
constant whose name collides with a tidy name gets the prefix `uc_`.

## Identifiers and metadata

| Column | Source | Notes |
|---|---|---|
| `source_file` | file name | |
| `instrument_id` | `<Console s/n>/<Head s/n>`, `--instrument-col`, or `unknown:<file>` | |
| `instrument_id_source` | | how `instrument_id` was set |
| `curve_id` | `--curve-cols` joined with " - " | only with `--curve-cols` |
| `file_type` | | `Excel` or `plaintext` |
| `console_serial`, `head_serial` | file header | NA when the header lacks them |
| `obs` | obs | observation number in the log |
| `time_unix_s` | time | seconds since 1970-01-01 |
| `elapsed_s` | elapsed | s since the log was opened |
| `date_local`, `hhmmss` | date, hhmmss | instrument clock, as text |
| `user_remark` | latest user remark | |

## Gas exchange

| Column | LI-6800 column | Units |
|---|---|---|
| `a_umol_m2_s` | A | umol m-2 s-1 |
| `ci_umol_mol` | Ci | umol mol-1 |
| `ca_umol_mol` | Ca | umol mol-1 |
| `gsw_mol_m2_s` | gsw | mol m-2 s-1 |
| `gbw_mol_m2_s` | gbw | mol m-2 s-1 |
| `gtc_mol_m2_s` | gtc | mol m-2 s-1 |
| `e_mol_m2_s` | E | mol m-2 s-1 |
| `tleaf_c` | TleafCnd (Tleaf when TleafCnd is absent) | C |
| `tleaf_thermocouple_c` | Tleaf | C |
| `tair_c` | Tair | C |
| `qin_umol_m2_s` | Qin (incident PPFD) | umol m-2 s-1 |
| `qabs_umol_m2_s` | Qabs (absorbed, with assumed absorptance) | umol m-2 s-1 |
| `leaf_absorptance_assumed` | alpha | fraction |
| `co2_r_umol_mol`, `co2_s_umol_mol` | CO2_r, CO2_s | umol mol-1 |
| `co2_r_sp_umol_mol` | CO2_r_sp (setpoint) | umol mol-1 |
| `delta_co2_umol_mol` | CO2_r - CO2_s | umol mol-1 |
| `h2o_r_mmol_mol`, `h2o_s_mmol_mol` | H2O_r, H2O_s | mmol mol-1 |
| `rh_cham_pct` | RHcham | % |
| `vpd_leaf_kpa` | VPDleaf (leaf to chamber air) | kPa |
| `svp_leaf_kpa` | SVPleaf | kPa |
| `pa_kpa`, `delta_pcham_kpa` | Pa, DeltaPcham | kPa |
| `oxygen_pct` | Oxygen | % |
| `flow_umol_s` | Flow | umol s-1 |
| `leaf_area_cm2` | S | cm2 |
| `stomatal_ratio_k` | K | |
| `stable_n`, `stable_total` | Stable, Total | stability criteria met / defined |
| `a_ok`, `gsw_ok` | A:OK, gsw:OK | 1 = stable |

## QC flags

| Column | TRUE when | Kind |
|---|---|---|
| `qc_missing_values` | A, Ci or gsw missing | hard |
| `qc_ci_negative` | Ci < 0 | hard |
| `qc_ci_gt_ca` | Ci > Ca while A > 0 | hard |
| `qc_gsw_negative` | gsw < 0 | hard |
| `qc_tleaf_out_of_range` | leaf temperature outside `--tleaf-min`..`--tleaf-max` (5-45 C) | hard |
| `qc_duplicate` | same instrument and timestamp as a row read earlier from another file (NA when the instrument is unknown) | hard |
| `qc_unstable` | Stable < Total (or A:OK / gsw:OK = 0) | caution |
| `qc_small_delta_co2` | abs(CO2_r - CO2_s) < `--min-delta-co2` (1 umol mol-1) | caution |
| `qc_hard_fail` | any hard flag | |
| `qc_flags` | text list of the flags that are TRUE | |

A flag is NA when its check could not be made because the column is absent.
