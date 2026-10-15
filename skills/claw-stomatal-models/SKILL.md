---
name: claw-stomatal-models
description: Fit the Ball-Berry (gsw = g0 + g1 A hs / Cs) and Medlyn USO (gsw = g0 + 1.6 (1 + g1 / sqrt(D)) A / Cs) stomatal conductance models per leaf or curve from LI-COR LI-6800 gas exchange with PhotoGEA, computing leaf-surface humidity, CO2 and VPD from the logged data, and report intrinsic water-use efficiency (iWUE = A / gsw). Stops on unit mistakes (VPD in Pa or hPa, gsw in mmol, relative humidity as a fraction, CO2 in Pa) and flags implausible g1, negative g0 and narrow data ranges. Use when the user asks for stomatal slope g1, g0, Ball-Berry or Medlyn parameters, stomatal sensitivity, iWUE or water-use efficiency from gas-exchange curves, for crop models or genotype comparisons.
license: MIT
compatibility: R >= 4.2 with PhotoGEA (tested with 1.4.0). No internet needed.
metadata:
  version: "1.0.0"
  domain: physiology
  language: R
  status: tested
  repository: https://github.com/bgtamang/AgriClaw
---

# Stomatal conductance models (Ball-Berry, Medlyn) and iWUE

## What it does

Computes leaf-surface relative humidity (hs), CO2 (Cs) and vapour pressure deficit (D) from LI-6800 gas exchange with PhotoGEA `calculate_gas_properties`, fits the Ball-Berry model by linear regression (`fit_ball_berry`) and the Medlyn model by nonlinear least squares (`fit_medlyn`) for each curve, and reports g0 and g1 with standard errors, fit statistics and intrinsic water-use efficiency.

## When to use / when not to use

- Use for "Ball-Berry curves": steady-state points on one leaf where light, CO2 or humidity was stepped so that A, hs, Cs or D varied, read by `claw-licor-6800`.
- Raw LI-6800 files: run `claw-licor-6800` first.
- A-Ci or light-response parameters: `claw-aci-fit`, `claw-light-response`.
- Transpiration responses to VPD with breakpoints (whole plants): `claw-vpd-transpiration` (planned).
- Stomatal kinetics (speed of opening and closing): not covered.

## Inputs

The tidy CSV from `claw-licor-6800`, or any CSV with:

| Column | Units | Needed for |
|---|---|---|
| `a_umol_m2_s`, `gsw_mol_m2_s` | umol m-2 s-1, mol m-2 s-1 (water vapour) | required |
| curve columns | | required, named with `--curve-cols` |
| `e_mol_m2_s`, `gbw_mol_m2_s`, `h2o_s_mmol_mol`, `ca_umol_mol`, `tleaf_c`, `pa_kpa`, `delta_pcham_kpa` | mol m-2 s-1, mmol/mol, umol/mol, C, kPa | to compute hs, Cs and D (LI-6800 data) |
| or `rh_leaf_pct`, `cs_umol_mol`, `vpd_leaf_surface_kpa` | %, umol/mol, kPa | precomputed leaf-surface values (other instruments) |
| `vpd_leaf_kpa` | kPa | `--vpd-source licor` |
| `qc_hard_fail`, `qc_unstable` | | exclusions |

Ask the user, never guess:

- **Which columns identify one curve or leaf** (`--curve-cols`).
- **Which D they want for Medlyn g1** if they compare with published values: leaf-surface VPD (default) or the LI-6800 `VPDleaf`, which is leaf to chamber air and higher.
- **Whether every point was at stomatal steady state** (wait time per step); stomata often need 10-30 minutes.
- **Units** whenever columns come from another instrument or a spreadsheet: gsw in mol (not mmol) m-2 s-1, D in kPa, hs in %, Cs in umol/mol.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/stomatal_models.R" --demo
Rscript "${CLAUDE_SKILL_DIR}/scripts/stomatal_models.R" --input results/licor/tables/licor_tidy.csv --curve-cols genotype,plot --out results/stomata
Rscript "${CLAUDE_SKILL_DIR}/scripts/stomatal_models.R" --input licor_tidy.csv --curve-cols plot --models medlyn --vpd-source licor
```

The demo fits four synthetic light series generated with the Medlyn model and known g1 (`demo/README.md`). Fits take a fraction of a second per curve.

## Outputs

- `tables/stomatal_parameters.csv`, one row per curve: curve metadata, `n_points_used`, then
  - Ball-Berry: `bb_g0_mol_m2_s`, `bb_g0_se`, `bb_g1` (dimensionless), `bb_g1_se`, `bb_r2`, `bb_p_value`, `bb_rmse_mol_m2_s`, `bb_index_cv`;
  - Medlyn: `medlyn_g0_mol_m2_s`, `medlyn_g0_se`, `medlyn_g1_kpa05` (kPa^0.5), `medlyn_g1_se`, `medlyn_r2`, `medlyn_rmse_mol_m2_s`;
  - `iwue_mean_umol_mol`, `iwue_sd_umol_mol`, `a_mean_umol_m2_s`, `gsw_mean_mol_m2_s`, `vpd_mean_kpa`, `vpd_source`, `qc_flags`, `fit_error`.
- `tables/stomatal_points.csv`: per point hs, Cs, D, Ball-Berry and Medlyn indices, fitted gsw from each model, iWUE, and whether the point was used.
- `figures/stomatal_fits.png`: gsw against each model's index, one colour per curve (open circles excluded).
- `report.md`: warnings first, results, method, limitations.

## Expert decisions

- **Models.** Ball-Berry with an intercept: gsw = g0 + g1 A hs / Cs (Ball, Woodrow & Berry 1987; intercept form as in Collatz et al. 1991), g1 dimensionless with hs as a fraction. Medlyn: gsw = g0 + 1.6 (1 + g1 / sqrt(D)) A / Cs (Medlyn et al. 2011, with the 2012 corrigendum), gsw for water vapour, D in kPa, g1 in kPa^0.5. Both fitted by PhotoGEA (`lm`; `nls` started at g1 = 4, g0 = 0.005).
- **Leaf-surface values.** hs, Cs and the leaf-surface water vapour come from PhotoGEA `calculate_gas_properties` (LI-6800 manual Appendix C; von Caemmerer & Farquhar 1981), using gbw as logged. D defaults to the leaf-surface VPD, SVP(Tleaf) x (1 - hs), consistent with hs and Cs. The LI-6800 `VPDleaf` is leaf to chamber air and was 0.03-0.29 kPa higher in PhotoGEA's example data; `--vpd-source licor` uses it, and the report gives the mean difference. g1 values computed with different D definitions are not comparable.
- **Steady state.** Points logged before the stability criteria were met are excluded by default (PhotoGEA's Ball-Berry article does the same); `--keep-unstable` keeps them.
- **iWUE** = A / gsw (umol CO2 per mol H2O; Leakey et al. 2019), per point and as the curve mean of used points. It is meaningful only at stomatal steady state; within an A-Ci or light curve usually only the first point qualifies.
- **Unit checks (stop before fitting).** gsw above 3 mol m-2 s-1 (rule of thumb: mmol units); D above 10 kPa (physically impossible for leaves: SVP is about 9.6 kPa at 45 C, so D is in Pa or hPa); median D below 0.02 kPa; hs at most 1 (a fraction, not %); median Cs below 50 umol/mol (Pa); computed hs outside 0-100% (E in mmol, or wrong gbw or H2O units).
- **Flags (rules of thumb):** Ball-Berry g1 outside 0-30; Medlyn g1 outside 0-20 kPa^0.5 (compare with the compilation by Lin et al. 2015 and the review by Miner et al. 2017); negative g0; Ball-Berry index CV below 10% (slope poorly constrained); Ball-Berry R2 below 0.5; fewer than 4 points (`--min-points`) means no fit.

## Validation

- Reproduces the Ball-Berry intercepts, slopes and R2 of the five curves in PhotoGEA's article "Analyzing Ball-Berry Data" (version 1.4.0) from PhotoGEA's example files to 1e-6, through this skill's tidy-CSV path.
- Simulated LI-6800 states (gsw solved so the model holds at the leaf surface): hs, Cs and D recovered from E, gbw, H2O_s, Ca, leaf temperature and pressure to 1e-9; near-noise-free curves give g0 and g1 to 1e-4.
- 10 curves per model (8 points, gsw scatter sd 0.01 mol m-2 s-1): Ball-Berry g1 within 6.6% and Medlyn g1 within 11% of the truth, every error within 3.5 standard errors; g0 within 0.02 mol m-2 s-1.
- Unit mistakes caught: D in Pa and hPa, gsw in mmol, hs as a fraction, Cs in Pa, E in mmol, LI-6800 VPD in Pa.
- Not validated: Medlyn fits with g0 fixed at 0 (not offered), other stomatal models (Leuning, optimisation-based), and data from instruments other than the LI-6800 beyond the column contract.

## Checks before reporting

- Read the Warnings in `report.md` and look at `figures/stomatal_fits.png`: points should rise roughly linearly with each index.
- **Stop and tell the user** when: the run stops on units; `bb_g1_implausible` or `medlyn_g1_implausible` appears; `narrow_index_range` or `low_r2` affects curves being compared (g1 is not reliable); many points were excluded as unstable; groups differ in mean D or leaf temperature (g1 comparisons across very different D are confounded with the model's assumptions).
- `negative_g0` is common and means g0 is an extrapolation; do not interpret it as night-time or minimum conductance.
- When quoting g1, state the model, the D definition and that g0 was fitted.

## References

- Ball JT, Woodrow IE, Berry JA (1987) A model predicting stomatal conductance and its contribution to the control of photosynthesis under different environmental conditions. In: Biggins J (ed) Progress in Photosynthesis Research, vol 4, pp 221-224. Martinus Nijhoff, Dordrecht. https://doi.org/10.1007/978-94-017-0519-6_48
- Collatz GJ, Ball JT, Grivet C, Berry JA (1991) Physiological and environmental regulation of stomatal conductance, photosynthesis and transpiration: a model that includes a laminar boundary layer. Agricultural and Forest Meteorology 54: 107-136. https://doi.org/10.1016/0168-1923(91)90002-8
- Medlyn BE, Duursma RA, Eamus D, Ellsworth DS, Prentice IC, Barton CVM, Crous KY, de Angelis P, Freeman M, Wingate L (2011) Reconciling the optimal and empirical approaches to modelling stomatal conductance. Global Change Biology 17: 2134-2144. https://doi.org/10.1111/j.1365-2486.2010.02375.x
- Medlyn BE, Duursma RA, Eamus D, Ellsworth DS, Prentice IC, Barton CVM, Crous KY, de Angelis P, Freeman M, Wingate L (2012) Reconciling the optimal and empirical approaches to modelling stomatal conductance (corrigendum). Global Change Biology 18: 3476. https://doi.org/10.1111/j.1365-2486.2012.02790.x
- Lin YS, Medlyn BE, Duursma RA, et al. (2015) Optimal stomatal behaviour around the world. Nature Climate Change 5: 459-464. https://doi.org/10.1038/nclimate2550
- Miner GL, Bauerle WL, Baldocchi DD (2017) Estimating the sensitivity of stomatal conductance to photosynthesis: a review. Plant, Cell & Environment 40: 1214-1238. https://doi.org/10.1111/pce.12871
- Leakey ADB, Ferguson JN, Pignon CP, Wu A, Jin Z, Hammer GL, Lobell DB (2019) Water use efficiency as a constraint and target for improving the resilience and productivity of C3 and C4 crops. Annual Review of Plant Biology 70: 781-808. https://doi.org/10.1146/annurev-arplant-042817-040305
- von Caemmerer S, Farquhar GD (1981) Some relationships between the biochemistry of photosynthesis and the gas exchange of leaves. Planta 153: 376-387. https://doi.org/10.1007/BF00384257
- Lochocki EB, Salesse-Smith CE, McGrath JM (2025) PhotoGEA: an R package for closer fitting of photosynthetic gas exchange data with non-Gaussian confidence interval estimation. Plant, Cell & Environment 48: 5104-5119. https://doi.org/10.1111/pce.15501
