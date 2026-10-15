---
name: claw-aci-fit
description: Fit C3 and C4 CO2-response (A-Ci) curves from LI-COR LI-6800 gas exchange with PhotoGEA to estimate Vcmax, J, Jmax, TPU and RL (respiration in the light) at 25 C and at leaf temperature with confidence intervals, and Vcmax, Vpmax and RL for C4 crops such as maize and sorghum (von Caemmerer model). Flags curves with too few low- or high-Ci points, limitations supported by too few points, a drop in A when CO2 returns to the reference level (incomplete induction or non-steady state, which can masquerade as down-regulation), stomatal closure at high CO2, negative respiration and poor fits. Use when the user has A-Ci or CO2-response curves and asks for Vcmax, Jmax, photosynthetic capacity, Farquhar / FvCB parameters, biochemical limitation or acclimation to elevated CO2.
license: MIT
compatibility: R >= 4.2 with PhotoGEA (tested with 1.4.0; it installs DEoptim and dfoptim). Base R parallel for --cores. No internet needed.
metadata:
  version: "1.0.0"
  domain: physiology
  language: R
  status: tested
  summary: "Fit C3 and C4 A-Ci curves: Vcmax, J, TPU, Vpmax (PhotoGEA)"
  repository: https://github.com/bgtamang/AgriClaw
---

# A-Ci curve fitting (C3 FvCB and C4)

## What it does

Fits each A-Ci curve with PhotoGEA's maximum-likelihood fitters (`fit_c3_aci`, `fit_c4_aci`) and reports parameters normalized to 25 C and at leaf temperature with relative-likelihood confidence intervals, sets parameters PhotoGEA judges unreliable to NA, computes Jmax from J, and runs curve QC. C3 parameters are **apparent**: the fits assume infinite mesophyll conductance (Cc = Ci).

## When to use / when not to use

- Use for A-Ci (CO2-response) curves already read and QC'd by `claw-licor-6800`.
- Raw LI-6800 files: run `claw-licor-6800` first.
- Light-response (A-Q) curves: `claw-light-response`. Stomatal conductance models: `claw-stomatal-models`.
- Mesophyll conductance, C3 variable-J fits or TDL isotope data: use PhotoGEA directly (`fit_c3_variable_j`, `calculate_gm_busch`); not covered here.
- Comparing parameters between genotypes or treatments: take `tables/aci_parameters.csv` to `claw-trial-analysis`, with the curve (leaf) as the experimental unit and instrument as a block if it varied.

## Inputs

The tidy CSV from `claw-licor-6800` (`tables/licor_tidy.csv`), or any CSV with these columns:

| Column | Units | Needed for |
|---|---|---|
| `a_umol_m2_s`, `ci_umol_mol`, `ca_umol_mol` | umol m-2 s-1, umol/mol | required |
| `tleaf_c` | C | required (temperature responses) |
| `pa_kpa`, `delta_pcham_kpa` | kPa | required (partial pressures); or `--pressure-kpa` |
| `oxygen_pct` | % | required; or `--oxygen-pct 21` if measured in air |
| curve columns | | required, named with `--curve-cols` |
| `co2_r_sp_umol_mol` | umol/mol | recommended: detects return-to-reference points (falls back to `co2_r_umol_mol`, then `ca_umol_mol`) |
| `gsw_mol_m2_s` | mol m-2 s-1 | stomatal-closure check |
| `qin_umol_m2_s` | umol m-2 s-1 | Jmax (C3) |
| `obs`, `elapsed_s` or `time_unix_s` | | measurement order (else row order) |
| `qc_hard_fail`, `qc_unstable` | | exclusions |

Ask the user, never guess:

- **C3 or C4** (`--pathway`): soybean, wheat, rice, cotton are C3; maize, sorghum, sugarcane, millets are C4. Ask if the species is not stated.
- **Which columns identify one curve** (`--curve-cols`).
- **Measurement protocol**: did each curve start at a reference CO2, go down, return, then go up? The induction check needs the return; one-direction curves get no induction ratio.
- **Growth CO2** if they want A at the operating point (`--ca-atmospheric`, e.g. 420 ambient or 550-600 in FACE).

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/aci_fit.R" --demo --cores 4
Rscript "${CLAUDE_SKILL_DIR}/scripts/aci_fit.R" --demo --pathway C4
Rscript "${CLAUDE_SKILL_DIR}/scripts/aci_fit.R" --input results/licor/tables/licor_tidy.csv --curve-cols genotype,plot --cores 4 --out results/aci
Rscript "${CLAUDE_SKILL_DIR}/scripts/aci_fit.R" --input maize_tidy.csv --curve-cols plot --pathway C4 --ca-atmospheric 420
```

The default optimizer (DEoptim, 200 generations) takes roughly 5-15 s per curve; `--cores N` fits curves in parallel with identical results. The demo fits four synthetic C3 curves generated with known parameters (one with a loss of activation after the low-CO2 steps; see `demo/README.md`). `--help` lists all options.

## Outputs

- `tables/aci_parameters.csv`, one row per curve: curve metadata (genotype, plot, instrument...), then
  - C3: `vcmax_25_umol_m2_s`, `j_25_umol_m2_s`, `jmax_25_umol_m2_s`, `tp_25_umol_m2_s`, `rl_25_umol_m2_s`, each with `_lower`/`_upper` confidence limits and `_tleaf` (value at the curve's mean leaf temperature); `*_trust`; `alpha_old`; `n_ac_limiting`, `n_aj_limiting`, `n_ap_limiting`; the 25 C constants used (`gamma_star_25_umol_mol`, `kc_25_umol_mol`, `ko_25_mmol_mol`).
  - C4: `vcmax_25_umol_m2_s`, `vpmax_25_umol_m2_s`, `rl_25_umol_m2_s` with limits and `_tleaf`; `gmc_25_mol_m2_s_bar` and `vpr_umol_m2_s` as fixed.
  - Fit and QC: `rmse_umol_m2_s`, `aic`, `seed`, `n_points_used`, `n_ci_below_low`, `n_ci_above_high`, `max_ci_umol_mol`, `induction_ratio` (A at the last return point / A at the starting point), `gsw_high_co2_ratio`, `qc_flags`, `fit_error`.
- `tables/aci_fitted_points.csv` (measured and fitted A, residual, limiting process per point), `tables/point_use.csv` (which points were used and why others were excluded).
- `figures/aci_fits_NN.png` (9 curves per page: points, fitted A, dashed Ac/Aj/Ap or C4 rates; open circles excluded, triangles recovery points), `figures/residuals.png`.
- `report.md`: warnings first, then the parameter table, method and limitations.

NA in a parameter column means PhotoGEA found the estimate unreliable: the process never limited, or the upper confidence limit was infinite. That is normal for J or TPU on many curves; report it, do not fill it.

## Expert decisions

- **Model and fitter.** PhotoGEA's FvCB implementation (Farquhar et al. 1980; equations as in von Caemmerer 2000, with RuBP-regeneration coefficients 4 and 8) and its C4 model (von Caemmerer 2000, 2021), fitted by maximum likelihood with DEoptim (PhotoGEA's default, chosen because the min() of limiting rates defeats gradient optimizers). On PhotoGEA's own example curves the fast `nmkb` optimizer stopped in worse local minima on 3 of 6 curves (e.g. Vcmax 130.6 instead of 103.8, higher RMSE), so it is offered only for quick looks.
- **Defaults fitted (C3):** Vcmax, J, TPU, RL and `alpha_old` (PhotoGEA defaults); gm infinite, so parameters are apparent. **C4:** Vcmax, Vpmax and RL, with J and Vpr fixed at 1000 (Rubisco-limited plateau assumed), gm 1 mol m-2 s-1 bar-1 at 25 C, bundle-sheath conductance 0.003, Rm_frac 0.5 (PhotoGEA defaults from von Caemmerer 2000).
- **Temperature responses (C3), `--temp-params`.** Default `bernacchi`: Bernacchi et al. (2001) Kc, Ko, Gamma*, Vcmax and RL, derived in vivo on an intercellular-CO2 basis, which matches Cc = Ci; J from Bernacchi et al. (2003); TPU from Sharkey et al. (2007). `sharkey` uses chloroplastic constants (Bernacchi et al. 2002, as in Sharkey et al. 2007) intended for fits that estimate gm; PhotoGEA's C3 article uses it. Parameters from different sets are not comparable. C4 uses `c4_temperature_param_vc` (von Caemmerer 2021, estimated for Setaria viridis).
- **Confidence intervals:** PhotoGEA relative-likelihood intervals with threshold 0.147 (about 95%). Unreliable estimates are set to NA (`remove_unreliable_param = 2`).
- **Jmax** from J with PhotoGEA `calculate_jmax` and Bernacchi et al. (2003) alpha_J (0.293) and theta_J (0.979) at the curve's mean incident PPFD. Jmax depends on these assumed values and is NA when J >= alpha_J x PPFD.
- **Recovery points.** Points logged after CO2 returned to the starting level are excluded from the fit (PhotoGEA articles do the same) but used for the induction ratio. `--keep-recovery-points` keeps them.
- **Seeds.** Curve i (sorted IDs) uses seed + i - 1, so results do not depend on `--cores` or on which other curves are in the file.
- **QC rules of thumb (adjustable):** at least 3 points below Ci 200 and above Ci 600 umol/mol for C3 (100 and 400 for C4); at least 3 points limited by each process; induction ratio outside 0.9-1.1; gsw at the highest CO2 below 0.5 x its starting value; RMSE above 5% of the curve's maximum A; negative RL. Each sets a flag, never a deletion.
- **Point exclusions:** `qc_hard_fail` points are excluded (from `claw-licor-6800`); unstable points are kept unless `--exclude-unstable` (PhotoGEA advises fitting first and removing points only when they cause problems).

## Validation

- **Published values:** reproduces the parameter tables of PhotoGEA's articles "Analyzing C3 A-Ci Curves" and "Analyzing C4 A-Ci Curves" (version 1.4.0) from PhotoGEA's example files: all six C3 Vcmax values within 0.5% (the first exactly), confidence limits within 1%, J, Jmax and TPU including which curves give NA; all five C4 Vpmax values within 0.5% and the species-mean Vcmax and RL (`tests/test_aci_fit_photogea.R`).
- **Simulation:** 20 C3 curves simulated with PhotoGEA's own model (Vcmax 60-140, leaf temperature 24-33 C, A noise sd 0.4, Ci noise sd 2) recovered Vcmax within 3.1% (median about 1%), J within 0.9% and TPU within 4.2% where those processes limited at least 3 points (`tests/test_aci_fit_recovery.R`); noise-free C3 and C4 curves are recovered exactly.
- **Not validated:** C4 recovery from noisy curves (C4 Vpmax and RL are sensitive to errors in the few low-Ci points: Ci noise of 1.5 umol/mol biased Vpmax by about 10% in a check run), `alpha_old`, and the operating-point estimate.

## Checks before reporting

- Read the Warnings in `report.md` and look at every page of `figures/aci_fits_NN.png`: fitted lines should pass through the points and residuals should show no trend.
- **Stop and tell the user** when:
  - `induction_drop` or `induction_rise` appears, especially if it differs between treatments: a lower J or TPU in one group may be an artifact of loss of activation, not down-regulation;
  - J (or TPU) is NA for most curves but the user wants to compare J or Jmax: the curves do not support it (more high-Ci points, or measurements at saturating light, are needed);
  - `few_low_ci_points` affects curves used for Vcmax comparisons;
  - `negative_rl`, `poor_fit_rmse` or `fit_failed` appear;
  - groups differ in mean leaf temperature by more than a few degrees: the 25 C values then rest heavily on the assumed temperature responses;
  - the `claw-licor-6800` report showed treatment confounded with instrument.
- When comparing with literature, state that values are apparent (Ci-based), the temperature-response set and the PhotoGEA version.

## References

- Farquhar GD, von Caemmerer S, Berry JA (1980) A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species. Planta 149: 78-90. https://doi.org/10.1007/BF00386231
- von Caemmerer S (2000) Biochemical Models of Leaf Photosynthesis. CSIRO Publishing, Collingwood. https://doi.org/10.1071/9780643103405
- von Caemmerer S (2021) Updating the steady-state model of C4 photosynthesis. Journal of Experimental Botany 72: 6003-6017. https://doi.org/10.1093/jxb/erab266
- Sharkey TD, Bernacchi CJ, Farquhar GD, Singsaas EL (2007) Fitting photosynthetic carbon dioxide response curves for C3 leaves. Plant, Cell & Environment 30: 1035-1040. https://doi.org/10.1111/j.1365-3040.2007.01710.x
- Bernacchi CJ, Singsaas EL, Pimentel C, Portis AR Jr, Long SP (2001) Improved temperature response functions for models of Rubisco-limited photosynthesis. Plant, Cell & Environment 24: 253-259. https://doi.org/10.1111/j.1365-3040.2001.00668.x
- Bernacchi CJ, Portis AR, Nakano H, von Caemmerer S, Long SP (2002) Temperature response of mesophyll conductance. Implications for the determination of Rubisco enzyme kinetics and for limitations to photosynthesis in vivo. Plant Physiology 130: 1992-1998. https://doi.org/10.1104/pp.008250
- Bernacchi CJ, Pimentel C, Long SP (2003) In vivo temperature response functions of parameters required to model RuBP-limited photosynthesis. Plant, Cell & Environment 26: 1419-1430. https://doi.org/10.1046/j.0016-8025.2003.01050.x
- Long SP, Bernacchi CJ (2003) Gas exchange measurements, what can they tell us about the underlying limitations to photosynthesis? Procedures and sources of error. Journal of Experimental Botany 54: 2393-2401. https://doi.org/10.1093/jxb/erg262
- Lochocki EB, Salesse-Smith CE, McGrath JM (2025) PhotoGEA: an R package for closer fitting of photosynthetic gas exchange data with non-Gaussian confidence interval estimation. Plant, Cell & Environment 48: 5104-5119. https://doi.org/10.1111/pce.15501
