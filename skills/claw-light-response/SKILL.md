---
name: claw-light-response
description: Fit photosynthetic light-response curves (A-Q, A vs PPFD or PAR) from LI-COR LI-6800 gas exchange with the non-rectangular hyperbola (minpack.lm) to estimate light-saturated net photosynthesis (Asat), apparent quantum yield, curvature (theta), dark respiration (Rd), light compensation point and a light saturation point defined as the PPFD where A reaches 90% of Asat. Flags too few low-light points for the quantum yield, photoinhibition at high light, curvature at its bounds and extrapolated saturation points; reports the quantum yield on an incident-light basis unless leaf absorptance was measured. Use when the user has light curves and asks for Amax, Asat, quantum efficiency, light compensation or saturation point, or compares light responses between genotypes or treatments.
license: MIT
compatibility: R >= 4.2 with minpack.lm. No internet needed.
metadata:
  version: "1.0.0"
  domain: physiology
  language: R
  status: tested
  summary: "Light-response curves: Asat, quantum yield, compensation and saturation points"
  repository: https://github.com/bgtamang/AgriClaw
---

# Light-response (A-Q) curve fitting

## What it does

Fits each light-response curve with the non-rectangular hyperbola by Levenberg-Marquardt least squares (`minpack.lm::nlsLM`, three starting curvatures, best fit kept) and reports Asat, the apparent quantum yield, curvature, Rd, the light compensation point (LCP) and a saturation point, with standard errors and Wald 95% intervals. LCP and the saturation point are solved exactly from the fitted curve.

## When to use / when not to use

- Use for A vs incident PPFD curves read and QC'd by `claw-licor-6800` (or a CSV with the same columns).
- Raw LI-6800 files: run `claw-licor-6800` first.
- CO2-response (A-Ci) curves: `claw-aci-fit`. Stomatal models from light series: `claw-stomatal-models`.
- Respiration in the light (RL): this skill's Rd is an extrapolation through the low-light points and includes the Kok effect; for RL use a Laisk or Kok analysis (PhotoGEA `fit_laisk`).
- Curves with strong photoinhibition: the hyperbola cannot decline; see the photoinhibition flag below.

## Inputs

The tidy CSV from `claw-licor-6800`, or any CSV with:

| Column | Units | Needed for |
|---|---|---|
| `qin_umol_m2_s` | umol m-2 s-1 (incident PPFD) | required |
| `a_umol_m2_s` | umol m-2 s-1 | required |
| curve columns | | required, named with `--curve-cols` |
| `tleaf_c` | C | reported |
| `qc_hard_fail`, `qc_unstable` | | exclusions |

Ask the user, never guess:

- **Which columns identify one curve** (`--curve-cols`).
- **Was leaf absorptance measured** (integrating sphere or similar)? Only then pass `--absorptance`; the LI-6800 `Qabs` and `alpha` columns use an assumed absorptance and must not be treated as measured.
- **How long each light step lasted**, and whether the curve ran from high to low light. Short steps after high light give non-steady-state points.
- **Which definition of the saturation point** they need, if they compare with a paper (default: 90% of Asat; `--sat-fraction`).

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/light_response.R" --demo
Rscript "${CLAUDE_SKILL_DIR}/scripts/light_response.R" --input results/licor/tables/licor_tidy.csv --curve-cols genotype,plot --out results/light
Rscript "${CLAUDE_SKILL_DIR}/scripts/light_response.R" --input licor_tidy.csv --curve-cols plot --absorptance 0.86 --sat-fraction 0.75
Rscript "${CLAUDE_SKILL_DIR}/scripts/light_response.R" --input licor_tidy.csv --curve-cols plot --theta 0.7
```

The demo fits four synthetic curves with known parameters, one with photoinhibition added (`demo/README.md`). Each curve takes well under a second.

## Outputs

- `tables/light_parameters.csv`, one row per curve: curve metadata, then
  - `asat_umol_m2_s` (light-saturated net A = Pmax - Rd, the asymptote) with `_se`, `_lower`, `_upper`;
  - `phi_incident_mol_mol` (mol CO2 per mol incident photons) with SE and limits, and `phi_absorbed_mol_mol` when `--absorptance` is given;
  - `theta`, `theta_se`, `theta_fitted`; `rd_umol_m2_s` (positive = CO2 release) with SE and limits; `pmax_gross_umol_m2_s`;
  - `lcp_umol_m2_s`; `qsat_umol_m2_s` (PPFD where A = `sat_fraction` x Asat);
  - `a_model_at_qmax_umol_m2_s`, `a_max_observed_umol_m2_s`, `q_at_a_max_observed`, `q_max_umol_m2_s`;
  - `n_points`, `n_low_light_points`, `rmse_umol_m2_s`, `tleaf_mean_c`, `qc_flags`, `fit_error`.
- `tables/light_fitted_points.csv`: measured and fitted A and residuals.
- `figures/light_fits_NN.png`: 9 curves per page with the fit, LCP (x), saturation point (dashed) and Asat (dotted).
- `report.md`: warnings first, results, method, limitations.

## Expert decisions

- **Model.** Non-rectangular hyperbola (Prioul & Chartier 1977; Marshall & Biscoe 1980): A = [phi Q + Pmax - sqrt((phi Q + Pmax)^2 - 4 theta phi Q Pmax)] / (2 theta) - Rd, with phi the initial slope, Pmax the gross asymptote and 0 < theta <= 1 (theta = 1 is the Blackman response, theta -> 0 the rectangular hyperbola). It is the usual empirical model for leaf light responses; Ögren & Evans (1993) analyse what controls the curvature, and Lobo et al. (2013) compare it with other light-response models.
- **Asat** is defined as the asymptote Pmax - Rd. With low theta it is approached only far above the measured range, so the model A at the highest PPFD and the maximum observed A are also reported. Say which one you quote.
- **LCP and saturation point** are exact: gross rate g = Rd (LCP) or g = Rd + f (Pmax - Rd), and Q = g (Pmax - theta g) / (phi (Pmax - g)). The saturation point depends on the chosen fraction f; definitions vary between studies (Lobo et al. 2013), so 0.9 is an adjustable rule of thumb and is recorded in the output.
- **Quantum yield basis.** phi is per mol incident photons unless `--absorptance` is given (phi_absorbed = phi / absorptance). The theoretical maximum is 0.125 mol CO2 per mol absorbed photons (8 photons per CO2); measured C3 values are lower (Long, Postl & Bolhár-Nordenkampf 1993). A fitted phi above 0.125 is flagged.
- **Fitting.** Bounds phi 1e-4 to 0.3, theta 0.001 to 1, Pmax up to 3 x the initial guess + 50, Rd -10 to 50; starting theta 0.3, 0.7 and 0.95, best residual sum of squares kept. Standard errors come from the nlsLM covariance matrix; intervals are Wald intervals (approximate, especially for theta and for parameters near a bound).
- **QC rules of thumb (adjustable):** at least 3 points between 0 and 150 umol m-2 s-1 for phi (`--low-light`, `--min-low-points`); a point below 50 for Rd and LCP; photoinhibition when A at the highest PPFD is more than 10% below the maximum observed A at a lower PPFD (`--photoinhibition-drop`); theta at a bound; RMSE above 5% of the maximum A; saturation point beyond the measured range.
- **Exclusions:** `qc_hard_fail` points excluded; unstable points kept unless `--exclude-unstable`.

## Validation

- LCP and saturation point equal the numerical roots of the fitted curve to 1e-6 (25 random parameter sets); a noise-free curve returns phi, theta, Asat and Rd to 1e-4.
- 20 simulated curves (phi 0.04-0.075, theta 0.55-0.9, A noise sd 0.3): every estimate within 3.5 standard errors of the truth and at least 75% within 2; largest errors phi 16%, Asat 3.8%, theta 0.15, Rd 0.45, LCP 24% (phi and theta trade off); true phi inside the Wald interval on at least 75% of curves.
- The demo's three regular curves give Asat within 5% of the generating values; the photoinhibited one is flagged.
- Not validated: behaviour on curves with Kok-effect breaks at low light, and a photoinhibition model (none is fitted).

## Checks before reporting

- Read the Warnings in `report.md` and look at `figures/light_fits_NN.png`.
- **Stop and tell the user** when: `photoinhibition` appears (offer a refit without the points beyond the peak, and say Asat is biased); `few_low_light_points` or `no_dark_point` affects curves used to compare phi, Rd or LCP; `theta_at_bound` appears (suggest `--theta`, e.g. 0.7, and report it as fixed); `phi_above_theoretical_max` appears (unit or low-light measurement problem); `qsat_beyond_measured_range` affects a quoted saturation point.
- Always state the basis of phi (incident or absorbed) and the definition of Asat and of the saturation point.

## References

- Prioul JL, Chartier P (1977) Partitioning of transfer and carboxylation components of intracellular resistance to photosynthetic CO2 fixation: a critical analysis of the methods used. Annals of Botany 41: 789-800. https://doi.org/10.1093/oxfordjournals.aob.a085354
- Marshall B, Biscoe PV (1980) A model for C3 leaves describing the dependence of net photosynthesis on irradiance. Journal of Experimental Botany 31: 29-39. https://doi.org/10.1093/jxb/31.1.29
- Ögren E, Evans JR (1993) Photosynthetic light-response curves. Planta 189: 182-190. https://doi.org/10.1007/BF00195075
- Long SP, Postl WF, Bolhár-Nordenkampf HR (1993) Quantum yields for uptake of carbon dioxide in C3 vascular plants of contrasting habitats and taxonomic groupings. Planta 189: 226-234. https://doi.org/10.1007/BF00195081
- Lobo F de A et al. (2013) Fitting net photosynthetic light-response curves with Microsoft Excel - a critical look at the models. Photosynthetica 51: 445-456. https://doi.org/10.1007/s11099-013-0045-y
- Elzhov TV, Mullen KM, Spiess A-N, Bolker B. minpack.lm: R interface to the Levenberg-Marquardt nonlinear least-squares algorithm found in MINPACK. R package. https://CRAN.R-project.org/package=minpack.lm
