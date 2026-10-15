---
name: claw-selection-index
description: Predict response to selection and build multi-trait selection indices for plant breeding. Computes selection intensity from the proportion selected (with an exact small-population correction), the breeder's equation R = i h2 sigma_P and genetic gain per year i r sigma_A / L, the Smith-Hazel optimum index from phenotypic (P) and genetic (G) covariance matrices and economic weights, the Pesek-Baker desired-gains index, and the Kempthorne-Nordskog restricted index, with expected gain per trait and index accuracy r_IH. Checks that G is positive semi-definite and offers bending. Use when the user asks about selection intensity, expected genetic gain, response to selection, genetic gain per year or cycle, the breeder's equation, a selection index, economic weights, index coefficients, desired gains, restricting a trait (for example holding plant height or maturity constant), or comparing breeding strategies.
license: MIT
compatibility: R >= 4.2, base R only (openssl optional, for SHA-256 checksums). No internet needed.
metadata:
  version: "1.0.0"
  domain: genomics
  language: R
  status: tested
  summary: "Response to selection, selection intensity, Smith-Hazel and desired-gains indices"
  repository: https://github.com/bgtamang/AgriClaw
---

# Selection response and selection indices

## What it does

Turns variance components into breeding decisions. For one trait it gives the selection intensity and the expected response per cycle and per year (breeder's equation). For several traits it builds the Smith-Hazel index (b = P^-1 G w), the Pesek-Baker desired-gains index and the Kempthorne-Nordskog restricted index, and reports index coefficients, expected genetic gain per trait, the accuracy r_IH and the response in aggregate economic merit.

## When to use / when not to use

- Use to compare breeding strategies (intensity, accuracy, cycle length), to weight traits for selection, or to hold a trait constant while improving others.
- Variance components, heritabilities, P and G come from other skills: single trials from `claw-trial-analysis`, multi-environment trials from `claw-met-gxe`. This skill does not estimate them.
- For the accuracy r of genomic selection, run `claw-genomic-prediction` first and pass its accuracy here.
- Don't use for selection on unreplicated single plants with unknown heritability: ask for estimates first.

## Inputs

**Selection index**: two square CSV files with the same traits in the same order:

```
trait,yield_kg_ha,protein_pct,oil_pct
yield_kg_ha,90000,-90,36
protein_pct,-90,1.44,-0.432
oil_pct,36,-0.432,0.64
```

- `--p-matrix`: phenotypic (co)variances of the unit you select on, in trait units squared.
- `--g-matrix`: genetic (co)variances (additive, or among lines for inbred-line selection) on the same basis.
- Or, from a multi-environment analysis, `--e-matrix` (residual, plot level) with `--ge-matrix` (genotype x environment), `--n-env` and `--n-rep`: P for entry means is built as G + GE/e + E/(re) and saved to `tables/P_entry_means.csv`.
- `--weights`: economic weight per trait unit, `trait=value,...` or a CSV with `trait,weight`. Every trait needs a value; use 0 for a trait that has no economic value but is in the index.
- `--desired`: desired gain per cycle in trait units (`trait=value,...` or CSV `trait,desired_gain`).
- `--restrict`: traits whose genetic change is held at zero (needs `--weights`).
- `--proportion` (default 0.05) or `--intensity`; `--n-candidates` for small populations; `--cycle-years` for gain per year.

**Breeder's equation**: `--scenarios` CSV, one row per strategy. Columns: `scenario`, `proportion_selected`, then `h2` + `sigma_p` (phenotypic selection), or variance components `var_g`, `var_e`, `n_rep` (+ `var_ge`, `n_env` for several environments) from which entry-mean h2 and sigma_P are computed, and/or `accuracy` + `sigma_a`; optional `n_candidates`, `cycle_years`. Example: `demo/scenarios.csv`.

**Ask the user** (never guess):

- Do P and G come from the **same population and the same analysis**, on the **same basis** (plot, entry mean over r reps and e environments)? Heritability and P change with the number of replicates; G does not.
- What are the **units** of each trait and of each economic weight (value per kg/ha, per %)?
- What **proportion** will actually be selected, from how many candidates, and how many **years per cycle**?
- For gain per year: is the accuracy `r` from phenotypic selection (r = sqrt(h2) on the selection basis) or from genomic prediction (from `claw-genomic-prediction`)?

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/selection_index.R" --demo
Rscript "${CLAUDE_SKILL_DIR}/scripts/selection_index.R" --p-matrix P.csv --g-matrix G.csv \
    --weights "yield_kg_ha=1,protein_pct=60,oil_pct=40" --proportion 0.05 --out results/
Rscript "${CLAUDE_SKILL_DIR}/scripts/selection_index.R" --p-matrix P.csv --g-matrix G.csv \
    --weights w.csv --restrict height_cm --desired "yield_kg_ha=60,protein_pct=0.3,oil_pct=0,height_cm=0"
Rscript "${CLAUDE_SKILL_DIR}/scripts/selection_index.R" --scenarios scenarios.csv --out strategies/
Rscript "${CLAUDE_SKILL_DIR}/scripts/selection_index.R" --g-matrix G.csv --ge-matrix GE.csv --e-matrix E.csv \
    --n-env 4 --n-rep 2 --weights w.csv --proportion 0.1
```

## Outputs

- `report.md`: warnings first, then scenarios, input checks and index results.
- `tables/breeders_equation.csv`: per scenario `intensity`, `n_selected`, `response_per_cycle_phenotypic` (i h2 sigma_P), `response_per_cycle_accuracy` (i r sigma_A), `gain_per_year`.
- `tables/index_summary.csv`: per index `sigma_I`, `r_IH`, `response_H` (gain in economic merit per cycle), `intensity`.
- `tables/index_gains.csv`: per index and trait: coefficient `b`, `expected_gain` (trait units per cycle), `gain_per_year`, `gain_in_sd_g`, `gain_in_sd_p`.
- `tables/input_checks.csv`: heritabilities, eigenvalues of P, G and P - G, canonical heritabilities; `tables/genetic_correlations.csv`.
- `tables/bending_changes.csv`, `tables/G_bent.csv`: only if `--bend` was needed; `tables/P_entry_means.csv`: only with `--e-matrix`.
- `figures/expected_gains.png` (gains in genetic SD units by index), `figures/breeders_equation.png`.
- `reproducibility/`: `commands.sh`, `sessionInfo.txt`, checksums.

Read `b` only within an index: coefficients depend on trait units. Compare indices with `r_IH`, `response_H` and the gains per trait.

## Expert decisions

- **Selection intensity** for truncation selection is i = phi(z)/p, z the standard normal truncation point (Falconer & Mackay 1996). Tested values: p = 0.10, 0.05, 0.01, 0.005, 0.001 give 1.755, 2.063, 2.665, 2.892, 3.367 (Cobb et al. 2019, Table 1).
- **Small populations** (`--n-candidates`, or `n_candidates` in scenarios): when the best k = round(pN) of N are kept, the expected intensity is lower than the table value (Burrows 1972). The script computes the exact expected mean of the top k of N normal values by numerical integration, e.g. 1.674 for the best 3 of 30 versus 1.755 for p = 0.10.
- **Breeder's equation**: R = i h2 sigma_P per cycle (Falconer & Mackay 1996) and Delta G = i r sigma_A / L per year (Cobb et al. 2019). With r = h and sigma_A = h sigma_P the two agree. h2, sigma_P and r must refer to the unit actually selected (e.g. entry means).
- **Entry means from variance components**: var_P = var_G + var_GE/e + var_e/(r e) and h2 = var_G / var_P (Holland et al. 2002); the same formula on covariance matrices gives the entry-mean P. G x E is never assumed absent: with more than one environment the GE matrix is required.
- **Smith-Hazel index**: b = P^-1 G w maximises the correlation r_IH between the index and the aggregate genotype H = w'g (Smith 1936; Hazel 1943). sigma_I = sqrt(b'Pb); expected gain per trait = i G b / sigma_I; r_IH = sigma_I / sigma_H; response in H = i sigma_I. The base index (b = w) is reported for comparison.
- **Desired gains**: b = P^-1 G (G P^-1 G)^-1 d (Itoh & Yamada 1986), equal to the Pesek & Baker (1969) b = G^-1 d when G is invertible. Gains are proportional to d; the report says how many cycles reach d. It ignores economic weights and does not maximise r_IH.
- **Restricted index**: b_R = [I - P^-1 C (C'P^-1 C)^-1 C'] P^-1 G w with C the G columns of restricted traits (Kempthorne & Nordskog 1959); their expected change is exactly zero, at a cost in r_IH.
- **Checks on P and G**: P must be positive definite (it is inverted). G must be positive semi-definite; if not, the run stops unless `--bend` is given. Canonical heritabilities (eigenvalues of P^-1 G) must lie in [0, 1]; a value above 1 means P - G is not a covariance matrix, so P and G are inconsistent.
- **Bending** (`--bend`, opt-in): eigenvalues of G below `--bend-eps` (default 1e-4, as in Jorjani et al. 2003 and the mbend package; Nilforooshan 2020) are set to that value and G is rebuilt. Bending was introduced for selection indices by Hayes & Hill (1981). The changes are written to `bending_changes.csv`; tell the user G was modified.
- **Default proportion selected 0.05**: an adjustable rule of thumb for the demo; always use the user's real proportion.
- **|r_G| > 0.95 warning**: adjustable rule of thumb; near-collinear traits make b unstable.
- Gains smaller than 1e-10 times the largest gain are printed as 0 (floating-point noise).

**What the tests validate**: intensities against Cobb et al. (2019, Table 1) to 3 decimals; exact small-population values (best 1 of 2 = 1/sqrt(pi), best 1 of 3 = 3/(2 sqrt(pi))); a hand-computed two-trait Smith-Hazel index; the published worked examples of Ceron-Rojas & Crossa (2018): Smith-Hazel b = (0.555, -1.063, 1.087), sigma_I = 9.312, response 19.21, gains (9.51, -5.48, 4.22), r_IH = 0.894, base index 0.875 (section 2.8.2), and the restricted index b_R = (-0.35, -0.41, 0.59, 0.89), response 6.87, gains (0, -2.2, 2.03, 2.66), r_IH = 0.35 (section 3.1.4); Smith-Hazel beats 5000 random indices; desired gains are reached in the requested proportions; bending reproduces Nilforooshan (2020, Table 1) to 4 decimals; entry-mean P and h2 from variance components by hand.

## Checks before reporting

- Read the Warnings section of `report.md` first and pass every warning on to the user.
- **Stop and tell the user** if: G was not positive semi-definite (do not bend without their agreement); a canonical heritability exceeds 1 or a heritability exceeds 1; P and G come from different populations, years or analyses; trait units of P, G and weights disagree.
- Check `input_checks.csv`: heritabilities and genetic correlations should be plausible for the crop and trait. Implausible values (h2 near 0 or 1, |r_G| near 1) usually mean estimation problems.
- A Smith-Hazel index with a large weight on a trait of zero economic value is normal (it is used as an indirect predictor); explain it rather than "fixing" it.
- Expected gains assume P and G are known; with small trials the realised gain is lower. Say so when reporting gains.
- If the desired-gains report says many cycles are needed, the desired gains are ambitious for the stated intensity; do not present them as achievable in one cycle.

## References

- Falconer DS, Mackay TFC (1996) Introduction to Quantitative Genetics, 4th edn. Longman, Harlow. ISBN 978-0-582-24302-6.
- Cobb JN, Juma RU, Biswas PS, Arbelaez JD, Rutkoski J, Atlin G, Hagen T, Quinn M, Ng EH (2019) Enhancing the rate of genetic gain in public-sector plant breeding programs: lessons from the breeder's equation. Theor Appl Genet 132:627-645. https://doi.org/10.1007/s00122-019-03317-0
- Burrows PM (1972) Expected selection differentials for directional selection. Biometrics 28:1091-1100. https://doi.org/10.2307/2528642
- Smith HF (1936) A discriminant function for plant selection. Ann Eugen 7:240-250. https://doi.org/10.1111/j.1469-1809.1936.tb02143.x
- Hazel LN (1943) The genetic basis for constructing selection indexes. Genetics 28:476-490. https://doi.org/10.1093/genetics/28.6.476
- Pesek J, Baker RJ (1969) Desired improvement in relation to selection indices. Can J Plant Sci 49:803-804. https://doi.org/10.4141/cjps69-137
- Itoh Y, Yamada Y (1986) Re-examination of selection index for desired gains. Genet Sel Evol 18:499-504. https://doi.org/10.1186/1297-9686-18-4-499
- Kempthorne O, Nordskog AW (1959) Restricted selection indices. Biometrics 15:10-19. https://doi.org/10.2307/2527598
- Holland JB, Nyquist WE, Cervantes-Martinez CT (2002) Estimating and interpreting heritability for plant breeding: an update. Plant Breed Rev 22:9-112. https://doi.org/10.1002/9780470650202.ch2
- Ceron-Rojas JJ, Crossa J (2018) Linear Selection Indices in Modern Plant Breeding. Springer, Cham. https://doi.org/10.1007/978-3-319-91223-3 (worked examples in sections 2.8.2 and 3.1.4 are reproduced by the tests)
- Hayes JF, Hill WG (1981) Modification of estimates of parameters in the construction of genetic selection indices ('bending'). Biometrics 37:483. https://doi.org/10.2307/2530561
- Jorjani H, Klei L, Emanuelson U (2003) A simple method for weighted bending of genetic (co)variance matrices. J Dairy Sci 86:677-679. https://doi.org/10.3168/jds.S0022-0302(03)73646-7
- Nilforooshan MA (2020) mbend: an R package for bending non-positive-definite symmetric matrices to positive-definite. BMC Genetics 21:97. https://doi.org/10.1186/s12863-020-00881-z
