---
name: claw-growth-analysis
description: Classical and functional plant growth analysis from sequential destructive harvests. Computes relative growth rate (RGR) from the mean of ln dry mass with standard errors and confidence intervals (avoiding the ln-of-means bias), net assimilation rate (NAR, unit leaf rate), leaf area ratio (LAR), specific leaf area (SLA) and leaf mass fraction (LMF), plus a spline-based RGR(t) curve with confidence bands and RGR versus plant size. Use when the user has dry mass, leaf area and leaf mass from repeated harvests of plants (genotypes or treatments in pots, growth chambers or the field) and asks for growth rates, RGR, NAR, ULR, LAR, SLA, LMF/LWR, growth curves, Hunt or Poorter style growth analysis, or comparing growth between genotypes or treatments.
license: MIT
compatibility: R >= 4.2 with mgcv (ships with R); openssl optional for SHA-256 checksums. No internet needed.
metadata:
  version: "1.0.0"
  domain: physiology
  language: R
  status: tested
  summary: "Growth analysis from sequential harvests: RGR, NAR, LAR, SLA, LMF"
  repository: https://github.com/bgtamang/AgriClaw
---

# Plant growth analysis (RGR, NAR, LAR, SLA, LMF)

## What it does

Takes one row per harvested plant and returns, for each genotype or treatment: harvest means; classical RGR and NAR between consecutive harvests with confidence intervals; SLA, LAR and LMF at each harvest; the mean RGR over the experiment; and a functional RGR(t) curve from a penalized spline of ln mass on time. It uses the unbiased RGR estimator of Hoffmann & Poorter (2002) and treats harvests as unpaired samples of different plants.

## When to use / when not to use

- Use for destructive sequential harvests (each plant weighed once) with total dry mass, and optionally leaf area and leaf dry mass.
- Not for repeated non-destructive measurements of the same plants (imaging, scanner area): those are paired, so per-plant growth curves are better; the unpaired SEs here would be conservative. For image-based time series see `claw-vi-timeseries`.
- Not for canopy-level crop growth rates per unit ground area (g m-2 d-1 of field plots) or yield components.
- Not for whole-plant water use: use `claw-vpd-transpiration`.

## Inputs

CSV, one row per harvested plant:

| Column (option) | Meaning |
|---|---|
| `genotype` (`--group`; `''` for one group) | genotype or treatment |
| `days` (`--time-col`) | harvest time in days (after sowing, transplanting or treatment start), or ISO dates (converted to days after the first harvest) |
| `total_dry_mass_g` (`--mass-col`) | total plant dry mass, g, > 0 |
| `leaf_area_cm2` (`--leaf-area-col`) | optional, one-sided leaf area, cm2 (for NAR, LAR, SLA) |
| `leaf_dry_mass_g` (`--leaf-mass-col`) | optional, leaf dry mass, g (for SLA, LMF) |

Missing total mass is an error (remove those plants); missing leaf area or leaf mass is allowed and counted in the report.

**Ask the user before running:**

1. Is the mass whole-plant (shoot + root)? If shoot only, every result is a shoot growth rate: say so in the report.
2. Were the harvests destructive, with plants chosen at random at each harvest? Never pair plants across harvests.
3. What time zero means (sowing, transplanting, treatment start) and that times are in days.
4. Leaf area units (cm2), one-sided, and whether petioles are included in leaf area and leaf mass.
5. How many plants per harvest, and whether some plants died or were lost (do not fill them in).

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/growth_analysis.R" --demo --out demo_growth
Rscript "${CLAUDE_SKILL_DIR}/scripts/growth_analysis.R" --input harvests.csv --group genotype \
  --time-col days --mass-col total_dry_mass_g --out results/
Rscript "${CLAUDE_SKILL_DIR}/scripts/growth_analysis.R" --input shoots.csv --group treatment \
  --mass-col shoot_dry_mass_g --leaf-mass-col "" --out results/   # no leaf mass column
```

Options: `--k-max` (largest spline basis, default 10), `--n-boot` (NAR bootstrap resamples, default 2000), `--min-n` (warn below this many plants per harvest, default 4), `--seed`.

## Outputs

- `report.md`: warnings first, then harvest summary, classical intervals, mean RGR, functional RGR, figures, limitations.
- `tables/harvest_summary.csv`: per group and harvest: `n`, `mass_mean_g`, `ln_mass_mean`, `ln_mass_sd`, `mass_geomean_g`, `leaf_area_mean_cm2`, `sla_m2_kg`, `lar_m2_kg`, `lmf_g_g` (each ratio with `_se`).
- `tables/classical_intervals.csv`: per consecutive harvest pair: `rgr_g_g_d` with `rgr_se`, `rgr_df`, `rgr_lo`, `rgr_hi`; `rgr_ln_of_means_g_g_d` (the biased estimator, only for comparison; never report it as the result); `nar_g_m2_d` with bootstrap `nar_se`, `nar_lo`, `nar_hi`.
- `tables/overall_rgr.csv`: slope of ln mass on time per group with CI, `r2`, and `p_curvature` (quadratic term; p < 0.05 means RGR changed over time).
- `tables/functional_rgr.csv`: RGR(t) on a time grid and at each harvest (`at_harvest`): `rgr_g_g_d`, `rgr_lo`, `rgr_hi`, fitted `ln_mass_fit` and `mass_fit_g`.
- `figures/ln_mass.png`, `figures/rgr.png` (functional curve with band; dashed segments and bars = classical interval RGR and CI), `figures/rgr_vs_size.png`, `figures/allocation.png`.
- `reproducibility/`: `commands.sh`, `sessionInfo.txt`, checksums.

Units: RGR g g-1 d-1 (multiply by 1000 for mg g-1 d-1); NAR g m-2 d-1; SLA and LAR m2 kg-1 (= cm2 g-1 / 10); LMF g g-1.

## Expert decisions

- **RGR** between harvests = (mean ln W2 - mean ln W1) / (t2 - t1). Using ln of the mean mass instead is biased whenever the variance of ln W changes between harvests (Hoffmann & Poorter 2002); that version is printed only to show the size of the bias.
- **Unpaired SE**: SE = sqrt(s1^2/n1 + s2^2/n2) / (t2 - t1) from the two harvests' ln-mass variances, with a Welch-Satterthwaite t interval (variances may differ between harvests).
- **NAR** (unit leaf rate) = (W2 - W1)/(t2 - t1) x (ln LA2 - ln LA1)/(LA2 - LA1) on harvest means (Williams 1946; Radford 1967). This is exact only if mass is linear in leaf area within the interval (test 2); otherwise it is an approximation. Its CI is a percentile bootstrap that resamples plants within each harvest independently (2000 resamples, adjustable).
- **SLA, LAR, LMF**: per-plant ratios averaged over the plants of a harvest (SE across plants), so per plant LAR = SLA x LMF exactly; RGR = NAR x LAR holds only instantaneously, not for interval means.
- **Mean RGR over the experiment**: OLS slope of ln W on time using every plant (harvests are independent samples). A significant quadratic term (`p_curvature` < 0.05) means RGR drifted, so the slope is only an average.
- **Functional RGR(t)**: penalized cubic regression spline `s(time, bs = "cr")` of ln W on time in mgcv, smoothing chosen by REML (Wood 2011, 2017), basis size k = min(`--k-max`, number of harvest times). The penalty shrinks toward a straight line, i.e. constant RGR (exponential growth). RGR(t) is the derivative of the fitted curve; its 95 % band uses the covariance corrected for smoothing-parameter uncertainty. These are Bayesian intervals with good average, not exact pointwise, coverage (Marra & Wood 2012). This is the "form-free" functional approach (Shipley & Hunt 1996; Poorter 1989; Hunt 1990) rather than a fitted parametric growth function (Paine et al. 2012). With k capped at 5, simulated logistic growth over 9 harvests had pointwise coverage down to 72 %, so the cap is 10.
- **Thresholds that are rules of thumb**: fewer than 4 plants per harvest (`--min-n`) gives unreliable SEs; RGR above 0.5 g g-1 d-1 for whole plants suggests wrong time units or labels; groups differing by > 30 % in mass at the first harvest trigger the size-dependence warning.

## Checks before reporting

Stop and tell the user when:

- any warning appears at the top of `report.md` (pass every one on);
- an interval has RGR <= 0, or RGR > 0.5 g g-1 d-1: check harvest labels, time units, senescence or lost plants;
- groups differ in size at the start: RGR falls as plants grow (ontogenetic drift), so faster growth of a smaller group may be a size effect. Compare RGR at a common mass (`figures/rgr_vs_size.png`; Rees et al. 2010) and interpret allocation (LMF) against size, not only time (Poorter & Sack 2012);
- `p_curvature` < 0.05: do not report the overall slope as "the" RGR; report interval or functional RGR;
- NAR changes sign or its CI is very wide: leaf area and mass may not be linearly related, or harvests are too far apart;
- leaf area or mass look like the wrong units (rule of thumb: SLA outside about 5-60 m2 kg-1 for crop leaves is suspicious; ask).

Always state that harvests are unpaired, the mass basis (whole plant or shoot), and the units.

## Validation

`tests/test_growth_analysis.R` (40 checks): hand-calculated RGR (ln 4 / 10), its SE, NAR (17.655 g m-2 d-1), SLA, LAR and LMF on a two-plant example; NAR equal to the exact interval mean (4 g m-2 d-1) when mass is linear in leaf area; simulated exponential growth (RGR 0.12, 5 harvests x 8 plants, 80 data sets): classical mean 0.1198 with 95 % CI coverage 95 %, overall slope unbiased, functional band coverage 100 % (25 data sets); simulated logistic growth (9 harvests x 6 plants, 20 data sets): functional RGR(t) mean absolute error 0.011 g g-1 d-1 with 98 % band coverage, classical intervals unbiased for the exact interval means; Hoffmann & Poorter bias (SD of ln W 0.1 then 0.6, 600 data sets): the skill's estimator averaged 0.0989 (truth 0.1) while ln of mean mass overestimated by 0.0213 (delta-method prediction 0.0220); bad inputs; and the demo run. Not validated: the NAR bootstrap interval's coverage, and data with plant mortality or very unequal harvest sizes.

## References

- Hoffmann WA, Poorter H (2002) Avoiding bias in calculations of relative growth rate. Annals of Botany 90:37-42. https://doi.org/10.1093/aob/mcf140
- Hunt R (1990) Basic Growth Analysis. Springer Netherlands, Dordrecht. https://doi.org/10.1007/978-94-010-9117-6
- Hunt R, Causton DR, Shipley B, Askew AP (2002) A modern tool for classical plant growth analysis. Annals of Botany 90:485-488. https://doi.org/10.1093/aob/mcf214
- Marra G, Wood SN (2012) Coverage properties of confidence intervals for generalized additive model components. Scandinavian Journal of Statistics 39:53-74. https://doi.org/10.1111/j.1467-9469.2011.00760.x
- Paine CET, Marthews TR, Vogt DR, Purves D, Rees M, Hector A, Turnbull LA (2012) How to fit nonlinear plant growth models and calculate growth rates: an update for ecologists. Methods in Ecology and Evolution 3:245-256. https://doi.org/10.1111/j.2041-210X.2011.00155.x
- Poorter H (1989) Plant growth analysis: towards a synthesis of the classical and the functional approach. Physiologia Plantarum 75:237-244. https://doi.org/10.1111/j.1399-3054.1989.tb06175.x
- Poorter H, Garnier E (1996) Plant growth analysis: an evaluation of experimental design and computational methods. Journal of Experimental Botany 47:1343-1351. https://doi.org/10.1093/jxb/47.9.1343
- Poorter H, Sack L (2012) Pitfalls and possibilities in the analysis of biomass allocation patterns in plants. Frontiers in Plant Science 3:259. https://doi.org/10.3389/fpls.2012.00259
- Radford PJ (1967) Growth analysis formulae: their use and abuse. Crop Science 7:171-175. https://doi.org/10.2135/cropsci1967.0011183X000700030001x
- Rees M, Osborne CP, Woodward FI, Hulme SP, Turnbull LA, Taylor SH (2010) Partitioning the components of relative growth rate: how important is plant size variation? The American Naturalist 176:E152-E161. https://doi.org/10.1086/657037
- Shipley B, Hunt R (1996) Regression smoothers for estimating parameters of growth analyses. Annals of Botany 78:569-576. https://doi.org/10.1006/anbo.1996.0162
- Williams RF (1946) The physiology of plant growth with special reference to the concept of net assimilation rate. Annals of Botany 10:41-72. https://doi.org/10.1093/oxfordjournals.aob.a083119
- Wood SN (2011) Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. Journal of the Royal Statistical Society Series B 73:3-36. https://doi.org/10.1111/j.1467-9868.2010.00749.x
- Wood SN (2017) Generalized Additive Models: An Introduction with R, 2nd edn. Chapman and Hall/CRC. https://doi.org/10.1201/9781315370279
