---
name: claw-ld-decay
description: Linkage disequilibrium (LD) decay analysis for crop SNP data. Computes r2 between marker pairs within chromosomes and fits the Hill and Weir (1988) decay curve by nonlinear least squares, as in Remington et al. (2001), alongside a model-free binned curve. Reports the distance at which r2 falls to 0.2 and to half its maximum, per chromosome and genome-wide, plus the background r2 of unlinked markers. Writes an LD window table that sets GWAS QTL-region and candidate-gene search windows. Use when the user asks how fast LD decays, for the LD decay or half-decay distance, the extent of LD in a panel or population, r2-versus-distance plots, how wide a GWAS candidate window should be, or what marker density association mapping needs.
license: MIT
compatibility: R >= 4.2 with the vcfR package (openssl optional, for SHA-256 checksums). No internet needed.
metadata:
  version: "1.0.0"
  domain: genomics
  language: R
  status: tested
  repository: https://github.com/bgtamang/AgriClaw
---

# LD decay and candidate-window distance

## What it does

For each chromosome, the skill samples marker pairs and computes r², the squared correlation of allele counts. It then fits the Hill-Weir expectation of r² against physical distance and computes a model-free binned curve. It reports the distance at which LD falls to r² = 0.2 and to half its maximum, plus the background r² between markers on different chromosomes. The per-chromosome distance is written to `tables/ld_window.csv`, which `claw-gwas-crop` (`--ld-table`) and `claw-gene-lookup` (`--ld-table`) read.

## When to use / when not to use

- Use before GWAS candidate-gene work, to set how far from a significant SNP to look, instead of a distance remembered from a paper on another panel.
- Use to compare LD extent between chromosomes, or between panels (run once per panel or subpopulation file).
- First run `claw-marker-qc` if the genotypes are not filtered (MAF, call rate, duplicates).
- Don't use for LD pruning, tagging SNP selection, haplotype blocks or phasing: use PLINK or Beagle-type tools.
- Structure-corrected LD (r² adjusted for kinship; Mangin et al. 2012) is not implemented. The skill reports the background r², so that structure-inflated LD is visible.

## Inputs

- `--geno`: the 0/1/2 marker table from `claw-marker-qc` (`marker, chr, pos`, then one column per line) or a VCF. Positions must be physical positions in base pairs. Genetic positions in cM give meaningless "bp" distances.

**Ask the user:**

- Which genome assembly the positions use, to record with the window.
- Whether the panel contains distinct subpopulations (for example, maturity groups or breeding programmes). Pooled panels show inflated LD; consider one run per subpopulation.
- Inbred or outbred material. For inbred lines, r² is the haplotype r²; for unphased heterozygous genotypes, it is the genotypic (composite) r² (Rogers & Huff 2009).

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
Rscript "${CLAUDE_SKILL_DIR}/scripts/ld_decay.R" --demo
Rscript "${CLAUDE_SKILL_DIR}/scripts/ld_decay.R" --geno qc/tables/genotypes_filtered_012.csv --out ld/
# short-range LD (e.g. outcrossing maize): restrict pair distances to resolve the decay
Rscript "${CLAUDE_SKILL_DIR}/scripts/ld_decay.R" --geno geno.csv --max-dist 1000000 --max-pairs 50000 --out ld_1Mb/
```

The demo has 150 doubled-haploid lines from a forward-simulated population with three chromosomes. Their recombination rates are known: 0.5, 1 and 2 cM/Mb. The reported decay distance shrinks from chromosome 1 to chromosome 3, as it should.

Other options:
- `--min-maf` (default 0.05)
- `--r2-threshold` (default 0.2)
- `--window-basis hw|binned`
- `--n-bins`
- `--n-unlinked`
- `--seed`

## Outputs

- `report.md`: warnings, the genome-wide decay distances, and a per-chromosome table:
  - r² = 0.2 distance by Hill-Weir fit and by binned means;
  - half-maximum distance, binned and Hill-Weir;
  - background r².
- `tables/ld_decay_summary.csv`: one row per chromosome plus `all` (pooled). Columns:
  - `n_markers`, `n_pairs`, `n_samples`, `max_pair_dist_bp`;
  - `rho_per_bp` (fitted 4Nc per bp), `hw_fit_at_boundary`;
  - `dist_r2_threshold_hw_bp`, `r2_max_hw`, `dist_half_max_hw_bp`;
  - `dist_r2_threshold_binned_bp`, `r2_max_binned`, `dist_half_max_binned_bp`;
  - `background_r2_p95`, `dist_background_hw_bp`.
- `tables/ld_window.csv`: `chr, ld_distance_bp, basis`, with an `all` row. This is the window file for the GWAS and gene-lookup skills.
- `tables/ld_binned.csv`: binned mean r² by distance, per chromosome and pooled.
- `figures/ld_decay_genome.png` and `figures/ld_decay_by_chr.png`. Each shows the binned means (orange), the Hill-Weir fit (blue), the r² threshold, the background r² and the decay distance.

## Expert decisions

- **r² rather than D′.** r² (Hill & Robertson 1968) measures how well one marker predicts another, which is what limits association-mapping resolution. Flint-Garcia et al. (2003) review LD measures in plants.
- **Pair sampling.** Up to 20,000 pairs per chromosome are sampled uniformly at random from all pairs (or all pairs closer than `--max-dist`), with a recorded seed. This is an adjustable rule of thumb for speed. When there are fewer pairs, all of them are used.
- **MAF ≥ 0.05** (adjustable rule of thumb): the possible range of r² shrinks when allele frequencies differ, so rare markers depress r².
- **Hill-Weir fit.** E(r²) = [(10 + C)/((2 + C)(11 + C))] · [1 + (3 + C)(12 + 12C + C²)/(n(2 + C)(11 + C))], with C = ρ·d (Hill & Weir 1988, as written by Remington et al. 2001).
  - n is the median number of samples per pair.
  - ρ is fitted by least squares on individual pairs: a grid over log10 ρ, then `optimize()`.
  - The curve's maximum is fixed by the model at about 0.45 + 0.74/n, so its half-maximum distance is close to its r² ≈ 0.23 distance. The **binned half-maximum** is the data-driven one.
- **Thresholds.** Adjustable conventions rather than biological constants:
  - r² = 0.2 is a common fixed threshold;
  - half of the maximum is a relative one;
  - the background level is the parametric 95th percentile of √r² for markers on different chromosomes, squared (Breseghello & Sorrells 2006).
- **Window.** `ld_window.csv` uses the Hill-Weir distance at the r² threshold, but only when all of these hold:
  - the fit did not hit the edge of its search range;
  - the distance is not an extrapolation below the closest sampled pair;
  - the binned data actually reach the threshold. The model's maximum is fixed near 0.45, so the curve "crosses" 0.2 even when no pair group does.

  Otherwise the window falls back to the binned value. If neither estimate is usable, the window is **NA**, with the reason in the `basis` column. `claw-gwas-crop` and `claw-gene-lookup` then use the `all` row, or stop and ask for an explicit window. `--window-basis binned` prefers binned values throughout.
- **Validation in the tests.**
  - Chromosomes with 0.5, 1 and 2 cM/Mb (3 founders, 20 generations of random mating, 200 lines) gave decay distances in the correct order by both estimators in all 10 seeds tried.
  - The binned half-maximum distance was 0.5 to 1.3 times the closed-form value (1 − θ)^(2G) = ½ per chromosome, and 0.7 to 1.0 times averaged over chromosomes.
  - Drift in a finite population makes single-chromosome estimates noisy.

## Checks before reporting

- Read the Warnings in `report.md` first.
- **Stop and tell the user** when:
  - Hill-Weir and binned distances differ more than 2-fold: the model fits poorly; report both and prefer the binned value;
  - a window is NA, or LD between linked markers is no higher than between unlinked ones: there is no measurable decay at this marker spacing. Never substitute a remembered value; ask the user, or try a lower `--r2-threshold` and say so;
  - the background r² is above about 0.1 (rule of thumb): structure or relatedness inflates LD, so decay distances are overestimated; suggest one run per subpopulation;
  - the decay distance exceeds half the longest pair distance (decay not resolved), or few pairs lie below it (short-range decay poorly sampled): re-run with a different `--max-dist`, as the warning says;
  - a chromosome has fewer than 50 markers, or there are fewer than 50 samples (rules of thumb).
- Check plausibility. Mating system strongly shapes LD (Flint-Garcia et al. 2003), and LD is predicted to be more extensive in selfing species (Nordborg 2000). A self-pollinated crop panel (soybean, rice, wheat) that decays within a few kb, or an outcrossing maize diversity panel with Mb-scale LD, suggests a map, unit or structure problem: ask the user before reporting.
- Decay varies along chromosomes (pericentromeric regions decay slowly), so a per-chromosome value is an average. Mention this when giving candidate windows.

## References

- Breseghello F, Sorrells ME (2006) Association mapping of kernel size and milling quality in wheat (*Triticum aestivum* L.) cultivars. Genetics 172:1165-1177. https://doi.org/10.1534/genetics.105.044586
- Flint-Garcia SA, Thornsberry JM, Buckler ES (2003) Structure of linkage disequilibrium in plants. Annual Review of Plant Biology 54:357-374. https://doi.org/10.1146/annurev.arplant.54.031902.134907
- Hill WG, Robertson A (1968) Linkage disequilibrium in finite populations. Theoretical and Applied Genetics 38:226-231. https://doi.org/10.1007/BF01245622
- Hill WG, Weir BS (1988) Variances and covariances of squared linkage disequilibria in finite populations. Theoretical Population Biology 33:54-78. https://doi.org/10.1016/0040-5809(88)90004-4
- Mangin B, Siberchicot A, Nicolas S, Doligez A, This P, Cierco-Ayrolles C (2012) Novel measures of linkage disequilibrium that correct the bias due to population structure and relatedness. Heredity 108:285-291. https://doi.org/10.1038/hdy.2011.73
- Nordborg M (2000) Linkage disequilibrium, gene trees and selfing: an ancestral recombination graph with partial self-fertilization. Genetics 154:923-929. https://doi.org/10.1093/genetics/154.2.923
- Remington DL, Thornsberry JM, Matsuoka Y, Wilson LM, Whitt SR, Doebley J, Kresovich S, Goodman MM, Buckler ES (2001) Structure of linkage disequilibrium and phenotypic associations in the maize genome. PNAS 98:11479-11484. https://doi.org/10.1073/pnas.201394398
- Rogers AR, Huff C (2009) Linkage disequilibrium between loci with unknown phase. Genetics 182:839-844. https://doi.org/10.1534/genetics.108.093153
