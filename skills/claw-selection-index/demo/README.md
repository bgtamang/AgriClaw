# Demo data

All files are synthetic and illustrative. They are written by `make_demo.R` (no random numbers; rerun with `Rscript demo/make_demo.R`). They are not estimates from any real population.

| File | Content |
|---|---|
| `P_matrix.csv` | phenotypic covariance matrix of four soybean-like traits (`yield_kg_ha`, `protein_pct`, `oil_pct`, `height_cm`), built from invented phenotypic SDs (300, 1.2, 0.8, 8) and phenotypic correlations |
| `G_matrix.csv` | genetic covariance matrix from invented heritabilities (0.45, 0.70, 0.65, 0.60) and genetic correlations (for example yield-protein -0.40, protein-oil -0.60) |
| `weights.csv` | invented relative economic weights per trait unit (1, 60, 40, 0) |
| `desired_gains.csv` | desired gains per cycle (60 kg/ha yield, +0.3 % protein, no change in oil or height) |
| `scenarios.csv` | three breeder's-equation scenarios: phenotypic selection on a 4-year cycle (infinite and 30-candidate populations) and genomic selection on a 2-year cycle |

The generator checks that P is positive definite and that G and P - G are positive semi-definite, so every implied heritability lies between 0 and 1.

The tests also use two published worked examples, typed in from Ceron-Rojas & Crossa (2018, Linear Selection Indices in Modern Plant Breeding, Springer, open access, https://doi.org/10.1007/978-3-319-91223-3, sections 2.8.2 and 3.1.4), and the 5 x 5 example matrix from Nilforooshan (2020, BMC Genetics 21:97, https://doi.org/10.1186/s12863-020-00881-z).
