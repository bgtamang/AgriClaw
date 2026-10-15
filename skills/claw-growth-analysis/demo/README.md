# Demo data

`harvests.csv` is synthetic, written by `make_demo.R` (fixed seed 20261015). Re-create it with `Rscript demo/make_demo.R`.

Three genotypes (`A`, `B`, `C`), five destructive harvests (days 7, 14, 21, 28 and 35 after transplanting), six plants per harvest: 90 rows with `total_dry_mass_g`, `leaf_dry_mass_g` and `leaf_area_cm2`.

Total dry mass follows logistic growth, W(t) = K / (1 + ((K - W0) / W0) exp(-r t)) with W0 = 0.08 g, so the true RGR(t) = r (1 - W(t) / K):

| genotype | r (d-1) | K (g) |
|---|---|---|
| A | 0.20 | 40 |
| B | 0.17 | 40 |
| C | 0.22 | 25 |

Plant-to-plant variation is log-normal; the SD of ln mass rises from 0.15 at day 7 to 0.35 at day 35. Leaf mass fraction falls from about 0.55 to 0.41 and SLA from about 30 to 21 m2 kg-1 (genotype B higher, C lower).
