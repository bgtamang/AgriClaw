# Demo data

All files are synthetic, written by `make_demo.R` with a fixed seed (`Rscript demo/make_demo.R` regenerates them exactly). The same script provides the population simulators used by the tests.

A population of 160 recombinant inbred lines (RIL by selfing) was simulated from two inbred parents, A and B, by single-seed descent to F8 with an explicit model of meiosis (5 chromosomes of 100 cM, Haldane map, no crossover interference). Residual heterozygous calls are written as missing, and 1% of calls are set missing at random.

| File | Content |
|---|---|
| `geno.csv` | 160 lines x 150 markers (about every 3.4 cM), codes `A`, `B`, `-` (missing) |
| `gmap.csv` | `marker,chr,pos`: genetic map in cM |
| `pmap.csv` | `marker,chr,pos`: invented physical positions in Mb (about 0.45 Mb per cM, varying) |
| `pheno.csv` | `id`, `days_to_flowering` (3 values missing), `plant_height_cm` |
| `truth.csv` | the simulated QTL: chromosome, position (cM) and additive effect (half the BB - AA difference) |

Simulated QTL: `days_to_flowering` on chr 1 at 45.3 cM (a = +2.0 days) and chr 3 at 70.8 cM (a = -1.4 days), residual SD 2.2 days; `plant_height_cm` on chr 2 at 21.6 cM (a = +3.0 cm), residual SD 8 cm.

The tests also use the `iron` example data distributed with the qtl2 package (`system.file("extdata", "iron.zip", package = "qtl2")`; Grant et al. 2006, Hepatology 44:174-185), loaded from the installed package rather than copied.
