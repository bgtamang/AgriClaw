# Demo data

All files are synthetic, written by `make_demo.R` with a fixed seed (`Rscript demo/make_demo.R` regenerates them exactly).

Eight inbred founders with random SNP alleles were crossed in seven biparental families (F1-F7). Recombinant inbred lines were derived by single-seed descent to F6 with a forward simulation of meiosis (5 chromosomes of 100 cM, Haldane map, no crossover interference). Genetic values come from 60 additive QTL placed on randomly chosen SNPs (genetic SD 250 kg/ha).

| File | Content |
|---|---|
| `genotypes.csv` | 270 lines x 500 SNPs, coded 0/1/2 (count of the allele coded 1); a little residual heterozygosity remains at F6 |
| `phenotypes.csv` | `id`, `family`, `yield_kg_ha`: entry-mean BLUEs for families F1-F6 (240 lines), mean 3500 kg/ha, true heritability 0.5 |
| `phenotypes_multienv.csv` | the same 240 lines in three environments (E1-E3: means 3000, 3500, 4100 kg/ha), with G x E (genetic correlation 0.8 between environments) and error SD 300 kg/ha |
| `truth.csv` | true genetic value of every line, for checking predictions only (never an input) |

Family F7 (30 lines) has genotypes but no phenotypes: these are the selection candidates whose predictions appear in `predictions_all_lines.csv`.
