# Demo data

`fluorescence.csv` is synthetic, written by `make_demo.R` (fixed seed 20261015). Re-create it with `Rscript demo/make_demo.R`.

Three genotypes (`G1`-`G3`) x 4 plants. Each plant has one dark-adapted Fo and Fm (repeated on its rows) and light-adapted flashes at PPFD 100, 300, 600, 1000 and 1500 umol m-2 s-1: 60 rows with `fo`, `fm`, `fs`, `fm_prime`, `fo_prime` and `ppfd_umol_m2_s`.

Signals are built from known parameters: Fm' = Fm / (1 + NPQ), Fo' = Fo / (Fv/Fm + Fo/Fm') and Fs = Fm' - qP (Fm' - Fo'), with NPQ rising and qP falling with light (G3 is the most quenched and has the lowest Fv/Fm, about 0.79), plus 1 % multiplicative noise. `fo_prime` is "measured" for G1 and G2 and missing for G3, so G3's Fo' is estimated.

Planted problems for the quality control: plant `G1-4` has an underestimated Fo (Fv/Fm about 0.90), and plant `G2-3`'s Fm was recorded before full dark adaptation (70 % of the true value), so its Fv/Fm is about 0.75 and Fm' exceeds Fm at the lowest light.

LI-6800 example files are not copied here: the tests read the eight files that ship with the PhotoGEA package through `system.file("extdata", package = "PhotoGEA")`.
