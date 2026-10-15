# Demo data

All four files are synthetic, written by `make_demo.R` (fixed seed 20261015). Re-create them with `Rscript demo/make_demo.R`.

- `weights.csv`: 12 sealed pots (`P01`-`P12`) logged every minute from 04:00 to 14:00 on one day; `weight_g` with 0.03 g balance noise.
- `plants.csv`: genotype, growth chamber (`A` or `B`) and leaf area (cm2) of each pot.
- `env.csv`: air temperature and RH every minute from two sensors in each chamber.
- `steps.csv`: one night window (04:10-05:00, dark, 20 degC) and seven 1-h daytime VPD steps (0.9 to 3.9 kPa at 30 degC); each window starts 10 min after the set-point change so plants have settled.

True responses (TR in mg H2O m-2 s-1): `G1` two segments with a breakpoint at 2.0 kPa (slopes 32, then 6); `G2` breakpoint at 2.6 kPa (slopes 30, then 10); `G3` linear (slope 28). Each pot's TR is scaled by a random factor (SD 6 %) to mimic leaf-area error, and responds to VPD changes with a 3-min lag.

Planted problems for the quality control to find: an irrigation of +150 g on `P07` at 10:20, a one-reading +25 g spike on `P02` at 08:41, and a 20-min logger dropout on `P10` from 12:30.

Real data from the published gravimetric platform (Tamang & Sadok 2018) are in https://github.com/bgtamang/Gravimetric-phenotyping-of-water-saving-traits; that repository has no licence file, so its data are not copied here.
