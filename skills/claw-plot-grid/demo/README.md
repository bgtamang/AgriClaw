# Demo data

Both files are synthetic and were written by `make_demo.R` (seed 2026; re-run with `Rscript demo/make_demo.R` to regenerate them identically).

- `corners.csv`: the four outer corners (lon, lat, WGS84) of an invented block of 6 ranges x 8 columns of 2-row plots (5.2 m long, 2 x 0.76 m wide, 0.9 m alleys) planted at a true bearing of 12 degrees with columns to the right. Each corner carries about 1 cm of simulated RTK noise, so the corner check in the demo report shows small, realistic misfits. The location (near 41.30 N, 95.80 W, UTM zone 15N) is arbitrary farmland; no real field is represented.
- `fieldbook.csv`: an invented RCBD field book with 12 entries x 4 reps = 48 plots, numbered serpentine along planter passes from corner c1. `row` and `col` give each plot's position (range and pass counted from c1); the demo uses them with `--fb-coords row,col` to cross-check the numbering.
