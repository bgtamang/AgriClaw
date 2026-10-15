# Demo data

`plot_yields.csv` is synthetic. `make_demo.py` generates it with a fixed seed (`python demo/make_demo.py`). It holds 12 soybean plots (4 entries × 3 replicates) with:

- `plot_weight`: harvested grain weight in lb at harvest moisture, about 8.3 ± 0.6 lb;
- `harvest_moisture_pct`: 11.5–16.5% wet basis.

The demo command treats each plot as 2 harvested rows × 30 in × 17.5 ft (87.5 ft²) and reports yields at the soybean standard of 13% moisture with a 60 lb bushel, which gives about 60–74 bu/ac. The values only exercise the skill; they are not data from a real trial.
