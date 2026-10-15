"""Generate demo/plot_yields.csv: 12 synthetic soybean plots (fixed seed).

Plots are harvested as 2 rows x 30 in x 17.5 ft (the CLI flags in --demo); weights are
in lb at harvest moisture. Run from the skill folder: python demo/make_demo.py
"""
import csv
import random
from pathlib import Path

random.seed(20251015)
out = Path(__file__).resolve().parent / "plot_yields.csv"
entries = ["LD11", "Williams82", "NIL-n", "NIL-b"]
with out.open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["plot", "entry", "rep", "plot_weight", "harvest_moisture_pct"])
    plot = 101
    for rep in (1, 2, 3):
        for e in entries:
            # 8.3 lb from 87.5 ft2 of harvested area is about 69 bu/ac (60 lb/bu, before moisture adjustment)
            weight = random.gauss(8.3, 0.6)
            moisture = random.uniform(11.5, 16.5)
            w.writerow([plot, e, rep, f"{weight:.2f}", f"{moisture:.1f}"])
            plot += 1
print(f"wrote {out}")
