# Demo data

The demo data are synthetic and are generated on the fly by `make_demo.R` (no image files are committed). `Rscript scripts/canopy_cover.R --demo` writes them to `<out>/demo_input/`; `Rscript demo/make_demo.R DIR` writes them to `DIR`.

- `rgb_ortho.tif`: an 8-bit RGB orthomosaic with 2 cm pixels (about 290 KB), 5.5 x 9.5 m, EPSG:32615. Nine plots (1.5 x 2.5 m) hold two crop rows of green disks whose size differs between plots, on brown soil, with a brightness gradient across the field, pixel noise and anti-aliased leaf edges (4 x 4 supersampling). Plot P07 has straw-colored residue patches between the rows, P03 a small saturated (glare) patch and P05 a dark shadow strip, so the QC outputs have something to show.
- `plots.gpkg`: the nine plot polygons with `plot_id`, `plot` and `true_cover_pct`, the exact vegetation area fraction of each plot in the synthetic scene, for comparison with the estimate.

The scene is drawn from a fixed seed (1), so every run produces the same files. Nothing comes from real imagery.
