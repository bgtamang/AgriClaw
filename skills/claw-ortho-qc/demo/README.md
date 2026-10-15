# Demo data

The demo data are synthetic and are generated on the fly by `make_demo.R` (no image files are committed). `Rscript scripts/ortho_qc.R --demo` writes them to `<out>/demo_input/`; `Rscript demo/make_demo.R DIR` writes them to `DIR`. Together they are about 200 KB.

- `ortho_reflectance.tif`: 4-band surface reflectance (blue, green, red, NIR) stored as reflectance x 10000 in 16-bit integers with NoData = 0; 5 cm pixels, 5.5 x 10.5 m, EPSG:32615. The demo runs with `--expect reflectance --scale 0.0001 --expect-bands 4`.
- `plots.gpkg`: 12 plots of 1.5 x 2 m in 3 columns and 4 ranges, each with two crop rows.
- `dsm.tif`: surface model (float32, 5 cm): ground near 300 m on a gentle slope plus 0.2-0.9 m of canopy on the crop rows, with 1 cm noise.
- `dtm.tif`: bare-ground terrain model (float32, 10 cm).

Problems injected on purpose, which the demo report flags as warnings (and nothing else):

1. A NoData hole over 10% of plot P06 (`plot_coverage`).
2. Plot P12 runs past the edge of the image (`extent`, `plot_coverage`, `dsm_coverage`).
3. 0.3% of DSM pixels inside plots P01-P11 carry +3 m spikes (`dsm_spikes`).

The data are drawn from a fixed seed (7), so every run produces the same files. Nothing comes from real imagery.
