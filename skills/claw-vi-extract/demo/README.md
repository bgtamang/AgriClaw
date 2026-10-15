# Demo data

The demo data are synthetic and are generated on the fly by `make_demo.R` (no image files are committed). `Rscript scripts/vi_extract.R --demo` writes them to `<out>/demo_input/`; `Rscript demo/make_demo.R DIR` writes them to `DIR`.

- `ortho_multispectral.tif`: a 5-band surface-reflectance image stored as 16-bit integers (reflectance x 10000, about 310 KB), so the demo also exercises `--scale 0.0001`. Layers are named blue, green, red, rededge and nir, in wavelength order as in a stacked MicaSense export. 5 cm pixels, EPSG:32615, 6.5 x 16.25 m.
- `plots.gpkg`: 12 plots (3 columns x 4 ranges of 1.5 x 3 m with 0.75 m alleys) with `plot_id`, `plot`, `rep` and `genotype`.

Each plot has two crop rows whose canopy width differs between plots (25-95% cover) and its own canopy spectrum (darker or paler green); soil fills the gaps between rows, the alleys and the margins. Every band carries 3% multiplicative noise. The demo runs with an Otsu soil mask and a -0.25 m buffer, so the report shows the mask threshold and the vegetation fraction per plot.

The scene is drawn from a fixed seed (42), so every run produces the same files. Nothing comes from real imagery.
