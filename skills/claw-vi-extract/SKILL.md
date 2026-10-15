---
name: claw-vi-extract
description: Extract per-plot vegetation indices (NDVI, NDRE, GNDVI, OSAVI, EVI2, MCARI2, CIred-edge) from multispectral drone/UAV orthomosaics (MicaSense RedEdge-MX or Altum, DJI P4 Multispectral) or satellite images (Sentinel-2 L2A, PlanetScope) with plot polygons. Gives coverage-weighted plot mean, median, p10, p90, SD and pixel counts (exactextractr), handles reflectance scale and offset (including the Sentinel-2 -1000 offset), refuses uncalibrated digital numbers, offers a soil/shadow mask (NDVI threshold or Otsu) and inward buffers, and warns about mixed pixels. Use when the user asks for NDVI or other vegetation indices per plot, plot-level reflectance, or zonal statistics of drone or satellite imagery over field-trial plots.
license: MIT
compatibility: R >= 4.2 with terra, sf and exactextractr (GDAL, GEOS, PROJ). No internet needed.
metadata:
  version: "1.0.0"
  domain: remote-sensing
  language: R
  status: tested
  repository: https://github.com/bgtamang/AgriClaw
---

# Per-plot vegetation indices

## What it does

Converts a multispectral raster to surface reflectance (scale, offset, no-data), computes vegetation indices per pixel, optionally removes soil and shadow pixels, and summarises each plot polygon with coverage-weighted statistics from exactextractr. It refuses imagery that still holds digital numbers and flags plots too small for the pixel size.

## When to use / when not to use

- Use after `claw-ortho-qc` (image checks) and `claw-plot-grid` (polygons), once per flight or satellite scene. Give `--date` and stack the per-date `tables/plot_indices.csv` files for `claw-vi-timeseries`; plot means also feed `claw-yield-prediction`.
- RGB-only imagery (no NIR): use `claw-canopy-cover` instead; NIR-based indices cannot be computed from RGB.
- Thermal bands and hyperspectral narrow-band indices are not covered: bands are treated as broad multispectral bands.
- Not a radiometric calibration tool: the input must already be surface reflectance (processing advice: `claw-ortho-qc` `references/processing-guide.md`).

## Inputs

- `--input`: a multi-band raster of **surface reflectance** (0-1, or scaled integers with a known scale), with a CRS, ideally projected in metres.
- `--plots`: plot polygons (GeoPackage, shapefile, GeoJSON) with a unique id column (`--id-col`, default `plot_id`, else `plot`). Any CRS; they are transformed to the raster's.
- `--bands`: which raster layer holds blue, green, red, red edge and NIR, by number or layer name.

Sensor presets (`--sensor`) suggest a band order and scaling. They are only a starting point: **band order in stacked orthomosaics depends on the processing software**.

| Preset | Suggested bands | Scaling |
|---|---|---|
| `micasense-rededge-mx`, `micasense-altum` | 1 blue, 2 green, 3 red, 4 red edge, 5 NIR (wavelength order, as Pix4D/Metashape write it; the camera itself numbers NIR 4 and red edge 5) | none (float reflectance) |
| `dji-p4m` | 1 blue, 2 green, 3 red, 4 red edge, 5 NIR | none |
| `sentinel2-l2a` | layers named B02, B03, B04, B05, B08 | (DN - 1000) / 10000, DN 0 = no data |
| `sentinel2-l2a-harmonized` | as above | DN / 10000 |
| `planetscope-4band` / `-8band` | B, G, R, NIR / 2 B, 4 G, 6 R, 7 RE, 8 NIR | DN / 10000, 0 = no data |

**Ask the user, never guess:**

1. How the orthomosaic was made (Pix4D, Metashape, OpenDroneMap, the MicaSense library) and the band order it wrote. Check layer names with `terra::rast(f)` or `gdalinfo`.
2. Whether it is reflectance: were reflectance panels and/or the sunlight sensor used? If values are integers, what is the scale factor?
3. For Sentinel-2: where the data came from. ESA or Copernicus products from processing baseline 04.00 (25 January 2022) carry the -1000 offset (`sentinel2-l2a`); Google Earth Engine's `S2_SR_HARMONIZED` and older baselines do not (`sentinel2-l2a-harmonized`). Applying it twice makes red negative.
4. The flight or acquisition date (`--date`).
5. Whether plot values should include soil (whole-plot NDVI) or only vegetation (`--mask`); this changes what the numbers mean.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
# Demo: synthetic 5-band reflectance (16-bit x 10000), 12 plots, Otsu soil mask, -0.25 m buffer
Rscript "${CLAUDE_SKILL_DIR}/scripts/vi_extract.R" --demo

# UAV reflectance orthomosaic, band order confirmed with the user, inner plots from claw-plot-grid
Rscript "${CLAUDE_SKILL_DIR}/scripts/vi_extract.R" --input ortho_reflectance.tif \
  --plots grid/spatial/plots_inner.gpkg --bands "blue=1,green=2,red=3,rededge=4,nir=5" \
  --date 2026-07-15 --out results/vi_2026-07-15

# Canopy-only indices: drop soil with one Otsu NDVI threshold for the whole flight, trim 0.2 m
Rscript "${CLAUDE_SKILL_DIR}/scripts/vi_extract.R" --input ortho_reflectance.tif --plots plots.gpkg \
  --bands "blue=1,green=2,red=3,rededge=4,nir=5" --mask otsu --buffer -0.2 --date 2026-07-15 --out results/vi

# Sentinel-2 L2A from ESA (baseline >= 04.00), layers named B02...B12
Rscript "${CLAUDE_SKILL_DIR}/scripts/vi_extract.R" --input S2_L2A_stack.tif --plots fields.gpkg \
  --sensor sentinel2-l2a --date 2026-07-12 --out results/s2
```

## Outputs

- `report.md`: warnings first; inputs, band mapping and its source, a reflectance check per band, the mask threshold, a summary per index.
- `tables/plot_indices.csv`: one row per plot: id, `date`, the plot attributes, one column per index (`ndvi`, `ndre`, `gndvi`, `osavi`, `evi2`, `mcari2`, `cire`: coverage-weighted means), `n_pixels` (effective pixels used), `n_pixels_all` (before masking), `veg_fraction` (kept / all), `n_pure_pixels` (pixels wholly inside), `plot_area_m2`, `pixels_per_plot`, `flag`.
- `tables/plot_index_stats.csv`: one row per plot x index: `mean`, `median`, `p10`, `p90`, `sd` (population, coverage-weighted), `n_pixels`.
- `tables/plot_band_reflectance.csv`: mean reflectance per band (after masking), to compute indices of mean reflectance if wanted.
- `tables/image_bands.csv`: raw and reflectance quantiles per band inside the plots.
- `spatial/plot_indices.gpkg`, `figures/index_map.png`, `figures/ndvi_histogram.png`, `figures/plot_means.png`, `reproducibility/`.

## Expert decisions

| Index | Formula (reflectance) | Source |
|---|---|---|
| NDVI | (NIR - R) / (NIR + R) | Rouse et al. 1974 |
| NDRE | (NIR - RE) / (NIR + RE) | red-edge normalised difference (750/705 nm) of Gitelson & Merzlyak 1994 |
| GNDVI | (NIR - G) / (NIR + G) | Gitelson et al. 1996 |
| OSAVI | (NIR - R) / (NIR + R + 0.16) | Rondeaux et al. 1996 (some software multiplies by 1.16) |
| EVI2 | 2.5 (NIR - R) / (NIR + 2.4 R + 1) | Jiang et al. 2008, eq. 15 |
| MCARI2 | 1.5 [2.5 (NIR - R) - 1.3 (NIR - G)] / sqrt((2 NIR + 1)^2 - (6 NIR - 5 sqrt(R)) - 0.5) | Haboudane et al. 2004, eq. 12 (defined for 800, 670, 550 nm) |
| CIred-edge (`cire`) | NIR / RE - 1 | Gitelson et al. 2003, 2005 |

- **Reflectance only.** Indices with additive constants (OSAVI, EVI2, MCARI2) are meaningless on digital numbers, and even ratio indices are wrong because each band has its own gain. The run stops when the 99th percentile of any band exceeds 1.5 after scaling (rule of thumb: reflectance rarely exceeds 1); `--confirm-calibrated` overrides it only when the user confirms the values are reflectance.
- **Scaling.** reflectance = value x scale + offset. Sentinel-2 L2A from baseline 04.00: (DN + BOA_ADD_OFFSET) / QUANTIFICATION_VALUE = (DN - 1000) / 10000, DN 0 = no data (ESA SentiWiki, Sentinel-2 products). PlanetScope surface reflectance: DN / 10000.
- **Statistics.** Indices are computed per pixel, then averaged with each pixel weighted by the fraction inside the plot (exactextractr); median and quantiles are coverage-weighted too (they equal R's type-7 quantiles when all pixels are whole). Intermediate rasters are float32 (about 7 significant digits).
- **Soil and shadow mask** (off by default, because whole-plot NDVI is also a valid trait). `--mask otsu` sets one NDVI threshold for the whole image with Otsu's method on pixels inside the plots, never per plot (a plot that is all canopy or all soil would be split anyway). A threshold outside 0.1-0.8, or an Otsu separability (between-class / total variance) below 0.8, is flagged as a non-bimodal histogram (rule of thumb; a single normal mode gives 2/pi = 0.64). `--mask ndvi` uses a fixed threshold: soil is usually below about 0.3 and green canopy above about 0.5 in UAV reflectance, so 0.3-0.4 is a common starting point (rule of thumb; check the histogram). `--shadow-nir` drops pixels darker than a NIR reflectance (rule of thumb 0.10-0.15).
- **Buffers** are applied in metres (a projected CRS), never in degrees. A negative buffer of at least one to two pixels (rule of thumb, often 0.2-0.5 m) keeps alley, edge and neighbour pixels out.
- **Mixed pixels.** Each plot reports effective and whole-pixel counts. Plots with fewer than `--min-pixels` (10, rule of thumb) effective pixels or no whole pixel are flagged; a median below 10 pixels per plot raises a warning (for example Sentinel-2 10 m pixels over 1.5 x 5 m plots: 0.075 pixels per plot).
- **Band-order checks.** A layer named, say, `NIR` mapped to red stops the run. Pixels with NDVI below -0.3 at the 5th percentile suggest swapped red and NIR; a 95th percentile below 0.2 suggests no vegetation signal.

## Checks before reporting

- Read the Warnings section first and pass each warning on.
- Confirm the band mapping line in `report.md` with the user unless it says `user (--bands)`.
- Look at `figures/ndvi_histogram.png` when masking: two clear modes, the threshold between them.
- **Stop and tell the user** if: the run refused digital numbers; more than 5% of pixels are negative in any band (offset applied twice?); any band exceeds reflectance 1 widely; NDVI looks swapped or flat; plots are flagged for few pixels or incomplete coverage; or median pixels per plot is below 10.

## References

- Rouse JW, Haas RH, Schell JA, Deering DW (1974) Monitoring vegetation systems in the Great Plains with ERTS. Third ERTS-1 Symposium, NASA Goddard Space Flight Center, vol. 1, sect. A. https://ntrs.nasa.gov/citations/19740022614
- Gitelson A, Merzlyak MN (1994) Spectral reflectance changes associated with autumn senescence of Aesculus hippocastanum L. and Acer platanoides L. leaves. Journal of Plant Physiology 143:286-292. https://doi.org/10.1016/S0176-1617(11)81633-0
- Gitelson AA, Kaufman YJ, Merzlyak MN (1996) Use of a green channel in remote sensing of global vegetation from EOS-MODIS. Remote Sensing of Environment 58:289-298. https://doi.org/10.1016/S0034-4257(96)00072-7
- Rondeaux G, Steven M, Baret F (1996) Optimization of soil-adjusted vegetation indices. Remote Sensing of Environment 55:95-107. https://doi.org/10.1016/0034-4257(95)00186-7
- Jiang Z, Huete AR, Didan K, Miura T (2008) Development of a two-band enhanced vegetation index without a blue band. Remote Sensing of Environment 112:3833-3845. https://doi.org/10.1016/j.rse.2008.06.006
- Haboudane D, Miller JR, Pattey E, Zarco-Tejada PJ, Strachan IB (2004) Hyperspectral vegetation indices and novel algorithms for predicting green LAI of crop canopies. Remote Sensing of Environment 90:337-352. https://doi.org/10.1016/j.rse.2003.12.013
- Gitelson AA, Viña A, Arkebauer TJ, Rundquist DC, Keydan G, Leavitt B (2003) Remote estimation of leaf area index and green leaf biomass in maize canopies. Geophysical Research Letters 30(5). https://doi.org/10.1029/2002GL016450
- Gitelson AA, Viña A, Ciganda V, Rundquist DC, Arkebauer TJ (2005) Remote estimation of canopy chlorophyll content in crops. Geophysical Research Letters 32(8). https://doi.org/10.1029/2005GL022688
- Otsu N (1979) A threshold selection method from gray-level histograms. IEEE Transactions on Systems, Man, and Cybernetics 9:62-66. https://doi.org/10.1109/TSMC.1979.4310076
- Drusch M, Del Bello U, Carlier S, et al. (2012) Sentinel-2: ESA's optical high-resolution mission for GMES operational services. Remote Sensing of Environment 120:25-36. https://doi.org/10.1016/j.rse.2011.11.026
- ESA SentiWiki, Sentinel-2 products (Level-2A radiometric offset from processing baseline 04.00). https://sentiwiki.copernicus.eu/web/s2-products
- Planet, PlanetScope product documentation (band order). https://docs.planet.com/data/imagery/planetscope/
- Baston D. exactextractr: fast extraction from raster datasets using polygons. R package. https://CRAN.R-project.org/package=exactextractr
