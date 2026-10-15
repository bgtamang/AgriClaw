---
name: claw-ortho-qc
description: Quality-check a drone/UAV orthomosaic (and optional DSM/DTM) before per-plot extraction - CRS and projection, ground sampling distance (pixel size) against plot size, rotation, extent versus the plot grid, band count and data type, reflectance versus digital numbers, declared NoData and zero collars, NoData holes and plot coverage, saturation, DSM elevation sanity, spikes and canopy height - each check pass, warn or fail against stated thresholds. Includes an OpenDroneMap processing guide (radiometric calibration with reflectance panels, ground control points, overlap and GSD trade-offs). Use when the user asks to check, validate or QC an orthomosaic, orthophoto, DSM or drone image, asks why plot values look wrong, or asks how to process drone imagery for plot trials.
license: MIT
compatibility: R >= 4.2 with terra, sf and exactextractr (GDAL, GEOS, PROJ). No internet needed. The processing guide is advice for OpenDroneMap; nothing in it is executed.
metadata:
  version: "1.0.0"
  domain: remote-sensing
  language: R
  status: tested
  summary: "Check a drone orthomosaic and DSM before plot extraction; processing guide"
  repository: https://github.com/bgtamang/AgriClaw
---

# Orthomosaic and DSM quality checks

## What it does

Runs a fixed list of checks on an orthomosaic (and optionally the plot grid, a DSM and a DTM) and reports each as pass, warn or fail against a stated threshold, with a per-plot table of coverage, saturation and canopy height. It catches the problems that silently corrupt plot values: wrong or missing CRS, digital numbers instead of reflectance, NoData holes, glare, plots off the image, and DSM artefacts. `references/processing-guide.md` gives OpenDroneMap processing advice for plot trials.

## When to use / when not to use

- Use on every new orthomosaic before `claw-vi-extract` or `claw-canopy-cover`, with the plot polygons from `claw-plot-grid` so the plot-level checks run.
- Use the processing guide when the user asks how to fly or process imagery (overlap, pixel size, GCPs, reflectance panels, ODM options).
- It does not make orthomosaics (use OpenDroneMap or commercial software), fix problems, or measure geolocation accuracy (that needs independent check points, see the guide).

## Inputs

- `--input`: the orthomosaic (GeoTIFF or any GDAL raster).
- `--plots` (optional, recommended): plot polygons; enables extent, coverage, saturation-per-plot and pixel-size-versus-plot checks.
- `--dsm`, `--dtm` (optional): surface and bare-ground terrain models in metres.
- `--expect`: what the values should be: `reflectance` (multispectral reflectance, 0-1 after `--scale`), `rgb`, `dn` (raw digital numbers) or `any` (default: describe only).

**Ask the user, never guess:**

1. What the file should contain (reflectance, RGB, raw DN) and, for integer reflectance, its scale factor (for example 0.0001).
2. How many bands it should have and their order (needed later by `claw-vi-extract`).
3. The vertical units and datum of the DSM/DTM (metres above the ellipsoid or sea level?), and how the DTM was made.
4. The tallest plausible crop height for `--max-height-m` (default 4 m suits maize; use about 1.5 m for soybean).

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
# Demo: synthetic reflectance mosaic, plots, DSM and DTM with deliberately injected problems
Rscript "${CLAUDE_SKILL_DIR}/scripts/ortho_qc.R" --demo

# Multispectral reflectance mosaic (float 0-1) with the plot grid
Rscript "${CLAUDE_SKILL_DIR}/scripts/ortho_qc.R" --input ortho_reflectance.tif \
  --plots grid/spatial/plots.gpkg --expect reflectance --expect-bands 5 --out results/qc_0715

# RGB mosaic with surface and terrain models; stop a pipeline on any failed check
Rscript "${CLAUDE_SKILL_DIR}/scripts/ortho_qc.R" --input rgb.tif --plots plots.gpkg --dsm dsm.tif \
  --dtm dtm_bare_soil.tif --expect rgb --max-height-m 1.5 --strict --out results/qc_rgb
```

## Outputs

- `report.md`: the QC result line (counts of fail, warn, pass), warnings, a table of every check (value, threshold, note), band statistics and the worst plots.
- `tables/qc_checks.csv`: `check`, `status` (pass, warn, fail, info), `value`, `threshold`, `message`.
- `tables/band_stats.csv`: per band: share of NoData, min, 1%, median, 99%, max, share negative, above 1 and zero (inside the plots when given).
- `tables/plot_qc.csv`: per plot: `inside_extent`, `coverage` (share of the plot with data), `saturated_fraction`, and with a DSM `dsm_coverage`, `dsm_spike_fraction`, `chm_median_m`, `chm_p95_m`.
- `figures/overview.png` (plots outlined by coverage status), `figures/band_histograms.png`, `figures/dsm.png`; `reproducibility/`.
- With `--strict`, the script exits with an error after writing everything if any check fails.

## Expert decisions

| Check | Pass / warn / fail rule | Why |
|---|---|---|
| `rotation` | fail if the geotransform has rotation terms | most tools mis-place rotated pixels; warp to north-up first |
| `crs` | fail if missing; warn if lon/lat or not metres | plots, buffers and areas need metres |
| `pixel_shape` | warn if x and y pixel sizes differ by > 1% | anisotropic resampling |
| `gsd` | warn above `--max-gsd-cm` | user's requirement |
| `gsd_vs_plot` | warn < 5, fail < 1 pixel across the narrowest plot side | mixed pixels (rule of thumb) |
| `band_count` | fail if not `--expect-bands`; RGB warns unless 3 or 4 | wrong file or stack |
| `datatype` | warn for 8-bit data when reflectance is expected | 8-bit reflectance is a quantised display product |
| `nodata_declared` | warn if NoData equals the saturation level (255, 65535) | glare would read as missing data |
| `zero_as_nodata` | warn if > 1% of pixels are 0 in every band and no NoData is declared | collars read as data (rule of thumb) |
| `extent` | fail if a plot is entirely outside the image; warn if partly outside | wrong field, CRS mismatch, truncated mosaic |
| `value_range` (reflectance) | fail if any band's 99th percentile > 1.5; warn if > 1% negative or > 1% above 1 | digital numbers or a missing scale; offset applied twice; glare (rules of thumb) |
| `plot_coverage` | warn < 0.95, fail < 0.75 of a plot's area with data | holes and edges bias plot means (rules of thumb) |
| `saturation` | warn > 0.5%, fail > 5% of pixels at the saturation level (worst plot) | clipped colour and reflectance (rules of thumb) |
| `dsm_values` | fail for 8-bit DSMs or elevations outside -500 to 9000 m | not elevations, wrong units, or unflagged NoData |
| `dsm_coverage` | as `plot_coverage` | holes in the surface model |
| `dsm_spikes` | warn > 0.1%, fail > 1% of plot pixels more than `--spike-m` (1 m) from their 5 x 5 median | matching errors, birds, dust, moving leaves (rules of thumb) |
| `dsm_relief` | warn if the 1st-99th percentile range over plots exceeds 20 m | artefacts or wrong units in a trial field (rule of thumb) |
| `canopy_height` | warn if a plot's median DSM - DTM < -0.2 m or its 95th percentile > `--max-height-m` | DSM/DTM vertical offset or misregistration; spikes or the wrong DTM (rules of thumb) |

- Value statistics use a seeded random sample of pixels at least half inside the plots (`--sample-size`, `--seed`), or a regular grid over the whole image without plots. Coverage, saturation and spikes per plot use every pixel, weighted by the fraction inside the plot (exactextractr).
- A band 4 holding only 0 and the data-type maximum is treated as an alpha mask (`--alpha-band 0` turns this off).
- terra and GDAL write 255 (8-bit) or 65535 (16-bit) as NoData by default; in a camera product that value is also saturation. This is the reason for `nodata_declared`.

## Checks before reporting

- Read the QC result line and the Warnings section; pass on every warn and fail with its message.
- Open `figures/overview.png`: plot outlines on the plots, orange or red outlines where coverage is low.
- **Stop and tell the user** before extracting values if any check fails, if `value_range` says the data are digital numbers, if plots fall outside the image, or if `canopy_height` warns (the DSM and DTM do not line up). Point to `references/processing-guide.md` for the likely processing fix.

## References

- OpenDroneMap documentation (flying, ground control points, multispectral, options). https://docs.opendronemap.org/
- James MR, Chandler JH, Eltner A, et al. (2019) Guidelines on the use of structure-from-motion photogrammetry in geomorphic research. Earth Surface Processes and Landforms 44:2081-2084. https://doi.org/10.1002/esp.4637
- Sanz-Ablanedo E, Chandler JH, Rodríguez-Pérez JR, Ordóñez C (2018) Accuracy of unmanned aerial vehicle (UAV) and SfM photogrammetry survey as a function of the number and location of ground control points used. Remote Sensing 10:1606. https://doi.org/10.3390/rs10101606
- Assmann JJ, Kerby JT, Cunliffe AM, Myers-Smith IH (2019) Vegetation monitoring using multispectral sensors - best practices and lessons learned from high latitudes. Journal of Unmanned Vehicle Systems 7:54-75. https://doi.org/10.1139/juvs-2018-0018
- Aasen H, Honkavaara E, Lucieer A, Zarco-Tejada PJ (2018) Quantitative remote sensing at ultra-high resolution with UAV spectroscopy. Remote Sensing 10:1091. https://doi.org/10.3390/rs10071091
- Smith GM, Milton EJ (1999) The use of the empirical line method to calibrate remotely sensed data to reflectance. International Journal of Remote Sensing 20:2653-2662. https://doi.org/10.1080/014311699211994
- Holman F, Riche A, Michalski A, Castle M, Wooster M, Hawkesford M (2016) High throughput field phenotyping of wheat plant height and growth rate in field plot trials using UAV based remote sensing. Remote Sensing 8:1031. https://doi.org/10.3390/rs8121031
