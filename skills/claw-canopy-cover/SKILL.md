---
name: claw-canopy-cover
description: Estimate per-plot green canopy cover (fractional vegetation cover, percent ground cover) from RGB drone/UAV orthomosaics with the excess-green index (ExG) and one Otsu threshold per flight, or ExG minus ExR with a zero threshold. Counts and excludes saturated (glare), dark (shadow) and no-data pixels, cross-checks the two methods, and writes a vegetation mask plus QC images. Use when the user asks for canopy cover, ground cover, green cover percentage, vegetation fraction, canopy closure or early vigour per plot from RGB or drone images.
license: MIT
compatibility: R >= 4.2 with terra, sf and exactextractr (GDAL, GEOS, PROJ). No internet needed.
metadata:
  version: "1.0.0"
  domain: remote-sensing
  language: R
  status: tested
  repository: https://github.com/bgtamang/AgriClaw
---

# Canopy cover from RGB imagery

## What it does

Classifies every pixel of an RGB orthomosaic as green vegetation or background with a colour index on chromatic coordinates, using one threshold for the whole flight, and reports the canopy cover of each plot (coverage-weighted share of usable pixels). Saturated, very dark and no-data pixels are counted and left out, and a second method is computed as a cross-check.

## When to use / when not to use

- Use for green cover from RGB (or the red, green and blue bands of a multispectral camera), typically from emergence to canopy closure, once per flight. Polygons come from `claw-plot-grid`; check the image first with `claw-ortho-qc`. Stack per-date tables for `claw-vi-timeseries`.
- For NDVI or other indices from multispectral reflectance, use `claw-vi-extract` (its `--mask` gives an NDVI-based vegetation fraction).
- Not for counting plants, measuring height (DSM), or senescence scoring: green-index methods treat yellow and brown leaves as background.
- Near full canopy (cover above about 90%) colour thresholds saturate; use height or NDVI-based traits instead.

## Inputs

- `--input`: RGB orthomosaic (8-bit, 16-bit or float GeoTIFF) with a CRS. A fourth band holding only 0 and the maximum is used as an alpha (transparency) mask.
- `--plots`: plot polygons with a unique id (`--id-col`, default `plot_id`, else `plot`).
- `--bands`: layer order (default `red=1,green=2,blue=3`, the usual RGB order).

**Ask the user, never guess:**

1. The band order if the file is not a plain RGB export (multispectral stacks, BGR exports).
2. The flight date (`--date`) and crop stage; with ground cover above about 90% the result is not informative.
3. Whether the field has green weeds, green crop residue or algae: these count as canopy.
4. The ground sampling distance; early-stage cover needs pixels of 1-2 cm or finer.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
# Demo: synthetic 8-bit RGB orthomosaic, 9 plots with known cover, residue, glare and shadow
Rscript "${CLAUDE_SKILL_DIR}/scripts/canopy_cover.R" --demo

# RGB orthomosaic, ExG with one Otsu threshold for the flight, plots trimmed by 0.1 m
Rscript "${CLAUDE_SKILL_DIR}/scripts/canopy_cover.R" --input rgb_ortho.tif \
  --plots grid/spatial/plots_inner.gpkg --buffer -0.1 --date 2026-06-20 --out results/cover_0620

# Residue-covered field: ExG - ExR with its built-in zero threshold
Rscript "${CLAUDE_SKILL_DIR}/scripts/canopy_cover.R" --input rgb_ortho.tif --plots plots.gpkg \
  --index exgr --out results/cover_exgr
```

## Outputs

- `report.md`: warnings first; inputs, threshold (value, method, Otsu separability), cross-check agreement, a table per plot.
- `tables/canopy_cover.csv`: one row per plot: id, `date`, plot attributes, `canopy_cover_pct` (primary method), `cover_exgr_pct` or `cover_exg_pct` (cross-check), `n_pixels` (usable, coverage-weighted), `valid_fraction` (plot area with image data), `saturated_fraction`, `dark_fraction`, `plot_area_m2`, `flag`.
- `spatial/veg_mask.tif`: 1 = vegetation, 0 = background, 2 = left out (saturated or dark), no data elsewhere.
- `figures/index_histogram.png` (index histogram with the threshold), `figures/cover_map.png` (mask over the orthomosaic with plot labels), `figures/qc_plots.png` (RGB and mask for plots spanning the cover range and flagged plots).
- `reproducibility/`.

## Expert decisions

- **Index.** Chromatic coordinates r = R/(R+G+B), g = G/(R+G+B), b = B/(R+G+B) remove overall brightness. ExG = 2g - r - b (Woebbecke et al. 1995). ExR = 1.4r - g and ExGR = ExG - ExR = 3g - 2.4r - b, classified with a fixed zero threshold (Meyer & Neto 2008, who found it worked especially well on wheat-straw residue).
- **One threshold per flight.** The Otsu (1979) threshold is computed once from pixels pooled across all plots, never per plot: a plot that is all canopy or all soil has a one-mode histogram that Otsu would split anyway (the tests show a per-plot threshold calling about half of a full canopy green).
- **Bimodality check.** Otsu's separability (between-class / total variance) is 1 for two spikes, 0.75 for a uniform and 2/pi = 0.64 for a single normal histogram; below 0.8 (rule of thumb) the report warns that the image may lack one of the two classes.
- **Cross-check.** Both ExG-Otsu and ExGR > 0 are reported; a median disagreement over 10 percentage points (rule of thumb) is flagged, because it points to residue, a colour cast (white balance) or a poor threshold.
- **Saturated pixels** (any channel at the data-type maximum, 255 or 65535, or `--max-value`) have distorted colour and are left out. Files that declare NoData = 255 hide glare as missing data; the report says so.
- **Dark pixels** (mean of R, G, B below 5% of full scale, `--dark-threshold`, rule of thumb) are left out: with 8-bit data a one-level change in G then moves ExG by about 0.1, so their class is noise. Deep shadow therefore lowers `n_pixels` instead of being guessed.
- **Plot flags** (rules of thumb): saturated > 1%, dark > 10%, image coverage < 95%, no usable pixels. A median plot under 100 pixels raises a warning.
- **Statistics.** Cover is the coverage-weighted mean of the 0/1 vegetation layer over usable pixels (exactextractr), so pixels cut by the plot edge count by the fraction inside.

## Checks before reporting

- Open `figures/qc_plots.png` and `figures/cover_map.png`: magenta should sit on leaves only, not on soil, residue or shadows. This visual check is essential.
- Check `figures/index_histogram.png`: two modes with the threshold between them.
- **Stop and tell the user** if: the separability warning appears; the two methods disagree by a median of more than 10 points; plots are flagged for glare, shadow or incomplete coverage; the NoData = 255 warning appears; or the pixel size is coarser than about 2 cm for early-stage plots.

## References

- Woebbecke DM, Meyer GE, Von Bargen K, Mortensen DA (1995) Color indices for weed identification under various soil, residue, and lighting conditions. Transactions of the ASAE 38:259-269. https://doi.org/10.13031/2013.27838
- Meyer GE, Neto JC (2008) Verification of color vegetation indices for automated crop imaging applications. Computers and Electronics in Agriculture 63:282-293. https://doi.org/10.1016/j.compag.2008.03.009
- Otsu N (1979) A threshold selection method from gray-level histograms. IEEE Transactions on Systems, Man, and Cybernetics 9:62-66. https://doi.org/10.1109/TSMC.1979.4310076
- Hamuda E, Glavin M, Jones E (2016) A survey of image processing techniques for plant extraction and segmentation in the field. Computers and Electronics in Agriculture 125:184-199. https://doi.org/10.1016/j.compag.2016.04.024
- Baston D. exactextractr: fast extraction from raster datasets using polygons. R package. https://CRAN.R-project.org/package=exactextractr
