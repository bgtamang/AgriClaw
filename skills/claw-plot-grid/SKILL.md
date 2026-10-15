---
name: claw-plot-grid
description: Build plot polygons (GeoPackage) for a field trial from the four field corners, or one corner plus a bearing, and the layout - ranges x columns (rows), plot length, plot width (rows x row spacing), alleys - with serpentine or row-wise plot numbering from any start corner. Converts lon/lat to the local UTM zone so plots, buffers and areas are in meters, trims plot ends and sides, joins a claw-trial-design field book, and draws a numbered plot map. Use when the user asks for a plot grid, plot boundaries, plot shapefile or polygons, a field map, or needs plots to extract drone/UAV or satellite data per plot or to match plot numbers to an orthomosaic.
license: MIT
compatibility: R >= 4.2 with sf (GDAL, GEOS, PROJ); the tests also use terra. No internet needed.
metadata:
  version: "1.0.0"
  domain: remote-sensing
  language: R
  status: tested
  summary: "Plot polygons from field corners and layout, joined to the field book"
  repository: https://github.com/bgtamang/AgriClaw
---

# Plot grid from field corners and layout

## What it does

Builds one rectangle per plot from the field geometry (four corners, or one corner plus the planting bearing) and the trial layout, numbers the plots the way the planter did (serpentine or row-wise, from any corner), and writes them as a GeoPackage in meters (the local UTM zone by default). Optional inner polygons trim plot ends and sides for clean extraction, and a field book can be joined so each polygon carries its entry, rep and genotype.

## When to use / when not to use

- Use before `claw-vi-extract` and `claw-canopy-cover`, once per field (the same polygons serve every flight). Run `claw-ortho-qc` with the polygons to check they sit on the orthomosaic.
- Don't use it to randomize a trial: use `claw-trial-design`, then join its field book here.
- Don't use it for irregular plots (curved rows, missing blocks, plots of different sizes in one block). Digitize those by hand in QGIS, or run this skill once per regular block.
- If the user already has plot polygons (QGIS, the planter's as-planted map, or FIELDimageR, which builds a plot shapefile after cropping and rotating the orthomosaic; Matias et al. 2020), use those and only check them with `claw-ortho-qc`. This skill builds polygons from surveyed geometry instead, so the same grid fits every flight.

## Inputs

Layout terms used here (see `references/layout-conventions.md` for diagrams):

- **Range**: a tier of plots along the planting direction; ranges are separated by alleys. Range 1 is at corner c1.
- **Column** (a "row" or "pass" in many field books): plots side by side across the planting direction. Column 1 is at c1.
- **Corners**: c1 = outer corner of the plot at range 1, column 1; c2 = far end of column 1 along the planting direction; c3 = diagonal from c1; c4 = across from c1. Corners are the outer edges of the plot block, not stakes in the border or headland.

| Need | Options | Notes |
|---|---|---|
| Field geometry | `--corners` (4 points, CSV or string) or `--origin` + `--bearing` | lon,lat order (x first) for EPSG:4326 |
| Coordinate CRS | `--crs` (default EPSG:4326) | `local` for unreferenced meters |
| Layout | `--ranges --columns --plot-length --alley` and `--rows-per-plot --row-spacing` (or `--plot-width`) | all lengths in `--units` |
| Units | `--units m` or `--units ft` | required, never assumed |
| Numbering | `--numbering --sweep --start-corner --first-plot` | defaults are flagged in the report |
| Field book | `--fieldbook --join-by [--fb-coords]` | CSV, one row per plot |

**Ask the user, never guess:**

1. The units of every length (ft or m), and whether plot length excludes the alley (it must: give the alley separately).
2. Coordinate order and CRS: Google Maps gives "lat, lon"; this skill wants lon,lat.
3. For a bearing: is it true, grid (measured on the UTM orthomosaic) or magnetic? Magnetic bearings must first be corrected by the local declination (NOAA calculator); this skill does not do that.
4. Where plot 1 is, which way the planter traveled, and whether numbering snakes (serpentine) or restarts at the same end (row-wise). Show `figures/plot_map.png` and get a yes.
5. Whether the corners are the outer edges of the plot block. If they are stakes outside the block, use `--origin` with `--offset-along/--offset-across`.

## How to run

Commands use `${CLAUDE_SKILL_DIR}`, this skill's folder: Claude Code fills it in; in other tools, replace it with the folder that contains this SKILL.md.

```bash
# Demo: 6 ranges x 8 columns of 2-row soybean plots from RTK corners, RCBD field book joined
Rscript "${CLAUDE_SKILL_DIR}/scripts/plot_grid.R" --demo

# Four RTK corners (lon/lat), metric layout, field book from claw-trial-design
# (its field rows are ranges here and its serpentine numbering runs across columns)
Rscript "${CLAUDE_SKILL_DIR}/scripts/plot_grid.R" --corners corners.csv --ranges 20 --columns 12 \
  --plot-length 5.2 --alley 0.9 --rows-per-plot 2 --row-spacing 0.76 --units m \
  --numbering serpentine --sweep across --start-corner 1 --trim-ends 0.5 --trim-sides 0.15 \
  --fieldbook field_book.csv --join-by "plot=plot" --fb-coords "row,col" --out results/grid

# One corner and a true bearing, layout in feet (17.5 ft plots, 2.5 ft alleys, 30 in rows)
Rscript "${CLAUDE_SKILL_DIR}/scripts/plot_grid.R" --origin "-88.2251,40.0623" --bearing 3.5 \
  --bearing-ref true --across right --ranges 10 --columns 16 --plot-length 17.5 --alley 2.5 \
  --rows-per-plot 2 --row-spacing 2.5 --units ft --out results/grid

# Corners digitized on the orthomosaic (CSV with x,y in UTM 15N), layout stretched to fit them
Rscript "${CLAUDE_SKILL_DIR}/scripts/plot_grid.R" --corners corners_utm.csv --crs EPSG:32615 \
  --fit-corners --ranges 20 --columns 12 --plot-length 5.2 --alley 0.9 --rows-per-plot 2 \
  --row-spacing 0.76 --units m --out results/grid
```

## Outputs

- `report.md`: warnings first, then geometry (CRS, bearing true and grid, scale factor, block size), the corner check, the overlap check and the field-book join summary.
- `spatial/plots.gpkg` (layer `plots`): one polygon per plot with `plot`, `range`, `column`, `plot_id` (e.g. `R03C12`), `length_m`, `width_m`, `area_m2` (ground area), centroid `x`, `y` (output CRS), `lon`, `lat`, the field-book columns and `matched` (has a field-book row).
- `spatial/plots_inner.gpkg`: trimmed polygons (only with `--trim-ends`/`--trim-sides`); use these for extraction.
- `tables/plots.csv`: the attributes without geometry. `tables/corner_check.csv`: given versus fitted corners.
- `figures/plot_map.png`: plots with numbers, the numbering path from plot 1 (red dot), corners c1-c4, grid north and a scale bar.
- `reproducibility/`: command, session info, checksums.

## Expert decisions

- **Metric CRS.** Lon/lat input is projected to the WGS84 UTM zone of the field (EPSG 326xx north, 327xx south, with the south-west Norway and Svalbard exceptions). `--out-crs` must be projected in meters; plots are never built or buffered in degrees. Use the orthomosaic's CRS when it is not UTM.
- **Ground versus map lengths.** Plot dimensions are ground lengths. A UTM map stretches them by the point scale factor k (0.9996 on the central meridian, up to about 1.001 at zone edges), so polygons are built with k x length; without this a 100 m block drifts by up to 10 cm. k is computed from the projection at the field (it matches Snyder 1987, eq. 8-11, to 1e-9). `--scale-factor 1` turns it off, for example when dimensions were measured on the orthomosaic itself.
- **True versus grid north.** A true bearing is converted with the meridian convergence at the field (up to about 3 deg in UTM); the report gives both bearings. A grid bearing is relative to the output CRS.
- **Corner fit.** By default the layout is placed on the corners by a rigid least-squares fit (rotation and shift, no stretching), which tolerates a few centimeters of corner error and reports each corner's misfit. Misfit over 0.10 m is a warning; over `--corner-tol` (0.5 m) the run stops, because the layout and corners disagree (usually ft versus m, the number of ranges, or the alley). Both thresholds are adjustable rules of thumb. `--fit-corners` instead stretches the layout to the four corners by bilinear interpolation (for corners clicked on the orthomosaic when the dimensions are approximate); a stretch over 2% is flagged.
- **Plot width** = rows per plot x row spacing: plot edges lie half a row spacing outside the outer rows, so side-by-side plots tile without gaps. Use `--col-gap` only for extra unplanted space between columns.
- **Trimming.** Ends and sides are trimmed separately because alley-edge effects (ends) and neighbor-canopy overlap (sides) differ. Typical values (rules of thumb): 0.3-0.5 m off each end; on the sides, enough to clear canopy leaning in from the next plot (often a quarter to half a row spacing).
- **Dimension sanity** (rules of thumb, warnings only): plot length 1-30 m, width 0.3-20 m, row spacing 0.1-1.6 m, alley up to 6 m. Values outside usually mean a ft/m mix-up.
- **Field-book join.** Duplicate keys and field-book rows that match no plot stop the run; grid plots without a row (fillers, borders) are counted as a warning. `--fb-coords` checks that the field book's own position columns map one-to-one onto ranges and columns; if not, the plot numbers were assigned with a different serpentine direction or start corner, and the run stops.

## Checks before reporting

- Open `figures/plot_map.png`: plot 1 at the right corner, the numbering path matching the planter, columns growing toward the right side. Then overlay `plots.gpkg` on the orthomosaic (QGIS, or `claw-ortho-qc --plots`). Polygons sitting on alleys or shifted by a plot are the most common error.
- Read the Warnings section first and pass every warning on.
- **Stop and tell the user** if: the corner misfit exceeds 0.10 m; `--fit-corners` stretches more than 2%; any dimension warning appears; the field book has unmatched plots; the report says numbering defaults were used and the user has not confirmed the scheme; or the true bearing differs from what the user expects by more than a degree or two (possible magnetic bearing or swapped corners).

## References

- Matias FI, Caraza-Harter MV, Endelman JB (2020) FIELDimageR: an R package to analyze orthomosaic images from agricultural field trials. The Plant Phenome Journal 3:e20005. https://doi.org/10.1002/ppj2.20005
- Snyder JP (1987) Map projections: a working manual. USGS Professional Paper 1395. https://doi.org/10.3133/pp1395
- Pebesma E (2018) Simple features for R: standardized support for spatial vector data. The R Journal 10:439-446. https://doi.org/10.32614/RJ-2018-009
