# Processing guide: orthomosaics and DSMs for plot trials (OpenDroneMap)

**This guide is advice, not code.** `claw-ortho-qc` does not run OpenDroneMap and cannot check most of what follows; it only checks the products afterwards. Numbers come from the cited sources; everything marked *rule of thumb* is a starting point to adapt to the crop, sensor and site.

## 1. Plan the flight

- **Pixel size (ground sampling distance, GSD).** GSD = flying height x sensor pixel size / focal length. Read the camera's GSD at a given height from its specification sheet. Choose the height so that the narrowest plot, after trimming its edges, still spans enough pixels: *rule of thumb* at least 5-10 pixels for plot-mean vegetation indices, and pixels of 1-2 cm or finer for canopy cover of young plants.
- **Overlap.** OpenDroneMap recommends 70-80% overlap, slightly off-nadir (5-10 degrees), for 2D and 2.5D products, and closer to 80-83% for complex vegetation (ODM, Flying tips). Uniform canopies are hard to match, so stay at the high end over crops.
- **Time of day and weather.** Fly close to solar noon and at the same time of day on every date (*rule of thumb*: within about two hours of solar noon), under stable light (clear sky or uniform overcast, not passing clouds) and light wind: moving leaves blur the mosaic and spoil surface models. Assmann et al. (2019) give a standardised protocol for flight planning, weather, geolocation and radiometric calibration with multispectral sensors.
- **Repeat flights identically.** Same flight plan, height, overlap and camera settings on every date, so that products are comparable and one plot grid (`claw-plot-grid`) fits every flight.

## 2. Georeferencing with ground control points (GCPs)

- **How many and where.** A minimum of 5 GCPs works for most jobs and 8-10 suffice for larger projects; spread them evenly, some near the corners and some near the centre, covering high and low ground; space them further apart than one image footprint; keep them 10-30 m inside the project perimeter; make each GCP visible in at least 5 images (ODM, Ground control points).
- **File format.** `gcp_list.txt` starts with the projection (EPSG code, PROJ string or `WGS84 UTM 15N`), followed by one line per observation: `geo_x geo_y geo_z im_x im_y image_name [gcp_name]` (ODM, Ground control points). Use 0.0, not NaN, when an elevation is unknown.
- **Independent check points.** Keep some surveyed points out of the adjustment and report their errors: accuracy judged only from the control points is overestimated, and accuracy improves with the number and even distribution of GCPs (Sanz-Ablanedo et al. 2018). James et al. (2019) set out how to design, report and assess structure-from-motion surveys.
- **Permanent markers.** Fixed or RTK-surveyed GCPs used on every date keep the flights co-registered. Check co-registration by overlaying the plot grid on each date's mosaic.

## 3. Radiometric calibration (multispectral cameras)

- **ODM option.** `--radiometric-calibration` takes `none` (default; the output stays in digital numbers), `camera` (black level, vignetting, row gradient and exposure compensation from the EXIF tags) or `camera+sun` (experimental: also the downwelling light sensor, accounting for sun angle) (ODM, radiometric-calibration). For reflectance, `none` is wrong.
- **Reflectance panels.** The ODM documentation does not describe calibration with reflectance panels. Photograph a calibrated reflectance panel before and after each flight anyway (*rule of thumb*: panel level, unshaded, filling a good part of the frame) so that an empirical-line correction (Smith & Milton 1999) can be applied in software that supports panels, or checked by hand against the panel's known reflectance.
- **Keep radiometry intact.** Use `--texturing-skip-global-seam-leveling` (ODM: useful for radiometric data), because global seam levelling alters pixel values to hide seams.
- **Bands.** Process all bands together, not one folder per band (ODM, Multispectral). `--primary-band` picks the band used for reconstruction (sharp and in focus); use `--skip-band-alignment` only when the bands are already aligned (ODM, arguments). Supported multispectral cameras include the MicaSense RedEdge-MX and Altum, Sentera 6X, DJI Phantom 4 Multispectral and DJI Mavic 3 Multispectral (ODM, Multispectral).
- **Check the output.** Run `claw-ortho-qc --expect reflectance`: values should lie between 0 and 1. Note the band order ODM wrote, and pass it to `claw-vi-extract --bands`.
- Aasen et al. (2018) review sensor technology, measurement procedures and data correction for quantitative UAV spectral data; Tmušić et al. (2020) summarise current practice for UAS environmental monitoring.

## 4. Orthomosaic and surface-model settings (ODM)

| Option | What it does (ODM documentation) | Advice for plot trials |
|---|---|---|
| `--orthophoto-resolution` | orthophoto cm/pixel, default 5, capped by the estimated GSD | set to the GSD |
| `--dsm`, `--dtm` | build a surface model (ground plus objects) and a terrain model (ground only); default off | turn on for plant height |
| `--dem-resolution` | DSM/DTM cm/pixel, default 5, capped by the GSD | match the GSD for height work |
| `--pc-quality`, `--feature-quality` | point-cloud and feature quality (defaults medium, high); higher costs time and memory | raise for height work |
| `--fast-orthophoto` | skips dense reconstruction; orthophoto from the sparse model | quick checks only; no usable DSM |
| `--auto-boundary` | limits reconstruction to the area under the camera positions | removes distant background |
| `--gcp` | ground control file | see section 2 |

## 5. Plant height from surface models

- Height = DSM - DTM, where the DTM is bare ground: a flight before emergence or after harvest, or ground points in the alleys. Use the same GCPs for both so that their vertical datums match; `claw-ortho-qc --dsm --dtm` flags negative or implausibly tall plot heights.
- Crop surface models from UAV RGB images have been used for barley biomass (Bendig et al. 2014) and for wheat height and growth rate in plot trials (Holman et al. 2016); validate against ruler measurements in a subset of plots.
- Wind and mixed pixels at canopy edges bias plot-mean heights low; upper percentiles (for example the 95th) are less affected (*rule of thumb*).

## 6. Trade-offs (rules of thumb)

| Choice | Gain | Cost |
|---|---|---|
| Lower flying height | finer pixels | more images, longer flights, more processing, more light change during the flight |
| Higher overlap | more robust matching, better DSMs | more images and processing time |
| More GCPs | better and verifiable accuracy | field time; markers must stay visible |
| Higher ODM quality settings | better point clouds and DSMs | time and memory |

## 7. Before extracting plot values

1. `claw-ortho-qc` with the plot grid (and DSM/DTM): no failed checks; warnings understood.
2. Overlay the plot grid on the mosaic and confirm it sits on the plots on every date.
3. Note the band order, the reflectance scale and the flight date for `claw-vi-extract` or `claw-canopy-cover`.

## References

- OpenDroneMap documentation: Flying tips, https://docs.opendronemap.org/flying/; Ground control points, https://docs.opendronemap.org/gcp/; Multispectral and thermal, https://docs.opendronemap.org/multispectral/; Options and flags, https://docs.opendronemap.org/arguments/ (radiometric calibration: https://docs.opendronemap.org/arguments/radiometric-calibration/). Accessed 2026-10-15.
- Assmann JJ, Kerby JT, Cunliffe AM, Myers-Smith IH (2019) Vegetation monitoring using multispectral sensors - best practices and lessons learned from high latitudes. Journal of Unmanned Vehicle Systems 7:54-75. https://doi.org/10.1139/juvs-2018-0018
- Aasen H, Honkavaara E, Lucieer A, Zarco-Tejada PJ (2018) Quantitative remote sensing at ultra-high resolution with UAV spectroscopy: a review of sensor technology, measurement procedures, and data correction workflows. Remote Sensing 10:1091. https://doi.org/10.3390/rs10071091
- Tmušić G, Manfreda S, Aasen H, et al. (2020) Current practices in UAS-based environmental monitoring. Remote Sensing 12:1001. https://doi.org/10.3390/rs12061001
- Smith GM, Milton EJ (1999) The use of the empirical line method to calibrate remotely sensed data to reflectance. International Journal of Remote Sensing 20:2653-2662. https://doi.org/10.1080/014311699211994
- Sanz-Ablanedo E, Chandler JH, Rodríguez-Pérez JR, Ordóñez C (2018) Accuracy of unmanned aerial vehicle (UAV) and SfM photogrammetry survey as a function of the number and location of ground control points used. Remote Sensing 10:1606. https://doi.org/10.3390/rs10101606
- James MR, Chandler JH, Eltner A, et al. (2019) Guidelines on the use of structure-from-motion photogrammetry in geomorphic research. Earth Surface Processes and Landforms 44:2081-2084. https://doi.org/10.1002/esp.4637
- Bendig J, Bolten A, Bennertz S, Broscheit J, Eichfuss S, Bareth G (2014) Estimating biomass of barley using crop surface models (CSMs) derived from UAV-based RGB imaging. Remote Sensing 6:10395-10412. https://doi.org/10.3390/rs61110395
- Holman F, Riche A, Michalski A, Castle M, Wooster M, Hawkesford M (2016) High throughput field phenotyping of wheat plant height and growth rate in field plot trials using UAV based remote sensing. Remote Sensing 8:1031. https://doi.org/10.3390/rs8121031
