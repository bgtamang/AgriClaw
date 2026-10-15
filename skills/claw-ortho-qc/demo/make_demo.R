# Synthetic demo inputs for claw-ortho-qc (fixed seed). The main script calls
# make_demo() for --demo and writes the files to <out>/demo_input; run this file
# directly (Rscript demo/make_demo.R [DIR]) to write them to DIR (default: a
# temporary folder). Nothing here comes from real imagery.
#
#   ortho_reflectance.tif  4-band reflectance (blue, green, red, NIR) x 10000 as
#                          16-bit integers, NoData = 0, 5 cm pixels, UTM 15N
#   plots.gpkg             12 plots (3 columns x 4 ranges, 1.5 x 2 m)
#   dsm.tif                surface model (float32, 5 cm): ground about 300 m with a
#                          gentle slope plus 0.2-0.9 m of crop canopy
#   dtm.tif                bare-ground terrain model (float32, 10 cm)
#
# Deliberate problems, so the demo report shows warnings: a NoData hole over 10%
# of plot P06, plot P12 reaching past the edge of the image, and 0.3% of DSM
# pixels inside the plots replaced by spikes of +3 m.

make_demo <- function(dir, seed = 7L) {
  set.seed(seed)
  gsd <- 0.05
  x0 <- 500000; y0 <- 4570000
  n_cols <- 3L; n_ranges <- 4L; W <- 1.5; L <- 2; A <- 0.5; margin <- 0.5
  width <- n_cols * W + 2 * margin
  height <- n_ranges * L + (n_ranges - 1) * A + 2 * margin
  nc <- round(width / gsd); nr <- round(height / gsd)
  ext <- terra::ext(x0, x0 + nc * gsd, y0, y0 + nr * gsd)
  X <- matrix(x0 + (seq_len(nc) - 0.5) * gsd, nr, nc, byrow = TRUE)
  Y <- matrix(y0 + nr * gsd - (seq_len(nr) - 0.5) * gsd, nr, nc)
  soil <- c(0.09, 0.12, 0.16, 0.23); leaf <- c(0.03, 0.08, 0.04, 0.45)
  vegm <- matrix(FALSE, nr, nc)
  height_m <- matrix(0, nr, nc)
  plots <- list(); k <- 0L
  for (cl in seq_len(n_cols)) for (rg in seq_len(n_ranges)) {
    k <- k + 1L
    xa <- x0 + margin + (cl - 1) * W; ya <- y0 + margin + (rg - 1) * (L + A)
    xb <- xa + W; yb <- ya + L
    if (k == 12L) xb <- x0 + nc * gsd + 0.4                       # P12 runs off the image edge
    plots[[k]] <- sf::st_polygon(list(rbind(c(xa, ya), c(xb, ya), c(xb, yb), c(xa, yb), c(xa, ya))))
    inside <- X > xa & X < xa + W & Y > ya & Y < yb
    v <- inside & (abs(X - (xa + 0.375)) < 0.26 | abs(X - (xa + 1.125)) < 0.26)
    vegm <- vegm | v
    height_m[v] <- stats::runif(1, 0.2, 0.9)
  }
  bands <- lapply(1:4, function(b) {
    m <- ifelse(vegm, leaf[b], soil[b]) * (1 + stats::rnorm(nr * nc, 0, 0.03))
    pmax(round(m * 10000), 1)
  })
  # NoData hole over 10% of plot P06 (column 2, range 2)
  p6 <- sf::st_bbox(plots[[6]])
  hole <- X > p6[["xmin"]] & X < p6[["xmin"]] + 0.15 & Y > p6[["ymin"]] & Y < p6[["ymax"]]
  for (b in 1:4) bands[[b]][hole] <- 0
  ortho <- terra::rast(lapply(bands, function(m) terra::rast(m, extent = ext, crs = "EPSG:32615")))
  names(ortho) <- c("blue", "green", "red", "nir")

  ground <- 300 + 0.02 * (X - x0) + 0.01 * (Y - y0)
  dsm <- ground + height_m + stats::rnorm(nr * nc, 0, 0.01)
  inplots <- Reduce(`|`, lapply(plots[1:11], function(p) {
    b <- sf::st_bbox(p); X > b[["xmin"]] & X < b[["xmax"]] & Y > b[["ymin"]] & Y < b[["ymax"]]
  }))
  spikes <- sample(which(inplots), round(0.003 * sum(inplots)))
  dsm[spikes] <- dsm[spikes] + 3
  dsm_r <- terra::rast(dsm, extent = ext, crs = "EPSG:32615")
  dtm_r <- terra::aggregate(terra::rast(ground, extent = ext, crs = "EPSG:32615"), 2, fun = "mean")

  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  paths <- list(raster = file.path(dir, "ortho_reflectance.tif"), plots = file.path(dir, "plots.gpkg"),
                dsm = file.path(dir, "dsm.tif"), dtm = file.path(dir, "dtm.tif"))
  terra::writeRaster(ortho, paths$raster, datatype = "INT2U", NAflag = 0, overwrite = TRUE,
                     gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  terra::writeRaster(dsm_r, paths$dsm, datatype = "FLT4S", overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
  terra::writeRaster(dtm_r, paths$dtm, datatype = "FLT4S", overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
  p <- sf::st_sf(plot_id = sprintf("P%02d", seq_along(plots)), plot = seq_along(plots),
                 geometry = sf::st_sfc(plots, crs = 32615))
  if (file.exists(paths$plots)) unlink(paths$plots)
  sf::st_write(p, paths$plots, layer = "plots", quiet = TRUE)
  invisible(paths)
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  d <- if (length(a)) a[[1]] else file.path(tempdir(), "claw-ortho-qc-demo")
  p <- make_demo(d)
  cat("Wrote demo inputs to", d, "\n")
}
