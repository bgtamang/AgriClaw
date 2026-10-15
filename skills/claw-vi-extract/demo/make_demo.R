# Synthetic demo inputs for claw-vi-extract (fixed seed). The main script calls
# make_demo() for --demo and writes the files to <out>/demo_input; run this file
# directly (Rscript demo/make_demo.R [DIR]) to write them to DIR (default: a
# temporary folder). Nothing here comes from real imagery.
#
#   ortho_multispectral.tif  5-band surface reflectance stored as 16-bit integers
#                            (reflectance x 10000, so the demo runs with --scale 0.0001),
#                            layers named blue, green, red, rededge, nir in wavelength
#                            order like a stacked MicaSense export, 5 cm pixels, UTM 15N
#   plots.gpkg               12 plots (3 columns x 4 ranges, 1.5 x 3 m, 0.75 m alleys)
#                            with plot, plot_id, rep and genotype
#
# Each plot has two crop rows whose canopy width differs between plots (cover
# 25-95%) and its own canopy spectrum; soil fills the inter-row gaps and alleys.

make_demo <- function(dir, seed = 42L) {
  set.seed(seed)
  gsd <- 0.05
  x0 <- 500000; y0 <- 4570000                     # SW corner of the image (UTM 15N)
  n_cols <- 3L; n_ranges <- 4L; W <- 1.5; L <- 3; A <- 0.75; margin <- 1
  width <- n_cols * W + 2 * margin
  height <- n_ranges * L + (n_ranges - 1) * A + 2 * margin
  nc <- round(width / gsd); nr <- round(height / gsd)
  xs <- x0 + (seq_len(nc) - 0.5) * gsd
  ys <- y0 + height - (seq_len(nr) - 0.5) * gsd     # row 1 is the top of the image
  X <- matrix(xs, nr, nc, byrow = TRUE)
  Y <- matrix(ys, nr, nc)

  soil <- c(blue = 0.09, green = 0.12, red = 0.16, rededge = 0.19, nir = 0.23)
  bands <- lapply(names(soil), function(b) matrix(soil[[b]], nr, nc))
  names(bands) <- names(soil)
  n_plots <- n_cols * n_ranges
  half_w <- seq(0.09, 0.36, length.out = n_plots)[sample(n_plots)]   # canopy half-width per row (m)
  chl <- stats::runif(n_plots)                                       # 0 = pale, 1 = dark green
  plots <- vector("list", n_plots)
  info <- data.frame(plot = seq_len(n_plots), rep = 0L, genotype = "", stringsAsFactors = FALSE)
  k <- 0L
  for (cl in seq_len(n_cols)) for (rg in seq_len(n_ranges)) {
    k <- k + 1L
    xa <- x0 + margin + (cl - 1) * W; xb <- xa + W
    ya <- y0 + margin + (rg - 1) * (L + A); yb <- ya + L
    plots[[k]] <- sf::st_polygon(list(rbind(c(xa, ya), c(xb, ya), c(xb, yb), c(xa, yb), c(xa, ya))))
    rows <- xa + c(0.375, 1.125)
    inside <- X > xa & X < xb & Y > ya & Y < yb
    veg <- inside & (abs(X - rows[1]) < half_w[k] | abs(X - rows[2]) < half_w[k])
    spec <- c(blue = 0.035 - 0.01 * chl[k], green = 0.09 - 0.03 * chl[k], red = 0.06 - 0.03 * chl[k],
              rededge = 0.28 - 0.10 * chl[k], nir = 0.40 + 0.10 * chl[k])
    for (b in names(spec)) bands[[b]][veg] <- spec[[b]]
    info$rep[k] <- cl
    info$genotype[k] <- sprintf("G%d", rg)
  }
  for (b in names(bands)) {
    bands[[b]] <- bands[[b]] * (1 + stats::rnorm(nr * nc, 0, 0.03))
    bands[[b]] <- round(pmax(bands[[b]], 0.001), 4)
  }
  r <- terra::rast(lapply(names(bands), function(b)
    terra::rast(bands[[b]], extent = terra::ext(x0, x0 + nc * gsd, y0, y0 + nr * gsd), crs = "EPSG:32615")))
  names(r) <- names(bands)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  rpath <- file.path(dir, "ortho_multispectral.tif")
  terra::writeRaster(round(r * 10000), rpath, datatype = "INT2U", overwrite = TRUE,
                     gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  info$plot_id <- sprintf("P%02d", info$plot)
  p <- sf::st_sf(info[, c("plot_id", "plot", "rep", "genotype")], geometry = sf::st_sfc(plots, crs = 32615))
  ppath <- file.path(dir, "plots.gpkg")
  if (file.exists(ppath)) unlink(ppath)
  sf::st_write(p, ppath, layer = "plots", quiet = TRUE)
  invisible(list(raster = rpath, plots = ppath))
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  d <- if (length(a)) a[[1]] else file.path(tempdir(), "claw-vi-extract-demo")
  p <- make_demo(d)
  cat("Wrote", p$raster, "and", p$plots, "\n")
}
