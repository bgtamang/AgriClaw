# Synthetic demo inputs for claw-canopy-cover (fixed seed). The main script calls
# make_demo() for --demo and writes the files to <out>/demo_input; run this file
# directly (Rscript demo/make_demo.R [DIR]) to write them to DIR (default: a
# temporary folder). Nothing here comes from real imagery.
#
#   rgb_ortho.tif  8-bit RGB orthomosaic, 2 cm pixels, UTM 15N, about 5.5 x 9.5 m
#   plots.gpkg     9 plots (3 columns x 3 ranges, 1.5 x 2.5 m, 0.5 m alleys) with
#                  plot_id, plot and the true canopy cover of the synthetic scene
#                  (true_cover_pct, exact area fraction from 4 x 4 supersampling)
#
# Plants are green discs along two rows per plot (0.75 m apart) whose size differs
# between plots, on brown soil, with a brightness gradient across the field, pixel
# noise and anti-aliased leaf edges. Plot P07 has straw residue patches, P03 a
# small specular (saturated) patch and P05 a dark shadow strip.

# Discs painted into a plot-sized matrix (xl, yl: its sub-pixel centres); only
# each disc's bounding box is scanned. Vegetation outside plots is cut anyway.
paint_discs <- function(xl, yl, cx, cy, r) {
  loc <- matrix(0, length(yl), length(xl))
  for (j in seq_along(cx)) {
    if (r[j] <= 0) next
    jx <- which(xl >= cx[j] - r[j] & xl <= cx[j] + r[j])
    jy <- which(yl >= cy[j] - r[j] & yl <= cy[j] + r[j])
    if (!length(jx) || !length(jy)) next
    d2 <- outer((yl[jy] - cy[j])^2, (xl[jx] - cx[j])^2, "+")
    blk <- loc[jy, jx, drop = FALSE]
    blk[d2 <= r[j]^2] <- 1
    loc[jy, jx] <- blk
  }
  loc
}

synth_rgb_scene <- function(seed = 1L, gsd = 0.02, x0 = 500000, y0 = 4570000, artifacts = TRUE,
                            radius = c(0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.22, 0.27, 0.34)) {
  set.seed(seed)
  n_cols <- 3L; n_ranges <- 3L; W <- 1.5; L <- 2.5; A <- 0.5; margin <- 0.5
  width <- n_cols * W + 2 * margin
  height <- n_ranges * L + (n_ranges - 1) * A + 2 * margin
  nc <- round(width / gsd); nr <- round(height / gsd)
  ss <- 4L                                          # supersampling per pixel side
  sub <- (seq_len(ss) - 0.5) / ss
  xs <- as.vector(outer(sub, 0:(nc - 1), "+")) * gsd           # sub-pixel x centres
  ys <- height - as.vector(outer(sub, 0:(nr - 1), "+")) * gsd  # top to bottom
  vegf <- matrix(0, nr * ss, nc * ss)              # vegetation indicator at sub-pixels
  resid <- matrix(0, nr * ss, nc * ss)
  radius <- radius[sample(9)]
  plots <- list(); info <- data.frame(plot = 1:9, plot_id = sprintf("P%02d", 1:9))
  k <- 0L
  for (cl in seq_len(n_cols)) for (rg in seq_len(n_ranges)) {
    k <- k + 1L
    xa <- margin + (cl - 1) * W; ya <- margin + (rg - 1) * (L + A)
    plots[[k]] <- c(xa, ya, xa + W, ya + L)
    centres <- expand.grid(x = xa + c(0.375, 1.125), y = seq(ya + 0.1, ya + L - 0.1, by = 0.15))
    centres$x <- centres$x + stats::rnorm(nrow(centres), 0, 0.02)
    rr <- radius[k] * stats::runif(nrow(centres), 0.85, 1.15)
    ix <- which(xs >= xa & xs <= xa + W)
    iy <- which(ys >= ya & ys <= ya + L)
    vegf[iy, ix] <- paint_discs(xs[ix], ys[iy], centres$x, centres$y, rr)
    if (artifacts && k == 7L) {                     # straw residue patches between the rows
      cx <- xa + 0.75 + stats::runif(6, -0.15, 0.15); cy <- ya + 0.3 + (0:5) * 0.38
      resid[iy, ix] <- paint_discs(xs[ix], ys[iy], cx, cy, rep(0.12, 6))
    }
  }
  # Clip vegetation to plots (alleys and margins are bare soil).
  inplot <- matrix(FALSE, nr * ss, nc * ss)
  for (p in plots) inplot[ys >= p[2] & ys <= p[4], xs >= p[1] & xs <= p[3]] <- TRUE
  vegf[!inplot] <- 0
  resid[vegf == 1] <- 0
  agg <- function(m) {                              # sub-pixels -> pixel fractions (block means)
    a <- rowsum(m, rep(seq_len(nr), each = ss), reorder = FALSE)
    t(rowsum(t(a), rep(seq_len(nc), each = ss), reorder = FALSE)) / ss^2
  }
  fv <- agg(vegf); fr <- agg(resid)
  soil <- c(140, 110, 85); leaf <- c(55, 105, 38); straw <- c(196, 178, 128)
  gx <- matrix(seq(0.85, 1.15, length.out = nc), nr, nc, byrow = TRUE)     # illumination gradient
  shade <- matrix(stats::runif(nr * nc, 0.9, 1.1), nr, nc)
  bands <- lapply(1:3, function(b) {
    v <- (fv * leaf[b] + fr * straw[b] + (1 - fv - fr) * soil[b]) * gx * shade
    v + stats::rnorm(nr * nc, 0, 4)
  })
  # P05: shadow strip (dark); P03: specular patch (saturated)
  px <- function(x) floor(x / gsd) + 1; py <- function(y) floor((height - y) / gsd) + 1
  if (artifacts) {
    p5 <- plots[[5]]; rows5 <- py(p5[2] + 1.2):py(p5[2] + 0.9); cols5 <- px(p5[1]):(px(p5[3]) - 1)
    for (b in 1:3) bands[[b]][rows5, cols5] <- bands[[b]][rows5, cols5] * 0.08
    p3 <- plots[[3]]; rows3 <- py(p3[2] + 1.0):py(p3[2] + 0.9); cols3 <- px(p3[1] + 0.3):px(p3[1] + 0.5)
    for (b in 1:3) bands[[b]][rows3, cols3] <- 255
  }
  bands <- lapply(bands, function(v) pmin(pmax(round(v), 0), 255))
  r <- terra::rast(lapply(bands, function(m) terra::rast(m, extent = terra::ext(x0, x0 + nc * gsd, y0, y0 + nr * gsd),
                                                         crs = "EPSG:32615")))
  names(r) <- c("red", "green", "blue")
  truth <- vapply(plots, function(p) {
    sel <- inplot[ys >= p[2] & ys <= p[4], xs >= p[1] & xs <= p[3]]
    100 * mean(vegf[ys >= p[2] & ys <= p[4], xs >= p[1] & xs <= p[3]][sel])
  }, numeric(1))
  info$true_cover_pct <- round(truth, 2)
  geoms <- lapply(plots, function(p) sf::st_polygon(list(rbind(c(p[1], p[2]), c(p[3], p[2]), c(p[3], p[4]),
                                                               c(p[1], p[4]), c(p[1], p[2])) + rep(c(x0, y0), each = 5))))
  list(raster = r, plots = sf::st_sf(info, geometry = sf::st_sfc(geoms, crs = 32615)), veg_fraction = fv)
}

make_demo <- function(dir, seed = 1L) {
  sc <- synth_rgb_scene(seed)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  rpath <- file.path(dir, "rgb_ortho.tif")
  # NAflag = NA: terra would otherwise declare 255 as NoData and hide saturated pixels.
  terra::writeRaster(sc$raster, rpath, datatype = "INT1U", overwrite = TRUE, NAflag = NA,
                     gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  ppath <- file.path(dir, "plots.gpkg")
  if (file.exists(ppath)) unlink(ppath)
  sf::st_write(sc$plots, ppath, layer = "plots", quiet = TRUE)
  invisible(list(raster = rpath, plots = ppath))
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  d <- if (length(a)) a[[1]] else file.path(tempdir(), "claw-canopy-cover-demo")
  p <- make_demo(d)
  cat("Wrote", p$raster, "and", p$plots, "\n")
}
