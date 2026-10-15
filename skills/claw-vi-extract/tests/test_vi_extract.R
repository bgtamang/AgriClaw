# Tests for claw-vi-extract. Run from anywhere: Rscript tests/test_vi_extract.R
# Known answers: hand-computed index values, synthetic rasters with known
# per-plot reflectance and exact soil fractions, coverage weights of plots that
# cut pixels in half, Sentinel-2 scaling with the -1000 offset, DN refusal, a
# hand-checked Otsu example, and an end-to-end demo run.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "vi_extract.R"))  # loads functions; main() does not run

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
expect_error <- function(expr, pattern, what) {
  msg <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  check(!is.na(msg) && grepl(pattern, msg, ignore.case = TRUE), paste0(what, " (got: ", msg, ")"))
}
run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "vi_extract.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
tmp <- file.path(tempdir(), "claw-vi-extract-tests")
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
CRS <- "EPSG:32616"
X0 <- 500000; Y0 <- 4400000                       # local test coordinates sit on the zone's central meridian
sq <- function(xa, ya, xb, yb, x0 = X0, y0 = Y0) {
  xa <- xa + x0; xb <- xb + x0; ya <- ya + y0; yb <- yb + y0
  sf::st_polygon(list(rbind(c(xa, ya), c(xb, ya), c(xb, yb), c(xa, yb), c(xa, ya))))
}
write_plots <- function(geoms, name, crs = 32616, extra = NULL) {
  df <- data.frame(plot_id = sprintf("P%d", seq_along(geoms)), stringsAsFactors = FALSE)
  if (!is.null(extra)) df <- cbind(df, extra)
  p <- sf::st_sf(df, geometry = sf::st_sfc(geoms, crs = crs))
  path <- file.path(tmp, paste0(name, ".gpkg"))
  if (file.exists(path)) unlink(path)
  sf::st_write(p, path, quiet = TRUE)
  path
}
# Raster from a list of band matrices (rows = top to bottom).
write_raster <- function(bands, res, name, crs = CRS, datatype = "FLT8S", x0 = X0, y0 = Y0, nms = names(bands)) {
  nr <- nrow(bands[[1]]); nc <- ncol(bands[[1]])
  r <- terra::rast(lapply(bands, function(m) terra::rast(m, extent = terra::ext(x0, x0 + nc * res, y0, y0 + nr * res), crs = crs)))
  names(r) <- nms
  path <- file.path(tmp, paste0(name, ".tif"))
  terra::writeRaster(r, path, datatype = datatype, overwrite = TRUE)
  path
}
opt_for <- function(input, plots, ...) {
  o <- DEFAULTS
  o$input <- input; o$plots <- plots; o$bands <- "blue=1,green=2,red=3,rededge=4,nir=5"
  mods <- list(...)
  for (k in names(mods)) o[[k]] <- mods[[k]]
  o
}
roles <- c("blue", "green", "red", "rededge", "nir")
TOL <- 1e-6   # terra stores computed rasters as float32 (about 7 significant digits)

# 1. Index formulas against values computed by hand (independently, in Python)
b <- list(blue = 0.04, green = 0.08, red = 0.05, rededge = 0.20, nir = 0.50)
im <- index_matrix(b, names(INDEX_DEFS))
hand <- c(ndvi = 0.8181818182, ndre = 0.4285714286, gndvi = 0.7241379310, osavi = 0.6338028169,
          evi2 = 0.6944444444, mcari2 = 0.6827724716, cire = 1.5)
for (ix in names(hand)) check(abs(im[1, ix] - hand[[ix]]) < 1e-9, paste("hand-computed", ix))
z <- index_matrix(list(blue = 0, green = 0, red = 0, rededge = 0, nir = 0), names(INDEX_DEFS))
check(all(is.na(z[1, c("ndvi", "ndre", "gndvi", "cire")])), "0/0 and x/0 give NA, not Inf")

# 2. Otsu: hand-checked histogram, plateau midpoint, symmetric mixture, brute force
check(abs(otsu_from_hist(c(3, 5, 2, 1, 1, 6, 3), 0:6) - 3.5) < 1e-12, "Otsu on the hand-checked histogram is 3.5")
check(abs(otsu_from_hist(c(5, 5, 0, 0, 0, 5, 5), 0:6) - 3) < 1e-12, "Otsu takes the middle of an empty gap")
set.seed(3)
check(abs(attr(otsu_threshold(rnorm(2e5)), "eta") - 2 / pi) < 0.01, "Otsu separability of one normal mode is 2/pi")
x <- c(rnorm(5000, 0.2, 0.05), rnorm(5000, 0.8, 0.05))
check(abs(otsu_threshold(x) - 0.5) < 0.02, sprintf("Otsu on a symmetric mixture is 0.5 (got %.3f)", otsu_threshold(x)))
y <- c(rnorm(3000, 0.15, 0.06), rnorm(7000, 0.7, 0.1))
t_o <- otsu_threshold(y, 256)
cand <- seq(min(y), max(y), length.out = 2001)
sbv <- vapply(cand, function(t) { a <- y <= t; if (all(a) || !any(a)) return(NA_real_)
  mean(a) * mean(!a) * (mean(y[a]) - mean(y[!a]))^2 }, numeric(1))
check(abs(t_o - cand[which.max(sbv)]) < 2 * diff(range(y)) / 256,
      "Otsu threshold maximizes the between-class variance of the raw values")

# 3. Known per-plot reflectance, plots aligned to the pixel grid
res <- 0.1
spectra <- rbind(c(0.03, 0.07, 0.04, 0.22, 0.45), c(0.05, 0.10, 0.09, 0.20, 0.30), c(0.02, 0.06, 0.03, 0.30, 0.55),
                 c(0.04, 0.09, 0.06, 0.25, 0.40), c(0.06, 0.11, 0.12, 0.18, 0.25), c(0.03, 0.08, 0.05, 0.21, 0.48))
colnames(spectra) <- roles
bands <- lapply(roles, function(b) matrix(0.1, 40, 60))
names(bands) <- roles
geoms <- list(); k <- 0
for (rr in 1:2) for (cc in 1:3) {
  k <- k + 1
  rows <- if (rr == 1) 21:40 else 1:20          # rr 1 = bottom half
  cols <- (cc - 1) * 20 + 1:20
  for (b in roles) bands[[b]][rows, cols] <- spectra[k, b]
  ya <- if (rr == 1) 0 else 2
  geoms[[k]] <- sq((cc - 1) * 2, ya, cc * 2, ya + 2)
}
r1 <- write_raster(bands, res, "uniform")
p1 <- write_plots(geoms, "uniform")
out1 <- run_extraction(opt_for(r1, p1))
truth <- index_matrix(lapply(stats::setNames(roles, roles), function(b) spectra[, b]), names(INDEX_DEFS))
for (ix in names(INDEX_DEFS))
  check(max(abs(out1$wide[[ix]] - truth[, ix])) < TOL, paste("uniform plots: mean", ix, "equals the truth"))
lg <- out1$long
check(max(abs(lg$median - lg$mean), abs(lg$p10 - lg$mean), abs(lg$p90 - lg$mean)) < TOL && max(lg$sd) < TOL,
      "uniform plots: median, p10, p90 equal the mean and sd is 0")
check(max(abs(out1$wide$n_pixels - 400)) < TOL && all(out1$wide$n_pure_pixels == 400), "each 2 x 2 m plot has 400 whole 10-cm pixels")
check(all(abs(out1$bands$nir_reflectance - spectra[, "nir"]) < TOL), "band means equal the true reflectance")
check(!any(grepl("Band order", out1$warnings)), "explicit --bands gives no band-order warning")

# 4. Soil pixels in known fractions: unmasked, NDVI-threshold and Otsu masks
veg <- c(0.03, 0.08, 0.04, 0.25, 0.50); soil <- c(0.09, 0.12, 0.16, 0.19, 0.23)
names(veg) <- names(soil) <- roles
fr <- c(0, 0.25, 0.5, 0.75, 0.1, 0.9)            # soil fraction per plot
bands <- lapply(roles, function(b) matrix(soil[[b]], 40, 60)); names(bands) <- roles
k <- 0
for (rr in 1:2) for (cc in 1:3) {
  k <- k + 1
  rows <- if (rr == 1) 21:40 else 1:20
  cols <- (cc - 1) * 20 + 1:20
  vcols <- cols[seq_len(20) > round(20 * fr[k])]
  for (b in roles) bands[[b]][rows, vcols] <- veg[[b]]
}
r2 <- write_raster(bands, res, "soil")
nd_v <- (veg[["nir"]] - veg[["red"]]) / (veg[["nir"]] + veg[["red"]])
nd_s <- (soil[["nir"]] - soil[["red"]]) / (soil[["nir"]] + soil[["red"]])
o_none <- run_extraction(opt_for(r2, p1))
check(max(abs(o_none$wide$ndvi - (fr * nd_s + (1 - fr) * nd_v))) < TOL, "no mask: NDVI mean = soil/vegetation mixture")
o_thr <- run_extraction(opt_for(r2, p1, mask = "ndvi", mask_threshold = 0.5))
check(max(abs(o_thr$wide$ndvi - nd_v)) < TOL, "NDVI-threshold mask: plot NDVI = vegetation NDVI")
check(max(abs(o_thr$wide$veg_fraction - (1 - fr))) < TOL, "veg_fraction = true vegetation fraction")
check(max(abs(o_thr$wide$n_pixels - 400 * (1 - fr))) < TOL && max(abs(o_thr$wide$n_pixels_all - 400)) < TOL,
      "masked pixel counts")
check(max(abs(o_thr$bands$red_reflectance - veg[["red"]])) < TOL, "masked band means are the vegetation spectrum")
o_ot <- run_extraction(opt_for(r2, p1, mask = "otsu"))
check(o_ot$threshold > nd_s && o_ot$threshold < nd_v, sprintf("Otsu threshold %.3f lies between soil and vegetation", o_ot$threshold))
check(max(abs(o_ot$wide$ndvi - nd_v)) < TOL, "Otsu mask: plot NDVI = vegetation NDVI")
o_sh <- run_extraction(opt_for(r2, p1, shadow_nir = 0.3))
check(max(abs(o_sh$wide$ndvi - nd_v)) < TOL, "shadow mask on NIR removes the darker soil pixels")
expect_error(run_extraction(opt_for(r2, p1, mask = "ndvi")), "mask-threshold", "--mask ndvi needs a threshold")

# 5. Coverage weighting when plot edges cut pixels in half
nirm <- matrix(rep(0.3 + 0.01 * (1:10), each = 10), 10, 10)
bands <- list(blue = matrix(0.04, 10, 10), green = matrix(0.08, 10, 10), red = matrix(0.05, 10, 10),
              rededge = matrix(0.2, 10, 10), nir = nirm)
r3 <- write_raster(bands, 0.1, "gradient")
p3 <- write_plots(list(sq(0.05, 0.2, 0.35, 0.8)), "halfpix")
o3 <- run_extraction(opt_for(r3, p3))
w <- c(0.5, 1, 1, 0.5); nv <- 0.3 + 0.01 * (1:4)
check(abs(o3$bands$nir_reflectance - sum(w * nv) / sum(w)) < TOL, "coverage-weighted band mean (0.325)")
nd <- (nv - 0.05) / (nv + 0.05)
check(abs(o3$wide$ndvi - sum(w * nd) / sum(w)) < TOL, "coverage-weighted NDVI mean")
check(abs(o3$wide$n_pixels - 18) < 1e-9 && o3$wide$n_pure_pixels == 12, "effective (18) and whole (12) pixel counts")

# 6. Buffers in meters, and plots in another CRS (lon/lat)
o4 <- run_extraction(opt_for(r1, p1, buffer = -0.5))
check(all(abs(o4$wide$plot_area_m2 - 1) < 1e-9) && all(abs(o4$wide$n_pixels - 100) < 1e-9), "-0.5 m buffer leaves 1 x 1 m (100 pixels)")
expect_error(run_extraction(opt_for(r1, p1, buffer = -1.1)), "removes", "a buffer that deletes plots is refused")
pll <- sf::st_transform(sf::st_read(p1, quiet = TRUE), 4326)
p1ll <- file.path(tmp, "uniform_ll.gpkg"); if (file.exists(p1ll)) unlink(p1ll)
sf::st_write(pll, p1ll, quiet = TRUE)
o5 <- run_extraction(opt_for(r1, p1ll, buffer = -0.5))
check(max(abs(o5$wide$ndvi - o4$wide$ndvi)) < 1e-9 && max(abs(o5$wide$n_pixels - 100)) < 0.05,
      "lon/lat plots are buffered in meters and give the same values")

# 7. Sentinel-2 L2A scaling with the -1000 offset; DN 0 = no data; mixed pixels
refl <- c(B02 = 0.0312, B03 = 0.0655, B04 = 0.0421, B05 = 0.1288, B08 = 0.4106, B11 = 0.2000)
bands <- lapply(refl, function(v) matrix(round(v * 10000) + 1000, 10, 10))
bands$B04[1, 1] <- 0                               # one no-data pixel (outside the plots below)
r5 <- write_raster(bands, 10, "s2", datatype = "INT2U", x0 = 600000, y0 = 4400000)
p5 <- write_plots(list(sq(600010, 4400010, 600040, 4400040, 0, 0), sq(600050, 4400050, 600090, 4400090, 0, 0)), "s2plots")
o6 <- run_extraction(opt_for(r5, p5, bands = "", sensor = "sentinel2-l2a"))
check(max(abs(o6$bands$red_reflectance - refl[["B04"]])) < TOL && max(abs(o6$bands$nir_reflectance - refl[["B08"]])) < TOL,
      "Sentinel-2 L2A: (DN - 1000) / 10000 recovers reflectance")
check(identical(unname(o6$idx), c(1L, 2L, 3L, 4L, 5L)) && any(grepl("sensor preset", o6$warnings)),
      "preset bands matched by layer name and flagged as unconfirmed")
o7 <- run_extraction(opt_for(r5, p5, bands = "", sensor = "sentinel2-l2a-harmonized"))
check(max(abs(o7$bands$red_reflectance - refl[["B04"]] - 0.1)) < TOL, "without the offset reflectance is 0.1 too high")
check(is.na(terra::values(to_reflectance(terra::rast(r5)[[3]], 0, 1e-4, -0.1))[1]), "DN 0 becomes NA, not -0.1")
p6 <- write_plots(list(sq(600012, 4400012, 600013.5, 4400017, 0, 0), sq(600031, 4400031, 600032.5, 4400036, 0, 0)), "tiny")
o8 <- run_extraction(opt_for(r5, p6, bands = "", sensor = "sentinel2-l2a"))
check(any(grepl("Mixed pixels", o8$warnings)) && all(grepl("few pixels", o8$wide$flag)) &&
        all(grepl("no whole pixel", o8$wide$flag)), "1.5 x 5 m plots on 10 m pixels are flagged as mixed")

# 8. Digital numbers are refused unless the user confirms calibration
set.seed(1)
dn <- lapply(roles, function(b) matrix(round(runif(400, 15000, 45000)), 20, 20)); names(dn) <- roles
r6 <- write_raster(dn, 0.1, "dn", datatype = "INT2U")
p7 <- write_plots(list(sq(0.2, 0.2, 1.8, 1.8)), "dnplot")
expect_error(run_extraction(opt_for(r6, p7)), "digital numbers", "uint16 DN imagery is refused")
o9 <- run_extraction(opt_for(r6, p7, confirm_calibrated = TRUE))
check(any(grepl("confirm-calibrated", o9$warnings)), "--confirm-calibrated proceeds with a warning")
o10 <- run_extraction(opt_for(r6, p7, scale = 1e-5))
check(!any(grepl("digital numbers", o10$warnings)), "a correct --scale clears the DN check")

# 9. Band-order protection
expect_error(run_extraction(opt_for(r1, p1, bands = "blue=1,green=2,red=5,rededge=4,nir=3")), "conflict",
             "mapping contradicts the layer names")
r1x <- write_raster(lapply(roles, function(b) terra::as.matrix(terra::rast(r1)[[b]], wide = TRUE)), res, "unnamed",
                    nms = paste0("band", 1:5))
o11 <- run_extraction(opt_for(r1x, p1, bands = "blue=1,green=2,red=5,rededge=4,nir=3"))
check(any(grepl("swapped", o11$warnings)), "red/NIR swap is detected from negative NDVI")
expect_error(run_extraction(opt_for(r1x, p1, bands = "")), "Band order is unknown", "no band information is refused")
check(identical(infer_bands_from_names(terra::rast(r1)), c(blue = "1", green = "2", red = "3", rededge = "4", nir = "5")),
      "band roles read from layer names")
expect_error(run_extraction(opt_for(r1, p1, bands = "red=3,nir=5", indices = "ndre")), "need band", "NDRE without red edge is refused")

# 10. End to end on the demo data
out <- file.path(tempdir(), "claw-vi-extract-demo-test")
unlink(out, recursive = TRUE)
cli <- run_cli(c("--demo", "--out", shQuote(out)))
check(cli$status == 0L, paste("demo run exits cleanly:", paste(cli$output, collapse = "\n")))
for (f in c("report.md", "tables/plot_indices.csv", "tables/plot_index_stats.csv", "tables/plot_band_reflectance.csv",
            "tables/image_bands.csv", "spatial/plot_indices.gpkg", "figures/index_map.png", "figures/ndvi_histogram.png",
            "figures/plot_means.png", "reproducibility/commands.sh", "reproducibility/sessionInfo.txt"))
  check(file.exists(file.path(out, f)), paste("demo writes", f))
wd <- utils::read.csv(file.path(out, "tables", "plot_indices.csv"))
check(nrow(wd) == 12L && all(c("plot_id", "date", "ndvi", "ndre", "mcari2", "cire", "flag") %in% names(wd)),
      "demo table has 12 plots with date and index columns")
check(all(wd$ndvi > 0.7 & wd$ndvi < 0.95) && all(is.na(wd$flag) | wd$flag == ""), "demo plot NDVI plausible and unflagged")
check(file.size(file.path(out, "demo_input", "ortho_multispectral.tif")) < 500000, "demo raster is under 500 KB")
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
expect_error(main(c("--demo", "--color", "red")), "Unknown option", "unknown option is rejected")
expect_error(main(c("--input", r1, "--out", file.path(tmp, "noplots"))), "--plots", "missing --plots is rejected")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
