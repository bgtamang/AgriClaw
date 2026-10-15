# Tests for claw-ortho-qc. Run from anywhere: Rscript tests/test_ortho_qc.R
# A clean synthetic orthomosaic + DSM + DTM must pass every check; then each
# known problem is injected into a copy and must be flagged by the right check
# with the right status (and, where it has one, the exact value).

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "ortho_qc.R"))  # loads functions; main() does not run

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
  script <- file.path(SKILL_TEST_DIR, "scripts", "ortho_qc.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
tmp <- file.path(tempdir(), "claw-ortho-qc-tests")
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
X0 <- 500000; Y0 <- 4570000; CRS <- "EPSG:32615"
mk <- function(layers, name, res = 0.05, crs = CRS, datatype = "FLT4S", NAflag = NULL, x0 = X0, y0 = Y0, yres = res) {
  nr <- nrow(layers[[1]]); nc <- ncol(layers[[1]])
  r <- terra::rast(lapply(layers, function(m) terra::rast(m, extent = terra::ext(x0, x0 + nc * res, y0, y0 + nr * yres), crs = crs)))
  names(r) <- names(layers)
  path <- file.path(tmp, paste0(name, ".tif"))
  args <- list(r, path, datatype = datatype, overwrite = TRUE)
  if (!is.null(NAflag)) args$NAflag <- NAflag
  do.call(terra::writeRaster, args)
  path
}
sq <- function(xa, ya, xb, yb) sf::st_polygon(list(rbind(c(xa, ya), c(xb, ya), c(xb, yb), c(xa, yb), c(xa, ya)) +
                                                     rep(c(X0, Y0), each = 5)))
mkplots <- function(geoms, name, crs = 32615) {
  p <- sf::st_sf(plot_id = sprintf("P%d", seq_along(geoms)), geometry = sf::st_sfc(geoms, crs = crs))
  path <- file.path(tmp, paste0(name, ".gpkg"))
  if (file.exists(path)) unlink(path)
  sf::st_write(p, path, quiet = TRUE)
  path
}
qc_run <- function(input, ...) {
  o <- DEFAULTS
  o$input <- input
  mods <- list(...)
  for (k in names(mods)) o[[k]] <- mods[[k]]
  run_qc(o)
}
status_of <- function(res, id) res$qc$status[res$qc$check == id]

# Clean 4-band reflectance image (120 x 100 px of 5 cm), 4 plots, DSM and DTM
nr <- 120; nc <- 100
set.seed(1)
refl <- lapply(c(blue = 0.05, green = 0.09, red = 0.08, nir = 0.35), function(v) matrix(v * (1 + rnorm(nr * nc, 0, 0.05)), nr, nc))
clean <- mk(refl, "clean")
plots4 <- mkplots(list(sq(0.5, 0.5, 2, 2.5), sq(2.5, 0.5, 4, 2.5), sq(0.5, 3, 2, 5), sq(2.5, 3, 4, 5)), "plots4")
ground <- matrix(300, nr, nc)
canopy <- matrix(rep(c(0, 0.6), length.out = nc), nr, nc, byrow = TRUE)
dsm <- mk(list(dsm = ground + canopy + matrix(rnorm(nr * nc, 0, 0.005), nr, nc)), "dsm")
dtm <- mk(list(dtm = ground), "dtm")
rc <- qc_run(clean, plots = plots4, dsm = dsm, dtm = dtm, expect = "reflectance", expect_bands = 4)
check(!any(rc$qc$status %in% c("warn", "fail")),
      paste("clean data pass every check; flagged:", paste(rc$qc$check[rc$qc$status %in% c("warn", "fail")], collapse = ", ")))
check(all(c("rotation", "crs", "pixel_shape", "band_count", "datatype", "nodata_declared", "zero_as_nodata", "extent",
            "gsd_vs_plot", "value_range", "plot_coverage", "saturation", "dsm_crs", "dsm_values", "dsm_coverage",
            "dsm_spikes", "dsm_relief", "canopy_height") %in% rc$qc$check), "all checks are reported")
check(all(abs(rc$plot_table$coverage - 1) < 1e-9), "clean plots are fully covered")
check(all(abs(rc$plot_table$chm_p95_m - 0.6) < 0.02), "canopy height p95 = 0.6 m")

# 1. CRS problems, pixel shape, rotation
nocrs <- mk(refl, "nocrs", crs = "")
check(status_of(qc_run(nocrs), "crs") == "fail", "missing CRS fails")
geo <- mk(refl, "geo", res = 5e-7, crs = "EPSG:4326", x0 = -93.5, y0 = 41.3)
rg <- qc_run(geo)
check(status_of(rg, "crs") == "warn", "lon/lat CRS warns")
check(abs(rg$gsd[2] - 5e-7 * 110574) < 1e-6, "pixel size of a lon/lat raster converted to metres")
aniso <- mk(refl, "aniso", res = 0.05, yres = 0.08)
check(status_of(qc_run(aniso), "pixel_shape") == "warn", "non-square pixels warn")
b <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10, crs = CRS, vals = 1:100)
terra::writeRaster(b, file.path(tmp, "base.tif"), overwrite = TRUE)
writeLines(c('<VRTDataset rasterXSize="10" rasterYSize="10">', '  <SRS>EPSG:32615</SRS>',
             '  <GeoTransform>500000, 0.0707, 0.0707, 4570000, 0.0707, -0.0707</GeoTransform>',
             '  <VRTRasterBand dataType="Float32" band="1">',
             '    <SimpleSource><SourceFilename relativeToVRT="1">base.tif</SourceFilename><SourceBand>1</SourceBand></SimpleSource>',
             '  </VRTRasterBand>', '</VRTDataset>'), file.path(tmp, "rotated.vrt"))
check(status_of(qc_run(file.path(tmp, "rotated.vrt")), "rotation") == "fail", "rotated geotransform fails")

# 2. Pixel size against plot size
coarse <- mk(lapply(refl, function(m) m[1:12, 1:10]), "coarse", res = 0.5)
rco <- qc_run(coarse, plots = plots4)
check(status_of(rco, "gsd_vs_plot") == "warn" && grepl("min 3.0 px", rco$qc$value[rco$qc$check == "gsd_vs_plot"]),
      "50 cm pixels on 1.5 m plots (3 px across) warn")
vcoarse <- mk(lapply(refl, function(m) m[1:3, 1:3]), "vcoarse", res = 2)
check(status_of(qc_run(vcoarse, plots = plots4), "gsd_vs_plot") == "fail", "2 m pixels on 1.5 m plots fail")
check(status_of(qc_run(clean, max_gsd_cm = 3), "gsd") == "warn", "--max-gsd-cm 3 flags 5 cm pixels")

# 3. Values: digital numbers, offset applied twice, band count
dn <- mk(lapply(refl, function(m) round(m * 1e5)), "dn", datatype = "INT2U")
check(status_of(qc_run(dn, plots = plots4, expect = "reflectance"), "value_range") == "fail", "DN imagery fails the reflectance check")
check(status_of(qc_run(dn, plots = plots4, expect = "reflectance", scale = 1e-5), "value_range") == "pass", "the right --scale passes")
shifted <- mk(lapply(refl, function(m) m - 0.1), "shifted")
rs <- qc_run(shifted, plots = plots4, expect = "reflectance")
check(status_of(rs, "value_range") == "warn" && grepl("negative", rs$qc$message[rs$qc$check == "value_range"]),
      "negative reflectance (offset applied twice) warns")
check(status_of(qc_run(clean, expect_bands = 5), "band_count") == "fail", "wrong band count fails")
check(status_of(qc_run(mk(lapply(refl, function(m) round(m * 255)), "r8", datatype = "INT1U"), expect = "reflectance"),
                "datatype") == "warn", "8-bit reflectance warns")

# 4. NoData handling
rgb <- lapply(c(red = 150, green = 120, blue = 90), function(v) matrix(round(v + rnorm(nr * nc, 0, 10)), nr, nc))
rgb255 <- mk(rgb, "rgb255", datatype = "INT1U")                     # terra declares NoData = 255 by default
check(status_of(qc_run(rgb255, expect = "rgb"), "nodata_declared") == "warn", "NoData = 255 on 8-bit data warns")
collar <- lapply(rgb, function(m) { m[, 1:15] <- 0; m })
zero <- mk(collar, "zero", datatype = "INT1U", NAflag = NA)
rz <- qc_run(zero, expect = "rgb")
check(status_of(rz, "zero_as_nodata") == "warn" && abs(as.numeric(rz$qc$value[rz$qc$check == "zero_as_nodata"]) - 0.15) < 0.01,
      "15% zero collar without NoData warns")
zero_nd <- mk(collar, "zero_nd", datatype = "INT1U", NAflag = 0)
check(status_of(qc_run(zero_nd, expect = "rgb"), "zero_as_nodata") == "pass", "declared NoData = 0 passes")

# 5. Plot coverage, extent and saturation with exact values
holes <- lapply(refl, function(m) {
  m[(nr - 50 + 1):(nr - 10), 11:19] <- NA           # plot P1 (x 0.5-2 m, y 0.5-2.5 m): 9 of 30 columns = 30%
  m[(nr - 50 + 1):(nr - 10), 51:53] <- NA           # plot P2 (x 2.5-4 m): 3 of 30 columns = 10%
  m
})
rh <- qc_run(mk(holes, "holes"), plots = plots4)
check(abs(rh$plot_table$coverage[1] - 0.7) < 1e-6 && abs(rh$plot_table$coverage[2] - 0.9) < 1e-6,
      sprintf("plot coverage 0.700 and 0.900 (got %.3f, %.3f)", rh$plot_table$coverage[1], rh$plot_table$coverage[2]))
check(status_of(rh, "plot_coverage") == "fail", "30% NoData in a plot fails")
rh2 <- qc_run(mk(lapply(holes, function(m) { m[(nr - 50 + 1):(nr - 10), 11:19] <- 0.1; m }), "hole3"), plots = plots4)
check(status_of(rh2, "plot_coverage") == "warn", "10% NoData in a plot warns")
part <- mkplots(list(sq(0.5, 0.5, 2, 2.5), sq(4, 0.5, 5.5, 2.5)), "part")      # image is 5 m wide
rp <- qc_run(clean, plots = part)
check(status_of(rp, "extent") == "warn" && abs(rp$plot_table$inside_extent[2] - 2 / 3) < 1e-6 &&
        abs(rp$plot_table$coverage[2] - 2 / 3) < 1e-6, "plot two-thirds inside the image: extent warns, coverage 0.667")
out_pl <- mkplots(list(sq(0.5, 0.5, 2, 2.5), sq(20, 20, 21.5, 22)), "outside")
check(status_of(qc_run(clean, plots = out_pl), "extent") == "fail", "plot entirely outside the image fails")
sat <- lapply(refl, function(m) round(m * 1e4))
sat$red[(nr - 50 + 1):(nr - 10), 11:11] <- 65535                    # 1 of 30 columns of P1 = 3.3%
rsat <- qc_run(mk(sat, "sat", datatype = "INT2U", NAflag = 0), plots = plots4, expect = "reflectance", scale = 1e-4)
check(abs(rsat$plot_table$saturated_fraction[1] - 1 / 30) < 1e-6, "saturated fraction 1/30 in P1")
check(status_of(rsat, "saturation") == "warn", "3% saturated pixels warn")
sat$red[(nr - 50 + 1):(nr - 10), 12:14] <- 65535                    # 4 of 30 columns = 13%
check(status_of(qc_run(mk(sat, "sat2", datatype = "INT2U", NAflag = 0), plots = plots4, expect = "reflectance", scale = 1e-4),
                "saturation") == "fail", "13% saturated pixels fail")
rgba <- c(rgb, list(alpha = matrix(255, nr, nc)))
rgba$alpha[(nr - 50 + 1):(nr - 10), 11:25] <- 0                    # half of P1 transparent
ra <- qc_run(mk(rgba, "rgba", datatype = "INT1U", NAflag = NA), plots = plots4, expect = "rgb")
check(ra$alpha == 4L && abs(ra$plot_table$coverage[1] - 0.5) < 1e-6, "alpha band detected and honoured (coverage 0.5)")

# 6. DSM and DTM problems
check(status_of(qc_run(clean, plots = plots4, dsm = mk(list(d = matrix(120, nr, nc)), "dsm8", datatype = "INT1U", NAflag = NA)),
                "dsm_values") == "fail", "8-bit DSM fails")
check(status_of(qc_run(clean, plots = plots4, dsm = mk(list(d = ground + 99999), "dsmbig")), "dsm_values") == "fail",
      "implausible elevations fail")
spk <- ground + canopy
set.seed(2)
cells_in_p1 <- as.vector(outer((nr - 50 + 1):(nr - 10), (11:40 - 1) * nr, "+"))
hit <- sample(cells_in_p1, 12)                                       # 12 of 1200 pixels of P1 = 1%; 0.25% of all plots
spk[hit] <- spk[hit] + 5
rk <- qc_run(clean, plots = plots4, dsm = mk(list(d = spk), "spikes"))
check(abs(rk$plot_table$dsm_spike_fraction[1] - 12 / 1200) < 1e-6 && status_of(rk, "dsm_spikes") == "warn",
      sprintf("DSM spikes: fraction %.4f in P1, warn overall", rk$plot_table$dsm_spike_fraction[1]))
rn <- qc_run(clean, plots = plots4, dsm = dsm, dtm = mk(list(d = ground + 1), "dtm_high"))
check(status_of(rn, "canopy_height") == "warn" && all(rn$plot_table$chm_median_m < -0.2), "DTM above the DSM: negative heights warn")
rt <- qc_run(clean, plots = plots4, dsm = mk(list(d = ground + canopy * 10), "tall"), dtm = dtm)
check(status_of(rt, "canopy_height") == "warn" && all(rt$plot_table$chm_p95_m > 4), "6 m canopy exceeds --max-height-m 4")
dsm3857 <- file.path(tmp, "dsm3857.tif")
terra::writeRaster(terra::project(terra::rast(dsm), "EPSG:3857", res = 0.065), dsm3857, overwrite = TRUE)
rw <- qc_run(clean, plots = plots4, dsm = dsm3857)
check(status_of(rw, "dsm_crs") == "warn" && status_of(rw, "dsm_coverage") == "pass", "DSM in another CRS warns but is still matched")
expect_error(qc_run(clean, dtm = dtm), "--dtm needs --dsm", "a DTM without a DSM is refused")

# 7. Command line: --strict, demo end to end
expect_error(main(c("--input", dn, "--expect", "reflectance", "--strict", "--out", file.path(tmp, "strict"))),
             "QC failed", "--strict stops on a failed check")
expect_error(main(c("--demo", "--colour", "red")), "Unknown option", "unknown option is rejected")
out <- file.path(tempdir(), "claw-ortho-qc-demo-test")
unlink(out, recursive = TRUE)
cli <- run_cli(c("--demo", "--out", shQuote(out)))
check(cli$status == 0L, paste("demo run exits cleanly:", paste(cli$output, collapse = "\n")))
for (f in c("report.md", "tables/qc_checks.csv", "tables/band_stats.csv", "tables/plot_qc.csv", "figures/overview.png",
            "figures/band_histograms.png", "figures/dsm.png", "reproducibility/commands.sh", "reproducibility/sessionInfo.txt"))
  check(file.exists(file.path(out, f)), paste("demo writes", f))
dq <- utils::read.csv(file.path(out, "tables", "qc_checks.csv"))
check(!any(dq$status == "fail") && setequal(dq$check[dq$status == "warn"], c("extent", "plot_coverage", "dsm_coverage", "dsm_spikes")),
      "demo flags exactly the injected problems (extent, coverage, DSM coverage, spikes)")
dp <- utils::read.csv(file.path(out, "tables", "plot_qc.csv"))
check(abs(dp$coverage[dp$plot_id == "P06"] - 0.9) < 1e-6, "demo hole covers 10% of P06")
check(sum(file.size(Sys.glob(file.path(out, "demo_input", "*.tif")))) < 500000, "demo rasters are under 500 KB")
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
