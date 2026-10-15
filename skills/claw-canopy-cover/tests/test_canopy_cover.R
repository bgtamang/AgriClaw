# Tests for claw-canopy-cover. Run from anywhere: Rscript tests/test_canopy_cover.R
# Known answers: hand-computed ExG/ExGR, a hand-checked Otsu histogram and Otsu's
# separability for a normal distribution (2/pi), and synthetic RGB scenes whose
# true canopy cover per plot is known exactly from the generator (anti-aliased
# leaf edges, noise, illumination gradient, residue, glare, shadow).

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "canopy_cover.R"))  # loads functions; main() does not run
source(file.path(SKILL_TEST_DIR, "demo", "make_demo.R"))        # synth_rgb_scene()

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
  script <- file.path(SKILL_TEST_DIR, "scripts", "canopy_cover.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
tmp <- file.path(tempdir(), "claw-canopy-cover-tests")
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
save_scene <- function(sc, name, datatype = "INT1U", NAflag = NA, scale = 1, alpha = NULL) {
  r <- sc$raster * scale
  if (!is.null(alpha)) r <- c(r, alpha)
  rp <- file.path(tmp, paste0(name, ".tif"))
  terra::writeRaster(r, rp, datatype = datatype, overwrite = TRUE, NAflag = NAflag)
  pp <- file.path(tmp, paste0(name, ".gpkg"))
  if (file.exists(pp)) unlink(pp)
  sf::st_write(sc$plots, pp, quiet = TRUE)
  list(raster = rp, plots = pp)
}
opt_for <- function(paths, ...) {
  o <- DEFAULTS
  o$input <- paths$raster; o$plots <- paths$plots
  mods <- list(...)
  for (k in names(mods)) o[[k]] <- mods[[k]]
  o
}

# 1. Colour indices by hand
ci <- color_indices(60, 120, 20)                  # r = 0.3, g = 0.6, b = 0.1
check(abs(ci$exg - 0.8) < 1e-12, "ExG = 2g - r - b = 0.8")
check(abs(ci$exgr - 0.98) < 1e-12, "ExGR = ExG - (1.4r - g) = 0.98")
cs <- color_indices(140, 110, 85)
check(abs(cs$exg - (2 * 110 - 140 - 85) / 335) < 1e-12 && cs$exgr < 0, "soil: ExG near 0, ExGR negative")
check(is.na(color_indices(0, 0, 0)$exg), "black pixels have no chromatic coordinates")
check(abs(color_indices(120, 240, 40)$exg - 0.8) < 1e-12, "ExG ignores overall brightness")

# 2. Otsu: hand-checked histogram, plateau, separability of a normal distribution
check(abs(otsu_from_hist(c(3, 5, 2, 1, 1, 6, 3), 0:6) - 3.5) < 1e-12, "Otsu on the hand-checked histogram is 3.5")
check(abs(otsu_from_hist(c(5, 5, 0, 0, 0, 5, 5), 0:6) - 3) < 1e-12, "Otsu takes the middle of an empty gap")
check(abs(attr(otsu_from_hist(c(10, 0, 0, 10), 1:4), "eta") - 1) < 1e-12, "two spikes: separability 1")
set.seed(11)
eta_n <- attr(otsu_threshold(rnorm(2e5)), "eta")
check(abs(eta_n - 2 / pi) < 0.01, sprintf("one normal mode: separability %.3f = 2/pi", eta_n))

# 3. Artifact-free scene with known cover, including bare (0%) and full (100%) plots
sc <- synth_rgb_scene(seed = 5L, artifacts = FALSE, radius = c(0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.25, 0.3, 0.6))
truth <- sc$plots$true_cover_pct
check(min(truth) == 0 && max(truth) == 100, "scene has a bare plot and a full-canopy plot")
p8 <- save_scene(sc, "clean")
res <- run_cover(opt_for(p8))
err <- res$table$canopy_cover_pct - truth
check(max(abs(err)) <= 2, sprintf("ExG + Otsu recovers cover within 2 percentage points (max error %.2f)", max(abs(err))))
check(abs(mean(err)) <= 1, sprintf("ExG + Otsu mean bias within 1 pp (%.2f)", mean(err)))
err2 <- res$table$cover_exgr_pct - truth
check(max(abs(err2)) <= 2, sprintf("ExGR > 0 cross-check within 2 pp (max error %.2f)", max(abs(err2))))
full <- which(truth == 100); bare <- which(truth == 0)
check(res$table$canopy_cover_pct[full] > 98 && res$table$canopy_cover_pct[bare] < 1,
      "one threshold per flight: full plot near 100%, bare plot near 0%")
# Why not per plot: Otsu on the full-canopy plot alone splits the canopy
dfv <- exactextractr::exact_extract(res$px[["exg"]], res$plots[full, ], progress = FALSE)[[1]]
vals <- dfv$value[dfv$coverage_fraction >= 0.99 & !is.na(dfv$value)]
own <- otsu_threshold(vals)
check(mean(vals > own) > 0.2 && mean(vals > own) < 0.8,
      sprintf("a per-plot Otsu threshold would call only %.0f%% of a full canopy green", 100 * mean(vals > own)))
check(res$threshold$method == "Otsu" && res$threshold$eta > 0.8, "clean scene is clearly bimodal")
check(all(res$table$valid_fraction > 0.999) && all(res$table$saturated_fraction == 0) && all(res$table$dark_fraction == 0),
      "clean scene: full coverage, no saturated or dark pixels")
r_exgr <- run_cover(opt_for(p8, index = "exgr"))
check(r_exgr$threshold$value == 0 && max(abs(r_exgr$table$canopy_cover_pct - truth)) <= 2,
      "--index exgr uses the fixed zero threshold and recovers cover")
check(max(abs(r_exgr$table$canopy_cover_pct - res$table$cover_exgr_pct)) < 1e-9, "cross-check column equals the ExGR run")

# 4. Same scene as 16-bit data: identical cover (chromatic coordinates are scale-free)
p16 <- save_scene(sc, "clean16", datatype = "INT2U", scale = 257)
r16 <- run_cover(opt_for(p16))
check(max(abs(r16$table$canopy_cover_pct - res$table$canopy_cover_pct)) < 1e-6 && r16$max_value == 65535,
      "16-bit copy gives the same cover; saturation level 65535")

# 5. Scene with residue, glare and shadow
sa <- synth_rgb_scene(seed = 5L, artifacts = TRUE)
pa <- save_scene(sa, "artifacts")
ra <- run_cover(opt_for(pa))
tb <- ra$table
r8 <- terra::rast(pa$raster)
count_in <- function(i, fun) {
  m <- terra::mask(terra::crop(r8, terra::vect(sa$plots[i, ])), terra::vect(sa$plots[i, ]))
  v <- terra::values(m)
  c(hit = sum(fun(v), na.rm = TRUE), n = sum(stats::complete.cases(v)))
}
s3 <- count_in(3, function(v) v[, 1] >= 255 | v[, 2] >= 255 | v[, 3] >= 255)
check(abs(tb$saturated_fraction[3] - s3[["hit"]] / s3[["n"]]) < 1e-9 && s3[["hit"]] > 0,
      sprintf("saturated fraction in P03 = %d / %d pixels", s3[["hit"]], s3[["n"]]))
d5 <- count_in(5, function(v) (v[, 1] + v[, 2] + v[, 3]) / 765 < 0.05 & !(v[, 1] >= 255 | v[, 2] >= 255 | v[, 3] >= 255))
check(abs(tb$dark_fraction[5] - d5[["hit"]] / d5[["n"]]) < 1e-9 && grepl("dark", tb$flag[5]),
      sprintf("dark fraction in P05 = %d / %d pixels, flagged", d5[["hit"]], d5[["n"]]))
check(all(tb$dark_fraction[-5] == 0), "no dark pixels outside the shadow strip")
check(abs(tb$canopy_cover_pct[7] - sa$plots$true_cover_pct[7]) <= 2 && abs(tb$cover_exgr_pct[7] - sa$plots$true_cover_pct[7]) <= 2,
      "straw residue is not counted as canopy (ExG and ExGR)")
check(max(abs(tb$canopy_cover_pct[-5] - sa$plots$true_cover_pct[-5])) <= 2, "other plots still within 2 pp")
check(tb$n_pixels[3] < tb$n_pixels[1] && tb$n_pixels[5] < tb$n_pixels[1], "glare and shadow pixels are left out of the count")

# 6. NoData declared at the saturation level hides glare; alpha bands are honoured
pn <- save_scene(sa, "nodata255", NAflag = 255)
rn <- run_cover(opt_for(pn))
check(any(grepl("NoData = 255", rn$warnings)) && rn$table$valid_fraction[3] < 1 && rn$table$saturated_fraction[3] == 0,
      "NoData = 255 is reported: glare becomes missing coverage, not saturation")
alpha <- terra::rast(sc$raster[[1]]); terra::values(alpha) <- 255
e9 <- terra::ext(terra::vect(sc$plots[9, ]))
cells <- terra::cells(alpha, terra::ext(e9$xmin, e9$xmin + 0.5, e9$ymin, e9$ymax))
alpha[cells] <- 0
pal <- save_scene(sc, "alpha", alpha = alpha)
ral <- run_cover(opt_for(pal))
v9 <- terra::vect(sc$plots[9, ])
a9 <- terra::values(terra::mask(terra::crop(alpha, v9), v9))
exp_valid <- mean(a9 > 0, na.rm = TRUE)
check(!is.na(ral$alpha_idx) && abs(ral$table$valid_fraction[9] - exp_valid) < 1e-6 && grepl("incomplete", ral$table$flag[9]),
      sprintf("alpha band detected; plot 9 valid fraction %.3f (expected %.3f)", ral$table$valid_fraction[9], exp_valid))
expect_error(run_cover(opt_for(p8, bands = "red=1,green=1,blue=3")), "different layer", "duplicate band numbers refused")
expect_error(run_cover(opt_for(p8, threshold = "abc")), "threshold", "bad --threshold refused")

# 7. End to end on the demo data
out <- file.path(tempdir(), "claw-canopy-cover-demo-test")
unlink(out, recursive = TRUE)
cli <- run_cli(c("--demo", "--out", shQuote(out)))
check(cli$status == 0L, paste("demo run exits cleanly:", paste(cli$output, collapse = "\n")))
for (f in c("report.md", "tables/canopy_cover.csv", "spatial/veg_mask.tif", "figures/index_histogram.png",
            "figures/cover_map.png", "figures/qc_plots.png", "reproducibility/commands.sh",
            "reproducibility/sessionInfo.txt"))
  check(file.exists(file.path(out, f)), paste("demo writes", f))
dt <- utils::read.csv(file.path(out, "tables", "canopy_cover.csv"))
check(nrow(dt) == 9L && all(c("plot_id", "date", "canopy_cover_pct", "cover_exgr_pct", "flag") %in% names(dt)),
      "demo table has 9 plots")
check(max(abs(dt$canopy_cover_pct - dt$true_cover_pct)[dt$plot_id != "P05"]) <= 2,
      "demo cover within 2 pp of the true cover (except the shadowed plot)")
check(file.size(file.path(out, "demo_input", "rgb_ortho.tif")) < 500000, "demo raster is under 500 KB")
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
expect_error(main(c("--demo", "--colour", "red")), "Unknown option", "unknown option is rejected")
expect_error(main(c("--input", p8$raster, "--out", file.path(tmp, "x"))), "--plots", "missing --plots is rejected")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
