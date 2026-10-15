# Tests for claw-stomatal-models. Run from anywhere: Rscript tests/test_stomatal_models.R
# Known answers: (1) leaf-surface hs, Cs and D recovered exactly from simulated
# LI-6800 columns; (2) known g0 and g1 recovered from noisy curves simulated
# with each model; (3) unit mistakes caught; (4) iWUE = A / gsw; (5) the
# Ball-Berry table published in PhotoGEA's "Analyzing Ball-Berry Data" article
# reproduced from PhotoGEA's example files; (6) an end-to-end --demo run.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "stomatal_models.R"))   # loads functions; main() does not run
source(file.path(SKILL_TEST_DIR, "demo", "make_demo.R"))            # simulate_point(), simulate_curve()

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "stomatal_models.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
tmp <- file.path(tempdir(), "claw-stomatal-models-test")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE)
opt <- DEFAULTS
opt$curve_cols <- "curve"

# Curves vary light (A 3-35) and chamber humidity (leaf-to-air VPD about 0.8-2.6 kPa).
sim_set <- function(model, truth, n_pts = 8L, noise = 0.01) {
  do.call(rbind, lapply(seq_len(nrow(truth)), function(i) {
    a <- stats::runif(n_pts, 3, 35)
    h2o <- stats::runif(n_pts, 9.7, 26.8)
    data.frame(curve = truth$curve[i], simulate_curve(a, h2o, model, truth$g0[i], truth$g1[i], noise = noise))
  }))
}

# 1. Leaf-surface quantities recovered exactly from the logged columns ----------------
set.seed(11)
one <- sim_set("medlyn", data.frame(curve = "x", g0 = 0.03, g1 = 4), noise = 0)
one$used_in_fit <- TRUE
one$curve_id <- one$curve
x <- build_exdf(one, opt, precomputed = FALSE)
check(isTRUE(all.equal(x[, "RHleaf"] / 100, one$true_hs, tolerance = 1e-9)), "hs from calculate_gas_properties = simulated hs")
check(isTRUE(all.equal(x[, "Csurface"], one$true_cs_umol_mol, tolerance = 1e-9)), "Cs = simulated Cs")
check(isTRUE(all.equal(x[, "VPDsurf"], one$true_d_kpa, tolerance = 1e-9)), "leaf-surface VPD = simulated D")
check(all(one$vpd_leaf_kpa > one$true_d_kpa), "LI-6800 VPDleaf (to chamber air) exceeds leaf-surface VPD")
# nls (inside fit_medlyn) cannot converge on exactly zero-residual data (see ?nls),
# so the near-exact recovery check uses scatter of sd 1e-6 mol m-2 s-1.
set.seed(12)
tiny <- sim_set("medlyn", data.frame(curve = "x", g0 = 0.03, g1 = 4), noise = 1e-6)
r1 <- analyse_stomatal(tiny, opt)$params
check(abs(r1$medlyn_g1_kpa05 - 4) / 4 < 1e-4 && abs(r1$medlyn_g0_mol_m2_s - 0.03) < 1e-5,
      "near-noise-free Medlyn curve: g1 and g0 recovered to 1e-4")
set.seed(13)
tiny_bb <- sim_set("ball-berry", data.frame(curve = "x", g0 = 0.02, g1 = 9), noise = 1e-6)
r1b <- analyse_stomatal(tiny_bb, modifyList(opt, list(models = "ball-berry")))$params
check(abs(r1b$bb_g1 - 9) / 9 < 1e-4 && abs(r1b$bb_g0_mol_m2_s - 0.02) < 1e-5,
      "near-noise-free Ball-Berry curve: g1 and g0 recovered to 1e-4")

# 2. Recovery of known g0 and g1 from noisy curves (gsw scatter sd 0.01) ---------------
set.seed(21)
tb <- data.frame(curve = sprintf("bb%02d", 1:10), g0 = stats::runif(10, 0.01, 0.05), g1 = stats::runif(10, 6, 12))
rb <- analyse_stomatal(sim_set("ball-berry", tb), modifyList(opt, list(models = "ball-berry")))$params
rb <- rb[match(tb$curve, rb$curve_id), ]
e_g1 <- (rb$bb_g1 - tb$g1) / tb$g1
e_g0 <- rb$bb_g0_mol_m2_s - tb$g0
cat(sprintf("Ball-Berry: max |g1 error| %.1f%%, max |g0 error| %.4f\n", 100 * max(abs(e_g1)), max(abs(e_g0))))
# Errors are judged against the fit's own standard errors (z = error / SE), plus absolute bounds.
zb <- (rb$bb_g1 - tb$g1) / rb$bb_g1_se
check(all(abs(zb) <= 3.5) && mean(abs(zb) <= 2) >= 0.75,
      "Ball-Berry g1 errors within 3.5 SE on every curve and within 2 SE on >= 75%")
check(max(abs(e_g1)) <= 0.15 && stats::median(abs(e_g1)) <= 0.05, "Ball-Berry g1 within 15% (median within 5%)")
check(max(abs(e_g0)) <= 0.02, "Ball-Berry g0 within 0.02 mol m-2 s-1")

set.seed(31)
tm <- data.frame(curve = sprintf("md%02d", 1:10), g0 = stats::runif(10, 0.01, 0.05), g1 = stats::runif(10, 2, 6))
rm_ <- analyse_stomatal(sim_set("medlyn", tm), modifyList(opt, list(models = "medlyn")))$params
rm_ <- rm_[match(tm$curve, rm_$curve_id), ]
e_m1 <- (rm_$medlyn_g1_kpa05 - tm$g1) / tm$g1
e_m0 <- rm_$medlyn_g0_mol_m2_s - tm$g0
cat(sprintf("Medlyn: max |g1 error| %.1f%%, max |g0 error| %.4f\n", 100 * max(abs(e_m1)), max(abs(e_m0))))
zm <- (rm_$medlyn_g1_kpa05 - tm$g1) / rm_$medlyn_g1_se
check(all(abs(zm) <= 3.5) && mean(abs(zm) <= 2) >= 0.75,
      "Medlyn g1 errors within 3.5 SE on every curve and within 2 SE on >= 75%")
check(max(abs(e_m1)) <= 0.15 && stats::median(abs(e_m1)) <= 0.05, "Medlyn g1 within 15% (median within 5%)")
check(max(abs(e_m0)) <= 0.02, "Medlyn g0 within 0.02 mol m-2 s-1")
check(all(rm_$qc_flags == "" | rm_$qc_flags == "negative_g0"), "clean simulated curves raise no unit flags")

# Precomputed leaf-surface columns give the same fit as the raw LI-6800 columns.
r1 <- analyse_stomatal(tiny, opt)$params
pre <- tiny
pre$rh_leaf_pct <- 100 * pre$true_hs
pre$cs_umol_mol <- pre$true_cs_umol_mol
pre$vpd_leaf_surface_kpa <- pre$true_d_kpa
pre <- pre[, c("curve", "a_umol_m2_s", "gsw_mol_m2_s", "rh_leaf_pct", "cs_umol_mol", "vpd_leaf_surface_kpa")]
r_pre <- analyse_stomatal(pre, opt)$params
check(isTRUE(all.equal(r_pre$medlyn_g1_kpa05, r1$medlyn_g1_kpa05, tolerance = 1e-7)) &&
      isTRUE(all.equal(r_pre$bb_g1, r1$bb_g1, tolerance = 1e-7)), "precomputed and raw inputs give identical fits")

# 3. Unit mistakes are caught before fitting ------------------------------------------------
err <- function(d, o = opt) tryCatch({ analyse_stomatal(d, o); "no error" }, error = function(e) conditionMessage(e))
check(grepl("Pa", err(transform(pre, vpd_leaf_surface_kpa = vpd_leaf_surface_kpa * 1000))), "VPD in Pa is caught")
check(grepl("hPa", err(transform(pre, vpd_leaf_surface_kpa = vpd_leaf_surface_kpa * 10))), "VPD in hPa is caught")
check(grepl("mmol", err(transform(pre, gsw_mol_m2_s = gsw_mol_m2_s * 1000))), "gsw in mmol m-2 s-1 is caught")
check(grepl("fraction", err(transform(pre, rh_leaf_pct = rh_leaf_pct / 100))), "hs as a fraction is caught")
check(grepl("Pa", err(transform(pre, cs_umol_mol = cs_umol_mol * 0.0996))), "Cs in Pa is caught")
raw <- tiny
check(grepl("0-100%", err(transform(raw, e_mol_m2_s = e_mol_m2_s * 1000))), "E in mmol (impossible humidity) is caught")
lic <- raw
lic$vpd_leaf_kpa <- lic$vpd_leaf_kpa * 1000
check(grepl("Pa", err(lic, modifyList(opt, list(vpd_source = "licor")))), "LI-6800 VPDleaf in Pa is caught")

gap <- raw
gap$e_mol_m2_s[3] <- NA
rg <- analyse_stomatal(gap, opt)
check(rg$params$n_points_used == nrow(raw) - 1L && grepl("missing_inputs", rg$data$excluded_reason[3]),
      "a row with missing E is excluded with a reason, not dropped silently")

# 4. iWUE = A / gsw ---------------------------------------------------------------------------
pts <- analyse_stomatal(raw, opt)$points
check(isTRUE(all.equal(pts$iwue_umol_mol, raw$a_umol_m2_s / raw$gsw_mol_m2_s)), "iWUE = A / gsw at every point")

# 5. PhotoGEA "Analyzing Ball-Berry Data" article ------------------------------------------------
suppressPackageStartupMessages(library(PhotoGEA))
gl <- lapply(c("ball_berry_1.xlsx", "ball_berry_2.xlsx"),
             function(f) read_gasex_file(PhotoGEA_example_file_path(f), "time"))
cc <- do.call(identify_common_columns, gl)
g <- do.call(rbind, lapply(gl, function(z) z[, cc, TRUE]))$main_data
bbt <- data.frame(instrument = g$instrument, species = g$species, plot = g$plot,
                  a_umol_m2_s = g$A, gsw_mol_m2_s = g$gsw, e_mol_m2_s = g$E, gbw_mol_m2_s = g$gbw,
                  h2o_s_mmol_mol = g$H2O_s, ca_umol_mol = g$Ca, tleaf_c = g$TleafCnd, pa_kpa = g$Pa,
                  delta_pcham_kpa = g$DeltaPcham, vpd_leaf_kpa = g$VPDleaf,
                  qc_unstable = g$Stable < g$Total, stringsAsFactors = FALSE)
# The article drops three curves it judged abnormal and excludes unstable points.
drop <- paste(bbt$instrument, bbt$species, bbt$plot) %in% c("ripe4 soybean 1a", "ripe4 soybean 1b", "ripe4 soybean 5")
ob <- modifyList(opt, list(curve_cols = "instrument,species,plot", models = "ball-berry"))
rv <- analyse_stomatal(bbt[!drop, ], ob)$params
check(identical(rv$curve_id, c("ripe4 - tobacco - 2", "ripe5 - soybean - 1", "ripe5 - soybean - 5a",
                               "ripe5 - soybean - 5c", "ripe5 - tobacco - 5")), "article: five curves")
check(isTRUE(all.equal(rv$bb_g0_mol_m2_s, c(0.08778255, 0.21437406, 0.10722148, 0.11118324, 0.13565311), tolerance = 1e-6)),
      "article: Ball-Berry intercepts reproduced")
check(isTRUE(all.equal(rv$bb_g1, c(3.864316, 4.462244, 4.823671, 4.750743, 4.097481), tolerance = 1e-6)),
      "article: Ball-Berry slopes reproduced")
check(isTRUE(all.equal(rv$bb_r2, c(0.6452714, 0.9004803, 0.8866346, 0.7726284, 0.9849869), tolerance = 1e-6)),
      "article: R-squared reproduced")

# 6. End to end on the demo data -----------------------------------------------------------------
out <- file.path(tmp, "demo")
r <- run_cli(c("--demo", "--out", shQuote(out)))
check(r$status == 0L, paste("demo run exits cleanly:", paste(r$output, collapse = "\n")))
for (f in c("report.md", "tables/stomatal_parameters.csv", "tables/stomatal_points.csv", "figures/stomatal_fits.png",
            "reproducibility/commands.sh", "reproducibility/sessionInfo.txt"))
  check(file.exists(file.path(out, f)), paste("demo writes", f))
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
dp <- utils::read.csv(file.path(out, "tables", "stomatal_parameters.csv"), stringsAsFactors = FALSE)
check(nrow(dp) == 4L && all(abs(dp$medlyn_g1_kpa05 - c(3.5, 3.5, 4.5, 4.5)) / c(3.5, 3.5, 4.5, 4.5) <= 0.10),
      "demo: Medlyn g1 within 10% of the generating values")
check(identical(unique(dp$vpd_source), "surface"), "demo records the VPD definition used")
bad_run <- run_cli(c("--demo", "--vpd-source", "air"))
check(bad_run$status != 0L, "invalid --vpd-source exits with an error")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
