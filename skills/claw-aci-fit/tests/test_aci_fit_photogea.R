# Known published answers for claw-aci-fit: reproduce the parameter tables in
# PhotoGEA's "Analyzing C3 A-Ci Curves" and "Analyzing C4 A-Ci Curves" articles
# (PhotoGEA 1.4.0, https://eloch216.github.io/PhotoGEA/articles/) from the
# example files shipped with PhotoGEA. Data preparation follows the articles
# exactly; fitting, consolidation and Jmax go through this skill's run_fits()
# and collect_results(). Curve i uses seed 1234 + i - 1, so the first curve
# reproduces the article's random stream exactly and the others converge to
# the same optimum (DEoptim, 200 generations).
# Run from anywhere: Rscript tests/test_aci_fit_photogea.R

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "aci_fit.R"))  # loads functions; main() does not run
suppressPackageStartupMessages(library(PhotoGEA))

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
close_to <- function(x, y, rel) length(x) == length(y) && all(is.na(x) == is.na(y)) &&
  all(abs(x[!is.na(y)] - y[!is.na(y)]) <= rel * abs(y[!is.na(y)]))
cores <- max(1L, min(6L, parallel::detectCores(), na.rm = TRUE))

load_example <- function(files) {
  gl <- lapply(files, function(f) read_gasex_file(PhotoGEA_example_file_path(f), "time"))
  cc <- do.call(identify_common_columns, gl)
  g <- do.call(rbind, lapply(gl, function(x) x[, cc, TRUE]))
  g[, "curve_identifier"] <- paste(g[, "instrument"], "-", g[, "species"], "-", g[, "plot"])
  g
}

# 1. C3: "Analyzing C3 A-Ci Curves" -------------------------------------------
g <- load_example(c("c3_aci_1.xlsx", "c3_aci_2.xlsx"))
check_response_curve_data(g, "curve_identifier", 16, "CO2_r_sp")
g <- organize_response_curve_data(g, "curve_identifier", c(9, 10), "Ci",
                                  columns_to_average = c("TleafCnd", "Qin"), print_information = FALSE)
g <- g[g[, "Stable"] == 2, , TRUE]
npts <- by(g, g[, "curve_identifier"], nrow)
g <- g[g[, "curve_identifier"] %in% names(npts[npts > 2]), , TRUE]
g <- remove_points(g, list(instrument = "ripe1", CO2_r_sp = 1800), method = "exclude")
g <- calculate_temperature_response(g, c3_temperature_param_sharkey)
g <- calculate_total_pressure(g)
g[, "curve_id"] <- g[, "curve_identifier"]
keep <- c("curve_id", "species", "A", "Ca", "Ci", "TleafCnd", "Oxygen", "total_pressure",
          "include_when_fitting", "TleafCnd_avg", "Qin_avg", "Gamma_star_at_25", "Gamma_star_norm",
          "Kc_at_25", "Kc_norm", "Ko_at_25", "Ko_norm", "J_norm", "RL_norm", "Tp_norm",
          "Vcmax_norm", "gmc_norm")
g <- g[, keep, TRUE]
opt3 <- DEFAULTS
opt3$pathway <- "C3"; opt3$ca_atmospheric <- 420; opt3$cores <- cores
t0 <- Sys.time()
c3 <- collect_results(run_fits(g, opt3), opt3)$params
cat(sprintf("C3: 6 curves fitted in %.0f s on %d core(s)\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")), cores))
ids <- c("ripe1 - soybean - 5a", "ripe1 - tobacco - 1", "ripe1 - tobacco - 2",
         "ripe2 - soybean - 1", "ripe2 - soybean - 5b", "ripe2 - tobacco - 4")
check(identical(c3$curve_id, ids), "C3: six curves in the article's order")
# Values printed in the article (PhotoGEA 1.4.0).
vc_art <- c(103.76133, 133.10357, 113.58840, 110.45928, 124.55297, 89.21564)
j_art <- c(NA, 236.1173, 208.5933, NA, 160.8934, NA)
jmax_art <- c(NA, 239.5027, 211.0427, NA, 162.1871, NA)
vc_lo_art <- c(96.90014, 126.18030, 104.97289, 103.52513, 107.49445, 85.58188)
vc_up_art <- c(111.42324, 140.02684, 121.19946, 121.96362, 141.61104, 92.84993)
check(isTRUE(all.equal(c3$vcmax_25_umol_m2_s[1], vc_art[1], tolerance = 1e-6)),
      "C3: first curve Vcmax_at_25 identical to the article (same random stream)")
check(close_to(c3$vcmax_25_umol_m2_s, vc_art, 0.005), "C3: all six Vcmax_at_25 within 0.5% of the article")
check(close_to(c3$vcmax_25_lower_umol_m2_s, vc_lo_art, 0.01) && close_to(c3$vcmax_25_upper_umol_m2_s, vc_up_art, 0.01),
      "C3: Vcmax confidence limits within 1% of the article")
check(close_to(c3$j_25_umol_m2_s, j_art, 0.005), "C3: J_at_25 (and which curves give NA) match the article")
check(close_to(c3$jmax_25_umol_m2_s, jmax_art, 0.005), "C3: Jmax_at_25 matches the article")
# Species means of TPU from the article's basic_stats table: soybean 9.909124, tobacco 11.168970,
# each from the single curve with a reliable TPU estimate.
check(isTRUE(all.equal(sort(c3$tp_25_umol_m2_s[!is.na(c3$tp_25_umol_m2_s)]), c(9.909124, 11.168970), tolerance = 0.005)),
      "C3: the two reliable TPU estimates match the article")
check(close_to(c3$rl_25_umol_m2_s[c(2, 3, 5)], c(2.2118752, 1.9404546, 2.7843997), 0.02),
      "C3: RL_at_25 matches the reference run for the J-limited curves")

# 2. C4: "Analyzing C4 A-Ci Curves", mechanistic fit assuming Rubisco limitation
g4 <- load_example(c("c4_aci_1.xlsx", "c4_aci_2.xlsx"))
g4 <- organize_response_curve_data(g4, "curve_identifier", c(9, 10), "Ci",
                                   columns_to_average = c("TleafCnd", "Qin"), print_information = FALSE)
g4 <- remove_points(g4, list(instrument = "ripe1", CO2_r_sp = 1800), method = "exclude")
g4 <- calculate_temperature_response(g4, c4_temperature_param_vc)
g4 <- calculate_total_pressure(g4)
g4[, "curve_id"] <- g4[, "curve_identifier"]
g4 <- g4[, c("curve_id", "species", "A", "Ca", "Ci", "TleafCnd", "Oxygen", "total_pressure",
             "include_when_fitting", "TleafCnd_avg", "Qin_avg", "Vcmax_norm", "Vpmax_norm", "RL_norm",
             "Kc", "Ko", "Kp", "gamma_star", "ao", "gmc_norm", "J_norm"), TRUE]
opt4 <- opt3
opt4$pathway <- "C4"
t0 <- Sys.time()
c4 <- collect_results(run_fits(g4, opt4), opt4)$params
cat(sprintf("C4: 5 curves fitted in %.0f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
check(nrow(c4) == 5L, "C4: five curves")
vp_art <- c(158.2040, 149.4453, 124.2466, 154.8675, 131.1353)
check(isTRUE(all.equal(c4$vpmax_25_umol_m2_s[1], vp_art[1], tolerance = 1e-6)),
      "C4: first curve Vpmax_at_25 identical to the article")
check(close_to(c4$vpmax_25_umol_m2_s, vp_art, 0.005), "C4: all five Vpmax_at_25 within 0.5% of the article")
check(close_to(c4$vpmax_25_lower_umol_m2_s, c(133.8529, 123.4999, 115.9223, 135.7686, 110.0679), 0.01) &&
      close_to(c4$vpmax_25_upper_umol_m2_s, c(187.1139, 179.5492, 133.2167, 176.9160, 156.6364), 0.01),
      "C4: Vpmax confidence limits within 1% of the article")
sp <- ifelse(grepl("maize", c4$curve_id), "maize", "sorghum")
vc_mean <- as.numeric(tapply(c4$vcmax_25_umol_m2_s, sp, mean)[c("maize", "sorghum")])
check(isTRUE(all.equal(vc_mean, c(36.53579, 38.39997), tolerance = 0.005)),
      "C4: species mean Vcmax_at_25 (maize 36.54, sorghum 38.40) matches the article")
rl_mean <- as.numeric(tapply(c4$rl_25_umol_m2_s, sp, mean)[c("maize", "sorghum")])
check(isTRUE(all.equal(rl_mean, c(-0.1333101, -0.9351559), tolerance = 0.05)),
      "C4: species mean RL_at_25 matches the article (negative values, as published)")

cat(sprintf("PASS: %d checks (%s, PhotoGEA articles)\n", n_checks, SKILL))
