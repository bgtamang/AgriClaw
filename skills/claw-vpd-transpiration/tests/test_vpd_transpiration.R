# Tests for claw-vpd-transpiration. Run from anywhere: Rscript tests/test_vpd_transpiration.R
# Known-answer checks: VPD against hand calculations and FAO-56 worked examples,
# TR unit conversion, irrigation/spike removal, AICc against a GraphPad Prism
# output, breakpoint recovery and false-breakpoint rates on simulated data,
# bad-input errors, and end-to-end command-line runs.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "vpd_transpiration.R"))  # loads functions; main() does not run

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
near <- function(a, b, tol) isTRUE(all(abs(a - b) <= tol))

run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "vpd_transpiration.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
errmsg <- function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))

# 1. VPD against hand calculations (computed independently in Python) ----------
check(near(vpd_kpa(30, 50), 2.130572026, 1e-8), "VPD at 30 degC, 50 % RH (Buck 1981 moist-air form)")
check(near(vpd_kpa(20, 70), 0.704098649, 1e-8), "VPD at 20 degC, 70 % RH")
check(near(svp_kpa(25, "buck1981"), 3.167034531, 1e-8), "Buck (1981) pure-water SVP at 25 degC")
check(near(vpd_kpa(25, 100), 0, 1e-12) && near(vpd_kpa(25, 0), svp_kpa(25), 1e-12),
      "VPD is 0 at RH 100 % and equals SVP at RH 0 %")
check(near(vpd_kpa(25, 101), 0, 1e-12), "RH 100-103 % is treated as saturation")
# FAO-56 worked examples (Allen et al. 1998, chapter 3, examples 3-5), to 3 decimals
fao_t <- c(24.5, 15, 18, 25, 19.5)
fao_es <- c(3.075, 1.705, 2.064, 3.168, 2.267)
check(near(svp_kpa(fao_t, "fao56"), fao_es, 0.0006), "FAO-56 e0(T) matches the FAO-56 worked examples")
tt <- seq(5, 45, by = 0.1)
check(max(abs(svp_kpa(tt, "buck1981") / svp_kpa(tt, "fao56") - 1)) < 0.003,
      "Buck (1981) and FAO-56 SVP agree within 0.3 % from 5 to 45 degC")
check(near(svp_kpa(tt, "buck1981_moist") / svp_kpa(tt, "buck1981"), 0.61375 / 0.61121, 1e-12),
      "moist-air form = Buck pure-water form x 1.0042")

# 2. TR unit conversion on a noiseless pot ---------------------------------------
t0 <- as.POSIXct("2026-07-01 10:00:00", tz = "UTC")
tm <- t0 + 60 * (0:60)
w <- data.frame(plant_id = "X", time = tm, weight_g = 2000 - 0.12 * (0:60))
wq <- qc_weights(w)
win <- data.frame(step = "S1", start = t0, end = t0 + 3600, phase = "day",
                  chamber = NA_character_, stringsAsFactors = FALSE)
st <- step_tr(wq, data.frame(plant_id = "X", leaf_area_cm2 = 500), win, 0, 10)
check(near(st$tr_g_h, 7.2, 1e-9), "0.12 g/min = 7.2 g/h")
check(near(st$tr_mg_m2_s, 40, 1e-9), "0.12 g/min over 500 cm2 = 40 mg m-2 s-1 (hand calculation)")
check(near(st$tr_mmol_m2_s, 40 / 18.015, 1e-9), "mg to mmol conversion uses 18.015 mg/mmol")

# 3. Irrigation and spike removal ---------------------------------------------------
set.seed(11)
mins <- 0:179
true_rate <- 0.1                                     # g/min
wt <- 2500 - true_rate * mins + rnorm(length(mins), 0, 0.02)
wt[mins >= 50] <- wt[mins >= 50] + 100               # irrigation at minute 50
wt[mins == 120] <- wt[mins == 120] + 30              # one-reading spike at minute 120
pot <- data.frame(plant_id = "P", time = t0 + 60 * mins, weight_g = wt)
q <- qc_weights(pot, irrigation_g = 5, settle_min = 30)
check(identical(which(q$irrigation), 51L), "irrigation detected exactly at minute 50")
check(q$outlier[mins == 120] && !q$irrigation[mins == 120], "one-reading spike flagged as an outlier, not an irrigation")
check(all(!q$use[mins >= 50 & mins < 80]) && mean(q$use[mins >= 80 & mins != 120]) >= 0.95,
      "readings in the 30-min settling window are excluded and later ones kept")
check(sum(q$outlier) <= 3, "few false spike flags on clean data")
r <- step_tr(q, data.frame(plant_id = "P", leaf_area_cm2 = 1e4),
             data.frame(step = "all", start = t0, end = t0 + 180 * 60, phase = "day",
                        chamber = NA_character_, stringsAsFactors = FALSE), 0, 10)
check(near(r$tr_g_h / 60, true_rate, 0.002), "TR across an irrigation and a spike within 2 % of the truth")
naive <- -mean(diff(wt))                             # mean of minute-by-minute losses
check(naive < 0, "without screening, the mean of minute differences is wrong (negative TR)")
noisy <- data.frame(plant_id = "N", time = t0 + 60 * (0:599),
                    weight_g = 2500 - 0.2 * (0:599) + rnorm(600, 0, 0.05))
qn <- qc_weights(noisy)
check(!any(qn$irrigation | qn$disturbance) && mean(qn$outlier) < 0.02,
      "clean series: no jumps and < 2 % chance spike flags")
drop <- noisy; drop$weight_g[301:600] <- drop$weight_g[301:600] - 40
qd <- qc_weights(drop)
check(identical(which(qd$disturbance), 301L) && !any(qd$irrigation), "a 40 g drop is flagged as a disturbance")

# 4. AICc reproduces a GraphPad Prism comparison (Tamang & Sadok repository,
#    'AICc method/Nonlin fit of WBDC.016-1.csv': n = 36, SS line 444.2,
#    SS segmental 345.6, Prism: difference in AICc 3.788, P(segmental) 86.92 %)
d_aicc <- aicc_rss(444.2, 36, 3) - aicc_rss(345.6, 36, 5)
check(near(d_aicc, 3.788, 0.01), "AICc difference matches Prism (3.788)")
check(near(1 / (1 + exp(-0.5 * d_aicc)), 0.8692, 0.001), "Akaike weight matches Prism (86.92 %)")
check(is.na(aicc_rss(1, 6, 5)), "AICc is undefined when n <= k + 1")

# 5. Simulations with known truth ----------------------------------------------------
# Design: 7 VPD steps (0.9-3.9 kPa) x 6 plants, residual SD 4 mg m-2 s-1.
vpd_steps <- c(0.9, 1.4, 1.9, 2.4, 2.9, 3.4, 3.9)
sim <- function(linear, bp = 2.2, s1 = 30, s2 = 8, sd = 4) {
  v <- rep(vpd_steps, 6) + rnorm(42, 0, 0.05)
  mu <- if (linear) -10 + s1 * v else -10 + s1 * v + (s2 - s1) * pmax(v - bp, 0)
  data.frame(vpd = v, tr = mu + rnorm(42, 0, sd))
}
set.seed(2026)
n_sim <- 30
seg_rows <- do.call(rbind, lapply(seq_len(n_sim), function(i) {
  d <- sim(FALSE); fit_tr_vpd(d$vpd, d$tr)$row }))
hit <- mean(abs(seg_rows$breakpoint_kpa - 2.2) <= 0.2)
hit3 <- mean(abs(seg_rows$breakpoint_kpa - 2.2) <= 0.3)
cls <- mean(seg_rows$model == "segmented")
check(length(unique(seg_rows$breakpoint_kpa)) == n_sim, "simulated data sets are independent (segmented does not reset the RNG)")
check(hit >= 0.75 && hit3 >= 0.9,
      sprintf("breakpoint within +/-0.2 kPa of 2.2 in >= 75%% (got %.0f%%) and within +/-0.3 kPa in >= 90%% (got %.0f%%)", 100 * hit, 100 * hit3))
check(cls >= 0.9, sprintf("two-segment truth classified as segmented in >= 90%% (got %.0f%%)", 100 * cls))
check(abs(stats::median(seg_rows$slope1) - 30) < 2 && abs(stats::median(seg_rows$slope2) - 8) < 3,
      "median slopes close to the true 30 and 8")
cover <- mean(seg_rows$breakpoint_lo <= 2.2 & seg_rows$breakpoint_hi >= 2.2)
check(cover >= 0.85, sprintf("score-based breakpoint 95%% CI covers the truth in >= 85%% (got %.0f%%)", 100 * cover))
red <- seg_rows$slope_reduction_pct
check(near(red, 100 * (seg_rows$slope1 - seg_rows$slope2) / seg_rows$slope1, 1e-9),
      "slope reduction % = 100 (slope1 - slope2) / slope1")
set.seed(7)
lin_sets <- lapply(seq_len(n_sim), function(i) sim(TRUE))
lin_rows <- do.call(rbind, lapply(lin_sets, function(d) fit_tr_vpd(d$vpd, d$tr)$row))
# AICc-only rule applied to the same fits: breakpoint if valid and lower AICc
aicc_only_linear <- !(is.finite(lin_rows$delta_aicc) & lin_rows$delta_aicc > 0 &
                        lin_rows$n_below_bp >= 3 & lin_rows$n_above_bp >= 3 &
                        is.finite(lin_rows$slope1) & is.finite(lin_rows$slope2))
direct <- sapply(lin_sets[1:8], function(d) fit_tr_vpd(d$vpd, d$tr, selection = "aicc")$row$model == "linear")
check(identical(unname(direct), aicc_only_linear[1:8]), "--selection aicc applies the AICc-only rule")
check(mean(lin_rows$model == "linear") >= 0.9,
      sprintf("linear truth classified linear in >= 90%% with AICc + Davies (got %.0f%%)", 100 * mean(lin_rows$model == "linear")))
check(mean(aicc_only_linear) >= 0.65,
      sprintf("linear truth classified linear in most simulations with AICc alone (got %.0f%%)", 100 * mean(aicc_only_linear)))
few <- fit_tr_vpd(vpd_steps, -10 + 30 * vpd_steps)$row
check(few$model == "linear" && grepl("not tried", few$note), "fewer than --min-points: no breakpoint model")

# 6. Bad input fails loudly -------------------------------------------------------------
check(grepl("fraction", errmsg(check_temp_rh(25, 0.5, "x"))), "RH given as a fraction is rejected")
check(grepl("Fahrenheit", errmsg(check_temp_rh(86, 50, "x"))), "temperature in Fahrenheit is rejected")
check(grepl("svp", errmsg(svp_kpa(20, "magnus"))), "unknown SVP formula is rejected")
check(grepl("timestamp", errmsg(parse_time("01/07/2026 10:00", "weights"))), "unreadable timestamps are rejected")

# 7. End to end: demo (raw logger files) and precomputed input --------------------------
out <- file.path(tempdir(), "claw-vpd-test")
unlink(out, recursive = TRUE)
res <- run_cli(c("--demo", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "tables/tr_vpd_fits.csv", "tables/step_tr_vpd.csv", "tables/weight_flags.csv",
            "tables/night_tr.csv", "figures/tr_vpd_fits_1.png", "figures/weight_qc.png",
            "reproducibility/commands.sh", "reproducibility/sessionInfo.txt")) {
  check(file.exists(file.path(out, f)), paste("demo writes", f))
}
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
fits <- utils::read.csv(file.path(out, "tables", "tr_vpd_fits.csv"), stringsAsFactors = FALSE)
check(identical(fits$model, c("segmented", "segmented", "linear")), "demo: G1 and G2 segmented, G3 linear")
check(near(fits$breakpoint_kpa[1], 2.0, 0.15) && near(fits$breakpoint_kpa[2], 2.6, 0.2),
      "demo: breakpoints near the true 2.0 and 2.6 kPa")
fl <- utils::read.csv(file.path(out, "tables", "weight_flags.csv"), stringsAsFactors = FALSE)
check(any(fl$plant_id == "P07" & fl$irrigation & fl$timestamp == "2026-07-01 10:20:00"),
      "demo: the planted irrigation on P07 at 10:20 is found")
check(any(fl$plant_id == "P02" & fl$outlier & fl$timestamp == "2026-07-01 08:41:00"),
      "demo: the planted spike on P02 at 08:41 is found")

st <- utils::read.csv(file.path(out, "tables", "step_tr_vpd.csv"), stringsAsFactors = FALSE)
# Each pot gets VPD from its own chamber's sensors (mean of sensor means), by hand:
env <- utils::read.csv(file.path(SKILL_TEST_DIR, "demo", "env.csv"), stringsAsFactors = FALSE)
env$vpd <- vpd_kpa(env$air_temp_c, env$rh_pct)
in_s4 <- env$timestamp >= "2026-07-01 10:10:00" & env$timestamp < "2026-07-01 11:00:00"
hand <- sapply(c("A", "B"), function(ch) {
  e <- env[in_s4 & env$chamber == ch, ]
  mean(tapply(e$vpd, e$sensor, mean))
})
tab_a <- st$vpd_kpa[st$plant_id == "P01" & st$step == "S4"]   # chamber A
tab_b <- st$vpd_kpa[st$plant_id == "P07" & st$step == "S4"]   # chamber B
check(near(tab_a, hand[["A"]], 1e-9) && near(tab_b, hand[["B"]], 1e-9) && abs(hand[["A"]] - hand[["B"]]) > 0.01,
      "step VPD comes from the pot's own chamber (chamber A and B differ)")
pre <- st[st$phase == "day", c("genotype", "plant_id", "vpd_kpa", "tr_mg_m2_s")]
pre_file <- file.path(tempdir(), "claw-vpd-precomputed.csv")
utils::write.csv(pre, pre_file, row.names = FALSE)
out2 <- file.path(tempdir(), "claw-vpd-test-pre")
unlink(out2, recursive = TRUE)
res2 <- run_cli(c("--input", shQuote(pre_file), "--out", shQuote(out2)))
check(res2$status == 0L, paste("precomputed run exits cleanly:", paste(res2$output, collapse = "\n")))
fits2 <- utils::read.csv(file.path(out2, "tables", "tr_vpd_fits.csv"), stringsAsFactors = FALSE)
check(identical(fits2$model, fits$model) &&
        isTRUE(all.equal(fits2$breakpoint_kpa, fits$breakpoint_kpa, tolerance = 1e-4)),
      "precomputed input reproduces the raw-data fits")

bad <- run_cli(c("--demo", "--colour", "red"))
check(bad$status != 0L, "unknown option exits with an error")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
