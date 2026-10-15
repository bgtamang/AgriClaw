# Tests for claw-vi-timeseries. Run from anywhere: Rscript tests/test_vi_timeseries.R
# Known-answer checks: model formulas against phenofit, metric definitions on an
# analytic curve, recovery of simulated double-logistic curves with known
# parameters, outlier handling, sparse-series refusal, axis alignment, input
# validation, and an end-to-end run of the command line on the demo data.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "vi_timeseries.R"))  # loads functions; main() does not run

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}

run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "vi_timeseries.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}

err_msg <- function(expr) tryCatch({ force(expr); "" }, error = function(e) conditionMessage(e))

opt0 <- DEFAULTS
with_opt <- function(...) { o <- opt0; v <- list(...); o[names(v)] <- v; o }

# Simulated crop seasons: Beck curves with random but plausible parameters.
sim_params <- function(i) {
  set.seed(1000 + i)
  c(vmin = runif(1, 0.15, 0.25), vmax = runif(1, 0.75, 0.92), sos = runif(1, 32, 45),
    rsp = runif(1, 0.10, 0.22), eos = runif(1, 85, 105), rau = runif(1, 0.08, 0.20))
}
sim_dates <- function(i) {  # 12-16 irregular dates, DAP 12-128 (UAV-like)
  set.seed(2000 + i)
  n <- sample(12:16, 1)
  sort(unique(round(seq(15, 125, length.out = n) + runif(n, -3, 3))))
}
true_metrics <- function(f, t) {
  grid <- seq(min(t), max(t), length.out = GRID_N)
  curve_metrics(grid, f(grid))
}
metric_errors <- function(a, f, t) {
  tm <- true_metrics(f, t)
  m <- a$metrics
  c(peak_d = m$peak_t - tm$peak_t, peak_v = m$peak_value - tm$peak_value,
    gu = m$greenup_t - tm$greenup_t, sen = m$senescence_t - tm$senescence_t,
    auc_pct = 100 * (m$auc - tm$auc) / tm$auc)
}

# 1. Model formulas equal phenofit's reference implementations -------------------
tt <- seq(0, 150, by = 0.5)
pb <- c(vmin = 0.2, vmax = 0.85, sos = 40, rsp = 0.15, eos = 95, rau = 0.12)
pe <- c(base = 0.2, amp = 0.75, sos = 40, rsp = 0.15, eos = 95, rau = 0.12, m7 = 0.001)
if (requireNamespace("phenofit", quietly = TRUE)) {
  check(max(abs(dl_beck(tt, pb) - phenofit::doubleLog.Beck(unname(pb), tt))) < 1e-12,
        "dl_beck equals phenofit::doubleLog.Beck (Beck et al. 2006)")
  check(max(abs(dl_elmore(tt, pe) - phenofit::doubleLog.Elmore(unname(pe), tt))) < 1e-12,
        "dl_elmore equals phenofit::doubleLog.Elmore (Elmore et al. 2012)")
} else {
  cat("note: phenofit not installed; formula cross-check skipped\n")
}
# reparameterised fitting forms reproduce the published forms
u <- (tt - 10) / 120
check(max(abs(beck_u(u, 0.2, 0.65, (40 - 10) / 120, 0.15 * 120, 55 / 120, 0.12 * 120) -
                dl_beck(tt, pb))) < 1e-12, "beck_u on the unit axis equals dl_beck")

# 2. Metric definitions on an analytic curve -----------------------------------------
# f(t) = 0.2 + 0.6 exp(-(t - 50)^2 / (2 * 15^2)) on [0, 100]
g <- function(t) 0.2 + 0.6 * exp(-(t - 50)^2 / 450)
grid <- seq(0, 100, length.out = 20001)
m <- curve_metrics(grid, g(grid), threshold = 0.5)
base <- g(0)
thr <- base + 0.5 * (0.8 - base)
cross <- sqrt(-450 * log((thr - 0.2) / 0.6))
check(abs(m$peak_t - 50) < 1e-6 && abs(m$peak_value - 0.8) < 1e-9, "peak of the analytic curve")
check(abs(m$greenup_t - (50 - cross)) < 1e-3, "green-up = analytic 50% crossing")
check(abs(m$senescence_t - (50 + cross)) < 1e-3, "senescence = analytic 50% crossing")
check(abs(m$season_length - 2 * cross) < 2e-3, "season length = senescence - green-up")
max_slope <- 0.6 / 15 * exp(-0.5)            # at t = 50 -/+ 15
check(abs(m$greenup_rate - max_slope) < 1e-5 && abs(m$senescence_rate - max_slope) < 1e-5,
      "green-up and senescence rates = analytic maximum slopes")
auc_true <- 0.2 * 100 + 0.6 * 15 * sqrt(2 * pi) * (pnorm(50 / 15) - pnorm(-50 / 15))
check(abs(m$auc - auc_true) / auc_true < 1e-6, "AUC = analytic integral")
np <- sqrt(-450 * log(1 - 0.05 * (0.8 - base) / 0.6))
check(abs(m$near_peak_end_t - m$near_peak_start_t - 2 * np) < 0.02, "near-peak window width")

# A window that ends before senescence: no senescence value is invented.
grid2 <- seq(0, 60, length.out = 6001)
m2 <- curve_metrics(grid2, g(grid2))
check(is.na(m2$senescence_t) && is.na(m2$senescence_rate) && is.finite(m2$greenup_t),
      "truncated window: senescence NA, green-up still reported")
check(any(grepl("senescence not reached", m2$notes)), "truncated window explains the NA")

# 3. Double-logistic recovery on 20 simulated seasons --------------------------------
errs <- t(vapply(1:20, function(i) {
  p <- sim_params(i)
  t <- sim_dates(i)
  set.seed(3000 + i)
  y <- dl_beck(t, p) + rnorm(length(t), 0, 0.02)
  a <- analyse_plot(t, y, rep(1, length(t)), with_opt(method = "dl-beck"))
  metric_errors(a, function(x) dl_beck(x, p), t)
}, numeric(5)))
mx <- apply(abs(errs), 2, max)
check(mx[["peak_d"]] <= 4, sprintf("dl-beck peak date within 4 d (max %.2f)", mx[["peak_d"]]))
check(mx[["peak_v"]] <= 0.02, sprintf("dl-beck peak value within 0.02 (max %.3f)", mx[["peak_v"]]))
check(mx[["gu"]] <= 2.5, sprintf("dl-beck green-up within 2.5 d (max %.2f)", mx[["gu"]]))
check(mx[["sen"]] <= 2.5, sprintf("dl-beck senescence within 2.5 d (max %.2f)", mx[["sen"]]))
check(mx[["auc_pct"]] <= 2.5, sprintf("dl-beck AUC within 2.5%% (max %.2f)", mx[["auc_pct"]]))

errs_e <- t(vapply(1:10, function(i) {
  p <- sim_params(i)
  t <- sim_dates(i)
  pe <- c(base = p[["vmin"]], amp = 1.15 * (p[["vmax"]] - p[["vmin"]]) + 0.09, sos = p[["sos"]],
          rsp = p[["rsp"]], eos = p[["eos"]], rau = p[["rau"]], m7 = 0.0015)
  set.seed(3000 + i)
  y <- dl_elmore(t, pe) + rnorm(length(t), 0, 0.02)
  a <- analyse_plot(t, y, rep(1, length(t)), with_opt(method = "dl-elmore"))
  metric_errors(a, function(x) dl_elmore(x, pe), t)
}, numeric(5)))
mx <- apply(abs(errs_e), 2, max)
check(mx[["gu"]] <= 2.5 && mx[["sen"]] <= 2.5,
      sprintf("dl-elmore green-up/senescence within 2.5 d (max %.2f / %.2f)", mx[["gu"]], mx[["sen"]]))
check(mx[["auc_pct"]] <= 2.5, sprintf("dl-elmore AUC within 2.5%% (max %.2f)", mx[["auc_pct"]]))

for (meth in c("loess", "sg")) {
  e <- t(vapply(1:10, function(i) {
    p <- sim_params(i)
    t <- sim_dates(i)
    set.seed(3000 + i)
    y <- dl_beck(t, p) + rnorm(length(t), 0, 0.02)
    metric_errors(analyse_plot(t, y, rep(1, length(t)), with_opt(method = meth)),
                  function(x) dl_beck(x, p), t)
  }, numeric(5)))
  mx <- apply(abs(e), 2, max)
  check(mx[["gu"]] <= 3.5 && mx[["sen"]] <= 3.5,
        sprintf("%s green-up/senescence within 3.5 d (max %.2f / %.2f)", meth, mx[["gu"]], mx[["sen"]]))
  check(mx[["auc_pct"]] <= 3, sprintf("%s AUC within 3%% (max %.2f)", meth, mx[["auc_pct"]]))
}

# 4. Sparse and one-sided series are refused for the double logistic -----------------
p <- sim_params(1)
t5 <- c(20, 45, 65, 90, 110)
a <- analyse_plot(t5, dl_beck(t5, p), rep(1, 5), with_opt(method = "dl-beck"))
check(a$status == "sparse" && is.null(a$metrics), "5 dates: double logistic refused")
check(a$simple$n_dates == 5 && abs(a$simple$max_obs_value - max(dl_beck(t5, p))) < 1e-12,
      "5 dates: simple metrics still reported")
t6 <- c(20, 40, 55, 70, 90, 110)
check(analyse_plot(t6, dl_beck(t6, p), rep(1, 6), with_opt(method = "dl-beck"))$status == "sparse",
      "6 dates: 6-parameter Beck model refused (needs 7)")
check(analyse_plot(t6, dl_beck(t6, p), rep(1, 6), with_opt(method = "loess"))$status == "fitted",
      "6 dates: LOESS allowed")
check(grepl("at least 7", err_msg(min_dates_for("dl-beck", 5))), "--min-dates below p + 1 rejected")
t_rise <- seq(15, 60, by = 5)
a <- analyse_plot(t_rise, dl_beck(t_rise, p), rep(1, length(t_rise)), with_opt(method = "dl-beck"))
check(a$status == "limb_missing", "rise only: double logistic refused (limb missing)")

# 5. Cloud-like drops on satellite-like series: upper envelope removes the bias -------
sat <- function(i, weighting) {
  p <- sim_params(i)
  t <- seq(10, 130, by = 5)
  set.seed(4000 + i)
  y <- dl_beck(t, p) + rnorm(length(t), 0, 0.015)
  k <- sample(3:(length(t) - 2), 6)
  y[k] <- y[k] - runif(6, 0.1, 0.4)
  a <- analyse_plot(t, y, rep(1, length(t)), with_opt(weighting = weighting))
  c(metric_errors(a, function(x) dl_beck(x, p), t), caught = mean(a$weights[k] < 0.2),
    clean_w = mean(a$weights[-k]))
}
e_none <- t(vapply(1:10, sat, numeric(7), weighting = "none"))
e_up <- t(vapply(1:10, sat, numeric(7), weighting = "upper"))
check(mean(e_none[, "auc_pct"]) < -5, sprintf("unweighted fit is biased low by clouds (mean AUC error %.1f%%)", mean(e_none[, "auc_pct"])))
check(mean(abs(e_up[, "auc_pct"])) < 2 && max(abs(e_up[, "auc_pct"])) < 6,
      sprintf("upper envelope: mean |AUC error| %.2f%%, max %.2f%%", mean(abs(e_up[, "auc_pct"])), max(abs(e_up[, "auc_pct"]))))
check(mean(abs(e_up[, "peak_v"])) < 0.02, "upper envelope: peak value recovered")
check(mean(e_up[, "caught"]) >= 0.8 && mean(e_up[, "clean_w"]) >= 0.9,
      sprintf("upper envelope: %.0f%% of cloud points weighted < 0.2; clear dates keep mean weight %.2f",
              100 * mean(e_up[, "caught"]), mean(e_up[, "clean_w"])))

# 6. Upward outlier on the plateau: upper trusts it, robust rejects it ---------------
up <- function(i, weighting) {
  p <- sim_params(i)
  t <- round(seq(15, 125, length.out = 15))
  set.seed(5000 + i)
  y <- dl_beck(t, p) + rnorm(length(t), 0, 0.015)
  k <- which.min(abs(t - (p[["sos"]] + p[["eos"]]) / 2))
  y[k] <- y[k] + 0.15
  a <- analyse_plot(t, y, rep(1, length(t)), with_opt(weighting = weighting))
  c(metric_errors(a, function(x) dl_beck(x, p), t), w = a$weights[k])
}
u_up <- t(vapply(1:10, up, numeric(6), weighting = "upper"))
u_rb <- t(vapply(1:10, up, numeric(6), weighting = "robust"))
check(all(u_up[, "w"] > 0.99), "upper envelope keeps full weight on an upward outlier (documented)")
check(all(u_rb[, "w"] < 0.05), "robust weighting rejects the upward outlier")
check(max(abs(u_rb[, "peak_v"])) <= 0.03 && mean(u_up[, "peak_v"]) > 0.05,
      "robust recovers the peak value; upper is pulled up by the outlier")

# 7. Through the full pipeline: flags, date diagnostics, axes ------------------------
make_obs <- function(curves, planting = "2025-05-01", shift_date = NULL, shift = 0) {
  rows <- lapply(seq_along(curves), function(i) {
    cv <- curves[[i]]
    data.frame(plot_id = sprintf("P%02d", i), date = format(as.Date(planting) + cv$t),
               ndvi = cv$y, stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, rows)
  if (!is.null(shift_date)) d$ndvi[d$date == shift_date] <- d$ndvi[d$date == shift_date] + shift
  d
}
pipeline <- function(dat, ...) {
  o <- with_opt(...)
  prep <- prepare_series(dat, o, if (nzchar(o$gdd)) read_csv_checked(o$gdd) else NULL)
  run_timeseries(prep, o)
}
tdates <- round(seq(15, 125, length.out = 15))
curves <- lapply(1:8, function(i) {
  p <- sim_params(i)
  set.seed(6000 + i)
  list(t = tdates, y = dl_beck(tdates, p) + rnorm(length(tdates), 0, 0.015))
})
k_peak <- which.min(abs(tdates - 70))
curves[[2]]$y[k_peak] <- curves[[2]]$y[k_peak] + 0.15   # upward error
curves[[5]]$y[k_peak] <- curves[[5]]$y[k_peak] - 0.25   # shadow
dat <- make_obs(curves)
for (w in c("none", "upper")) {
  r <- pipeline(dat, planting_date = "2025-05-01", weighting = w)
  ob <- r$observations
  hit_hi <- ob$outlier[ob$plot_id == "P02" & ob$dap == tdates[k_peak]]
  hit_lo <- ob$outlier[ob$plot_id == "P05" & ob$dap == tdates[k_peak]]
  check(hit_hi == "high" && hit_lo == "low",
        sprintf("weighting %s: upward error flagged high and shadow flagged low", w))
  check(sum(ob$outlier != "") <= 4, sprintf("weighting %s: few false flags (%d flags)", w, sum(ob$outlier != "")))
  check(abs(r$noise_sd - 0.015) < 0.005, sprintf("noise SD estimate %.4f near the true 0.015", r$noise_sd))
}
# a whole-flight calibration shift is detected as a date problem
shift_day <- format(as.Date("2025-05-01") + tdates[6])
r <- pipeline(make_obs(curves, shift_date = shift_day, shift = -0.08), planting_date = "2025-05-01")
dd <- date_diagnostics(r$observations, r$noise_sd)
check(identical(dd$date[dd$flag], shift_day), "whole-flight shift flagged on the right date only")

# never extrapolate: curves stay inside each plot's observed window
cv <- r$curves
rng <- tapply(r$observations$dap, r$observations$plot_id, range)
check(all(vapply(names(rng), function(p) {
  x <- cv$dap[cv$plot_id == p]
  min(x) >= rng[[p]][1] && max(x) <= rng[[p]][2]
}, TRUE)), "fitted curves never leave the observed window")
r_fa <- pipeline(dat, planting_date = "2025-05-01", features_at = "10,60,130")
check(all(is.na(r_fa$metrics$ndvi_at_dap10)) && all(is.na(r_fa$metrics$ndvi_at_dap130)) &&
        all(is.finite(r_fa$metrics$ndvi_at_dap60)), "--features-at gives NA outside the window")

# planting-date alignment: same season shape, planted three weeks apart
one <- curves[1]
d_early <- make_obs(one, planting = "2025-05-01")
d_late <- make_obs(one, planting = "2025-05-22")
d_late$plot_id <- "LATE"
two <- rbind(transform(d_early, planting_date = "2025-05-01"),
             transform(d_late, planting_date = "2025-05-22"))
r_dap <- pipeline(two)
md <- r_dap$metrics
check(abs(diff(md$peak_dap)) < 1e-8 && abs(diff(md$greenup_dap)) < 1e-8,
      "DAP axis: identical metrics for identical seasons planted on different dates")
check(abs(as.numeric(diff(as.Date(md$peak_date)))) == 21, "calendar peak dates differ by the 21-day planting gap")
r_doy <- pipeline(two, axis = "doy")
check(abs(abs(diff(r_doy$metrics$peak_doy)) - 21) < 1e-6, "DOY axis shifts the same season by 21 days")

# GDD axis: 10 GDD every day makes GDD = 10 x DAP, so metrics scale by 10
gdd_file <- tempfile(fileext = ".csv")
utils::write.csv(data.frame(date = format(seq(as.Date("2025-04-25"), as.Date("2025-10-01"), by = "day")),
                            gdd_c = 10), gdd_file, row.names = FALSE)
r_g <- pipeline(dat, planting_date = "2025-05-01", axis = "gdd", gdd = gdd_file)
r_d <- pipeline(dat, planting_date = "2025-05-01")
check(max(abs(r_g$metrics$peak_gdd - 10 * r_d$metrics$peak_dap)) < 1e-6 &&
        max(abs(r_g$metrics$senescence_gdd - 10 * r_d$metrics$senescence_dap)) < 1e-6,
      "GDD axis: metrics = 10 x DAP metrics for 10 GDD per day")
check(identical(r_g$metrics$peak_date, r_d$metrics$peak_date), "GDD axis maps back to the same calendar dates")
gap_file <- tempfile(fileext = ".csv")
gd <- utils::read.csv(gdd_file)
utils::write.csv(gd[-40, ], gap_file, row.names = FALSE)
check(grepl("missing day", err_msg(pipeline(dat, planting_date = "2025-05-01", axis = "gdd", gdd = gap_file))),
      "a gap in the GDD series is an error, not silently filled")

# 8. Input validation fails loudly ----------------------------------------------------
check(grepl("never guess", err_msg(pipeline(dat))), "missing planting date: asks for it")
dup <- rbind(dat, dat[1, ])
check(grepl("repeated plot x date", err_msg(pipeline(dup, planting_date = "2025-05-01"))), "duplicate plot x date rejected")
scaled <- dat
scaled$ndvi <- scaled$ndvi * 10000
check(grepl("scaled", err_msg(pipeline(scaled, planting_date = "2025-05-01"))), "NDVI x 10000 rejected")
check(grepl("before planting", err_msg(pipeline(dat, planting_date = "2025-06-01"))), "observations before planting rejected")
check(grepl("not found", err_msg(pipeline(dat, planting_date = "2025-05-01", vi = "evi"))), "missing VI column named")

# 9. End to end on the demo data ------------------------------------------------------
out <- file.path(tempdir(), "claw-vi-timeseries-test")
unlink(out, recursive = TRUE)
res <- run_cli(c("--demo", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "tables/metrics.csv", "tables/observations.csv", "tables/fitted_curves.csv",
            "tables/dl_parameters.csv", "tables/date_diagnostics.csv", "tables/feature_time.csv",
            "figures/curves_examples.png", "figures/curves_all.png", "figures/metric_distributions.png",
            "reproducibility/commands.sh", "reproducibility/sessionInfo.txt")) {
  check(file.exists(file.path(out, f)), paste("demo writes", f))
}
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
met <- utils::read.csv(file.path(out, "tables", "metrics.csv"))
check(nrow(met) == 24L && met$status[met$plot_id == "E2-12"] == "sparse", "demo: 24 plots, E2-12 sparse")
ob <- utils::read.csv(file.path(out, "tables", "observations.csv"))
check(ob$outlier[ob$plot_id == "E1-03" & ob$dap == 52] == "low" &&
        ob$outlier[ob$plot_id == "E2-07" & ob$dap == 71] == "low", "demo: both planted shadow drops flagged")
rep_txt <- readLines(file.path(out, "report.md"))
check(any(grepl("## Warnings", rep_txt)) && any(grepl("E2-12", rep_txt)), "report lists the sparse plot under Warnings")
out_g <- file.path(tempdir(), "claw-vi-timeseries-test-gdd")
unlink(out_g, recursive = TRUE)
res_g <- run_cli(c("--demo", "--axis", "gdd", "--method", "loess", "--out", shQuote(out_g)))
check(res_g$status == 0L && file.exists(file.path(out_g, "tables", "metrics.csv")),
      "demo on the GDD axis with LOESS runs")
check(run_cli(c("--demo", "--colour", "red"))$status != 0L, "unknown option exits with an error")
check(any(grepl("--weighting", run_cli("--help")$output)), "--help prints usage")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
