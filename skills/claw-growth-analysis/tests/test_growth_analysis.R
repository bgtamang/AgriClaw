# Tests for claw-growth-analysis. Run from anywhere: Rscript tests/test_growth_analysis.R
# Known-answer checks: hand-calculated RGR, NAR, SLA, LAR and LMF; exact NAR
# when mass is linear in leaf area; recovery of known RGR from simulated
# exponential and logistic growth (classical and functional); the
# ln-of-means bias described by Hoffmann & Poorter (2002); bad input; and an
# end-to-end run of the command line on the demo data.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "growth_analysis.R"))  # loads functions; main() does not run

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
near <- function(a, b, tol) isTRUE(all(abs(a - b) <= tol))

run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "growth_analysis.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
errmsg <- function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))
sim_frame <- function(time, mass) data.frame(group = "S", time = time, mass = mass, la = NA_real_, lm = NA_real_)

# 1. Hand calculations (two plants per harvest, values worked out in Python) ----
raw <- data.frame(genotype = "A", days = c(0, 0, 10, 10), total_dry_mass_g = c(1, 2, 4, 8),
                  leaf_area_cm2 = c(100, 200, 300, 500), leaf_dry_mass_g = c(0.5, 1, 2, 3))
d <- prepare_data(raw)
cl <- classical_intervals(d, n_boot = 0)
check(near(cl$rgr_g_g_d, 0.138629436112, 1e-10), "RGR = (mean ln W2 - mean ln W1)/dt = ln(4)/10")
check(near(cl$rgr_se, 0.0490129071734, 1e-10), "RGR SE = sqrt(s1^2/n1 + s2^2/n2)/dt")
check(near(cl$nar_g_m2_d, 17.6549265542, 1e-8), "NAR = (W2-W1)/dt x (ln LA2 - ln LA1)/(LA2 - LA1), in g m-2 d-1")
hs <- harvest_summary(d)
check(near(hs$sla_m2_kg, c(20, 15.8333333333), 1e-8), "SLA (m2/kg) = mean of per-plant leaf area / leaf mass x 0.1")
check(near(hs$lar_m2_kg, c(10, 6.875), 1e-10), "LAR (m2/kg) = mean of per-plant leaf area / total mass x 0.1")
check(near(hs$lmf_g_g, c(0.5, 0.4375), 1e-12), "LMF = mean of per-plant leaf mass / total mass")
check(near(hs$mass_geomean_g, c(sqrt(2), sqrt(32)), 1e-12), "geometric mean mass = exp(mean ln W)")
check(near(0.1 * d$la / d$mass, (0.1 * d$la / d$lm) * (d$lm / d$mass), 1e-12), "per plant LAR = SLA x LMF")

# 2. NAR is exact when mass is linear in leaf area -------------------------------
# LA(t) = 50 exp(0.1 t) cm2, W = 0.2 + 0.004 LA g: instantaneous NAR = (dW/dt)/LA,
# its mean over [t1, t2] = 0.004 x 0.1 g cm-2 d-1 = 4 g m-2 d-1.
la_t <- function(t) 50 * exp(0.1 * t)
w_t <- function(t) 0.2 + 0.004 * la_t(t)
raw2 <- data.frame(genotype = "N", days = rep(c(5, 12), each = 3),
                   total_dry_mass_g = w_t(rep(c(5, 12), each = 3)), leaf_area_cm2 = la_t(rep(c(5, 12), each = 3)))
cl2 <- classical_intervals(prepare_data(raw2, leaf_mass_col = ""), n_boot = 0)
inst <- function(t) 0.004 * 0.1 * la_t(t) / la_t(t)
num <- stats::integrate(inst, 5, 12)$value / 7 * 1e4
check(near(cl2$nar_g_m2_d, 4, 1e-9) && near(num, 4, 1e-9), "NAR equals the exact interval mean (4 g m-2 d-1)")

# 3. Exponential growth with known RGR ----------------------------------------------
set.seed(1)
r <- 0.12
days <- c(0, 7, 14, 21, 28)
est <- cover <- f_err <- f_cov <- ov <- c()
for (s in 1:80) {                      # 80 data sets; the spline fit on the first 25
  dd <- sim_frame(rep(days, each = 8), 0)
  dd$mass <- 0.05 * exp(r * dd$time + rnorm(nrow(dd), 0, 0.25))
  ci <- classical_intervals(dd, n_boot = 0)
  est <- c(est, ci$rgr_g_g_d); cover <- c(cover, ci$rgr_lo <= r & ci$rgr_hi >= r)
  ov <- c(ov, overall_rgr(dd)$rgr_g_g_d)
  if (s > 25) next
  fr <- functional_rgr(dd)$curve
  fh <- fr[fr$at_harvest, ]
  f_err <- c(f_err, fh$rgr_g_g_d - r); f_cov <- c(f_cov, fh$rgr_lo <= r & fh$rgr_hi >= r)
}
check(near(mean(est), r, 0.003), sprintf("classical RGR unbiased for exponential growth (mean %.4f vs 0.12)", mean(est)))
check(mean(cover) >= 0.9, sprintf("classical 95%% CI covers the true RGR in >= 90%% of intervals (got %.0f%%)", 100 * mean(cover)))
check(max(abs(tapply(f_err, rep(seq_along(days), 25), mean))) < 0.006, "functional RGR(t) unbiased (within 0.006) at every harvest time")
check(mean(f_cov) >= 0.85, sprintf("functional 95%% band covers the true RGR >= 85%% (got %.0f%%)", 100 * mean(f_cov)))
check(near(mean(ov), r, 0.003), "overall slope of ln W on time recovers the constant RGR")

# 4. Logistic growth: RGR(t) = r (1 - W/K) declines with time --------------------
set.seed(2)
K <- 20; w0 <- 0.05; rr <- 0.25
W <- function(t) K / (1 + ((K - w0) / w0) * exp(-rr * t))
true_rgr <- function(t) rr * (1 - W(t) / K)
days_l <- seq(0, 40, 5)
abs_err <- covl <- c(); cl_err <- cl_cov <- c()
for (s in 1:20) {
  dd <- sim_frame(rep(days_l, each = 6), 0)
  dd$mass <- W(dd$time) * exp(rnorm(nrow(dd), 0, 0.2))
  fr <- functional_rgr(dd)$curve
  fh <- fr[fr$at_harvest, ]
  abs_err <- c(abs_err, abs(fh$rgr_g_g_d - true_rgr(fh$time)))
  covl <- c(covl, fh$rgr_lo <= true_rgr(fh$time) & fh$rgr_hi >= true_rgr(fh$time))
  ci <- classical_intervals(dd, n_boot = 0)
  truth_int <- (log(W(ci$t2)) - log(W(ci$t1))) / (ci$t2 - ci$t1)   # exact interval-mean RGR
  cl_err <- c(cl_err, ci$rgr_g_g_d - truth_int); cl_cov <- c(cl_cov, ci$rgr_lo <= truth_int & ci$rgr_hi >= truth_int)
}
check(mean(abs_err) < 0.015, sprintf("functional RGR(t) within 0.015 g/g/d of the logistic truth on average (got %.4f)", mean(abs_err)))
check(mean(covl) >= 0.85, sprintf("functional band covers the logistic RGR(t) >= 85%% (got %.0f%%)", 100 * mean(covl)))
check(abs(mean(cl_err)) < 0.003 && mean(cl_cov) >= 0.9, "classical interval RGR unbiased for the logistic interval means")
check(true_rgr(0) > true_rgr(40) && all(diff(true_rgr(days_l)) < 0), "logistic truth has declining RGR (sanity)")

# 5. Hoffmann & Poorter (2002): ln of mean mass is biased when var(ln W) changes --
# True RGR 0.1 over 7 days; SD of ln W 0.1 at harvest 1 and 0.6 at harvest 2; n = 10.
set.seed(3)
ub <- bi <- numeric(600)
for (s in seq_along(ub)) {
  dd <- sim_frame(rep(c(0, 7), each = 10), c(exp(rnorm(10, 0, 0.1)), exp(0.7 + rnorm(10, 0, 0.6))))
  ci <- classical_intervals(dd, n_boot = 0)
  ub[s] <- ci$rgr_g_g_d; bi[s] <- ci$rgr_ln_of_means_g_g_d
}
pred <- ((0.36 / 2 - (exp(0.36) - 1) / 20) - (0.01 / 2 - (exp(0.01) - 1) / 20)) / 7   # delta-method bias
check(near(mean(ub), 0.1, 0.004), sprintf("the skill's estimator (mean of ln W) is unbiased (mean %.4f)", mean(ub)))
check(mean(bi) - 0.1 > 0.015, sprintf("ln of mean W overestimates RGR by > 0.015 (got %.4f)", mean(bi) - 0.1))
check(near(mean(bi) - 0.1, pred, 0.004), sprintf("that bias matches the predicted %.4f g/g/d", pred))

# 6. Bad input fails loudly ------------------------------------------------------------
check(grepl("> 0", errmsg(prepare_data(transform(raw, total_dry_mass_g = c(0, 2, 4, 8))))), "zero mass is rejected")
check(grepl("biomass_g", errmsg(prepare_data(raw, mass_col = "biomass_g"))), "missing column is named")
check(grepl("two harvest", errmsg(prepare_data(raw[raw$days == 0, ]))), "a single harvest is rejected")
check(grepl("Leaf dry mass exceeds", errmsg(prepare_data(transform(raw, leaf_dry_mass_g = c(5, 1, 2, 3))))),
      "leaf mass above total mass is rejected")

# 7. End to end on the demo data -----------------------------------------------------------
out <- file.path(tempdir(), "claw-growth-test")
unlink(out, recursive = TRUE)
res <- run_cli(c("--demo", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "tables/harvest_summary.csv", "tables/classical_intervals.csv", "tables/overall_rgr.csv",
            "tables/functional_rgr.csv", "figures/ln_mass.png", "figures/rgr.png", "figures/rgr_vs_size.png",
            "figures/allocation.png", "reproducibility/commands.sh", "reproducibility/sessionInfo.txt")) {
  check(file.exists(file.path(out, f)), paste("demo writes", f))
}
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
fr <- utils::read.csv(file.path(out, "tables", "functional_rgr.csv"))
fa <- fr[fr$group == "A" & fr$at_harvest, ]
truth_a <- 0.20 * (1 - (40 / (1 + ((40 - 0.08) / 0.08) * exp(-0.20 * fa$time))) / 40)
check(mean(fa$rgr_lo <= truth_a & fa$rgr_hi >= truth_a) >= 0.8, "demo: genotype A bands cover the true logistic RGR(t)")
bad <- run_cli(c("--demo", "--colour", "red"))
check(bad$status != 0L, "unknown option exits with an error")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
