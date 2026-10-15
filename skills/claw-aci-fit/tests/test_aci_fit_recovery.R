# Parameter recovery for claw-aci-fit: 20 C3 curves simulated with PhotoGEA's
# own FvCB model (calculate_c3_assimilation, Bernacchi et al. temperature
# responses, Cc = Ci) with known Vcmax, J, TPU and RL, leaf temperature 24-33 C,
# Gaussian noise (A sd 0.4 umol m-2 s-1, Ci sd 2 umol/mol), a 16-point LI-6800
# sequence with two recovery points, fitted with the skill's defaults
# (DEoptim 200 generations, alpha_old fitted) through the tidy-CSV path.
# Tolerances: Vcmax within 10% on every curve and median error <= 5%; J and
# TPU within 10% where the true model has >= 3 points limited by them; RL
# within 1 umol m-2 s-1.
# Run from anywhere: Rscript tests/test_aci_fit_recovery.R

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
cores <- max(1L, min(8L, parallel::detectCores(), na.rm = TRUE))
tmp <- file.path(tempdir(), "claw-aci-fit-recovery")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE)

SP <- c(400, 300, 200, 150, 100, 75, 50, 20, 400, 400, 600, 800, 1000, 1200, 1500, 1800)
CI <- c(290, 215, 145, 110, 76, 60, 45, 25, 290, 290, 430, 575, 720, 865, 1080, 1300)

sim_c3 <- function(id, Vcmax, J, Tp, RL, tleaf, noise_a = 0.4, noise_ci = 2) {
  d <- exdf(data.frame(Cc = CI, TleafCnd = tleaf, total_pressure = 0.996, Oxygen = 21))
  d <- document_variables(d, c("", "Cc", "micromol mol^(-1)"), c("", "TleafCnd", "degrees C"),
                          c("", "total_pressure", "bar"), c("", "Oxygen", "percent"))
  d <- calculate_temperature_response(d, c3_temperature_param_bernacchi, "TleafCnd")
  m <- calculate_c3_assimilation(d, 0, 0, 0, 0, "", J, "", "", RL, Tp, Vcmax)
  w <- cbind(m[, "Wc"], m[, "Wj"], m[, "Wp"])
  w[!is.finite(w)] <- Inf
  lim <- c("Ac", "Aj", "Ap")[apply(w, 1L, which.min)]
  data.frame(curve = id, obs = seq_along(SP), co2_r_sp_umol_mol = SP, ca_umol_mol = SP * 0.97,
             ci_umol_mol = CI + stats::rnorm(length(CI), 0, noise_ci),
             a_umol_m2_s = m[, "An"] + stats::rnorm(length(CI), 0, noise_a),
             gsw_mol_m2_s = 0.3, tleaf_c = tleaf, qin_umol_m2_s = 1800, pa_kpa = 99.5,
             delta_pcham_kpa = 0.1, oxygen_pct = 21, true_limit = lim, stringsAsFactors = FALSE)
}

set.seed(42)
n_curves <- 20L
truth <- data.frame(curve = sprintf("c%02d", seq_len(n_curves)),
                    Vcmax = stats::runif(n_curves, 60, 140), tleaf = stats::runif(n_curves, 24, 33),
                    RL = stats::runif(n_curves, 0.8, 2.0), stringsAsFactors = FALSE)
truth$J <- truth$Vcmax * stats::runif(n_curves, 1.5, 2.1)
# Odd curves: TPU limits at high Ci; even curves: TPU never limits.
truth$Tp <- truth$J * ifelse(seq_len(n_curves) %% 2 == 1, stats::runif(n_curves, 0.060, 0.068), 0.09)
sim <- do.call(rbind, lapply(seq_len(n_curves), function(i)
  with(truth[i, ], sim_c3(curve, Vcmax, J, Tp, RL, tleaf))))
csv <- file.path(tmp, "sim_c3.csv")
utils::write.csv(sim, csv, row.names = FALSE)

opt <- DEFAULTS
opt$curve_cols <- "curve"
opt$cores <- cores
t0 <- Sys.time()
res <- analyse_aci(read_tidy(csv), opt)
cat(sprintf("20 simulated curves fitted in %.0f s on %d core(s)\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")), cores))
p <- res$params[match(truth$curve, res$params$curve_id), ]
used <- sim[!sim$obs %in% c(9, 10), ]  # recovery points are excluded from the fits
n_lim <- function(l) as.integer(tapply(used$true_limit == l, used$curve, sum)[truth$curve])
true_j_pts <- n_lim("Aj")
true_tp_pts <- n_lim("Ap")
rel <- function(est, tru) (est - tru) / tru

e_v <- rel(p$vcmax_25_umol_m2_s, truth$Vcmax)
cat(sprintf("Vcmax errors (%%): %s\n", paste(sprintf("%.1f", 100 * e_v), collapse = " ")))
check(all(!is.na(e_v)), "Vcmax estimated for all 20 curves")
check(max(abs(e_v)) <= 0.10, sprintf("Vcmax within 10%% on every curve (max error %.1f%%)", 100 * max(abs(e_v))))
check(stats::median(abs(e_v)) <= 0.05, sprintf("median Vcmax error <= 5%% (%.1f%%)", 100 * stats::median(abs(e_v))))
cover <- mean(truth$Vcmax >= p$vcmax_25_lower_umol_m2_s & truth$Vcmax <= p$vcmax_25_upper_umol_m2_s)
check(cover >= 0.75, sprintf("true Vcmax inside the ~95%% interval for >= 75%% of curves (%.0f%%)", 100 * cover))

jt <- true_j_pts >= 3
e_j <- rel(p$j_25_umol_m2_s, truth$J)[jt]
cat(sprintf("J errors (%%) where testable: %s\n", paste(sprintf("%.1f", 100 * e_j), collapse = " ")))
check(sum(jt) >= 8 && mean(!is.na(e_j)) >= 0.8,
      sprintf("J estimated for >= 80%% of the %d curves with >= 3 truly J-limited points", sum(jt)))
check(max(abs(e_j), na.rm = TRUE) <= 0.10, sprintf("J within 10%% (max error %.1f%%)", 100 * max(abs(e_j), na.rm = TRUE)))
tt <- true_tp_pts >= 3
e_t <- rel(p$tp_25_umol_m2_s, truth$Tp)[tt]
cat(sprintf("TPU errors (%%) where testable: %s\n", paste(sprintf("%.1f", 100 * e_t), collapse = " ")))
check(sum(tt) >= 5 && mean(!is.na(e_t)) >= 0.8,
      sprintf("TPU estimated for >= 80%% of the %d curves with >= 3 truly TPU-limited points", sum(tt)))
check(max(abs(e_t), na.rm = TRUE) <= 0.10, sprintf("TPU within 10%% (max error %.1f%%)", 100 * max(abs(e_t), na.rm = TRUE)))
check(all(is.na(p$tp_25_umol_m2_s[true_tp_pts == 0]) | grepl("tp_few", p$qc_flags[true_tp_pts == 0])),
      "TPU is NA or flagged on curves where TPU never limits")
e_rl <- p$rl_25_umol_m2_s - truth$RL
check(max(abs(e_rl)) <= 1.0, sprintf("RL within 1 umol m-2 s-1 (max error %.2f)", max(abs(e_rl))))
check(all(p$n_recovery_points == 2L) && all(p$n_points_used == 14L), "2 recovery points per curve detected and excluded")
check(all(abs(p$tleaf_mean_c - truth$tleaf) < 1e-8), "mean leaf temperature reported per curve")
check(all(res$params$seed == 1234 + seq_len(n_curves) - 1L), "per-curve seeds recorded")

cat(sprintf("PASS: %d checks (%s, parameter recovery)\n", n_checks, SKILL))
