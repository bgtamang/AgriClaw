# Tests for claw-aci-fit. Run from anywhere: Rscript tests/test_aci_fit.R
# Known answers: recovery-point detection and the induction ratio on
# hand-built curves (exact values), limiting-process logic, input validation,
# exact recovery from a noise-free curve simulated with PhotoGEA's own model,
# and an end-to-end --demo run whose generating parameters are known.
# Other test files: test_aci_fit_recovery.R (20 noisy simulated curves) and
# test_aci_fit_photogea.R (PhotoGEA's published article tables).

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

run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "aci_fit.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
cores <- max(1L, min(4L, parallel::detectCores(), na.rm = TRUE))
tmp <- file.path(tempdir(), "claw-aci-fit-test")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE)
rel <- function(est, tru) (est - tru) / tru

SP <- c(400, 300, 200, 150, 100, 75, 50, 20, 400, 400, 600, 800, 1000, 1200, 1500, 1800)
CI <- c(290, 215, 145, 110, 76, 60, 45, 25, 290, 290, 430, 575, 720, 865, 1080, 1300)

# 1. Recovery points and the induction check on hand-built curves ----------------
mk_curve <- function(a_return, gsw_top = 0.25) {
  data.frame(curve_id = "x", co2_r_sp_umol_mol = SP, ca_umol_mol = SP,
             a_umol_m2_s = c(20, 16, 10, 7, 3, 1, -1, -3, a_return - 1, a_return, 30, 35, 38, 39, 39, 39),
             ci_umol_mol = CI, gsw_mol_m2_s = c(0.40, rep(0.42, 14), gsw_top))
}
m <- mark_recovery(mk_curve(16))
check(identical(which(m$recovery), c(9L, 10L)), "recovery points = the repeated 400 points (9 and 10)")
check(isTRUE(all.equal(m$ratio, 0.8)), "induction ratio = A(last recovery) / A(initial) = 16/20")
fake_params <- data.frame(curve_id = "x", rmse_umol_m2_s = 0.5, rl_25_umol_m2_s = 1, n_ac_limiting = 6,
                          n_aj_limiting = 5, n_ap_limiting = 3, fit_error = "", stringsAsFactors = FALSE)
qc_for <- function(a_return, gsw_top = 0.25) {
  d <- mk_curve(a_return, gsw_top)
  mm <- mark_recovery(d)
  d$used_in_fit <- !mm$recovery
  ind <- data.frame(curve_id = "x", induction_ratio = mm$ratio, stringsAsFactors = FALSE)
  curve_qc(d, fake_params, ind, DEFAULTS)$qc_flags
}
check(grepl("induction_drop", qc_for(16)), "A falling 20% on return to reference CO2 is flagged")
check(grepl("induction_rise", qc_for(24)), "A rising 20% on return is flagged as incomplete initial induction")
check(!grepl("induction", qc_for(21)), "a 5% change is within tolerance")
check(grepl("stomatal_closure_high_co2", qc_for(20, gsw_top = 0.15)), "gsw at the highest CO2 < 0.5 x initial is flagged")
check(!grepl("stomatal_closure", qc_for(20, gsw_top = 0.30)), "moderate gsw decline is not flagged")
no_sp <- mk_curve(16)
no_sp$co2_r_sp_umol_mol <- NULL
check(identical(which(mark_recovery(no_sp)$recovery), c(9L, 10L)),
      "recovery detection falls back to Ca when no setpoint column")
mono <- mk_curve(16)
mono$co2_r_sp_umol_mol <- sort(SP)
check(!any(mark_recovery(mono)$recovery) && is.na(mark_recovery(mono)$ratio),
      "a one-direction curve has no recovery points and no induction ratio")

# Limiting process must come from the carboxylation rates, not Ac/Aj/Ap.
lp <- limiting_process(data.frame(Wc = c(5, 30), Wj = c(8, 20), Wp = c(Inf, 25)))
check(identical(lp, c("Ac", "Aj")), "limiting process = smallest of Wc, Wj, Wp")

# 2. Exact recovery from a noise-free curve (PhotoGEA's own model) ------------------
sim_c3 <- function(id, Vcmax, J, Tp, RL, tleaf) {
  d <- exdf(data.frame(Cc = CI, TleafCnd = tleaf, total_pressure = 0.996, Oxygen = 21))
  d <- document_variables(d, c("", "Cc", "micromol mol^(-1)"), c("", "TleafCnd", "degrees C"),
                          c("", "total_pressure", "bar"), c("", "Oxygen", "percent"))
  d <- calculate_temperature_response(d, c3_temperature_param_bernacchi, "TleafCnd")
  a <- calculate_c3_assimilation(d, 0, 0, 0, 0, "", J, "", "", RL, Tp, Vcmax)[, "An"]
  data.frame(curve = id, obs = seq_along(SP), co2_r_sp_umol_mol = SP, ca_umol_mol = SP * 0.97,
             ci_umol_mol = CI, a_umol_m2_s = a, gsw_mol_m2_s = 0.3, tleaf_c = tleaf,
             qin_umol_m2_s = 1800, pa_kpa = 99.5, delta_pcham_kpa = 0.1, oxygen_pct = 21)
}
# nf_tpu: J far above need (never limits); nf_j: TPU far above need (never limits).
nf <- rbind(sim_c3("nf_tpu", Vcmax = 90, J = 250, Tp = 11, RL = 1.5, tleaf = 30),
            sim_c3("nf_j", Vcmax = 80, J = 140, Tp = 20, RL = 1.2, tleaf = 27))
extra <- nf[nf$curve == "nf_tpu", ][16, ]   # a 17th logged point with A missing
extra$obs <- 17
extra$a_umol_m2_s <- NA
nf <- rbind(nf, extra)
o <- DEFAULTS
o$curve_cols <- "curve"
o$cores <- min(2L, cores)
r_nf <- analyze_aci(nf, o)$params
a1 <- r_nf[r_nf$curve_id == "nf_tpu", ]
a2 <- r_nf[r_nf$curve_id == "nf_j", ]
check(a1$n_points == 17L && a1$n_points_used == 14L,
      "a point with missing A is excluded (and reported), as are the 2 recovery points")
check(abs(rel(a1$vcmax_25_umol_m2_s, 90)) < 0.002 && abs(rel(a1$tp_25_umol_m2_s, 11)) < 0.002 &&
      abs(a1$rl_25_umol_m2_s - 1.5) < 0.01 && is.na(a1$j_25_umol_m2_s),
      "noise-free TPU-limited curve: Vcmax, TPU, RL within 0.2%; J NA because it never limits")
check(abs(rel(a2$vcmax_25_umol_m2_s, 80)) < 0.002 && abs(rel(a2$j_25_umol_m2_s, 140)) < 0.002 &&
      abs(a2$rl_25_umol_m2_s - 1.2) < 0.01 && is.na(a2$tp_25_umol_m2_s),
      "noise-free J-limited curve: Vcmax, J, RL within 0.2%; TPU NA because it never limits")
vnorm <- calculate_temperature_response(
  exdf(data.frame(T = 30), units = data.frame(T = "degrees C")), c3_temperature_param_bernacchi, "T")[, "Vcmax_norm"]
check(isTRUE(all.equal(a1$vcmax_tleaf_umol_m2_s / a1$vcmax_25_umol_m2_s, vnorm)),
      "Vcmax at leaf temperature / Vcmax at 25 C = Bernacchi et al. (2001) factor at 30 C")
check(is.finite(a2$jmax_25_umol_m2_s) && a2$jmax_25_umol_m2_s > a2$j_25_umol_m2_s,
      "Jmax computed from J (Jmax > J at 1800 umol m-2 s-1)")
nf <- nf[nf$curve == "nf_tpu", ]

# 3. Input validation fails loudly -----------------------------------------------------
e1 <- tryCatch(validate_input(transform(nf, ci_umol_mol = ci_umol_mol / 10), "curve", DEFAULTS),
               error = function(e) conditionMessage(e))
check(grepl("Pa", e1), "Ci given in Pa is caught")
e2 <- tryCatch(validate_input(nf[, setdiff(names(nf), "oxygen_pct")], "curve", DEFAULTS),
               error = function(e) conditionMessage(e))
check(grepl("--oxygen-pct", e2), "missing oxygen column asks for --oxygen-pct")
e3 <- tryCatch(validate_input(nf, character(), DEFAULTS), error = function(e) conditionMessage(e))
check(grepl("--curve-cols", e3), "missing curve columns ask which columns identify a curve")
e4 <- tryCatch(read_tidy(system.file("extdata", "c3_aci_1.xlsx", package = "PhotoGEA")),
               error = function(e) conditionMessage(e))
check(grepl("claw-licor-6800", e4), "a raw LI-6800 file is redirected to claw-licor-6800")

# 4. End to end on the demo data -----------------------------------------------------------
out <- file.path(tmp, "demo")
r <- run_cli(c("--demo", "--cores", cores, "--out", shQuote(out)))
check(r$status == 0L, paste("demo run exits cleanly:", paste(r$output, collapse = "\n")))
for (f in c("report.md", "tables/aci_parameters.csv", "tables/aci_fitted_points.csv", "tables/point_use.csv",
            "figures/aci_fits_01.png", "figures/residuals.png", "reproducibility/commands.sh",
            "reproducibility/sessionInfo.txt")) check(file.exists(file.path(out, f)), paste("demo writes", f))
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
dp <- utils::read.csv(file.path(out, "tables", "aci_parameters.csv"), stringsAsFactors = FALSE)
check(identical(dp$curve_id, c("G1 - 101", "G1 - 102", "G2 - 201", "G2 - 202")), "demo: four curves")
check(identical(grepl("induction_drop", dp$qc_flags), c(FALSE, FALSE, FALSE, TRUE)),
      "demo: only curve G2 - 202 (simulated loss of activation) is flagged for induction")
check(all(abs(rel(dp$vcmax_25_umol_m2_s[1:3], c(95, 105, 75))) <= 0.10),
      "demo: Vcmax of the three steady curves within 10% of the generating values")
check(all(dp$genotype == c("G1", "G1", "G2", "G2")) && all(nzchar(dp$instrument_id)),
      "curve metadata (genotype, instrument) carried into the parameter table")
rep_txt <- readLines(file.path(out, "report.md"))
check(any(grepl("apparent", rep_txt)) && any(grepl("masquerade as down-regulation", rep_txt)),
      "report states parameters are apparent and explains the induction flag")
bad_run <- run_cli(c("--demo", "--pathway", "C5"))
check(bad_run$status != 0L, "invalid pathway exits with an error")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
