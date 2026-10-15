# Tests for claw-light-response. Run from anywhere: Rscript tests/test_light_response.R
# Known answers: the light compensation and saturation points against
# numerical roots of the fitted model; exact recovery from a noise-free curve;
# recovery of known parameters from 20 noisy simulated curves; quantum yield
# on an absorbed basis; QC flags on curves built to trigger them; input
# validation; an end-to-end --demo run.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "light_response.R"))  # loads functions; main() does not run

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "light_response.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
tmp <- file.path(tempdir(), "claw-light-response-test")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE)
Q <- c(2000, 1500, 1200, 1000, 800, 600, 400, 300, 200, 150, 100, 75, 50, 25, 0)
rel <- function(est, tru) (est - tru) / tru
mk <- function(id, phi, theta, pmax, rd, noise = 0, q = Q) {
  data.frame(curve = id, qin_umol_m2_s = q,
             a_umol_m2_s = nrh_net(q, phi, theta, pmax, rd) + stats::rnorm(length(q), 0, noise))
}
opt <- DEFAULTS
opt$curve_cols <- "curve"

# 1. Light compensation and saturation points equal the numerical roots ----------
set.seed(1)
for (k in 1:25) {
  phi <- stats::runif(1, 0.03, 0.09); th <- stats::runif(1, 0.05, 0.99)
  pm <- stats::runif(1, 10, 50); rd <- stats::runif(1, 0.3, 3); f <- stats::runif(1, 0.5, 0.95)
  lcp <- lcp_of(phi, th, pm, rd)
  root <- stats::uniroot(function(q) nrh_net(q, phi, th, pm, rd), c(0, 1000), tol = 1e-12)$root
  qs <- qsat_of(phi, th, pm, rd, f)
  if (!isTRUE(abs(lcp - root) <= 1e-6 * root) ||
      !isTRUE(abs(nrh_net(qs, phi, th, pm, rd) - f * (pm - rd)) <= 1e-9 * pm)) {
    stop("FAIL: analytical LCP/saturation point, parameter set ", k, call. = FALSE)
  }
}
n_checks <- n_checks + 1L  # 25 parameter sets, one check
check(isTRUE(all.equal(nrh_gross(1e6, 0.05, 0.7, 30), 30, tolerance = 1e-3)), "gross rate tends to Pmax at high light")
check(isTRUE(all.equal(nrh_gross(1e-3, 0.05, 0.7, 30) / 1e-3, 0.05, tolerance = 1e-4)), "initial slope equals phi")
check(isTRUE(all.equal(nrh_gross(500, 0.05, 1, 30), min(0.05 * 500, 30))), "theta = 1 gives the Blackman response")
check(isTRUE(all.equal(nrh_gross(500, 0.05, 0, 30), 0.05 * 500 * 30 / (0.05 * 500 + 30))),
      "theta = 0 gives the rectangular hyperbola")

# 2. Exact recovery from a noise-free curve ----------------------------------------
r0 <- analyze_light(mk("nf", 0.055, 0.75, 30, 1.4), opt)$params
check(abs(rel(r0$phi_incident_mol_mol, 0.055)) < 1e-4 && abs(r0$theta - 0.75) < 1e-4 &&
      abs(rel(r0$asat_umol_m2_s, 28.6)) < 1e-4 && abs(r0$rd_umol_m2_s - 1.4) < 1e-4,
      "noise-free curve: phi, theta, Asat, Rd recovered to 1e-4")
check(abs(rel(r0$lcp_umol_m2_s, lcp_of(0.055, 0.75, 30, 1.4))) < 1e-4, "noise-free LCP equals the analytical value")
check(nchar(r0$qc_flags) == 0L, "a clean curve raises no flags")

# 3. Recovery from 20 noisy curves (A noise sd 0.3 umol m-2 s-1) ---------------------
set.seed(7)
n <- 20L
tr <- data.frame(curve = sprintf("n%02d", seq_len(n)), phi = stats::runif(n, 0.04, 0.075),
                 theta = stats::runif(n, 0.55, 0.9), pmax = stats::runif(n, 20, 45), rd = stats::runif(n, 0.8, 2.5))
sim <- do.call(rbind, lapply(seq_len(n), function(i) with(tr[i, ], mk(curve, phi, theta, pmax, rd, noise = 0.3))))
rs <- analyze_light(sim, opt)$params
rs <- rs[match(tr$curve, rs$curve_id), ]
e_phi <- rel(rs$phi_incident_mol_mol, tr$phi)
e_as <- rel(rs$asat_umol_m2_s, tr$pmax - tr$rd)
e_th <- rs$theta - tr$theta
e_rd <- rs$rd_umol_m2_s - tr$rd
e_lcp <- rel(rs$lcp_umol_m2_s, lcp_of(tr$phi, tr$theta, tr$pmax, tr$rd))
cat(sprintf("max |error|: phi %.1f%%, Asat %.1f%%, theta %.3f, Rd %.2f, LCP %.1f%%\n",
            100 * max(abs(e_phi)), 100 * max(abs(e_as)), max(abs(e_th)), max(abs(e_rd)), 100 * max(abs(e_lcp))))
# Errors are judged against the fit's own standard errors (z = error / SE):
# an unbiased fit with honest SEs keeps every |z| small and about 95% within 2.
z <- cbind(phi = e_phi * tr$phi / rs$phi_se, theta = e_th / rs$theta_se, rd = e_rd / rs$rd_se,
           asat = e_as * (tr$pmax - tr$rd) / rs$asat_se)
cat(sprintf("max |z|: %s\n", paste(sprintf("%s %.2f", colnames(z), apply(abs(z), 2, max)), collapse = ", ")))
check(all(abs(z) <= 3.5), "every estimate of phi, theta, Rd and Asat within 3.5 standard errors of the truth")
check(all(colMeans(abs(z) <= 2) >= 0.75), "at least 75% of estimates within 2 standard errors, for each parameter")
check(max(abs(e_phi)) <= 0.20 && stats::median(abs(e_phi)) <= 0.05, "phi within 20% (median within 5%)")
check(max(abs(e_as)) <= 0.05 && stats::median(abs(e_as)) <= 0.02, "Asat within 5% (median within 2%)")
check(max(abs(e_th)) <= 0.2, "theta within 0.2")
check(max(abs(e_rd)) <= 0.5, "Rd within 0.5 umol m-2 s-1")
check(max(abs(e_lcp)) <= 0.30, "light compensation point within 30%")
cov <- mean(tr$phi >= rs$phi_lower_mol_mol & tr$phi <= rs$phi_upper_mol_mol)
check(cov >= 0.75, sprintf("true phi inside the 95%% Wald interval for >= 75%% of curves (%.0f%%)", 100 * cov))

# 4. Absorbed-light basis -------------------------------------------------------------
oa <- opt
oa$absorptance <- 0.85
ra <- analyze_light(mk("nf", 0.055, 0.75, 30, 1.4), oa)$params
check(isTRUE(all.equal(ra$phi_absorbed_mol_mol, ra$phi_incident_mol_mol / 0.85)), "phi(absorbed) = phi(incident) / absorptance")

# 5. QC flags --------------------------------------------------------------------------
pi_curve <- mk("pi", 0.055, 0.75, 30, 1.4)
pi_curve$a_umol_m2_s <- pi_curve$a_umol_m2_s - 0.006 * pmax(pi_curve$qin_umol_m2_s - 1000, 0)
check(grepl("photoinhibition", analyze_light(pi_curve, opt)$params$qc_flags), "decline at high light flagged as photoinhibition")
sparse <- mk("sp", 0.055, 0.75, 30, 1.4, q = c(2000, 1500, 1200, 1000, 800, 600, 400, 300, 200, 100, 0))
check(grepl("few_low_light_points", analyze_light(sparse, opt)$params$qc_flags),
      "only one point between 0 and 150 flags the quantum yield")
nodark <- mk("nd", 0.055, 0.75, 30, 1.4, q = c(2000, 1500, 1200, 1000, 800, 600, 400, 300, 200, 150, 120, 100, 75))
check(grepl("no_dark_point", analyze_light(nodark, opt)$params$qc_flags), "no point below 50 flags extrapolated Rd")
of <- opt
of$theta <- 0.7
rf <- analyze_light(mk("fx", 0.055, 0.75, 30, 1.4), of)$params
check(rf$theta == 0.7 && !rf$theta_fitted, "--theta fixes the curvature")

# 6. Input validation ------------------------------------------------------------------------
bad <- mk("u", 0.055, 0.75, 30, 1.4)
bad$qin_umol_m2_s <- bad$qin_umol_m2_s / 1e6
e1 <- tryCatch(validate_input(bad, "curve", DEFAULTS), error = function(e) conditionMessage(e))
check(grepl("umol m-2 s-1", e1), "PPFD in mol m-2 s-1 is caught")
e2 <- tryCatch(validate_input(mk("u", 0.055, 0.75, 30, 1.4), character(), DEFAULTS), error = function(e) conditionMessage(e))
check(grepl("--curve-cols", e2), "missing curve columns ask which columns identify a curve")

# 7. End to end on the demo data -----------------------------------------------------------
out <- file.path(tmp, "demo")
r <- run_cli(c("--demo", "--out", shQuote(out)))
check(r$status == 0L, paste("demo run exits cleanly:", paste(r$output, collapse = "\n")))
for (f in c("report.md", "tables/light_parameters.csv", "tables/light_fitted_points.csv", "figures/light_fits_01.png",
            "reproducibility/commands.sh", "reproducibility/sessionInfo.txt"))
  check(file.exists(file.path(out, f)), paste("demo writes", f))
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
dp <- utils::read.csv(file.path(out, "tables", "light_parameters.csv"), stringsAsFactors = FALSE)
check(nrow(dp) == 4L && identical(grepl("photoinhibition", dp$qc_flags), c(FALSE, FALSE, FALSE, TRUE)),
      "demo: only the curve generated with photoinhibition is flagged")
check(all(abs(rel(dp$asat_umol_m2_s[1:3], c(30.5, 33.3, 22.8))) <= 0.05),
      "demo: Asat of the three regular curves within 5% of the generating values")
rep_txt <- readLines(file.path(out, "report.md"))
check(any(grepl("incident-light basis", rep_txt)), "report says the quantum yield is on an incident basis")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
