# Tests for claw-fluorescence. Run from anywhere: Rscript tests/test_fluorescence.R
# Known-answer checks: a hand-worked example, algebraic identities between the
# parameters, agreement with the LI-6800's own calculations in the example files
# shipped with PhotoGEA, QC flags, bad input, and end-to-end command-line runs.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "fluorescence.R"))  # loads functions; main() does not run

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
near <- function(a, b, tol) isTRUE(all(abs(a - b) <= tol))

run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "fluorescence.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
errmsg <- function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))

# 1. Worked example, by hand: Fo 300, Fm 1500, Fs 600, Fm' 1000, Fo' 250, PPFD 1000
p <- fluor_params(fo = 300, fm = 1500, fs = 600, fm_prime = 1000, fo_prime = 250, ppfd = 1000)
check(near(p$fv_fm, 0.8, 1e-12), "Fv/Fm = (1500 - 300) / 1500 = 0.8")
check(near(p$phi_psii, 0.4, 1e-12), "PhiPSII = (1000 - 600) / 1000 = 0.4")
check(near(p$fvp_fmp, 0.75, 1e-12), "Fv'/Fm' = (1000 - 250) / 1000 = 0.75")
check(near(p$qp, 400 / 750, 1e-12), "qP = (1000 - 600) / (1000 - 250) = 0.5333")
check(near(p$ql, (400 / 750) * 250 / 600, 1e-12), "qL = qP x Fo'/Fs = 0.2222")
check(near(p$npq, 0.5, 1e-12), "NPQ = (1500 - 1000) / 1000 = 0.5")
check(near(p$etr_umol_m2_s, 168, 1e-9), "ETR = 0.4 x 1000 x 0.84 x 0.5 = 168")
check(p$fo_prime_source == "measured", "measured Fo' is used when present")
e <- fluor_params(fo = 300, fm = 1500, fs = 600, fm_prime = 1000, fo_prime = NA, ppfd = 1000)
check(near(e$fo_prime_used, 300 / 1.1, 1e-9) && e$fo_prime_source == "estimated",
      "Fo' estimated as Fo / (Fv/Fm + Fo/Fm') = 300 / 1.1 = 272.73")
check(near(e$qp, 0.55, 1e-12) && near(e$ql, 0.25, 1e-12) && near(e$fvp_fmp, 8 / 11, 1e-12),
      "with estimated Fo': qP = 0.55, qL = 0.25, Fv'/Fm' = 0.7273")
check(near(fluor_params(300, 1500, 600, 1000, 250, 1000, absorptance = 0.9, f_psii = 0.4)$etr_umol_m2_s,
           0.4 * 1000 * 0.9 * 0.4, 1e-9), "ETR uses the given absorptance and PSII fraction")
check(near(fluor_params(300, 1500, 600, 1000, 250, ppfd = 1000, qabs = 840)$etr_umol_m2_s, 168, 1e-9),
      "ETR from absorbed PPFD (LI-6800 Qabs_fs) x PS2 fraction")

# 2. Identities on random valid data -------------------------------------------------
set.seed(5)
n <- 2000
fo <- runif(n, 200, 400); fm <- fo / runif(n, 0.15, 0.25)
fm_p <- fm / (1 + runif(n, 0, 3))
fo_p <- fo / ((fm - fo) / fm + fo / fm_p)
fs <- fm_p - runif(n, 0.05, 1) * (fm_p - fo_p)
r <- fluor_params(fo, fm, fs, fm_p, fo_p, ppfd = 1000)
check(near(r$phi_psii, r$qp * r$fvp_fmp, 1e-12), "PhiPSII = qP x Fv'/Fm' (Genty et al. 1989)")
check(all(r$ql <= r$qp + 1e-12) && all(fs >= fo_p), "qL <= qP whenever Fs >= Fo'")
check(near(r$ql, r$qp * r$fo_prime_used / fs, 1e-12), "qL = qP x Fo'/Fs")
est <- fluor_params(fo, fm, fs, fm_p, NA, 1000)$fo_prime_used
check(near(1 / est - 1 / fm_p, 1 / fo - 1 / fm, 1e-12),
      "estimated Fo' satisfies 1/Fo' - 1/Fm' = 1/Fo - 1/Fm (Oxborough & Baker 1997)")
same <- fluor_params(fo, fm, fs, fm, NA, 1000)
check(near(same$fo_prime_used, fo, 1e-9) && near(same$npq, 0, 1e-12), "with Fm' = Fm: estimated Fo' = Fo and NPQ = 0")

# 3. LI-6800 example files in PhotoGEA: match the instrument's own values --------------
files <- c("ball_berry_1.xlsx", "ball_berry_2.xlsx", "c3_aci_1.xlsx", "c3_aci_2.xlsx", "c4_aci_1.xlsx",
           "c4_aci_2.xlsx", "plaintext_licor_file", "plaintext_licor_file_v2")
n_rows <- 0L; worst <- 0
for (f in files) {
  path <- system.file("extdata", f, package = "PhotoGEA")
  check(nzchar(path) && detect_format(path) == "licor", paste("PhotoGEA example found and detected as LI-6800:", f))
  res <- analyze(standardize(read_licor(path), "licor", NA_real_, NA_real_), dark_adapted = "yes")
  for (v in names(LICOR_OUTPUTS)) {
    a <- res[[v]]; b <- res[[paste0("licor_", v)]]
    d <- max(abs(a - b) / abs(b))
    worst <- max(worst, d)
    check(is.finite(d) && d < 1e-9, sprintf("%s: %s matches the LI-6800 value (max rel. diff %.1e)", f, LICOR_OUTPUTS[[v]], d))
  }
  n_rows <- n_rows + nrow(res)
  check(all(grepl("Fv/Fm > 0.86", res$qc_check)), paste(f, "- Fv/Fm 0.925 in these files is flagged as implausible"))
}
check(n_rows == 424L, "all 424 example rows compared")
td <- licor_id_time("RECT-12243-20210724-05_20_30"); tl <- licor_id_time("MPF-13020-20210823-09_30_24")
check(near(as.numeric(difftime(tl, td, units = "days")), 30.17, 0.01), "LI-6800 flash IDs decoded to timestamps")

# 4. QC flags ------------------------------------------------------------------------------
q <- analyze(data.frame(fo = c(300, 300, 300, 300), fm = c(1500, 900, 1500, 1500),
                        fs = c(600, 600, 1000, 200), fm_prime = c(1000, 1000, 1000, 1000),
                        fo_prime = c(250, 250, 250, 250), ppfd = 1000, qabs = NA_real_,
                        absorptance = 0.84, f_psii = 0.5))
check(!nzchar(q$qc_check[1]) && !nzchar(q$qc_invalid[1]), "a valid row has no flags")
check(grepl("Fm' > Fm", q$qc_check[2]) && q$npq[2] < 0, "Fm' > Fm flagged (negative NPQ)")
check(grepl("Fv/Fm < 0.75", q$qc_check[2]), "low Fv/Fm flagged")
check(grepl("Fs >= Fm'", q$qc_invalid[3]), "Fs >= Fm' is invalid")
check(grepl("Fs < Fo'", q$qc_check[4]) && q$qp[4] > 1, "Fs < Fo' flagged (qP > 1)")
nd <- analyze(data.frame(fo = 300, fm = 1500, fs = 600, fm_prime = 1000, fo_prime = NA_real_, ppfd = 1000,
                         qabs = NA_real_, absorptance = 0.84, f_psii = 0.5), dark_adapted = "no")
check(is.na(nd$fv_fm) && is.na(nd$npq) && is.na(nd$qp) && near(nd$phi_psii, 0.4, 1e-12),
      "not dark-adapted: Fv/Fm, NPQ and estimated-Fo' parameters withheld; PhiPSII kept")

lo <- analyze(standardize(data.frame(fs = 600, fm_prime = 1000, ppfd_umol_m2_s = 1000), "csv", NA_real_, NA_real_))
check(near(lo$phi_psii, 0.4, 1e-12) && near(lo$etr_umol_m2_s, 168, 1e-9) && is.na(lo$fv_fm) && is.na(lo$qp),
      "light-adapted data only: PhiPSII and ETR computed, Fv/Fm and qP left empty")

# 5. Bad input -----------------------------------------------------------------------------
check(grepl("fm_prime", errmsg(standardize(data.frame(fo = 1, fm = 2, fs = 1), "csv", NA, NA))),
      "missing fm_prime is named in the error")
check(grepl("numeric", errmsg(standardize(data.frame(fs = "a", fm_prime = 2), "csv", NA, NA))),
      "non-numeric signal is rejected")

# 6. End to end --------------------------------------------------------------------------------
out <- file.path(tempdir(), "claw-fluorescence-test")
unlink(out, recursive = TRUE)
res <- run_cli(c("--demo", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "tables/fluorescence_parameters.csv", "tables/group_summary.csv",
            "figures/fv_fm.png", "figures/phipsii_etr_npq.png", "reproducibility/commands.sh",
            "reproducibility/sessionInfo.txt")) {
  check(file.exists(file.path(out, f)), paste("demo writes", f))
}
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
fp <- utils::read.csv(file.path(out, "tables", "fluorescence_parameters.csv"), check.names = FALSE, stringsAsFactors = FALSE)
check(all(grepl("Fv/Fm > 0.86", fp$qc_check[fp$plant_id == "G1-4"])), "demo: underestimated Fo on G1-4 is flagged")
check(any(grepl("Fm' > Fm", fp$qc_check[fp$plant_id == "G2-3"])), "demo: Fm' > Fm on G2-3 is flagged")
check(sum(nzchar(fp$qc_check[!fp$plant_id %in% c("G1-4", "G2-3")])) == 0, "demo: no other plant is flagged")

lic <- system.file("extdata", "c3_aci_1.xlsx", package = "PhotoGEA")
out2 <- file.path(tempdir(), "claw-fluorescence-licor")
unlink(out2, recursive = TRUE)
res2 <- run_cli(c("--input", shQuote(lic), "--dark-adapted", "yes", "--out", shQuote(out2)))
check(res2$status == 0L, paste("LI-6800 run exits cleanly:", paste(res2$output, collapse = "\n")))
rep2 <- readLines(file.path(out2, "report.md"), encoding = "UTF-8")
check(any(grepl("recorded more than 12 h before", rep2)), "LI-6800 run warns that the dark-adapted Fo/Fm is from another day")
check(any(grepl("Agreement with the LI-6800", rep2)), "LI-6800 run reports agreement with the instrument")
bad <- run_cli(c("--demo", "--color", "red"))
check(bad$status != 0L, "unknown option exits with an error")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
