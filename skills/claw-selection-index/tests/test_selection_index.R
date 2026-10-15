# Tests for claw-selection-index. Run from anywhere: Rscript tests/test_selection_index.R
# Known answers: published selection-intensity tables, analytic order-statistic
# values, a hand-computed two-trait Smith-Hazel index, the worked examples in
# Ceron-Rojas & Crossa (2018, sections 2.8.2 and 3.1.4), and the bending example
# in Nilforooshan (2020). Design properties: Smith-Hazel maximises r_IH, desired
# gains are achieved in the requested proportions, restricted traits do not change.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "selection_index.R"))  # loads functions only

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
near <- function(x, y, tol) isTRUE(all(abs(x - y) <= tol))

run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "selection_index.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}

# 1. Selection intensity, infinite population: Cobb et al. (2019) Table 1
#    (p = 0.05, 0.01, 0.005, 0.001) and the standard 10% value 1.755;
#    p = 0.5 has the analytic value 2 * dnorm(0).
tab <- c(`0.1` = 1.755, `0.05` = 2.063, `0.01` = 2.665, `0.005` = 2.892, `0.001` = 3.367)
got <- intensity_infinite(as.numeric(names(tab)))
check(near(got, tab, 5e-4), "selection intensities match the published table to 3 decimals")
check(near(intensity_infinite(0.5), 2 * dnorm(0), 1e-12), "i(p = 0.5) = 2 * phi(0)")

# 2. Small populations: exact expected mean of the best k of N.
#    E[max of 2 N(0,1)] = 1/sqrt(pi); E[max of 3] = 3/(2 sqrt(pi)); the middle of
#    3 has expectation 0, so the mean of the best 2 of 3 is 3/(4 sqrt(pi)).
check(near(intensity_finite(1, 2), 1 / sqrt(pi), 1e-8), "best 1 of 2 = 1/sqrt(pi)")
check(near(intensity_finite(1, 3), 3 / (2 * sqrt(pi)), 1e-8), "best 1 of 3 = 3/(2 sqrt(pi))")
check(near(intensity_finite(2, 3), 3 / (4 * sqrt(pi)), 1e-8), "best 2 of 3 = 3/(4 sqrt(pi))")
i_fin <- intensity_finite(100, 1000)
check(i_fin < intensity_infinite(0.1) && intensity_infinite(0.1) - i_fin < 0.005,
      "finite-population intensity is just below the infinite value for N = 1000")
check(intensity_finite(3, 30) < intensity_finite(30, 300),
      "intensity increases with population size at a fixed proportion")

# 3. Breeder's equation by hand: i = 1.754983 (p = 0.1), h2 = 0.45, sigma_P = 300.
i10 <- intensity_infinite(0.1)
be <- breeders_equation(i10, h2 = 0.45, sigma_p = 300, cycle_years = 4)
check(near(be$response_per_cycle_phenotypic, i10 * 0.45 * 300, 1e-10), "R = i h2 sigma_P")
check(near(be$gain_per_year, i10 * 0.45 * 300 / 4, 1e-10), "gain per year = R / L")
be2 <- breeders_equation(i10, accuracy = sqrt(0.45), sigma_a = sqrt(0.45) * 300)
check(near(be2$response_per_cycle_accuracy, be$response_per_cycle_phenotypic, 1e-9),
      "i r sigma_A equals i h2 sigma_P when r = h and sigma_A = h sigma_P")
err <- tryCatch(breeders_equation(i10, h2 = 1.3, sigma_p = 1), error = conditionMessage)
check(grepl("h2", err), "h2 > 1 is rejected")

# 3b. Entry means from variance components (Holland et al. 2002):
#     var_P = var_G + var_GE/e + var_e/(r e). Hand example: 40 + 20/4 + 80/8 = 55.
hm <- entry_mean_h2(var_g = 40, var_e = 80, n_rep = 2, var_ge = 20, n_env = 4)
check(near(hm[["h2"]], 40 / 55, 1e-12) && near(hm[["sigma_p"]], sqrt(55), 1e-12),
      "entry-mean h2 and sigma_P from variance components")
scv <- run_scenarios(data.frame(scenario = "vc", proportion_selected = 0.1, var_g = 40,
                                var_ge = 20, var_e = 80, n_env = 4, n_rep = 2))
check(near(scv$response_per_cycle_phenotypic, i10 * (40 / 55) * sqrt(55), 1e-10),
      "scenario built from variance components")
err <- tryCatch(run_scenarios(data.frame(scenario = "bad", proportion_selected = 0.1, var_g = 40,
                                         var_ge = 20, var_e = 80, n_rep = 2)), error = conditionMessage)
check(grepl("var_ge and n_env together", err), "var_ge without n_env is rejected")
Gm <- matrix(c(2, 0.5, 0.5, 3), 2, dimnames = list(c("a", "b"), c("a", "b")))
GEm <- diag(2); dimnames(GEm) <- dimnames(Gm)
Em <- matrix(c(4, 1, 1, 8), 2, dimnames = dimnames(Gm))
check(near(entry_mean_P(Gm, Em, GEm, n_env = 2, n_rep = 2), matrix(c(3.5, 0.75, 0.75, 5.5), 2), 1e-12),
      "entry-mean P = G + GE/e + E/(re) (hand computed)")
wm <- function(M, f) utils::write.csv(data.frame(trait = rownames(M), M, check.names = FALSE), f, row.names = FALSE)
tdir <- tempfile("si")
dir.create(tdir)
wm(Gm, file.path(tdir, "G.csv")); wm(GEm, file.path(tdir, "GE.csv")); wm(Em, file.path(tdir, "E.csv"))
invisible(utils::capture.output(main(c("--g-matrix", file.path(tdir, "G.csv"), "--e-matrix", file.path(tdir, "E.csv"),
                                       "--ge-matrix", file.path(tdir, "GE.csv"), "--n-env", "2", "--n-rep", "2",
                                       "--weights", "a=1,b=1", "--out", file.path(tdir, "out")))))
pe <- utils::read.csv(file.path(tdir, "out", "tables", "P_entry_means.csv"))
check(near(as.matrix(pe[, -1]), matrix(c(3.5, 0.75, 0.75, 5.5), 2), 1e-12),
      "command line builds the entry-mean P from G, GE and E")
err <- tryCatch(main(c("--g-matrix", file.path(tdir, "G.csv"), "--e-matrix", file.path(tdir, "E.csv"),
                       "--n-env", "3", "--weights", "a=1,b=1", "--out", file.path(tdir, "o2"))),
                error = conditionMessage)
check(grepl("ge-matrix", err), "several environments without a GE matrix are rejected")

# 4. Hand-computed two-trait Smith-Hazel index.
#    P = [4 1; 1 9], G = [2 0.5; 0.5 3], w = (1, 1):
#    Gw = (2.5, 3.5); P^-1 = [9 -1; -1 4] / 35; b = (19/35, 23/70);
#    b'Pb = b'Gw = 175.5/70; sigma_H^2 = 6; Gb = (1.25, 88/70).
P2 <- matrix(c(4, 1, 1, 9), 2, dimnames = list(c("a", "b"), c("a", "b")))
G2 <- matrix(c(2, 0.5, 0.5, 3), 2, dimnames = dimnames(P2))
w2 <- c(a = 1, b = 1)
b2 <- smith_hazel(P2, G2, w2)
check(near(b2, c(19 / 35, 23 / 70), 1e-12), "Smith-Hazel b = P^-1 G w (hand computed)")
m2 <- index_metrics(b2, P2, G2, i10, w2)
check(near(m2$sigma_I, sqrt(175.5 / 70), 1e-12), "sigma_I (hand computed)")
check(near(m2$r_IH, sqrt((175.5 / 70) / 6), 1e-12), "r_IH = sigma_I / sigma_H for the optimum index")
check(near(m2$response_H, i10 * sqrt(175.5 / 70), 1e-12), "response in H = i sigma_I")
check(near(m2$gains, i10 * c(1.25, 88 / 70) / sqrt(175.5 / 70), 1e-12),
      "expected gain per trait = i G b / sigma_I (hand computed)")
check(near(sum(w2 * m2$gains), m2$response_H, 1e-12), "sum of w * gains = response in H")

# 5. Published worked example: Ceron-Rojas & Crossa (2018), section 2.8.2.
#    Inputs are printed to 2 decimals, so results agree to about 3 significant digits.
P3 <- matrix(c(62.50, -12.74, 8.53, -12.74, 17.52, -3.38, 8.53, -3.38, 12.31), 3, 3)
G3 <- matrix(c(36.21, -12.93, 8.35, -12.93, 13.04, -3.40, 8.35, -3.40, 9.96), 3, 3)
dimnames(P3) <- dimnames(G3) <- list(paste0("T", 1:3), paste0("T", 1:3))
w3 <- c(T1 = 1, T2 = -1, T3 = 1)
b3 <- smith_hazel(P3, G3, w3)
check(near(b3, c(0.555, -1.063, 1.087), 1e-3), "published LPSI coefficients (0.555, -1.063, 1.087)")
m3 <- index_metrics(b3, P3, G3, 2.063, w3)
check(near(m3$sigma_I^2, 86.72, 0.01), "published index variance 86.72")
check(near(m3$sigma_I, 9.312, 1e-3), "published sigma_I 9.312")
check(near(m3$response_H, 19.21, 5e-3), "published selection response 19.21 (k = 2.063)")
check(near(m3$gains, c(9.51, -5.48, 4.22), 0.01), "published expected gains (9.51, -5.48, 4.22)")
check(near(m3$r_IH, 0.894, 1e-3), "published r_IH 0.894")
check(near(index_metrics(w3, P3, G3, 2.063, w3)$r_IH, 0.875, 1e-3),
      "published base-index r_IH 0.875")

# 6. Published restricted index: Ceron-Rojas & Crossa (2018), section 3.1.4,
#    T1 restricted, k = 1.755 (10%).
P4 <- matrix(c(62.50, -12.74, 8.53, 2.73, -12.74, 17.52, -3.38, -2.28,
               8.53, -3.38, 12.31, 0.16, 2.73, -2.28, 0.16, 7.27), 4, 4, byrow = TRUE)
G4 <- matrix(c(36.21, -12.93, 8.35, 2.74, -12.93, 13.04, -3.40, -2.24,
               8.35, -3.40, 9.96, 0.16, 2.74, -2.24, 0.16, 6.64), 4, 4, byrow = TRUE)
dimnames(P4) <- dimnames(G4) <- list(paste0("T", 1:4), paste0("T", 1:4))
w4 <- c(T1 = 1, T2 = -1, T3 = 1, T4 = 1)
check(near(smith_hazel(P4, G4, w4), c(0.55, -1.05, 1.09, 1.06), 0.005),
      "published 4-trait LPSI coefficients")
br <- restricted_index(P4, G4, w4, "T1")
check(near(br, c(-0.35, -0.41, 0.59, 0.89), 0.005), "published RLPSI coefficients")
mr <- index_metrics(br, P4, G4, 1.755, w4)
check(near(mr$response_H, 6.87, 0.01), "published RLPSI response 6.87")
check(mr$gains[["T1"]] == 0, "restricted trait has zero expected gain")
check(near(mr$gains[-1], c(-2.2, 2.03, 2.66), c(0.05, 0.01, 0.01)),
      "published RLPSI gains (-2.2, 2.03, 2.66)")
check(near(mr$r_IH, 0.35, 0.005), "published RLPSI r_IH 0.35")

# 7. Smith-Hazel maximises r_IH: no random weighting beats it (demo matrices).
dd <- file.path(SKILL_TEST_DIR, "demo")
Pd <- read_matrix(file.path(dd, "P_matrix.csv"), "P")
Gd <- read_matrix(file.path(dd, "G_matrix.csv"), "G")
wd <- read_trait_values(file.path(dd, "weights.csv"), "weight", rownames(Pd), "Weights")
r_sh <- index_metrics(smith_hazel(Pd, Gd, wd), Pd, Gd, 1, wd)$r_IH
set.seed(11)
scale_b <- 1 / sqrt(diag(Pd))
r_rand <- vapply(seq_len(5000), function(k)
  index_metrics(stats::rnorm(4) * scale_b, Pd, Gd, 1, wd)$r_IH, 0)
check(max(r_rand) <= r_sh + 1e-12, "no random index has a higher r_IH than Smith-Hazel")
check(r_sh - max(r_rand) < 0.05, "best random index approaches the optimum (search is not degenerate)")

# 8. Desired gains are achieved in the requested proportions (Pesek & Baker 1969).
dg <- read_trait_values(file.path(dd, "desired_gains.csv"), "desired_gain", rownames(Pd), "D")
md <- index_metrics(desired_gains_index(Pd, Gd, dg), Pd, Gd, 2.063)
nz <- dg != 0
ratio <- md$gains[nz] / dg[nz]
check(max(abs(ratio / ratio[1] - 1)) < 1e-10, "expected gains proportional to desired gains")
check(all(md$gains[!nz] == 0), "traits with desired gain 0 do not change")
d_mix <- c(yield_kg_ha = -30, protein_pct = 0.5, oil_pct = -0.1, height_cm = 2)
mm <- index_metrics(desired_gains_index(Pd, Gd, d_mix), Pd, Gd, 2.063)
check(max(abs(mm$gains / d_mix / (mm$gains[1] / d_mix[1]) - 1)) < 1e-10 && mm$gains[1] / d_mix[1] > 0,
      "proportional gains with mixed-sign desired gains")
check(near(desired_gains_index(Pd, Gd, d_mix), solve(Gd, d_mix), 1e-8),
      "general form equals the Pesek-Baker b = G^-1 d when G is invertible")

# 9. Restricted index on the demo: zero change in height, lower r_IH than Smith-Hazel.
mrd <- index_metrics(restricted_index(Pd, Gd, wd, "height_cm"), Pd, Gd, 1, wd)
check(mrd$gains[["height_cm"]] == 0 && mrd$r_IH < r_sh, "restriction costs accuracy and holds height")

# 10. Bending: Nilforooshan (2020) Table 1, unweighted method of Jorjani et al.
#     (2003) on the 5 x 5 non-PD matrix V.
V <- matrix(c(100, 95, 80, 40, 40, 95, 100, 95, 80, 40, 80, 95, 100, 95, 80,
              40, 80, 95, 100, 95, 40, 40, 80, 95, 100), 5, 5, byrow = TRUE)
check(near(round(eigen(V)$values, 2), c(399.48, 98.52, 23.65, -3.12, -18.52), 1e-9),
      "published eigenvalues of the non-PD example matrix")
bs <- bend_stats(V, bend_matrix(V, 1e-4))
check(near(unlist(bs), c(-5.9296, 6.4973, 0.7235, 3.3727, 0.9856, 3.9275), 1e-4),
      "published bending statistics (eps = 1e-4)")
bs2 <- bend_stats(V, bend_matrix(V, 1e-2))
check(near(unlist(bs2), c(-5.9320, 6.5016, 0.7241, 3.3754, 0.9856, 3.9300), 1e-4),
      "published bending statistics (eps = 1e-2)")
check(min(eigen(bend_matrix(V, 1e-4))$values) >= 1e-4 - 1e-9, "bent matrix is positive definite")

# 11. Bad inputs fail loudly.
Gbad <- Gd
Gbad["yield_kg_ha", "protein_pct"] <- Gbad["protein_pct", "yield_kg_ha"] <- 300
err <- tryCatch(check_matrices(Pd, Gbad), error = conditionMessage)
check(grepl("positive semi-definite", err) && grepl("--bend", err),
      "non-PSD G stops with an explanation and the --bend option")
ok <- check_matrices(Pd, Gbad, bend = TRUE)
check(min(ok$ev_g) > 0 && !is.null(ok$bent), "--bend repairs G and records the change")
err <- tryCatch(read_trait_values("yield_kg_ha=1", "weight", rownames(Pd), "Weights"),
                error = conditionMessage)
check(grepl("no value for", err), "missing economic weights are not filled silently")

# 12. End to end on the demo data.
out <- file.path(tempdir(), "claw-selection-index-test")
unlink(out, recursive = TRUE)
res <- run_cli(c("--demo", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "tables/index_gains.csv", "tables/index_summary.csv",
            "tables/breeders_equation.csv", "tables/input_checks.csv",
            "figures/expected_gains.png", "figures/breeders_equation.png",
            "reproducibility/commands.sh", "reproducibility/sessionInfo.txt")) {
  check(file.exists(file.path(out, f)), paste("demo writes", f))
}
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L,
      "demo writes a checksum file")
s <- utils::read.csv(file.path(out, "tables", "index_summary.csv"))
check(nrow(s) == 4L, "demo computes four indices")
g <- utils::read.csv(file.path(out, "tables", "index_gains.csv"))
check(g$expected_gain[g$index == "restricted_height_cm" & g$trait == "height_cm"] == 0,
      "demo restricted index holds height")
sc <- utils::read.csv(file.path(out, "tables", "breeders_equation.csv"))
check(near(sc$intensity[1], intensity_infinite(0.1), 1e-9), "demo scenario intensity")

# 13. Unknown options are rejected.
bad <- run_cli(c("--demo", "--colour", "red"))
check(bad$status != 0L, "unknown option exits with an error")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
