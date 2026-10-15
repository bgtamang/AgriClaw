# Generate the synthetic demo data for claw-stomatal-models (fixed seed).
#   Rscript demo/make_demo.R
# Writes demo/stomatal_demo.csv in the tidy layout of claw-licor-6800.
#
# Each point is a self-consistent leaf gas-exchange state: given A, chamber
# CO2 (Ca) and water vapour (H2O_s), leaf temperature, boundary-layer
# conductance and pressure, gsw is solved by fixed-point iteration so that the
# stomatal model holds at the leaf surface, where hs, Cs and D are computed
# exactly as PhotoGEA::calculate_gas_properties will compute them back from
# the logged columns (LI-6800 manual Appendix C). The demo uses the Medlyn
# model (g0 = 0.02 mol m-2 s-1; g1 = 3.5 kPa^0.5 for G1 and 4.5 for G2) with
# Gaussian scatter (sd 0.01 mol m-2 s-1) on gsw, along a light series.

svp_kpa <- function(t) 0.6135 * exp(17.502 * t / (240.97 + t))  # LI-6800 eq. C-70

# One gas-exchange state. model: "medlyn" or "ball-berry"; eps: scatter on gsw.
simulate_point <- function(a, ca, tleaf, h2o_s, gbw, p_kpa, model, g0, g1, eps = 0) {
  svp <- svp_kpa(tleaf)
  wi <- 1000 * svp / p_kpa
  if (h2o_s >= wi) stop("chamber air is saturated relative to the leaf")
  gsw <- 0.3
  for (k in seq_len(500)) {
    f <- function(w) gbw * (w - h2o_s) / (1000 - (w + h2o_s) / 2) - gsw * (wi - w) / (1000 - (wi + w) / 2)
    w <- stats::uniroot(f, c(h2o_s, wi), tol = 1e-13)$root
    e <- gbw * (w - h2o_s) / (1000 - (w + h2o_s) / 2)
    gbc <- gbw / 1.37
    cs <- ((gbc - e / 2) * ca - a) / (gbc + e / 2)
    hs <- 0.1 * w * p_kpa / svp / 100
    d <- svp - w * p_kpa / 1000
    g_new <- if (model == "ball-berry") g0 + g1 * a * hs / cs else g0 + 1.6 * (1 + g1 / sqrt(d)) * a / cs
    g_new <- g_new + eps
    if (abs(g_new - gsw) < 1e-14) break
    gsw <- g_new
  }
  c(gsw = gsw, e = e, h2o_surf = w, cs = cs, hs = hs, d = d)
}

# A curve of points with given A and chamber humidity values.
simulate_curve <- function(a, h2o_s, model, g0, g1, tleaf = 27, ca = 400, gbw = 2.9,
                           pa = 99.5, dp = 0.1, noise = 0.01) {
  p <- pa + dp
  rows <- lapply(seq_along(a), function(i) {
    s <- simulate_point(a[i], ca, tleaf, h2o_s[i], gbw, p, model, g0, g1, stats::rnorm(1, 0, noise))
    data.frame(a_umol_m2_s = a[i], gsw_mol_m2_s = s[["gsw"]], e_mol_m2_s = s[["e"]],
               gbw_mol_m2_s = gbw, h2o_s_mmol_mol = h2o_s[i], ca_umol_mol = ca, tleaf_c = tleaf,
               pa_kpa = pa, delta_pcham_kpa = dp,
               vpd_leaf_kpa = svp_kpa(tleaf) - h2o_s[i] * p / 1000,  # LI-6800 VPDleaf: leaf to chamber air
               true_hs = s[["hs"]], true_cs_umol_mol = s[["cs"]], true_d_kpa = s[["d"]])
  })
  do.call(rbind, rows)
}

TRUTH <- data.frame(genotype = c("G1", "G1", "G2", "G2"), plot = c("101", "102", "201", "202"),
                    g0 = 0.02, g1 = c(3.5, 3.5, 4.5, 4.5), amax = c(30, 33, 26, 28),
                    stringsAsFactors = FALSE)
Q_SERIES <- c(2000, 1500, 1100, 800, 500, 350, 200)

make_demo <- function() {
  parts <- lapply(seq_len(nrow(TRUTH)), function(i) {
    t <- TRUTH[i, ]
    a <- t$amax * (1 - exp(-0.0022 * Q_SERIES)) + stats::rnorm(length(Q_SERIES), 0, 0.2)
    h2o <- 0.6 * 1000 * svp_kpa(27) / 99.6 + stats::rnorm(length(Q_SERIES), 0, 0.3)
    d <- simulate_curve(a, h2o, "medlyn", t$g0, t$g1)
    data.frame(source_file = "demo_synthetic", instrument_id = "68C-000101/68H-000101",
               genotype = t$genotype, plot = t$plot, obs = seq_along(a), qin_umol_m2_s = Q_SERIES,
               d[, setdiff(names(d), c("true_hs", "true_cs_umol_mol", "true_d_kpa"))],
               stable_n = 2, stable_total = 2, qc_hard_fail = FALSE, qc_unstable = FALSE,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, parts)
  num <- vapply(out, is.numeric, logical(1))
  out[num] <- lapply(out[num], function(v) signif(v, 8))
  out
}

main <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else "."
  set.seed(20261015)
  utils::write.csv(make_demo(), file.path(here, "stomatal_demo.csv"), row.names = FALSE)
  cat("Wrote stomatal_demo.csv to", here, "\n")
}

if (sys.nframe() == 0L) main()
