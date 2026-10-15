# Generate the synthetic demo curves for claw-aci-fit (fixed seed).
#   Rscript demo/make_demo.R
# Writes demo/aci_demo_c3.csv and demo/aci_demo_c4.csv in the tidy layout of
# claw-licor-6800. True parameter values are listed in demo/README.md.
#
# Model: PhotoGEA's own FvCB (calculate_c3_assimilation, Bernacchi et al.
# temperature responses, Cc = Ci) and von Caemmerer C4 model
# (calculate_c4_assimilation, PhotoGEA defaults gmc = 1, gbs = 0.003).
# Ci is solved so that A = gsw / 1.6 * (Ca - Ci) at every point, and
# Ca = CO2 setpoint - 1.2 * A (6 cm2 leaf, 500 umol/s flow).

suppressPackageStartupMessages(library(PhotoGEA))

SETPOINTS <- c(400, 300, 200, 150, 100, 75, 50, 20, 400, 400, 600, 800, 1000, 1200, 1500, 1800)

# Relative stomatal conductance: opens at low CO2, closes at high CO2.
gsw_rel <- function(ca) ifelse(ca < 400, 1 + 0.15 * (400 - ca) / 380, 1 - 0.25 * (ca - 400) / 1400)

one_row_exdf <- function(cc, tleaf, pathway) {
  d <- exdf(data.frame(Cc = cc, PCm = cc, TleafCnd = tleaf, total_pressure = 0.996, Oxygen = 21))
  d <- document_variables(d, c("", "Cc", "micromol mol^(-1)"), c("", "PCm", "microbar"),
                          c("", "TleafCnd", "degrees C"), c("", "total_pressure", "bar"),
                          c("", "Oxygen", "percent"))
  tp <- if (pathway == "C3") c3_temperature_param_bernacchi else c4_temperature_param_vc
  calculate_temperature_response(d, tp, "TleafCnd")
}

# Net assimilation at intercellular CO2 `ci` (C3: Cc = Ci; C4: PCm = Ci*P - A/gmc).
model_a <- function(ci, tleaf, p, pathway) {
  if (pathway == "C3") {
    d <- one_row_exdf(ci, tleaf, "C3")
    return(calculate_c3_assimilation(d, 0, 0, 0, 0, "", p$J, "", "", p$RL, p$Tp, p$Vcmax)[, "An"])
  }
  d <- one_row_exdf(ci, tleaf, "C4")
  gm <- 1 * d[, "gmc_norm"]
  f <- function(a) {
    d[, "PCm"] <- ci * 0.996 - a / gm
    calculate_c4_assimilation(d, 0, 0.003, 1000, p$RL, 0.5, p$Vcmax, p$Vpmax, 1000)[, "An"] - a
  }
  stats::uniroot(f, c(-30, 120), tol = 1e-10)$root
}

solve_point <- function(sp, gs0, tleaf, p, pathway) {
  ca <- sp
  for (k in 1:3) {
    g <- gs0 * gsw_rel(ca)
    f <- function(ci) model_a(ci, tleaf, p, pathway) - g / 1.6 * (ca - ci)
    ci <- stats::uniroot(f, c(1, ca + 300), tol = 1e-8)$root
    a <- g / 1.6 * (ca - ci)
    ca <- sp - 1.2 * a
  }
  c(ca = ca, ci = ci, a = a, gsw = g)
}

# One curve; `after` scales Vcmax and J for points after the low-CO2 phase
# (loss of activation: the induction problem).
sim_curve <- function(p, tleaf, gs0, pathway, after = 1, noise_a = 0.3, noise_ci = 1.5) {
  rows <- lapply(seq_along(SETPOINTS), function(i) {
    pp <- p
    if (i >= 9) { pp$Vcmax <- p$Vcmax * after; if (!is.null(p$J)) pp$J <- p$J * after }
    tl <- tleaf + stats::rnorm(1, 0, 0.15)
    s <- solve_point(SETPOINTS[i], gs0, tl, pp, pathway)
    data.frame(obs = i, elapsed_s = 150 * (i - 1), co2_r_sp_umol_mol = SETPOINTS[i],
               co2_r_umol_mol = SETPOINTS[i] + stats::rnorm(1, 0, 0.3),
               ca_umol_mol = s[["ca"]], ci_umol_mol = s[["ci"]] + stats::rnorm(1, 0, noise_ci),
               a_umol_m2_s = s[["a"]] + stats::rnorm(1, 0, noise_a), gsw_mol_m2_s = s[["gsw"]],
               tleaf_c = tl, qin_umol_m2_s = 1800 + stats::rnorm(1, 0, 1), pa_kpa = 99.5,
               delta_pcham_kpa = 0.1, oxygen_pct = 21, stable_n = 2, stable_total = 2)
  })
  d <- do.call(rbind, rows)
  d$co2_s_umol_mol <- d$ca_umol_mol
  d$qc_hard_fail <- FALSE
  d$qc_unstable <- FALSE
  d
}

C3_TRUTH <- data.frame(
  genotype = c("G1", "G1", "G2", "G2"), plot = c("101", "102", "201", "202"),
  instrument_id = c("68C-000101/68H-000101", "68C-000202/68H-000202",
                    "68C-000101/68H-000101", "68C-000202/68H-000202"),
  Vcmax = c(95, 105, 75, 80), J = c(175, 185, 150, 155), Tp = c(11.5, 12.5, 13, 13),
  RL = c(1.4, 1.6, 1.2, 1.3), tleaf = c(28, 28.5, 28, 28.5), gs0 = c(0.35, 0.40, 0.30, 0.32),
  after = c(1, 1, 1, 0.78), stringsAsFactors = FALSE)

C4_TRUTH <- data.frame(
  genotype = c("maize_A", "maize_B"), plot = c("301", "302"),
  instrument_id = c("68C-000101/68H-000101", "68C-000202/68H-000202"),
  Vcmax = c(40, 34), Vpmax = c(130, 110), RL = c(1.0, 0.8), tleaf = c(30, 30),
  gs0 = c(0.20, 0.18), after = c(1, 1), stringsAsFactors = FALSE)

make_set <- function(truth, pathway) {
  parts <- lapply(seq_len(nrow(truth)), function(i) {
    t <- truth[i, ]
    p <- if (pathway == "C3") list(Vcmax = t$Vcmax, J = t$J, Tp = t$Tp, RL = t$RL) else
      list(Vcmax = t$Vcmax, Vpmax = t$Vpmax, RL = t$RL)
    d <- sim_curve(p, t$tleaf, t$gs0, pathway, after = t$after)
    data.frame(source_file = "demo_synthetic", instrument_id = t$instrument_id,
               genotype = t$genotype, plot = t$plot, d, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, parts)
  num <- vapply(out, is.numeric, logical(1))
  out[num] <- lapply(out[num], function(v) round(v, 4))
  out
}

main <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else "."
  set.seed(20261015)
  utils::write.csv(make_set(C3_TRUTH, "C3"), file.path(here, "aci_demo_c3.csv"), row.names = FALSE)
  utils::write.csv(make_set(C4_TRUTH, "C4"), file.path(here, "aci_demo_c4.csv"), row.names = FALSE)
  cat("Wrote aci_demo_c3.csv and aci_demo_c4.csv to", here, "\n")
}

if (sys.nframe() == 0L) main()
