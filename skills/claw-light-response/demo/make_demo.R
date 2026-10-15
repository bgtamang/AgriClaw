# Generate the synthetic demo curves for claw-light-response (fixed seed).
#   Rscript demo/make_demo.R
# Writes demo/light_demo.csv in the tidy layout of claw-licor-6800. Curves
# follow the non-rectangular hyperbola with the parameters below plus
# Gaussian noise (sd 0.25 umol m-2 s-1). Curve G2 - 202 has photoinhibition
# added above 1000 umol m-2 s-1 (A reduced by 0.005 x (Q - 1000)).

Q_SEQ <- c(2000, 1500, 1200, 1000, 800, 600, 400, 300, 200, 150, 100, 75, 50, 25, 0)

nrh_net <- function(q, phi, theta, pmax, rd) {
  x <- phi * q + pmax
  (x - sqrt(x^2 - 4 * theta * phi * q * pmax)) / (2 * theta) - rd
}

TRUTH <- data.frame(
  genotype = c("G1", "G1", "G2", "G2"), plot = c("101", "102", "201", "202"),
  instrument_id = c("68C-000101/68H-000101", "68C-000202/68H-000202",
                    "68C-000101/68H-000101", "68C-000202/68H-000202"),
  phi = c(0.060, 0.058, 0.052, 0.055), theta = c(0.75, 0.80, 0.70, 0.72),
  pmax = c(32, 35, 24, 27), rd = c(1.5, 1.7, 1.2, 1.3), photoinhibition = c(0, 0, 0, 0.005),
  stringsAsFactors = FALSE)

make_demo <- function() {
  parts <- lapply(seq_len(nrow(TRUTH)), function(i) {
    t <- TRUTH[i, ]
    q <- Q_SEQ + stats::rnorm(length(Q_SEQ), 0, 0.5)
    q[Q_SEQ == 0] <- abs(q[Q_SEQ == 0])
    a <- nrh_net(q, t$phi, t$theta, t$pmax, t$rd) - t$photoinhibition * pmax(q - 1000, 0) +
      stats::rnorm(length(q), 0, 0.25)
    data.frame(source_file = "demo_synthetic", instrument_id = t$instrument_id, genotype = t$genotype,
               plot = t$plot, obs = seq_along(q), elapsed_s = 180 * (seq_along(q) - 1),
               qin_umol_m2_s = round(q, 2), a_umol_m2_s = round(a, 4),
               tleaf_c = round(25 + stats::rnorm(length(q), 0, 0.15), 3), ca_umol_mol = 400,
               stable_n = 2, stable_total = 2, qc_hard_fail = FALSE, qc_unstable = FALSE,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, parts)
}

main <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else "."
  set.seed(20261015)
  utils::write.csv(make_demo(), file.path(here, "light_demo.csv"), row.names = FALSE)
  cat("Wrote light_demo.csv to", here, "\n")
}

if (sys.nframe() == 0L) main()
