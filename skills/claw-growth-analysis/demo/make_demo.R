# Generates the synthetic demo data for claw-growth-analysis (fixed seed).
# Run from anywhere: Rscript demo/make_demo.R
#
# Three genotypes, five destructive harvests (days 7-35 after transplanting),
# six plants per harvest. Total dry mass follows logistic growth
#   W(t) = K / (1 + ((K - W0) / W0) exp(-r t)),  W0 = 0.08 g at t = 0,
# so the true RGR(t) = r (1 - W(t) / K) declines with time:
#   A: r = 0.20 d-1, K = 40 g;  B: r = 0.17 d-1, K = 40 g;  C: r = 0.22 d-1, K = 25 g.
# Plant-to-plant variation is log-normal, with SD of ln mass rising from 0.15
# to 0.35 (as in real harvests), which biases RGR computed from ln(mean mass).
# Leaf mass fraction falls from about 0.55 to 0.41 and SLA from about 30 to
# 21 m2 kg-1 over the experiment.

set.seed(20261015)
out_dir <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else "."
})

logistic <- function(t, r, K, w0 = 0.08) K / (1 + ((K - w0) / w0) * exp(-r * t))
pars <- list(A = c(r = 0.20, K = 40, sla_shift = 0), B = c(r = 0.17, K = 40, sla_shift = 25),
             C = c(r = 0.22, K = 25, sla_shift = -20))
days <- c(7, 14, 21, 28, 35)
sd_ln <- stats::setNames(seq(0.15, 0.35, length.out = length(days)), days)

rows <- list()
for (g in names(pars)) {
  p <- pars[[g]]
  for (t in days) {
    n <- 6
    w <- logistic(t, p["r"], p["K"]) * exp(rnorm(n, 0, sd_ln[as.character(t)]))
    lmf <- pmin(pmax(0.57 - 0.0045 * t + rnorm(n, 0, 0.02), 0.2), 0.8)
    sla_cm2_g <- (320 - 3.2 * t + p["sla_shift"]) * exp(rnorm(n, 0, 0.06))
    leaf_mass <- lmf * w
    rows[[length(rows) + 1]] <- data.frame(
      genotype = g, plant_id = sprintf("%s-d%02d-%d", g, t, seq_len(n)), days = t,
      total_dry_mass_g = round(w, 4), leaf_dry_mass_g = round(leaf_mass, 4),
      leaf_area_cm2 = round(leaf_mass * sla_cm2_g, 2))
  }
}
harvests <- do.call(rbind, rows)
utils::write.csv(harvests, file.path(out_dir, "harvests.csv"), row.names = FALSE)
cat("Wrote harvests.csv to", out_dir, "\n")
