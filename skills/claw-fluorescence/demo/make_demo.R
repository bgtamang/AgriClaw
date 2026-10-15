# Generates the synthetic demo data for claw-fluorescence (fixed seed).
# Run from anywhere: Rscript demo/make_demo.R
#
# Three genotypes x 4 plants. Each plant has one dark-adapted Fo/Fm and
# light-adapted flashes at five PPFD levels. Signals are built from known
# parameters so the true values are exactly defined:
#   Fm' = Fm / (1 + NPQ),  Fo' = Fo / (Fv/Fm + Fo/Fm')  (Oxborough & Baker 1997),
#   Fs = Fm' - qP (Fm' - Fo'),
# with NPQ rising and qP falling with light; 1 % multiplicative noise is added.
# fo_prime is "measured" (with noise) for G1 and G2 and missing for G3.
# Planted problems: plant G1-4 has an underestimated Fo (Fv/Fm about 0.89), and
# plant G2-3's "dark-adapted" Fm was taken before full dark adaptation, so
# Fm' > Fm at low light (negative NPQ).

set.seed(20261015)
out_dir <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else "."
})

ppfd <- c(100, 300, 600, 1000, 1500)
geno <- list(G1 = c(fvfm = 0.83, npq_max = 2.2, qp_drop = 0.45),
             G2 = c(fvfm = 0.82, npq_max = 1.6, qp_drop = 0.40),
             G3 = c(fvfm = 0.79, npq_max = 2.8, qp_drop = 0.55))
noise <- function(x) x * exp(rnorm(length(x), 0, 0.01))

rows <- list()
for (g in names(geno)) {
  p <- geno[[g]]
  for (i in 1:4) {
    fo <- runif(1, 280, 340)
    fm <- fo / (1 - (p[["fvfm"]] + rnorm(1, 0, 0.008)))
    npq <- 0.05 + p[["npq_max"]] * ppfd / (ppfd + 500)
    qp <- 0.97 - p[["qp_drop"]] * ppfd / 1500
    fm_p <- fm / (1 + npq)
    fo_p <- fo / ((fm - fo) / fm + fo / fm_p)
    fs <- fm_p - qp * (fm_p - fo_p)
    fo_rec <- fo; fm_rec <- fm
    if (g == "G1" && i == 4) fo_rec <- fo * 0.62                 # Fo underestimated
    if (g == "G2" && i == 3) fm_rec <- fm * 0.70                 # Fm not fully dark-adapted
    rows[[length(rows) + 1]] <- data.frame(
      genotype = g, plant_id = paste0(g, "-", i), ppfd_umol_m2_s = ppfd,
      fo = round(noise(fo_rec)), fm = round(noise(fm_rec)),
      fs = round(noise(fs), 1), fm_prime = round(noise(fm_p), 1),
      fo_prime = if (g == "G3") NA_real_ else round(noise(fo_p), 1))
  }
}
dat <- do.call(rbind, rows)
# One dark-adapted Fo/Fm per plant: the same values on every row of that plant
for (pid in unique(dat$plant_id)) {
  k <- dat$plant_id == pid
  dat$fo[k] <- dat$fo[k][1]; dat$fm[k] <- dat$fm[k][1]
}
utils::write.csv(dat, file.path(out_dir, "fluorescence.csv"), row.names = FALSE)
cat("Wrote fluorescence.csv to", out_dir, "\n")
