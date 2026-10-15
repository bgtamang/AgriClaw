# Generates the synthetic demo data for claw-yield-prediction (fixed seed).
# Run from the skill folder:  Rscript demo/make_demo.R
# Writes demo/rs_yield.csv and demo/feature_time.csv. Every value is invented.

set.seed(20251015)
out_dir <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else "demo"
})

# Four environments (2 sites x 2 years). Each has its own mean yield, a
# different sensor calibration offset (shifts every NDVI value of that
# environment) and its own season weather. The offsets and weather totals are
# unrelated to the environment means, as with 4 environments in practice: they
# identify an environment (a random split can memorize its mean) but do not
# transfer to a new one.
envs <- data.frame(env = c("IL2024", "IA2024", "IL2025", "IA2025"),
                   year = c(2024, 2024, 2025, 2025),
                   mean_yield = c(3900, 3300, 4400, 3600),
                   ndvi_offset = c(0.02, -0.02, -0.03, 0.03),
                   precip_season_mm = c(430, 380, 405, 455),
                   gdd_season_c = c(1480, 1420, 1395, 1450),
                   dap_shift = c(0, 3, -2, 2), stringsAsFactors = FALSE)

# Sixty genotypes. Maturity and early vigor raise yield and are visible to the
# sensors (senescence date, early NDVI); "other" is yield potential they miss.
n_geno <- 60
geno <- data.frame(genotype = sprintf("G%03d", seq_len(n_geno)),
                   maturity = rnorm(n_geno), vigor = rnorm(n_geno), other = rnorm(n_geno),
                   stringsAsFactors = FALSE)

rows <- vector("list", nrow(envs) * n_geno)
k <- 0L
for (e in seq_len(nrow(envs))) for (g in seq_len(n_geno)) {
  k <- k + 1L
  ev <- envs[e, ]
  gv <- geno[g, ]
  yield <- ev$mean_yield + 180 * gv$maturity + 120 * gv$vigor + 150 * gv$other +
    rnorm(1, 0, 120) + rnorm(1, 0, 220)                       # G x E + plot error
  sen <- 96 + 5 * gv$maturity + ev$dap_shift + rnorm(1, 0, 1.5)
  gu <- 38 - 1.5 * gv$vigor + ev$dap_shift + rnorm(1, 0, 1.2)
  peak <- 0.84 + 0.012 * gv$vigor + ev$ndvi_offset + rnorm(1, 0, 0.008)
  rows[[k]] <- data.frame(
    plot_id = sprintf("%s-%03d", ev$env, g), env = ev$env, year = ev$year,
    genotype = gv$genotype, yield_kg_ha = round(yield),
    greenup_dap = round(gu, 1),
    ndvi_at_dap50 = round(0.62 + 0.04 * gv$vigor + ev$ndvi_offset - 0.01 * ev$dap_shift + rnorm(1, 0, 0.02), 4),
    ndvi_at_dap60 = round(0.78 + 0.02 * gv$vigor + ev$ndvi_offset + rnorm(1, 0, 0.015), 4),
    ndvi_peak = round(peak, 4),
    senescence_dap = round(sen, 1),
    season_length_days = round(sen - gu, 1),
    auc_vi_days = round(0.95 * peak * (sen - gu) + rnorm(1, 0, 1), 2),
    precip_season_mm = ev$precip_season_mm, gdd_season_c = ev$gdd_season_c,
    stringsAsFactors = FALSE)
}
d <- do.call(rbind, rows)
utils::write.csv(d, file.path(out_dir, "rs_yield.csv"), row.names = FALSE)

# When each feature is known, in days after planting. Green-up and the early
# NDVI values come from flights up to DAP 60 (claw-vi-timeseries --cutoff 60);
# season-long metrics and weather totals are known only at the end.
ft <- data.frame(feature = c("greenup_dap", "ndvi_at_dap50", "ndvi_at_dap60", "ndvi_peak",
                             "senescence_dap", "season_length_days", "auc_vi_days",
                             "precip_season_mm", "gdd_season_c"),
                 available_t = c(60, 50, 60, 90, 125, 125, 125, 130, 130),
                 stringsAsFactors = FALSE)
utils::write.csv(ft, file.path(out_dir, "feature_time.csv"), row.names = FALSE)
cat("Wrote", file.path(out_dir, "rs_yield.csv"), "and", file.path(out_dir, "feature_time.csv"), "\n")
