# Generates the synthetic demo data for claw-vi-timeseries (fixed seed).
# Run from the skill folder:  Rscript demo/make_demo.R
# Writes demo/vi_long.csv and demo/gdd_daily.csv. Every value is invented.

set.seed(20251015)
out_dir <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else "demo"
})

beck <- function(t, vmin, vmax, sos, rsp, eos, rau)
  vmin + (vmax - vmin) * (1 / (1 + exp(-rsp * (t - sos))) + 1 / (1 + exp(rau * (t - eos))) - 1)

# Two environments planted three weeks apart, so the same genotype reaches the
# same stage on different calendar dates but similar days after planting.
envs <- list(
  E1 = list(planting = as.Date("2025-05-08"),
            flights_dap = c(17, 24, 31, 38, 46, 52, 59, 66, 74, 80, 88, 95, 102, 109, 117, 124)),
  E2 = list(planting = as.Date("2025-05-29"),
            flights_dap = c(15, 22, 30, 36, 43, 50, 58, 64, 71, 79, 86, 93, 100, 107, 115)))

n_geno <- 12
geno <- data.frame(genotype = sprintf("G%02d", seq_len(n_geno)),
                   sos = 38 + rnorm(n_geno, 0, 2.5),          # green-up midpoint, DAP
                   eos = 96 + seq(-9, 9, length.out = n_geno), # maturity group spread, DAP
                   vmax = 0.86 + rnorm(n_geno, 0, 0.02),
                   stringsAsFactors = FALSE)

rows <- list()
for (e in names(envs)) {
  ev <- envs[[e]]
  shift <- if (e == "E2") -3 else 0   # E2 is planted later into warmer soil: faster
  for (g in seq_len(n_geno)) {
    t <- ev$flights_dap
    y <- beck(t, vmin = 0.18, vmax = geno$vmax[g], sos = geno$sos[g] + shift, rsp = 0.16,
              eos = geno$eos[g] + shift, rau = 0.13) + rnorm(length(t), 0, 0.015)
    rows[[length(rows) + 1L]] <- data.frame(
      plot_id = sprintf("%s-%02d", e, g), env = e, genotype = geno$genotype[g],
      planting_date = format(ev$planting), date = format(ev$planting + t),
      ndvi = round(y, 4), stringsAsFactors = FALSE)
  }
}
d <- do.call(rbind, rows)

# Two shadow-like drops (low outliers) and one plot that missed most flights.
i1 <- which(d$plot_id == "E1-03")[6]
d$ndvi[i1] <- round(d$ndvi[i1] - 0.25, 4)
i2 <- which(d$plot_id == "E2-07")[9]
d$ndvi[i2] <- round(d$ndvi[i2] - 0.18, 4)
sparse <- which(d$plot_id == "E2-12")
d <- d[-sparse[-c(2, 5, 8, 11, 14)], ]
utils::write.csv(d, file.path(out_dir, "vi_long.csv"), row.names = FALSE)

# Daily growing degree days (base 10 C, ceiling 30 C, deg C days) from a
# synthetic seasonal temperature curve; E2 runs 1 C warmer.
gdd_rows <- list()
for (e in names(envs)) {
  ev <- envs[[e]]
  days <- seq(ev$planting, ev$planting + max(ev$flights_dap) + 5, by = "day")
  doy <- as.numeric(format(days, "%j"))
  tmean <- 22 + 5 * sin(2 * pi * (doy - 110) / 365) + (e == "E2") + rnorm(length(days), 0, 2)
  tmax <- tmean + 6
  tmin <- tmean - 6
  gdd <- (pmin(tmax, 30) + pmax(pmin(tmin, 30), 10)) / 2 - 10
  gdd_rows[[e]] <- data.frame(env = e, date = format(days), gdd_c = round(pmax(gdd, 0), 2),
                              stringsAsFactors = FALSE)
}
utils::write.csv(do.call(rbind, gdd_rows), file.path(out_dir, "gdd_daily.csv"), row.names = FALSE)
cat("Wrote", file.path(out_dir, "vi_long.csv"), "and", file.path(out_dir, "gdd_daily.csv"), "\n")
