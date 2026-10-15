# Generates the synthetic demo data for claw-vpd-transpiration (fixed seed).
# Run from anywhere: Rscript demo/make_demo.R
#
# 12 sealed pots in 2 growth chambers, logged every minute from 04:00 to 14:00:
#   04:00-05:00  night step (dark, 20 degC)
#   05:00-07:00  lights on, warm-up (not analyzed)
#   07:00-14:00  seven 1-h daytime VPD steps at 30 degC (0.9 to 3.9 kPa)
# True responses (TR in mg H2O m-2 s-1, VPD in kPa):
#   G1  two segments, breakpoint 2.0 kPa, slopes 32 then 6
#   G2  two segments, breakpoint 2.6 kPa, slopes 30 then 10
#   G3  linear, slope 28
# Planted problems: P07 irrigated (+150 g) at 10:20, a one-reading +25 g spike
# on P02 at 08:41, and a 20-min logger dropout on P10 from 12:30.

set.seed(20261015)
out_dir <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else "."
})

svp <- function(t) 0.61375 * exp(17.502 * t / (240.97 + t))
tr_true <- function(g, vpd) {
  switch(g,
         G1 = -8 + 32 * vpd + (6 - 32) * pmax(vpd - 2.0, 0),
         G2 = -8 + 30 * vpd + (10 - 30) * pmax(vpd - 2.6, 0),
         G3 = -8 + 28 * vpd)
}

t0 <- as.POSIXct("2026-07-01 04:00:00", tz = "UTC")
minutes <- 0:599
times <- t0 + 60 * minutes
hour <- 4 + minutes / 60

# Chamber set points (VPD steps change on the hour; 5-min first-order transition)
vpd_steps <- c(0.9, 1.4, 1.9, 2.4, 2.9, 3.4, 3.9)
set_t <- ifelse(hour < 5, 20, ifelse(hour < 6, 20 + 10 * (hour - 5), 30))
set_vpd <- ifelse(hour < 5, 0.70, ifelse(hour < 7, 0.9, vpd_steps[pmax(pmin(floor(hour - 7) + 1, 7), 1)]))
smooth <- function(x, tau = 5) {
  y <- x
  for (i in 2:length(x)) y[i] <- y[i - 1] + (x[i] - y[i - 1]) / tau
  y
}

env <- list()
chamber_air <- list()
for (ch in c("A", "B")) {
  dT <- if (ch == "A") 0 else 0.3
  temp <- smooth(set_t) + dT
  rh_set <- 100 * (1 - set_vpd / svp(set_t))
  rh <- smooth(rh_set)
  chamber_air[[ch]] <- data.frame(temp = temp, rh = rh, vpd = svp(temp) * (1 - rh / 100))
  for (s in 1:2) {
    env[[length(env) + 1]] <- data.frame(
      timestamp = format(times, "%Y-%m-%d %H:%M:%S"), chamber = ch, sensor = paste0(ch, s),
      air_temp_c = round(temp + rnorm(length(times), 0, 0.15), 2),
      rh_pct = round(pmin(rh + rnorm(length(times), 0, 0.8), 100), 1))
  }
}
env <- do.call(rbind, env)

plants <- data.frame(plant_id = sprintf("P%02d", 1:12),
                     genotype = rep(c("G1", "G2", "G3"), times = 4),
                     chamber = rep(c("A", "B"), each = 6),
                     leaf_area_cm2 = round(runif(12, 350, 650)))
plants <- plants[order(plants$plant_id), ]

weights <- list()
for (i in seq_len(nrow(plants))) {
  p <- plants[i, ]
  air <- chamber_air[[p$chamber]]
  factor_p <- rnorm(1, 1, 0.06)
  tr <- ifelse(hour < 5, 7, tr_true(p$genotype, air$vpd)) * factor_p   # mg m-2 s-1
  tr <- ifelse(hour >= 5 & hour < 5.5, tr * (hour - 5) * 2, tr)        # stomata opening
  tr <- smooth(tr, tau = 3)                                            # stomatal lag
  loss_g_min <- tr * 60 / 1000 * p$leaf_area_cm2 / 10000
  w <- 2600 + rnorm(1, 0, 150) - cumsum(loss_g_min)
  # minutes after 04:00: 10:20 = 380, 08:41 = 281, 12:30-12:49 = 510-529
  if (p$plant_id == "P07") w[minutes >= 380] <- w[minutes >= 380] + 150
  w <- w + rnorm(length(w), 0, 0.03)
  if (p$plant_id == "P02") w[minutes == 281] <- w[minutes == 281] + 25
  keep <- !(p$plant_id == "P10" & minutes >= 510 & minutes < 530)
  weights[[i]] <- data.frame(plant_id = p$plant_id,
                             timestamp = format(times[keep], "%Y-%m-%d %H:%M:%S"),
                             weight_g = round(w[keep], 2))
}
weights <- do.call(rbind, weights)

steps <- data.frame(step = c("night", paste0("S", 1:7)),
                    phase = c("night", rep("day", 7)),
                    start = format(c(t0 + 10 * 60, t0 + (3 * 60 + 10 + 60 * 0:6) * 60), "%Y-%m-%d %H:%M:%S"),
                    end = format(c(t0 + 60 * 60, t0 + (4 * 60 + 60 * 0:6) * 60), "%Y-%m-%d %H:%M:%S"))

utils::write.csv(weights, file.path(out_dir, "weights.csv"), row.names = FALSE)
utils::write.csv(plants, file.path(out_dir, "plants.csv"), row.names = FALSE)
utils::write.csv(env, file.path(out_dir, "env.csv"), row.names = FALSE)
utils::write.csv(steps, file.path(out_dir, "steps.csv"), row.names = FALSE)
cat("Wrote weights.csv, plants.csv, env.csv, steps.csv to", out_dir, "\n")
