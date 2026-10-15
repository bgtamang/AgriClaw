# Generates the synthetic demo inputs for claw-plot-grid (fixed seed). Run from
# anywhere with: Rscript demo/make_demo.R
#
#   corners.csv   the 4 outer corners (lon, lat; WGS84) of a block of 6 ranges x
#                 8 columns of 2-row plots (5.2 m long, 2 x 0.76 m wide, 0.9 m
#                 alleys), planted at a true bearing of 12 degrees with columns to
#                 the right, plus 1 cm of simulated RTK noise per coordinate.
#   fieldbook.csv an RCBD field book: 12 entries x 4 reps = 48 plots, numbered
#                 serpentine along planter passes from corner c1, with the field
#                 position of each plot (row = range from c1, col = pass from c1).
#
# The location (near 41.30 N, 95.80 W, UTM zone 15N) is arbitrary farmland; the
# field and the entries are invented.

make_demo <- function(dir, seed = 2026L) {
  set.seed(seed)
  n_ranges <- 6L; n_cols <- 8L
  L <- 5.2; A <- 0.9; W <- 2 * 0.76
  Ltot <- n_ranges * L + (n_ranges - 1) * A
  Wtot <- n_cols * W
  origin <- c(-95.8000, 41.3000)
  bearing_true <- 12

  # Project the origin to UTM 15N and find grid north and the scale factor there.
  utm <- sf::st_crs(32615)
  to_utm <- function(lon, lat) {
    p <- sf::st_sfc(lapply(seq_along(lon), function(i) sf::st_point(c(lon[i], lat[i]))), crs = 4326)
    sf::st_coordinates(sf::st_transform(p, utm))[, 1:2, drop = FALSE]
  }
  d <- 1e-4
  q <- to_utm(origin[1] + c(-d, d, 0, 0), origin[2] + c(0, 0, -d, d))
  north_grid <- atan2(q[4, 1] - q[3, 1], q[4, 2] - q[3, 2]) * 180 / pi
  phi <- origin[2] * pi / 180
  a <- 6378137; f <- 1 / 298.257223563; e2 <- f * (2 - f)
  w <- sqrt(1 - e2 * sin(phi)^2)
  k_e <- sqrt(sum((q[2, ] - q[1, ])^2)) / (a / w * cos(phi) * 2 * d * pi / 180)
  k_n <- sqrt(sum((q[4, ] - q[3, ])^2)) / (a * (1 - e2) / w^3 * 2 * d * pi / 180)
  k <- (k_e + k_n) / 2

  b <- (bearing_true + north_grid) * pi / 180
  along <- c(sin(b), cos(b))
  right <- c(cos(b), -sin(b))
  O <- to_utm(origin[1], origin[2])[1, ]
  corners <- rbind(O, O + k * Ltot * along, O + k * Ltot * along + k * Wtot * right,
                   O + k * Wtot * right)
  corners <- corners + matrix(stats::rnorm(8, 0, 0.01), 4, 2)
  ll <- sf::st_coordinates(sf::st_transform(
    sf::st_sfc(lapply(1:4, function(i) sf::st_point(corners[i, ])), crs = utm), 4326))
  cdf <- data.frame(corner = 1:4, lon = round(ll[, 1], 9), lat = round(ll[, 2], 9))

  # RCBD: 12 entries randomised within each rep; rep r occupies passes 2r-1 and 2r.
  plot <- seq_len(n_ranges * n_cols)
  pass <- (plot - 1L) %/% n_ranges + 1L
  pos <- (plot - 1L) %% n_ranges + 1L
  row <- ifelse(pass %% 2L == 1L, pos, n_ranges + 1L - pos)
  rep <- (pass + 1L) %/% 2L
  entry <- unlist(lapply(1:4, function(r) sample(1:12)))
  fb <- data.frame(plot = plot, rep = rep, entry = entry,
                   genotype = sprintf("G%02d", entry), row = row, col = pass)

  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(cdf, file.path(dir, "corners.csv"), row.names = FALSE)
  utils::write.csv(fb, file.path(dir, "fieldbook.csv"), row.names = FALSE)
  invisible(list(corners = cdf, fieldbook = fb))
}

if (sys.nframe() == 0L) {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else getwd()
  make_demo(here)
  cat("Wrote corners.csv and fieldbook.csv to", here, "\n")
}
