# Tests for claw-plot-grid. Run from anywhere: Rscript tests/test_plot_grid.R
# Known answers: hand-built numbering matrices, exact plot geometry in a local
# metric frame, UTM zones of known places, the transverse Mercator series for
# meridian convergence and for the point scale factor (Snyder 1987, eq. 8-11),
# ellipsoidal ground distances (terra uses GeographicLib), corner fits of known
# layouts, and field-book joins.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "plot_grid.R"))  # loads functions; main() does not run

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
expect_error <- function(expr, pattern, what) {
  msg <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  check(!is.na(msg) && grepl(pattern, msg, ignore.case = TRUE),
        paste0(what, " (got: ", msg, ")"))
}
run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "plot_grid.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
base_opt <- function(...) {
  o <- DEFAULTS
  o$ranges <- 4; o$columns <- 3; o$plot_length <- 5; o$alley <- 1
  o$plot_width <- 1.5; o$col_gap <- 0.25; o$units <- "m"
  mods <- list(...)
  for (k in names(mods)) o[[k]] <- mods[[k]]
  o
}
ring_xy <- function(g, i) sf::st_coordinates(g[i])[1:4, 1:2]

# 1. UTM zones of known places, including the Norway and Svalbard exceptions
check(utm_epsg(-88.2, 40.1) == 32616L, "Urbana IL is UTM 16N")
check(utm_epsg(-93.7, 41.7) == 32615L, "central Iowa is UTM 15N")
check(utm_epsg(151.2, -33.9) == 32756L, "Sydney is UTM 56S")
check(utm_epsg(-0.5, 51.5) == 32630L, "London is UTM 30N")
check(utm_epsg(5.3, 60.4) == 32632L, "Bergen uses the 32V exception")
check(utm_epsg(8.0, 78.0) == 32631L, "Svalbard 0-9E is zone 31X")
check(utm_epsg(15.6, 78.2) == 32633L, "Longyearbyen is zone 33X")
expect_error(utm_epsg(10, 85), "80 deg S and 84 deg N", "UTM refuses polar latitudes")

# 2. Plot numbering: hand-built matrices ([range, column]) for 4 ranges x 3 columns
m <- plot_numbers(4, 3, "serpentine", "along", 1, 1)
check(identical(m, matrix(c(1:4, 8:5, 9:12), 4, 3)), "serpentine along a pass from c1")
m <- plot_numbers(4, 3, "rowwise", "along", 1, 1)
check(identical(m, matrix(1:12, 4, 3)), "row-wise along a pass from c1")
m <- plot_numbers(4, 3, "serpentine", "across", 1, 1)
check(identical(m, rbind(1:3, 6:4, 7:9, 12:10)), "serpentine across ranges from c1")
m <- plot_numbers(4, 3, "serpentine", "along", 3, 1)
check(identical(m, matrix(c(12:9, 5:8, 4:1), 4, 3)), "serpentine along a pass from c3")
m <- plot_numbers(4, 3, "rowwise", "across", 2, 101)
check(identical(m, rbind(110:112, 107:109, 104:106, 101:103)), "row-wise across from c2, first plot 101")
check(all(sort(as.vector(plot_numbers(7, 5, "serpentine", "along", 4, 1))) == 1:35),
      "every plot number used exactly once")

# 3. Exact geometry in a local metric frame (no projection): bearing 30, right
opt <- base_opt(crs = "local", origin = "1000,2000", bearing = 30, across = "right")
res <- build_grid(opt)
g <- res$geom
at <- res$attrs
check(isTRUE(all.equal(as.numeric(sf::st_area(g)), rep(7.5, 12), tolerance = 1e-12)),
      "plot areas equal length x width")
a <- c(sin(pi / 6), cos(pi / 6)); cv <- c(cos(pi / 6), -sin(pi / 6))
exp_x <- 1000 + ((at$range - 1) * 6 + 2.5) * a[1] + ((at$column - 1) * 1.75 + 0.75) * cv[1]
exp_y <- 2000 + ((at$range - 1) * 6 + 2.5) * a[2] + ((at$column - 1) * 1.75 + 0.75) * cv[2]
check(max(abs(at$x - exp_x), abs(at$y - exp_y)) < 1e-9, "centroids at the exact layout positions")
c11 <- at[at$range == 1 & at$column == 1, ]; c21 <- at[at$range == 2 & at$column == 1, ]
c12 <- at[at$range == 1 & at$column == 2, ]
check(abs(sqrt((c21$x - c11$x)^2 + (c21$y - c11$y)^2) - 6) < 1e-9, "range spacing = length + alley")
check(abs(sqrt((c12$x - c11$x)^2 + (c12$y - c11$y)^2) - 1.75) < 1e-9, "column spacing = width + gap")
r1 <- ring_xy(g, 1)
edge_bearing <- atan2(r1[2, 1] - r1[1, 1], r1[2, 2] - r1[1, 2]) * 180 / pi
check(abs(edge_bearing - 30) < 1e-9, "long plot edges follow the bearing")
check(abs(sqrt(sum((r1[2, ] - r1[1, ])^2)) - 5) < 1e-9 && abs(sqrt(sum((r1[3, ] - r1[2, ])^2)) - 1.5) < 1e-9,
      "plot edges are 5 m x 1.5 m")
check(res$overlap$area == 0, "no plot overlap")
d12 <- c(c12$x - c11$x, c12$y - c11$y)
check(a[1] * d12[2] - a[2] * d12[1] < 0, "columns increase to the right of the planting direction")
resL <- build_grid(base_opt(crs = "local", origin = "0,0", bearing = 30, across = "left"))
l11 <- resL$attrs[resL$attrs$range == 1 & resL$attrs$column == 1, ]
l12 <- resL$attrs[resL$attrs$range == 1 & resL$attrs$column == 2, ]
check(a[1] * (l12$y - l11$y) - a[2] * (l12$x - l11$x) > 0, "--across left puts columns on the left")
res0 <- build_grid(base_opt(crs = "local", origin = "0,0", bearing = 0, alley = 0, col_gap = 0))
check(res0$overlap$area == 0 && res0$overlap$pairs == 0L, "touching plots (no alleys) do not overlap")
check(abs(as.numeric(sf::st_area(sf::st_union(res0$geom))) - 12 * 7.5) < 1e-9,
      "touching plots tile the block exactly")
resT <- build_grid(base_opt(crs = "local", origin = "0,0", bearing = 0, trim_ends = 0.5, trim_sides = 0.2))
check(isTRUE(all.equal(as.numeric(sf::st_area(resT$inner)), rep(4 * 1.1, 12))), "trimmed plots are 4 x 1.1 m")
check(all(sf::st_within(resT$inner, resT$geom, sparse = FALSE)[cbind(1:12, 1:12)]), "trimmed plots lie inside their plots")
resO <- build_grid(base_opt(crs = "local", origin = "0,0", bearing = 0, offset_along = 2, offset_across = 1))
o11 <- resO$attrs[resO$attrs$plot == 1, ]
check(abs(o11$x - (1 + 0.75)) < 1e-9 && abs(o11$y - (2 + 2.5)) < 1e-9, "offsets shift the first plot from the origin")
resF <- build_grid(base_opt(crs = "local", origin = "0,0", bearing = 0, units = "ft"))
check(isTRUE(all.equal(as.numeric(sf::st_area(resF$geom))[1], 5 * 1.5 * 0.3048^2)), "feet are converted exactly (0.3048 m)")

# 4. Lon/lat input: UTM zone, convergence, scale factor, ground distances, round trip
lon0 <- -95.8; lat0 <- 41.3
optG <- base_opt(origin = sprintf("%.6f,%.6f", lon0, lat0), bearing = 12, across = "right",
                 ranges = 6, columns = 8, plot_length = 5.2, alley = 0.9, plot_width = 1.52, col_gap = 0)
resG <- build_grid(optG)
check(resG$info$crs_out$epsg == 32615L, "lon/lat input goes to the local UTM zone (15N)")
phi <- lat0 * pi / 180; dl <- (lon0 - (-93)) * pi / 180
e2 <- WGS84_F * (2 - WGS84_F); ep2 <- e2 / (1 - e2); eta2 <- ep2 * cos(phi)^2
gamma <- dl * sin(phi) * (1 + dl^2 * cos(phi)^2 / 3 * (1 + 3 * eta2 + 2 * eta2^2))   # TM convergence series
check(abs(resG$info$gf$north_bearing - (-gamma * 180 / pi)) < 1e-4,
      sprintf("grid bearing of true north (%.5f) matches the convergence series (%.5f)",
              resG$info$gf$north_bearing, -gamma * 180 / pi))
A <- dl * cos(phi); Tt <- tan(phi)^2
k_snyder <- 0.9996 * (1 + (1 + eta2) * A^2 / 2 + (5 - 4 * Tt + 42 * eta2 + 13 * eta2^2 - 28 * ep2) * A^4 / 24)
check(abs(resG$info$k - k_snyder) < 2e-8, sprintf("scale factor %.9f matches Snyder eq. 8-11 (%.9f)", resG$info$k, k_snyder))
check(abs(geo_factors(-93, 40, sf::st_crs(32615))$k - 0.9996) < 1e-9, "scale factor is 0.9996 on the central meridian")
atG <- resG$attrs
pairs_d <- function(i, j) terra::distance(cbind(atG$lon[i], atG$lat[i]), cbind(atG$lon[j], atG$lat[j]), lonlat = TRUE)
i11 <- which(atG$range == 1 & atG$column == 1); i61 <- which(atG$range == 6 & atG$column == 1)
i21 <- which(atG$range == 2 & atG$column == 1); i18 <- which(atG$range == 1 & atG$column == 8)
check(abs(pairs_d(i11, i21) - 6.1) < 1e-3, sprintf("ground range spacing 6.100 m (got %.4f)", pairs_d(i11, i21)))
check(abs(pairs_d(i11, i61) - 5 * 6.1) < 1e-3, sprintf("ground distance over 5 ranges 30.500 m (got %.4f)", pairs_d(i11, i61)))
check(abs(pairs_d(i11, i18) - 7 * 1.52) < 1e-3, sprintf("ground distance over 7 columns 10.640 m (got %.4f)", pairs_d(i11, i18)))
M <- WGS84_A * (1 - e2) / (1 - e2 * sin(phi)^2)^1.5; N <- WGS84_A / sqrt(1 - e2 * sin(phi)^2)
az <- atan2((atG$lon[i61] - atG$lon[i11]) * pi / 180 * N * cos(phi), (atG$lat[i61] - atG$lat[i11]) * pi / 180 * M) * 180 / pi
check(abs(az - 12) < 1e-3, sprintf("true bearing of the planting direction is 12 deg (got %.5f)", az))
ll <- sf::st_transform(resG$geom, 4326)
back <- sf::st_transform(ll, 32615)
check(max(abs(sf::st_coordinates(back)[, 1:2] - sf::st_coordinates(resG$geom)[, 1:2])) < 1e-6,
      "UTM -> lon/lat -> UTM round trip within 1 micrometre")
optGrid <- optG; optGrid$bearing <- resG$info$bearing_grid; optGrid$bearing_ref <- "grid"
check(max(abs(build_grid(optGrid)$attrs$x - atG$x)) < 1e-6, "--bearing-ref grid with the grid bearing gives the same plots")

# 5. Corners mode
blk <- function(res) {           # outer block corners c1..c4 of a result
  at <- res$attrs; nr <- max(at$range); nc <- max(at$column)
  pick <- function(r, cl, k) ring_xy(res$geom, which(at$range == r & at$column == cl))[k, ]
  rbind(pick(1, 1, 1), pick(nr, 1, 2), pick(nr, nc, 3), pick(1, nc, 4))
}
P <- blk(resG)
Pll <- transform_xy(P, sf::st_crs(32615), sf::st_crs(4326))
corner_str <- paste(apply(Pll, 1, function(z) sprintf("%.10f,%.10f", z[1], z[2])), collapse = ";")
optC <- optG; optC$origin <- ""; optC$corners <- corner_str
resC <- build_grid(optC)
check(max(resC$corner_df$misfit_m) < 1e-4, "exact corners fit with no misfit")
check(abs(resC$info$bearing_true - 12) < 1e-5, "bearing recovered from the corners")
check(max(abs(resC$attrs$x - atG$x), abs(resC$attrs$y - atG$y)) < 1e-4, "corner mode reproduces origin-mode plots")
check(resC$info$across == "right", "handedness read from the corner order")
set.seed(7)
Pn <- P + matrix(rnorm(8, 0, 0.01), 4, 2)
optN <- optC; optN$corners <- paste(apply(Pn, 1, function(z) sprintf("%.4f,%.4f", z[1], z[2])), collapse = ";")
optN$crs <- "EPSG:32615"
resN <- build_grid(optN)
check(max(resN$corner_df$misfit_m) < 0.05 && abs(resN$info$bearing_true - 12) < 0.05,
      "1 cm corner noise: small misfit, bearing within 0.05 deg")
optU <- optC; optU$units <- "ft"   # metre values read as feet
expect_error(build_grid(optU), "units", "layout/corner mismatch (ft vs m) stops with a units hint")
optR <- optC; optR$ranges <- 7
expect_error(build_grid(optR), "disagree", "wrong number of ranges is caught by the corners")
optX <- optC; optX$corners <- paste(strsplit(corner_str, ";")[[1]][c(1, 3, 2, 4)], collapse = ";")
expect_error(build_grid(optX), "convex quadrilateral", "corners out of order are refused")
Q <- rbind(c(0, 0), c(0, 30), c(14, 31), c(13, 0))          # trapezoid, local metres
optB <- base_opt(crs = "local", corners = paste(apply(Q, 1, paste, collapse = ","), collapse = ";"),
                 ranges = 5, columns = 4, plot_length = 5.2, alley = 1, plot_width = 3, col_gap = 0.3,
                 fit_corners = TRUE)
resB <- build_grid(optB)
check(max(abs(blk(resB) - Q)) < 1e-9, "--fit-corners puts the outer plot corners on the given corners")
check(resB$overlap$area == 0, "bilinear plots do not overlap")

# 6. Field book join
fb <- data.frame(plot = 1:12, genotype = sprintf("G%02d", c(3, 1, 2, 4, 6, 5, 9, 7, 8, 12, 10, 11)),
                 range = 99, stringsAsFactors = FALSE)
j <- join_fieldbook(res$attrs, fb, "plot")
check(identical(j$data$genotype[match(1:12, j$data$plot)], fb$genotype), "genotypes joined to the right plots")
check("fb_range" %in% names(j$data) && j$renamed == "range", "clashing field-book column renamed to fb_range")
check(j$n_unmatched == 0L, "all plots matched")
j2 <- join_fieldbook(res$attrs, fb[-(1:2), ], "plot")
check(j2$n_unmatched == 2L && sum(!j2$data$matched) == 2L, "plots without a field-book row are counted")
expect_error(join_fieldbook(res$attrs, rbind(fb, fb[1, ]), "plot"), "duplicate", "duplicate plots refused")
fbx <- fb; fbx$plot[12] <- 99
expect_error(join_fieldbook(res$attrs, fbx, "plot"), "match no plot", "unknown plot numbers refused")
fbrc <- data.frame(ROW = rep(1:4, 3), COL = rep(1:3, each = 4), entry = 1:12)
j3 <- join_fieldbook(res$attrs, fbrc, "ROW=range,COL=column")
check(all(j3$data$entry == (j3$data$column - 1) * 4 + j3$data$range), "join on range and column")
fbp <- data.frame(plot = 1:12, row = 0, col = 0)
for (p in 1:12) { w <- which(plot_numbers(4, 3, "serpentine", "along", 1, 1) == p, arr.ind = TRUE)
  fbp$row[p] <- w[1]; fbp$col[p] <- w[2] }
check(check_fb_coords(join_fieldbook(res$attrs, fbp, "plot")$data, c("row", "col"))$ok,
      "field-book positions agree with the matching numbering")
resRW <- build_grid(base_opt(crs = "local", origin = "0,0", bearing = 0, numbering = "rowwise"))
check(!check_fb_coords(join_fieldbook(resRW$attrs, fbp, "plot")$data, c("row", "col"))$ok,
      "field-book positions expose a serpentine/row-wise mismatch")

# 7. Input validation
expect_error(build_grid(base_opt(crs = "local", origin = "0,0", bearing = 0, units = "")), "units is required",
             "missing --units is refused")
expect_error(build_grid(base_opt(crs = "local", origin = "0,0", bearing = 0, alley = -1)), "negative",
             "negative alley refused")
expect_error(build_grid(base_opt(crs = "local", origin = "0,0", bearing = 0, trim_ends = 3)), "whole plot length",
             "trimming the whole plot refused")
expect_error(build_grid(base_opt(origin = "41.3,-95.8", bearing = 0)), "lat,lon|out of range",
             "swapped lat,lon refused")
expect_error(build_grid(base_opt(origin = "-95.8,41.3", bearing = 0, out_crs = "EPSG:4326")), "metres|degrees",
             "lon/lat output CRS refused")
check(length(dimension_warnings(list(L = 5, W = 1.5, A = 1, row_spacing = 30), "m")) == 1L,
      "a 30 m row spacing triggers a units warning")

# 8. End to end on the demo data
out <- file.path(tempdir(), "claw-plot-grid-test")
unlink(out, recursive = TRUE)
cli <- run_cli(c("--demo", "--out", shQuote(out)))
check(cli$status == 0L, paste("demo run exits cleanly:", paste(cli$output, collapse = "\n")))
for (f in c("report.md", "spatial/plots.gpkg", "spatial/plots_inner.gpkg", "tables/plots.csv",
            "tables/corner_check.csv", "figures/plot_map.png", "reproducibility/commands.sh",
            "reproducibility/sessionInfo.txt"))
  check(file.exists(file.path(out, f)), paste("demo writes", f))
gp <- sf::st_read(file.path(out, "spatial", "plots.gpkg"), quiet = TRUE)
check(nrow(gp) == 48L && sf::st_crs(gp)$epsg == 32615L, "demo GeoPackage has 48 plots in UTM 15N")
check(all(c("plot", "range", "column", "genotype", "rep") %in% names(gp)), "demo GeoPackage carries field-book columns")
fbd <- utils::read.csv(file.path(SKILL_TEST_DIR, "demo", "fieldbook.csv"))
check(identical(gp$genotype[order(gp$plot)], fbd$genotype[order(fbd$plot)]), "demo genotypes joined by plot")
cc <- utils::read.csv(file.path(out, "tables", "corner_check.csv"))
check(max(cc$misfit_m) < 0.05, "demo corners (1 cm noise) fit within 5 cm")
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
expect_error(main(c("--demo", "--colour", "red")), "Unknown option", "unknown option is rejected")
expect_error(main(c("--origin", "0,0", "--crs", "local", "--bearing", "0", "--ranges", "2", "--columns", "2",
                    "--plot-length", "5", "--alley", "1", "--plot-width", "1.5", "--out", out)),
             "units is required", "command line without --units is rejected")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
