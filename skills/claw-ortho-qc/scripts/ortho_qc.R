#!/usr/bin/env Rscript
# claw-ortho-qc: quality checks for an orthomosaic (and optional DSM/DTM) before
# per-plot extraction: CRS, pixel size, rotation, extent against the plot grid,
# band count and data type, value ranges (reflectance or digital numbers),
# declared NoData, NoData stored as zeros, plot coverage, saturation, DSM value
# sanity, spikes and canopy height. Each check is pass, warn or fail against a
# stated threshold. Writes a report; --strict exits with an error on any fail.
#
# Usage:
#   Rscript scripts/ortho_qc.R --demo
#   Rscript scripts/ortho_qc.R --input ortho.tif --plots plots.gpkg --expect reflectance
#   Rscript scripts/ortho_qc.R --input rgb.tif --plots plots.gpkg --dsm dsm.tif --dtm dtm.tif --expect rgb
#   Rscript scripts/ortho_qc.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT license.

SKILL   <- "claw-ortho-qc"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/ortho_qc.R [--demo] --input RASTER [--plots POLYGONS] [--dsm DSM [--dtm DTM]] [options]

  --demo          run on synthetic demo data with injected problems (written to <out>/demo_input)
  --input         orthomosaic (GeoTIFF etc.)
  --plots         plot polygons, e.g. from claw-plot-grid (enables extent, coverage and GSD-vs-plot checks)
  --layer         layer name inside --plots (default: first layer)
  --id-col        plot id column (default: plot_id, else plot, else row number)
  --dsm, --dtm    digital surface model and bare-ground terrain model (meters)
  --expect        reflectance | rgb | dn | any (default any): what the values should be
  --expect-bands  number of bands the file should have
  --scale         multiply values by this before reflectance checks (e.g. 0.0001; default 1)
  --max-value     saturation level (default: data-type maximum for integer data)
  --alpha-band    alpha band number (default: auto-detect a 0/255 band 4; 0 = none)
  --max-gsd-cm    warn if the pixel size is coarser than this (cm)
  --max-height-m  warn if canopy height (DSM - DTM) exceeds this in a plot (default 4)
  --spike-m       a DSM pixel this far from its 5 x 5 median counts as a spike (default 1)
  --sample-size   pixels sampled for value statistics (default 200000)
  --seed          random seed for that sample (default 42)
  --strict        exit with an error if any check fails
  --out           output folder (default agriclaw_output/claw-ortho-qc)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, input = "", plots = "", layer = "", id_col = "", dsm = "",
                 dtm = "", expect = "any", expect_bands = NA_real_, scale = 1,
                 max_value = NA_real_, alpha_band = NA_real_, max_gsd_cm = NA_real_,
                 max_height_m = 4, spike_m = 1, sample_size = 200000, seed = 42,
                 strict = FALSE, out = "")

# ---------------------------------------------------------------------------
# Helpers shared by all AgriClaw R scripts (copy unchanged)
# ---------------------------------------------------------------------------

parse_args <- function(argv, defaults = DEFAULTS) {
  opt <- defaults
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (!startsWith(a, "--"))
      stop("Unexpected argument '", a, "'. See --help.", call. = FALSE)
    key <- gsub("-", "_", substring(a, 3L), fixed = TRUE)
    if (!key %in% names(defaults))
      stop("Unknown option '", a, "'. See --help.", call. = FALSE)
    if (is.logical(defaults[[key]])) {
      opt[[key]] <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(argv)) stop("Option '", a, "' needs a value.", call. = FALSE)
    val <- argv[[i + 1L]]
    if (is.numeric(defaults[[key]])) {
      num <- suppressWarnings(as.numeric(val))
      if (is.na(num))
        stop("Option '", a, "' needs a number, got '", val, "'.", call. = FALSE)
      val <- num
    }
    opt[[key]] <- val
    i <- i + 2L
  }
  opt
}

# Path of this script, whether it is run with Rscript or loaded with source().
this_script <- function() {
  for (i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if (!is.null(f)) return(normalizePath(f, winslash = "/"))
  }
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) return(normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/"))
  NA_character_
}
SCRIPT_FILE <- this_script()
SKILL_DIR <- if (is.na(SCRIPT_FILE)) getwd() else dirname(dirname(SCRIPT_FILE))

write_utf8 <- function(lines, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con))
  writeLines(lines, con)
}

save_table <- function(df, out_dir, name) {
  d <- file.path(out_dir, "tables")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(d, name)
  utils::write.csv(df, path, row.names = FALSE)
  path
}

# Open a PNG device in <out>/figures; the caller draws and then calls dev.off().
open_png <- function(out_dir, name, width = 7, height = 5) {
  d <- file.path(out_dir, "figures")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  grDevices::png(file.path(d, name), width = width, height = height,
                 units = "in", res = 150)
}

file_hash <- function(path, algo) {
  if (algo == "sha256") {
    con <- file(path, open = "rb")
    on.exit(close(con))
    return(as.character(openssl::sha256(con)))
  }
  unname(tools::md5sum(path))
}

# commands.sh, sessionInfo.txt and checksums of inputs + outputs.
write_repro_bundle <- function(out_dir, argv, inputs = character()) {
  rd <- file.path(out_dir, "reproducibility")
  dir.create(rd, recursive = TRUE, showWarnings = FALSE)
  script <- if (is.na(SCRIPT_FILE)) paste0("scripts/", SKILL, ".R") else SCRIPT_FILE
  cmd <- paste(c("Rscript", shQuote(script, type = "sh"),
                 vapply(argv, shQuote, "", type = "sh")), collapse = " ")
  write_utf8(c("#!/usr/bin/env bash",
               sprintf("# %s %s, run %s", SKILL, VERSION,
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
               cmd), file.path(rd, "commands.sh"))
  write_utf8(utils::capture.output(utils::sessionInfo()),
             file.path(rd, "sessionInfo.txt"))
  outs <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  outs <- outs[!startsWith(normalizePath(outs, winslash = "/"),
                           normalizePath(rd, winslash = "/"))]
  ins <- inputs[nzchar(inputs) & file.exists(inputs)]
  files <- c(ins, outs)
  algo <- if (requireNamespace("openssl", quietly = TRUE)) "sha256" else "md5"
  hashes <- vapply(files, file_hash, "", algo = algo)
  base <- normalizePath(out_dir, winslash = "/")
  labels <- normalizePath(files, winslash = "/")
  inside <- startsWith(labels, paste0(base, "/"))
  labels[inside] <- substring(labels[inside], nchar(base) + 2L)
  write_utf8(paste0(hashes, "  ", labels),
             file.path(rd, paste0("checksums.", algo)))
}

# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

# Thresholds (rules of thumb unless stated; see SKILL.md, Expert decisions).
TH <- list(pix_shape = 0.01, px_across_fail = 1, px_across_warn = 5, refl_p99_fail = 1.5,
           refl_neg_warn = 0.01, refl_over1_warn = 0.01, zero_warn = 0.01, cov_warn = 0.95,
           cov_fail = 0.75, sat_warn = 0.005, sat_fail = 0.05, dsm_min = -500, dsm_max = 9000,
           spike_warn = 0.001, spike_fail = 0.01, relief_warn = 20, chm_neg = -0.2)
DTYPE_MAX <- c(INT1U = 255, INT2U = 65535, INT1S = 127, INT2S = 32767, INT4U = 4294967295, INT4S = 2147483647)

need_pkgs <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop("Missing R package(s): ", paste(miss, collapse = ", "), ". Install with ",
         "install.packages(c(", paste0('"', miss, '"', collapse = ", "), ")).", call. = FALSE)
  invisible(TRUE)
}

qc_add <- function(qc, check, status, value, threshold, message) {
  rbind(qc, data.frame(check = check, status = status, value = paste(value, collapse = "; "),
                       threshold = threshold, message = message, stringsAsFactors = FALSE))
}

# Read a raster, recording GDAL's warning about rotated geotransforms.
read_raster <- function(path, what) {
  if (!file.exists(path)) stop(what, " not found: ", path, call. = FALSE)
  msgs <- character()
  r <- withCallingHandlers(terra::rast(path), warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  attr(r, "read_warnings") <- msgs
  r
}

# NoData values declared in the file (GDAL metadata); NaN on integer data means none.
declared_nodata <- function(path) {
  d <- tryCatch(terra::describe(path), error = function(e) character())
  v <- regmatches(d, regexpr("NoData Value=[^ ]+", d))
  v <- unique(suppressWarnings(as.numeric(sub("NoData Value=", "", v))))
  v[is.finite(v)]
}

# Number of pixels a fully covered plot would hold, in the raster's own units.
expected_pixels <- function(plots_r, r) {
  if (is_ll(r)) return(as.numeric(sf::st_area(plots_r)) / prod(gsd_meters(r)))
  as.numeric(sf::st_area(plots_r)) / prod(terra::res(r))
}

is_ll <- function(r) isTRUE(suppressWarnings(terra::is.lonlat(r)))

gsd_meters <- function(r) {
  rs <- terra::res(r)
  if (is_ll(r)) {
    lat <- mean(as.vector(terra::ext(r))[3:4]) * pi / 180
    return(c(rs[1] * 111320 * cos(lat), rs[2] * 110574))
  }
  lu <- terra::linearUnits(r)
  if (!is.finite(lu) || lu == 0) lu <- 1
  rs * lu
}

metric_crs <- function(p, r) {
  if (!isTRUE(sf::st_is_longlat(p))) return(sf::st_crs(p))
  if (!is_ll(r) && nzchar(terra::crs(r))) return(sf::st_crs(terra::crs(r)))
  cc <- sf::st_coordinates(sf::st_transform(sf::st_centroid(sf::st_as_sfc(sf::st_bbox(p))), 4326))
  sf::st_crs((if (cc[2] >= 0) 32600 else 32700) + floor((cc[1] + 180) / 6) + 1)
}

read_plots <- function(path, layer, id_col) {
  if (!file.exists(path)) stop("Plot file not found: ", path, call. = FALSE)
  p <- if (nzchar(layer)) sf::st_read(path, layer = layer, quiet = TRUE) else sf::st_read(path, quiet = TRUE)
  if (is.na(sf::st_crs(p))) stop("The plot polygons have no CRS; assign the correct one first.", call. = FALSE)
  if (!nzchar(id_col)) id_col <- if ("plot_id" %in% names(p)) "plot_id" else if ("plot" %in% names(p)) "plot" else ""
  ids <- if (nzchar(id_col)) {
    if (!id_col %in% names(p)) stop("Plot id column '", id_col, "' not found.", call. = FALSE)
    as.character(p[[id_col]])
  } else as.character(seq_len(nrow(p)))
  list(plots = p, ids = ids)
}

# Values of all bands at sampled pixels: inside plots (cells at least half inside)
# or, without plots, on a regular grid over the whole image.
sample_pixels <- function(r, plots, n, seed) {
  if (is.null(plots)) {
    v <- terra::spatSample(r, min(n, terra::ncell(r)), method = "regular")
    return(as.data.frame(v))
  }
  per <- max(50L, ceiling(n / nrow(plots)))
  set.seed(seed)
  v <- exactextractr::exact_extract(r, plots, function(values, coverage_fraction) {
    values <- as.data.frame(values)
    ok <- which(coverage_fraction >= 0.5)
    if (length(ok) > per) ok <- sample(ok, per)
    values[ok, , drop = FALSE]
  }, progress = FALSE)
  if (is.data.frame(v)) v else do.call(rbind, v)
}

band_stats <- function(v, scale) {
  do.call(rbind, lapply(seq_along(v), function(i) {
    x <- v[[i]]
    ok <- x[!is.na(x)]
    s <- ok * scale
    q <- if (length(s)) stats::quantile(s, c(0.01, 0.5, 0.99), names = FALSE) else rep(NA_real_, 3)
    data.frame(band = names(v)[i], n = length(x), na_fraction = mean(is.na(x)),
               min = if (length(s)) min(s) else NA, p01 = q[1], median = q[2], p99 = q[3],
               max = if (length(s)) max(s) else NA, frac_negative = if (length(s)) mean(s < 0) else NA,
               frac_above_1 = if (length(s)) mean(s > 1) else NA, frac_zero = if (length(ok)) mean(ok == 0) else NA)
  }))
}

detect_alpha <- function(r, given) {
  if (!is.na(given)) {
    if (given == 0) return(NA_integer_)
    if (given < 1 || given > terra::nlyr(r)) stop("--alpha-band is not a layer of the raster.", call. = FALSE)
    return(as.integer(given))
  }
  if (terra::nlyr(r) != 4L || !terra::datatype(r)[4] %in% c("INT1U", "INT2U")) return(NA_integer_)
  v <- unique(stats::na.omit(terra::spatSample(r[[4]], min(20000, terra::ncell(r)), method = "regular")[[1]]))
  if (length(v) && all(v %in% c(0, DTYPE_MAX[[terra::datatype(r)[4]]]))) 4L else NA_integer_
}

dsm_checks <- function(qc, opt, r, plots, pm, ids, ptab) {
  dsm <- read_raster(opt$dsm, "DSM")
  if (!nzchar(terra::crs(dsm))) {
    qc <- qc_add(qc, "dsm_crs", "fail", "none", "CRS required", "The DSM has no CRS; it cannot be matched to the plots.")
    return(list(qc = qc, ptab = ptab, dsm = NULL, chm = NULL))
  }
  same <- terra::same.crs(dsm, r)
  qc <- qc_add(qc, "dsm_crs", if (same) "pass" else "warn", terra::crs(dsm, describe = TRUE)$name, "same as orthomosaic",
               if (same) "DSM and orthomosaic share a CRS." else "DSM CRS differs from the orthomosaic: plots are transformed separately; check that the two line up.")
  dt <- terra::datatype(dsm)[1]
  gd <- gsd_meters(dsm)
  qc <- qc_add(qc, "dsm_resolution", "info", sprintf("%.2f cm", 100 * gd[1]), "-",
               sprintf("DSM pixel %.2f cm (orthomosaic %.2f cm).", 100 * gd[1], 100 * gsd_meters(r)[1]))
  if (identical(dt, "INT1U")) {
    qc <- qc_add(qc, "dsm_values", "fail", dt, "elevations in meters (float or 16/32-bit)",
                 "The DSM is 8-bit: it looks like an image or hillshade, not elevations.")
    return(list(qc = qc, ptab = ptab, dsm = dsm, chm = NULL))
  }
  pd <- if (!is.null(plots)) sf::st_transform(pm, sf::st_crs(terra::crs(dsm))) else NULL
  zone <- if (!is.null(pd)) terra::ext(terra::vect(pd)) + 5 * max(terra::res(dsm)) else terra::ext(dsm)
  e1 <- as.vector(terra::ext(dsm)); e2 <- as.vector(zone)
  if (!(e2[1] < e1[2] && e2[2] > e1[1] && e2[3] < e1[4] && e2[4] > e1[3])) {
    qc <- qc_add(qc, "dsm_coverage", "fail", "no overlap", "DSM covers the plots", "The DSM does not overlap the plots.")
    return(list(qc = qc, ptab = ptab, dsm = NULL, chm = NULL))
  }
  dsm_c <- terra::crop(dsm, terra::intersect(zone, terra::ext(dsm)))
  sv <- terra::spatSample(dsm_c, min(opt$sample_size, terra::ncell(dsm_c)), method = "regular", na.rm = TRUE)[[1]]
  if (!length(sv)) {
    qc <- qc_add(qc, "dsm_values", "fail", "no data", "values present", "The DSM has no values over the plots.")
    return(list(qc = qc, ptab = ptab, dsm = dsm, chm = NULL))
  }
  rng <- range(sv)
  bad <- rng[1] < TH$dsm_min || rng[2] > TH$dsm_max
  qc <- qc_add(qc, "dsm_values", if (bad) "fail" else "pass", sprintf("%.2f to %.2f", rng[1], rng[2]),
               sprintf("within %g to %g m", TH$dsm_min, TH$dsm_max),
               if (bad) "Elevations outside any plausible range in meters: wrong units (feet?), a no-data value read as data, or not a DSM." else "Elevations plausible.")
  if (!is.null(pd)) {
    vs <- exactextractr::exact_extract(!is.na(dsm_c), pd, "sum", progress = FALSE)
    ptab$dsm_coverage <- pmin(1, vs / expected_pixels(pd, dsm_c))
    w <- which.min(ptab$dsm_coverage)
    st <- if (min(ptab$dsm_coverage) < TH$cov_fail) "fail" else if (min(ptab$dsm_coverage) < TH$cov_warn) "warn" else "pass"
    qc <- qc_add(qc, "dsm_coverage", st, sprintf("min %.3f (%s)", ptab$dsm_coverage[w], ids[w]),
                 sprintf("warn < %.2f, fail < %.2f", TH$cov_warn, TH$cov_fail),
                 sprintf("%d plot(s) with DSM coverage below %.2f.", sum(ptab$dsm_coverage < TH$cov_warn), TH$cov_warn))
    med <- terra::focal(dsm_c, w = 5, fun = "median", na.rm = TRUE)
    spike <- abs(dsm_c - med) > opt$spike_m
    sf_spk <- exactextractr::exact_extract(spike, pd, "mean", progress = FALSE)
    ptab$dsm_spike_fraction <- sf_spk
    u <- sf::st_sf(geometry = sf::st_union(sf::st_geometry(pd)))
    tot <- exactextractr::exact_extract(spike, u, "mean", progress = FALSE)
    st <- if (is.na(tot)) "warn" else if (tot > TH$spike_fail) "fail" else if (tot > TH$spike_warn) "warn" else "pass"
    qc <- qc_add(qc, "dsm_spikes", st, sprintf("%.4f", tot), sprintf("warn > %.3f, fail > %.2f (|z - 5x5 median| > %g m)",
                                                                     TH$spike_warn, TH$spike_fail, opt$spike_m),
                 sprintf("%.2f%% of DSM pixels in plots are spikes or pits (matching errors, birds, dust, moving leaves).", 100 * tot))
    q <- exactextractr::exact_extract(dsm, u, "quantile", quantiles = c(0.01, 0.99), progress = FALSE)
    relief <- as.numeric(q[[2]] - q[[1]])
    qc <- qc_add(qc, "dsm_relief", if (relief > TH$relief_warn) "warn" else "pass", sprintf("%.2f m", relief),
                 sprintf("warn > %g m (1st-99th percentile over plots)", TH$relief_warn),
                 if (relief > TH$relief_warn) "Large elevation range over the plots: artifacts, a very steep field, or wrong units." else "Elevation range over plots plausible.")
  }
  chm <- NULL
  if (nzchar(opt$dtm)) {
    dtm <- read_raster(opt$dtm, "DTM")
    if (!nzchar(terra::crs(dtm))) {
      qc <- qc_add(qc, "canopy_height", "fail", "DTM has no CRS", "CRS required", "The DTM has no CRS.")
    } else {
      dtm2 <- if (terra::same.crs(dtm, dsm_c)) terra::resample(dtm, dsm_c, method = "bilinear") else
        terra::project(dtm, dsm_c, method = "bilinear")
      chm <- dsm_c - dtm2
      if (!is.null(pd)) {
        hs <- exactextractr::exact_extract(chm, pd, c("median", "quantile"), quantiles = 0.95, progress = FALSE)
        ptab$chm_median_m <- hs$median
        ptab$chm_p95_m <- hs$q95
        neg <- which(ptab$chm_median_m < TH$chm_neg)
        tall <- which(ptab$chm_p95_m > opt$max_height_m)
        st <- if (length(neg) || length(tall)) "warn" else "pass"
        msg <- c(if (length(neg)) sprintf("%d plot(s) have median canopy height below %.1f m: DSM and DTM are vertically offset or misaligned.", length(neg), TH$chm_neg),
                 if (length(tall)) sprintf("%d plot(s) have 95th-percentile height above %g m: spikes, objects, or the wrong DTM.", length(tall), opt$max_height_m))
        qc <- qc_add(qc, "canopy_height", st, sprintf("plot medians %.2f to %.2f m; p95 max %.2f m",
                                                      min(ptab$chm_median_m, na.rm = TRUE), max(ptab$chm_median_m, na.rm = TRUE),
                                                      max(ptab$chm_p95_m, na.rm = TRUE)),
                     sprintf("median >= %.1f m and p95 <= %g m", TH$chm_neg, opt$max_height_m),
                     if (length(msg)) paste(msg, collapse = " ") else "Canopy heights plausible.")
      }
    }
  }
  list(qc = qc, ptab = ptab, dsm = dsm_c, chm = chm)
}

run_qc <- function(opt) {
  need_pkgs(c("terra", "sf", "exactextractr"))
  if (!opt$expect %in% c("reflectance", "rgb", "dn", "any"))
    stop("--expect must be reflectance, rgb, dn or any.", call. = FALSE)
  if (nzchar(opt$dtm) && !nzchar(opt$dsm)) stop("--dtm needs --dsm.", call. = FALSE)
  qc <- data.frame(check = character(), status = character(), value = character(), threshold = character(),
                   message = character(), stringsAsFactors = FALSE)
  r <- read_raster(opt$input, "Raster")
  rot <- isTRUE(tryCatch(terra::is.rotated(r), error = function(e) FALSE))
  qc <- qc_add(qc, "rotation", if (rot) "fail" else "pass", rot, "no rotation terms",
               if (rot) "The geotransform is rotated: most GIS tools mis-place the pixels. Warp it to a north-up grid (gdalwarp)." else "North-up grid.")

  if (rot) {
    qc <- qc_add(qc, "other_checks", "info", "skipped", "-",
                 "Values of a rotated raster cannot be read reliably; warp it to a north-up grid and run the QC again.")
    return(list(qc = qc, band_stats = data.frame(), plot_table = NULL, raster = r, plots = NULL, alpha = NA_integer_,
                data_bands = seq_len(terra::nlyr(r)), sample = data.frame(), dsm = NULL, chm = NULL,
                gsd = abs(terra::res(r)), max_value = NA_real_, read_warnings = attr(r, "read_warnings")))
  }

  # CRS and pixel size
  cr <- terra::crs(r)
  if (!nzchar(cr)) {
    qc <- qc_add(qc, "crs", "fail", "none", "CRS defined", "No CRS: plots cannot be placed. Re-export with the CRS from the processing software.")
  } else if (is_ll(r)) {
    qc <- qc_add(qc, "crs", "warn", terra::crs(r, describe = TRUE)$name, "projected, meters",
                 "Geographic CRS (degrees): pixel size and buffers are not in meters. Export in the local UTM zone, or reproject once with nearest-neighbor resampling.")
  } else {
    lu <- terra::linearUnits(r)
    ok <- is.finite(lu) && abs(lu - 1) < 1e-9
    qc <- qc_add(qc, "crs", if (ok) "pass" else "warn", terra::crs(r, describe = TRUE)$name, "projected, meters",
                 if (ok) "Projected CRS in meters." else sprintf("Projected CRS but units are not meters (1 unit = %.6f m).", lu))
  }
  gs <- gsd_meters(r)
  shape <- abs(gs[1] - gs[2]) / max(gs)
  qc <- qc_add(qc, "pixel_shape", if (shape > TH$pix_shape) "warn" else "pass", sprintf("%.2f x %.2f cm", 100 * gs[1], 100 * gs[2]),
               sprintf("x and y within %.0f%%", 100 * TH$pix_shape),
               if (shape > TH$pix_shape) "Pixels are not square: the mosaic was resampled anisotropically." else "Square pixels.")
  if (!is.na(opt$max_gsd_cm))
    qc <- qc_add(qc, "gsd", if (100 * max(gs) > opt$max_gsd_cm) "warn" else "pass", sprintf("%.2f cm", 100 * max(gs)),
                 sprintf("<= %g cm (user)", opt$max_gsd_cm), "Ground sampling distance against the requested maximum.")

  # Bands and data type
  nb <- terra::nlyr(r)
  dts <- unique(terra::datatype(r))
  alpha <- detect_alpha(r, opt$alpha_band)
  data_bands <- setdiff(seq_len(nb), alpha)
  bc_status <- if (!is.na(opt$expect_bands) && nb != opt$expect_bands) "fail" else
    if (opt$expect == "rgb" && !nb %in% 3:4) "warn" else "pass"
  qc <- qc_add(qc, "band_count", bc_status, nb, if (!is.na(opt$expect_bands)) sprintf("%d (user)", as.integer(opt$expect_bands)) else if (opt$expect == "rgb") "3 or 4 for RGB(A)" else "-",
               sprintf("%d band(s)%s.", nb, if (!is.na(alpha)) sprintf("; band %d is an alpha mask", alpha) else ""))
  dt_status <- if (opt$expect == "reflectance" && any(dts == "INT1U")) "warn" else "pass"
  qc <- qc_add(qc, "datatype", dt_status, paste(dts, collapse = ","), if (opt$expect == "reflectance") "not 8-bit for reflectance" else "-",
               if (dt_status == "warn") "8-bit reflectance is a quantized display product, not for quantitative indices." else "Data type recorded.")
  nd <- declared_nodata(opt$input)
  int_max <- if (all(terra::datatype(r)[data_bands] %in% names(DTYPE_MAX))) DTYPE_MAX[[terra::datatype(r)[data_bands[1]]]] else NA_real_
  max_value <- if (!is.na(opt$max_value)) opt$max_value else int_max
  nd_bad <- length(nd) && !is.na(max_value) && any(nd == max_value, na.rm = TRUE)
  qc <- qc_add(qc, "nodata_declared", if (nd_bad) "warn" else "pass", if (length(nd)) paste(nd, collapse = ",") else "none",
               "NoData differs from the saturation level",
               if (nd_bad) sprintf("NoData = %s equals the saturation level: saturated (glare) pixels read as missing data.", format(max_value)) else
                 if (length(nd)) "Declared NoData value recorded." else "No NoData value declared.")

  # Whole-image sample: NoData fraction and zeros used as NoData
  whole <- as.data.frame(terra::spatSample(r[[data_bands]], min(opt$sample_size, terra::ncell(r)), method = "regular"))
  na_frac <- mean(!stats::complete.cases(whole))
  if (!is.na(alpha)) {
    av <- terra::spatSample(r[[alpha]], min(opt$sample_size, terra::ncell(r)), method = "regular")[[1]]
    na_frac <- mean(!stats::complete.cases(whole) | av == 0)
  }
  qc <- qc_add(qc, "nodata_fraction", "info", sprintf("%.3f", na_frac), "-",
               sprintf("%.1f%% of the image is NoData (collars around a mosaic are normal).", 100 * na_frac))
  zero_frac <- mean(stats::complete.cases(whole) & rowSums(whole == 0, na.rm = TRUE) == length(data_bands))
  zbad <- !length(nd) && is.na(alpha) && zero_frac > TH$zero_warn
  qc <- qc_add(qc, "zero_as_nodata", if (zbad) "warn" else "pass", sprintf("%.3f", zero_frac),
               sprintf("<= %.2f all-zero pixels without a NoData value", TH$zero_warn),
               if (zbad) "Many pixels are 0 in every band but no NoData value is declared: zeros are probably the mosaic collar and will be read as data. Declare NoData (gdal_translate -a_nodata 0)." else "No sign of undeclared NoData.")

  # Plots: extent, pixels across plots
  plots <- NULL; pm <- NULL; ids <- NULL; ptab <- NULL
  if (nzchar(opt$plots) && !nzchar(cr)) {
    qc <- qc_add(qc, "extent", "fail", "not checked", "raster CRS needed",
                 "The raster has no CRS, so the plots cannot be placed on it; plot-level checks were skipped.")
  } else if (nzchar(opt$plots)) {
    rp <- read_plots(opt$plots, opt$layer, opt$id_col)
    ids <- rp$ids
    pm <- sf::st_transform(rp$plots, metric_crs(rp$plots, r))
    plots <- sf::st_transform(rp$plots, sf::st_crs(cr))
    ext_poly <- sf::st_as_sfc(sf::st_bbox(c(xmin = terra::xmin(r), ymin = terra::ymin(r), xmax = terra::xmax(r),
                                            ymax = terra::ymax(r)), crs = sf::st_crs(cr)))
    inside <- vapply(seq_len(nrow(plots)), function(i) {
      a <- as.numeric(sf::st_area(sf::st_geometry(plots)[i]))
      x <- suppressWarnings(sf::st_intersection(sf::st_geometry(plots)[i], ext_poly))
      if (!length(x)) 0 else sum(as.numeric(sf::st_area(x))) / a
    }, numeric(1))
    st <- if (all(inside == 0)) "fail" else if (any(inside < 1e-6)) "fail" else if (any(inside < 1 - 1e-6)) "warn" else "pass"
    msg <- if (all(inside == 0)) "No plot overlaps the image: CRS mismatch or the wrong field?" else
      sprintf("%d plot(s) entirely outside and %d partly outside the image extent.", sum(inside < 1e-6), sum(inside >= 1e-6 & inside < 1 - 1e-6))
    qc <- qc_add(qc, "extent", st, sprintf("%d of %d plots fully inside", sum(inside >= 1 - 1e-6), length(inside)),
                 "all plots inside", msg)
    mrr <- sf::st_minimum_rotated_rectangle(sf::st_geometry(pm))
    side <- vapply(seq_along(mrr), function(i) {
      xy <- sf::st_coordinates(mrr[i])[1:3, 1:2]
      min(sqrt(sum((xy[2, ] - xy[1, ])^2)), sqrt(sum((xy[3, ] - xy[2, ])^2)))
    }, numeric(1))
    across <- side / max(gs)
    st <- if (min(across) < TH$px_across_fail) "fail" else if (min(across) < TH$px_across_warn) "warn" else "pass"
    qc <- qc_add(qc, "gsd_vs_plot", st, sprintf("min %.1f px across (narrowest plot %.2f m)", min(across), min(side)),
                 sprintf("warn < %d, fail < %d pixels across the narrowest plot side", TH$px_across_warn, TH$px_across_fail),
                 if (st == "pass") "Enough pixels across each plot." else "Pixels are coarse relative to the plots: values mix neighbors, alleys and soil (mixed pixels).")
    ptab <- data.frame(id = ids, inside_extent = inside, area_m2 = as.numeric(sf::st_area(pm)), stringsAsFactors = FALSE)
  }

  # Value ranges
  v <- if (!is.null(plots)) sample_pixels(r[[data_bands]], plots, opt$sample_size, opt$seed) else whole
  names(v) <- names(r)[data_bands]
  if (!is.null(plots) && !is.na(alpha)) {
    va <- sample_pixels(r[[c(data_bands, alpha)]], plots, opt$sample_size, opt$seed)
    v <- va[va[[ncol(va)]] > 0, seq_along(data_bands), drop = FALSE]
    names(v) <- names(r)[data_bands]
  }
  bs <- band_stats(v, if (opt$expect == "reflectance") opt$scale else 1)
  dn_like <- any(bs$p99 > TH$refl_p99_fail, na.rm = TRUE)
  if (opt$expect == "reflectance") {
    neg <- bs$band[bs$frac_negative > TH$refl_neg_warn]
    over <- bs$band[bs$frac_above_1 > TH$refl_over1_warn]
    st <- if (dn_like) "fail" else if (length(neg) || length(over)) "warn" else "pass"
    msg <- if (dn_like) "Values look like digital numbers (99th percentile above 1.5): not reflectance, or a scale factor is missing (--scale)." else
      paste(c(if (length(neg)) sprintf("More than 1%% negative values in %s (offset applied twice, or calibration error).", paste(neg, collapse = ", ")),
              if (length(over)) sprintf("More than 1%% of values above 1 in %s (glare or calibration error).", paste(over, collapse = ", "))), collapse = " ")
    qc <- qc_add(qc, "value_range", st, sprintf("p99 %s", paste(sprintf("%s %.3g", bs$band, bs$p99), collapse = ", ")),
                 "reflectance 0-1: fail if p99 > 1.5; warn if > 1% negative or > 1% above 1",
                 if (nzchar(msg)) msg else "Values look like surface reflectance.")
  } else {
    kind <- if (dn_like) "integer-scaled or digital numbers" else "0-1 (reflectance-like)"
    st <- if (opt$expect == "dn" && !dn_like) "warn" else if (opt$expect == "rgb" && !dn_like && !any(dts %in% c("FLT4S", "FLT8S"))) "warn" else "info"
    qc <- qc_add(qc, "value_range", st, sprintf("p99 %s", paste(sprintf("%s %.3g", bs$band, bs$p99), collapse = ", ")),
                 sprintf("expected: %s", opt$expect), sprintf("Values look %s.", kind))
  }

  # Plot coverage and saturation
  valid <- sum(is.na(r[[data_bands]])) == 0
  if (!is.na(alpha)) valid <- valid & (r[[alpha]] > 0)
  sat <- if (!is.na(max_value)) valid & (sum(r[[data_bands]] >= max_value) > 0) else NULL
  if (!is.null(plots)) {
    stk <- if (is.null(sat)) valid else c(valid, sat)
    names(stk) <- c("valid", if (!is.null(sat)) "sat")
    s <- exactextractr::exact_extract(stk, plots, "sum", progress = FALSE, force_df = TRUE, full_colnames = TRUE)
    ptab$coverage <- pmin(1, s$sum.valid / expected_pixels(plots, r))
    w <- which.min(ptab$coverage)
    st <- if (min(ptab$coverage) < TH$cov_fail) "fail" else if (min(ptab$coverage) < TH$cov_warn) "warn" else "pass"
    qc <- qc_add(qc, "plot_coverage", st, sprintf("min %.3f (%s)", ptab$coverage[w], ids[w]),
                 sprintf("warn < %.2f, fail < %.2f of plot area with data", TH$cov_warn, TH$cov_fail),
                 sprintf("%d plot(s) below %.2f coverage (outside the image, NoData holes, or masked by alpha).",
                         sum(ptab$coverage < TH$cov_warn), TH$cov_warn))
    if (!is.null(sat)) ptab$saturated_fraction <- ifelse(s$sum.valid > 0, s$sum.sat / s$sum.valid, NA)
  }
  if (is.null(sat)) {
    qc <- qc_add(qc, "saturation", "info", "not checked", "-", "Floating-point data: give --max-value to check saturation.")
  } else {
    sfrac <- if (!is.null(plots)) max(ptab$saturated_fraction, na.rm = TRUE) else
      mean(apply(as.matrix(whole) >= max_value, 1, any), na.rm = TRUE)
    st <- if (sfrac > TH$sat_fail) "fail" else if (sfrac > TH$sat_warn) "warn" else "pass"
    qc <- qc_add(qc, "saturation", st, sprintf("%.4f%s", sfrac, if (!is.null(plots)) " (worst plot)" else ""),
                 sprintf("warn > %.3f, fail > %.2f of pixels at %s", TH$sat_warn, TH$sat_fail, format(max_value)),
                 if (st == "pass") "Little or no saturation." else "Saturated pixels (glare, overexposure): their color and reflectance are clipped.")
  }

  dres <- NULL
  if (nzchar(opt$dsm)) {
    dres <- dsm_checks(qc, opt, r, plots, pm, ids, ptab)
    qc <- dres$qc
    ptab <- dres$ptab
  }
  if (!is.null(ptab)) names(ptab)[1] <- "plot_id"
  list(qc = qc, band_stats = bs, plot_table = ptab, raster = r, plots = plots, alpha = alpha,
       data_bands = data_bands, sample = v, dsm = if (!is.null(dres)) dres$dsm else NULL,
       chm = if (!is.null(dres)) dres$chm else NULL, gsd = gs, max_value = max_value,
       read_warnings = attr(r, "read_warnings"))
}

# ---------------------------------------------------------------------------
# Figures, report, main
# ---------------------------------------------------------------------------

status_col <- c(pass = "#2e7d32", warn = "#ef6c00", fail = "#c62828", info = "grey40")

# Layers to show as red, green, blue: by layer name when the names say so, else the first three.
composite_bands <- function(r, db) {
  nm <- gsub("[^a-z]", "", tolower(names(r)))
  hit <- match(c("red", "green", "blue"), nm)
  if (!anyNA(hit) && all(hit %in% db)) return(list(idx = hit, label = "red, green, blue layers"))
  list(idx = db[1:3], label = sprintf("layers %s shown as red, green, blue", paste(db[1:3], collapse = ", ")))
}

draw_figures <- function(res, out_dir) {
  r <- res$raster
  db <- res$data_bands
  open_png(out_dir, "overview.png", width = 7, height = 7)
  ok <- tryCatch({
    if (length(db) >= 3) {
      cb <- composite_bands(r, db)
      terra::plotRGB(r[[cb$idx]], stretch = "lin", mar = c(1, 1, 3, 1), axes = FALSE,
                     main = sprintf("Overview (%s, linear stretch)
plot outlines: yellow ok, orange/red low coverage", cb$label),
                     cex.main = 0.85)
    } else terra::plot(r[[db[1]]], main = "Overview (band 1)", col = grDevices::gray.colors(50))
    TRUE
  }, error = function(e) FALSE)
  if (!ok) graphics::plot.new()
  if (!is.null(res$plots) && ok) {
    cov <- res$plot_table$coverage
    col <- ifelse(cov < TH$cov_fail, status_col[["fail"]], ifelse(cov < TH$cov_warn, status_col[["warn"]], "yellow"))
    terra::plot(terra::vect(res$plots), add = TRUE, border = col, lwd = 2)
  }
  grDevices::dev.off()
  v <- res$sample
  if (NROW(v) && NCOL(v)) {
    nb <- ncol(v)
    open_png(out_dir, "band_histograms.png", width = 7, height = max(3, 1.8 * ceiling(nb / 2)))
    op <- graphics::par(mfrow = c(ceiling(nb / 2), min(2, nb)), mar = c(3.5, 3, 2, 1), mgp = c(2, 0.6, 0))
    for (j in seq_len(nb)) {
      x <- v[[j]][!is.na(v[[j]])]
      if (length(x)) graphics::hist(x, breaks = 60, col = "grey75", border = NA, main = names(v)[j], xlab = "value")
      else graphics::plot.new()
    }
    graphics::par(op)
    grDevices::dev.off()
  }
  if (!is.null(res$dsm)) {
    open_png(out_dir, "dsm.png", width = if (is.null(res$chm)) 7 else 11, height = 6)
    op <- graphics::par(mfrow = c(1, if (is.null(res$chm)) 1 else 2))
    terra::plot(res$dsm, main = "DSM (m)", col = grDevices::terrain.colors(50))
    if (!is.null(res$plots)) terra::plot(terra::vect(sf::st_transform(res$plots, sf::st_crs(terra::crs(res$dsm)))), add = TRUE)
    if (!is.null(res$chm)) {
      terra::plot(res$chm, main = "Canopy height, DSM - DTM (m)", col = grDevices::hcl.colors(50, "YlGn", rev = TRUE))
      if (!is.null(res$plots)) terra::plot(terra::vect(sf::st_transform(res$plots, sf::st_crs(terra::crs(res$chm)))), add = TRUE)
    }
    graphics::par(op)
    grDevices::dev.off()
  }
}

write_report <- function(res, opt, out_dir) {
  qc <- res$qc
  n <- table(factor(qc$status, levels = c("fail", "warn", "pass", "info")))
  ord <- order(match(qc$status, c("fail", "warn", "pass", "info")))
  bad <- qc[qc$status %in% c("fail", "warn"), ]
  r <- res$raster
  lines <- c(sprintf("# %s report", SKILL), "",
             sprintf("**QC result: %d fail, %d warn, %d pass, %d info.**", n[["fail"]], n[["warn"]], n[["pass"]], n[["info"]]), "")
  if (opt$demo) lines <- c(lines, "Demo data: the problems below were injected on purpose (see demo/README.md).", "")
  if (nrow(bad)) lines <- c(lines, "## Warnings", "", sprintf("- **%s** (%s): %s", bad$check[order(bad$status != "fail")],
                                                             bad$status[order(bad$status != "fail")], bad$message[order(bad$status != "fail")]), "")
  if (length(res$read_warnings)) lines <- c(lines, "GDAL/terra messages when reading:", paste("-", unique(res$read_warnings)), "")
  lines <- c(lines, "## Checks", "", "| check | status | value | threshold | note |", "|---|---|---|---|---|",
             sprintf("| %s | %s | %s | %s | %s |", qc$check[ord], qc$status[ord], qc$value[ord], qc$threshold[ord], qc$message[ord]), "",
             "## Image", "",
             sprintf("- File: `%s`; %d x %d pixels, %d band(s), %s.", opt$input, terra::ncol(r), terra::nrow(r), terra::nlyr(r),
                     paste(unique(terra::datatype(r)), collapse = ",")),
             sprintf("- Pixel size: %.2f x %.2f cm.", 100 * res$gsd[1], 100 * res$gsd[2]),
             sprintf("- Expected values: %s%s.", opt$expect, if (opt$expect == "reflectance" && opt$scale != 1) sprintf(" (after x %g)", opt$scale) else ""), "",
             "| band | min | 1% | median | 99% | max | negative | above 1 | zero |", "|---|---|---|---|---|---|---|---|---|",
             sprintf("| %s | %.4g | %.4g | %.4g | %.4g | %.4g | %.3f | %.3f | %.3f |", res$band_stats$band, res$band_stats$min,
                     res$band_stats$p01, res$band_stats$median, res$band_stats$p99, res$band_stats$max,
                     res$band_stats$frac_negative, res$band_stats$frac_above_1, res$band_stats$frac_zero), "",
             "![Overview](figures/overview.png)", "", "![Band histograms](figures/band_histograms.png)", "")
  if (!is.null(res$plot_table)) {
    pt <- res$plot_table
    worst <- utils::head(pt[order(pt$coverage), ], 10)
    lines <- c(lines, "## Plots", "", sprintf("%d plots; median coverage %.3f. Lowest coverage:", nrow(pt), stats::median(pt$coverage)), "",
               "| plot | inside extent | coverage | saturated |", "|---|---|---|---|",
               sprintf("| %s | %.3f | %.3f | %s |", worst$plot_id, worst$inside_extent, worst$coverage,
                       if (!is.null(worst$saturated_fraction)) sprintf("%.4f", worst$saturated_fraction) else "-"), "")
  }
  if (!is.null(res$dsm)) lines <- c(lines, "![DSM](figures/dsm.png)", "")
  lines <- c(lines, "## Outputs", "",
             "- `tables/qc_checks.csv`: one row per check: status, value, threshold, message.",
             "- `tables/band_stats.csv`: value distribution per band in the plots (or the whole image).",
             if (!is.null(res$plot_table)) "- `tables/plot_qc.csv`: per plot: inside_extent, coverage, saturated_fraction and, with a DSM, DSM coverage, spike fraction and canopy height." else NULL,
             "", "## Processing advice", "",
             "See `references/processing-guide.md` in this skill for OpenDroneMap settings, reflectance panels, ground control points, overlap and pixel size. That guide is advice, not something this script checks.",
             "", "## References", "",
             "- OpenDroneMap documentation. https://docs.opendronemap.org/",
             "- James MR, Chandler JH, Eltner A, et al. (2019) Guidelines on the use of structure-from-motion photogrammetry in geomorphic research. Earth Surface Processes and Landforms 44:2081-2084. https://doi.org/10.1002/esp.4637",
             "- Aasen H, Honkavaara E, Lucieer A, Zarco-Tejada PJ (2018) Quantitative remote sensing at ultra-high resolution with UAV spectroscopy. Remote Sensing 10:1091. https://doi.org/10.3390/rs10071091",
             "", sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(lines, file.path(out_dir, "report.md"))
}

demo_options <- function(opt) {
  src <- file.path(SKILL_DIR, "demo", "make_demo.R")
  if (!file.exists(src)) stop("Demo generator not found: ", src, call. = FALSE)
  e <- new.env()
  sys.source(src, envir = e)
  p <- e$make_demo(file.path(opt$out, "demo_input"))
  opt$input <- p$raster; opt$plots <- p$plots; opt$dsm <- p$dsm; opt$dtm <- p$dtm
  opt$expect <- "reflectance"; opt$scale <- 1e-4; opt$expect_bands <- 4
  opt
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  if (!opt$demo && !nzchar(opt$input)) stop("Provide --input RASTER, or --demo. See --help.", call. = FALSE)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  if (opt$demo) opt <- demo_options(opt)

  res <- run_qc(opt)
  save_table(res$qc, opt$out, "qc_checks.csv")
  save_table(res$band_stats, opt$out, "band_stats.csv")
  if (!is.null(res$plot_table)) save_table(res$plot_table, opt$out, "plot_qc.csv")
  draw_figures(res, opt$out)
  write_report(res, opt, opt$out)
  write_repro_bundle(opt$out, argv, inputs = c(opt$input, opt$plots, opt$dsm, opt$dtm))
  n <- table(factor(res$qc$status, levels = c("fail", "warn", "pass", "info")))
  cat(sprintf("QC: %d fail, %d warn, %d pass. Report: %s\n", n[["fail"]], n[["warn"]], n[["pass"]],
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  bad <- res$qc[res$qc$status %in% c("fail", "warn"), ]
  if (nrow(bad)) cat(sprintf("%s: %s - %s\n", toupper(bad$status), bad$check, bad$message), sep = "")
  if (opt$strict && n[["fail"]] > 0)
    stop(sprintf("QC failed (%s). See the report.", paste(res$qc$check[res$qc$status == "fail"], collapse = ", ")), call. = FALSE)
  invisible(res)
}

if (sys.nframe() == 0L) main()
