#!/usr/bin/env Rscript
# claw-canopy-cover: per-plot green canopy cover from an RGB orthomosaic with the
# excess-green index (ExG = 2g - r - b on chromatic coordinates, Woebbecke et al.
# 1995) and one Otsu threshold per flight, or ExG - ExR with a zero threshold
# (Meyer & Neto 2008). Saturated, dark and no-data pixels are counted and left
# out; QC images show the vegetation mask.
#
# Usage:
#   Rscript scripts/canopy_cover.R --demo
#   Rscript scripts/canopy_cover.R --input rgb_ortho.tif --plots plots.gpkg
#   Rscript scripts/canopy_cover.R --input rgb_ortho.tif --plots plots.gpkg --index exgr
#   Rscript scripts/canopy_cover.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-canopy-cover"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/canopy_cover.R [--demo] --input RGB_RASTER --plots POLYGONS [options]

  --demo            run on a synthetic RGB orthomosaic (written to <out>/demo_input)
  --input           RGB orthomosaic (GeoTIFF etc.), 8-bit, 16-bit or float
  --plots           plot polygons (GeoPackage, shapefile, GeoJSON), e.g. from claw-plot-grid
  --layer           layer name inside --plots (default: first layer)
  --id-col          unique plot id column (default: plot_id, else plot)
  --bands           band order (default 'red=1,green=2,blue=3')
  --alpha-band      alpha (transparency) band number; 0 = none; default: auto-detect band 4
  --index           exg (default) or exgr
  --threshold       auto (default: Otsu for exg, 0 for exgr), otsu, zero, or a number
  --max-value       saturation level (default: 255 for 8-bit, 65535 for 16-bit, 1 for float)
  --dark-threshold  pixels darker than this fraction of full brightness are left out (default 0.05)
  --buffer          buffer in metres applied to the plots; negative shrinks (default 0)
  --date            flight date YYYY-MM-DD, copied to the output table
  --qc-plots        number of plots shown in the QC image (default 6)
  --sample-size     pixels sampled for the histogram and Otsu threshold (default 200000)
  --seed            random seed for that sample (default 42)
  --out             output folder (default agriclaw_output/claw-canopy-cover)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, input = "", plots = "", layer = "", id_col = "",
                 bands = "red=1,green=2,blue=3", alpha_band = NA_real_, index = "exg",
                 threshold = "auto", max_value = NA_real_, dark_threshold = 0.05,
                 buffer = 0, date = "", qc_plots = 6, sample_size = 200000, seed = 42,
                 out = "")

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

need_pkgs <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop("Missing R package(s): ", paste(miss, collapse = ", "), ". Install with ",
         "install.packages(c(", paste0('"', miss, '"', collapse = ", "), ")).", call. = FALSE)
  invisible(TRUE)
}

# Chromatic coordinates r, g, b = R/(R+G+B) ... and the two colour indices.
# ExG = 2g - r - b (Woebbecke et al. 1995); ExR = 1.4r - g and
# ExGR = ExG - ExR = 3g - 2.4r - b (Meyer & Neto 2008).
color_indices <- function(R, G, B) {
  S <- R + G + B
  r <- R / S; g <- G / S; b <- B / S
  exg <- 2 * g - r - b
  exgr <- exg - (1.4 * r - g)
  bad <- !is.finite(exg) | S <= 0
  exg[bad] <- NA
  exgr[bad] <- NA
  list(exg = exg, exgr = exgr)
}

# Otsu (1979) threshold from a histogram: maximise the between-class variance
# sigma_B^2(k) = (mu_T w(k) - mu(k))^2 / (w(k) (1 - w(k))). Class 1 holds bins
# 1..k; the threshold lies between bin k and bin k + 1. Empty bins between the
# two modes give a plateau of equal maxima; the middle of the plateau is used.
# Attribute 'eta' = sigma_B^2 / sigma_T^2, Otsu's separability (1 = two spikes,
# 0.75 = uniform, 2/pi = 0.64 = one normal distribution).
otsu_from_hist <- function(counts, mids) {
  p <- counts / sum(counts)
  w <- cumsum(p)
  mu <- cumsum(p * mids)
  mu_t <- sum(p * mids)
  sb2 <- (mu_t * w - mu)^2 / (w * (1 - w))
  sb2[!is.finite(sb2)] <- NA
  sb2[length(sb2)] <- NA
  if (all(is.na(sb2))) return(NA_real_)
  k1 <- which.max(sb2)
  top <- !is.na(sb2) & sb2 >= sb2[k1] * (1 - 1e-12)
  k2 <- k1
  while (k2 < length(top) && top[k2 + 1L]) k2 <- k2 + 1L
  t1 <- (mids[k1] + mids[k1 + 1L]) / 2
  t2 <- (mids[k2] + mids[k2 + 1L]) / 2
  st2 <- sum(p * (mids - mu_t)^2)
  structure((t1 + t2) / 2, eta = sb2[k1] / st2)
}

otsu_threshold <- function(x, nbins = 256L) {
  x <- x[is.finite(x)]
  if (length(x) < 2L || diff(range(x)) == 0) return(NA_real_)
  br <- seq(min(x), max(x), length.out = nbins + 1L)
  h <- graphics::hist(x, breaks = br, plot = FALSE)
  otsu_from_hist(h$counts, h$mids)
}

parse_rgb_map <- function(s) {
  parts <- trimws(strsplit(s, ",", fixed = TRUE)[[1]])
  kv <- strsplit(parts[nzchar(parts)], "=", fixed = TRUE)
  if (length(kv) != 3L || any(lengths(kv) != 2L))
    stop("--bands must look like 'red=1,green=2,blue=3'; got '", s, "'.", call. = FALSE)
  m <- stats::setNames(suppressWarnings(as.integer(vapply(kv, `[`, "", 2))), tolower(trimws(vapply(kv, `[`, "", 1))))
  if (!setequal(names(m), c("red", "green", "blue")) || anyNA(m) || anyDuplicated(m))
    stop("--bands needs red, green and blue, each a different layer number; got '", s, "'.", call. = FALSE)
  m[c("red", "green", "blue")]
}

# Plot polygons: checks, id column, buffer in metres, then the raster CRS.
prepare_plots <- function(path, layer, id_col, r, buffer_m) {
  if (!file.exists(path)) stop("Plot file not found: ", path, call. = FALSE)
  p <- if (nzchar(layer)) sf::st_read(path, layer = layer, quiet = TRUE) else sf::st_read(path, quiet = TRUE)
  if (is.na(sf::st_crs(p))) stop("The plot polygons have no CRS; assign the correct one first.", call. = FALSE)
  gt <- unique(as.character(sf::st_geometry_type(p)))
  if (!all(gt %in% c("POLYGON", "MULTIPOLYGON")))
    stop("--plots must contain polygons (found ", paste(gt, collapse = ", "), ").", call. = FALSE)
  if (!nzchar(id_col)) id_col <- if ("plot_id" %in% names(p)) "plot_id" else if ("plot" %in% names(p)) "plot" else ""
  if (!nzchar(id_col) || !id_col %in% names(p))
    stop("Plot id column not found. Give --id-col; columns: ",
         paste(setdiff(names(p), attr(p, "sf_column")), collapse = ", "), ".", call. = FALSE)
  ids <- p[[id_col]]
  if (anyNA(ids) || anyDuplicated(ids))
    stop("Plot id column '", id_col, "' has missing or duplicated values.", call. = FALSE)
  rcrs <- sf::st_crs(terra::crs(r))
  work <- if (!isTRUE(sf::st_is_longlat(p))) sf::st_crs(p) else if (!terra::is.lonlat(r)) rcrs else {
    cc <- sf::st_coordinates(sf::st_transform(sf::st_centroid(sf::st_as_sfc(sf::st_bbox(p))), 4326))
    zone <- floor((cc[1] + 180) / 6) + 1
    sf::st_crs((if (cc[2] >= 0) 32600 else 32700) + zone)
  }
  pm <- sf::st_transform(p, work)
  if (buffer_m != 0) {
    pm <- sf::st_buffer(pm, dist = buffer_m, joinStyle = "MITRE", mitreLimit = 5)
    gone <- sf::st_is_empty(pm)
    if (any(gone))
      stop(sprintf("--buffer %.3f m removes %d plot(s) entirely (e.g. %s). Use a smaller inward buffer.",
                   buffer_m, sum(gone), paste(utils::head(ids[gone], 5), collapse = ", ")), call. = FALSE)
  }
  list(plots = sf::st_transform(pm, rcrs), id_col = id_col, area_m2 = as.numeric(sf::st_area(pm)))
}

pixel_area_m2 <- function(r) {
  rs <- terra::res(r)
  if (terra::is.lonlat(r)) {
    lat <- mean(as.vector(terra::ext(r))[3:4]) * pi / 180
    return(rs[1] * 111320 * cos(lat) * rs[2] * 110574)
  }
  lu <- terra::linearUnits(r)
  if (!is.finite(lu) || lu == 0) lu <- 1
  prod(rs) * lu^2
}

# NoData values declared in the file (GDAL metadata), if any.
declared_nodata <- function(path) {
  d <- tryCatch(terra::describe(path), error = function(e) character())
  v <- regmatches(d, regexpr("NoData Value=[^ ]+", d))
  v <- unique(suppressWarnings(as.numeric(sub("NoData Value=", "", v))))
  v[is.finite(v)]
}

# Saturation level from the data type unless given.
saturation_level <- function(r, given, observed_max) {
  if (!is.na(given)) return(given)
  dt <- terra::datatype(r)[1]
  if (identical(dt, "INT1U")) return(255)
  if (identical(dt, "INT2U")) return(65535)
  if (is.finite(observed_max) && observed_max <= 1) return(1)
  NA_real_
}

# Alpha band: explicit number, 0 = none, NA = auto (band 4 holding only 0 and the maximum).
detect_alpha <- function(r, given, rgb) {
  if (!is.na(given)) {
    if (given == 0) return(NA_integer_)
    if (given < 1 || given > terra::nlyr(r) || given %in% rgb)
      stop("--alpha-band must be a layer number other than the RGB bands (or 0 for none).", call. = FALSE)
    return(as.integer(given))
  }
  if (terra::nlyr(r) < 4L || 4L %in% rgb) return(NA_integer_)
  v <- unique(stats::na.omit(terra::spatSample(r[[4]], min(10000, terra::ncell(r)), method = "regular")[[1]]))
  mx <- if (identical(terra::datatype(r)[4], "INT2U")) 65535 else 255
  if (length(v) && all(v %in% c(0, mx))) 4L else NA_integer_
}

# Per-pixel layers: indices (NA if unusable) and 0/1 indicators for counting.
# `bright_ref` is full brightness (the saturation level, or the image maximum).
pixel_layers <- function(rgb, alpha, max_value, dark_threshold, bright_ref = max_value) {
  f <- function(R, G, B, A) {
    valid <- !is.na(R) & !is.na(G) & !is.na(B) & (R + G + B) > 0
    if (!is.null(A)) valid <- valid & !is.na(A) & A > 0
    sat <- valid & (if (is.na(max_value)) FALSE else (R >= max_value | G >= max_value | B >= max_value))
    dark <- valid & !sat & ((R + G + B) / (3 * bright_ref)) < dark_threshold
    use <- valid & !sat & !dark
    ci <- color_indices(R, G, B)
    exg <- ifelse(use, ci$exg, NA)
    exgr <- ifelse(use, ci$exgr, NA)
    cbind(exg = exg, exgr = exgr, valid = as.numeric(valid), sat = as.numeric(sat),
          dark = as.numeric(dark), use = as.numeric(use))
  }
  if (is.null(alpha)) terra::lapp(rgb, function(R, G, B) f(R, G, B, NULL))
  else terra::lapp(c(rgb, alpha), function(R, G, B, A) f(R, G, B, A))
}

sample_values <- function(x, plots, n_total, seed) {
  per <- max(50L, ceiling(n_total / nrow(plots)))
  set.seed(seed)
  v <- exactextractr::exact_extract(x, plots, function(values, coverage_fraction) {
    ok <- which(!is.na(values) & coverage_fraction >= 0.5)
    if (length(ok) > per) ok <- sample(ok, per)
    data.frame(v = values[ok])
  }, progress = FALSE)
  if (is.data.frame(v)) v$v else unlist(v)
}

resolve_threshold <- function(spec, index, sample) {
  spec <- tolower(trimws(spec))
  if (spec == "auto") spec <- if (index == "exg") "otsu" else "zero"
  if (spec == "otsu") {
    t <- otsu_threshold(sample)
    if (is.na(t)) stop("Otsu threshold could not be computed (no usable pixels or no spread).", call. = FALSE)
    return(list(value = as.numeric(t), method = "Otsu", eta = attr(t, "eta")))
  }
  if (spec == "zero") return(list(value = 0, method = "fixed 0", eta = NA_real_))
  v <- suppressWarnings(as.numeric(spec))
  if (is.na(v)) stop("--threshold must be auto, otsu, zero or a number; got '", spec, "'.", call. = FALSE)
  list(value = v, method = "user value", eta = NA_real_)
}

run_cover <- function(opt) {
  need_pkgs(c("terra", "sf", "exactextractr"))
  if (!file.exists(opt$input)) stop("Raster not found: ", opt$input, call. = FALSE)
  if (!opt$index %in% c("exg", "exgr")) stop("--index must be exg or exgr.", call. = FALSE)
  if (opt$dark_threshold < 0 || opt$dark_threshold >= 1) stop("--dark-threshold must be in [0, 1).", call. = FALSE)
  if (nzchar(opt$date) && is.na(as.Date(opt$date, format = "%Y-%m-%d")))
    stop("--date must be YYYY-MM-DD, got '", opt$date, "'.", call. = FALSE)
  warn <- character()
  notes <- character()
  r0 <- terra::rast(opt$input)
  if (!nzchar(terra::crs(r0))) stop("The raster has no CRS; it cannot be matched to the plots.", call. = FALSE)
  if (terra::nlyr(r0) < 3L) stop("An RGB raster needs at least 3 layers (found ", terra::nlyr(r0), ").", call. = FALSE)
  rgb_idx <- parse_rgb_map(opt$bands)
  if (any(rgb_idx > terra::nlyr(r0))) stop("--bands refers to a layer the raster does not have.", call. = FALSE)
  alpha_idx <- detect_alpha(r0, opt$alpha_band, rgb_idx)
  if (!is.na(alpha_idx) && is.na(opt$alpha_band)) notes <- c(notes, "Band 4 holds only 0 and the maximum value: used as an alpha (transparency) mask.")
  if (terra::nlyr(r0) > 3L && is.na(alpha_idx))
    notes <- c(notes, sprintf("The raster has %d layers; only %s are used.", terra::nlyr(r0),
                              paste(sprintf("%s = %d", names(rgb_idx), rgb_idx), collapse = ", ")))
  if (identical(opt$bands, DEFAULTS$bands))
    notes <- c(notes, "Band order assumed red, green, blue (layers 1, 2, 3), the usual order of RGB orthomosaics.")
  pp <- prepare_plots(opt$plots, opt$layer, opt$id_col, r0, opt$buffer)
  plots <- pp$plots
  ids <- plots[[pp$id_col]]
  pv <- terra::vect(plots)
  e1 <- as.vector(terra::ext(r0)); e2 <- as.vector(terra::ext(pv))
  if (!(e2[1] < e1[2] && e2[2] > e1[1] && e2[3] < e1[4] && e2[4] > e1[3]))
    stop("The plots do not overlap the raster. Check that both are georeferenced to the same place.", call. = FALSE)
  keep_layers <- c(rgb_idx, if (!is.na(alpha_idx)) alpha_idx)
  r <- terra::crop(r0[[keep_layers]], terra::ext(pv) + 3 * max(terra::res(r0)))
  rgb <- r[[1:3]]
  names(rgb) <- c("red", "green", "blue")
  alpha <- if (!is.na(alpha_idx)) r[[4]] else NULL

  obs_max <- max(terra::global(rgb, "max", na.rm = TRUE)$max)
  max_value <- saturation_level(r0, opt$max_value, obs_max)
  if (is.na(max_value))
    warn <- c(warn, "Saturation level unknown for this data type: saturated pixels were not detected. Give --max-value.")
  nd <- declared_nodata(opt$input)
  if (!is.na(max_value) && any(nd == max_value, na.rm = TRUE))
    warn <- c(warn, sprintf(paste("The file declares NoData = %s, the saturation level: saturated pixels read as no data,",
                                  "so glare shows up as missing coverage (valid_fraction) instead of saturation."), format(max_value)))
  if (!is.na(max_value) && identical(terra::datatype(r0)[1], "INT2U") && obs_max < 4096)
    notes <- c(notes, "16-bit data use at most 12 bits here; saturation is judged at 65535 unless --max-value is given.")

  px <- pixel_layers(rgb, alpha, max_value, opt$dark_threshold,
                     bright_ref = if (is.na(max_value)) obs_max else max_value)
  smp <- sample_values(px[[opt$index]], plots, opt$sample_size, opt$seed)
  if (length(smp) < 100L) stop("Fewer than 100 usable pixels inside the plots.", call. = FALSE)
  th <- resolve_threshold(opt$threshold, opt$index, smp)
  if (th$method == "Otsu" && is.finite(th$eta) && th$eta < 0.8)
    warn <- c(warn, sprintf(paste("The %s histogram is not clearly bimodal (Otsu separability %.2f; 0.64 for a single",
                                  "normal mode): the image may be nearly all canopy or all soil, and Otsu splits it",
                                  "anyway. Check figures/index_histogram.png; consider --index exgr."), toupper(opt$index), th$eta))
  other <- if (opt$index == "exg") "exgr" else "exg"
  smp_o <- sample_values(px[[other]], plots, opt$sample_size, opt$seed)
  th_o <- if (other == "exgr") list(value = 0, method = "fixed 0") else {
    t <- otsu_threshold(smp_o); list(value = as.numeric(t), method = "Otsu", eta = attr(t, "eta"))
  }

  veg <- terra::ifel(px[[opt$index]] > th$value, 1, 0)
  veg2 <- terra::ifel(px[[other]] > th_o$value, 1, 0)
  names(veg) <- "veg"; names(veg2) <- "veg2"
  cov <- exactextractr::exact_extract(c(veg, veg2), plots, "mean", progress = FALSE, force_df = TRUE,
                                      full_colnames = TRUE)
  cnt <- exactextractr::exact_extract(px[[c("valid", "sat", "dark", "use")]], plots, "sum", progress = FALSE,
                                      force_df = TRUE, full_colnames = TRUE)
  pa <- pixel_area_m2(r)
  expected <- pp$area_m2 / pa
  res <- data.frame(id = ids, stringsAsFactors = FALSE)
  res$canopy_cover_pct <- 100 * cov$mean.veg
  res[[sprintf("cover_%s_pct", other)]] <- 100 * cov$mean.veg2
  res$n_pixels <- cnt$sum.use
  res$valid_fraction <- pmin(1, cnt$sum.valid / expected)
  res$saturated_fraction <- ifelse(cnt$sum.valid > 0, cnt$sum.sat / cnt$sum.valid, NA)
  res$dark_fraction <- ifelse(cnt$sum.valid > 0, cnt$sum.dark / cnt$sum.valid, NA)
  res$plot_area_m2 <- pp$area_m2
  flag <- rep("", nrow(res))
  addf <- function(cond, txt) ifelse(!is.na(cond) & cond, ifelse(nzchar(flag), paste(flag, txt, sep = "; "), txt), flag)
  flag <- addf(res$n_pixels < 1, "no usable pixels")
  flag <- addf(res$valid_fraction < 0.95, "incomplete image coverage")
  flag <- addf(res$saturated_fraction > 0.01, "saturated pixels")
  flag <- addf(res$dark_fraction > 0.10, "dark/shadow pixels")
  res$flag <- flag

  diff_pp <- abs(res$canopy_cover_pct - res[[sprintf("cover_%s_pct", other)]])
  if (stats::median(diff_pp, na.rm = TRUE) > 10)
    warn <- c(warn, sprintf(paste("ExG (Otsu) and ExGR (> 0) disagree by a median of %.1f percentage points: residue,",
                                  "colour cast or a poor threshold. Inspect figures/qc_plots.png before using either."),
                            stats::median(diff_pp, na.rm = TRUE)))
  nf <- sum(nzchar(flag))
  if (nf) warn <- c(warn, sprintf("%d of %d plot(s) flagged (%s). See the flag column.", nf, nrow(res),
                                  paste(names(table(unlist(strsplit(flag[nzchar(flag)], "; ")))), collapse = ", ")))
  if (stats::median(expected) < 100)
    warn <- c(warn, sprintf("A median plot has only %.0f pixels: too coarse for canopy cover (rule of thumb: >= 100, ideally 1-2 cm pixels).",
                            stats::median(expected)))

  attrs <- sf::st_drop_geometry(plots)
  extra <- attrs[, setdiff(names(attrs), pp$id_col), drop = FALSE]
  extra <- extra[, !names(extra) %in% names(res), drop = FALSE]
  res <- cbind(res[1], extra, res[-1])
  names(res)[1] <- pp$id_col
  if (nzchar(opt$date)) res <- cbind(res[1], date = opt$date, res[-1])
  list(table = res, threshold = th, other = other, threshold_other = th_o, max_value = max_value,
       sample = smp, sample_other = smp_o, warnings = warn, notes = notes, plots = plots, rgb = rgb,
       px = px, veg = veg, id_col = pp$id_col, pixel_area = pa, raster = r0, alpha_idx = alpha_idx,
       rgb_idx = rgb_idx, diff_pp = diff_pp)
}

# ---------------------------------------------------------------------------
# Figures, report, main
# ---------------------------------------------------------------------------

stretch_rgb <- function(rgb, max_value) {
  mv <- if (is.na(max_value)) max(terra::global(rgb, "max", na.rm = TRUE)$max) else max_value
  terra::clamp(rgb / mv * 255, 0, 255)
}

draw_figures <- function(res, opt, out_dir) {
  idx <- toupper(opt$index)
  open_png(out_dir, "index_histogram.png", width = 7, height = 4.5)
  graphics::hist(res$sample, breaks = 120, col = "grey80", border = NA, xlab = idx,
                 main = sprintf("%s of usable plot pixels (sample)", idx))
  graphics::abline(v = res$threshold$value, col = "#cc0000", lwd = 2)
  graphics::legend("topleft", bty = "n", col = "#cc0000", lwd = 2,
                   legend = sprintf("threshold %.3f (%s, one per flight)", res$threshold$value, res$threshold$method))
  grDevices::dev.off()

  rgb8 <- stretch_rgb(res$rgb, res$max_value)
  mask <- terra::ifel(res$veg == 1, 1, NA)
  open_png(out_dir, "cover_map.png", width = 7, height = 7)
  terra::plotRGB(rgb8, r = 1, g = 2, b = 3, scale = 255, mar = c(1, 1, 3, 1), axes = FALSE,
                 main = sprintf("Vegetation mask (%s > %.3f) over the orthomosaic", idx, res$threshold$value))
  terra::plot(mask, col = "#ff00ff", alpha = 0.45, legend = FALSE, add = TRUE)
  terra::plot(terra::vect(res$plots), add = TRUE, border = "yellow", lwd = 1.2)
  cc <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(res$plots)))
  graphics::text(cc[, 1], cc[, 2], labels = sprintf("%.0f%%", res$table$canopy_cover_pct), col = "white", cex = 0.7, font = 2)
  grDevices::dev.off()

  tab <- res$table
  o <- order(tab$canopy_cover_pct)
  n <- min(opt$qc_plots, nrow(tab))
  flagged <- utils::head(which(nzchar(tab$flag)), floor(n / 2))
  spread <- o[round(seq(1, length(o), length.out = n))]
  pick <- utils::head(unique(c(flagged, spread)), n)
  pick <- pick[order(tab$canopy_cover_pct[pick])]
  open_png(out_dir, "qc_plots.png", width = 5.5, height = 2.2 * length(pick) + 0.5)
  op <- graphics::par(mfrow = c(length(pick), 2), mar = c(0.5, 0.5, 2, 0.5))
  for (i in pick) {
    g <- terra::vect(res$plots[i, ])
    e <- terra::ext(g) + 2 * max(terra::res(res$rgb))
    sub <- terra::crop(rgb8, e)
    terra::plotRGB(sub, r = 1, g = 2, b = 3, scale = 255, mar = c(0.5, 0.5, 2, 0.5), axes = FALSE,
                   main = sprintf("%s: RGB", tab[[res$id_col]][i]))
    terra::plot(g, add = TRUE, border = "yellow", lwd = 1.5)
    terra::plotRGB(sub, r = 1, g = 2, b = 3, scale = 255, mar = c(0.5, 0.5, 2, 0.5), axes = FALSE,
                   main = sprintf("mask: %.1f%% cover%s", tab$canopy_cover_pct[i], if (nzchar(tab$flag[i])) paste0(" (", tab$flag[i], ")") else ""))
    terra::plot(terra::crop(mask, e), col = "#ff00ff", alpha = 0.6, legend = FALSE, add = TRUE)
    terra::plot(g, add = TRUE, border = "yellow", lwd = 1.5)
  }
  graphics::par(op)
  grDevices::dev.off()
}

fmt <- function(x, d = 1) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))

write_report <- function(res, opt, out_dir) {
  t <- res$table
  r <- res$raster
  oc <- sprintf("cover_%s_pct", res$other)
  lines <- c(sprintf("# %s report", SKILL), "")
  if (length(res$warnings)) lines <- c(lines, "## Warnings", "", paste("-", res$warnings), "")
  lines <- c(lines, "## Inputs", "",
             sprintf("- Raster: `%s` (%d layers, %s, %s, pixel %.4f m2 = %.2f cm on a side).", opt$input, terra::nlyr(r),
                     terra::datatype(r)[1], terra::crs(r, describe = TRUE)$name, res$pixel_area, 100 * sqrt(res$pixel_area)),
             sprintf("- Bands: %s%s.", paste(sprintf("%s = layer %d", names(res$rgb_idx), res$rgb_idx), collapse = ", "),
                     if (!is.na(res$alpha_idx)) sprintf(", alpha = layer %d", res$alpha_idx) else ""),
             sprintf("- Plots: `%s`, %d polygons, id column `%s`%s.", opt$plots, nrow(t), res$id_col,
                     if (opt$buffer != 0) sprintf(", buffered by %.3f m", opt$buffer) else ""),
             sprintf("- Saturation level %s; pixels darker than %.0f%% of full brightness left out.",
                     if (is.na(res$max_value)) "unknown" else format(res$max_value), 100 * opt$dark_threshold),
             if (nzchar(opt$date)) sprintf("- Date: %s.", opt$date) else NULL,
             if (length(res$notes)) paste("-", res$notes) else NULL, "",
             "## Method", "",
             sprintf("- Primary: %s with threshold %.4f (%s%s), one threshold for the whole flight.", toupper(opt$index),
                     res$threshold$value, res$threshold$method,
                     if (is.finite(res$threshold$eta)) sprintf(", separability %.2f", res$threshold$eta) else ""),
             sprintf("- Cross-check: %s with threshold %.4f (%s): median absolute difference %.1f percentage points.",
                     toupper(res$other), res$threshold_other$value, res$threshold_other$method,
                     stats::median(res$diff_pp, na.rm = TRUE)), "",
             "## Results", "",
             sprintf("Canopy cover over %d plots: median %.1f%% (range %.1f-%.1f%%).", nrow(t),
                     stats::median(t$canopy_cover_pct, na.rm = TRUE), min(t$canopy_cover_pct, na.rm = TRUE),
                     max(t$canopy_cover_pct, na.rm = TRUE)), "",
             "| plot | cover (%) | cross-check (%) | usable pixels | saturated | dark | flag |", "|---|---|---|---|---|---|---|",
             sprintf("| %s | %s | %s | %.0f | %s | %s | %s |", t[[res$id_col]], fmt(t$canopy_cover_pct), fmt(t[[oc]]),
                     t$n_pixels, fmt(100 * t$saturated_fraction, 2), fmt(100 * t$dark_fraction, 2), t$flag)[
                       seq_len(min(nrow(t), 60))],
             if (nrow(t) > 60) sprintf("| ... | %d more plots in tables/canopy_cover.csv | | | | | |", nrow(t) - 60) else NULL,
             "", "![Histogram](figures/index_histogram.png)", "", "![Cover map](figures/cover_map.png)", "",
             "![QC plots](figures/qc_plots.png)", "",
             "## Outputs", "",
             "- `tables/canopy_cover.csv`: one row per plot: `canopy_cover_pct` (primary method), the cross-check column, `n_pixels` (usable, coverage-weighted), `valid_fraction`, `saturated_fraction`, `dark_fraction`, `plot_area_m2`, `flag`.",
             "- `spatial/veg_mask.tif`: 1 = vegetation, 0 = background, 2 = left out (saturated or dark), no data elsewhere.",
             "- `figures/`: index histogram with the threshold, cover map, QC panels for selected plots.",
             "", "## Limitations", "",
             "- Green canopy only: senescent, chlorotic or purple leaves count as background; residue with a green tint, algae and weeds count as canopy.",
             "- Cover is the fraction of usable pixels; saturated and dark pixels are left out, so heavy shadow or glare biases it.",
             "- Mixed pixels at leaf edges are classified whole; coarse pixels bias cover towards the dominant class.",
             "", "## References", "",
             "- Woebbecke DM, Meyer GE, Von Bargen K, Mortensen DA (1995) Color indices for weed identification under various soil, residue, and lighting conditions. Transactions of the ASAE 38:259-269. https://doi.org/10.13031/2013.27838",
             "- Meyer GE, Neto JC (2008) Verification of color vegetation indices for automated crop imaging applications. Computers and Electronics in Agriculture 63:282-293. https://doi.org/10.1016/j.compag.2008.03.009",
             "- Otsu N (1979) A threshold selection method from gray-level histograms. IEEE Transactions on Systems, Man, and Cybernetics 9:62-66. https://doi.org/10.1109/TSMC.1979.4310076",
             "", sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(lines, file.path(out_dir, "report.md"))
}

demo_options <- function(opt) {
  src <- file.path(SKILL_DIR, "demo", "make_demo.R")
  if (!file.exists(src)) stop("Demo generator not found: ", src, call. = FALSE)
  e <- new.env()
  sys.source(src, envir = e)
  paths <- e$make_demo(file.path(opt$out, "demo_input"))
  opt$input <- paths$raster
  opt$plots <- paths$plots
  opt$date <- "2026-06-20"
  opt
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  if (!opt$demo && (!nzchar(opt$input) || !nzchar(opt$plots)))
    stop("Provide --input RGB_RASTER and --plots POLYGONS, or --demo. See --help.", call. = FALSE)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  if (opt$demo) opt <- demo_options(opt)

  res <- run_cover(opt)
  save_table(res$table, opt$out, "canopy_cover.csv")
  dir.create(file.path(opt$out, "spatial"), showWarnings = FALSE)
  code <- terra::ifel(res$px[["valid"]] == 1, terra::ifel(res$px[["use"]] == 1, res$veg, 2), NA)
  terra::writeRaster(code, file.path(opt$out, "spatial", "veg_mask.tif"), datatype = "INT1U", overwrite = TRUE,
                     gdal = c("COMPRESS=DEFLATE"))
  draw_figures(res, opt, opt$out)
  write_report(res, opt, opt$out)
  write_repro_bundle(opt$out, argv, inputs = c(opt$input, opt$plots))
  cat(sprintf("Done: %d plots, median cover %.1f%%. Report: %s\n", nrow(res$table),
              stats::median(res$table$canopy_cover_pct, na.rm = TRUE),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  if (length(res$warnings)) cat(paste0("Warning: ", res$warnings, "\n"), sep = "")
  invisible(res)
}

if (sys.nframe() == 0L) main()
