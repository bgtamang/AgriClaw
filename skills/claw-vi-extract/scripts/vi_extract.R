#!/usr/bin/env Rscript
# claw-vi-extract: per-plot vegetation indices (NDVI, NDRE, GNDVI, OSAVI, EVI2,
# MCARI2, CIred-edge) from a multispectral orthomosaic or satellite image, with
# coverage-weighted plot statistics from exactextractr, optional soil/shadow
# masking, negative buffers and mixed-pixel warnings. Refuses digital numbers.
#
# Usage:
#   Rscript scripts/vi_extract.R --demo
#   Rscript scripts/vi_extract.R --input ortho.tif --plots plots.gpkg \
#       --bands "blue=1,green=2,red=3,rededge=4,nir=5" --date 2026-07-15
#   Rscript scripts/vi_extract.R --input s2.tif --plots plots.gpkg --sensor sentinel2-l2a
#   Rscript scripts/vi_extract.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT license.

SKILL   <- "claw-vi-extract"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/vi_extract.R [--demo] --input RASTER --plots POLYGONS
                                [--bands MAP | --sensor PRESET] [options]

  --demo        run on synthetic demo data (written to <out>/demo_input)
  --input       multispectral raster (GeoTIFF etc.), surface reflectance
  --plots       plot polygons (GeoPackage, shapefile, GeoJSON), e.g. from claw-plot-grid
  --layer       layer name inside --plots (default: first layer)
  --id-col      unique plot id column (default: plot_id, else plot)
  --bands       band roles to raster layers, e.g. 'blue=1,green=2,red=3,rededge=4,nir=5'
                (numbers or layer names); confirm the order with the user
  --sensor      preset (suggestion only): micasense-rededge-mx, micasense-altum, dji-p4m,
                sentinel2-l2a, sentinel2-l2a-harmonized, planetscope-4band, planetscope-8band
  --scale, --offset   reflectance = value x scale + offset (default 1 and 0, or the preset)
  --nodata      raw value meaning no data (default: the file's own; 0 for satellite presets)
  --indices     comma list or 'auto' (all possible): ndvi,ndre,gndvi,osavi,evi2,mcari2,cire
  --mask        none | ndvi | otsu: drop soil pixels below an NDVI threshold (default none)
  --mask-threshold  NDVI threshold for --mask ndvi
  --shadow-nir  also drop pixels with NIR reflectance below this value (default: off)
  --buffer      buffer in meters applied to plots first; negative shrinks (default 0)
  --min-pixels  flag plots with fewer effective pixels than this (default 10)
  --date        flight or acquisition date YYYY-MM-DD, copied to the output tables
  --confirm-calibrated  proceed although values look like digital numbers
  --sample-size pixels sampled for histograms and the Otsu threshold (default 200000)
  --seed        random seed for that sample (default 42)
  --out         output folder (default agriclaw_output/claw-vi-extract)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, input = "", plots = "", layer = "", id_col = "",
                 bands = "", sensor = "", scale = NA_real_, offset = NA_real_,
                 nodata = NA_real_, indices = "auto", mask = "none",
                 mask_threshold = NA_real_, shadow_nir = NA_real_, buffer = 0,
                 min_pixels = 10, date = "", confirm_calibrated = FALSE,
                 sample_size = 200000, seed = 42, out = "")

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

ROLES <- c("blue", "green", "red", "rededge", "nir")

# Sensor presets: a suggested band order and reflectance scaling. Band order in
# stacked orthomosaics depends on the processing software: always confirm.
SENSORS <- list(
  "micasense-rededge-mx" = list(
    bands = c(blue = "1", green = "2", red = "3", rededge = "4", nir = "5"), scale = 1, offset = 0, nodata = NA,
    note = paste("MicaSense RedEdge-MX camera bands are 1 blue (475 nm), 2 green (560), 3 red (668), 4 NIR (842),",
                 "5 red edge (717). Pix4D and Metashape write stacked reflectance in wavelength order",
                 "(blue, green, red, red edge, NIR), assumed here; the MicaSense image-processing library keeps",
                 "camera order (red edge last).")),
  "micasense-altum" = list(
    bands = c(blue = "1", green = "2", red = "3", rededge = "4", nir = "5"), scale = 1, offset = 0, nodata = NA,
    note = paste("MicaSense Altum camera bands are 1 blue, 2 green, 3 red, 4 NIR, 5 red edge, 6 thermal.",
                 "Wavelength order (blue, green, red, red edge, NIR) is assumed here; check where the thermal",
                 "band and red edge ended up in the stack.")),
  "dji-p4m" = list(
    bands = c(blue = "1", green = "2", red = "3", rededge = "4", nir = "5"), scale = 1, offset = 0, nodata = NA,
    note = "DJI P4 Multispectral bands: blue 450, green 560, red 650, red edge 730, NIR 840 nm (wavelength order assumed)."),
  "sentinel2-l2a" = list(
    bands = c(blue = "B02", green = "B03", red = "B04", rededge = "B05", nir = "B08"), scale = 1e-4, offset = -0.1, nodata = 0,
    note = paste("ESA Sentinel-2 L2A from processing baseline 04.00 (January 25, 2022) on: reflectance =",
                 "(DN + BOA_ADD_OFFSET) / QUANTIFICATION_VALUE = (DN - 1000) / 10000; DN 0 = no data.",
                 "Layers are matched by name (B02 ...); red edge = B05 (705 nm, 20 m).")),
  "sentinel2-l2a-harmonized" = list(
    bands = c(blue = "B02", green = "B03", red = "B04", rededge = "B05", nir = "B08"), scale = 1e-4, offset = 0, nodata = 0,
    note = paste("Sentinel-2 L2A without the offset: products before baseline 04.00, or collections that",
                 "already removed it (Google Earth Engine COPERNICUS/S2_SR_HARMONIZED): reflectance = DN / 10000.")),
  "planetscope-4band" = list(
    bands = c(blue = "1", green = "2", red = "3", nir = "4"), scale = 1e-4, offset = 0, nodata = 0,
    note = "PlanetScope 4-band surface reflectance: blue, green, red, NIR; reflectance x 10000 in 16-bit integers."),
  "planetscope-8band" = list(
    bands = c(blue = "2", green = "4", red = "6", rededge = "7", nir = "8"), scale = 1e-4, offset = 0, nodata = 0,
    note = paste("PlanetScope 8-band surface reflectance: 1 coastal blue, 2 blue, 3 green I, 4 green, 5 yellow,",
                 "6 red, 7 red edge, 8 NIR; reflectance x 10000 in 16-bit integers."))
)

# Index definitions on surface reflectance (numbers or SpatRasters).
INDEX_DEFS <- list(
  ndvi   = list(needs = c("nir", "red"), label = "NDVI",
                f = function(b) (b$nir - b$red) / (b$nir + b$red)),
  ndre   = list(needs = c("nir", "rededge"), label = "NDRE",
                f = function(b) (b$nir - b$rededge) / (b$nir + b$rededge)),
  gndvi  = list(needs = c("nir", "green"), label = "GNDVI",
                f = function(b) (b$nir - b$green) / (b$nir + b$green)),
  osavi  = list(needs = c("nir", "red"), label = "OSAVI",
                f = function(b) (b$nir - b$red) / (b$nir + b$red + 0.16)),
  evi2   = list(needs = c("nir", "red"), label = "EVI2",
                f = function(b) 2.5 * (b$nir - b$red) / (b$nir + 2.4 * b$red + 1)),
  mcari2 = list(needs = c("nir", "red", "green"), label = "MCARI2",
                f = function(b) 1.5 * (2.5 * (b$nir - b$red) - 1.3 * (b$nir - b$green)) /
                  sqrt((2 * b$nir + 1)^2 - (6 * b$nir - 5 * sqrt(b$red)) - 0.5)),
  cire   = list(needs = c("nir", "rededge"), label = "CIred-edge",
                f = function(b) b$nir / b$rededge - 1)
)

need_pkgs <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop("Missing R package(s): ", paste(miss, collapse = ", "), ". Install with ",
         "install.packages(c(", paste0('"', miss, '"', collapse = ", "), ")).", call. = FALSE)
  invisible(TRUE)
}

norm_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))

role_alias <- function(x) {
  n <- norm_name(x)
  out <- rep(NA_character_, length(n))
  out[n %in% c("blue", "b")] <- "blue"
  out[n %in% c("green", "g")] <- "green"
  out[n %in% c("red", "r")] <- "red"
  out[n %in% c("rededge", "re", "redge")] <- "rededge"
  out[n %in% c("nir", "nearinfrared", "nearir")] <- "nir"
  out
}

# Role implied by a raster layer name (NA if the name says nothing).
role_from_layer_name <- function(nm) {
  n <- norm_name(nm)
  out <- rep(NA_character_, length(n))
  out[grepl("rededge|^re$|redge", n) | n %in% c("b05", "b06", "b07", "b5", "b6", "b7")] <- "rededge"
  out[is.na(out) & (grepl("nir|nearinfrared", n) | n %in% c("b08", "b8a", "b8"))] <- "nir"
  out[is.na(out) & (grepl("^blue", n) | n %in% c("b02", "b2"))] <- "blue"
  out[is.na(out) & (n %in% c("green", "greenreflectance") | n %in% c("b03", "b3"))] <- "green"
  out[is.na(out) & (n %in% c("red", "redreflectance") | n %in% c("b04", "b4"))] <- "red"
  out
}

parse_band_map <- function(s) {
  parts <- trimws(strsplit(s, ",", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  kv <- strsplit(parts, "=", fixed = TRUE)
  if (!length(kv) || any(lengths(kv) != 2L))
    stop("--bands must look like 'blue=1,green=2,red=3,rededge=4,nir=5'; got '", s, "'.", call. = FALSE)
  roles <- role_alias(vapply(kv, `[`, "", 1))
  if (anyNA(roles))
    stop("--bands: unknown band role(s) ", paste(vapply(kv, `[`, "", 1)[is.na(roles)], collapse = ", "),
         ". Use blue, green, red, rededge, nir.", call. = FALSE)
  if (anyDuplicated(roles)) stop("--bands names a band role twice.", call. = FALSE)
  stats::setNames(trimws(vapply(kv, `[`, "", 2)), roles)
}

# Map band roles to layer indices of raster r (numbers or layer names).
resolve_bands <- function(r, map) {
  nl <- terra::nlyr(r)
  nms <- names(r)
  idx <- vapply(names(map), function(role) {
    v <- map[[role]]
    if (grepl("^[0-9]+$", v)) {
      i <- as.integer(v)
      if (i < 1L || i > nl)
        stop("--bands: ", role, "=", v, " but the raster has ", nl, " layer(s).", call. = FALSE)
      return(i)
    }
    i <- match(norm_name(v), norm_name(nms))
    if (is.na(i))
      stop("--bands: no raster layer named '", v, "' (layers: ", paste(nms, collapse = ", "),
           "). Give layer numbers instead.", call. = FALSE)
    i
  }, integer(1))
  if (anyDuplicated(idx)) stop("--bands maps two band roles to the same raster layer.", call. = FALSE)
  idx
}

# Stop if a layer's own name contradicts the role it is mapped to.
check_name_conflicts <- function(r, idx) {
  ln <- role_from_layer_name(names(r)[idx])
  bad <- !is.na(ln) & ln != names(idx)
  if (any(bad))
    stop("Band order conflict: ", paste(sprintf("layer %d is named '%s' but is mapped to %s",
                                               idx[bad], names(r)[idx[bad]], names(idx)[bad]), collapse = "; "),
         ". Fix --bands (and confirm the order with whoever processed the imagery).", call. = FALSE)
  invisible(TRUE)
}

infer_bands_from_names <- function(r) {
  rl <- role_from_layer_name(names(r))
  keep <- !is.na(rl) & !duplicated(rl) & !(rl %in% rl[duplicated(rl)])
  if (!all(c("red", "nir") %in% rl[keep])) return(NULL)
  stats::setNames(as.character(which(keep)), rl[keep])
}

# Choose the band map: user > sensor preset > layer names; record the source.
choose_bands <- function(r, opt) {
  if (nzchar(opt$bands)) return(list(map = parse_band_map(opt$bands), source = "user (--bands)"))
  if (nzchar(opt$sensor)) {
    p <- SENSORS[[opt$sensor]]
    return(list(map = p$bands, source = paste0("sensor preset '", opt$sensor, "' (unconfirmed)")))
  }
  m <- infer_bands_from_names(r)
  if (!is.null(m)) return(list(map = m, source = "raster layer names"))
  stop("Band order is unknown: give --bands (e.g. 'blue=1,green=2,red=3,rededge=4,nir=5') after ",
       "confirming the order with the user, or a --sensor preset as a starting point. Layers: ",
       paste(names(r), collapse = ", "), ".", call. = FALSE)
}

# Otsu (1979) threshold from a histogram: maximize the between-class variance
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

# Reflectance from raw values: no-data value -> NA, then value x scale + offset.
to_reflectance <- function(r, nodata = NA, scale = 1, offset = 0) {
  if (!is.na(nodata)) r <- terra::classify(r, cbind(nodata, NA))
  if (scale != 1 || offset != 0) r <- r * scale + offset
  r
}

# Per-pixel indices (matrix: pixels x indices); non-finite values become NA.
index_matrix <- function(b, indices) {
  out <- do.call(cbind, lapply(indices, function(ix) INDEX_DEFS[[ix]]$f(b)))
  if (is.null(dim(out))) out <- matrix(out, ncol = length(indices))
  colnames(out) <- indices
  out[!is.finite(out)] <- NA
  out
}

index_raster <- function(refl, indices) {
  roles <- names(refl)
  terra::lapp(refl, function(...) {
    b <- list(...)
    names(b) <- roles
    index_matrix(b, indices)
  })
}

pick_indices <- function(spec, roles) {
  possible <- names(INDEX_DEFS)[vapply(INDEX_DEFS, function(d) all(d$needs %in% roles), logical(1))]
  if (tolower(spec) == "auto") {
    if (!length(possible)) stop("The bands given cannot make any index (need at least red and NIR).", call. = FALSE)
    return(possible)
  }
  want <- tolower(trimws(strsplit(spec, ",", fixed = TRUE)[[1]]))
  unknown <- setdiff(want, names(INDEX_DEFS))
  if (length(unknown))
    stop("Unknown index: ", paste(unknown, collapse = ", "), ". Choose from ",
         paste(names(INDEX_DEFS), collapse = ", "), ".", call. = FALSE)
  missing <- setdiff(want, possible)
  if (length(missing))
    stop("Index(es) ", paste(missing, collapse = ", "), " need band(s) not in --bands (",
         paste(unique(unlist(lapply(INDEX_DEFS[missing], `[[`, "needs"))), collapse = ", "), ").", call. = FALSE)
  want
}

# Plot polygons: checks, id column, buffer in meters, then the raster CRS.
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
  # Buffer in meters: the plots' own CRS if projected, else the raster's, else UTM.
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
  area <- as.numeric(sf::st_area(pm))
  pr <- sf::st_transform(pm, rcrs)
  list(plots = pr, id_col = id_col, area_m2 = area, crs_work = work)
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

# Coverage-weighted summaries of value rasters over the union of all plots.
band_summary <- function(x, plots) {
  u <- sf::st_sf(geometry = sf::st_union(sf::st_geometry(plots)))
  q <- c(0.01, 0.05, 0.5, 0.95, 0.99)
  s <- exactextractr::exact_extract(x, u, c("min", "max", "quantile", "count"), quantiles = q,
                                    progress = FALSE, force_df = TRUE, full_colnames = TRUE)
  do.call(rbind, lapply(names(x), function(b) {
    g <- function(f) s[[paste0(f, ".", b)]]
    data.frame(band = b, min = g("min"), q01 = g("q01"), q05 = g("q05"), median = g("q50"),
               q95 = g("q95"), q99 = g("q99"), max = g("max"), n_pixels = g("count"))
  }))
}

# Random sample of pixel values (cells at least half inside a plot).
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

run_extraction <- function(opt) {
  need_pkgs(c("terra", "sf", "exactextractr"))
  if (!file.exists(opt$input)) stop("Raster not found: ", opt$input, call. = FALSE)
  if (nzchar(opt$sensor) && is.null(SENSORS[[opt$sensor]]))
    stop("Unknown --sensor '", opt$sensor, "'. Choose from ", paste(names(SENSORS), collapse = ", "), ".",
         call. = FALSE)
  if (!opt$mask %in% c("none", "ndvi", "otsu")) stop("--mask must be none, ndvi or otsu.", call. = FALSE)
  if (nzchar(opt$date) && is.na(as.Date(opt$date, format = "%Y-%m-%d")))
    stop("--date must be YYYY-MM-DD, got '", opt$date, "'.", call. = FALSE)
  warn <- character()
  notes <- character()

  r0 <- terra::rast(opt$input)
  if (!nzchar(terra::crs(r0))) stop("The raster has no CRS; it cannot be matched to the plots.", call. = FALSE)
  if (terra::is.lonlat(r0))
    warn <- c(warn, "The raster is in lon/lat degrees; plots were buffered in meters (UTM) and transformed to it. Prefer a projected raster.")
  pp <- prepare_plots(opt$plots, opt$layer, opt$id_col, r0, opt$buffer)
  plots <- pp$plots
  ids <- plots[[pp$id_col]]

  bm <- choose_bands(r0, opt)
  idx <- resolve_bands(r0, bm$map)
  check_name_conflicts(r0, idx)
  if (!grepl("^user", bm$source))
    warn <- c(warn, sprintf("Band order taken from %s: %s. Confirm with the person who processed the imagery; stacked orthomosaics are often reordered.",
                            bm$source, paste(sprintf("%s = layer %d", names(idx), idx), collapse = ", ")))
  preset <- if (nzchar(opt$sensor)) SENSORS[[opt$sensor]] else NULL
  if (!is.null(preset)) notes <- c(notes, preset$note)
  scale <- if (!is.na(opt$scale)) opt$scale else if (!is.null(preset)) preset$scale else 1
  offset <- if (!is.na(opt$offset)) opt$offset else if (!is.null(preset)) preset$offset else 0
  nodata <- if (!is.na(opt$nodata)) opt$nodata else if (!is.null(preset)) preset$nodata else NA
  indices <- pick_indices(opt$indices, names(idx))

  pv <- terra::vect(plots)
  e1 <- as.vector(terra::ext(r0)); e2 <- as.vector(terra::ext(pv))
  if (!(e2[1] < e1[2] && e2[2] > e1[1] && e2[3] < e1[4] && e2[4] > e1[3]))
    stop("The plots do not overlap the raster (raster CRS ", terra::crs(r0, describe = TRUE)$code,
         "). Check that both are georeferenced to the same place.", call. = FALSE)
  r <- terra::crop(r0[[idx]], terra::ext(pv) + 3 * max(terra::res(r0)))
  names(r) <- names(idx)

  raw <- band_summary(r, plots)
  refl <- to_reflectance(r, nodata, scale, offset)
  names(refl) <- names(idx)
  rs <- band_summary(refl, plots)
  dn_like <- any(rs$q99 > 1.5, na.rm = TRUE)
  if (dn_like) {
    msg <- sprintf(paste("Values look like digital numbers, not surface reflectance: 99th percentile after scaling is %s",
                         "(%s). Indices need reflectance (0-1). If the file stores reflectance x 10000 (Sentinel-2 L2A,",
                         "PlanetScope SR) give --scale 0.0001 (Sentinel-2 from baseline 04.00 also --offset -0.1, or use",
                         "--sensor); reflectance in percent needs --scale 0.01. Raw camera DN must be radiometrically",
                         "calibrated first (see the claw-ortho-qc processing guide)."),
                   paste(signif(rs$q99, 4), collapse = ", "), paste(rs$band, collapse = ", "))
    if (!opt$confirm_calibrated) stop(msg, " Use --confirm-calibrated only if the user confirms the values are calibrated reflectance.", call. = FALSE)
    warn <- c(warn, paste("Proceeding on --confirm-calibrated although", sub("^Values", "values", msg)))
  }
  if (any(rs$q05 < 0, na.rm = TRUE))
    warn <- c(warn, sprintf("More than 5%% of pixels in plots are negative in band(s) %s after scaling: an offset applied twice, or a calibration problem?",
                            paste(rs$band[rs$q05 < 0], collapse = ", ")))
  if (any(rs$q99 > 1 & rs$q99 <= 1.5, na.rm = TRUE))
    warn <- c(warn, sprintf("More than 1%% of pixels exceed reflectance 1 in band(s) %s: specular glare or a calibration problem.",
                            paste(rs$band[rs$q99 > 1 & rs$q99 <= 1.5], collapse = ", ")))
  if ("nir" %in% rs$band && rs$q99[rs$band == "nir"] < 0.05)
    warn <- c(warn, "NIR reflectance is below 0.05 almost everywhere: scale applied twice, or not reflectance?")

  idx_r <- index_raster(refl, indices)
  have_ndvi <- all(c("red", "nir") %in% names(refl))
  ndvi_r <- if (have_ndvi) (if ("ndvi" %in% indices) idx_r[["ndvi"]] else index_raster(refl, "ndvi")) else NULL
  ndvi_s <- if (have_ndvi) sample_values(ndvi_r, plots, opt$sample_size, opt$seed) else numeric()
  if (length(ndvi_s)) {
    q <- stats::quantile(ndvi_s, c(0.05, 0.95), names = FALSE)
    if (q[1] < -0.3)
      warn <- c(warn, sprintf("5%% of plot pixels have NDVI below %.2f: red and NIR may be swapped. Check the band order.", q[1]))
    if (q[2] < 0.2)
      warn <- c(warn, sprintf("95th percentile NDVI is only %.2f: no vegetation signal (bare soil, wrong bands, or not reflectance).", q[2]))
  }

  thr <- NA_real_
  eta <- NA_real_
  keep <- NULL
  if (opt$mask != "none" || !is.na(opt$shadow_nir)) {
    if (opt$mask != "none" && !have_ndvi) stop("--mask needs red and NIR bands.", call. = FALSE)
    if (opt$mask == "ndvi") {
      if (is.na(opt$mask_threshold)) stop("--mask ndvi needs --mask-threshold (e.g. 0.3; see SKILL.md).", call. = FALSE)
      thr <- opt$mask_threshold
    } else if (opt$mask == "otsu") {
      thr <- otsu_threshold(ndvi_s)
      if (is.na(thr)) stop("Otsu threshold could not be computed (NDVI has no spread).", call. = FALSE)
      eta <- attr(thr, "eta")
      thr <- as.numeric(thr)
      if (thr < 0.1 || thr > 0.8 || eta < 0.8)
        warn <- c(warn, sprintf(paste("Otsu NDVI threshold %.3f (separability %.2f; 0.64 for a single normal mode): the",
                                      "histogram is probably not bimodal (closed canopy or bare soil). Prefer --mask none",
                                      "or a fixed threshold."), thr, eta))
    }
    keep <- if (!is.na(thr)) ndvi_r >= thr else NULL
    if (!is.na(opt$shadow_nir)) {
      if (!"nir" %in% names(refl)) stop("--shadow-nir needs a NIR band.", call. = FALSE)
      sh <- refl[["nir"]] >= opt$shadow_nir
      keep <- if (is.null(keep)) sh else keep & sh
    }
  }
  idx_m <- if (is.null(keep)) idx_r else terra::mask(idx_r, keep, maskvalues = c(0, NA))
  refl_m <- if (is.null(keep)) refl else terra::mask(refl, keep, maskvalues = c(0, NA))

  st <- exactextractr::exact_extract(idx_m, plots, c("mean", "median", "quantile", "stdev", "count"),
                                     quantiles = c(0.1, 0.9), progress = FALSE, force_df = TRUE,
                                     full_colnames = TRUE)
  n_all <- exactextractr::exact_extract(idx_r[[1]], plots, "count", progress = FALSE)
  n_pure <- exactextractr::exact_extract(refl[[1]], plots, function(values, coverage_fraction)
    sum(coverage_fraction >= 0.999999 & !is.na(values)), progress = FALSE)
  bands_mean <- exactextractr::exact_extract(refl_m, plots, "mean", progress = FALSE, force_df = TRUE,
                                             full_colnames = TRUE)

  long <- do.call(rbind, lapply(indices, function(ix) {
    g <- function(f) st[[paste0(f, ".", ix)]]
    data.frame(id = ids, index = ix, mean = g("mean"), median = g("median"), p10 = g("q10"),
               p90 = g("q90"), sd = g("stdev"), n_pixels = g("count"), stringsAsFactors = FALSE)
  }))
  n_eff <- st[[paste0("count.", indices[1])]]
  pa <- pixel_area_m2(r)
  wide <- data.frame(id = ids, stringsAsFactors = FALSE)
  for (ix in indices) wide[[ix]] <- st[[paste0("mean.", ix)]]
  wide$n_pixels <- n_eff
  wide$n_pixels_all <- n_all
  wide$veg_fraction <- ifelse(n_all > 0, n_eff / n_all, NA_real_)
  wide$n_pure_pixels <- n_pure
  wide$plot_area_m2 <- pp$area_m2
  wide$pixels_per_plot <- pp$area_m2 / pa
  flag <- rep("", nrow(wide))
  addf <- function(cond, txt) ifelse(cond, ifelse(nzchar(flag), paste(flag, txt, sep = "; "), txt), flag)
  flag <- addf(n_all == 0, "no data")
  flag <- addf(n_all > 0 & n_eff < opt$min_pixels, "few pixels")
  flag <- addf(n_all > 0 & n_pure == 0, "no whole pixel inside plot")
  flag <- addf(n_all > 0 & n_all * pa < 0.95 * pp$area_m2, "incomplete image coverage")
  wide$flag <- flag
  bandtab <- data.frame(id = ids, bands_mean[paste0("mean.", names(refl))])
  names(bandtab) <- c("id", paste0(names(refl), "_reflectance"))

  attrs <- sf::st_drop_geometry(plots)
  extra <- attrs[, setdiff(names(attrs), pp$id_col), drop = FALSE]
  extra <- extra[, !tolower(names(extra)) %in% c(indices, "date", "id", "flag", "n_pixels"), drop = FALSE]
  finalize <- function(df) {
    names(df)[1] <- pp$id_col
    if (nzchar(opt$date)) df <- cbind(df[1], date = opt$date, df[-1])
    df
  }
  long <- finalize(long)
  bandtab <- finalize(bandtab)
  wide <- finalize(cbind(wide[1], extra, wide[-1]))

  med_px <- stats::median(pp$area_m2 / pa)
  if (med_px < 10)
    warn <- c(warn, sprintf(paste("Mixed pixels: a median plot covers only %.1f pixel(s) of %.2f m2. Plot values mix",
                                  "neighboring plots, alleys and soil; treat them as unreliable (rule of thumb: at least",
                                  "10 whole pixels inside each plot)."), med_px, pa))
  nf <- sum(nzchar(flag))
  if (nf) warn <- c(warn, sprintf("%d of %d plot(s) flagged (%s). See the flag column.", nf, nrow(wide),
                                  paste(names(table(unlist(strsplit(flag[nzchar(flag)], "; ")))), collapse = ", ")))
  list(wide = wide, long = long, bands = bandtab, raw = raw, refl_stats = rs, plots = plots,
       indices = indices, idx = idx, band_source = bm$source, scale = scale, offset = offset,
       nodata = nodata, threshold = thr, eta = eta, warnings = warn, notes = notes, ndvi_sample = ndvi_s,
       idx_r = idx_r, keep = keep, pixel_area = pa, id_col = pp$id_col, raster = r0, crs_work = pp$crs_work)
}

# ---------------------------------------------------------------------------
# Figures, report, main
# ---------------------------------------------------------------------------

draw_figures <- function(res, opt, out_dir) {
  first <- if ("ndvi" %in% res$indices) "ndvi" else res$indices[1]
  open_png(out_dir, "index_map.png", width = 7, height = 6)
  ras <- res$idx_r[[first]]
  if (!is.null(res$keep)) ras <- terra::mask(ras, res$keep, maskvalues = c(0, NA))
  terra::plot(ras, main = sprintf("%s%s", toupper(first), if (!is.null(res$keep)) " (masked pixels removed)" else ""),
              col = grDevices::hcl.colors(50, "RdYlGn"), mar = c(3, 3, 3, 5))
  terra::plot(terra::vect(res$plots), add = TRUE, border = "black", lwd = 0.8)
  grDevices::dev.off()
  if (length(res$ndvi_sample)) {
    open_png(out_dir, "ndvi_histogram.png", width = 7, height = 4.5)
    graphics::hist(res$ndvi_sample, breaks = 100, col = "grey80", border = NA, main = "NDVI of plot pixels (sample)",
                   xlab = "NDVI")
    if (!is.na(res$threshold)) {
      graphics::abline(v = res$threshold, col = "#cc0000", lwd = 2)
      graphics::legend("topleft", bty = "n", legend = sprintf("mask threshold %.3f (%s)", res$threshold, opt$mask),
                       col = "#cc0000", lwd = 2)
    }
    grDevices::dev.off()
  }
  w <- res$wide
  open_png(out_dir, "plot_means.png", width = 7, height = max(4, min(10, 0.12 * nrow(w) + 2)))
  o <- order(w[[first]])
  graphics::dotchart(w[[first]][o], labels = if (nrow(w) <= 80) w[[res$id_col]][o] else NULL,
                     xlab = sprintf("plot mean %s", toupper(first)), pch = 19, cex = 0.7,
                     col = ifelse(nzchar(w$flag[o]), "#cc0000", "black"),
                     main = "Plot means (red = flagged)")
  grDevices::dev.off()
}

fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))

write_report <- function(res, opt, out_dir) {
  w <- res$wide
  r <- res$raster
  lines <- c(sprintf("# %s report", SKILL), "")
  if (length(res$warnings)) lines <- c(lines, "## Warnings", "", paste("-", res$warnings), "")
  lines <- c(lines, "## Inputs", "",
             sprintf("- Raster: `%s` (%d layers, %s, pixel %s x %s CRS units, %.4f m2).", opt$input, terra::nlyr(r),
                     terra::crs(r, describe = TRUE)$name, signif(terra::res(r)[1], 4), signif(terra::res(r)[2], 4), res$pixel_area),
             sprintf("- Plots: `%s`, %d polygons, id column `%s`%s.", opt$plots, nrow(w), res$id_col,
                     if (opt$buffer != 0) sprintf(", buffered by %.3f m", opt$buffer) else ""),
             sprintf("- Bands (%s): %s.", res$band_source, paste(sprintf("%s = layer %d (%s)", names(res$idx), res$idx,
                                                                         names(r)[res$idx]), collapse = "; ")),
             sprintf("- Reflectance = value x %s + %s; no-data value %s.", format(res$scale), format(res$offset),
                     if (is.na(res$nodata)) "from the file" else format(res$nodata)),
             if (nzchar(opt$date)) sprintf("- Date: %s.", opt$date) else NULL,
             if (length(res$notes)) paste("- Sensor note:", res$notes) else NULL, "",
             "## Reflectance check (pixels inside plots)", "",
             "| band | min | 1% | median | 99% | max |", "|---|---|---|---|---|---|",
             sprintf("| %s | %s | %s | %s | %s | %s |", res$refl_stats$band, fmt(res$refl_stats$min, 4),
                     fmt(res$refl_stats$q01, 4), fmt(res$refl_stats$median, 4), fmt(res$refl_stats$q99, 4),
                     fmt(res$refl_stats$max, 4)), "")
  if (!is.na(res$threshold) || !is.null(res$keep))
    lines <- c(lines, "## Mask", "",
               sprintf("- Method: %s%s%s%s.", opt$mask, if (!is.na(res$threshold)) sprintf(", NDVI threshold %.3f (one threshold for the whole image)", res$threshold) else "",
                       if (!is.na(res$eta)) sprintf(", Otsu separability %.2f", res$eta) else "",
                       if (!is.na(opt$shadow_nir)) sprintf("; shadow mask NIR < %.3f", opt$shadow_nir) else ""),
               sprintf("- Fraction of plot pixels kept: median %.2f (range %.2f-%.2f).", stats::median(w$veg_fraction, na.rm = TRUE),
                       min(w$veg_fraction, na.rm = TRUE), max(w$veg_fraction, na.rm = TRUE)),
               "", "![NDVI histogram](figures/ndvi_histogram.png)", "")
  lines <- c(lines, "## Plot results", "",
             "| index | min | median | max | plots with values |", "|---|---|---|---|---|",
             vapply(res$indices, function(ix) sprintf("| %s | %s | %s | %s | %d |", INDEX_DEFS[[ix]]$label,
                                                      fmt(min(w[[ix]], na.rm = TRUE)), fmt(stats::median(w[[ix]], na.rm = TRUE)),
                                                      fmt(max(w[[ix]], na.rm = TRUE)), sum(!is.na(w[[ix]]))), ""),
             "", sprintf("Pixels per plot: median %.1f effective (coverage-weighted), %.0f whole pixels inside.",
                         stats::median(w$n_pixels, na.rm = TRUE), stats::median(w$n_pure_pixels, na.rm = TRUE)),
             "", "![Index map](figures/index_map.png)", "", "![Plot means](figures/plot_means.png)", "",
             "## Outputs", "",
             "- `tables/plot_indices.csv`: one row per plot: plot attributes, the coverage-weighted mean of each index (columns named after the index), pixel counts, `veg_fraction` and `flag`.",
             "- `tables/plot_index_stats.csv`: one row per plot and index: mean, median, p10, p90, sd, n_pixels.",
             "- `tables/plot_band_reflectance.csv`: mean reflectance per band and plot (after masking).",
             "- `tables/image_bands.csv`: raw and reflectance distributions per band inside the plots.",
             "- `spatial/plot_indices.gpkg`: the plot polygons used (after buffering) with the plot means.",
             "", "## Method notes", "",
             "- Statistics are means of per-pixel index values, weighted by the fraction of each pixel inside the plot (exactextractr). The index of plot-mean reflectance can differ for nonlinear indices; band means are given for that purpose.",
             "- MCARI2 was defined for narrow bands at 800, 670 and 550 nm; broad multispectral bands approximate it.",
             "- OSAVI is computed as (NIR - red)/(NIR + red + 0.16) as in Rondeaux et al. (1996); some software multiplies by 1.16.",
             "", "## References", "",
             "- Rouse JW, Haas RH, Schell JA, Deering DW (1974) Monitoring vegetation systems in the Great Plains with ERTS. Third ERTS-1 Symposium, NASA Goddard Space Flight Center. https://ntrs.nasa.gov/citations/19740022614",
             "- Gitelson A, Merzlyak MN (1994) J Plant Physiol 143:286-292. https://doi.org/10.1016/S0176-1617(11)81633-0",
             "- Gitelson AA, Kaufman YJ, Merzlyak MN (1996) Remote Sens Environ 58:289-298. https://doi.org/10.1016/S0034-4257(96)00072-7",
             "- Rondeaux G, Steven M, Baret F (1996) Remote Sens Environ 55:95-107. https://doi.org/10.1016/0034-4257(95)00186-7",
             "- Jiang Z, Huete AR, Didan K, Miura T (2008) Remote Sens Environ 112:3833-3845. https://doi.org/10.1016/j.rse.2008.06.006",
             "- Haboudane D, Miller JR, Pattey E, Zarco-Tejada PJ, Strachan IB (2004) Remote Sens Environ 90:337-352. https://doi.org/10.1016/j.rse.2003.12.013",
             "- Gitelson AA, Viña A, Arkebauer TJ, Rundquist DC, Keydan G, Leavitt B (2003) Geophys Res Lett 30(5). https://doi.org/10.1029/2002GL016450",
             "- Otsu N (1979) IEEE Trans Syst Man Cybern 9:62-66. https://doi.org/10.1109/TSMC.1979.4310076",
             "", sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(lines, file.path(out_dir, "report.md"))
}

demo_options <- function(opt) {
  src <- file.path(SKILL_DIR, "demo", "make_demo.R")
  if (!file.exists(src)) stop("Demo generator not found: ", src, call. = FALSE)
  e <- new.env()
  sys.source(src, envir = e)
  d <- file.path(opt$out, "demo_input")
  paths <- e$make_demo(d)
  opt$input <- paths$raster
  opt$plots <- paths$plots
  opt$bands <- "blue=1,green=2,red=3,rededge=4,nir=5"
  opt$scale <- 1e-4
  opt$mask <- "otsu"
  opt$buffer <- -0.25
  opt$date <- "2026-07-15"
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
    stop("Provide --input RASTER and --plots POLYGONS, or --demo. See --help.", call. = FALSE)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  if (opt$demo) opt <- demo_options(opt)

  res <- run_extraction(opt)
  save_table(res$wide, opt$out, "plot_indices.csv")
  save_table(res$long, opt$out, "plot_index_stats.csv")
  save_table(res$bands, opt$out, "plot_band_reflectance.csv")
  img <- merge(transform(res$raw, stage = "raw"), transform(res$refl_stats, stage = "reflectance"), all = TRUE)
  save_table(img[order(img$stage, img$band), c("stage", setdiff(names(img), "stage"))], opt$out, "image_bands.csv")
  gp <- sf::st_sf(res$wide, geometry = sf::st_geometry(res$plots))
  dir.create(file.path(opt$out, "spatial"), showWarnings = FALSE)
  gpath <- file.path(opt$out, "spatial", "plot_indices.gpkg")
  if (file.exists(gpath)) unlink(gpath)
  sf::st_write(gp, gpath, layer = "plots", quiet = TRUE)
  draw_figures(res, opt, opt$out)
  write_report(res, opt, opt$out)
  write_repro_bundle(opt$out, argv, inputs = c(opt$input, opt$plots))
  cat(sprintf("Done: %d plots, %d indices. Report: %s\n", nrow(res$wide), length(res$indices),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  if (length(res$warnings)) cat(paste0("Warning: ", res$warnings, "\n"), sep = "")
  invisible(res)
}

if (sys.nframe() == 0L) main()
