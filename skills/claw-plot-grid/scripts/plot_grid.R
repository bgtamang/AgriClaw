#!/usr/bin/env Rscript
# claw-plot-grid: plot polygons for a field trial from the four field corners
# (or one origin corner and a bearing) and the layout (ranges x columns, plot
# length, plot width, alleys). Polygons are built in metres in a projected CRS
# (the local UTM zone by default), never in degrees.
#
# Usage:
#   Rscript scripts/plot_grid.R --demo
#   Rscript scripts/plot_grid.R --corners corners.csv --ranges 6 --columns 8 \
#       --plot-length 5.2 --alley 0.9 --rows-per-plot 2 --row-spacing 0.76 --units m
#   Rscript scripts/plot_grid.R --origin "-95.80,41.30" --bearing 12 --across right ...
#   Rscript scripts/plot_grid.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-plot-grid"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/plot_grid.R [--demo] (--corners C | --origin X,Y --bearing DEG)
          --ranges N --columns N --plot-length L --alley A
          (--rows-per-plot N --row-spacing S | --plot-width W) --units m|ft [options]

Field geometry (give one):
  --corners      4 plot-block corners in order c1;c2;c3;c4 as 'x1,y1;x2,y2;x3,y3;x4,y4'
                 or a CSV with columns x,y (or lon,lat). c1 = corner of range 1 / column 1,
                 c2 = far end along the planting direction, c3 = diagonal, c4 = across from c1
  --origin       corner c1 as 'x,y' (lon,lat order for EPSG:4326), with --bearing
  --bearing      direction of planting (range 1 -> range N), degrees clockwise from north
  --bearing-ref  true | grid | auto (default auto: true for lon/lat input, grid otherwise)
  --across       right | left: side of the planting direction where columns increase
                 (origin mode; default right)
  --offset-along, --offset-across  distance from the origin point to the first plot
                 corner (origin mode; default 0)
  --crs          CRS of the input coordinates (default EPSG:4326; 'local' = unreferenced metres)
  --out-crs      CRS of the output polygons: auto (local UTM zone, default) or e.g. EPSG:32615
  --fit-corners  stretch the layout to the 4 corners (bilinear) instead of a rigid fit
  --corner-tol   max corner misfit in metres before stopping (default 0.5)
  --scale-factor auto (default: map-projection point scale factor) or a number such as 1

Layout (lengths in --units; --units is required):
  --ranges, --columns          number of plots along / across the planting direction
  --plot-length                plot length along the planting direction
  --alley                      gap between consecutive ranges (give 0 if none)
  --rows-per-plot, --row-spacing   plot width = rows x spacing (or give --plot-width)
  --col-gap                    extra gap between adjacent columns of plots (default 0)
  --units                      m or ft

Numbering:
  --numbering    serpentine | rowwise (default serpentine)
  --sweep        along | across: consecutive plots run along a planter pass (through the
                 ranges of one column) or across a range (default along)
  --start-corner 1-4: corner where plot --first-plot is (default 1 = c1)
  --first-plot   number of the first plot (default 1)

Trimming and field book:
  --trim-ends, --trim-sides   inner polygons trimmed by these distances (per end / side)
  --fieldbook    CSV field book to join (e.g. from claw-trial-design)
  --join-by      'plot' (default) or 'FBCOL=gridcol,...' with gridcol in plot, range, column
  --fb-coords    'ROWCOL,COLCOL': field-book position columns to cross-check the numbering
  --out          output folder (default agriclaw_output/claw-plot-grid)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, corners = "", origin = "", crs = "EPSG:4326",
                 out_crs = "auto", bearing = NA_real_, bearing_ref = "auto",
                 across = "right", offset_along = 0, offset_across = 0,
                 ranges = NA_real_, columns = NA_real_, plot_length = NA_real_,
                 alley = NA_real_, rows_per_plot = NA_real_, row_spacing = NA_real_,
                 plot_width = NA_real_, col_gap = 0, units = "",
                 numbering = "serpentine", sweep = "along", start_corner = 1,
                 first_plot = 1, trim_ends = 0, trim_sides = 0, fit_corners = FALSE,
                 corner_tol = 0.5, scale_factor = "auto", fieldbook = "",
                 join_by = "plot", fb_coords = "", out = "")

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

FT_TO_M <- 0.3048            # international foot, exact
WGS84_A <- 6378137           # WGS84 semi-major axis (m)
WGS84_F <- 1 / 298.257223563 # WGS84 flattening
GRID_COLS <- c("plot", "range", "column")
# Names the output already uses; field-book columns with these names (in any
# letter case, since GeoPackage column names are case-insensitive) get 'fb_'.
RESERVED <- c("plot", "range", "column", "plot_id", "length_m", "width_m", "area_m2",
              "x", "y", "lon", "lat", "fid", "geom", "geometry", "matched")

need_pkgs <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop("Missing R package(s): ", paste(miss, collapse = ", "), ". Install with ",
         "install.packages(c(", paste0('"', miss, '"', collapse = ", "), ")).", call. = FALSE)
  invisible(TRUE)
}

# Options the user typed (to tell explicit choices from silent defaults).
given_options <- function(argv) {
  k <- argv[startsWith(argv, "--")]
  gsub("-", "_", substring(k, 3L), fixed = TRUE)
}

# EPSG code of the WGS84 UTM zone containing (lon, lat), with the standard
# exceptions for south-west Norway (zone 32V) and Svalbard (31X/33X/35X/37X).
utm_epsg <- function(lon, lat) {
  if (!is.finite(lon) || !is.finite(lat))
    stop("Longitude/latitude must be finite numbers.", call. = FALSE)
  if (lat < -80 || lat > 84)
    stop("UTM is defined only between 80 deg S and 84 deg N; give --out-crs ",
         "(for example a polar stereographic EPSG code).", call. = FALSE)
  lon <- ((lon + 180) %% 360) - 180
  zone <- floor((lon + 180) / 6) + 1
  if (lat >= 56 && lat < 64 && lon >= 3 && lon < 12) zone <- 32
  if (lat >= 72 && lat < 84) {
    if (lon >= 0 && lon < 9) zone <- 31
    else if (lon >= 9 && lon < 21) zone <- 33
    else if (lon >= 21 && lon < 33) zone <- 35
    else if (lon >= 33 && lon < 42) zone <- 37
  }
  as.integer((if (lat >= 0) 32600L else 32700L) + zone)
}

parse_xy <- function(s, what = "--origin") {
  v <- suppressWarnings(as.numeric(strsplit(trimws(s), "[,[:space:]]+")[[1]]))
  if (length(v) != 2L || anyNA(v))
    stop(what, " needs two numbers 'x,y' (lon,lat order for EPSG:4326); got '", s, "'.",
         call. = FALSE)
  v
}

# Four corners from 'x1,y1;x2,y2;x3,y3;x4,y4' or a CSV (columns x,y or lon,lat;
# an optional 'corner' column 1-4 sets the order).
read_corners <- function(s) {
  if (grepl("\\.csv$", s, ignore.case = TRUE)) {
    if (!file.exists(s)) stop("Corners file not found: ", s, call. = FALSE)
    d <- utils::read.csv(s, check.names = FALSE, stringsAsFactors = FALSE)
    nm <- tolower(trimws(names(d)))
    if ("corner" %in% nm) d <- d[order(d[[match("corner", nm)]]), , drop = FALSE]
    pick <- function(cands) { i <- match(cands, nm); i <- i[!is.na(i)]; if (length(i)) i[1] else NA }
    ix <- pick(c("x", "lon", "longitude", "easting"))
    iy <- pick(c("y", "lat", "latitude", "northing"))
    if (is.na(ix) || is.na(iy))
      stop("Corners CSV needs columns x,y (or lon,lat or easting,northing); found: ",
           paste(names(d), collapse = ", "), call. = FALSE)
    m <- cbind(suppressWarnings(as.numeric(d[[ix]])), suppressWarnings(as.numeric(d[[iy]])))
  } else {
    parts <- strsplit(trimws(s), ";", fixed = TRUE)[[1]]
    m <- do.call(rbind, lapply(parts, parse_xy, what = "--corners"))
  }
  if (nrow(m) != 4L || anyNA(m))
    stop("--corners needs exactly 4 corners in order c1;c2;c3;c4 (got ", nrow(m),
         " usable rows).", call. = FALSE)
  m
}

# sf CRS object from a string ('EPSG:4326', '32615', WKT, PROJ); NULL = local.
parse_crs <- function(s, what) {
  if (tolower(trimws(s)) %in% c("local", "none")) return(NULL)
  if (grepl("^[0-9]+$", trimws(s))) s <- as.integer(s)
  cr <- tryCatch(sf::st_crs(s), error = function(e) sf::st_crs(NA))
  if (is.na(cr))
    stop("Unknown CRS for ", what, ": '", s, "'. Use e.g. EPSG:4326 (lon/lat) or ",
         "EPSG:32615 (UTM 15N).", call. = FALSE)
  cr
}

is_geographic <- function(cr) !is.null(cr) && isTRUE(sf::st_is_longlat(cr))

transform_xy <- function(m, from, to) {
  m <- matrix(m, ncol = 2)
  pts <- sf::st_sfc(lapply(seq_len(nrow(m)), function(i) sf::st_point(m[i, ])), crs = from)
  unname(sf::st_coordinates(sf::st_transform(pts, to))[, 1:2, drop = FALSE])
}

# Point scale factor k (map distance / ground distance on the WGS84 ellipsoid)
# and the grid bearing of true north for CRS `cr` at (lon, lat). Ground lengths
# of small steps use the ellipsoid's radii of curvature; map lengths come from
# PROJ. For a conformal projection such as UTM, k is the same in every direction.
geo_factors <- function(lon, lat, cr, step_deg = 1e-4) {
  phi <- lat * pi / 180
  e2 <- WGS84_F * (2 - WGS84_F)
  w <- sqrt(1 - e2 * sin(phi)^2)
  N <- WGS84_A / w
  M <- WGS84_A * (1 - e2) / w^3
  pts <- rbind(c(lon - step_deg, lat), c(lon + step_deg, lat),
               c(lon, lat - step_deg), c(lon, lat + step_deg))
  q <- transform_xy(pts, sf::st_crs(4326), cr)
  de <- q[2, ] - q[1, ]
  dn <- q[4, ] - q[3, ]
  step_rad <- 2 * step_deg * pi / 180
  k_e <- sqrt(sum(de^2)) / (N * cos(phi) * step_rad)
  k_n <- sqrt(sum(dn^2)) / (M * step_rad)
  list(k = (k_e + k_n) / 2, k_east = k_e, k_north = k_n,
       north_bearing = atan2(dn[1], dn[2]) * 180 / pi)
}

unit_along <- function(bearing) c(sin(bearing * pi / 180), cos(bearing * pi / 180))
unit_across <- function(bearing, across = "right") {
  a <- unit_along(bearing)
  if (across == "right") c(a[2], -a[1]) else c(-a[2], a[1])
}
norm360 <- function(x) x %% 360

# Plot numbers as a [range, column] matrix. Corners: c1 = (range 1, column 1),
# c2 = (range N, column 1), c3 = (range N, column M), c4 = (range 1, column M).
plot_numbers <- function(n_ranges, n_cols, numbering = "serpentine", sweep = "along",
                         start_corner = 1, first_plot = 1) {
  rseq <- if (start_corner %in% c(1, 4)) seq_len(n_ranges) else rev(seq_len(n_ranges))
  cseq <- if (start_corner %in% c(1, 2)) seq_len(n_cols) else rev(seq_len(n_cols))
  num <- matrix(NA_integer_, n_ranges, n_cols)
  k <- as.integer(first_plot)
  if (sweep == "along") {
    for (li in seq_along(cseq)) {
      rr <- if (numbering == "serpentine" && li %% 2L == 0L) rev(rseq) else rseq
      for (r in rr) { num[r, cseq[li]] <- k; k <- k + 1L }
    }
  } else {
    for (li in seq_along(rseq)) {
      cc <- if (numbering == "serpentine" && li %% 2L == 0L) rev(cseq) else cseq
      for (cl in cc) { num[rseq[li], cl] <- k; k <- k + 1L }
    }
  }
  num
}

# One row per plot with its extent in local layout coordinates (metres):
# s along the planting direction from c1, t across from c1.
plot_table <- function(n_ranges, n_cols, L, A, W, G, num) {
  g <- expand.grid(range = seq_len(n_ranges), column = seq_len(n_cols))
  g$plot <- num[cbind(g$range, g$column)]
  g$s0 <- (g$range - 1) * (L + A)
  g$s1 <- g$s0 + L
  g$t0 <- (g$column - 1) * (W + G)
  g$t1 <- g$t0 + W
  g <- g[order(g$plot), c("plot", "range", "column", "s0", "s1", "t0", "t1")]
  wr <- nchar(n_ranges)
  wc <- nchar(n_cols)
  g$plot_id <- sprintf("R%0*dC%0*d", wr, g$range, wc, g$column)
  rownames(g) <- NULL
  g
}

# Least-squares rotation + translation (no scaling) taking Q onto P (2-D Procrustes).
fit_rigid <- function(P, Q) {
  pc <- colMeans(P)
  qc <- colMeans(Q)
  p <- sweep(P, 2, pc)
  q <- sweep(Q, 2, qc)
  phi <- atan2(sum(q[, 1] * p[, 2] - q[, 2] * p[, 1]), sum(q[, 1] * p[, 1] + q[, 2] * p[, 2]))
  R <- matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), 2)
  shift <- pc - as.vector(R %*% qc)
  fitted <- t(R %*% t(Q)) + matrix(shift, nrow(Q), 2, byrow = TRUE)
  list(phi = phi, fitted = fitted, resid = P - fitted)
}

bilinear_xy <- function(u, v, P) {
  w1 <- (1 - u) * (1 - v); w2 <- u * (1 - v); w3 <- u * v; w4 <- (1 - u) * v
  cbind(w1 * P[1, 1] + w2 * P[2, 1] + w3 * P[3, 1] + w4 * P[4, 1],
        w1 * P[1, 2] + w2 * P[2, 2] + w3 * P[3, 2] + w4 * P[4, 2])
}

# Is the quadrilateral c1-c2-c3-c4 simple and convex?
is_convex_quad <- function(P) {
  z <- vapply(1:4, function(i) {
    a <- P[i, ]; b <- P[(i %% 4) + 1, ]; c <- P[((i + 1) %% 4) + 1, ]
    (b[1] - a[1]) * (c[2] - b[2]) - (b[2] - a[2]) * (c[1] - b[1])
  }, numeric(1))
  all(z > 0) || all(z < 0)
}

make_polys <- function(tab, mapper, crs, trim_s = 0, trim_t = 0) {
  geoms <- lapply(seq_len(nrow(tab)), function(i) {
    s0 <- tab$s0[i] + trim_s; s1 <- tab$s1[i] - trim_s
    t0 <- tab$t0[i] + trim_t; t1 <- tab$t1[i] - trim_t
    ring <- mapper(c(s0, s1, s1, s0, s0), c(t0, t0, t1, t1, t0))
    ring[5, ] <- ring[1, ]
    sf::st_polygon(list(ring))
  })
  sf::st_sfc(geoms, crs = if (is.null(crs)) sf::NA_crs_ else crs)
}

# Total area where plot interiors overlap (should be 0) and the number of pairs.
overlap_check <- function(geom) {
  rel <- sf::st_relate(geom, geom, pattern = "2********")
  tot <- 0
  pairs <- 0L
  for (i in seq_along(rel)) for (j in rel[[i]]) if (j > i) {
    a <- sum(as.numeric(sf::st_area(sf::st_intersection(geom[i], geom[j]))))
    tot <- tot + a
    pairs <- pairs + 1L
  }
  list(area = tot, pairs = pairs)
}

parse_join <- function(s) {
  parts <- trimws(strsplit(s, ",", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  if (!length(parts)) stop("--join-by is empty.", call. = FALSE)
  fb <- character(); grid <- character()
  for (p in parts) {
    kv <- trimws(strsplit(p, "=", fixed = TRUE)[[1]])
    if (length(kv) == 1L) kv <- c(kv, kv)
    if (length(kv) != 2L || !nzchar(kv[1]) || !nzchar(kv[2]))
      stop("--join-by must look like 'plot' or 'PLOT=plot,ROW=range'; got '", s, "'.", call. = FALSE)
    if (!tolower(kv[2]) %in% GRID_COLS)
      stop("--join-by: grid column '", kv[2], "' must be one of plot, range, column.", call. = FALSE)
    fb <- c(fb, kv[1]); grid <- c(grid, tolower(kv[2]))
  }
  list(fb = fb, grid = grid)
}

norm_key <- function(x) {
  x <- trimws(as.character(x))
  num <- suppressWarnings(as.numeric(x))
  ifelse(!is.na(num), format(num, scientific = FALSE, trim = TRUE, digits = 15), x)
}

# Join a field book onto the grid table. Returns the joined table, the number
# of grid plots without a field-book row, and messages.
join_fieldbook <- function(grid_df, fb, join_by = "plot") {
  jb <- parse_join(join_by)
  miss <- setdiff(jb$fb, names(fb))
  if (length(miss))
    stop("Field book has no column(s) ", paste(miss, collapse = ", "), ". Columns: ",
         paste(names(fb), collapse = ", "), ". Set --join-by, e.g. 'PLOT=plot'.", call. = FALSE)
  fkey <- do.call(paste, c(lapply(jb$fb, function(cn) norm_key(fb[[cn]])), sep = "|"))
  gkey <- do.call(paste, c(lapply(jb$grid, function(cn) norm_key(grid_df[[cn]])), sep = "|"))
  if (anyNA(fb[jb$fb]) || any(!nzchar(trimws(unlist(fb[jb$fb])))))
    stop("Field book has empty values in the join column(s) ", paste(jb$fb, collapse = ", "),
         ".", call. = FALSE)
  dup <- unique(fkey[duplicated(fkey)])
  if (length(dup))
    stop("Field book has duplicate join keys (", paste(utils::head(dup, 10), collapse = ", "),
         "): each plot must appear once.", call. = FALSE)
  extra <- setdiff(fkey, gkey)
  if (length(extra))
    stop(length(extra), " field-book row(s) match no plot in the grid (e.g. ",
         paste(utils::head(extra, 10), collapse = ", "), "). Check --ranges, --columns, ",
         "--first-plot, the numbering options and --join-by.", call. = FALSE)
  fbx <- fb
  ren <- tolower(names(fbx)) %in% RESERVED & !(names(fbx) %in% jb$fb[jb$fb == jb$grid])
  names(fbx)[ren] <- paste0("fb_", names(fbx)[ren])
  keep <- setdiff(names(fbx), jb$fb[jb$fb == jb$grid])
  idx <- match(gkey, fkey)
  out <- grid_df
  for (cn in keep) out[[cn]] <- fbx[[cn]][idx]
  out$matched <- !is.na(idx)
  list(data = out, n_unmatched = sum(is.na(idx)), renamed = names(fb)[ren])
}

# The field book's own position columns (e.g. ROW, COLUMN) must map one-to-one
# onto the grid's range and column indices; otherwise the plot numbers were
# assigned with a different numbering scheme or start corner.
check_fb_coords <- function(joined, fb_cols) {
  d <- joined[joined$matched, , drop = FALSE]
  f1 <- d[[fb_cols[1]]]; f2 <- d[[fb_cols[2]]]
  is_bij <- function(g, f) {
    m1 <- tapply(f, g, function(z) length(unique(z)))
    m2 <- tapply(g, f, function(z) length(unique(z)))
    all(m1 == 1L) && all(m2 == 1L)
  }
  if (is_bij(d$range, f1) && is_bij(d$column, f2))
    return(list(ok = TRUE, mapping = sprintf("%s <-> range, %s <-> column", fb_cols[1], fb_cols[2])))
  if (is_bij(d$column, f1) && is_bij(d$range, f2))
    return(list(ok = TRUE, mapping = sprintf("%s <-> column, %s <-> range", fb_cols[1], fb_cols[2])))
  list(ok = FALSE, mapping = NA_character_)
}

to_m <- function(x, units) if (units == "ft") x * FT_TO_M else x

# Validate options and convert lengths to metres.
layout_from_opt <- function(opt) {
  if (!opt$units %in% c("m", "ft"))
    stop("--units is required: m or ft (all lengths use it). Ask the user which unit ",
         "the plot dimensions are in; never guess.", call. = FALSE)
  for (k in c("ranges", "columns", "plot_length", "alley"))
    if (is.na(opt[[k]])) stop("--", gsub("_", "-", k), " is required. See --help.", call. = FALSE)
  if (opt$ranges < 1 || opt$columns < 1 || opt$ranges != round(opt$ranges) ||
      opt$columns != round(opt$columns))
    stop("--ranges and --columns must be whole numbers >= 1.", call. = FALSE)
  W <- if (!is.na(opt$plot_width)) opt$plot_width else {
    if (is.na(opt$rows_per_plot) || is.na(opt$row_spacing))
      stop("Give --rows-per-plot and --row-spacing (plot width = rows x spacing), or --plot-width.",
           call. = FALSE)
    opt$rows_per_plot * opt$row_spacing
  }
  if (!is.na(opt$plot_width) && !is.na(opt$rows_per_plot) && !is.na(opt$row_spacing) &&
      abs(opt$plot_width - opt$rows_per_plot * opt$row_spacing) > 1e-9)
    stop("--plot-width disagrees with --rows-per-plot x --row-spacing; give one or the other.",
         call. = FALSE)
  L <- opt$plot_length; A <- opt$alley; G <- opt$col_gap
  if (L <= 0 || W <= 0) stop("Plot length and width must be positive.", call. = FALSE)
  if (A < 0 || G < 0) stop("--alley and --col-gap cannot be negative.", call. = FALSE)
  te <- opt$trim_ends; ts <- opt$trim_sides
  if (te < 0 || ts < 0) stop("--trim-ends and --trim-sides cannot be negative.", call. = FALSE)
  if (2 * te >= L) stop("--trim-ends removes the whole plot length.", call. = FALSE)
  if (2 * ts >= W) stop("--trim-sides removes the whole plot width.", call. = FALSE)
  u <- opt$units
  list(n_ranges = as.integer(opt$ranges), n_cols = as.integer(opt$columns),
       L = to_m(L, u), A = to_m(A, u), W = to_m(W, u), G = to_m(G, u),
       row_spacing = if (is.na(opt$row_spacing)) NA_real_ else to_m(opt$row_spacing, u),
       te = to_m(te, u), ts = to_m(ts, u),
       off_s = to_m(opt$offset_along, u), off_t = to_m(opt$offset_across, u))
}

# Plausibility of plot dimensions in metres (rule-of-thumb ranges for row-crop
# and small-grain trials), used only to warn about a likely ft/m mix-up.
dimension_warnings <- function(lay, units) {
  w <- character()
  if (lay$L < 1 || lay$L > 30)
    w <- c(w, sprintf("Plot length is %.2f m, outside the usual 1-30 m. Check --units (%s).", lay$L, units))
  if (lay$W < 0.3 || lay$W > 20)
    w <- c(w, sprintf("Plot width is %.2f m, outside the usual 0.3-20 m. Check --units (%s).", lay$W, units))
  if (!is.na(lay$row_spacing) && (lay$row_spacing < 0.1 || lay$row_spacing > 1.6))
    w <- c(w, sprintf("Row spacing is %.3f m, outside the usual 0.1-1.6 m. Check --units (%s).",
                      lay$row_spacing, units))
  if (lay$A > 6)
    w <- c(w, sprintf("Alley is %.2f m, wider than the usual 0-6 m. Check --units (%s).", lay$A, units))
  w
}

# Build the grid. Returns plots (sf), inner plots (sf or NULL), corner check,
# summary info and warnings.
build_grid <- function(opt) {
  need_pkgs("sf")
  lay <- layout_from_opt(opt)
  warn <- dimension_warnings(lay, opt$units)
  if (!opt$numbering %in% c("serpentine", "rowwise"))
    stop("--numbering must be serpentine or rowwise.", call. = FALSE)
  if (!opt$sweep %in% c("along", "across"))
    stop("--sweep must be along or across.", call. = FALSE)
  if (!opt$start_corner %in% 1:4) stop("--start-corner must be 1, 2, 3 or 4.", call. = FALSE)
  if (!opt$across %in% c("right", "left")) stop("--across must be right or left.", call. = FALSE)
  if (!opt$bearing_ref %in% c("auto", "true", "grid"))
    stop("--bearing-ref must be true, grid or auto.", call. = FALSE)
  have_c <- nzchar(opt$corners)
  have_o <- nzchar(opt$origin)
  if (have_c == have_o)
    stop("Give either --corners (4 corners) or --origin with --bearing, not both or neither.",
         call. = FALSE)

  crs_in <- parse_crs(opt$crs, "--crs")
  local <- is.null(crs_in)
  if (local && !tolower(opt$out_crs) %in% c("auto", "local", "none"))
    stop("Input coordinates are local (--crs local); the output cannot be georeferenced.",
         call. = FALSE)

  pts_in <- if (have_c) read_corners(opt$corners) else matrix(parse_xy(opt$origin), 1)
  if (is_geographic(crs_in)) {
    if (any(abs(pts_in[, 2]) > 90) || any(abs(pts_in[, 1]) > 180))
      stop("Coordinates are out of range for ", opt$crs, " (lon,lat). If they are ",
           "lat,lon swap them; if they are projected give --crs.", call. = FALSE)
  }
  lonlat <- if (local) NULL else transform_xy(pts_in, crs_in, sf::st_crs(4326))
  ref_ll <- if (local) NULL else colMeans(lonlat)

  if (local) {
    crs_out <- NULL
    out_label <- "local (no CRS)"
  } else if (tolower(opt$out_crs) == "auto") {
    crs_out <- sf::st_crs(utm_epsg(ref_ll[1], ref_ll[2]))
    out_label <- sprintf("EPSG:%d (%s, chosen automatically)", crs_out$epsg, crs_out$Name)
  } else {
    crs_out <- parse_crs(opt$out_crs, "--out-crs")
    if (is_geographic(crs_out))
      stop("--out-crs must be a projected CRS in metres, never lon/lat: buffers and areas ",
           "would be in degrees.", call. = FALSE)
    out_label <- if (!is.na(crs_out$epsg)) sprintf("EPSG:%d (%s)", crs_out$epsg, crs_out$Name) else crs_out$Name
  }
  if (!local) {
    u <- crs_out$units_gdal
    if (is.null(u) || is.na(u) || !u %in% c("metre", "meter", "m"))
      stop("--out-crs must have metre units (got '", u, "').", call. = FALSE)
  }

  gf <- if (local) list(k = 1, k_east = 1, k_north = 1, north_bearing = 0) else
    geo_factors(lonlat[1, 1], lonlat[1, 2], crs_out)
  if (abs(gf$k_east - gf$k_north) > 1e-4 * gf$k)
    warn <- c(warn, sprintf(paste("The output CRS is not conformal here (scale %.6f east-west,",
                                  "%.6f north-south): plot shapes are distorted. Prefer the UTM zone."),
                            gf$k_east, gf$k_north))
  k <- if (tolower(opt$scale_factor) == "auto") gf$k else suppressWarnings(as.numeric(opt$scale_factor))
  if (is.na(k) || k <= 0 || abs(k - 1) > 0.01)
    stop("--scale-factor must be 'auto' or a number within 1% of 1 (got '", opt$scale_factor, "').",
         call. = FALSE)

  P <- if (local) pts_in else transform_xy(pts_in, crs_in, crs_out)
  Ltot <- lay$n_ranges * lay$L + (lay$n_ranges - 1) * lay$A
  Wtot <- lay$n_cols * lay$W + (lay$n_cols - 1) * lay$G
  corner_df <- NULL
  info <- list(mode = if (have_c) "corners" else "origin", k = k, gf = gf,
               out_label = out_label, crs_out = crs_out, local = local, Ltot = Ltot,
               Wtot = Wtot, lay = lay, ref_ll = ref_ll, fit = "rigid")

  if (have_o) {
    if (is.na(opt$bearing)) stop("--origin needs --bearing (degrees clockwise from north).",
                                 call. = FALSE)
    ref <- if (opt$bearing_ref == "auto") (if (is_geographic(crs_in)) "true" else "grid") else opt$bearing_ref
    if (local && ref == "true")
      stop("--bearing-ref true needs georeferenced coordinates (not --crs local).", call. = FALSE)
    b_grid <- norm360(if (ref == "true") opt$bearing + gf$north_bearing else opt$bearing)
    across <- opt$across
    O <- P[1, ]
    off_s <- lay$off_s; off_t <- lay$off_t
    info$bearing_ref <- ref
  } else {
    if (!is_convex_quad(P))
      stop("The 4 corners do not form a simple convex quadrilateral in the order given. ",
           "List them in order around the block: c1 (range 1, column 1), c2 (far end along ",
           "the planting direction), c3 (diagonal), c4 (across from c1).", call. = FALSE)
    if (lay$off_s != 0 || lay$off_t != 0)
      stop("--offset-along/--offset-across apply only with --origin; in corners mode the ",
           "corners must be the outer corners of the plot block.", call. = FALSE)
    d12 <- P[2, ] - P[1, ]; d14 <- P[4, ] - P[1, ]
    cross <- d12[1] * d14[2] - d12[2] * d14[1]
    across <- if (cross < 0) "right" else "left"
    sgn <- if (across == "left") 1 else -1
    Q <- k * rbind(c(0, 0), c(Ltot, 0), c(Ltot, sgn * Wtot), c(0, sgn * Wtot))
    fr <- fit_rigid(P, Q)
    b_grid <- norm360(90 - fr$phi * 180 / pi)
    O <- fr$fitted[1, ]
    off_s <- 0; off_t <- 0
    ground <- function(a, b) sqrt(sum((a - b)^2)) / k
    corner_df <- data.frame(
      corner = paste0("c", 1:4), x_given = P[, 1], y_given = P[, 2],
      x_fitted = fr$fitted[, 1], y_fitted = fr$fitted[, 2],
      misfit_m = sqrt(rowSums(fr$resid^2)) / k)
    info$side_lengths <- c(c1c2 = ground(P[1, ], P[2, ]), c4c3 = ground(P[4, ], P[3, ]),
                           c1c4 = ground(P[1, ], P[4, ]), c2c3 = ground(P[2, ], P[3, ]))
    v1 <- d12 / sqrt(sum(d12^2)); v2 <- d14 / sqrt(sum(d14^2))
    info$corner_angle <- acos(max(-1, min(1, sum(v1 * v2)))) * 180 / pi
    if (opt$fit_corners) {
      info$fit <- "bilinear"
      sc <- c(along = mean(info$side_lengths[c("c1c2", "c4c3")]) / Ltot,
              across = mean(info$side_lengths[c("c1c4", "c2c3")]) / Wtot)
      info$fit_scale <- sc
      if (any(abs(sc - 1) > 0.02))
        warn <- c(warn, sprintf(paste("Fitting to the corners stretches the layout by %.1f%% along",
                                      "and %.1f%% across (more than 2%%): check the layout",
                                      "dimensions, --units and the corner positions."),
                                100 * (sc[1] - 1), 100 * (sc[2] - 1)))
      if (abs(info$corner_angle - 90) > 2)
        warn <- c(warn, sprintf("The corner angle at c1 is %.1f deg, not 90: check the corners.",
                                info$corner_angle))
    } else {
      worst <- max(corner_df$misfit_m)
      if (worst > opt$corner_tol)
        stop(sprintf(paste0("The corners and the layout disagree: worst corner misfit %.2f m ",
                            "(tolerance %.2f m). Corners span %.2f x %.2f m (c1-c2 x c1-c4) but the ",
                            "layout needs %.2f x %.2f m. Check --units (ft vs m), --ranges/--columns, ",
                            "--alley, the plot width and the corner order; or use --fit-corners if ",
                            "the corners are right and the dimensions approximate."),
                     worst, opt$corner_tol, info$side_lengths[["c1c2"]], info$side_lengths[["c1c4"]],
                     Ltot, Wtot), call. = FALSE)
      if (worst > 0.10)
        warn <- c(warn, sprintf(paste("Worst corner misfit is %.2f m (more than 0.10 m): plot edges may",
                                      "be offset by up to this much; check the corners."), worst))
    }
    info$bearing_ref <- "from corners"
  }
  info$bearing_grid <- b_grid
  info$bearing_true <- norm360(b_grid - gf$north_bearing)
  info$across <- across

  a <- unit_along(b_grid)
  cvec <- unit_across(b_grid, across)
  mapper <- if (info$fit == "bilinear") {
    function(s, t) bilinear_xy(s / Ltot, t / Wtot, P)
  } else {
    function(s, t) cbind(O[1] + k * ((s + off_s) * a[1] + (t + off_t) * cvec[1]),
                         O[2] + k * ((s + off_s) * a[2] + (t + off_t) * cvec[2]))
  }
  num <- plot_numbers(lay$n_ranges, lay$n_cols, opt$numbering, opt$sweep,
                      opt$start_corner, opt$first_plot)
  tab <- plot_table(lay$n_ranges, lay$n_cols, lay$L, lay$A, lay$W, lay$G, num)
  geom <- make_polys(tab, mapper, crs_out)
  inner <- if (lay$te > 0 || lay$ts > 0) make_polys(tab, mapper, crs_out, lay$te, lay$ts) else NULL

  ov <- overlap_check(geom)
  if (ov$area > 1e-6)
    stop(sprintf("Plots overlap (%d pairs, %.3f m2). Check --alley and --col-gap.", ov$pairs, ov$area),
         call. = FALSE)
  cen <- sf::st_coordinates(sf::st_centroid(geom))
  attrs <- data.frame(plot = tab$plot, range = tab$range, column = tab$column,
                      plot_id = tab$plot_id, length_m = lay$L, width_m = lay$W,
                      area_m2 = lay$L * lay$W, x = cen[, 1], y = cen[, 2])
  if (info$fit == "bilinear") {
    attrs$area_m2 <- as.numeric(sf::st_area(geom)) / k^2
    attrs$length_m <- NA_real_
    attrs$width_m <- NA_real_
  }
  if (!local) {
    ll <- sf::st_coordinates(sf::st_transform(sf::st_centroid(geom), 4326))
    attrs$lon <- ll[, 1]
    attrs$lat <- ll[, 2]
  }
  list(attrs = attrs, geom = geom, inner = inner, corner_df = corner_df, info = info,
       warnings = warn, overlap = ov, tab = tab, num = num)
}

# ---------------------------------------------------------------------------
# Map
# ---------------------------------------------------------------------------

nice_length <- function(x) {
  cand <- c(1, 2, 5) * rep(10^(-1:4), each = 3)
  cand[max(which(cand <= x))]
}

draw_map <- function(res, out_dir, title) {
  geom <- res$geom
  at <- res$attrs
  bb <- sf::st_bbox(geom)
  wd <- as.numeric(bb["xmax"] - bb["xmin"]); ht <- as.numeric(bb["ymax"] - bb["ymin"])
  fig_h <- max(4.5, min(11, 7.5 * ht / wd + 1.5))
  open_png(out_dir, "plot_map.png", width = 8, height = fig_h)
  on.exit(grDevices::dev.off())
  op <- graphics::par(mar = c(3.2, 3.2, 3.2, 1), mgp = c(2, 0.6, 0), cex.axis = 0.7)
  on.exit(graphics::par(op), add = TRUE, after = FALSE)
  fill <- if (!is.null(at$matched)) ifelse(at$matched, "#d9ead3", "#eeeeee") else "#d9ead3"
  xl <- c(bb["xmin"], bb["xmax"]) + c(-0.06, 0.06) * wd
  yl <- c(bb["ymin"], bb["ymax"]) + c(-0.08, 0.06) * ht
  graphics::plot(sf::st_geometry(geom), col = fill, border = "grey35", axes = TRUE,
                 xlim = xl, ylim = yl, main = title, cex.main = 0.85)
  if (!is.null(res$inner))
    graphics::plot(sf::st_geometry(res$inner), border = "#6aa84f", lty = 2, add = TRUE)
  o <- order(at$plot)
  graphics::lines(at$x[o], at$y[o], col = "#3d85c6", lwd = 0.8)
  step <- max(1L, floor(length(o) / 12))
  for (i in seq(1L, length(o) - 1L, by = step))
    graphics::arrows(at$x[o[i]], at$y[o[i]],
                     at$x[o[i]] + 0.45 * (at$x[o[i + 1]] - at$x[o[i]]),
                     at$y[o[i]] + 0.45 * (at$y[o[i + 1]] - at$y[o[i]]),
                     length = 0.07, col = "#3d85c6", lwd = 1.2)
  cex_lab <- max(0.3, min(0.8, 9 / sqrt(nrow(at))))
  graphics::text(at$x, at$y, labels = at$plot, cex = cex_lab, pos = 3, offset = 0.15)
  graphics::points(at$x[o[1]], at$y[o[1]], pch = 19, col = "#cc0000", cex = 1.2)
  # Corners of the block
  cn <- sf::st_coordinates(geom)
  pick <- function(r, cl) at$range == r & at$column == cl
  nr <- max(at$range); nc <- max(at$column)
  corner_plot <- c(which(pick(1, 1)), which(pick(nr, 1)), which(pick(nr, nc)), which(pick(1, nc)))
  ring <- function(i) sf::st_coordinates(geom[i])[1:4, 1:2, drop = FALSE]
  far <- function(i, ref) { r <- ring(i); r[which.max(rowSums((r - matrix(ref, 4, 2, byrow = TRUE))^2)), ] }
  ctr <- c(mean(at$x), mean(at$y))
  for (j in 1:4) {
    xy <- far(corner_plot[j], ctr)
    graphics::points(xy[1], xy[2], pch = 15, col = "black", cex = 0.8)
    graphics::text(xy[1], xy[2], labels = paste0("c", j), pos = 1, cex = 0.75, font = 2)
  }
  # North arrow (grid north) and scale bar
  ax <- xl[1] + 0.04 * diff(xl); ay <- yl[2] - 0.16 * diff(yl)
  graphics::arrows(ax, ay, ax, ay + 0.1 * diff(yl), length = 0.08, lwd = 2)
  graphics::text(ax, ay + 0.1 * diff(yl), "N (grid)", pos = 4, cex = 0.7)
  sl <- nice_length(0.25 * wd)
  sx <- xl[1] + 0.04 * diff(xl); sy <- yl[1] + 0.04 * diff(yl)
  graphics::segments(sx, sy, sx + sl, sy, lwd = 3)
  graphics::text(sx + sl / 2, sy, paste(sl, "m"), pos = 3, cex = 0.7)
  leg <- data.frame(label = c(sprintf("plot %d (start)", at$plot[o[1]]), "numbering path"),
                    col = c("#cc0000", "#3d85c6"), pch = c(19, NA), lty = c(NA, 1), stringsAsFactors = FALSE)
  if (!is.null(res$inner))
    leg <- rbind(leg, data.frame(label = "trimmed (inner) plot", col = "#6aa84f", pch = NA, lty = 2))
  if (!is.null(at$matched) && any(!at$matched))
    leg <- rbind(leg, data.frame(label = "no field-book row", col = "grey60", pch = 15, lty = NA))
  graphics::legend("bottomright", bty = "n", cex = 0.7, legend = leg$label, col = leg$col,
                   pch = leg$pch, lty = leg$lty)
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Report and main
# ---------------------------------------------------------------------------

write_gpkg <- function(attrs, geom, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  obj <- sf::st_sf(attrs, geometry = geom)
  if (file.exists(path)) unlink(path)
  sf::st_write(obj, path, layer = "plots", quiet = TRUE)
  path
}

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)

write_report <- function(res, opt, fb_info, out_dir, warn) {
  inf <- res$info
  lay <- inf$lay
  u <- opt$units
  lines <- c(sprintf("# %s report", SKILL), "")
  if (length(warn)) lines <- c(lines, "## Warnings", "", paste("-", warn), "")
  lines <- c(lines, "## Geometry", "",
             sprintf("- Input: %s mode, coordinates in %s.", inf$mode, opt$crs),
             sprintf("- Output CRS: %s.", inf$out_label))
  if (!inf$local) {
    lines <- c(lines,
               sprintf("- Field reference point: lon %.7f, lat %.7f ([map](https://www.openstreetmap.org/?mlat=%.7f&mlon=%.7f#map=18/%.7f/%.7f)).",
                       inf$ref_ll[1], inf$ref_ll[2], inf$ref_ll[2], inf$ref_ll[1], inf$ref_ll[2], inf$ref_ll[1]),
               sprintf("- Planting direction (range 1 to range %d): %.3f deg from true north = %.3f deg from grid north (bearing source: %s; grid north lies %.3f deg %s of true north here).",
                       lay$n_ranges, inf$bearing_true, inf$bearing_grid, inf$bearing_ref,
                       abs(inf$gf$north_bearing), if (inf$gf$north_bearing <= 0) "east" else "west"),
               sprintf("- Columns increase to the %s of the planting direction.", inf$across),
               sprintf("- Map scale factor applied to ground lengths: k = %.7f (%s). Polygon areas in map units are k^2 times the ground areas reported here.",
                       inf$k, if (tolower(opt$scale_factor) == "auto") "projection point scale at the field" else "user value"))
  } else {
    lines <- c(lines, sprintf("- Planting direction: %.3f deg from local grid north; columns to the %s.",
                              inf$bearing_grid, inf$across))
  }
  lines <- c(lines,
             sprintf("- Layout: %d ranges x %d columns = %d plots; plot %.3f m long x %.3f m wide; alley %.3f m; column gap %.3f m (input in %s).",
                     lay$n_ranges, lay$n_cols, lay$n_ranges * lay$n_cols, lay$L, lay$W, lay$A, lay$G, u),
             sprintf("- Plot block: %.3f m along x %.3f m across.", inf$Ltot, inf$Wtot),
             sprintf("- Numbering: %s, sweep %s, starting at corner c%d with plot %d.",
                     opt$numbering, opt$sweep, as.integer(opt$start_corner), as.integer(opt$first_plot)),
             if (!is.null(res$inner)) sprintf("- Inner (trimmed) plots: %.3f m off each end and %.3f m off each side, leaving %.3f x %.3f m.",
                                              lay$te, lay$ts, lay$L - 2 * lay$te, lay$W - 2 * lay$ts) else NULL,
             "")
  if (!is.null(res$corner_df)) {
    cd <- res$corner_df
    lines <- c(lines, "## Corner check", "",
               if (inf$fit == "bilinear") sprintf("Layout stretched to the corners (bilinear). Implied scale: %.4f along, %.4f across (1 = layout dimensions match the corners exactly).",
                                                inf$fit_scale[1], inf$fit_scale[2]) else
                 "Rigid least-squares fit of the layout to the 4 corners (rotation and shift only).",
               "", "| corner | given x | given y | fitted x | fitted y | misfit (m) |", "|---|---|---|---|---|---|",
               sprintf("| %s | %s | %s | %s | %s | %s |", cd$corner, fmt(cd$x_given), fmt(cd$y_given),
                       fmt(cd$x_fitted), fmt(cd$y_fitted), fmt(cd$misfit_m)),
               "",
               sprintf("Corner sides (ground m): c1-c2 %.3f, c4-c3 %.3f, c1-c4 %.3f, c2-c3 %.3f; angle at c1 %.2f deg.",
                       inf$side_lengths[1], inf$side_lengths[2], inf$side_lengths[3], inf$side_lengths[4],
                       inf$corner_angle), "")
  }
  lines <- c(lines, "## Checks", "",
             sprintf("- Plot interiors overlapping: %d pairs, %.6f m2 (must be 0).", res$overlap$pairs, res$overlap$area),
             sprintf("- Plot ground area: %s m2 each.", if (inf$fit == "bilinear")
               paste(fmt(range(res$attrs$area_m2), 3), collapse = " to ") else fmt(lay$L * lay$W, 3)),
             "")
  if (!is.null(fb_info)) {
    lines <- c(lines, "## Field book", "",
               sprintf("- File: `%s`; joined on %s.", opt$fieldbook, opt$join_by),
               sprintf("- Plots with a field-book row: %d of %d.", sum(res$attrs$matched), nrow(res$attrs)),
               if (length(fb_info$renamed)) sprintf("- Field-book columns renamed with prefix fb_: %s.",
                                                     paste(fb_info$renamed, collapse = ", ")) else NULL,
               if (!is.null(fb_info$coords)) sprintf("- Position cross-check (%s): consistent (%s).",
                                                     opt$fb_coords, fb_info$coords$mapping) else NULL,
               "")
  }
  lines <- c(lines, "## Map", "", "![Plot map](figures/plot_map.png)", "",
             "Check the map against the planting plan: plot 1 position, numbering direction and the side the columns grow towards.", "",
             "## Outputs", "",
             "- `spatial/plots.gpkg` (layer `plots`): plot polygons with plot, range, column, plot_id and field-book columns.",
             if (!is.null(res$inner)) "- `spatial/plots_inner.gpkg` (layer `plots`): trimmed polygons for extraction." else NULL,
             "- `tables/plots.csv`: the same attributes with centroid coordinates.",
             if (!is.null(res$corner_df)) "- `tables/corner_check.csv`: corner misfits." else NULL,
             "", "## Limitations", "",
             "- Plots are rectangles of the nominal size; real plots can drift from the planter path. Overlay the polygons on the orthomosaic before extracting values.",
             "- Lengths are ground distances on the WGS84 ellipsoid; the height of the field above the ellipsoid is ignored (error below 0.01% under about 600 m).",
             "", "## References", "",
             "- Matias FI, Caraza-Harter MV, Endelman JB (2020) FIELDimageR: an R package to analyze orthomosaic images from agricultural field trials. The Plant Phenome Journal 3:e20005. https://doi.org/10.1002/ppj2.20005",
             "- Snyder JP (1987) Map projections: a working manual. USGS Professional Paper 1395. https://doi.org/10.3133/pp1395",
             "- Pebesma E (2018) Simple features for R: standardized support for spatial vector data. The R Journal 10:439-446. https://doi.org/10.32614/RJ-2018-009",
             "", sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(lines, file.path(out_dir, "report.md"))
}

demo_options <- function(opt) {
  d <- file.path(SKILL_DIR, "demo")
  opt$corners <- file.path(d, "corners.csv")
  opt$crs <- "EPSG:4326"
  opt$ranges <- 6; opt$columns <- 8
  opt$plot_length <- 5.2; opt$alley <- 0.9
  opt$rows_per_plot <- 2; opt$row_spacing <- 0.76
  opt$units <- "m"
  opt$numbering <- "serpentine"; opt$sweep <- "along"; opt$start_corner <- 1
  opt$trim_ends <- 0.5; opt$trim_sides <- 0.15
  opt$fieldbook <- file.path(d, "fieldbook.csv")
  opt$join_by <- "plot"; opt$fb_coords <- "row,col"
  opt
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  given <- given_options(argv)
  if (opt$demo) opt <- demo_options(opt)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  need_pkgs("sf")

  res <- build_grid(opt)
  warn <- res$warnings
  if (!opt$demo) {
    dflt <- setdiff(c("numbering", "sweep", "start_corner"), given)
    if (length(dflt))
      warn <- c(warn, sprintf(paste("Plot numbering uses default value(s) for %s (serpentine, along a",
                                    "planter pass, plot 1 at corner c1). Confirm against the planting",
                                    "map before joining data."), paste0("--", gsub("_", "-", dflt), collapse = ", ")))
    if (res$info$mode == "origin" && !"across" %in% given)
      warn <- c(warn, "--across not given: columns assumed to increase to the right of the planting direction.")
    if (res$info$mode == "origin" && !"bearing_ref" %in% given)
      warn <- c(warn, sprintf("--bearing-ref not given: the bearing was taken as %s north.", res$info$bearing_ref))
  }

  attrs <- res$attrs
  fb_info <- NULL
  if (nzchar(opt$fieldbook)) {
    if (!file.exists(opt$fieldbook)) stop("Field book not found: ", opt$fieldbook, call. = FALSE)
    fb <- utils::read.csv(opt$fieldbook, check.names = FALSE, stringsAsFactors = FALSE)
    j <- join_fieldbook(attrs, fb, opt$join_by)
    attrs <- j$data
    fb_info <- list(renamed = j$renamed)
    if (j$n_unmatched > 0)
      warn <- c(warn, sprintf("%d plot(s) have no field-book row (fillers, borders or a layout mismatch?).",
                              j$n_unmatched))
    if (nzchar(opt$fb_coords)) {
      fc <- trimws(strsplit(opt$fb_coords, ",", fixed = TRUE)[[1]])
      fc_names <- ifelse(tolower(fc) %in% RESERVED, paste0("fb_", fc), fc)
      if (length(fc) != 2L || !all(fc_names %in% names(attrs)))
        stop("--fb-coords needs two field-book column names, e.g. 'ROW,COLUMN'.", call. = FALSE)
      chk <- check_fb_coords(attrs, fc_names)
      if (!chk$ok)
        stop("The field book's position columns (", opt$fb_coords, ") do not map one-to-one onto ",
             "the grid's ranges and columns for the joined plots: the field book numbers plots with a ",
             "different scheme. Adjust --numbering, --sweep and --start-corner until they agree.",
             call. = FALSE)
      fb_info$coords <- chk
    }
  }
  res$attrs <- attrs

  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  write_gpkg(attrs, res$geom, file.path(opt$out, "spatial", "plots.gpkg"))
  if (!is.null(res$inner)) {
    ia <- attrs
    ia$length_m <- ia$length_m - 2 * res$info$lay$te
    ia$width_m <- ia$width_m - 2 * res$info$lay$ts
    ia$area_m2 <- as.numeric(sf::st_area(res$inner)) / res$info$k^2
    write_gpkg(ia, res$inner, file.path(opt$out, "spatial", "plots_inner.gpkg"))
  }
  save_table(attrs, opt$out, "plots.csv")
  if (!is.null(res$corner_df)) save_table(res$corner_df, opt$out, "corner_check.csv")
  draw_map(res, opt$out, sprintf("%d plots, %s numbering (%s), start c%d", nrow(attrs),
                                 opt$numbering, opt$sweep, as.integer(opt$start_corner)))
  write_report(res, opt, fb_info, opt$out, warn)
  inputs <- c(if (grepl("\\.csv$", opt$corners, ignore.case = TRUE)) opt$corners else "", opt$fieldbook)
  write_repro_bundle(opt$out, argv, inputs = inputs)
  cat(sprintf("Done: %d plots. Report: %s\n", nrow(attrs),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  if (length(warn)) cat(paste0("Warning: ", warn, "\n"), sep = "")
  invisible(res)
}

if (sys.nframe() == 0L) main()
