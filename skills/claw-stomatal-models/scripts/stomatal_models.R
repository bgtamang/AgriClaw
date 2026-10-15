#!/usr/bin/env Rscript
# claw-stomatal-models: fit the Ball-Berry (gsw = g0 + g1 A hs / Cs) and
# Medlyn USO (gsw = g0 + 1.6 (1 + g1 / sqrt(D)) A / Cs) stomatal models per
# curve with PhotoGEA (fit_ball_berry, fit_medlyn), computing leaf-surface
# humidity, CO2 and VPD from LI-6800 gas exchange with
# PhotoGEA::calculate_gas_properties, and report intrinsic water-use
# efficiency (A / gsw). Implausible units stop the run with an explanation.
#
# Input: the tidy CSV written by claw-licor-6800 (or any CSV with the same
# column names and units, or with precomputed leaf-surface columns).
#
# Usage:
#   Rscript scripts/stomatal_models.R --demo
#   Rscript scripts/stomatal_models.R --input licor_tidy.csv --curve-cols genotype,plot
#   Rscript scripts/stomatal_models.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT license.

SKILL   <- "claw-stomatal-models"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/stomatal_models.R [--demo] [--input FILE.csv] [--curve-cols COLS]
                                       [--models ball-berry,medlyn] [--vpd-source surface|licor]
                                       [--keep-unstable] [--include-flagged] [--min-points N]
                                       [--out DIR]

  --demo             fit the bundled synthetic curves (demo/stomatal_demo.csv)
  --input            tidy CSV from claw-licor-6800 (one row per logged point)
  --curve-cols       column(s) that together identify one curve, comma-separated
  --models           ball-berry, medlyn or both (default: ball-berry,medlyn)
  --vpd-source       D for the Medlyn model: surface (default; leaf-surface VPD from
                     calculate_gas_properties) or licor (LI-6800 VPDleaf, leaf to chamber air)
  --keep-unstable    keep points logged before stability (default: excluded, because
                     stomata are slow to reach steady state)
  --include-flagged  keep points with qc_hard_fail = TRUE (not advised)
  --min-points       minimum points per curve for a fit (default 4; rule of thumb)
  --out              output folder (default: agriclaw_output/claw-stomatal-models)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, input = "", curve_cols = "", models = "ball-berry,medlyn",
                 vpd_source = "surface", keep_unstable = FALSE, include_flagged = FALSE,
                 min_points = 4, out = "")

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

U_A <- "micromol m^(-2) s^(-1)"
U_G <- "mol m^(-2) s^(-1)"
RAW_COLS <- c("e_mol_m2_s", "gbw_mol_m2_s", "h2o_s_mmol_mol", "ca_umol_mol", "tleaf_c", "pa_kpa")
SURFACE_COLS <- c("rh_leaf_pct", "cs_umol_mol", "vpd_leaf_surface_kpa")

# Rules of thumb for plausible values (see SKILL.md).
PLAUSIBLE <- list(gsw_max = 3, vpd_max = 10, vpd_min = 0.02, cs_min = 50, cs_max = 5000,
                  bb_g1 = c(0, 30), medlyn_g1 = c(0, 20))

split_cols <- function(s) {
  if (!nzchar(s)) return(character())
  trimws(strsplit(s, ",", fixed = TRUE)[[1L]])
}

read_tidy <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
  if (tolower(tools::file_ext(path)) %in% c("xlsx", "", "txt"))
    stop("'", basename(path), "' looks like a raw LI-6800 log. Run claw-licor-6800 first and ",
         "pass its tables/licor_tidy.csv here.", call. = FALSE)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

# Unit plausibility on the inputs; stops with the likely cause.
check_input_units <- function(df) {
  g <- df$gsw_mol_m2_s
  if (max(g, na.rm = TRUE) > PLAUSIBLE$gsw_max)
    stop("gsw reaches ", signif(max(g, na.rm = TRUE), 3), " mol m-2 s-1; leaf values above ",
         PLAUSIBLE$gsw_max, " are implausible. gsw looks like mmol m-2 s-1: divide by 1000.", call. = FALSE)
  if (max(df$a_umol_m2_s, na.rm = TRUE) > 100)
    stop("A exceeds 100 umol m-2 s-1: check the units of a_umol_m2_s.", call. = FALSE)
  if ("rh_leaf_pct" %in% names(df)) {
    rh <- df$rh_leaf_pct
    if (max(rh, na.rm = TRUE) <= 1)
      stop("rh_leaf_pct is at most 1: it looks like a fraction. Give relative humidity in % (0-100).", call. = FALSE)
    if (any(rh < 0 | rh > 100, na.rm = TRUE)) stop("rh_leaf_pct outside 0-100%.", call. = FALSE)
  }
  if ("cs_umol_mol" %in% names(df)) check_cs(df$cs_umol_mol, "cs_umol_mol")
  if ("vpd_leaf_surface_kpa" %in% names(df)) check_vpd(df$vpd_leaf_surface_kpa, "vpd_leaf_surface_kpa")
  invisible(TRUE)
}

check_cs <- function(cs, what) {
  med <- stats::median(cs, na.rm = TRUE)
  if (med < PLAUSIBLE$cs_min)
    stop(what, " has median ", signif(med, 3), ": CO2 at the leaf surface should be in umol/mol ",
         "(about 400 in air). Values near 40 are Pa; convert with umol/mol = Pa / kPa total pressure * 1000.",
         call. = FALSE)
  if (med > PLAUSIBLE$cs_max) stop(what, " has median ", signif(med, 3), " umol/mol: check units.", call. = FALSE)
}

check_vpd <- function(vpd, what) {
  mx <- max(vpd, na.rm = TRUE)
  if (mx > PLAUSIBLE$vpd_max)
    stop(what, " reaches ", signif(mx, 4), ". Leaf VPD cannot exceed the saturation vapor pressure ",
         "(about 9.6 kPa at 45 C), so D is not in kPa: values in the hundreds or thousands are Pa ",
         "(divide by 1000), tens are hPa or mbar (divide by 10).", call. = FALSE)
  if (stats::median(vpd, na.rm = TRUE) < PLAUSIBLE$vpd_min)
    stop(what, " has median ", signif(stats::median(vpd, na.rm = TRUE), 3), " kPa: this looks like MPa ",
         "or a fraction. Give D in kPa.", call. = FALSE)
  if (any(vpd <= 0, na.rm = TRUE)) stop(what, " has zero or negative values; the Medlyn model needs D > 0.", call. = FALSE)
}

validate_input <- function(df, curve_cols, opt) {
  miss <- setdiff(c("a_umol_m2_s", "gsw_mol_m2_s"), names(df))
  if (length(miss))
    stop("Missing column(s): ", paste(miss, collapse = ", "), ". The input must be the tidy CSV from ",
         "claw-licor-6800 (A in umol m-2 s-1, gsw in mol m-2 s-1).", call. = FALSE)
  if (!length(curve_cols))
    stop("--curve-cols is required: which column(s) identify one stomatal-response curve or leaf?", call. = FALSE)
  miss <- setdiff(curve_cols, names(df))
  if (length(miss)) stop("--curve-cols column(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
  has_surface <- all(SURFACE_COLS[1:2] %in% names(df))
  has_raw <- all(RAW_COLS %in% names(df))
  if (!has_surface && !has_raw)
    stop("Need either the LI-6800 columns ", paste(RAW_COLS, collapse = ", "), " (to compute leaf-surface ",
         "humidity and CO2), or precomputed ", paste(SURFACE_COLS, collapse = ", "), ".", call. = FALSE)
  if (opt$vpd_source == "licor" && !"vpd_leaf_kpa" %in% names(df))
    stop("--vpd-source licor needs the vpd_leaf_kpa column (LI-6800 VPDleaf).", call. = FALSE)
  if (opt$vpd_source == "surface" && has_surface && !has_raw && !"vpd_leaf_surface_kpa" %in% names(df) &&
      "medlyn" %in% split_cols(opt$models))
    stop("The Medlyn model needs vpd_leaf_surface_kpa (or the raw LI-6800 columns to compute it).", call. = FALSE)
  check_input_units(df)
  invisible(has_surface && !has_raw)
}

curve_key <- function(df, curve_cols) {
  do.call(paste, c(lapply(curve_cols, function(cn) as.character(df[[cn]])), sep = " - "))
}

# Build a PhotoGEA exdf with leaf-surface RH (RHleaf, %), CO2 (Csurface),
# leaf-surface VPD (VPDsurf, kPa) and the Ball-Berry index.
build_exdf <- function(df, opt, precomputed) {
  base <- data.frame(curve_id = df$curve_id, A = df$a_umol_m2_s, gsw = df$gsw_mol_m2_s,
                     include_when_fitting = df$used_in_fit, stringsAsFactors = FALSE)
  if (precomputed) {
    base$RHleaf <- df$rh_leaf_pct
    base$Csurface <- df$cs_umol_mol
    if ("vpd_leaf_surface_kpa" %in% names(df)) base$VPDsurf <- df$vpd_leaf_surface_kpa
    x <- PhotoGEA::exdf(base, units = data.frame(A = U_A, gsw = U_G, RHleaf = "%", Csurface = "micromol mol^(-1)",
                                                 VPDsurf = "kPa")[, intersect(c("A", "gsw", "RHleaf", "Csurface", "VPDsurf"), names(base))])
  } else {
    base$E <- df$e_mol_m2_s
    base$gbw <- df$gbw_mol_m2_s
    base$H2O_s <- df$h2o_s_mmol_mol
    base$Ca <- df$ca_umol_mol
    base$TleafCnd <- df$tleaf_c
    base$Pa <- df$pa_kpa
    base$DeltaPcham <- if ("delta_pcham_kpa" %in% names(df)) df$delta_pcham_kpa else 0
    x <- PhotoGEA::exdf(base, units = data.frame(A = U_A, gsw = U_G, E = U_G, gbw = U_G, H2O_s = "mmol mol^(-1)",
                                                 Ca = "micromol mol^(-1)", TleafCnd = "degrees C", Pa = "kPa",
                                                 DeltaPcham = "kPa"))
    x <- PhotoGEA::calculate_total_pressure(x)
    x <- PhotoGEA::calculate_gas_properties(x)
    x[, "VPDsurf"] <- x[, "SVPleaf"] * (1 - x[, "RHleaf"] / 100)
    x <- PhotoGEA::document_variables(x, c("claw-stomatal-models", "VPDsurf", "kPa"))
    rh <- x[, "RHleaf"]
    if (any(rh < 0 | rh > 100, na.rm = TRUE))
      stop("Leaf-surface relative humidity computed from E, gbw and H2O_s falls outside 0-100%: check the ",
           "units of e_mol_m2_s (mol, not mmol), gbw_mol_m2_s and h2o_s_mmol_mol.", call. = FALSE)
    check_cs(x[, "Csurface"], "Leaf-surface CO2 computed from Ca")
  }
  if (opt$vpd_source == "licor") {
    x[, "VPDuse"] <- df$vpd_leaf_kpa
  } else if ("VPDsurf" %in% colnames(x)) {
    x[, "VPDuse"] <- x[, "VPDsurf"]
  }
  if ("VPDuse" %in% colnames(x)) {
    x <- PhotoGEA::document_variables(x, c("claw-stomatal-models", "VPDuse", "kPa"))
    check_vpd(x[, "VPDuse"][df$used_in_fit], if (opt$vpd_source == "licor") "vpd_leaf_kpa" else "Leaf-surface VPD")
  }
  PhotoGEA::calculate_ball_berry_index(x)
}

r_squared <- function(y, f) {
  ok <- is.finite(y) & is.finite(f)
  1 - sum((y[ok] - f[ok])^2) / sum((y[ok] - mean(y[ok]))^2)
}

fit_one <- function(xc, models, min_points) {
  use <- xc[, "include_when_fitting"]
  n <- sum(use)
  out <- data.frame(curve_id = xc[1, "curve_id"], n_points_used = n, stringsAsFactors = FALSE)
  pts <- data.frame(curve_id = xc[, "curve_id"], used_in_fit = use, a_umol_m2_s = xc[, "A"],
                    gsw_mol_m2_s = xc[, "gsw"], rh_leaf_pct = xc[, "RHleaf"],
                    cs_umol_mol = xc[, "Csurface"], bb_index_mol_m2_s = xc[, "bb_index"],
                    stringsAsFactors = FALSE)
  if ("VPDuse" %in% colnames(xc)) {
    pts$vpd_kpa <- xc[, "VPDuse"]
    pts$medlyn_index <- xc[, "A"] / (xc[, "Csurface"] * sqrt(xc[, "VPDuse"]))
  }
  pts$iwue_umol_mol <- xc[, "A"] / xc[, "gsw"]
  errs <- character()
  if (n < min_points) {
    out$fit_error <- sprintf("only %d usable points (need %d)", n, min_points)
    return(list(params = out, points = pts))
  }
  if ("ball-berry" %in% models) {
    bb <- tryCatch(PhotoGEA::fit_ball_berry(xc), error = function(e) e)
    if (inherits(bb, "error")) errs <- c(errs, paste("Ball-Berry:", conditionMessage(bb))) else {
      p <- bb$parameters
      out$bb_g0_mol_m2_s <- p[, "bb_intercept"]
      out$bb_g0_se <- p[, "bb_intercept_err"]
      out$bb_g1 <- p[, "bb_slope"]
      out$bb_g1_se <- p[, "bb_slope_err"]
      out$bb_r2 <- p[, "r_squared"]
      out$bb_p_value <- p[, "p_value"]
      out$bb_rmse_mol_m2_s <- p[, "RMSE"]
      pts$gsw_fit_bb <- bb$fits[, "gsw_fit"]
      ix <- xc[use, "bb_index"]
      out$bb_index_cv <- stats::sd(ix) / mean(ix)
    }
  }
  if ("medlyn" %in% models && "VPDuse" %in% colnames(xc)) {
    md <- tryCatch(PhotoGEA::fit_medlyn(xc, vpdleaf_column_name = "VPDuse"), error = function(e) e)
    if (inherits(md, "error")) errs <- c(errs, paste("Medlyn:", conditionMessage(md))) else {
      p <- md$parameters
      out$medlyn_g0_mol_m2_s <- p[, "medlyn_g0"]
      out$medlyn_g0_se <- p[, "medlyn_g0_err"]
      out$medlyn_g1_kpa05 <- p[, "medlyn_g1"]
      out$medlyn_g1_se <- p[, "medlyn_g1_err"]
      out$medlyn_rmse_mol_m2_s <- p[, "RMSE"]
      pts$gsw_fit_medlyn <- md$fits[, "gsw_fit"]
      out$medlyn_r2 <- r_squared(xc[use, "gsw"], pts$gsw_fit_medlyn[use])
    }
  }
  out$iwue_mean_umol_mol <- mean(pts$iwue_umol_mol[use])
  out$iwue_sd_umol_mol <- stats::sd(pts$iwue_umol_mol[use])
  out$a_mean_umol_m2_s <- mean(xc[use, "A"])
  out$gsw_mean_mol_m2_s <- mean(xc[use, "gsw"])
  if ("VPDuse" %in% colnames(xc)) out$vpd_mean_kpa <- mean(xc[use, "VPDuse"])
  out$fit_error <- paste(errs, collapse = "; ")
  list(params = out, points = pts)
}

stomatal_flags <- function(p) {
  f <- character()
  g <- function(cn) if (cn %in% names(p)) p[[cn]] else NA
  if (nzchar(g("fit_error") %||% "")) f <- c(f, "fit_failed")
  if (!is.na(g("bb_g1")) && (g("bb_g1") < PLAUSIBLE$bb_g1[1] || g("bb_g1") > PLAUSIBLE$bb_g1[2]))
    f <- c(f, "bb_g1_implausible")
  if (!is.na(g("medlyn_g1_kpa05")) && (g("medlyn_g1_kpa05") < PLAUSIBLE$medlyn_g1[1] ||
                                       g("medlyn_g1_kpa05") > PLAUSIBLE$medlyn_g1[2]))
    f <- c(f, "medlyn_g1_implausible")
  if ((!is.na(g("bb_g0_mol_m2_s")) && g("bb_g0_mol_m2_s") < 0) ||
      (!is.na(g("medlyn_g0_mol_m2_s")) && g("medlyn_g0_mol_m2_s") < 0)) f <- c(f, "negative_g0")
  if (!is.na(g("bb_index_cv")) && g("bb_index_cv") < 0.1) f <- c(f, "narrow_index_range")
  if (!is.na(g("bb_r2")) && g("bb_r2") < 0.5) f <- c(f, "low_r2")
  paste(f, collapse = ";")
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || (length(a) == 1L && is.na(a))) b else a

rbind_fill <- function(dfs) {
  cols <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function(d) {
    for (cn in setdiff(cols, names(d))) d[[cn]] <- NA
    d[, cols, drop = FALSE]
  })
  out <- do.call(rbind, dfs)
  rownames(out) <- NULL
  out
}

curve_metadata <- function(points, curve_cols) {
  cand <- setdiff(names(points), c("curve_id", "used_in_fit", "excluded_reason"))
  cand <- cand[!startsWith(cand, "qc_")]
  cand <- cand[!vapply(points[cand], function(v) is.numeric(v) || is.logical(v), logical(1)) |
                 cand %in% curve_cols]
  keep <- cand[vapply(cand, function(cn) all(tapply(points[[cn]], points$curve_id,
                                                    function(v) length(unique(v)) == 1L)), logical(1))]
  keep <- setdiff(keep, c("date_local", "hhmmss", "user_remark", "file_type"))
  m <- points[!duplicated(points$curve_id), c("curve_id", keep), drop = FALSE]
  rownames(m) <- NULL
  m
}

analyze_stomatal <- function(df, opt) {
  curve_cols <- split_cols(opt$curve_cols)
  models <- tolower(split_cols(opt$models))
  if (!length(models) || !all(models %in% c("ball-berry", "medlyn")))
    stop("--models must be ball-berry, medlyn or both.", call. = FALSE)
  if (!opt$vpd_source %in% c("surface", "licor")) stop("--vpd-source must be surface or licor.", call. = FALSE)
  precomputed <- validate_input(df, curve_cols, opt)
  df$curve_id <- curve_key(df, curve_cols)
  reason <- rep("", nrow(df))
  add <- function(cond, why) {
    cond <- cond %in% TRUE
    reason[cond] <<- ifelse(nzchar(reason[cond]), paste(reason[cond], why, sep = ";"), why)
  }
  if (!opt$include_flagged && "qc_hard_fail" %in% names(df)) add(df$qc_hard_fail, "qc_hard_fail")
  if (!opt$keep_unstable && "qc_unstable" %in% names(df)) add(df$qc_unstable, "unstable")
  add(is.na(df$a_umol_m2_s) | is.na(df$gsw_mol_m2_s), "missing_A_or_gsw")
  # Missing inputs to hs, Cs or D would otherwise be dropped silently by lm/nls.
  req <- if (precomputed) intersect(SURFACE_COLS, names(df)) else intersect(c(RAW_COLS, "delta_pcham_kpa"), names(df))
  if (opt$vpd_source == "licor") req <- c(req, "vpd_leaf_kpa")
  if (length(req)) add(Reduce(`|`, lapply(req, function(cn) is.na(df[[cn]]))), "missing_inputs")
  df$excluded_reason <- reason
  df$used_in_fit <- !nzchar(reason)
  x <- build_exdf(df, opt, precomputed)
  ids <- sort(unique(df$curve_id))
  fits <- lapply(ids, function(id) fit_one(x[x[, "curve_id"] == id, , TRUE], models, opt$min_points))
  params <- rbind_fill(lapply(fits, `[[`, "params"))
  params$qc_flags <- vapply(seq_len(nrow(params)), function(i) stomatal_flags(params[i, , drop = FALSE]), "")
  params$vpd_source <- if ("medlyn" %in% models) opt$vpd_source else NA_character_
  meta <- curve_metadata(df, curve_cols)
  params <- merge(meta, params, by = "curve_id", all.y = TRUE, sort = FALSE)
  params <- params[order(params$curve_id), , drop = FALSE]
  rownames(params) <- NULL
  pts <- rbind_fill(lapply(fits, `[[`, "points"))
  vpd_diff <- NA_real_
  if (!precomputed && "vpd_leaf_kpa" %in% names(df))
    vpd_diff <- mean(df$vpd_leaf_kpa[df$used_in_fit] - x[, "VPDsurf"][df$used_in_fit], na.rm = TRUE)
  list(params = params, points = pts, models = models, vpd_diff = vpd_diff, data = df)
}

plot_fits <- function(res, out_dir) {
  pts <- res$points
  ids <- sort(unique(pts$curve_id))
  cols <- grDevices::hcl.colors(max(2L, length(ids)), "Dark 3")
  nm <- "stomatal_fits.png"
  n_pan <- length(res$models)
  open_png(out_dir, nm, width = 5 * n_pan, height = 4.6)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(1, n_pan), mar = c(4.2, 4.4, 2.5, 1))
  draw <- function(xcol, fcol, xlab, main) {
    u <- pts$used_in_fit
    graphics::plot(pts[[xcol]], pts$gsw_mol_m2_s, type = "n", xlab = xlab, ylab = "gsw (mol m-2 s-1)", main = main)
    for (k in seq_along(ids)) {
      d <- pts[pts$curve_id == ids[k], , drop = FALSE]
      graphics::points(d[[xcol]][d$used_in_fit], d$gsw_mol_m2_s[d$used_in_fit], pch = 16, col = cols[k])
      graphics::points(d[[xcol]][!d$used_in_fit], d$gsw_mol_m2_s[!d$used_in_fit], pch = 1, col = cols[k])
      if (fcol %in% names(d)) {
        o <- order(d[[xcol]])
        graphics::lines(d[[xcol]][o], d[[fcol]][o], col = cols[k])
      }
    }
    graphics::legend("topleft", legend = substr(ids, 1L, 25L), col = cols[seq_along(ids)], pch = 16, bty = "n", cex = 0.65)
  }
  if ("ball-berry" %in% res$models)
    draw("bb_index_mol_m2_s", "gsw_fit_bb", "Ball-Berry index A hs / Cs (mol m-2 s-1)", "Ball-Berry")
  if ("medlyn" %in% res$models && "medlyn_index" %in% names(pts))
    draw("medlyn_index", "gsw_fit_medlyn", "Medlyn index A / (Cs sqrt(D)) (mol m-2 s-1 kPa^-0.5)", "Medlyn")
  nm
}

fmt_table <- function(df, digits = 3) {
  if (is.null(df) || !nrow(df)) return(character())
  f <- vapply(df, function(v) if (is.numeric(v)) formatC(v, digits = digits, format = "g") else as.character(v),
              character(nrow(df)))
  if (is.null(dim(f))) f <- matrix(f, nrow = 1L)
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(f, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

FLAG_TEXT <- c(
  fit_failed = "a fit failed (see fit_error)",
  bb_g1_implausible = "Ball-Berry g1 outside 0-30: check the units of A, gsw, hs and Cs",
  medlyn_g1_implausible = "Medlyn g1 outside 0-20 kPa^0.5: check that D is in kPa and Cs in umol/mol",
  negative_g0 = "negative g0: the intercept is extrapolated and poorly constrained; do not interpret it as night conductance",
  narrow_index_range = "the Ball-Berry index varies by less than 10% (CV): the slope is poorly constrained",
  low_r2 = "Ball-Berry R2 below 0.5")

write_report <- function(res, opt, input, out_dir, fig) {
  p <- res$params
  pts <- res$data
  w <- character()
  n_un <- sum(grepl("unstable", pts$excluded_reason))
  if (n_un) w <- c(w, sprintf("%d point(s) logged before stability were excluded (stomata are slow to reach steady state; --keep-unstable keeps them).", n_un))
  for (fl in names(FLAG_TEXT)) {
    hit <- p$curve_id[grepl(paste0("(^|;)", fl, "($|;)"), p$qc_flags)]
    if (length(hit)) w <- c(w, sprintf("%s%s (%d curve(s): %s).", toupper(substr(FLAG_TEXT[[fl]], 1, 1)),
                                       substring(FLAG_TEXT[[fl]], 2), length(hit), paste(hit, collapse = "; ")))
  }
  if (!is.na(res$vpd_diff))
    w <- c(w, sprintf("The LI-6800 VPDleaf (leaf to chamber air) exceeds the leaf-surface VPD by %.2f kPa on average; Medlyn g1 used the %s value.",
                      res$vpd_diff, if (opt$vpd_source == "licor") "LI-6800" else "leaf-surface"))
  cols <- intersect(c("curve_id", "n_points_used", "bb_g0_mol_m2_s", "bb_g1", "bb_r2", "medlyn_g0_mol_m2_s",
                      "medlyn_g1_kpa05", "medlyn_r2", "iwue_mean_umol_mol", "qc_flags"), names(p))
  show <- p[, cols]
  lab <- c(curve_id = "curve", n_points_used = "n", bb_g0_mol_m2_s = "BB g0", bb_g1 = "BB g1", bb_r2 = "BB R2",
           medlyn_g0_mol_m2_s = "Medlyn g0", medlyn_g1_kpa05 = "Medlyn g1", medlyn_r2 = "Medlyn R2",
           iwue_mean_umol_mol = "iWUE", qc_flags = "flags")
  names(show) <- lab[cols]
  lines <- c(sprintf("# %s report", SKILL), "", "## Warnings", "",
             if (length(w)) paste("-", w) else "None.", "",
             "## Results", "",
             sprintf("Input: `%s`. %d curve(s); %d of %d points used.", input, nrow(p), sum(pts$used_in_fit), nrow(pts)), "",
             "Units: g0 in mol m-2 s-1 (water vapor); Ball-Berry g1 dimensionless (hs as a fraction, Cs in umol/mol); Medlyn g1 in kPa^0.5; iWUE = A / gsw in umol CO2 per mol H2O. Standard errors are in `tables/stomatal_parameters.csv`.", "",
             fmt_table(show), "",
             paste0("![Stomatal model fits](figures/", fig, ")"), "",
             "## Method", "",
             paste0("Leaf-surface relative humidity (hs), CO2 (Cs) and water vapor were computed from E, gbw, H2O_s, Ca, leaf temperature and pressure with PhotoGEA ", as.character(utils::packageVersion("PhotoGEA")),
                    " `calculate_gas_properties` (LI-6800 manual Appendix C; von Caemmerer & Farquhar 1981). Leaf-surface VPD = SVPleaf x (1 - hs). ",
                    "Ball-Berry: gsw = g0 + g1 A hs / Cs, fitted by linear regression (`fit_ball_berry`). Medlyn: gsw = g0 + 1.6 (1 + g1 / sqrt(D)) A / Cs, fitted by nonlinear least squares (`fit_medlyn`) with D = ",
                    if (opt$vpd_source == "licor") "the LI-6800 VPDleaf (leaf to chamber air)." else "leaf-surface VPD."), "",
             "## Limitations", "",
             "- Both models are empirical or semi-empirical and assume each point is at stomatal steady state.",
             "- g0 and g1 trade off; g0 is an extrapolation to A = 0 and is often poorly constrained.",
             "- g1 values from different studies are comparable only if the same model, D definition and g0 treatment were used.", "",
             "## References", "",
             "- Ball JT, Woodrow IE, Berry JA (1987) A model predicting stomatal conductance and its contribution to the control of photosynthesis under different environmental conditions. In: Biggins J (ed) Progress in Photosynthesis Research, vol 4, pp 221-224. Martinus Nijhoff, Dordrecht. https://doi.org/10.1007/978-94-017-0519-6_48",
             "- Medlyn BE, Duursma RA, Eamus D, Ellsworth DS, Prentice IC, Barton CVM, Crous KY, de Angelis P, Freeman M, Wingate L (2011) Reconciling the optimal and empirical approaches to modelling stomatal conductance. Global Change Biology 17: 2134-2144. https://doi.org/10.1111/j.1365-2486.2010.02375.x",
             "- Medlyn BE et al. (2012) Corrigendum to Medlyn et al. (2011). Global Change Biology 18: 3476. https://doi.org/10.1111/j.1365-2486.2012.02790.x",
             "- Leakey ADB, Ferguson JN, Pignon CP, Wu A, Jin Z, Hammer GL, Lobell DB (2019) Water use efficiency as a constraint and target for improving the resilience and productivity of C3 and C4 crops. Annual Review of Plant Biology 70: 781-808. https://doi.org/10.1146/annurev-arplant-042817-040305",
             "- von Caemmerer S, Farquhar GD (1981) Some relationships between the biochemistry of photosynthesis and the gas exchange of leaves. Planta 153: 376-387. https://doi.org/10.1007/BF00384257",
             "- Lochocki EB, Salesse-Smith CE, McGrath JM (2025) PhotoGEA: an R package for closer fitting of photosynthetic gas exchange data with non-Gaussian confidence interval estimation. Plant, Cell & Environment 48: 5104-5119. https://doi.org/10.1111/pce.15501",
             "", sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(lines, file.path(out_dir, "report.md"))
  invisible(w)
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (!requireNamespace("PhotoGEA", quietly = TRUE))
    stop("The R package PhotoGEA is required: install.packages('PhotoGEA').", call. = FALSE)
  if (opt$demo) {
    opt$input <- file.path(SKILL_DIR, "demo", "stomatal_demo.csv")
    if (!nzchar(opt$curve_cols)) opt$curve_cols <- "genotype,plot"
  }
  if (!nzchar(opt$input)) stop("Provide --input FILE.csv (claw-licor-6800 output), or --demo. See --help.", call. = FALSE)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  res <- analyze_stomatal(read_tidy(opt$input), opt)
  save_table(res$params, opt$out, "stomatal_parameters.csv")
  save_table(res$points, opt$out, "stomatal_points.csv")
  fig <- plot_fits(res, opt$out)
  w <- write_report(res, opt, opt$input, opt$out, fig)
  write_repro_bundle(opt$out, argv, inputs = opt$input)
  cat(sprintf("Done: %d curve(s), %d warning(s). Report: %s\n", nrow(res$params), length(w),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(res)
}

if (sys.nframe() == 0L) main()
