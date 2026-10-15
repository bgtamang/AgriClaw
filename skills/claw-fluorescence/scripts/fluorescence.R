#!/usr/bin/env Rscript
# claw-fluorescence: PAM chlorophyll-fluorescence parameters (Fv/Fm, PhiPSII,
# Fv'/Fm', qP, qL, NPQ, estimated Fo', ETR) from a generic CSV or LI-COR
# LI-6800 files, with quality-control flags (dark adaptation, Fv/Fm range,
# Fm' > Fm, impossible signal orders).
#
# Usage:
#   Rscript scripts/fluorescence.R --demo
#   Rscript scripts/fluorescence.R --input fluor.csv --group genotype --dark-adapted yes
#   Rscript scripts/fluorescence.R --input licor.xlsx --format licor --dark-adapted yes
#   Rscript scripts/fluorescence.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-fluorescence"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/fluorescence.R [--demo] [--input FILE] [options]

  --input FILE        CSV (columns fo, fm, fs, fm_prime [, fo_prime, ppfd_umol_m2_s])
                      or an LI-6800 Excel/plaintext file
  --format F          auto (default), csv or licor
  --group COLS        comma-separated grouping columns for the summary (optional)
  --by-ppfd           also group the summary by PPFD (rounded to 10 umol m-2 s-1)
  --exclude-flagged   leave rows with QC flags out of the summary (they stay in
                      the per-row table)
  --dark-adapted A    were fo and fm measured on dark-adapted leaves? yes, no or
                      unknown (default: unknown; Fv/Fm and NPQ need 'yes')
  --absorptance A     leaf absorptance for ETR (default: 0.84; LI-6800 files use
                      the instrument's absorbed PPFD, Qabs_fs, unless this is set)
  --f-psii F          fraction of absorbed light reaching PSII (default: 0.5;
                      LI-6800 files use PS2/1 unless this is set)
  --fvfm-high X       flag Fv/Fm above this (default: 0.86)
  --fvfm-low X        flag Fv/Fm below this (default: 0.75)
  --out DIR           output folder (default: agriclaw_output/claw-fluorescence)
  --demo              run on the bundled synthetic data
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, input = "", format = "auto", group = "", by_ppfd = FALSE,
                 exclude_flagged = FALSE,
                 dark_adapted = "unknown", absorptance = NA_real_, f_psii = NA_real_,
                 fvfm_high = 0.86, fvfm_low = 0.75, out = "")

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

DEFAULT_ABSORPTANCE <- 0.84  # leaf absorptance often used for ETR (see SKILL.md)
DEFAULT_F_PSII <- 0.5        # share of absorbed photons reaching PSII (assumption)

# LI-6800 column names (verified against the example files shipped with PhotoGEA).
LICOR_COLS <- c(fo = "Fo", fm = "Fm", fs = "Fs", fm_prime = "Fm'", fo_prime = "Fo'",
                ppfd = "Qin", qabs = "Qabs_fs", f_psii = "PS2/1")
LICOR_OUTPUTS <- c(fv_fm = "Fv/Fm", phi_psii = "PhiPS2", fvp_fmp = "Fv'/Fm'", qp = "qP",
                   ql = "qL", npq = "NPQ", etr_umol_m2_s = "ETR", fo_prime_est = "alt. Fo'")
LICOR_IDS <- c("DarkAdaptedID", "LightAdaptedID", "DarkPulseID")

# Fluorescence parameters, row by row. Inputs are fluorescence signals in any
# consistent unit; fo, fm from a dark-adapted leaf; fs, fm_prime (and fo_prime
# if measured) in the light.
#   Fv/Fm   = (Fm - Fo) / Fm                         Kitajima & Butler (1975)
#   PhiPSII = (Fm' - Fs) / Fm'                       Genty et al. (1989)
#   Fo'     = Fo / (Fv/Fm + Fo/Fm')  if not measured Oxborough & Baker (1997)
#   Fv'/Fm' = (Fm' - Fo') / Fm'
#   qP      = (Fm' - Fs) / (Fm' - Fo')               (puddle model)
#   qL      = qP x Fo' / Fs                          Kramer et al. (2004) (lake model)
#   NPQ     = (Fm - Fm') / Fm'                       Bilger & Bjorkman (1990)
#   ETR     = PhiPSII x absorbed PPFD x f_psii, absorbed PPFD = PPFD x absorptance
fluor_params <- function(fo, fm, fs, fm_prime, fo_prime = NA_real_, ppfd = NA_real_,
                         absorptance = DEFAULT_ABSORPTANCE, f_psii = DEFAULT_F_PSII,
                         qabs = NA_real_) {
  n <- max(length(fo), length(fm), length(fs), length(fm_prime))
  rn <- function(x) if (length(x) == 1L) rep(x, n) else x
  fo <- rn(fo); fm <- rn(fm); fs <- rn(fs); fm_prime <- rn(fm_prime)
  fo_prime <- rn(fo_prime); ppfd <- rn(ppfd); absorptance <- rn(absorptance)
  f_psii <- rn(f_psii); qabs <- rn(qabs)
  fv_fm <- (fm - fo) / fm
  phi <- (fm_prime - fs) / fm_prime
  fo_est <- fo / (fv_fm + fo / fm_prime)
  measured <- is.finite(fo_prime) & fo_prime > 0
  fo_used <- ifelse(measured, fo_prime, fo_est)
  qp <- (fm_prime - fs) / (fm_prime - fo_used)
  q_abs <- ifelse(is.finite(qabs), qabs, ppfd * absorptance)
  data.frame(fv_fm = fv_fm, phi_psii = phi, fo_prime_used = fo_used,
             fo_prime_source = ifelse(measured, "measured", ifelse(is.finite(fo_est), "estimated", NA_character_)),
             fo_prime_est = fo_est, fvp_fmp = (fm_prime - fo_used) / fm_prime,
             qp = qp, ql = qp * fo_used / fs, npq = (fm - fm_prime) / fm_prime,
             qabs_umol_m2_s = q_abs, etr_umol_m2_s = phi * q_abs * f_psii,
             stringsAsFactors = FALSE)
}

# Quality-control flags. 'invalid' rows break the definitions (excluded from
# summaries); the others are kept but must be checked.
qc_flags <- function(fo, fm, fs, fm_prime, fo_used, fv_fm, fvfm_high = 0.86, fvfm_low = 0.75) {
  f <- function(cond, label) ifelse(!is.na(cond) & cond, label, "")
  invalid <- cbind(f(fo <= 0 | fm <= 0 | fs <= 0 | fm_prime <= 0, "non-positive signal"),
                   f(fo >= fm, "Fo >= Fm"), f(fs >= fm_prime, "Fs >= Fm'"),
                   f(fo_used >= fm_prime, "Fo' >= Fm'"))
  warn <- cbind(f(fv_fm > fvfm_high, sprintf("Fv/Fm > %g", fvfm_high)),
                f(fv_fm < fvfm_low, sprintf("Fv/Fm < %g", fvfm_low)),
                f(fm_prime > fm, "Fm' > Fm (NPQ < 0)"), f(fs < fo_used, "Fs < Fo' (qP > 1)"))
  join <- function(m) apply(m, 1, function(r) paste(r[nzchar(r)], collapse = "; "))
  data.frame(qc_invalid = join(invalid), qc_check = join(warn), stringsAsFactors = FALSE)
}

# Dark-adapted (Fo, Fm) and light-adapted IDs carry timestamps in LI-6800 files,
# e.g. "RECT-12243-20210724-05_20_30". A dark-adapted measurement taken more than
# 12 h before (or after) the light-adapted one is almost certainly from another
# leaf or session: the instrument keeps the last Fo/Fm and reuses it.
licor_id_time <- function(id) {
  m <- regmatches(id, regexpr("[0-9]{8}-[0-9]{2}_[0-9]{2}_[0-9]{2}$", id))
  out <- rep(NA_real_, length(id))
  has <- grepl("[0-9]{8}-[0-9]{2}_[0-9]{2}_[0-9]{2}$", id)
  out[has] <- as.numeric(as.POSIXct(m, format = "%Y%m%d-%H_%M_%S", tz = "UTC"))
  as.POSIXct(out, origin = "1970-01-01", tz = "UTC")
}

detect_format <- function(path) {
  if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) return("licor")
  first <- readLines(path, n = 1L, warn = FALSE)
  if (length(first) && grepl("^\\[Header\\]|^\\[Data\\]|^Open ", first)) return("licor")
  "csv"
}

read_licor <- function(path) {
  if (!requireNamespace("PhotoGEA", quietly = TRUE))
    stop("Reading LI-6800 files needs the PhotoGEA package: install.packages('PhotoGEA').", call. = FALSE)
  x <- PhotoGEA::read_gasex_file(path)
  dat <- x$main_data
  miss <- setdiff(LICOR_COLS[c("fs", "fm_prime")], names(dat))
  if (length(miss))
    stop("LI-6800 file has no fluorometer column(s) ", paste(miss, collapse = ", "),
         ": was the fluorometer head used and were flashes logged?", call. = FALSE)
  dat
}

# Standardise either input to columns fo, fm, fs, fm_prime, fo_prime, ppfd,
# qabs, f_psii (+ any other columns, kept for grouping).
standardise <- function(dat, format, absorptance, f_psii) {
  if (format == "licor") {
    get <- function(k) if (LICOR_COLS[[k]] %in% names(dat)) suppressWarnings(as.numeric(dat[[LICOR_COLS[[k]]]])) else rep(NA_real_, nrow(dat))
    out <- data.frame(fo = get("fo"), fm = get("fm"), fs = get("fs"), fm_prime = get("fm_prime"),
                      fo_prime = get("fo_prime"), ppfd = get("ppfd"),
                      qabs = if (is.finite(absorptance)) NA_real_ else get("qabs"),
                      f_psii = if (is.finite(f_psii)) f_psii else get("f_psii"))
    out$f_psii[!is.finite(out$f_psii)] <- DEFAULT_F_PSII
    out$absorptance <- if (is.finite(absorptance)) absorptance else DEFAULT_ABSORPTANCE
    keep <- setdiff(names(dat), c(unname(LICOR_COLS), unname(LICOR_OUTPUTS), names(out)))
    extra <- dat[, keep, drop = FALSE]
    inst <- as.data.frame(lapply(LICOR_OUTPUTS, function(cl) if (cl %in% names(dat)) suppressWarnings(as.numeric(dat[[cl]])) else NA_real_))
    names(inst) <- paste0("licor_", names(LICOR_OUTPUTS))
    return(cbind(extra, out, inst))
  }
  need <- c("fs", "fm_prime")
  miss <- setdiff(need, names(dat))
  if (length(miss))
    stop("Column(s) not found: ", paste(miss, collapse = ", "), ". A CSV needs fs and fm_prime ",
         "(light-adapted); fo and fm (dark-adapted), fo_prime and ppfd_umol_m2_s are optional. ",
         "For an LI-6800 file use --format licor.", call. = FALSE)
  for (cl in intersect(c("fo", "fm", need, "fo_prime", "ppfd_umol_m2_s", "absorptance"), names(dat)))
    if (!is.numeric(dat[[cl]])) stop("Column '", cl, "' must be numeric.", call. = FALSE)
  out <- dat
  if (!"fo" %in% names(dat)) out$fo <- NA_real_
  if (!"fm" %in% names(dat)) out$fm <- NA_real_
  out$fo_prime <- if ("fo_prime" %in% names(dat)) dat$fo_prime else NA_real_
  out$ppfd <- if ("ppfd_umol_m2_s" %in% names(dat)) dat$ppfd_umol_m2_s else NA_real_
  out$qabs <- NA_real_
  out$absorptance <- if (is.finite(absorptance)) absorptance else if ("absorptance" %in% names(dat)) dat$absorptance else DEFAULT_ABSORPTANCE
  out$f_psii <- if (is.finite(f_psii)) f_psii else DEFAULT_F_PSII
  out
}

analyse <- function(s, dark_adapted = "unknown", fvfm_high = 0.86, fvfm_low = 0.75) {
  fo <- s$fo; fm <- s$fm
  if (dark_adapted == "no") fo <- fm <- rep(NA_real_, nrow(s))  # not dark-adapted: not used
  p <- fluor_params(fo, fm, s$fs, s$fm_prime, s$fo_prime, s$ppfd, s$absorptance, s$f_psii, s$qabs)
  q <- qc_flags(fo, fm, s$fs, s$fm_prime, p$fo_prime_used, p$fv_fm, fvfm_high, fvfm_low)
  cbind(s, p, q)
}

# Group label per row ("all" without grouping columns; missing values become "NA").
group_key <- function(res, group_cols) {
  if (!length(group_cols)) return(rep("all", nrow(res)))
  parts <- lapply(res[group_cols], function(x) ifelse(is.na(x), "NA", as.character(x)))
  do.call(paste, c(parts, sep = ", "))
}

# Group labels in a natural order (numeric grouping columns sorted as numbers).
group_levels <- function(res, group_cols) {
  if (!length(group_cols)) return("all")
  o <- do.call(order, unname(as.list(res[group_cols])))
  unique(group_key(res, group_cols)[o])
}

summarise_groups <- function(res, group_cols, exclude_flagged = FALSE) {
  vars <- c("fv_fm", "phi_psii", "fvp_fmp", "qp", "ql", "npq", "etr_umol_m2_s")
  ok <- !nzchar(res$qc_invalid) & (!exclude_flagged | !nzchar(res$qc_check))
  key <- factor(group_key(res, group_cols), levels = group_levels(res, group_cols))
  rows <- lapply(levels(key), function(k) {
    sel <- key == k
    s <- res[sel & ok, , drop = FALSE]
    out <- data.frame(group = k, n_rows = sum(sel), n_invalid = sum(sel & nzchar(res$qc_invalid)),
                      n_flagged = sum(sel & nzchar(res$qc_check)), n_used = sum(sel & ok),
                      stringsAsFactors = FALSE)
    for (v in vars) {
      x <- s[[v]][is.finite(s[[v]])]
      out[[paste0(v, "_mean")]] <- if (length(x)) mean(x) else NA_real_
      out[[paste0(v, "_se")]] <- if (length(x) > 1) stats::sd(x) / sqrt(length(x)) else NA_real_
      out[[paste0(v, "_n")]] <- length(x)
    }
    out
  })
  do.call(rbind, rows)
}

# --- Figures -------------------------------------------------------------------

plot_fluor <- function(res, group_cols, out_dir, fvfm_high, fvfm_low) {
  key <- group_key(res, group_cols)
  groups <- group_levels(res, group_cols)
  cols <- stats::setNames(if (length(groups) == 1) "#1f78b4" else grDevices::hcl.colors(length(groups), "Dark 3"), groups)
  files <- character()
  # one point per dark-adapted measurement (rows sharing Fo and Fm are one leaf)
  u <- !duplicated(data.frame(key, res$fo, res$fm)) & is.finite(res$fv_fm)
  if (any(u)) {
    fv <- res$fv_fm[u]; k <- key[u]
    open_png(out_dir, "fv_fm.png")
    gi <- match(k, groups)
    graphics::plot(jitter(gi, amount = 0.15), fv, xaxt = "n", xlim = c(0.5, length(groups) + 0.5),
                   ylim = range(c(fv, fvfm_low - 0.02, fvfm_high + 0.02), na.rm = TRUE),
                   pch = ifelse(fv > fvfm_high | fv < fvfm_low, 4, 16), cex = 1.2,
                   col = grDevices::adjustcolor(cols[k], 0.8), xlab = "", ylab = "Fv/Fm (dark-adapted)")
    graphics::axis(1, at = seq_along(groups), labels = groups, las = if (length(groups) > 4) 2 else 1, cex.axis = 0.8)
    graphics::abline(h = 0.83, lty = 2, col = "grey40")
    graphics::abline(h = c(fvfm_low, fvfm_high), lty = 3, col = "#d7301f")
    grDevices::dev.off(); files <- c(files, "fv_fm.png")
  }
  if (any(is.finite(res$ppfd)) && length(unique(round(res$ppfd[is.finite(res$ppfd)]))) > 1) {
    open_png(out_dir, "phipsii_etr_npq.png", width = 10, height = 3.6)
    op <- graphics::par(mfrow = c(1, 3), mar = c(4.2, 4.2, 1, 1))
    for (v in list(c("phi_psii", "PhiPSII"), c("etr_umol_m2_s", "ETR (umol m-2 s-1)"), c("npq", "NPQ"))) {
      if (!any(is.finite(res[[v[1]]]))) {
        graphics::plot.new(); graphics::title(ylab = v[2])
        graphics::text(0.5, 0.5, "not available", col = "grey40")
        next
      }
      graphics::plot(res$ppfd, res[[v[1]]], pch = 16, col = grDevices::adjustcolor(cols[key], 0.7),
                     xlab = "PPFD (umol m-2 s-1)", ylab = v[2])
    }
    if (length(groups) > 1) graphics::legend("topleft", legend = groups, col = cols, pch = 16, bty = "n", cex = 0.8)
    graphics::par(op); grDevices::dev.off(); files <- c(files, "phipsii_etr_npq.png")
  }
  files
}

# --- Main -------------------------------------------------------------------------

fmt <- function(x, d = 3) ifelse(is.finite(x), formatC(x, format = "f", digits = d), "NA")

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (opt$demo) {
    opt$input <- file.path(SKILL_DIR, "demo", "fluorescence.csv")
    if (!"--group" %in% argv) opt$group <- "genotype"
    if (!"--dark-adapted" %in% argv) opt$dark_adapted <- "yes"
    opt$by_ppfd <- TRUE
  }
  if (!nzchar(opt$input)) stop("Provide --input FILE, or --demo. See --help.", call. = FALSE)
  if (!file.exists(opt$input)) stop("Input file not found: ", opt$input, call. = FALSE)
  if (!opt$dark_adapted %in% c("yes", "no", "unknown"))
    stop("--dark-adapted must be yes, no or unknown.", call. = FALSE)
  if (!opt$format %in% c("auto", "csv", "licor")) stop("--format must be auto, csv or licor.", call. = FALSE)
  if (is.finite(opt$absorptance) && (opt$absorptance <= 0 || opt$absorptance > 1))
    stop("--absorptance must be between 0 and 1.", call. = FALSE)
  if (is.finite(opt$f_psii) && (opt$f_psii <= 0 || opt$f_psii > 1))
    stop("--f-psii must be between 0 and 1.", call. = FALSE)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

  warnings_list <- character()
  warn <- function(msg) warnings_list <<- c(warnings_list, msg)
  fmt_in <- if (opt$format == "auto") detect_format(opt$input) else opt$format
  dat <- if (fmt_in == "licor") read_licor(opt$input) else
    utils::read.csv(opt$input, check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
  s <- standardise(dat, fmt_in, opt$absorptance, opt$f_psii)
  group_cols <- if (nzchar(opt$group)) trimws(strsplit(opt$group, ",")[[1]]) else character()
  miss <- setdiff(group_cols, names(s))
  if (length(miss)) stop("Grouping column(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
  na_rows <- !stats::complete.cases(s[, c("fs", "fm_prime")])
  if (all(na_rows)) stop("No row has both fs and fm_prime.", call. = FALSE)
  if (any(na_rows)) {
    warn(sprintf("%d row(s) lack fs or fm_prime and were dropped.", sum(na_rows)))
    s <- s[!na_rows, , drop = FALSE]
  }
  if (!any(is.finite(s$fo) & is.finite(s$fm)))
    warn("No dark-adapted fo and fm: Fv/Fm, NPQ and estimated Fo' were not computed (qP and qL need a measured fo_prime).")
  res <- analyse(s, opt$dark_adapted, opt$fvfm_high, opt$fvfm_low)
  if (opt$by_ppfd && any(is.finite(res$ppfd))) {
    res$ppfd_level <- round(res$ppfd, -1)
    group_cols <- c(group_cols, "ppfd_level")
  }

  # Dark-adaptation status and stale dark-adapted values
  has_dark <- any(is.finite(res$fo) & is.finite(res$fm))
  if (opt$dark_adapted == "unknown" && has_dark)
    warn("Dark adaptation not confirmed (--dark-adapted unknown): Fv/Fm, NPQ and the estimated Fo' are valid only if fo and fm came from leaves dark-adapted for long enough (commonly 20-30 min). Ask the user and re-run with --dark-adapted yes or no.")
  if (opt$dark_adapted == "no")
    warn("fo/fm were not dark-adapted: Fv/Fm, NPQ and estimated Fo' (and qP, qL that depend on it) were not computed.")
  if (all(c("DarkAdaptedID", "LightAdaptedID") %in% names(res))) {
    td <- licor_id_time(as.character(res$DarkAdaptedID)); tl <- licor_id_time(as.character(res$LightAdaptedID))
    hrs <- as.numeric(difftime(tl, td, units = "hours"))
    res$dark_to_light_h <- hrs
    stale <- is.finite(hrs) & (hrs > 12 | hrs < 0)
    res$qc_check[stale] <- ifelse(nzchar(res$qc_check[stale]), paste0(res$qc_check[stale], "; dark-adapted Fo/Fm from another session"),
                                  "dark-adapted Fo/Fm from another session")
    if (any(stale))
      warn(sprintf("%d row(s) use dark-adapted Fo/Fm recorded more than 12 h before (or after) the light-adapted flash (DarkAdaptedID %s): the instrument reused stored values, so Fv/Fm, NPQ and the estimated Fo' for these rows are not valid for this leaf.",
                   sum(stale), paste(utils::head(unique(as.character(res$DarkAdaptedID[stale])), 3), collapse = ", ")))
    if (length(unique(res$DarkAdaptedID)) == 1 && length(unique(res$LightAdaptedID)) > 10)
      warn("Every row shares one DarkAdaptedID: confirm that all rows are from the same leaf; otherwise NPQ uses another leaf's Fm.")
  }
  n_inv <- sum(nzchar(res$qc_invalid))
  if (n_inv) warn(sprintf("%d row(s) are physically impossible (%s) and were excluded from summaries.", n_inv,
                          paste(unique(unlist(strsplit(res$qc_invalid[nzchar(res$qc_invalid)], "; "))), collapse = ", ")))
  for (lab in c(sprintf("Fv/Fm > %g", opt$fvfm_high), sprintf("Fv/Fm < %g", opt$fvfm_low), "Fm' > Fm (NPQ < 0)", "Fs < Fo' (qP > 1)")) {
    k <- sum(grepl(lab, res$qc_check, fixed = TRUE))
    if (!k) next
    why <- switch(lab,
      "Fm' > Fm (NPQ < 0)" = "Fm was not from a fully dark-adapted leaf, came from another leaf, or Fm and Fm' used different flash settings",
      "Fs < Fo' (qP > 1)" = "Fo' is overestimated or Fs mis-measured",
      if (startsWith(lab, "Fv/Fm >")) "implausibly high for healthy leaves (about 0.83): Fo may be underestimated (measuring light, background signal)" else
        "photoinhibition, stress, or incomplete dark adaptation")
    warn(sprintf("%d row(s) with %s: %s.", k, lab, why))
  }
  est <- sum(res$fo_prime_source == "estimated", na.rm = TRUE)
  if (est) warn(sprintf("Fo' was estimated (Oxborough & Baker 1997) for %d row(s); qP, qL and Fv'/Fm' for those rows assume dark-adapted Fo and Fm from the same leaf.", est))
  etr_basis <- if (fmt_in == "licor" && !is.finite(opt$absorptance))
    "absorbed PPFD from the LI-6800 (Qabs_fs, from its own absorptance for the light source) x PS2/1" else
      sprintf("PPFD x absorptance (%s) x f_PSII (%s)", if (is.finite(opt$absorptance)) opt$absorptance else
        if ("absorptance" %in% names(dat)) "per-row column" else DEFAULT_ABSORPTANCE,
        if (is.finite(opt$f_psii)) opt$f_psii else if (fmt_in == "licor") "PS2/1" else DEFAULT_F_PSII)
  if (!any(is.finite(res$ppfd)) && !any(is.finite(res$qabs)))
    warn("No PPFD column (ppfd_umol_m2_s): ETR was not computed.")

  # Agreement with the instrument's own calculations (LI-6800 files)
  cmp_lines <- NULL
  if (fmt_in == "licor") {
    cmp <- sapply(names(LICOR_OUTPUTS), function(v) {
      a <- res[[v]]; b <- res[[paste0("licor_", v)]]
      ok <- is.finite(a) & is.finite(b)
      if (!any(ok)) return(NA_real_)
      max(abs(a[ok] - b[ok]) / pmax(abs(b[ok]), 1e-12))
    })
    cmp_lines <- c("## Agreement with the LI-6800's own values", "",
                   "| parameter | max relative difference |", "|---|---|",
                   sprintf("| %s | %s |", names(cmp), ifelse(is.finite(cmp), formatC(cmp, format = "g", digits = 2), "not in file")), "")
    if (any(cmp > 1e-6, na.rm = TRUE))
      warn(sprintf("Values differ from the LI-6800's for: %s. Check the constants (absorptance, PS2/1) and whether the file was recomputed.",
                   paste(names(cmp)[which(cmp > 1e-6)], collapse = ", ")))
  }

  id_cols <- if (fmt_in == "licor") c("obs", "date", "hhmmss", LICOR_IDS, "dark_to_light_h") else
    setdiff(names(dat), c("fo", "fm", "fs", "fm_prime", "fo_prime", "ppfd_umol_m2_s"))
  out_cols <- c(intersect(c(group_cols, "ppfd_level", id_cols), names(res)),
                "fo", "fm", "fs", "fm_prime", "fo_prime", "ppfd", "fv_fm", "phi_psii", "fo_prime_used", "fo_prime_source",
                "fo_prime_est", "fvp_fmp", "qp", "ql", "npq", "qabs_umol_m2_s", "etr_umol_m2_s",
                grep("^licor_", names(res), value = TRUE), "qc_invalid", "qc_check")
  out_tab <- res[, unique(out_cols[out_cols %in% names(res)]), drop = FALSE]
  names(out_tab)[names(out_tab) == "ppfd"] <- "ppfd_umol_m2_s"
  save_table(out_tab, opt$out, "fluorescence_parameters.csv")
  sm <- summarise_groups(res, group_cols, opt$exclude_flagged)
  save_table(sm, opt$out, "group_summary.csv")
  figs <- plot_fluor(res, setdiff(group_cols, "ppfd_level"), opt$out, opt$fvfm_high, opt$fvfm_low)

  sm_lines <- c("| group | n | flagged | Fv/Fm | PhiPSII | Fv'/Fm' | qP | qL | NPQ | ETR |",
                "|---|---|---|---|---|---|---|---|---|---|",
                sprintf("| %s | %d | %d | %s | %s | %s | %s | %s | %s | %s |", sm$group, sm$n_used, sm$n_flagged,
                        fmt(sm$fv_fm_mean), fmt(sm$phi_psii_mean), fmt(sm$fvp_fmp_mean), fmt(sm$qp_mean),
                        fmt(sm$ql_mean), fmt(sm$npq_mean, 2), fmt(sm$etr_umol_m2_s_mean, 1)))
  report <- c(sprintf("# %s report", SKILL), "",
              if (length(warnings_list)) c("## Warnings", "", paste("-", warnings_list), "") else NULL,
              "## What was done", "",
              sprintf("Input: `%s` (read as %s); %d row(s). Dark-adapted Fo/Fm: %s.", opt$input,
                      if (fmt_in == "licor") "an LI-6800 file via PhotoGEA" else "CSV", nrow(res), opt$dark_adapted), "",
              "| parameter | formula | source |", "|---|---|---|",
              "| Fv/Fm | (Fm - Fo) / Fm | Kitajima & Butler 1975 |",
              "| PhiPSII | (Fm' - Fs) / Fm' | Genty et al. 1989 |",
              "| Fo' (if not measured) | Fo / (Fv/Fm + Fo/Fm') | Oxborough & Baker 1997 |",
              "| Fv'/Fm' | (Fm' - Fo') / Fm' | |",
              "| qP | (Fm' - Fs) / (Fm' - Fo') | |",
              "| qL | qP x Fo' / Fs | Kramer et al. 2004 |",
              "| NPQ | (Fm - Fm') / Fm' | Bilger & Bjorkman 1990 |",
              sprintf("| ETR (umol m-2 s-1) | PhiPSII x %s | |", etr_basis), "",
              "## Summary", "",
              sprintf("Means over valid rows%s%s. Per-row values and flags: `tables/fluorescence_parameters.csv`; SEs: `tables/group_summary.csv`.",
                      if (length(group_cols)) sprintf(" by %s", paste(group_cols, collapse = ", ")) else "",
                      if (opt$exclude_flagged) ", leaving out flagged rows" else ", including flagged rows (column 'flagged'; use --exclude-flagged to leave them out)"), "",
              sm_lines, "",
              cmp_lines,
              if (length(figs)) c(sprintf("![%s](figures/%s)", sub("\\.png$", "", figs), figs), "") else NULL,
              "## Limitations", "",
              "- ETR assumes a leaf absorptance and an equal PSII/PSI share; measured absorptance varies widely among species (Stemke & Santiago 2011), so compare ETR across species or treatments with care.",
              "- qP assumes separate PSII antennae (puddle model); qL assumes connected ones (lake model) and is preferred for land plants (Kramer et al. 2004).",
              "- NPQ compares Fm' with a dark-adapted Fm from the same leaf; it is not comparable between leaves with different Fm.",
              "",
              "## References", "",
              "- Genty B, Briantais J-M, Baker NR (1989) Biochim Biophys Acta 990:87-92. https://doi.org/10.1016/S0304-4165(89)80016-9",
              "- Kramer DM, Johnson G, Kiirats O, Edwards GE (2004) Photosynth Res 79:209-218. https://doi.org/10.1023/B:PRES.0000015391.99477.0d",
              "- Oxborough K, Baker NR (1997) Photosynth Res 54:135-142. https://doi.org/10.1023/A:1005936823310",
              "- Bilger W, Bjorkman O (1990) Photosynth Res 25:173-185. https://doi.org/10.1007/BF00033159",
              "- Kitajima M, Butler WL (1975) Biochim Biophys Acta 376:105-115. https://doi.org/10.1016/0005-2728(75)90209-1",
              "- Murchie EH, Lawson T (2013) J Exp Bot 64:3983-3998. https://doi.org/10.1093/jxb/ert208",
              "- Stemke JA, Santiago LS (2011) Photosynthetica 49:195-200. https://doi.org/10.1007/s11099-011-0026-y",
              "",
              sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(report, file.path(opt$out, "report.md"))
  write_repro_bundle(opt$out, argv, inputs = opt$input)
  cat(sprintf("Done: %d row(s), %d flagged, %d invalid. Report: %s\n", nrow(res), sum(nzchar(res$qc_check)), n_inv,
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(res)
}

if (sys.nframe() == 0L) main()
