#!/usr/bin/env Rscript
# claw-fieldbook: move trial data between the design and the analysis.
#   export: a claw-trial-design field book -> Field Book app (PhenoApps) field
#           import CSV (unique ID + primary and secondary order columns).
#   import: Field Book "database" or "table" exports, or any wide CSV, merged
#           back onto the design, with QC (unknown or duplicated plots, missing
#           data, out-of-range values, robust outliers, date problems) and a
#           clean, analysis-ready table for claw-trial-analysis.
#
# Usage:
#   Rscript scripts/fieldbook.R --demo
#   Rscript scripts/fieldbook.R --mode export --design field_book.csv
#   Rscript scripts/fieldbook.R --mode import --design field_book.csv --data export.csv --limits limits.csv
#   Rscript scripts/fieldbook.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-fieldbook"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/fieldbook.R --mode export|import --design FIELD_BOOK.csv [options]

  --mode M            export (design -> Field Book import file) or import (data -> QC + clean table)
  --design FILE       field book from claw-trial-design (one row per plot)
  --design-id COL     unique plot ID column in the design (default plot_id)
  --primary COL       primary order column for Field Book (default row)
  --secondary COL     secondary order column for Field Book (default col)
Import:
  --data FILES        collected data, comma-separated: Field Book 'database' exports
                      (trait/value columns), 'table' exports, or any wide CSV
  --id-col COL        plot ID column in the data (default: design-id if present,
                      else the column whose values match the design IDs)
  --traits LIST       traits to keep, comma-separated (default: all)
  --limits SPEC       CSV with columns trait,min,max, or 'trait=min:max,trait2=min:max'
  --keep-out-of-range keep values outside --limits (default: set to NA and flag)
  --mad-threshold Z   flag |robust z| above Z (default 3.5; 0 turns it off)
  --duplicates P      repeated plot x trait records: flag (default; conflicting -> NA),
                      first, last (by timestamp), or mean
  --planting-date D   YYYY-MM-DD; dates before it are flagged; adds <trait>_dap columns
  --season-start D    YYYY-MM-DD; observations earlier are flagged
  --season-end D      YYYY-MM-DD; observations later are flagged
  --date-format F     format of date values (default %Y-%m-%d)
  --today D           reference date for 'future date' checks (default: today)
General:
  --demo              export + import on the demo design and simulated Field Book data
  --out DIR           output folder (default: agriclaw_output/claw-fieldbook)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, mode = "", design = "", design_id = "plot_id",
                 primary = "row", secondary = "col", data = "", id_col = "",
                 traits = "", limits = "", keep_out_of_range = FALSE,
                 mad_threshold = 3.5, duplicates = "flag", planting_date = "",
                 season_start = "", season_end = "", date_format = "%Y-%m-%d",
                 today = "", out = "")

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
# Reading and export
# ---------------------------------------------------------------------------

fail <- function(...) stop(paste0(...), call. = FALSE)

# Field Book "database" export columns that are metadata, not plot attributes
# (Field Book source, CSVWriter.writeDatabaseFormat).
FB_META <- c("trait", "value", "timestamp", "person", "location", "number",
             "attached_photo", "attached_video", "attached_audio", "device_name")
MISSING_CODES <- c("", "NA", "N/A", "n/a", "na", "NaN", "NULL", "null", ".", "-")

read_csv_chr <- function(path) {
  if (!file.exists(path)) fail("File not found: ", path)
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
                       na.strings = character(), strip.white = TRUE, fileEncoding = "UTF-8-BOM")
  names(x) <- trimws(names(x))
  if (anyDuplicated(names(x))) fail("Duplicated column names in ", path, ": ",
                                    paste(unique(names(x)[duplicated(names(x))]), collapse = ", "), ".")
  x
}

read_design <- function(path, id) {
  if (!file.exists(path)) fail("Design file not found: ", path)
  d <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  names(d) <- trimws(names(d))
  if (!id %in% names(d))
    fail("The design has no '", id, "' column (the unique plot ID). Columns: ", paste(names(d), collapse = ", "),
         ". Use the field book from claw-trial-design or set --design-id.")
  ids <- trimws(as.character(d[[id]]))
  if (any(is.na(ids) | ids == "")) fail("The design has blank plot IDs in '", id, "'.")
  if (anyDuplicated(ids)) fail("Plot IDs in the design are not unique: ", paste(utils::head(unique(ids[duplicated(ids)]), 10), collapse = ", "), ".")
  d[[id]] <- ids
  if (!"is_filler" %in% names(d)) d$is_filler <- FALSE
  d$is_filler <- toupper(as.character(d$is_filler)) %in% c("TRUE", "T", "1", "YES")
  d
}

# Column names Field Book accepts: it strips [ ] ` " ' from headers; anything
# other than letters, digits, _ . - is replaced here as well.
fb_safe_names <- function(nm) {
  x <- gsub("[][`\"']", "", nm)
  x <- gsub("[^A-Za-z0-9_.-]", "_", x)
  make.unique(x, sep = "_")
}

export_fieldbook <- function(design, id = "plot_id", primary = "row", secondary = "col") {
  for (cl in c(id, primary, secondary))
    if (!cl %in% names(design)) fail("Column '", cl, "' not found in the design. Columns: ", paste(names(design), collapse = ", "), ".")
  ids <- as.character(design[[id]])
  bad <- grepl("[/\\\\]", ids)
  if (any(bad)) fail("Field Book rejects '/' and '\\' in the unique ID; offending IDs: ", paste(utils::head(ids[bad], 5), collapse = ", "), ".")
  for (cl in c(primary, secondary))
    if (anyNA(design[[cl]]) || any(trimws(as.character(design[[cl]])) == ""))
      fail("Order column '", cl, "' has missing values; every plot needs a primary and secondary order.")
  if (anyDuplicated(paste(design[[primary]], design[[secondary]])))
    fail("Some plots share the same ", primary, " and ", secondary, "; Field Book navigation needs a unique position per plot.")
  front <- c(id, primary, secondary)
  out <- design[order(design[[primary]], design[[secondary]]), c(front, setdiff(names(design), front)), drop = FALSE]
  names(out) <- fb_safe_names(names(out))
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Import: data -> long records
# ---------------------------------------------------------------------------

detect_id_col <- function(df, design_ids, id_col, design_id, exclude = character()) {
  if (nzchar(id_col)) {
    if (!id_col %in% names(df)) fail("--id-col '", id_col, "' not found in the data. Columns: ", paste(names(df), collapse = ", "), ".")
    return(id_col)
  }
  if (design_id %in% names(df)) return(design_id)
  cand <- setdiff(names(df), c(exclude, FB_META))
  share <- vapply(cand, function(cl) mean(trimws(df[[cl]]) %in% design_ids), 0)
  if (length(share) && max(share) >= 0.5) return(cand[which.max(share)])
  fail("Cannot tell which data column holds the plot ID (none matches the design's '", design_id,
       "' values). Give --id-col. Columns: ", paste(names(df), collapse = ", "), ".")
}

to_long <- function(df, src, design, opt) {
  design_ids <- design[[opt$design_id]]
  fmt <- if (all(c("trait", "value") %in% names(df))) "database" else "wide"
  idc <- detect_id_col(df, design_ids, opt$id_col, opt$design_id)
  if (fmt == "database") {
    out <- data.frame(plot_id = trimws(df[[idc]]), trait = trimws(df$trait), value = trimws(df$value),
                      timestamp = if ("timestamp" %in% names(df)) trimws(df$timestamp) else NA_character_,
                      source = src, line = seq_len(nrow(df)) + 1L, stringsAsFactors = FALSE)
  } else {
    meta <- c(idc, names(design), FB_META, "source", "line")
    tr <- setdiff(names(df), meta)
    if (!length(tr)) fail("No trait columns found in ", src, " (every column is a design or ID column).")
    out <- do.call(rbind, lapply(tr, function(t)
      data.frame(plot_id = trimws(df[[idc]]), trait = t, value = trimws(df[[t]]), timestamp = NA_character_,
                 source = src, line = seq_len(nrow(df)) + 1L, stringsAsFactors = FALSE)))
  }
  list(long = out, format = fmt, id_col = idc)
}

# ---------------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------------

parse_date <- function(x, fmt = "%Y-%m-%d") {
  d <- suppressWarnings(as.Date(x, format = fmt))
  iso <- is.na(d) & grepl("^\\d{4}-\\d{2}-\\d{2}", x)
  d[iso] <- suppressWarnings(as.Date(substr(x[iso], 1, 10), format = "%Y-%m-%d"))
  d
}

as_num <- function(x) suppressWarnings(as.numeric(x))

parse_limits <- function(spec) {
  if (!nzchar(spec)) return(data.frame(trait = character(), min = numeric(), max = numeric()))
  if (file.exists(spec)) {
    l <- utils::read.csv(spec, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
    names(l) <- tolower(trimws(names(l)))
    if (!all(c("trait", "min", "max") %in% names(l))) fail("The limits file needs columns trait, min, max.")
  } else {
    parts <- trimws(strsplit(spec, ",", fixed = TRUE)[[1]])
    ok <- grepl("^[^=]+=[^:]*:[^:]*$", parts)
    if (!all(ok)) fail("--limits must be a CSV file or 'trait=min:max,trait2=min:max'; could not read: ", paste(parts[!ok], collapse = ", "))
    l <- data.frame(trait = sub("=.*$", "", parts), min = sub("^[^=]*=([^:]*):.*$", "\\1", parts),
                    max = sub("^.*:", "", parts), stringsAsFactors = FALSE)
  }
  l$trait <- trimws(l$trait)
  l$min <- as_num(l$min); l$max <- as_num(l$max)
  if (any(!is.na(l$min) & !is.na(l$max) & l$min > l$max)) fail("A --limits entry has min > max.")
  l
}

opt_date <- function(x, what) {
  if (!nzchar(x)) return(as.Date(NA))
  d <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  if (is.na(d)) fail(what, " must be a date in YYYY-MM-DD form, got '", x, "'.")
  d
}

robust_z <- function(x) {
  m <- stats::median(x, na.rm = TRUE)
  s <- stats::mad(x, center = m, na.rm = TRUE)   # 1.4826 * median(|x - m|)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

flag_row <- function(plot_id, trait, issue, value, detail, action) {
  if (!length(plot_id)) return(NULL)
  data.frame(plot_id = plot_id, trait = trait, issue = issue, value = as.character(value),
             detail = detail, action = action, stringsAsFactors = FALSE)
}

import_fieldbook <- function(design, files, opt) {
  id <- opt$design_id
  pieces <- lapply(files, function(f) to_long(read_csv_chr(f), basename(f), design, opt))
  long <- do.call(rbind, lapply(pieces, `[[`, "long"))
  formats <- vapply(pieces, `[[`, "", "format")
  if (nzchar(opt$traits)) {
    keep <- trimws(strsplit(opt$traits, ",", fixed = TRUE)[[1]])
    miss <- setdiff(keep, unique(long$trait))
    if (length(miss)) fail("Trait(s) not found in the data: ", paste(miss, collapse = ", "), ". Found: ", paste(unique(long$trait), collapse = ", "), ".")
    long <- long[long$trait %in% keep, ]
  }
  traits <- unique(long$trait)
  if (!length(traits)) fail("The data contain no trait values.")
  clash <- intersect(traits, names(design))
  if (length(clash)) fail("Trait name(s) clash with design columns: ", paste(clash, collapse = ", "), ". Rename them.")
  today <- if (nzchar(opt$today)) opt_date(opt$today, "--today") else Sys.Date()
  plant <- opt_date(opt$planting_date, "--planting-date")
  s0 <- opt_date(opt$season_start, "--season-start")
  s1 <- opt_date(opt$season_end, "--season-end")
  dup_policy <- match.arg(opt$duplicates, c("flag", "first", "last", "mean"))
  flags <- list()
  addf <- function(x) if (!is.null(x)) flags[[length(flags) + 1L]] <<- x

  # 1. Unknown plot IDs and data recorded on filler plots
  unknown <- !long$plot_id %in% design[[id]]
  if (any(unknown)) {
    u <- long[unknown, ]
    agg <- stats::aggregate(line ~ plot_id + source, data = u, FUN = function(z) paste(utils::head(z, 5), collapse = ";"))
    addf(flag_row(agg$plot_id, "", "unknown_plot_id", "", paste0("not in the design; ", agg$source, " line(s) ", agg$line), "excluded"))
  }
  long <- long[!unknown, ]
  fill_ids <- design[[id]][design$is_filler]
  onfill <- long$plot_id %in% fill_ids & !long$value %in% MISSING_CODES
  if (any(onfill))
    addf(flag_row(long$plot_id[onfill], long$trait[onfill], "data_on_filler", long$value[onfill], "value recorded on a FILLER plot", "excluded"))
  long <- long[!long$plot_id %in% fill_ids, ]

  # 2. Timestamps
  long$obs_date <- parse_date(long$timestamp)
  has_ts <- !is.na(long$timestamp) & !long$timestamp %in% MISSING_CODES
  bad_ts <- has_ts & is.na(long$obs_date)
  addf(flag_row(long$plot_id[bad_ts], long$trait[bad_ts], "bad_timestamp", long$timestamp[bad_ts], "timestamp is not a date", "kept"))
  chk_dates <- function(d, sel, what, impossible_action) {
    fut <- sel & !is.na(d) & d > today
    addf(flag_row(long$plot_id[fut], long$trait[fut], "future_date", format(d[fut]), paste(what, "after", format(today)), impossible_action))
    out <- fut
    if (!is.na(plant)) {
      bp <- sel & !is.na(d) & d < plant
      addf(flag_row(long$plot_id[bp], long$trait[bp], "before_planting", format(d[bp]), paste(what, "before planting date", format(plant)), impossible_action))
      out <- out | bp
    }
    if (!is.na(s0) || !is.na(s1)) {
      os <- sel & !is.na(d) & ((!is.na(s0) & d < s0) | (!is.na(s1) & d > s1)) & !out
      addf(flag_row(long$plot_id[os], long$trait[os], "outside_season", format(d[os]),
                    sprintf("%s outside season window %s to %s", what, format(s0), format(s1)), "kept"))
    }
    out
  }
  invisible(chk_dates(long$obs_date, has_ts & !bad_ts, "observation timestamp", "kept"))

  # 3. Missing codes and trait types
  long$value[long$value %in% MISSING_CODES] <- NA
  types <- stats::setNames(character(length(traits)), traits)
  for (t in traits) {
    v <- long$value[long$trait == t & !is.na(long$value)]
    if (!length(v)) { types[t] <- "empty"; next }
    pn <- mean(!is.na(as_num(v)))
    pd <- mean(!is.na(parse_date(v, opt$date_format)))
    types[t] <- if (pn > 0.5) "numeric" else if (pd > 0.5) "date" else "text"
  }
  long$num <- NA_real_
  long$date <- as.Date(NA)
  for (t in traits) {
    i <- which(long$trait == t & !is.na(long$value))
    if (types[t] == "numeric") {
      x <- as_num(long$value[i])
      bad <- is.na(x)
      det <- ifelse(grepl("^-?[0-9]+,[0-9]+$", long$value[i][bad]), "not a number (decimal comma?)", "not a number")
      addf(flag_row(long$plot_id[i][bad], t, "non_numeric_value", long$value[i][bad], det, "set to NA"))
      long$num[i] <- x
      long$value[i][bad] <- NA
    } else if (types[t] == "date") {
      x <- parse_date(long$value[i], opt$date_format)
      bad <- is.na(x)
      addf(flag_row(long$plot_id[i][bad], t, "unparseable_date", long$value[i][bad],
                    paste("not a date in format", opt$date_format), "set to NA"))
      long$date[i] <- x
      long$value[i][bad] <- NA
    }
  }

  # 4. Duplicates: more than one record for a plot x trait
  key <- paste(long$plot_id, long$trait, sep = "\r")
  dk <- unique(key[duplicated(key)])
  long$keep <- TRUE
  resolved <- list()
  for (k in dk) {
    i <- which(key == k)
    vals <- long$value[i]
    pid <- long$plot_id[i[1]]; t <- long$trait[i[1]]
    nn <- vals[!is.na(vals)]
    if (length(unique(nn)) <= 1) {
      addf(flag_row(pid, t, "duplicate_identical", if (length(nn)) nn[1] else NA, sprintf("%d identical records collapsed", length(i)), "collapsed"))
      long$keep[i[-1]] <- FALSE
      if (length(nn)) { j <- i[!is.na(vals)][1]; long$keep[i] <- FALSE; long$keep[j] <- TRUE }
      next
    }
    ord <- i[order(long$obs_date[i], long$line[i], na.last = TRUE)]
    detail <- sprintf("%d records: %s", length(i), paste(vals, collapse = " | "))
    long$keep[i] <- FALSE
    if (dup_policy == "flag") {
      j <- i[1]; long$keep[j] <- TRUE
      long$value[j] <- NA; long$num[j] <- NA; long$date[j] <- as.Date(NA)
      addf(flag_row(pid, t, "duplicate_conflict", paste(vals, collapse = " | "), detail, "set to NA"))
    } else if (dup_policy %in% c("first", "last")) {
      cand <- ord[!is.na(long$value[ord])]
      j <- if (dup_policy == "first") cand[1] else cand[length(cand)]
      long$keep[j] <- TRUE
      addf(flag_row(pid, t, "duplicate_conflict", paste(vals, collapse = " | "), detail, paste("kept", dup_policy, "record")))
    } else {
      j <- i[1]; long$keep[j] <- TRUE
      if (types[t] == "numeric") {
        long$num[j] <- mean(long$num[i], na.rm = TRUE); long$value[j] <- format(long$num[j])
        addf(flag_row(pid, t, "duplicate_conflict", paste(vals, collapse = " | "), detail, "kept mean"))
      } else {
        long$value[j] <- NA; long$date[j] <- as.Date(NA)
        addf(flag_row(pid, t, "duplicate_conflict", paste(vals, collapse = " | "), paste(detail, "(mean impossible for non-numeric trait)"), "set to NA"))
      }
    }
  }
  long <- long[long$keep, ]

  # 5. Range limits, robust outliers and date-trait checks
  lim <- parse_limits(opt$limits)
  unk_lim <- setdiff(lim$trait, traits)
  long$set_na <- FALSE
  for (t in traits) {
    i <- which(long$trait == t & !is.na(long$num))
    if (types[t] == "numeric" && t %in% lim$trait && length(i)) {
      lo <- lim$min[lim$trait == t]; hi <- lim$max[lim$trait == t]
      out <- (!is.na(lo) & long$num[i] < lo) | (!is.na(hi) & long$num[i] > hi)
      act <- if (opt$keep_out_of_range) "kept" else "set to NA"
      addf(flag_row(long$plot_id[i][out], t, "out_of_range", long$value[i][out], sprintf("outside [%s, %s]", lo, hi), act))
      if (!opt$keep_out_of_range) long$set_na[i[out]] <- TRUE
    }
    if (types[t] == "numeric" && opt$mad_threshold > 0) {
      i <- which(long$trait == t & !is.na(long$num) & !long$set_na)
      if (length(i) >= 10) {
        z <- robust_z(long$num[i])
        o <- !is.na(z) & abs(z) > opt$mad_threshold
        addf(flag_row(long$plot_id[i][o], t, "mad_outlier", long$value[i][o], sprintf("robust z = %.1f (|z| > %g)", z[o], opt$mad_threshold), "kept"))
      }
    }
    if (types[t] == "date") {
      i <- long$trait == t & !is.na(long$date)
      bad <- chk_dates(long$date, i, "date value", if (opt$keep_out_of_range) "kept" else "set to NA")
      if (!opt$keep_out_of_range) long$set_na[bad] <- TRUE
    }
  }
  long$value[long$set_na] <- NA
  long$num[long$set_na] <- NA
  long$date[long$set_na] <- as.Date(NA)

  # 6. Missing values per plot and trait (filler plots are not expected to have data)
  real <- design[!design$is_filler, ]
  seen <- unique(long$plot_id)
  notobs <- setdiff(real[[id]], seen)
  addf(flag_row(notobs, "", "plot_not_observed", "", "no record for any trait", "kept (all traits NA)"))
  for (t in traits) {
    have <- long$plot_id[long$trait == t & !is.na(long$value)]
    miss <- setdiff(real[[id]], c(have, notobs))
    if (length(miss) == nrow(real)) {
      addf(flag_row("", t, "trait_not_observed", "", "no valid value on any plot", "info"))
    } else {
      addf(flag_row(miss, t, "missing_value", "", "no valid value for this plot", "kept as NA"))
    }
  }

  # 7. Clean wide table: one row per (non-filler) design plot
  clean <- real
  for (t in traits) {
    li <- long[long$trait == t, ]
    m <- match(clean[[id]], li$plot_id)
    clean[[t]] <- switch(types[t], numeric = li$num[m], date = format(li$date[m]), li$value[m])
    if (types[t] == "date") {
      clean[[t]][clean[[t]] == "NA"] <- NA
      if (!is.na(plant)) clean[[paste0(t, "_dap")]] <- as.integer(li$date[m] - plant)
    }
  }
  fl <- do.call(rbind, flags)
  if (is.null(fl)) fl <- flag_row(character(), character(), character(), character(), character(), character())
  if (is.null(fl)) fl <- data.frame(plot_id = character(), trait = character(), issue = character(),
                                    value = character(), detail = character(), action = character())
  serious <- fl[!fl$issue %in% c("missing_value", "plot_not_observed", "trait_not_observed") & nzchar(fl$plot_id), ]
  clean$qc_issues <- vapply(clean[[id]], function(p) {
    s <- serious[serious$plot_id == p, ]
    if (!nrow(s)) "" else paste(unique(paste0(s$trait, ":", s$issue)), collapse = "; ")
  }, "")

  # 8. Per-trait summary
  cnt <- function(t, iss) sum(fl$trait == t & fl$issue == iss)
  summ <- do.call(rbind, lapply(traits, function(t) {
    x <- if (types[t] == "numeric") clean[[t]] else NULL
    data.frame(trait = t, type = types[t], plots = nrow(clean), observed = sum(!is.na(clean[[t]])),
               missing = sum(is.na(clean[[t]])), out_of_range = cnt(t, "out_of_range"),
               mad_outliers = cnt(t, "mad_outlier"), duplicate_conflicts = cnt(t, "duplicate_conflict"),
               non_numeric = cnt(t, "non_numeric_value"), date_problems = sum(fl$trait == t & fl$issue %in% c("future_date", "before_planting", "outside_season", "unparseable_date")),
               min = if (!is.null(x) && any(!is.na(x))) min(x, na.rm = TRUE) else NA,
               median = if (!is.null(x) && any(!is.na(x))) stats::median(x, na.rm = TRUE) else NA,
               max = if (!is.null(x) && any(!is.na(x))) max(x, na.rm = TRUE) else NA,
               stringsAsFactors = FALSE)
  }))
  list(clean = clean, flags = fl, summary = summ, types = types, formats = formats,
       id_cols = vapply(pieces, `[[`, "", "id_col"), today = today, unknown_limit_traits = unk_lim,
       long = long)
}

# ---------------------------------------------------------------------------
# Figures and report
# ---------------------------------------------------------------------------

INK <- "#0b0b0b"; INK2 <- "#52514e"; SURFACE <- "#fcfcfb"
SEQ <- c("#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b")

draw_qc_maps <- function(res, design, opt, out_dir) {
  num <- names(res$types)[res$types == "numeric"]
  if (!length(num) || !all(c("row", "col") %in% names(design))) return(FALSE)
  num <- utils::head(num, 6)
  cl <- res$clean
  R <- max(design$row, na.rm = TRUE); C <- max(design$col, na.rm = TRUE)
  np <- length(num); nc <- min(2, np); nr <- ceiling(np / nc)
  open_png(out_dir, "qc_field_maps.png", width = min(16, 1.2 + nc * max(4, C * 0.32)), height = min(16, 0.6 + nr * max(3, R * 0.36 + 1)))
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(nr, nc), mar = c(1.2, 2.4, 3, 0.6), bg = SURFACE, fg = INK2, col.axis = INK2)
  for (t in num) {
    z <- matrix(NA_real_, R, C)
    z[cbind(cl$row, cl$col)] <- cl[[t]]
    graphics::plot.new()
    graphics::plot.window(xlim = c(0.5, C + 0.5), ylim = c(R + 0.5, 0.5), xaxs = "i", yaxs = "i")
    rng <- range(z, na.rm = TRUE)
    if (all(is.finite(rng))) {
      br <- seq(rng[1], rng[2] + 1e-9 * max(1, abs(rng[2])), length.out = length(SEQ) + 1)
      cols <- SEQ[pmax(1, findInterval(z, br, rightmost.closed = TRUE))]
      cols[is.na(z)] <- "#ffffff"
      rr <- row(z); cc <- col(z)
      graphics::rect(cc - 0.5, rr + 0.5, cc + 0.5, rr - 0.5, col = cols, border = SURFACE, lwd = 0.6)
    }
    miss <- which(is.na(z), arr.ind = TRUE)
    if (nrow(miss)) graphics::points(miss[, 2], miss[, 1], pch = 4, cex = 0.6, col = INK2)
    f <- res$flags[res$flags$trait == t & res$flags$issue %in% c("out_of_range", "mad_outlier", "duplicate_conflict", "non_numeric_value"), ]
    if (nrow(f)) {
      m <- match(f$plot_id, cl[[opt$design_id]])
      m <- m[!is.na(m)]
      graphics::rect(cl$col[m] - 0.45, cl$row[m] + 0.45, cl$col[m] + 0.45, cl$row[m] - 0.45, border = "#eb6834", lwd = 2.2)
    }
    graphics::axis(2, at = unique(round(pretty(seq_len(R)))), las = 1, tick = FALSE, cex.axis = 0.7, line = -0.8)
    graphics::mtext(t, side = 3, line = 1.2, adj = 0, cex = 0.85, font = 2, col = INK)
    graphics::mtext(sprintf("light = low (%.4g), dark = high (%.4g); x = missing; orange box = flagged", rng[1], rng[2]),
                    side = 3, line = 0.3, adj = 0, cex = 0.55, col = INK2)
  }
  TRUE
}

md_table <- function(df) {
  df[] <- lapply(df, function(x) { x <- as.character(x); x[is.na(x)] <- ""; gsub("\\|", "/", x) })
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

REFERENCES <- c(
  "- Rife TW, Poland JA (2014). Field Book: An open-source application for field data collection on Android. Crop Science 54(4):1624-1627. https://doi.org/10.2135/cropsci2013.08.0579",
  "- Field Book documentation (PhenoApps): field import and export formats. https://field-book.readthedocs.io/",
  "- Leys C, Ley C, Klein O, Bernard P, Licata L (2013). Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median. Journal of Experimental Social Psychology 49(4):764-766. https://doi.org/10.1016/j.jesp.2013.03.013",
  "- Iglewicz B, Hoaglin DC (1993). How to Detect and Handle Outliers. ASQC Basic References in Quality Control: Statistical Techniques, vol. 16. Modified z-score cut-off 3.5 as summarised in the NIST/SEMATECH e-Handbook of Statistical Methods, section 1.3.5.17: https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm")

import_warnings <- function(res, opt) {
  fl <- res$flags
  n <- function(iss) sum(fl$issue == iss)
  w <- character()
  if (n("unknown_plot_id")) w <- c(w, sprintf("%d plot ID(s) in the data are not in the design and were excluded: %s. Check for typos or the wrong trial.", n("unknown_plot_id"), paste(utils::head(unique(fl$plot_id[fl$issue == "unknown_plot_id"]), 6), collapse = ", ")))
  if (n("duplicate_conflict")) w <- c(w, sprintf("%d plot x trait combination(s) have conflicting repeated records (policy: %s). Resolve them with the data collector.", n("duplicate_conflict"), opt$duplicates))
  if (n("out_of_range")) w <- c(w, sprintf("%d value(s) outside the --limits ranges (%s).", n("out_of_range"), if (opt$keep_out_of_range) "kept" else "set to NA"))
  if (n("non_numeric_value")) w <- c(w, sprintf("%d non-numeric value(s) in numeric traits were set to NA (decimal commas or notes typed into a number field?).", n("non_numeric_value")))
  if (n("mad_outlier")) w <- c(w, sprintf("%d value(s) flagged as robust outliers (|z| > %g); they are kept. Check them against field notes before analysis.", n("mad_outlier"), opt$mad_threshold))
  dp <- sum(fl$issue %in% c("future_date", "before_planting", "unparseable_date", "bad_timestamp"))
  if (dp) w <- c(w, sprintf("%d date problem(s) (future, before planting or unreadable dates).", dp))
  if (n("outside_season")) w <- c(w, sprintf("%d observation(s) dated outside the season window.", n("outside_season")))
  if (n("data_on_filler")) w <- c(w, sprintf("%d value(s) were recorded on FILLER plots and excluded.", n("data_on_filler")))
  if (n("plot_not_observed")) w <- c(w, sprintf("%d design plot(s) have no data at all: %s.", n("plot_not_observed"), paste(utils::head(fl$plot_id[fl$issue == "plot_not_observed"], 8), collapse = ", ")))
  if (n("trait_not_observed")) w <- c(w, sprintf("Trait(s) with no valid value: %s.", paste(fl$trait[fl$issue == "trait_not_observed"], collapse = ", ")))
  if (length(res$unknown_limit_traits)) w <- c(w, sprintf("--limits names trait(s) not in the data: %s.", paste(res$unknown_limit_traits, collapse = ", ")))
  num <- names(res$types)[res$types == "numeric"]
  nolim <- setdiff(num, parse_limits(opt$limits)$trait)
  if (length(nolim)) w <- c(w, sprintf("No --limits for: %s. Only the robust-outlier screen was applied to them.", paste(nolim, collapse = ", ")))
  w
}

write_report <- function(out_dir, mode, opt, exp = NULL, imp = NULL, warn = character(), maps = FALSE) {
  L <- c(sprintf("# %s report", SKILL), "")
  if (length(warn)) L <- c(L, "## Warnings", "", paste("-", warn), "")
  if (!is.null(exp)) {
    L <- c(L, "## Field Book import file", "",
           sprintf("`tables/fieldbook_import.csv`: %d plots, %d columns. Copy it to Field Book's `field_import` folder (or import it from storage).", nrow(exp), ncol(exp)), "",
           sprintf("At import choose **%s** as the unique identifier, **%s** as the primary order and **%s** as the secondary order.",
                   names(exp)[1], names(exp)[2], names(exp)[3]), "",
           "- Unique IDs are checked to be non-empty, unique and free of '/' and '\\', which Field Book rejects.",
           "- Column names are stripped of [ ] ` \" ' (Field Book removes them) and other symbols are replaced by '_'.",
           "- Every field you import must have IDs that are unique across all fields on the device.", "")
  }
  if (!is.null(imp)) {
    L <- c(L, "## Data import and QC", "",
           sprintf("Data: %s (format: %s; plot ID column: %s). Traits: %s.", paste(basename(strsplit(opt$data, ",")[[1]]), collapse = ", "),
                   paste(unique(imp$formats), collapse = ", "), paste(unique(imp$id_cols), collapse = ", "), paste(names(imp$types), collapse = ", ")), "",
           sprintf("Reference date for future-date checks: %s. Duplicate policy: %s. Robust-outlier threshold: |z| > %g.", format(imp$today), opt$duplicates, opt$mad_threshold), "",
           md_table(imp$summary), "")
    if (nrow(imp$flags)) {
      tab <- as.data.frame(table(issue = imp$flags$issue), stringsAsFactors = FALSE)
      names(tab) <- c("issue", "count")
      L <- c(L, "### Flags by issue", "", md_table(tab), "")
    } else {
      L <- c(L, "No QC flags: every design plot has a valid value for every trait.", "")
    }
    ser <- imp$flags[!imp$flags$issue %in% c("missing_value"), ]
    if (nrow(ser)) L <- c(L, "### Flagged values (first 40, missing values omitted)", "", md_table(utils::head(ser, 40)), "")
    if (maps) L <- c(L, "![QC field maps](figures/qc_field_maps.png)", "")
    L <- c(L, "## Files", "",
           "- `tables/clean_data.csv`: one row per design plot (FILLER plots dropped), design columns + one column per trait; numeric traits are numbers, date traits ISO dates (plus `<trait>_dap` when a planting date is given); `qc_issues` lists flags per plot.",
           "- `tables/qc_flags.csv`: every flag with plot, trait, issue, value, detail and the action taken.",
           "- `tables/qc_summary.csv`: per-trait counts.",
           "- Next step: `claw-trial-analysis --input tables/clean_data.csv --traits <trait>`.", "")
  }
  L <- c(L, "## References", "", REFERENCES, "", sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(L, file.path(out_dir, "report.md"))
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (opt$demo) {
    dd <- file.path(SKILL_DIR, "demo")
    opt$mode <- "both"
    opt$design <- file.path(dd, "design_field_book.csv")
    opt$data <- file.path(dd, "fieldbook_database_export.csv")
    if (!nzchar(opt$limits)) opt$limits <- file.path(dd, "limits.csv")
    if (!nzchar(opt$planting_date)) opt$planting_date <- "2026-05-15"
    if (!nzchar(opt$season_start)) opt$season_start <- "2026-05-01"
    if (!nzchar(opt$season_end)) opt$season_end <- "2026-11-30"
  }
  mode <- tolower(opt$mode)
  if (!mode %in% c("export", "import", "both")) fail("Give --mode export or --mode import (or --demo). See --help.")
  if (!nzchar(opt$design)) fail("Give --design FIELD_BOOK.csv (the field book from claw-trial-design).")
  if (mode %in% c("import", "both") && !nzchar(opt$data)) fail("Give --data with the collected data file(s).")
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  design <- read_design(opt$design, opt$design_id)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  exp <- imp <- NULL
  warn <- character()
  maps <- FALSE
  inputs <- opt$design
  if (mode %in% c("export", "both")) {
    exp <- export_fieldbook(design, opt$design_id, opt$primary, opt$secondary)
    d <- file.path(opt$out, "tables"); dir.create(d, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(exp, file.path(d, "fieldbook_import.csv"), row.names = FALSE, na = "")
  }
  if (mode %in% c("import", "both")) {
    if (!nzchar(opt$data)) fail("Give --data with the collected data file(s).")
    files <- trimws(strsplit(opt$data, ",", fixed = TRUE)[[1]])
    imp <- import_fieldbook(design, files, opt)
    save_table(imp$clean, opt$out, "clean_data.csv")
    save_table(imp$flags, opt$out, "qc_flags.csv")
    save_table(imp$summary, opt$out, "qc_summary.csv")
    maps <- draw_qc_maps(imp, design, opt, opt$out)
    warn <- import_warnings(imp, opt)
    inputs <- c(inputs, files, if (file.exists(opt$limits)) opt$limits)
  }
  write_report(opt$out, mode, opt, exp, imp, warn, maps)
  write_repro_bundle(opt$out, argv, inputs = inputs)
  cat(sprintf("Done (%s). Report: %s\n", mode, normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  if (length(warn)) cat(paste("Warning:", warn), sep = "\n")
  invisible(list(export = exp, import = imp))
}

if (sys.nframe() == 0L) main()
