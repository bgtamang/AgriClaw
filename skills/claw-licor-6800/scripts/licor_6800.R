#!/usr/bin/env Rscript
# claw-licor-6800: read LI-COR LI-6800 log files (plaintext or Excel) with
# PhotoGEA, stop on Excel files whose formulas were never evaluated, combine
# files with a file and instrument ID, flag questionable points, check whether
# treatments were confounded with instruments, and write one tidy CSV that
# claw-aci-fit, claw-light-response and claw-stomatal-models read.
#
# Usage:
#   Rscript scripts/licor_6800.R --demo
#   Rscript scripts/licor_6800.R --input "day1.xlsx,day2" --curve-cols genotype,plot --design-cols genotype
#   Rscript scripts/licor_6800.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT license.

SKILL   <- "claw-licor-6800"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/licor_6800.R [--demo] [--input FILES] [--input-dir DIR]
                                  [--curve-cols COLS] [--design-cols COLS]
                                  [--instrument-col COL] [--tleaf-min C] [--tleaf-max C]
                                  [--min-delta-co2 PPM] [--tleaf-drift C] [--out DIR]

  --demo            run on two LI-6800 Excel files shipped with PhotoGEA
  --input           LI-6800 log file(s): .xlsx or plaintext (no extension), comma-separated
  --input-dir       folder of LI-6800 log files (all .xlsx and plaintext logs inside)
  --curve-cols      user-constant column(s) that together identify one curve or leaf,
                    comma-separated (e.g. genotype,plot); enables per-curve checks
  --design-cols     treatment/genotype column(s) to check for confounding with instrument
  --instrument-col  column holding the instrument ID; default: serial numbers from the file
  --tleaf-min       lowest plausible leaf temperature, C (default 5; rule of thumb)
  --tleaf-max       highest plausible leaf temperature, C (default 45; rule of thumb)
  --min-delta-co2   flag |CO2_r - CO2_s| below this, umol/mol (default 1; rule of thumb)
  --tleaf-drift     flag curves whose leaf temperature range exceeds this, C (default 2)
  --out             output folder (default: agriclaw_output/claw-licor-6800)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, input = "", input_dir = "", curve_cols = "",
                 design_cols = "", instrument_col = "", tleaf_min = 5,
                 tleaf_max = 45, min_delta_co2 = 1, tleaf_drift = 2, out = "")

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

# Tidy column contract shared with claw-aci-fit, claw-light-response and
# claw-stomatal-models. `licor` is the LI-6800 column name as read (and
# standardized) by PhotoGEA; `units` is what PhotoGEA must report for it.
TIDY_MAP <- data.frame(
  tidy = c("obs", "time_unix_s", "elapsed_s", "date_local", "hhmmss",
           "a_umol_m2_s", "ci_umol_mol", "ca_umol_mol", "gsw_mol_m2_s",
           "gbw_mol_m2_s", "gtc_mol_m2_s", "e_mol_m2_s", "tleaf_c",
           "tleaf_thermocouple_c", "tair_c", "qin_umol_m2_s", "qabs_umol_m2_s",
           "leaf_absorptance_assumed", "co2_r_umol_mol", "co2_s_umol_mol",
           "co2_r_sp_umol_mol", "h2o_r_mmol_mol", "h2o_s_mmol_mol",
           "rh_cham_pct", "vpd_leaf_kpa", "svp_leaf_kpa", "pa_kpa",
           "delta_pcham_kpa", "oxygen_pct", "flow_umol_s", "leaf_area_cm2",
           "stomatal_ratio_k", "stable_n", "stable_total", "a_ok", "gsw_ok"),
  licor = c("obs", "time", "elapsed", "date", "hhmmss", "A", "Ci", "Ca", "gsw",
            "gbw", "gtc", "E", "TleafCnd", "Tleaf", "Tair", "Qin", "Qabs",
            "alpha", "CO2_r", "CO2_s", "CO2_r_sp", "H2O_r", "H2O_s", "RHcham",
            "VPDleaf", "SVPleaf", "Pa", "DeltaPcham", "Oxygen", "Flow", "S",
            "K", "Stable", "Total", "A:OK", "gsw:OK"),
  units = c(NA, "s", "s", NA, NA, "micromol m^(-2) s^(-1)", "micromol mol^(-1)",
            "micromol mol^(-1)", "mol m^(-2) s^(-1)", "mol m^(-2) s^(-1)",
            "mol m^(-2) s^(-1)", "mol m^(-2) s^(-1)", "degrees C", "degrees C",
            "degrees C", "micromol m^(-2) s^(-1)", "micromol m^(-2) s^(-1)", NA,
            "micromol mol^(-1)", "micromol mol^(-1)", "micromol mol^(-1)",
            "mmol mol^(-1)", "mmol mol^(-1)", "%", "kPa", "kPa", "kPa", "kPa",
            "percent", "micromol s^(-1)", "cm^2", NA, NA, NA, NA, NA),
  type = c("num", "num", "num", "chr", "chr", rep("num", 31)),
  stringsAsFactors = FALSE)

# Computed (formula) columns that must not be all zero/empty, and raw sensor
# columns that are stored as values even in unevaluated Excel files.
CORE_COMPUTED <- c("A", "gsw", "Ci")
RAW_SENSOR <- c("CO2_r", "CO2_s", "H2O_r", "H2O_s", "Flow")
QC_FLAG_COLS <- c("qc_missing_values", "qc_ci_negative", "qc_ci_gt_ca",
                  "qc_gsw_negative", "qc_tleaf_out_of_range", "qc_duplicate",
                  "qc_unstable", "qc_small_delta_co2")
HARD_FLAGS <- c("qc_missing_values", "qc_ci_negative", "qc_ci_gt_ca",
                "qc_gsw_negative", "qc_tleaf_out_of_range", "qc_duplicate")

unevaluated_message <- function(path, cols) {
  paste0(
    "LI-6800 Excel file '", basename(path), "': the computed column(s) ",
    paste(cols, collapse = ", "), " are all zero or empty, although the raw ",
    "CO2/H2O readings are present. The LI-6800 stores A, gsw, Ci and the other ",
    "computed variables as Excel formulas, and they are only evaluated when the ",
    "file is opened and saved in Excel; a program reading the file sees 0 or ",
    "blank instead (documented in ?PhotoGEA::read_licor_6800_Excel).\n",
    "Fix: use the plaintext log the LI-6800 always writes next to the .xlsx ",
    "(same name, no extension), which stores the computed values. Or open the ",
    ".xlsx in Excel, choose Formulas > Calculate Now (F9), save, and run again. ",
    "Do not type values in by hand.")
}

dead_column <- function(v) {
  v <- suppressWarnings(as.numeric(v))
  length(v) == 0L || all(is.na(v) | v == 0)
}

# Stop with an explanation when the computed gas-exchange columns of an Excel
# file are all zero or empty (formulas never evaluated). Returns the names of
# any other computed columns that are all zero or empty.
check_unevaluated <- function(md, path, file_type) {
  if (nrow(md) == 0L) stop("'", basename(path), "' contains no logged observations.", call. = FALSE)
  core <- intersect(CORE_COMPUTED, names(md))
  if (!length(core))
    stop("'", basename(path), "' has none of the columns A, gsw, Ci: it does not look like an ",
         "LI-6800 log. Use the unmodified instrument file.", call. = FALSE)
  dead <- vapply(core, function(cn) dead_column(md[[cn]]), logical(1))
  raw <- intersect(RAW_SENSOR, names(md))
  raw_alive <- any(!vapply(raw, function(cn) dead_column(md[[cn]]), logical(1)))
  if (all(dead)) {
    if (raw_alive && file_type == "Excel") stop(unevaluated_message(path, core), call. = FALSE)
    stop("'", basename(path), "': columns ", paste(core, collapse = ", "),
         " are all zero or empty and so are the raw CO2/H2O readings. The file holds no usable ",
         "measurements.", call. = FALSE)
  }
  extra <- intersect(c("E", "Ca", "gtc", "gbw"), names(md))
  c(core[dead], extra[vapply(extra, function(cn) dead_column(md[[cn]]), logical(1))])
}

file_type_of <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "xlsx") return("Excel")
  if (ext %in% c("", "txt")) return("plaintext")
  if (ext == "xls")
    stop("'", basename(path), "' is a legacy .xls file. The LI-6800 writes .xlsx and plaintext ",
         "logs; use one of those.", call. = FALSE)
  if (ext == "csv")
    stop("'", basename(path), "' is a CSV file. This skill reads raw LI-6800 logs (.xlsx or ",
         "plaintext). If it is already an AgriClaw tidy CSV, pass it to claw-aci-fit, ",
         "claw-light-response or claw-stomatal-models directly.", call. = FALSE)
  stop("'", basename(path), "' has extension '.", ext, "': not an LI-6800 log (.xlsx or ",
       "plaintext without extension).", call. = FALSE)
}

read_one_licor <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
  ftype <- file_type_of(path)
  x <- tryCatch(
    PhotoGEA::read_gasex_file(path, timestamp_colname = NA, file_type = ftype,
                              instrument_type = "Licor LI-6800", check_for_zero = c()),
    error = function(e) stop("Could not read '", basename(path), "' as an LI-6800 ", ftype,
                             " log: ", conditionMessage(e), "\nUse the file as the instrument ",
                             "wrote it (do not delete header rows or rename columns).",
                             call. = FALSE))
  dead_other <- check_unevaluated(x$main_data, path, ftype)
  list(x = x, file_type = ftype, dead_columns = dead_other)
}

remark_value <- function(x, key) {
  r <- x$remarks
  if (is.data.frame(r) && key %in% names(r) && nrow(r) >= 1L) {
    v <- trimws(as.character(r[[key]][1L]))
    if (!is.na(v) && nzchar(v)) return(v)
  }
  NA_character_
}

safe_name <- function(nm) {
  nm <- gsub("[^A-Za-z0-9_.]+", "_", trimws(nm))
  ifelse(nzchar(nm), nm, "unnamed")
}

# Convert one PhotoGEA exdf to the tidy layout. Units of mapped columns are
# checked against the contract so a unit change fails loudly.
licor_to_tidy <- function(x, path, file_type) {
  md <- x$main_data
  n <- nrow(md)
  out <- data.frame(source_file = rep(basename(path), n), file_type = file_type,
                    console_serial = remark_value(x, "Console s/n"),
                    head_serial = remark_value(x, "Head s/n"),
                    stringsAsFactors = FALSE)
  units <- x$units
  for (i in seq_len(nrow(TIDY_MAP))) {
    src <- TIDY_MAP$licor[i]
    if (TIDY_MAP$tidy[i] == "tleaf_c" && !src %in% names(md) && "Tleaf" %in% names(md)) src <- "Tleaf"
    if (!src %in% names(md)) next
    want <- TIDY_MAP$units[i]
    got <- if (src %in% names(units)) as.character(units[[src]][1L]) else NA_character_
    if (!is.na(want) && !identical(got, want))
      stop("'", basename(path), "': column ", src, " has units '", got, "' but '", want,
           "' was expected. Check the instrument's unit settings before continuing.", call. = FALSE)
    v <- md[[src]]
    out[[TIDY_MAP$tidy[i]]] <- if (TIDY_MAP$type[i] == "num") suppressWarnings(as.numeric(v)) else as.character(v)
  }
  cats <- vapply(names(md), function(cn) as.character(x$categories[[cn]][1L]), "")
  for (u in names(md)[cats %in% "UserDefCon"]) {
    nm <- safe_name(u)
    if (nm %in% names(out)) nm <- paste0("uc_", nm)
    out[[nm]] <- md[[u]]
  }
  if ("user_remark" %in% names(md)) out$user_remark <- as.character(md$user_remark)
  out
}

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

# Instrument ID: an explicit column if given, else "<console s/n>/<head s/n>"
# from the file header, else "unknown:<file>".
assign_instrument <- function(df, instrument_col = "") {
  if (nzchar(instrument_col)) {
    if (!instrument_col %in% names(df))
      stop("--instrument-col '", instrument_col, "' is not a column. Available user constants: ",
           paste(setdiff(names(df), c(TIDY_MAP$tidy, "source_file", "file_type",
                                      "console_serial", "head_serial", "user_remark")),
                 collapse = ", "), call. = FALSE)
    id <- trimws(as.character(df[[instrument_col]]))
    if (any(is.na(id) | !nzchar(id)))
      stop("Column '", instrument_col, "' has empty instrument IDs in ",
           sum(is.na(id) | !nzchar(id)), " row(s). Fill them in before running.", call. = FALSE)
    df$instrument_id <- id
    df$instrument_id_source <- paste0("column ", instrument_col)
    return(df)
  }
  cs <- df$console_serial
  hs <- df$head_serial
  id <- ifelse(!is.na(cs) & !is.na(hs), paste0(cs, "/", hs), ifelse(!is.na(cs), cs, hs))
  src <- ifelse(is.na(id), "unknown", "file header serial numbers")
  id[is.na(id)] <- paste0("unknown:", df$source_file[is.na(id)])
  df$instrument_id <- id
  df$instrument_id_source <- src
  df
}

# The same observation (same instrument, same timestamp) read from more than
# one file, e.g. a log's .xlsx and plaintext versions: later copies (in input
# order) are TRUE. Checked only when the instrument is known, because two
# different instruments can log in the same second.
duplicate_observations <- function(df) {
  if (!all(c("instrument_id", "time_unix_s", "source_file") %in% names(df))) return(rep(NA, nrow(df)))
  known <- !startsWith(df$instrument_id, "unknown:") & !is.na(df$time_unix_s)
  key <- paste(df$instrument_id, format(df$time_unix_s, digits = 15))
  dup <- known & duplicated(key)
  first_file <- df$source_file[match(key, key)]
  out <- dup & first_file != df$source_file
  out[!known] <- NA
  out
}

# Per-row QC flags. NA means the check could not be made (column absent).
add_qc_flags <- function(df, tleaf_min = 5, tleaf_max = 45, min_delta_co2 = 1) {
  n <- nrow(df)
  has <- function(cols) all(cols %in% names(df))
  na_vec <- rep(NA, n)
  core <- intersect(c("a_umol_m2_s", "ci_umol_mol", "gsw_mol_m2_s"), names(df))
  df$qc_missing_values <- if (length(core))
    Reduce(`|`, lapply(core, function(cn) is.na(df[[cn]]))) else na_vec
  df$qc_ci_negative <- if (has("ci_umol_mol")) df$ci_umol_mol < 0 else na_vec
  df$qc_ci_gt_ca <- if (has(c("ci_umol_mol", "ca_umol_mol", "a_umol_m2_s")))
    df$ci_umol_mol > df$ca_umol_mol & df$a_umol_m2_s > 0 else na_vec
  df$qc_gsw_negative <- if (has("gsw_mol_m2_s")) df$gsw_mol_m2_s < 0 else na_vec
  df$qc_tleaf_out_of_range <- if (has("tleaf_c"))
    df$tleaf_c < tleaf_min | df$tleaf_c > tleaf_max else na_vec
  if (has(c("stable_n", "stable_total"))) {
    df$qc_unstable <- ifelse(df$stable_total > 0, df$stable_n < df$stable_total, NA)
  } else if (has(c("a_ok", "gsw_ok"))) {
    df$qc_unstable <- df$a_ok == 0 | df$gsw_ok == 0
  } else {
    df$qc_unstable <- na_vec
  }
  if (has(c("co2_r_umol_mol", "co2_s_umol_mol"))) {
    df$delta_co2_umol_mol <- df$co2_r_umol_mol - df$co2_s_umol_mol
    df$qc_small_delta_co2 <- abs(df$delta_co2_umol_mol) < min_delta_co2
  } else {
    df$qc_small_delta_co2 <- na_vec
  }
  df$qc_duplicate <- duplicate_observations(df)
  df$qc_hard_fail <- Reduce(`|`, lapply(HARD_FLAGS, function(cn) df[[cn]] %in% TRUE))
  m <- matrix(unlist(lapply(QC_FLAG_COLS, function(cn) df[[cn]] %in% TRUE)), nrow = n)
  labels <- sub("^qc_", "", QC_FLAG_COLS)
  df$qc_flags <- vapply(seq_len(n), function(i) paste(labels[m[i, ]], collapse = ";"), "")
  df
}

split_cols <- function(s) {
  if (!nzchar(s)) return(character())
  trimws(strsplit(s, ",", fixed = TRUE)[[1L]])
}

check_cols_exist <- function(df, cols, what) {
  miss <- setdiff(cols, names(df))
  if (length(miss))
    stop(what, " column(s) not found: ", paste(miss, collapse = ", "), ". Available columns ",
         "include: ", paste(setdiff(names(df), c(TIDY_MAP$tidy, QC_FLAG_COLS)), collapse = ", "),
         call. = FALSE)
  invisible(TRUE)
}

curve_key <- function(df, curve_cols) {
  if (!length(curve_cols)) return(rep(NA_character_, nrow(df)))
  do.call(paste, c(lapply(curve_cols, function(cn) as.character(df[[cn]])), sep = " - "))
}

# One row per curve: size, files, instruments, temperature range and flags.
curve_summary <- function(df, curve_cols, tleaf_drift = 2) {
  key <- curve_key(df, curve_cols)
  ids <- unique(key)
  rows <- lapply(ids, function(k) {
    d <- df[key == k, , drop = FALSE]
    tl <- if ("tleaf_c" %in% names(d)) d$tleaf_c else NA_real_
    rng <- if (all(is.na(tl))) NA_real_ else diff(range(tl, na.rm = TRUE))
    sp <- if ("co2_r_sp_umol_mol" %in% names(d)) paste(round(d$co2_r_sp_umol_mol), collapse = " ") else ""
    data.frame(curve = k, n_points = nrow(d),
               source_files = paste(unique(d$source_file), collapse = ";"),
               instrument_ids = paste(unique(d$instrument_id), collapse = ";"),
               tleaf_min_c = if (all(is.na(tl))) NA else min(tl, na.rm = TRUE),
               tleaf_max_c = if (all(is.na(tl))) NA else max(tl, na.rm = TRUE),
               n_hard_fail = sum(d$qc_hard_fail), n_unstable = sum(d$qc_unstable %in% TRUE),
               n_small_delta_co2 = sum(d$qc_small_delta_co2 %in% TRUE),
               spans_files = length(unique(d$source_file)) > 1L,
               spans_instruments = length(unique(d$instrument_id)) > 1L,
               tleaf_drift = !is.na(rng) && rng > tleaf_drift,
               co2_setpoints = sp, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Confounding of a design factor with instrument, counted in curves (or rows
# when no curve columns are given). Status categories:
#   single instrument / single level: nothing to confound;
#   completely confounded: no instrument measured more than one level;
#   partially confounded: some level was never measured on some instrument;
#   crossed: every level measured on every instrument (Cramer's V reported).
confounding_check <- function(df, design_cols, curve_cols = character(), v_warn = 0.5) {
  if (!length(design_cols)) return(list(summary = NULL, counts = NULL))
  unit <- if (length(curve_cols)) curve_key(df, curve_cols) else as.character(seq_len(nrow(df)))
  sums <- list()
  counts <- list()
  for (dc in design_cols) {
    lev <- as.character(df[[dc]])
    u <- unique(data.frame(unit = unit, level = lev, instrument = df$instrument_id,
                           stringsAsFactors = FALSE))
    tab <- table(u$level, u$instrument)
    counts[[dc]] <- data.frame(design_col = dc, as.data.frame(tab, stringsAsFactors = FALSE),
                               stringsAsFactors = FALSE)
    names(counts[[dc]]) <- c("design_col", "level", "instrument_id", "n_units")
    nl <- nrow(tab)
    ni <- ncol(tab)
    v <- NA_real_
    if (nl >= 2L && ni >= 2L) {
      chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE)$statistic)
      v <- as.numeric(sqrt(chi / (sum(tab) * (min(nl, ni) - 1L))))
    }
    status <- if (ni < 2L) "single instrument" else if (nl < 2L) "single level" else {
      levels_per_inst <- colSums(tab > 0)
      if (all(levels_per_inst <= 1L)) "completely confounded"
      else if (any(tab == 0)) "partially confounded"
      else "crossed"
    }
    msg <- switch(status,
      "single instrument" = "All data come from one instrument: no instrument confounding.",
      "single level" = "Only one level present: nothing to compare.",
      "completely confounded" = paste0(
        dc, " is completely confounded with instrument: no instrument measured more than one ",
        "level, so instrument differences cannot be separated from ", dc, " effects. No ",
        "statistical model can fix this; re-measure with instruments balanced across levels, ",
        "or measure the same leaves on both instruments to estimate the instrument offset."),
      "partially confounded" = paste0(
        dc, " is partially confounded with instrument: some levels were never measured on some ",
        "instruments. Include instrument as a blocking factor (fixed block or random effect) or ",
        "a covariate in the downstream model, and say that comparisons between levels measured ",
        "on different instruments rely on that adjustment."),
      "crossed" = if (!is.na(v) && v >= v_warn) paste0(
        dc, " is crossed with instrument but unevenly (Cramer's V = ", sprintf("%.2f", v),
        "). Include instrument as a blocking factor or covariate.") else paste0(
        dc, " is crossed with instrument (every level measured on every instrument; Cramer's V = ",
        sprintf("%.2f", v), "). Instrument can still be used as a blocking factor."))
    warn <- status %in% c("completely confounded", "partially confounded") ||
      (status == "crossed" && !is.na(v) && v >= v_warn)
    sums[[dc]] <- data.frame(design_col = dc, n_levels = nl, n_instruments = ni,
                             status = status, cramers_v = v, warn = warn, message = msg,
                             stringsAsFactors = FALSE)
  }
  list(summary = do.call(rbind, sums), counts = do.call(rbind, counts))
}

list_input_files <- function(input, input_dir) {
  files <- split_cols(input)
  if (nzchar(input_dir)) {
    if (!dir.exists(input_dir)) stop("--input-dir not found: ", input_dir, call. = FALSE)
    cand <- list.files(input_dir, full.names = TRUE)
    cand <- cand[!dir.exists(cand) & !startsWith(basename(cand), "~$")]
    keep <- vapply(cand, function(f) {
      ext <- tolower(tools::file_ext(f))
      if (ext == "xlsx") return(TRUE)
      if (!ext %in% c("", "txt")) return(FALSE)
      first <- tryCatch(readLines(f, n = 5L, warn = FALSE), error = function(e) character())
      any(grepl("^\\[Head", first))
    }, logical(1))
    files <- c(files, cand[keep])
  }
  files <- unique(files)
  if (!length(files))
    stop("No input files. Give --input FILE[,FILE] or --input-dir DIR (LI-6800 .xlsx or ",
         "plaintext logs), or --demo. See --help.", call. = FALSE)
  files
}

# An LI-6800 log is often saved twice: plaintext `X` and Excel `X.xlsx`.
# Keep the plaintext version (it stores computed values); report the skip.
drop_xlsx_twins <- function(files) {
  stem <- sub("[.](xlsx|txt)$", "", files, ignore.case = TRUE)
  is_x <- grepl("[.]xlsx$", files, ignore.case = TRUE)
  twin <- is_x & (stem %in% stem[!is_x])
  structure(files[!twin], skipped = files[twin])
}

read_licor_files <- function(files) {
  parts <- list()
  meta <- list()
  for (f in files) {
    r <- read_one_licor(f)
    td <- licor_to_tidy(r$x, f, r$file_type)
    parts[[f]] <- td
    meta[[f]] <- data.frame(
      source_file = basename(f), file_type = r$file_type, n_rows = nrow(td),
      console_serial = remark_value(r$x, "Console s/n"),
      head_serial = remark_value(r$x, "Head s/n"),
      file_opened = remark_value(r$x, "File opened"),
      console_version = remark_value(r$x, "Console ver"),
      chamber_type = remark_value(r$x, "Chamber type"),
      all_zero_columns = paste(r$dead_columns, collapse = ";"),
      stringsAsFactors = FALSE)
  }
  list(data = rbind_fill(parts), files = do.call(rbind, meta))
}

process_licor <- function(files, curve_cols = character(), design_cols = character(),
                          instrument_col = "", tleaf_min = 5, tleaf_max = 45,
                          min_delta_co2 = 1, tleaf_drift = 2) {
  rd <- read_licor_files(files)
  df <- rd$data
  check_cols_exist(df, curve_cols, "--curve-cols")
  check_cols_exist(df, design_cols, "--design-cols")
  df <- assign_instrument(df, instrument_col)
  df <- add_qc_flags(df, tleaf_min, tleaf_max, min_delta_co2)
  if (length(curve_cols)) df$curve_id <- curve_key(df, curve_cols)
  first <- c("source_file", "instrument_id", if (length(curve_cols)) "curve_id", curve_cols, design_cols)
  df <- df[, c(unique(first), setdiff(names(df), first)), drop = FALSE]
  curves <- if (length(curve_cols)) curve_summary(df, curve_cols, tleaf_drift) else NULL
  conf <- confounding_check(df, design_cols, curve_cols)
  rd$files$instrument_id <- vapply(rd$files$source_file, function(f)
    paste(unique(df$instrument_id[df$source_file == f]), collapse = ";"), "")
  missing_uc <- list()
  uc <- setdiff(names(df), c(TIDY_MAP$tidy, "source_file", "file_type", "console_serial",
                             "head_serial", "user_remark", "instrument_id",
                             "instrument_id_source", "curve_id", "delta_co2_umol_mol",
                             QC_FLAG_COLS, "qc_hard_fail", "qc_flags"))
  for (u in uc) {
    absent <- unique(df$source_file[is.na(df[[u]])])
    if (length(absent) && length(absent) < length(unique(df$source_file)))
      missing_uc[[u]] <- absent
  }
  list(data = df, files = rd$files, curves = curves, confounding = conf,
       missing_user_constants = missing_uc)
}

flag_table <- function(df) {
  by_file <- lapply(split(df, df$source_file), function(d) {
    data.frame(source_file = d$source_file[1L], instrument_id = paste(unique(d$instrument_id), collapse = ";"),
               n_rows = nrow(d),
               t(vapply(c(QC_FLAG_COLS, "qc_hard_fail"), function(cn) sum(d[[cn]] %in% TRUE), 0)),
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  do.call(rbind, by_file)
}

plot_qc <- function(df, out_dir) {
  open_png(out_dir, "qc_overview.png", width = 10, height = 4.8)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  inst <- factor(df$instrument_id)
  cols <- grDevices::hcl.colors(max(2L, nlevels(inst)), "Dark 3")[as.integer(inst)]
  if (all(c("ci_umol_mol", "a_umol_m2_s") %in% names(df))) {
    ok <- is.finite(df$ci_umol_mol) & is.finite(df$a_umol_m2_s)
    graphics::plot(df$ci_umol_mol[ok], df$a_umol_m2_s[ok], col = cols[ok], pch = 16, cex = 0.8,
                   xlab = "Ci (umol/mol)", ylab = "A (umol m-2 s-1)",
                   main = "All points (x = hard QC failure)")
    bad <- ok & df$qc_hard_fail
    if (any(bad)) graphics::points(df$ci_umol_mol[bad], df$a_umol_m2_s[bad], pch = 4, cex = 1.6, lwd = 2)
    graphics::legend("bottomright", legend = levels(inst), col = grDevices::hcl.colors(max(2L, nlevels(inst)), "Dark 3")[seq_len(nlevels(inst))],
                     pch = 16, bty = "n", cex = 0.7, title = "instrument")
  } else graphics::plot.new()
  if ("tleaf_c" %in% names(df)) {
    graphics::plot(seq_len(nrow(df)), df$tleaf_c, col = cols, pch = 16, cex = 0.7,
                   xlab = "Row (files in input order)", ylab = "Leaf temperature (C)",
                   main = "Leaf temperature by row")
  } else graphics::plot.new()
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

build_warnings <- function(res, opt) {
  df <- res$data
  w <- character()
  dead <- res$files$all_zero_columns
  if (any(nzchar(dead)))
    w <- c(w, sprintf("File %s: computed column(s) %s are all zero or empty; check them before use.",
                      res$files$source_file[nzchar(dead)], dead[nzchar(dead)]))
  unk <- unique(df$source_file[startsWith(df$instrument_id, "unknown:")])
  if (length(unk))
    w <- c(w, sprintf(paste0("No serial number in %s (Remarks sheet missing?): instrument ID set to ",
                             "'unknown:<file>', and duplicate observations cannot be checked. Pass ",
                             "--instrument-col if a user constant records the instrument."),
                      paste(unk, collapse = ", ")))
  if (length(res$skipped_twins))
    w <- c(w, sprintf("Skipped %s: the plaintext version of the same log was also given (it stores the computed values).",
                      paste(basename(res$skipped_twins), collapse = ", ")))
  dup <- df$qc_duplicate %in% TRUE
  if (any(dup))
    w <- c(w, sprintf(paste0("%d row(s) in %s repeat observations already read from another file (same ",
                             "instrument and timestamp): the same log was probably read twice. They are ",
                             "flagged qc_duplicate (a hard failure), so fits use the first copy only."),
                      sum(dup), paste(unique(df$source_file[dup]), collapse = ", ")))
  if (!is.null(res$confounding$summary)) {
    cs <- res$confounding$summary
    w <- c(w, cs$message[cs$warn])
  }
  n_hard <- sum(df$qc_hard_fail)
  if (n_hard) w <- c(w, sprintf("%d of %d points fail a hard QC check (missing values, Ci < 0, Ci > Ca with A > 0, gsw < 0, leaf temperature outside %g-%g C, or a duplicate of an observation in another file). They are kept in the CSV with qc_hard_fail = TRUE; downstream skills exclude them from fits by default.",
                                n_hard, nrow(df), opt$tleaf_min, opt$tleaf_max))
  n_un <- sum(df$qc_unstable %in% TRUE)
  if (n_un) w <- c(w, sprintf("%d point(s) were logged before the instrument's stability criteria were met (qc_unstable).", n_un))
  if (!is.null(res$curves)) {
    cv <- res$curves
    if (any(cv$spans_files)) w <- c(w, sprintf("Curve ID(s) %s span more than one file: the curve columns do not identify single curves. Add source_file or a date to --curve-cols.",
                                               paste(cv$curve[cv$spans_files], collapse = "; ")))
    if (any(cv$tleaf_drift)) w <- c(w, sprintf("Leaf temperature drifted by more than %g C within curve(s): %s.",
                                               opt$tleaf_drift, paste(cv$curve[cv$tleaf_drift], collapse = "; ")))
  }
  for (u in names(res$missing_user_constants))
    w <- c(w, sprintf("User constant '%s' is missing from file(s) %s (left empty in the CSV).", u,
                      paste(res$missing_user_constants[[u]], collapse = ", ")))
  if ("oxygen_pct" %in% names(df) && any(abs(df$oxygen_pct - 21) > 0.5, na.rm = TRUE))
    w <- c(w, sprintf("Oxygen is not 21%% in %d row(s) (values: %s). A-Ci fits use this column; confirm it matches the gas actually supplied.",
                      sum(abs(df$oxygen_pct - 21) > 0.5, na.rm = TRUE),
                      paste(sort(unique(df$oxygen_pct)), collapse = ", ")))
  w
}

write_report <- function(res, opt, files, out_dir) {
  df <- res$data
  w <- build_warnings(res, opt)
  ft <- res$files[, c("source_file", "file_type", "n_rows", "instrument_id", "console_serial",
                      "head_serial", "file_opened")]
  fl <- flag_table(df)
  lines <- c(sprintf("# %s report", SKILL), "",
             if (length(w)) c("## Warnings", "", paste("-", w), "") else c("## Warnings", "", "None.", ""),
             "## Files read", "",
             sprintf("%d file(s), %d rows. Instrument IDs from: %s.", nrow(res$files), nrow(df),
                     paste(unique(df$instrument_id_source), collapse = "; ")), "",
             fmt_table(ft), "",
             "## QC flags (number of rows)", "",
             "Hard failures (`qc_hard_fail`): missing A/Ci/gsw, Ci < 0, Ci > Ca while A > 0, gsw < 0, leaf temperature out of range, duplicate of an observation already read from another file. Cautions: `qc_unstable`, `qc_small_delta_co2`. Thresholds are listed below; all are adjustable rules of thumb.", "",
             fmt_table(fl), "",
             sprintf("Thresholds: leaf temperature %g-%g C; |CO2_r - CO2_s| < %g umol/mol; curve temperature drift > %g C.",
                     opt$tleaf_min, opt$tleaf_max, opt$min_delta_co2, opt$tleaf_drift), "")
  if (!is.null(res$confounding$summary)) {
    lines <- c(lines, "## Instrument confounding", "",
               fmt_table(res$confounding$summary[, c("design_col", "n_levels", "n_instruments", "status", "cramers_v")]), "",
               "Counts (curves, or rows if no --curve-cols) per level and instrument are in `tables/confounding_counts.csv`.", "")
  }
  if (!is.null(res$curves)) {
    cv <- res$curves
    lines <- c(lines, "## Curves", "",
               sprintf("%d curve(s) identified by %s. Per-curve details: `tables/curve_summary.csv`.",
                       nrow(cv), opt$curve_cols), "",
               fmt_table(cv[, c("curve", "n_points", "instrument_ids", "tleaf_min_c", "tleaf_max_c",
                                "n_hard_fail", "n_unstable")]), "")
  }
  lines <- c(lines, "## Output", "",
             "- `tables/licor_tidy.csv`: one row per logged point; units are in the column names. Read by claw-aci-fit, claw-light-response and claw-stomatal-models.",
             "- `tables/file_summary.csv`, `tables/qc_flag_counts.csv`, `tables/curve_summary.csv`, `tables/confounding_*.csv`.",
             "- `figures/qc_overview.png`: A vs Ci (hard failures marked) and leaf temperature.", "",
             "## Method", "",
             sprintf("Files were read with PhotoGEA %s (`read_gasex_file`), which handles Excel and plaintext LI-6800 logs, files closed and reopened, and changes of oxygen during a log. Computed columns of Excel files were checked for unevaluated formulas. Instrument IDs combine the console and head serial numbers from the file header unless --instrument-col is given.",
                     as.character(utils::packageVersion("PhotoGEA"))), "",
             "## References", "",
             "- Lochocki EB, Salesse-Smith CE, McGrath JM (2025) PhotoGEA: an R package for closer fitting of photosynthetic gas exchange data with non-Gaussian confidence interval estimation. Plant, Cell & Environment 48: 5104-5119. https://doi.org/10.1111/pce.15501",
             "- Long SP, Bernacchi CJ (2003) Gas exchange measurements, what can they tell us about the underlying limitations to photosynthesis? Procedures and sources of error. Journal of Experimental Botany 54: 2393-2401. https://doi.org/10.1093/jxb/erg262",
             "- von Caemmerer S, Farquhar GD (1981) Some relationships between the biochemistry of photosynthesis and the gas exchange of leaves. Planta 153: 376-387. https://doi.org/10.1007/BF00384257",
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
    files <- vapply(c("c3_aci_1.xlsx", "c3_aci_2.xlsx"), function(f)
      system.file("extdata", f, package = "PhotoGEA"), "")
    if (any(!nzchar(files))) stop("PhotoGEA example files not found.", call. = FALSE)
    if (!nzchar(opt$curve_cols)) opt$curve_cols <- "instrument,species,plot"
    if (!nzchar(opt$design_cols)) opt$design_cols <- "species"
    if (!nzchar(opt$instrument_col)) opt$instrument_col <- "instrument"
  } else {
    files <- drop_xlsx_twins(list_input_files(opt$input, opt$input_dir))
  }
  skipped <- attr(files, "skipped")
  files <- as.character(files)
  if (opt$tleaf_min >= opt$tleaf_max) stop("--tleaf-min must be below --tleaf-max.", call. = FALSE)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

  res <- process_licor(files, split_cols(opt$curve_cols), split_cols(opt$design_cols),
                       opt$instrument_col, opt$tleaf_min, opt$tleaf_max, opt$min_delta_co2,
                       opt$tleaf_drift)
  res$skipped_twins <- skipped
  save_table(res$data, opt$out, "licor_tidy.csv")
  save_table(res$files, opt$out, "file_summary.csv")
  save_table(flag_table(res$data), opt$out, "qc_flag_counts.csv")
  if (!is.null(res$curves)) save_table(res$curves, opt$out, "curve_summary.csv")
  if (!is.null(res$confounding$summary)) {
    save_table(res$confounding$summary, opt$out, "confounding_summary.csv")
    save_table(res$confounding$counts, opt$out, "confounding_counts.csv")
  }
  plot_qc(res$data, opt$out)
  w <- write_report(res, opt, files, opt$out)
  write_repro_bundle(opt$out, argv, inputs = files)
  cat(sprintf("Done: %d file(s), %d rows, %d warning(s). Report: %s\n", length(files),
              nrow(res$data), length(w),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(res)
}

if (sys.nframe() == 0L) main()
