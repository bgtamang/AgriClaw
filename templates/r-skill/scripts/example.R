#!/usr/bin/env Rscript
# claw-example: summarize a numeric trait by group (n, mean, SD, CV) and flag
# groups with too few observations. Template for AgriClaw R skills: copy this
# file, keep the helper section unchanged, and replace the analysis section.
#
# Usage:
#   Rscript scripts/example.R --demo
#   Rscript scripts/example.R --input data.csv --trait yield_kg_ha --group genotype
#   Rscript scripts/example.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT license.

SKILL   <- "claw-example"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/example.R [--demo] [--input FILE.csv] [--trait COL]
                               [--group COL] [--min-n N] [--out DIR]

  --demo     run on demo/example.csv
  --input    CSV file, one row per observation
  --trait    numeric column to summarize (default: yield_kg_ha)
  --group    grouping column (default: genotype)
  --min-n    flag groups with fewer observations than this (default: 3)
  --out      output folder (default: agriclaw_output/claw-example)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, input = "", trait = "yield_kg_ha",
                 group = "genotype", min_n = 3, out = "")

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
# Analysis (replace in a real skill)
# ---------------------------------------------------------------------------

read_input <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

check_columns <- function(dat, trait, group) {
  missing <- setdiff(c(trait, group), names(dat))
  if (length(missing))
    stop("Column(s) not found: ", paste(missing, collapse = ", "),
         ". Available: ", paste(names(dat), collapse = ", "), call. = FALSE)
  if (!is.numeric(dat[[trait]]))
    stop("Column '", trait, "' must be numeric.", call. = FALSE)
  invisible(TRUE)
}

summarize_trait <- function(dat, trait, group, min_n = 3) {
  check_columns(dat, trait, group)
  groups <- sort(unique(dat[[group]]))
  rows <- lapply(groups, function(g) {
    x <- dat[[trait]][dat[[group]] == g]
    x <- x[!is.na(x)]
    n <- length(x)
    m <- if (n > 0) mean(x) else NA_real_
    s <- if (n > 1) stats::sd(x) else NA_real_
    data.frame(group = g, n = n, mean = m, sd = s, cv_pct = 100 * s / m,
               low_n = n < min_n, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, rows)
  names(res)[1] <- group
  res
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (opt$demo) opt$input <- file.path(SKILL_DIR, "demo", "example.csv")
  if (!nzchar(opt$input)) stop("Provide --input FILE.csv, or --demo. See --help.", call. = FALSE)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

  dat <- read_input(opt$input)
  check_columns(dat, opt$trait, opt$group)
  n_na <- sum(is.na(dat[[opt$trait]]))
  res <- summarize_trait(dat, opt$trait, opt$group, opt$min_n)
  save_table(res, opt$out, "summary.csv")

  open_png(opt$out, "means.png")
  ord <- order(res$mean)
  graphics::dotchart(res$mean[ord], labels = res[[opt$group]][ord],
                     xlab = opt$trait, main = "Group means")
  grDevices::dev.off()

  warn <- character()
  if (n_na > 0) warn <- c(warn, sprintf("%d missing value(s) in '%s' were excluded.", n_na, opt$trait))
  if (any(res$low_n)) warn <- c(warn, sprintf("%d group(s) have fewer than %d observations.",
                                              sum(res$low_n), opt$min_n))
  report <- c(sprintf("# %s report", SKILL), "",
              if (length(warn)) c("## Warnings", "", paste("-", warn), "") else NULL,
              "## Summary", "",
              sprintf("Input: `%s`  ", opt$input),
              sprintf("Groups: %d; observations used: %d", nrow(res), sum(res$n)), "",
              "| group | n | mean | sd | cv (%) |", "|---|---|---|---|---|",
              sprintf("| %s | %d | %.3g | %.3g | %.3g |", res[[opt$group]], res$n,
                      res$mean, res$sd, res$cv_pct),
              "", "![Group means](figures/means.png)", "",
              sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(report, file.path(opt$out, "report.md"))
  write_repro_bundle(opt$out, argv, inputs = opt$input)
  cat(sprintf("Done: %d groups. Report: %s\n", nrow(res),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(res)
}

if (sys.nframe() == 0L) main()
