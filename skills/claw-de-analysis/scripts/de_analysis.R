#!/usr/bin/env Rscript
# claw-de-analysis: differential expression with DESeq2. Design formula with batch,
# explicit reference levels, Wald contrasts or likelihood-ratio tests (time courses,
# multi-level factors), independent filtering, BH-adjusted p-values, and apeglm
# log2 fold-change shrinkage for ranking and plots only. Refuses designs with
# unreplicated groups or batch completely confounded with condition, and inputs
# that are TPM/FPKM instead of counts.
#
# Usage:
#   Rscript scripts/de_analysis.R --demo
#   Rscript scripts/de_analysis.R --txi txi.rds --samples samples.csv \
#       --design "~ batch + condition" --reference condition=control
#   Rscript scripts/de_analysis.R --counts counts.csv --samples samples.csv \
#       --design "~ genotype + time + genotype:time" --test LRT --reduced "~ genotype + time"
#   Rscript scripts/de_analysis.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-de-analysis"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/de_analysis.R [--demo] (--txi TXI.rds | --counts COUNTS.csv [--lengths LENGTHS.csv])
         --samples SAMPLES.csv --design \"~ batch + condition\" [options]

  --demo              run on the bundled simulated counts (demo/)
  --txi               tximport gene-level list saved by claw-rnaseq-import (txi.rds)
  --counts            CSV of raw integer counts: first column gene_id, one column per sample
  --lengths           CSV of average transcript lengths matching --counts (tximport offsets)
  --samples           sample sheet CSV with a 'sample' column and the design variables
  --design            design formula, e.g. \"~ batch + condition\" (variable of interest last)
  --reference         reference levels, e.g. condition=control,genotype=WT
  --contrast          factor,numerator,denominator (e.g. condition,drought,control) or a
                      coefficient name; several separated by ';'. Default: every level of
                      the last design variable against its reference level
  --test              Wald or LRT (default: Wald)
  --reduced           reduced formula for --test LRT, e.g. \"~ batch\"
  --alpha             FDR level for calling genes and for independent filtering (default: 0.05)
  --lfc-threshold     test |log2 fold change| > this value instead of 0 (default: 0)
  --shrink            apeglm or none: LFC shrinkage for ranking and plots (default: apeglm)
  --min-count         pre-filter: keep genes with at least this count... (default: 10)
  --min-samples       ...in at least this many samples (default: 0 = smallest group size)
  --min-reps          warn when a compared group has fewer replicates (default: 3)
  --exclude           comma-separated sample names to leave out (reported)
  --as-factor         comma-separated numeric columns to treat as factors (e.g. time)
  --no-length-offset  with --txi: ignore transcript-length offsets (3' tag-seq libraries)
  --allow-noninteger  round non-integer counts (only for tximport scaledTPM/lengthScaledTPM counts)
  --out               output folder (default: agriclaw_output/claw-de-analysis)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, txi = "", counts = "", lengths = "", samples = "",
                 design = "", reference = "", contrast = "", test = "Wald", reduced = "",
                 alpha = 0.05, lfc_threshold = 0, shrink = "apeglm", min_count = 10,
                 min_samples = 0, min_reps = 3, exclude = "", as_factor = "",
                 no_length_offset = FALSE, allow_noninteger = FALSE, out = "")

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

fail <- function(...) stop(paste0(...), call. = FALSE)
split_list <- function(x, sep = ",") {
  if (!nzchar(x)) return(character())
  v <- trimws(strsplit(x, sep, fixed = TRUE)[[1L]])
  v[nzchar(v)]
}

read_matrix_csv <- function(path, what) {
  if (!file.exists(path)) fail(what, " file not found: ", path)
  d <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(d) < 2L) fail(what, " file needs a gene_id column plus one column per sample.")
  ids <- as.character(d[[1L]])
  if (anyDuplicated(ids)) fail("Duplicated gene IDs in ", path, ", e.g. ", ids[duplicated(ids)][1])
  m <- as.matrix(d[, -1L, drop = FALSE])
  if (!is.numeric(m))
    fail(what, " columns 2..n must all be numeric. Non-numeric column(s): ",
         paste(names(d)[-1L][!vapply(d[-1L], is.numeric, TRUE)], collapse = ", "))
  rownames(m) <- ids
  m
}

# Read counts (and optional length offsets). Returns list(counts, lengths, source, cfa).
read_counts_input <- function(txi = "", counts = "", lengths = "", no_length_offset = FALSE) {
  if (nzchar(txi) && nzchar(counts)) fail("Give either --txi or --counts, not both.")
  if (nzchar(txi)) {
    if (!file.exists(txi)) fail("txi file not found: ", txi)
    x <- readRDS(txi)
    if (!is.list(x) || is.null(x$counts))
      fail("--txi must be the tximport list saved by claw-rnaseq-import (txi.rds).")
    cfa <- x$countsFromAbundance %||% "no"
    len <- if (identical(cfa, "no") && !no_length_offset) x$length else NULL
    if (identical(cfa, "no") && !no_length_offset && is.null(len))
      fail("The txi object has no 'length' matrix; re-run claw-rnaseq-import.")
    return(list(counts = as.matrix(x$counts), lengths = len, source = "tximport",
                cfa = cfa, noninteger_ok = TRUE))
  }
  if (!nzchar(counts)) fail("Give --txi TXI.rds or --counts COUNTS.csv (or --demo). See --help.")
  m <- read_matrix_csv(counts, "Counts")
  len <- NULL
  if (nzchar(lengths)) {
    len <- read_matrix_csv(lengths, "Lengths")
    if (!identical(dim(len), dim(m)) || !identical(rownames(len), rownames(m)) ||
        !identical(colnames(len), colnames(m)))
      fail("--lengths must have the same genes and samples, in the same order, as --counts.")
    if (any(!is.finite(len)) || any(len <= 0)) fail("--lengths must be positive numbers.")
  }
  list(counts = m, lengths = len, source = if (is.null(len)) "count matrix" else
         "count matrix + length offsets", cfa = "no", noninteger_ok = !is.null(len))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Reject TPM/FPKM/CPM/log values: DESeq2 needs raw (estimated) counts.
check_count_matrix <- function(m, noninteger_ok = FALSE, allow_noninteger = FALSE) {
  if (anyNA(m)) fail("The count matrix has missing values (NA). DESeq2 needs complete counts; ",
                     "a gene not detected in a sample has count 0, not NA.")
  if (any(m < 0)) fail("The count matrix has negative values: these are not counts ",
                       "(log-transformed or batch-corrected data?). Provide raw counts.")
  cs <- colSums(m)
  frac_nonint <- mean(abs(m - round(m)) > 1e-8)
  per_million <- all(abs(cs - 1e6) / 1e6 < 0.01)
  if (per_million && frac_nonint > 0)
    fail("Every sample sums to about 1,000,000: these look like TPM (or CPM) values, not counts. ",
         "DESeq2 must be given raw counts; TPM removes the information about sequencing depth ",
         "that the count model needs. Use tximport estimated counts (txi.rds from ",
         "claw-rnaseq-import) or featureCounts/HTSeq counts.")
  if (per_million)
    return("Every sample has about 1,000,000 counts: if the counts were rarefied (downsampled), use the original counts; DESeq2 normalises library size itself.")
  if (frac_nonint > 0 && !noninteger_ok && !allow_noninteger) {
    logish <- max(m) < 50
    fail(sprintf("%.0f%% of values are not whole numbers. ", 100 * frac_nonint),
         "DESeq2 needs raw counts, not FPKM/RPKM, TPM, CPM or log values",
         if (logish) " (values below 50 suggest a log scale)" else "", ". ",
         "If these are tximport estimated counts, pass --txi txi.rds (or --counts with ",
         "--lengths gene_lengths.csv) so the length offsets are kept; if they are tximport ",
         "scaledTPM/lengthScaledTPM counts, add --allow-noninteger.")
  }
  if (max(cs) / max(min(cs), 1) > 20)
    return(sprintf("Library sizes differ more than 20-fold (%.3g to %.3g counts).", min(cs), max(cs)))
  character()
}

read_samples <- function(path, count_samples, exclude = character()) {
  if (!nzchar(path)) fail("Give --samples SAMPLES.csv with a 'sample' column and the design variables.")
  if (!file.exists(path)) fail("Sample sheet not found: ", path)
  s <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
  if (!"sample" %in% names(s)) fail("The sample sheet needs a 'sample' column.")
  s$sample <- as.character(s$sample)
  if (anyDuplicated(s$sample)) fail("Duplicated sample names in the sample sheet.")
  bad_ex <- setdiff(exclude, s$sample)
  if (length(bad_ex)) fail("--exclude names not in the sample sheet: ", paste(bad_ex, collapse = ", "))
  s <- s[!s$sample %in% exclude, , drop = FALSE]
  count_samples <- setdiff(count_samples, exclude)
  miss <- setdiff(s$sample, count_samples)
  if (length(miss)) fail("Samples in the sheet but not in the counts: ", paste(miss, collapse = ", "))
  extra <- setdiff(count_samples, s$sample)
  if (length(extra))
    fail("Count columns with no row in the sample sheet: ", paste(extra, collapse = ", "),
         ". Add them to the sheet, or leave them out explicitly with --exclude.")
  rownames(s) <- s$sample
  s
}

parse_formula <- function(x, what = "--design") {
  if (!nzchar(x)) fail(what, " is required, e.g. \"~ batch + condition\" (variable of interest last).")
  f <- tryCatch(stats::as.formula(x), error = function(e) NULL)
  if (is.null(f) || length(f) != 2L) fail(what, " must be a one-sided formula such as \"~ batch + condition\".")
  f
}

# Make design variables factors/numerics, set reference levels. Returns list(samples, notes, warn).
prepare_design <- function(samples, design, reference = list(), as_factor = character()) {
  vars <- all.vars(design)
  miss <- setdiff(vars, names(samples))
  if (length(miss)) fail("Design variable(s) not in the sample sheet: ", paste(miss, collapse = ", "),
                         ". Columns: ", paste(names(samples), collapse = ", "))
  warn <- character()
  ref_used <- character()
  for (v in vars) {
    x <- samples[[v]]
    if (anyNA(x) || any(as.character(x) == ""))
      fail("Design variable '", v, "' has missing values. Fill them in or drop those samples with --exclude.")
    if (v %in% as_factor || is.character(x) || is.logical(x) || is.factor(x)) {
      samples[[v]] <- factor(as.character(x))
    } else if (is.numeric(x) && length(unique(x)) <= 6) {
      warn <- c(warn, sprintf(paste0("'%s' is numeric with %d distinct values, so it is modelled as a ",
                                     "continuous (linear) effect. If it is a group label (time point, ",
                                     "block), re-run with --as-factor %s."), v, length(unique(x)), v))
    }
  }
  for (v in names(reference)) {
    if (!v %in% vars) fail("--reference names '", v, "', which is not in the design.")
    if (!is.factor(samples[[v]])) fail("--reference: '", v, "' is numeric, not a factor.")
    if (!reference[[v]] %in% levels(samples[[v]]))
      fail("--reference ", v, "=", reference[[v]], ": level not found. Levels: ",
           paste(levels(samples[[v]]), collapse = ", "))
    samples[[v]] <- stats::relevel(samples[[v]], ref = reference[[v]])
  }
  for (v in vars) if (is.factor(samples[[v]]))
    ref_used <- c(ref_used, sprintf("%s = %s (%s)", v, levels(samples[[v]])[1],
                                    if (v %in% names(reference)) "given" else "alphabetical"))
  list(samples = samples, warn = warn, references = ref_used)
}

parse_reference <- function(x) {
  items <- split_list(x, ",")
  out <- list()
  for (it in items) {
    kv <- strsplit(it, "=", fixed = TRUE)[[1L]]
    if (length(kv) != 2L) fail("--reference entries look like condition=control; got '", it, "'.")
    out[[trimws(kv[1])]] <- trimws(kv[2])
  }
  out
}

# Contrasts: list of list(type = "triple"|"name", factor, num, den, label)
parse_contrasts <- function(x) {
  lapply(split_list(x, ";"), function(it) {
    p <- trimws(strsplit(it, ",", fixed = TRUE)[[1L]])
    if (length(p) == 3L) return(list(type = "triple", factor = p[1], num = p[2], den = p[3],
                                     label = paste0(p[1], "_", p[2], "_vs_", p[3])))
    if (length(p) == 1L) return(list(type = "name", name = p[1], label = p[1]))
    fail("--contrast entries are factor,numerator,denominator or one coefficient name; got '", it, "'.")
  })
}

term_vars <- function(design) {
  tl <- attr(stats::terms(design), "term.labels")
  stats::setNames(lapply(tl, function(t) all.vars(stats::as.formula(paste("~", t)))), tl)
}

default_contrasts <- function(samples, design, reference) {
  tl <- attr(stats::terms(design), "term.labels")
  if (!length(tl)) fail("The design has no terms to test.")
  last <- tl[length(tl)]
  if (grepl(":", last, fixed = TRUE))
    fail("The last design term is an interaction (", last, "). Give the comparison explicitly with ",
         "--contrast (a coefficient name) or use --test LRT with --reduced.")
  v <- last
  x <- samples[[v]]
  if (!is.factor(x)) return(list(list(type = "name", name = v, label = v)))
  if (!v %in% names(reference))
    fail("Set the reference (control) level of '", v, "' explicitly, e.g. --reference ", v, "=",
         levels(x)[1], ", or give --contrast ", v, ",TREATMENT,CONTROL. Without it R uses ",
         "alphabetical order, which can silently flip the direction of every fold change.")
  lapply(levels(x)[-1], function(l) list(type = "triple", factor = v, num = l, den = levels(x)[1],
                                         label = paste0(v, "_", l, "_vs_", levels(x)[1])))
}

# Which factor variables are tested (for replicate checks and pre-filter group size).
tested_vars <- function(contrasts, design, reduced, test, samples) {
  if (test == "LRT") {
    full_t <- term_vars(design)
    red_t <- names(term_vars(reduced))
    dropped <- setdiff(names(full_t), red_t)
    if (!length(dropped)) fail("--reduced must drop at least one term of --design.")
    return(full_t[dropped])
  }
  out <- list()
  for (ct in contrasts) {
    if (ct$type == "triple") {
      out[[ct$label]] <- structure(ct$factor, levels = c(ct$num, ct$den))  # only these levels are compared
    } else {
      mm <- stats::model.matrix(design, samples)
      cn <- make.names(colnames(mm))
      tl <- attr(stats::terms(design), "term.labels")
      a <- attr(mm, "assign")
      j <- which(cn == make.names(ct$name) | colnames(mm) == ct$name)
      out[[ct$label]] <- if (length(j) && a[j[1]] > 0) term_vars(design)[[tl[a[j[1]]]]] else character()
    }
  }
  out
}

# Replicates, rank, confounding. Returns warnings; stops on fatal problems.
check_design <- function(samples, design, tested, min_reps = 3) {
  warn <- character()
  n <- nrow(samples)
  vars <- all.vars(design)
  fac <- vars[vapply(vars, function(v) is.factor(samples[[v]]), TRUE)]
  for (v in fac) {
    samples[[v]] <- droplevels(samples[[v]])
    if (nlevels(samples[[v]]) < 2)
      fail("Design factor '", v, "' has a single level (", levels(samples[[v]]), "); remove it from --design.")
  }
  # tested groups need biological replicates
  for (nm in names(tested)) {
    tv <- intersect(tested[[nm]], fac)
    if (!length(tv)) next
    cell <- interaction(samples[tv], drop = TRUE, sep = ":")
    tab <- table(cell)
    lv <- attr(tested[[nm]], "levels")
    if (!is.null(lv) && length(tv) == 1L) {
      other <- tab[!names(tab) %in% lv & tab < 2]
      if (length(other))
        warn <- c(warn, sprintf("%s: level(s) %s of '%s' have a single sample; they add nothing to the variance estimate.",
                                nm, paste(names(other), collapse = ", "), tv))
      tab <- tab[names(tab) %in% lv]
    }
    if (any(tab < 2))
      fail("No biological replicates for ", paste(tv, collapse = ":"), " = ",
           paste(names(tab)[tab < 2], collapse = ", "), " (n = 1) in ", nm, ". DESeq2 cannot ",
           "estimate within-group variability for an unreplicated group; at least 2 replicates ",
           "are required and 3 or more are recommended.")
    low <- tab[tab < min_reps]
    if (length(low))
      warn <- c(warn, sprintf(paste0("%s: group(s) %s have fewer than %d replicates. With 3 replicates ",
                                     "most DE tools find only 20-40%% of the DE genes found with 42 ",
                                     "(Schurch et al. 2016), who recommend at least 6."),
                              nm, paste(sprintf("%s (n = %d)", names(low), as.integer(low)), collapse = ", "),
                              min_reps))
  }
  # confounding
  mm <- stats::model.matrix(design, samples)
  q <- qr(mm)
  if (q$rank < ncol(mm)) {
    why <- character()
    if (length(fac) >= 2) for (i in seq_along(fac)) for (j in seq_along(fac)) {
      if (j <= i) next
      tab <- table(samples[[fac[i]]], samples[[fac[j]]])
      row_one <- all(rowSums(tab > 0) == 1)
      col_one <- all(colSums(tab > 0) == 1)
      if (row_one && col_one) {
        why <- c(why, sprintf(paste0("'%s' and '%s' are completely confounded: each level of one occurs ",
                                     "with exactly one level of the other, so their effects cannot be ",
                                     "separated. A batch that coincides with the treatment makes the ",
                                     "treatment effect untestable; no statistical fix exists."),
                              fac[i], fac[j]))
      } else if (row_one || col_one) {
        inner <- if (row_one) fac[i] else fac[j]
        outer <- if (row_one) fac[j] else fac[i]
        why <- c(why, sprintf(paste0("'%s' is nested within '%s' (each %s level belongs to one %s level). ",
                                     "See the DESeq2 vignette section on individuals nested within groups."),
                              inner, outer, inner, outer))
      }
    }
    tl <- attr(stats::terms(design), "term.labels")
    for (t in tl[grepl(":", tl, fixed = TRUE)]) {
      tv <- intersect(term_vars(design)[[t]], fac)
      full <- prod(vapply(tv, function(v) nlevels(samples[[v]]), 1))
      have <- nlevels(interaction(samples[tv], drop = TRUE))
      if (have < full) why <- c(why, sprintf("Interaction %s has %d of %d level combinations with no samples.",
                                             t, full - have, full))
    }
    alias <- colnames(mm)[q$pivot[(q$rank + 1L):ncol(mm)]]
    fail("The design matrix is not full rank, so the model cannot be fitted as specified ",
         "(aliased column(s): ", paste(alias, collapse = ", "), ").\n  ",
         if (length(why)) paste(why, collapse = "\n  ") else
           "Some design columns are linear combinations of others.")
  }
  df_resid <- n - q$rank
  if (df_resid < 1)
    fail("No residual degrees of freedom (", n, " samples, ", q$rank, " coefficients): there are ",
         "no replicates left to estimate dispersion. Simplify the design or add replicates.")
  # partial confounding: a level of one factor missing within another
  if (length(fac) >= 2) for (i in seq_along(fac)) for (j in seq_along(fac)) {
    if (j <= i) next
    tab <- table(samples[[fac[i]]], samples[[fac[j]]])
    if (any(tab == 0))
      warn <- c(warn, sprintf(paste0("'%s' and '%s' are partly confounded (some combinations have no ",
                                     "samples): estimates rely on the other combinations and lose precision."),
                              fac[i], fac[j]))
  }
  list(warn = warn, df_resid = df_resid)
}

smallest_group <- function(samples, tested) {
  sizes <- unlist(lapply(tested, function(tv) {
    lv <- attr(tv, "levels")
    tv <- tv[vapply(tv, function(v) is.factor(samples[[v]]), TRUE)]
    if (!length(tv)) return(NULL)
    tab <- table(interaction(samples[tv], drop = TRUE))
    if (!is.null(lv) && length(tv) == 1L) tab <- tab[names(tab) %in% lv]
    as.integer(tab)
  }))
  if (length(sizes)) min(sizes) else 2L
}

build_dds <- function(counts, lengths, samples, design) {
  cd <- samples
  if (!is.null(lengths)) {
    txi <- list(counts = counts, length = lengths, abundance = NULL, countsFromAbundance = "no")
    dds <- suppressMessages(DESeq2::DESeqDataSetFromTximport(txi, cd, design))
  } else {
    m <- round(counts)
    storage.mode(m) <- "integer"
    dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(m, cd, design))
  }
  dds
}

# apeglm needs a coefficient: relevel the factor so the denominator is the
# reference and re-run nbinomWaldTest (DESeq2 vignette), then shrink.
shrink_apeglm <- function(dds, ct) {
  n_warn <- 0L
  muffle <- function(w) {
    if (grepl("line search|optim|converge", conditionMessage(w))) {
      n_warn <<- n_warn + 1L
      invokeRestart("muffleWarning")
    }
  }
  if (ct$type == "triple") {
    f <- ct$factor
    d2 <- dds
    if (levels(d2[[f]])[1] != ct$den) {
      d2[[f]] <- stats::relevel(d2[[f]], ref = ct$den)
      d2 <- suppressMessages(DESeq2::nbinomWaldTest(d2, quiet = TRUE))
    }
    coef <- paste0(f, "_", make.names(ct$num), "_vs_", make.names(ct$den))
  } else {
    d2 <- dds
    coef <- ct$name
  }
  rn <- DESeq2::resultsNames(d2)
  if (!coef %in% rn) coef <- rn[make.names(rn) == make.names(coef)][1]
  if (is.na(coef)) return(NULL)
  s <- withCallingHandlers(
    suppressMessages(DESeq2::lfcShrink(d2, coef = coef, type = "apeglm", quiet = TRUE)),
    warning = muffle)
  list(lfc = s$log2FoldChange, se = s$lfcSE, n_warn = n_warn)
}

gene_status <- function(res) {
  st <- rep("tested", nrow(res))
  st[is.na(res$padj) & !is.na(res$pvalue)] <- "filtered_low_mean (independent filtering)"
  st[is.na(res$pvalue) & !is.na(res$baseMean) & res$baseMean > 0] <- "count_outlier (Cook's distance)"
  st[!is.na(res$baseMean) & res$baseMean == 0] <- "all_zero"
  st
}

# Core analysis. Returns list(dds, tables (named list of data.frames), summary, warn, info).
run_de <- function(counts, samples, design, lengths = NULL, contrasts = NULL, test = "Wald",
                   reduced = NULL, alpha = 0.05, lfc_threshold = 0, shrink = "apeglm",
                   min_count = 10, min_samples = 0, min_reps = 3, reference = list()) {
  if (!test %in% c("Wald", "LRT")) fail("--test must be Wald or LRT.")
  if (test == "LRT" && is.null(reduced)) fail("--test LRT needs --reduced, e.g. \"~ batch\".")
  if (test == "LRT" && length(contrasts))
    fail("--contrast is for Wald tests; the LRT tests the terms dropped in --reduced.")
  if (alpha <= 0 || alpha >= 1) fail("--alpha must be between 0 and 1.")
  if (test == "Wald" && !length(contrasts)) contrasts <- default_contrasts(samples, design, reference)
  for (ct in contrasts) if (ct$type == "triple") {
    if (!ct$factor %in% all.vars(design)) fail("--contrast factor '", ct$factor, "' is not in the design.")
    lv <- levels(samples[[ct$factor]])
    bad <- setdiff(c(ct$num, ct$den), lv)
    if (length(bad)) fail("--contrast level(s) ", paste(bad, collapse = ", "), " not in '", ct$factor,
                          "'. Levels: ", paste(lv, collapse = ", "))
  }
  if (test == "LRT") {
    extra <- setdiff(all.vars(reduced), all.vars(design))
    if (length(extra)) fail("--reduced uses variables not in --design: ", paste(extra, collapse = ", "))
  }
  tested <- tested_vars(contrasts, design, reduced, test, samples)
  cd <- check_design(samples, design, tested, min_reps)
  warn <- cd$warn
  if (test == "LRT" && smallest_group(samples, tested) <= 3)
    warn <- c(warn, paste0("With 3 or fewer replicates per group, DESeq2's likelihood-ratio test was anti-conservative ",
                           "in AgriClaw simulations (empirical FDR about twice the nominal level; nominal 0.025 gave ",
                           "about 0.06). Consider --alpha 0.025, and treat genes near the cut-off with caution."))
  if (min_samples <= 0) min_samples <- smallest_group(samples, tested)
  keep <- rowSums(counts >= min_count) >= min_samples
  if (sum(keep) < 50) fail("Only ", sum(keep), " genes pass the pre-filter (>= ", min_count, " counts in >= ",
                           min_samples, " samples). Check that the input holds counts.")
  dds <- build_dds(counts[keep, , drop = FALSE], if (is.null(lengths)) NULL else lengths[keep, , drop = FALSE],
                   samples, design)
  dds <- if (test == "LRT") suppressMessages(DESeq2::DESeq(dds, test = "LRT", reduced = reduced, quiet = TRUE))
         else suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
  tables <- list()
  summ <- list()
  items <- if (test == "LRT") list(list(type = "lrt", label = "LRT")) else contrasts
  for (ct in items) {
    res <- if (ct$type == "lrt") DESeq2::results(dds, alpha = alpha)
           else if (ct$type == "triple") DESeq2::results(dds, contrast = c(ct$factor, ct$num, ct$den),
                                                         alpha = alpha, lfcThreshold = lfc_threshold)
           else {
             rn <- DESeq2::resultsNames(dds)
             if (!ct$name %in% rn) fail("--contrast coefficient '", ct$name, "' not found. Available: ",
                                        paste(rn, collapse = ", "))
             DESeq2::results(dds, name = ct$name, alpha = alpha, lfcThreshold = lfc_threshold)
           }
    df <- data.frame(gene_id = rownames(res), baseMean = res$baseMean,
                     log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE, stat = res$stat,
                     pvalue = res$pvalue, padj = res$padj, stringsAsFactors = FALSE)
    shr <- NULL
    if (shrink == "apeglm" && ct$type != "lrt") shr <- shrink_apeglm(dds, ct)
    df$log2FoldChange_shrunk <- if (is.null(shr)) NA_real_ else shr$lfc
    df$lfcSE_shrunk <- if (is.null(shr)) NA_real_ else shr$se
    df$status <- gene_status(res)
    df$significant <- !is.na(df$padj) & df$padj < alpha
    df$direction <- ifelse(!df$significant, "ns", ifelse(df$log2FoldChange > 0, "up", "down"))
    pre <- rownames(counts)[!keep]
    if (length(pre)) {
      add <- data.frame(gene_id = pre, baseMean = NA_real_, log2FoldChange = NA_real_, lfcSE = NA_real_,
                        stat = NA_real_, pvalue = NA_real_, padj = NA_real_, log2FoldChange_shrunk = NA_real_,
                        lfcSE_shrunk = NA_real_, status = "prefiltered (low counts)", significant = FALSE,
                        direction = "ns", stringsAsFactors = FALSE)
      df <- rbind(df, add)
    }
    df <- df[order(df$padj, df$pvalue, na.last = TRUE), , drop = FALSE]
    md <- S4Vectors::metadata(res)
    lfc_name <- if (ct$type == "lrt") utils::tail(DESeq2::resultsNames(dds), 1) else ct$label
    tables[[ct$label]] <- df
    summ[[ct$label]] <- data.frame(
      comparison = ct$label, test = test, n_genes_input = nrow(counts), n_prefiltered = sum(!keep),
      n_tested = sum(!is.na(df$padj)), n_low_mean_filtered = sum(grepl("^filtered", df$status)),
      n_count_outliers = sum(grepl("^count_outlier", df$status)), filter_threshold_baseMean =
        if (is.null(md$filterThreshold)) NA_real_ else unname(md$filterThreshold),
      n_up = sum(df$direction == "up"), n_down = sum(df$direction == "down"), alpha = alpha,
      lfc_threshold = lfc_threshold, lfc_column = lfc_name,
      apeglm_optimizer_warnings = if (is.null(shr)) NA_integer_ else shr$n_warn, stringsAsFactors = FALSE)
  }
  list(dds = dds, tables = tables, summary = do.call(rbind, summ), warn = warn,
       info = list(min_samples = min_samples, n_kept = sum(keep), df_resid = cd$df_resid,
                   tested = tested, contrasts = contrasts))
}

# --- Figures ------------------------------------------------------------------

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

plot_ma <- function(df, label, out, alpha) {
  d <- df[!is.na(df$baseMean) & df$baseMean > 0, ]
  y <- if (all(is.na(d$log2FoldChange_shrunk))) d$log2FoldChange else d$log2FoldChange_shrunk
  ylab <- if (all(is.na(d$log2FoldChange_shrunk))) "log2 fold change (MLE)" else "log2 fold change (apeglm-shrunken)"
  lim <- max(2, stats::quantile(abs(y), 0.999, na.rm = TRUE))
  yc <- pmax(pmin(y, lim), -lim)
  open_png(out, paste0("ma_", safe_name(label), ".png"))
  op <- graphics::par(mar = c(4.5, 4.5, 3, 1))
  on.exit({graphics::par(op); grDevices::dev.off()})
  graphics::plot(d$baseMean, yc, log = "x", pch = ifelse(abs(y) > lim, 2, 20), cex = 0.5,
                 col = ifelse(d$significant, "#C0392B", "#95A5A6"), ylim = c(-lim, lim),
                 xlab = "Mean of normalised counts", ylab = ylab,
                 main = sprintf("%s: %d up, %d down (padj < %g)", label, sum(d$direction == "up"),
                                sum(d$direction == "down"), alpha))
  graphics::abline(h = 0, col = "grey40")
}

plot_volcano <- function(df, label, out, alpha) {
  d <- df[!is.na(df$pvalue), ]
  x <- if (all(is.na(d$log2FoldChange_shrunk))) d$log2FoldChange else d$log2FoldChange_shrunk
  y <- -log10(pmax(d$pvalue, .Machine$double.xmin))
  open_png(out, paste0("volcano_", safe_name(label), ".png"))
  op <- graphics::par(mar = c(4.5, 4.5, 3, 1))
  on.exit({graphics::par(op); grDevices::dev.off()})
  graphics::plot(x, y, pch = 20, cex = 0.5, col = ifelse(d$significant, "#C0392B", "#95A5A6"),
                 xlab = if (all(is.na(d$log2FoldChange_shrunk))) "log2 fold change (MLE)" else
                   "log2 fold change (apeglm-shrunken)", ylab = "-log10 p-value", main = label)
  if (any(d$significant)) graphics::abline(h = -log10(max(d$pvalue[d$significant])), lty = 2)
  graphics::legend("topleft", c(sprintf("padj < %g", alpha), "not significant"), pch = 20,
                   col = c("#C0392B", "#95A5A6"), bty = "n", cex = 0.8)
}

plot_pvalue_hist <- function(df, label, out) {
  p <- df$pvalue[!is.na(df$padj)]
  open_png(out, paste0("pvalue_hist_", safe_name(label), ".png"), width = 6, height = 4.5)
  op <- graphics::par(mar = c(4.5, 4.5, 3, 1))
  on.exit({graphics::par(op); grDevices::dev.off()})
  graphics::hist(p, breaks = seq(0, 1, 0.05), col = "#BDC3C7", border = "white",
                 main = paste("p-values of tested genes:", label), xlab = "p-value")
}

transform_dds <- function(dds) {
  n_ok <- sum(rowMeans(DESeq2::counts(dds, normalized = TRUE)) > 5)
  tr <- if (n_ok >= 1000) DESeq2::vst(dds, blind = FALSE) else
    DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
  SummarizedExperiment::assay(tr)
}

plot_pca <- function(vmat, samples, colour_var, shape_var, out) {
  rv <- apply(vmat, 1, stats::var)
  sel <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
  p <- stats::prcomp(t(vmat[sel, , drop = FALSE]))
  ve <- 100 * p$sdev^2 / sum(p$sdev^2)
  g <- if (length(colour_var)) factor(samples[colnames(vmat), colour_var]) else factor(rep("all", ncol(vmat)))
  sh <- if (length(shape_var)) factor(samples[colnames(vmat), shape_var]) else factor(rep("all", ncol(vmat)))
  pal <- grDevices::hcl.colors(max(3, nlevels(g)), "Dark 3")[seq_len(nlevels(g))]
  pchs <- c(19, 17, 15, 18, 1, 2, 0, 5)[(seq_len(nlevels(sh)) - 1) %% 8 + 1]
  open_png(out, "pca.png", width = 7, height = 5.5)
  op <- graphics::par(mar = c(4.5, 4.5, 3, 9), xpd = NA)
  on.exit({graphics::par(op); grDevices::dev.off()})
  graphics::plot(p$x[, 1], p$x[, 2], col = pal[as.integer(g)], pch = pchs[as.integer(sh)], cex = 1.4,
                 xlab = sprintf("PC1 (%.1f%%)", ve[1]), ylab = sprintf("PC2 (%.1f%%)", ve[2]),
                 main = sprintf("PCA, VST (blind = FALSE), %d genes", length(sel)))
  graphics::text(p$x[, 1], p$x[, 2], colnames(vmat), pos = 3, cex = 0.65)
  usr <- graphics::par("usr")
  leg_x <- usr[2] + 0.03 * (usr[2] - usr[1])
  graphics::legend(leg_x, usr[4], levels(g), col = pal, pch = 19, bty = "n", cex = 0.8,
                   title = if (length(colour_var)) colour_var else NULL)
  if (length(shape_var))
    graphics::legend(leg_x, usr[4] - 0.45 * (usr[4] - usr[3]), levels(sh), pch = pchs, bty = "n",
                     cex = 0.8, title = shape_var)
}

# --- Main -------------------------------------------------------------------

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (opt$demo) {
    opt$counts <- file.path(SKILL_DIR, "demo", "counts.csv")
    opt$samples <- file.path(SKILL_DIR, "demo", "samples.csv")
    if (!nzchar(opt$design)) opt$design <- "~ batch + condition"
    if (!nzchar(opt$reference)) opt$reference <- "condition=control"
  }
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  if (!opt$shrink %in% c("apeglm", "none")) fail("--shrink must be apeglm or none.")
  inp <- read_counts_input(opt$txi, opt$counts, opt$lengths, opt$no_length_offset)
  warn <- check_count_matrix(inp$counts, inp$noninteger_ok, opt$allow_noninteger)
  notes <- character()
  if (opt$allow_noninteger) notes <- c(notes, "Non-integer counts were rounded (--allow-noninteger).")
  if (opt$no_length_offset) notes <- c(notes, "Transcript-length offsets were not used (--no-length-offset).")
  if (!identical(inp$cfa, "no"))
    notes <- c(notes, sprintf("tximport counts were made with countsFromAbundance = '%s': no length offsets.", inp$cfa))
  samples <- read_samples(opt$samples, colnames(inp$counts), split_list(opt$exclude))
  counts <- inp$counts[, samples$sample, drop = FALSE]
  lengths <- if (is.null(inp$lengths)) NULL else inp$lengths[, samples$sample, drop = FALSE]
  design <- parse_formula(opt$design)
  reduced <- if (nzchar(opt$reduced)) parse_formula(opt$reduced, "--reduced") else NULL
  reference <- parse_reference(opt$reference)
  pd <- prepare_design(samples, design, reference, split_list(opt$as_factor))
  samples <- pd$samples
  warn <- c(warn, pd$warn)
  contrasts <- parse_contrasts(opt$contrast)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

  res <- run_de(counts, samples, design, lengths, contrasts, opt$test, reduced, opt$alpha,
                opt$lfc_threshold, opt$shrink, opt$min_count, opt$min_samples, opt$min_reps, reference)
  warn <- c(warn, res$warn)
  dds <- res$dds
  for (nm in names(res$tables)) {
    df <- res$tables[[nm]]
    save_table(df, opt$out, paste0("de_", safe_name(nm), ".csv"))
    plot_ma(df, nm, opt$out, opt$alpha)
    plot_volcano(df, nm, opt$out, opt$alpha)
    plot_pvalue_hist(df, nm, opt$out)
  }
  save_table(res$summary, opt$out, "de_summary.csv")
  ncnt <- DESeq2::counts(dds, normalized = TRUE)
  save_table(data.frame(gene_id = rownames(ncnt), round(ncnt, 3), check.names = FALSE), opt$out,
             "normalized_counts.csv")
  vmat <- transform_dds(dds)
  save_table(data.frame(gene_id = rownames(vmat), round(vmat, 5), check.names = FALSE), opt$out, "vst.csv")
  nf <- if (is.null(DESeq2::normalizationFactors(dds))) DESeq2::sizeFactors(dds) else
    exp(colMeans(log(DESeq2::normalizationFactors(dds))))
  save_table(data.frame(sample = colnames(dds), size_factor = round(nf, 5)), opt$out, "size_factors.csv")
  save_table(samples, opt$out, "samples_used.csv")
  fac_tested <- unique(unlist(res$info$tested))
  fac_tested <- fac_tested[vapply(fac_tested, function(v) is.factor(samples[[v]]), TRUE)]
  other <- setdiff(all.vars(design), fac_tested)
  other <- other[vapply(other, function(v) is.factor(samples[[v]]), TRUE)]
  plot_pca(vmat, samples, utils::head(fac_tested, 1), utils::head(other, 1), opt$out)
  open_png(opt$out, "dispersion.png", width = 6.5, height = 5)
  DESeq2::plotDispEsts(dds, main = "Dispersion estimates")
  grDevices::dev.off()

  s <- res$summary
  big_out <- s$n_count_outliers > 0.01 * s$n_tested
  if (any(big_out)) warn <- c(warn, sprintf("%s: %d genes (> 1%% of tested) have count outliers by Cook's distance; check for an outlier sample.",
                                            s$comparison[big_out], s$n_count_outliers[big_out]))
  none <- s$n_up + s$n_down == 0
  if (any(none)) warn <- c(warn, sprintf("%s: no genes at padj < %g.", s$comparison[none], opt$alpha))

  grp_lines <- character()
  for (v in all.vars(design)) if (is.factor(samples[[v]])) {
    tb <- table(samples[[v]])
    grp_lines <- c(grp_lines, sprintf("- %s: %s", v, paste(sprintf("%s (n = %d)", names(tb), as.integer(tb)), collapse = ", ")))
  }
  rep <- c(sprintf("# %s report", SKILL), "",
           if (length(warn)) c("## Warnings", "", paste("-", warn), "") else NULL,
           "## Design", "",
           sprintf("- Input: %s (%d genes, %d samples)%s.", inp$source, nrow(counts), ncol(counts),
                   if (nzchar(opt$exclude)) paste0("; excluded: ", opt$exclude) else ""),
           sprintf("- Design: `%s`; test: %s%s.", opt$design, opt$test,
                   if (opt$test == "LRT") sprintf(" (reduced `%s`)", opt$reduced) else ""),
           sprintf("- Reference levels: %s.", paste(pd$references, collapse = "; ")),
           grp_lines,
           sprintf("- Pre-filter: kept %d of %d genes with >= %g counts in >= %d samples; residual df = %d.",
                   res$info$n_kept, nrow(counts), opt$min_count, res$info$min_samples, res$info$df_resid),
           if (length(notes)) paste("-", notes) else NULL, "",
           "## Results", "",
           sprintf("FDR (Benjamini-Hochberg) < %g; independent filtering optimised for alpha = %g%s.",
                   opt$alpha, opt$alpha, if (opt$lfc_threshold > 0)
                     sprintf("; tests H0: |log2FC| <= %g", opt$lfc_threshold) else ""), "",
           "| comparison | tested | up | down | low-mean filtered | count outliers | filter threshold (baseMean) |",
           "|---|---|---|---|---|---|---|",
           sprintf("| %s | %d | %d | %d | %d | %d | %.2f |", s$comparison, s$n_tested, s$n_up, s$n_down,
                   s$n_low_mean_filtered, s$n_count_outliers, s$filter_threshold_baseMean), "")
  if (opt$test == "LRT")
    rep <- c(rep, sprintf(paste0("The LRT p-value tests all terms dropped from the full model at once. ",
                                 "The `log2FoldChange` column is the estimate for one coefficient only (%s); ",
                                 "do not read it as the effect being tested."), s$lfc_column[1]), "")
  for (nm in names(res$tables)) {
    top <- utils::head(res$tables[[nm]][res$tables[[nm]]$significant, ], 10)
    rep <- c(rep, sprintf("### %s", nm), "")
    if (nrow(top)) {
      rep <- c(rep, "| gene | baseMean | log2FC (MLE) | log2FC (shrunken) | padj |", "|---|---|---|---|---|",
               sprintf("| %s | %.1f | %.2f | %s | %.2e |", top$gene_id, top$baseMean, top$log2FoldChange,
                       ifelse(is.na(top$log2FoldChange_shrunk), "NA", sprintf("%.2f", top$log2FoldChange_shrunk)),
                       top$padj), "")
    } else rep <- c(rep, "No significant genes.", "")
    rep <- c(rep, sprintf("![MA plot](figures/ma_%s.png) ![Volcano](figures/volcano_%s.png) ![p-values](figures/pvalue_hist_%s.png)",
                          safe_name(nm), safe_name(nm), safe_name(nm)), "")
  }
  rep <- c(rep, "![PCA](figures/pca.png) ![Dispersion](figures/dispersion.png)", "",
           "## How to read the tables", "",
           "- `log2FoldChange`/`pvalue`/`padj`: Wald (or LRT) test on the maximum-likelihood estimate; `significant` = padj < alpha.",
           "- `log2FoldChange_shrunk`: apeglm-shrunken estimate for ranking and plots. Its p-values were not used.",
           "- `status`: tested; filtered_low_mean (independent filtering, padj = NA); count_outlier (Cook's distance, p = NA); prefiltered (low counts, not tested).",
           "- Genes with non-missing `padj` are the tested universe to use as background in claw-enrichment.",
           "- `tables/vst.csv` (VST, blind = FALSE) is for clustering, PCA and claw-coexpression, never for DE tests.", "",
           "## Methods", "",
           paste0(sprintf("Differential expression was tested with DESeq2 %s (%s test, design `%s`). ",
                          as.character(utils::packageVersion("DESeq2")), opt$test, opt$design),
                  if (!is.null(lengths)) "Gene counts and average transcript lengths from tximport were used as offsets (DESeqDataSetFromTximport). " else "",
                  sprintf("Genes with fewer than %g counts in %d samples were removed before fitting. ", opt$min_count, res$info$min_samples),
                  sprintf("P-values were adjusted with Benjamini-Hochberg; independent filtering on the mean count used alpha = %g. ", opt$alpha),
                  if (opt$shrink == "apeglm" && opt$test == "Wald") sprintf("Log2 fold changes for ranking and plots were shrunk with apeglm %s.", as.character(utils::packageVersion("apeglm"))) else ""),
           "", "## References", "",
           "- Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15:550. https://doi.org/10.1186/s13059-014-0550-8",
           "- Zhu A, Ibrahim JG, Love MI (2019). Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics 35:2084-2092. https://doi.org/10.1093/bioinformatics/bty895",
           "- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society B 57:289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x",
           "- Bourgon R, Gentleman R, Huber W (2010). Independent filtering increases detection power for high-throughput experiments. PNAS 107:9546-9551. https://doi.org/10.1073/pnas.0914005107",
           "- Schurch NJ, Schofield P, Gierlinski M, et al. (2016). How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? RNA 22:839-851. https://doi.org/10.1261/rna.053959.115",
           "", sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(rep, file.path(opt$out, "report.md"))
  write_repro_bundle(opt$out, argv, inputs = c(opt$txi, opt$counts, opt$lengths, opt$samples))
  cat(sprintf("Done: %s. Report: %s\n", paste(sprintf("%s %d up / %d down", s$comparison, s$n_up, s$n_down),
                                               collapse = "; "),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(res)
}

if (sys.nframe() == 0L) main()
