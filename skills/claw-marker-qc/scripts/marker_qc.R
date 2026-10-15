#!/usr/bin/env Rscript
# claw-marker-qc: quality control and format conversion for crop genotype data
# (VCF, HapMap or a 0/1/2 matrix). Per marker: call rate, minor allele
# frequency, heterozygosity, monomorphic markers. Per sample: call rate,
# heterozygosity (contamination or outcrossing in inbred lines) and
# near-identical duplicates by identity-by-state. Optional checks for
# biparental populations: segregation distortion (RIL/DH 1:1, F2 1:2:1,
# BC 1:1) and consistency with the expected parents. Writes a filtered 0/1/2
# matrix of ALT-allele counts plus QC tables and plots.
#
# Usage:
#   Rscript scripts/marker_qc.R --demo
#   Rscript scripts/marker_qc.R --input panel.vcf.gz --germplasm inbred
#   Rscript scripts/marker_qc.R --input rils.hmp.txt --germplasm inbred \
#           --population ril --parents P1,P2
#   Rscript scripts/marker_qc.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT license.

SKILL   <- "claw-marker-qc"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/marker_qc.R [--demo] --input FILE --germplasm inbred|outbred
                                 [options] [--out DIR]

Input
  --demo                   run on demo/demo_panel.vcf.gz (inbred panel with planted problems)
  --input FILE             VCF (.vcf/.vcf.gz), HapMap (.hmp.txt) or 0/1/2 matrix (.csv/.tsv, .gz ok)
  --format F               auto (default), vcf, hapmap or matrix
  --orientation O          matrix only: auto (default), markers-in-rows or samples-in-rows
  --map FILE               matrix only: CSV with marker, chr, pos (used to place markers
                           and to detect the orientation)
  --germplasm G            REQUIRED: inbred (lines, RILs, DH, pure-line cultivars) or
                           outbred (open-pollinated, hybrids, clones, F2 individuals)

Thresholds (adjustable rules of thumb, see SKILL.md)
  --min-marker-call-rate X drop markers called in fewer than X of samples (default 0.9)
  --min-sample-call-rate X drop samples called at fewer than X of markers (default 0.8)
  --min-maf X              drop markers with minor allele frequency below X (default 0.05)
  --max-het X              inbred: flag samples with heterozygosity above X (default 0.10)
  --het-sd X               outbred: flag samples outside mean +/- X SD heterozygosity (default 3)
  --ibs-dup X              flag sample pairs with identity-by-state >= X as duplicates (default 0.98)
  --ibs-markers N          markers used for identity-by-state (random subset; default 5000)
  --min-shared N           minimum markers called in both samples to judge a pair (default 100)
  --drop-duplicates        keep only the best-called sample of each duplicate group
  --drop-het-outliers      remove samples flagged for heterozygosity

Biparental populations (optional)
  --population P           ril, dh, f2 or bc1 (backcross to the first parent)
  --parents P1,P2          sample IDs of the two parents (required for bc1)
  --seg-alpha X            family-wise alpha for segregation distortion, Bonferroni (default 0.05)
  --drop-distorted         remove markers with significant segregation distortion
  --max-nonparental X      flag progeny with more than X non-parental calls (default 0.05)
  --drop-nonparental       remove progeny flagged as inconsistent with the parents

Output
  --write-vcf              also write the filtered genotypes as VCF (needs REF/ALT alleles)
  --seed N                 random seed for the identity-by-state marker subset (default 2026)
  --out DIR                output folder (default: agriclaw_output/claw-marker-qc)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, input = "", format = "auto", orientation = "auto",
                 map = "", germplasm = "", min_marker_call_rate = 0.9,
                 min_sample_call_rate = 0.8, min_maf = 0.05, max_het = 0.10,
                 het_sd = 3, ibs_dup = 0.98, ibs_markers = 5000, min_shared = 100,
                 drop_duplicates = FALSE, drop_het_outliers = FALSE,
                 population = "", parents = "", seg_alpha = 0.05,
                 drop_distorted = FALSE, max_nonparental = 0.05,
                 drop_nonparental = FALSE, write_vcf = FALSE, seed = 2026, out = "")

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
# Internal representation: `geno` is an integer matrix with markers in rows and
# samples in columns, coded as the number of copies of the ALT (counted)
# allele: 0, 1, 2 or NA. `map` has one row per marker: marker, chr, pos, ref, alt.

COL_BLUE <- "#2a78d6"; COL_ORANGE <- "#eb6834"; COL_INK <- "#0b0b0b"
COL_INK2 <- "#52514e"; COL_MUTED <- "#898781"; COL_FILL <- "#86b6ef"

open_text <- function(path) {
  if (!nzchar(path) || !file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
  gzfile(path, open = "rt")  # reads plain text and .gz alike
}

first_line <- function(path) {
  con <- open_text(path)
  on.exit(close(con))
  l <- readLines(con, n = 1L, warn = FALSE)
  if (!length(l)) stop("Input file is empty: ", path, call. = FALSE)
  l
}

detect_format <- function(path) {
  l <- first_line(path)
  if (startsWith(l, "##fileformat=VCF")) return("vcf")
  h <- tolower(strsplit(l, "\t", fixed = TRUE)[[1]])
  if (length(h) >= 12L && h[1] %in% c("rs#", "rs") && h[2] == "alleles") return("hapmap")
  "matrix"
}

make_map <- function(marker, chr = NA, pos = NA, ref = NA, alt = NA) {
  data.frame(marker = as.character(marker), chr = as.character(chr),
             pos = suppressWarnings(as.numeric(pos)), ref = as.character(ref),
             alt = as.character(alt), stringsAsFactors = FALSE)
}

unique_ids <- function(ids, what, notes) {
  dup <- duplicated(ids)
  if (any(dup)) {
    notes <- c(notes, sprintf("%d duplicated %s ID(s) were made unique with suffixes (e.g. %s).",
                              sum(dup), what, ids[dup][1]))
    ids <- make.unique(ids, sep = "_dup")
  }
  list(ids = ids, notes = notes)
}

# --- VCF --------------------------------------------------------------------

# Convert a character matrix of VCF GT strings to ALT-allele counts.
# Biallelic diploid calls only: 0/0, 0/1, 1/0, 1/1 (phased '|' accepted).
gt_to_dosage <- function(gt) {
  g <- gsub("|", "/", gt, fixed = TRUE)
  g[g %in% c(".", "./.")] <- NA
  called <- !is.na(g)
  n_sep <- nchar(g) - nchar(gsub("/", "", g, fixed = TRUE))
  bad_ploidy <- called & n_sep != 1L
  if (any(bad_ploidy)) {
    stop("Only diploid genotype calls (GT like 0/1) are supported; found ",
         paste(utils::head(unique(g[bad_ploidy]), 3), collapse = ", "),
         ". Convert haploid or polyploid calls before running claw-marker-qc.", call. = FALSE)
  }
  a1 <- sub("/.*$", "", g)
  a2 <- sub("^.*/", "", g)
  half <- called & (a1 == "." | a2 == ".")
  ok <- called & !half & a1 %in% c("0", "1") & a2 %in% c("0", "1")
  bad <- called & !half & !ok
  if (any(bad)) {
    stop("Unexpected GT values at biallelic sites: ",
         paste(utils::head(unique(g[bad]), 3), collapse = ", "),
         ". Only 0/0, 0/1, 1/1 and missing calls are allowed.", call. = FALSE)
  }
  d <- matrix(NA_integer_, nrow(g), ncol(g), dimnames = dimnames(gt))
  d[ok] <- (a1[ok] == "1") + (a2[ok] == "1")
  attr(d, "n_half_missing") <- sum(half)
  d
}

read_vcf_geno <- function(path) {
  vcf <- vcfR::read.vcfR(path, verbose = FALSE)
  fix <- vcfR::getFIX(vcf)
  if (!nrow(fix)) stop("The VCF has no variant records: ", path, call. = FALSE)
  if (ncol(vcf@gt) < 2L) stop("The VCF has no sample (genotype) columns: ", path, call. = FALSE)
  notes <- character()
  gt <- vcfR::extract.gt(vcf, element = "GT", IDtoRowNames = FALSE)
  alt <- fix[, "ALT"]
  multi <- !is.na(alt) & grepl(",", alt, fixed = TRUE)
  if (any(multi)) {
    notes <- c(notes, sprintf("%d multi-allelic site(s) were excluded (ALT-allele counts need biallelic markers).",
                              sum(multi)))
  }
  keep <- !multi
  if (!any(keep)) stop("No biallelic sites in the VCF.", call. = FALSE)
  fix <- fix[keep, , drop = FALSE]
  gt <- gt[keep, , drop = FALSE]
  ids <- fix[, "ID"]
  noid <- is.na(ids) | ids == "."
  ids[noid] <- paste0(fix[noid, "CHROM"], "_", fix[noid, "POS"])
  u <- unique_ids(ids, "marker", notes); ids <- u$ids; notes <- u$notes
  phased <- sum(grepl("|", gt, fixed = TRUE), na.rm = TRUE)
  if (phased > 0) notes <- c(notes, "Phased calls (0|1) were read as unphased genotypes; phase is not kept.")
  d <- gt_to_dosage(gt)
  nh <- attr(d, "n_half_missing")
  if (nh > 0) notes <- c(notes, sprintf("%d half-missing call(s) such as ./1 were set to missing.", nh))
  attr(d, "n_half_missing") <- NULL
  samples <- colnames(vcf@gt)[-1]
  u <- unique_ids(samples, "sample", notes); samples <- u$ids; notes <- u$notes
  dimnames(d) <- list(ids, samples)
  map <- make_map(ids, fix[, "CHROM"], fix[, "POS"], fix[, "REF"],
                  ifelse(is.na(alt[keep]), ".", alt[keep]))
  list(geno = d, map = map, notes = notes, format = "vcf", orientation = "VCF (markers in rows)")
}

# --- HapMap -------------------------------------------------------------------

IUPAC <- c(A = "AA", C = "CC", G = "GG", T = "TT", R = "AG", Y = "CT", S = "CG",
           W = "AT", K = "GT", M = "AC", "+" = "++", "-" = "--")

# HapMap calls (two-letter 'AG' or single-letter IUPAC 'R') to copies of the
# second allele listed in the 'alleles' column.
hapmap_dosage <- function(calls, a1, a2) {
  calls <- toupper(trimws(calls))
  miss <- is.na(calls) | calls %in% c("", "N", "NN", "NA", ".", "./.")
  two <- calls
  one <- !miss & nchar(calls) == 1L
  two[one] <- unname(IUPAC[calls[one]])
  two <- gsub("[/|:]", "", two)
  bad <- !miss & (is.na(two) | nchar(two) != 2L)
  x1 <- substr(two, 1L, 1L)
  x2 <- substr(two, 2L, 2L)
  A1 <- matrix(a1, nrow(calls), ncol(calls))
  A2 <- matrix(a2, nrow(calls), ncol(calls))
  inset <- (x1 == A1 | x1 == A2) & (x2 == A1 | x2 == A2)
  bad <- bad | (!miss & !bad & !inset)
  d <- matrix(NA_integer_, nrow(calls), ncol(calls), dimnames = dimnames(calls))
  okc <- !miss & !bad
  d[okc] <- (x1[okc] == A2[okc]) + (x2[okc] == A2[okc])
  attr(d, "n_bad") <- sum(bad)
  d
}

read_hapmap_geno <- function(path) {
  con <- open_text(path)
  df <- tryCatch(utils::read.delim(con, colClasses = "character", check.names = FALSE,
                                   comment.char = "", quote = "", na.strings = character()),
                 finally = close(con))
  if (ncol(df) < 12L) stop("A HapMap file needs 11 annotation columns plus sample columns.", call. = FALSE)
  notes <- character()
  al <- strsplit(df[[2]], "/", fixed = TRUE)
  n_al <- lengths(al)
  keep <- n_al == 2L
  if (any(!keep)) notes <- c(notes, sprintf("%d marker(s) without exactly two alleles in the 'alleles' column were excluded.",
                                            sum(!keep)))
  if (!any(keep)) stop("No biallelic markers in the HapMap file.", call. = FALSE)
  df <- df[keep, , drop = FALSE]
  al <- al[keep]
  a1 <- toupper(vapply(al, `[`, "", 1L))
  a2 <- toupper(vapply(al, `[`, "", 2L))
  calls <- as.matrix(df[, -(1:11), drop = FALSE])
  u <- unique_ids(df[[1]], "marker", notes); ids <- u$ids; notes <- u$notes
  u <- unique_ids(colnames(calls), "sample", notes); samples <- u$ids; notes <- u$notes
  dimnames(calls) <- list(ids, samples)
  d <- hapmap_dosage(calls, a1, a2)
  nb <- attr(d, "n_bad")
  attr(d, "n_bad") <- NULL
  if (nb > 0) notes <- c(notes, sprintf("%d call(s) did not match the marker's two alleles and were set to missing.", nb))
  notes <- c(notes, "HapMap coding: 0/1/2 counts copies of the SECOND allele in the 'alleles' column (A/G counts G), reported as ALT.")
  map <- make_map(ids, df[[3]], df[[4]], a1, a2)
  list(geno = d, map = map, notes = notes, format = "hapmap", orientation = "HapMap (markers in rows)")
}

# --- 0/1/2 matrix ---------------------------------------------------------------

ALIASES <- list(marker = c("marker", "snp", "rs#", "marker_id", "snp_id", "id", "name"),
                chr = c("chr", "chrom", "chromosome"),
                pos = c("pos", "position", "bp", "pos_bp", "position_bp"),
                ref = c("ref", "ref_allele"), alt = c("alt", "alt_allele"))

find_col <- function(nm, aliases) {
  i <- match(TRUE, tolower(nm) %in% aliases)
  if (is.na(i)) NA_integer_ else i
}

MARKER_ID_RE <- "(^ss[0-9]+$)|(^rs[0-9]+$)|((^|[_.-])(chr|gm|s)?[0-9]{1,2}[a-z]?[_:.-][0-9]{3,}($|[_.-]))"

marker_like <- function(ids) mean(grepl(MARKER_ID_RE, ids, ignore.case = TRUE, perl = TRUE))

# Decide whether the first column of a matrix holds marker or sample IDs.
decide_orientation <- function(row_ids, col_ids, orientation, map = NULL) {
  if (orientation %in% c("markers-in-rows", "samples-in-rows"))
    return(list(orientation = orientation, basis = "given with --orientation"))
  if (orientation != "auto")
    stop("--orientation must be auto, markers-in-rows or samples-in-rows.", call. = FALSE)
  if (!is.null(map)) {
    in_rows <- mean(row_ids %in% map$marker)
    in_cols <- mean(col_ids %in% map$marker)
    if (in_rows >= 0.5 && in_rows > in_cols)
      return(list(orientation = "markers-in-rows", basis = sprintf("%.0f%% of row IDs are in the --map file", 100 * in_rows)))
    if (in_cols >= 0.5 && in_cols > in_rows)
      return(list(orientation = "samples-in-rows", basis = sprintf("%.0f%% of column IDs are in the --map file", 100 * in_cols)))
  }
  fr <- marker_like(row_ids)
  fc <- marker_like(col_ids)
  if (fr >= 0.9 && fc <= 0.1)
    return(list(orientation = "markers-in-rows", basis = sprintf("%.0f%% of row IDs look like marker names (chromosome_position)", 100 * fr)))
  if (fc >= 0.9 && fr <= 0.1)
    return(list(orientation = "samples-in-rows", basis = sprintf("%.0f%% of column IDs look like marker names (chromosome_position)", 100 * fc)))
  stop("Cannot tell whether markers are in rows or in columns of this matrix (first row IDs: ",
       paste(utils::head(row_ids, 3), collapse = ", "), "; first column IDs: ",
       paste(utils::head(col_ids, 3), collapse = ", "), "). Ask the user, then re-run with ",
       "--orientation markers-in-rows or --orientation samples-in-rows, or give --map FILE ",
       "with the marker IDs.", call. = FALSE)
}

# Character matrix of 0/1/2 (or -1/0/1) codes to integer ALT counts.
codes_to_dosage <- function(vals, notes) {
  v <- trimws(vals)
  v[v %in% c("", "NA", "N", "NN", ".", "./.", "NaN")] <- NA
  num <- suppressWarnings(as.numeric(v))
  bad <- !is.na(v) & is.na(num)
  if (any(bad)) stop("Genotype codes must be 0/1/2 (or -1/0/1); found '",
                     paste(utils::head(unique(v[bad]), 3), collapse = "', '"),
                     "'. Use --format vcf or hapmap for allele calls.", call. = FALSE)
  u <- sort(unique(num[!is.na(num)]))
  if (any(u %% 1 != 0)) stop("Genotype codes must be hard calls (0/1/2); found non-integer values such as ",
                             u[u %% 1 != 0][1], ". Imputed dosages cannot be used for call-rate and ",
                             "heterozygosity QC.", call. = FALSE)
  if (any(u == -1)) {
    if (!all(u %in% c(-1, 0, 1)))
      stop("Genotype codes mix -1 and 2; use either 0/1/2 or -1/0/1 coding.", call. = FALSE)
    num <- num + 1
    notes <- c(notes, "Codes -1/0/1 (rrBLUP style) were converted to 0/1/2 by adding 1; the counted allele is the one coded 1.")
  } else if (!all(u %in% c(0, 1, 2))) {
    stop("Unexpected genotype codes: ", paste(setdiff(u, c(0, 1, 2)), collapse = ", "),
         ". Allowed: 0, 1, 2 and missing (NA or empty).", call. = FALSE)
  }
  d <- matrix(as.integer(num), nrow(vals), ncol(vals), dimnames = dimnames(vals))
  list(geno = d, notes = notes)
}

read_map_file <- function(path) {
  m <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
  im <- find_col(names(m), ALIASES$marker); ic <- find_col(names(m), ALIASES$chr)
  ip <- find_col(names(m), ALIASES$pos)
  if (is.na(im) || is.na(ic) || is.na(ip))
    stop("The --map file needs columns marker, chr and pos. Found: ", paste(names(m), collapse = ", "), call. = FALSE)
  make_map(m[[im]], m[[ic]], m[[ip]])
}

read_matrix_geno <- function(path, orientation = "auto", map_path = "") {
  l1 <- first_line(path)
  n_tab <- sum(gregexpr("\t", l1, fixed = TRUE)[[1]] > 0)
  n_com <- sum(gregexpr(",", l1, fixed = TRUE)[[1]] > 0)
  sep <- if (n_tab > n_com) "\t" else ","
  con <- open_text(path)
  df <- tryCatch(utils::read.table(con, header = TRUE, sep = sep, colClasses = "character",
                                   check.names = FALSE, quote = "\"", comment.char = "",
                                   na.strings = character(), stringsAsFactors = FALSE),
                 finally = close(con))
  notes <- character()
  nm <- names(df)
  ic <- find_col(nm, ALIASES$chr); ip <- find_col(nm, ALIASES$pos)
  map_file <- if (nzchar(map_path)) read_map_file(map_path) else NULL
  if (!is.na(ic) && !is.na(ip)) {
    im <- find_col(nm, ALIASES$marker)
    if (is.na(im)) im <- 1L
    ir <- find_col(nm, ALIASES$ref); ia <- find_col(nm, ALIASES$alt)
    meta <- unique(stats::na.omit(c(im, ic, ip, ir, ia)))
    vals <- as.matrix(df[, -meta, drop = FALSE])
    u <- unique_ids(df[[im]], "marker", notes); ids <- u$ids; notes <- u$notes
    rownames(vals) <- ids
    map <- make_map(ids, df[[ic]], df[[ip]], if (!is.na(ir)) df[[ir]] else NA,
                    if (!is.na(ia)) df[[ia]] else NA)
    orient <- list(orientation = "markers-in-rows", basis = "the file has chr and pos columns")
  } else {
    row_ids <- df[[1]]
    vals <- as.matrix(df[, -1, drop = FALSE])
    orient <- decide_orientation(row_ids, colnames(vals), orientation, map_file)
    if (orient$orientation == "samples-in-rows") {
      vals <- t(vals)
      colnames(vals) <- row_ids
    } else {
      rownames(vals) <- row_ids
    }
    u <- unique_ids(rownames(vals), "marker", notes); rownames(vals) <- u$ids; notes <- u$notes
    map <- make_map(rownames(vals))
    if (!is.null(map_file)) {
      i <- match(map$marker, map_file$marker)
      if (anyNA(i)) notes <- c(notes, sprintf("%d marker(s) are not in the --map file (chr/pos missing).", sum(is.na(i))))
      map$chr <- map_file$chr[i]
      map$pos <- map_file$pos[i]
    }
  }
  u <- unique_ids(colnames(vals), "sample", notes); colnames(vals) <- u$ids; notes <- u$notes
  cv <- codes_to_dosage(vals, notes)
  list(geno = cv$geno, map = map, notes = cv$notes, format = "matrix",
       orientation = sprintf("%s (%s)", orient$orientation, orient$basis))
}

read_genotypes <- function(path, format = "auto", orientation = "auto", map_path = "") {
  fmt <- if (format == "auto") detect_format(path) else format
  res <- switch(fmt,
                vcf = read_vcf_geno(path),
                hapmap = read_hapmap_geno(path),
                matrix = read_matrix_geno(path, orientation, map_path),
                stop("--format must be auto, vcf, hapmap or matrix.", call. = FALSE))
  if (nrow(res$geno) < 2L || ncol(res$geno) < 2L)
    stop("Need at least 2 markers and 2 samples; got ", nrow(res$geno), " x ", ncol(res$geno), ".", call. = FALSE)
  res
}

# --- QC statistics ------------------------------------------------------------------

marker_stats <- function(G) {
  called <- !is.na(G)
  n_called <- rowSums(called)
  s <- rowSums(G, na.rm = TRUE)
  het <- rowSums(G == 1L, na.rm = TRUE)
  alt_freq <- ifelse(n_called > 0, s / (2 * n_called), NA_real_)
  maf <- pmin(alt_freq, 1 - alt_freq)
  data.frame(n_called = n_called, call_rate = n_called / ncol(G), alt_freq = alt_freq,
             maf = maf, het_obs = ifelse(n_called > 0, het / n_called, NA_real_),
             monomorphic = !is.na(maf) & maf == 0, row.names = NULL)
}

sample_stats <- function(G) {
  called <- !is.na(G)
  n_called <- colSums(called)
  het <- colSums(G == 1L, na.rm = TRUE)
  data.frame(sample = colnames(G), n_called = n_called, call_rate = n_called / nrow(G),
             het_obs = ifelse(n_called > 0, het / n_called, NA_real_), row.names = NULL,
             stringsAsFactors = FALSE)
}

# Heterozygosity flags. Inbred: absolute ceiling. Outbred: mean +/- k SD
# (Anderson et al. 2010).
flag_het <- function(het, germplasm, max_het = 0.10, het_sd = 3) {
  if (germplasm == "inbred") return(!is.na(het) & het > max_het)
  mu <- mean(het, na.rm = TRUE)
  s <- stats::sd(het, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(FALSE, length(het)))
  !is.na(het) & abs(het - mu) > het_sd * s
}

# Identity-by-state similarity between samples (rows of X, coded 0/1/2/NA):
# 1 - mean |g_i - g_j| / 2 over markers called in both samples.
ibs_similarity <- function(X) {
  ind <- function(v) { a <- X == v; a[is.na(a)] <- FALSE; a + 0 }
  A0 <- ind(0L); A1 <- ind(1L); A2 <- ind(2L)
  N <- tcrossprod(A0 + A1 + A2)
  t01 <- tcrossprod(A0, A1); t12 <- tcrossprod(A1, A2); t02 <- tcrossprod(A0, A2)
  D <- (t01 + t(t01) + t12 + t(t12)) + 2 * (t02 + t(t02))
  ibs <- 1 - D / (2 * N)
  ibs[N == 0] <- NA
  dimnames(ibs) <- dimnames(N) <- list(rownames(X), rownames(X))
  list(ibs = ibs, n_shared = N)
}

duplicate_pairs <- function(ibs, n_shared, threshold = 0.98, min_shared = 100) {
  up <- upper.tri(ibs) & !is.na(ibs) & ibs >= threshold & n_shared >= min_shared
  w <- which(up, arr.ind = TRUE)
  if (!nrow(w)) return(data.frame(sample_1 = character(), sample_2 = character(),
                                  ibs = numeric(), n_shared_markers = numeric()))
  ids <- rownames(ibs)
  out <- data.frame(sample_1 = ids[w[, 1]], sample_2 = ids[w[, 2]], ibs = ibs[w],
                    n_shared_markers = n_shared[w], stringsAsFactors = FALSE)
  out[order(-out$ibs), , drop = FALSE]
}

# Connected groups of duplicate pairs; returns a group label per sample (NA = none).
duplicate_groups <- function(samples, pairs) {
  grp <- stats::setNames(rep(NA_integer_, length(samples)), samples)
  if (!nrow(pairs)) return(grp)
  lab <- stats::setNames(seq_along(samples), samples)
  repeat {
    changed <- FALSE
    for (k in seq_len(nrow(pairs))) {
      a <- lab[[pairs$sample_1[k]]]; b <- lab[[pairs$sample_2[k]]]
      if (a != b) { lab[lab == max(a, b)] <- min(a, b); changed <- TRUE }
    }
    if (!changed) break
  }
  inv <- unique(c(pairs$sample_1, pairs$sample_2))
  grp[inv] <- as.integer(factor(lab[inv], levels = unique(lab[inv])))
  grp
}

# Segregation-distortion chi-square tests for biparental populations.
# Classes: A = homozygous for parent 1 allele, B = parent 2, H = heterozygous.
# Without parents, A = 0 and B = 2 copies of the ALT allele.
SEG_RATIOS <- list(ril = c(A = 0.5, B = 0.5), dh = c(A = 0.5, B = 0.5),
                   f2 = c(A = 0.25, H = 0.5, B = 0.25), bc1 = c(A = 0.5, H = 0.5))

segregation_test <- function(G, population, p1 = NULL, p2 = NULL, alpha = 0.05, min_n = 10) {
  population <- tolower(population)
  if (!population %in% names(SEG_RATIOS))
    stop("--population must be ril, dh, f2 or bc1.", call. = FALSE)
  if (population == "bc1" && is.null(p1))
    stop("A backcross test needs the parents: give --parents RECURRENT,DONOR.", call. = FALSE)
  if (!is.null(p1)) {
    informative <- !is.na(p1) & !is.na(p2) & p1 %in% c(0L, 2L) & p2 %in% c(0L, 2L) & p1 != p2
    nA <- rowSums(G == p1, na.rm = TRUE)
    nB <- rowSums(G == p2, na.rm = TRUE)
  } else {
    nA <- rowSums(G == 0L, na.rm = TRUE)
    nB <- rowSums(G == 2L, na.rm = TRUE)
  }
  nH <- rowSums(G == 1L, na.rm = TRUE)
  # Without parents, a marker monomorphic in the progeny cannot be told apart
  # from one that is monomorphic between the parents, so it is not tested.
  if (is.null(p1)) informative <- ((nA > 0) + (nB > 0) + (nH > 0)) >= 2L
  ratio <- SEG_RATIOS[[population]]
  obs <- switch(population, ril = cbind(A = nA, B = nB), dh = cbind(A = nA, B = nB),
                f2 = cbind(A = nA, H = nH, B = nB), bc1 = cbind(A = nA, H = nH))
  n <- rowSums(obs)
  E <- outer(n, ratio)
  chisq <- rowSums((obs - E)^2 / E)
  df <- length(ratio) - 1L
  tested <- informative & n >= min_n
  p <- ifelse(tested, stats::pchisq(chisq, df, lower.tail = FALSE), NA_real_)
  n_tested <- sum(tested)
  thr <- if (n_tested > 0) alpha / n_tested else NA_real_
  data.frame(seg_n_A = nA, seg_n_H = nH, seg_n_B = nB,
             seg_expected = paste(names(ratio), ratio, sep = "=", collapse = ";"),
             seg_chisq = ifelse(tested, chisq, NA_real_), seg_df = df, seg_p = p,
             seg_tested = tested, distorted = tested & !is.na(p) & p < thr, row.names = NULL)
}

# Non-parental calls: at markers where both parents are homozygous for the same
# allele, every progeny genotype must equal the parental genotype.
nonparental_check <- function(G, p1, p2) {
  mono <- !is.na(p1) & !is.na(p2) & p1 %in% c(0L, 2L) & p1 == p2
  Gm <- G[mono, , drop = FALSE]
  called <- !is.na(Gm)
  bad <- called & (Gm != p1[mono])
  list(n_markers = sum(mono),
       sample = data.frame(sample = colnames(G), nonparental_n_checked = colSums(called),
                           nonparental_n = colSums(bad),
                           nonparental_rate = ifelse(colSums(called) > 0, colSums(bad) / colSums(called), NA_real_),
                           stringsAsFactors = FALSE, row.names = NULL),
       marker_rate = { r <- rep(NA_real_, nrow(G)); r[mono] <- ifelse(rowSums(called) > 0, rowSums(bad) / rowSums(called), NA_real_); r })
}

# --- Writers ----------------------------------------------------------------------

write_marker_table <- function(G, map, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  pos <- ifelse(is.na(map$pos), NA_character_, format(map$pos, scientific = FALSE, trim = TRUE))
  df <- data.frame(marker = map$marker, chr = map$chr, pos = pos, ref = map$ref,
                   alt = map$alt, stringsAsFactors = FALSE, check.names = FALSE)
  df <- cbind(df, as.data.frame(G, optional = TRUE, stringsAsFactors = FALSE))
  utils::write.csv(df, path, row.names = FALSE, na = "")
  invisible(path)
}

write_vcf <- function(G, map, path, source = paste(SKILL, VERSION)) {
  if (anyNA(map$ref) || anyNA(map$alt) || any(!nzchar(map$ref)) || any(!nzchar(map$alt)))
    stop("Writing a VCF needs REF and ALT alleles for every marker; the input did not provide them.", call. = FALSE)
  if (anyNA(map$chr) || anyNA(map$pos))
    stop("Writing a VCF needs chromosome and position for every marker.", call. = FALSE)
  gt <- matrix("./.", nrow(G), ncol(G))
  gt[!is.na(G) & G == 0L] <- "0/0"
  gt[!is.na(G) & G == 1L] <- "0/1"
  gt[!is.na(G) & G == 2L] <- "1/1"
  header <- c("##fileformat=VCFv4.2", paste0("##source=", source),
              "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
              paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                      colnames(G)), collapse = "\t"))
  body <- paste(map$chr, format(map$pos, scientific = FALSE, trim = TRUE), map$marker,
                map$ref, map$alt, ".", "PASS", ".", "GT",
                apply(gt, 1, paste, collapse = "\t"), sep = "\t")
  con <- if (grepl("\\.gz$", path)) gzfile(path, "w") else file(path, "w")
  on.exit(close(con))
  writeLines(c(header, body), con)
  invisible(path)
}

# --- Plots -----------------------------------------------------------------------------

plot_marker_qc <- function(ms, out, min_cr, min_maf) {
  open_png(out, "marker_qc.png", width = 9, height = 4)
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), col.axis = COL_INK2, fg = COL_INK2)
  graphics::hist(ms$call_rate, breaks = 40, col = COL_FILL, border = "white",
                 main = "Marker call rate", xlab = "Call rate", ylab = "Markers", col.main = COL_INK)
  graphics::abline(v = min_cr, col = COL_ORANGE, lwd = 2, lty = 2)
  graphics::hist(ms$maf, breaks = seq(0, 0.5, by = 0.0125), col = COL_FILL, border = "white",
                 main = "Minor allele frequency", xlab = "MAF", ylab = "Markers", col.main = COL_INK)
  graphics::abline(v = min_maf, col = COL_ORANGE, lwd = 2, lty = 2)
  graphics::par(op)
  grDevices::dev.off()
}

plot_sample_qc <- function(ss, out, min_cr, germplasm, max_het) {
  open_png(out, "sample_qc.png", width = 7, height = 5)
  op <- graphics::par(mar = c(4.5, 4.5, 3, 1), col.axis = COL_INK2, fg = COL_INK2)
  flag <- ss$fail_call_rate | ss$het_flag | !is.na(ss$duplicate_group) |
    (if ("nonparental_flag" %in% names(ss)) ss$nonparental_flag else FALSE)
  graphics::plot(ss$call_rate, ss$het_obs, pch = 19, cex = 0.9,
                 col = ifelse(flag, COL_ORANGE, COL_BLUE), xlab = "Sample call rate",
                 ylab = "Observed heterozygosity", main = "Sample call rate and heterozygosity",
                 col.main = COL_INK, xlim = range(c(ss$call_rate, min_cr), na.rm = TRUE),
                 ylim = range(c(0, ss$het_obs, if (germplasm == "inbred") max_het), na.rm = TRUE))
  graphics::abline(v = min_cr, col = COL_MUTED, lty = 2)
  if (germplasm == "inbred") graphics::abline(h = max_het, col = COL_MUTED, lty = 2)
  graphics::legend("topleft", bty = "n", pch = 19, col = c(COL_BLUE, COL_ORANGE),
                   legend = c("passes", "flagged (call rate, heterozygosity, duplicate or parents)"),
                   text.col = COL_INK2)
  graphics::par(op)
  grDevices::dev.off()
}

plot_ibs <- function(ibs, out, thr) {
  v <- ibs[upper.tri(ibs)]
  v <- v[!is.na(v)]
  if (!length(v)) return(invisible(NULL))
  open_png(out, "ibs_pairs.png", width = 7, height = 4.5)
  op <- graphics::par(mar = c(4.5, 4.5, 3, 1), col.axis = COL_INK2, fg = COL_INK2)
  graphics::hist(v, breaks = 60, col = COL_FILL, border = "white", main = "Pairwise identity-by-state",
                 xlab = "IBS similarity", ylab = "Sample pairs", col.main = COL_INK, xlim = c(min(v), 1))
  graphics::abline(v = thr, col = COL_ORANGE, lwd = 2, lty = 2)
  graphics::par(op)
  grDevices::dev.off()
}

plot_segregation <- function(ms, out) {
  ok <- !is.na(ms$seg_p)
  if (!any(ok)) return(invisible(NULL))
  open_png(out, "segregation.png", width = 9, height = 4)
  op <- graphics::par(mar = c(4.5, 4.5, 3, 1), col.axis = COL_INK2, fg = COL_INK2)
  y <- -log10(pmax(ms$seg_p, 1e-300))
  x <- seq_len(nrow(ms))
  graphics::plot(x[ok], y[ok], pch = 19, cex = 0.6, col = ifelse(ms$distorted[ok], COL_ORANGE, COL_BLUE),
                 xlab = "Marker (file order)", ylab = expression(-log[10](p)),
                 main = "Segregation distortion (chi-square)", col.main = COL_INK)
  thr <- 0.05 / sum(ok)
  graphics::abline(h = -log10(thr), col = COL_MUTED, lty = 2)
  graphics::legend("topright", bty = "n", pch = 19, col = c(COL_BLUE, COL_ORANGE),
                   legend = c("consistent with expected ratio", "distorted (Bonferroni)"), text.col = COL_INK2)
  graphics::par(op)
  grDevices::dev.off()
}

# --- Main QC pipeline --------------------------------------------------------------------

check_options <- function(opt) {
  if (!opt$germplasm %in% c("inbred", "outbred"))
    stop("Tell the skill what the samples are: --germplasm inbred (inbred lines, RILs, DH, ",
         "pure-line cultivars) or --germplasm outbred (open-pollinated, hybrids, clones, F2 ",
         "individuals). Ask the user if unsure; heterozygosity QC depends on it.", call. = FALSE)
  for (k in c("min_marker_call_rate", "min_sample_call_rate", "min_maf", "max_het", "ibs_dup",
              "max_nonparental", "seg_alpha")) {
    if (is.na(opt[[k]]) || opt[[k]] < 0 || opt[[k]] > 1)
      stop("--", gsub("_", "-", k), " must be between 0 and 1.", call. = FALSE)
  }
  if (opt$min_maf > 0.5) stop("--min-maf must be at most 0.5.", call. = FALSE)
  if (nzchar(opt$population) && !tolower(opt$population) %in% names(SEG_RATIOS))
    stop("--population must be ril, dh, f2 or bc1.", call. = FALSE)
  invisible(TRUE)
}

run_qc <- function(geno, map, opt) {
  G <- geno
  samples <- colnames(G)
  parents <- if (nzchar(opt$parents)) trimws(strsplit(opt$parents, ",", fixed = TRUE)[[1]]) else character()
  if (length(parents) && length(parents) != 2L)
    stop("--parents needs exactly two sample IDs separated by a comma.", call. = FALSE)
  miss_par <- setdiff(parents, samples)
  if (length(miss_par)) stop("Parent sample(s) not found in the genotype file: ",
                             paste(miss_par, collapse = ", "), call. = FALSE)

  # 1. Sample call rate on all markers.
  ss <- sample_stats(G)
  ss$fail_call_rate <- ss$call_rate < opt$min_sample_call_rate
  keep_s <- !ss$fail_call_rate | ss$sample %in% parents
  if (sum(keep_s) < 2L) stop("Fewer than 2 samples pass the sample call-rate threshold.", call. = FALSE)

  # 2. Marker statistics on samples that pass call rate (parents excluded from frequencies).
  stat_s <- keep_s & !ss$sample %in% parents
  ms <- cbind(map, marker_stats(G[, stat_s, drop = FALSE]))
  ms$fail_call_rate <- ms$call_rate < opt$min_marker_call_rate
  ms$fail_maf <- is.na(ms$maf) | ms$maf < opt$min_maf
  pass_m <- !ms$fail_call_rate & !ms$fail_maf

  # 3. Sample heterozygosity on markers that pass.
  if (sum(pass_m) < 1L) stop("No markers pass the call-rate and MAF thresholds.", call. = FALSE)
  Gp <- G[pass_m, , drop = FALSE]
  het <- sample_stats(Gp)
  ss$het_obs <- het$het_obs
  ss$het_flag <- FALSE
  ss$het_flag[stat_s] <- flag_het(ss$het_obs[stat_s], opt$germplasm, opt$max_het, opt$het_sd)

  # 4. Duplicates by identity-by-state on a marker subset.
  set.seed(opt$seed)
  idx <- which(pass_m)
  if (length(idx) > opt$ibs_markers) idx <- sort(sample(idx, opt$ibs_markers))
  ib <- ibs_similarity(t(G[idx, keep_s, drop = FALSE]))
  pairs <- duplicate_pairs(ib$ibs, ib$n_shared, opt$ibs_dup, min(opt$min_shared, length(idx)))
  ss$duplicate_group <- unname(duplicate_groups(ss$sample, pairs))

  # 5. Biparental checks.
  pop <- tolower(opt$population)
  np <- NULL
  if (length(parents)) {
    p1 <- G[, parents[1]]; p2 <- G[, parents[2]]
    prog <- ss$sample[keep_s & !ss$sample %in% parents]
    np <- nonparental_check(G[, prog, drop = FALSE], p1, p2)
    ss <- merge(ss, np$sample, by = "sample", all.x = TRUE, sort = FALSE)
    ss <- ss[match(samples, ss$sample), , drop = FALSE]
    ss$nonparental_flag <- !is.na(ss$nonparental_rate) & ss$nonparental_rate > opt$max_nonparental
    ms$nonparental_call_rate <- np$marker_rate
    ms$parent1 <- p1; ms$parent2 <- p2
    ms$parents_polymorphic <- !is.na(p1) & !is.na(p2) & p1 %in% c(0L, 2L) & p2 %in% c(0L, 2L) & p1 != p2
  }
  if (nzchar(pop)) {
    prog <- ss$sample[keep_s & !ss$sample %in% parents]
    if (length(parents) && "nonparental_flag" %in% names(ss))
      prog <- setdiff(prog, ss$sample[ss$nonparental_flag])
    st <- segregation_test(G[, prog, drop = FALSE], pop,
                           if (length(parents)) G[, parents[1]] else NULL,
                           if (length(parents)) G[, parents[2]] else NULL, opt$seg_alpha)
    ms <- cbind(ms, st)
  }

  # 6. Final keep decisions.
  drop_dist <- if (nzchar(pop) && opt$drop_distorted) ms$distorted else rep(FALSE, nrow(ms))
  pass_m <- pass_m & !drop_dist
  ms$pass <- pass_m
  ms$fail_reason <- gsub("\\s+", " ", trimws(paste(
    ifelse(ms$fail_call_rate, "call_rate", ""),
    ifelse(ms$fail_maf, ifelse(ms$monomorphic, "monomorphic", "maf"), ""),
    ifelse(drop_dist, "distorted", ""))))
  ss$keep <- !ss$fail_call_rate
  if (opt$drop_het_outliers) ss$keep <- ss$keep & !ss$het_flag
  if (opt$drop_duplicates && nrow(pairs)) {
    for (g in unique(stats::na.omit(ss$duplicate_group))) {
      mem <- which(ss$duplicate_group == g)
      best <- mem[which.max(ss$call_rate[mem])]
      ss$keep[setdiff(mem, best)] <- FALSE
    }
  }
  if (opt$drop_nonparental && "nonparental_flag" %in% names(ss)) ss$keep <- ss$keep & !ss$nonparental_flag
  ss$keep[ss$sample %in% parents] <- TRUE
  ss$fail_reason <- gsub("\\s+", " ", trimws(paste(
    ifelse(ss$fail_call_rate, "call_rate", ""),
    ifelse(ss$het_flag, "heterozygosity", ""),
    ifelse(!is.na(ss$duplicate_group), "duplicate", ""),
    if ("nonparental_flag" %in% names(ss)) ifelse(ss$nonparental_flag, "non_parental", "") else "")))
  list(markers = ms, samples = ss, pairs = pairs, ibs = ib$ibs, n_ibs_markers = length(idx),
       parents = parents, population = pop, n_parental_mono = if (!is.null(np)) np$n_markers else NA)
}

fmt_list <- function(x, n = 10) {
  if (!length(x)) return("none")
  s <- paste(utils::head(x, n), collapse = ", ")
  if (length(x) > n) s <- paste0(s, sprintf(" and %d more", length(x) - n))
  s
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  options(scipen = 100)  # positions as 1000000, not 1e+06, in CSV output
  if (opt$demo) {
    opt$input <- file.path(SKILL_DIR, "demo", "demo_panel.vcf.gz")
    if (!nzchar(opt$germplasm)) opt$germplasm <- "inbred"
  }
  if (!nzchar(opt$input)) stop("Provide --input FILE (VCF, HapMap or 0/1/2 matrix), or --demo. See --help.", call. = FALSE)
  check_options(opt)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

  gd <- read_genotypes(opt$input, opt$format, opt$orientation, opt$map)
  res <- run_qc(gd$geno, gd$map, opt)
  ms <- res$markers; ss <- res$samples

  keep_m <- ms$pass; keep_s <- ss$keep
  Gf <- gd$geno[keep_m, keep_s, drop = FALSE]
  mapf <- gd$map[keep_m, , drop = FALSE]
  write_marker_table(Gf, mapf, file.path(opt$out, "tables", "genotypes_filtered_012.csv"))
  vcf_note <- NULL
  if (opt$write_vcf) {
    vcf_note <- tryCatch({
      write_vcf(Gf, mapf, file.path(opt$out, "genotypes_filtered.vcf"))
      "genotypes_filtered.vcf written."
    }, error = function(e) paste("VCF not written:", conditionMessage(e)))
  }
  save_table(ms, opt$out, "marker_qc.csv")
  save_table(ss, opt$out, "sample_qc.csv")
  save_table(res$pairs, opt$out, "duplicate_pairs.csv")
  summ <- data.frame(
    item = c("samples_in", "markers_in", "samples_failing_call_rate", "samples_het_flagged",
             "duplicate_pairs", "markers_failing_call_rate", "markers_monomorphic",
             "markers_low_maf", "samples_out", "markers_out"),
    value = c(nrow(ss), nrow(ms), sum(ss$fail_call_rate), sum(ss$het_flag), nrow(res$pairs),
              sum(ms$fail_call_rate), sum(ms$monomorphic & !ms$fail_call_rate),
              sum(ms$fail_maf & !ms$monomorphic & !ms$fail_call_rate), sum(keep_s), sum(keep_m)))
  if (nzchar(res$population)) summ <- rbind(summ, data.frame(item = "markers_distorted", value = sum(ms$distorted)))
  if ("nonparental_flag" %in% names(ss)) summ <- rbind(summ, data.frame(item = "samples_nonparental", value = sum(ss$nonparental_flag)))
  save_table(summ, opt$out, "qc_summary.csv")

  plot_marker_qc(ms, opt$out, opt$min_marker_call_rate, opt$min_maf)
  plot_sample_qc(ss, opt$out, opt$min_sample_call_rate, opt$germplasm, opt$max_het)
  plot_ibs(res$ibs, opt$out, opt$ibs_dup)
  if (nzchar(res$population)) plot_segregation(ms, opt$out)

  # Report
  warn <- character()
  if (any(ss$fail_call_rate)) warn <- c(warn, sprintf("%d sample(s) have call rate below %.2f and were removed: %s.",
                                                      sum(ss$fail_call_rate), opt$min_sample_call_rate, fmt_list(ss$sample[ss$fail_call_rate])))
  if (any(ss$het_flag)) warn <- c(warn, sprintf(
    "%d sample(s) flagged for heterozygosity (%s): %s. %s",
    sum(ss$het_flag),
    if (opt$germplasm == "inbred") sprintf("above %.2f in inbred material: possible DNA contamination, seed mix or outcrossing", opt$max_het)
    else sprintf("outside mean +/- %g SD", opt$het_sd),
    fmt_list(ss$sample[ss$het_flag]),
    if (opt$drop_het_outliers) "They were removed." else "They were KEPT in the filtered matrix; ask the user whether to drop them (--drop-het-outliers)."))
  if (nrow(res$pairs)) warn <- c(warn, sprintf(
    "%d sample pair(s) are near-identical (IBS >= %.2f): %s. %s",
    nrow(res$pairs), opt$ibs_dup, fmt_list(paste(res$pairs$sample_1, res$pairs$sample_2, sep = "~"), 8),
    if (opt$drop_duplicates) "Only the best-called sample of each group was kept." else
      "All were KEPT; confirm with the user whether these are true duplicates, sister lines or mislabels."))
  if (nzchar(res$population) && any(ms$distorted)) warn <- c(warn, sprintf(
    "%d marker(s) show segregation distortion (Bonferroni, family-wise alpha %.2f). %s",
    sum(ms$distorted), opt$seg_alpha,
    if (opt$drop_distorted) "They were removed." else "They were KEPT; distortion can be real biology, so ask before dropping."))
  if ("nonparental_flag" %in% names(ss) && any(ss$nonparental_flag)) warn <- c(warn, sprintf(
    "%d progeny have more than %.0f%% non-parental calls (likely not from this cross): %s.",
    sum(ss$nonparental_flag), 100 * opt$max_nonparental, fmt_list(ss$sample[ss$nonparental_flag])))
  if (all(is.na(gd$map$chr)) || all(is.na(gd$map$pos)))
    warn <- c(warn, "No chromosome/position information: give --map so that GWAS and LD skills can use the output.")
  if (sum(keep_m) < 100) warn <- c(warn, sprintf("Only %d markers remain after QC.", sum(keep_m)))
  if (!is.null(vcf_note) && startsWith(vcf_note, "VCF not written")) warn <- c(warn, vcf_note)

  thr <- data.frame(
    threshold = c("min marker call rate", "min sample call rate", "min MAF",
                  if (opt$germplasm == "inbred") "max heterozygosity (inbred)" else "heterozygosity SD limit (outbred)",
                  "IBS duplicate threshold"),
    value = c(opt$min_marker_call_rate, opt$min_sample_call_rate, opt$min_maf,
              if (opt$germplasm == "inbred") opt$max_het else opt$het_sd, opt$ibs_dup))
  report <- c(
    sprintf("# %s report", SKILL), "",
    if (length(warn)) c("## Warnings", "", paste("-", warn), "") else NULL,
    "## Summary", "",
    sprintf("- Input: `%s` (format: %s; layout: %s)", opt$input, gd$format, gd$orientation),
    sprintf("- Germplasm: %s", opt$germplasm),
    sprintf("- Samples: %d in, %d out. Markers: %d in, %d out.", nrow(ss), sum(keep_s), nrow(ms), sum(keep_m)),
    sprintf("- Genotype coding in all outputs: number of copies of the ALT (counted) allele, 0/1/2, empty = missing."),
    if (length(gd$notes)) paste("- Input note:", gd$notes) else NULL, "",
    "## Thresholds used", "", "| threshold | value |", "|---|---|",
    sprintf("| %s | %s |", thr$threshold, format(thr$value)), "",
    "These are adjustable rules of thumb (see SKILL.md), not universal standards.", "",
    "## Marker QC", "",
    sprintf("- Failing call rate: %d", sum(ms$fail_call_rate)),
    sprintf("- Monomorphic (among samples that pass): %d", sum(ms$monomorphic & !ms$fail_call_rate)),
    sprintf("- Below MAF %.2f (polymorphic): %d", opt$min_maf, sum(ms$fail_maf & !ms$monomorphic & !ms$fail_call_rate)),
    sprintf("- Median marker call rate: %.3f; median MAF: %.3f", stats::median(ms$call_rate), stats::median(ms$maf, na.rm = TRUE)), "",
    "![Marker QC](figures/marker_qc.png)", "",
    "## Sample QC", "",
    sprintf("- Failing call rate: %s", fmt_list(ss$sample[ss$fail_call_rate])),
    sprintf("- Heterozygosity flagged: %s", fmt_list(ss$sample[ss$het_flag])),
    sprintf("- Median sample heterozygosity: %.4f", stats::median(ss$het_obs, na.rm = TRUE)),
    sprintf("- Duplicate check: %d markers, %d pair(s) with IBS >= %.2f", res$n_ibs_markers, nrow(res$pairs), opt$ibs_dup), "",
    "![Sample QC](figures/sample_qc.png)", "", "![IBS](figures/ibs_pairs.png)", "",
    if (length(res$parents) || nzchar(res$population)) c(
      "## Biparental population checks", "",
      if (length(res$parents)) c(
        sprintf("- Parents: %s (P1) and %s (P2)", res$parents[1], res$parents[2]),
        sprintf("- Markers polymorphic between parents: %d; monomorphic between parents (used for the non-parental check): %d",
                sum(ms$parents_polymorphic), res$n_parental_mono),
        sprintf("- Progeny with > %.0f%% non-parental calls: %s", 100 * opt$max_nonparental,
                fmt_list(ss$sample[ss$nonparental_flag]))) else NULL,
      if (nzchar(res$population)) c(
        sprintf("- Population: %s; expected ratio %s; markers tested: %d; distorted (Bonferroni): %d",
                res$population, ms$seg_expected[1], sum(ms$seg_tested), sum(ms$distorted)),
        "", "![Segregation](figures/segregation.png)") else NULL, "") else NULL,
    "## Outputs", "",
    "- `tables/genotypes_filtered_012.csv`: markers in rows (marker, chr, pos, ref, alt), one column per sample. Input for claw-gwas-crop and claw-ld-decay.",
    "- `tables/marker_qc.csv`, `tables/sample_qc.csv`, `tables/duplicate_pairs.csv`, `tables/qc_summary.csv`",
    if (!is.null(vcf_note)) paste("-", vcf_note) else NULL, "",
    "## Limitations", "",
    "- Diploid, biallelic markers only; multi-allelic sites and polyploid dosage calls are not handled.",
    "- Missing calls are left missing (no imputation). Impute before GWAS if missingness is substantial.",
    "- Thresholds are rules of thumb; flagged samples are reported for review, not proof of contamination or duplication.",
    "- Hardy-Weinberg tests are not applied: they are not meaningful for inbred lines or structured crop panels.", "",
    "## References", "",
    "- Anderson CA et al. (2010) Data quality control in genetic case-control association studies. Nat Protoc 5:1564-1573. https://doi.org/10.1038/nprot.2010.116",
    "- Broman KW, Sen S (2009) A Guide to QTL Mapping with R/qtl. Springer. https://doi.org/10.1007/978-0-387-92125-9",
    "- Knaus BJ, Grunwald NJ (2017) vcfR: a package to manipulate and visualize variant call format data in R. Mol Ecol Resour 17:44-53. https://doi.org/10.1111/1755-0998.12549", "",
    sprintf("*%s %s; seed %d*", SKILL, VERSION, opt$seed))
  write_utf8(report, file.path(opt$out, "report.md"))
  write_repro_bundle(opt$out, argv, inputs = c(opt$input, opt$map))
  cat(sprintf("Done: %d of %d samples and %d of %d markers kept. Report: %s\n",
              sum(keep_s), nrow(ss), sum(keep_m), nrow(ms),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(res)
}

if (sys.nframe() == 0L) main()
