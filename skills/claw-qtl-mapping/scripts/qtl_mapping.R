#!/usr/bin/env Rscript
# claw-qtl-mapping: single-QTL genome scans for biparental and multiparent
# populations with R/qtl2: genotype probabilities (HMM with genotyping error),
# Haley-Knott regression or a LOCO linear mixed model, permutation thresholds,
# LOD support or Bayes credible intervals, QTL effects and variance explained.
#
# Usage:
#   Rscript scripts/qtl_mapping.R --demo
#   Rscript scripts/qtl_mapping.R --geno geno.csv --map gmap.csv --pheno blues.csv \
#       --cross-type riself
#   Rscript scripts/qtl_mapping.R --control magic.yaml --kinship loco
#   Rscript scripts/qtl_mapping.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT license.

SKILL   <- "claw-qtl-mapping"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/qtl_mapping.R [--demo] [--control FILE.yaml|.json|.zip]
            [--geno FILE --map FILE --pheno FILE --cross-type TYPE] [--pmap FILE]
            [--geno-codes CODES] [--covar FILE] [--addcovar COLS] [--id-col COL]
            [--map-units cM|Mb|bp] [--cm-per-mb X] [--traits COLS]
            [--error-prob P] [--step CM] [--map-function haldane|kosambi]
            [--kinship none|loco] [--n-perm N] [--interval lod|bayes]
            [--lod-drop D] [--bayes-prob P] [--peakdrop D] [--seed S]
            [--cores N] [--out DIR]

  --demo          run on the simulated RIL population in demo/ (CSV input)
  --control       R/qtl2 control file or zip (any cross type, including multiparent)
  --geno          genotype CSV: id + one column per marker (biparental crosses)
  --map           genetic map CSV: marker, chr, pos (cM unless --map-units)
  --pmap          optional physical map CSV: marker, chr, pos (Mb), to report intervals in Mb
  --pheno         phenotype CSV: id + one column per trait (e.g. BLUEs)
  --cross-type    riself, risib, f2, bc, dh or haploid (CSV input)
  --geno-codes    genotype codes, e.g. 'A=1,B=2' (riself, risib, dh, haploid; default),
                  'A=1,H=2,B=3' (f2; default), 'A=1,H=2' (bc; default)
  --covar         covariate CSV: id + columns (CSV input); --addcovar names columns to use
  --id-col        ID column name in --pheno and --covar (default: id)
  --map-units     units of --map positions: cM (default), Mb or bp
  --cm-per-mb     genome-average cM per Mb to convert a physical map (rough; see SKILL.md)
  --traits        comma-separated traits (default: all)
  --error-prob    genotyping error probability for calc_genoprob (default 0.002)
  --step          pseudomarker spacing in cM; 0 = markers only (default 1)
  --map-function  haldane (default) or kosambi
  --kinship       none = Haley-Knott regression (default); loco = LMM with LOCO kinship
  --n-perm        permutations for genome-wide thresholds (default 1000)
  --interval      lod (LOD support interval, default) or bayes (Bayes credible interval)
  --lod-drop      LOD drop for support intervals (default 1.5)
  --bayes-prob    coverage for Bayes credible intervals (default 0.95)
  --peakdrop      LOD drop required between two peaks on one chromosome (default Inf:
                  one peak per chromosome and trait)
  --seed          random seed for permutations (default 20261015)
  --cores         CPU cores for qtl2 (default 1)
  --out           output folder (default: agriclaw_output/claw-qtl-mapping)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, control = "", geno = "", map = "", pmap = "", pheno = "",
                 cross_type = "", geno_codes = "", covar = "", addcovar = "",
                 id_col = "id", map_units = "cM", cm_per_mb = 0, traits = "",
                 error_prob = 0.002, step = 1, map_function = "haldane",
                 kinship = "none", n_perm = 1000, interval = "lod", lod_drop = 1.5,
                 bayes_prob = 0.95, peakdrop = Inf, seed = 20261015, cores = 1, out = "")

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
# Building the cross object
# ---------------------------------------------------------------------------

BIPARENTAL <- c("riself", "risib", "f2", "bc", "dh", "haploid")
DEFAULT_CODES <- list(riself = c(A = 1L, B = 2L), risib = c(A = 1L, B = 2L),
                      dh = c(A = 1L, B = 2L), haploid = c(A = 1L, B = 2L),
                      f2 = c(A = 1L, H = 2L, B = 3L), bc = c(A = 1L, H = 2L))
MISSING_CODES <- c("", "NA", "-", ".", "N", "--", "NN")

parse_codes <- function(spec, cross_type) {
  if (!nzchar(spec)) return(DEFAULT_CODES[[cross_type]])
  parts <- strsplit(strsplit(spec, ",", fixed = TRUE)[[1]], "=", fixed = TRUE)
  if (any(lengths(parts) != 2L))
    stop("--geno-codes must look like 'A=1,H=2,B=3'.", call. = FALSE)
  v <- suppressWarnings(as.integer(vapply(parts, `[`, "", 2L)))
  if (anyNA(v)) stop("--geno-codes values must be integers.", call. = FALSE)
  stats::setNames(v, trimws(vapply(parts, `[`, "", 1L)))
}

read_map_csv <- function(path, what) {
  if (!file.exists(path)) stop(what, " file not found: ", path, call. = FALSE)
  m <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                       colClasses = c("character", "character", "numeric"))
  names(m)[1:3] <- c("marker", "chr", "pos")
  if (anyNA(m$pos)) stop(what, " has missing or non-numeric positions.", call. = FALSE)
  if (anyDuplicated(m$marker))
    stop(what, " lists marker(s) twice: ", paste(utils::head(m$marker[duplicated(m$marker)], 5),
                                                  collapse = ", "), call. = FALSE)
  m
}

# Convert user CSV files into the qtl2 file format, write a control file in
# work_dir and read it with qtl2::read_cross2. Biparental cross types only.
build_cross_from_csv <- function(geno_file, map_file, pheno_file, cross_type,
                                 geno_codes = "", pmap_file = "", covar_file = "",
                                 id_col = "id", map_units = "cM", cm_per_mb = 0,
                                 work_dir = tempfile("qtl2_input")) {
  notes <- character()
  if (!cross_type %in% BIPARENTAL)
    stop("--cross-type must be one of ", paste(BIPARENTAL, collapse = ", "),
         " for CSV input. Multiparent populations (MAGIC, 4-/8-way RILs, DO, ...) need ",
         "founder genotypes and cross information: give an R/qtl2 control file with --control.",
         call. = FALSE)
  codes <- parse_codes(geno_codes, cross_type)
  if (!file.exists(geno_file)) stop("Genotype file not found: ", geno_file, call. = FALSE)
  g <- utils::read.csv(geno_file, check.names = FALSE, colClasses = "character",
                       na.strings = character())
  ids <- trimws(g[[1L]])
  if (anyDuplicated(ids)) stop("Duplicated individual IDs in the genotype file.", call. = FALSE)
  G <- as.matrix(g[, -1L, drop = FALSE])
  G <- trimws(G)
  G[G %in% MISSING_CODES] <- "-"
  bad <- setdiff(unique(as.vector(G)), c(names(codes), "-"))
  if (length(bad))
    stop("Genotype code(s) not covered by --geno-codes (", paste(names(codes), codes, sep = "=",
                                                                 collapse = ","), "): ",
         paste(utils::head(bad, 10), collapse = ", "),
         if (cross_type %in% c("riself", "risib", "dh", "haploid") && any(bad %in% c("H", "AB", "BA", "1")))
           ". Heterozygous calls are not allowed for this cross type: code residual heterozygotes as missing ('-')."
         else ". Give the right codes with --geno-codes, or recode the file.", call. = FALSE)
  map <- read_map_csv(map_file, "Map")
  if (map_units != "cM") {
    if (!map_units %in% c("Mb", "bp"))
      stop("--map-units must be cM, Mb or bp.", call. = FALSE)
    if (cm_per_mb <= 0)
      stop("The map is in ", map_units, ". QTL mapping needs genetic distances (cM): ask the ",
           "user for a genetic map of this population (or a published linkage map). If none ",
           "exists, --cm-per-mb converts with one genome-average rate, which is only a rough ",
           "approximation (recombination varies along chromosomes).", call. = FALSE)
    pmap_from_map <- map
    if (map_units == "bp") pmap_from_map$pos <- pmap_from_map$pos / 1e6
    map$pos <- pmap_from_map$pos * cm_per_mb
    notes <- c(notes, sprintf(paste0(
      "Genetic positions were approximated from the physical map with one genome-average ",
      "rate of %g cM/Mb. QTL positions and intervals in cM, and the HMM genotype probabilities, ",
      "are only as good as that approximation: obtain a genetic map for final results."), cm_per_mb))
  } else {
    if (max(map$pos) > 1000)
      stop("Map positions reach ", signif(max(map$pos), 4), ", which looks like base pairs or Mb, ",
           "not cM. Give --map-units bp or Mb (and a genetic map, or --cm-per-mb).", call. = FALSE)
    pmap_from_map <- NULL
  }
  miss_map <- setdiff(colnames(G), map$marker)
  miss_geno <- setdiff(map$marker, colnames(G))
  if (length(miss_map))
    notes <- c(notes, sprintf("%d genotyped marker(s) have no map position and were not used (e.g. %s).",
                              length(miss_map), paste(utils::head(miss_map, 3), collapse = ", ")))
  if (length(miss_geno))
    notes <- c(notes, sprintf("%d mapped marker(s) are not in the genotype file and were ignored.",
                              length(miss_geno)))
  map <- map[map$marker %in% colnames(G), , drop = FALSE]
  if (nrow(map) < 2L) stop("Fewer than two markers have both genotypes and map positions.", call. = FALSE)
  chr_order <- unique(map$chr)
  num <- suppressWarnings(as.numeric(chr_order))
  if (!anyNA(num)) chr_order <- chr_order[order(num)]
  map <- map[order(match(map$chr, chr_order), map$pos), , drop = FALSE]
  G <- G[, map$marker, drop = FALSE]

  ph <- utils::read.csv(pheno_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!id_col %in% names(ph))
    stop("Phenotype file has no '", id_col, "' column (set --id-col).", call. = FALSE)
  traits <- setdiff(names(ph), id_col)
  nonnum <- traits[!vapply(ph[traits], is.numeric, TRUE)]
  if (length(nonnum))
    stop("Non-numeric phenotype column(s): ", paste(nonnum, collapse = ", "),
         ". Keep only numeric traits (and the ID column) in --pheno.", call. = FALSE)
  ph <- data.frame(id = as.character(ph[[id_col]]), ph[traits], check.names = FALSE)
  if (!any(ph$id %in% ids)) stop("No phenotype ID matches a genotype ID.", call. = FALSE)
  no_g <- setdiff(ph$id, ids)
  if (length(no_g))
    notes <- c(notes, sprintf("%d phenotyped individual(s) have no genotypes and were not used.", length(no_g)))

  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(data.frame(id = ids, G, check.names = FALSE), file.path(work_dir, "geno.csv"),
                   row.names = FALSE, quote = FALSE)
  utils::write.csv(map, file.path(work_dir, "gmap.csv"), row.names = FALSE, quote = FALSE)
  utils::write.csv(ph, file.path(work_dir, "pheno.csv"), row.names = FALSE, na = "NA")
  pm_file <- NULL
  pmap <- if (nzchar(pmap_file)) read_map_csv(pmap_file, "Physical map") else pmap_from_map
  if (nzchar(pmap_file) && max(pmap$pos) > 1e4)
    stop("Physical map positions reach ", signif(max(pmap$pos), 4), ": --pmap must be in Mb ",
         "(divide base pairs by 1e6).", call. = FALSE)
  if (!is.null(pmap)) {
    pmap <- pmap[match(map$marker, pmap$marker), , drop = FALSE]
    if (anyNA(pmap$marker)) stop("The physical map lacks some mapped markers.", call. = FALSE)
    pmap$chr <- map$chr
    utils::write.csv(pmap, file.path(work_dir, "pmap.csv"), row.names = FALSE, quote = FALSE)
    pm_file <- "pmap.csv"
  }
  cv_file <- NULL
  if (nzchar(covar_file)) {
    cv <- utils::read.csv(covar_file, check.names = FALSE, stringsAsFactors = FALSE)
    if (!id_col %in% names(cv)) stop("Covariate file has no '", id_col, "' column.", call. = FALSE)
    names(cv)[names(cv) == id_col] <- "id"
    cv <- cv[, c("id", setdiff(names(cv), "id")), drop = FALSE]
    utils::write.csv(cv, file.path(work_dir, "covar.csv"), row.names = FALSE, na = "NA")
    cv_file <- "covar.csv"
  }
  alleles <- c("A", "B")
  ctrl <- file.path(work_dir, "cross.yaml")
  qtl2::write_control_file(ctrl, crosstype = cross_type, geno_file = "geno.csv",
                           gmap_file = "gmap.csv", pmap_file = pm_file, pheno_file = "pheno.csv",
                           covar_file = cv_file, geno_codes = codes, alleles = alleles,
                           na.strings = c("-", "NA"), overwrite = TRUE,
                           description = "Written by claw-qtl-mapping from CSV input")
  cross <- qtl2::read_cross2(ctrl, quiet = TRUE)
  list(cross = cross, notes = notes, control = ctrl)
}

read_cross_control <- function(path) {
  if (!file.exists(path)) stop("Control file not found: ", path, call. = FALSE)
  list(cross = qtl2::read_cross2(path, quiet = TRUE), notes = character(), control = path)
}

# ---------------------------------------------------------------------------
# Quality control
# ---------------------------------------------------------------------------

expected_genotype_freq <- function(cross_type) {
  switch(cross_type,
         riself = c(0.5, 0.5), risib = c(0.5, 0.5), dh = c(0.5, 0.5), haploid = c(0.5, 0.5),
         bc = c(0.5, 0.5), f2 = c(0.25, 0.5, 0.25), NULL)
}

# Marker and individual missingness, segregation distortion (chi-square against
# the Mendelian ratio of the cross type, autosomes of biparental crosses only)
# and near-duplicate individuals (qtl2::compare_geno).
qc_cross <- function(cross) {
  is_x <- cross$is_x_chr
  exp_f <- expected_genotype_freq(cross$crosstype)
  mk <- do.call(rbind, lapply(names(cross$geno), function(ch) {
    g <- cross$geno[[ch]]
    typed <- colSums(g > 0L)
    p <- rep(NA_real_, ncol(g))
    if (!is.null(exp_f) && !isTRUE(is_x[[ch]])) {
      k <- length(exp_f)
      counts <- matrix(vapply(seq_len(k), function(v) colSums(g == v), numeric(ncol(g))),
                       ncol = k)
      tot <- rowSums(counts)
      E <- outer(tot, exp_f)
      chi <- rowSums((counts - E)^2 / E)
      p <- stats::pchisq(chi, df = k - 1L, lower.tail = FALSE)
      p[tot < 10 | tot != typed] <- NA_real_
    }
    data.frame(marker = colnames(g), chr = ch, pos = unname(cross$gmap[[ch]][colnames(g)]),
               n_typed = unname(typed), missing_rate = unname(1 - typed / nrow(g)),
               freq_geno1 = unname(colSums(g == 1L) / pmax(typed, 1)), seg_dist_p = unname(p),
               stringsAsFactors = FALSE)
  }))
  mk$distorted <- !is.na(mk$seg_dist_p) & mk$seg_dist_p < 0.05 / nrow(mk)
  gall <- do.call(cbind, cross$geno)
  ind <- data.frame(id = rownames(gall), missing_rate = rowMeans(gall == 0L),
                    stringsAsFactors = FALSE)
  dup <- NULL
  if (nrow(gall) >= 2L) {
    cg <- qtl2::compare_geno(cross, quiet = TRUE)
    s <- summary(cg, threshold = 0.95)
    if (nrow(s)) dup <- as.data.frame(s)
  }
  list(markers = mk, individuals = ind, duplicates = dup)
}

skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L) return(NA_real_)
  mean((x - mean(x))^3) / stats::sd(x)^3
}

# ---------------------------------------------------------------------------
# Genome scan, permutations, peaks, effects
# ---------------------------------------------------------------------------

covariate_matrix <- function(cross, names_csv) {
  if (!nzchar(names_csv)) return(NULL)
  nm <- trimws(strsplit(names_csv, ",", fixed = TRUE)[[1]])
  cv <- cross$covar
  if (is.null(cv)) stop("--addcovar given but the cross has no covariates (--covar).", call. = FALSE)
  miss <- setdiff(nm, colnames(cv))
  if (length(miss)) stop("Covariate(s) not found: ", paste(miss, collapse = ", "),
                         ". Available: ", paste(colnames(cv), collapse = ", "), call. = FALSE)
  parts <- lapply(nm, function(v) {
    x <- cv[[v]]
    if (is.numeric(x)) return(matrix(x, ncol = 1L, dimnames = list(NULL, v)))
    f <- factor(x)
    mm <- stats::model.matrix(~ f, data.frame(f = f), na.action = stats::na.pass)[, -1L, drop = FALSE]
    colnames(mm) <- paste0(v, "_", levels(f)[-1L])
    mm
  })
  m <- do.call(cbind, parts)
  rownames(m) <- rownames(cv)
  m
}

run_scan <- function(cross, traits, step = 1, error_prob = 0.002, map_function = "haldane",
                     kinship = "none", addcovar = NULL, cores = 1) {
  map <- if (step > 0) qtl2::insert_pseudomarkers(cross$gmap, step = step) else cross$gmap
  pr <- qtl2::calc_genoprob(cross, map, error_prob = error_prob,
                            map_function = map_function, cores = cores)
  multiparent <- !cross$crosstype %in% BIPARENTAL
  probs <- if (multiparent) qtl2::genoprob_to_alleleprob(pr, cores = cores) else pr
  K <- if (kinship == "loco") qtl2::calc_kinship(probs, "loco", cores = cores) else NULL
  Xcovar <- if (any(cross$is_x_chr)) qtl2::get_x_covar(cross) else NULL
  pheno <- cross$pheno[, traits, drop = FALSE]
  out <- qtl2::scan1(probs, pheno, kinship = K, addcovar = addcovar, Xcovar = Xcovar,
                     cores = cores)
  list(map = map, probs = probs, K = K, Xcovar = Xcovar, addcovar = addcovar,
       out = out, pheno = pheno, multiparent = multiparent,
       is_x = attr(pr, "is_x_chr"), crosstype = cross$crosstype)
}

# Genome-wide thresholds by permutation (Churchill & Doerge 1994), with separate
# autosome and X thresholds when the cross has an X chromosome (Broman et al.
# 2006). Traits are permuted one at a time: qtl2 is much slower when traits
# with different missing values are permuted together. The same seed is used
# for every trait, so each trait sees the same permutation sequence.
run_perm <- function(sc, n_perm = 1000, seed = 1, cores = 1, alpha = c(0.05, 0.10),
                     x_specific = TRUE) {
  has_x <- any(sc$is_x) && x_specific
  per_trait <- lapply(colnames(sc$pheno), function(tr) {
    set.seed(seed)
    op <- qtl2::scan1perm(sc$probs, sc$pheno[, tr, drop = FALSE], kinship = sc$K,
                          addcovar = sc$addcovar, Xcovar = sc$Xcovar, n_perm = n_perm,
                          perm_Xsp = has_x,
                          chr_lengths = if (has_x) qtl2::chr_lengths(sc$map) else NULL,
                          cores = cores)
    list(operm = op, summary = summary(op, alpha = alpha))
  })
  names(per_trait) <- colnames(sc$pheno)
  pick <- function(part) {
    m <- do.call(cbind, lapply(per_trait, function(p)
      if (has_x) unclass(p$summary[[part]])[, 1L] else unclass(p$summary)[, 1L]))
    dimnames(m) <- list(as.character(alpha), colnames(sc$pheno))
    m
  }
  thrA <- pick("A")
  thrX <- if (has_x) pick("X") else NULL
  operm <- lapply(per_trait, `[[`, "operm")
  tab <- do.call(rbind, lapply(colnames(thrA), function(tr)
    data.frame(trait = tr, alpha = alpha, threshold_autosomes = thrA[, tr],
               threshold_X = if (has_x) thrX[, tr] else NA_real_,
               n_perm = n_perm, stringsAsFactors = FALSE)))
  list(operm = operm, table = tab, thrA = thrA, thrX = thrX)
}

# find_peaks with per-trait thresholds; LOD support (drop) or Bayes (prob) intervals.
get_peaks <- function(sc, threshold, thresholdX = NULL, interval = "lod", drop = 1.5,
                      prob = 0.95, peakdrop = Inf) {
  args <- list(scan1_output = sc$out, map = sc$map, threshold = threshold,
               peakdrop = peakdrop, expand2markers = TRUE)
  if (!is.null(thresholdX)) args$thresholdX <- thresholdX
  if (interval == "lod") args$drop <- drop else args$prob <- prob
  pk <- do.call(qtl2::find_peaks, args)
  if (is.null(pk) || !nrow(pk)) return(NULL)
  pk
}

# Proportion of phenotypic variance explained by a single QTL, from its LOD:
# for least squares (Haley-Knott), LOD = (n/2) log10(RSS0/RSS1), so
# R2 = 1 - RSS1/RSS0 = 1 - 10^(-2 LOD / n). Approximate for the LMM.
pve_from_lod <- function(lod, n) 1 - 10^(-2 * lod / n)

effect_contrasts <- function(crosstype, n_geno) {
  if (crosstype == "f2" && n_geno == 3L)
    return(cbind(mu = c(1, 1, 1), a = c(-1, 0, 1), d = c(0, 1, 0)))
  if (crosstype %in% c("riself", "risib", "dh", "haploid") && n_geno == 2L)
    return(cbind(mu = c(1, 1), a = c(-1, 1)))
  if (crosstype == "bc" && n_geno == 2L)
    return(cbind(mu = c(1, 1), AB_minus_AA = c(0, 1)))
  NULL
}

# Effects at one peak with qtl2::fit1 (proper standard errors). Biparental:
# additive a = (BB - AA)/2 (and dominance d = AB - (AA + BB)/2 for F2), plus
# genotype means. Multiparent: founder allele effects centered to sum to zero.
peak_effects <- function(sc, trait, chr, pos) {
  pr_pk <- qtl2::pull_genoprobpos(sc$probs, sc$map, chr, pos)
  y <- sc$pheno[, trait]
  names(y) <- rownames(sc$pheno)
  K <- if (!is.null(sc$K)) sc$K[[chr]] else NULL
  isx <- isTRUE(sc$is_x[[chr]])
  nullcov <- if (isx) sc$Xcovar else NULL
  fit_m <- qtl2::fit1(pr_pk, y, kinship = K, addcovar = sc$addcovar, nullcovar = nullcov,
                      zerosum = sc$multiparent, se = TRUE)
  keep_m <- names(fit_m$coef) %in% c(colnames(pr_pk), "intercept")   # drop covariate terms
  eff <- data.frame(trait = trait, chr = chr, pos = pos,
                    effect = names(fit_m$coef)[keep_m], estimate = unname(fit_m$coef[keep_m]),
                    se = unname(fit_m$SE[keep_m]),
                    type = if (sc$multiparent) "allele_effect" else "genotype_mean",
                    stringsAsFactors = FALSE)
  contr <- if (!sc$multiparent && !isx) effect_contrasts(sc$crosstype, ncol(pr_pk)) else NULL
  if (!is.null(contr)) {
    fit_c <- qtl2::fit1(pr_pk, y, kinship = K, addcovar = sc$addcovar, contrasts = contr, se = TRUE)
    keep <- setdiff(colnames(contr), "mu")
    eff <- rbind(eff, data.frame(trait = trait, chr = chr, pos = pos, effect = keep,
                                 estimate = unname(fit_c$coef[keep]), se = unname(fit_c$SE[keep]),
                                 type = "contrast", stringsAsFactors = FALSE))
  }
  attr(eff, "lod") <- fit_m$lod
  eff
}

annotate_peaks <- function(pk, sc, cross, perm, interval) {
  n_all <- attr(sc$out, "sample_size")
  rows <- lapply(seq_len(nrow(pk)), function(k) {
    p <- pk[k, ]
    tr <- as.character(p$lodcolumn)
    isx <- isTRUE(sc$is_x[[as.character(p$chr)]])
    thr <- if (isx && !is.null(perm$thrX)) perm$thrX[, tr] else perm$thrA[, tr]
    n <- unname(n_all[tr])
    nearest <- qtl2::find_marker(cross$gmap, p$chr, p$pos)
    lo_m <- qtl2::find_marker(cross$gmap, p$chr, p$ci_lo)
    hi_m <- qtl2::find_marker(cross$gmap, p$chr, p$ci_hi)
    out <- data.frame(trait = tr, chr = as.character(p$chr), pos_cM = p$pos, lod = p$lod,
                      ci_lo_cM = p$ci_lo, ci_hi_cM = p$ci_hi,
                      interval = if (interval == "lod") "LOD support" else "Bayes credible",
                      nearest_marker = nearest, ci_lo_marker = lo_m, ci_hi_marker = hi_m,
                      significance = ifelse(p$lod >= thr[["0.05"]], "genome-wide 5%", "suggestive 10%"),
                      threshold_5pct = thr[["0.05"]], threshold_10pct = thr[["0.1"]],
                      n = n, pve_pct = 100 * pve_from_lod(p$lod, n), stringsAsFactors = FALSE)
    if (!is.null(cross$pmap) && as.character(p$chr) %in% names(cross$pmap)) {
      gm <- cross$gmap[[as.character(p$chr)]]
      pm <- cross$pmap[[as.character(p$chr)]][names(gm)]
      # Physical span of every marker inside the genetic interval (robust to local
      # order differences between the genetic and physical maps); the peak itself
      # is interpolated.
      inside <- gm >= p$ci_lo - 1e-8 & gm <= p$ci_hi + 1e-8
      conv <- function(x) stats::approx(gm, pm, xout = x, rule = 2, ties = mean)$y
      out$ci_lo_Mb <- if (any(inside)) min(pm[inside]) else conv(p$ci_lo)
      out$pos_Mb <- conv(p$pos)
      out$ci_hi_Mb <- if (any(inside)) max(pm[inside]) else conv(p$ci_hi)
    }
    out
  })
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

fmt <- function(x, d = 2) ifelse(is.na(x), "", sprintf(paste0("%.", d, "f"), x))

md_table <- function(df) {
  hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  body <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  c(hdr, sep, body)
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (opt$demo) {
    dd <- file.path(SKILL_DIR, "demo")
    opt$geno <- file.path(dd, "geno.csv")
    opt$map <- file.path(dd, "gmap.csv")
    opt$pmap <- file.path(dd, "pmap.csv")
    opt$pheno <- file.path(dd, "pheno.csv")
    opt$cross_type <- "riself"
  }
  csv_mode <- nzchar(opt$geno) || nzchar(opt$map) || nzchar(opt$pheno)
  if (!nzchar(opt$control) && !csv_mode)
    stop("Provide --control, or --geno, --map, --pheno and --cross-type, or --demo. See --help.",
         call. = FALSE)
  if (nzchar(opt$control) && csv_mode)
    stop("Give either --control or the CSV files, not both.", call. = FALSE)
  if (csv_mode && !(nzchar(opt$geno) && nzchar(opt$map) && nzchar(opt$pheno) && nzchar(opt$cross_type)))
    stop("CSV input needs --geno, --map, --pheno and --cross-type.", call. = FALSE)
  if (!opt$kinship %in% c("none", "loco")) stop("--kinship must be none or loco.", call. = FALSE)
  if (!opt$interval %in% c("lod", "bayes")) stop("--interval must be lod or bayes.", call. = FALSE)
  if (!opt$map_function %in% c("haldane", "kosambi"))
    stop("--map-function must be haldane or kosambi.", call. = FALSE)
  if (opt$error_prob <= 0 || opt$error_prob >= 0.5) stop("--error-prob must be in (0, 0.5).", call. = FALSE)
  if (opt$n_perm < 1) stop("--n-perm must be at least 1.", call. = FALSE)
  if (opt$interval == "lod" && opt$lod_drop > opt$peakdrop)
    stop("--lod-drop must not exceed --peakdrop.", call. = FALSE)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  warn <- character()

  # --- data -----------------------------------------------------------------
  input <- if (csv_mode) {
    build_cross_from_csv(opt$geno, opt$map, opt$pheno, opt$cross_type, opt$geno_codes,
                         opt$pmap, opt$covar, opt$id_col, opt$map_units, opt$cm_per_mb,
                         work_dir = file.path(opt$out, "qtl2_input"))
  } else read_cross_control(opt$control)
  cross <- input$cross
  warn <- c(warn, input$notes)
  ok <- suppressWarnings(qtl2::check_cross2(cross))
  if (!isTRUE(ok)) stop("qtl2::check_cross2 found problems in the cross: ",
                        paste(attr(ok, "problems") %||% "see qtl2 warnings", collapse = "; "),
                        call. = FALSE)
  traits <- if (nzchar(opt$traits)) trimws(strsplit(opt$traits, ",", fixed = TRUE)[[1]]) else
    colnames(cross$pheno)
  miss <- setdiff(traits, colnames(cross$pheno))
  if (length(miss)) stop("Trait(s) not found: ", paste(miss, collapse = ", "), ". Available: ",
                         paste(colnames(cross$pheno), collapse = ", "), call. = FALSE)
  multiparent <- !cross$crosstype %in% BIPARENTAL
  addcovar <- covariate_matrix(cross, opt$addcovar)

  qc <- qc_cross(cross)
  save_table(qc$markers, opt$out, "marker_qc.csv")
  save_table(qc$individuals, opt$out, "individual_qc.csv")
  n_dist <- sum(qc$markers$distorted)
  if (n_dist)
    warn <- c(warn, sprintf(paste0("%d marker(s) show segregation distortion (chi-square, Bonferroni 5%%; ",
                                   "see tables/marker_qc.csv). Isolated distorted markers suggest genotyping ",
                                   "errors; clustered ones suggest selection or a wrong cross type."), n_dist))
  hi_m <- sum(qc$markers$missing_rate > 0.2)
  if (hi_m) warn <- c(warn, sprintf("%d marker(s) have more than 20%% missing genotypes.", hi_m))
  hi_i <- sum(qc$individuals$missing_rate > 0.2)
  if (hi_i) warn <- c(warn, sprintf("%d individual(s) have more than 20%% missing genotypes.", hi_i))
  if (!is.null(qc$duplicates)) {
    save_table(qc$duplicates, opt$out, "duplicate_pairs.csv")
    warn <- c(warn, sprintf(paste0("%d pair(s) of individuals share > 95%% of genotypes (possible duplicates ",
                                   "or sample mix-ups; see tables/duplicate_pairs.csv)."), nrow(qc$duplicates)))
  }
  ph_rows <- lapply(traits, function(tr) {
    y <- cross$pheno[, tr]
    data.frame(trait = tr, n = sum(is.finite(y)), mean = mean(y, na.rm = TRUE),
               sd = stats::sd(y, na.rm = TRUE), skewness = skewness(y),
               n_beyond_4sd = sum(abs(y - mean(y, na.rm = TRUE)) > 4 * stats::sd(y, na.rm = TRUE), na.rm = TRUE))
  })
  phs <- do.call(rbind, ph_rows)
  save_table(phs, opt$out, "phenotype_summary.csv")
  for (k in seq_len(nrow(phs))) {
    if (phs$n[k] < 100)
      warn <- c(warn, sprintf(paste0("Trait %s has only %d individuals: power is low and the effects and ",
                                     "variance explained of detected QTL are overestimated (Beavis effect)."),
                              phs$trait[k], phs$n[k]))
    if (is.finite(phs$skewness[k]) && abs(phs$skewness[k]) > 1)
      warn <- c(warn, sprintf("Trait %s is skewed (skewness %.2f): consider a transformation; Haley-Knott assumes roughly normal residuals.",
                              phs$trait[k], phs$skewness[k]))
    if (phs$n_beyond_4sd[k] > 0)
      warn <- c(warn, sprintf("Trait %s has %d value(s) beyond 4 SD: check for data errors.",
                              phs$trait[k], phs$n_beyond_4sd[k]))
  }
  if (opt$n_perm < 1000)
    warn <- c(warn, sprintf("Only %d permutations: thresholds are imprecise; use at least 1000 for reporting.",
                            as.integer(opt$n_perm)))

  # --- scan, thresholds, peaks ----------------------------------------------
  t0 <- Sys.time()
  sc <- run_scan(cross, traits, opt$step, opt$error_prob, opt$map_function, opt$kinship,
                 addcovar, opt$cores)
  perm <- run_perm(sc, opt$n_perm, as.integer(opt$seed), opt$cores)
  save_table(perm$table, opt$out, "thresholds.csv")
  thr10 <- perm$thrA["0.1", traits]
  thr10X <- if (!is.null(perm$thrX)) perm$thrX["0.1", traits] else NULL
  pk <- get_peaks(sc, thr10, thr10X, opt$interval, opt$lod_drop, opt$bayes_prob, opt$peakdrop)
  peaks <- NULL
  effects <- NULL
  if (!is.null(pk)) {
    peaks <- annotate_peaks(pk, sc, cross, perm, opt$interval)
    effects <- do.call(rbind, lapply(seq_len(nrow(peaks)), function(k)
      peak_effects(sc, peaks$trait[k], peaks$chr[k], peaks$pos_cM[k])))
    if (!sc$multiparent) {
      for (nm in c("a", "d", "AB_minus_AA")) {
        sel <- effects$type == "contrast" & effects$effect == nm
        if (any(sel)) {
          key <- paste(effects$trait[sel], effects$chr[sel], effects$pos[sel])
          kk <- match(paste(peaks$trait, peaks$chr, peaks$pos_cM), key)
          peaks[[paste0(nm, "_effect")]] <- effects$estimate[sel][kk]
          peaks[[paste0(nm, "_se")]] <- effects$se[sel][kk]
        }
      }
    }
    save_table(peaks, opt$out, "peaks.csv")
    save_table(effects, opt$out, "peak_effects.csv")
    wide <- peaks$ci_hi_cM - peaks$ci_lo_cM > 20
    if (any(wide))
      warn <- c(warn, sprintf("%d QTL interval(s) are wider than 20 cM: positional resolution is low.", sum(wide)))
    if (any(table(paste(peaks$trait, peaks$chr)) > 1))
      warn <- c(warn, "Several peaks for one trait on one chromosome: a single-QTL scan cannot separate linked QTL; confirm with a multiple-QTL model.")
  }
  lod_tab <- data.frame(position = rownames(sc$out),
                        chr = rep(names(sc$map), lengths(sc$map)),
                        pos_cM = unlist(sc$map, use.names = FALSE), unclass(sc$out),
                        check.names = FALSE, row.names = NULL)
  save_table(lod_tab, opt$out, "lod_scores.csv")
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # --- figures -----------------------------------------------------------------
  for (tr in traits) {
    open_png(opt$out, sprintf("lod_%s.png", gsub("[^A-Za-z0-9_.-]", "_", tr)), width = 9, height = 4.5)
    ymax <- max(c(qtl2::maxlod(qtl2::subset_scan1(sc$out, lodcolumn = tr)), perm$thrA[, tr],
                  if (!is.null(perm$thrX)) perm$thrX[, tr]), na.rm = TRUE)
    graphics::par(mar = c(4.5, 4.5, 2.5, 1))
    qtl2::plot_scan1(sc$out, sc$map, lodcolumn = tr, ylim = c(0, ymax * 1.08),
                     col = "slateblue", main = tr)
    graphics::abline(h = perm$thrA["0.05", tr], col = "red", lty = 2)
    graphics::abline(h = perm$thrA["0.1", tr], col = "darkorange", lty = 3)
    leg <- c("genome-wide 5% (autosomes)", "suggestive 10% (autosomes)")
    if (!is.null(perm$thrX))
      leg <- c(leg, sprintf("X chromosome: 5%% = %.2f, 10%% = %.2f", perm$thrX["0.05", tr],
                            perm$thrX["0.1", tr]))
    graphics::legend("topright", leg, col = c("red", "darkorange", NA)[seq_along(leg)],
                     lty = c(2, 3, NA)[seq_along(leg)], bty = "n", cex = 0.8)
    grDevices::dev.off()
  }
  if (!is.null(peaks)) {
    for (k in seq_len(nrow(peaks))) {
      ch <- peaks$chr[k]
      tr <- peaks$trait[k]
      y <- sc$pheno[, tr]
      K <- if (!is.null(sc$K)) sc$K[[ch]] else NULL
      contr <- if (!sc$multiparent && !isTRUE(sc$is_x[[ch]]))
        effect_contrasts(sc$crosstype, dim(sc$probs[[ch]])[2]) else NULL
      co <- qtl2::scan1coef(sc$probs[, ch], y, kinship = K, addcovar = sc$addcovar,
                            contrasts = contr, zerosum = sc$multiparent)
      cols <- if (is.null(contr)) seq_len(dim(sc$probs[[ch]])[2]) else 2:ncol(contr)
      open_png(opt$out, sprintf("effects_%s_chr%s_%03d.png", gsub("[^A-Za-z0-9_.-]", "_", tr), ch,
                                as.integer(round(peaks$pos_cM[k]))), width = 8, height = 6)
      qtl2::plot_coef(co, sc$map[ch], columns = cols,
                      scan1_output = qtl2::subset_scan1(sc$out, lodcolumn = tr),
                      main = sprintf("%s, chr %s", tr, ch), legend = "topright",
                      ylab = "Effect (trait units)",
                      col = grDevices::hcl.colors(length(cols), "Dark 3"))
      grDevices::dev.off()
    }
  }

  # --- report ---------------------------------------------------------------------
  nmk <- vapply(cross$geno, ncol, 1L)
  len <- vapply(cross$gmap, function(m) diff(range(m)), 0)
  body <- c("## Data", "",
            sprintf("Cross type `%s`%s; %d individuals genotyped; %d markers on %d chromosomes (%.0f cM in total, median spacing %.2f cM).",
                    cross$crosstype, if (multiparent) " (multiparent)" else "",
                    nrow(cross$geno[[1]]), sum(nmk), length(nmk), sum(len),
                    stats::median(unlist(lapply(cross$gmap, diff)))),
            sprintf("Traits: %s.", paste(traits, collapse = ", ")), "",
            md_table(data.frame(trait = phs$trait, n = phs$n, mean = fmt(phs$mean), sd = fmt(phs$sd),
                                skewness = fmt(phs$skewness))), "",
            "## Method", "",
            sprintf("- Genotype probabilities: qtl2 hidden Markov model, genotyping error probability %g, %s map function, pseudomarkers every %g cM%s.",
                    opt$error_prob, opt$map_function, opt$step, if (opt$step == 0) " (none)" else ""),
            sprintf("- Genome scan: %s%s.",
                    if (opt$kinship == "loco") "linear mixed model with leave-one-chromosome-out (LOCO) kinship" else "Haley-Knott regression",
                    if (multiparent) " on founder allele probabilities (additive allele model)" else ""),
            sprintf("- Thresholds: %d permutations (seed %s), genome-wide alpha = 0.05 (significant) and 0.10 (suggestive)%s.",
                    as.integer(opt$n_perm), format(opt$seed, scientific = FALSE),
                    if (!is.null(perm$thrX)) ", separate autosome and X thresholds" else ""),
            sprintf("- Intervals: %s; peaks above the 10%% threshold, %s.",
                    if (opt$interval == "lod") sprintf("%.1f-LOD support intervals expanded to flanking markers", opt$lod_drop)
                    else sprintf("%.0f%% Bayes credible intervals expanded to flanking markers", 100 * opt$bayes_prob),
                    if (is.infinite(opt$peakdrop)) "at most one per chromosome and trait" else sprintf("peakdrop %.1f", opt$peakdrop)),
            "- Variance explained: PVE = 1 - 10^(-2 LOD / n) for each QTL on its own (exact for Haley-Knott, approximate for the LMM).",
            "",
            "## Thresholds", "",
            md_table(data.frame(trait = perm$table$trait, alpha = perm$table$alpha,
                                autosomes = fmt(perm$table$threshold_autosomes),
                                X = fmt(perm$table$threshold_X))), "",
            "## QTL", "")
  if (is.null(peaks)) {
    body <- c(body, "No peak exceeded the suggestive (10%) genome-wide threshold for any trait.", "")
  } else {
    pt <- data.frame(trait = peaks$trait, chr = peaks$chr, pos_cM = fmt(peaks$pos_cM, 1),
                     interval_cM = paste0(fmt(peaks$ci_lo_cM, 1), "-", fmt(peaks$ci_hi_cM, 1)),
                     LOD = fmt(peaks$lod), significance = peaks$significance,
                     PVE_pct = fmt(peaks$pve_pct, 1), marker = peaks$nearest_marker)
    if (!is.null(peaks$pos_Mb))
      pt$interval_Mb <- paste0(fmt(peaks$ci_lo_Mb, 3), "-", fmt(peaks$ci_hi_Mb, 3))
    if (!is.null(peaks$a_effect)) pt$additive <- paste0(fmt(peaks$a_effect), " (", fmt(peaks$a_se), ")")
    if (!is.null(peaks$d_effect)) pt$dominance <- paste0(fmt(peaks$d_effect), " (", fmt(peaks$d_se), ")")
    al <- if (length(cross$alleles) == 2L) cross$alleles else c("A", "B")
    body <- c(body, md_table(pt), "",
              if (!sc$multiparent)
                sprintf(paste0("Additive effect a = (%s%s - %s%s)/2 in trait units (SE): positive means the %s ",
                               "(second parent) allele increases the trait. Dominance d = %s%s - (%s%s + %s%s)/2 ",
                               "(F2 only). Genotype means are in `tables/peak_effects.csv`."),
                        al[2], al[2], al[1], al[1], al[2], al[1], al[2], al[1], al[1], al[2], al[2]) else
                "Founder allele effects (centered to sum to zero) are in `tables/peak_effects.csv`.",
              "")
  }
  body <- c(body, "## Figures", "",
            paste0("![LOD ", traits, "](figures/lod_", gsub("[^A-Za-z0-9_.-]", "_", traits), ".png)"), "",
            sprintf("Run time of scans and permutations: %.0f s.", secs), "")
  limits <- c("## Limitations", "",
              "- Single-QTL genome scans: each position is tested on its own. Linked QTL can merge into one peak or create a ghost peak between them, and effects are not adjusted for other QTL. Confirm multiple QTL with a multiple-QTL model before reporting them as separate loci.",
              "- Effects and PVE of detected QTL are biased upward, strongly so with fewer than a few hundred individuals (Beavis effect; Xu 2003). PVE values of different QTL are not additive.",
              "- Intervals are approximate: LOD support and Bayes intervals assume a single QTL on the chromosome.",
              "")
  refs <- c("## References", "",
            "- Broman KW et al. (2019) Genetics 211:495-502. https://doi.org/10.1534/genetics.118.301595",
            "- Haley CS, Knott SA (1992) Heredity 69:315-324. https://doi.org/10.1038/hdy.1992.131",
            "- Lander ES, Botstein D (1989) Genetics 121:185-199. https://doi.org/10.1093/genetics/121.1.185",
            "- Churchill GA, Doerge RW (1994) Genetics 138:963-971. https://doi.org/10.1093/genetics/138.3.963",
            "- Broman KW et al. (2006) Genetics 174:2151-2158. https://doi.org/10.1534/genetics.106.061176",
            "- Manichaikul A et al. (2006) Genetics 174:481-489. https://doi.org/10.1534/genetics.106.061549",
            "- Xu S (2003) Genetics 165:2259-2268. https://doi.org/10.1093/genetics/165.4.2259",
            "- Broman KW, Sen S (2009) A Guide to QTL Mapping with R/qtl. Springer. https://doi.org/10.1007/978-0-387-92125-9")
  report <- c(sprintf("# %s report", SKILL), "",
              if (length(warn)) c("## Warnings", "", paste("-", warn), "") else NULL,
              body, limits, refs, "", sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(report, file.path(opt$out, "report.md"))
  inputs <- if (csv_mode) c(opt$geno, opt$map, opt$pmap, opt$pheno, opt$covar) else opt$control
  write_repro_bundle(opt$out, argv, inputs = inputs)
  cat(sprintf("Done: %d QTL peak(s). Report: %s\n", if (is.null(peaks)) 0L else nrow(peaks),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(list(cross = cross, scan = sc, perm = perm, peaks = peaks, effects = effects,
                 qc = qc, warnings = warn))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

if (sys.nframe() == 0L) main()
