#!/usr/bin/env Rscript
# claw-gwas-crop: mixed-model genome-wide association study (GWAS) for crop
# panels with statgenGWAS: kinship (Astle-Balding, VanRaden or IBS), optional
# leave-one-chromosome-out (LOCO) kinship, covariates and principal components,
# MAF filtering, Bonferroni or Benjamini-Hochberg FDR thresholds, genomic-control
# lambda, QQ and Manhattan plots, and grouping of significant SNPs into QTL
# regions with an LD distance.
#
# Usage:
#   Rscript scripts/gwas_crop.R --demo
#   Rscript scripts/gwas_crop.R --geno genotypes_filtered_012.csv --pheno blues.csv \
#           --traits yield_kg_ha,height_cm --ld-distance 250000
#   Rscript scripts/gwas_crop.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-gwas-crop"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/gwas_crop.R [--demo] --geno FILE --pheno FILE --traits A,B
                                [options] [--out DIR]

Input
  --demo              run on the bundled synthetic structured panel
  --geno FILE         0/1/2 marker table from claw-marker-qc (marker, chr, pos, then one
                      column per line; .csv/.tsv/.gz) or a VCF (.vcf/.vcf.gz)
  --pheno FILE        CSV with one row per line (BLUEs), an ID column and trait columns
  --id-col NAME       line ID column in --pheno (default: genotype)
  --traits A,B        trait columns to analyse (required)
  --covar A,B         extra covariate columns in --pheno (numeric or categorical)

Model
  --kinship M         astle (default), vanRaden or IBS
  --loco              chromosome-specific kinship (leave one chromosome out)
  --n-pcs K           add the first K marker principal components as covariates (default 0)
  --maf X             test only markers with minor allele frequency >= X (default 0.05)
  --max-missing X     drop markers with more than X missing calls before mean
                      imputation (default 0.2)
  --threshold T       bonferroni (default) or fdr (Benjamini-Hochberg)
  --alpha X           family-wise error (bonferroni) or false discovery rate (fdr) (default 0.05)
  --compare-naive     also fit a model without kinship to show the inflation it causes
  --n-cores N         cores used by statgenGWAS (default 1)

QTL regions
  --ld-distance BP    join significant SNPs closer than BP and extend regions by BP
  --ld-table FILE     per-chromosome distances (chr, ld_distance_bp) from claw-ld-decay

Output
  --out DIR           output folder (default: agriclaw_output/claw-gwas-crop)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, geno = "", pheno = "", id_col = "genotype", traits = "",
                 covar = "", kinship = "astle", loco = FALSE, n_pcs = 0, maf = 0.05,
                 max_missing = 0.2, threshold = "bonferroni", alpha = 0.05,
                 compare_naive = FALSE, n_cores = 1, ld_distance = NA_real_,
                 ld_table = "", out = "")

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

COL_BLUE <- "#2a78d6"; COL_ORANGE <- "#eb6834"; COL_INK <- "#0b0b0b"
COL_INK2 <- "#52514e"; COL_MUTED <- "#898781"; COL_BAND <- "#e1e0d9"; COL_GREY <- "#a3a29b"

open_text <- function(path) {
  if (!nzchar(path) || !file.exists(path)) stop("File not found: ", path, call. = FALSE)
  gzfile(path, open = "rt")
}

first_line <- function(path) {
  con <- open_text(path)
  on.exit(close(con))
  l <- readLines(con, n = 1L, warn = FALSE)
  if (!length(l)) stop("File is empty: ", path, call. = FALSE)
  l
}

ALIASES <- list(marker = c("marker", "snp", "rs#", "marker_id", "snp_id", "id", "name"),
                chr = c("chr", "chrom", "chromosome"),
                pos = c("pos", "position", "bp", "pos_bp", "position_bp"),
                ref = c("ref", "ref_allele"), alt = c("alt", "alt_allele"))

find_col <- function(nm, aliases) {
  i <- match(TRUE, tolower(nm) %in% aliases)
  if (is.na(i)) NA_integer_ else i
}

# Biallelic diploid GT strings to ALT-allele counts (half-missing -> NA).
gt_to_dosage <- function(gt) {
  g <- gsub("|", "/", gt, fixed = TRUE)
  g[g %in% c(".", "./.")] <- NA
  a1 <- sub("/.*$", "", g)
  a2 <- sub("^.*/", "", g)
  called <- !is.na(g)
  if (any(called & (nchar(g) - nchar(gsub("/", "", g, fixed = TRUE))) != 1L))
    stop("Only diploid GT calls are supported. Convert the VCF with claw-marker-qc first.", call. = FALSE)
  ok <- called & a1 %in% c("0", "1") & a2 %in% c("0", "1")
  d <- matrix(NA_real_, nrow(g), ncol(g))
  d[ok] <- (a1[ok] == "1") + (a2[ok] == "1")
  d
}

read_vcf_dosage <- function(path) {
  vcf <- vcfR::read.vcfR(path, verbose = FALSE)
  fix <- vcfR::getFIX(vcf)
  gt <- vcfR::extract.gt(vcf, element = "GT", IDtoRowNames = FALSE)
  multi <- !is.na(fix[, "ALT"]) & grepl(",", fix[, "ALT"], fixed = TRUE)
  fix <- fix[!multi, , drop = FALSE]
  gt <- gt[!multi, , drop = FALSE]
  ids <- fix[, "ID"]
  noid <- is.na(ids) | ids == "."
  ids[noid] <- paste0(fix[noid, "CHROM"], "_", fix[noid, "POS"])
  G <- gt_to_dosage(gt)
  dimnames(G) <- list(make.unique(ids), colnames(vcf@gt)[-1])
  notes <- if (any(multi)) sprintf("%d multi-allelic VCF site(s) excluded.", sum(multi)) else character()
  list(G = G, map = data.frame(marker = rownames(G), chr = fix[, "CHROM"],
                               pos = as.numeric(fix[, "POS"]), stringsAsFactors = FALSE),
       notes = notes)
}

# AgriClaw 0/1/2 marker table (claw-marker-qc output) or VCF. Returns a numeric
# matrix with markers in rows and lines in columns, plus the map.
read_genotypes <- function(path) {
  l1 <- first_line(path)
  if (startsWith(l1, "##fileformat=VCF")) return(read_vcf_dosage(path))
  n_tab <- sum(gregexpr("\t", l1, fixed = TRUE)[[1]] > 0)
  n_com <- sum(gregexpr(",", l1, fixed = TRUE)[[1]] > 0)
  con <- open_text(path)
  df <- tryCatch(utils::read.table(con, header = TRUE, sep = if (n_tab > n_com) "\t" else ",",
                                   colClasses = "character", check.names = FALSE, quote = "\"",
                                   comment.char = "", na.strings = c("", "NA"), stringsAsFactors = FALSE),
                 finally = close(con))
  nm <- names(df)
  im <- find_col(nm, ALIASES$marker); ic <- find_col(nm, ALIASES$chr); ip <- find_col(nm, ALIASES$pos)
  if (is.na(ic) || is.na(ip)) {
    stop("The genotype file needs marker, chr and pos columns followed by one column per line ",
         "(the table written by claw-marker-qc). Convert HapMap files, VCFs with unusual coding, ",
         "or matrices with lines in rows using claw-marker-qc first.", call. = FALSE)
  }
  if (is.na(im)) im <- 1L
  meta <- unique(stats::na.omit(c(im, ic, ip, find_col(nm, ALIASES$ref), find_col(nm, ALIASES$alt))))
  vals <- as.matrix(df[, -meta, drop = FALSE])
  num <- suppressWarnings(matrix(as.numeric(vals), nrow(vals), ncol(vals)))
  bad <- !is.na(vals) & is.na(num)
  if (any(bad)) stop("Genotype values must be numeric ALT-allele counts (0/1/2); found '",
                     vals[bad][1], "'.", call. = FALSE)
  if (any(num < 0 | num > 2, na.rm = TRUE))
    stop("Genotype values must lie between 0 and 2 (ALT-allele counts). Recode -1/0/1 data with claw-marker-qc.", call. = FALSE)
  dimnames(num) <- list(make.unique(df[[im]]), colnames(vals))
  list(G = num, map = data.frame(marker = rownames(num), chr = df[[ic]],
                                 pos = suppressWarnings(as.numeric(df[[ip]])), stringsAsFactors = FALSE),
       notes = character())
}

read_pheno <- function(path, id_col, traits, covar) {
  if (!file.exists(path)) stop("Phenotype file not found: ", path, call. = FALSE)
  ph <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!id_col %in% names(ph))
    stop("ID column '", id_col, "' not found in the phenotype file. Use --id-col. Columns: ",
         paste(names(ph), collapse = ", "), call. = FALSE)
  if (!length(traits)) {
    num <- setdiff(names(ph)[vapply(ph, is.numeric, NA)], c(id_col, covar))
    stop("Name the trait column(s) with --traits. Numeric columns available: ",
         paste(num, collapse = ", "), call. = FALSE)
  }
  miss <- setdiff(c(traits, covar), names(ph))
  if (length(miss)) stop("Column(s) not found in the phenotype file: ", paste(miss, collapse = ", "),
                         ". Available: ", paste(names(ph), collapse = ", "), call. = FALSE)
  for (t in traits) if (!is.numeric(ph[[t]]))
    stop("Trait '", t, "' is not numeric. Check for text such as 'NA*' or units in the cells.", call. = FALSE)
  ph[[id_col]] <- as.character(ph[[id_col]])
  dup <- unique(ph[[id_col]][duplicated(ph[[id_col]])])
  if (length(dup)) stop("The phenotype file has ", length(dup), " repeated line ID(s) (e.g. ", dup[1],
                        "). GWAS needs one value per line: compute BLUEs across plots and ",
                        "environments first (claw-trial-analysis, claw-met-gxe).", call. = FALSE)
  ph
}

# Mean imputation per marker (ALT-allele count = 2p), after dropping markers
# with too much missing data.
impute_mean <- function(G, max_missing = 0.2) {
  miss_frac <- rowMeans(is.na(G))
  keep <- miss_frac <= max_missing
  G <- G[keep, , drop = FALSE]
  n_imp <- sum(is.na(G))
  if (n_imp) {
    mu <- rowMeans(G, na.rm = TRUE)
    idx <- which(is.na(G), arr.ind = TRUE)
    G[idx] <- mu[idx[, 1]]
  }
  list(G = G, keep = keep, n_imputed = n_imp, frac_imputed = n_imp / length(G))
}

# Principal components of the standardised marker matrix (lines x markers).
marker_pcs <- function(X, k = 10) {
  p <- colMeans(X) / 2
  poly <- p > 0 & p < 1
  Z <- sweep(X[, poly, drop = FALSE], 2, 2 * p[poly])
  Z <- sweep(Z, 2, sqrt(2 * p[poly] * (1 - p[poly])), "/")
  k <- min(k, nrow(X) - 1L, ncol(Z))
  pc <- stats::prcomp(Z, center = FALSE, scale. = FALSE, rank. = k)
  ve <- pc$sdev^2 / sum(pc$sdev^2)
  sc <- pc$x[, seq_len(k), drop = FALSE]
  colnames(sc) <- paste0("PC", seq_len(k))
  list(scores = sc, var_explained = ve[seq_len(k)])
}

# Genomic-control lambda from p-values: median chi-square (1 df) / 0.4549
# (Devlin and Roeder 1999).
lambda_gc <- function(p) {
  p <- p[!is.na(p)]
  stats::median(stats::qchisq(p, 1, lower.tail = FALSE)) / stats::qchisq(0.5, 1)
}

# Smallest QTN R2 (fraction of phenotypic variance) detectable with the given
# power at a per-test level alpha_per_test, for a 1-df test with non-centrality
# n * R2 / (1 - R2) (Sham and Purcell 2014). Ignores loss of power from kinship.
min_detectable_r2 <- function(n, alpha_per_test, power = 0.8) {
  crit <- stats::qchisq(alpha_per_test, 1, lower.tail = FALSE)
  f <- function(ncp) stats::pchisq(crit, 1, ncp = ncp, lower.tail = FALSE) - power
  ncp <- stats::uniroot(f, c(1e-6, 1e4), tol = 1e-10)$root
  ncp / (n + ncp)
}

natural_chr_order <- function(chr) {
  u <- unique(as.character(chr))
  num <- suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", u)))
  num[!grepl("[0-9]", u)] <- NA
  u[order(is.na(num), num, u)]
}

# Per-chromosome LD distances: a scalar, or a table with chr and ld_distance_bp
# (an 'all' row is the default for chromosomes not listed).
ld_lookup <- function(chrs, ld_distance = NA_real_, ld_table = "") {
  out <- stats::setNames(rep(NA_real_, length(chrs)), chrs)
  if (nzchar(ld_table)) {
    lt <- utils::read.csv(ld_table, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("chr", "ld_distance_bp") %in% names(lt)))
      stop("--ld-table needs columns chr and ld_distance_bp (tables/ld_window.csv from claw-ld-decay).", call. = FALSE)
    lt$chr <- as.character(lt$chr)
    lt$ld_distance_bp <- suppressWarnings(as.numeric(lt$ld_distance_bp))
    dflt <- lt$ld_distance_bp[tolower(lt$chr) %in% c("all", "genome", "genome-wide")]
    out[] <- if (length(dflt)) dflt[1] else NA_real_
    hit <- intersect(chrs, lt$chr)
    v <- lt$ld_distance_bp[match(hit, lt$chr)]
    out[hit[!is.na(v)]] <- v[!is.na(v)]   # NA rows (not estimable) fall back to the 'all' row
    if (anyNA(out)) stop("--ld-table has no usable distance for chromosome(s) ",
                         paste(names(out)[is.na(out)], collapse = ", "),
                         " (missing or NA, and no usable 'all' row; see the basis column and the claw-ld-decay report). ",
                         "Ask the user for a window and pass --ld-distance instead.", call. = FALSE)
  } else if (!is.na(ld_distance)) {
    if (ld_distance < 0) stop("--ld-distance must be >= 0 (base pairs).", call. = FALSE)
    out[] <- ld_distance
  }
  out
}

# Group significant SNPs (data frame with snp, chr, pos, p_value) into regions:
# consecutive significant SNPs on a chromosome closer than the LD distance are
# joined; each region is extended by the LD distance on both sides (start >= 1).
group_regions <- function(sig, ld) {
  empty <- data.frame(region_id = character(), chr = character(), start_bp = numeric(),
                      end_bp = numeric(), lead_snp = character(), lead_pos_bp = numeric(),
                      lead_p = numeric(), n_significant_snps = integer(), ld_distance_bp = numeric(),
                      stringsAsFactors = FALSE)
  if (!nrow(sig)) return(empty)
  sig <- sig[order(match(sig$chr, natural_chr_order(sig$chr)), sig$pos), , drop = FALSE]
  out <- list()
  for (ch in unique(sig$chr)) {
    s <- sig[sig$chr == ch, , drop = FALSE]
    d <- if (ch %in% names(ld) && !is.na(ld[[ch]])) ld[[ch]] else 0
    grp <- cumsum(c(1, diff(s$pos) > d))
    for (g in unique(grp)) {
      r <- s[grp == g, , drop = FALSE]
      lead <- which.min(r$p_value)
      out[[length(out) + 1L]] <- data.frame(
        region_id = NA_character_, chr = as.character(ch), start_bp = max(1, min(r$pos) - d),
        end_bp = max(r$pos) + d, lead_snp = r$snp[lead], lead_pos_bp = r$pos[lead],
        lead_p = r$p_value[lead], n_significant_snps = nrow(r), ld_distance_bp = d,
        stringsAsFactors = FALSE)
    }
  }
  res <- do.call(rbind, out)
  res$region_id <- sprintf("R%02d", seq_len(nrow(res)))
  res
}

build_gdata <- function(X, map, ph, id_col, traits, covar_df) {
  map_sg <- data.frame(chr = map$chr, pos = map$pos, row.names = map$marker)
  pheno <- data.frame(genotype = ph[[id_col]], ph[, traits, drop = FALSE],
                      check.names = FALSE, stringsAsFactors = FALSE)
  # createGData names the trial with deparse(substitute(pheno)), so it must be
  # called directly with a variable (not through do.call with a data frame).
  if (!is.null(covar_df) && ncol(covar_df) > 0) {
    suppressMessages(suppressWarnings(statgenGWAS::createGData(geno = X, map = map_sg, pheno = pheno,
                                                               covar = covar_df)))
  } else {
    suppressMessages(suppressWarnings(statgenGWAS::createGData(geno = X, map = map_sg, pheno = pheno)))
  }
}

run_statgen <- function(gdat, traits, covar_names, kinship, loco, maf, alpha, n_cores, kin = NULL) {
  covs <- if (length(covar_names)) covar_names else NULL
  suppressMessages(statgenGWAS::runSingleTraitGwas(
    gData = gdat, traits = traits, covar = covs, kin = kin, kinshipMethod = kinship,
    GLSMethod = if (loco) "multi" else "single", useMAF = TRUE, MAF = maf,
    thrType = "bonf", alpha = alpha, nCores = n_cores))
}

# Tidy statgenGWAS results for one trait and apply the chosen threshold.
tidy_results <- function(gw, trait, threshold, alpha) {
  r <- as.data.frame(gw$GWAResult[[1]])
  r <- r[r$trait == trait, , drop = FALSE]
  tested <- !is.na(r$pValue)
  m <- sum(tested)
  q <- rep(NA_real_, nrow(r))
  q[tested] <- stats::p.adjust(r$pValue[tested], method = "BH")
  sig <- if (threshold == "fdr") tested & q <= alpha else tested & r$pValue <= alpha / m
  out <- data.frame(trait = trait, snp = r$snp, chr = as.character(r$chr), pos = r$pos,
                    alt_allele_freq = r$allFreq, maf = pmin(r$allFreq, 1 - r$allFreq),
                    p_value = r$pValue, minus_log10_p = -log10(r$pValue),
                    effect_per_alt_allele = r$effect, effect_se = r$effectSe,
                    r2_lr = r$RLR2, q_value_bh = q, significant = sig,
                    stringsAsFactors = FALSE)
  out[order(match(out$chr, natural_chr_order(out$chr)), out$pos), , drop = FALSE]
}

# --- Plots -------------------------------------------------------------------

plot_manhattan <- function(res, out, file, trait, line_p = NA, line_label = "") {
  ok <- !is.na(res$p_value)
  res <- res[ok, , drop = FALSE]
  chrs <- natural_chr_order(res$chr)
  len <- tapply(res$pos, factor(res$chr, levels = chrs), max)
  offs <- c(0, cumsum(as.numeric(len))[-length(len)])
  names(offs) <- chrs
  x <- res$pos + offs[res$chr]
  y <- res$minus_log10_p
  ci <- match(res$chr, chrs)
  col <- ifelse(ci %% 2 == 1, COL_BLUE, COL_GREY)
  col[res$significant] <- COL_ORANGE
  open_png(out, file, width = 10, height = 4.5)
  op <- graphics::par(mar = c(4.5, 4.5, 3.5, 1), col.axis = COL_INK2, fg = COL_INK2)
  ymax <- max(c(y, if (!is.na(line_p)) -log10(line_p)), na.rm = TRUE) * 1.05
  graphics::plot(x, y, pch = 19, cex = ifelse(res$significant, 0.9, 0.5), col = col, xaxt = "n",
                 xlab = "Chromosome", ylab = expression(-log[10](p)), ylim = c(0, ymax), bty = "l")
  graphics::title(main = paste("Manhattan plot:", trait), adj = 0, col.main = COL_INK)
  graphics::axis(1, at = offs + len / 2, labels = chrs, tick = FALSE, cex.axis = 0.8)
  if (!is.na(line_p)) graphics::abline(h = -log10(line_p), col = COL_INK2, lty = 2)
  graphics::legend("topright", inset = c(0, -0.16), xpd = NA, horiz = TRUE, bty = "n",
                   pch = c(19, NA), lty = c(NA, 2), col = c(COL_ORANGE, COL_INK2),
                   legend = c("significant SNP", line_label), text.col = COL_INK2, cex = 0.85)
  graphics::par(op)
  grDevices::dev.off()
}

plot_qq <- function(p, out, file, trait, lambda) {
  p <- sort(p[!is.na(p)])
  n <- length(p)
  e <- -log10(stats::ppoints(n))
  o <- -log10(p)
  i <- seq_len(n)
  lo <- -log10(stats::qbeta(0.975, i, n - i + 1))
  hi <- -log10(stats::qbeta(0.025, i, n - i + 1))
  open_png(out, file, width = 5.5, height = 5.5)
  op <- graphics::par(mar = c(4.5, 4.5, 3, 1), col.axis = COL_INK2, fg = COL_INK2)
  lim <- c(0, max(c(o, hi), na.rm = TRUE))
  graphics::plot(NA, xlim = c(0, max(e)), ylim = lim, xlab = expression(Expected ~ -log[10](p)),
                 ylab = expression(Observed ~ -log[10](p)), main = paste("QQ plot:", trait),
                 col.main = COL_INK, bty = "l")
  graphics::polygon(c(e, rev(e)), c(lo, rev(hi)), col = COL_BAND, border = NA)
  graphics::abline(0, 1, col = COL_MUTED, lty = 2)
  graphics::points(e, o, pch = 19, cex = 0.5, col = COL_BLUE)
  graphics::legend("topleft", bty = "n", legend = c(sprintf("lambda = %.3f", lambda), "95% band under the null"),
                   pch = c(19, 15), col = c(COL_BLUE, COL_BAND), text.col = COL_INK2, cex = 0.85)
  graphics::par(op)
  grDevices::dev.off()
}

plot_pca <- function(pcs, out) {
  if (ncol(pcs$scores) < 2L) return(invisible(NULL))
  open_png(out, "pca_structure.png", width = 9, height = 4.5)
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), col.axis = COL_INK2, fg = COL_INK2)
  graphics::plot(pcs$scores[, 1], pcs$scores[, 2], pch = 19, cex = 0.6, col = COL_BLUE,
                 xlab = sprintf("PC1 (%.1f%%)", 100 * pcs$var_explained[1]),
                 ylab = sprintf("PC2 (%.1f%%)", 100 * pcs$var_explained[2]),
                 main = "Population structure", col.main = COL_INK, bty = "l")
  k <- length(pcs$var_explained)
  graphics::plot(seq_len(k), 100 * pcs$var_explained, type = "b", pch = 19, col = COL_BLUE,
                 xlab = "Principal component", ylab = "Variance explained (%)", main = "Scree plot",
                 col.main = COL_INK, bty = "l")
  graphics::par(op)
  grDevices::dev.off()
}

# --- Main ---------------------------------------------------------------------

split_list <- function(x) if (nzchar(x)) trimws(strsplit(x, ",", fixed = TRUE)[[1]]) else character()

check_options <- function(opt) {
  if (!opt$kinship %in% c("astle", "vanRaden", "IBS"))
    stop("--kinship must be astle, vanRaden or IBS.", call. = FALSE)
  if (!opt$threshold %in% c("bonferroni", "fdr"))
    stop("--threshold must be bonferroni or fdr.", call. = FALSE)
  if (opt$alpha <= 0 || opt$alpha >= 1) stop("--alpha must be between 0 and 1.", call. = FALSE)
  if (opt$maf < 0 || opt$maf >= 0.5) stop("--maf must be in [0, 0.5).", call. = FALSE)
  if (opt$max_missing < 0 || opt$max_missing >= 1) stop("--max-missing must be in [0, 1).", call. = FALSE)
  if (opt$n_pcs < 0 || opt$n_pcs != round(opt$n_pcs)) stop("--n-pcs must be a whole number >= 0.", call. = FALSE)
  invisible(TRUE)
}

fmt_p <- function(p) ifelse(is.na(p), "NA", formatC(p, format = "e", digits = 2))

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  options(scipen = 100)
  if (opt$demo) {
    opt$geno <- file.path(SKILL_DIR, "demo", "demo_geno_012.csv.gz")
    opt$pheno <- file.path(SKILL_DIR, "demo", "demo_pheno_blues.csv")
    if (!nzchar(opt$traits)) opt$traits <- "height_cm,protein_pct"
    if (is.na(opt$ld_distance) && !nzchar(opt$ld_table)) opt$ld_distance <- 2e6
    opt$compare_naive <- TRUE
  }
  if (!nzchar(opt$geno) || !nzchar(opt$pheno))
    stop("Provide --geno FILE and --pheno FILE, or --demo. See --help.", call. = FALSE)
  check_options(opt)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  traits <- split_list(opt$traits)
  covar <- split_list(opt$covar)
  warn <- character()

  # Data
  gd <- read_genotypes(opt$geno)
  ph <- read_pheno(opt$pheno, opt$id_col, traits, covar)
  ph <- ph[rowSums(!is.na(ph[, traits, drop = FALSE])) > 0, , drop = FALSE]
  no_geno <- setdiff(ph[[opt$id_col]], colnames(gd$G))
  ids <- intersect(ph[[opt$id_col]], colnames(gd$G))
  if (length(no_geno)) warn <- c(warn, sprintf("%d phenotyped line(s) have no genotypes and were excluded (e.g. %s). Check ID spelling.",
                                               length(no_geno), paste(utils::head(no_geno, 3), collapse = ", ")))
  if (length(ids) < 30L) stop("Only ", length(ids), " lines have both genotypes and phenotypes; GWAS needs many more ",
                              "(about 200 or more for useful power). Check that IDs match between files.", call. = FALSE)
  ph <- ph[match(ids, ph[[opt$id_col]]), , drop = FALSE]
  map <- gd$map
  okmap <- !is.na(map$chr) & nzchar(map$chr) & !is.na(map$pos)
  if (any(!okmap)) warn <- c(warn, sprintf("%d marker(s) without chromosome or position were excluded.", sum(!okmap)))
  G <- gd$G[okmap, ids, drop = FALSE]
  map <- map[okmap, , drop = FALSE]
  imp <- impute_mean(G, opt$max_missing)
  if (any(!imp$keep)) warn <- c(warn, sprintf("%d marker(s) with more than %.0f%% missing calls were excluded.",
                                              sum(!imp$keep), 100 * opt$max_missing))
  G <- imp$G
  map <- map[imp$keep, , drop = FALSE]
  if (imp$n_imputed > 0) {
    msg <- sprintf("%.2f%% of genotype calls (%d) were missing and were mean-imputed (2 x ALT frequency).",
                   100 * imp$frac_imputed, imp$n_imputed)
    if (imp$frac_imputed > 0.10) msg <- paste(msg, "That is a lot: impute with a dedicated tool (e.g. Beagle) before GWAS.")
    warn <- c(warn, msg)
  }
  X <- t(G)  # lines x markers for statgenGWAS
  n <- nrow(X)
  if (opt$loco && length(unique(map$chr)) < 2L) stop("--loco needs markers on at least two chromosomes.", call. = FALSE)

  # Structure diagnostics and covariates
  pcs <- marker_pcs(X, k = max(10, opt$n_pcs))
  plot_pca(pcs, opt$out)
  covar_df <- data.frame(row.names = ids)
  covar_names <- character()
  if (opt$n_pcs > 0) {
    covar_df <- cbind(covar_df, pcs$scores[, seq_len(opt$n_pcs), drop = FALSE])
    covar_names <- colnames(pcs$scores)[seq_len(opt$n_pcs)]
  }
  for (cv in covar) {
    v <- ph[[cv]]
    covar_df[[cv]] <- if (is.numeric(v)) v else factor(v)
    covar_names <- c(covar_names, cv)
  }
  gdat <- build_gdata(X, map, ph, opt$id_col, traits, covar_df)
  gw <- run_statgen(gdat, traits, covar_names, opt$kinship, opt$loco, opt$maf, opt$alpha, opt$n_cores)
  gw0 <- NULL
  if (opt$compare_naive) {
    K0 <- diag(n); dimnames(K0) <- list(ids, ids)
    gw0 <- run_statgen(gdat, traits, covar_names, opt$kinship, FALSE, opt$maf, opt$alpha, opt$n_cores, kin = K0)
  }

  ld <- ld_lookup(unique(map$chr), opt$ld_distance, opt$ld_table)
  if (all(is.na(ld))) warn <- c(warn, paste("No LD distance given: QTL regions are not extended beyond the significant SNPs.",
                                            "Run claw-ld-decay and pass --ld-table (or --ld-distance) for candidate windows."))
  diag_rows <- list(); all_sig <- list(); all_reg <- list(); trait_lines <- character()
  for (tr in traits) {
    res <- tidy_results(gw, tr, opt$threshold, opt$alpha)
    m_t <- sum(!is.na(res$p_value))
    n_t <- sum(!is.na(ph[[tr]]))
    lam <- unname(gw$GWASInfo$inflationFactor[[1]][tr])
    lam0 <- if (!is.null(gw0)) unname(gw0$GWASInfo$inflationFactor[[1]][tr]) else NA_real_
    vc <- if (!opt$loco) gw$GWASInfo$varComp[[1]][[tr]] else c(NA, NA)
    sig <- res[res$significant, , drop = FALSE]
    reg <- group_regions(sig[, c("snp", "chr", "pos", "p_value")], ld)
    if (nrow(reg)) {
      reg <- cbind(trait = tr, reg, stringsAsFactors = FALSE)
      reg$lead_effect_per_alt_allele <- res$effect_per_alt_allele[match(reg$lead_snp, res$snp)]
      reg$lead_alt_allele_freq <- res$alt_allele_freq[match(reg$lead_snp, res$snp)]
      all_reg[[tr]] <- reg
    }
    all_sig[[tr]] <- sig
    save_table(res, opt$out, sprintf("gwas_%s.csv", tr))
    if (opt$threshold == "fdr") {
      line_p <- if (nrow(sig)) max(sig$p_value) else opt$alpha / m_t
      line_lab <- if (nrow(sig)) sprintf("BH FDR %.2f", opt$alpha) else sprintf("Bonferroni %.2f (no FDR hits)", opt$alpha)
    } else {
      line_p <- opt$alpha / m_t
      line_lab <- sprintf("Bonferroni %.2f", opt$alpha)
    }
    plot_manhattan(res, opt$out, sprintf("manhattan_%s.png", tr), tr, line_p, line_lab)
    plot_qq(res$p_value, opt$out, sprintf("qq_%s.png", tr), tr, lam)
    r2min <- min_detectable_r2(n_t, opt$alpha / m_t)
    diag_rows[[tr]] <- data.frame(trait = tr, n_lines = n_t, n_markers_tested = m_t,
                                  lambda_gc = lam, lambda_naive = lam0, vg = vc[1], ve = vc[2],
                                  threshold = opt$threshold, alpha = opt$alpha,
                                  bonferroni_p = opt$alpha / m_t, n_significant_snps = nrow(sig),
                                  n_regions = nrow(reg), min_r2_80pct_power = r2min,
                                  stringsAsFactors = FALSE)
    if (lam > 1.1) warn <- c(warn, sprintf(
      "%s: lambda = %.2f (> 1.1): p-values look inflated. Investigate before interpreting hits: add principal components (--n-pcs 2-5)%s, check for outlier lines and duplicated genotypes.",
      tr, lam, if (opt$loco) " (LOCO without PCs can be inflated in strongly structured panels)" else ""))
    if (lam < 0.9) warn <- c(warn, sprintf(
      "%s: lambda = %.2f (< 0.9): p-values are deflated (conservative); check marker QC, the MAF filter and the sample size.", tr, lam))
    if (n_t < 200) warn <- c(warn, sprintf(
      "%s: only %d lines. Power is low: at this threshold a QTN must explain about %.0f%% of the phenotypic variance to be found with 80%% power.",
      tr, n_t, 100 * r2min))
    z <- (ph[[tr]] - mean(ph[[tr]], na.rm = TRUE)) / stats::sd(ph[[tr]], na.rm = TRUE)
    if (any(abs(z) > 4, na.rm = TRUE)) warn <- c(warn, sprintf(
      "%s: %d line(s) are more than 4 SD from the mean (%s). Check for data-entry or unit errors.",
      tr, sum(abs(z) > 4, na.rm = TRUE), paste(utils::head(ph[[opt$id_col]][which(abs(z) > 4)], 5), collapse = ", ")))
    trait_lines <- c(trait_lines, sprintf("| %s | %d | %d | %.3f | %s | %d | %d | %s |", tr, n_t, m_t, lam,
                                          if (is.na(lam0)) "-" else sprintf("%.2f", lam0), nrow(sig), nrow(reg),
                                          sprintf("%.1f%%", 100 * r2min)))
  }
  diagn <- do.call(rbind, diag_rows)
  sig_all <- do.call(rbind, all_sig)
  reg_all <- if (length(all_reg)) do.call(rbind, all_reg) else
    data.frame(trait = character(), group_regions(sig_all[0, c("snp", "chr", "pos", "p_value")], ld))
  save_table(diagn, opt$out, "diagnostics.csv")
  save_table(sig_all, opt$out, "significant_snps.csv")
  save_table(reg_all, opt$out, "qtl_regions.csv")
  save_table(data.frame(line = ids, pcs$scores, check.names = FALSE), opt$out, "principal_components.csv")

  model_desc <- sprintf("Mixed model y = covariates + marker + g + e, with polygenic effect g ~ N(0, Vg K); kinship %s (%s), %s; %s; MAF >= %.2f; variance components by EMMA REML; marker effects by generalised least squares (statgenGWAS %s).",
                        opt$kinship, if (opt$loco) "chromosome-specific, leave one chromosome out" else "one genome-wide matrix",
                        if (opt$n_pcs > 0) sprintf("%d PCs as fixed covariates", opt$n_pcs) else "no PCs",
                        if (length(covar)) paste("covariates:", paste(covar, collapse = ", ")) else "no other covariates",
                        opt$maf, as.character(utils::packageVersion("statgenGWAS")))
  reg_lines <- if (nrow(reg_all)) c("| trait | region | chr | start (bp) | end (bp) | lead SNP | lead p | SNPs |",
                                    "|---|---|---|---|---|---|---|---|",
                                    sprintf("| %s | %s | %s | %.0f | %.0f | %s | %s | %d |", reg_all$trait, reg_all$region_id,
                                            reg_all$chr, reg_all$start_bp, reg_all$end_bp, reg_all$lead_snp,
                                            fmt_p(reg_all$lead_p), reg_all$n_significant_snps)) else "No significant SNPs."
  report <- c(
    sprintf("# %s report", SKILL), "",
    if (length(warn)) c("## Warnings", "", paste("-", warn), "") else NULL,
    "## Summary", "",
    sprintf("- Genotypes: `%s`; phenotypes: `%s`", opt$geno, opt$pheno),
    sprintf("- Lines analysed: %d; markers after missing-data filter: %d", n, nrow(G)),
    sprintf("- Model: %s", model_desc),
    sprintf("- Threshold: %s", if (opt$threshold == "fdr") sprintf("Benjamini-Hochberg FDR <= %.2f over tested markers", opt$alpha)
            else sprintf("Bonferroni, family-wise alpha %.2f over tested markers", opt$alpha)),
    "- Effects are per copy of the ALT allele, in trait units (half the homozygote difference in inbred lines).", "",
    "| trait | lines | markers tested | lambda | lambda without kinship | significant SNPs | regions | QTN R2 for 80% power |",
    "|---|---|---|---|---|---|---|---|", trait_lines, "",
    "lambda is the genomic-control inflation factor (Devlin and Roeder 1999) from statgenGWAS. Values near 1 mean the test statistics are calibrated; 'lambda without kinship' is an ordinary regression on the same data, shown to make the effect of the kinship correction visible.", "",
    "## QTL regions", "",
    if (!all(is.na(ld))) sprintf("Significant SNPs closer than the LD distance were joined, and regions were extended by it (%s).",
                                 if (nzchar(opt$ld_table)) paste("per chromosome from", opt$ld_table) else sprintf("%.0f bp", opt$ld_distance)) else NULL,
    "", reg_lines, "",
    "## Figures", "",
    paste0("![Manhattan ", traits, "](figures/manhattan_", traits, ".png)"),
    paste0("![QQ ", traits, "](figures/qq_", traits, ".png)"),
    "![Structure](figures/pca_structure.png)", "",
    "## Limitations", "",
    "- Single-locus tests: each marker is tested on its own; large-effect loci can mask nearby signals.",
    "- Missing genotypes are mean-imputed; heavy missingness lowers power.",
    "- Associations are statistical: a lead SNP marks a region in LD with a causal variant, not the causal gene itself.",
    "- The power figure ignores the loss of power when a QTN is correlated with population structure.", "",
    "## References", "",
    "- Kang HM et al. (2008) Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723. https://doi.org/10.1534/genetics.107.080101",
    "- Astle W, Balding DJ (2009) Population structure and cryptic relatedness in genetic association studies. Stat Sci 24:451-471. https://doi.org/10.1214/09-STS307",
    "- VanRaden PM (2008) Efficient methods to compute genomic predictions. J Dairy Sci 91:4414-4423. https://doi.org/10.3168/jds.2007-0980",
    "- Rincent R et al. (2014) Recovering power in association mapping panels with variable levels of linkage disequilibrium. Genetics 197:375-387. https://doi.org/10.1534/genetics.113.159731",
    "- Devlin B, Roeder K (1999) Genomic control for association studies. Biometrics 55:997-1004. https://doi.org/10.1111/j.0006-341X.1999.00997.x",
    "- Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate. J R Stat Soc B 57:289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x", "",
    sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(report, file.path(opt$out, "report.md"))
  write_repro_bundle(opt$out, argv, inputs = c(opt$geno, opt$pheno, opt$ld_table))
  cat(sprintf("Done: %d trait(s), %d lines, %d markers; %d QTL region(s). Report: %s\n", length(traits), n,
              nrow(G), nrow(reg_all), normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(list(diagnostics = diagn, regions = reg_all, significant = sig_all))
}

if (sys.nframe() == 0L) main()
