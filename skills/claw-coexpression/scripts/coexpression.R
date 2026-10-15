#!/usr/bin/env Rscript
# claw-coexpression: weighted gene co-expression network analysis (WGCNA) of
# variance-stabilised RNA-seq data. Filters low-expression and low-variance genes,
# flags outlier samples, optionally removes batch effects (ComBat), picks the
# soft-threshold power for a signed network, detects modules (blockwiseModules),
# computes module eigengenes, module-trait correlations with BH correction, hub
# genes by module membership (kME), and exports module edges for Cytoscape.
#
# Usage:
#   Rscript scripts/coexpression.R --demo
#   Rscript scripts/coexpression.R --expr vst.csv --traits traits.csv
#   Rscript scripts/coexpression.R --counts counts.csv --traits traits.csv --batch batch
#   Rscript scripts/coexpression.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-coexpression"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/coexpression.R [--demo] (--expr VST.csv | --counts COUNTS.csv)
         [--traits TRAITS.csv] [options]

  --demo              run on the bundled simulated expression data (demo/)
  --expr              genes x samples CSV of variance-stabilised (DESeq2 VST) or log-CPM values;
                      first column gene_id (e.g. claw-de-analysis tables/vst.csv)
  --counts            genes x samples CSV of raw counts: filtered and VST-transformed here
  --traits            CSV with a 'sample' column and trait columns (numeric; text columns with
                      2 levels become 0/1, more levels become one 0/1 column per level)
  --batch             column of --traits holding a batch to remove with ComBat (sva)
  --protect           comma-separated --traits columns to protect during batch removal
  --n-top-var         keep this many most variable genes (default: 5000; 0 = all)
  --min-count         --counts only: keep genes with at least this count... (default: 10)
  --min-fraction      ...in at least this fraction of samples (default: 0.1)
  --min-mean          --expr only: drop genes with mean expression below this (default: none)
  --power             soft-threshold power (default: 0 = choose by scale-free fit)
  --r2-cut            scale-free topology fit target, signed R^2 (default: 0.8)
  --cor-type          bicor or pearson (default: bicor, with --max-p-outliers 0.05)
  --max-p-outliers    bicor: maximum share of samples treated as outliers (default: 0.05)
  --min-module-size   smallest module (default: 30)
  --merge-cut-height  merge modules whose eigengenes correlate above 1 - this (default: 0.25)
  --deep-split        tree-cut sensitivity 0-4 (default: 2)
  --max-block-size    genes per network block (default: 5000)
  --outlier-z         flag samples with standardised connectivity below this (default: -2.5)
  --remove-outliers   remove flagged samples before building the network (reported)
  --allow-small-n     run with fewer than 15 samples (not recommended; results are noise-prone)
  --hub-n             hub genes listed per module (default: 10)
  --hub-kme           kME a gene needs to be called a hub (default: 0.8)
  --edge-threshold    smallest topological overlap exported as an edge (default: 0.1)
  --max-edges         strongest edges kept per module for Cytoscape (default: 2000)
  --alpha             FDR for module-trait associations (default: 0.05)
  --seed              random seed (default: 12345)
  --out               output folder (default: agriclaw_output/claw-coexpression)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, expr = "", counts = "", traits = "", batch = "", protect = "",
                 n_top_var = 5000, min_count = 10, min_fraction = 0.1, min_mean = -Inf,
                 power = 0, r2_cut = 0.8, cor_type = "bicor", max_p_outliers = 0.05,
                 min_module_size = 30, merge_cut_height = 0.25, deep_split = 2,
                 max_block_size = 5000, outlier_z = -2.5, remove_outliers = FALSE,
                 allow_small_n = FALSE, hub_n = 10, hub_kme = 0.8, edge_threshold = 0.1,
                 max_edges = 2000, alpha = 0.05, seed = 12345, out = "")

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

read_gene_matrix <- function(path, what) {
  if (!file.exists(path)) fail(what, " file not found: ", path)
  d <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(d) < 3L) fail(what, " needs a gene_id column plus at least two sample columns.")
  ids <- as.character(d[[1L]])
  if (anyDuplicated(ids)) fail("Duplicated gene IDs in ", path)
  m <- as.matrix(d[, -1L, drop = FALSE])
  if (!is.numeric(m)) fail(what, ": sample columns must be numeric.")
  rownames(m) <- ids
  m
}

# Low-count filter (WGCNA FAQ) and DESeq2 VST for raw counts. Returns genes x samples.
vst_from_counts <- function(counts, min_count = 10, min_fraction = 0.1) {
  if (anyNA(counts) || any(counts < 0)) fail("--counts must be non-negative counts without NA.")
  if (any(abs(counts - round(counts)) > 1e-8))
    fail("--counts must be raw integer counts; for VST or log-CPM values use --expr.")
  need <- ceiling(min_fraction * ncol(counts))
  keep <- rowSums(counts >= min_count) >= need
  if (sum(keep) < 50) fail("Only ", sum(keep), " genes have >= ", min_count, " counts in >= ", need, " samples.")
  m <- round(counts[keep, , drop = FALSE])
  storage.mode(m) <- "integer"
  cd <- data.frame(sample = colnames(m), row.names = colnames(m))
  dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(m, cd, ~ 1))
  dds <- suppressMessages(DESeq2::estimateSizeFactors(dds))
  n_ok <- sum(rowMeans(DESeq2::counts(dds, normalized = TRUE)) > 5)
  tr <- if (n_ok >= 1000) suppressMessages(DESeq2::vst(dds, blind = TRUE)) else
    suppressMessages(DESeq2::varianceStabilizingTransformation(dds, blind = TRUE))
  list(expr = SummarizedExperiment::assay(tr), n_in = nrow(counts), n_kept = sum(keep), need = need)
}

check_log_scale <- function(m) {
  if (anyNA(m)) fail("--expr has missing values; WGCNA here needs complete data.")
  if (min(m) >= 0 && max(m) > 100)
    fail("--expr values reach ", signif(max(m), 3), ": they look like counts or normalised counts, ",
         "not a log-scale (VST/log-CPM) matrix. Use --counts for raw counts, or give VST values ",
         "(claw-de-analysis tables/vst.csv).")
  invisible(TRUE)
}

# Traits: numeric columns kept; text with 2 levels -> 0/1; more levels -> one-hot.
read_traits <- function(path, samples, exclude = character()) {
  if (!nzchar(path)) return(list(traits = NULL, raw = NULL, notes = character()))
  if (!file.exists(path)) fail("Traits file not found: ", path)
  d <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
  if (!"sample" %in% names(d)) fail("The traits file needs a 'sample' column.")
  if (anyDuplicated(d$sample)) fail("Duplicated samples in the traits file.")
  miss <- setdiff(samples, d$sample)
  if (length(miss)) fail("Samples with expression data but no row in the traits file: ", paste(miss, collapse = ", "))
  d <- d[match(samples, d$sample), , drop = FALSE]
  rownames(d) <- d$sample
  notes <- character()
  out <- list()
  for (v in setdiff(names(d), c("sample", exclude))) {
    x <- d[[v]]
    if (is.numeric(x)) {
      out[[v]] <- x
    } else {
      lv <- sort(unique(stats::na.omit(as.character(x))))
      if (length(lv) < 2) next
      if (length(lv) == 2) {
        out[[paste0(v, "_", lv[2])]] <- ifelse(is.na(x), NA, as.numeric(x == lv[2]))
        notes <- c(notes, sprintf("Trait '%s' coded 0 = %s, 1 = %s.", v, lv[1], lv[2]))
      } else {
        for (l in lv) out[[paste0(v, "_", l)]] <- ifelse(is.na(x), NA, as.numeric(x == l))
        notes <- c(notes, sprintf("Trait '%s' (%d levels) coded as one 0/1 column per level.", v, length(lv)))
      }
    }
  }
  tr <- if (length(out)) as.data.frame(out, check.names = FALSE, row.names = samples) else NULL
  list(traits = tr, raw = d, notes = notes)
}

# ComBat (sva) on a samples x genes matrix; 'protect' holds covariates whose effects
# are kept. Stops when the batch is completely confounded with a protected factor.
remove_batch <- function(dat, batch, protect = NULL, batch_name = "batch") {
  b <- factor(batch)
  if (nlevels(b) < 2) fail("The batch column has a single level; nothing to remove.")
  if (any(table(b) < 2)) fail("Every batch needs at least 2 samples for ComBat.")
  mod <- NULL
  if (!is.null(protect) && ncol(protect)) {
    for (p in names(protect)) {
      if (!is.numeric(protect[[p]])) {
        tab <- table(b, protect[[p]])
        if (all(rowSums(tab > 0) == 1) && all(colSums(tab > 0) == 1))
          fail("Batch '", batch_name, "' is completely confounded with '", p, "': removing the batch would ",
               "remove the biology. Analyse without batch removal and report the confounding.")
      }
    }
    mod <- stats::model.matrix(stats::as.formula(paste("~", paste(names(protect), collapse = " + "))), data = protect)
  }
  t(suppressMessages(sva::ComBat(dat = t(dat), batch = b, mod = mod)))
}

# Standardised sample connectivity (Oldham et al. 2012): signed adjacency
# ((1 + r) / 2)^2 between samples, k = row sums, Z.k = (k - mean) / sd.
sample_zk <- function(dat_expr) {
  a <- ((1 + stats::cor(t(dat_expr))) / 2)^2
  k <- rowSums(a) - 1
  (k - mean(k)) / stats::sd(k)
}

# Lowest power whose signed scale-free fit reaches r2_cut; otherwise the WGCNA FAQ
# table for signed networks (by sample number).
select_power <- function(fit, n_samples, r2_cut = 0.8) {
  signed_r2 <- -sign(fit$slope) * fit$SFT.R.sq
  ok <- which(signed_r2 >= r2_cut)
  if (length(ok)) return(list(power = fit$Power[min(ok)], fallback = FALSE, signed_r2 = signed_r2))
  pw <- if (n_samples < 20) 18 else if (n_samples <= 30) 16 else if (n_samples <= 40) 14 else 12
  list(power = pw, fallback = TRUE, signed_r2 = signed_r2)
}

# Adjusted Rand index (Hubert & Arabie 1985), for comparing partitions in tests.
adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  c2 <- function(x) x * (x - 1) / 2
  s_ij <- sum(c2(tab))
  s_a <- sum(c2(rowSums(tab)))
  s_b <- sum(c2(colSums(tab)))
  e <- s_a * s_b / c2(length(a))
  mx <- (s_a + s_b) / 2
  if (mx == e) return(1)
  (s_ij - e) / (mx - e)
}

module_trait <- function(MEs, traits, alpha = 0.05) {
  MEs <- MEs[, colnames(MEs) != "ME0", drop = FALSE]
  if (is.null(traits) || !ncol(MEs)) return(NULL)
  r <- stats::cor(MEs, traits, use = "pairwise.complete.obs")
  ok <- crossprod(!is.na(as.matrix(MEs)), !is.na(as.matrix(traits)))
  p <- WGCNA::corPvalueStudent(r, ok)
  out <- data.frame(module = rep(rownames(r), ncol(r)), trait = rep(colnames(r), each = nrow(r)),
                    r = as.vector(r), n = as.vector(ok), pvalue = as.vector(p), stringsAsFactors = FALSE)
  out$padj <- stats::p.adjust(out$pvalue, method = "BH")
  out$significant <- out$padj < alpha
  out[order(out$pvalue), , drop = FALSE]
}

# Smallest |r| reaching p < alpha (two-sided, unadjusted) with n samples.
min_detectable_r <- function(n, alpha = 0.05) {
  t <- stats::qt(1 - alpha / 2, n - 2)
  t / sqrt(n - 2 + t^2)
}

# Element of a 'dist' object of size n for pairs (i, j), i < j.
dist_value <- function(d, n, i, j) {
  lo <- pmin(i, j)
  hi <- pmax(i, j)
  d[n * (lo - 1) - lo * (lo - 1) / 2 + hi - lo]
}

# Within-module TOM edges from the TOM files saved by blockwiseModules.
module_edges <- function(net, gene_ids, labels, threshold = 0.1, max_edges = 2000) {
  edges <- list()
  for (b in seq_along(net$TOMFiles)) {
    env <- new.env()
    load(net$TOMFiles[[b]], envir = env)
    tom <- env$TOM
    bg <- net$blockGenes[[b]]
    nb <- length(bg)
    for (m in setdiff(unique(labels[bg]), 0)) {
      li <- which(labels[bg] == m)
      if (length(li) < 2) next
      pr <- utils::combn(li, 2)
      w <- dist_value(tom, nb, pr[1, ], pr[2, ])
      keep <- w >= threshold
      if (!any(keep)) next
      e <- data.frame(fromNode = gene_ids[bg[pr[1, keep]]], toNode = gene_ids[bg[pr[2, keep]]],
                      weight = w[keep], module = m, stringsAsFactors = FALSE)
      e <- e[order(-e$weight), , drop = FALSE]
      edges[[length(edges) + 1L]] <- utils::head(e, max_edges)
    }
  }
  if (!length(edges)) return(data.frame(fromNode = character(), toNode = character(), weight = numeric(),
                                        module = integer()))
  do.call(rbind, edges)
}

# Core pipeline on a samples x genes matrix. Returns a list of results.
run_wgcna <- function(dat_expr, traits = NULL, power = 0, r2_cut = 0.8, cor_type = "bicor",
                      max_p_outliers = 0.05, min_module_size = 30, merge_cut_height = 0.25,
                      deep_split = 2, max_block_size = 5000, seed = 12345, alpha = 0.05,
                      tom_dir = tempfile("claw_tom_")) {
  # blockwiseModules looks correlation functions up by name, so WGCNA must be attached.
  if (!"package:WGCNA" %in% search()) suppressPackageStartupMessages(library(WGCNA))
  n <- nrow(dat_expr)
  powers <- c(1:10, seq(12, 30, by = 2))
  cor_fnc <- if (cor_type == "bicor") WGCNA::bicor else stats::cor
  cor_opt <- if (cor_type == "bicor") list(use = "p", maxPOutliers = max_p_outliers) else list(use = "p")
  utils::capture.output(sft <- suppressWarnings(WGCNA::pickSoftThreshold(
    dat_expr, powerVector = powers, networkType = "signed", corFnc = cor_fnc, corOptions = cor_opt,
    RsquaredCut = r2_cut, verbose = 0)))
  fit <- sft$fitIndices
  sel <- select_power(fit, n, r2_cut)
  pw <- if (power > 0) power else sel$power
  dir.create(tom_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)
  net <- WGCNA::blockwiseModules(dat_expr, power = pw, networkType = "signed", TOMType = "signed",
                                 corType = cor_type, maxPOutliers = max_p_outliers,
                                 minModuleSize = min_module_size, mergeCutHeight = merge_cut_height,
                                 deepSplit = deep_split, maxBlockSize = max_block_size,
                                 numericLabels = TRUE, pamRespectsDendro = FALSE, reassignThreshold = 0,
                                 randomSeed = seed, saveTOMs = TRUE,
                                 saveTOMFileBase = file.path(tom_dir, "TOM"), verbose = 0)
  labels <- as.integer(net$colors)
  MEs <- WGCNA::orderMEs(net$MEs)
  kme <- WGCNA::signedKME(dat_expr, MEs[, colnames(MEs) != "ME0", drop = FALSE], outputColumnName = "kME")
  mt <- module_trait(MEs, traits, alpha)
  list(sft = fit, power = pw, power_auto = sel, net = net, labels = labels, MEs = MEs, kme = kme,
       module_trait = mt, n = n)
}

hub_table <- function(labels, kme, gene_ids, hub_n = 10, hub_kme = 0.8) {
  rows <- list()
  for (m in sort(setdiff(unique(labels), 0))) {
    col <- paste0("kME", m)
    idx <- which(labels == m)
    v <- kme[idx, col]
    o <- order(-v)
    top <- utils::head(o, hub_n)
    rows[[length(rows) + 1L]] <- data.frame(module = m, module_color = WGCNA::labels2colors(m),
                                            gene_id = gene_ids[idx[top]], kME = round(v[top], 4),
                                            rank = seq_along(top), hub = v[top] >= hub_kme,
                                            stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# --- Figures ------------------------------------------------------------------

plot_samples <- function(dat_expr, zk, z_cut, out) {
  open_png(out, "sample_outliers.png", width = 9, height = 5)
  op <- graphics::par(mfrow = c(1, 2), mar = c(6, 4.5, 3, 1))
  on.exit({graphics::par(op); grDevices::dev.off()})
  hc <- stats::hclust(stats::dist(dat_expr), method = "average")
  graphics::plot(hc, main = "Sample clustering (Euclidean, average)", xlab = "", sub = "", cex = 0.6)
  col <- ifelse(zk < z_cut, "#C0392B", "#7F8C8D")
  graphics::barplot(zk, names.arg = rownames(dat_expr), las = 2, col = col, cex.names = 0.6,
                    ylab = "Standardised connectivity Z.k", main = "Sample connectivity")
  graphics::abline(h = z_cut, lty = 2)
}

plot_soft_threshold <- function(fit, sel, power, r2_cut, out) {
  open_png(out, "soft_threshold.png", width = 9, height = 4.5)
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  on.exit({graphics::par(op); grDevices::dev.off()})
  graphics::plot(fit$Power, sel$signed_r2, type = "n", xlab = "Soft-threshold power",
                 ylab = "Scale-free fit, signed R^2", main = "Scale independence", ylim = c(min(0, sel$signed_r2), 1))
  graphics::text(fit$Power, sel$signed_r2, fit$Power, col = ifelse(fit$Power == power, "#C0392B", "black"), cex = 0.8)
  graphics::abline(h = r2_cut, lty = 2, col = "#C0392B")
  graphics::plot(fit$Power, fit$mean.k., type = "n", xlab = "Soft-threshold power", ylab = "Mean connectivity",
                 main = "Mean connectivity", log = "y")
  graphics::text(fit$Power, fit$mean.k., fit$Power, col = ifelse(fit$Power == power, "#C0392B", "black"), cex = 0.8)
}

plot_dendro <- function(net, out) {
  open_png(out, "gene_dendrogram.png", width = 9, height = 5)
  on.exit(grDevices::dev.off())
  bg <- net$blockGenes[[1]]
  WGCNA::plotDendroAndColors(net$dendrograms[[1]], WGCNA::labels2colors(net$colors[bg]), "Module",
                             dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                             main = "Gene dendrogram and modules (block 1)")
}

plot_module_trait <- function(mt, MEs, out, alpha) {
  mods <- setdiff(colnames(MEs), "ME0")
  trs <- unique(mt$trait)
  r <- matrix(mt$r[match(paste(rep(mods, length(trs)), rep(trs, each = length(mods))), paste(mt$module, mt$trait))],
              length(mods), length(trs), dimnames = list(mods, trs))
  pa <- matrix(mt$padj[match(paste(rep(mods, length(trs)), rep(trs, each = length(mods))), paste(mt$module, mt$trait))],
               length(mods), length(trs))
  txt <- matrix(sprintf("%.2f\n(%s)", r, ifelse(pa < 0.001, "<0.001", sprintf("%.3f", pa))), nrow(r))
  open_png(out, "module_trait.png", width = 3 + 1.3 * length(trs), height = 3.5 + 0.45 * length(mods))
  op <- graphics::par(mar = c(12, 10, 3, 2))
  on.exit({graphics::par(op); grDevices::dev.off()})
  ylab <- paste0(mods, " (", WGCNA::labels2colors(as.integer(sub("ME", "", mods))), ")")
  WGCNA::labeledHeatmap(Matrix = r, xLabels = trs, yLabels = ylab, colors = WGCNA::blueWhiteRed(50),
                        textMatrix = txt, setStdMargins = FALSE, cex.text = 0.6, zlim = c(-1, 1),
                        main = "Module-trait correlation r (BH-adjusted p)")
}

plot_kme_gs <- function(dat_expr, labels, kme, traits, pair, out) {
  m <- as.integer(sub("ME", "", pair$module))
  idx <- which(labels == m)
  gs <- stats::cor(dat_expr[, idx, drop = FALSE], traits[[pair$trait]], use = "p")
  open_png(out, "kme_vs_gene_significance.png", width = 6, height = 5)
  op <- graphics::par(mar = c(4.5, 4.5, 3, 1))
  on.exit({graphics::par(op); grDevices::dev.off()})
  graphics::plot(kme[idx, paste0("kME", m)], abs(gs), pch = 19, cex = 0.6, col = WGCNA::labels2colors(m),
                 xlab = sprintf("Module membership kME (%s)", pair$module),
                 ylab = sprintf("|Gene significance| for %s", pair$trait),
                 main = sprintf("%s vs %s: cor = %.2f", pair$module, pair$trait,
                                stats::cor(kme[idx, paste0("kME", m)], abs(gs), use = "p")))
}

# --- Main -------------------------------------------------------------------

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (opt$demo) {
    opt$expr <- file.path(SKILL_DIR, "demo", "expr_vst.csv")
    opt$traits <- file.path(SKILL_DIR, "demo", "traits.csv")
  }
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  if (!opt$cor_type %in% c("bicor", "pearson")) fail("--cor-type must be bicor or pearson.")
  if (nzchar(opt$expr) == nzchar(opt$counts)) fail("Give exactly one of --expr VST.csv or --counts COUNTS.csv (or --demo).")
  warn <- character()
  notes <- character()
  if (nzchar(opt$counts)) {
    raw <- read_gene_matrix(opt$counts, "--counts")
    v <- vst_from_counts(raw, opt$min_count, opt$min_fraction)
    expr <- v$expr
    notes <- c(notes, sprintf("Counts: kept %d of %d genes with >= %g counts in >= %d samples (WGCNA FAQ rule), then DESeq2 VST (blind).",
                              v$n_kept, v$n_in, opt$min_count, v$need))
  } else {
    expr <- read_gene_matrix(opt$expr, "--expr")
    check_log_scale(expr)
    if (is.finite(opt$min_mean)) {
      keep <- rowMeans(expr) >= opt$min_mean
      notes <- c(notes, sprintf("Mean-expression filter (>= %g): kept %d of %d genes.", opt$min_mean, sum(keep), nrow(expr)))
      expr <- expr[keep, , drop = FALSE]
    }
  }
  n_samples <- ncol(expr)
  if (n_samples < 15 && !opt$allow_small_n)
    fail("Only ", n_samples, " samples. The WGCNA authors do not recommend WGCNA with fewer than 15 samples ",
         "(correlations are too noisy) and suggest at least 20. Use claw-de-analysis for group comparisons, ",
         "or re-run with --allow-small-n and treat modules as exploratory.")
  if (n_samples < 15) warn <- c(warn, sprintf("Only %d samples (--allow-small-n): modules and module-trait correlations are unreliable (WGCNA FAQ: at least 15, preferably 20).", n_samples))
  else if (n_samples < 20) warn <- c(warn, sprintf("%d samples: below the 20 the WGCNA FAQ recommends; interpret modules cautiously.", n_samples))

  # variance filter
  rv <- apply(expr, 1, stats::var)
  if (any(rv == 0)) {
    notes <- c(notes, sprintf("%d genes with zero variance removed.", sum(rv == 0)))
    expr <- expr[rv > 0, , drop = FALSE]
    rv <- rv[rv > 0]
  }
  if (opt$n_top_var > 0 && nrow(expr) > opt$n_top_var) {
    notes <- c(notes, sprintf("Kept the %d most variable of %d genes.", opt$n_top_var, nrow(expr)))
    expr <- expr[order(-rv)[seq_len(opt$n_top_var)], , drop = FALSE]
  }
  dat <- t(expr)
  gsg <- WGCNA::goodSamplesGenes(dat, verbose = 0)
  if (!gsg$allOK) {
    notes <- c(notes, sprintf("goodSamplesGenes removed %d genes and %d samples with too many missing values.",
                              sum(!gsg$goodGenes), sum(!gsg$goodSamples)))
    dat <- dat[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  }

  # traits and batch
  tr <- read_traits(opt$traits, rownames(dat), exclude = c(opt$batch))
  notes <- c(notes, tr$notes)
  if (nzchar(opt$batch)) {
    if (is.null(tr$raw) || !opt$batch %in% names(tr$raw)) fail("--batch column '", opt$batch, "' not in --traits.")
    prot <- split_list(opt$protect)
    miss <- setdiff(prot, names(tr$raw))
    if (length(miss)) fail("--protect column(s) not in --traits: ", paste(miss, collapse = ", "))
    dat <- remove_batch(dat, tr$raw[[opt$batch]], tr$raw[, prot, drop = FALSE], opt$batch)
    notes <- c(notes, sprintf("Batch '%s' (%d levels) removed with ComBat%s.", opt$batch,
                              length(unique(tr$raw[[opt$batch]])),
                              if (length(prot)) paste0(", protecting ", paste(prot, collapse = ", ")) else ""))
  }

  # sample outliers
  zk <- sample_zk(dat)
  flagged <- names(zk)[zk < opt$outlier_z]
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  plot_samples(dat, zk, opt$outlier_z, opt$out)
  save_table(data.frame(sample = names(zk), Zk = round(zk, 3), flagged = zk < opt$outlier_z,
                        removed = opt$remove_outliers & zk < opt$outlier_z), opt$out, "sample_connectivity.csv")
  if (length(flagged)) {
    if (opt$remove_outliers) {
      dat <- dat[!rownames(dat) %in% flagged, , drop = FALSE]
      notes <- c(notes, sprintf("Removed outlier sample(s) %s (Z.k < %g).", paste(flagged, collapse = ", "), opt$outlier_z))
      if (!is.null(tr$traits)) tr$traits <- tr$traits[rownames(dat), , drop = FALSE]
      if (nrow(dat) < 15 && !opt$allow_small_n) fail("Fewer than 15 samples remain after removing outliers.")
    } else {
      warn <- c(warn, sprintf("Sample(s) %s have standardised connectivity Z.k < %g (possible outliers); they were kept. Inspect figures/sample_outliers.png and consider --remove-outliers.",
                              paste(flagged, collapse = ", "), opt$outlier_z))
    }
  }

  res <- run_wgcna(dat, tr$traits, opt$power, opt$r2_cut, opt$cor_type, opt$max_p_outliers,
                   opt$min_module_size, opt$merge_cut_height, opt$deep_split, opt$max_block_size,
                   opt$seed, opt$alpha)
  if (opt$power == 0 && res$power_auto$fallback)
    warn <- c(warn, sprintf(paste0("No power up to 30 reached a signed scale-free fit R^2 >= %g. Power %d was taken ",
                                   "from the WGCNA FAQ table for %d samples. This usually means a strong driver ",
                                   "(batch, outlier samples or a large treatment effect) dominates the data: check ",
                                   "the sample clustering before interpreting modules."), opt$r2_cut, res$power, nrow(dat)))
  labels <- res$labels
  gene_ids <- colnames(dat)
  sizes <- table(labels)
  n_mod <- sum(names(sizes) != "0")
  if (n_mod == 0) warn <- c(warn, "No modules found: every gene is grey (unassigned).")
  grey_frac <- mean(labels == 0)
  if (grey_frac > 0.6) warn <- c(warn, sprintf("%.0f%% of genes are unassigned (grey).", 100 * grey_frac))

  # tables
  kme <- res$kme
  own <- vapply(seq_along(labels), function(i) if (labels[i] == 0) NA_real_ else kme[i, paste0("kME", labels[i])], 0)
  save_table(data.frame(gene_id = gene_ids, module = labels, module_color = WGCNA::labels2colors(labels),
                        kME_own = round(own, 4)), opt$out, "gene_modules.csv")
  save_table(data.frame(gene_id = gene_ids, round(kme, 4), check.names = FALSE), opt$out, "kME_all_modules.csv")
  save_table(data.frame(module = as.integer(names(sizes)), module_color = WGCNA::labels2colors(as.integer(names(sizes))),
                        n_genes = as.integer(sizes)), opt$out, "module_sizes.csv")
  save_table(data.frame(sample = rownames(dat), round(res$MEs, 5), check.names = FALSE), opt$out, "module_eigengenes.csv")
  fit <- res$sft
  save_table(data.frame(power = fit$Power, signed_R2 = round(res$power_auto$signed_r2, 4), slope = round(fit$slope, 4),
                        mean_k = round(fit$mean.k., 3), median_k = round(fit$median.k., 3)), opt$out, "soft_threshold.csv")
  hubs <- hub_table(labels, kme, gene_ids, opt$hub_n, opt$hub_kme)
  if (!is.null(hubs)) save_table(hubs, opt$out, "hub_genes.csv")
  edges <- module_edges(res$net, gene_ids, labels, opt$edge_threshold, opt$max_edges)
  save_table(edges, opt$out, "cytoscape_edges.csv")
  nodes <- data.frame(gene_id = gene_ids, module = labels, module_color = WGCNA::labels2colors(labels),
                      kME_own = round(own, 4), hub = gene_ids %in% hubs$gene_id[hubs$hub])
  save_table(nodes[nodes$gene_id %in% c(edges$fromNode, edges$toNode), , drop = FALSE], opt$out, "cytoscape_nodes.csv")
  mt <- res$module_trait
  if (!is.null(mt)) {
    mt$module_color <- WGCNA::labels2colors(as.integer(sub("ME", "", mt$module)))
    save_table(transform(mt, r = round(r, 4)), opt$out, "module_trait_correlations.csv")
  }

  # figures
  plot_soft_threshold(fit, res$power_auto, res$power, opt$r2_cut, opt$out)
  if (n_mod > 0) plot_dendro(res$net, opt$out)
  if (!is.null(mt) && n_mod > 0) {
    plot_module_trait(mt, res$MEs, opt$out, opt$alpha)
    plot_kme_gs(dat, labels, kme, tr$traits, mt[1, ], opt$out)
  }
  r_min <- min_detectable_r(nrow(dat), opt$alpha)
  if (!is.null(mt)) {
    small <- mt$significant & abs(mt$r) < 0.5
    if (nrow(dat) < 30)
      warn <- c(warn, sprintf(paste0("With %d samples only correlations with |r| >= %.2f reach p < %g before multiple-testing ",
                                     "correction, and estimates of r are imprecise: treat module-trait links as hypotheses."),
                              nrow(dat), r_min, opt$alpha))
    if (any(small)) warn <- c(warn, sprintf("%d significant module-trait correlation(s) have |r| < 0.5.", sum(small)))
  }

  # report
  sel <- res$power_auto
  rep <- c(sprintf("# %s report", SKILL), "",
           if (length(warn)) c("## Warnings", "", paste("-", warn), "") else NULL,
           "## Data", "",
           sprintf("- Input: `%s`; %d samples x %d genes used for the network.",
                   basename(if (nzchar(opt$expr)) opt$expr else opt$counts), nrow(dat), ncol(dat)),
           paste("-", notes),
           sprintf("- Sample outliers: Z.k < %g flagged: %s.", opt$outlier_z, if (length(flagged)) paste(flagged, collapse = ", ") else "none"),
           "", "## Network", "",
           sprintf("- Signed network, %s correlation%s; soft-threshold power %d (%s).", opt$cor_type,
                   if (opt$cor_type == "bicor") sprintf(" (maxPOutliers %g)", opt$max_p_outliers) else "", res$power,
                   if (opt$power > 0) "set by the user" else if (sel$fallback) "WGCNA FAQ fallback: no power reached the fit target"
                   else sprintf("lowest power with signed scale-free R^2 >= %g", opt$r2_cut)),
           sprintf("- Signed R^2 at the chosen power: %.2f; mean connectivity %.1f.",
                   sel$signed_r2[match(res$power, fit$Power)], fit$mean.k.[match(res$power, fit$Power)]),
           sprintf("- blockwiseModules: minModuleSize %g, mergeCutHeight %g, deepSplit %g, maxBlockSize %g (%d block(s)).",
                   opt$min_module_size, opt$merge_cut_height, opt$deep_split, opt$max_block_size, length(res$net$TOMFiles)),
           sprintf("- %d modules; %d of %d genes (%.0f%%) unassigned (grey, module 0).", n_mod, sum(labels == 0),
                   length(labels), 100 * grey_frac), "",
           "| module | colour | genes | top hub (kME) |", "|---|---|---|---|")
  for (m in sort(setdiff(as.integer(names(sizes)), 0))) {
    h <- hubs[hubs$module == m, ][1, ]
    rep <- c(rep, sprintf("| %d | %s | %d | %s (%.2f) |", m, WGCNA::labels2colors(m), sizes[[as.character(m)]], h$gene_id, h$kME))
  }
  rep <- c(rep, "", "![Soft threshold](figures/soft_threshold.png)", "", "![Samples](figures/sample_outliers.png)", "",
           if (n_mod > 0) c("![Dendrogram](figures/gene_dendrogram.png)", "") else NULL)
  if (!is.null(mt)) {
    sig <- mt[mt$significant, , drop = FALSE]
    rep <- c(rep, "## Module-trait associations", "",
             sprintf("Pearson correlation of module eigengenes with traits; Student p-values; BH correction over all %d module x trait tests (grey excluded). Minimum |r| for unadjusted p < %g with n = %d: %.2f.",
                     nrow(mt), opt$alpha, nrow(dat), r_min), "")
    if (nrow(sig)) {
      rep <- c(rep, "| module | colour | trait | r | padj |", "|---|---|---|---|---|",
               sprintf("| %s | %s | %s | %.2f | %.2e |", sig$module, sig$module_color, sig$trait, sig$r, sig$padj), "")
    } else rep <- c(rep, sprintf("No module-trait correlation at padj < %g.", opt$alpha), "")
    rep <- c(rep, "![Module-trait](figures/module_trait.png)", "", "![kME vs GS](figures/kme_vs_gene_significance.png)", "")
  }
  rep <- c(rep, "## Outputs for follow-up", "",
           sprintf("- `tables/cytoscape_edges.csv`: %d within-module edges with topological overlap >= %g (strongest %d per module); import with `cytoscape_nodes.csv` as node table.",
                   nrow(edges), opt$edge_threshold, opt$max_edges),
           "- `tables/hub_genes.csv`: top genes per module by kME; `hub` = kME >= the hub threshold.",
           "- Module gene lists (`tables/gene_modules.csv`) can go to claw-enrichment with the network genes as the universe.", "",
           "## Methods", "",
           paste0(sprintf("Co-expression modules were built with WGCNA %s: a signed network ", as.character(utils::packageVersion("WGCNA"))),
                  sprintf("(%s correlation, soft-threshold power %d) and signed topological overlap; modules were detected with blockwiseModules ",
                          opt$cor_type, res$power),
                  sprintf("(dynamic tree cut, minimum size %g, merging eigengenes correlated above %.2f). ", opt$min_module_size, 1 - opt$merge_cut_height),
                  "Module eigengenes were correlated with traits and p-values adjusted with Benjamini-Hochberg; hub genes were ranked by kME."),
           "", "## References", "",
           "- Langfelder P, Horvath S (2008). WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9:559. https://doi.org/10.1186/1471-2105-9-559",
           "- Zhang B, Horvath S (2005). A general framework for weighted gene co-expression network analysis. Statistical Applications in Genetics and Molecular Biology 4:Article17. https://doi.org/10.2202/1544-6115.1128",
           "- Langfelder P, Zhang B, Horvath S (2008). Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R. Bioinformatics 24:719-720. https://doi.org/10.1093/bioinformatics/btm563",
           "- Oldham MC, Langfelder P, Horvath S (2012). Network methods for describing sample relationships in genomic datasets: application to Huntington's disease. BMC Systems Biology 6:63. https://doi.org/10.1186/1752-0509-6-63",
           "- Langfelder P, Horvath S. WGCNA package: Frequently Asked Questions. https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html",
           "", sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(rep, file.path(opt$out, "report.md"))
  write_repro_bundle(opt$out, argv, inputs = c(opt$expr, opt$counts, opt$traits))
  cat(sprintf("Done: %d modules (power %d) from %d genes x %d samples. Report: %s\n", n_mod, res$power, ncol(dat),
              nrow(dat), normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(res)
}

if (sys.nframe() == 0L) main()
