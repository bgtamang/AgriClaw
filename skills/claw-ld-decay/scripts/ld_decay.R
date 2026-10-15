#!/usr/bin/env Rscript
# claw-ld-decay: linkage disequilibrium (LD) decay with physical distance.
# Computes r2 between marker pairs within chromosomes (uniformly subsampled
# pairs for speed), fits the Hill and Weir (1988) expectation of r2 by nonlinear
# least squares as in Remington et al. (2001), adds a model-free binned curve,
# and reports the distance at which r2 falls to 0.2 and to half its maximum,
# per chromosome and genome-wide, plus the background r2 of unlinked markers.
# The distances define candidate windows for claw-gwas-crop and claw-gene-lookup.
#
# Usage:
#   Rscript scripts/ld_decay.R --demo
#   Rscript scripts/ld_decay.R --geno genotypes_filtered_012.csv
#   Rscript scripts/ld_decay.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-ld-decay"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/ld_decay.R [--demo] --geno FILE [options] [--out DIR]

  --demo              run on the bundled synthetic population (3 chromosomes with
                      known recombination rates)
  --geno FILE         0/1/2 marker table from claw-marker-qc (marker, chr, pos, then one
                      column per line; .csv/.tsv/.gz) or a VCF (.vcf/.vcf.gz)
  --min-maf X         use markers with minor allele frequency >= X (default 0.05)
  --max-dist BP       only pairs closer than BP (default: whole chromosome)
  --max-pairs N       random marker pairs per chromosome (default 20000)
  --n-unlinked N      random pairs on different chromosomes for background r2 (default 5000)
  --r2-threshold X    r2 level for the decay distance and LD window (default 0.2)
  --window-basis B    distance written to ld_window.csv: hw (Hill-Weir fit, default)
                      or binned (model-free binned means)
  --n-bins N          distance bins for the binned curve (default 40, log-spaced)
  --seed N            random seed for pair sampling (default 2026)
  --out DIR           output folder (default: agriclaw_output/claw-ld-decay)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, geno = "", min_maf = 0.05, max_dist = NA_real_,
                 max_pairs = 20000, n_unlinked = 5000, r2_threshold = 0.2,
                 window_basis = "hw", n_bins = 40, seed = 2026, out = "")

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
COL_INK2 <- "#52514e"; COL_MUTED <- "#898781"; COL_POINT <- "#c3c2b733"

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

# Marker table (claw-marker-qc output) or VCF -> markers x samples matrix + map.
read_genotypes <- function(path) {
  l1 <- first_line(path)
  if (startsWith(l1, "##fileformat=VCF")) {
    vcf <- vcfR::read.vcfR(path, verbose = FALSE)
    fix <- vcfR::getFIX(vcf)
    gt <- vcfR::extract.gt(vcf, element = "GT", IDtoRowNames = FALSE)
    keep <- !(!is.na(fix[, "ALT"]) & grepl(",", fix[, "ALT"], fixed = TRUE))
    ids <- fix[keep, "ID"]
    noid <- is.na(ids) | ids == "."
    ids[noid] <- paste0(fix[keep, "CHROM"][noid], "_", fix[keep, "POS"][noid])
    G <- gt_to_dosage(gt[keep, , drop = FALSE])
    dimnames(G) <- list(make.unique(ids), colnames(vcf@gt)[-1])
    return(list(G = G, map = data.frame(marker = rownames(G), chr = fix[keep, "CHROM"],
                                        pos = as.numeric(fix[keep, "POS"]), stringsAsFactors = FALSE)))
  }
  n_tab <- sum(gregexpr("\t", l1, fixed = TRUE)[[1]] > 0)
  n_com <- sum(gregexpr(",", l1, fixed = TRUE)[[1]] > 0)
  con <- open_text(path)
  df <- tryCatch(utils::read.table(con, header = TRUE, sep = if (n_tab > n_com) "\t" else ",",
                                   colClasses = "character", check.names = FALSE, quote = "\"",
                                   comment.char = "", na.strings = c("", "NA"), stringsAsFactors = FALSE),
                 finally = close(con))
  nm <- names(df)
  im <- find_col(nm, ALIASES$marker); ic <- find_col(nm, ALIASES$chr); ip <- find_col(nm, ALIASES$pos)
  if (is.na(ic) || is.na(ip))
    stop("The genotype file needs marker, chr and pos columns followed by one column per line ",
         "(the table written by claw-marker-qc). Convert other layouts with claw-marker-qc first.", call. = FALSE)
  if (is.na(im)) im <- 1L
  meta <- unique(stats::na.omit(c(im, ic, ip, find_col(nm, ALIASES$ref), find_col(nm, ALIASES$alt))))
  vals <- as.matrix(df[, -meta, drop = FALSE])
  num <- suppressWarnings(matrix(as.numeric(vals), nrow(vals), ncol(vals)))
  if (any(!is.na(vals) & is.na(num))) stop("Genotype values must be numeric 0/1/2.", call. = FALSE)
  if (any(num < 0 | num > 2, na.rm = TRUE)) stop("Genotype values must lie between 0 and 2.", call. = FALSE)
  dimnames(num) <- list(make.unique(df[[im]]), colnames(vals))
  list(G = num, map = data.frame(marker = rownames(num), chr = df[[ic]],
                                 pos = suppressWarnings(as.numeric(df[[ip]])), stringsAsFactors = FALSE))
}

natural_chr_order <- function(chr) {
  u <- unique(as.character(chr))
  num <- suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", u)))
  num[!grepl("[0-9]", u)] <- NA
  u[order(is.na(num), num, u)]
}

# Uniform random sample (without replacement) of marker pairs i < j with
# 0 < pos_j - pos_i <= max_dist; all pairs when there are at most max_pairs.
# Returns indices into the original (unsorted) vector.
sample_pairs <- function(pos, max_dist = Inf, max_pairs = 20000) {
  o <- order(pos)
  ps <- pos[o]
  m <- length(ps)
  if (m < 2L) return(data.frame(i = integer(), j = integer()))
  hi <- findInterval(ps + max_dist, ps)
  cnt <- pmax(hi - seq_len(m), 0L)
  total <- sum(as.numeric(cnt))
  if (total == 0) return(data.frame(i = integer(), j = integer()))
  if (total <= max_pairs) {
    i <- rep.int(seq_len(m), cnt)
    j <- i + sequence(cnt)
  } else {
    k <- sample.int(total, max_pairs)
    cs <- cumsum(as.numeric(cnt))
    i <- findInterval(k - 1, cs) + 1L
    j <- i + as.integer(k - c(0, cs)[i])
  }
  keep <- ps[j] > ps[i]   # pairs at the same position carry no distance information
  data.frame(i = o[i[keep]], j = o[j[keep]])
}

# Squared Pearson correlation between columns I and J of X (samples x markers),
# using samples called at both markers. For inbred lines (0/2) this is the
# haplotype r2; for unphased outbred genotypes it is the genotypic (composite)
# r2 (Rogers and Huff 2009).
r2_pairs <- function(X, I, J, chunk = 5000L) {
  r2 <- n_obs <- numeric(length(I))
  for (s in seq(1L, length(I), by = chunk)) {
    k <- s:min(length(I), s + chunk - 1L)
    A <- X[, I[k], drop = FALSE]; B <- X[, J[k], drop = FALSE]
    M <- !is.na(A) & !is.na(B)
    A[!M] <- 0; B[!M] <- 0
    n <- colSums(M)
    sa <- colSums(A); sb <- colSums(B)
    saa <- colSums(A * A); sbb <- colSums(B * B); sab <- colSums(A * B)
    va <- n * saa - sa^2; vb <- n * sbb - sb^2
    r <- (n * sab - sa * sb) / sqrt(va * vb)
    r[va <= 1e-12 | vb <= 1e-12 | n < 3] <- NA
    r2[k] <- r^2
    n_obs[k] <- n
  }
  list(r2 = r2, n_obs = n_obs)
}

# Hill and Weir (1988) expectation of r2 with low mutation and a sample of n,
# as written by Remington et al. (2001), with C = rho * distance and rho the
# population recombination parameter (4 N c) per base pair.
hw_expect <- function(d, rho, n) {
  C <- rho * d
  ((10 + C) / ((2 + C) * (11 + C))) * (1 + ((3 + C) * (12 + 12 * C + C^2)) / (n * (2 + C) * (11 + C)))
}

# Least-squares fit of rho: grid search over log10(rho) then optimize().
fit_hw <- function(d, r2, n) {
  ok <- !is.na(r2) & !is.na(d) & d > 0
  d <- d[ok]; r2 <- r2[ok]
  if (length(d) < 10L) return(list(rho = NA_real_, rss = NA_real_, at_boundary = NA, n_pairs = length(d)))
  rss <- function(l) sum((r2 - hw_expect(d, 10^l, n))^2)
  grid <- seq(-12, 0, by = 0.1)
  v <- vapply(grid, rss, 0)
  b <- which.min(v)
  if (b == 1L || b == length(grid)) {
    return(list(rho = 10^grid[b], rss = v[b], at_boundary = TRUE, n_pairs = length(d)))
  }
  o <- stats::optimize(rss, c(grid[b - 1L], grid[b + 1L]), tol = 1e-10)
  list(rho = 10^o$minimum, rss = o$objective, at_boundary = FALSE, n_pairs = length(d))
}

# Distance (bp) at which the fitted Hill-Weir curve falls to `target`.
hw_distance <- function(target, rho, n) {
  if (is.na(rho) || is.na(target)) return(NA_real_)
  f <- function(C) hw_expect(C, 1, n) - target
  if (f(0) <= 0) return(0)
  hiC <- 1e9
  if (f(hiC) >= 0) return(NA_real_)   # the curve never falls that low (floor is about 1/n)
  stats::uniroot(f, c(0, hiC), tol = 1e-10)$root / rho
}

# Model-free curve: mean r2 in log-spaced distance bins (bins with >= min_n pairs).
bin_r2 <- function(d, r2, n_bins = 40, min_n = 20) {
  ok <- !is.na(r2) & d > 0
  d <- d[ok]; r2 <- r2[ok]
  empty <- data.frame(bin_start_bp = numeric(), bin_end_bp = numeric(), mean_dist_bp = numeric(),
                      mean_r2 = numeric(), n_pairs = integer())
  if (length(d) < min_n) return(empty)
  br <- exp(seq(log(min(d)), log(max(d)), length.out = n_bins + 1L))
  br[1] <- min(d) * (1 - 1e-9)                 # guard against exp(log(x)) rounding
  br[length(br)] <- max(d) * (1 + 1e-9)
  b <- findInterval(d, br, rightmost.closed = TRUE)
  out <- data.frame(bin_start_bp = br[-length(br)], bin_end_bp = br[-1],
                    mean_dist_bp = as.numeric(tapply(d, factor(b, levels = seq_len(n_bins)), mean)),
                    mean_r2 = as.numeric(tapply(r2, factor(b, levels = seq_len(n_bins)), mean)),
                    n_pairs = as.integer(table(factor(b, levels = seq_len(n_bins)))))
  out[out$n_pairs >= min_n, , drop = FALSE]
}

# First distance at which the binned curve falls below `target` (linear
# interpolation between bin mean distances). NA if it never does or if the
# first bin is already below the target.
binned_distance <- function(bins, target) {
  if (!nrow(bins) || is.na(target)) return(NA_real_)
  y <- bins$mean_r2; x <- bins$mean_dist_bp
  if (y[1] < target) return(NA_real_)
  k <- which(y < target)[1]
  if (is.na(k)) return(NA_real_)
  x[k - 1] + (target - y[k - 1]) * (x[k] - x[k - 1]) / (y[k] - y[k - 1])
}

# Background r2 from markers on different chromosomes: parametric 95th
# percentile of sqrt(r2), squared (Breseghello and Sorrells 2006).
background_r2 <- function(X, chr, n_pairs = 5000) {
  if (length(unique(chr)) < 2L) return(list(r2_p95 = NA_real_, n_pairs = 0L))
  i <- sample.int(length(chr), n_pairs * 3L, replace = TRUE)
  j <- sample.int(length(chr), n_pairs * 3L, replace = TRUE)
  keep <- chr[i] != chr[j]
  i <- i[keep][seq_len(min(n_pairs, sum(keep)))]
  j <- j[keep][seq_len(min(n_pairs, sum(keep)))]
  r <- r2_pairs(X, i, j)$r2
  s <- sqrt(r[!is.na(r)])
  list(r2_p95 = (mean(s) + stats::qnorm(0.95) * stats::sd(s))^2, n_pairs = length(s))
}

decay_summary <- function(label, d, r2, n, n_markers, target, bins_n, bg) {
  fit <- fit_hw(d, r2, n)
  bins <- bin_r2(d, r2, bins_n)
  r2max_fit <- hw_expect(0, 1, n)
  bmax <- if (nrow(bins)) max(bins$mean_r2) else NA_real_
  dmin <- if (any(d > 0)) min(d[d > 0]) else NA_real_
  dh <- hw_distance(target, fit$rho, n)
  data.frame(chr = label, n_markers = n_markers, n_pairs = fit$n_pairs, n_samples = n,
             min_pair_dist_bp = dmin, max_pair_dist_bp = if (length(d)) max(d) else NA_real_,
             rho_per_bp = fit$rho, hw_fit_at_boundary = fit$at_boundary,
             r2_threshold = target,
             dist_r2_threshold_hw_bp = dh,
             hw_extrapolated = !is.na(dh) && !is.na(dmin) && dh < dmin,
             r2_max_hw = r2max_fit,
             dist_half_max_hw_bp = hw_distance(r2max_fit / 2, fit$rho, n),
             dist_r2_threshold_binned_bp = binned_distance(bins, target),
             r2_max_binned = bmax,
             dist_half_max_binned_bp = binned_distance(bins, bmax / 2),
             background_r2_p95 = bg,
             dist_background_hw_bp = hw_distance(bg, fit$rho, n),
             stringsAsFactors = FALSE)
}

# Candidate-window distance per chromosome. A Hill-Weir distance is usable only
# if the fit did not hit the search edge, the distance is not an extrapolation
# below the shortest sampled pair, and the binned data actually reach the
# threshold (the model's maximum is fixed near 0.45, so it "crosses" 0.2 even
# when the data never do). Otherwise the other estimate is used, and if neither
# is usable the window is NA with the reason given.
ld_window <- function(s, basis = "hw", target = 0.2) {
  hw_ok <- !is.na(s$dist_r2_threshold_hw_bp) & !s$hw_extrapolated & !(s$hw_fit_at_boundary %in% TRUE) &
    !is.na(s$r2_max_binned) & s$r2_max_binned >= target
  bin_ok <- !is.na(s$dist_r2_threshold_binned_bp)
  first <- if (basis == "hw") hw_ok else bin_ok
  second <- if (basis == "hw") bin_ok else hw_ok
  v1 <- if (basis == "hw") s$dist_r2_threshold_hw_bp else s$dist_r2_threshold_binned_bp
  v2 <- if (basis == "hw") s$dist_r2_threshold_binned_bp else s$dist_r2_threshold_hw_bp
  n1 <- if (basis == "hw") "Hill-Weir fit" else "binned means"
  n2 <- if (basis == "hw") "binned means" else "Hill-Weir fit"
  d <- ifelse(first, v1, ifelse(second, v2, NA_real_))
  why <- ifelse(!is.na(s$r2_max_binned) & s$r2_max_binned < target,
                sprintf("not estimable: binned r2 stays below %g at all sampled distances (max %.3f)", target, s$r2_max_binned),
                "not estimable: Hill-Weir estimate is an extrapolation or failed, and the binned curve never crosses")
  b <- ifelse(first, sprintf("r2 = %g, %s", target, n1),
              ifelse(second, sprintf("r2 = %g, %s (fallback)", target, n2), why))
  data.frame(chr = s$chr, ld_distance_bp = round(d), basis = b, stringsAsFactors = FALSE)
}

# --- Plots ---------------------------------------------------------------------

dist_unit <- function(maxd) if (maxd >= 2e6) list(div = 1e6, lab = "Mb") else list(div = 1e3, lab = "kb")

draw_decay <- function(d, r2, bins, rho, n, target, dist_t, bg, main, points = TRUE, cex_main = 1) {
  u <- dist_unit(max(d, na.rm = TRUE))
  graphics::plot(NA, xlim = c(0, max(d, na.rm = TRUE) / u$div), ylim = c(0, 1),
                 xlab = sprintf("Distance (%s)", u$lab), ylab = expression(r^2), bty = "l",
                 main = main, col.main = COL_INK, cex.main = cex_main)
  if (points) graphics::points(d / u$div, r2, pch = 16, cex = 0.3, col = COL_POINT)
  if (nrow(bins)) graphics::points(bins$mean_dist_bp / u$div, bins$mean_r2, pch = 19, cex = 0.7, col = COL_ORANGE)
  if (!is.na(rho)) {
    xs <- seq(0, max(d, na.rm = TRUE), length.out = 400)
    graphics::lines(xs / u$div, hw_expect(xs, rho, n), col = COL_BLUE, lwd = 2)
  }
  graphics::abline(h = target, col = COL_INK2, lty = 2)
  if (!is.na(bg)) graphics::abline(h = bg, col = COL_MUTED, lty = 3)
  if (!is.na(dist_t)) graphics::abline(v = dist_t / u$div, col = COL_BLUE, lty = 3)
}

plot_genome <- function(d, r2, bins, s, out, target) {
  open_png(out, "ld_decay_genome.png", width = 8, height = 5)
  op <- graphics::par(mar = c(4.5, 4.5, 3, 1), col.axis = COL_INK2, fg = COL_INK2)
  draw_decay(d, r2, bins, s$rho_per_bp, s$n_samples, target, s$dist_r2_threshold_hw_bp,
             s$background_r2_p95, "LD decay, all chromosomes")
  graphics::legend("topright", bty = "n", cex = 0.85, text.col = COL_INK2,
                   pch = c(16, 19, NA, NA, NA, NA), lty = c(NA, NA, 1, 2, 3, 3), lwd = c(NA, NA, 2, 1, 1, 1),
                   col = c("#a3a29b", COL_ORANGE, COL_BLUE, COL_INK2, COL_MUTED, COL_BLUE),
                   legend = c("marker pairs", "binned mean", "Hill-Weir fit", sprintf("r2 = %g", target),
                              "background r2 (unlinked, 95th pct)", "decay distance (Hill-Weir)"))
  graphics::par(op)
  grDevices::dev.off()
}

plot_by_chr <- function(per_chr, pairs, out, target, n) {
  chrs <- per_chr$chr
  k <- length(chrs)
  nc <- ceiling(sqrt(k)); nr <- ceiling(k / nc)
  open_png(out, "ld_decay_by_chr.png", width = min(12, 3.2 * nc), height = min(12, 2.8 * nr))
  op <- graphics::par(mfrow = c(nr, nc), mar = c(4, 4, 2, 0.5), col.axis = COL_INK2, fg = COL_INK2)
  for (ch in chrs) {
    p <- pairs[pairs$chr == ch, , drop = FALSE]
    s <- per_chr[per_chr$chr == ch, ]
    draw_decay(p$dist_bp, p$r2, bin_r2(p$dist_bp, p$r2), s$rho_per_bp, n, target,
               s$dist_r2_threshold_hw_bp, s$background_r2_p95, paste("Chromosome", ch),
               points = FALSE, cex_main = 0.9)
  }
  graphics::par(op)
  grDevices::dev.off()
}

# --- Main ---------------------------------------------------------------------

fmt_bp <- function(x) ifelse(is.na(x), "NA", ifelse(x >= 1e6, sprintf("%.2f Mb", x / 1e6), sprintf("%.1f kb", x / 1e3)))

run_ld <- function(G, map, opt) {
  maf <- {
    af <- rowMeans(G, na.rm = TRUE) / 2
    pmin(af, 1 - af)
  }
  ncall <- rowSums(!is.na(G))
  keep <- !is.na(maf) & maf >= opt$min_maf & ncall >= 20 & !is.na(map$chr) & !is.na(map$pos)
  if (sum(keep) < 10L) stop("Fewer than 10 markers pass the MAF and missing-data filters.", call. = FALSE)
  G <- G[keep, , drop = FALSE]
  map <- map[keep, , drop = FALSE]
  X <- t(G)
  n <- nrow(X)
  set.seed(opt$seed)
  chrs <- natural_chr_order(map$chr)
  maxd <- if (is.na(opt$max_dist)) Inf else opt$max_dist
  pair_list <- list()
  for (ch in chrs) {
    idx <- which(map$chr == ch)
    sp <- sample_pairs(map$pos[idx], maxd, opt$max_pairs)
    if (!nrow(sp)) next
    I <- idx[sp$i]; J <- idx[sp$j]
    rr <- r2_pairs(X, I, J)
    pair_list[[ch]] <- data.frame(chr = ch, marker_1 = map$marker[I], marker_2 = map$marker[J],
                                  dist_bp = abs(map$pos[J] - map$pos[I]), r2 = rr$r2, n_obs = rr$n_obs,
                                  stringsAsFactors = FALSE)
  }
  if (!length(pair_list)) stop("No marker pairs within chromosomes: check chr and pos.", call. = FALSE)
  pairs <- do.call(rbind, pair_list)
  pairs <- pairs[!is.na(pairs$r2), , drop = FALSE]
  n_eff <- stats::median(pairs$n_obs)
  bg <- background_r2(X, map$chr, opt$n_unlinked)
  per_chr <- do.call(rbind, lapply(names(pair_list), function(ch) {
    p <- pairs[pairs$chr == ch, , drop = FALSE]
    decay_summary(ch, p$dist_bp, p$r2, n_eff, sum(map$chr == ch), opt$r2_threshold, opt$n_bins, bg$r2_p95)
  }))
  genome <- decay_summary("all", pairs$dist_bp, pairs$r2, n_eff, nrow(map), opt$r2_threshold, opt$n_bins, bg$r2_p95)
  list(pairs = pairs, per_chr = per_chr, genome = genome, n_samples = n, n_eff = n_eff,
       n_markers = nrow(map), n_dropped = sum(!keep), background = bg)
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  options(scipen = 100)
  if (opt$demo) opt$geno <- file.path(SKILL_DIR, "demo", "demo_ld_012.csv.gz")
  if (!nzchar(opt$geno)) stop("Provide --geno FILE, or --demo. See --help.", call. = FALSE)
  if (opt$r2_threshold <= 0 || opt$r2_threshold >= 1) stop("--r2-threshold must be between 0 and 1.", call. = FALSE)
  if (!opt$window_basis %in% c("hw", "binned")) stop("--window-basis must be hw or binned.", call. = FALSE)
  if (opt$min_maf < 0 || opt$min_maf >= 0.5) stop("--min-maf must be in [0, 0.5).", call. = FALSE)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

  gd <- read_genotypes(opt$geno)
  res <- run_ld(gd$G, gd$map, opt)
  s_all <- rbind(res$per_chr, res$genome)
  save_table(s_all, opt$out, "ld_decay_summary.csv")
  bins_all <- bin_r2(res$pairs$dist_bp, res$pairs$r2, opt$n_bins)
  bins_chr <- do.call(rbind, lapply(unique(res$pairs$chr), function(ch) {
    p <- res$pairs[res$pairs$chr == ch, ]
    b <- bin_r2(p$dist_bp, p$r2, opt$n_bins)
    if (nrow(b)) cbind(chr = ch, b) else NULL
  }))
  save_table(rbind(cbind(chr = "all", bins_all), bins_chr), opt$out, "ld_binned.csv")
  win <- ld_window(s_all, opt$window_basis, opt$r2_threshold)
  save_table(win, opt$out, "ld_window.csv")

  plot_genome(res$pairs$dist_bp, res$pairs$r2, bins_all, res$genome, opt$out, opt$r2_threshold)
  plot_by_chr(res$per_chr, res$pairs, opt$out, opt$r2_threshold, res$n_eff)

  # Warnings
  warn <- character()
  for (k in seq_len(nrow(s_all))) {
    s <- s_all[k, ]
    lab <- if (s$chr == "all") "Genome-wide" else paste("Chromosome", s$chr)
    if (isTRUE(s$hw_fit_at_boundary)) warn <- c(warn, sprintf("%s: the Hill-Weir fit hit the edge of the search range; use the binned estimate.", lab))
    dh <- s$dist_r2_threshold_hw_bp; db <- s$dist_r2_threshold_binned_bp
    if (!is.na(dh) && !is.na(db) && max(dh, db) / max(min(dh, db), 1) > 2)
      warn <- c(warn, sprintf("%s: Hill-Weir (%s) and binned (%s) distances differ more than 2-fold; the Hill-Weir model fits poorly here, so report both and prefer the binned value.",
                              lab, fmt_bp(dh), fmt_bp(db)))
    dd <- win$ld_distance_bp[k]
    if (!is.na(dd) && dd > 0.5 * s$max_pair_dist_bp)
      warn <- c(warn, sprintf("%s: the decay distance (%s) is more than half of the longest pair distance; LD decay is not resolved. Increase --max-dist or check the map.", lab, fmt_bp(dd)))
    if (s$chr != "all" && !is.na(dd)) {
      p <- res$pairs[res$pairs$chr == s$chr, ]
      if (mean(p$dist_bp < dd) < 0.05)
        warn <- c(warn, sprintf("%s: under 5%% of sampled pairs are closer than the decay distance (%s); short-range LD is poorly sampled. Re-run with --max-dist about 10 times that distance.", lab, fmt_bp(dd)))
    }
    if (is.na(dh) && is.na(db)) warn <- c(warn, sprintf("%s: r2 never falls to %g (or starts below it); no decay distance.", lab, opt$r2_threshold))
    if (!is.na(dh) && !isTRUE(s$hw_extrapolated) && !is.na(s$r2_max_binned) && s$r2_max_binned < opt$r2_threshold)
      warn <- c(warn, sprintf(
        "%s: the Hill-Weir curve crosses r2 = %g at %s only because its maximum is fixed near 0.45; the binned data never exceed %.3f, so that distance was not used.",
        lab, opt$r2_threshold, fmt_bp(dh), s$r2_max_binned))
    if (isTRUE(s$hw_extrapolated)) warn <- c(warn, sprintf(
      "%s: the Hill-Weir distance (%s) is shorter than the closest sampled pair (%s): it is an extrapolation and was not used for the LD window.",
      lab, fmt_bp(dh), fmt_bp(s$min_pair_dist_bp)))
    if (!is.na(s$r2_max_binned) && !is.na(s$background_r2_p95) && s$r2_max_binned <= s$background_r2_p95)
      warn <- c(warn, sprintf("%s: LD between linked markers (max binned mean r2 %.3f) is no higher than between unlinked markers (%.3f): no LD decay is detectable at the sampled distances. Markers may be too sparse for the LD extent, or structure dominates.",
                              lab, s$r2_max_binned, s$background_r2_p95))
  }
  na_win <- win$chr[is.na(win$ld_distance_bp)]
  if (length(na_win)) warn <- c(warn, sprintf(
    "No LD window could be estimated for %s (see the basis column of ld_window.csv). Do not invent one: try a lower --r2-threshold, denser markers, or ask the user for a window.",
    paste(ifelse(na_win == "all", "the genome-wide row", paste("chromosome", na_win)), collapse = ", ")))
  if (any(win$chr != "all" & is.na(win$ld_distance_bp)) && !any(win$chr == "all" & is.na(win$ld_distance_bp)))
    warn <- c(warn, "Chromosomes without a window fall back to the genome-wide 'all' row in claw-gwas-crop and claw-gene-lookup.")
  if (!is.na(res$background$r2_p95) && res$background$r2_p95 > 0.1)
    warn <- c(warn, sprintf("Background r2 between unlinked markers is %.3f (> 0.1): population structure or relatedness inflates r2, so decay distances are overestimated. Structure-corrected LD (r2 with kinship, Mangin et al. 2012) is not implemented here.",
                            res$background$r2_p95))
  if (res$n_samples < 50) warn <- c(warn, sprintf("Only %d samples: r2 estimates are noisy (rule of thumb: 50 or more).", res$n_samples))
  small <- res$per_chr$chr[res$per_chr$n_markers < 50]
  if (length(small)) warn <- c(warn, sprintf("Chromosome(s) with fewer than 50 markers: %s. Their estimates are unreliable.", paste(small, collapse = ", ")))

  tab <- c("| chr | markers | pairs | r2 = t (Hill-Weir) | r2 = t (binned) | half max (binned) | half max (Hill-Weir) |",
           "|---|---|---|---|---|---|---|",
           sprintf("| %s | %d | %d | %s | %s | %s | %s |", s_all$chr, s_all$n_markers, s_all$n_pairs,
                   fmt_bp(s_all$dist_r2_threshold_hw_bp), fmt_bp(s_all$dist_r2_threshold_binned_bp),
                   fmt_bp(s_all$dist_half_max_binned_bp), fmt_bp(s_all$dist_half_max_hw_bp)))
  g <- res$genome
  report <- c(
    sprintf("# %s report", SKILL), "",
    if (length(warn)) c("## Warnings", "", paste("-", warn), "") else NULL,
    "## Summary", "",
    sprintf("- Input: `%s`", opt$geno),
    sprintf("- Samples: %d; markers used: %d (MAF >= %.2f; %d dropped)", res$n_samples, res$n_markers, opt$min_maf, res$n_dropped),
    sprintf("- Pairs: up to %d random pairs per chromosome%s; %d pairs in total", opt$max_pairs,
            if (is.na(opt$max_dist)) " (any distance)" else sprintf(" closer than %s", fmt_bp(opt$max_dist)), nrow(res$pairs)),
    sprintf("- Genome-wide distance at r2 = %g: **%s** (Hill-Weir fit%s), %s (binned means); LD window: %s (%s)", opt$r2_threshold,
            fmt_bp(g$dist_r2_threshold_hw_bp), if (isTRUE(g$hw_extrapolated)) ", extrapolated below the closest pair: not used" else "",
            fmt_bp(g$dist_r2_threshold_binned_bp), fmt_bp(win$ld_distance_bp[win$chr == "all"]), win$basis[win$chr == "all"]),
    sprintf("- Genome-wide half-maximum distance (binned): %s (maximum binned mean r2 %.2f)", fmt_bp(g$dist_half_max_binned_bp), g$r2_max_binned),
    sprintf("- Background r2 (unlinked markers, 95th percentile of sqrt(r2), Breseghello and Sorrells 2006): %.3f; %s",
            g$background_r2_p95,
            if (is.na(g$dist_background_hw_bp) || g$dist_background_hw_bp > g$max_pair_dist_bp)
              "the Hill-Weir curve does not reach it within the observed distances"
            else sprintf("the Hill-Weir curve reaches it at %s", fmt_bp(g$dist_background_hw_bp))),
    sprintf("- LD window written to `tables/ld_window.csv` (basis: r2 = %g, %s) for claw-gwas-crop --ld-table and claw-gene-lookup --ld-table.",
            opt$r2_threshold, if (opt$window_basis == "hw") "Hill-Weir fit" else "binned means"), "",
    sprintf("t = %g. The Hill-Weir maximum is fixed by the model at about %.2f, so its half-maximum distance is close to its r2 = 0.2 distance; the binned half-maximum is data-driven.",
            opt$r2_threshold, g$r2_max_hw), "",
    tab, "",
    "![Genome-wide LD decay](figures/ld_decay_genome.png)", "", "![By chromosome](figures/ld_decay_by_chr.png)", "",
    "## Method", "",
    sprintf("r2 is the squared correlation of allele counts between two markers over samples called at both. The Hill and Weir (1988) expectation E(r2) = [(10 + C)/((2 + C)(11 + C))][1 + (3 + C)(12 + 12C + C^2)/(n(2 + C)(11 + C))], with C = rho x distance and n = %.0f samples, was fitted by least squares as in Remington et al. (2001). Binned means use %d log-spaced distance bins with at least 20 pairs each. Random seed: %d.",
            res$n_eff, opt$n_bins, opt$seed), "",
    "## Limitations", "",
    "- The Hill-Weir model assumes drift-recombination equilibrium; bottlenecks, selection, structure and few founders violate it. Treat the fitted curve as a summary.",
    "- r2 from unphased heterozygous genotypes is the genotypic (composite) r2.",
    "- Structure and relatedness inflate r2 at all distances; this skill does not correct for them.",
    "- Decay distances vary along chromosomes (centromeres, recombination hot and cold spots); one value per chromosome is an average.", "",
    "## References", "",
    "- Hill WG, Weir BS (1988) Variances and covariances of squared linkage disequilibria in finite populations. Theor Popul Biol 33:54-78. https://doi.org/10.1016/0040-5809(88)90004-4",
    "- Remington DL et al. (2001) Structure of linkage disequilibrium and phenotypic associations in the maize genome. PNAS 98:11479-11484. https://doi.org/10.1073/pnas.201394398",
    "- Breseghello F, Sorrells ME (2006) Association mapping of kernel size and milling quality in wheat (Triticum aestivum L.) cultivars. Genetics 172:1165-1177. https://doi.org/10.1534/genetics.105.044586", "",
    sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(report, file.path(opt$out, "report.md"))
  write_repro_bundle(opt$out, argv, inputs = opt$geno)
  cat(sprintf("Done: %d chromosomes, %d pairs; genome-wide LD window %s (%s). Report: %s\n",
              nrow(res$per_chr), nrow(res$pairs), fmt_bp(win$ld_distance_bp[win$chr == "all"]),
              win$basis[win$chr == "all"], normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(res)
}

if (sys.nframe() == 0L) main()
