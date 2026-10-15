# Tests for claw-de-analysis, part 2: LRT time course, tximport offsets and the command
# line. Run from anywhere: Rscript tests/test_de_lrt_cli.R
# Known answers: planted genotype-by-time interactions (empirical FDR and power at the
# calibrated levels), isoform-switching genes with and without tximport length offsets,
# and an end-to-end demo run checked against the planted truth. Part 1 is test_de_analysis.R.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "de_analysis.R"))  # loads functions only

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
expect_error <- function(expr, pattern, what) {
  msg <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  check(!is.na(msg) && grepl(pattern, msg), paste0(what, " (got: ", msg, ")"))
}
run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "de_analysis.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
T0 <- Sys.time()
stamp <- function(what) cat(sprintf("  [%5.1fs] %s\n", as.numeric(difftime(Sys.time(), T0, units = "secs")), what))
tmp <- file.path(tempdir(), "claw-de-lrt-cli-test")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE)

# Negative-binomial counts: two conditions x n_per replicates in two balanced batches.
sim <- function(seed, n_genes = 1200, n_per = 4, prop_de = 0.1, lfc = c(1, 2)) {
  set.seed(seed)
  base <- stats::rlnorm(n_genes, log(150), 1.5)
  disp <- (0.05 + 1 / base) * exp(stats::rnorm(n_genes, 0, 0.3))
  de <- seq_len(n_genes) %in% sample(n_genes, round(prop_de * n_genes))
  true <- ifelse(de, sample(c(-1, 1), n_genes, TRUE) * stats::runif(n_genes, lfc[1], lfc[2]), 0)
  s <- data.frame(sample = sprintf("s%02d", seq_len(2 * n_per)),
                  condition = rep(c("control", "treated"), each = n_per),
                  batch = rep(rep(c("b1", "b2"), length.out = n_per), 2), stringsAsFactors = FALSE)
  rownames(s) <- s$sample
  bt <- stats::rnorm(n_genes, 0, 0.3)
  sf <- exp(stats::rnorm(nrow(s), 0, 0.25))
  m <- sapply(seq_len(nrow(s)), function(j) stats::rnbinom(n_genes, size = 1 / disp,
         mu = base * sf[j] * 2^(true * (s$condition[j] == "treated") + bt * (s$batch[j] == "b2"))))
  dimnames(m) <- list(sprintf("g%05d", seq_len(n_genes)), s$sample)
  list(counts = m, samples = s, de = de, true = true, base = base, disp = disp, sf = sf, bt = bt)
}
design <- ~ batch + condition
ref <- list(condition = "control")
prep <- function(x) prepare_design(x$samples, design, ref)$samples


stamp("6. LRT time course")
# 6. LRT time course: genotype-by-time interaction (2 seeds) -------------------------
# AgriClaw calibration (10 seeds, this design): DESeq2's LRT p-values are anti-conservative
# with 3 replicates per cell (null P(p < 0.01) = 0.020); empirical FDR 0.098 at nominal 0.05
# (10 seeds) and 0.058 at nominal 0.025 (8 seeds). The checks below encode those measured levels.
lrt_sim <- function(seed, ng = 1000) {
  set.seed(seed)
  st <- expand.grid(rep = 1:3, time = c("t0", "t6", "t24"), genotype = c("WT", "mut"), stringsAsFactors = FALSE)
  st$sample <- sprintf("x%02d", seq_len(nrow(st)))
  rownames(st) <- st$sample
  base <- stats::rlnorm(ng, log(200), 1.3)
  disp <- (0.05 + 1 / base) * exp(stats::rnorm(ng, 0, 0.3))
  idx <- sample(ng)
  inter <- seq_len(ng) %in% idx[1:100]          # genotype-specific time response
  time_only <- seq_len(ng) %in% idx[101:270]    # same time response in both genotypes
  geno_only <- seq_len(ng) %in% idx[271:370]    # constant genotype difference
  eff <- function(j) {
    t6 <- st$time[j] == "t6"; t24 <- st$time[j] == "t24"; mt <- st$genotype[j] == "mut"
    1.5 * time_only * (t6 + t24) + 1.2 * geno_only * mt + inter * mt * (t6 + 2 * t24)  # 2x at 6 h, 4x at 24 h
  }
  m <- sapply(seq_len(nrow(st)), function(j) stats::rnbinom(ng, size = 1 / disp, mu = base * 2^eff(j)))
  dimnames(m) <- list(sprintf("g%05d", seq_len(ng)), st$sample)
  list(m = m, st = st, inter = inter, main = time_only | geno_only)
}
d_full <- ~ genotype + time + genotype:time
lrt <- lapply(7:8, function(seed) {
  x <- lrt_sim(seed)
  st2 <- prepare_design(x$st, d_full, list(genotype = "WT", time = "t0"))$samples
  r <- run_de(x$m, st2, d_full, test = "LRT", reduced = ~ genotype + time, shrink = "none")
  padj <- stats::setNames(r$tables[[1]]$padj, r$tables[[1]]$gene_id)[rownames(x$m)]
  at <- function(a) { call <- !is.na(padj) & padj < a; c(fdr = mean(!x$inter[call]), power = mean(call[x$inter]),
                                                         main = mean(call[x$main])) }
  list(a05 = at(0.05), a025 = at(0.025), warn = r$warn)
})
f05 <- mean(sapply(lrt, function(z) z$a05[["fdr"]]))
f025 <- mean(sapply(lrt, function(z) z$a025[["fdr"]]))
pw <- mean(sapply(lrt, function(z) z$a05[["power"]]))
cat(sprintf("  LRT (2 seeds): FDR %.3f at nominal 0.05, %.3f at 0.025; power %.2f\n", f05, f025, pw))
check(f025 <= 0.05 + 0.03, sprintf("LRT at nominal 0.025: empirical FDR %.3f <= 0.08", f025))
check(f05 <= 0.13, sprintf("LRT at nominal 0.05: empirical FDR %.3f <= 0.13 (known anti-conservative; see SKILL.md)", f05))
check(pw >= 0.6, sprintf(paste0("LRT power %.2f >= 0.6 for interactions of 2-fold at 6 h ",
      "and 4-fold at 24 h (n = 3 per genotype and time)"), pw))
check(mean(sapply(lrt, function(z) z$a05[["main"]])) <= 0.05,
      "LRT does not call genes with only main effects (they are nulls for the interaction)")
check(any(grepl("anti-conservative", lrt[[1]]$warn)), "LRT with 3 replicates per cell carries the calibration warning")

stamp("7. tximport length offsets")
# 7. tximport length offsets remove false calls caused by isoform switching ----------
set.seed(11)
xo <- sim(11, n_genes = 1000, n_per = 3, prop_de = 0.1)
switch <- !xo$de & seq_len(1000) %in% sample(which(!xo$de), 100)
lenm <- matrix(1500, 1000, ncol(xo$counts), dimnames = dimnames(xo$counts))
trt <- xo$samples$condition == "treated"
lenm[switch, trt] <- 3000                     # longer isoform dominates: 2x reads, same abundance
cnt_sw <- xo$counts
for (j in which(trt)) {                       # re-draw with doubled mean (same dispersion)
  mu <- 2 * xo$base[switch] * xo$sf[j] * 2^(xo$bt[switch] * (xo$samples$batch[j] == "b2"))
  cnt_sw[switch, j] <- stats::rnbinom(sum(switch), size = 1 / xo$disp[switch], mu = mu)
}
so <- prepare_design(xo$samples, design, ref)$samples
with_off <- run_de(cnt_sw, so, design, lengths = lenm, shrink = "none", reference = ref)$tables[[1]]
no_off <- run_de(cnt_sw, so, design, shrink = "none", reference = ref)$tables[[1]]
f_with <- mean(stats::setNames(with_off$significant, with_off$gene_id)[rownames(cnt_sw)][switch])
f_no <- mean(stats::setNames(no_off$significant, no_off$gene_id)[rownames(cnt_sw)][switch])
check(f_with <= 0.1 && f_no >= 0.5, sprintf(paste0("isoform-switch genes called DE: %.2f with tximport ",
      "length offsets (<= 0.10) vs %.2f without (>= 0.50)"), f_with, f_no))

stamp("8. End to end")
# 8. End to end on the demo ---------------------------------------------------------
out <- file.path(tmp, "demo_out")
res <- run_cli(c("--demo", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(utils::tail(res$output, 5), collapse = "\n")))
for (f in c("report.md", "tables/de_condition_drought_vs_control.csv", "tables/de_summary.csv",
            "tables/vst.csv", "tables/normalized_counts.csv", "tables/size_factors.csv",
            "tables/samples_used.csv", "figures/ma_condition_drought_vs_control.png",
            "figures/volcano_condition_drought_vs_control.png",
            "figures/pvalue_hist_condition_drought_vs_control.png", "figures/pca.png", "figures/dispersion.png",
            "reproducibility/commands.sh", "reproducibility/sessionInfo.txt"))
  check(file.exists(file.path(out, f)), paste("demo writes", f))
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
dd <- utils::read.csv(file.path(out, "tables", "de_condition_drought_vs_control.csv"))
truth <- utils::read.csv(file.path(SKILL_TEST_DIR, "demo", "truth.csv"))
m <- merge(dd, truth, by = "gene_id")
check(nrow(m) == 1200L, "demo table lists every input gene (tested, filtered or pre-filtered)")
dfdr <- mean(!m$is_de[m$significant])
dpow <- mean(m$significant[m$is_de])
check(dfdr <= 0.1 && dpow >= 0.7,
      sprintf("demo: FDR %.3f <= 0.10 and power %.2f >= 0.70 against the planted truth", dfdr, dpow))
check(all(sign(m$log2FoldChange[m$significant & m$is_de]) == sign(m$true_log2fc[m$significant & m$is_de])),
      "demo: every true call has the right direction (drought vs control, not reversed)")
sm <- utils::read.csv(file.path(out, "tables", "de_summary.csv"))
check(sm$n_tested == sum(!is.na(dd$padj)), "summary n_tested = genes with non-missing padj (enrichment universe)")
vst <- utils::read.csv(file.path(out, "tables", "vst.csv"), check.names = FALSE)
check(ncol(vst) == 9L && all(is.finite(as.matrix(vst[, -1]))), "VST table: gene_id plus 8 finite sample columns")

# 9. Command-line refusals (these stop before DESeq2 is fitted) ---------------------
cnt <- as.matrix(utils::read.csv(file.path(SKILL_TEST_DIR, "demo", "counts.csv"), row.names = 1, check.names = FALSE))
tpm <- t(t(cnt / 1500) / colSums(cnt / 1500)) * 1e6
tpm_file <- file.path(tmp, "tpm.csv")
utils::write.csv(data.frame(gene_id = rownames(tpm), tpm, check.names = FALSE), tpm_file, row.names = FALSE)
demo_samples <- file.path(SKILL_TEST_DIR, "demo", "samples.csv")
bad <- run_cli(c("--counts", shQuote(tpm_file), "--samples", shQuote(demo_samples), "--design", shQuote("~ condition"),
                 "--reference", "condition=control", "--out", shQuote(file.path(tmp, "bad"))))
check(bad$status != 0L && any(grepl("TPM", bad$output)), "CLI refuses TPM input")
conf <- utils::read.csv(demo_samples)
conf$batch <- ifelse(conf$condition == "control", "b1", "b2")
conf_file <- file.path(tmp, "confounded.csv")
utils::write.csv(conf, conf_file, row.names = FALSE)
bad2 <- run_cli(c("--counts", shQuote(file.path(SKILL_TEST_DIR, "demo", "counts.csv")), "--samples", shQuote(conf_file),
                  "--design", shQuote("~ batch + condition"), "--reference", "condition=control",
                  "--out", shQuote(file.path(tmp, "bad2"))))
check(bad2$status != 0L && any(grepl("completely confounded", bad2$output)), "CLI refuses batch confounded with condition")
bad3 <- run_cli(c("--counts", shQuote(file.path(SKILL_TEST_DIR, "demo", "counts.csv")), "--samples", shQuote(demo_samples),
                  "--design", shQuote("~ batch + condition"), "--out", shQuote(file.path(tmp, "bad3"))))
check(bad3$status != 0L && any(grepl("reference", bad3$output)), "CLI asks for an explicit reference level")
check(run_cli(c("--demo", "--color", "red"))$status != 0L, "unknown option exits with an error")

cat(sprintf("PASS: %d checks (%s, part 2)\n", n_checks, SKILL))
