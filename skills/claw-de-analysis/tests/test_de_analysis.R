# Tests for claw-de-analysis, part 1: Wald tests. Run: Rscript tests/test_de_analysis.R
# Known answers from negative-binomial simulations with planted DE genes: empirical
# FDR and power over several seeds, apeglm accuracy, contrast direction; refusal of
# TPM/FPKM, unreplicated groups and confounded designs. Part 2 (LRT time course,
# tximport offsets, command line) is tests/test_de_lrt_cli.R.

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
T0 <- Sys.time()
stamp <- function(what) cat(sprintf("  [%5.1fs] %s\n", as.numeric(difftime(Sys.time(), T0, units = "secs")), what))
tmp <- file.path(tempdir(), "claw-de-analysis-test")
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

stamp("1. Empirical FDR")
# 1. Empirical FDR and power over 5 seeds (Wald test, alpha 0.05) -------------------
fdr <- power <- numeric()
for (seed in 1:5) {
  x <- sim(seed)
  r <- run_de(x$counts, prep(x), design, shrink = "none", reference = ref)
  tab <- r$tables[[1]]
  if (seed == 5) rn <- tab                      # reused below: same data without shrinkage
  call <- stats::setNames(tab$significant, tab$gene_id)[rownames(x$counts)]
  fdr[seed] <- if (any(call)) mean(!x$de[call]) else 0
  power[seed] <- mean(call[x$de])
}
cat(sprintf("  FDR by seed: %s (mean %.3f); power: %s (mean %.3f)\n",
            paste(sprintf("%.3f", fdr), collapse = " "), mean(fdr),
            paste(sprintf("%.2f", power), collapse = " "), mean(power)))
check(mean(fdr) <= 0.05 + 0.03, sprintf("mean empirical FDR %.3f <= nominal 0.05 + 0.03 tolerance", mean(fdr)))
check(mean(fdr) >= 0.015, sprintf("mean empirical FDR %.3f is not degenerate (>= 0.015; about 0.065 in a 15-seed calibration)", mean(fdr)))
check(mean(power) >= 0.6, sprintf("mean power %.2f >= 0.6 for 2- to 4-fold changes, n = 4 per group", mean(power)))

stamp("2. apeglm shrinkage")
# 2. apeglm shrinkage: accurate for strong genes, pulls noise to 0, does not change tests
x <- sim(5)
s1 <- prep(x)
ra <- run_de(x$counts, s1, design, shrink = "apeglm", reference = ref)$tables[[1]]
ra <- ra[match(rownames(x$counts), ra$gene_id), ]
rn <- rn[match(rownames(x$counts), rn$gene_id), ]
check(isTRUE(all.equal(ra$padj, rn$padj)), "shrinkage leaves the Wald p-values and padj unchanged")
strong <- x$de & !is.na(ra$baseMean) & ra$baseMean > 100
err <- stats::median(abs(ra$log2FoldChange_shrunk[strong] - x$true[strong]))
check(err < 0.25, sprintf("apeglm LFC within 0.25 of truth for DE genes with baseMean > 100 (median error %.3f)", err))
nul <- !x$de & !is.na(ra$log2FoldChange_shrunk)
check(mean(abs(ra$log2FoldChange_shrunk[nul]) <= abs(ra$log2FoldChange[nul]) + 1e-6) > 0.95,
      "apeglm shrinks null-gene LFCs toward 0 (|shrunken| <= |MLE| for > 95%)")
err_mle <- stats::median(abs(ra$log2FoldChange[strong] - x$true[strong]))
check(err <= err_mle + 0.02, sprintf("apeglm preserves large effects (median error %.3f vs MLE %.3f)", err, err_mle))
shrink_ratio <- function(b) stats::median(abs(ra$log2FoldChange_shrunk[b])) / stats::median(abs(ra$log2FoldChange[b]))
noisy <- nul & ra$lfcSE >= 0.5
precise <- nul & ra$lfcSE < 0.25
check(sum(noisy) >= 10 && shrink_ratio(noisy) < 0.5 && shrink_ratio(noisy) < shrink_ratio(precise),
      sprintf(paste0("apeglm shrinks noisy null estimates (lfcSE >= 0.5) by more than half (ratio %.2f) and ",
                     "more than precise ones (lfcSE < 0.25, ratio %.2f)"), shrink_ratio(noisy), shrink_ratio(precise)))

stamp("3. Contrast direction")
# 3. Contrast direction: reversed contrast flips signs; apeglm via relevel agrees
rr <- run_de(x$counts, s1, design, contrasts = parse_contrasts("condition,control,treated"),
             shrink = "apeglm", reference = ref)$tables[[1]]
rr <- rr[match(rownames(x$counts), rr$gene_id), ]
ok <- !is.na(rr$log2FoldChange) & !is.na(ra$log2FoldChange)
check(max(abs(rr$log2FoldChange[ok] + ra$log2FoldChange[ok])) < 1e-6 &&
        isTRUE(all.equal(rr$pvalue, ra$pvalue)), "reversed contrast: MLE LFC negated, same p-values")
oks <- ok & !is.na(rr$log2FoldChange_shrunk)
check(stats::cor(rr$log2FoldChange_shrunk[oks], -ra$log2FoldChange_shrunk[oks]) > 0.999,
      "apeglm on a non-reference denominator (relevel + nbinomWaldTest) mirrors the forward contrast")

stamp("4. Inputs that are not counts")
# 4. Inputs that are not counts are refused --------------------------------------
cnt <- x$counts
len <- matrix(stats::runif(nrow(cnt), 500, 3000), nrow(cnt), ncol(cnt), dimnames = dimnames(cnt))
tpm <- t(t(cnt / len) / colSums(cnt / len)) * 1e6
fpkm <- t(t(cnt / (len / 1e3)) / (colSums(cnt) / 1e6))
expect_error(check_count_matrix(tpm), "TPM", "TPM input rejected")
expect_error(check_count_matrix(fpkm), "FPKM", "FPKM input rejected")
expect_error(check_count_matrix(log2(cnt + 0.5)), "negative|not whole numbers", "log-scale input rejected")
cna <- cnt; cna[1, 1] <- NA
expect_error(check_count_matrix(cna), "NA", "missing counts rejected")
check(length(check_count_matrix(cnt)) == 0L, "integer counts accepted")
rare <- round(t(t(cnt) / colSums(cnt)) * 1e6)
check(grepl("rarefied", check_count_matrix(rare)), "integer counts summing to 1e6 (rarefied) warn instead of failing as TPM")
check(length(check_count_matrix(cnt + 0.37, noninteger_ok = TRUE)) == 0L,
      "non-integer tximport estimated counts accepted when length offsets are supplied")

stamp("5. Designs that cannot")
# 5. Designs that cannot be tested are refused ------------------------------------
s_unrep <- data.frame(sample = colnames(cnt)[1:4], condition = c("control", "control", "control", "treated"),
                      row.names = colnames(cnt)[1:4])
s_unrep <- prepare_design(s_unrep, ~ condition, ref)$samples
expect_error(run_de(cnt[, 1:4], s_unrep, ~ condition, shrink = "none", reference = ref),
             "No biological replicates", "group with one sample (no replicate) rejected")
s_conf <- data.frame(sample = colnames(cnt), condition = rep(c("control", "treated"), each = 4),
                     batch = rep(c("b1", "b2"), each = 4), row.names = colnames(cnt))
s_conf <- prepare_design(s_conf, design, ref)$samples
expect_error(run_de(cnt, s_conf, design, shrink = "none", reference = ref), "completely confounded",
             "batch completely confounded with condition rejected")
s_nest <- data.frame(sample = colnames(cnt), genotype = rep(c("WT", "mut"), each = 4),
                     plant = rep(c("p1", "p2", "p3", "p4"), each = 2), row.names = colnames(cnt))
s_nest <- prepare_design(s_nest, ~ plant + genotype, list(genotype = "WT"))$samples
expect_error(run_de(cnt, s_nest, ~ plant + genotype, shrink = "none", reference = list(genotype = "WT")),
             "nested", "plant nested within genotype detected")
expect_error(default_contrasts(s1, design, list()), "reference", "missing reference level refused")
small <- cnt[1:400, ]                       # warning-only checks: a gene subset is enough
s_three <- data.frame(sample = colnames(cnt), condition = c(rep("control", 3), rep("treated", 4), "heat"),
                      row.names = colnames(cnt))
s_three <- prepare_design(s_three, ~ condition, ref)$samples
r3 <- run_de(small, s_three, ~ condition, contrasts = parse_contrasts("condition,treated,control"),
             shrink = "none", reference = ref)
check(nrow(r3$tables[[1]]) == nrow(small) && any(grepl("single sample", r3$warn)),
      "an unreplicated level outside the contrast gives a warning, not a refusal")
expect_error(run_de(small, s_three, ~ condition, contrasts = parse_contrasts("condition,heat,control"),
                    shrink = "none", reference = ref), "No biological replicates",
             "comparing the unreplicated level itself is refused")
s_part <- data.frame(sample = colnames(cnt), condition = rep(c("control", "treated"), each = 4),
                     batch = c("b1", "b1", "b2", "b2", "b2", "b2", "b2", "b2"), row.names = colnames(cnt))
s_part <- prepare_design(s_part, design, ref)$samples
rp <- run_de(small, s_part, design, shrink = "none", reference = ref)
check(any(grepl("partly confounded", rp$warn)), "partial batch-condition confounding gives a warning")
s_two <- data.frame(sample = colnames(cnt)[c(1, 2, 5, 6)], condition = rep(c("control", "treated"), each = 2),
                    row.names = colnames(cnt)[c(1, 2, 5, 6)])
s_two <- prepare_design(s_two, ~ condition, ref)$samples
r2 <- run_de(small[, c(1, 2, 5, 6)], s_two, ~ condition, shrink = "none", reference = ref)
check(any(grepl("Schurch", r2$warn)), "fewer than 3 replicates gives the replicate warning")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
