# Tests for claw-gwas-crop. Run from anywhere: Rscript tests/test_gwas_crop.R
# Known answers: closed-form power, genomic-control lambda on simulated
# p-values, hand-worked QTL-region grouping, and a simulated structured panel
# of inbred lines with known QTNs and a trait with no marker effects.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "gwas_crop.R"))  # loads functions; main() does not run
source(file.path(SKILL_TEST_DIR, "demo", "make_demo.R"))     # simulate_panel(), simulate_traits()

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "gwas_crop.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
tmp <- file.path(tempdir(), "claw-gwas-crop-test")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE)

# 1. Power: non-centrality for 80% power at p = 5e-8 equals (z_{1-a/2} + z_{0.8})^2 -------
r2 <- min_detectable_r2(1000, 5e-8, 0.8)
ncp <- 1000 * r2 / (1 - r2)
ncp_closed <- (stats::qnorm(1 - 5e-8 / 2) + stats::qnorm(0.8))^2
check(abs(ncp - ncp_closed) / ncp_closed < 1e-3, sprintf("power NCP %.3f matches closed form %.3f", ncp, ncp_closed))

# 2. Genomic-control lambda ------------------------------------------------------------------
set.seed(1)
check(abs(lambda_gc(stats::runif(1e5)) - 1) < 0.02, "lambda of uniform p-values is 1")
p_infl <- stats::pchisq(1.5 * stats::rchisq(1e5, 1), 1, lower.tail = FALSE)
check(abs(lambda_gc(p_infl) - 1.5) < 0.05, "lambda recovers a 1.5-fold inflation of chi-square statistics")

# 3. QTL-region grouping (hand-worked) ----------------------------------------------------------
sig <- data.frame(snp = c("a", "b", "c", "d"), chr = c("1", "1", "1", "2"),
                  pos = c(100, 150, 400, 5000), p_value = c(1e-6, 1e-8, 1e-7, 1e-9))
reg <- group_regions(sig, c("1" = 100, "2" = 100))
check(nrow(reg) == 3L, "three regions: SNPs 50 bp apart join, 250 bp apart do not")
check(identical(reg$start_bp, c(1, 300, 4900)) && identical(reg$end_bp, c(250, 500, 5100)),
      "regions extend by the LD distance and start at bp 1 at the chromosome end")
check(identical(reg$lead_snp, c("b", "c", "d")) && identical(reg$n_significant_snps, c(2L, 1L, 1L)),
      "lead SNP is the smallest p-value in each region")
check(nrow(group_regions(sig[0, ], c("1" = 100))) == 0L, "no significant SNPs give an empty region table")
f_ld <- file.path(tmp, "ld_window.csv")
utils::write.csv(data.frame(chr = c("all", "2"), ld_distance_bp = c(500, 80)), f_ld, row.names = FALSE)
ld <- ld_lookup(c("1", "2", "3"), ld_table = f_ld)
check(identical(unname(ld), c(500, 80, 500)), "per-chromosome LD table with an 'all' default row")
f_ld_na <- file.path(tmp, "ld_window_na.csv")
utils::write.csv(data.frame(chr = c("all", "2"), ld_distance_bp = c(500, NA)), f_ld_na, row.names = FALSE)
check(identical(unname(ld_lookup(c("1", "2"), ld_table = f_ld_na)), c(500, 500)),
      "a chromosome whose LD window is NA (not estimable) falls back to the 'all' row")
utils::write.csv(data.frame(chr = c("all", "2"), ld_distance_bp = c(NA, NA)), f_ld_na, row.names = FALSE)
err <- tryCatch(ld_lookup(c("1", "2"), ld_table = f_ld_na), error = conditionMessage)
check(is.character(err) && grepl("no usable distance", err) && grepl("--ld-distance", err),
      "an LD table without any usable distance stops and asks for --ld-distance")

# 4. Structured panel with known QTNs (seed 5; see SKILL.md for a 12-seed summary) ------------
panel <- simulate_panel(seed = 5)
tr <- simulate_traits(panel, seed = 5)
X <- panel$G; storage.mode(X) <- "double"
ids <- rownames(X); map <- panel$map; n <- nrow(X)
ph <- data.frame(genotype = ids, y_qtn = tr$y_qtn, y_null = tr$y_null, stringsAsFactors = FALSE)
pcs <- marker_pcs(X, 10)
check(pcs$var_explained[1] > 0.05 && pcs$var_explained[2] > 0.03 && pcs$var_explained[3] < pcs$var_explained[2] / 2,
      "PCA shows the simulated three-subpopulation structure (two large PCs)")
gd <- build_gdata(X, map, ph, "genotype", c("y_qtn", "y_null"), data.frame(pcs$scores[, 1:2], row.names = ids))
gK <- run_statgen(gd, c("y_qtn", "y_null"), character(), "astle", FALSE, 0.05, 0.05, 1)
rq <- tidy_results(gK, "y_qtn", "bonferroni", 0.05)
regq <- group_regions(rq[rq$significant, c("snp", "chr", "pos", "p_value")], stats::setNames(rep(2e6, 5), as.character(1:5)))
qchr <- map$chr[match(tr$qtn, map$marker)]; qpos <- map$pos[match(tr$qtn, map$marker)]
in_reg <- vapply(seq_along(qpos), function(i) any(regq$chr == qchr[i] & regq$start_bp <= qpos[i] & regq$end_bp >= qpos[i]), NA)
check(all(in_reg), "kinship model: every QTN lies in a Bonferroni-significant region")
extra <- sum(!vapply(seq_len(nrow(regq)), function(r) any(regq$chr[r] == qchr & regq$start_bp[r] <= qpos & regq$end_bp[r] >= qpos), NA))
check(extra == 0L, "kinship model: no significant region away from the QTNs")
true_eff <- sqrt(0.10) / apply(X[, tr$qtn], 2, stats::sd)   # trait units per ALT-allele copy
est <- rq$effect_per_alt_allele[match(tr$qtn, rq$snp)]
se <- rq$effect_se[match(tr$qtn, rq$snp)]
check(all(abs(est - true_eff) < 3 * se), "QTN effects per ALT allele are within 3 SE of the simulated values")
check(all(est > 0), "effect sign follows the ALT (counted) allele")

lamK <- gK$GWASInfo$inflationFactor[[1]][["y_null"]]
rn <- tidy_results(gK, "y_null", "bonferroni", 0.05)
check(lamK > 0.9 && lamK < 1.1, sprintf("null trait, kinship model: lambda %.3f is within 0.9-1.1", lamK))
check(abs(lambda_gc(rn$p_value) - lamK) / lamK < 0.03, "statgenGWAS lambda agrees with the chi-square median lambda")
check(sum(rn$significant) == 0L, "null trait: no Bonferroni-significant SNP")
check(sum(tidy_results(gK, "y_null", "fdr", 0.05)$significant) == 0L, "null trait: no BH-FDR-significant SNP")
K0 <- diag(n); dimnames(K0) <- list(ids, ids)
g0 <- run_statgen(gd, "y_null", character(), "astle", FALSE, 0.05, 0.05, 1, kin = K0)
lam0 <- g0$GWASInfo$inflationFactor[[1]][["y_null"]]
r0 <- tidy_results(g0, "y_null", "bonferroni", 0.05)
check(lam0 > 2 && lam0 > lamK + 1, sprintf("null trait without kinship is visibly inflated (lambda %.2f vs %.2f)", lam0, lamK))
check(sum(r0$significant) > 0L, "null trait without kinship gives false-positive hits")
lamL <- run_statgen(gd, "y_null", character(), "astle", TRUE, 0.05, 0.05, 1)$GWASInfo$inflationFactor[[1]][["y_null"]]
lamLP <- run_statgen(gd, "y_null", c("PC1", "PC2"), "astle", TRUE, 0.05, 0.05, 1)$GWASInfo$inflationFactor[[1]][["y_null"]]
check(lamL > 1.1 && lamLP < 1.1 && lamLP > 0.9,
      sprintf("LOCO alone is inflated in this structured panel (%.2f); LOCO + 2 PCs is calibrated (%.2f)", lamL, lamLP))

# 5. Input checks --------------------------------------------------------------------------------
f_dup <- file.path(tmp, "dup.csv")
utils::write.csv(data.frame(genotype = c("L001", "L001", "L002"), y = 1:3), f_dup, row.names = FALSE)
err <- tryCatch(read_pheno(f_dup, "genotype", "y", character()), error = conditionMessage)
check(is.character(err) && grepl("BLUEs", err) && grepl("claw-trial-analysis", err),
      "repeated phenotype IDs stop with advice to compute BLUEs")
err <- tryCatch(read_pheno(f_dup, "genotype", character(), character()), error = conditionMessage)
check(is.character(err) && grepl("--traits", err), "missing --traits lists the available columns")
vcf <- c("##fileformat=VCFv4.2",
         paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "A", "B", sep = "\t"),
         paste("1", "100", "m1", "A", "G", ".", "PASS", ".", "GT", "0/0", "1|1", sep = "\t"),
         paste("1", "900", ".", "C", "T", ".", "PASS", ".", "GT", "./.", "0/1", sep = "\t"),
         paste("2", "50", "m3", "G", "A,C", ".", "PASS", ".", "GT", "0/0", "1/2", sep = "\t"))
f_vcf <- file.path(tmp, "tiny.vcf"); writeLines(vcf, f_vcf)
v <- read_genotypes(f_vcf)
check(identical(unname(v$G), matrix(c(0, NA, 2, 1), 2)) && identical(v$map$marker, c("m1", "1_900")),
      "VCF input: ALT counts, CHROM_POS IDs, multi-allelic site excluded")
im <- impute_mean(matrix(c(0, 2, NA, 2, NA, NA, 0, 0, 2), 3, byrow = TRUE), max_missing = 0.5)
check(identical(im$keep, c(TRUE, FALSE, TRUE)) && im$G[1, 3] == 1 && im$n_imputed == 1L,
      "mean imputation fills 2p and drops markers above --max-missing")

# 6. End to end on the demo data ----------------------------------------------------------------
out <- file.path(tmp, "demo")
res <- run_cli(c("--demo", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "tables/gwas_height_cm.csv", "tables/gwas_protein_pct.csv", "tables/diagnostics.csv",
            "tables/significant_snps.csv", "tables/qtl_regions.csv", "tables/principal_components.csv",
            "figures/manhattan_height_cm.png", "figures/qq_height_cm.png", "figures/manhattan_protein_pct.png",
            "figures/qq_protein_pct.png", "figures/pca_structure.png",
            "reproducibility/commands.sh", "reproducibility/sessionInfo.txt")) {
  check(file.exists(file.path(out, f)), paste("demo writes", f))
}
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
truth <- utils::read.csv(file.path(SKILL_TEST_DIR, "demo", "demo_truth.csv"), stringsAsFactors = FALSE)
qr <- utils::read.csv(file.path(out, "tables", "qtl_regions.csv"), stringsAsFactors = FALSE)
dg <- utils::read.csv(file.path(out, "tables", "diagnostics.csv"), stringsAsFactors = FALSE)
check(all(vapply(seq_len(nrow(truth)), function(i) any(qr$trait == "height_cm" & qr$chr == truth$chr[i] &
                                                         qr$start_bp <= truth$pos[i] & qr$end_bp >= truth$pos[i]), NA)),
      "demo: each of the 3 height QTNs is inside a reported region")
check(sum(qr$trait == "protein_pct") == 0L, "demo: no region for the trait without marker effects")
lp <- dg$lambda_gc[dg$trait == "protein_pct"]
check(lp > 0.9 && lp < 1.1 && dg$lambda_naive[dg$trait == "protein_pct"] > 2,
      "demo: lambda near 1 with kinship, inflated without")
rep_txt <- readLines(file.path(out, "report.md"))
check(!any(grepl("^## Warnings", rep_txt)) || !any(grepl("lambda = .*> 1.1", rep_txt)), "demo report has no inflation warning")

out2 <- file.path(tmp, "demo_fdr")
res <- run_cli(c("--demo", "--threshold", "fdr", "--ld-table", shQuote(f_ld), "--n-pcs", "2", "--loco",
                 "--out", shQuote(out2)))
check(res$status == 0L, paste("FDR + LD table + PCs + LOCO run exits cleanly:", paste(res$output, collapse = "\n")))
qr2 <- utils::read.csv(file.path(out2, "tables", "qtl_regions.csv"), stringsAsFactors = FALSE)
check(all(qr2$ld_distance_bp == ifelse(qr2$chr == "2", 80, 500)) && nrow(qr2) >= 3L,
      "LD distances come from --ld-table (chromosome 2 row, 'all' row for the rest)")
g2 <- utils::read.csv(file.path(out2, "tables", "gwas_height_cm.csv"), stringsAsFactors = FALSE)
check(all(g2$significant == (!is.na(g2$q_value_bh) & g2$q_value_bh <= 0.05)), "FDR significance is BH q <= alpha")

bad <- run_cli(c("--demo", "--threshold", "lenient"))
check(bad$status != 0L, "invalid --threshold exits with an error")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
