# Tests for claw-marker-qc. Run from anywhere: Rscript tests/test_marker_qc.R
# Known-answer checks: hand-coded genotype conversions, identity-by-state
# against brute force, chi-square against stats::chisq.test, and simulated
# panels with planted problems that must all be detected exactly (no misses,
# no false flags). Ends with VCF round trips and command-line runs.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "marker_qc.R"))  # loads functions; main() does not run
source(file.path(SKILL_TEST_DIR, "demo", "make_demo.R"))     # simulate_qc_panel()

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
same_set <- function(a, b) setequal(as.character(a), as.character(b)) && length(unique(a)) == length(unique(b))

run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "marker_qc.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
tmp <- file.path(tempdir(), "claw-marker-qc-test")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE)
opts <- function(...) { o <- DEFAULTS; o$germplasm <- "inbred"; m <- list(...); o[names(m)] <- m; o }

# 1. Genotype conversions (hand-coded answers) ------------------------------------
gt <- matrix(c("0/0", "0|1", "1/1", "./.", "1|0", ".", "./1", "0/0"), 2, 4)
d <- gt_to_dosage(gt)
check(identical(as.vector(d), c(0L, 1L, 2L, NA, 1L, NA, NA, 0L)), "VCF GT strings -> ALT counts")
check(attr(d, "n_half_missing") == 1L, "half-missing call ./1 is counted")
err <- tryCatch(gt_to_dosage(matrix(c("0", "1"), 1)), error = conditionMessage)
check(is.character(err) && grepl("diploid", err), "haploid GT is rejected with a clear message")
err <- tryCatch(gt_to_dosage(matrix(c("0/2", "0/0"), 1)), error = conditionMessage)
check(is.character(err) && grepl("Unexpected GT", err), "allele index 2 at a biallelic site is rejected")

calls <- matrix(c("GG", "AG", "AA", "NN", "R", "G", "A", "N", "AT"), 1)
h <- hapmap_dosage(calls, "G", "A")   # alleles column G/A: counts copies of A
check(identical(as.vector(h), c(0L, 1L, 2L, NA, 1L, 0L, 2L, NA, NA)), "HapMap two-letter and IUPAC calls")
check(attr(h, "n_bad") == 1L, "a call with an allele outside the alleles column is counted as bad")

cv <- codes_to_dosage(matrix(c("-1", "0", "1", "NA"), 2), character())
check(identical(as.vector(cv$geno), c(0L, 1L, 2L, NA)) && length(cv$notes) == 1L, "-1/0/1 coding converted with a note")
err <- tryCatch(codes_to_dosage(matrix(c("0", "0.37"), 1), character()), error = conditionMessage)
check(is.character(err) && grepl("hard calls", err), "imputed dosages are rejected")

# 2. Identity-by-state equals brute force -----------------------------------------------
set.seed(1)
X <- matrix(sample(c(0L, 1L, 2L, NA), 300, TRUE, prob = c(0.4, 0.1, 0.4, 0.1)), 10, 30)
rownames(X) <- paste0("s", 1:10)
ib <- ibs_similarity(X)
bf <- outer(1:10, 1:10, Vectorize(function(i, j) {
  ok <- !is.na(X[i, ]) & !is.na(X[j, ]); 1 - mean(abs(X[i, ok] - X[j, ok])) / 2 }))
nb <- outer(1:10, 1:10, Vectorize(function(i, j) sum(!is.na(X[i, ]) & !is.na(X[j, ]))))
check(max(abs(ib$ibs - bf)) < 1e-12, "IBS matrix equals brute-force computation")
check(all(ib$n_shared == nb), "shared-marker counts equal brute force")

# 3. Chi-square segregation statistic equals stats::chisq.test ----------------------------
Gs <- matrix(c(rep(0L, 70), rep(2L, 30), rep(0L, 25), rep(1L, 50), rep(2L, 25)), 2, byrow = TRUE)
st_ril <- segregation_test(Gs[1, , drop = FALSE], "ril")
ref <- suppressWarnings(stats::chisq.test(c(70, 30), p = c(0.5, 0.5)))
check(abs(st_ril$seg_chisq - unname(ref$statistic)) < 1e-10 && abs(st_ril$seg_p - ref$p.value) < 1e-12,
      "RIL 1:1 chi-square matches chisq.test")
st_f2 <- segregation_test(Gs[2, , drop = FALSE], "f2")
check(abs(st_f2$seg_chisq) < 1e-12 && st_f2$seg_df == 2L, "F2 1:2:1 counts 25:50:25 give chi-square 0 on 2 df")
st_bc <- segregation_test(matrix(c(rep(2L, 60), rep(1L, 40)), 1), "bc1", p1 = 2L, p2 = 0L)
ref <- stats::chisq.test(c(60, 40), p = c(0.5, 0.5))
check(abs(st_bc$seg_chisq - unname(ref$statistic)) < 1e-10, "BC1 1:1 chi-square uses the recurrent-parent class")

# 4. Planted problems in a simulated inbred panel are detected exactly -------------------
for (seed in c(7, 99)) {
  sim <- simulate_qc_panel(seed = seed)
  r <- run_qc(sim$geno, sim$map, opts())
  tr <- split(sim$truth$id, sim$truth$problem)
  ms <- r$markers; ss <- r$samples
  check(same_set(ms$marker[ms$fail_call_rate], tr$marker_call_rate), paste("seed", seed, "high-missing markers exactly"))
  check(same_set(ms$marker[ms$monomorphic & !ms$fail_call_rate], tr$marker_monomorphic), paste("seed", seed, "monomorphic markers exactly"))
  check(same_set(ms$marker[ms$fail_maf & !ms$monomorphic & !ms$fail_call_rate], tr$marker_low_maf), paste("seed", seed, "low-MAF markers exactly"))
  check(same_set(ss$sample[ss$fail_call_rate], tr$sample_call_rate), paste("seed", seed, "low-call-rate samples exactly"))
  check(same_set(ss$sample[ss$het_flag], tr$sample_heterozygosity), paste("seed", seed, "contaminated samples exactly"))
  pr <- apply(r$pairs[, c("sample_1", "sample_2")], 1, function(x) paste(sort(x), collapse = "~"))
  check(same_set(pr, tr$duplicate_pair), paste("seed", seed, "duplicate pairs exactly"))
  check(all(r$pairs$ibs > 0.98) && max(r$ibs[upper.tri(r$ibs)][!r$ibs[upper.tri(r$ibs)] %in% r$pairs$ibs], na.rm = TRUE) < 0.9,
        paste("seed", seed, "duplicates are well separated from other pairs"))
}
r2 <- run_qc(sim$geno, sim$map, opts(drop_duplicates = TRUE, drop_het_outliers = TRUE))
check(sum(!r2$samples$keep) == 4L + 3L + 3L, "drop options remove low-call, contaminated and one sample per duplicate pair")

# 5. Outbred germplasm: heterozygosity outliers by mean +/- 3 SD -------------------------
set.seed(3)
m <- 2000; n <- 100
p <- stats::runif(m, 0.2, 0.8)
Go <- matrix(stats::rbinom(m * n, 2, p), m, n, dimnames = list(sprintf("S1_%d", 1:m * 100), sprintf("I%03d", 1:n)))
mix <- 1:3; selfed <- 4:5
for (j in mix) { b <- stats::rbinom(m, 2, p); Go[, j] <- ifelse(Go[, j] == b & Go[, j] != 1L, Go[, j], 1L) }
for (j in selfed) { h <- Go[, j] == 1L; Go[h, j] <- sample(c(0L, 1L, 1L, 2L), sum(h), TRUE) }   # one selfing generation halves heterozygosity
mapo <- make_map(rownames(Go), "1", 1:m * 100)
ro <- run_qc(Go, mapo, opts(germplasm = "outbred"))
check(same_set(ro$samples$sample[ro$samples$het_flag], colnames(Go)[c(mix, selfed)]),
      "outbred: excess (mixed) and deficit (selfed) heterozygosity flagged exactly")

# 6. RIL population with parents: distortion and non-parental lines ----------------------
set.seed(11)
m <- 400; n <- 200
poly <- 1:300; mono <- 301:400
p1 <- c(sample(c(0L, 2L), 300, TRUE), sample(c(0L, 2L), 100, TRUE))
p2 <- c(2L - p1[poly], p1[mono])
q1 <- rep(0.5, m); dist <- sample(poly, 10); q1[dist] <- 0.8      # probability of the P1 allele
Gr <- matrix(0L, m, n)
for (j in 1:n) {
  fromP1 <- stats::runif(m) < q1
  Gr[, j] <- ifelse(fromP1, p1, p2)
}
het <- matrix(stats::runif(m * n) < 0.01, m, n); het[mono, ] <- FALSE; Gr[het] <- 1L
offtype <- 1:3
for (j in offtype) Gr[, j] <- 2L * stats::rbinom(m, 1, 0.5)
Gr[matrix(stats::runif(m * n) < 0.01, m, n)] <- NA
Gr <- cbind(Gr, p1, p2)
dimnames(Gr) <- list(sprintf("S%d_%d", rep(1:4, each = 100), rep(1:100, 4) * 1000), c(sprintf("RIL%03d", 1:n), "P1", "P2"))
mapr <- make_map(rownames(Gr), rep(1:4, each = 100), rep(1:100, 4) * 1000)
rr <- run_qc(Gr, mapr, opts(population = "ril", parents = "P1,P2", min_maf = 0))
check(same_set(rr$markers$marker[rr$markers$distorted], rownames(Gr)[dist]), "RIL: the 10 distorted markers are found exactly")
check(!any(rr$markers$seg_tested[mono]), "RIL: markers monomorphic between parents are not tested")
check(same_set(rr$samples$sample[rr$samples$nonparental_flag], sprintf("RIL%03d", offtype)),
      "RIL: off-type lines with non-parental alleles are flagged exactly")
check(rr$n_parental_mono == 100L, "RIL: 100 markers monomorphic between parents used for the parent check")

# 7. F2 without parents: 1:2:1 distortion ------------------------------------------------
set.seed(5)
m <- 200; n <- 300
Gf <- matrix(sample(0:2, m * n, TRUE, prob = c(0.25, 0.5, 0.25)), m, n)
dist <- sample.int(m, 5)
for (i in dist) Gf[i, ] <- sample(0:2, n, TRUE, prob = c(0.45, 0.40, 0.15))
dimnames(Gf) <- list(sprintf("S1_%d", 1:m * 1000), sprintf("F2_%03d", 1:n))
rf <- run_qc(Gf, make_map(rownames(Gf), "1", 1:m * 1000), opts(germplasm = "outbred", population = "f2"))
check(same_set(rf$markers$marker[rf$markers$distorted], rownames(Gf)[dist]), "F2: the 5 distorted markers are found exactly")

# 8. Orientation detection for 0/1/2 matrices --------------------------------------------
Gsmall <- sim$geno[1:30, 1:12]
Gsmall[is.na(Gsmall)] <- 0L
f_rows <- file.path(tmp, "rows.csv"); f_cols <- file.path(tmp, "cols.csv")
utils::write.csv(data.frame(id = rownames(Gsmall), Gsmall, check.names = FALSE), f_rows, row.names = FALSE)
utils::write.csv(data.frame(id = colnames(Gsmall), t(Gsmall), check.names = FALSE), f_cols, row.names = FALSE)
a <- read_genotypes(f_rows); b <- read_genotypes(f_cols)
check(identical(a$geno, b$geno) && identical(dimnames(a$geno), dimnames(Gsmall)),
      "markers-in-rows and samples-in-rows files give the same matrix")
check(grepl("look like marker names", a$orientation), "orientation detected from marker-like IDs")
anon <- Gsmall; dimnames(anon) <- list(paste0("m", 1:30), paste0("line", 1:12))
f_anon <- file.path(tmp, "anon.csv")
utils::write.csv(data.frame(id = rownames(anon), anon, check.names = FALSE), f_anon, row.names = FALSE)
err <- tryCatch(read_genotypes(f_anon), error = conditionMessage)
check(is.character(err) && grepl("--orientation", err), "ambiguous orientation stops and asks for --orientation")
f_map <- file.path(tmp, "map.csv")
utils::write.csv(data.frame(marker = rownames(anon), chr = 1, pos = 1:30), f_map, row.names = FALSE)
c2 <- read_genotypes(f_anon, map_path = f_map)
check(identical(unname(c2$geno), unname(anon)) && all(c2$map$chr == "1"), "orientation detected from --map IDs")
f_emb <- file.path(tmp, "embedded.csv")
write_marker_table(Gsmall, sim$map[1:30, ], f_emb)
e <- read_genotypes(f_emb)
check(identical(e$geno, Gsmall) && identical(e$map$pos, as.numeric(sim$map$pos[1:30])),
      "AgriClaw marker table (map columns embedded) round-trips")

# 9. VCF round trips ---------------------------------------------------------------------
vcf_txt <- c("##fileformat=VCFv4.2",
             "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
             paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "A", "B", "C", sep = "\t"),
             paste("1", "100", "snp1", "A", "G", ".", "PASS", ".", "GT", "0/0", "0|1", "1/1", sep = "\t"),
             paste("1", "200", ".", "C", "T", ".", "PASS", ".", "GT", "./.", "1|1", "0/1", sep = "\t"),
             paste("1", "300", "snp3", "G", "A,C", ".", "PASS", ".", "GT", "0/1", "1/2", "0/0", sep = "\t"),
             paste("2", "1500000", "snp4", "T", "C", ".", "PASS", ".", "GT:DP", "1/1:12", "./1:3", "0/0:20", sep = "\t"))
f_vcf <- file.path(tmp, "tiny.vcf")
writeLines(vcf_txt, f_vcf)
v <- read_genotypes(f_vcf)
check(identical(unname(v$geno), matrix(c(0L, NA, 2L, 1L, 2L, NA, 2L, 1L, 0L), 3)), "VCF genotypes read as ALT counts")
check(identical(rownames(v$geno), c("snp1", "1_200", "snp4")), "missing VCF IDs become CHROM_POS")
check(any(grepl("multi-allelic", v$notes)) && any(grepl("half-missing", v$notes)) && any(grepl("Phased", v$notes)),
      "multi-allelic, half-missing and phased calls are reported, not silently changed")
for (ext in c(".vcf", ".vcf.gz")) {
  f_out <- file.path(tmp, paste0("roundtrip", ext))
  write_vcf(v$geno, v$map, f_out)
  w <- read_genotypes(f_out)
  check(identical(w$geno, v$geno) && identical(w$map[, c("marker", "chr", "pos", "ref", "alt")],
                                               v$map[, c("marker", "chr", "pos", "ref", "alt")]),
        paste("VCF write/read round trip is exact for", ext))
}
err <- tryCatch(write_vcf(a$geno, a$map, file.path(tmp, "noalleles.vcf")), error = conditionMessage)
check(is.character(err) && grepl("REF and ALT", err), "VCF writing refuses markers without alleles")

hmp <- c(paste(c("rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center", "protLSID",
                 "assayLSID", "panelLSID", "QCcode", "L1", "L2", "L3"), collapse = "\t"),
         paste(c("m1", "A/G", "1", "1000", "+", rep("NA", 6), "AA", "R", "GG"), collapse = "\t"),
         paste(c("m2", "C/T", "1", "2000", "+", rep("NA", 6), "T", "N", "CC"), collapse = "\t"),
         paste(c("m3", "A/C/G", "2", "500", "+", rep("NA", 6), "AA", "CC", "GG"), collapse = "\t"))
f_hmp <- file.path(tmp, "tiny.hmp.txt")
writeLines(hmp, f_hmp)
hm <- read_genotypes(f_hmp)
check(hm$format == "hapmap" && identical(unname(hm$geno), matrix(c(0L, 2L, 1L, NA, 2L, 0L), 2)) &&
        identical(hm$map$alt, c("G", "T")), "HapMap file read; tri-allelic marker excluded")

# 10. End to end on the demo data ----------------------------------------------------------
out <- file.path(tmp, "demo")
res <- run_cli(c("--demo", "--write-vcf", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "tables/genotypes_filtered_012.csv", "tables/marker_qc.csv", "tables/sample_qc.csv",
            "tables/duplicate_pairs.csv", "tables/qc_summary.csv", "figures/marker_qc.png", "figures/sample_qc.png",
            "figures/ibs_pairs.png", "genotypes_filtered.vcf", "reproducibility/commands.sh",
            "reproducibility/sessionInfo.txt")) {
  check(file.exists(file.path(out, f)), paste("demo writes", f))
}
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
truth <- utils::read.csv(file.path(SKILL_TEST_DIR, "demo", "demo_truth.csv"), stringsAsFactors = FALSE)
tr <- split(truth$id, truth$problem)
mq <- utils::read.csv(file.path(out, "tables", "marker_qc.csv"), stringsAsFactors = FALSE)
sq <- utils::read.csv(file.path(out, "tables", "sample_qc.csv"), stringsAsFactors = FALSE)
dp <- utils::read.csv(file.path(out, "tables", "duplicate_pairs.csv"), stringsAsFactors = FALSE)
check(same_set(mq$marker[mq$fail_call_rate], tr$marker_call_rate) &&
        same_set(mq$marker[mq$monomorphic & !mq$fail_call_rate], tr$marker_monomorphic) &&
        same_set(sq$sample[sq$fail_call_rate], tr$sample_call_rate) &&
        same_set(sq$sample[sq$het_flag], tr$sample_heterozygosity) &&
        same_set(apply(dp[, 1:2], 1, function(x) paste(sort(x), collapse = "~")), tr$duplicate_pair),
      "demo: every planted problem is reported, and nothing else")
fm <- read_genotypes(file.path(out, "tables", "genotypes_filtered_012.csv"))
fv <- read_genotypes(file.path(out, "genotypes_filtered.vcf"))
check(identical(fm$geno, fv$geno) && nrow(fm$geno) == 560L && ncol(fm$geno) == 116L,
      "demo: filtered CSV and filtered VCF hold the same 560 x 116 genotypes")

# 11. Bad input fails loudly ------------------------------------------------------------------
bad <- run_cli(c("--input", shQuote(f_vcf), "--out", shQuote(file.path(tmp, "bad"))))
check(bad$status != 0L && any(grepl("germplasm", bad$output)), "missing --germplasm stops with a message")
bad <- run_cli(c("--demo", "--color", "red"))
check(bad$status != 0L, "unknown option exits with an error")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
