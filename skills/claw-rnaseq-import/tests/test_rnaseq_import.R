# Tests for claw-rnaseq-import. Run from anywhere: Rscript tests/test_rnaseq_import.R
# Known answers: gene-level sums of synthetic Salmon files computed by hand, a planted
# low-mapping sample, ID version mismatches, GFF3/GTF parsing on hand-written files,
# and an end-to-end run of the command line on the demo data.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "rnaseq_import.R"))  # loads functions only

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
expect_error <- function(expr, pattern, what) {
  msg <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  check(!is.na(msg) && grepl(pattern, msg), paste0(what, " (got: ", msg, ")"))
  invisible(msg)
}
run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "rnaseq_import.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
tmp <- file.path(tempdir(), "claw-rnaseq-import-test")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE)

# Synthetic Salmon output: n_samples x transcripts, values written with fixed decimals.
write_salmon <- function(root, tx_ids, len, numreads, eff, pct_mapped = NULL) {
  files <- character()
  for (j in seq_len(ncol(numreads))) {
    s <- colnames(numreads)[j]
    d <- file.path(root, s)
    dir.create(file.path(d, "aux_info"), recursive = TRUE, showWarnings = FALSE)
    rate <- numreads[, j] / eff[, j]
    tpm <- 1e6 * rate / sum(rate)
    q <- data.frame(Name = tx_ids, Length = len, EffectiveLength = sprintf("%.3f", eff[, j]),
                    TPM = sprintf("%.6f", tpm), NumReads = sprintf("%.3f", numreads[, j]))
    utils::write.table(q, file.path(d, "quant.sf"), sep = "\t", quote = FALSE, row.names = FALSE)
    if (!is.null(pct_mapped))
      writeLines(sprintf('{"num_processed": 1000000, "num_mapped": %d, "percent_mapped": %.3f, "library_types": ["ISR"]}',
                         round(1e4 * pct_mapped[j]), pct_mapped[j]),
                 file.path(d, "aux_info", "meta_info.json"))
    files[s] <- file.path(d, "quant.sf")
  }
  files
}

set.seed(1)
n_g <- 120
iso <- sample(1:3, n_g, replace = TRUE)
genes <- sprintf("Gm.%02dG%06d", rep(1:6, each = 20), seq_len(n_g) * 100)
tx_gene <- rep(genes, iso)
tx_ids <- paste0(tx_gene, ".", unlist(lapply(iso, seq_len)))
n_tx <- length(tx_ids)
len <- sample(400:4000, n_tx, replace = TRUE)
samp <- sprintf("S%d", 1:6)
nr <- matrix(round(stats::rgamma(n_tx * 6, 0.6, 1 / 800), 3), n_tx, 6, dimnames = list(tx_ids, samp))
nr[5, ] <- 0                                   # a transcript never seen
zero_gene <- genes[iso == 1][3]                # a single-isoform gene with zero reads in S2
nr[tx_gene == zero_gene, "S2"] <- 0
eff <- sapply(seq_len(6), function(j) round(pmax(len - 180 - j, 30), 3))
dimnames(eff) <- dimnames(nr)
files <- write_salmon(file.path(tmp, "salmon"), tx_ids, len, nr, eff, pct_mapped = c(88, 90, 87, 89, 91, 31))
tx2gene <- data.frame(transcript_id = tx_ids, gene_id = tx_gene)

# 1. Gene-level sums equal hand-computed values ------------------------------
txi_tx <- import_transcripts(files, "salmon")
sg <- summarise_to_genes(txi_tx, tx2gene)
txi <- sg$txi
read_back <- function(col) sapply(files, function(f) as.numeric(utils::read.delim(f)[[col]]))
nr_file <- read_back("NumReads"); tpm_file <- read_back("TPM"); eff_file <- read_back("EffectiveLength")
exp_counts <- rowsum(nr_file, tx_gene)
exp_tpm <- rowsum(tpm_file, tx_gene)
check(identical(sort(rownames(txi$counts)), sort(genes)), "one row per gene")
check(max(abs(txi$counts[rownames(exp_counts), ] - exp_counts)) < 1e-8,
      "gene counts = sum of transcript NumReads")
check(max(abs(txi$abundance[rownames(exp_tpm), ] - exp_tpm)) < 1e-8, "gene TPM = sum of transcript TPM")
exp_len <- rowsum(tpm_file * eff_file, tx_gene) / exp_tpm
ok <- is.finite(exp_len)
check(max(abs(txi$length[rownames(exp_len), ][ok] - exp_len[ok])) < 1e-6,
      "gene length = TPM-weighted mean of transcript effective lengths")
zl <- txi$length[zero_gene, ]
check(abs(zl[["S2"]] - exp(mean(log(zl[-2])))) < 1e-8,
      "zero-abundance gene gets the geometric mean length of the other samples (tximport rule)")
check(identical(txi$countsFromAbundance, "no"), "countsFromAbundance = 'no' (original counts kept)")
check(sg$match_frac == 1 && length(sg$unmatched) == 0, "all transcripts matched")

# DESeq2 offsets: DESeqDataSetFromTximport accepts the object and stores avgTxLength
dds <- suppressMessages(DESeq2::DESeqDataSetFromTximport(
  txi, data.frame(sample = samp, row.names = samp), ~ 1))
check("avgTxLength" %in% SummarizedExperiment::assayNames(dds),
      "txi feeds DESeqDataSetFromTximport with length offsets")

# kallisto: abundance.tsv + run_info.json in per-sample folders give the same gene sums
kdir <- file.path(tmp, "kallisto")
for (j in 1:3) {
  d <- file.path(kdir, samp[j])
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  rate <- nr[, j] / eff[, j]
  utils::write.table(data.frame(target_id = tx_ids, length = len, eff_length = sprintf("%.3f", eff[, j]),
                                est_counts = sprintf("%.3f", nr[, j]), tpm = sprintf("%.6f", 1e6 * rate / sum(rate))),
                     file.path(d, "abundance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines('{"n_processed": 2000000, "n_pseudoaligned": 1700000, "p_pseudoaligned": 85.0}',
             file.path(d, "run_info.json"))
}
rk <- resolve_quant_files(data.frame(sample = samp[1:3], path = samp[1:3]), kdir)
check(rk$type == "kallisto", "kallisto detected from abundance.tsv inside sample folders")
sk <- summarise_to_genes(import_transcripts(rk$files, "kallisto"), tx2gene)
check(max(abs(sk$txi$counts[rownames(exp_counts), ] - exp_counts[, 1:3])) < 1e-8,
      "kallisto: gene counts = hand sums of est_counts")
check(read_mapping_info(rk$files[[1]], "kallisto")$mapping_rate_pct == 85, "kallisto mapping rate read from run_info.json")

# 2. Planted low-mapping sample is flagged, the others are not ---------------
map <- do.call(rbind, lapply(names(files), function(s) {
  m <- read_mapping_info(files[[s]], "salmon")
  data.frame(sample = s, n_processed = m$n_processed, n_mapped = m$n_mapped,
             mapping_rate_pct = m$mapping_rate_pct)
}))
check(all(abs(map$mapping_rate_pct - c(88, 90, 87, 89, 91, 31)) < 1e-9),
      "mapping rate read from aux_info/meta_info.json")
tm <- qc_transform(txi$counts, txi$length)
sq <- sample_qc(txi$counts, map, tm, detect_min = 10, min_map = 50, k = 3)
check(grepl("mapping rate < 50%", sq$qc$flags[6]) && grepl("MAD", sq$qc$flags[6]),
      "planted S6 flagged by the absolute and the MAD mapping-rate rules")
check(!any(grepl("mapping rate", sq$qc$flags[1:5])), "no other sample gets a mapping-rate flag")

# MAD rule on a hand-computed example: median 10, MAD 1.4826 -> cut-off 5.55
check(identical(mad_low(c(10, 11, 9, 10, 12, 2), 3), c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)),
      "MAD outlier rule")

# 3. ID version mismatches are caught -----------------------------------------
ens <- sprintf("TXA%04d", seq_len(n_tx))
files_v <- write_salmon(file.path(tmp, "salmon_ver"), paste0(ens, ".2"), len, nr, eff)
t2g_nov <- data.frame(transcript_id = ens, gene_id = tx_gene)
txv <- import_transcripts(files_v, "salmon")
msg <- expect_error(summarise_to_genes(txv, t2g_nov), "trailing '\\.N'",
                    "versioned quant IDs vs unversioned tx2gene stops with a diagnosis")
check(grepl("--strip-tx-version", msg), "the error suggests --strip-tx-version")
sgv <- summarise_to_genes(txv, t2g_nov, strip_version = TRUE)
check(max(abs(sgv$txi$counts[rownames(exp_counts), ] - exp_counts)) < 1e-8,
      "after --strip-tx-version the gene sums are correct")
# the reverse direction: annotation versioned, quantification not
expect_error(summarise_to_genes(import_transcripts(write_salmon(file.path(tmp, "salmon_nov"),
             ens, len, nr, eff), "salmon"), data.frame(transcript_id = paste0(ens, ".1"),
             gene_id = tx_gene)), "annotation IDs carry a trailing", "unversioned quant vs versioned tx2gene")
# plant isoform numbers are not versions: stripping would merge isoforms
expect_error(summarise_to_genes(txi_tx, tx2gene, strip_version = TRUE), "isoform number",
             "--strip-tx-version refuses to merge plant isoforms (Gm.01G000100.1/.2)")
# tximport's own ignoreTxVersion cuts at the first dot: show why it is unsafe here
check(all(sub("\\..*", "", tx_ids) == "Gm"),
      "tximport-style first-dot stripping would collapse every soybean-style ID")
# soybean Wm82.a1 vs a2+ nomenclature is diagnosed
a1 <- sprintf("Glyma%02dg%05d.1", rep(1:6, each = 20), seq_len(120) * 10)
a2 <- sprintf("Glyma.%02dG%06d.1", rep(1:6, each = 20), seq_len(120) * 100)
h <- diagnose_id_mismatch(a2, a1)
check(any(grepl("Wm82.a1", h)) && any(grepl("Wm82.a2", h)), "a1 vs a2 soybean IDs are named in the diagnosis")
# tx2gene table that maps one transcript to two genes is rejected
bad_t2g <- file.path(tmp, "bad_t2g.csv")
utils::write.csv(data.frame(transcript_id = c("t1", "t1"), gene_id = c("g1", "g2")), bad_t2g, row.names = FALSE)
expect_error(read_tx2gene_table(bad_t2g), "more than one gene", "ambiguous tx2gene rejected")

# 4. GFF3 / GTF parsing ---------------------------------------------------------
phyto <- file.path(tmp, "phyto.gff3")
writeLines(c("##gff-version 3",
  "Gm01\tphytozome\tgene\t100\t900\t.\t+\t.\tID=Glyma.01G000100.Wm82.a4.v1;Name=Glyma.01G000100",
  "Gm01\tphytozome\tmRNA\t100\t900\t.\t+\t.\tID=Glyma.01G000100.1.Wm82.a4.v1;Name=Glyma.01G000100.1;Parent=Glyma.01G000100.Wm82.a4.v1",
  "Gm01\tphytozome\texon\t100\t900\t.\t+\t.\tID=Glyma.01G000100.1.Wm82.a4.v1.exon.1;Parent=Glyma.01G000100.1.Wm82.a4.v1",
  "Gm01\tphytozome\tmRNA\t100\t700\t.\t+\t.\tID=Glyma.01G000100.2.Wm82.a4.v1;Name=Glyma.01G000100.2;Parent=Glyma.01G000100.Wm82.a4.v1",
  "Gm01\tphytozome\tgene\t2000\t2900\t.\t-\t.\tID=Glyma.01G000200.Wm82.a4.v1;Name=Glyma.01G000200",
  "Gm01\tphytozome\tmRNA\t2000\t2900\t.\t-\t.\tID=Glyma.01G000200.1.Wm82.a4.v1;Name=Glyma.01G000200.1;Parent=Glyma.01G000200.Wm82.a4.v1",
  "Gm01\tphytozome\tCDS\t2000\t2900\t.\t-\t0\tParent=Glyma.01G000200.1.Wm82.a4.v1",
  "##FASTA", ">Gm01", "ACGT"), phyto)
pa <- parse_annotation(phyto)
ch <- choose_tx2gene(pa$candidates, c("Glyma.01G000100.1", "Glyma.01G000100.2", "Glyma.01G000200.1"))
check(ch$style == "Name" && ch$match[["Name"]] == 1 && ch$match[["ID"]] == 0,
      "Phytozome-style GFF3: Name attribute chosen automatically")
check(identical(ch$tx2gene$gene_id[match("Glyma.01G000100.2", ch$tx2gene$transcript_id)], "Glyma.01G000100"),
      "Phytozome-style GFF3: isoform 2 maps to the plain gene ID")
check(nrow(ch$tx2gene) == 3L, "exon/CDS lines and the ##FASTA section are ignored")

ensg <- file.path(tmp, "ensembl.gff3.gz")
con <- gzfile(ensg, "w")
writeLines(c("##gff-version 3",
  "1\tgramene\tgene\t1\t500\t.\t+\t.\tID=gene:Zm00001eb000010;biotype=protein_coding;gene_id=Zm00001eb000010",
  "1\tgramene\tmRNA\t1\t500\t.\t+\t.\tID=transcript:Zm00001eb000010_T001;Parent=gene:Zm00001eb000010;transcript_id=Zm00001eb000010_T001",
  "1\tgramene\tmRNA\t1\t400\t.\t+\t.\tID=transcript:Zm00001eb000010_T002;Parent=gene:Zm00001eb000010;transcript_id=Zm00001eb000010_T002",
  "1\tgramene\tncRNA_gene\t900\t990\t.\t+\t.\tID=gene:Zm00001eb000020;Name=my%3Bnc",
  "1\tgramene\tlnc_RNA\t900\t990\t.\t+\t.\tID=transcript:Zm00001eb000020_T001;Parent=gene:Zm00001eb000020",
  "1\tgramene\texon\t1\t500\t.\t+\t.\tParent=transcript:Zm00001eb000010_T001"), con)
close(con)
pe <- parse_annotation(ensg)
che <- choose_tx2gene(pe$candidates, c("Zm00001eb000010_T001", "Zm00001eb000010_T002", "Zm00001eb000020_T001"))
check(che$match[[che$style]] == 1 && all(che$tx2gene$gene_id == c("Zm00001eb000010", "Zm00001eb000010", "Zm00001eb000020")),
      "Ensembl-style gzipped GFF3: 'transcript:'/'gene:' prefixes removed, ncRNA genes kept")
check(identical(attr_value("ID=a%3Bb;Name=x", "ID"), "a;b"), "GFF3 %XX escapes are decoded")

gtf <- file.path(tmp, "a.gtf")
writeLines(c('1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id "OsG1"; transcript_id "OsG1.1";',
             '1\tsrc\texon\t1\t100\t.\t+\t.\tgene_id "OsG1"; transcript_id "OsG1.1";',
             '1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id "OsG1"; transcript_id "OsG1.2";'), gtf)
pg <- parse_annotation(gtf)
check(pg$format == "GTF" && nrow(pg$candidates$transcript_id) == 2L &&
        all(pg$candidates$transcript_id$gene_id == "OsG1"), "GTF transcript_id/gene_id parsed")

# 5. featureCounts input ---------------------------------------------------------
fcf <- file.path(tmp, "fc.txt")
writeLines(c("# Program:featureCounts v2.0.6; Command:\"featureCounts\" ...",
             "Geneid\tChr\tStart\tEnd\tStrand\tLength\t/data/bam/S1.sorted.bam\t/data/bam/S2.sorted.bam",
             "g1\t1\t1\t100\t+\t100\t5\t7", "g2\t1\t200\t400\t-\t201\t0\t12"), fcf)
fc <- read_featurecounts(fcf)
cm <- match_fc_columns(data.frame(sample = c("S2", "S1")), colnames(fc$counts))
check(identical(unname(fc$counts[, cm[["S2"]]]), c(7, 12)), "featureCounts columns matched to samples by BAM name")

# 6. End to end on the demo data ---------------------------------------------------
out <- file.path(tmp, "demo_out")
res <- run_cli(c("--demo", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "txi.rds", "tables/gene_counts.csv", "tables/gene_lengths.csv",
            "tables/gene_tpm.csv", "tables/sample_qc.csv", "tables/samples.csv", "tables/tx2gene.csv",
            "tables/unmatched_transcripts.csv", "figures/qc_overview.png", "figures/pca.png",
            "figures/sample_correlation.png", "reproducibility/commands.sh",
            "reproducibility/sessionInfo.txt")) check(file.exists(file.path(out, f)), paste("demo writes", f))
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
qc <- utils::read.csv(file.path(out, "tables", "sample_qc.csv"), stringsAsFactors = FALSE)
check(identical(qc$sample[qc$flagged], "S6"), "demo: only the planted low-mapping sample S6 is flagged")
um <- utils::read.csv(file.path(out, "tables", "unmatched_transcripts.csv"))
check(setequal(um$transcript_id, c("ChrC_psbA", "ChrC_rbcL")), "demo: the 2 organelle transcripts are reported as unmatched")
rep_txt <- readLines(file.path(out, "report.md"))
check(any(grepl("transcript IDs taken from `Name`", rep_txt, fixed = TRUE)), "demo: GFF3 Name attribute chosen")
check(grep("## Warnings", rep_txt) < grep("## Summary", rep_txt), "warnings come before the summary")
dtxi <- readRDS(file.path(out, "txi.rds"))
demo_files <- Sys.glob(file.path(SKILL_TEST_DIR, "demo", "quant", "S*", "quant.sf"))
names(demo_files) <- basename(dirname(demo_files))
d_nr <- sapply(demo_files, function(f) { q <- utils::read.delim(f); stats::setNames(q$NumReads, q$Name) })
d_t2g <- utils::read.csv(file.path(out, "tables", "tx2gene.csv"))
d_exp <- rowsum(d_nr[d_t2g$transcript_id, ], d_t2g$gene_id)
check(max(abs(dtxi$counts[rownames(d_exp), colnames(d_exp)] - d_exp)) < 1e-6,
      "demo: gene counts equal hand sums of quant.sf NumReads")

three <- main(c("--demo", "--three-prime", "--out", file.path(tmp, "three_prime")))
tp_rep <- readLines(file.path(tmp, "three_prime", "report.md"))
check(any(grepl("--no-length-offset", tp_rep, fixed = TRUE)), "--three-prime: the report tells claw-de-analysis to drop length offsets")

# 7. Bad inputs fail loudly ----------------------------------------------------------
bad <- run_cli(c("--demo", "--colour", "red"))
check(bad$status != 0L, "unknown option exits with an error")
a1_t2g <- file.path(tmp, "a1_tx2gene.csv")
utils::write.csv(data.frame(transcript_id = sub("^Demo\\.0([1-4])G", "Glyma0\\1g", d_t2g$transcript_id),
                            gene_id = d_t2g$gene_id), a1_t2g, row.names = FALSE)
bad2 <- run_cli(c("--samples", shQuote(file.path(SKILL_TEST_DIR, "demo", "samples.csv")),
                  "--tx2gene", shQuote(a1_t2g), "--out", shQuote(file.path(tmp, "bad2"))))
check(bad2$status != 0L && any(grepl("are in tx2gene", bad2$output)),
      "CLI stops when the tx2gene IDs do not match the quantification")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
