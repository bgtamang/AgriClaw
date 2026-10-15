# Generate the synthetic claw-rnaseq-import demo (fixed seed). Run from anywhere:
#   Rscript demo/make_demo.R
# Writes demo/annotation.gff3, demo/samples.csv and demo/quant/S1..S6/ (quant.sf and
# aux_info/meta_info.json, as written by Salmon). All values are simulated.

demo_dir <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else getwd()
})
set.seed(20261015)

n_genes <- 200
chrom <- sprintf("Dm%02d", rep(1:4, each = n_genes / 4))
gene_num <- ave(seq_len(n_genes), chrom, FUN = seq_along)
gene_id <- sprintf("Demo.%sG%04d00", substring(chrom, 3), gene_num)
suffix <- ".Dm1.v1"                        # Phytozome-style annotation-version suffix in ID=
n_iso <- sample(1:3, n_genes, replace = TRUE, prob = c(0.6, 0.3, 0.1))
is_nc <- seq_len(n_genes) %in% c(17, 88, 151)  # three long non-coding RNA genes
n_iso[is_nc] <- 1L

# --- GFF3 -------------------------------------------------------------------
gff <- c("##gff-version 3", sprintf("##sequence-region %s 1 2000000", unique(chrom)))
tx <- data.frame(transcript_id = character(), gene_id = character(), length = integer(),
                 stringsAsFactors = FALSE)
pos <- stats::setNames(rep(1000L, 4), unique(chrom))
for (g in seq_len(n_genes)) {
  ch <- chrom[g]
  start <- pos[[ch]]
  glen <- sample(1500:6000, 1)
  end <- start + glen
  strand <- sample(c("+", "-"), 1)
  gff <- c(gff, sprintf("%s\tdemo\tgene\t%d\t%d\t.\t%s\t.\tID=%s%s;Name=%s",
                        ch, start, end, strand, gene_id[g], suffix, gene_id[g]))
  for (k in seq_len(n_iso[g])) {
    tid <- sprintf("%s.%d", gene_id[g], k)
    ttype <- if (is_nc[g]) "lnc_RNA" else "mRNA"
    gff <- c(gff, sprintf("%s\tdemo\t%s\t%d\t%d\t.\t%s\t.\tID=%s%s;Name=%s;Parent=%s%s",
                          ch, ttype, start, end, strand, tid, suffix, tid, gene_id[g], suffix))
    n_ex <- sample(2:3, 1)
    cuts <- sort(sample(seq(start + 100, end - 100), 2 * n_ex - 2))
    ex_s <- c(start, cuts[seq(2, length(cuts), 2)])
    ex_e <- c(cuts[seq(1, length(cuts), 2)], end)
    for (e in seq_len(n_ex)) {
      gff <- c(gff, sprintf("%s\tdemo\texon\t%d\t%d\t.\t%s\t.\tParent=%s%s",
                            ch, ex_s[e], ex_e[e], strand, tid, suffix))
      if (!is_nc[g] && k == 1L)
        gff <- c(gff, sprintf("%s\tdemo\tCDS\t%d\t%d\t.\t%s\t0\tParent=%s%s",
                              ch, ex_s[e], ex_e[e], strand, tid, suffix))
    }
    tx[nrow(tx) + 1L, ] <- list(tid, gene_id[g], sum(ex_e - ex_s + 1L))
  }
  pos[[ch]] <- end + sample(2000:20000, 1)
}
writeLines(gff, file.path(demo_dir, "annotation.gff3"))

# Two chloroplast transcripts are in the (simulated) Salmon index but not in the
# nuclear GFF3, as happens when organelle sequences are added to the index.
tx_all <- rbind(tx, data.frame(transcript_id = c("ChrC_psbA", "ChrC_rbcL"), gene_id = NA,
                               length = c(1062L, 1428L)))

# --- Expression ---------------------------------------------------------------
samples <- data.frame(sample = sprintf("S%d", 1:6),
                      condition = rep(c("control", "drought"), each = 3),
                      replicate = rep(1:3, 2), stringsAsFactors = FALSE)
n_processed <- round(c(12.1, 10.8, 13.0, 11.6, 12.4, 14.2) * 1e6)
pct_mapped <- c(88.4, 89.9, 87.1, 90.2, 88.8, 32.5)   # S6: planted low mapping rate
n_mapped <- round(n_processed * pct_mapped / 100)

g_mean <- stats::rlnorm(n_genes, meanlog = log(400), sdlog = 1.6)   # relative abundance
lfc <- rep(0, n_genes)
de <- sample(which(!is_nc), 40)
lfc[de] <- sample(c(-1, 1), 40, replace = TRUE) * stats::runif(40, 1, 2.5)
iso_frac <- lapply(n_iso, function(k) { w <- stats::rgamma(k, 2); w / sum(w) })
switch_genes <- utils::head(which(n_iso >= 2 & !seq_len(n_genes) %in% de), 6)  # isoform switching
organelle_share <- 0.04                                                 # of mapped reads

dir.create(file.path(demo_dir, "quant"), showWarnings = FALSE)
for (j in seq_len(nrow(samples))) {
  drought <- samples$condition[j] == "drought"
  mu <- g_mean * 2^(lfc * drought)
  mu <- mu / sum(mu) * n_mapped[j] * (1 - organelle_share)
  gcount <- stats::rnbinom(n_genes, mu = mu, size = 1 / 0.04)   # dispersion 0.04
  numreads <- numeric(nrow(tx))
  for (g in seq_len(n_genes)) {
    fr <- iso_frac[[g]]
    if (drought && g %in% switch_genes) fr <- rev(fr)
    idx <- which(tx$gene_id == gene_id[g])
    numreads[idx] <- round(gcount[g] * fr, 3)
  }
  org <- round(n_mapped[j] * organelle_share * c(0.7, 0.3), 3)
  nr <- c(numreads, org)
  frag <- stats::runif(1, 185, 215)
  eff <- round(pmax(tx_all$length - frag, 50), 3)
  rate <- nr / eff
  tpm <- 1e6 * rate / sum(rate)
  qd <- file.path(demo_dir, "quant", samples$sample[j])
  dir.create(file.path(qd, "aux_info"), recursive = TRUE, showWarnings = FALSE)
  q <- data.frame(Name = tx_all$transcript_id, Length = tx_all$length,
                  EffectiveLength = sprintf("%.3f", eff), TPM = sprintf("%.6f", tpm),
                  NumReads = sprintf("%.3f", nr))
  utils::write.table(q, file.path(qd, "quant.sf"), sep = "\t", quote = FALSE, row.names = FALSE)
  meta <- c("{",
            '  "salmon_version": "1.10.1",',
            '  "note": "synthetic AgriClaw demo file, not produced by Salmon",',
            '  "library_types": ["ISR"],',
            sprintf('  "num_processed": %d,', n_processed[j]),
            sprintf('  "num_mapped": %d,', n_mapped[j]),
            sprintf('  "percent_mapped": %.4f', 100 * n_mapped[j] / n_processed[j]),
            "}")
  writeLines(meta, file.path(qd, "aux_info", "meta_info.json"))
}
samples$path <- file.path("quant", samples$sample)
utils::write.csv(samples[, c("sample", "path", "condition", "replicate")],
                 file.path(demo_dir, "samples.csv"), row.names = FALSE)
cat("Demo written to", demo_dir, "\n")
