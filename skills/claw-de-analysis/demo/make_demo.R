# Generate the synthetic claw-de-analysis demo (fixed seed). Run from anywhere:
#   Rscript demo/make_demo.R
# Writes demo/counts.csv (raw counts), demo/samples.csv and demo/truth.csv.

demo_dir <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else getwd()
})
set.seed(42)

n_genes <- 1200
genes <- sprintf("Demo.%02dG%04d00", rep(1:10, each = n_genes / 10), rep(seq_len(n_genes / 10), 10))
samples <- data.frame(sample = sprintf("P%02d", 1:8),
                      condition = rep(c("control", "drought"), each = 4),
                      batch = rep(c("b1", "b1", "b2", "b2"), 2), stringsAsFactors = FALSE)

base_mean <- stats::rlnorm(n_genes, meanlog = log(120), sdlog = 1.7)        # counts at size factor 1
disp <- (0.05 + 1 / base_mean) * exp(stats::rnorm(n_genes, 0, 0.3))         # NB dispersion trend + scatter
is_de <- seq_len(n_genes) %in% sample(n_genes, 120)
lfc <- ifelse(is_de, sample(c(-1, 1), n_genes, replace = TRUE) * stats::runif(n_genes, 0.8, 2.5), 0)
batch_eff <- stats::rnorm(n_genes, 0, 0.3)                                   # log2 shift in batch b2
size_factor <- c(0.8, 1.2, 1.0, 1.4, 0.9, 1.1, 0.7, 1.3)

counts <- matrix(0L, n_genes, nrow(samples), dimnames = list(genes, samples$sample))
for (j in seq_len(nrow(samples))) {
  mu <- base_mean * size_factor[j] * 2^(lfc * (samples$condition[j] == "drought") +
                                          batch_eff * (samples$batch[j] == "b2"))
  counts[, j] <- stats::rnbinom(n_genes, mu = mu, size = 1 / disp)
}
utils::write.csv(data.frame(gene_id = genes, counts, check.names = FALSE),
                 file.path(demo_dir, "counts.csv"), row.names = FALSE)
utils::write.csv(samples, file.path(demo_dir, "samples.csv"), row.names = FALSE)
utils::write.csv(data.frame(gene_id = genes, true_log2fc = round(lfc, 4), is_de = is_de,
                            base_mean = round(base_mean, 2), dispersion = round(disp, 4)),
                 file.path(demo_dir, "truth.csv"), row.names = FALSE)
cat("Demo written to", demo_dir, "\n")
