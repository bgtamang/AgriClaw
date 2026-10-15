# Generate the synthetic claw-coexpression demo (fixed seed). Run from anywhere:
#   Rscript demo/make_demo.R
# Writes demo/expr_vst.csv (600 genes x 24 samples on a VST-like log2 scale),
# demo/traits.csv and demo/true_modules.csv (the planted module of every gene).
# Expression is simulated with WGCNA::simulateDatExpr: 4 modules + unrelated genes.

demo_dir <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else getwd()
})
set.seed(7)

n <- 24
samples <- sprintf("S%02d", seq_len(n))
treatment <- rep(c("control", "drought"), each = n / 2)
me <- data.frame(ME1 = as.numeric(scale(1.5 * (treatment == "drought") + stats::rnorm(n))),  # drought module
                 ME2 = stats::rnorm(n), ME3 = stats::rnorm(n), ME4 = stats::rnorm(n))
n_genes <- 600
sim <- WGCNA::simulateDatExpr(me, nGenes = n_genes, modProportions = c(0.18, 0.14, 0.11, 0.08, 0.49),
                              minCor = 0.4, maxCor = 0.95, signed = TRUE, backgroundNoise = 0.2, verbose = 0)
x <- as.matrix(sim$datExpr)
gene_sd <- stats::runif(n_genes, 0.3, 1.2)
gene_mean <- stats::runif(n_genes, 6, 12)
x <- sweep(sweep(scale(x), 2, gene_sd, "*"), 2, gene_mean, "+")     # VST-like log2 scale
genes <- sprintf("Demo.%02dG%04d00", rep(1:6, each = 100), rep(1:100, 6))
colnames(x) <- genes
rownames(x) <- samples

utils::write.csv(data.frame(gene_id = genes, round(t(x), 4), check.names = FALSE),
                 file.path(demo_dir, "expr_vst.csv"), row.names = FALSE)
traits <- data.frame(sample = samples, treatment = treatment,
                     stomatal_conductance_mol_m2_s = round(0.35 - 0.06 * me$ME1 + stats::rnorm(n, 0, 0.03), 4),
                     leaf_n_mg_g = round(32 + 2.5 * me$ME3 + stats::rnorm(n, 0, 1.5), 2),
                     plant_height_cm = round(85 + stats::rnorm(n, 0, 6), 1))
utils::write.csv(traits, file.path(demo_dir, "traits.csv"), row.names = FALSE)
utils::write.csv(data.frame(gene_id = genes, true_module = sim$allLabels),
                 file.path(demo_dir, "true_modules.csv"), row.names = FALSE)
cat("Demo written to", demo_dir, "\n")
