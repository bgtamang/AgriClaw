# Generate the synthetic claw-enrichment demo (fixed seed). Run from anywhere:
#   Rscript demo/make_demo.R
# Writes demo/de_results.csv (a DE table like claw-de-analysis output) and
# demo/go_annotation.tsv (one row per gene, comma-separated GO IDs, most specific
# terms only, like the GO column of a Phytozome annotation_info file).
# Gene IDs and annotations are invented; the GO IDs are real (checked in GO.db).

demo_dir <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else getwd()
})
set.seed(2026)

n_genes <- 3000
genes <- sprintf("Demo.%02dG%04d00", rep(1:10, each = 300), rep(1:300, 10))
tested <- seq_len(n_genes) %in% sample(n_genes, 2400)        # 600 genes filtered before testing
up <- seq_len(n_genes) %in% sample(which(tested), 150)
down <- seq_len(n_genes) %in% sample(which(tested & !up), 150)

planted_up <- c("GO:0009414", "GO:0042631", "GO:0009737", "GO:0009269")   # drought / ABA responses
planted_down <- c("GO:0009765", "GO:0009768", "GO:0016168", "GO:0009535", "GO:0009523")  # light harvesting
null_terms <- c("GO:0009408", "GO:0006412", "GO:0006457", "GO:0055085", "GO:0006355", "GO:0006468",
                "GO:0005975", "GO:0006979", "GO:0009809", "GO:0010200", "GO:0006952", "GO:0009733",
                "GO:0009723", "GO:0006970", "GO:0009651", "GO:0006511", "GO:0007165", "GO:0006629",
                "GO:0006364", "GO:0006260", "GO:0007049", "GO:0009644", "GO:0009631", "GO:0016042",
                "GO:0010286", "GO:0005524", "GO:0003700", "GO:0004672", "GO:0016491", "GO:0005515",
                "GO:0003735", "GO:0046872", "GO:0005634", "GO:0005886", "GO:0005840", "GO:0005737")
null_w <- stats::rgamma(length(null_terms), 1.5)                 # some terms are common, some rare
stopifnot(all(c(planted_up, planted_down, null_terms) %in% AnnotationDbi::keys(GO.db::GO.db)))

ann <- vector("list", n_genes)
for (g in seq_len(n_genes)) {
  terms <- character()
  if (stats::runif(1) < 0.72)                                    # about 28% of genes unannotated
    terms <- sample(null_terms, sample(1:3, 1), prob = null_w)
  p_up <- if (up[g]) 0.30 else 0.015
  p_dn <- if (down[g]) 0.30 else 0.015
  terms <- c(terms, planted_up[stats::runif(length(planted_up)) < p_up],
             planted_down[stats::runif(length(planted_down)) < p_dn])
  ann[[g]] <- unique(terms)
}
obs <- sample(which(lengths(ann) > 0), 5)
ann[obs] <- lapply(ann[obs], function(t) c(t, "GO:0018298"))  # obsolete in GO.db 3.23
go_col <- vapply(ann, paste, "", collapse = ",")
utils::write.table(data.frame(gene_id = genes, GO = go_col), file.path(demo_dir, "go_annotation.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

lfc <- ifelse(up, stats::runif(n_genes, 1, 3), ifelse(down, -stats::runif(n_genes, 1, 3),
                                                      stats::rnorm(n_genes, 0, 0.3)))
padj <- ifelse(up | down, 10^-stats::runif(n_genes, 2, 12), stats::runif(n_genes, 0.06, 1))
padj[!tested] <- NA
de <- data.frame(gene_id = genes, baseMean = round(stats::rlnorm(n_genes, log(200), 1.2), 2),
                 log2FoldChange = round(lfc, 4), padj = signif(padj, 4))
de$baseMean[!tested] <- round(stats::runif(sum(!tested), 0.5, 8), 2)
de <- de[order(de$padj, na.last = TRUE), ]
utils::write.csv(de, file.path(demo_dir, "de_results.csv"), row.names = FALSE)
cat("Demo written to", demo_dir, "\n")
