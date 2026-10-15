# Install every R package used by AgriClaw skills (CRAN + Bioconductor).
# Safe to re-run: only missing packages are installed.
#
#   Rscript tools/install_r_deps.R
#
# On Windows and macOS, binary packages are used where available. GO.db is a
# data-only package, so it installs from source without a compiler.

options(repos = c(CRAN = "https://cloud.r-project.org"),
        Ncpus = max(1L, parallel::detectCores() - 1L))
binary_ok <- .Platform$OS.type == "windows" || Sys.info()[["sysname"]] == "Darwin"
if (binary_ok) options(install.packages.compile.from.source = "never")

cran <- c(
  # trials
  "agricolae", "FielDHub", "SpATS", "statgenSTA", "statgenGxE", "lme4", "lmerTest", "nlme", "emmeans", "desplot",
  # genomics
  "statgenGWAS", "qtl2", "rrBLUP", "BGLR", "vcfR",
  # remote sensing
  "terra", "sf", "exactextractr", "ranger", "glmnet", "signal", "phenofit",
  # physiology
  "PhotoGEA", "segmented", "minpack.lm",
  # field / data access
  "soilDB", "jsonlite", "httr2", "curl", "readxl", "openxlsx", "data.table",
  # co-expression and utilities
  "WGCNA", "fastcluster", "dynamicTreeCut", "readr", "openssl", "testthat"
)
bioc <- c("DESeq2", "tximport", "edgeR", "limma", "apeglm", "sva",
          "AnnotationDbi", "GO.db", "impute", "preprocessCore", "clusterProfiler")

have <- function() rownames(installed.packages())
if (!"BiocManager" %in% have()) install.packages("BiocManager")

need_b <- setdiff(bioc, have())
if ("GO.db" %in% need_b) {
  BiocManager::install("GO.db", update = FALSE, ask = FALSE, type = "source")
  need_b <- setdiff(need_b, "GO.db")
}
if (length(need_b)) BiocManager::install(need_b, update = FALSE, ask = FALSE)

need_c <- setdiff(cran, have())
if (length(need_c)) install.packages(need_c, dependencies = c("Depends", "Imports", "LinkingTo"))

failed <- character()
for (p in c(cran, bioc)) {
  if (!suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))) failed <- c(failed, p)
}
if (length(failed)) {
  cat("Could not load:", paste(failed, collapse = ", "), "\n")
  quit(status = 1)
}
cat(sprintf("All %d AgriClaw R packages load.\n", length(cran) + length(bioc)))
