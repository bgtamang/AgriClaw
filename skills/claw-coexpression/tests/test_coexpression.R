# Tests for claw-coexpression. Run from anywhere: Rscript tests/test_coexpression.R
# Known answers: planted modules from WGCNA::simulateDatExpr recovered (adjusted Rand
# index, core-gene accuracy), a planted module-trait association detected and a null
# trait not, soft-threshold selection rules, sample-outlier detection, exported TOM
# edges equal WGCNA's TOM, ComBat batch removal, refusals, and the demo end to end.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "coexpression.R"))  # loads functions only

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
  script <- file.path(SKILL_TEST_DIR, "scripts", "coexpression.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
tmp <- file.path(tempdir(), "claw-coexpression-test")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE)

# 1. Adjusted Rand index against hand-computed values --------------------------------
check(adjusted_rand(c(1, 1, 2, 2), c(2, 2, 1, 1)) == 1, "ARI = 1 for the same partition with permuted labels")
check(abs(adjusted_rand(c(1, 1, 2, 2), c(1, 1, 1, 2)) - 0) < 1e-12, "ARI = 0 for a hand-computed chance-level example")
check(abs(adjusted_rand(c(1, 1, 1, 2, 2, 2), c(1, 1, 2, 2, 3, 3)) - 0.8 / 3.3) < 1e-12, "ARI = 0.8/3.3 by hand (Hubert & Arabie formula)")

# 2. Planted modules are recovered (5 modules, 800 genes, 40 samples) ------------------
set.seed(1)
me <- data.frame(matrix(stats::rnorm(40 * 5), 40, 5))
names(me) <- paste0("ME", 1:5)
sim <- WGCNA::simulateDatExpr(me, nGenes = 800, modProportions = c(0.15, 0.12, 0.10, 0.08, 0.06, 0.49),
                              minCor = 0.4, maxCor = 0.95, signed = TRUE, backgroundNoise = 0.2, verbose = 0)
dat <- as.matrix(sim$datExpr)
dimnames(dat) <- list(sprintf("s%02d", 1:40), sprintf("g%04d", 1:800))
truth <- sim$allLabels
assoc_mod <- sim$labelOrder[2]                                     # true label seeded by ME2
traits <- data.frame(t_assoc = 0.8 * me$ME2 + stats::rnorm(40, 0, 0.6), t_null = stats::rnorm(40),
                     row.names = rownames(dat))
res <- run_wgcna(dat, traits, seed = 1)
lab <- res$labels
ari_mod <- adjusted_rand(truth[truth > 0], lab[truth > 0])
ari_all <- adjusted_rand(truth, lab)
cat(sprintf("  modules: %d found; ARI module genes %.3f, all genes %.3f; power %d\n",
            length(setdiff(unique(lab), 0)), ari_mod, ari_all, res$power))
check(length(setdiff(unique(lab), 0)) == 5L, "all 5 planted modules found, no extra modules")
check(ari_mod >= 0.85, sprintf("ARI among planted-module genes %.3f >= 0.85", ari_mod))
check(ari_all >= 0.6, sprintf("ARI over all genes including unassigned %.3f >= 0.6", ari_all))
tab <- table(truth, lab)
best <- apply(tab[rownames(tab) != "0", , drop = FALSE], 1, function(r) colnames(tab)[which.max(r)])
true_kme <- vapply(seq_len(ncol(dat)), function(i) if (truth[i] == 0) NA_real_ else
  stats::cor(dat[, i], me[[which(sim$labelOrder == truth[i])]]), 0)
core <- !is.na(true_kme) & true_kme >= 0.7
acc <- mean(as.character(lab[core]) == best[as.character(truth[core])])
check(acc >= 0.95, sprintf("genes with true kME >= 0.7 assigned to their module: %.3f >= 0.95", acc))
check(length(unique(best)) == 5L, "each planted module maps to a different detected module")

# 3. Module-trait association: planted trait detected, null trait not ------------------
mt <- res$module_trait
target <- paste0("ME", best[as.character(assoc_mod)])
top <- mt[mt$trait == "t_assoc", ][1, ]
check(top$module == target && top$padj < 0.05 && top$r > 0.5,
      sprintf("t_assoc: strongest module is the planted one (%s, r = %.2f, padj = %.1e)", top$module, top$r, top$padj))
check(!any(mt$significant[mt$trait == "t_null"]), "null trait: no module at padj < 0.05")
check(isTRUE(all.equal(mt$padj, stats::p.adjust(mt$pvalue, "BH"))), "BH applied across all module x trait tests")
check(!"ME0" %in% mt$module, "grey (ME0) excluded from module-trait tests")
check(abs(min_detectable_r(24) - 0.4044) < 1e-3, "minimum detectable |r| at n = 24, alpha 0.05 is 0.404")

# 4. Hub genes and kME ------------------------------------------------------------------
hubs <- hub_table(lab, res$kme, colnames(dat), hub_n = 10, hub_kme = 0.8)
h1 <- hubs[hubs$module == as.integer(best[1]), ]
check(all(diff(h1$kME) <= 0) && nrow(h1) == 10L, "hub genes ranked by decreasing kME, 10 per module")
check(mean(truth[match(hubs$gene_id, colnames(dat))] > 0) == 1, "all hub genes are planted module members")

# 5. Soft-threshold selection -------------------------------------------------------------
fit <- data.frame(Power = c(1:10, seq(12, 30, 2)), SFT.R.sq = c(0.1, 0.3, 0.6, 0.85, 0.9, rep(0.92, 15)),
                  slope = c(1, rep(-1.5, 19)))
check(select_power(fit, 30)$power == 4L, "lowest power with signed R^2 >= 0.8 is chosen")
fit_pos <- fit
fit_pos$slope[4] <- 1.2                                            # positive slope: not scale-free
check(select_power(fit_pos, 30)$power == 5L, "a high R^2 with positive slope does not count")
low <- transform(fit, SFT.R.sq = 0.5)
check(select_power(low, 18)$power == 18 && select_power(low, 25)$power == 16 &&
        select_power(low, 35)$power == 14 && select_power(low, 50)$power == 12 && select_power(low, 25)$fallback,
      "no power reaches 0.8: WGCNA FAQ table for signed networks (18/16/14/12)")

# 6. Sample outliers by standardised connectivity ----------------------------------------
# Realistic VST scale: genes differ in mean (6-12) and SD, so samples correlate highly.
d_real <- sweep(sweep(scale(dat), 2, stats::runif(ncol(dat), 0.3, 1.2), "*"), 2, stats::runif(ncol(dat), 6, 12), "+")
d_out <- d_real
d_out[7, ] <- d_out[7, ] + stats::rnorm(ncol(dat), 0, 1.5)           # a degraded, noisy sample
zk <- sample_zk(d_out)
cat(sprintf("  Z.k: degraded sample %.2f, lowest other %.2f; clean data minimum %.2f\n",
            zk[["s07"]], min(zk[names(zk) != "s07"]), min(sample_zk(d_real))))
check(identical(names(zk)[zk < -2.5], "s07"), sprintf("noisy sample s07 flagged (Z.k = %.1f), no other sample", zk[7]))
check(zk[["s07"]] < min(zk[names(zk) != "s07"]) - 3, "the degraded sample sits > 3 SD below every other sample")

# 7. Exported edges equal WGCNA's TOM -------------------------------------------------------
edges <- module_edges(res$net, colnames(dat), lab, threshold = 0.1, max_edges = 500)
tom <- WGCNA::TOMsimilarityFromExpr(dat, corType = "bicor", maxPOutliers = 0.05, networkType = "signed",
                                    power = res$power, TOMType = "signed", verbose = 0)
dimnames(tom) <- list(colnames(dat), colnames(dat))
check(max(abs(edges$weight - tom[cbind(edges$fromNode, edges$toNode)])) < 1e-5,
      "edge weights equal TOMsimilarityFromExpr (dist-index extraction is correct)")
check(all(lab[match(edges$fromNode, colnames(dat))] == lab[match(edges$toNode, colnames(dat))]) &&
        all(edges$weight >= 0.1) && all(edges$fromNode != edges$toNode) &&
        all(table(edges$module) <= 500), "edges are within-module, above threshold, capped per module, no self-loops")

# 8. Batch removal and input checks ----------------------------------------------------------
bdat <- dat
batch <- rep(c("A", "B"), 20)
bdat[batch == "B", ] <- bdat[batch == "B", ] + 2
cb <- remove_batch(bdat, batch)
shift_after <- mean(colMeans(cb[batch == "B", ]) - colMeans(cb[batch == "A", ]))
shift_before <- mean(colMeans(bdat[batch == "B", ]) - colMeans(bdat[batch == "A", ]))
check(abs(shift_after) < 0.1 && abs(shift_before - 2) < 0.2,
      sprintf("ComBat removes a +2 batch shift (mean shift %.2f before, %.3f after)", shift_before, shift_after))
expect_error(remove_batch(bdat, batch, data.frame(trt = ifelse(batch == "A", "c", "d"))),
             "completely confounded", "batch confounded with a protected factor is refused")
expect_error(check_log_scale(matrix(c(0, 5000, 20, 300), 2)), "counts", "count-scale --expr refused")
small <- file.path(tmp, "small.csv")
utils::write.csv(data.frame(gene_id = colnames(dat), t(dat[1:12, ]), check.names = FALSE), small, row.names = FALSE)
expect_error(main(c("--expr", small, "--out", file.path(tmp, "small_out"))), "15 samples",
             "fewer than 15 samples refused (WGCNA FAQ)")

# 9. End to end on the demo --------------------------------------------------------------------
out <- file.path(tmp, "demo_out")
r <- run_cli(c("--demo", "--out", shQuote(out)))
check(r$status == 0L, paste("demo run exits cleanly:", paste(utils::tail(r$output, 5), collapse = "\n")))
for (f in c("report.md", "tables/gene_modules.csv", "tables/module_eigengenes.csv", "tables/module_trait_correlations.csv",
            "tables/hub_genes.csv", "tables/cytoscape_edges.csv", "tables/cytoscape_nodes.csv",
            "tables/soft_threshold.csv", "tables/sample_connectivity.csv", "figures/soft_threshold.png",
            "figures/gene_dendrogram.png", "figures/module_trait.png", "figures/sample_outliers.png",
            "reproducibility/commands.sh", "reproducibility/sessionInfo.txt"))
  check(file.exists(file.path(out, f)), paste("demo writes", f))
gm <- utils::read.csv(file.path(out, "tables", "gene_modules.csv"))
tm <- utils::read.csv(file.path(SKILL_TEST_DIR, "demo", "true_modules.csv"))
m <- merge(gm, tm, by = "gene_id")
check(length(setdiff(unique(gm$module), 0)) == 4L, "demo: the 4 planted modules are found")
check(adjusted_rand(m$true_module[m$true_module > 0], m$module[m$true_module > 0]) >= 0.8,
      "demo: ARI among planted-module genes >= 0.8")
dmt <- utils::read.csv(file.path(out, "tables", "module_trait_correlations.csv"))
sig <- dmt[dmt$significant, ]
check(any(sig$trait == "stomatal_conductance_mol_m2_s" & sig$r < 0) && any(sig$trait == "leaf_n_mg_g" & sig$r > 0),
      "demo: planted stomatal-conductance (negative) and leaf-N (positive) module links found")
check(!any(sig$trait == "plant_height_cm"), "demo: the null trait (plant height) has no significant module")
check(run_cli(c("--demo", "--colour", "red"))$status != 0L, "unknown option exits with an error")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
