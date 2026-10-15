# Tests for claw-genomic-prediction. Run from anywhere: Rscript tests/test_genomic_prediction.R
# Known answers: the Daetwyler et al. (2008) expectation for GBLUP accuracy on
# simulated unlinked loci, equivalence of GBLUP and RR-BLUP, the A.mat
# definition A = WW'/c, true genetic values of the simulated demo, and the
# expected ranking of Bayesian models under a sparse architecture.
# Design properties: every CV scheme keeps test and training apart as defined.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "genomic_prediction.R"))  # loads functions only

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}

run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "genomic_prediction.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}

dd <- file.path(SKILL_TEST_DIR, "demo")

# 1. Marker coding.
X3 <- matrix(c(0, 1, 2, 2, 0, NA), 2, dimnames = list(c("a", "b"), c("m1", "m2", "m3")))
check(identical(code_markers(X3, "auto")$coding, "012"), "auto-detects 0/1/2 dosage")
check(isTRUE(all.equal(code_markers(X3, "012")$M, X3 - 1)), "0/1/2 recoded to -1/0/1")
check(isTRUE(all.equal(code_markers(X3 - 1)$M, X3 - 1)), "-1/0/1 kept as is")
X01 <- matrix(c(0, 1, 1, 0), 2)
err <- tryCatch(code_markers(X01, "auto"), error = conditionMessage)
check(grepl("ambiguous", err), "0/1-only matrix with auto coding fails and asks for --geno-coding")
check(isTRUE(all.equal(code_markers(X01, "01")$M, 2 * X01 - 1)), "0/1 inbred coding recoded to -1/1")
err <- tryCatch(code_markers(X3 + 5, "012"), error = conditionMessage)
check(grepl("outside", err), "out-of-range genotype values are rejected")

# 2. Relationship matrix: A.mat equals W W' / c with c = 2 sum p(1-p)
#    (VanRaden 2008; rrBLUP documentation), monomorphic markers dropped,
#    missing calls mean-imputed and counted.
set.seed(1)
Msim <- matrix(sample(c(-1, 1), 80 * 60, replace = TRUE), 80, 60,
               dimnames = list(sprintf("L%02d", 1:80), sprintf("m%02d", 1:60)))
Msim[, 3] <- 1
Msim[sample(length(Msim), 50)] <- NA
rl <- build_relationship(Msim)
imp <- Msim[, -3]
for (j in seq_len(ncol(imp))) imp[is.na(imp[, j]), j] <- mean(imp[, j], na.rm = TRUE)
p <- (colMeans(imp) + 1) / 2
W_hand <- sweep(imp, 2, colMeans(imp))
check(max(abs(rl$A - tcrossprod(W_hand) / (2 * sum(p * (1 - p))))) < 1e-10,
      "A = WW'/c from hand-imputed markers")
check(max(abs(rl$W - W_hand)) < 1e-12, "W = mean-imputed, centred markers")
check(rl$n_markers_used == 59L && rl$n_imputed == 50L, "monomorphic marker removed, 50 imputed calls counted")

# 3. Demo data, and GBLUP == RR-BLUP (equivalent models).
X <- read_genotypes(file.path(dd, "genotypes.csv"))
rel <- build_relationship(code_markers(X)$M)
ph <- read_phenotypes(file.path(dd, "phenotypes.csv"), "id", "yield_kg_ha", "family")$data
truth <- utils::read.csv(file.path(dd, "truth.csv"))
g_true <- stats::setNames(truth$true_genetic_value, truth$id)
y <- stats::setNames(ph$y, ph$id)
set.seed(2)
te <- sample(names(y), 48)
tr <- setdiff(names(y), te)
pg <- fit_predict("gblup", y[tr], te, rel)
pr <- fit_predict("rrblup", y[tr], te, rel)
check(stats::cor(pg, pr) > 0.99999 && max(abs(pg - pr)) < 1e-4 * stats::sd(y),
      "GBLUP and RR-BLUP predictions are identical")
check(abs(attr(pg, "vc")[["Vg"]] / (rel$c * attr(pr, "vc")[["Vu"]]) - 1) < 1e-4,
      "REML genetic variance: Vg(GBLUP) = c * Vu(RR-BLUP)")

# 4. Fold designs.
set.seed(3)
lab <- assign_kfold(103, 5)
check(all(table(lab) %in% c(20L, 21L)), "k-fold: balanced fold sizes")
fam <- rep(c("A", "B", "C", "D"), c(30, 25, 12, 7))
lab <- assign_family_stratified(fam, 5)
check(all(apply(table(fam, lab), 1, function(x) max(x) - min(x) <= 1)),
      "family-stratified: every family spread evenly over folds")
me <- utils::read.csv(file.path(dd, "phenotypes_multienv.csv"))
lab2 <- assign_cv2(me$id, 5)
check(all(tapply(lab2, me$id, function(x) !anyDuplicated(x))),
      "CV2: records of one line fall in different folds")
lines <- unique(me$id)
lf <- stats::setNames(assign_kfold(length(lines), 5), lines)
check(all(tapply(lf[me$id], me$id, function(x) length(unique(x)) == 1L)),
      "CV1: all records of a line share one fold")

# 5. Leave-one-family-out keeps families apart and lowers test-training relatedness.
cvd <- run_cv_single(ph, rel, "rrblup", c("kfold", "lofo"), 5, 1, 11L, list())
lo <- cvd$preds[cvd$preds$scheme == "lofo", ]
check(all(tapply(lo$family, lo$fold, function(x) length(unique(x)) == 1L)) &&
        all(lo$family == lo$fold), "LOFO: each test fold is exactly one family")
fd <- cvd$folds
check(stats::median(fd$median_max_relatedness[fd$scheme == "lofo"]) <
        stats::median(fd$median_max_relatedness[fd$scheme == "kfold"]),
      "LOFO test lines are less related to training lines than random k-fold")
acc_true <- function(d) mean(vapply(split(d, d$fold), function(z) stats::cor(g_true[z$id], z$predicted), 0))
a_k <- acc_true(cvd$preds[cvd$preds$scheme == "kfold", ])
a_l <- acc_true(lo)
check(a_k > a_l, sprintf("true accuracy is higher with relatives in training (%.2f vs %.2f)", a_k, a_l))
cat(sprintf("Demo true accuracy: k-fold %.3f, LOFO %.3f\n", a_k, a_l))
est_k <- mean(fd$r[fd$scheme == "kfold"]) / sqrt(0.5)
check(abs(est_k - a_k) < 0.1,
      sprintf("demo: r/sqrt(h2) = %.3f is within 0.1 of the true accuracy %.3f", est_k, a_k))
check(all(table(cvd$preds$id[cvd$preds$scheme == "kfold"]) == 1L) &&
        length(unique(cvd$preds$id[cvd$preds$scheme == "kfold"])) == nrow(ph),
      "k-fold: every line is predicted exactly once per repetition (never from its own record)")

# 6. Reproducibility with a seed.
cva <- run_cv_single(ph, rel, "rrblup", "kfold", 3, 1, 99L, list())
cvb <- run_cv_single(ph, rel, "rrblup", "kfold", 3, 1, 99L, list())
cvc <- run_cv_single(ph, rel, "rrblup", "kfold", 3, 1, 100L, list())
check(identical(cva$folds, cvb$folds) && identical(cva$preds, cvb$preds), "same seed, same result")
fold_of <- function(cv) cv$preds$fold[order(cv$preds$id)]
check(!identical(fold_of(cva), fold_of(cvc)), "different seed, different fold assignment")

# 7. Daetwyler et al. (2008), eq. 1: expected accuracy r1 = sqrt(L h2 / (L h2 + 1)),
#    L = N / M, for M independent loci. The authors call it conservative; their
#    second-order correction (Appendix S1) raises it by a fraction of about
#    r1^2 / (2 L). GBLUP accuracy on unlinked loci (all markers causal, h2 = 0.5,
#    M = 400; computed with the equivalent RR-BLUP solver, see section 3) for 200
#    unrelated test lines, 3 replicates per training size, must fall
#    between r1 - 0.05 and r1 * (1 + r1^2 / (2 L)) + 0.05 (0.05 = Monte Carlo margin,
#    about 2 SE of a 3-replicate mean).
sim_acc <- function(N, M = 400, h2 = 0.5, n_test = 200, seed = 1) {
  set.seed(seed)
  n <- N + n_test
  pf <- stats::runif(M, 0.1, 0.9)
  Xs <- matrix(ifelse(stats::runif(n * M) < rep(pf, each = n), 1, -1), n, M,
               dimnames = list(sprintf("S%04d", seq_len(n)), sprintf("q%03d", seq_len(M))))
  g <- drop(scale(Xs, scale = FALSE) %*% stats::rnorm(M))
  g <- g / stats::sd(g) * sqrt(h2)
  ys <- stats::setNames(g + stats::rnorm(n, 0, sqrt(1 - h2)), rownames(Xs))
  rs <- build_relationship(Xs)
  trs <- rownames(Xs)[seq_len(N)]
  tes <- rownames(Xs)[N + seq_len(n_test)]
  pp <- fit_predict("rrblup", ys[trs], tes, rs)
  c(acc = stats::cor(g[tes], pp), est = stats::cor(ys[tes], pp) / sqrt(h2))
}
sizes <- c(120, 240, 480)
res <- lapply(sizes, function(N) t(vapply(1:3, function(s) sim_acc(N, seed = 100 * N + s), c(acc = 0, est = 0))))
acc <- vapply(res, function(r) mean(r[, "acc"]), 0)
est <- vapply(res, function(r) mean(r[, "est"]), 0)
lam <- sizes / 400
r1 <- sqrt(lam * 0.5 / (lam * 0.5 + 1))
r2 <- r1 * (1 + r1^2 / (2 * lam))
for (j in seq_along(sizes))
  check(acc[j] > r1[j] - 0.05 && acc[j] < r2[j] + 0.05,
        sprintf("N = %d: GBLUP accuracy %.3f within Daetwyler first/second-order range [%.3f, %.3f]",
                sizes[j], acc[j], r1[j], r2[j]))
check(all(diff(acc) > 0), "accuracy increases with training size")
d_est <- unlist(lapply(res, function(r) r[, "est"] - r[, "acc"]))
check(abs(mean(d_est)) < 0.06,
      sprintf("r(y, y_hat)/sqrt(h2) is an unbiased estimate of accuracy (mean error %.3f over 9 runs)",
              mean(d_est)))
cat(sprintf("Daetwyler check: N = %s; GBLUP accuracy %s; first-order %s; second-order %s; mean error of r/sqrt(h2) %.3f\n",
            paste(sizes, collapse = "/"), paste(sprintf("%.3f", acc), collapse = "/"),
            paste(sprintf("%.3f", r1), collapse = "/"), paste(sprintf("%.3f", r2), collapse = "/"), mean(d_est)))

# 8. Bayesian models (BGLR, 1500 iterations, 500 burn-in). Bayesian ridge
#    regression is the same model as RR-BLUP, so predictions must agree. For a
#    sparse trait (5 QTL among 500 markers, h2 = 0.5), variable-selection priors
#    must beat GBLUP (Meuwissen et al. 2001; Daetwyler et al. 2010): mean
#    accuracy over two train/test splits, BayesB at least 0.03 above GBLUP
#    (development runs: about 0.1 above), BayesC above GBLUP.
set.seed(101)
W <- rel$W
q <- sample(ncol(W), 5)
g5 <- drop(W[, q] %*% stats::rnorm(5))
y5 <- stats::setNames(g5 + stats::rnorm(length(g5), 0, stats::sd(g5)), rownames(W))
bg <- list(n_iter = 1500, burn_in = 500, thin = 5, dir = tempfile("bglr"))
dir.create(bg$dir)
acc5 <- NULL
for (s in 1:2) {
  te5 <- names(y5)[sample(length(y5), 70)]
  tr5 <- setdiff(names(y5), te5)
  p_gb <- fit_predict("gblup", y5[tr5], te5, rel)
  p_bb <- fit_predict("bayesb", y5[tr5], te5, rel, bg, seed = s)
  p_bc <- fit_predict("bayesc", y5[tr5], te5, rel, bg, seed = s)
  acc5 <- rbind(acc5, c(gblup = stats::cor(g5[te5], p_gb), bayesb = stats::cor(g5[te5], p_bb),
                        bayesc = stats::cor(g5[te5], p_bc)))
  if (s == 1L) {
    p_rr <- fit_predict("rrblup", y5[tr5], te5, rel)
    p_brr <- fit_predict("brr", y5[tr5], te5, rel, bg, seed = s)
    check(stats::cor(p_brr, p_rr) > 0.99, "BRR (BGLR) agrees with RR-BLUP")
    dg <- bglr_diagnostics(attr(p_bb, "prefix"), 500, 5)
    check(length(dg$trace) == 300L && is.finite(dg$ess),
          "BGLR residual-variance chain read and summarised")
  }
}
m5 <- colMeans(acc5)
check(m5[["bayesb"]] > m5[["gblup"]] + 0.03, "BayesB beats GBLUP for a sparse (5-QTL) trait")
check(m5[["bayesc"]] > m5[["gblup"]], "BayesC beats GBLUP for a sparse (5-QTL) trait")
cat(sprintf("5-QTL trait accuracy: GBLUP %.3f, BayesB %.3f, BayesC %.3f; cor(BRR, RR-BLUP) %.4f\n",
            m5[["gblup"]], m5[["bayesb"]], m5[["bayesc"]], stats::cor(p_brr, p_rr)))
unlink(bg$dir, recursive = TRUE)

# 9. Multi-environment CV1 / CV2 on a subset (3 families, 3 environments).
mesub <- me[me$family %in% c("F1", "F2", "F3"), ]
names(mesub)[names(mesub) == "yield_kg_ha"] <- "y"
cvm <- run_cv_multienv(mesub, rel, c("cv1", "cv2"), 3, 1, 21L)
fm <- cvm$folds
check(all(fm$test_lines_seen_in_training[fm$scheme == "cv1"] == 0), "CV1: test lines never in training")
check(all(fm$test_lines_seen_in_training[fm$scheme == "cv2"] == 1), "CV2: test lines observed in other environments")
check(mean(fm$r[fm$scheme == "cv2"]) > mean(fm$r[fm$scheme == "cv1"]),
      "CV2 predictive ability exceeds CV1 when environments are genetically correlated")
cat(sprintf("CV1 r = %.3f, CV2 r = %.3f\n", mean(fm$r[fm$scheme == "cv1"]), mean(fm$r[fm$scheme == "cv2"])))

# 10. Bad inputs fail loudly.
err <- tryCatch(main(c("--geno", file.path(dd, "genotypes.csv"), "--pheno",
                       file.path(dd, "phenotypes.csv"), "--cv", "lofo", "--out", tempfile())),
                error = conditionMessage)
check(grepl("family-col", err), "LOFO without --family-col is rejected")
err <- tryCatch(read_phenotypes(file.path(dd, "phenotypes_multienv.csv"), "id", "yield_kg_ha"),
                error = conditionMessage)
check(grepl("Duplicated line IDs", err), "repeated IDs without --env-col are rejected")
tg <- tempfile(fileext = ".csv")
utils::write.csv(data.frame(marker = colnames(X)[1:20], t(X[1:30, 1:20]), check.names = FALSE),
                 tg, row.names = FALSE)
err <- tryCatch(main(c("--geno", tg, "--pheno", file.path(dd, "phenotypes.csv"),
                       "--trait", "yield_kg_ha", "--out", tempfile())), error = conditionMessage)
check(grepl("geno-transposed", err), "transposed genotype file is detected")

# 11. End to end on the demo data (one repetition to keep the test fast).
out <- file.path(tempdir(), "claw-genomic-prediction-test")
unlink(out, recursive = TRUE)
res <- run_cli(c("--demo", "--folds", "3", "--reps", "1", "--cv", "kfold", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "tables/cv_folds.csv", "tables/cv_summary.csv",
            "tables/cv_predictions.csv", "tables/predictions_all_lines.csv",
            "tables/variance_components.csv", "figures/cv_predictive_ability.png",
            "figures/cv_observed_vs_predicted.png", "reproducibility/commands.sh",
            "reproducibility/sessionInfo.txt")) {
  check(file.exists(file.path(out, f)), paste("demo writes", f))
}
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L,
      "demo writes a checksum file")
s <- utils::read.csv(file.path(out, "tables", "cv_summary.csv"))
check(nrow(s) == 2L && all(s$predictive_ability_mean > 0.3),
      "demo summary: 2 models, useful k-fold predictive ability")
a <- utils::read.csv(file.path(out, "tables", "predictions_all_lines.csv"))
check(sum(!a$phenotyped) == 30L && all(is.finite(a$predicted_gblup)),
      "predictions for the 30 unphenotyped candidates")
rep_txt <- readLines(file.path(out, "report.md"))
check(any(grepl("random k-fold", rep_txt)) && any(grepl("never reported as accuracy", rep_txt)),
      "report explains the scheme and states that training-set fit is not an accuracy")

# 12. Unknown options are rejected.
bad <- run_cli(c("--demo", "--colour", "red"))
check(bad$status != 0L, "unknown option exits with an error")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
