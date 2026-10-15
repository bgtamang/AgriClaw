# Tests for claw-yield-prediction. Run from anywhere: Rscript tests/test_yield_prediction.R
# Known-answer checks: metric definitions, fold design properties, no test-fold
# information in training, simulated leakage (random CV must beat grouped CV
# when genotype or environment identity leaks), recovery of a real signal under
# grouped CV, a pure-noise target, seeds, temporal leakage guard, input
# validation, and an end-to-end run of the command line on the demo data.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "yield_prediction.R"))  # loads functions; main() does not run

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}

run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "yield_prediction.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}

err_msg <- function(expr) tryCatch({ force(expr); "" }, error = function(e) conditionMessage(e))

# Light settings keep the file under ~2 minutes; the methods are unchanged.
opt0 <- DEFAULTS
opt0$num_trees <- 100
opt0$alpha_grid <- "0.5,1"
opt0$inner_folds <- 3
opt0$perm_repeats <- 2
opt0$folds <- 3
with_opt <- function(...) { o <- opt0; v <- list(...); o[names(v)] <- v; o }
metric <- function(cv, scheme, model, col) cv$metrics[[col]][cv$metrics$scheme == scheme & cv$metrics$model == model]

# 1. Metric definitions, by hand ---------------------------------------------------------
m <- pred_metrics(c(1, 2, 3, 4), c(2, 2, 2, 2))
check(is.na(m[["r"]]) && abs(m[["rmse"]] - sqrt(1.5)) < 1e-12 && abs(m[["bias"]] + 0.5) < 1e-12,
      "constant prediction: r undefined, RMSE and bias by hand")
check(abs(m[["r2_sse"]] - (1 - 6 / 5)) < 1e-12, "R2 = 1 - SSE/SST is -0.2 for a constant off-centre prediction")
check(abs(m[["rrmse_pct"]] - 100 * sqrt(1.5) / 2.5) < 1e-9, "relative RMSE = RMSE / mean observed")
m <- pred_metrics(1:4, 2 * (1:4))
check(abs(m[["r"]] - 1) < 1e-12 && abs(m[["r2_cor"]] - 1) < 1e-12 && abs(m[["r2_sse"]] - (1 - 30 / 5)) < 1e-12,
      "perfectly correlated but scaled predictions: r^2 = 1 while 1 - SSE/SST = -5")
check(abs(mean(within_group_r(c(1, 2, 3, 10, 11, 12), c(3, 2, 1, 10, 11, 12), rep(c("a", "b"), each = 3))) - 0) < 1e-12,
      "within-group r averages r = -1 and r = +1 to 0")

# 2. Fold design properties ----------------------------------------------------------------
design <- expand.grid(genotype = sprintf("G%02d", 1:12), site = c("IL", "IA"), year = c(2024, 2025),
                      stringsAsFactors = FALSE)
design$env <- paste0(design$site, design$year)
design$plot_id <- seq_len(nrow(design))
o <- with_opt(year_col = "year")
for (s in c("random", "leave-env-out", "leave-genotype-out", "leave-year-out", "leave-genotype-env-out")) {
  f <- make_folds(design, s, o)
  tested <- sort(unlist(lapply(f, `[[`, "test")))
  check(identical(tested, seq_len(nrow(design))), paste(s, ": every row is tested exactly once"))
  check(all(vapply(f, function(x) !length(intersect(x$train, x$test)), TRUE)), paste(s, ": train and test never overlap"))
  if (s %in% c("leave-env-out", "leave-genotype-env-out"))
    check(all(vapply(f, function(x) !length(intersect(design$env[x$train], design$env[x$test])), TRUE)),
          paste(s, ": no environment in both train and test"))
  if (s %in% c("leave-genotype-out", "leave-genotype-env-out"))
    check(all(vapply(f, function(x) !length(intersect(design$genotype[x$train], design$genotype[x$test])), TRUE)),
          paste(s, ": no genotype in both train and test"))
  if (s == "leave-year-out")
    check(all(vapply(f, function(x) !length(intersect(design$year[x$train], design$year[x$test])), TRUE)),
          "leave-year-out: no year in both train and test")
}
ids <- inner_foldid(30, rep(c("a", "b", "c", "d", "e", "f"), each = 5), 3, 1)
check(all(tapply(ids$id, rep(1:6, each = 5), function(v) length(unique(v))) == 1) && grepl("grouped", ids$type),
      "inner folds keep each group together")
check(grepl("random", inner_foldid(10, rep(c("a", "b"), each = 5), 3, 1)$type), "inner folds fall back to random with < 3 groups")
check(identical(available_schemes(design, with_opt(year_col = "")),
                c("leave-genotype-env-out", "leave-env-out", "leave-genotype-out", "random")),
      "available schemes follow the grouping columns present")

# 3. Simulated data ---------------------------------------------------------------------------
# Environment leakage: yield = environment mean + genotype + error; the only
# features are environment and genotype "fingerprints" that carry identity but
# no transferable information.
sim_env_leak <- function(seed, n_env = 6, n_geno = 30) {
  set.seed(seed)
  env <- rep(sprintf("E%d", seq_len(n_env)), each = n_geno)
  geno <- rep(sprintf("G%02d", seq_len(n_geno)), times = n_env)
  ie <- match(env, unique(env))
  ig <- match(geno, unique(geno))
  y <- 4000 + rnorm(n_env, 0, 800)[ie] + rnorm(n_geno, 0, 300)[ig] + rnorm(length(ie), 0, 250)
  x <- cbind(matrix(rnorm(n_env * 3), n_env)[ie, ] + rnorm(length(ie) * 3, 0, 0.05),
             matrix(rnorm(n_geno * 3), n_geno)[ig, ] + rnorm(length(ie) * 3, 0, 0.05))
  colnames(x) <- c(paste0("env_fp", 1:3), paste0("geno_fp", 1:3))
  data.frame(plot_id = seq_along(y), env = env, genotype = geno, yield_kg_ha = y, x)
}
# Genotype leakage: replicated genotypes, heritable traits that identify the
# genotype but do not cause yield.
sim_geno_leak <- function(seed, n_geno = 40, reps = 4, dims = 6) {
  set.seed(seed)
  geno <- rep(sprintf("G%02d", seq_len(n_geno)), times = reps)
  ig <- match(geno, unique(geno))
  n <- length(ig)
  y <- 4000 + rnorm(n_geno, 0, 400)[ig] + rnorm(n, 0, 200)
  x <- matrix(rnorm(n_geno * dims), n_geno)[ig, ] + rnorm(n * dims, 0, 0.1)
  colnames(x) <- paste0("trait", seq_len(dims))
  data.frame(plot_id = seq_len(n), genotype = geno, yield_kg_ha = y, x)
}
# Real signal: yield depends on two plot-level features the same way in every
# environment and genotype; eight pure-noise features.
sim_signal <- function(seed, n_env = 6, n_geno = 30, n_noise = 8) {
  set.seed(seed)
  env <- rep(sprintf("E%d", seq_len(n_env)), each = n_geno)
  geno <- rep(sprintf("G%02d", seq_len(n_geno)), times = n_env)
  ie <- match(env, unique(env))
  ig <- match(geno, unique(geno))
  n <- length(ie)
  x1 <- rnorm(n_geno)[ig] + rnorm(n_env, 0, 0.7)[ie] + rnorm(n, 0, 0.4)
  x2 <- 0.5 * rnorm(n_geno)[ig] + rnorm(n_env)[ie] + rnorm(n, 0, 0.4)
  signal <- 400 * x1 + 300 * x2
  y <- 4000 + signal + rnorm(n, 0, 250)
  noise <- matrix(rnorm(n * n_noise), n)
  colnames(noise) <- paste0("noise", seq_len(n_noise))
  list(d = data.frame(plot_id = seq_len(n), env = env, genotype = geno, yield_kg_ha = y,
                      x1 = x1, x2 = x2, noise), oracle_r = stats::cor(y, signal))
}

# 4. Leakage: random CV reports clearly higher accuracy than grouped CV ----------------------------
d <- sim_env_leak(1)
cv <- run_cv(prepare_data(d, opt0), with_opt(cv = "random,leave-env-out", models = "rf"))
r_rand <- metric(cv, "random", "rf", "r")
r_loeo <- metric(cv, "leave-env-out", "rf", "r")
check(r_rand >= 0.7, sprintf("environment fingerprints: random CV r = %.2f looks good", r_rand))
r2_rand <- metric(cv, "random", "rf", "r2_sse")
r2_loeo <- metric(cv, "leave-env-out", "rf", "r2_sse")
check(r_rand - r_loeo >= 0.5 && r2_rand - r2_loeo >= 0.5 && r2_loeo < 0.2,
      sprintf("environment fingerprints: leave-env-out exposes it (r %.2f vs %.2f, R2 %.2f vs %.2f)", r_loeo, r_rand, r2_loeo, r2_rand))
w <- cv_warnings(cv, prepare_data(d, opt0), with_opt(cv = "random,leave-env-out"))
check(any(grepl("Random k-fold overstates accuracy", w)), "report warns that random k-fold overstates accuracy")
check(any(grepl("GWAS", w)), "negative R2 triggers the do-not-use-as-phenotype warning")

d <- sim_geno_leak(1)
cv <- run_cv(prepare_data(d, with_opt(env_col = "")), with_opt(env_col = "", cv = "random,leave-genotype-out", models = "rf"))
r_rand <- metric(cv, "random", "rf", "r")
r_lgo <- metric(cv, "leave-genotype-out", "rf", "r")
check(r_rand >= 0.5 && r_lgo <= 0.2 && r_rand - r_lgo >= 0.4,
      sprintf("genotype fingerprints: random r %.2f vs leave-genotype-out r %.2f", r_rand, r_lgo))
check(metric(cv, "leave-genotype-out", "rf", "r2_sse") < 0, "genotype fingerprints: new genotypes are not predictable (R2 < 0)")

# 5. A real signal survives grouped CV -----------------------------------------------------------
sg <- sim_signal(1)
# Tolerances are relative to the oracle (the true signal itself): r within
# 0.05 and R2 within 0.1 of oracle r^2 for the elastic net (the true model is
# linear); R2 within 0.25 for the random forest.
oracle_r2 <- sg$oracle_r^2
cv <- run_cv(prepare_data(sg$d, opt0), with_opt(cv = "leave-env-out"))
s <- "leave-env-out"
check(metric(cv, s, "enet", "r") >= sg$oracle_r - 0.05 && metric(cv, s, "enet", "r2_sse") >= oracle_r2 - 0.1,
      sprintf("%s: elastic net r %.3f, R2 %.2f vs oracle r %.3f, R2 %.2f", s, metric(cv, s, "enet", "r"),
              metric(cv, s, "enet", "r2_sse"), sg$oracle_r, oracle_r2))
check(metric(cv, s, "rf", "r2_sse") >= oracle_r2 - 0.25,
      sprintf("%s: random forest R2 %.2f vs oracle %.2f", s, metric(cv, s, "rf", "r2_sse"), oracle_r2))
cv2 <- run_cv(prepare_data(sg$d, opt0), with_opt(cv = "leave-genotype-env-out", models = "enet"))
s <- "leave-genotype-env-out"
check(metric(cv2, s, "enet", "r") >= sg$oracle_r - 0.05 && metric(cv2, s, "enet", "r2_sse") >= oracle_r2 - 0.1,
      sprintf("%s: elastic net r %.3f, R2 %.2f vs oracle r %.3f, R2 %.2f", s, metric(cv2, s, "enet", "r"),
              metric(cv2, s, "enet", "r2_sse"), sg$oracle_r, oracle_r2))
imp <- cv$importance[cv$importance$scheme == cv$primary, ]
for (mm in c("enet", "rf"))
  check(setequal(imp$feature[imp$model == mm & imp$rank <= 2], c("x1", "x2")),
        paste(mm, ": permutation importance ranks the two informative features first"))

# 6. Pure-noise target: R2 (1 - SSE/SST) is zero or below ---------------------------------------
d <- sim_env_leak(2)
set.seed(99)
d$yield_kg_ha <- rnorm(nrow(d), 4000, 500)
cv <- run_cv(prepare_data(d, opt0), with_opt(cv = "random,leave-env-out"))
check(all(cv$metrics$r2_sse <= 0), sprintf("pure noise: every R2 <= 0 (max %.3f)", max(cv$metrics$r2_sse)))
check(all(cv$metrics$flag[cv$metrics$model != "train_mean"] == "worse than predicting the mean"),
      "pure noise: models flagged as worse than predicting the mean")

# 7. No test-fold information reaches training (scaling, tuning, fitting) -----------------------
sg2 <- sim_signal(2, n_env = 4, n_geno = 20)
d1 <- sg2$d
f1 <- make_folds(d1, "leave-env-out", opt0)[[1]]
d2 <- d1
feat <- setdiff(names(d1), c("plot_id", "env", "genotype", "yield_kg_ha"))
d2[f1$test, feat] <- d2[f1$test, feat] * 100 + 50
d2$yield_kg_ha[f1$test] <- d2$yield_kg_ha[f1$test] + 1e5
o <- with_opt(perm_repeats = 0)
k1 <- run_scheme(d1, feat, "leave-env-out", o, keep_models = TRUE)$models
k2 <- run_scheme(d2, feat, "leave-env-out", o, keep_models = TRUE)$models
lab <- paste("enet", f1$label)
check(identical(k1[[lab]]$coef, k2[[lab]]$coef), "elastic net fitted without the held-out rows (coefficients unchanged)")
probe <- as.matrix(d1[f1$train, feat])
check(identical(k1[[paste("rf", f1$label)]]$predict(probe), k2[[paste("rf", f1$label)]]$predict(probe)),
      "random forest fitted without the held-out rows (predictions unchanged)")
check(!identical(k1[[paste("enet", make_folds(d1, "leave-env-out", opt0)[[2]]$label)]]$coef,
                 k2[[paste("enet", make_folds(d1, "leave-env-out", opt0)[[2]]$label)]]$coef),
      "control: folds that train on the changed rows do change")

# 8. Seeds make results reproducible ---------------------------------------------------------------
p1 <- run_cv(prepare_data(d1, opt0), with_opt(cv = "random", seed = 7))$pred
p2 <- run_cv(prepare_data(d1, opt0), with_opt(cv = "random", seed = 7))$pred
p3 <- run_cv(prepare_data(d1, opt0), with_opt(cv = "random", seed = 8))$pred
check(identical(p1, p2), "same seed: identical predictions")
check(!identical(p1$fold, p3$fold) || !identical(p1$predicted, p3$predicted), "different seed: different random folds")

# 9. Temporal leakage guard ---------------------------------------------------------------------------
ft <- data.frame(feature = c("x1", "x2", paste0("noise", 1:8)), available_t = c(45, 110, rep(60, 8)))
pr <- prepare_data(d1, with_opt(cutoff = 60), ft)
check(!"x2" %in% pr$features && "x1" %in% pr$features, "--cutoff 60 drops a feature known only at DAP 110")
check(grepl("cutoff", pr$feature_log$reason[pr$feature_log$feature == "x2"]), "the dropped feature is logged with its reason")
check(grepl("needs --feature-time", err_msg(prepare_data(d1, with_opt(cutoff = 60)))), "--cutoff without feature times is refused")
check(grepl("No available_t", err_msg(prepare_data(d1, with_opt(cutoff = 60, features = "x1,x2"), ft[-1, ]))),
      "a listed feature without an available time is refused")
pr <- prepare_data(d1, with_opt(cutoff = 60), ft[-1, ])
check(!"x1" %in% pr$features && any(grepl("x1", pr$notes)), "columns absent from the feature-time table are reported, not used")

# 10. Input validation and warnings -------------------------------------------------------------------
check(grepl("not found", err_msg(prepare_data(d1, with_opt(target = "yield_t_ha")))), "missing target column is named")
dna <- d1
dna$x1[3] <- NA
check(grepl("missing feature values", err_msg(prepare_data(dna, opt0))), "missing features stop the run (nothing imputed)")
pr <- prepare_data(dna, with_opt(drop_na_rows = TRUE))
check(nrow(pr$dat) == nrow(d1) - 1 && any(grepl("dropped", pr$notes)), "--drop-na-rows drops and reports the row")
ddup <- d1
ddup$plot_id[2] <- ddup$plot_id[1]
check(grepl("repeated identifiers", err_msg(prepare_data(ddup, opt0))), "duplicate plot ids rejected")
check(grepl("not possible", err_msg(run_cv(prepare_data(d1, opt0), with_opt(cv = "leave-year-out")))),
      "a scheme without its grouping column is refused")
dw <- d1
dw$row <- rep(1:10, length.out = nrow(dw))
dw$yield_copy <- dw$yield_kg_ha + rnorm(nrow(dw), 0, 1)
dw$rain_mm <- c(300, 350, 400, 450)[match(dw$env, unique(dw$env))]
names(dw)[names(dw) == "yield_kg_ha"] <- "yield"
ww <- data_warnings(prepare_data(dw, with_opt(target = "yield")), with_opt(target = "yield"))
check(any(grepl("field-position", ww)), "warns about a field-position feature")
check(any(grepl("target leakage", ww)), "warns about a feature that copies the target")
check(any(grepl("constant within each environment", ww)), "warns about environment-level features")
check(any(grepl("no units", ww)), "warns when the target has no units")

# 11. End to end on the demo data ---------------------------------------------------------------------------
out <- file.path(tempdir(), "claw-yield-prediction-test")
unlink(out, recursive = TRUE)
res <- run_cli(c("--demo", "--cutoff", "60", "--cv", "random,leave-env-out", "--num-trees", "60",
                 "--alpha-grid", "0.5,1", "--inner-folds", "3", "--perm-repeats", "2", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "tables/cv_metrics.csv", "tables/cv_predictions.csv", "tables/fold_metrics.csv",
            "tables/tuning.csv", "tables/permutation_importance.csv", "tables/features.csv",
            "figures/accuracy_by_scheme.png", "figures/obs_vs_pred.png", "figures/importance.png",
            "reproducibility/commands.sh", "reproducibility/sessionInfo.txt")) {
  check(file.exists(file.path(out, f)), paste("demo writes", f))
}
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
met <- utils::read.csv(file.path(out, "tables", "cv_metrics.csv"))
check(nrow(met) == 6L && setequal(met$scheme, c("random", "leave-env-out")), "demo: 2 schemes x (enet, rf, train_mean)")
pred <- utils::read.csv(file.path(out, "tables", "cv_predictions.csv"))
check(nrow(pred) == 240L * 6L, "demo: every plot predicted once per scheme and model")
fl <- utils::read.csv(file.path(out, "tables", "features.csv"))
check(setequal(fl$feature[fl$used], c("greenup_dap", "ndvi_at_dap50", "ndvi_at_dap60")),
      "demo --cutoff 60 keeps only features known by DAP 60")
rep_txt <- readLines(file.path(out, "report.md"))
check(any(grepl("## Warnings", rep_txt)) && any(grepl("seed 1", rep_txt)), "report has warnings and records the seed")
check(run_cli(c("--demo", "--colour", "red"))$status != 0L, "unknown option exits with an error")
check(any(grepl("--feature-time", run_cli("--help")$output)), "--help prints usage")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
