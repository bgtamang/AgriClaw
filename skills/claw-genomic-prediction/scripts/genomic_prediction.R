#!/usr/bin/env Rscript
# claw-genomic-prediction: genomic prediction (GBLUP, RR-BLUP with rrBLUP;
# optional BayesB, BayesC and Bayesian ridge regression with BGLR) evaluated by
# leakage-aware cross-validation: repeated k-fold, family-stratified k-fold,
# leave-one-family-out, and CV1/CV2 for multi-environment data.
#
# Usage:
#   Rscript scripts/genomic_prediction.R --demo
#   Rscript scripts/genomic_prediction.R --geno geno.csv --pheno blues.csv \
#       --trait yield_kg_ha --family-col family --cv kfold,lofo --h2 0.6
#   Rscript scripts/genomic_prediction.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-genomic-prediction"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/genomic_prediction.R [--demo | --demo-multienv]
            --geno FILE.csv --pheno FILE.csv [--trait COL] [--id-col COL]
            [--family-col COL] [--env-col COL] [--geno-coding auto|012|-101|01]
            [--geno-transposed] [--models gblup,rrblup,bayesb,bayesc,brr]
            [--cv kfold,family,lofo,cv1,cv2] [--folds K] [--reps R] [--h2 H2]
            [--n-iter N] [--burn-in B] [--thin T] [--seed S] [--out DIR]

  --demo            run on demo/ (single-environment BLUEs, kfold + family + lofo)
  --demo-multienv   run on demo/ multi-environment BLUEs (cv1 + cv2, GBLUP)
  --geno            genotype CSV: first column = line ID, one column per marker
                    (use --geno-transposed if markers are rows and lines columns)
  --geno-coding     auto (default), 012 (allele dosage), -101 (rrBLUP), 01 (inbred/dominant)
  --pheno           phenotype CSV: one row per line (BLUEs) or, with --env-col, one
                    row per line x environment
  --trait           phenotype column (default: the only numeric column)
  --id-col          line ID column in --pheno (default: id)
  --family-col      family / cross column (needed for --cv family and lofo)
  --env-col         environment column (needed for --cv cv1 and cv2)
  --models          comma list (default gblup,rrblup); bayesb, bayesc, brr use BGLR
  --cv              comma list of schemes (default kfold)
  --folds           folds for k-fold schemes (default 5)
  --reps            repetitions of k-fold schemes (default 5)
  --h2              entry-mean heritability of the phenotypes, for accuracy = r/sqrt(h2)
  --n-iter, --burn-in, --thin   BGLR sampler settings (default 12000, 2000, 5)
  --seed            random seed (default 2026)
  --out             output folder (default: agriclaw_output/claw-genomic-prediction)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, demo_multienv = FALSE, geno = "", pheno = "",
                 trait = "", id_col = "id", family_col = "", env_col = "",
                 geno_coding = "auto", geno_transposed = FALSE,
                 models = "gblup,rrblup", cv = "kfold", folds = 5, reps = 5,
                 h2 = 0, n_iter = 12000, burn_in = 2000, thin = 5, seed = 2026,
                 out = "")

# ---------------------------------------------------------------------------
# Helpers shared by all AgriClaw R scripts (copy unchanged)
# ---------------------------------------------------------------------------

parse_args <- function(argv, defaults = DEFAULTS) {
  opt <- defaults
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (!startsWith(a, "--"))
      stop("Unexpected argument '", a, "'. See --help.", call. = FALSE)
    key <- gsub("-", "_", substring(a, 3L), fixed = TRUE)
    if (!key %in% names(defaults))
      stop("Unknown option '", a, "'. See --help.", call. = FALSE)
    if (is.logical(defaults[[key]])) {
      opt[[key]] <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(argv)) stop("Option '", a, "' needs a value.", call. = FALSE)
    val <- argv[[i + 1L]]
    if (is.numeric(defaults[[key]])) {
      num <- suppressWarnings(as.numeric(val))
      if (is.na(num))
        stop("Option '", a, "' needs a number, got '", val, "'.", call. = FALSE)
      val <- num
    }
    opt[[key]] <- val
    i <- i + 2L
  }
  opt
}

# Path of this script, whether it is run with Rscript or loaded with source().
this_script <- function() {
  for (i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if (!is.null(f)) return(normalizePath(f, winslash = "/"))
  }
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) return(normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/"))
  NA_character_
}
SCRIPT_FILE <- this_script()
SKILL_DIR <- if (is.na(SCRIPT_FILE)) getwd() else dirname(dirname(SCRIPT_FILE))

write_utf8 <- function(lines, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con))
  writeLines(lines, con)
}

save_table <- function(df, out_dir, name) {
  d <- file.path(out_dir, "tables")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(d, name)
  utils::write.csv(df, path, row.names = FALSE)
  path
}

# Open a PNG device in <out>/figures; the caller draws and then calls dev.off().
open_png <- function(out_dir, name, width = 7, height = 5) {
  d <- file.path(out_dir, "figures")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  grDevices::png(file.path(d, name), width = width, height = height,
                 units = "in", res = 150)
}

file_hash <- function(path, algo) {
  if (algo == "sha256") {
    con <- file(path, open = "rb")
    on.exit(close(con))
    return(as.character(openssl::sha256(con)))
  }
  unname(tools::md5sum(path))
}

# commands.sh, sessionInfo.txt and checksums of inputs + outputs.
write_repro_bundle <- function(out_dir, argv, inputs = character()) {
  rd <- file.path(out_dir, "reproducibility")
  dir.create(rd, recursive = TRUE, showWarnings = FALSE)
  script <- if (is.na(SCRIPT_FILE)) paste0("scripts/", SKILL, ".R") else SCRIPT_FILE
  cmd <- paste(c("Rscript", shQuote(script, type = "sh"),
                 vapply(argv, shQuote, "", type = "sh")), collapse = " ")
  write_utf8(c("#!/usr/bin/env bash",
               sprintf("# %s %s, run %s", SKILL, VERSION,
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
               cmd), file.path(rd, "commands.sh"))
  write_utf8(utils::capture.output(utils::sessionInfo()),
             file.path(rd, "sessionInfo.txt"))
  outs <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  outs <- outs[!startsWith(normalizePath(outs, winslash = "/"),
                           normalizePath(rd, winslash = "/"))]
  ins <- inputs[nzchar(inputs) & file.exists(inputs)]
  files <- c(ins, outs)
  algo <- if (requireNamespace("openssl", quietly = TRUE)) "sha256" else "md5"
  hashes <- vapply(files, file_hash, "", algo = algo)
  base <- normalizePath(out_dir, winslash = "/")
  labels <- normalizePath(files, winslash = "/")
  inside <- startsWith(labels, paste0(base, "/"))
  labels[inside] <- substring(labels[inside], nchar(base) + 2L)
  write_utf8(paste0(hashes, "  ", labels),
             file.path(rd, paste0("checksums.", algo)))
}

# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

MODELS_RRBLUP <- c("gblup", "rrblup")
MODELS_BGLR <- c(bayesb = "BayesB", bayesc = "BayesC", brr = "BRR")
CV_SCHEMES <- c("kfold", "family", "lofo", "cv1", "cv2")

read_csv_fast <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  if (requireNamespace("data.table", quietly = TRUE))
    return(as.data.frame(data.table::fread(path, data.table = FALSE,
                                           check.names = FALSE)))
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

# Genotype file -> numeric matrix, lines x markers, row names = line IDs.
read_genotypes <- function(path, transposed = FALSE) {
  df <- read_csv_fast(path)
  if (ncol(df) < 3L) stop("Genotype file needs an ID column and at least two markers.", call. = FALSE)
  ids <- as.character(df[[1L]])
  X <- as.matrix(df[, -1L, drop = FALSE])
  if (!is.numeric(X)) {
    X <- suppressWarnings(apply(X, 2, as.numeric))
    if (all(is.na(X)))
      stop("Genotype values are not numeric. Recode them as 0/1/2 or -1/0/1 ",
           "(for example with claw-marker-qc).", call. = FALSE)
  }
  dimnames(X) <- list(ids, colnames(df)[-1L])
  if (transposed) X <- t(X)
  if (anyDuplicated(rownames(X)))
    stop("Duplicated line IDs in the genotype file: ",
         paste(utils::head(unique(rownames(X)[duplicated(rownames(X))]), 5), collapse = ", "),
         call. = FALSE)
  X
}

# Recode to rrBLUP's -1/0/1 (fractional values allowed).
code_markers <- function(X, coding = "auto") {
  v <- X[!is.na(X)]
  if (!length(v)) stop("The genotype matrix has no observed values.", call. = FALSE)
  if (coding == "auto") {
    coding <- if (any(v < 0)) "-101" else if (any(v > 1)) "012" else
      stop("All genotype values lie between 0 and 1, so the coding is ambiguous. ",
           "Give --geno-coding 01 (two homozygous classes coded 0/1), 012 or -101.",
           call. = FALSE)
  }
  rng <- switch(coding, "012" = c(0, 2), "-101" = c(-1, 1), "01" = c(0, 1),
                stop("--geno-coding must be auto, 012, -101 or 01.", call. = FALSE))
  if (any(v < rng[1] - 1e-8 | v > rng[2] + 1e-8))
    stop(sprintf("Genotype values outside [%g, %g] for coding %s: check --geno-coding.",
                 rng[1], rng[2], coding), call. = FALSE)
  M <- switch(coding, "012" = X - 1, "-101" = X, "01" = 2 * X - 1)
  list(M = M, coding = coding)
}

# A.mat (rrBLUP): mean imputation, monomorphic markers removed; W = centred
# imputed markers so that A = W W' / c and RR-BLUP on W equals GBLUP on A.
build_relationship <- function(M) {
  n_na <- sum(is.na(M))
  am <- rrBLUP::A.mat(M, impute.method = "mean", return.imputed = TRUE)
  imp <- am$imputed
  W <- sweep(imp, 2, colMeans(imp))
  p <- (colMeans(imp) + 1) / 2                     # frequency of the allele coded 1
  list(A = am$A, W = W, c = 2 * sum(p * (1 - p)), n_imputed = n_na, n_markers_in = ncol(M),
       n_markers_used = ncol(imp), missing_rate = n_na / length(M),
       marker_missing = colMeans(is.na(M)))
}

read_phenotypes <- function(path, id_col, trait, family_col = "", env_col = "") {
  df <- read_csv_fast(path)
  need <- c(id_col, if (nzchar(family_col)) family_col, if (nzchar(env_col)) env_col)
  miss <- setdiff(need, names(df))
  if (length(miss))
    stop("Phenotype file is missing column(s): ", paste(miss, collapse = ", "),
         ". Available: ", paste(names(df), collapse = ", "), call. = FALSE)
  if (!nzchar(trait)) {
    cand <- setdiff(names(df)[vapply(df, is.numeric, TRUE)], need)
    if (length(cand) != 1L)
      stop("Give --trait; numeric columns: ", paste(cand, collapse = ", "), call. = FALSE)
    trait <- cand
  }
  if (!trait %in% names(df))
    stop("Trait column '", trait, "' not found. Available: ",
         paste(names(df), collapse = ", "), call. = FALSE)
  if (!is.numeric(df[[trait]])) stop("Trait '", trait, "' must be numeric.", call. = FALSE)
  out <- data.frame(id = as.character(df[[id_col]]), y = df[[trait]],
                    family = if (nzchar(family_col)) as.character(df[[family_col]]) else NA_character_,
                    env = if (nzchar(env_col)) as.character(df[[env_col]]) else NA_character_,
                    stringsAsFactors = FALSE)
  n_na <- sum(is.na(out$y))
  out <- out[!is.na(out$y), , drop = FALSE]
  key <- if (nzchar(env_col)) paste(out$id, out$env) else out$id
  if (anyDuplicated(key))
    stop(if (nzchar(env_col)) "Duplicated line x environment records. " else
           "Duplicated line IDs. Give one BLUE per line, or --env-col for multi-environment data. ",
         "Examples: ", paste(utils::head(unique(key[duplicated(key)]), 5), collapse = ", "),
         call. = FALSE)
  list(data = out, trait = trait, n_missing = n_na)
}

# ---------------------------------------------------------------------------
# Cross-validation folds
# ---------------------------------------------------------------------------

# Random balanced k folds; returns an integer fold label per element.
assign_kfold <- function(n, k) sample(rep_len(seq_len(k), n))

# Folds balanced within each family (each fold gets about 1/k of every family).
assign_family_stratified <- function(family, k) {
  f <- integer(length(family))
  for (fa in unique(family)) {
    idx <- which(family == fa)
    idx <- idx[sample.int(length(idx))]
    f[idx] <- rep_len(sample.int(k), length(idx))
  }
  f
}

# CV2 (Burgueno et al. 2012): records of a line are spread over different folds,
# so each test record has records of the same line in other environments in
# the training set.
assign_cv2 <- function(line, k) {
  f <- integer(length(line))
  for (l in unique(line)) {
    idx <- which(line == l)
    idx <- idx[sample.int(length(idx))]
    f[idx] <- rep_len(sample.int(k), length(idx))
  }
  f
}

# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------

# Train on y (named by line ID) and predict lines in pred_ids. Returns a named
# vector of predicted phenotypes (intercept + genetic value).
fit_predict <- function(model, y, pred_ids, rel, bglr = list(), seed = 1L) {
  ids <- names(y)
  if (model == "gblup") {
    kb <- rrBLUP::kin.blup(data.frame(id = ids, y = unname(y)), geno = "id",
                           pheno = "y", K = rel$A)
    p <- kb$pred[pred_ids]
    attr(p, "vc") <- c(Vg = kb$Vg, Ve = kb$Ve)
    return(p)
  }
  if (model == "rrblup") {
    ms <- rrBLUP::mixed.solve(unname(y), Z = rel$W[ids, , drop = FALSE])
    p <- ms$beta[1] + drop(rel$W[pred_ids, , drop = FALSE] %*% ms$u)
    names(p) <- pred_ids
    attr(p, "vc") <- c(Vu = ms$Vu, Ve = ms$Ve)
    return(p)
  }
  if (model %in% names(MODELS_BGLR)) {
    all_ids <- union(ids, pred_ids)
    yy <- stats::setNames(rep(NA_real_, length(all_ids)), all_ids)
    yy[ids] <- y
    pre <- file.path(bglr$dir, sprintf("%s_%d_", model, seed))
    set.seed(seed)
    fm <- suppressMessages(BGLR::BGLR(
      y = unname(yy), ETA = list(list(X = rel$W[all_ids, , drop = FALSE],
                                      model = MODELS_BGLR[[model]])),
      nIter = bglr$n_iter, burnIn = bglr$burn_in, thin = bglr$thin,
      saveAt = pre, verbose = FALSE))
    p <- stats::setNames(fm$yHat, all_ids)[pred_ids]
    attr(p, "prefix") <- pre
    return(p)
  }
  stop("Unknown model '", model, "'. Use gblup, rrblup, bayesb, bayesc or brr.", call. = FALSE)
}

# Multi-environment main-effects GBLUP: y_ij = env_j + g_i + e_ij, g ~ N(0, A Vg),
# fitted with rrBLUP::mixed.solve (environments fixed, lines random).
fit_me <- function(rec, rel) {
  lev <- sort(unique(rec$env))
  X <- stats::model.matrix(~ factor(env, levels = lev), rec)
  lines <- rownames(rel$A)
  Z <- matrix(0, nrow(rec), length(lines))
  Z[cbind(seq_len(nrow(rec)), match(rec$id, lines))] <- 1
  ms <- rrBLUP::mixed.solve(rec$y, Z = Z, K = rel$A, X = X)
  env_means <- drop(stats::model.matrix(~ factor(env, levels = lev),
                                        data.frame(env = lev)) %*% ms$beta)
  list(u = stats::setNames(ms$u, lines), env_means = stats::setNames(env_means, lev),
       Vg = ms$Vu, Ve = ms$Ve)
}

fit_predict_me <- function(rec_train, rec_pred, rel) {
  fit <- fit_me(rec_train, rel)
  if (!all(rec_pred$env %in% names(fit$env_means)))
    stop("A test environment has no training records.", call. = FALSE)
  unname(fit$env_means[rec_pred$env] + fit$u[rec_pred$id])
}

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

# Predictive ability r(y, yhat), bias = slope of y on yhat (1 = unbiased),
# RMSE. Never computed on training records.
fold_metrics <- function(y, yhat) {
  n <- length(y)
  if (n < 3L || stats::sd(yhat) == 0 || stats::sd(y) == 0)
    return(c(n_test = n, r = NA_real_, slope = NA_real_,
             rmse = sqrt(mean((y - yhat)^2))))
  c(n_test = n, r = stats::cor(y, yhat), slope = stats::cov(y, yhat) / stats::var(yhat),
    rmse = sqrt(mean((y - yhat)^2)))
}

# For each test line, the largest genomic relationship (scaled as a correlation)
# with any training line: a direct measure of how related test and training are.
max_relatedness <- function(Kn, test_ids, train_ids) {
  s <- Kn[test_ids, train_ids, drop = FALSE]
  stats::median(apply(s, 1, max))
}

# ---------------------------------------------------------------------------
# Cross-validation drivers
# ---------------------------------------------------------------------------

run_cv_single <- function(ph, rel, models, schemes, k, reps, seed, bglr) {
  Kn <- stats::cov2cor(rel$A)
  folds_out <- list()
  preds_out <- list()
  for (s in schemes) {
    nrep <- if (s == "lofo") 1L else reps
    for (r in seq_len(nrep)) {
      set.seed(seed + 1000L * match(s, CV_SCHEMES) + r)
      lab <- switch(s,
                    kfold = assign_kfold(nrow(ph), k),
                    family = assign_family_stratified(ph$family, k),
                    lofo = as.integer(factor(ph$family)))
      fold_names <- if (s == "lofo") levels(factor(ph$family)) else as.character(seq_len(k))
      for (f in sort(unique(lab))) {
        te <- ph$id[lab == f]
        tr <- ph$id[lab != f]
        y_tr <- stats::setNames(ph$y[lab != f], tr)
        rel_max <- max_relatedness(Kn, te, tr)
        for (m in models) {
          fs <- seed + 100000L * match(m, c(MODELS_RRBLUP, names(MODELS_BGLR))) +
            1000L * match(s, CV_SCHEMES) + 10L * r + f
          p <- fit_predict(m, y_tr, te, rel, bglr, seed = fs)
          y_te <- ph$y[lab == f]
          met <- fold_metrics(y_te, unname(p))
          folds_out[[length(folds_out) + 1L]] <- data.frame(
            scheme = s, model = m, rep = r, fold = fold_names[f],
            n_train = length(tr), n_test = met[["n_test"]], r = met[["r"]],
            slope = met[["slope"]], rmse = met[["rmse"]],
            median_max_relatedness = rel_max, stringsAsFactors = FALSE)
          preds_out[[length(preds_out) + 1L]] <- data.frame(
            scheme = s, model = m, rep = r, fold = fold_names[f], id = te,
            family = ph$family[lab == f], observed = y_te, predicted = unname(p),
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  list(folds = do.call(rbind, folds_out), preds = do.call(rbind, preds_out))
}

run_cv_multienv <- function(ph, rel, schemes, k, reps, seed) {
  folds_out <- list()
  preds_out <- list()
  for (s in schemes) {
    for (r in seq_len(reps)) {
      set.seed(seed + 1000L * match(s, CV_SCHEMES) + r)
      if (s == "cv1") {
        lines <- unique(ph$id)
        lf <- stats::setNames(assign_kfold(length(lines), k), lines)
        lab <- unname(lf[ph$id])
      } else {
        lab <- assign_cv2(ph$id, k)
      }
      for (f in sort(unique(lab))) {
        tr <- ph[lab != f, , drop = FALSE]
        te <- ph[lab == f, , drop = FALSE]
        p <- fit_predict_me(tr, te, rel)
        for (e in sort(unique(te$env))) {
          sel <- te$env == e
          met <- fold_metrics(te$y[sel], p[sel])
          folds_out[[length(folds_out) + 1L]] <- data.frame(
            scheme = s, model = "gblup", rep = r, fold = as.character(f), env = e,
            n_train = nrow(tr), n_test = met[["n_test"]], r = met[["r"]],
            slope = met[["slope"]], rmse = met[["rmse"]],
            test_lines_seen_in_training = mean(te$id[sel] %in% tr$id),
            stringsAsFactors = FALSE)
        }
        preds_out[[length(preds_out) + 1L]] <- data.frame(
          scheme = s, model = "gblup", rep = r, fold = as.character(f), id = te$id,
          env = te$env, observed = te$y, predicted = p, stringsAsFactors = FALSE)
      }
    }
  }
  list(folds = do.call(rbind, folds_out), preds = do.call(rbind, preds_out))
}

summarise_cv <- function(folds, h2 = 0) {
  key <- unique(folds[, c("scheme", "model")])
  rows <- lapply(seq_len(nrow(key)), function(j) {
    d <- folds[folds$scheme == key$scheme[j] & folds$model == key$model[j], ]
    mr <- mean(d$r, na.rm = TRUE)
    data.frame(scheme = key$scheme[j], model = key$model[j], n_folds = sum(!is.na(d$r)),
               mean_n_train = mean(d$n_train), mean_n_test = mean(d$n_test),
               predictive_ability_mean = mr,
               predictive_ability_sd = stats::sd(d$r, na.rm = TRUE),
               accuracy = if (h2 > 0) mr / sqrt(h2) else NA_real_,
               slope_mean = mean(d$slope, na.rm = TRUE),
               slope_sd = stats::sd(d$slope, na.rm = TRUE),
               median_max_relatedness = if ("median_max_relatedness" %in% names(d))
                 stats::median(d$median_max_relatedness) else NA_real_,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Convergence summary of the residual variance chain written by BGLR.
bglr_diagnostics <- function(prefix, burn_in, thin) {
  f <- paste0(prefix, "varE.dat")
  if (!file.exists(f)) return(NULL)
  v <- scan(f, quiet = TRUE)
  post <- v[-seq_len(min(length(v) - 2L, floor(burn_in / thin)))]
  out <- list(trace = v, post = post, geweke_z = NA_real_, ess = NA_real_)
  if (requireNamespace("coda", quietly = TRUE) && length(post) > 20L) {
    ch <- coda::mcmc(post)
    out$geweke_z <- unname(coda::geweke.diag(ch)$z)
    out$ess <- unname(coda::effectiveSize(ch))
  }
  out
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

fmt <- function(x, d = 3) ifelse(is.na(x), "", sprintf(paste0("%.", d, "f"), x))

md_table <- function(df) {
  hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  body <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  c(hdr, sep, body)
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  dd <- file.path(SKILL_DIR, "demo")
  if (opt$demo) {
    opt$geno <- file.path(dd, "genotypes.csv")
    opt$pheno <- file.path(dd, "phenotypes.csv")
    opt$trait <- "yield_kg_ha"
    opt$family_col <- "family"
    if (!"--cv" %in% argv) opt$cv <- "kfold,family,lofo"
    if (!"--h2" %in% argv) opt$h2 <- 0.5
  }
  if (opt$demo_multienv) {
    opt$geno <- file.path(dd, "genotypes.csv")
    opt$pheno <- file.path(dd, "phenotypes_multienv.csv")
    opt$trait <- "yield_kg_ha"
    opt$family_col <- "family"
    opt$env_col <- "env"
    if (!"--cv" %in% argv) opt$cv <- "cv1,cv2"
    if (!"--models" %in% argv) opt$models <- "gblup"
  }
  if (!nzchar(opt$geno) || !nzchar(opt$pheno))
    stop("Provide --geno and --pheno, or --demo. See --help.", call. = FALSE)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  models <- unique(trimws(strsplit(tolower(opt$models), ",", fixed = TRUE)[[1]]))
  schemes <- unique(trimws(strsplit(tolower(opt$cv), ",", fixed = TRUE)[[1]]))
  bad <- setdiff(models, c(MODELS_RRBLUP, names(MODELS_BGLR)))
  if (length(bad)) stop("Unknown model(s): ", paste(bad, collapse = ", "),
                        ". Use gblup, rrblup, bayesb, bayesc, brr.", call. = FALSE)
  bad <- setdiff(schemes, CV_SCHEMES)
  if (length(bad)) stop("Unknown CV scheme(s): ", paste(bad, collapse = ", "),
                        ". Use kfold, family, lofo, cv1, cv2.", call. = FALSE)
  multienv <- nzchar(opt$env_col)
  if (multienv && any(!schemes %in% c("cv1", "cv2")))
    stop("With --env-col use --cv cv1 and/or cv2 (multi-environment schemes).", call. = FALSE)
  if (!multienv && any(schemes %in% c("cv1", "cv2")))
    stop("--cv cv1/cv2 need multi-environment data: give --env-col.", call. = FALSE)
  if (any(schemes %in% c("family", "lofo")) && !nzchar(opt$family_col))
    stop("--cv family and lofo need --family-col.", call. = FALSE)
  if (opt$folds < 2 || opt$reps < 1) stop("--folds must be >= 2 and --reps >= 1.", call. = FALSE)
  if (opt$h2 < 0 || opt$h2 > 1) stop("--h2 must be between 0 and 1.", call. = FALSE)
  if (opt$burn_in >= opt$n_iter) stop("--burn-in must be smaller than --n-iter.", call. = FALSE)
  k <- as.integer(opt$folds)
  reps <- as.integer(opt$reps)
  seed <- as.integer(opt$seed)
  warn <- character()

  # --- data ---------------------------------------------------------------
  X <- read_genotypes(opt$geno, opt$geno_transposed)
  cm <- code_markers(X, opt$geno_coding)
  phr <- read_phenotypes(opt$pheno, opt$id_col, opt$trait, opt$family_col, opt$env_col)
  ph <- phr$data
  if (!any(ph$id %in% rownames(X))) {
    hint <- if (any(ph$id %in% colnames(X))) " The genotype file looks transposed: add --geno-transposed." else ""
    stop("No phenotyped line ID matches a genotype ID.", hint, call. = FALSE)
  }
  no_geno <- unique(ph$id[!ph$id %in% rownames(X)])
  if (length(no_geno))
    warn <- c(warn, sprintf("%d phenotyped line(s) have no genotypes and were excluded (e.g. %s).",
                            length(no_geno), paste(utils::head(no_geno, 3), collapse = ", ")))
  ph <- ph[ph$id %in% rownames(X), , drop = FALSE]
  if (phr$n_missing > 0)
    warn <- c(warn, sprintf("%d record(s) with a missing '%s' were excluded.", phr$n_missing, phr$trait))
  rel <- build_relationship(cm$M)
  if (rel$missing_rate > 0)
    warn <- c(warn, sprintf(paste0("%.2f%% of genotype calls were missing and were imputed with ",
                                   "the marker mean (rrBLUP::A.mat)."), 100 * rel$missing_rate))
  if (rel$missing_rate > 0.10 || any(rel$marker_missing > 0.5))
    warn <- c(warn, "Genotype missingness is high (>10% overall or markers >50% missing): run claw-marker-qc first.")
  if (rel$n_markers_used < rel$n_markers_in)
    warn <- c(warn, sprintf("%d monomorphic or fully missing marker(s) removed by A.mat.",
                            rel$n_markers_in - rel$n_markers_used))
  if (grepl("blup", opt$trait, ignore.case = TRUE) || grepl("blup", basename(opt$pheno), ignore.case = TRUE))
    warn <- c(warn, paste0("The phenotype looks like BLUPs. Training on BLUPs shrinks twice and biases ",
                           "the slope; use BLUEs (fixed genotype effects) unless the user confirms otherwise."))
  n_lines <- length(unique(ph$id))
  if (n_lines < 100)
    warn <- c(warn, sprintf("Only %d phenotyped lines: cross-validation estimates will be very noisy.", n_lines))
  z <- (ph$y - mean(ph$y)) / stats::sd(ph$y)
  if (!multienv && any(abs(z) > 4))
    warn <- c(warn, sprintf("%d phenotype value(s) are more than 4 SD from the mean; check for data errors.",
                            sum(abs(z) > 4)))
  if (any(schemes %in% c("family", "lofo"))) {
    if (anyNA(ph$family) || any(ph$family == ""))
      stop("Some phenotyped lines have no family label.", call. = FALSE)
    fam_n <- table(ph$family[!duplicated(ph$id)])
    if ("lofo" %in% schemes && length(fam_n) < 3)
      stop("Leave-one-family-out needs at least 3 families.", call. = FALSE)
    if (any(fam_n < 10))
      warn <- c(warn, sprintf("Families with fewer than 10 lines (%s): their leave-one-family-out r is unreliable.",
                              paste(names(fam_n)[fam_n < 10], collapse = ", ")))
  }
  Kn <- stats::cov2cor(rel$A)
  dup <- which(Kn > 0.95 & upper.tri(Kn), arr.ind = TRUE)
  if (nrow(dup))
    warn <- c(warn, sprintf(paste0("%d pair(s) of lines are near-identical (genomic correlation > 0.95), ",
                                   "e.g. %s. Duplicates split across folds inflate cross-validation."),
                            nrow(dup), paste(rownames(Kn)[dup[1, 1]], rownames(Kn)[dup[1, 2]], sep = " & ")))
  if (multienv && any(models != "gblup")) {
    warn <- c(warn, "Multi-environment CV uses the main-effects GBLUP model only; other models were skipped.")
    models <- "gblup"
  }

  # --- cross-validation --------------------------------------------------------
  bglr <- list(n_iter = opt$n_iter, burn_in = opt$burn_in, thin = opt$thin,
               dir = file.path(tempdir(), paste0("bglr_", Sys.getpid())))
  dir.create(bglr$dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(bglr$dir, recursive = TRUE), add = TRUE)
  t0 <- Sys.time()
  cv <- if (multienv) run_cv_multienv(ph, rel, schemes, k, reps, seed) else
    run_cv_single(ph, rel, models, schemes, k, reps, seed, bglr)
  summ <- summarise_cv(cv$folds, opt$h2)
  if (multienv) {
    env_summ <- stats::aggregate(r ~ scheme + env, cv$folds, function(x) mean(x, na.rm = TRUE))
    names(env_summ)[3] <- "predictive_ability_mean"
    save_table(env_summ, opt$out, "cv_by_environment.csv")
  }
  save_table(cv$folds, opt$out, "cv_folds.csv")
  save_table(summ, opt$out, "cv_summary.csv")
  save_table(cv$preds, opt$out, "cv_predictions.csv")
  if (opt$h2 == 0)
    warn <- c(warn, paste0("Accuracy (r / sqrt(h2)) was not computed: give --h2, the entry-mean heritability ",
                           "of the phenotypes used (e.g. from claw-trial-analysis)."))
  if (opt$h2 > 0 && any(summ$accuracy > 1, na.rm = TRUE))
    warn <- c(warn, "Estimated accuracy exceeds 1: --h2 is probably too low or not on the entry-mean basis.")
  if (any(abs(summ$slope_mean - 1) > 0.2, na.rm = TRUE))
    warn <- c(warn, paste0("Mean slope of observed on predicted differs from 1 by more than 0.2 for: ",
                           paste(unique(paste(summ$scheme, summ$model)[abs(summ$slope_mean - 1) > 0.2]),
                                 collapse = ", "), " (biased or over/under-dispersed predictions)."))
  if (all(c("kfold", "lofo") %in% schemes)) {
    gap <- mean(summ$predictive_ability_mean[summ$scheme == "kfold"]) -
      mean(summ$predictive_ability_mean[summ$scheme == "lofo"])
    if (is.finite(gap) && gap > 0.1)
      warn <- c(warn, sprintf(paste0("Random k-fold predictive ability is %.2f higher than leave-one-family-out: ",
                                     "random folds keep relatives in training. Report the scheme that matches ",
                                     "the real use (new families vs. new lines from known families)."), gap))
  }

  # --- final model on all phenotyped lines: predictions for every genotyped line
  gebv <- data.frame(id = rownames(rel$A), phenotyped = rownames(rel$A) %in% ph$id,
                     stringsAsFactors = FALSE)
  diag_rows <- list()
  vc <- NULL
  if (!multienv) {
    y_all <- stats::setNames(ph$y, ph$id)
    gebv$observed <- unname(y_all[gebv$id])
    for (m in models) {
      p <- fit_predict(m, y_all, gebv$id, rel, bglr, seed = seed)
      gebv[[paste0("predicted_", m)]] <- unname(p)
      if (m == "gblup") vc <- attr(p, "vc")
      if (m %in% names(MODELS_BGLR)) {
        dg <- bglr_diagnostics(attr(p, "prefix"), opt$burn_in, opt$thin)
        if (!is.null(dg)) {
          diag_rows[[m]] <- data.frame(model = m, n_iter = opt$n_iter, burn_in = opt$burn_in,
                                       thin = opt$thin, varE_mean = mean(dg$post),
                                       geweke_z = dg$geweke_z, ess = dg$ess)
          open_png(opt$out, sprintf("bglr_trace_%s.png", m))
          graphics::plot(seq_along(dg$trace) * opt$thin, dg$trace, type = "l",
                         xlab = "Iteration", ylab = "Residual variance",
                         main = sprintf("%s: residual variance trace", MODELS_BGLR[[m]]))
          graphics::abline(v = opt$burn_in, lty = 2, col = "red")
          grDevices::dev.off()
          if ((!is.na(dg$geweke_z) && abs(dg$geweke_z) > 2) || (!is.na(dg$ess) && dg$ess < 100))
            warn <- c(warn, sprintf(paste0("%s chain may not have converged (Geweke z = %.2f, effective ",
                                           "sample size = %.0f); increase --n-iter and --burn-in."),
                                    m, dg$geweke_z, dg$ess))
        }
      }
    }
  } else {
    fit <- fit_me(ph, rel)
    gebv$genetic_value_gblup <- unname(fit$u[gebv$id])
    for (e in names(fit$env_means))
      gebv[[paste0("predicted_", e)]] <- unname(fit$env_means[[e]] + fit$u[gebv$id])
    vc <- c(Vg = fit$Vg, Ve = fit$Ve)
  }
  save_table(gebv, opt$out, "predictions_all_lines.csv")
  if (length(diag_rows)) save_table(do.call(rbind, diag_rows), opt$out, "bglr_convergence.csv")
  vc_tab <- NULL
  if (!is.null(vc)) {
    pid <- unique(ph$id)
    Ap <- rel$A[pid, pid]
    dk <- mean(diag(Ap)) - mean(Ap)
    vc_tab <- data.frame(Vg = vc[["Vg"]], Ve = vc[["Ve"]], scaling_D = dk,
                         genomic_h2 = vc[["Vg"]] * dk / (vc[["Vg"]] * dk + vc[["Ve"]]))
    save_table(vc_tab, opt$out, "variance_components.csv")
  }
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # --- figures --------------------------------------------------------------
  open_png(opt$out, "cv_predictive_ability.png", width = 8, height = 5)
  grp <- interaction(cv$folds$model, cv$folds$scheme, sep = "\n", lex.order = TRUE, drop = TRUE)
  graphics::par(mar = c(6, 5, 3, 1))
  graphics::boxplot(cv$folds$r ~ grp, las = 1, col = "lightsteelblue", xlab = "",
                    ylab = "Predictive ability r(y, y_hat) per fold",
                    main = "Cross-validated predictive ability", cex.axis = 0.8)
  graphics::abline(h = 0, lty = 3)
  grDevices::dev.off()
  first <- cv$preds[cv$preds$scheme == schemes[1] & cv$preds$model == models[1] &
                      cv$preds$rep == 1L, ]
  open_png(opt$out, "cv_observed_vs_predicted.png", width = 6, height = 6)
  graphics::plot(first$predicted, first$observed, pch = 16, col = grDevices::adjustcolor("steelblue", 0.6),
                 xlab = sprintf("Cross-validated prediction (%s, %s, rep 1)", models[1], schemes[1]),
                 ylab = sprintf("Observed %s", phr$trait), main = "Out-of-fold predictions")
  graphics::abline(stats::lm(observed ~ predicted, first), col = "red", lwd = 2)
  graphics::abline(0, 1, lty = 2)
  graphics::legend("topleft", c("regression of y on prediction", "1:1 line"),
                   col = c("red", "black"), lty = c(1, 2), lwd = c(2, 1), bty = "n")
  grDevices::dev.off()

  # --- report ---------------------------------------------------------------
  agree <- NULL
  if (all(c("gblup", "rrblup") %in% models)) {
    a <- cv$preds[cv$preds$model == "gblup", ]
    b <- cv$preds[cv$preds$model == "rrblup", ]
    agree <- stats::cor(a$predicted, b$predicted)
  }
  scheme_text <- c(
    kfold = "random k-fold: lines assigned to folds at random; relatives of test lines stay in training",
    family = "family-stratified k-fold: each fold holds about 1/k of every family (prediction of new lines from known families)",
    lofo = "leave-one-family-out: each family is predicted from the others (prediction of new families)",
    cv1 = "CV1: whole lines held out (new lines never tested in any environment)",
    cv2 = "CV2: line x environment records held out; the same line is observed in other environments (sparse testing)")
  body <- c("## Data", "",
            sprintf("Genotypes: %d lines, %d markers used of %d (coding %s, recoded to -1/0/1).  ",
                    nrow(rel$A), rel$n_markers_used, rel$n_markers_in, cm$coding),
            sprintf("Phenotypes: trait `%s`, %d records on %d genotyped lines%s.  ",
                    phr$trait, nrow(ph), n_lines,
                    if (multienv) sprintf(" in %d environments", length(unique(ph$env))) else ""),
            sprintf("Selection candidates without phenotypes: %d.", sum(!gebv$phenotyped)), "",
            "## Cross-validation", "",
            paste0("- ", scheme_text[schemes]),
            sprintf("- %d folds x %d repetitions for k-fold schemes; seed %d. Training-set fit is never reported as accuracy.",
                    k, reps, seed), "",
            md_table(data.frame(scheme = summ$scheme, model = summ$model, folds = summ$n_folds,
                                n_train = fmt(summ$mean_n_train, 0),
                                r_mean = fmt(summ$predictive_ability_mean),
                                r_sd = fmt(summ$predictive_ability_sd),
                                accuracy = fmt(summ$accuracy),
                                slope = fmt(summ$slope_mean),
                                relatedness = fmt(summ$median_max_relatedness))),
            "",
            paste0("`r_mean` and `r_sd`: predictive ability r(y, y_hat) across folds. `accuracy` = r_mean / sqrt(h2), ",
                   if (opt$h2 > 0) sprintf("with h2 = %.3f supplied by --h2 (entry-mean heritability of the phenotypes)", opt$h2)
                   else "not computed because --h2 was not given",
                   ". `slope`: regression of observed on predicted (1 = unbiased; < 1 = predictions too dispersed). ",
                   "`relatedness`: median over folds of the median largest genomic correlation between a test line and any training line."),
            "")
  if (multienv)
    body <- c(body, "Correlations are computed within each environment and then averaged (pooling environments would inflate r with environment means); see `tables/cv_by_environment.csv`.", "")
  if (!is.null(agree))
    body <- c(body, sprintf("GBLUP and RR-BLUP out-of-fold predictions correlate at %.6f (they are equivalent models).", agree), "")
  if (!is.null(vc_tab))
    body <- c(body, sprintf(paste0("Full-data GBLUP (REML): Vg = %.4g, Ve = %.4g, genomic heritability ",
                                   "= %.2f on the basis of %s (Legarra 2016 scaling; a sanity check, not the h2 used for accuracy)."),
                            vc_tab$Vg, vc_tab$Ve, vc_tab$genomic_h2,
                            if (multienv) "single line x environment records" else "the phenotypes supplied"), "")
  if (length(diag_rows))
    body <- c(body, "BGLR convergence (residual variance chain after burn-in):", "",
              md_table(data.frame(model = names(diag_rows),
                                  geweke_z = fmt(vapply(diag_rows, function(d) d$geweke_z, 0)),
                                  ess = fmt(vapply(diag_rows, function(d) d$ess, 0), 0))), "")
  body <- c(body, "![Predictive ability](figures/cv_predictive_ability.png)", "",
            "![Observed vs predicted](figures/cv_observed_vs_predicted.png)", "",
            "## Predictions for all lines", "",
            "`tables/predictions_all_lines.csv` has predictions from models trained on all phenotyped lines, for every genotyped line. For unphenotyped candidates these are genomic predictions; for phenotyped lines they are fitted values and say nothing about accuracy.",
            sprintf("Run time: %.0f s.", secs), "")
  limits <- c("## Limitations", "",
              "- Cross-validated accuracy applies to lines like the test folds: random folds within a breeding cycle overestimate accuracy for new families, new cycles or new environments.",
              "- Accuracy = r / sqrt(h2) assumes h2 is the entry-mean heritability of the phenotypes used, estimated without bias; errors in h2 propagate directly.",
              "- The multi-environment model has no G x E term; predictions differ between environments only by the environment mean.",
              "")
  refs <- c("## References", "",
            "- Endelman JB (2011) Plant Genome 4:250-255. https://doi.org/10.3835/plantgenome2011.08.0024",
            "- VanRaden PM (2008) J Dairy Sci 91:4414-4423. https://doi.org/10.3168/jds.2007-0980",
            "- Meuwissen THE, Hayes BJ, Goddard ME (2001) Genetics 157:1819-1829. https://doi.org/10.1093/genetics/157.4.1819",
            "- Perez P, de los Campos G (2014) Genetics 198:483-495. https://doi.org/10.1534/genetics.114.164442",
            "- Burgueno J et al. (2012) Crop Sci 52:707-719. https://doi.org/10.2135/cropsci2011.06.0299",
            "- Daetwyler HD, Villanueva B, Woolliams JA (2008) PLoS ONE 3:e3395. https://doi.org/10.1371/journal.pone.0003395",
            "- Habier D, Fernando RL, Dekkers JCM (2007) Genetics 177:2389-2397. https://doi.org/10.1534/genetics.107.081190",
            "- Legarra A (2016) Theor Popul Biol 107:26-30. https://doi.org/10.1016/j.tpb.2015.08.005")
  report <- c(sprintf("# %s report", SKILL), "",
              if (length(warn)) c("## Warnings", "", paste("-", warn), "") else NULL,
              body, limits, refs, "", sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(report, file.path(opt$out, "report.md"))
  write_repro_bundle(opt$out, argv, inputs = c(opt$geno, opt$pheno))
  cat(sprintf("Done: %d folds evaluated. Report: %s\n", nrow(cv$folds),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(list(cv = cv, summary = summ, predictions = gebv, warnings = warn))
}

if (sys.nframe() == 0L) main()
