#!/usr/bin/env Rscript
# claw-yield-prediction: predict yield from remote-sensing features (VI
# time-series metrics, per-date VIs, weather summaries) with elastic net
# (glmnet) and random forest (ranger), and estimate accuracy with nested,
# leakage-safe cross-validation: random k-fold (for comparison only),
# leave-environment-out, leave-genotype-out, leave-year-out and
# leave-genotype-and-environment-out.
#
# Usage:
#   Rscript scripts/yield_prediction.R --demo
#   Rscript scripts/yield_prediction.R --input plots.csv --target yield_kg_ha --env-col env --genotype-col genotype
#   Rscript scripts/yield_prediction.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-yield-prediction"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/yield_prediction.R [--demo] [--input FILE.csv] [options]

Input
  --input          CSV, one row per plot (or per genotype x environment mean/BLUE)
  --target         yield column, units in its name (default: yield_kg_ha)
  --features       comma-separated feature columns (default: the features in
                   --feature-time, else every numeric column that is not the
                   target or an id/grouping column)
  --id-col         plot identifier column (default: plot_id; optional)
  --env-col        environment column (default: env)
  --genotype-col   genotype column (default: genotype)
  --year-col       year column (default: year)
  --feature-time   CSV with columns feature, available_t: the time (e.g. DAP)
                   at which each feature is known, such as
                   tables/feature_time.csv from claw-vi-timeseries
  --cutoff         prediction time, same units as available_t: features known
                   later are excluded (prevents temporal leakage)
  --drop-na-rows   drop rows with missing features (counted in the report)
Validation
  --cv             comma-separated schemes: random, leave-env-out,
                   leave-genotype-out, leave-year-out, leave-genotype-env-out
                   (default: random plus every grouped scheme the columns allow)
  --primary        scheme to headline (default: the strictest one run)
  --folds          folds for random and leave-genotype-out (default 5)
  --inner-folds    inner folds for tuning (default 5)
  --seed           random seed (default 1)
Models
  --models         enet, rf, or enet,rf (default)
  --alpha-grid     elastic-net mixing values (default 0,0.25,0.5,0.75,1)
  --enet-lambda    1se (default, glmnet's default) or min
  --num-trees      trees per random forest (default 500)
  --mtry-grid      comma-separated mtry values (default: sqrt(p), p/3, 2p/3)
  --min-node-grid  comma-separated min.node.size values (default 5)
  --perm-repeats   permutations per feature for importance (default 5)
Output
  --out            output folder (default: agriclaw_output/claw-yield-prediction)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, input = "", target = "yield_kg_ha", features = "",
                 id_col = "plot_id", env_col = "env", genotype_col = "genotype",
                 year_col = "year", feature_time = "", cutoff = NA_real_,
                 drop_na_rows = FALSE, cv = "", primary = "", folds = 5,
                 inner_folds = 5, seed = 1, models = "enet,rf",
                 alpha_grid = "0,0.25,0.5,0.75,1", enet_lambda = "1se",
                 num_trees = 500, mtry_grid = "", min_node_grid = "5",
                 perm_repeats = 5, out = "")

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
# Analysis
# ---------------------------------------------------------------------------

# Strictest first: the default headline is the first scheme run from this list.
SCHEMES <- c("leave-genotype-env-out", "leave-year-out", "leave-env-out",
             "leave-genotype-out", "random")
MODELS <- c("enet", "rf")
MIN_ROWS        <- 10   # hard minimum to attempt cross-validation
SMALL_N         <- 50   # warn below this many rows (rule of thumb)
FEW_GROUPS      <- 4    # warn when a grouped scheme has fewer test groups (rule of thumb)
MAX_P_RATIO     <- 0.2  # warn when features > 0.2 x rows (rule of thumb)
TARGET_COPY_R   <- 0.95 # |r| between a feature and the target that suggests a copy (rule of thumb)
LEAKAGE_GAP_R   <- 0.1  # report a random-vs-grouped gap in r above this (rule of thumb)
UNIT_TOKENS <- c("kg_ha", "t_ha", "mg_ha", "bu_ac", "bu_a", "g_m2", "g_plot", "kg_plot",
                 "lb_ac", "g_plant", "kg", "g", "t", "bu", "lb", "pct", "percent")
POSITION_NAMES <- c("row", "col", "column", "range", "rep", "replicate", "block", "plot",
                    "x", "y", "lat", "lon", "latitude", "longitude", "easting", "northing",
                    "pass", "tier", "bay")
PAL <- c(enet = "#2a78d6", rf = "#eb6834", train_mean = "#8a8984")

split_list <- function(x) {
  x <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  x[nzchar(x)]
}

num_list <- function(x, what) {
  v <- suppressWarnings(as.numeric(split_list(x)))
  if (!length(v) || anyNA(v)) stop(what, " needs comma-separated numbers.", call. = FALSE)
  v
}

read_csv_checked <- function(path, what = "Input") {
  if (!nzchar(path) || !file.exists(path))
    stop(what, " file not found: '", path, "'.", call. = FALSE)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

has_col <- function(dat, col) nzchar(col) && col %in% names(dat)

# --- data preparation -------------------------------------------------------

# Validate the table, choose the features and drop the ones not yet known at
# the cutoff. Returns the analysis data plus a feature log and notes.
prepare_data <- function(dat, opt, feature_time = NULL) {
  notes <- character()
  if (!opt$target %in% names(dat))
    stop("Target column '", opt$target, "' not found. Available: ",
         paste(names(dat), collapse = ", "), ". Name it with --target.", call. = FALSE)
  if (!is.numeric(dat[[opt$target]]))
    stop("Target column '", opt$target, "' must be numeric.", call. = FALSE)
  n_na_y <- sum(is.na(dat[[opt$target]]))
  if (n_na_y) {
    notes <- c(notes, sprintf("%d row(s) with missing '%s' were excluded.", n_na_y, opt$target))
    dat <- dat[!is.na(dat[[opt$target]]), , drop = FALSE]
  }
  if (has_col(dat, opt$id_col) && anyDuplicated(dat[[opt$id_col]]))
    stop("Column '", opt$id_col, "' has repeated identifiers, e.g. '",
         dat[[opt$id_col]][duplicated(dat[[opt$id_col]])][1],
         "'. Give one row per plot (or per genotype x environment).", call. = FALSE)
  group_cols <- c(opt$env_col, opt$genotype_col, opt$year_col)
  group_cols <- group_cols[vapply(group_cols, function(k) has_col(dat, k), TRUE)]
  for (k in group_cols)
    if (anyNA(dat[[k]]) || any(!nzchar(as.character(dat[[k]]))))
      stop("Grouping column '", k, "' has missing values: fill them or drop those rows.", call. = FALSE)
  reserved <- c(opt$target, if (has_col(dat, opt$id_col)) opt$id_col, group_cols)

  # candidate features
  if (nzchar(opt$features)) {
    feats <- split_list(opt$features)
    miss <- setdiff(feats, names(dat))
    if (length(miss)) stop("Feature column(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
    bad <- feats[!vapply(feats, function(f) is.numeric(dat[[f]]), TRUE)]
    if (length(bad)) stop("Feature column(s) must be numeric: ", paste(bad, collapse = ", "), call. = FALSE)
    if (length(intersect(feats, reserved)))
      stop("The target and id/grouping columns cannot be features: ",
           paste(intersect(feats, reserved), collapse = ", "), call. = FALSE)
    source_note <- "--features"
  } else if (!is.null(feature_time)) {
    feats <- intersect(feature_time$feature, names(dat))
    feats <- feats[vapply(feats, function(f) is.numeric(dat[[f]]), TRUE)]
    feats <- setdiff(feats, reserved)
    other <- setdiff(names(dat)[vapply(dat, is.numeric, TRUE)], c(feats, reserved))
    if (length(other))
      notes <- c(notes, sprintf("Numeric column(s) not listed in --feature-time were not used as features: %s. Add them to that table (with their available time) or name features with --features.",
                                paste(other, collapse = ", ")))
    source_note <- "features listed in --feature-time"
  } else {
    feats <- names(dat)[vapply(dat, is.numeric, TRUE)]
    feats <- setdiff(feats, reserved)
    source_note <- "every numeric column except the target and id/grouping columns"
  }
  log <- data.frame(feature = feats, used = TRUE, reason = "", available_t = NA_real_,
                    stringsAsFactors = FALSE)

  # temporal leakage: drop features not yet known at the cutoff
  if (!is.na(opt$cutoff)) {
    if (is.null(feature_time))
      stop("--cutoff needs --feature-time (a CSV with feature, available_t) so that features known only after the cutoff can be excluded.", call. = FALSE)
    log$available_t <- feature_time$available_t[match(log$feature, feature_time$feature)]
    unk <- log$feature[is.na(log$available_t)]
    if (length(unk))
      stop("No available_t for feature(s): ", paste(unk, collapse = ", "),
           ". Add them to the --feature-time table or leave them out with --features.", call. = FALSE)
    late <- log$available_t > opt$cutoff
    log$used[late] <- FALSE
    log$reason[late] <- sprintf("known only at %g > cutoff %g", log$available_t[late], opt$cutoff)
  } else if (!is.null(feature_time)) {
    log$available_t <- feature_time$available_t[match(log$feature, feature_time$feature)]
  }

  # missing values
  used <- log$feature[log$used]
  na_rows <- if (length(used)) !stats::complete.cases(dat[, used, drop = FALSE]) else rep(FALSE, nrow(dat))
  if (any(na_rows)) {
    if (!opt$drop_na_rows) {
      cnt <- colSums(is.na(dat[, used, drop = FALSE]))
      cnt <- cnt[cnt > 0]
      stop(sprintf("%d row(s) have missing feature values (%s). Drop those features with --features, or the rows with --drop-na-rows; nothing is imputed.",
                   sum(na_rows), paste(sprintf("%s: %d", names(cnt), cnt), collapse = ", ")), call. = FALSE)
    }
    notes <- c(notes, sprintf("%d row(s) with missing feature values were dropped (--drop-na-rows).", sum(na_rows)))
    dat <- dat[!na_rows, , drop = FALSE]
  }
  # constant features carry no information
  const <- vapply(used, function(f) stats::sd(dat[[f]]) == 0, TRUE)
  if (any(const)) {
    log$used[log$feature %in% used[const]] <- FALSE
    log$reason[log$feature %in% used[const]] <- "constant"
  }
  used <- log$feature[log$used]
  if (!length(used)) stop("No usable features remain. See the feature list and --cutoff.", call. = FALSE)
  if (nrow(dat) < MIN_ROWS)
    stop(sprintf("Only %d rows: cross-validated accuracy needs at least %d (and many more to be useful).",
                 nrow(dat), MIN_ROWS), call. = FALSE)
  list(dat = dat, features = used, feature_log = log, notes = notes, source = source_note,
       group_cols = group_cols)
}

# Leakage and design warnings that do not stop the analysis.
data_warnings <- function(prep, opt) {
  dat <- prep$dat
  f <- prep$features
  w <- character()
  y <- dat[[opt$target]]
  low <- tolower(opt$target)
  if (!any(vapply(UNIT_TOKENS, function(u) grepl(paste0("(^|_)", u, "$"), low), TRUE)))
    w <- c(w, sprintf("Target '%s' has no units in its name (e.g. yield_kg_ha). Confirm the units and whether these are plot values or BLUEs before reporting RMSE.", opt$target))
  pos <- f[tolower(f) %in% POSITION_NAMES | grepl("^(row|col|range|rep|block)(_|$)", tolower(f))]
  if (length(pos))
    w <- c(w, sprintf("Feature(s) look like field-position or design columns: %s. Position predicts yield through spatial trends and makes the model a field map, not a sensor model. Remove them unless that is intended.", paste(pos, collapse = ", ")))
  rr <- vapply(f, function(k) abs(stats::cor(dat[[k]], y)), 0)
  cp <- f[rr > TARGET_COPY_R]
  if (length(cp))
    w <- c(w, sprintf("Feature(s) correlate with the target at |r| > %g: %s. Check that they are not yield itself, harvest data or derived from the target (target leakage).", TARGET_COPY_R, paste(cp, collapse = ", ")))
  if (has_col(dat, opt$env_col)) {
    env <- as.character(dat[[opt$env_col]])
    ne <- length(unique(env))
    envlev <- f[vapply(f, function(k) all(tapply(dat[[k]], env, function(v) length(unique(v)) == 1)), TRUE)]
    if (length(envlev) && ne > 1)
      w <- c(w, sprintf("%d feature(s) are constant within each environment (e.g. weather summaries: %s). With %d environments they act as environment labels: random CV rewards them, leave-env-out tests them honestly.",
                        length(envlev), paste(utils::head(envlev, 5), collapse = ", "), ne))
  }
  if (nrow(dat) < SMALL_N)
    w <- c(w, sprintf("Only %d rows: cross-validated accuracy is very uncertain at this size. Treat differences between models as noise unless they are large.", nrow(dat)))
  if (length(f) > MAX_P_RATIO * nrow(dat))
    w <- c(w, sprintf("%d features for %d rows: high risk of overfitting; nested CV keeps the estimate honest but it will be noisy.", length(f), nrow(dat)))
  w
}

# --- folds -----------------------------------------------------------------------

scheme_group_col <- function(scheme, opt) {
  switch(scheme, random = "", "leave-env-out" = opt$env_col,
         "leave-genotype-out" = opt$genotype_col, "leave-year-out" = opt$year_col,
         "leave-genotype-env-out" = opt$env_col)
}

available_schemes <- function(dat, opt) {
  lv <- function(k) if (has_col(dat, k)) length(unique(dat[[k]])) else 0L
  ok <- c("leave-genotype-env-out" = lv(opt$env_col) >= 2 && lv(opt$genotype_col) >= 2,
          "leave-year-out" = lv(opt$year_col) >= 2,
          "leave-env-out" = lv(opt$env_col) >= 2,
          "leave-genotype-out" = lv(opt$genotype_col) >= 2,
          random = TRUE)
  names(ok)[ok]
}

scheme_seed <- function(seed, scheme) as.integer(seed + 1000L * match(scheme, rev(SCHEMES)))

# Outer folds: list of list(train, test, label). Every row is tested exactly once.
make_folds <- function(dat, scheme, opt) {
  n <- nrow(dat)
  set.seed(scheme_seed(opt$seed, scheme))
  if (scheme == "random") {
    k <- min(opt$folds, n)
    f <- sample(rep_len(seq_len(k), n))
    return(lapply(seq_len(k), function(j) list(test = which(f == j), train = which(f != j),
                                               label = paste0("fold", j))))
  }
  if (scheme %in% c("leave-env-out", "leave-year-out")) {
    g <- as.character(dat[[scheme_group_col(scheme, opt)]])
    return(lapply(sort(unique(g)), function(v) list(test = which(g == v), train = which(g != v),
                                                    label = v)))
  }
  geno <- as.character(dat[[opt$genotype_col]])
  gl <- sort(unique(geno))
  k <- min(opt$folds, length(gl))
  gfold <- stats::setNames(sample(rep_len(seq_len(k), length(gl))), gl)
  gf <- unname(gfold[geno])
  if (scheme == "leave-genotype-out")
    return(lapply(seq_len(k), function(j) {
      te <- gf == j
      list(test = which(te), train = which(!te), label = paste0("genotype-fold", j))
    }))
  # leave-genotype-env-out: test = one environment x one genotype fold; train =
  # other environments AND other genotype folds (new lines in a new environment)
  env <- as.character(dat[[opt$env_col]])
  out <- list()
  for (e in sort(unique(env))) for (j in seq_len(k)) {
    te <- env == e & gf == j
    if (!any(te)) next
    tr <- env != e & gf != j
    out[[length(out) + 1L]] <- list(test = which(te), train = which(tr),
                                    label = sprintf("%s|genotype-fold%d", e, j))
  }
  out
}

# Inner fold ids for tuning, grouped like the outer scheme when the training
# set has at least 3 groups; otherwise random (and say so).
inner_foldid <- function(n, groups, k, seed) {
  set.seed(seed)
  if (!is.null(groups)) {
    gl <- unique(groups)
    if (length(gl) >= 3) {
      kk <- min(k, length(gl))
      gf <- stats::setNames(sample(rep_len(seq_len(kk), length(gl))), gl)
      return(list(id = unname(gf[groups]), type = sprintf("grouped (%d folds)", kk)))
    }
  }
  kk <- max(3L, min(k, n))
  list(id = sample(rep_len(seq_len(kk), n)),
       type = if (is.null(groups)) sprintf("random (%d folds)", kk) else
         sprintf("random (%d folds; fewer than 3 groups in training)", kk))
}

# --- models ----------------------------------------------------------------------

pad_x <- function(x) if (ncol(x) == 1L) cbind(x, .pad = 0) else x  # glmnet needs >= 2 columns

# Elastic net (Zou & Hastie 2005) via glmnet (Friedman et al. 2010). Features
# are standardised inside every fit from that fit's training rows only
# (glmnet's standardize = TRUE); alpha and lambda are tuned by inner CV.
fit_enet <- function(x, y, foldid, alphas, rule = "1se") {
  x <- pad_x(x)
  best <- NULL
  for (a in alphas) {
    cv <- glmnet::cv.glmnet(x, y, alpha = a, foldid = foldid, family = "gaussian",
                            standardize = TRUE, type.measure = "mse")
    err <- min(cv$cvm)
    if (is.null(best) || err < best$err) best <- list(err = err, alpha = a, cv = cv)
  }
  lam <- if (rule == "min") best$cv$lambda.min else best$cv$lambda.1se
  fit <- best$cv$glmnet.fit
  list(predict = function(xn) as.numeric(stats::predict(fit, newx = pad_x(xn), s = lam)),
       params = c(alpha = best$alpha, lambda = lam, inner_cv_mse = best$err),
       coef = as.numeric(stats::coef(fit, s = lam)))
}

# Random forest (Breiman 2001) via ranger (Wright & Ziegler 2017); mtry and
# min.node.size tuned by inner CV. Trees need no scaling.
fit_rf <- function(x, y, foldid, mtry_grid, node_grid, num_trees, seed) {
  xd <- as.data.frame(x)
  grid <- expand.grid(mtry = mtry_grid, node = node_grid)
  cv_mse <- vapply(seq_len(nrow(grid)), function(g) {
    se <- 0
    for (j in sort(unique(foldid))) {
      tr <- foldid != j
      m <- ranger::ranger(x = xd[tr, , drop = FALSE], y = y[tr], num.trees = num_trees,
                          mtry = grid$mtry[g], min.node.size = grid$node[g], seed = seed,
                          verbose = FALSE)
      se <- se + sum((stats::predict(m, xd[!tr, , drop = FALSE], num.threads = 1)$predictions - y[!tr])^2)
    }
    se / length(y)
  }, 0)
  b <- which.min(cv_mse)
  m <- ranger::ranger(x = xd, y = y, num.trees = num_trees, mtry = grid$mtry[b],
                      min.node.size = grid$node[b], seed = seed, verbose = FALSE)
  list(predict = function(xn) stats::predict(m, as.data.frame(xn), num.threads = 1)$predictions,
       params = c(mtry = grid$mtry[b], min_node_size = grid$node[b], inner_cv_mse = cv_mse[b]))
}

default_mtry <- function(p) sort(unique(pmin(p, pmax(1L, round(c(sqrt(p), p / 3, 2 * p / 3))))))

# --- metrics ---------------------------------------------------------------------

# r, RMSE, relative RMSE, bias and both R2 definitions (Kvalseth 1985):
# r2_cor = squared correlation (ignores bias and scale); r2_sse = 1 - SSE/SST
# (can be negative: worse than predicting the mean of the observations).
pred_metrics <- function(obs, pred) {
  n <- length(obs)
  r <- if (n > 2 && stats::sd(pred) > 0 && stats::sd(obs) > 0) stats::cor(obs, pred) else NA_real_
  rmse <- sqrt(mean((pred - obs)^2))
  sst <- sum((obs - mean(obs))^2)
  c(n = n, r = r, rmse = rmse,
    rrmse_pct = if (mean(obs) > 0) 100 * rmse / mean(obs) else NA_real_,
    bias = mean(pred - obs), r2_cor = r^2,
    r2_sse = if (sst > 0) 1 - sum((obs - pred)^2) / sst else NA_real_)
}

within_group_r <- function(obs, pred, group) {
  vapply(split(seq_along(obs), group), function(i) {
    if (length(i) < 3 || stats::sd(pred[i]) == 0 || stats::sd(obs[i]) == 0) return(NA_real_)
    stats::cor(obs[i], pred[i])
  }, 0)
}

rmse <- function(obs, pred) sqrt(mean((pred - obs)^2))

# Permutation importance on held-out rows (Breiman 2001): the increase in test
# RMSE when one feature is shuffled within the test fold, averaged over
# `repeats` shuffles. All shuffled copies are predicted in one call.
perm_importance <- function(predict_fun, x, y, repeats, seed) {
  n <- nrow(x)
  p <- ncol(x)
  base <- rmse(y, predict_fun(x))
  set.seed(seed)
  big <- x[rep(seq_len(n), p * repeats), , drop = FALSE]
  for (j in seq_len(p)) for (b in seq_len(repeats)) {
    rows <- ((j - 1L) * repeats + b - 1L) * n + seq_len(n)
    big[rows, j] <- x[sample.int(n), j]
  }
  pr <- predict_fun(big)
  vapply(seq_len(p), function(j) {
    mean(vapply(seq_len(repeats), function(b) {
      rmse(y, pr[((j - 1L) * repeats + b - 1L) * n + seq_len(n)])
    }, 0)) - base
  }, 0)
}

# --- cross-validation ------------------------------------------------------------

run_scheme <- function(dat, features, scheme, opt, keep_models = FALSE, importance = TRUE) {
  x <- as.matrix(dat[, features, drop = FALSE])
  y <- dat[[opt$target]]
  folds <- make_folds(dat, scheme, opt)
  gcol <- scheme_group_col(scheme, opt)
  models <- split_list(opt$models)
  alphas <- num_list(opt$alpha_grid, "--alpha-grid")
  mtry <- if (nzchar(opt$mtry_grid)) num_list(opt$mtry_grid, "--mtry-grid") else default_mtry(length(features))
  mtry <- sort(unique(pmin(length(features), pmax(1, round(mtry)))))
  nodes <- num_list(opt$min_node_grid, "--min-node-grid")
  base_seed <- scheme_seed(opt$seed, scheme)
  pred <- list()
  tuning <- list()
  imp <- list()
  kept <- list()
  skipped <- character()
  for (fi in seq_along(folds)) {
    fo <- folds[[fi]]
    if (length(fo$train) < 5) {
      skipped <- c(skipped, sprintf("%s (%d training rows)", fo$label, length(fo$train)))
      next
    }
    xtr <- x[fo$train, , drop = FALSE]
    ytr <- y[fo$train]
    xte <- x[fo$test, , drop = FALSE]
    yte <- y[fo$test]
    groups <- if (nzchar(gcol)) as.character(dat[[gcol]][fo$train]) else NULL
    inner <- inner_foldid(length(ytr), groups, opt$inner_folds, base_seed + fi)
    fold_pred <- list(train_mean = rep(mean(ytr), length(yte)))
    for (m in models) {
      fit <- if (m == "enet") fit_enet(xtr, ytr, inner$id, alphas, opt$enet_lambda) else
        fit_rf(xtr, ytr, inner$id, mtry, nodes, opt$num_trees, base_seed + fi)
      fold_pred[[m]] <- fit$predict(xte)
      tuning[[length(tuning) + 1L]] <- data.frame(
        scheme = scheme, model = m, fold = fo$label, n_train = length(ytr), n_test = length(yte),
        inner_cv = inner$type, parameter = names(fit$params), value = unname(fit$params),
        stringsAsFactors = FALSE)
      if (importance && opt$perm_repeats > 0) {
        im <- perm_importance(fit$predict, xte, yte, opt$perm_repeats, base_seed + 7L * fi)
        imp[[length(imp) + 1L]] <- data.frame(scheme = scheme, model = m, fold = fo$label,
                                              n_test = length(yte), feature = colnames(xte),
                                              rmse_increase = im, stringsAsFactors = FALSE)
      }
      if (keep_models) kept[[paste(m, fo$label)]] <- fit
    }
    for (m in names(fold_pred))
      pred[[length(pred) + 1L]] <- data.frame(row = fo$test, scheme = scheme, model = m,
                                              fold = fo$label, observed = yte,
                                              predicted = fold_pred[[m]], stringsAsFactors = FALSE)
  }
  list(pred = do.call(rbind, pred), tuning = do.call(rbind, tuning), imp = do.call(rbind, imp),
       n_folds = length(folds), models = kept, skipped = skipped)
}

summarise_cv <- function(pred, dat, opt) {
  env <- if (has_col(dat, opt$env_col)) as.character(dat[[opt$env_col]]) else NULL
  keys <- unique(pred[, c("scheme", "model")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    s <- pred[pred$scheme == keys$scheme[i] & pred$model == keys$model[i], , drop = FALSE]
    m <- pred_metrics(s$observed, s$predicted)
    wr <- if (!is.null(env)) within_group_r(s$observed, s$predicted, env[s$row]) else NA_real_
    data.frame(scheme = keys$scheme[i], model = keys$model[i], n = m[["n"]],
               n_folds = length(unique(s$fold)), r = m[["r"]], rmse = m[["rmse"]],
               rrmse_pct = m[["rrmse_pct"]], bias = m[["bias"]], r2_cor = m[["r2_cor"]],
               r2_sse = m[["r2_sse"]],
               r_within_env = if (all(is.na(wr))) NA_real_ else mean(wr, na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$flag <- ifelse(is.finite(out$r2_sse) & out$r2_sse < 0 & out$model != "train_mean",
                     "worse than predicting the mean", "")
  out
}

fold_metrics <- function(pred) {
  keys <- unique(pred[, c("scheme", "model", "fold")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    s <- pred[pred$scheme == keys$scheme[i] & pred$model == keys$model[i] &
                pred$fold == keys$fold[i], , drop = FALSE]
    m <- pred_metrics(s$observed, s$predicted)
    data.frame(keys[i, ], n_test = m[["n"]], r = m[["r"]], rmse = m[["rmse"]],
               bias = m[["bias"]], r2_sse = m[["r2_sse"]], stringsAsFactors = FALSE)
  }))
}

summarise_importance <- function(imp) {
  if (is.null(imp)) return(NULL)
  keys <- unique(imp[, c("scheme", "model", "feature")])
  out <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    s <- imp[imp$scheme == keys$scheme[i] & imp$model == keys$model[i] &
               imp$feature == keys$feature[i], , drop = FALSE]
    data.frame(keys[i, ], rmse_increase = stats::weighted.mean(s$rmse_increase, s$n_test),
               sd_over_folds = if (nrow(s) > 1) stats::sd(s$rmse_increase) else NA_real_,
               stringsAsFactors = FALSE)
  }))
  out <- out[order(out$scheme, out$model, -out$rmse_increase), , drop = FALSE]
  out$rank <- stats::ave(out$rmse_increase, out$scheme, out$model, FUN = function(v) seq_along(v))
  rownames(out) <- NULL
  out
}

run_cv <- function(prep, opt, keep_models = FALSE) {
  dat <- prep$dat
  avail <- available_schemes(dat, opt)
  schemes <- if (nzchar(opt$cv)) split_list(opt$cv) else avail
  bad <- setdiff(schemes, SCHEMES)
  if (length(bad))
    stop("Unknown --cv scheme(s): ", paste(bad, collapse = ", "), ". Choose from: ",
         paste(rev(SCHEMES), collapse = ", "), call. = FALSE)
  na <- setdiff(schemes, avail)
  if (length(na))
    stop("Scheme(s) not possible with these columns: ", paste(na, collapse = ", "),
         ". Each grouped scheme needs its grouping column (--env-col, --genotype-col, --year-col) with at least 2 levels.",
         call. = FALSE)
  schemes <- SCHEMES[SCHEMES %in% schemes]
  models <- split_list(opt$models)
  if (!length(models) || length(setdiff(models, MODELS)))
    stop("--models must be enet, rf or enet,rf.", call. = FALSE)
  if (!opt$enet_lambda %in% c("1se", "min")) stop("--enet-lambda must be 1se or min.", call. = FALSE)
  primary <- if (nzchar(opt$primary)) opt$primary else schemes[1]
  if (!primary %in% schemes) stop("--primary must be one of the schemes run: ", paste(schemes, collapse = ", "), call. = FALSE)
  res <- lapply(schemes, function(s) run_scheme(dat, prep$features, s, opt, keep_models,
                                                importance = s == primary))
  names(res) <- schemes
  pred <- do.call(rbind, lapply(res, `[[`, "pred"))
  rownames(pred) <- NULL
  list(schemes = schemes, primary = primary, pred = pred,
       metrics = summarise_cv(pred, dat, opt), folds = fold_metrics(pred),
       tuning = do.call(rbind, lapply(res, `[[`, "tuning")),
       importance = summarise_importance(do.call(rbind, lapply(res, `[[`, "imp"))),
       models = lapply(res, `[[`, "models"),
       skipped = unlist(lapply(names(res), function(s) if (length(res[[s]]$skipped))
         paste0(s, ": ", res[[s]]$skipped))))
}

# --- figures -----------------------------------------------------------------------

plot_accuracy <- function(cv, out_dir) {
  met <- cv$metrics
  sch <- rev(cv$schemes)
  models <- c(setdiff(unique(met$model), "train_mean"), "train_mean")
  panels <- c(r = "Pearson r (pooled)", r2_sse = "R2 = 1 - SSE/SST",
              r_within_env = "Mean r within environments")
  panels <- panels[vapply(names(panels), function(k) any(is.finite(met[[k]])), TRUE)]
  open_png(out_dir, "accuracy_by_scheme.png", width = 3.4 * length(panels) + 1.6, height = 0.5 * length(sch) + 2.4)
  op <- graphics::par(mfrow = c(1, length(panels)), mar = c(4, 1, 1, 0.8), oma = c(0, 11, 2.2, 0))
  for (k in names(panels)) {
    v <- met[[k]][met$model %in% models]
    xl <- range(c(0, v), na.rm = TRUE)
    xl <- xl + c(-0.05, 0.05) * diff(xl)
    graphics::plot(NA, xlim = xl, ylim = c(0.5, length(sch) + 0.5), yaxt = "n", xlab = panels[[k]],
                   ylab = "", main = "", bty = "n")
    graphics::abline(h = seq_along(sch), col = "grey92")
    graphics::abline(v = 0, col = "grey55", lty = 2)
    if (k == names(panels)[1]) graphics::axis(2, at = seq_along(sch), labels = sch, las = 1, tick = FALSE, outer = TRUE)
    off <- seq(-0.18, 0.18, length.out = length(models))
    for (mi in seq_along(models)) {
      s <- met[met$model == models[mi], , drop = FALSE]
      yy <- match(s$scheme, sch) + off[mi]
      graphics::points(s[[k]], yy, pch = if (models[mi] == "train_mean") 1 else 19,
                       col = PAL[[models[mi]]], cex = 1.4)
    }
  }
  graphics::par(op)
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  graphics::plot.new()
  graphics::legend("top", legend = models, pch = ifelse(models == "train_mean", 1, 19),
                   col = PAL[models], bty = "n", horiz = TRUE, cex = 0.9, pt.cex = 1.3)
  grDevices::dev.off()
}

plot_obs_pred <- function(cv, dat, opt, out_dir) {
  s <- cv$pred[cv$pred$scheme == cv$primary & cv$pred$model != "train_mean", , drop = FALSE]
  models <- unique(s$model)
  open_png(out_dir, "obs_vs_pred.png", width = 4.6 * length(models), height = 4.8)
  op <- graphics::par(mfrow = c(1, length(models)), mar = c(4.2, 4.2, 3, 1))
  lim <- range(c(s$observed, s$predicted))
  env <- if (has_col(dat, opt$env_col)) as.character(dat[[opt$env_col]]) else NULL
  for (m in models) {
    d <- s[s$model == m, , drop = FALSE]
    mt <- cv$metrics[cv$metrics$scheme == cv$primary & cv$metrics$model == m, ]
    graphics::plot(d$predicted, d$observed, xlim = lim, ylim = lim, pch = 19, cex = 0.7,
                   col = grDevices::adjustcolor(PAL[[m]], 0.45),
                   xlab = paste("predicted", opt$target), ylab = paste("observed", opt$target),
                   main = sprintf("%s, %s\nr = %.2f, R2 (1-SSE/SST) = %.2f", m, cv$primary, mt$r, mt$r2_sse),
                   cex.main = 0.9)
    graphics::abline(0, 1, col = "grey40", lty = 2)
    if (!is.null(env)) {  # environment centroids, labelled directly
      e <- env[d$row]
      cx <- tapply(d$predicted, e, mean)
      cy <- tapply(d$observed, e, mean)
      graphics::points(cx, cy, pch = 23, bg = "white", col = "grey15", cex = 1.6)
      graphics::text(cx, cy, names(cx), pos = 4, cex = 0.75, col = "grey15")
    }
  }
  graphics::par(op)
  grDevices::dev.off()
}

plot_importance <- function(cv, out_dir, top = 15) {
  imp <- cv$importance
  if (is.null(imp)) return(invisible(NULL))
  imp <- imp[imp$scheme == cv$primary, , drop = FALSE]
  models <- unique(imp$model)
  nshow <- min(top, max(table(imp$model)))
  open_png(out_dir, "importance.png", width = 4.8 * length(models), height = 0.32 * nshow + 1.8)
  op <- graphics::par(mfrow = c(1, length(models)), mar = c(4.2, 10, 3, 1))
  for (m in models) {
    s <- utils::head(imp[imp$model == m, , drop = FALSE], top)
    s <- s[rev(seq_len(nrow(s))), , drop = FALSE]
    graphics::barplot(s$rmse_increase, names.arg = s$feature, horiz = TRUE, las = 1,
                      col = PAL[[m]], border = NA, cex.names = 0.8,
                      xlab = "held-out RMSE increase",
                      main = sprintf("%s, %s", m, cv$primary), cex.main = 0.95)
    graphics::abline(v = 0, col = "grey40")
  }
  graphics::par(op)
  grDevices::dev.off()
}

# --- report --------------------------------------------------------------------------

fmt <- function(x, d = 3) ifelse(is.finite(x), as.character(signif(x, d)), "NA")

cv_warnings <- function(cv, prep, opt) {
  met <- cv$metrics
  w <- character()
  neg <- met[met$flag != "", , drop = FALSE]
  if (nrow(neg))
    w <- c(w, sprintf("Out-of-sample R2 (1 - SSE/SST) is negative for %s: those models predict worse than the mean of the observed yields. Do not use their predictions as phenotypes (e.g. in GWAS, selection or as a covariate); they carry no reliable yield information.",
                      paste(sprintf("%s under %s (%s)", neg$model, neg$scheme, fmt(neg$r2_sse)), collapse = "; ")))
  if ("random" %in% cv$schemes && cv$primary != "random") {
    for (m in setdiff(unique(met$model), "train_mean")) {
      a <- met[met$scheme == "random" & met$model == m, ]
      b <- met[met$scheme == cv$primary & met$model == m, ]
      if (nrow(a) && nrow(b) && is.finite(a$r) && (!is.finite(b$r) || a$r - b$r > LEAKAGE_GAP_R))
        w <- c(w, sprintf("Random k-fold overstates accuracy for %s: r = %s (R2 %s) vs %s under %s. Rows of the same genotype or environment in training and test sets leak information; report the grouped estimate.",
                          m, fmt(a$r), fmt(a$r2_sse), fmt(b$r), cv$primary))
    }
  }
  if (has_col(prep$dat, opt$env_col)) {
    pm <- met[met$scheme == cv$primary & met$model != "train_mean", , drop = FALSE]
    gap <- pm[is.finite(pm$r) & is.finite(pm$r_within_env) & pm$r - pm$r_within_env > 0.3, , drop = FALSE]
    if (nrow(gap))
      w <- c(w, sprintf("Pooled r is much higher than the mean within-environment r (%s): accuracy comes mostly from environment means, not from ranking genotypes within an environment.",
                        paste(sprintf("%s %s vs %s", gap$model, fmt(gap$r), fmt(gap$r_within_env)), collapse = "; ")))
    rk <- pm[is.finite(pm$r_within_env) & pm$r_within_env >= 0.3 & is.finite(pm$r2_sse) & pm$r2_sse < 0, , drop = FALSE]
    if (nrow(rk))
      w <- c(w, sprintf("Under %s the model still ranks genotypes within environments (mean within-environment r: %s) but misses the environment means (R2 < 0). Use its predictions only to rank entries within an environment, or remove environment means from yield and features before modelling.",
                        cv$primary, paste(sprintf("%s %s", rk$model, fmt(rk$r_within_env)), collapse = "; ")))
  }
  for (s in setdiff(cv$schemes, "random")) {
    nf <- length(unique(cv$pred$fold[cv$pred$scheme == s]))
    if (s %in% c("leave-env-out", "leave-year-out") && nf < FEW_GROUPS)
      w <- c(w, sprintf("%s rests on only %d held-out group(s): the estimate is anecdotal; more environments or years are needed for a stable number.", s, nf))
  }
  if (length(cv$skipped))
    w <- c(w, sprintf("Fold(s) skipped because fewer than 5 rows were left for training, so their rows have no prediction: %s.",
                      paste(cv$skipped, collapse = "; ")))
  if (!is.null(cv$tuning) && any(grepl("fewer than 3 groups", cv$tuning$inner_cv)))
    w <- c(w, "Some training sets had fewer than 3 groups, so their tuning folds were random rather than grouped (see tables/tuning.csv).")
  w
}

write_report <- function(cv, prep, opt, out_dir, warn) {
  met <- cv$metrics
  ord <- order(match(met$scheme, cv$schemes), match(met$model, c("enet", "rf", "train_mean")))
  met <- met[ord, , drop = FALSE]
  fl <- prep$feature_log
  excl <- fl[!fl$used, , drop = FALSE]
  imp <- cv$importance
  top <- if (!is.null(imp)) imp[imp$scheme == cv$primary & imp$rank <= 5, , drop = FALSE] else NULL
  lines <- c(sprintf("# %s report", SKILL), "",
             if (length(warn)) c("## Warnings", "", paste("-", warn), "") else NULL,
             "## Summary", "",
             sprintf("- Input: `%s`; target: `%s`; rows used: %d", opt$input, opt$target, nrow(prep$dat)),
             sprintf("- Features (%d, %s): %s", length(prep$features), prep$source, paste(prep$features, collapse = ", ")),
             if (nrow(excl)) sprintf("- Excluded features: %s", paste(sprintf("%s (%s)", excl$feature, excl$reason), collapse = "; ")),
             if (!is.na(opt$cutoff)) sprintf("- Prediction time (--cutoff): %g; only features known by then were used.", opt$cutoff),
             sprintf("- Models: %s; seed %g; outer schemes: %s; headline scheme: **%s**",
                     opt$models, opt$seed, paste(cv$schemes, collapse = ", "), cv$primary),
             sprintf("- Nested CV: tuning by inner %g-fold CV inside each outer training set, grouped like the outer scheme where possible; elastic-net lambda rule: %s; random forest: %g trees.",
                     opt$inner_folds, opt$enet_lambda, opt$num_trees),
             "", "## Cross-validated accuracy", "",
             "Every row is predicted once per scheme by a model that never saw it. `train_mean` predicts the training-fold mean: the baseline any useful model must beat.", "",
             "| scheme | model | n | r | RMSE | rRMSE % | bias | R2 (r^2) | R2 (1-SSE/SST) | mean r within env | flag |",
             "|---|---|---|---|---|---|---|---|---|---|---|",
             sprintf("| %s | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s |", met$scheme, met$model, met$n,
                     fmt(met$r), fmt(met$rmse, 4), fmt(met$rrmse_pct), fmt(met$bias), fmt(met$r2_cor),
                     fmt(met$r2_sse), fmt(met$r_within_env), met$flag),
             "",
             "- **r** and **R2 (r^2)** ignore bias and scale; **R2 (1-SSE/SST)** does not, and is negative when predictions are worse than the mean of the observations.",
             "- **Pooled r** across environments includes differences between environment means; **mean r within env** measures how well genotypes are ranked inside each environment.",
             "- Scheme meanings: random = rows shuffled (comparison only); leave-env-out = new environment, known genotypes; leave-year-out = new season; leave-genotype-out = new lines in known environments; leave-genotype-env-out = new lines in a new environment.", "",
             if (!is.null(top) && nrow(top)) c(sprintf("## Permutation importance (%s)", cv$primary), "",
                                               "| model | rank | feature | RMSE increase | SD over folds |", "|---|---|---|---|---|",
                                               sprintf("| %s | %d | %s | %s | %s |", top$model, as.integer(top$rank), top$feature,
                                                       fmt(top$rmse_increase), fmt(top$sd_over_folds)), "",
                                               "Correlated features share importance, and shuffling one of a correlated pair creates unrealistic rows; read importance as 'the model relies on', not 'yield is caused by'.", "") else NULL,
             "## Figures", "",
             "![Accuracy by scheme](figures/accuracy_by_scheme.png)", "",
             "![Observed vs predicted](figures/obs_vs_pred.png)", "",
             if (file.exists(file.path(out_dir, "figures", "importance.png"))) c("![Importance](figures/importance.png)", "") else NULL,
             "## Limitations", "",
             "- Accuracy applies to the scheme's scenario only: leave-env-out says nothing about new genotypes, leave-genotype-out nothing about new seasons.",
             "- Few environments or years make grouped estimates noisy; random forests cannot extrapolate beyond the training range of the features.",
             "- Predicted yields are shrunk towards the mean and inherit the features' genetics; they are not substitutes for measured yield in genetic analyses.", "",
             "## References", "",
             "- Zou H, Hastie T (2005) J R Stat Soc B 67:301-320. https://doi.org/10.1111/j.1467-9868.2005.00503.x",
             "- Friedman J, Hastie T, Tibshirani R (2010) J Stat Softw 33(1):1-22. https://doi.org/10.18637/jss.v033.i01",
             "- Breiman L (2001) Machine Learning 45:5-32. https://doi.org/10.1023/A:1010933404324",
             "- Wright MN, Ziegler A (2017) J Stat Softw 77(1):1-17. https://doi.org/10.18637/jss.v077.i01",
             "- Varma S, Simon R (2006) BMC Bioinformatics 7:91. https://doi.org/10.1186/1471-2105-7-91",
             "- Roberts DR et al. (2017) Ecography 40:913-929. https://doi.org/10.1111/ecog.02881",
             "- Kvalseth TO (1985) The American Statistician 39:279-285. https://doi.org/10.1080/00031305.1985.10479448", "",
             sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(lines, file.path(out_dir, "report.md"))
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (opt$demo) {
    opt$input <- file.path(SKILL_DIR, "demo", "rs_yield.csv")
    if (!nzchar(opt$feature_time)) opt$feature_time <- file.path(SKILL_DIR, "demo", "feature_time.csv")
  }
  if (!nzchar(opt$input)) stop("Provide --input FILE.csv, or --demo. See --help.", call. = FALSE)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  for (k in c("folds", "inner_folds", "num_trees"))
    if (opt[[k]] < 2 || opt[[k]] != round(opt[[k]]))
      stop("--", gsub("_", "-", k), " must be a whole number >= 2.", call. = FALSE)
  if (opt$perm_repeats < 0) stop("--perm-repeats must be >= 0.", call. = FALSE)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

  dat <- read_csv_checked(opt$input)
  ft <- NULL
  if (nzchar(opt$feature_time)) {
    ft <- read_csv_checked(opt$feature_time, "Feature-time")
    if (!all(c("feature", "available_t") %in% names(ft)) || !is.numeric(ft$available_t))
      stop("The --feature-time file needs columns 'feature' and numeric 'available_t'.", call. = FALSE)
  }
  prep <- prepare_data(dat, opt, ft)
  cv <- run_cv(prep, opt)
  warn <- c(prep$notes, data_warnings(prep, opt), cv_warnings(cv, prep, opt))

  ids <- data.frame(row = seq_len(nrow(prep$dat)))
  for (k in c(opt$id_col, opt$env_col, opt$genotype_col, opt$year_col))
    if (has_col(prep$dat, k)) ids[[k]] <- prep$dat[[k]]
  preds <- merge(ids, cv$pred, by = "row", sort = FALSE)
  preds <- preds[order(match(preds$scheme, cv$schemes), preds$model, preds$row), ]
  preds$row <- NULL
  save_table(cv$metrics, opt$out, "cv_metrics.csv")
  save_table(preds, opt$out, "cv_predictions.csv")
  save_table(cv$folds, opt$out, "fold_metrics.csv")
  save_table(cv$tuning, opt$out, "tuning.csv")
  if (!is.null(cv$importance)) save_table(cv$importance, opt$out, "permutation_importance.csv")
  save_table(prep$feature_log, opt$out, "features.csv")

  plot_accuracy(cv, opt$out)
  plot_obs_pred(cv, prep$dat, opt, opt$out)
  plot_importance(cv, opt$out)
  write_report(cv, prep, opt, opt$out, warn)
  write_repro_bundle(opt$out, argv, inputs = c(opt$input, opt$feature_time))
  hm <- cv$metrics[cv$metrics$scheme == cv$primary & cv$metrics$model != "train_mean", ]
  cat(sprintf("Done: %s; %s. %d warning(s). Report: %s\n", cv$primary,
              paste(sprintf("%s r = %.2f, R2 = %.2f", hm$model, hm$r, hm$r2_sse), collapse = "; "),
              length(warn), normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(cv)
}

if (sys.nframe() == 0L) main()
