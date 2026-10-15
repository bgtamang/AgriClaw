#!/usr/bin/env Rscript
# claw-trial-analysis: single-trial mixed-model analysis with statgenSTA.
# Fits a spatial P-spline model (SpATS) and the design-based model (lme4) with
# genotype fixed (BLUEs) and random (BLUPs, heritability), compares them, flags
# outliers on standardised residuals and writes per-genotype BLUEs with
# weights for multi-environment analysis (claw-met-gxe) or GWAS.
#
# Usage:
#   Rscript scripts/trial_analysis.R --demo
#   Rscript scripts/trial_analysis.R --input clean_data.csv --traits yield_kg_ha
#   Rscript scripts/trial_analysis.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-trial-analysis"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/trial_analysis.R --input DATA.csv --traits T1[,T2] [options]

  --input FILE       one row per plot (clean_data.csv from claw-fieldbook or similar)
  --traits LIST      numeric trait columns, comma-separated (units in the names)
  --genotype COL     genotype / entry column (default: entry, else genotype)
  --trial COL        trial column; each trial is analysed separately (default trial)
  --rep COL          replicate column (default rep)
  --block COL        incomplete-block column (default block)
  --row COL          field row coordinate (default row)
  --col COL          field column coordinate (default col)
  --check COL        TRUE/FALSE check column (default is_check)
  --design D         auto (default), rcbd, ibd, res.ibd, rowcol or res.rowcol
  --engine E         auto (default: SpATS and lme4 when coordinates exist), both,
                     SpATS or lme4
  --blues-from E     model whose BLUEs go to met_input.csv: auto (SpATS when fitted),
                     SpATS or lme4
  --outlier-limit Z  |standardised residual| limit (default: statgenSTA's
                     qnorm(1 - 0.5/residual df), kept between 2 and 4)
  --nseg C,R         P-spline segments for columns,rows (default: half the
                     number of columns and rows, statgenSTA's default)
  --demo             analyse demo/trial_clean.csv (yield_kg_ha, plant_height_cm)
  --out DIR          output folder (default: agriclaw_output/claw-trial-analysis)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, input = "", traits = "", genotype = "", trial = "trial",
                 rep = "rep", block = "block", row = "row", col = "col", check = "is_check",
                 design = "auto", engine = "auto", blues_from = "auto", outlier_limit = 0,
                 nseg = "", out = "")

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
# Input handling and design detection
# ---------------------------------------------------------------------------

fail <- function(...) stop(paste0(...), call. = FALSE)

STA_DESIGNS <- c("rcbd", "ibd", "res.ibd", "rowcol", "res.rowcol")
# statgenSTA renames columns to these names, so traits may not use them.
STA_RESERVED <- c("genotype", "trial", "loc", "year", "repId", "subBlock", "plotId", "rowCoord",
                  "colCoord", "rowId", "colId", "checkId")
# claw-trial-design design names -> statgenSTA model
DESIGN_MAP <- c(rcbd = "rcbd", alpha = "res.ibd", rowcol = "res.rowcol", augmented = "ibd",
                prep = "rowcol", diagonal = "rowcol", crd = "rowcol")

as_logical_col <- function(x) toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES", "Y")

read_trial_data <- function(path) {
  if (!file.exists(path)) fail("Input file not found: ", path)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
}

# Standardised working data for one data set; fails loudly on bad input.
prepare_data <- function(dat, opt) {
  gcol <- if (nzchar(opt$genotype)) opt$genotype else if ("entry" %in% names(dat)) "entry" else "genotype"
  if (!gcol %in% names(dat))
    fail("Genotype column '", gcol, "' not found. Use --genotype. Columns: ", paste(names(dat), collapse = ", "), ".")
  traits <- trimws(strsplit(opt$traits, ",", fixed = TRUE)[[1]])
  traits <- traits[nzchar(traits)]
  if (!length(traits)) fail("Give --traits (numeric trait columns, comma-separated).")
  miss <- setdiff(traits, names(dat))
  if (length(miss)) fail("Trait column(s) not found: ", paste(miss, collapse = ", "), ". Columns: ", paste(names(dat), collapse = ", "), ".")
  bad <- intersect(traits, STA_RESERVED)
  if (length(bad)) fail("Trait name(s) ", paste(bad, collapse = ", "), " clash with statgenSTA column names; rename them.")
  for (t in traits) if (!is.numeric(dat[[t]]))
    fail("Trait '", t, "' is not numeric (", paste(utils::head(unique(dat[[t]][is.na(suppressWarnings(as.numeric(dat[[t]])))]), 3), collapse = ", "),
         "). Clean it first, e.g. with claw-fieldbook.")
  has <- function(cl) nzchar(cl) && cl %in% names(dat) && any(!is.na(dat[[cl]]) & trimws(as.character(dat[[cl]])) != "")
  w <- data.frame(gen = trimws(as.character(dat[[gcol]])), stringsAsFactors = FALSE)
  if (any(is.na(w$gen) | w$gen == "")) fail("Blank genotype names in '", gcol, "'.")
  w$trialx <- if (has(opt$trial)) as.character(dat[[opt$trial]]) else "trial1"
  w$repx <- if (has(opt$rep)) as.character(dat[[opt$rep]]) else NA_character_
  w$blockx <- if (has(opt$block)) as.character(dat[[opt$block]]) else NA_character_
  w$rowx <- if (has(opt$row)) suppressWarnings(as.numeric(dat[[opt$row]])) else NA_real_
  w$colx <- if (has(opt$col)) suppressWarnings(as.numeric(dat[[opt$col]])) else NA_real_
  w$is_check <- if (has(opt$check)) as_logical_col(dat[[opt$check]]) else FALSE
  w$design_hint <- if ("design" %in% names(dat)) as.character(dat$design) else NA_character_
  w$plot_id <- if ("plot_id" %in% names(dat)) as.character(dat$plot_id) else as.character(seq_len(nrow(dat)))
  for (t in traits) w[[t]] <- dat[[t]]
  dropped <- 0L
  if ("is_filler" %in% names(dat)) {
    fil <- as_logical_col(dat$is_filler)
    dropped <- sum(fil)
    w <- w[!fil, ]
  }
  list(w = w, traits = traits, genotype_col = gcol, fillers_dropped = dropped)
}

detect_design <- function(w, requested) {
  if (requested != "auto") {
    if (!requested %in% STA_DESIGNS) fail("--design must be auto or one of: ", paste(STA_DESIGNS, collapse = ", "), ".")
    return(list(design = requested, reason = "given with --design"))
  }
  hint <- unique(stats::na.omit(w$design_hint))
  if (length(hint) == 1 && hint %in% names(DESIGN_MAP)) {
    d <- DESIGN_MAP[[hint]]
    if (hint == "latin") d <- if (length(unique(w$repx)) > 1) "res.rowcol" else "rowcol"
    return(list(design = d, reason = sprintf("design column says '%s'", hint)))
  }
  if (length(hint) == 1 && hint == "latin")
    return(list(design = if (length(unique(w$repx)) > 1) "res.rowcol" else "rowcol", reason = "design column says 'latin'"))
  if (length(hint) == 1 && hint %in% c("split", "strip", "factorial"))
    fail("This is a '", hint, "' treatment experiment. claw-trial-analysis fits one genotype factor; analyse factorial treatments with an ANOVA that has the right error strata.")
  has_rep <- !all(is.na(w$repx)) && length(unique(stats::na.omit(w$repx))) > 1
  has_blk <- !all(is.na(w$blockx))
  if (has_rep && has_blk) {
    # Blocks may be numbered across the trial or restart in each replicate;
    # statgenSTA fits them as repId:subBlock either way.
    if (length(unique(paste(w$repx, w$blockx))) > length(unique(w$repx)))
      return(list(design = "res.ibd", reason = "incomplete blocks within replicates"))
  }
  if (has_rep) return(list(design = "rcbd", reason = "replicates without incomplete blocks"))
  if (has_blk) return(list(design = "ibd", reason = "blocks without replicates"))
  if (!all(is.na(w$rowx))) return(list(design = "rowcol", reason = "no replicates or blocks: rows and columns used as blocking factors"))
  fail("No replicate, block or coordinate columns: the design cannot be modelled. Give --rep/--block/--row/--col.")
}

# "none", "complete", or an error message describing incomplete/duplicated positions.
coordinate_state <- function(w) {
  have <- !is.na(w$rowx) & !is.na(w$colx)
  if (!any(have)) return("none")
  if (!all(have))
    return(sprintf("Coordinates are incomplete: %d plot(s) have no row or column (e.g. %s). Spatial and row-column models need every plot's position; fix the data, or use --engine lme4 with a block design.",
                   sum(!have), paste(utils::head(w$plot_id[!have], 5), collapse = ", ")))
  dup <- duplicated(paste(w$rowx, w$colx))
  if (any(dup)) return(sprintf("%d plot(s) share a row/column position with another plot (e.g. %s). Each plot needs a unique position.",
                               sum(dup), paste(utils::head(w$plot_id[dup], 5), collapse = ", ")))
  "complete"
}

replication_info <- function(w, trait) {
  ok <- !is.na(w[[trait]])
  n <- table(w$gen[ok])
  chk <- unique(w$gen[w$is_check])
  test <- setdiff(names(n), chk)
  list(max_test = if (length(test)) max(n[test]) else 0, mean_test = if (length(test)) mean(n[test]) else 0,
       mean_check = if (length(chk)) mean(n[intersect(names(n), chk)]) else 0,
       n_check = length(chk), n_gen = length(n), n_obs = sum(ok))
}

# ---------------------------------------------------------------------------
# Heritability helpers (tested against closed forms)
# ---------------------------------------------------------------------------

# Cullis et al. (2006) generalised heritability for an lme4 model with genotype
# random: 1 - mean PEV of BLUP differences / (2 sigma2_g). The PEV matrix comes
# from the full mixed-model equations in lme4's spherical form
# [X'X, X'Z L; L'Z'X, L'Z'Z L + I], so it allows for the fixed effects.
cullis_h2 <- function(m, term = "genotype") {
  X <- lme4::getME(m, "X")
  Z <- as.matrix(lme4::getME(m, "Z"))
  L <- as.matrix(lme4::getME(m, "Lambda"))
  s2 <- stats::sigma(m)^2
  ZL <- Z %*% L
  C <- rbind(cbind(crossprod(X), crossprod(X, ZL)), cbind(crossprod(ZL, X), crossprod(ZL) + diag(ncol(ZL))))
  Ci <- solve(C)
  p <- ncol(X); q <- ncol(Z)
  Pu <- L %*% (s2 * Ci[p + seq_len(q), p + seq_len(q), drop = FALSE]) %*% t(L)
  gp <- lme4::getME(m, "Gp")
  k <- which(names(lme4::getME(m, "cnms")) == term)
  if (length(k) != 1) fail("Random term '", term, "' not found in the lme4 model.")
  idx <- (gp[k] + 1):gp[k + 1]
  P <- Pu[idx, idx, drop = FALSE]
  n <- length(idx)
  vbar <- 2 * (n * sum(diag(P)) - sum(P)) / (n * (n - 1))
  s2g <- as.numeric(lme4::VarCorr(m)[[term]][1, 1])
  if (s2g <= 1e-10) return(0)
  1 - vbar / (2 * s2g)
}

# SpATS generalised heritability (Oakey et al. 2006): effective dimension of the
# genotype term over its nominal dimension. Same as SpATS::getHeritability()
# without its rounding to two decimals.
spats_h2 <- function(m) {
  unname(m$eff.dim["genotype"] / m$dim.nom["genotype"])
}

# ---------------------------------------------------------------------------
# Model fitting for one trial
# ---------------------------------------------------------------------------

make_td <- function(w, use_check = FALSE) {
  x <- w
  x$checkx <- NA_character_
  if (use_check) {
    x$checkx <- ifelse(x$is_check, x$gen, "noCheck")
    x$gen[x$is_check] <- NA
  }
  # plot_id is kept as an ordinary column: statgenSTA prefixes the trial name
  # onto its own plotId column each time a model is refitted.
  args <- list(data = x, genotype = "gen", trial = "trialx")
  if (!all(is.na(x$repx))) args$repId <- "repx"
  if (!all(is.na(x$blockx))) args$subBlock <- "blockx"
  if (!all(is.na(x$rowx))) { args$rowCoord <- "rowx"; args$colCoord <- "colx" }
  if (use_check) args$checkId <- "checkx"
  do.call(statgenSTA::createTD, args)
}

quiet_fit <- function(expr) {
  msgs <- character()
  res <- withCallingHandlers(suppressMessages(expr),
                             warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(value = res, warnings = unique(msgs))
}

fit_one_trial <- function(w, traits, design, engines, opt) {
  ctrl <- NULL
  if (nzchar(opt$nseg)) {
    ns <- as.numeric(strsplit(opt$nseg, ",", fixed = TRUE)[[1]])
    if (length(ns) != 2 || anyNA(ns)) fail("--nseg must be two numbers, columns,rows (e.g. 10,8).")
    ctrl <- list(nSeg = ns)
  }
  tr <- unique(w$trialx)
  reps <- lapply(traits, function(t) replication_info(w, t)); names(reps) <- traits
  unrep <- vapply(reps, function(r) r$max_test <= 1, TRUE)
  use_check <- any(w$is_check) && all(vapply(reps, function(r) r$n_check > 0 && r$mean_check > r$mean_test, TRUE))
  coords <- !all(is.na(w$rowx))
  for (t in traits) if (unrep[t] && !coords && !use_check)
    fail(sprintf("Trial '%s', trait '%s': test entries are unreplicated and there are no replicated checks or field coordinates. Unreplicated trials need replicated checks or a spatial model (give --row and --col).", tr, t))
  td <- make_td(w)
  # Checks as fixed checkId effects in the genotype-random model: statgenSTA
  # codes check plots with genotype = NA. SpATS keeps those plots; lme4 drops
  # rows with an NA grouping factor (all check plots), so with lme4 the checks
  # stay in the model as random genotypes instead.
  td_r <- if (use_check) make_td(w, use_check = TRUE) else td
  fits <- list()
  for (eng in engines) {
    uc <- use_check && eng == "SpATS"
    ff <- quiet_fit(statgenSTA::fitTD(td, design = design, traits = traits, what = "fixed", engine = eng, control = ctrl))
    fr <- quiet_fit(statgenSTA::fitTD(if (uc) td_r else td, design = design, traits = traits, what = "random", engine = eng,
                                      useCheckId = uc, control = ctrl))
    fits[[eng]] <- list(fixed = ff$value, random = fr$value, use_check = uc,
                        warnings = unique(c(ff$warnings, fr$warnings)))
  }
  list(fits = fits, use_check = use_check, reps = reps, unreplicated = unrep, td = td)
}

# emmeans warns "Model has N prior weights, but we recovered N-1 rows of data"
# for lme4 fits with missing values; the weights are all 1, so it is harmless.
muffle_known <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("prior weights", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
  })
}

# Summaries, BLUEs, BLUPs and outliers for one engine and one trial.
summarise_engine <- function(fit, eng, w, traits, opt) {
  tr <- unique(w$trialx)
  sf <- fit$fixed; sr <- fit$random
  ef <- muffle_known(statgenSTA::extractSTA(sf, what = c("BLUEs", "seBLUEs", "CV", "rDfF")))[[tr]]
  er <- muffle_known(statgenSTA::extractSTA(sr, what = c("BLUPs", "seBLUPs", "varGen", "varErr")))[[tr]]
  rows <- list(); blues <- list(); blups <- list(); outl <- list()
  for (t in traits) {
    mf <- sf[[tr]]$mFix[[t]]; mr <- sr[[tr]]$mRand[[t]]
    if (is.null(mf) || is.null(mr)) {
      rows[[t]] <- data.frame(trial = tr, trait = t, model = eng, status = "model failed", stringsAsFactors = FALSE)
      next
    }
    b <- data.frame(trial = tr, genotype = as.character(ef$BLUEs$genotype), trait = t, model = eng,
                    blue = ef$BLUEs[[t]], se_blue = ef$seBLUEs[[t]][match(ef$BLUEs$genotype, ef$seBLUEs$genotype)],
                    stringsAsFactors = FALSE)
    b <- b[!is.na(b$blue), ]
    p <- data.frame(trial = tr, genotype = as.character(er$BLUPs$genotype), trait = t, model = eng,
                    blup = er$BLUPs[[t]], se_blup = er$seBLUPs[[t]][match(er$BLUPs$genotype, er$seBLUPs$genotype)],
                    stringsAsFactors = FALSE)
    p <- p[!is.na(p$blup) & !is.na(p$genotype), ]
    if (eng == "SpATS") {
      h2 <- spats_h2(mr)
      h2_def <- "generalised H2 (Oakey et al. 2006) = ED_genotype / nominal dimension, from SpATS"
      h2_statgen <- h2
    } else {
      h2 <- cullis_h2(mr)
      h2_def <- "generalised H2 (Cullis et al. 2006) = 1 - mean PEV of BLUP differences / (2 sigma2_g), from the lme4 model"
      h2_statgen <- muffle_known(statgenSTA::extractSTA(sr, what = "heritability"))[[t]]
    }
    rl <- if (opt$outlier_limit > 0) opt$outlier_limit else NULL
    o <- tryCatch(muffle_known(statgenSTA::outlierSTA(sf, traits = t, what = "fixed", rLimit = rl, verbose = FALSE))$outliers,
                  error = function(e) NULL)
    if (!is.null(o) && nrow(o)) {
      o <- o[o$outlier, , drop = FALSE]
      pid <- if ("plot_id" %in% names(o)) as.character(o$plot_id) else NA_character_
      outl[[t]] <- data.frame(trial = tr, trait = t, model = eng, plot_id = pid,
                              genotype = as.character(o$genotype),
                              row = if ("rowCoord" %in% names(o)) o$rowCoord else NA,
                              col = if ("colCoord" %in% names(o)) o$colCoord else NA,
                              value = o$value, std_residual = round(o$res, 3), stringsAsFactors = FALSE)
    }
    rdf <- ef$rDfF[[t]]
    lim <- if (!is.null(rl)) rl else min(4, max(2, stats::qnorm(1 - 0.5 / rdf)))
    rows[[t]] <- data.frame(trial = tr, trait = t, model = eng, status = "ok",
                            n_obs = sum(!is.na(w[[t]])), n_genotypes = length(unique(w$gen[!is.na(w[[t]])])),
                            mean = mean(w[[t]], na.rm = TRUE), cv_pct = round(ef$CV[[t]], 2),
                            var_genetic = signif(er$varGen[[t]], 4), var_residual = signif(er$varErr[[t]], 4),
                            h2 = round(h2, 3), h2_definition = h2_def,
                            h2_statgenSTA = round(h2_statgen, 3),
                            mean_se_blue = signif(mean(b$se_blue, na.rm = TRUE), 4),
                            residual_df = round(rdf), outlier_limit = round(lim, 2),
                            n_outliers = if (is.null(outl[[t]])) 0L else nrow(outl[[t]]),
                            stringsAsFactors = FALSE)
    blues[[t]] <- b; blups[[t]] <- p
  }
  list(summary = bind_rows_fill(rows), blues = do.call(rbind, blues), blups = do.call(rbind, blups),
       outliers = do.call(rbind, outl))
}

bind_rows_fill <- function(lst) {
  lst <- Filter(Negate(is.null), lst)
  if (!length(lst)) return(NULL)
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(x) { for (cl in setdiff(cols, names(x))) x[[cl]] <- NA; x[, cols] }))
}

analyse <- function(dat, opt) {
  prep <- prepare_data(dat, opt)
  w_all <- prep$w
  traits <- prep$traits
  eng_req <- tolower(opt$engine)
  if (!eng_req %in% c("auto", "both", "spats", "lme4")) fail("--engine must be auto, both, SpATS or lme4.")
  out <- list(summary = list(), blues = list(), blups = list(), outliers = list(), comparison = list(),
              warnings = character(), info = list(), fits = list())
  for (tr in unique(w_all$trialx)) {
    w <- w_all[w_all$trialx == tr, ]
    cs <- coordinate_state(w)
    dd <- detect_design(w, opt$design)
    needs_xy <- eng_req != "lme4" || dd$design %in% c("rowcol", "res.rowcol")
    if (!cs %in% c("none", "complete")) {
      if (needs_xy) fail(cs)
      out$warnings <- c(out$warnings, sprintf("Trial %s: %s They were ignored (the lme4 block model does not use them).", tr, cs))
      w$rowx <- NA_real_; w$colx <- NA_real_
      cs <- "none"
    }
    coords <- cs == "complete"
    if (dd$design %in% c("rcbd", "res.ibd", "res.rowcol") && all(is.na(w$repx))) fail("Design ", dd$design, " needs a replicate column (--rep).")
    if (dd$design %in% c("ibd", "res.ibd") && all(is.na(w$blockx))) fail("Design ", dd$design, " needs a block column (--block).")
    if (dd$design %in% c("rowcol", "res.rowcol") && !coords) fail("Design ", dd$design, " needs row and column coordinates.")
    engines <- switch(eng_req, auto = if (coords) c("SpATS", "lme4") else "lme4", both = c("SpATS", "lme4"),
                      spats = "SpATS", lme4 = "lme4")
    if ("SpATS" %in% engines && !coords) fail("SpATS (spatial) needs row and column coordinates for every plot; give --row/--col or use --engine lme4.")
    ft <- fit_one_trial(w, traits, dd$design, engines, opt)
    out$fits[[tr]] <- ft
    out$info[[tr]] <- list(design = dd$design, reason = dd$reason, engines = engines, use_check = ft$use_check,
                           checkid_engines = names(Filter(function(f) isTRUE(f$use_check), ft$fits)),
                           checks_replicated = vapply(ft$reps, function(r) r$n_check > 0 && r$mean_check > r$mean_test, TRUE),
                           unreplicated = ft$unreplicated, coords = coords, n_plots = nrow(w))
    for (eng in engines) {
      s <- summarise_engine(ft$fits[[eng]], eng, w, traits, opt)
      out$summary[[paste(tr, eng)]] <- s$summary
      out$blues[[paste(tr, eng)]] <- s$blues
      out$blups[[paste(tr, eng)]] <- s$blups
      out$outliers[[paste(tr, eng)]] <- s$outliers
      if (length(ft$fits[[eng]]$warnings))
        out$warnings <- c(out$warnings, sprintf("Trial %s, %s: %s", tr, eng, ft$fits[[eng]]$warnings))
    }
    # Spatial vs non-spatial comparison
    if (all(c("SpATS", "lme4") %in% engines)) {
      for (t in traits) {
        S <- out$summary[[paste(tr, "SpATS")]]; L <- out$summary[[paste(tr, "lme4")]]
        S <- S[S$trait == t, ]; L <- L[L$trait == t, ]
        bs <- out$blues[[paste(tr, "SpATS")]]; bl <- out$blues[[paste(tr, "lme4")]]
        if (is.null(bs) || is.null(bl) || !nrow(S) || !nrow(L) || S$status != "ok" || L$status != "ok") next
        bs <- bs[bs$trait == t, ]; bl <- bl[bl$trait == t, ]
        m <- match(bs$genotype, bl$genotype)
        out$comparison[[paste(tr, t)]] <- data.frame(
          trial = tr, trait = t, h2_spatial = S$h2, h2_nonspatial = L$h2,
          mean_se_blue_spatial = S$mean_se_blue, mean_se_blue_nonspatial = L$mean_se_blue,
          cor_blues = round(stats::cor(bs$blue, bl$blue[m], use = "complete.obs"), 4),
          outliers_spatial = S$n_outliers, outliers_nonspatial = L$n_outliers,
          more_precise = if (S$mean_se_blue < L$mean_se_blue) "spatial (SpATS)" else "non-spatial (lme4)",
          stringsAsFactors = FALSE)
      }
    }
    if (prep$fillers_dropped > 0) out$warnings <- c(out$warnings, sprintf("%d FILLER plot(s) were dropped before analysis.", prep$fillers_dropped))
  }
  out$summary <- bind_rows_fill(out$summary)
  out$blues <- bind_rows_fill(out$blues)
  out$blups <- bind_rows_fill(out$blups)
  out$outliers <- bind_rows_fill(out$outliers)
  out$comparison <- bind_rows_fill(out$comparison)
  out$traits <- traits
  out$prep <- prep
  # BLUEs for MET / GWAS
  pick <- tolower(opt$blues_from)
  if (!pick %in% c("auto", "spats", "lme4")) fail("--blues-from must be auto, SpATS or lme4.")
  met <- NULL
  for (tr in names(out$info)) {
    eng <- if (pick == "auto") out$info[[tr]]$engines[1] else if (pick == "spats") "SpATS" else "lme4"
    if (!eng %in% out$info[[tr]]$engines) fail("--blues-from ", eng, " was not fitted for trial ", tr, ".")
    out$info[[tr]]$blues_from <- eng
    b <- out$blues[out$blues$trial == tr & out$blues$model == eng, ]
    wide <- data.frame(trial = tr, genotype = sort(unique(b$genotype)), stringsAsFactors = FALSE)
    for (t in traits) {
      bt <- b[b$trait == t, ]
      m <- match(wide$genotype, bt$genotype)
      wide[[t]] <- bt$blue[m]
      wide[[paste0("se_", t)]] <- bt$se_blue[m]
      wide[[paste0("wt_", t)]] <- 1 / bt$se_blue[m]^2
    }
    met <- rbind(met, wide)
  }
  out$met <- met
  out
}

# ---------------------------------------------------------------------------
# Warnings, figures and report
# ---------------------------------------------------------------------------

analysis_warnings <- function(res) {
  w <- character()
  s <- res$summary
  for (i in seq_len(nrow(s))) {
    r <- s[i, ]
    if (r$status != "ok") { w <- c(w, sprintf("Trial %s, %s, %s: the model failed; see the console warnings.", r$trial, r$trait, r$model)); next }
    if (r$var_genetic <= 1e-8 * max(1, abs(r$var_residual)))
      w <- c(w, sprintf("Trial %s, %s (%s): genetic variance estimated as zero, so BLUPs are all equal and H2 = 0. Check the data before using this trial.", r$trial, r$trait, r$model))
    else if (r$h2 < 0.2)
      w <- c(w, sprintf("Trial %s, %s (%s): H2 = %.2f is low (rule of thumb < 0.2); genotype means are poorly separated.", r$trial, r$trait, r$model, r$h2))
    if (!is.na(r$cv_pct) && r$cv_pct > 20)
      w <- c(w, sprintf("Trial %s, %s (%s): CV = %.1f%% is high (rule of thumb > 20%% for yield); check units and outliers.", r$trial, r$trait, r$model, r$cv_pct))
    if (r$n_outliers > 0)
      w <- c(w, sprintf("Trial %s, %s (%s): %d outlier(s) with |standardised residual| > %.2f; listed in outliers.csv, not removed.", r$trial, r$trait, r$model, r$n_outliers, r$outlier_limit))
  }
  for (tr in names(res$info)) {
    inf <- res$info[[tr]]
    if (any(inf$unreplicated))
      w <- c(w, sprintf("Trial %s: test entries are unreplicated; genotype effects rest on %s.", tr,
                        paste(c(if (any(inf$checks_replicated)) "replicated checks" else NULL,
                                if (inf$coords) "the spatial model" else NULL), collapse = " and ")))
    if (isTRUE(inf$use_check) && "lme4" %in% inf$engines)
      w <- c(w, sprintf("Trial %s: in the lme4 BLUP/H2 model the checks are random genotypes (statgenSTA's checkId option drops all check plots with lme4); the SpATS model fits them as fixed checkId effects. BLUEs are unaffected.", tr))
    if (!inf$coords) w <- c(w, sprintf("Trial %s: no field coordinates, so no spatial model was fitted.", tr))
  }
  c(w, res$warnings)
}

DIV <- c("#2a78d6", "#5598e7", "#86b6ef", "#b7d3f6", "#f0efec", "#f5c2c1", "#ee9796", "#e86b6a", "#e34948")
SEQ <- c("#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b")
INK <- "#0b0b0b"; INK2 <- "#52514e"; SURFACE <- "#fcfcfb"

grid_panel <- function(x, y, z, pal, diverging, main, sub) {
  graphics::plot.new()
  graphics::plot.window(xlim = range(x) + c(-0.5, 0.5), ylim = rev(range(y) + c(-0.5, 0.5)), xaxs = "i", yaxs = "i")
  lim <- if (diverging) c(-1, 1) * max(abs(z), na.rm = TRUE) else range(z, na.rm = TRUE)
  if (!all(is.finite(lim)) || diff(lim) == 0) lim <- lim + c(-1, 1)
  br <- seq(lim[1], lim[2], length.out = length(pal) + 1)
  cl <- pal[pmax(1, pmin(length(pal), findInterval(z, br, rightmost.closed = TRUE)))]
  cl[is.na(z)] <- "#ffffff"
  graphics::rect(x - 0.5, y + 0.5, x + 0.5, y - 0.5, col = cl, border = NA)
  graphics::axis(2, las = 1, tick = FALSE, cex.axis = 0.65, line = -0.8, col.axis = INK2)
  graphics::axis(3, tick = FALSE, cex.axis = 0.65, line = -0.8, col.axis = INK2)
  graphics::mtext(main, side = 1, line = 0.5, adj = 0, cex = 0.8, font = 2, col = INK)
  graphics::mtext(sub, side = 1, line = 1.4, adj = 0, cex = 0.58, col = INK2)
}

draw_spatial <- function(res, out_dir) {
  files <- character()
  for (tr in names(res$fits)) {
    ft <- res$fits[[tr]]
    if (!"SpATS" %in% names(ft$fits)) next
    for (t in res$traits) {
      m <- ft$fits$SpATS$fixed[[tr]]$mFix[[t]]
      if (is.null(m)) next
      d <- m$data
      st <- SpATS::obtain.spatialtrend(m, grid = c(length(unique(d$colCoord)), length(unique(d$rowCoord))))
      tg <- expand.grid(row = st$row.p, col = st$col.p)
      tg$z <- as.vector(st$fit)
      rsd <- stats::residuals(m)
      o <- res$outliers
      o <- if (is.null(o)) NULL else o[o$trial == tr & o$trait == t & o$model == "SpATS", ]
      nm <- sprintf("spatial_%s_%s.png", gsub("[^A-Za-z0-9_.-]", "_", tr), gsub("[^A-Za-z0-9_.-]", "_", t))
      R <- length(unique(d$rowCoord)); C <- length(unique(d$colCoord))
      open_png(out_dir, nm, width = min(18, 3 * max(3.2, C * 0.22) + 0.8), height = min(10, max(3.2, R * 0.22) + 1.6))
      graphics::par(mfrow = c(1, 3), mar = c(2.6, 2.2, 1.6, 0.6), bg = SURFACE)
      grid_panel(d$colCoord, d$rowCoord, d[[t]], SEQ, FALSE, paste0(t, ": raw data"),
                 sprintf("light = low (%.4g), dark = high (%.4g)", min(d[[t]], na.rm = TRUE), max(d[[t]], na.rm = TRUE)))
      grid_panel(tg$col, tg$row, tg$z, DIV, TRUE, "fitted spatial trend (SpATS)",
                 sprintf("blue = below average, red = above; range %.3g to %.3g", min(tg$z), max(tg$z)))
      grid_panel(d$colCoord, d$rowCoord, rsd, DIV, TRUE, "residuals after genotype + trend",
                 "rings = outliers (outlierSTA); a pattern left here means the model missed structure")
      if (!is.null(o) && nrow(o)) graphics::points(o$col, o$row, pch = 1, cex = 2.2, lwd = 2, col = INK)
      grDevices::dev.off()
      files <- c(files, nm)
    }
  }
  files
}

draw_comparison <- function(res, out_dir) {
  if (is.null(res$comparison) || !nrow(res$comparison)) return(character())
  files <- character()
  for (i in seq_len(nrow(res$comparison))) {
    cmp <- res$comparison[i, ]
    b <- res$blues[res$blues$trial == cmp$trial & res$blues$trait == cmp$trait, ]
    s <- b[b$model == "SpATS", ]; l <- b[b$model == "lme4", ]
    m <- match(s$genotype, l$genotype)
    nm <- sprintf("blues_spatial_vs_nonspatial_%s_%s.png", gsub("[^A-Za-z0-9_.-]", "_", cmp$trial), gsub("[^A-Za-z0-9_.-]", "_", cmp$trait))
    open_png(out_dir, nm, width = 5.2, height = 5.2)
    graphics::par(mar = c(4.2, 4.2, 2.6, 1), bg = SURFACE, fg = INK2, col.axis = INK2, col.lab = INK2)
    rng <- range(c(s$blue, l$blue[m]), na.rm = TRUE)
    graphics::plot(l$blue[m], s$blue, xlim = rng, ylim = rng, pch = 16, col = grDevices::adjustcolor("#2a78d6", 0.75), cex = 0.9,
                   xlab = paste("BLUE, non-spatial (lme4):", cmp$trait), ylab = paste("BLUE, spatial (SpATS):", cmp$trait), bty = "l")
    graphics::abline(0, 1, col = INK2, lty = 2)
    graphics::mtext(sprintf("Trial %s: r = %.3f; H2 spatial %.2f vs non-spatial %.2f", cmp$trial, cmp$cor_blues, cmp$h2_spatial, cmp$h2_nonspatial),
                    side = 3, line = 0.8, adj = 0, cex = 0.8, col = INK, font = 2)
    grDevices::dev.off()
    files <- c(files, nm)
  }
  files
}

md_table <- function(df) {
  df[] <- lapply(df, function(x) { x <- as.character(x); x[is.na(x)] <- ""; gsub("\\|", "/", x) })
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

REFERENCES <- c(
  "- van Rossum B (2025). statgenSTA: Single Trial Analysis (STA) of Field Trials. R package version 1.0.15. https://doi.org/10.32614/CRAN.package.statgenSTA",
  "- Rodriguez-Alvarez MX, Boer MP, van Eeuwijk FA, Eilers PHC (2018). Correcting for spatial heterogeneity in plant breeding experiments with P-splines. Spatial Statistics 23:52-71. https://doi.org/10.1016/j.spasta.2017.10.003",
  "- Bates D, Maechler M, Bolker B, Walker S (2015). Fitting linear mixed-effects models using lme4. Journal of Statistical Software 67(1):1-48. https://doi.org/10.18637/jss.v067.i01",
  "- Cullis BR, Smith AB, Coombes NE (2006). On the design of early generation variety trials with correlated data. Journal of Agricultural, Biological, and Environmental Statistics 11(4):381-393. https://doi.org/10.1198/108571106X154443",
  "- Oakey H, Verbyla A, Pitchford W, Cullis B, Kuchel H (2006). Joint modeling of additive and non-additive genetic line effects in single field trials. Theoretical and Applied Genetics 113(5):809-819. https://doi.org/10.1007/s00122-006-0333-z",
  "- Piepho H-P, Mohring J (2007). Computing heritability and selection response from unbalanced plant breeding trials. Genetics 177(3):1881-1888. https://doi.org/10.1534/genetics.107.074229",
  "- Gilmour AR, Cullis BR, Verbyla AP (1997). Accounting for natural and extraneous variation in the analysis of field experiments. Journal of Agricultural, Biological, and Environmental Statistics 2(3):269-293. https://doi.org/10.2307/1400446",
  "- Smith AB, Cullis BR, Thompson R (2005). The analysis of crop cultivar breeding and evaluation trials: an overview of current mixed model approaches. Journal of Agricultural Science 143(6):449-462. https://doi.org/10.1017/S0021859605005587",
  "- Piepho H-P, Mohring J, Schulz-Streeck T, Ogutu JO (2012). A stage-wise approach for the analysis of multi-environment trials. Biometrical Journal 54(6):844-860. https://doi.org/10.1002/bimj.201100219",
  "- Garrick DJ, Taylor JF, Fernando RL (2009). Deregressing estimated breeding values and weighting information for genomic regression analyses. Genetics Selection Evolution 41:55. https://doi.org/10.1186/1297-9686-41-55")

write_report <- function(res, out_dir, warn, figs, opt) {
  L <- c(sprintf("# %s report", SKILL), "")
  if (length(warn)) L <- c(L, "## Warnings", "", paste("-", warn), "")
  L <- c(L, "## Trials and models", "")
  for (tr in names(res$info)) {
    inf <- res$info[[tr]]
    L <- c(L, sprintf("- **%s**: %d plots; design `%s` (%s); engines: %s; BLUEs exported from **%s**%s.", tr, inf$n_plots, inf$design, inf$reason,
                      paste(inf$engines, collapse = ", "), inf$blues_from,
                      if (length(inf$checkid_engines)) sprintf("; checks fitted as fixed checkId effects in the %s BLUP model", paste(inf$checkid_engines, collapse = "/")) else ""))
  }
  L <- c(L, "", "Models (statgenSTA::fitTD): SpATS = design terms + 2-D P-spline surface PSANOVA(col, row) + genotype; lme4 = design terms + genotype. Genotype is fixed for BLUEs and random for BLUPs and H2.", "",
         "## Results", "")
  show <- res$summary[, intersect(c("trial", "trait", "model", "n_obs", "n_genotypes", "mean", "cv_pct", "var_genetic", "var_residual",
                                    "h2", "h2_statgenSTA", "mean_se_blue", "n_outliers"), names(res$summary))]
  show$mean <- signif(show$mean, 5)
  L <- c(L, md_table(show), "",
         "- `h2` for SpATS: generalised heritability of Oakey et al. (2006), effective dimension of genotype / nominal dimension (SpATS::getHeritability, unrounded).",
         "- `h2` for lme4: generalised heritability of Cullis et al. (2006), 1 - mean PEV of BLUP differences / (2 sigma2_g), from the full mixed-model equations. `h2_statgenSTA` for lme4 is statgenSTA's sigma2_g / (sigma2_g + sigma2_e / r) (r = number of replicates; plot-level when there is no replicate column); it equals the Cullis value only for balanced complete-block data.",
         "- `cv_pct` = 100 x sqrt(residual variance) / mean (statgenSTA). For SpATS the residual is what is left after the spatial trend, so its CV is lower than an ANOVA CV.", "")
  if (!is.null(res$comparison) && nrow(res$comparison)) {
    L <- c(L, "## Spatial vs non-spatial", "", md_table(res$comparison), "",
           "Compared on precision (mean standard error of BLUEs; lower is better) and generalised heritability (higher is better). Likelihood criteria (AIC/BIC) are not comparable between SpATS and lme4 fits, so they are not used. A BLUE correlation near 1 means the choice changes little; a clear gain in H2 or standard errors favours the spatial model.", "")
  }
  if (!is.null(res$outliers) && nrow(res$outliers)) L <- c(L, "## Outliers", "", md_table(utils::head(res$outliers, 30)), "")
  if (length(figs)) L <- c(L, "## Figures", "", sprintf("![%s](figures/%s)", figs, figs), "")
  L <- c(L, "## Files", "",
         "- `tables/blues.csv`, `tables/blups.csv`: per trial, genotype, trait and model, with standard errors.",
         "- `tables/model_summary.csv`, `tables/model_comparison.csv`, `tables/outliers.csv`.",
         "- `tables/met_input.csv`: one row per trial x genotype with BLUE, `se_<trait>` and `wt_<trait>` = 1/se^2. Row-bind these files from several trials and pass them to `claw-met-gxe`; use BLUEs (not BLUPs) for GWAS too.", "",
         "## References", "", REFERENCES, "", sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(L, file.path(out_dir, "report.md"))
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (opt$demo) {
    opt$input <- file.path(SKILL_DIR, "demo", "trial_clean.csv")
    if (!nzchar(opt$traits)) opt$traits <- "yield_kg_ha,plant_height_cm"
  }
  if (!nzchar(opt$input)) fail("Provide --input DATA.csv and --traits, or --demo. See --help.")
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  dat <- read_trial_data(opt$input)
  res <- analyse(dat, opt)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  save_table(res$summary, opt$out, "model_summary.csv")
  save_table(res$blues, opt$out, "blues.csv")
  save_table(res$blups, opt$out, "blups.csv")
  save_table(if (is.null(res$outliers)) data.frame(trial = character(), trait = character(), model = character(), plot_id = character()) else res$outliers,
             opt$out, "outliers.csv")
  if (!is.null(res$comparison)) save_table(res$comparison, opt$out, "model_comparison.csv")
  save_table(res$met, opt$out, "met_input.csv")
  figs <- c(draw_spatial(res, opt$out), draw_comparison(res, opt$out))
  warn <- analysis_warnings(res)
  write_report(res, opt$out, warn, figs, opt)
  write_repro_bundle(opt$out, argv, inputs = opt$input)
  cat(sprintf("Done: %d trial(s), %d trait(s). Report: %s\n", length(res$info), length(res$traits),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  if (length(warn)) cat(paste("Warning:", warn), sep = "\n")
  invisible(res)
}

if (sys.nframe() == 0L) main()
