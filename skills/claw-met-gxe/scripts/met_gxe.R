#!/usr/bin/env Rscript
# claw-met-gxe: multi-environment trial (MET) analysis of per-environment
# genotype BLUEs with statgenGxE: variance components and heritability across
# environments, MET genotype BLUPs, AMMI and GGE biplots, Finlay-Wilkinson
# regression, stability measures (Lin-Binns superiority, Wricke's ecovalence,
# Shukla's stability variance, static stability) and mega-environments. When
# the BLUEs come with standard errors, a stage-2 model with known stage-1
# variances (nlme; Piepho et al. 2012) gives weighted BLUPs and splits GxE
# from stage-1 error.
#
# Usage:
#   Rscript scripts/met_gxe.R --demo
#   Rscript scripts/met_gxe.R --input met_input.csv --trait yield_kg_ha
#   Rscript scripts/met_gxe.R --input a/met_input.csv,b/met_input.csv --trait yield_kg_ha --best max
#   Rscript scripts/met_gxe.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT license.

SKILL   <- "claw-met-gxe"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/met_gxe.R --input FILE.csv[,FILE2.csv] --trait COL [options]

  --input FILES     one or more CSVs (comma-separated, stacked): one row per genotype
                    per environment, e.g. tables/met_input.csv from claw-trial-analysis
  --trait COL       trait column (per-environment BLUEs), e.g. yield_kg_ha
  --env COL         environment column (default trial)
  --genotype COL    genotype column (default genotype)
  --weight COL      weight column, 1/SE^2 (default wt_<trait>; else 1/se_<trait>^2)
  --no-weights      ignore weights / standard errors
  --nesting COL     trials nested in this column (e.g. region or scenario) for the
                    variance-component model
  --location COL    location column; with --year fits locations x years (each trial
  --year COL        must be one location-year)
  --npc N           principal components for AMMI and GGE (default 2)
  --best S          max (default) or min: the direction that is better for the trait
                    (superiority measure, mega-environment winners, rankings)
  --complete-only   use only genotypes observed in every environment for AMMI, GGE,
                    Finlay-Wilkinson and mega-environments (stability always does)
  --demo            synthetic 40 genotypes x 8 trials in demo/met_input.csv
  --out DIR         output folder (default: agriclaw_output/claw-met-gxe)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, input = "", trait = "", env = "trial", genotype = "genotype",
                 weight = "", no_weights = FALSE, nesting = "", location = "", year = "",
                 npc = 2, best = "max", complete_only = FALSE, out = "")

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

fail <- function(...) stop(paste0(...), call. = FALSE)

# Stack several CSVs; columns missing from a file are filled with NA.
read_inputs <- function(paths) {
  paths <- trimws(strsplit(paths, ",", fixed = TRUE)[[1]])
  paths <- paths[nzchar(paths)]
  if (!length(paths)) fail("Give --input FILE.csv (per-environment BLUEs), or --demo. See --help.")
  dl <- lapply(paths, function(p) {
    if (!file.exists(p)) fail("Input file not found: ", p)
    x <- utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
    names(x) <- trimws(names(x))
    x
  })
  cols <- unique(unlist(lapply(dl, names)))
  dl <- lapply(dl, function(x) { for (cl in setdiff(cols, names(x))) x[[cl]] <- NA; x[cols] })
  list(data = do.call(rbind, dl), paths = paths)
}

# Standardized MET table: env, genotype, y, optional se2 (stage-1 variance),
# optional location/year/nesting columns. Fails loudly on structural problems.
prepare_met <- function(raw, opt) {
  for (cl in c(opt$env, opt$genotype, opt$trait))
    if (!nzchar(cl) || !cl %in% names(raw))
      fail("Column '", cl, "' not found. Columns: ", paste(names(raw), collapse = ", "),
           ". Set --env, --genotype and --trait.")
  y <- raw[[opt$trait]]
  if (!is.numeric(y)) {
    yn <- suppressWarnings(as.numeric(as.character(y)))
    bad <- !is.na(y) & is.na(yn) & !trimws(as.character(y)) %in% c("", "NA")
    if (any(bad)) fail("Trait '", opt$trait, "' has non-numeric values, e.g. ", paste(utils::head(unique(y[bad]), 3), collapse = ", "), ".")
    y <- yn
  }
  d <- data.frame(env = trimws(as.character(raw[[opt$env]])), genotype = trimws(as.character(raw[[opt$genotype]])),
                  y = y, stringsAsFactors = FALSE)
  if (any(is.na(d$env) | d$env == "" | is.na(d$genotype) | d$genotype == ""))
    fail("Some rows have a blank environment or genotype.")
  notes <- character()
  # stage-1 variances for weights
  wcol <- if (nzchar(opt$weight)) opt$weight else paste0("wt_", opt$trait)
  secol <- paste0("se_", opt$trait)
  d$se2 <- NA_real_
  wsrc <- "none"
  if (!opt$no_weights) {
    if (wcol %in% names(raw)) {
      w <- suppressWarnings(as.numeric(raw[[wcol]]))
      d$se2 <- ifelse(is.finite(w) & w > 0, 1 / w, NA)
      wsrc <- wcol
    } else if (nzchar(opt$weight)) {
      fail("Weight column '", opt$weight, "' not found.")
    } else if (secol %in% names(raw)) {
      s <- suppressWarnings(as.numeric(raw[[secol]]))
      d$se2 <- ifelse(is.finite(s) & s > 0, s^2, NA)
      wsrc <- paste0("1/", secol, "^2")
    }
  }
  for (k in c("location", "year", "nesting")) {
    cl <- opt[[k]]
    if (nzchar(cl)) {
      if (!cl %in% names(raw)) fail("--", k, " column '", cl, "' not found.")
      d[[k]] <- trimws(as.character(raw[[cl]]))
    }
  }
  nmiss <- sum(is.na(d$y))
  if (nmiss) notes <- c(notes, sprintf("%d row(s) with a missing %s were removed before analysis.", nmiss, opt$trait))
  d <- d[!is.na(d$y), ]
  dup <- duplicated(d[c("env", "genotype")])
  if (any(dup))
    fail(sprintf("%d genotype x environment combination(s) occur more than once (e.g. %s in %s). The input must be one BLUE per genotype per environment (claw-trial-analysis output), not plot data.",
                 sum(dup), d$genotype[dup][1], d$env[dup][1]))
  ne <- length(unique(d$env)); ng <- length(unique(d$genotype))
  if (ne < 2) fail("Only ", ne, " environment: a MET analysis needs at least 2 (3 for AMMI/GGE, stability and mega-environments).")
  if (ng < 3) fail("Only ", ng, " genotypes; a MET analysis needs at least 3.")
  if (wsrc != "none" && anyNA(d$se2)) {
    notes <- c(notes, sprintf("%d row(s) lack a usable weight (%s); the weighted stage-2 model was skipped.", sum(is.na(d$se2)), wsrc))
    d$se2 <- NA_real_
    wsrc <- "none"
  }
  if (!is.null(d$nesting)) {
    k <- tapply(d$nesting, d$env, function(z) length(unique(z)))
    if (any(k > 1)) fail("--nesting column must hold one value per environment; ", names(k)[k > 1][1], " has several.")
  }
  if (xor(nzchar(opt$location), nzchar(opt$year))) fail("Give both --location and --year (or neither).")
  if (nzchar(opt$location)) {
    ly <- unique(d[c("env", "location", "year")])
    if (anyDuplicated(ly$env)) fail("Each environment must have a single location and year.")
    if (anyDuplicated(ly[c("location", "year")]))
      fail("Several environments share a location-year (e.g. water regimes); the location x year model needs one trial per location-year. Use --nesting instead.")
  }
  list(dat = d, notes = notes, weight_source = wsrc)
}

completeness <- function(d) {
  tab <- table(d$genotype, d$env)
  full <- rownames(tab)[rowSums(tab > 0) == ncol(tab)]
  list(cells = length(tab), missing = sum(tab == 0), complete_genotypes = full,
       envs_per_genotype = rowSums(tab > 0), table = tab)
}

make_td <- function(d, trait) {
  x <- data.frame(genotype = d$genotype, trial = d$env, stringsAsFactors = FALSE)
  x[[trait]] <- d$y
  args <- list(data = x, genotype = "genotype", trial = "trial")
  if (!is.null(d$location)) { x$loc <- d$location; x$year <- d$year; args$data <- x; args$loc <- "loc"; args$year <- "year" }
  if (!is.null(d$nesting)) { x$nest <- d$nesting; args$data <- x }
  do.call(statgenSTA::createTD, args)
}

# ---------------------------------------------------------------------------
# Analyses
# ---------------------------------------------------------------------------

quiet <- function(expr) {
  w <- character()
  val <- withCallingHandlers(suppressMessages(expr), warning = function(cw) {
    w <<- c(w, conditionMessage(cw))
    invokeRestart("muffleWarning")
  })
  attr(val, "warnings") <- unique(trimws(w))
  val
}

# Variance components, heritability and genotype BLUPs (statgenGxE::gxeVarComp, lme4).
run_varcomp <- function(td, d, trait) {
  args <- list(TD = td, trait = trait)
  model <- "trait = trial + genotype (random) + residual (= G x E + stage-1 error)"
  if (!is.null(d$nesting)) {
    args$nestingFactor <- "nest"
    model <- "trait = nest + nest:trial + genotype + genotype:nest (random) + residual"
  }
  if (!is.null(d$location)) {
    args$locationYear <- TRUE
    model <- "trait = year + loc + year:loc + genotype + genotype:year + genotype:loc (random) + residual"
  }
  vcm <- quiet(do.call(statgenGxE::gxeVarComp, args))
  comps <- statgenGxE::vc(vcm)
  full <- vcm$fullRandVC
  tab <- data.frame(component = rownames(comps), variance = comps$Component, stringsAsFactors = FALSE)
  tab$component[tab$component == "residuals"] <- "residual"
  tab$component <- sub("nest", "nesting", tab$component, fixed = TRUE)
  tot <- data.frame(term = sub("nest", "nesting", rownames(full), fixed = TRUE), variance = full$vcov,
                    percent_of_total = round(100 * full$vcovPerc, 1), stringsAsFactors = FALSE)
  tot$term[tot$term == "residuals"] <- "residual"
  h2 <- statgenGxE::herit(vcm)
  # predict.varComp calls emmeans, whose default small-sample df for lmer
  # (Kenward-Roger) is slow and does not change the predicted values;
  # asymptotic df gives identical predictions (mean trial effect + genotype BLUP).
  old <- emmeans::get_emm_option("lmer.df")
  emmeans::emm_options(lmer.df = "asymptotic")
  on.exit(emmeans::emm_options(lmer.df = old), add = TRUE)
  pred <- quiet(stats::predict(vcm, predictLevel = "genotype"))
  list(fit = vcm, components = tab, full_random = tot, h2 = h2, blup = pred, model = model,
       warnings = attr(vcm, "warnings"))
}

# Stage 2 with known stage-1 variances (Piepho et al. 2012): residual variance of
# each BLUE fixed at its squared SE (nlme: varFixed + sigma = 1), so the G x E
# variance is estimated separately from stage-1 error.
run_weighted_stage2 <- function(d) {
  x <- data.frame(y = d$y, env = factor(d$env), genotype = factor(d$genotype), se2 = d$se2)
  m <- nlme::lme(y ~ env, random = ~ 1 | genotype / env, data = x,
                 weights = nlme::varFixed(~ se2), control = nlme::lmeControl(sigma = 1, returnObject = TRUE))
  vcx <- nlme::VarCorr(m)
  s2g <- as.numeric(vcx[grep("^genotype", rownames(vcx)) + 1L, "Variance"][1])
  s2ge <- as.numeric(vcx[grep("^env", rownames(vcx)) + 1L, "Variance"][1])
  re <- nlme::ranef(m)$genotype
  mu <- mean(stats::predict(m, newdata = data.frame(env = levels(x$env)), level = 0))
  ehm <- length(unique(x$genotype)) / sum(1 / table(x$genotype))
  v1 <- mean(x$se2)
  list(sigma2_g = s2g, sigma2_ge = s2ge, mean_stage1_var = v1, env_harmonic_mean = ehm,
       h2 = s2g / (s2g + (s2ge + v1) / ehm),
       blup = data.frame(genotype = rownames(re), blup_weighted = mu + re[, 1], stringsAsFactors = FALSE))
}

ammi_like <- function(td, trait, npc, gge, exclude) {
  f <- if (gge) statgenGxE::gxeGGE else statgenGxE::gxeAmmi
  res <- tryCatch(quiet(f(TD = td, trait = trait, nPC = npc, excludeGeno = if (length(exclude)) exclude else NULL)),
                  error = function(e) structure(list(msg = conditionMessage(e)), class = "skipped"))
  res
}

scores_table <- function(fit, npc) {
  imp <- fit$importance
  ng <- nrow(fit$genoScores)
  d <- as.numeric(imp["Standard deviation", seq_len(npc)]) * sqrt(ng - 1)   # singular values
  gs <- sweep(fit$genoScores[, seq_len(npc), drop = FALSE], 2, sqrt(d), "/")   # symmetric scaling
  es <- sweep(fit$envScores[, seq_len(npc), drop = FALSE], 2, sqrt(d), "*")
  rbind(data.frame(type = "genotype", name = rownames(gs), gs, row.names = NULL, check.names = FALSE),
        data.frame(type = "environment", name = rownames(es), es, row.names = NULL, check.names = FALSE))
}

# Shukla's (1972) stability variance from Wricke's ecovalence W (p genotypes,
# q environments, complete two-way table of means). Equals agricolae's
# stability.par "Sigma-square" divided by its rep; its mean equals the G x E
# mean square.
shukla_variance <- function(W, q) {
  p <- length(W)
  if (p < 3 || q < 2) fail("Shukla's stability variance needs at least 3 genotypes and 2 environments.")
  p * W / ((p - 2) * (q - 1)) - sum(W) / ((p - 1) * (p - 2) * (q - 1))
}

run_stability <- function(td, d, trait, best, complete) {
  keep <- complete$complete_genotypes
  if (length(keep) < 3) return(structure(list(msg = "fewer than 3 genotypes are observed in every environment"), class = "skipped"))
  tdc <- make_td(d[d$genotype %in% keep, ], trait)
  st <- quiet(statgenGxE::gxeStability(TD = tdc, trait = trait, bestMethod = best, sorted = "none"))
  q <- length(unique(d$env))
  out <- data.frame(genotype = as.character(st$superiority$Genotype), mean = st$superiority$Mean,
                    superiority = st$superiority$Superiority, stringsAsFactors = FALSE)
  out$static_env_variance <- st$static$Static[match(out$genotype, as.character(st$static$Genotype))]
  out$wricke_ecovalence <- st$wricke$Wricke[match(out$genotype, as.character(st$wricke$Genotype))]
  out$shukla_variance <- shukla_variance(out$wricke_ecovalence, q)
  out$rank_mean <- rank(if (best == "max") -out$mean else out$mean, ties.method = "min")
  out$rank_superiority <- rank(out$superiority, ties.method = "min")
  out$rank_shukla <- rank(out$shukla_variance, ties.method = "min")
  out[order(out$rank_mean), ]
}

run_fw <- function(td, trait, exclude) {
  geno <- setdiff(unique(as.character(unlist(lapply(td, `[[`, "genotype")))), exclude)
  fw <- tryCatch(quiet(statgenGxE::gxeFw(TD = td, trait = trait, genotypes = if (length(exclude)) geno else NULL, sorted = "none")),
                 error = function(e) structure(list(msg = conditionMessage(e)), class = "skipped"))
  fw
}

run_megaenv <- function(td, trait, best, exclude, d) {
  tdx <- if (length(exclude)) make_td(d[!d$genotype %in% exclude, ], trait) else td
  tryCatch(quiet(statgenGxE::gxeMegaEnv(TD = tdx, trait = trait, method = best)),
           error = function(e) structure(list(msg = conditionMessage(e)), class = "skipped"))
}

analyze_met <- function(d, opt, weight_source = "none") {
  trait <- opt$trait
  best <- match.arg(tolower(opt$best), c("max", "min"))
  npc <- as.integer(opt$npc)
  if (npc < 1) fail("--npc must be at least 1.")
  comp <- completeness(d)
  ne <- ncol(comp$table); ng <- nrow(comp$table)
  td <- make_td(d, trait)
  res <- list(trait = trait, best = best, n_env = ne, n_geno = ng, completeness = comp, weight_source = weight_source)
  res$env_summary <- do.call(rbind, lapply(split(d, d$env), function(z)
    data.frame(environment = z$env[1], genotypes = nrow(z), mean = mean(z$y), sd = stats::sd(z$y),
               min = min(z$y), max = max(z$y),
               mean_se = if (all(is.na(z$se2))) NA else sqrt(mean(z$se2)), stringsAsFactors = FALSE)))
  res$varcomp <- run_varcomp(td, d, trait)
  if (weight_source != "none") {
    if (is.null(d$nesting) && is.null(d$location)) {
      res$weighted <- tryCatch(run_weighted_stage2(d), error = function(e) structure(list(msg = conditionMessage(e)), class = "skipped"))
    } else {
      res$weighted <- structure(list(msg = "the weighted stage-2 model is fitted only for the trial model (no --nesting or --location/--year)"), class = "skipped")
    }
  }
  # genotype table
  gm <- tapply(d$y, d$genotype, mean)
  bl <- res$varcomp$blup
  gt <- data.frame(genotype = names(gm), n_env = as.integer(comp$envs_per_genotype[names(gm)]),
                   raw_mean = as.numeric(gm), stringsAsFactors = FALSE)
  gt$blup <- bl$predictedValue[match(gt$genotype, as.character(bl$genotype))]
  if (!is.null(res$weighted) && !inherits(res$weighted, "skipped"))
    gt$blup_weighted <- res$weighted$blup$blup_weighted[match(gt$genotype, res$weighted$blup$genotype)]
  key <- if (!is.null(gt$blup_weighted)) gt$blup_weighted else gt$blup
  gt$rank <- rank(if (best == "max") -key else key, ties.method = "min")
  res$genotypes <- gt[order(gt$rank), ]
  exclude <- if (isTRUE(opt$complete_only)) setdiff(rownames(comp$table), comp$complete_genotypes) else character()
  res$excluded_incomplete <- exclude
  if (ne >= 3) {
    npc_use <- min(npc, ne - 1, ng - length(exclude) - 1)
    res$npc <- npc_use
    res$ammi <- ammi_like(td, trait, npc_use, FALSE, exclude)
    res$gge <- ammi_like(td, trait, npc_use, TRUE, exclude)
    res$stability <- run_stability(td, d, trait, best, comp)
    res$megaenv <- run_megaenv(td, trait, best, exclude, d)
  } else {
    sk <- structure(list(msg = "needs at least 3 environments"), class = "skipped")
    res$ammi <- res$gge <- res$stability <- res$megaenv <- sk
  }
  res$fw <- run_fw(td, trait, exclude)
  res
}

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

PAL <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e4e3df"; SURFACE <- "#fcfcfb"; MUTED <- "#a7a69f"

base_par <- function(...) graphics::par(bg = SURFACE, fg = INK2, col.axis = INK2, col.lab = INK, las = 1, mgp = c(2.3, 0.6, 0), tcl = -0.25, ...)

pc_label <- function(fit, k) sprintf("PC%d (%.1f%% of %s)", k, 100 * fit$importance["Proportion of Variance", k],
                                     if (isTRUE(attr(fit, "gge"))) "G + GxE" else "GxE")

biplot_base <- function(sc, fit, main, gge, best) {
  g <- sc[sc$type == "genotype", ]; e <- sc[sc$type == "environment", ]
  lim <- range(c(g$PC1, g$PC2, e$PC1, e$PC2)) * 1.4
  graphics::plot(NA, xlim = lim, ylim = lim, asp = 1, xlab = sprintf("PC1 (%.1f%% of %s)", 100 * fit$importance["Proportion of Variance", 1], if (gge) "G + GxE" else "GxE"),
                 ylab = sprintf("PC2 (%.1f%% of %s)", 100 * fit$importance["Proportion of Variance", 2], if (gge) "G + GxE" else "GxE"), main = "")
  graphics::abline(h = 0, v = 0, col = GRID)
  if (gge && nrow(g) >= 3) {
    h <- grDevices::chull(g$PC1, g$PC2)
    graphics::polygon(g$PC1[h], g$PC2[h], border = MUTED, lty = 1)
    hh <- c(h, h[1])
    for (k in seq_along(h)) {
      a <- c(g$PC1[hh[k]], g$PC2[hh[k]]); b <- c(g$PC1[hh[k + 1]], g$PC2[hh[k + 1]])
      dv <- b - a
      if (sum(dv^2) == 0) next
      foot <- a - sum(a * dv) / sum(dv^2) * dv           # foot of the perpendicular from the origin
      if (sqrt(sum(foot^2)) < 1e-12) next
      graphics::segments(0, 0, foot[1] * 50, foot[2] * 50, col = MUTED, lty = 2)
    }
    graphics::text(g$PC1[h], g$PC2[h], g$name[h], pos = 3, cex = 0.62, col = INK)
  }
  graphics::points(g$PC1, g$PC2, pch = 16, cex = 0.7, col = PAL[1])
  graphics::arrows(0, 0, e$PC1, e$PC2, length = 0.06, col = PAL[2], lwd = 1.6)
  graphics::text(e$PC1, e$PC2, e$name, pos = ifelse(e$PC1 >= 0, 4, 2), cex = 0.66, col = INK, xpd = TRUE)
  graphics::mtext(main, side = 3, line = 1.2, adj = 0, font = 2, cex = 0.9, col = INK)
  graphics::legend("bottomright", legend = c("genotype", "environment"), pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 1.6),
                   col = c(PAL[1], PAL[2]), bty = "n", cex = 0.7, text.col = INK2)
}

draw_figures <- function(res, out_dir) {
  made <- character()
  if (!inherits(res$ammi, "skipped") && res$npc >= 2) {
    open_png(out_dir, "ammi2_biplot.png", 7, 6.5); base_par(mar = c(4, 4, 3, 1))
    biplot_base(scores_table(res$ammi, 2), res$ammi, sprintf("AMMI2 biplot, %s (symmetric scaling)", res$trait), FALSE, res$best)
    grDevices::dev.off(); made <- c(made, "ammi2_biplot.png")
  }
  if (!inherits(res$gge, "skipped") && res$npc >= 2) {
    open_png(out_dir, "gge_biplot.png", 7, 6.5); base_par(mar = c(4, 4, 3, 1))
    biplot_base(scores_table(res$gge, 2), res$gge, sprintf("GGE biplot, which-won-where, %s", res$trait), TRUE, res$best)
    graphics::mtext(if (res$best == "max") "polygon vertices = highest-scoring genotypes for the environments in the sectors cut by the dashed rays"
                    else "lower is better: a vertex genotype is lowest for environments in the sector opposite its own",
                    side = 3, line = 0.3, adj = 0, cex = 0.62, col = INK2)
    grDevices::dev.off(); made <- c(made, "gge_biplot.png")
  }
  if (!inherits(res$fw, "skipped")) {
    est <- res$fw$estimates; env <- res$fw$envEffs
    open_png(out_dir, "finlay_wilkinson.png", 11, 5); base_par(mfrow = c(1, 2), mar = c(4, 5.6, 3, 1))
    xm <- env$EnvMean
    rng <- range(xm)
    xx <- seq(rng[1], rng[2], length.out = 50)
    eff <- (xx - mean(est$GenMean)) / mean(est$Sens)        # environment effect at each environment mean
    yl <- range(c(est$GenMean + outer(est$Sens, eff)), na.rm = TRUE)
    graphics::plot(NA, xlim = rng, ylim = yl, xlab = "environment mean", ylab = "", main = "")
    graphics::title(ylab = paste("fitted", res$trait), line = 4.2)
    graphics::grid(col = GRID, lty = 1)
    for (i in seq_len(nrow(est))) graphics::lines(xx, est$GenMean[i] + est$Sens[i] * eff, col = "#d6d5d0", lwd = 0.8)
    hi <- c(which.max(est$GenMean), which.max(est$Sens), which.min(est$Sens))
    lab <- c("highest mean", "most responsive", "least responsive")
    for (k in seq_along(hi)) graphics::lines(xx, est$GenMean[hi[k]] + est$Sens[hi[k]] * eff, col = PAL[k], lwd = 2)
    graphics::rug(xm, col = INK2)
    graphics::legend("topleft", legend = paste0(lab, ": ", est$Genotype[hi]), col = PAL[1:3], lwd = 2, bty = "n", cex = 0.7, text.col = INK2)
    graphics::mtext("Finlay-Wilkinson regressions (gray = all genotypes)", side = 3, line = 1, adj = 0, font = 2, cex = 0.85, col = INK)
    graphics::plot(est$GenMean, est$Sens, pch = 16, col = PAL[1], cex = 0.8, xlab = paste("genotype mean,", res$trait), ylab = "sensitivity (slope)", main = "")
    graphics::abline(h = 1, v = mean(est$GenMean), col = MUTED, lty = 2)
    graphics::mtext("sensitivity vs mean", side = 3, line = 1.6, adj = 0, font = 2, cex = 0.85, col = INK)
    graphics::mtext("slope > 1 = responds more to good environments", side = 3, line = 0.6, adj = 0, cex = 0.7, col = INK2)
    grDevices::dev.off(); made <- c(made, "finlay_wilkinson.png")
  }
  if (!inherits(res$stability, "skipped")) {
    st <- res$stability
    open_png(out_dir, "stability_vs_mean.png", 11, 5); base_par(mfrow = c(1, 2), mar = c(4, 6, 3, 1))
    graphics::plot(st$mean, st$shukla_variance, pch = 16, col = PAL[1], cex = 0.8, xlab = paste("genotype mean,", res$trait),
                   ylab = "", main = "")
    graphics::title(ylab = "Shukla stability variance", line = 4.6)
    graphics::abline(v = mean(st$mean), col = MUTED, lty = 2)
    top <- utils::head(order(st$rank_mean), 5)
    graphics::text(st$mean[top], st$shukla_variance[top], st$genotype[top], pos = 3, cex = 0.6, col = INK)
    graphics::mtext("stability is not performance: low variance = stable", side = 3, line = 1, adj = 0, font = 2, cex = 0.85, col = INK)
    graphics::plot(st$mean, st$superiority, pch = 16, col = PAL[3], cex = 0.8, xlab = paste("genotype mean,", res$trait),
                   ylab = "", main = "")
    graphics::title(ylab = "Lin-Binns superiority (low = near the best)", line = 4.6)
    graphics::text(st$mean[top], st$superiority[top], st$genotype[top], pos = 3, cex = 0.6, col = INK)
    graphics::mtext("superiority combines mean and consistency", side = 3, line = 1, adj = 0, font = 2, cex = 0.85, col = INK)
    grDevices::dev.off(); made <- c(made, "stability_vs_mean.png")
  }
  made
}

# ---------------------------------------------------------------------------
# Tables, warnings and report
# ---------------------------------------------------------------------------

fmt <- function(x, d = 3) ifelse(is.na(x), "", formatC(x, digits = d, format = "g"))

md_table <- function(df, digits = 4) {
  df[] <- lapply(df, function(x) if (is.numeric(x)) fmt(x, digits) else { x <- as.character(x); x[is.na(x)] <- ""; x })
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

anova_df <- function(a) {
  a <- as.data.frame(a)
  data.frame(term = rownames(a), df = a$Df, sum_sq = a[["Sum Sq"]], mean_sq = a[["Mean Sq"]],
             f_value = a[["F value"]], p_value = a[["Pr(>F)"]], stringsAsFactors = FALSE, row.names = NULL)
}

met_warnings <- function(res, notes) {
  w <- notes
  comp <- res$completeness
  if (res$n_env < 5)
    w <- c(w, sprintf("Only %d environments: G x E variance, stability and biplot patterns are weakly estimated (rule of thumb: at least about 5 environments). Treat them as descriptive.", res$n_env))
  if (comp$missing > 0)
    w <- c(w, if (length(res$excluded_incomplete))
      sprintf("%d of %d genotype x environment cells (%.1f%%) are missing. --complete-only dropped the %d incomplete genotypes from AMMI, GGE, Finlay-Wilkinson and mega-environments; stability also uses only the %d complete genotypes. Variance components and BLUPs use all data.",
              comp$missing, comp$cells, 100 * comp$missing / comp$cells, length(res$excluded_incomplete), length(comp$complete_genotypes))
    else
      sprintf("%d of %d genotype x environment cells (%.1f%%) are missing. AMMI, GGE and mega-environments fill them by iterative regression (statgenGxE's multMissing) before the PCA; --complete-only drops incomplete genotypes instead. Stability uses only the %d genotypes seen in every environment.",
              comp$missing, comp$cells, 100 * comp$missing / comp$cells, length(comp$complete_genotypes)))
  one <- names(comp$envs_per_genotype)[comp$envs_per_genotype == 1]
  if (length(one)) w <- c(w, sprintf("%d genotype(s) were tested in only one environment (%s): their BLUPs, sensitivities and scores rest on one trial.", length(one), paste(utils::head(one, 5), collapse = ", ")))
  for (nm in c("ammi", "gge", "fw", "stability", "megaenv"))
    if (inherits(res[[nm]], "skipped")) w <- c(w, sprintf("%s skipped: %s", c(ammi = "AMMI", gge = "GGE", fw = "Finlay-Wilkinson", stability = "Stability", megaenv = "Mega-environments")[[nm]], trimws(res[[nm]]$msg)))
  if (!is.null(res$weighted) && inherits(res$weighted, "skipped")) w <- c(w, paste("Weighted stage-2 model skipped:", res$weighted$msg))
  vw <- res$varcomp$warnings
  if (length(vw)) w <- c(w, paste("gxeVarComp:", vw))
  noise <- "essentially perfect fit"
  for (nm in c("ammi", "gge", "fw", "stability", "megaenv")) {
    pw <- attr(res[[nm]], "warnings")
    pw <- pw[!grepl(noise, pw, fixed = TRUE)]
    if (length(pw)) w <- c(w, paste0(c(ammi = "gxeAmmi", gge = "gxeGGE", fw = "gxeFw", stability = "gxeStability", megaenv = "gxeMegaEnv")[[nm]], ": ", pw))
  }
  vc <- res$varcomp$components
  if (any(vc$variance < 1e-8 * max(vc$variance)))
    w <- c(w, sprintf("Variance component(s) estimated at zero: %s.", paste(vc$component[vc$variance < 1e-8 * max(vc$variance)], collapse = ", ")))
  if (!inherits(res$megaenv, "skipped") && res$n_env < 8)
    w <- c(w, "Mega-environments from a single set of trials describe this data set only; they are mega-environments only if the grouping repeats across years (Yan et al. 2000).")
  w
}

REFERENCES <- c(
  "- van Rossum B (2025). statgenGxE: Genotype by Environment (GxE) Analysis. R package version 1.0.11. https://doi.org/10.32614/CRAN.package.statgenGxE",
  "- Pinheiro JC, Bates DM (2000). Mixed-Effects Models in S and S-PLUS. Springer, New York (R package nlme). https://doi.org/10.1007/b98882",
  "- Smith AB, Cullis BR, Thompson R (2005). The analysis of crop cultivar breeding and evaluation trials: an overview of current mixed model approaches. Journal of Agricultural Science 143(6):449-462. https://doi.org/10.1017/S0021859605005587",
  "- Finlay KW, Wilkinson GN (1963). The analysis of adaptation in a plant-breeding programme. Australian Journal of Agricultural Research 14(6):742-754. https://doi.org/10.1071/AR9630742",
  "- Gauch HG (1988). Model selection and validation for yield trials with interaction. Biometrics 44(3):705-715. https://doi.org/10.2307/2531585",
  "- Gauch HG (2006). Statistical analysis of yield trials by AMMI and GGE. Crop Science 46(4):1488-1500. https://doi.org/10.2135/cropsci2005.07-0193",
  "- Yan W, Hunt LA, Sheng Q, Szlavnics Z (2000). Cultivar evaluation and mega-environment investigation based on the GGE biplot. Crop Science 40(3):597-605. https://doi.org/10.2135/cropsci2000.403597x",
  "- Yan W, Tinker NA (2006). Biplot analysis of multi-environment trial data: Principles and applications. Canadian Journal of Plant Science 86(3):623-645. https://doi.org/10.4141/P05-169",
  "- Shukla GK (1972). Some statistical aspects of partitioning genotype-environmental components of variability. Heredity 29(2):237-245. https://doi.org/10.1038/hdy.1972.87",
  "- Lin CS, Binns MR (1988). A superiority measure of cultivar performance for cultivar x location data. Canadian Journal of Plant Science 68(1):193-198. https://doi.org/10.4141/cjps88-018",
  "- Becker HC, Leon J (1988). Stability analysis in plant breeding. Plant Breeding 101(1):1-23. https://doi.org/10.1111/j.1439-0523.1988.tb00261.x",
  "- Atlin GN, Baker RJ, McRae KB, Lu X (2000). Selection response in subdivided target regions. Crop Science 40(1):7-13. https://doi.org/10.2135/cropsci2000.4017",
  "- Holland JB, Nyquist WE, Cervantes-Martinez CT (2003). Estimating and interpreting heritability for plant breeding: an update. Plant Breeding Reviews 22:9-112. https://doi.org/10.1002/9780470650202.ch2",
  "- Piepho H-P, Mohring J, Schulz-Streeck T, Ogutu JO (2012). A stage-wise approach for the analysis of multi-environment trials. Biometrical Journal 54(6):844-860. https://doi.org/10.1002/bimj.201100219",
  "- Bates D, Machler M, Bolker B, Walker S (2015). Fitting linear mixed-effects models using lme4. Journal of Statistical Software 67(1):1-48. https://doi.org/10.18637/jss.v067.i01")

write_report <- function(res, opt, out_dir, warn, figs, inputs) {
  vc <- res$varcomp
  L <- c(sprintf("# %s report", SKILL), "")
  if (length(warn)) L <- c(L, "## Warnings", "", paste("-", warn), "")
  comp <- res$completeness
  L <- c(L, "## Data", "",
         sprintf("Input: %s. Trait: `%s`. %d environments x %d genotypes; %d of %d cells observed (%.1f%% missing); %d genotypes in every environment.",
                 paste(basename(inputs), collapse = ", "), res$trait, res$n_env, res$n_geno, comp$cells - comp$missing, comp$cells,
                 100 * comp$missing / comp$cells, length(comp$complete_genotypes)),
         sprintf("Weights: %s. Better direction: %s.", if (res$weight_source == "none") "none" else paste("stage-1 variances from", res$weight_source), res$best), "",
         md_table(res$env_summary), "",
         "## Variance components and heritability (statgenGxE::gxeVarComp, lme4, unweighted)", "",
         sprintf("Model: %s.", vc$model), "",
         md_table(vc$components), "",
         sprintf("Heritability across environments (statgenGxE::herit, entry-mean form of Atlin et al. 2000 with harmonic means for unbalanced data, Holland et al. 2003): **H2 = %.3f**.", vc$h2), "",
         "Share of total variance in the fully random model (environment as random):", "", md_table(vc$full_random), "")
  if (!is.null(res$weighted) && !inherits(res$weighted, "skipped")) {
    w <- res$weighted
    L <- c(L, "## Weighted stage-2 model (known stage-1 variances; Piepho et al. 2012)", "",
           sprintf("nlme::lme with residual variance of each BLUE fixed at its SE^2: genotype variance %.4g, G x E variance %.4g, mean stage-1 error variance %.4g. Entry-mean H2 = %.3f (same form as above, with G x E and stage-1 error separated; %.2f environments per genotype, harmonic mean).",
                   w$sigma2_g, w$sigma2_ge, w$mean_stage1_var, w$h2, w$env_harmonic_mean),
           sprintf("Correlation of weighted and unweighted genotype BLUPs: %.3f.", stats::cor(res$genotypes$blup, res$genotypes$blup_weighted, use = "complete.obs")), "")
  }
  L <- c(L, "## Genotype predictions across environments", "",
         sprintf("`tables/genotype_blups.csv` (top 15 shown; ranked by %s BLUP):", if (is.null(res$genotypes$blup_weighted)) "the" else "the weighted"), "",
         md_table(utils::head(res$genotypes, 15)), "")
  if (!inherits(res$ammi, "skipped")) {
    imp <- res$ammi$importance
    L <- c(L, "## AMMI", "",
           sprintf("PC1 explains %.1f%% and PC2 %.1f%% of the G x E sum of squares (statgenGxE::gxeAmmi, %d PCs).",
                   100 * imp["Proportion of Variance", 1], if (ncol(imp) > 1) 100 * imp["Proportion of Variance", 2] else NA, res$npc), "",
           md_table(anova_df(res$ammi$anova)), "")
  }
  if (!inherits(res$gge, "skipped")) {
    imp <- res$gge$importance
    L <- c(L, "## GGE", "", sprintf("PC1 explains %.1f%% and PC2 %.1f%% of G + G x E (%.1f%% together). The which-won-where reading of the biplot is only as good as this share.",
                                    100 * imp["Proportion of Variance", 1], 100 * imp["Proportion of Variance", 2],
                                    100 * imp["Cumulative Proportion", min(2, ncol(imp))]), "")
  }
  if (!inherits(res$fw, "skipped")) {
    est <- res$fw$estimates
    L <- c(L, "## Finlay-Wilkinson", "",
           sprintf("Sensitivities range %.2f to %.2f (mean %.2f; 1 = average response). Converged in %d iterations.",
                   min(est$Sens), max(est$Sens), mean(est$Sens), res$fw$iter), "",
           md_table(anova_df(res$fw$anova)), "")
  }
  if (!inherits(res$stability, "skipped")) {
    L <- c(L, "## Stability (genotypes observed in every environment)", "",
           "Superiority (Lin & Binns 1988): distance to the best genotype in each environment, low = good. Wricke's ecovalence and Shukla's (1972) stability variance: contribution to G x E, low = stable. Static stability (environmental variance; Becker & Leon 1988; labeled 'static' in statgenGxE): low = little change across environments, usually low-yielding genotypes.", "",
           md_table(utils::head(res$stability, 15)), "")
  }
  if (!inherits(res$megaenv, "skipped")) {
    L <- c(L, "## Mega-environments (statgenGxE::gxeMegaEnv: trials grouped by their AMMI2 winner)", "",
           md_table(res$megaenv$summTab), "")
  }
  if (length(figs)) L <- c(L, "## Figures", "", sprintf("![%s](figures/%s)", sub("\\.png$", "", figs), figs), "")
  L <- c(L, "## Files", "",
         "- `tables/genotype_blups.csv`: genotype, environments, raw mean, MET BLUP (unweighted), weighted BLUP (when SEs were given), rank.",
         "- `tables/variance_components.csv`, `tables/heritability.csv`, `tables/environment_summary.csv`.",
         "- `tables/ammi_anova.csv`, `tables/pc_importance.csv`, `tables/ammi_scores.csv`, `tables/gge_scores.csv` (symmetric scaling).",
         "- `tables/finlay_wilkinson.csv`, `tables/fw_anova.csv`, `tables/fw_environments.csv`, `tables/stability.csv`, `tables/mega_environments.csv`.", "",
         "## References", "", REFERENCES, "", sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(L, file.path(out_dir, "report.md"))
}

save_outputs <- function(res, out_dir) {
  save_table(res$env_summary, out_dir, "environment_summary.csv")
  vc <- res$varcomp
  save_table(rbind(data.frame(model = "gxeVarComp (lme4, unweighted)", component = vc$components$component, variance = vc$components$variance),
                   if (!is.null(res$weighted) && !inherits(res$weighted, "skipped"))
                     data.frame(model = "stage-2 with known stage-1 variances (nlme)", component = c("genotype", "genotype:trial (GxE)", "stage-1 error (mean SE^2, fixed)"),
                                variance = c(res$weighted$sigma2_g, res$weighted$sigma2_ge, res$weighted$mean_stage1_var))),
             out_dir, "variance_components.csv")
  save_table(data.frame(model = c("gxeVarComp + herit (unweighted)", if (!is.null(res$weighted) && !inherits(res$weighted, "skipped")) "stage-2 with known stage-1 variances"),
                        h2 = c(vc$h2, if (!is.null(res$weighted) && !inherits(res$weighted, "skipped")) res$weighted$h2),
                        definition = "sigma2_G / (sigma2_G + sum(sigma2_term / n_term) + sigma2_residual / (r n)), harmonic means (Atlin et al. 2000; Holland et al. 2003)"),
             out_dir, "heritability.csv")
  save_table(res$genotypes, out_dir, "genotype_blups.csv")
  imp <- NULL
  if (!inherits(res$ammi, "skipped")) {
    save_table(anova_df(res$ammi$anova), out_dir, "ammi_anova.csv")
    save_table(scores_table(res$ammi, res$npc), out_dir, "ammi_scores.csv")
    imp <- rbind(imp, data.frame(model = "AMMI", pc = colnames(res$ammi$importance), t(res$ammi$importance), check.names = FALSE, row.names = NULL))
  }
  if (!inherits(res$gge, "skipped")) {
    save_table(scores_table(res$gge, res$npc), out_dir, "gge_scores.csv")
    imp <- rbind(imp, data.frame(model = "GGE", pc = colnames(res$gge$importance), t(res$gge$importance), check.names = FALSE, row.names = NULL))
  }
  if (!is.null(imp)) save_table(imp, out_dir, "pc_importance.csv")
  if (!inherits(res$fw, "skipped")) {
    est <- res$fw$estimates
    save_table(data.frame(genotype = est$Genotype, mean = est$GenMean, se_mean = est$SE_GenMean, sensitivity = est$Sens,
                          se_sensitivity = est$SE_Sens, ms_deviation = est$MSdeviation), out_dir, "finlay_wilkinson.csv")
    save_table(anova_df(res$fw$anova), out_dir, "fw_anova.csv")
    save_table(res$fw$envEffs, out_dir, "fw_environments.csv")
  }
  if (!inherits(res$stability, "skipped")) save_table(res$stability, out_dir, "stability.csv")
  if (!inherits(res$megaenv, "skipped")) save_table(res$megaenv$summTab, out_dir, "mega_environments.csv")
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (opt$demo) {
    opt$input <- file.path(SKILL_DIR, "demo", "met_input.csv")
    if (!nzchar(opt$trait)) opt$trait <- "yield_kg_ha"
  }
  if (!nzchar(opt$trait)) fail("Give --trait COL (the per-environment BLUE column), or --demo. See --help.")
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  rd <- read_inputs(opt$input)
  prep <- prepare_met(rd$data, opt)
  res <- analyze_met(prep$dat, opt, prep$weight_source)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  save_outputs(res, opt$out)
  figs <- draw_figures(res, opt$out)
  warn <- met_warnings(res, prep$notes)
  write_report(res, opt, opt$out, warn, figs, rd$paths)
  write_repro_bundle(opt$out, argv, inputs = rd$paths)
  cat(sprintf("Done: %d environments x %d genotypes, H2 = %.2f. Report: %s\n", res$n_env, res$n_geno, res$varcomp$h2,
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  if (length(warn)) cat(paste("Warning:", warn), sep = "\n")
  invisible(res)
}

if (sys.nframe() == 0L) main()
