#!/usr/bin/env Rscript
# claw-growth-analysis: classical and functional plant growth analysis from
# sequential destructive harvests. Classical: relative growth rate (RGR) from
# the mean of ln(mass) per harvest (Hoffmann & Poorter 2002) with standard
# errors, net assimilation rate (NAR, unit leaf rate), leaf area ratio (LAR),
# specific leaf area (SLA) and leaf mass fraction (LMF). Functional: a
# penalized spline of ln(mass) on time (mgcv) giving RGR(t) with confidence
# bands.
#
# Usage:
#   Rscript scripts/growth_analysis.R --demo
#   Rscript scripts/growth_analysis.R --input harvests.csv --group genotype --time-col days
#   Rscript scripts/growth_analysis.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT license.

SKILL   <- "claw-growth-analysis"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/growth_analysis.R [--demo] [--input FILE.csv] [options]

  --input FILE        CSV, one row per harvested plant
  --group COL         treatment/genotype column (default: genotype; '' = one group)
  --time-col COL      harvest time in days, or ISO dates (default: days)
  --mass-col COL      total plant dry mass, g (default: total_dry_mass_g)
  --leaf-area-col COL leaf area, cm2; optional (default: leaf_area_cm2)
  --leaf-mass-col COL leaf dry mass, g; optional (default: leaf_dry_mass_g)
  --k-max N           maximum spline basis size for the functional fit (default: 10)
  --n-boot N          bootstrap resamples for NAR confidence limits (default: 2000)
  --min-n N           warn when a harvest has fewer plants than this (default: 4)
  --seed N            random seed for the bootstrap (default: 1)
  --out DIR           output folder (default: agriclaw_output/claw-growth-analysis)
  --demo              run on the bundled synthetic harvest data
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, input = "", out = "", group = "genotype",
                 time_col = "days", mass_col = "total_dry_mass_g",
                 leaf_area_col = "leaf_area_cm2", leaf_mass_col = "leaf_dry_mass_g",
                 k_max = 10, n_boot = 2000, min_n = 4, seed = 1)

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

read_input <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
}

# One row per harvested plant -> data frame with group, time, mass, la, lm.
prepare_data <- function(raw, group = "genotype", time_col = "days",
                         mass_col = "total_dry_mass_g", leaf_area_col = "leaf_area_cm2",
                         leaf_mass_col = "leaf_dry_mass_g", warn = function(m) NULL) {
  need <- c(if (nzchar(group)) group, time_col, mass_col)
  miss <- setdiff(need, names(raw))
  if (length(miss))
    stop("Column(s) not found: ", paste(miss, collapse = ", "), ". Available: ",
         paste(names(raw), collapse = ", "), call. = FALSE)
  tm <- raw[[time_col]]
  if (!is.numeric(tm)) {
    dt <- suppressWarnings(as.Date(as.character(tm), tryFormats = c("%Y-%m-%d", "%Y/%m/%d")))
    if (anyNA(dt))
      stop("'", time_col, "' must be numeric days or ISO dates (YYYY-MM-DD).", call. = FALSE)
    tm <- as.numeric(dt - min(dt))
    warn(sprintf("'%s' holds dates: converted to days after the first harvest (%s).", time_col, format(min(dt))))
  }
  if (anyNA(tm)) stop("'", time_col, "' has missing values.", call. = FALSE)
  w <- raw[[mass_col]]
  if (!is.numeric(w)) stop("'", mass_col, "' must be numeric (g).", call. = FALSE)
  if (anyNA(w))
    stop("'", mass_col, "' has ", sum(is.na(w)), " missing value(s). Remove those plants ",
         "first; nothing is filled in.", call. = FALSE)
  if (any(w <= 0))
    stop("'", mass_col, "' must be > 0 (ln mass is used); check rows ",
         paste(utils::head(which(w <= 0), 5), collapse = ", "), ".", call. = FALSE)
  opt_col <- function(cl, what) {
    if (!nzchar(cl) || !cl %in% names(raw)) return(rep(NA_real_, nrow(raw)))
    x <- raw[[cl]]
    if (!is.numeric(x)) stop("'", cl, "' must be numeric (", what, ").", call. = FALSE)
    if (any(x <= 0, na.rm = TRUE)) stop("'", cl, "' must be > 0 where present.", call. = FALSE)
    if (anyNA(x)) warn(sprintf("'%s' is missing for %d plant(s): ratios using it are based on the other plants.", cl, sum(is.na(x))))
    x
  }
  la <- opt_col(leaf_area_col, "cm2")
  lm <- opt_col(leaf_mass_col, "g")
  if (any(lm > w, na.rm = TRUE))
    stop("Leaf dry mass exceeds total dry mass for ", sum(lm > w, na.rm = TRUE),
         " plant(s): is the mass column total plant mass (g)?", call. = FALSE)
  d <- data.frame(group = if (nzchar(group)) as.character(raw[[group]]) else "all",
                  time = as.numeric(tm), mass = w, la = la, lm = lm,
                  stringsAsFactors = FALSE)
  nt <- tapply(d$time, d$group, function(x) length(unique(x)))
  if (any(nt < 2))
    stop("Every group needs at least two harvest times; only one in: ",
         paste(names(nt)[nt < 2], collapse = ", "), call. = FALSE)
  d[order(d$group, d$time), ]
}

se_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) > 1) stats::sd(x) / sqrt(length(x)) else NA_real_
}

# Per group and harvest: means of mass and ln mass, and per-plant ratios
#   SLA = leaf area / leaf mass, LAR = leaf area / total mass (both m2 kg-1;
#   1 cm2 g-1 = 0.1 m2 kg-1), LMF = leaf mass / total mass (g g-1),
# averaged over plants (mean of per-plant ratios, with SE).
harvest_summary <- function(d) {
  keys <- unique(d[, c("group", "time")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    s <- d[d$group == keys$group[i] & d$time == keys$time[i], ]
    sla <- 0.1 * s$la / s$lm
    lar <- 0.1 * s$la / s$mass
    lmf <- s$lm / s$mass
    data.frame(group = keys$group[i], time = keys$time[i], n = nrow(s),
               mass_mean_g = mean(s$mass), mass_sd_g = if (nrow(s) > 1) stats::sd(s$mass) else NA_real_,
               ln_mass_mean = mean(log(s$mass)),
               ln_mass_sd = if (nrow(s) > 1) stats::sd(log(s$mass)) else NA_real_,
               mass_geomean_g = exp(mean(log(s$mass))),
               leaf_area_mean_cm2 = if (any(is.finite(s$la))) mean(s$la, na.rm = TRUE) else NA_real_,
               sla_m2_kg = if (any(is.finite(sla))) mean(sla, na.rm = TRUE) else NA_real_,
               sla_se = se_mean(sla),
               lar_m2_kg = if (any(is.finite(lar))) mean(lar, na.rm = TRUE) else NA_real_,
               lar_se = se_mean(lar),
               lmf_g_g = if (any(is.finite(lmf))) mean(lmf, na.rm = TRUE) else NA_real_,
               lmf_se = se_mean(lmf), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Mean NAR (unit leaf rate) over a harvest interval, Williams (1946) and
# Radford (1967): NAR = (W2 - W1)/(t2 - t1) x (ln LA2 - ln LA1)/(LA2 - LA1),
# valid when total mass is linearly related to leaf area over the interval.
# Inputs g, cm2 and days; returns g m-2 d-1.
nar_interval <- function(w1, w2, la1, la2, dt) {
  out <- ifelse(abs(la2 - la1) <= 1e-12 * pmax(la1, la2),
                (w2 - w1) / (dt * la1),
                (w2 - w1) / dt * (log(la2) - log(la1)) / (la2 - la1))
  out * 1e4
}

# Classical analysis between consecutive harvests (plants are not paired).
#   RGR = (mean ln W2 - mean ln W1) / (t2 - t1)        Hoffmann & Poorter (2002)
#   SE  = sqrt(s1^2/n1 + s2^2/n2) / (t2 - t1), Welch-Satterthwaite df for the CI
# rgr_ln_of_means_g_g_d is the biased ln(mean W) version, kept only for comparison.
# NAR confidence limits: percentile bootstrap, resampling plants within each
# harvest independently.
classical_intervals <- function(d, n_boot = 2000, conf = 0.95) {
  rows <- list()
  for (g in sort(unique(d$group))) {
    dg <- d[d$group == g, ]
    times <- sort(unique(dg$time))
    for (i in seq_len(length(times) - 1L)) {
      a <- dg[dg$time == times[i], ]
      b <- dg[dg$time == times[i + 1L], ]
      dt <- times[i + 1L] - times[i]
      la <- log(a$mass); lb <- log(b$mass)
      n1 <- length(la); n2 <- length(lb)
      rgr <- (mean(lb) - mean(la)) / dt
      se <- df <- lo <- hi <- NA_real_
      if (n1 > 1 && n2 > 1) {
        v1 <- stats::var(la) / n1; v2 <- stats::var(lb) / n2
        se <- sqrt(v1 + v2) / dt
        df <- if (v1 + v2 > 0) (v1 + v2)^2 / (v1^2 / (n1 - 1) + v2^2 / (n2 - 1)) else n1 + n2 - 2
        q <- stats::qt(1 - (1 - conf) / 2, df)
        lo <- rgr - q * se; hi <- rgr + q * se
      }
      nar <- nar_se <- nar_lo <- nar_hi <- NA_real_
      ak <- a[is.finite(a$la), ]; bk <- b[is.finite(b$la), ]
      if (nrow(ak) > 0 && nrow(bk) > 0) {
        nar <- nar_interval(mean(ak$mass), mean(bk$mass), mean(ak$la), mean(bk$la), dt)
        if (nrow(ak) > 1 && nrow(bk) > 1 && n_boot > 0) {
          ia <- matrix(sample.int(nrow(ak), nrow(ak) * n_boot, replace = TRUE), nrow = n_boot)
          ib <- matrix(sample.int(nrow(bk), nrow(bk) * n_boot, replace = TRUE), nrow = n_boot)
          bw1 <- rowMeans(matrix(ak$mass[ia], nrow = n_boot)); bl1 <- rowMeans(matrix(ak$la[ia], nrow = n_boot))
          bw2 <- rowMeans(matrix(bk$mass[ib], nrow = n_boot)); bl2 <- rowMeans(matrix(bk$la[ib], nrow = n_boot))
          bn <- nar_interval(bw1, bw2, bl1, bl2, dt)
          nar_se <- stats::sd(bn)
          qq <- stats::quantile(bn, c((1 - conf) / 2, 1 - (1 - conf) / 2), names = FALSE)
          nar_lo <- qq[1]; nar_hi <- qq[2]
        }
      }
      rows[[length(rows) + 1L]] <- data.frame(
        group = g, t1 = times[i], t2 = times[i + 1L], n1 = n1, n2 = n2,
        rgr_g_g_d = rgr, rgr_se = se, rgr_df = df, rgr_lo = lo, rgr_hi = hi,
        rgr_ln_of_means_g_g_d = (log(mean(b$mass)) - log(mean(a$mass))) / dt,
        nar_g_m2_d = nar, nar_se = nar_se, nar_lo = nar_lo, nar_hi = nar_hi,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

# Mean RGR over the whole experiment: slope of ln W on time per group (valid as
# an average; a significant quadratic term means RGR changed with time).
overall_rgr <- function(d, conf = 0.95) {
  rows <- lapply(sort(unique(d$group)), function(g) {
    dg <- d[d$group == g, ]
    f <- stats::lm(log(mass) ~ time, data = dg)
    co <- summary(f)$coefficients
    ci <- stats::confint(f, "time", level = conf)
    p_quad <- NA_real_
    if (length(unique(dg$time)) >= 3) {
      f2 <- stats::lm(log(mass) ~ time + I(time^2), data = dg)
      p_quad <- summary(f2)$coefficients["I(time^2)", 4]
    }
    data.frame(group = g, n = nrow(dg), n_harvests = length(unique(dg$time)),
               rgr_g_g_d = co["time", 1], rgr_se = co["time", 2], rgr_lo = ci[1], rgr_hi = ci[2],
               r2 = summary(f)$r.squared, p_curvature = p_quad, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Functional approach: penalized cubic regression spline of ln W on time
# (mgcv, REML smoothing; the penalty shrinks toward a straight line, i.e.
# constant RGR). RGR(t) = d ln W / dt from the derivative of the basis; its SE
# uses the smoothing-parameter-corrected covariance (Wood 2017).
functional_rgr <- function(d, k_max = 10, n_grid = 101, conf = 0.95) {
  z <- stats::qnorm(1 - (1 - conf) / 2)
  out <- list(); info <- list()
  for (g in sort(unique(d$group))) {
    dg <- data.frame(lnw = log(d$mass[d$group == g]), time = d$time[d$group == g])
    nt <- length(unique(dg$time))
    if (nt < 3) {
      info[[g]] <- data.frame(group = g, k = NA, edf = NA, note = "fewer than 3 harvest times: no functional fit",
                              stringsAsFactors = FALSE)
      next
    }
    k <- as.integer(min(k_max, nt))
    fit <- mgcv::gam(lnw ~ s(time, k = k, bs = "cr"), data = dg, method = "REML")
    V <- tryCatch(stats::vcov(fit, unconditional = TRUE), error = function(e) fit$Vp)
    grid <- sort(unique(c(seq(min(dg$time), max(dg$time), length.out = n_grid), unique(dg$time))))
    eps <- 1e-4 * diff(range(dg$time))
    X <- stats::predict(fit, newdata = data.frame(time = grid), type = "lpmatrix")
    D <- (stats::predict(fit, newdata = data.frame(time = grid + eps), type = "lpmatrix") -
            stats::predict(fit, newdata = data.frame(time = grid - eps), type = "lpmatrix")) / (2 * eps)
    b <- stats::coef(fit)
    fit_ln <- drop(X %*% b); se_ln <- sqrt(rowSums((X %*% V) * X))
    rgr <- drop(D %*% b); se_r <- sqrt(rowSums((D %*% V) * D))
    out[[g]] <- data.frame(group = g, time = grid, at_harvest = grid %in% unique(dg$time),
                           ln_mass_fit = fit_ln, ln_mass_se = se_ln, mass_fit_g = exp(fit_ln),
                           rgr_g_g_d = rgr, rgr_se = se_r, rgr_lo = rgr - z * se_r, rgr_hi = rgr + z * se_r,
                           stringsAsFactors = FALSE)
    info[[g]] <- data.frame(group = g, k = k, edf = sum(fit$edf[-1]), note = "", stringsAsFactors = FALSE)
  }
  list(curve = if (length(out)) do.call(rbind, out) else NULL, info = do.call(rbind, info))
}

# --- Figures ---------------------------------------------------------------------

group_cols <- function(groups) {
  cols <- if (length(groups) <= 1) "#1f78b4" else grDevices::hcl.colors(length(groups), "Dark 3")
  stats::setNames(cols, groups)
}

plot_growth <- function(d, cls, fun, out_dir) {
  groups <- sort(unique(d$group)); cols <- group_cols(groups)
  files <- character()
  open_png(out_dir, "ln_mass.png")
  graphics::plot(d$time, log(d$mass), type = "n", xlab = "Time (d)", ylab = "ln(total dry mass, g)")
  for (g in groups) {
    s <- d[d$group == g, ]
    graphics::points(jitter(s$time, amount = 0.15), log(s$mass), pch = 16, cex = 0.7,
                     col = grDevices::adjustcolor(cols[g], 0.6))
    if (!is.null(fun$curve) && g %in% fun$curve$group) {
      cg <- fun$curve[fun$curve$group == g, ]
      graphics::polygon(c(cg$time, rev(cg$time)), c(cg$ln_mass_fit - 1.96 * cg$ln_mass_se,
                        rev(cg$ln_mass_fit + 1.96 * cg$ln_mass_se)), border = NA,
                        col = grDevices::adjustcolor(cols[g], 0.15))
      graphics::lines(cg$time, cg$ln_mass_fit, col = cols[g], lwd = 2)
    }
  }
  if (length(groups) > 1) graphics::legend("topleft", legend = groups, col = cols, pch = 16, lwd = 2, bty = "n", cex = 0.8)
  grDevices::dev.off(); files <- c(files, "ln_mass.png")

  open_png(out_dir, "rgr.png")
  yl <- range(c(cls$rgr_lo, cls$rgr_hi, cls$rgr_g_g_d, fun$curve$rgr_lo, fun$curve$rgr_hi), na.rm = TRUE)
  graphics::plot(range(d$time), yl, type = "n", xlab = "Time (d)",
                 ylab = expression(RGR ~ (g ~ g^-1 ~ d^-1)))
  graphics::abline(h = 0, col = "grey80")
  for (g in groups) {
    if (!is.null(fun$curve) && g %in% fun$curve$group) {
      cg <- fun$curve[fun$curve$group == g, ]
      graphics::polygon(c(cg$time, rev(cg$time)), c(cg$rgr_lo, rev(cg$rgr_hi)), border = NA,
                        col = grDevices::adjustcolor(cols[g], 0.15))
      graphics::lines(cg$time, cg$rgr_g_g_d, col = cols[g], lwd = 2)
    }
    cc <- cls[cls$group == g, ]
    dodge <- (match(g, groups) - (length(groups) + 1) / 2) * 0.012 * diff(range(d$time))
    mid <- (cc$t1 + cc$t2) / 2 + dodge
    graphics::segments(cc$t1, cc$rgr_g_g_d, cc$t2, cc$rgr_g_g_d, col = cols[g], lwd = 1.5, lty = 2)
    ok <- is.finite(cc$rgr_lo)
    if (any(ok)) graphics::arrows(mid[ok], cc$rgr_lo[ok], mid[ok], cc$rgr_hi[ok], angle = 90, code = 3,
                                  length = 0.03, col = cols[g])
  }
  if (length(groups) > 1) graphics::legend("topright", legend = groups, col = cols, lwd = 2, bty = "n", cex = 0.8)
  grDevices::dev.off(); files <- c(files, "rgr.png")

  if (!is.null(fun$curve)) {
    open_png(out_dir, "rgr_vs_size.png")
    graphics::plot(fun$curve$mass_fit_g, fun$curve$rgr_g_g_d, type = "n", log = "x",
                   xlab = "Fitted total dry mass (g, log scale)", ylab = expression(RGR ~ (g ~ g^-1 ~ d^-1)))
    for (g in intersect(groups, fun$curve$group)) {
      cg <- fun$curve[fun$curve$group == g, ]
      graphics::lines(cg$mass_fit_g, cg$rgr_g_g_d, col = cols[g], lwd = 2)
    }
    if (length(groups) > 1) graphics::legend("topright", legend = groups, col = cols, lwd = 2, bty = "n", cex = 0.8)
    grDevices::dev.off(); files <- c(files, "rgr_vs_size.png")
  }
  files
}

plot_allocation <- function(hs, out_dir) {
  vars <- list(c("sla_m2_kg", "sla_se", "SLA (m2/kg)"), c("lar_m2_kg", "lar_se", "LAR (m2/kg)"),
               c("lmf_g_g", "lmf_se", "LMF (g/g)"))
  vars <- Filter(function(v) any(is.finite(hs[[v[1]]])), vars)
  if (!length(vars)) return(character())
  groups <- sort(unique(hs$group)); cols <- group_cols(groups)
  open_png(out_dir, "allocation.png", width = 3.3 * length(vars) + 0.5, height = 3.6)
  op <- graphics::par(mfrow = c(1, length(vars)), mar = c(4, 4.2, 1, 1))
  for (v in vars) {
    yl <- range(c(hs[[v[1]]] - hs[[v[2]]], hs[[v[1]]] + hs[[v[2]]], hs[[v[1]]]), na.rm = TRUE)
    graphics::plot(range(hs$time), yl, type = "n", xlab = "Time (d)", ylab = v[3])
    for (g in groups) {
      s <- hs[hs$group == g, ]
      graphics::lines(s$time, s[[v[1]]], type = "b", pch = 16, col = cols[g])
      ok <- is.finite(s[[v[2]]]) & s[[v[2]]] > 0
      if (any(ok)) graphics::arrows(s$time[ok], s[[v[1]]][ok] - s[[v[2]]][ok], s$time[ok],
                                    s[[v[1]]][ok] + s[[v[2]]][ok], angle = 90, code = 3,
                                    length = 0.02, col = cols[g])
    }
  }
  if (length(groups) > 1) graphics::legend("topright", legend = groups, col = cols, pch = 16, bty = "n", cex = 0.8)
  graphics::par(op)
  grDevices::dev.off()
  "allocation.png"
}

# --- Main -------------------------------------------------------------------------

fmt <- function(x, d = 3) ifelse(is.finite(x), formatC(x, format = "f", digits = d), "NA")

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (opt$demo) opt$input <- file.path(SKILL_DIR, "demo", "harvests.csv")
  if (!nzchar(opt$input)) stop("Provide --input FILE.csv, or --demo. See --help.", call. = FALSE)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  if (opt$k_max < 3) stop("--k-max must be at least 3.", call. = FALSE)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  set.seed(opt$seed)

  warnings_list <- character()
  warn <- function(msg) warnings_list <<- c(warnings_list, msg)
  d <- prepare_data(read_input(opt$input), opt$group, opt$time_col, opt$mass_col,
                    opt$leaf_area_col, opt$leaf_mass_col, warn)
  has_la <- any(is.finite(d$la)); has_lm <- any(is.finite(d$lm))
  if (!has_la) warn(sprintf("No leaf area column ('%s'): NAR, LAR and SLA were not computed.", opt$leaf_area_col))
  if (!has_lm) warn(sprintf("No leaf mass column ('%s'): SLA and LMF were not computed.", opt$leaf_mass_col))

  hs <- harvest_summary(d)
  cls <- classical_intervals(d, n_boot = opt$n_boot)
  ov <- overall_rgr(d)
  fun <- functional_rgr(d, k_max = opt$k_max)
  save_table(hs, opt$out, "harvest_summary.csv")
  save_table(cls, opt$out, "classical_intervals.csv")
  save_table(ov, opt$out, "overall_rgr.csv")
  if (!is.null(fun$curve)) save_table(fun$curve, opt$out, "functional_rgr.csv")
  figs <- c(plot_growth(d, cls, fun, opt$out), if (has_la || has_lm) plot_allocation(hs, opt$out))

  # Warnings and red flags
  small <- hs[hs$n < opt$min_n, ]
  if (nrow(small)) warn(sprintf("%d harvest(s) have fewer than %d plants (e.g. %s at day %g, n = %d): RGR standard errors are unreliable.",
                                nrow(small), opt$min_n, small$group[1], small$time[1], small$n[1]))
  negs <- cls[cls$rgr_g_g_d <= 0, ]
  if (nrow(negs)) warn(sprintf("%d interval(s) have RGR <= 0 (e.g. %s, day %g-%g): check harvest labels, senescence or sampling of different plant sizes.",
                               nrow(negs), negs$group[1], negs$t1[1], negs$t2[1]))
  if (any(cls$rgr_g_g_d > 0.5) || any(ov$rgr_g_g_d > 0.5))
    warn("RGR above 0.5 g g-1 d-1 is unusually high for whole plants: confirm that time is in days and harvests are correctly labeled.")
  bias <- cls$rgr_ln_of_means_g_g_d - cls$rgr_g_g_d
  big <- which(abs(bias) > 0.1 * abs(cls$rgr_g_g_d) & abs(bias) > 0.005)
  if (length(big)) warn(sprintf("In %d interval(s) the variance of ln mass changed between harvests, so RGR from ln(mean mass) would differ by up to %.3f g g-1 d-1 (Hoffmann & Poorter 2002). This report uses the unbiased mean of ln mass.",
                                length(big), max(abs(bias[big]))))
  curved <- ov$group[is.finite(ov$p_curvature) & ov$p_curvature < 0.05]
  if (length(curved)) warn(sprintf("ln mass is curved over time in %s (quadratic p < 0.05): RGR changed during the experiment, so the overall slope is only an average. Use the interval or functional RGR.",
                                   paste(curved, collapse = ", ")))
  first <- hs[!duplicated(hs$group), ]
  if (nrow(first) > 1 && diff(range(first$ln_mass_mean)) > log(1.3))
    warn(sprintf("Groups differ by %.0f%% in mass at the first harvest. RGR falls as plants grow, so compare groups at a common size (figures/rgr_vs_size.png), not only at common dates.",
                 100 * (exp(diff(range(first$ln_mass_mean))) - 1)))
  if (any(!is.na(fun$info$note) & nzchar(fun$info$note)))
    warn(sprintf("No functional fit for: %s.", paste(fun$info$group[nzchar(fun$info$note)], collapse = ", ")))
  p_groups <- NA_real_
  if (length(unique(d$group)) > 1) {
    f0 <- stats::lm(log(mass) ~ time + group, data = d)
    f1 <- stats::lm(log(mass) ~ time * group, data = d)
    p_groups <- stats::anova(f0, f1)[2, "Pr(>F)"]
  }

  hs_lines <- c("| group | day | n | mean mass (g) | geometric mean (g) | sd ln mass | SLA (m2/kg) | LAR (m2/kg) | LMF (g/g) |",
                "|---|---|---|---|---|---|---|---|---|",
                sprintf("| %s | %g | %d | %s | %s | %s | %s | %s | %s |", hs$group, hs$time, hs$n,
                        fmt(hs$mass_mean_g), fmt(hs$mass_geomean_g), fmt(hs$ln_mass_sd, 2),
                        fmt(hs$sla_m2_kg, 1), fmt(hs$lar_m2_kg, 1), fmt(hs$lmf_g_g, 2)))
  cl_lines <- c("| group | interval (d) | RGR (g/g/d) | 95% CI | RGR from ln(mean) (biased) | NAR (g/m2/d) | NAR 95% CI |",
                "|---|---|---|---|---|---|---|",
                sprintf("| %s | %g-%g | %s | [%s, %s] | %s | %s | [%s, %s] |", cls$group, cls$t1, cls$t2,
                        fmt(cls$rgr_g_g_d, 4), fmt(cls$rgr_lo, 4), fmt(cls$rgr_hi, 4),
                        fmt(cls$rgr_ln_of_means_g_g_d, 4), fmt(cls$nar_g_m2_d, 2),
                        fmt(cls$nar_lo, 2), fmt(cls$nar_hi, 2)))
  ov_lines <- c("| group | harvests | RGR (g/g/d) | 95% CI | R2 | curvature p |", "|---|---|---|---|---|---|",
                sprintf("| %s | %d | %s | [%s, %s] | %s | %s |", ov$group, ov$n_harvests, fmt(ov$rgr_g_g_d, 4),
                        fmt(ov$rgr_lo, 4), fmt(ov$rgr_hi, 4), fmt(ov$r2, 3),
                        ifelse(is.finite(ov$p_curvature), formatC(ov$p_curvature, format = "g", digits = 3), "NA")))
  fun_lines <- NULL
  if (!is.null(fun$curve)) {
    fh <- fun$curve[fun$curve$at_harvest, ]
    fun_lines <- c("RGR(t) at the harvest times:", "",
                   "| group | day | RGR (g/g/d) | 95% band |", "|---|---|---|---|",
                   sprintf("| %s | %g | %s | [%s, %s] |", fh$group, fh$time, fmt(fh$rgr_g_g_d, 4),
                           fmt(fh$rgr_lo, 4), fmt(fh$rgr_hi, 4)), "",
                   sprintf("Spline basis size k and effective degrees of freedom: %s.",
                           paste(sprintf("%s k = %s, edf = %s", fun$info$group, fun$info$k, fmt(fun$info$edf, 2)), collapse = "; ")),
                   "")
  }
  report <- c(sprintf("# %s report", SKILL), "",
              if (length(warnings_list)) c("## Warnings", "", paste("-", warnings_list), "") else NULL,
              "## What was done", "",
              sprintf("Input: `%s`; %d plants, %d group(s), %d harvest time(s). Time column `%s` (days), mass `%s` (g)%s%s.",
                      opt$input, nrow(d), length(unique(d$group)), length(unique(d$time)), opt$time_col, opt$mass_col,
                      if (has_la) sprintf(", leaf area `%s` (cm2)", opt$leaf_area_col) else "",
                      if (has_lm) sprintf(", leaf mass `%s` (g)", opt$leaf_mass_col) else ""), "",
              "- Classical RGR between consecutive harvests from the mean of ln mass (Hoffmann & Poorter 2002); plants are not paired, so SE = sqrt(s1^2/n1 + s2^2/n2)/dt with Welch df.",
              sprintf("- NAR (unit leaf rate) from harvest means (Williams 1946; Radford 1967; assumes mass linear in leaf area within an interval), %g bootstrap resamples (seed %g) for its CI.", opt$n_boot, opt$seed),
              "- SLA, LAR and LMF as means of per-plant ratios at each harvest (SE across plants).",
              sprintf("- Functional RGR(t): penalized cubic regression spline of ln mass on time (mgcv, REML, k <= %g), derivative with 95%% bands.", opt$k_max),
              "",
              "## Harvest summary", "", hs_lines, "",
              "## Classical growth analysis", "", cl_lines, "",
              "## Mean RGR over the experiment", "",
              "Slope of ln mass on time (an average if RGR changed; see curvature p).", "", ov_lines, "",
              if (is.finite(p_groups)) c(sprintf("Equal-slope test across groups (time x group interaction): p = %s.", formatC(p_groups, format = "g", digits = 3)), "") else NULL,
              "## Functional growth analysis", "", fun_lines,
              sprintf("![%s](figures/%s)", sub("\\.png$", "", figs), figs), "",
              "## Limitations", "",
              "- Destructive harvests: each plant is measured once, so individual growth rates cannot be computed and plants must not be paired across harvests.",
              "- RGR declines as plants grow (ontogenetic drift); differences between groups of different size are partly size effects.",
              "- NAR assumes mass is linear in leaf area within each interval; the functional bands are Bayesian intervals with good average, not exact pointwise, coverage.",
              "",
              "## References", "",
              "- Hoffmann WA, Poorter H (2002) Ann Bot 90:37-42. https://doi.org/10.1093/aob/mcf140",
              "- Radford PJ (1967) Crop Sci 7:171-175. https://doi.org/10.2135/cropsci1967.0011183X000700030001x",
              "- Williams RF (1946) Ann Bot 10:41-72. https://doi.org/10.1093/oxfordjournals.aob.a083119",
              "- Hunt R, Causton DR, Shipley B, Askew AP (2002) Ann Bot 90:485-488. https://doi.org/10.1093/aob/mcf214",
              "- Wood SN (2017) Generalized Additive Models, 2nd edn. CRC. https://doi.org/10.1201/9781315370279",
              "",
              sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(report, file.path(opt$out, "report.md"))
  write_repro_bundle(opt$out, argv, inputs = opt$input)
  cat(sprintf("Done: %d group(s), %d interval(s). Report: %s\n", length(unique(d$group)), nrow(cls),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(list(harvests = hs, classical = cls, overall = ov, functional = fun))
}

if (sys.nframe() == 0L) main()
