#!/usr/bin/env Rscript
# claw-light-response: fit photosynthetic light-response (A-Q) curves with the
# non-rectangular hyperbola (Prioul & Chartier 1977; Marshall & Biscoe 1980)
# using minpack.lm::nlsLM, and report light-saturated net assimilation (Asat),
# apparent quantum yield, curvature, dark respiration, light compensation
# point and a saturation point defined as the PPFD where A reaches a set
# fraction of Asat. Includes low-light coverage and photoinhibition checks.
#
# Input: the tidy CSV written by claw-licor-6800 (or any CSV with the same
# column names and units).
#
# Usage:
#   Rscript scripts/light_response.R --demo
#   Rscript scripts/light_response.R --input licor_tidy.csv --curve-cols genotype,plot
#   Rscript scripts/light_response.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-light-response"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/light_response.R [--demo] [--input FILE.csv] [--curve-cols COLS]
                                     [--absorptance F] [--sat-fraction F] [--theta F]
                                     [--low-light PPFD] [--min-low-points N]
                                     [--photoinhibition-drop F] [--exclude-unstable]
                                     [--include-flagged] [--out DIR]

  --demo                  fit the bundled synthetic curves (demo/light_demo.csv)
  --input                 tidy CSV from claw-licor-6800 (one row per logged point)
  --curve-cols            column(s) that together identify one curve, comma-separated
  --absorptance           leaf absorptance (0-1), if measured; adds quantum yield on an
                          absorbed-light basis (default: incident basis only)
  --sat-fraction          saturation point = PPFD where A reaches this fraction of Asat
                          (default 0.9; rule of thumb)
  --theta                 fix the curvature theta (0-1) instead of fitting it
  --low-light             upper PPFD of the initial-slope region (default 150; rule of thumb)
  --min-low-points        points needed between 0 and --low-light for the quantum yield (default 3)
  --photoinhibition-drop  flag when A at the highest PPFD is this fraction below the
                          maximum observed A (default 0.10; rule of thumb)
  --exclude-unstable      exclude points logged before stability (qc_unstable)
  --include-flagged       keep points with qc_hard_fail = TRUE (not advised)
  --out                   output folder (default: agriclaw_output/claw-light-response)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, input = "", curve_cols = "", absorptance = NA_real_,
                 sat_fraction = 0.9, theta = NA_real_, low_light = 150, min_low_points = 3,
                 photoinhibition_drop = 0.10, exclude_unstable = FALSE,
                 include_flagged = FALSE, out = "")

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

# Non-rectangular hyperbola for gross assimilation (smaller root of
# theta*P^2 - (phi*Q + Pmax)*P + phi*Q*Pmax = 0); net A = gross - Rd.
# theta -> 0 gives the rectangular hyperbola, theta = 1 the Blackman response.
nrh_gross <- function(q, phi, theta, pmax) {
  x <- phi * q + pmax
  if (theta < 1e-6) return(phi * q * pmax / x)
  (x - sqrt(base::pmax(x^2 - 4 * theta * phi * q * pmax, 0))) / (2 * theta)
}
nrh_net <- function(q, phi, theta, pmax, rd) nrh_gross(q, phi, theta, pmax) - rd

# PPFD at which gross assimilation equals g (0 < g < pmax); exact inverse of
# the hyperbola: Q = g * (pmax - theta * g) / (phi * (pmax - g)).
q_for_gross <- function(g, phi, theta, pmax) {
  out <- g * (pmax - theta * g) / (phi * (pmax - g))
  out[!(g > 0 & g < pmax)] <- NA_real_
  out
}
# Light compensation point: net A = 0, i.e. gross = Rd.
lcp_of <- function(phi, theta, pmax, rd) q_for_gross(rd, phi, theta, pmax)
# Saturation point: net A = f * Asat, i.e. gross = Rd + f * (pmax - Rd).
qsat_of <- function(phi, theta, pmax, rd, f = 0.9) q_for_gross(rd + f * (pmax - rd), phi, theta, pmax)

split_cols <- function(s) {
  if (!nzchar(s)) return(character())
  trimws(strsplit(s, ",", fixed = TRUE)[[1L]])
}

read_tidy <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
  if (tolower(tools::file_ext(path)) %in% c("xlsx", "", "txt"))
    stop("'", basename(path), "' looks like a raw LI-6800 log. Run claw-licor-6800 first and ",
         "pass its tables/licor_tidy.csv here.", call. = FALSE)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

validate_input <- function(df, curve_cols, opt) {
  need <- c("qin_umol_m2_s", "a_umol_m2_s")
  miss <- setdiff(need, names(df))
  if (length(miss))
    stop("Missing column(s): ", paste(miss, collapse = ", "), ". The input must be the tidy CSV from ",
         "claw-licor-6800 (incident PPFD in umol m-2 s-1, A in umol m-2 s-1).", call. = FALSE)
  if (!length(curve_cols))
    stop("--curve-cols is required: which column(s) identify one light-response curve (e.g. genotype,plot)?",
         call. = FALSE)
  miss <- setdiff(curve_cols, names(df))
  if (length(miss)) stop("--curve-cols column(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
  q <- df$qin_umol_m2_s
  if (!is.numeric(q) || !is.numeric(df$a_umol_m2_s)) stop("qin_umol_m2_s and a_umol_m2_s must be numeric.", call. = FALSE)
  if (max(q, na.rm = TRUE) < 100)
    stop("Maximum PPFD is ", max(q, na.rm = TRUE), ": qin_umol_m2_s must be in umol m-2 s-1 ",
         "(not mol or W m-2), and a light curve should reach saturating light.", call. = FALSE)
  if (max(q, na.rm = TRUE) > 4000)
    stop("PPFD above 4000 umol m-2 s-1 (full sunlight is about 2000-2300): check the units of qin_umol_m2_s.",
         call. = FALSE)
  if (any(q < -5, na.rm = TRUE)) stop("Negative PPFD values: check qin_umol_m2_s.", call. = FALSE)
  if (max(df$a_umol_m2_s, na.rm = TRUE) > 100) stop("A exceeds 100 umol m-2 s-1: check units.", call. = FALSE)
  if (!is.na(opt$absorptance) && (opt$absorptance <= 0 || opt$absorptance > 1))
    stop("--absorptance must be between 0 and 1 (e.g. 0.85).", call. = FALSE)
  if (!is.na(opt$theta) && (opt$theta < 0 || opt$theta > 1)) stop("--theta must be between 0 and 1.", call. = FALSE)
  if (opt$sat_fraction <= 0 || opt$sat_fraction >= 1) stop("--sat-fraction must be between 0 and 1.", call. = FALSE)
  invisible(TRUE)
}

# Fit one curve with nlsLM from several starting curvatures; keep the best.
fit_nrh <- function(q, a, theta_fixed = NA_real_, low_light = 150) {
  ok <- is.finite(q) & is.finite(a)
  q <- q[ok]
  a <- a[ok]
  if (length(q) < 5L) stop("fewer than 5 usable points", call. = FALSE)
  low <- q <= low_light
  if (sum(low) >= 2L) {
    cf <- stats::coef(stats::lm(a[low] ~ q[low]))
    rd0 <- max(0.2, -cf[[1]])
    phi0 <- min(0.12, max(0.01, cf[[2]]))
  } else {
    rd0 <- max(0.2, -min(a))
    phi0 <- 0.05
  }
  pmax0 <- max(a) + rd0
  thetas <- if (is.na(theta_fixed)) c(0.3, 0.7, 0.95) else theta_fixed
  best <- NULL
  for (th in thetas) {
    fit <- tryCatch({
      if (is.na(theta_fixed)) {
        minpack.lm::nlsLM(a ~ nrh_net(q, phi, theta, pmax, rd),
                          start = list(phi = phi0, theta = th, pmax = pmax0, rd = rd0),
                          lower = c(1e-4, 1e-3, 1e-3, -10), upper = c(0.3, 1, 3 * pmax0 + 50, 50),
                          control = minpack.lm::nls.lm.control(maxiter = 1000, ftol = 1e-12, ptol = 1e-12))
      } else {
        minpack.lm::nlsLM(a ~ nrh_net(q, phi, theta_fixed, pmax, rd),
                          start = list(phi = phi0, pmax = pmax0, rd = rd0),
                          lower = c(1e-4, 1e-3, -10), upper = c(0.3, 3 * pmax0 + 50, 50),
                          control = minpack.lm::nls.lm.control(maxiter = 1000, ftol = 1e-12, ptol = 1e-12))
      }
    }, error = function(e) NULL)
    if (!is.null(fit) && (is.null(best) || stats::deviance(fit) < stats::deviance(best))) best <- fit
  }
  if (is.null(best)) stop("nlsLM did not converge from any starting value", call. = FALSE)
  best
}

# Parameter estimates, Wald SEs/CIs and derived quantities for one fit.
summarise_fit <- function(fit, q, a, opt) {
  est <- stats::coef(fit)
  theta_fit <- is.na(opt$theta)
  th <- if (theta_fit) est[["theta"]] else opt$theta
  vc <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  se <- function(p) if (!is.null(vc) && p %in% rownames(vc)) sqrt(max(vc[p, p], 0)) else NA_real_
  df_res <- stats::df.residual(fit)
  tq <- if (df_res > 0) stats::qt(0.975, df_res) else NA_real_
  phi <- est[["phi"]]
  pmax <- est[["pmax"]]
  rd <- est[["rd"]]
  asat <- pmax - rd
  asat_se <- if (!is.null(vc)) sqrt(max(vc["pmax", "pmax"] + vc["rd", "rd"] - 2 * vc["pmax", "rd"], 0)) else NA_real_
  res <- stats::residuals(fit)
  out <- data.frame(
    asat_umol_m2_s = asat, asat_se = asat_se,
    asat_lower_umol_m2_s = asat - tq * asat_se, asat_upper_umol_m2_s = asat + tq * asat_se,
    phi_incident_mol_mol = phi, phi_se = se("phi"),
    phi_lower_mol_mol = phi - tq * se("phi"), phi_upper_mol_mol = phi + tq * se("phi"),
    theta = th, theta_se = if (theta_fit) se("theta") else NA_real_, theta_fitted = theta_fit,
    rd_umol_m2_s = rd, rd_se = se("rd"),
    rd_lower_umol_m2_s = rd - tq * se("rd"), rd_upper_umol_m2_s = rd + tq * se("rd"),
    pmax_gross_umol_m2_s = pmax,
    lcp_umol_m2_s = lcp_of(phi, th, pmax, rd),
    qsat_umol_m2_s = qsat_of(phi, th, pmax, rd, opt$sat_fraction),
    sat_fraction = opt$sat_fraction,
    a_model_at_qmax_umol_m2_s = nrh_net(max(q), phi, th, pmax, rd),
    a_max_observed_umol_m2_s = max(a), q_at_a_max_observed = q[which.max(a)],
    q_max_umol_m2_s = max(q), n_points = length(q), dof = df_res,
    rmse_umol_m2_s = sqrt(mean(res^2)), stringsAsFactors = FALSE)
  if (!is.na(opt$absorptance)) {
    out$absorptance <- opt$absorptance
    out$phi_absorbed_mol_mol <- phi / opt$absorptance
  }
  out
}

# Curve QC; thresholds are rules of thumb (see SKILL.md).
light_qc <- function(q, a, s, opt) {
  flags <- character()
  n_low <- sum(q > 0 & q <= opt$low_light)
  if (n_low < opt$min_low_points) flags <- c(flags, "few_low_light_points")
  if (min(q) > 50) flags <- c(flags, "no_dark_point")
  qa <- q[which.max(a)]
  a_top <- mean(a[q == max(q)])
  if (qa < max(q) && a_top < (1 - opt$photoinhibition_drop) * max(a)) flags <- c(flags, "photoinhibition")
  if (!is.null(s)) {
    if (isTRUE(s$theta_fitted) && (s$theta >= 0.999 || s$theta <= 0.002)) flags <- c(flags, "theta_at_bound")
    if (s$phi_incident_mol_mol > 0.125) flags <- c(flags, "phi_above_theoretical_max")
    if (s$rd_umol_m2_s < 0) flags <- c(flags, "negative_rd")
    if (!is.na(s$qsat_umol_m2_s) && s$qsat_umol_m2_s > max(q)) flags <- c(flags, "qsat_beyond_measured_range")
    if (is.na(s$lcp_umol_m2_s)) flags <- c(flags, "lcp_undefined")
    if (s$rmse_umol_m2_s > 0.05 * max(abs(a))) flags <- c(flags, "poor_fit_rmse")
  }
  list(n_low = n_low, flags = paste(flags, collapse = ";"))
}

curve_key <- function(df, curve_cols) {
  do.call(paste, c(lapply(curve_cols, function(cn) as.character(df[[cn]])), sep = " - "))
}

analyse_light <- function(df, opt) {
  curve_cols <- split_cols(opt$curve_cols)
  validate_input(df, curve_cols, opt)
  df$curve_id <- curve_key(df, curve_cols)
  reason <- rep("", nrow(df))
  add <- function(cond, why) {
    cond <- cond %in% TRUE
    reason[cond] <<- ifelse(nzchar(reason[cond]), paste(reason[cond], why, sep = ";"), why)
  }
  if (!opt$include_flagged && "qc_hard_fail" %in% names(df)) add(df$qc_hard_fail, "qc_hard_fail")
  if (opt$exclude_unstable && "qc_unstable" %in% names(df)) add(df$qc_unstable, "unstable")
  add(is.na(df$a_umol_m2_s) | is.na(df$qin_umol_m2_s), "missing_A_or_Q")
  df$excluded_reason <- reason
  df$used_in_fit <- !nzchar(reason)
  ids <- sort(unique(df$curve_id))
  params <- list()
  fitted <- list()
  for (id in ids) {
    d <- df[df$curve_id == id & df$used_in_fit, , drop = FALSE]
    q <- d$qin_umol_m2_s
    a <- d$a_umol_m2_s
    fit <- tryCatch(fit_nrh(q, a, opt$theta, opt$low_light), error = function(e) e)
    if (inherits(fit, "error")) {
      qc <- if (length(q)) light_qc(q, a, NULL, opt) else list(n_low = 0L, flags = "")
      params[[id]] <- data.frame(curve_id = id, n_points = length(q), n_low_light_points = qc$n_low,
                                 fit_error = conditionMessage(fit),
                                 qc_flags = paste(c(qc$flags[nzchar(qc$flags)], "fit_failed"), collapse = ";"),
                                 stringsAsFactors = FALSE)
      next
    }
    s <- summarise_fit(fit, q, a, opt)
    qc <- light_qc(q, a, s, opt)
    params[[id]] <- data.frame(curve_id = id, s, n_low_light_points = qc$n_low,
                               tleaf_mean_c = if ("tleaf_c" %in% names(d)) mean(d$tleaf_c, na.rm = TRUE) else NA_real_,
                               fit_error = "", qc_flags = qc$flags, stringsAsFactors = FALSE)
    fitted[[id]] <- data.frame(curve_id = id, qin_umol_m2_s = q, a_umol_m2_s = a,
                               a_fit_umol_m2_s = stats::fitted(fit),
                               residual_umol_m2_s = stats::residuals(fit), stringsAsFactors = FALSE)
  }
  p <- rbind_fill(params)
  meta <- curve_metadata(df, curve_cols)
  p <- merge(meta, p, by = "curve_id", all.y = TRUE, sort = FALSE)
  p <- p[order(p$curve_id), , drop = FALSE]
  rownames(p) <- NULL
  list(params = p, fits = if (length(fitted)) do.call(rbind, fitted) else NULL, points = df)
}

rbind_fill <- function(dfs) {
  cols <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function(d) {
    for (cn in setdiff(cols, names(d))) d[[cn]] <- NA
    d[, cols, drop = FALSE]
  })
  out <- do.call(rbind, dfs)
  rownames(out) <- NULL
  out
}

curve_metadata <- function(points, curve_cols) {
  cand <- setdiff(names(points), c("curve_id", "used_in_fit", "excluded_reason"))
  cand <- cand[!startsWith(cand, "qc_")]
  cand <- cand[!vapply(points[cand], function(v) is.numeric(v) || is.logical(v), logical(1)) |
                 cand %in% curve_cols]
  keep <- cand[vapply(cand, function(cn) all(tapply(points[[cn]], points$curve_id,
                                                    function(v) length(unique(v)) == 1L)), logical(1))]
  keep <- setdiff(keep, c("date_local", "hhmmss", "user_remark", "file_type"))
  m <- points[!duplicated(points$curve_id), c("curve_id", keep), drop = FALSE]
  rownames(m) <- NULL
  m
}

plot_fits <- function(res, opt, out_dir) {
  p <- res$params
  ids <- p$curve_id
  pages <- split(ids, ceiling(seq_along(ids) / 9L))
  files <- character()
  for (pg in seq_along(pages)) {
    nm <- sprintf("light_fits_%02d.png", pg)
    open_png(out_dir, nm, width = 10, height = 9)
    graphics::par(mfrow = c(3, 3), mar = c(4, 4, 2.6, 0.8), cex = 0.75)
    for (id in pages[[pg]]) {
      d <- res$points[res$points$curve_id == id, , drop = FALSE]
      pr <- p[p$curve_id == id, , drop = FALSE]
      graphics::plot(d$qin_umol_m2_s, d$a_umol_m2_s, type = "n",
                     xlab = "Incident PPFD (umol m-2 s-1)", ylab = "A (umol m-2 s-1)",
                     main = substr(id, 1L, 40L))
      graphics::abline(h = 0, col = "grey80")
      if (is.finite(pr$asat_umol_m2_s %||% NA)) {
        qq <- seq(0, max(d$qin_umol_m2_s, na.rm = TRUE), length.out = 300)
        graphics::lines(qq, nrh_net(qq, pr$phi_incident_mol_mol, pr$theta, pr$pmax_gross_umol_m2_s, pr$rd_umol_m2_s), lwd = 2)
        graphics::abline(h = pr$asat_umol_m2_s, lty = 3, col = "grey40")
        if (is.finite(pr$lcp_umol_m2_s)) graphics::points(pr$lcp_umol_m2_s, 0, pch = 4, col = "#1f77b4", cex = 1.4, lwd = 2)
        if (is.finite(pr$qsat_umol_m2_s) && pr$qsat_umol_m2_s <= max(d$qin_umol_m2_s))
          graphics::abline(v = pr$qsat_umol_m2_s, lty = 2, col = "#d95f02")
      }
      u <- d$used_in_fit
      graphics::points(d$qin_umol_m2_s[u], d$a_umol_m2_s[u], pch = 16)
      graphics::points(d$qin_umol_m2_s[!u], d$a_umol_m2_s[!u], pch = 1, col = "grey40")
      if (nzchar(pr$qc_flags %||% "")) graphics::mtext(substr(pr$qc_flags, 1L, 60L), side = 3, line = 0.1,
                                                       cex = 0.55, col = "firebrick")
    }
    grDevices::dev.off()
    files <- c(files, nm)
  }
  files
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || (length(a) == 1L && is.na(a))) b else a

fmt_table <- function(df, digits = 3) {
  if (is.null(df) || !nrow(df)) return(character())
  f <- vapply(df, function(v) if (is.numeric(v)) formatC(v, digits = digits, format = "g") else as.character(v),
              character(nrow(df)))
  if (is.null(dim(f))) f <- matrix(f, nrow = 1L)
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(f, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

FLAG_TEXT <- c(
  few_low_light_points = "too few points between 0 and the low-light limit: the quantum yield is poorly constrained",
  no_dark_point = "no point below 50 umol m-2 s-1: Rd and the compensation point are extrapolated",
  photoinhibition = "A at the highest PPFD is well below the maximum observed A (photoinhibition or a non-steady state); the hyperbola cannot decline, so refit without the points beyond the peak or treat Asat with caution",
  theta_at_bound = "the curvature theta hit its bound (0 or 1): it is not identifiable from these points; consider fixing it with --theta",
  phi_above_theoretical_max = "the quantum yield exceeds the theoretical maximum of 0.125 mol CO2 per mol photons: check PPFD units and the low-light points",
  negative_rd = "negative Rd: respiration is poorly constrained",
  qsat_beyond_measured_range = "the saturation point lies above the highest PPFD measured, so it is extrapolated",
  lcp_undefined = "the light compensation point is undefined (Rd <= 0 or Rd >= gross Pmax)",
  poor_fit_rmse = "RMSE above 5% of the maximum A",
  fit_failed = "the fit failed (see fit_error)")

write_report <- function(res, opt, input, out_dir, figs) {
  p <- res$params
  w <- c(if (is.na(opt$absorptance))
    "The quantum yield is on an incident-light basis (PPFD reaching the leaf), not absorbed light; it is lower than the absorbed-basis yield by the leaf absorptance. Pass --absorptance only if absorptance was measured; the LI-6800 Qabs column uses an assumed absorptance.")
  for (fl in names(FLAG_TEXT)) {
    hit <- p$curve_id[grepl(paste0("(^|;)", fl, "($|;)"), p$qc_flags)]
    if (length(hit)) w <- c(w, sprintf("%s%s (%d curve(s): %s).", toupper(substr(FLAG_TEXT[[fl]], 1, 1)),
                                       substring(FLAG_TEXT[[fl]], 2), length(hit), paste(hit, collapse = "; ")))
  }
  show <- p[, intersect(c("curve_id", "asat_umol_m2_s", "phi_incident_mol_mol", "theta", "rd_umol_m2_s",
                          "lcp_umol_m2_s", "qsat_umol_m2_s", "rmse_umol_m2_s", "qc_flags"), names(p))]
  names(show) <- c("curve", "Asat", "phi (incident)", "theta", "Rd", "LCP", sprintf("Q%.0f", 100 * opt$sat_fraction),
                   "RMSE", "flags")[seq_len(ncol(show))]
  lines <- c(sprintf("# %s report", SKILL), "", "## Warnings", "",
             if (length(w)) paste("-", w) else "None.", "",
             "## Results", "",
             sprintf("Input: `%s`. %d curve(s); %d of %d points used.", input, nrow(p),
                     sum(res$points$used_in_fit), nrow(res$points)), "",
             "Asat, Rd, LCP and the saturation point in umol m-2 s-1; phi in mol CO2 per mol incident photons. Standard errors and 95% Wald intervals are in `tables/light_parameters.csv`.", "",
             fmt_table(show), "",
             paste0("![Light-response fits](figures/", figs[1L], ")"), "",
             "## Method", "",
             paste0("Each curve was fitted with the non-rectangular hyperbola A = [phi Q + Pmax - sqrt((phi Q + Pmax)^2 - 4 theta phi Q Pmax)] / (2 theta) - Rd by Levenberg-Marquardt least squares (minpack.lm::nlsLM, ",
                    if (is.na(opt$theta)) "theta fitted, best of three starting values" else sprintf("theta fixed at %g", opt$theta),
                    "). Asat = Pmax - Rd is the light-saturated net assimilation (the asymptote as PPFD rises without limit). ",
                    "The light compensation point solves A = 0 and the saturation point solves A = ", opt$sat_fraction,
                    " x Asat, both exactly: Q = g (Pmax - theta g) / (phi (Pmax - g)) with gross rate g = Rd or g = Rd + f (Pmax - Rd)."), "",
             sprintf("QC rules of thumb: at least %d points between 0 and %g umol m-2 s-1; a point below 50; photoinhibition when A at the highest PPFD is more than %.0f%% below the maximum observed A.",
                     opt$min_low_points, opt$low_light, 100 * opt$photoinhibition_drop), "",
             "## Limitations", "",
             "- The hyperbola has no photoinhibition term: declining A at high light biases Asat and theta.",
             "- Rd from a light curve is an extrapolation through the low-light points and includes the Kok effect; it is not day respiration (RL).",
             "- Each point is assumed to be at steady state; light curves stepped down quickly from high light often are not.", "",
             "## References", "",
             "- Prioul JL, Chartier P (1977) Partitioning of transfer and carboxylation components of intracellular resistance to photosynthetic CO2 fixation: a critical analysis of the methods used. Annals of Botany 41: 789-800. https://doi.org/10.1093/oxfordjournals.aob.a085354",
             "- Marshall B, Biscoe PV (1980) A model for C3 leaves describing the dependence of net photosynthesis on irradiance. Journal of Experimental Botany 31: 29-39. https://doi.org/10.1093/jxb/31.1.29",
             "- Ogren E, Evans JR (1993) Photosynthetic light-response curves. Planta 189: 182-190. https://doi.org/10.1007/BF00195075",
             "- Long SP, Postl WF, Bolhar-Nordenkampf HR (1993) Quantum yields for uptake of carbon dioxide in C3 vascular plants of contrasting habitats and taxonomic groupings. Planta 189: 226-234. https://doi.org/10.1007/BF00195081",
             "- Lobo F de A et al. (2013) Fitting net photosynthetic light-response curves with Microsoft Excel - a critical look at the models. Photosynthetica 51: 445-456. https://doi.org/10.1007/s11099-013-0045-y",
             "- Elzhov TV, Mullen KM, Spiess A-N, Bolker B. minpack.lm: R interface to the Levenberg-Marquardt nonlinear least-squares algorithm. https://CRAN.R-project.org/package=minpack.lm",
             "", sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(lines, file.path(out_dir, "report.md"))
  invisible(w)
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (!requireNamespace("minpack.lm", quietly = TRUE))
    stop("The R package minpack.lm is required: install.packages('minpack.lm').", call. = FALSE)
  if (opt$demo) {
    opt$input <- file.path(SKILL_DIR, "demo", "light_demo.csv")
    if (!nzchar(opt$curve_cols)) opt$curve_cols <- "genotype,plot"
  }
  if (!nzchar(opt$input)) stop("Provide --input FILE.csv (claw-licor-6800 output), or --demo. See --help.", call. = FALSE)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  res <- analyse_light(read_tidy(opt$input), opt)
  save_table(res$params, opt$out, "light_parameters.csv")
  if (!is.null(res$fits)) save_table(res$fits, opt$out, "light_fitted_points.csv")
  figs <- plot_fits(res, opt, opt$out)
  w <- write_report(res, opt, opt$input, opt$out, figs)
  write_repro_bundle(opt$out, argv, inputs = opt$input)
  cat(sprintf("Done: %d curve(s), %d warning(s). Report: %s\n", nrow(res$params), length(w),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(res)
}

if (sys.nframe() == 0L) main()
