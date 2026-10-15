#!/usr/bin/env Rscript
# claw-vpd-transpiration: whole-plant transpiration rate (TR) response to vapour
# pressure deficit (VPD) from gravimetric (pot-weight) phenotyping. Computes VPD
# from air temperature and RH, TR from pot-weight loss with irrigation and
# outlier screening, and compares a linear with a two-segment (breakpoint)
# TR-VPD model by AICc and the Davies test, to detect the "limited
# transpiration under high VPD" water-saving trait.
#
# Usage:
#   Rscript scripts/vpd_transpiration.R --demo
#   Rscript scripts/vpd_transpiration.R --weights weights.csv --plants plants.csv \
#       --env env.csv --steps steps.csv --out results/
#   Rscript scripts/vpd_transpiration.R --input tr_vpd.csv --tr-col tr_mg_m2_s
#   Rscript scripts/vpd_transpiration.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-vpd-transpiration"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/vpd_transpiration.R [--demo] [options]

Raw logger data (pot weights + air temperature/RH):
  --weights FILE      CSV: plant_id, timestamp, weight_g (one row per reading)
  --plants FILE       CSV: plant_id, genotype, leaf_area_cm2 [, chamber]
  --env FILE          CSV: timestamp, air_temp_c, rh_pct [, chamber]
  --steps FILE        CSV: step, start, end [, phase (day/night)] [, chamber];
                      without it, readings are binned into --bin-min windows
Precomputed data:
  --input FILE        CSV with --group, vpd_kpa (or air_temp_c + rh_pct) and --tr-col

Options:
  --group COL         grouping column for fits (default: genotype)
  --plant-col COL     plant identifier column (default: plant_id)
  --tr-col COL        TR column for --input (default: tr_mg_m2_s)
  --vpd-col COL       VPD column for --input (default: vpd_kpa)
  --level L           fit per 'genotype' (pooled plants, default) or per 'plant'
  --svp F             saturation vapour pressure: buck1981_moist (default),
                      buck1981 or fao56
  --selection S       'aicc+davies' (default) or 'aicc' (AICc only)
  --alpha A           Davies test level (default: 0.05)
  --ci-method M       breakpoint CI: score (default; Muggeo 2017) or delta
  --min-points N      fewest points to try a breakpoint model (default: 10)
  --min-side N        fewest points on each side of a breakpoint (default: 3)
  --irrigation-g G    weight rise per reading treated as irrigation (default: 5)
  --settle-min M      minutes excluded after an irrigation (default: 30)
  --outlier-k K       spike threshold in robust SDs (default: 3)
  --outlier-window N  running-median window, readings (default: 15)
  --gap-factor F      gaps longer than F x median interval split the series (default: 3)
  --bin-min M         window length when no --steps file (default: 60)
  --discard-first-min M  minutes dropped at the start of each step (default: 0)
  --env-offset-min M  minutes added to environment timestamps (default: 0)
  --min-readings N    fewest retained readings per plant and step (default: 10)
  --seed N            random seed for segmented's bootstrap restarts (default: 1)
  --out DIR           output folder (default: agriclaw_output/claw-vpd-transpiration)
  --demo              run on the bundled synthetic logger data
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, input = "", weights = "", plants = "", env = "",
                 steps = "", out = "", group = "genotype", plant_col = "plant_id",
                 tr_col = "tr_mg_m2_s", vpd_col = "vpd_kpa", level = "genotype",
                 svp = "buck1981_moist", selection = "aicc+davies", alpha = 0.05,
                 ci_method = "score", min_points = 10, min_side = 3,
                 irrigation_g = 5, settle_min = 30, outlier_k = 3,
                 outlier_window = 15, gap_factor = 3, bin_min = 60,
                 discard_first_min = 0, env_offset_min = 0, min_readings = 10,
                 seed = 1)

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

MW_WATER_MG_MMOL <- 18.015  # molar mass of water (mg per mmol)

# --- Input helpers -----------------------------------------------------------

read_csv_file <- function(path, what) {
  if (!nzchar(path)) stop("Provide --", what, " FILE.csv. See --help.", call. = FALSE)
  if (!file.exists(path)) stop(what, " file not found: ", path, call. = FALSE)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                  strip.white = TRUE)
}

need_cols <- function(dat, cols, what) {
  miss <- setdiff(cols, names(dat))
  if (length(miss))
    stop(what, ": column(s) not found: ", paste(miss, collapse = ", "),
         ". Available: ", paste(names(dat), collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

need_numeric <- function(dat, cols, what) {
  for (cl in cols) {
    if (!is.numeric(dat[[cl]]))
      stop(what, ": column '", cl, "' must be numeric (check for text such as ",
           "'NA ', units or decimal commas).", call. = FALSE)
    if (anyNA(dat[[cl]]))
      stop(what, ": column '", cl, "' has ", sum(is.na(dat[[cl]])),
           " missing value(s). Remove or fix those rows first; nothing is filled in.",
           call. = FALSE)
  }
  invisible(TRUE)
}

# Timestamps are read as clock times in UTC, so daylight-saving changes in the
# computer's time zone cannot shift them. Every logger must use the same clock.
parse_time <- function(x, what) {
  if (inherits(x, "POSIXct")) return(x)
  x <- trimws(as.character(x))
  fmts <- c("%Y-%m-%d %H:%M:%OS", "%Y-%m-%dT%H:%M:%OS", "%Y/%m/%d %H:%M:%OS",
            "%Y-%m-%d %H:%M", "%Y-%m-%dT%H:%M", "%Y/%m/%d %H:%M")
  out <- as.POSIXct(rep(NA_real_, length(x)), origin = "1970-01-01", tz = "UTC")
  for (f in fmts) {
    todo <- is.na(out)
    if (!any(todo)) break
    out[todo] <- as.POSIXct(x[todo], format = f, tz = "UTC")
  }
  if (anyNA(out))
    stop("Could not read ", what, " timestamp '", x[which(is.na(out))[1]],
         "'. Use ISO 8601 such as 2026-07-01 13:05:00.", call. = FALSE)
  out
}

# --- VPD ---------------------------------------------------------------------

# Saturation vapour pressure (kPa) over liquid water at air temperature t_c (degC).
#   buck1981       0.61121 exp(17.502 t / (240.97 + t))   Buck (1981), pure water vapour
#   buck1981_moist 0.61375 exp(17.502 t / (240.97 + t))   the same with an enhancement
#                  factor of ~1.0042 for moist air near 1000 hPa; the coefficient
#                  used in the published gravimetric workflow (Tamang & Sadok 2018)
#   fao56          0.6108 exp(17.27 t / (t + 237.3))      Allen et al. (1998), eq. 11
svp_kpa <- function(t_c, formula = "buck1981_moist") {
  switch(formula,
         buck1981 = 0.61121 * exp(17.502 * t_c / (240.97 + t_c)),
         buck1981_moist = 0.61375 * exp(17.502 * t_c / (240.97 + t_c)),
         fao56 = 0.6108 * exp(17.27 * t_c / (t_c + 237.3)),
         stop("Unknown --svp '", formula, "'. Use buck1981_moist, buck1981 or fao56.",
              call. = FALSE))
}

check_temp_rh <- function(t_c, rh_pct, where) {
  if (!is.numeric(t_c) || !is.numeric(rh_pct))
    stop(where, ": air_temp_c and rh_pct must be numeric.", call. = FALSE)
  if (anyNA(t_c) || anyNA(rh_pct))
    stop(where, ": ", sum(is.na(t_c) | is.na(rh_pct)), " reading(s) have missing ",
         "temperature or RH. Remove them first; nothing is filled in.", call. = FALSE)
  if (any(t_c < -10 | t_c > 60))
    stop(where, ": air_temp_c outside -10 to 60. Are these Fahrenheit or Kelvin? ",
         "Convert to degrees C.", call. = FALSE)
  if (max(rh_pct) <= 1)
    stop(where, ": rh_pct looks like a fraction (all values <= 1). Supply percent (0-100).",
         call. = FALSE)
  if (any(rh_pct < 0 | rh_pct > 103))
    stop(where, ": rh_pct outside 0-103 %. Check the sensor and the column.",
         call. = FALSE)
  invisible(TRUE)
}

# VPD (kPa) = SVP(T) x (1 - RH/100). RH readings between 100 and 103 % (sensor
# overshoot near saturation) are treated as 100 %; callers count and report them.
vpd_kpa <- function(t_c, rh_pct, formula = "buck1981_moist") {
  svp_kpa(t_c, formula) * (1 - pmin(rh_pct, 100) / 100)
}

# --- Pot-weight quality control -----------------------------------------------

# Flags the readings of one pot. Rules (each threshold is an option; h is half
# of outlier_window):
#  gap        a reading more than gap_factor x the median logging interval after
#             the previous one starts a new segment;
#  jump       a persistent level shift: the median of the next h readings minus
#             the median of the previous h readings, corrected for the local
#             downward trend (h x running median of the change per reading),
#             exceeds +irrigation_g (irrigation) or -irrigation_g (disturbance,
#             e.g. a pot knocked or a plant touched). A new segment starts at the
#             jump and readings in the next settle_min minutes are excluded.
#             Spikes of one to a few readings do not move the medians, so they
#             are not jumps;
#  outlier    within each jump-free segment, the series is detrended with the
#             running median of its loss rate, and readings whose residual from
#             a running median exceeds outlier_k x robust SD are spikes (Hampel
#             identifier; window outlier_window readings, robust SD = 1.4826 x
#             MAD of the residuals, floored at the balance resolution).
# Rates are never computed across segment boundaries.
odd_window <- function(k, n) {
  k <- min(as.integer(k), n)
  if (k %% 2L == 0L) k <- k - 1L
  k
}

qc_pot <- function(time, weight, irrigation_g = 5, settle_min = 30, outlier_k = 3,
                   outlier_window = 15, gap_factor = 3) {
  n <- length(weight)
  tmin <- as.numeric(difftime(time, time[1], units = "mins"))
  gap <- rep(FALSE, n)
  if (n > 1) {
    dt <- diff(tmin)
    if (any(dt <= 0)) stop("Duplicate or unsorted timestamps within a pot.", call. = FALSE)
    gap[-1] <- dt > gap_factor * stats::median(dt)
  }
  chunk <- cumsum(gap) + 1L
  h <- max(2L, as.integer(outlier_window) %/% 2L)
  jump_up <- jump_down <- rep(FALSE, n)
  for (ch in unique(chunk)) {
    idx <- which(chunk == ch)
    m <- length(idx)
    if (m < 2L * h + 1L) next
    x <- weight[idx]
    d <- diff(x)
    trend <- stats::runmed(d, odd_window(2L * h + 1L, length(d)), endrule = "median")
    med_from <- apply(stats::embed(x, h), 1, stats::median)  # median of x[r:(r+h-1)]
    i <- (h + 1L):(m - h + 1L)                               # i = first reading after a jump
    shift <- med_from[i] - med_from[i - h] - h * trend[i - 1L]
    cand <- i[abs(shift) > irrigation_g]
    if (!length(cand)) next
    for (run in split(cand, cumsum(c(1L, diff(cand) != 1L)))) {
      s <- shift[match(run, i)]
      sgn <- sign(s[which.max(abs(s))])
      # the jump is the largest single-reading step in that direction near the run
      p <- max(2L, min(run) - h):min(m, max(run) + h)
      j <- p[which.max(sgn * (x[p] - x[p - 1L]))]
      if (sgn > 0) jump_up[idx[j]] <- TRUE else jump_down[idx[j]] <- TRUE
    }
  }
  segment <- cumsum(gap | jump_up | jump_down) + 1L
  dw_all <- abs(diff(weight))
  resol <- if (any(dw_all > 0)) min(dw_all[dw_all > 0]) else 0
  outlier <- rep(FALSE, n)
  for (sg in unique(segment)) {
    idx <- which(segment == sg)
    k <- odd_window(outlier_window, length(idx) - 1L)
    if (k < 3L) next
    x <- weight[idx]
    tt <- tmin[idx]
    rate <- stats::runmed(diff(x) / diff(tt), k, endrule = "median")
    e <- x - c(0, cumsum(rate * diff(tt)))          # detrended series
    r <- e - stats::runmed(e, odd_window(outlier_window, length(e)), endrule = "median")
    s <- max(1.4826 * stats::median(abs(r - stats::median(r))), resol)
    outlier[idx] <- abs(r) > outlier_k * s
  }
  settling <- rep(FALSE, n)
  for (j in which(jump_up | jump_down))
    settling <- settling | (tmin >= tmin[j] & tmin < tmin[j] + settle_min)
  data.frame(gap_before = gap, outlier = outlier & !settling, irrigation = jump_up,
             disturbance = jump_down, settling = settling & !(jump_up | jump_down),
             segment = segment, use = !outlier & !settling & !jump_up & !jump_down)
}

qc_weights <- function(w, ...) {
  w <- w[order(w$plant_id, w$time), ]
  parts <- lapply(split(w, w$plant_id), function(p) cbind(p, qc_pot(p$time, p$weight_g, ...)))
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

# Water loss rate over a window: minus the common slope of weight on time,
# fitted to the retained readings with one intercept per segment, so an
# irrigation inside the window shifts the level but not the slope.
window_rate <- function(time, weight, segment) {
  tmin <- as.numeric(difftime(time, min(time), units = "mins"))
  seg <- factor(segment)
  f <- if (nlevels(seg) > 1) stats::lm(weight ~ tmin + seg) else stats::lm(weight ~ tmin)
  co <- suppressWarnings(summary(f))$coefficients
  if (!"tmin" %in% rownames(co)) return(c(rate_g_min = NA, se_g_min = NA))
  c(rate_g_min = -co["tmin", 1], se_g_min = co["tmin", 2])
}

# --- Step windows and per-step TR and VPD ---------------------------------------

make_windows <- function(steps, times, bin_min) {
  if (!is.null(steps)) {
    need_cols(steps, c("step", "start", "end"), "steps file")
    w <- data.frame(step = as.character(steps$step),
                    start = parse_time(steps$start, "steps start"),
                    end = parse_time(steps$end, "steps end"),
                    phase = if ("phase" %in% names(steps)) tolower(trimws(steps$phase)) else "day",
                    chamber = if ("chamber" %in% names(steps)) as.character(steps$chamber) else NA_character_,
                    stringsAsFactors = FALSE)
    if (any(w$end <= w$start))
      stop("steps file: 'end' must be after 'start' (step ", w$step[which(w$end <= w$start)[1]],
           ").", call. = FALSE)
    if (any(!w$phase %in% c("day", "night")))
      stop("steps file: 'phase' must be 'day' or 'night'.", call. = FALSE)
    if (anyDuplicated(w[, c("step", "chamber")]))
      stop("steps file: each step label must be unique (per chamber).", call. = FALSE)
    return(w)
  }
  if (bin_min <= 0) stop("--bin-min must be positive.", call. = FALSE)
  t0 <- as.POSIXct(floor(as.numeric(min(times)) / (bin_min * 60)) * bin_min * 60,
                   origin = "1970-01-01", tz = "UTC")
  starts <- seq(t0, max(times), by = bin_min * 60)
  data.frame(step = format(starts, "%Y-%m-%d %H:%M"), start = starts,
             end = starts + bin_min * 60, phase = "day", chamber = NA_character_,
             stringsAsFactors = FALSE)
}

# Mean VPD, temperature and RH per window (and chamber). With several sensors
# per chamber, sensor means are averaged; between-sensor SD is reported.
window_env <- function(env, windows, discard_first_min) {
  chambers <- if ("chamber" %in% names(env)) sort(unique(env$chamber)) else NA_character_
  rows <- list()
  for (i in seq_len(nrow(windows))) {
    chs <- if (is.na(windows$chamber[i])) chambers else windows$chamber[i]
    for (ch in chs) {
      s <- windows$start[i] + discard_first_min * 60
      sel <- env$time >= s & env$time < windows$end[i]
      if (!is.na(ch)) sel <- sel & env$chamber == ch
      e <- env[sel, ]
      if (!nrow(e)) {
        rows[[length(rows) + 1L]] <- data.frame(step = windows$step[i], chamber = ch,
          vpd_kpa = NA_real_, vpd_sd_kpa = NA_real_, sensor_sd_kpa = NA_real_,
          air_temp_c = NA_real_, rh_pct = NA_real_, n_env = 0L, stringsAsFactors = FALSE)
        next
      }
      sensor <- if ("sensor" %in% names(e)) as.character(e$sensor) else rep("1", nrow(e))
      sm <- tapply(e$vpd_kpa, sensor, mean)
      rows[[length(rows) + 1L]] <- data.frame(step = windows$step[i], chamber = ch,
        vpd_kpa = mean(sm), vpd_sd_kpa = if (nrow(e) > 1) stats::sd(e$vpd_kpa) else NA_real_,
        sensor_sd_kpa = if (length(sm) > 1) stats::sd(sm) else NA_real_,
        air_temp_c = mean(tapply(e$air_temp_c, sensor, mean)),
        rh_pct = mean(tapply(e$rh_pct, sensor, mean)), n_env = nrow(e),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

# TR per pot and window. Leaf-area-normalised TR in mg H2O m-2 s-1:
#   TR = rate (g min-1) x 1000 / 60 / (leaf area cm2 / 10000)
step_tr <- function(wq, plants, windows, discard_first_min = 0, min_readings = 10) {
  rows <- list()
  by_pot <- split(wq, wq$plant_id)
  for (p in names(by_pot)) {
    pw <- by_pot[[p]]
    info <- plants[plants$plant_id == p, , drop = FALSE]
    ch <- if ("chamber" %in% names(plants)) as.character(info$chamber[1]) else NA_character_
    win <- windows[is.na(windows$chamber) | (!is.na(ch) & windows$chamber == ch), , drop = FALSE]
    la <- if ("leaf_area_cm2" %in% names(plants)) info$leaf_area_cm2[1] else NA_real_
    for (i in seq_len(nrow(win))) {
      s <- win$start[i] + discard_first_min * 60
      inwin <- pw$time >= s & pw$time < win$end[i]
      q <- pw[inwin & pw$use, , drop = FALSE]
      rate <- c(rate_g_min = NA_real_, se_g_min = NA_real_)
      note <- ""
      if (nrow(q) >= min_readings) rate <- window_rate(q$time, q$weight_g, q$segment)
      else note <- sprintf("only %d retained readings (< %d)", nrow(q), min_readings)
      rows[[length(rows) + 1L]] <- data.frame(
        plant_id = p, step = win$step[i], phase = win$phase[i], chamber = ch,
        start = format(s, "%Y-%m-%d %H:%M:%S"),
        end = format(win$end[i], "%Y-%m-%d %H:%M:%S"),
        n_readings = sum(inwin), n_excluded = sum(inwin & !pw$use),
        tr_g_h = 60 * unname(rate["rate_g_min"]), tr_se_g_h = 60 * unname(rate["se_g_min"]),
        leaf_area_cm2 = la,
        tr_mg_m2_s = unname(rate["rate_g_min"]) * 1000 / 60 / (la / 10000),
        tr_se_mg_m2_s = unname(rate["se_g_min"]) * 1000 / 60 / (la / 10000),
        note = note, stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  out$tr_mmol_m2_s <- out$tr_mg_m2_s / MW_WATER_MG_MMOL
  out
}

# --- TR-VPD model comparison -----------------------------------------------------

# Small-sample AIC (Hurvich & Tsai 1989) from the residual sum of squares. k
# counts every estimated parameter including the residual variance: 3 for the
# line (intercept, slope, variance) and 5 for the two-segment model (intercept,
# slope 1, slope change, breakpoint, variance), as in GraphPad Prism.
aicc_rss <- function(rss, n, k) {
  if (n - k - 1 <= 0) return(NA_real_)
  n * log(rss / n) + 2 * k + 2 * k * (k + 1) / (n - k - 1)
}

# segmented's bootstrap restarts call set.seed(); keep the caller's random-number
# stream unchanged so repeated calls (and simulations) stay independent.
with_preserved_rng <- function(expr) {
  env <- globalenv()
  had <- exists(".Random.seed", envir = env, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = env, inherits = FALSE) else NULL
  on.exit(if (had) assign(".Random.seed", old, envir = env)
          else if (exists(".Random.seed", envir = env, inherits = FALSE)) rm(".Random.seed", envir = env))
  expr
}

empty_fit_row <- function(n) {
  data.frame(n = n, model = "not fitted", vpd_min_kpa = NA_real_, vpd_max_kpa = NA_real_,
             lin_intercept = NA_real_, lin_slope = NA_real_, lin_slope_se = NA_real_,
             lin_r2 = NA_real_, aicc_linear = NA_real_,
             seg_intercept = NA_real_, slope1 = NA_real_, slope1_se = NA_real_,
             slope1_lo = NA_real_, slope1_hi = NA_real_,
             breakpoint_kpa = NA_real_, breakpoint_se = NA_real_,
             breakpoint_lo = NA_real_, breakpoint_hi = NA_real_,
             slope2 = NA_real_, slope2_se = NA_real_, slope2_lo = NA_real_,
             slope2_hi = NA_real_, slope_reduction_pct = NA_real_, seg_r2 = NA_real_,
             aicc_segmented = NA_real_, delta_aicc = NA_real_,
             weight_segmented = NA_real_, davies_p = NA_real_,
             n_below_bp = NA_integer_, n_above_bp = NA_integer_, note = "",
             stringsAsFactors = FALSE)
}

# Compare TR = a + b VPD with the two-segment model
#   TR = a + b1 VPD + (b2 - b1) (VPD - psi)+    (segmented; Muggeo 2003)
# The breakpoint model is selected when its AICc is lower and, with selection
# "aicc+davies", the Davies (2002) test also rejects a constant slope at alpha.
fit_tr_vpd <- function(vpd, tr, selection = "aicc+davies", alpha = 0.05,
                       min_points = 10, min_side = 3, ci_method = "score", seed = 1) {
  if (!selection %in% c("aicc+davies", "aicc"))
    stop("--selection must be 'aicc+davies' or 'aicc'.", call. = FALSE)
  if (!ci_method %in% c("delta", "score"))
    stop("--ci-method must be 'delta' or 'score'.", call. = FALSE)
  ok <- is.finite(vpd) & is.finite(tr)
  d <- data.frame(vpd = vpd[ok], tr = tr[ok])
  n <- nrow(d)
  row <- empty_fit_row(n)
  notes <- character()
  if (sum(!ok)) notes <- c(notes, sprintf("%d point(s) without TR or VPD excluded", sum(!ok)))
  if (n < 3) {
    row$note <- paste(c(notes, "fewer than 3 points"), collapse = "; ")
    return(list(row = row, data = d, lin = NULL, seg = NULL))
  }
  row$vpd_min_kpa <- min(d$vpd)
  row$vpd_max_kpa <- max(d$vpd)
  lin <- stats::lm(tr ~ vpd, data = d)
  tss <- sum((d$tr - mean(d$tr))^2)
  rss_l <- sum(stats::residuals(lin)^2)
  cl <- suppressWarnings(summary(lin))$coefficients
  row$lin_intercept <- cl[1, 1]
  row$lin_slope <- cl[2, 1]
  row$lin_slope_se <- if (n > 2) cl[2, 2] else NA_real_
  row$lin_r2 <- 1 - rss_l / tss
  row$aicc_linear <- aicc_rss(rss_l, n, 3)
  row$model <- "linear"

  seg <- NULL
  if (n < min_points) {
    notes <- c(notes, sprintf("breakpoint model not tried: %d points (< %d)", n, min_points))
  } else if (length(unique(signif(d$vpd, 3))) < 4) {
    notes <- c(notes, "breakpoint model not tried: fewer than 4 distinct VPD levels")
  } else {
    seg <- with_preserved_rng(tryCatch(suppressWarnings(segmented::segmented(
      lin, seg.Z = ~vpd, npsi = 1, control = segmented::seg.control(seed = seed))),
      error = function(e) NULL))
    if (is.null(seg) || !inherits(seg, "segmented") || is.null(seg$psi)) {
      seg <- NULL
      notes <- c(notes, "segmented fit did not converge")
    }
  }
  if (!is.null(seg)) {
    psi <- unname(seg$psi[1, "Est."])
    rss_s <- sum(stats::residuals(seg)^2)
    sl <- tryCatch(segmented::slope(seg)$vpd, error = function(e) NULL)
    ci <- tryCatch(suppressWarnings(stats::confint(seg, method = ci_method)),
                   error = function(e) NULL)
    row$seg_intercept <- unname(stats::coef(seg)[["(Intercept)"]])
    row$breakpoint_kpa <- psi
    row$breakpoint_se <- unname(seg$psi[1, "St.Err"])
    if (!is.null(ci)) {
      row$breakpoint_lo <- ci[1, 2]
      row$breakpoint_hi <- ci[1, 3]
    }
    if (!is.null(sl)) {
      row$slope1 <- sl[1, 1]; row$slope1_se <- sl[1, 2]
      row$slope1_lo <- sl[1, 4]; row$slope1_hi <- sl[1, 5]
      row$slope2 <- sl[2, 1]; row$slope2_se <- sl[2, 2]
      row$slope2_lo <- sl[2, 4]; row$slope2_hi <- sl[2, 5]
      row$slope_reduction_pct <- 100 * (sl[1, 1] - sl[2, 1]) / sl[1, 1]
    }
    row$seg_r2 <- 1 - rss_s / tss
    row$aicc_segmented <- aicc_rss(rss_s, n, 5)
    row$n_below_bp <- sum(d$vpd <= psi)
    row$n_above_bp <- sum(d$vpd > psi)
    row$davies_p <- tryCatch(segmented::davies.test(lin, seg.Z = ~vpd)$p.value,
                             error = function(e) NA_real_)
    a <- c(row$aicc_linear, row$aicc_segmented)
    if (all(is.finite(a))) {
      row$delta_aicc <- a[1] - a[2]
      w <- exp(-0.5 * (a - min(a)))
      row$weight_segmented <- w[2] / sum(w)
    }
    valid <- row$n_below_bp >= min_side && row$n_above_bp >= min_side &&
      is.finite(row$slope1) && is.finite(row$slope2)
    if (!valid)
      notes <- c(notes, sprintf("breakpoint at the edge of the data (%s below, %s above; need %d each side)",
                                row$n_below_bp, row$n_above_bp, min_side))
    better <- is.finite(row$delta_aicc) && row$delta_aicc > 0
    davies_ok <- is.finite(row$davies_p) && row$davies_p < alpha
    if (valid && better && !davies_ok)
      notes <- c(notes, sprintf("AICc prefers a breakpoint but Davies p = %.3g", row$davies_p))
    if (valid && better && (selection == "aicc" || davies_ok)) {
      row$model <- if (row$slope2 < row$slope1) "segmented" else "segmented (slope increases)"
      if (is.finite(row$breakpoint_lo) &&
          (row$breakpoint_lo < row$vpd_min_kpa || row$breakpoint_hi > row$vpd_max_kpa))
        notes <- c(notes, "breakpoint CI extends beyond the measured VPD range")
    }
  }
  row$note <- paste(notes, collapse = "; ")
  list(row = row, data = d, lin = lin, seg = seg)
}

fit_groups <- function(dat, group_col, x_col, y_col, ...) {
  groups <- sort(unique(as.character(dat[[group_col]])))
  fits <- lapply(groups, function(g) {
    s <- dat[as.character(dat[[group_col]]) == g, , drop = FALSE]
    f <- fit_tr_vpd(s[[x_col]], s[[y_col]], ...)
    f$row <- cbind(data.frame(group = g, stringsAsFactors = FALSE), f$row)
    f
  })
  names(fits) <- groups
  tab <- do.call(rbind, lapply(fits, `[[`, "row"))
  names(tab)[1] <- group_col
  rownames(tab) <- NULL
  list(table = tab, fits = fits)
}

# --- Figures -------------------------------------------------------------------

plot_fits <- function(fits, out_dir, ylab, per_page = 9) {
  groups <- names(fits)
  pages <- split(groups, ceiling(seq_along(groups) / per_page))
  files <- character()
  for (pg in seq_along(pages)) {
    gs <- pages[[pg]]
    nc <- min(3L, length(gs)); nr <- ceiling(length(gs) / nc)
    name <- sprintf("tr_vpd_fits_%d.png", pg)
    open_png(out_dir, name, width = 3.2 * nc + 0.5, height = 3 * nr + 0.5)
    op <- graphics::par(mfrow = c(nr, nc), mar = c(4, 4.2, 3.2, 1))
    for (g in gs) {
      f <- fits[[g]]
      if (nrow(f$data) == 0) { graphics::plot.new(); graphics::title(g); next }
      graphics::plot(f$data$vpd, f$data$tr, pch = 16, col = grDevices::adjustcolor("grey30", 0.6),
                     xlab = "VPD (kPa)", ylab = ylab)
      graphics::title(main = g, line = 1.6, cex.main = 0.95)
      if (!is.null(f$lin)) graphics::abline(f$lin, lty = 2, col = "grey40")
      r <- f$row
      if (!is.null(f$seg) && grepl("^segmented", r$model)) {
        if (is.finite(r$breakpoint_lo)) {
          u <- graphics::par("usr")
          graphics::rect(r$breakpoint_lo, u[3], r$breakpoint_hi, u[4], border = NA,
                         col = grDevices::adjustcolor("#2b8cbe", 0.15))
        }
        xs <- seq(min(f$data$vpd), max(f$data$vpd), length.out = 200)
        graphics::lines(xs, stats::predict(f$seg, newdata = data.frame(vpd = xs)),
                        col = "#08589e", lwd = 2)
        graphics::abline(v = r$breakpoint_kpa, col = "#08589e", lty = 3)
        graphics::mtext(sprintf("BP %.2f kPa, -%.0f%% slope", r$breakpoint_kpa,
                                r$slope_reduction_pct), side = 3, line = 0.1, cex = 0.7)
      } else {
        graphics::mtext("linear", side = 3, line = 0.1, cex = 0.7)
      }
    }
    graphics::par(op)
    grDevices::dev.off()
    files <- c(files, name)
  }
  files
}

plot_weight_qc <- function(wq, windows, out_dir, max_pots = 12) {
  bad <- tapply(!wq$use, wq$plant_id, sum)
  pots <- names(sort(bad, decreasing = TRUE))[seq_len(min(max_pots, length(bad)))]
  pots <- sort(pots)
  nc <- min(4L, length(pots)); nr <- ceiling(length(pots) / nc)
  open_png(out_dir, "weight_qc.png", width = 3 * nc + 0.5, height = 2.6 * nr + 0.5)
  op <- graphics::par(mfrow = c(nr, nc), mar = c(3.5, 4, 2, 0.5), mgp = c(2.2, 0.7, 0))
  for (p in pots) {
    s <- wq[wq$plant_id == p, ]
    graphics::plot(s$time, s$weight_g, type = "l", col = "grey50", xlab = "time",
                   ylab = "weight (g)", main = p, cex.main = 0.9)
    graphics::abline(v = windows$start, col = "grey85", lty = 3)
    for (fl in list(list("settling", "#fdae61", 16, 0.5), list("outlier", "#d7301f", 4, 0.8),
                    list("irrigation", "#2c7bb6", 17, 1.1), list("disturbance", "#7b3294", 25, 1.1))) {
      x <- s[[fl[[1]]]]
      if (any(x)) graphics::points(s$time[x], s$weight_g[x], col = fl[[2]], pch = fl[[3]], cex = fl[[4]])
    }
  }
  graphics::par(op)
  grDevices::dev.off()
  "weight_qc.png"
}

# --- Main -------------------------------------------------------------------------

fmt <- function(x, d = 2) ifelse(is.finite(x), formatC(x, format = "f", digits = d), "NA")

prepare_raw <- function(opt, warn) {
  w <- read_csv_file(opt$weights, "weights")
  need_cols(w, c(opt$plant_col, "timestamp", "weight_g"), "weights file")
  wt <- suppressWarnings(as.numeric(as.character(w$weight_g)))
  if (all(is.na(wt))) stop("weights file: no numeric values in weight_g.", call. = FALSE)
  if (anyNA(wt))
    warn(sprintf("%d weight reading(s) were missing or not numeric (logger errors) and were dropped before analysis.", sum(is.na(wt))))
  w <- data.frame(plant_id = as.character(w[[opt$plant_col]]),
                  time = parse_time(w$timestamp, "weights"), weight_g = wt,
                  stringsAsFactors = FALSE)[!is.na(wt), ]
  dup <- duplicated(w[, c("plant_id", "time")])
  if (any(dup))
    stop("weights file: ", sum(dup), " duplicated plant/timestamp pair(s), e.g. ",
         w$plant_id[dup][1], " at ", format(w$time[dup][1]), ".", call. = FALSE)

  plants <- read_csv_file(opt$plants, "plants")
  need_cols(plants, c(opt$plant_col, if (opt$level == "genotype") opt$group), "plants file")
  names(plants)[names(plants) == opt$plant_col] <- "plant_id"
  plants$plant_id <- as.character(plants$plant_id)
  if (!opt$group %in% names(plants)) plants[[opt$group]] <- plants$plant_id
  if (anyDuplicated(plants$plant_id)) stop("plants file: plant_id values must be unique.", call. = FALSE)
  miss <- setdiff(unique(w$plant_id), plants$plant_id)
  if (length(miss))
    stop("plants file has no row for pot(s): ", paste(utils::head(miss, 10), collapse = ", "),
         ".", call. = FALSE)
  if ("leaf_area_cm2" %in% names(plants)) {
    need_numeric(plants, "leaf_area_cm2", "plants file")
    if (any(plants$leaf_area_cm2 <= 0)) stop("plants file: leaf_area_cm2 must be > 0.", call. = FALSE)
  } else if (opt$level != "plant") {
    stop("plants file needs leaf_area_cm2 to normalise TR before pooling plants. ",
         "Add it, or use --level plant to fit each pot on whole-plant TR (g/h).", call. = FALSE)
  }

  env <- read_csv_file(opt$env, "env")
  need_cols(env, c("timestamp", "air_temp_c", "rh_pct"), "env file")
  env$air_temp_c <- suppressWarnings(as.numeric(as.character(env$air_temp_c)))
  env$rh_pct <- suppressWarnings(as.numeric(as.character(env$rh_pct)))
  bad_env <- is.na(env$air_temp_c) | is.na(env$rh_pct)
  if (all(bad_env)) stop("env file: no numeric temperature/RH readings.", call. = FALSE)
  if (any(bad_env)) {
    warn(sprintf("%d environment reading(s) with missing or non-numeric temperature/RH were dropped.", sum(bad_env)))
    env <- env[!bad_env, ]
  }
  check_temp_rh(env$air_temp_c, env$rh_pct, "env file")
  env$time <- parse_time(env$timestamp, "env") + opt$env_offset_min * 60
  n_over <- sum(env$rh_pct > 100)
  if (n_over) warn(sprintf("%d RH reading(s) between 100 and 103 %% were treated as 100 %% (VPD = 0).", n_over))
  env$vpd_kpa <- vpd_kpa(env$air_temp_c, env$rh_pct, opt$svp)
  if ("chamber" %in% names(env)) {
    env$chamber <- as.character(env$chamber)
    if (!"chamber" %in% names(plants))
      stop("env file has a 'chamber' column, so the plants file needs one too, saying ",
           "which chamber each pot was in.", call. = FALSE)
    plants$chamber <- as.character(plants$chamber)
    bad <- setdiff(unique(plants$chamber), unique(env$chamber))
    if (length(bad)) stop("No env readings for chamber(s): ", paste(bad, collapse = ", "), call. = FALSE)
  } else if ("chamber" %in% names(plants)) {
    plants$chamber <- NULL
    warn("plants file has a 'chamber' column but the env file does not: one environment was used for all pots.")
  }
  if (opt$env_offset_min != 0)
    warn(sprintf("Environment timestamps were shifted by %+g min (--env-offset-min).", opt$env_offset_min))

  steps <- if (nzchar(opt$steps)) read_csv_file(opt$steps, "steps") else NULL
  if (is.null(steps))
    warn(sprintf("No --steps file: readings were binned into %g-min windows and all were treated as daytime. Remove dark periods and VPD transitions, or supply a steps file.", opt$bin_min))
  windows <- make_windows(steps, w$time, opt$bin_min)
  if (!"chamber" %in% names(plants) && any(!is.na(windows$chamber)))
    stop("steps file has a 'chamber' column but the pots have no chamber.", call. = FALSE)
  list(w = w, plants = plants, env = env, windows = windows)
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (opt$demo) {
    dd <- file.path(SKILL_DIR, "demo")
    opt$weights <- file.path(dd, "weights.csv"); opt$plants <- file.path(dd, "plants.csv")
    opt$env <- file.path(dd, "env.csv"); opt$steps <- file.path(dd, "steps.csv")
  }
  raw_mode <- nzchar(opt$weights)
  if (!raw_mode && !nzchar(opt$input))
    stop("Provide --weights/--plants/--env (raw logger data), --input (precomputed TR and VPD), or --demo. See --help.",
         call. = FALSE)
  if (raw_mode && nzchar(opt$input)) stop("Use either --weights ... or --input, not both.", call. = FALSE)
  if (!opt$level %in% c("genotype", "plant")) stop("--level must be 'genotype' or 'plant'.", call. = FALSE)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  set.seed(opt$seed)
  svp_kpa(20, opt$svp)  # validates --svp early

  warnings_list <- character()
  warn <- function(msg) warnings_list <<- c(warnings_list, msg)
  inputs <- c(opt$input, opt$weights, opt$plants, opt$env, opt$steps)
  method_lines <- character()
  qc_lines <- character()
  night <- NULL

  if (raw_mode) {
    r <- prepare_raw(opt, warn)
    wq <- qc_weights(r$w, irrigation_g = opt$irrigation_g, settle_min = opt$settle_min,
                     outlier_k = opt$outlier_k, outlier_window = opt$outlier_window,
                     gap_factor = opt$gap_factor)
    flags <- wq[!wq$use | wq$gap_before, ]
    flags$timestamp <- format(flags$time, "%Y-%m-%d %H:%M:%S")
    save_table(flags[, c("plant_id", "timestamp", "weight_g", "gap_before", "outlier",
                         "irrigation", "disturbance", "settling", "segment")],
               opt$out, "weight_flags.csv")
    st <- step_tr(wq, r$plants, r$windows, opt$discard_first_min, opt$min_readings)
    ev <- window_env(r$env, r$windows, opt$discard_first_min)
    st$chamber_key <- ifelse(is.na(st$chamber), "", st$chamber)
    ev$chamber_key <- ifelse(is.na(ev$chamber), "", ev$chamber)
    st <- merge(st, ev[, setdiff(names(ev), "chamber")], by = c("step", "chamber_key"),
                all.x = TRUE, sort = FALSE)
    st$chamber_key <- NULL
    st[[opt$group]] <- r$plants[[opt$group]][match(st$plant_id, r$plants$plant_id)]
    st <- st[order(st$plant_id, st$start), ]
    lead <- c("plant_id", opt$group)
    st <- st[, c(lead, setdiff(names(st), lead))]
    save_table(st, opt$out, "step_tr_vpd.csv")

    n_out <- sum(wq$outlier); n_irr <- sum(wq$irrigation); n_dis <- sum(wq$disturbance)
    n_set <- sum(wq$settling); n_gap <- sum(wq$gap_before)
    qc_lines <- c(
      sprintf("- Readings: %d from %d pots; excluded %d (%.2f%%).", nrow(wq), length(unique(wq$plant_id)),
              sum(!wq$use), 100 * mean(!wq$use)),
      sprintf("- Spikes (Hampel, k = %g, window = %g readings): %d; irrigation rises > %g g: %d; disturbance drops > %g g: %d; readings in the %g-min settling period after a jump: %d; gaps > %g x median interval: %d.",
              opt$outlier_k, opt$outlier_window, n_out, opt$irrigation_g, n_irr, opt$irrigation_g,
              n_dis, opt$settle_min, n_set, opt$gap_factor, n_gap))
    if (n_irr + n_dis > 0) {
      ev_tab <- wq[wq$irrigation | wq$disturbance, ]
      qc_lines <- c(qc_lines, sprintf("- Weight jump: pot %s at %s (%s).", ev_tab$plant_id,
                                      format(ev_tab$time, "%Y-%m-%d %H:%M"),
                                      ifelse(ev_tab$irrigation, "rise, irrigation", "drop, disturbance")))
    }
    med_loss <- stats::median(-diff(wq$weight_g[wq$use]), na.rm = TRUE)
    if (is.finite(med_loss) && med_loss > 0 && opt$irrigation_g < 5 * med_loss)
      warn(sprintf("--irrigation-g (%g g) is less than 5x the median loss per logging interval (%.2f g): normal transpiration may be flagged as jumps. Raise it.", opt$irrigation_g, med_loss))
    if (mean(!wq$use) > 0.2)
      warn(sprintf("%.0f%% of weight readings were excluded by quality control: inspect figures/weight_qc.png before trusting the TR values.", 100 * mean(!wq$use)))
    nfail <- sum(!is.finite(st$tr_g_h))
    if (nfail) warn(sprintf("%d pot x step combination(s) had too few retained readings for a TR estimate (see the note column of tables/step_tr_vpd.csv).", nfail))
    neg <- sum(st$phase == "day" & is.finite(st$tr_g_h) & st$tr_g_h <= 0)
    if (neg) warn(sprintf("%d daytime pot x step TR value(s) are zero or negative (the pot gained weight): check for condensation, irrigation or unsealed pots.", neg))
    unstable <- sum(st$phase == "day" & is.finite(st$vpd_sd_kpa) & st$vpd_sd_kpa > 0.2)
    if (unstable) warn(sprintf("%d pot x step combination(s) had a within-step VPD SD above 0.2 kPa: the VPD was not stable during those steps.", unstable))
    sens <- sum(is.finite(ev$sensor_sd_kpa) & ev$sensor_sd_kpa > 0.3)
    if (sens) warn(sprintf("%d step(s) had sensors in the same chamber disagreeing by SD > 0.3 kPa.", sens))
    if (any(st$phase == "night")) {
      night <- st[st$phase == "night", ]
      save_table(night, opt$out, "night_tr.csv")
    }
    dat <- st[st$phase == "day", ]
    y_col <- if (opt$level == "plant" && !"leaf_area_cm2" %in% names(r$plants)) "tr_g_h" else "tr_mg_m2_s"
    x_col <- "vpd_kpa"
    grp_col <- if (opt$level == "plant") "plant_id" else opt$group
    method_lines <- c(
      sprintf("- VPD = SVP(T) x (1 - RH/100) per environment reading (`--svp %s`), averaged over each step window%s.",
              opt$svp, if (opt$discard_first_min > 0) sprintf(" after dropping its first %g min", opt$discard_first_min) else ""),
      "- TR per pot and step = minus the slope of pot weight on time (least squares, one intercept per segment between weight jumps), from retained readings only.",
      "- TR (mg H2O m-2 s-1) = loss (g/min) x 1000 / 60 / (leaf area in m2).",
      sprintf("- Pots: %d; steps: %d day, %d night.", length(unique(st$plant_id)),
              length(unique(st$step[st$phase == "day"])), length(unique(st$step[st$phase == "night"]))))
    plot_weight_qc(wq, r$windows, opt$out)
  } else {
    dat <- read_csv_file(opt$input, "input")
    grp_col <- if (opt$level == "plant") opt$plant_col else opt$group
    need_cols(dat, c(grp_col, opt$tr_col), "input")
    need_numeric(dat, opt$tr_col, "input")
    if (!opt$vpd_col %in% names(dat)) {
      if (!all(c("air_temp_c", "rh_pct") %in% names(dat)))
        stop("input needs '", opt$vpd_col, "' or both air_temp_c and rh_pct.", call. = FALSE)
      check_temp_rh(dat$air_temp_c, dat$rh_pct, "input")
      n_over <- sum(dat$rh_pct > 100)
      if (n_over) warn(sprintf("%d RH value(s) between 100 and 103 %% were treated as 100 %%.", n_over))
      dat[[opt$vpd_col]] <- vpd_kpa(dat$air_temp_c, dat$rh_pct, opt$svp)
      method_lines <- c(method_lines, sprintf("- VPD computed from air_temp_c and rh_pct (`--svp %s`).", opt$svp))
    }
    need_numeric(dat, opt$vpd_col, "input")
    if (any(dat[[opt$vpd_col]] < 0 | dat[[opt$vpd_col]] > 8))
      stop("input: VPD outside 0-8 kPa. Is it in kPa?", call. = FALSE)
    x_col <- opt$vpd_col
    y_col <- opt$tr_col
    method_lines <- c(method_lines, sprintf("- Precomputed TR (`%s`) and VPD (`%s`) from `%s`.", y_col, x_col, opt$input))
  }

  ylab <- if (y_col == "tr_mg_m2_s") expression(TR ~ (mg ~ m^-2 ~ s^-1)) else y_col
  fg <- fit_groups(dat, grp_col, x_col, y_col, selection = opt$selection, alpha = opt$alpha,
                   min_points = opt$min_points, min_side = opt$min_side,
                   ci_method = opt$ci_method, seed = opt$seed)
  tab <- fg$table
  save_table(tab, opt$out, "tr_vpd_fits.csv")
  figs <- plot_fits(fg$fits, opt$out, ylab)

  # Red flags
  rng <- tab$vpd_max_kpa - tab$vpd_min_kpa
  if (any(is.finite(rng) & rng < 1.5))
    warn(sprintf("%d group(s) span less than 1.5 kPa of VPD: too narrow to locate a breakpoint reliably.", sum(is.finite(rng) & rng < 1.5)))
  inc <- tab[[grp_col]][tab$model == "segmented (slope increases)"]
  if (length(inc)) warn(sprintf("TR rises faster above the breakpoint in: %s. This is not a limited-transpiration response; check stomatal or temperature artefacts.", paste(inc, collapse = ", ")))
  wide <- tab[[grp_col]][grepl("^segmented", tab$model) & is.finite(tab$breakpoint_hi) & (tab$breakpoint_hi - tab$breakpoint_lo) > 1.5]
  if (length(wide)) warn(sprintf("Breakpoint 95%% CI wider than 1.5 kPa in: %s. Treat those thresholds as poorly determined.", paste(wide, collapse = ", ")))
  disagree <- tab[[grp_col]][grepl("Davies p", tab$note)]
  if (length(disagree)) warn(sprintf("AICc and the Davies test disagree for: %s (classified as %s).", paste(disagree, collapse = ", "), if (opt$selection == "aicc") "segmented" else "linear"))
  neg1 <- tab[[grp_col]][is.finite(tab$lin_slope) & tab$lin_slope <= 0]
  if (length(neg1)) warn(sprintf("TR does not increase with VPD in: %s. Check units, stomatal closure or water stress.", paste(neg1, collapse = ", ")))
  if (opt$level == "genotype" && raw_mode)
    warn("Points from the same pot are pooled as if independent (as in the published workflow): CIs are somewhat too narrow. Compare with --level plant when each pot has enough steps.")

  # Report
  res_tab <- c("| group | n | model | breakpoint kPa [95% CI] | slope 1 | slope 2 | slope reduction % | dAICc (lin - seg) | w(seg) | Davies p |",
               "|---|---|---|---|---|---|---|---|---|---|",
               sprintf("| %s | %d | %s | %s | %s | %s | %s | %s | %s | %s |",
                       tab[[grp_col]], tab$n, tab$model,
                       ifelse(grepl("^segmented", tab$model),
                              sprintf("%s [%s, %s]", fmt(tab$breakpoint_kpa), fmt(tab$breakpoint_lo), fmt(tab$breakpoint_hi)), "-"),
                       ifelse(grepl("^segmented", tab$model), fmt(tab$slope1, 1), fmt(tab$lin_slope, 1)),
                       ifelse(grepl("^segmented", tab$model), fmt(tab$slope2, 1), "-"),
                       ifelse(grepl("^segmented", tab$model), fmt(tab$slope_reduction_pct, 0), "-"),
                       fmt(tab$delta_aicc, 2), fmt(tab$weight_segmented, 2),
                       ifelse(is.finite(tab$davies_p), formatC(tab$davies_p, format = "g", digits = 3), "NA")))
  notes <- tab[nzchar(tab$note), ]
  night_lines <- NULL
  if (!is.null(night) && nrow(night)) {
    ns <- stats::aggregate(night$tr_mg_m2_s, list(g = night[[opt$group]]), function(x)
      c(mean = mean(x, na.rm = TRUE), se = stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))), n = sum(is.finite(x))))
    night_lines <- c("## Night-time TR", "",
                     "Night steps are excluded from the VPD-response fits. Mean over pot x night-step values:", "",
                     "| group | night TR (mg m-2 s-1) | SE | n |", "|---|---|---|---|",
                     sprintf("| %s | %s | %s | %d |", ns$g, fmt(ns$x[, "mean"], 2), fmt(ns$x[, "se"], 2), as.integer(ns$x[, "n"])), "")
  }
  report <- c(sprintf("# %s report", SKILL), "",
              if (length(warnings_list)) c("## Warnings", "", paste("-", warnings_list), "") else NULL,
              "## What was done", "",
              if (raw_mode) sprintf("Raw logger data: weights `%s`, plants `%s`, environment `%s`%s.", opt$weights, opt$plants, opt$env,
                                    if (nzchar(opt$steps)) sprintf(", steps `%s`", opt$steps) else "") else sprintf("Input: `%s`.", opt$input),
              "", method_lines,
              sprintf("- Fits per %s: linear vs two-segment regression (R package segmented, seed %g). Model choice: %s (alpha = %g); breakpoint CI by the %s method; breakpoint models need >= %g points and >= %g on each side.",
                      if (opt$level == "plant") "pot" else sprintf("`%s` (plants pooled)", grp_col), opt$seed,
                      if (opt$selection == "aicc") "lower AICc" else "lower AICc and Davies test p < alpha",
                      opt$alpha, opt$ci_method, opt$min_points, opt$min_side),
              "",
              if (length(qc_lines)) c("## Weight quality control", "", qc_lines, "",
                                      "Figure: pot weights (up to 12 pots with the most excluded readings); red x = spike, orange = settling after a jump, blue triangle = irrigation, purple = disturbance; dotted lines = step starts.", "",
                                      "![Weight QC](figures/weight_qc.png)", "") else NULL,
              "## TR-VPD response", "",
              sprintf("TR units: `%s`. Slopes in TR units per kPa. w(seg) = Akaike weight of the breakpoint model.", y_col), "",
              res_tab, "",
              if (nrow(notes)) c("Notes:", "", sprintf("- %s: %s", notes[[grp_col]], notes$note), "") else NULL,
              "Figure: dashed line = linear fit; solid line = two-segment fit (drawn only when selected); shaded band = 95% CI of the breakpoint.", "",
              sprintf("![TR-VPD fits](figures/%s)", figs), "",
              night_lines,
              "## Limitations", "",
              "- The breakpoint is a statistical summary of a curved, often gradual, stomatal response; its value depends on the VPD steps imposed and on the noise.",
              "- With few VPD steps the breakpoint can only be located between steps; delta-method CIs (--ci-method delta) assume a clear kink and are too narrow for gradual responses.",
              "- TR is per unit leaf area measured once; leaf-area error scales the whole curve of a pot.",
              "- Soil evaporation must be prevented (sealed pots); otherwise TR includes evaporation.",
              "",
              "## References", "",
              "- Tamang BG, Sadok W (2018) Environ Exp Bot 148:192-202. https://doi.org/10.1016/j.envexpbot.2017.11.016",
              "- Muggeo VMR (2003) Stat Med 22:3055-3071. https://doi.org/10.1002/sim.1545",
              "- Muggeo VMR (2017) Aust N Z J Stat 59:311-322. https://doi.org/10.1111/anzs.12200",
              "- Davies RB (2002) Biometrika 89:484-489. https://doi.org/10.1093/biomet/89.2.484",
              "- Hurvich CM, Tsai C-L (1989) Biometrika 76:297-307. https://doi.org/10.1093/biomet/76.2.297",
              "- Buck AL (1981) J Appl Meteorol 20:1527-1532. https://doi.org/10.1175/1520-0450(1981)020<1527:NEFCVP>2.0.CO;2",
              "- Sinclair TR et al. (2017) Plant Sci 260:109-118. https://doi.org/10.1016/j.plantsci.2017.04.007",
              "",
              sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(report, file.path(opt$out, "report.md"))
  write_repro_bundle(opt$out, argv, inputs = inputs)
  cat(sprintf("Done: %d group(s) fitted, %d with a breakpoint. Report: %s\n", nrow(tab),
              sum(grepl("^segmented", tab$model)),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(list(fits = tab, steps = if (raw_mode) st else NULL))
}

if (sys.nframe() == 0L) main()
