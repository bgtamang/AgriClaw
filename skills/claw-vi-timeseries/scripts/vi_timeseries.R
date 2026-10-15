#!/usr/bin/env Rscript
# claw-vi-timeseries: smooth per-plot vegetation-index (VI) time series and
# extract phenology metrics (peak, green-up, senescence, senescence rate, area
# under the curve, season length) on a days-after-planting, growing-degree-day
# or day-of-year axis.
#
# Usage:
#   Rscript scripts/vi_timeseries.R --demo
#   Rscript scripts/vi_timeseries.R --input vi_long.csv --vi ndvi --planting-date 2025-05-12
#   Rscript scripts/vi_timeseries.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-vi-timeseries"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/vi_timeseries.R [--demo] [--input FILE.csv] [options]

Input
  --input          long CSV, one row per plot x date (e.g. from claw-vi-extract)
  --vi             VI column to analyse (default: ndvi)
  --plot-col       plot identifier column (default: plot_id)
  --date-col       observation date column, YYYY-MM-DD (default: date)
  --planting-date  planting date YYYY-MM-DD for every plot; or give a
                   planting_date column in the input. Never guess it.
  --keep-cols      comma-separated per-plot columns copied to metrics.csv
                   (e.g. env,genotype,rep)
  --weight-col     optional prior weights in [0, 1] (e.g. from cloud/QA masks)
Axis
  --axis           dap (days after planting, default), gdd or doy
  --gdd            daily GDD CSV with a date column and --gdd-col (and
                   --env-col if it holds several series); needed for --axis gdd
  --gdd-col        daily GDD column in that file (default: gdd_c, as in
                   claw-gdd's tables/daily_gdd.csv; gdd_f for deg F days)
  --env-col        column that links plots to their GDD series
  --cutoff         use only observations with axis value <= this, e.g. for
                   early-season prediction features (default: all)
Fitting
  --method         dl-beck (default), dl-elmore, loess or sg
  --weighting      none (default), upper (upper envelope: down-weights points
                   below the curve, for cloud/shadow-affected satellite data)
                   or robust (bisquare: down-weights outliers on both sides)
  --min-dates      fewest dates needed to fit a curve (default 0 = auto:
                   7 dl-beck, 8 dl-elmore, 6 loess/sg)
  --span           LOESS span (default 0.4)
  --sg-order       Savitzky-Golay polynomial order (default 3)
  --sg-window      Savitzky-Golay window in axis units (default 0 = auto,
                   3 x the median sampling interval)
Metrics
  --threshold      amplitude fraction that defines green-up and senescence
                   (default 0.5)
  --min-limb-frac  report green-up/senescence only if that side of the curve
                   rises/falls by at least this fraction of the larger side
                   (default 0.5)
  --features-at    comma-separated axis values at which to report the fitted
                   VI, e.g. 40,50,60 (NA outside each plot's observed window)
Output
  --out            output folder (default: agriclaw_output/claw-vi-timeseries)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, input = "", vi = "ndvi", plot_col = "plot_id",
                 date_col = "date", planting_date = "", keep_cols = "",
                 weight_col = "", axis = "dap", gdd = "", gdd_col = "gdd_c",
                 env_col = "", cutoff = NA_real_, method = "dl-beck",
                 weighting = "none", min_dates = 0, span = 0.4, sg_order = 3,
                 sg_window = 0, threshold = 0.5, min_limb_frac = 0.5,
                 features_at = "", out = "")

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

METHODS    <- c("dl-beck", "dl-elmore", "loess", "sg")
WEIGHTINGS <- c("none", "upper", "robust")
AXES       <- c("dap", "gdd", "doy")
N_PAR      <- c("dl-beck" = 6L, "dl-elmore" = 7L)
BISQUARE_C     <- 6     # bisquare scale = 6 x median |residual| (Cleveland 1979)
MAX_ITER       <- 10L   # reweighting iterations (rule of thumb)
WEIGHT_TOL     <- 1e-3  # reweighting stops when no weight changes by more
OUTLIER_Z      <- 3.5   # flag |residual| > 3.5 noise SD (rule of thumb)
DATE_SHIFT_SD  <- 2     # flag a date whose median residual over plots > 2 noise SD (rule of thumb)
NEAR_PEAK_FRAC <- 0.05  # near-peak window: within 5% of the amplitude (rule of thumb)
GRID_N         <- 2001L # evaluation points across each plot's observed window
MIN_LOESS_PTS  <- 5     # each LOESS neighbourhood holds >= 5 points (rule of thumb)

axis_info <- function(axis) {
  switch(axis,
         dap = list(sfx = "dap", len = "days", per = "per_day",
                    label = "Days after planting"),
         doy = list(sfx = "doy", len = "days", per = "per_day",
                    label = "Day of year"),
         gdd = list(sfx = "gdd", len = "gdd", per = "per_gdd",
                    label = "Growing degree days after planting"))
}

min_dates_for <- function(method, min_dates = 0) {
  floor_n <- if (method %in% names(N_PAR)) N_PAR[[method]] + 1L else 6L
  if (min_dates <= 0) return(floor_n)
  if (method %in% names(N_PAR) && min_dates < floor_n)
    stop(sprintf("--min-dates %g is too low for %s: the model has %d parameters, so at least %d dates are needed.",
                 min_dates, method, N_PAR[[method]], floor_n), call. = FALSE)
  as.integer(min_dates)
}

split_list <- function(x) {
  x <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  x[nzchar(x)]
}

read_csv_checked <- function(path, what = "Input") {
  if (!nzchar(path) || !file.exists(path))
    stop(what, " file not found: '", path, "'.", call. = FALSE)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

parse_iso_date <- function(x, what) {
  s <- trimws(as.character(x))
  d <- as.Date(s, format = "%Y-%m-%d")
  bad <- is.na(d) & !is.na(x) & nzchar(s)
  if (any(bad))
    stop(sprintf("%s: %d value(s) are not dates in YYYY-MM-DD format, e.g. '%s'. Convert them first.",
                 what, sum(bad), s[which(bad)[1]]), call. = FALSE)
  d
}

check_vi_scale <- function(y, vi) {
  if (!grepl("(^|_)(g?ndvi|ndre|ndwi|osavi|savi|evi2?|gndvi)($|_)", tolower(vi))) return(NULL)
  if (any(abs(y) > 1.5, na.rm = TRUE))
    stop(sprintf("'%s' has values up to %g: normalised-difference indices lie in [-1, 1]. Is the column scaled (e.g. x 10000)? Rescale it before running.",
                 vi, max(abs(y), na.rm = TRUE)), call. = FALSE)
  if (any(abs(y) > 1, na.rm = TRUE))
    return(sprintf("%d value(s) of '%s' are outside [-1, 1]; check the index calculation.",
                   sum(abs(y) > 1, na.rm = TRUE), vi))
  NULL
}

# Cumulative GDD from each plot's planting date: the planting day contributes 0,
# so GDD = 0 at planting, matching DAP = 0.
gdd_lookup <- function(plots, gdd_tab, gdd_col, env_col) {
  need <- c("date", gdd_col, if (nzchar(env_col)) env_col)
  miss <- setdiff(need, names(gdd_tab))
  if (length(miss))
    stop("GDD file is missing column(s): ", paste(miss, collapse = ", "),
         ". It needs a date column and a daily GDD column (--gdd-col",
         if (nzchar(env_col)) ", plus the --env-col column" else "", ").", call. = FALSE)
  if (!is.numeric(gdd_tab[[gdd_col]]))
    stop("GDD column '", gdd_col, "' must be numeric (daily growing degree days).", call. = FALSE)
  gdd_tab$date <- parse_iso_date(gdd_tab$date, "GDD file, column 'date'")
  key <- if (nzchar(env_col)) as.character(gdd_tab[[env_col]]) else rep("all", nrow(gdd_tab))
  dup <- duplicated(paste(key, gdd_tab$date))
  if (any(dup))
    stop("GDD file has repeated dates (", sum(dup), "), e.g. ", format(gdd_tab$date[dup][1]),
         if (nzchar(env_col)) paste0(" for ", env_col, " = ", key[dup][1]) else "",
         ". Give one row per date", if (nzchar(env_col)) " per environment" else "",
         ".", call. = FALSE)
  out <- vector("list", nrow(plots))
  for (i in seq_len(nrow(plots))) {
    k <- if (nzchar(env_col)) as.character(plots$env[i]) else "all"
    days <- seq(plots$planting[i], plots$last_date[i], by = "day")
    sel <- key == k
    if (!any(sel))
      stop("No GDD series for ", env_col, " = '", k, "' (plot ", plots$plot[i], ").", call. = FALSE)
    g <- gdd_tab[[gdd_col]][sel][match(days[-1], gdd_tab$date[sel])]
    if (anyNA(g)) {
      gap <- days[-1][is.na(g)]
      stop(sprintf("GDD series%s has %d missing day(s) between planting (%s) and the last observation of plot %s, first %s. Fill the weather gaps (claw-gdd / claw-weather) before using --axis gdd.",
                   if (nzchar(env_col)) paste0(" for ", env_col, " = '", k, "'") else "",
                   length(gap), format(plots$planting[i]), plots$plot[i], format(gap[1])),
           call. = FALSE)
    }
    if (any(g < 0)) stop("Daily GDD values must be >= 0.", call. = FALSE)
    out[[i]] <- data.frame(plot = plots$plot[i], date = days, gdd = c(0, cumsum(g)),
                           stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

# Read, validate and put the input on the analysis axis.
prepare_series <- function(dat, opt, gdd_tab = NULL) {
  notes <- character()
  need <- c(opt$plot_col, opt$date_col, opt$vi,
            if (nzchar(opt$weight_col)) opt$weight_col,
            if (nzchar(opt$env_col)) opt$env_col)
  keep <- split_list(opt$keep_cols)
  need <- unique(c(need, keep))
  miss <- setdiff(need, names(dat))
  if (length(miss))
    stop("Column(s) not found: ", paste(miss, collapse = ", "), ". Available: ",
         paste(names(dat), collapse = ", "),
         ". Use --plot-col, --date-col and --vi to name the columns.", call. = FALSE)
  if (!is.numeric(dat[[opt$vi]]))
    stop("VI column '", opt$vi, "' must be numeric.", call. = FALSE)
  if (!opt$axis %in% AXES)
    stop("--axis must be one of: ", paste(AXES, collapse = ", "), call. = FALSE)

  n_na <- sum(is.na(dat[[opt$vi]]))
  if (n_na) notes <- c(notes, sprintf("%d row(s) with missing '%s' were excluded.", n_na, opt$vi))
  dat <- dat[!is.na(dat[[opt$vi]]), , drop = FALSE]
  if (!nrow(dat)) stop("No non-missing values in '", opt$vi, "'.", call. = FALSE)
  sc <- check_vi_scale(dat[[opt$vi]], opt$vi)
  if (!is.null(sc)) notes <- c(notes, sc)

  obs <- data.frame(plot = as.character(dat[[opt$plot_col]]),
                    date = parse_iso_date(dat[[opt$date_col]], paste0("Column '", opt$date_col, "'")),
                    y = dat[[opt$vi]], stringsAsFactors = FALSE)
  if (anyNA(obs$plot) || any(!nzchar(obs$plot)))
    stop("Column '", opt$plot_col, "' has empty plot identifiers.", call. = FALSE)
  if (anyNA(obs$date)) stop("Column '", opt$date_col, "' has missing dates.", call. = FALSE)
  dup <- duplicated(obs[, c("plot", "date")])
  if (any(dup))
    stop(sprintf("%d repeated plot x date row(s), e.g. plot '%s' on %s. Keep one value per plot and date (e.g. average repeated flights) before running; the script never averages silently.",
                 sum(dup), obs$plot[dup][1], format(obs$date[dup][1])), call. = FALSE)

  if (nzchar(opt$weight_col)) {
    w <- dat[[opt$weight_col]]
    if (!is.numeric(w) || anyNA(w) || any(w < 0 | w > 1))
      stop("Weight column '", opt$weight_col, "' must be numeric in [0, 1] with no missing values.", call. = FALSE)
    obs$w0 <- w
  } else obs$w0 <- 1

  # Planting date: from the option or a planting_date column, never guessed.
  has_col <- "planting_date" %in% names(dat)
  if (nzchar(opt$planting_date) && has_col)
    stop("Give the planting date either with --planting-date or as a planting_date column, not both.", call. = FALSE)
  if (nzchar(opt$planting_date)) {
    obs$planting <- parse_iso_date(rep(opt$planting_date, nrow(obs)), "--planting-date")
  } else if (has_col) {
    obs$planting <- parse_iso_date(dat$planting_date, "Column 'planting_date'")
    if (anyNA(obs$planting)) stop("Column 'planting_date' has missing values.", call. = FALSE)
  } else obs$planting <- as.Date(NA)
  if (opt$axis %in% c("dap", "gdd") && anyNA(obs$planting))
    stop("No planting date. Ask the user for it (never guess), then pass --planting-date YYYY-MM-DD or add a planting_date column. If it is truly unknown, --axis doy analyses calendar day of year, but plots planted on different dates are then not comparable.",
         call. = FALSE)
  if (!all(is.na(obs$planting))) {
    np <- tapply(obs$planting, obs$plot, function(x) length(unique(x)))
    if (any(np > 1))
      stop("Plot '", names(np)[np > 1][1], "' has more than one planting_date.", call. = FALSE)
    early <- obs$date < obs$planting
    if (any(early))
      stop(sprintf("%d observation(s) are dated before planting, e.g. plot '%s' on %s (planted %s). Check the dates.",
                   sum(early), obs$plot[early][1], format(obs$date[early][1]),
                   format(obs$planting[early][1])), call. = FALSE)
  }
  for (k in c(keep, if (nzchar(opt$env_col)) opt$env_col)) {
    nv <- tapply(dat[[k]], obs$plot, function(x) length(unique(x)))
    if (any(nv > 1))
      stop("Column '", k, "' changes within plot '", names(nv)[nv > 1][1],
           "'; --keep-cols and --env-col must be constant within each plot.", call. = FALSE)
    obs[[paste0("keep__", k)]] <- dat[[k]]
  }
  obs$env <- if (nzchar(opt$env_col)) as.character(dat[[opt$env_col]]) else NA_character_

  # Axis value t.
  daymap <- NULL
  if (opt$axis == "dap") {
    obs$t <- as.numeric(obs$date - obs$planting)
  } else if (opt$axis == "doy") {
    yrs <- tapply(format(obs$date, "%Y"), obs$plot, function(x) length(unique(x)))
    if (any(yrs > 1))
      stop("Plot '", names(yrs)[yrs > 1][1], "' spans two calendar years; --axis doy needs one season per calendar year. Use --axis dap.", call. = FALSE)
    obs$t <- as.numeric(format(obs$date, "%j"))
    if (length(unique(obs$planting[!is.na(obs$planting)])) > 1)
      notes <- c(notes, "Plots have different planting dates but the axis is day of year: metrics are not comparable between planting dates. Use --axis dap or gdd.")
  } else {
    if (is.null(gdd_tab))
      stop("--axis gdd needs --gdd FILE.csv with daily growing degree days (for example from claw-gdd).", call. = FALSE)
    plots <- unique(obs[, c("plot", "planting", "env")])
    plots$last_date <- as.Date(tapply(obs$date, obs$plot, max)[plots$plot])
    daymap <- gdd_lookup(plots, gdd_tab, opt$gdd_col, opt$env_col)
    obs$t <- daymap$gdd[match(paste(obs$plot, obs$date), paste(daymap$plot, daymap$date))]
  }

  if (!is.na(opt$cutoff)) {
    late <- obs$t > opt$cutoff
    if (any(late))
      notes <- c(notes, sprintf("%d observation(s) after the cutoff (%s %g) were excluded (--cutoff).",
                                sum(late), axis_info(opt$axis)$sfx, opt$cutoff))
    obs <- obs[!late, , drop = FALSE]
    if (!nrow(obs)) stop("No observations remain at or before --cutoff ", opt$cutoff, ".", call. = FALSE)
  }
  obs <- obs[order(obs$plot, obs$t), , drop = FALSE]
  rownames(obs) <- NULL
  list(obs = obs, notes = notes, daymap = daymap)
}

# --- curve models (Beck et al. 2006; Elmore et al. 2012) ---------------------

# Beck: vmin + (vmax - vmin) * (1/(1+exp(-rsp (t - sos))) + 1/(1+exp(rau (t - eos))) - 1)
dl_beck <- function(t, p) {
  p[["vmin"]] + (p[["vmax"]] - p[["vmin"]]) *
    (1 / (1 + exp(-p[["rsp"]] * (t - p[["sos"]]))) +
       1 / (1 + exp(p[["rau"]] * (t - p[["eos"]]))) - 1)
}

# Elmore: base + (amp - m7 t) * (1/(1+exp(-rsp (t - sos))) - 1/(1+exp(-rau (t - eos))))
dl_elmore <- function(t, p) {
  p[["base"]] + (p[["amp"]] - p[["m7"]] * t) *
    (1 / (1 + exp(-p[["rsp"]] * (t - p[["sos"]]))) -
       1 / (1 + exp(-p[["rau"]] * (t - p[["eos"]]))))
}

# The same models on u = (t - t0) / span in [0, 1], reparameterised so that
# amp >= 0 and eos = sos + gap >= sos (keeps the curve a single hump).
beck_u <- function(u, vmin, amp, sos, rsp, gap, rau) {
  vmin + amp * (1 / (1 + exp(-rsp * (u - sos))) + 1 / (1 + exp(rau * (u - sos - gap))) - 1)
}
elmore_u <- function(u, base, amp, sos, rsp, gap, rau, m7) {
  base + (amp - m7 * u) * (1 / (1 + exp(-rsp * (u - sos))) - 1 / (1 + exp(-rau * (u - sos - gap))))
}

# Fit a double logistic with minpack.lm::nlsLM (Levenberg-Marquardt with bounds)
# from several starting values; keep the converged fit with the lowest weighted
# residual sum of squares.
fit_dl <- function(t, y, w, model = c("dl-beck", "dl-elmore"), prev = NULL) {
  model <- match.arg(model)
  t0 <- min(t)
  span <- max(t) - t0
  u <- (t - t0) / span
  rng <- diff(range(y))
  if (rng <= 0) return(NULL)
  ys <- if (length(y) >= 5) stats::runmed(y, 3, endrule = "keep") else y
  u_max <- u[which.max(ys)]
  sos0 <- 0.5 * u_max                     # halfway from the first date to the highest value
  gap0 <- (u_max + 0.5 * (1 - u_max)) - sos0  # eos0 halfway from the highest value to the last date
  # A logistic rises from 10% to 90% in log(81)/rate. Rates are capped so that
  # this transition is no shorter than the median sampling interval: faster
  # changes cannot be resolved, and a capped rate is reported as "at a bound".
  r_max <- max(10, log(81) * span / stats::median(diff(sort(t))))
  if (model == "dl-beck") {
    lower <- c(vmin = min(y) - rng, amp = 0, sos = -0.5, rsp = 1, gap = 0, rau = 1)
    upper <- c(vmin = max(y), amp = 2 * rng, sos = 1.5, rsp = r_max, gap = 2, rau = r_max)
    base_start <- c(vmin = min(y), amp = rng, sos = sos0, rsp = 10, gap = gap0, rau = 10)
    form <- y ~ beck_u(u, vmin, amp, sos, rsp, gap, rau)
  } else {
    lower <- c(base = min(y) - rng, amp = 0, sos = -0.5, rsp = 1, gap = 0, rau = 1, m7 = -rng)
    upper <- c(base = max(y), amp = 3 * rng, sos = 1.5, rsp = r_max, gap = 2, rau = r_max, m7 = rng)
    base_start <- c(base = min(y), amp = rng, sos = sos0, rsp = 10, gap = gap0, rau = 10, m7 = 0)
    form <- y ~ elmore_u(u, base, amp, sos, rsp, gap, rau, m7)
  }
  starts <- list()
  if (!is.null(prev) && !is.null(prev$par_u)) starts <- list(prev$par_u)
  for (r in list(c(8, 8), c(25, 25), c(8, 25), c(25, 8))) {
    s <- base_start
    s[["rsp"]] <- r[1]
    s[["rau"]] <- r[2]
    starts[[length(starts) + 1L]] <- s
  }
  if (!is.null(prev) && !is.null(prev$par_u)) starts <- starts[1]  # warm start while reweighting
  best <- NULL
  dfit <- data.frame(u = u, y = y)
  for (s in starts) {
    s <- pmin(pmax(s, lower), upper)
    m <- tryCatch(
      suppressWarnings(minpack.lm::nlsLM(form, data = dfit, start = s, lower = lower,
                                         upper = upper, weights = w,
                                         control = minpack.lm::nls.lm.control(maxiter = 500))),
      error = function(e) NULL)
    if (is.null(m)) next
    rss <- sum(w * stats::residuals(m)^2)
    if (is.null(best) || rss < best$rss) best <- list(m = m, rss = rss)
  }
  if (is.null(best) && !is.null(prev)) return(fit_dl(t, y, w, model, prev = NULL))
  if (is.null(best)) return(NULL)
  pu <- stats::coef(best$m)
  at_bound <- names(pu)[abs(pu - lower[names(pu)]) < 1e-6 * pmax(1, abs(lower[names(pu)])) |
                          abs(pu - upper[names(pu)]) < 1e-6 * pmax(1, abs(upper[names(pu)]))]
  if (model == "dl-beck") {
    par <- c(vmin = pu[["vmin"]], vmax = pu[["vmin"]] + pu[["amp"]],
             sos = t0 + span * pu[["sos"]], rsp = pu[["rsp"]] / span,
             eos = t0 + span * (pu[["sos"]] + pu[["gap"]]), rau = pu[["rau"]] / span)
    f <- function(tt) dl_beck(tt, par)
  } else {
    m7 <- pu[["m7"]] / span
    par <- c(base = pu[["base"]], amp = pu[["amp"]] + m7 * t0,
             sos = t0 + span * pu[["sos"]], rsp = pu[["rsp"]] / span,
             eos = t0 + span * (pu[["sos"]] + pu[["gap"]]), rau = pu[["rau"]] / span, m7 = m7)
    f <- function(tt) dl_elmore(tt, par)
  }
  list(predict = f, fitted = f(t), par = par, par_u = pu, at_bound = at_bound,
       df_resid = sum(w > 0) - length(pu))
}

# LOESS (Cleveland & Devlin 1988), local quadratic. Each neighbourhood holds at
# least MIN_LOESS_PTS points.
fit_loess <- function(t, y, w, span = 0.4, prev = NULL) {
  span_eff <- min(1, max(span, MIN_LOESS_PTS / length(t)))
  m <- tryCatch(
    suppressWarnings(stats::loess(y ~ t, data = data.frame(t = t, y = y), weights = pmax(w, 1e-6),
                                  span = span_eff, degree = 2, family = "gaussian",
                                  control = stats::loess.control(surface = "direct"))),
    error = function(e) NULL)
  if (is.null(m)) return(NULL)
  f <- function(tt) as.numeric(stats::predict(m, newdata = data.frame(t = tt)))
  list(predict = f, fitted = f(t), par = c(span = span_eff), at_bound = character(),
       df_resid = m$one.delta)
}

# Savitzky-Golay (Savitzky & Golay 1964) through signal::sgolayfilt. The filter
# needs equal spacing, so observations are first interpolated linearly onto a
# regular grid (step = median sampling interval / 5). Weights cannot enter the
# filter itself; they act by replacement, y* = w y + (1 - w) yhat (the
# replacement step of Chen et al. 2004), from the second pass on.
fit_sg <- function(t, y, w, order = 3, window = 0, prev = NULL) {
  step <- stats::median(diff(t)) / 5
  if (!is.finite(step) || step <= 0) return(NULL)
  grid <- seq(min(t), max(t), by = step)
  if (utils::tail(grid, 1) < max(t)) grid <- c(grid, max(t))
  y_adj <- if (is.null(prev)) y else w * y + (1 - w) * prev$fitted
  yi <- stats::approx(t, y_adj, xout = grid)$y
  win_axis <- if (window > 0) window else 3 * stats::median(diff(t))
  n <- 2L * floor(win_axis / step / 2) + 1L
  n_min <- order + 3L - order %% 2L
  if (n < n_min) n <- n_min
  if (n > length(grid)) n <- length(grid) - (1L - length(grid) %% 2L)
  if (n < n_min || n <= order) return(NULL)
  smooth <- function(v) {
    s <- signal::sgolayfilt(stats::approx(t, v, xout = grid)$y, p = order, n = n)
    stats::approx(grid, s, xout = t)$y
  }
  s <- signal::sgolayfilt(yi, p = order, n = n)
  f <- function(tt) stats::approx(grid, s, xout = tt)$y
  # interpolation + filter is linear in y: build its hat matrix for the
  # residual degrees of freedom n - 2 tr(H) + tr(H'H)
  H <- vapply(seq_along(t), function(j) smooth(as.numeric(seq_along(t) == j)), numeric(length(t)))
  list(predict = f, fitted = f(t), par = c(order = order, window_points = n, step = step),
       at_bound = character(), df_resid = length(t) - 2 * sum(diag(H)) + sum(H^2))
}

# Bisquare robustness weights (Cleveland 1979): zero beyond 6 x median
# |residual| of the series. "upper" keeps full weight for points on or above
# the curve (upper envelope, as in TIMESAT, Jonsson & Eklundh 2002); "robust"
# down-weights large residuals on both sides. If `sigma` (a noise SD pooled
# over plots) is given, the cut-off is the same 6 x 0.6745 x sigma instead.
robust_weights <- function(r, weighting, sigma = NA_real_) {
  s <- if (is.finite(sigma) && sigma > 0) 0.6745 * sigma else stats::median(abs(r))
  if (!is.finite(s) || s <= 0) return(rep(1, length(r)))
  z <- r / (BISQUARE_C * s)
  b <- ifelse(abs(z) < 1, (1 - z^2)^2, 0)
  if (weighting == "upper") b[r >= 0] <- 1
  b
}

# Fit one series with optional iterative reweighting.
fit_series <- function(t, y, w0, method, weighting, opt, sigma = NA_real_) {
  fitter <- switch(method,
                   "dl-beck"   = function(w, prev) fit_dl(t, y, w, "dl-beck", prev),
                   "dl-elmore" = function(w, prev) fit_dl(t, y, w, "dl-elmore", prev),
                   "loess"     = function(w, prev) fit_loess(t, y, w, opt$span, prev),
                   "sg"        = function(w, prev) fit_sg(t, y, w, opt$sg_order, opt$sg_window, prev))
  n_need <- if (method %in% names(N_PAR)) N_PAR[[method]] + 1L else MIN_LOESS_PTS
  w <- w0
  fit <- fitter(w, NULL)
  if (is.null(fit)) return(NULL)
  if (method == "sg" && any(w0 < 1)) {
    fit2 <- fitter(w0, fit)
    if (!is.null(fit2)) fit <- fit2
  }
  iter <- 0L
  note <- character()
  if (weighting != "none") {
    for (iter in seq_len(MAX_ITER)) {
      w_new <- w0 * robust_weights(y - fit$fitted, weighting, sigma)
      if (max(abs(w_new - w)) < WEIGHT_TOL) break
      if (sum(w_new > 0) < n_need) {
        note <- "reweighting stopped: too few points kept"
        break
      }
      fit_new <- fitter(w_new, fit)
      if (is.null(fit_new)) {
        note <- "reweighting stopped: refit failed"
        break
      }
      w <- w_new
      fit <- fit_new
    }
  }
  fit$weights <- w
  fit$iterations <- iter
  fit$note <- note
  fit
}

crossing <- function(t1, t2, f1, f2, thr) t1 + (thr - f1) * (t2 - t1) / (f2 - f1)

# Phenology metrics from a curve evaluated on a fine grid inside the observed
# window. Green-up / senescence: the curve crosses base + threshold x amplitude,
# with base and amplitude taken separately on each side of the peak (White et
# al. 1997 midpoint idea; threshold method of TIMESAT / phenofit).
curve_metrics <- function(tt, ff, threshold = 0.5, min_limb_frac = 0.5) {
  ok <- is.finite(ff)
  tt <- tt[ok]
  ff <- ff[ok]
  n <- length(tt)
  out <- list(peak_value = NA_real_, peak_t = NA_real_, near_peak_start_t = NA_real_,
              near_peak_end_t = NA_real_, greenup_t = NA_real_, senescence_t = NA_real_,
              season_length = NA_real_, greenup_rate = NA_real_, senescence_rate = NA_real_,
              auc = NA_real_, base_pre = NA_real_, base_post = NA_real_,
              amp_pre = NA_real_, amp_post = NA_real_, notes = character())
  if (n < 3) {
    out$notes <- "curve could not be evaluated"
    return(out)
  }
  ip <- which.max(ff)
  peak <- ff[ip]
  base_pre <- min(ff[1:ip])
  base_post <- min(ff[ip:n])
  amp_pre <- peak - base_pre
  amp_post <- peak - base_post
  amp_max <- max(amp_pre, amp_post)
  out[c("peak_value", "peak_t", "base_pre", "base_post", "amp_pre", "amp_post")] <-
    list(peak, tt[ip], base_pre, base_post, amp_pre, amp_post)
  out$auc <- sum(diff(tt) * (ff[-1] + ff[-n]) / 2)
  notes <- character()
  if (ip == 1) notes <- c(notes, "peak at the first date: green-up not observed")
  if (ip == n) notes <- c(notes, "peak at the last date: senescence not observed")
  # near-peak window (curve within NEAR_PEAK_FRAC of the amplitude of the peak)
  if (amp_max > 0) {
    thr_np <- peak - NEAR_PEAK_FRAC * amp_max
    a <- ip
    while (a > 1 && ff[a - 1] >= thr_np) a <- a - 1
    b <- ip
    while (b < n && ff[b + 1] >= thr_np) b <- b + 1
    out$near_peak_start_t <- tt[a]
    out$near_peak_end_t <- tt[b]
  }
  d <- diff(ff) / diff(tt)
  # green-up: first upward crossing after the minimum that precedes the peak
  if (ip > 1 && amp_max > 0 && amp_pre >= min_limb_frac * amp_max) {
    i0 <- which.min(ff[1:ip])
    thr <- base_pre + threshold * amp_pre
    above <- ff[i0:ip] >= thr
    j <- i0 - 1 + which(above)[1]
    out$greenup_t <- crossing(tt[j - 1], tt[j], ff[j - 1], ff[j], thr)
    nc <- sum(diff(above) == 1)
    if (nc > 1) notes <- c(notes, sprintf("rising side crosses the green-up threshold %d times (curve not monotonic); the first crossing is reported", nc))
    seg <- i0:(ip - 1)
    k <- seg[which.max(d[seg])]
    if (k > 1) out$greenup_rate <- d[k] else
      notes <- c(notes, "steepest green-up at the first date: green-up rate not observed")
  } else if (ip > 1) {
    notes <- c(notes, sprintf("rise before the peak is only %.0f%% of the larger side: green-up not reached in the observed window",
                              100 * amp_pre / amp_max))
  }
  # senescence: first downward crossing after the peak (before the later minimum)
  if (ip < n && amp_max > 0 && amp_post >= min_limb_frac * amp_max) {
    i1 <- ip - 1 + which.min(ff[ip:n])
    thr <- base_post + threshold * amp_post
    below <- ff[ip:i1] < thr
    k <- ip - 1 + which(below)[1]
    out$senescence_t <- crossing(tt[k - 1], tt[k], ff[k - 1], ff[k], thr)
    nc <- sum(diff(below) == 1)
    if (nc > 1) notes <- c(notes, sprintf("falling side crosses the senescence threshold %d times (curve not monotonic); the first crossing is reported", nc))
    seg <- ip:(i1 - 1)
    kk <- seg[which.min(d[seg])]
    if (kk < n - 1) out$senescence_rate <- -d[kk] else
      notes <- c(notes, "steepest decline at the last date: senescence rate not observed")
  } else if (ip < n) {
    notes <- c(notes, sprintf("decline after the peak is only %.0f%% of the larger side: senescence not reached in the observed window",
                              100 * amp_post / amp_max))
  }
  out$season_length <- out$senescence_t - out$greenup_t
  out$notes <- notes
  out
}

# Model-free summaries from the observations themselves.
simple_metrics <- function(t, y) {
  o <- order(t)
  t <- t[o]
  y <- y[o]
  n <- length(t)
  list(n_dates = n, t_first = t[1], t_last = t[n], max_obs_value = max(y),
       max_obs_t = t[which.max(y)],
       auc_obs = if (n > 1) sum(diff(t) * (y[-1] + y[-n]) / 2) else NA_real_)
}

# Robust residual scale (1.4826 x MAD).
robust_scale <- function(r) {
  r <- r[is.finite(r)]
  if (length(r) < 3) return(NA_real_)
  1.4826 * stats::median(abs(r - stats::median(r)))
}

# Noise SD pooled over plots: robust scale of all residuals, rescaled by
# n / residual df because a fitted curve absorbs part of the noise. Pooled
# because noise is a property of the sensor and pipeline, and ~15 residuals per
# plot give too unstable a scale. `parts`: list of list(r = residuals, df = df).
noise_sd <- function(parts) {
  if (!length(parts)) return(NA_real_)
  r <- unlist(lapply(parts, `[[`, "r"))
  df <- sum(vapply(parts, function(p) max(p$df, 0), 0))
  s <- robust_scale(r)
  if (is.finite(s) && df > 0) s <- s * sqrt(sum(is.finite(r)) / df)
  s
}

# Flag observations whose residual from the reference fit exceeds OUTLIER_Z x
# the pooled noise SD. The reference is a bisquare fit of the same curve with
# the pooled noise SD as its fixed scale, so an outlier is judged the same way
# whatever weighting the reported curve uses: an unweighted fit dragged down by
# a cloud would otherwise flag the cloud's neighbours too, and an upper-envelope
# fit can absorb an upward error into the curve.
flag_outliers <- function(r, sigma) {
  out <- rep("", length(r))
  if (!is.finite(sigma) || sigma <= 0) return(out)
  lim <- OUTLIER_Z * sigma
  out[is.finite(r) & r < -lim] <- "low"
  out[is.finite(r) & r > lim] <- "high"
  out
}

# Analyse one plot. Returns status, metrics, per-observation results and the
# fitted curve on the evaluation grid.
analyse_plot <- function(t, y, w0, opt) {
  o <- order(t)
  t <- t[o]
  y <- y[o]
  w0 <- w0[o]
  n <- length(t)
  res <- list(status = "fitted", simple = simple_metrics(t, y), metrics = NULL,
              fitted = rep(NA_real_, n), weights = w0, fit = NULL,
              robust_resid = rep(NA_real_, n), notes = character(), order = o)
  need <- min_dates_for(opt$method, opt$min_dates)
  if (n < need) {
    res$status <- "sparse"
    res$notes <- sprintf("only %d date(s), fewer than the %d needed for %s: simple metrics only",
                         n, need, opt$method)
    return(res)
  }
  if (opt$method %in% names(N_PAR)) {
    ys <- if (n >= 5) stats::runmed(y, 3, endrule = "keep") else y
    tm <- t[which.max(ys)]
    if (sum(t < tm) < 2 || sum(t > tm) < 2) {
      res$status <- "limb_missing"
      res$notes <- sprintf("%d date(s) before and %d after the highest value: a double logistic needs both the rise and the decline sampled (at least 2 dates each). Simple metrics only; --method loess can still give peak and green-up.",
                           sum(t < tm), sum(t > tm))
      return(res)
    }
  }
  fit <- fit_series(t, y, w0, opt$method, opt$weighting, opt)
  if (is.null(fit)) {
    res$status <- "fit_failed"
    res$notes <- "curve fit did not converge: simple metrics only"
    return(res)
  }
  grid <- seq(t[1], t[n], length.out = GRID_N)
  ff <- fit$predict(grid)
  res$metrics <- curve_metrics(grid, ff, opt$threshold, opt$min_limb_frac)
  res$fitted <- fit$fitted
  res$weights <- fit$weights
  res$fit <- fit
  res$t <- t
  res$y <- y
  res$w0 <- w0
  res$notes <- c(fit$note, res$metrics$notes,
                 if (length(fit$at_bound)) paste("parameter(s) at a bound:", paste(fit$at_bound, collapse = ", ")))
  res$grid <- grid
  res$curve <- ff
  res
}

# Convert axis values back to calendar dates for one plot.
axis_to_date <- function(tv, axis, planting, year, daymap_plot) {
  out <- rep(as.Date(NA), length(tv))
  ok <- is.finite(tv)
  if (!any(ok)) return(out)
  if (axis == "dap") out[ok] <- planting + round(tv[ok])
  if (axis == "doy") out[ok] <- as.Date(sprintf("%s-01-01", year)) + round(tv[ok]) - 1
  if (axis == "gdd") {
    num <- stats::approx(daymap_plot$gdd, as.numeric(daymap_plot$date), xout = tv[ok],
                         ties = min, rule = 2)$y
    out[ok] <- as.Date(round(num), origin = "1970-01-01")
  }
  out
}

run_timeseries <- function(prep, opt) {
  obs <- prep$obs
  ax <- axis_info(opt$axis)
  feats_at <- suppressWarnings(as.numeric(split_list(opt$features_at)))
  if (anyNA(feats_at)) stop("--features-at needs comma-separated numbers, e.g. 40,50,60.", call. = FALSE)
  keep_cols <- grep("^keep__", names(obs), value = TRUE)
  plots <- unique(obs$plot)
  rows <- vector("list", length(plots))
  obs_out <- vector("list", length(plots))
  curves <- vector("list", length(plots))
  params <- vector("list", length(plots))
  # Pass 1: fit every plot. Pass 2: estimate the noise SD pooled over plots,
  # refit a bisquare reference curve with that fixed scale, flag outliers
  # against it, and assemble the tables.
  fits <- lapply(plots, function(p) {
    d <- obs[obs$plot == p, , drop = FALSE]
    a <- analyse_plot(d$t, d$y, d$w0, opt)
    list(d = d[a$order, , drop = FALSE], a = a)
  })
  ok <- vapply(fits, function(f) f$a$status == "fitted", TRUE)
  # starting noise SD from unweighted fits (weighted fits have skewed residuals)
  sigma0 <- noise_sd(lapply(fits[ok], function(f) {
    a <- f$a
    b <- if (opt$weighting == "none") a$fit else fit_series(a$t, a$y, a$w0, opt$method, "none", opt)
    if (is.null(b)) b <- a$fit
    list(r = a$y - b$fitted, df = b$df_resid)
  }))
  ref_df <- list()
  for (i in which(ok)) {
    a <- fits[[i]]$a
    ref <- fit_series(a$t, a$y, a$w0, opt$method, "robust", opt, sigma = sigma0)
    if (is.null(ref)) ref <- a$fit
    fits[[i]]$a$robust_resid <- a$y - ref$fitted
    ref_df[[length(ref_df) + 1L]] <- list(r = fits[[i]]$a$robust_resid, df = ref$df_resid)
  }
  sigma <- noise_sd(ref_df)
  for (i in seq_along(plots)) {
    p <- plots[i]
    d <- fits[[i]]$d
    a <- fits[[i]]$a
    a$flag <- flag_outliers(a$robust_resid, sigma)
    dm <- if (!is.null(prep$daymap)) prep$daymap[prep$daymap$plot == p, , drop = FALSE] else NULL
    yr <- format(min(d$date), "%Y")
    m <- a$metrics
    s <- a$simple
    to_date <- function(v) format(axis_to_date(v, opt$axis, d$planting[1], yr, dm))
    row <- data.frame(plot = p, stringsAsFactors = FALSE)
    for (k in keep_cols) row[[sub("^keep__", "", k)]] <- d[[k]][1]
    row$status <- a$status
    row$method <- opt$method
    row$weighting <- opt$weighting
    row$n_dates <- s$n_dates
    row[[paste0("t_first_", ax$sfx)]] <- s$t_first
    row[[paste0("t_last_", ax$sfx)]] <- s$t_last
    row$max_obs_value <- s$max_obs_value
    row[[paste0("max_obs_", ax$sfx)]] <- s$max_obs_t
    row[[paste0("auc_obs_vi_", ax$len)]] <- s$auc_obs
    g <- function(name) if (is.null(m)) NA_real_ else m[[name]]
    row$peak_value <- g("peak_value")
    row[[paste0("peak_", ax$sfx)]] <- g("peak_t")
    row$peak_date <- to_date(g("peak_t"))
    row[[paste0("near_peak_start_", ax$sfx)]] <- g("near_peak_start_t")
    row[[paste0("near_peak_end_", ax$sfx)]] <- g("near_peak_end_t")
    row[[paste0("greenup_", ax$sfx)]] <- g("greenup_t")
    row$greenup_date <- to_date(g("greenup_t"))
    row[[paste0("senescence_", ax$sfx)]] <- g("senescence_t")
    row$senescence_date <- to_date(g("senescence_t"))
    row[[paste0("season_length_", ax$len)]] <- g("season_length")
    row[[paste0("greenup_rate_", ax$per)]] <- g("greenup_rate")
    row[[paste0("senescence_rate_", ax$per)]] <- g("senescence_rate")
    row[[paste0("auc_vi_", ax$len)]] <- g("auc")
    row$base_pre_value <- g("base_pre")
    row$base_post_value <- g("base_post")
    row$rmse <- if (a$status == "fitted") sqrt(mean((d$y - a$fitted)^2)) else NA_real_
    row$n_downweighted <- if (a$status == "fitted") sum(a$weights < 0.5 * d$w0) else NA_integer_
    row$n_low_outliers <- sum(a$flag == "low")
    row$n_high_outliers <- sum(a$flag == "high")
    for (fa in feats_at) {
      v <- NA_real_
      if (a$status == "fitted" && fa >= s$t_first && fa <= s$t_last) v <- a$fit$predict(fa)
      row[[sprintf("%s_at_%s%g", opt$vi, ax$sfx, fa)]] <- v
    }
    row$notes <- paste(a$notes, collapse = "; ")
    rows[[i]] <- row
    obs_out[[i]] <- data.frame(plot = p, date = format(d$date), t = d$t, observed = d$y,
                               fitted = a$fitted, residual = d$y - a$fitted,
                               weight = a$weights, robust_residual = a$robust_resid,
                               outlier = a$flag, stringsAsFactors = FALSE)
    if (a$status == "fitted") {
      days <- seq(min(d$date), max(d$date), by = "day")
      tv <- switch(opt$axis,
                   dap = as.numeric(days - d$planting[1]),
                   doy = as.numeric(format(days, "%j")),
                   gdd = dm$gdd[match(days, dm$date)])
      curves[[i]] <- data.frame(plot = p, date = format(days), t = tv,
                                fitted = a$fit$predict(tv), stringsAsFactors = FALSE)
      if (opt$method %in% names(N_PAR))
        params[[i]] <- data.frame(plot = p, model = opt$method, parameter = names(a$fit$par),
                                  estimate = unname(a$fit$par),
                                  stringsAsFactors = FALSE)
    }
  }
  metrics <- do.call(rbind, rows)
  observations <- do.call(rbind, obs_out)
  names(observations)[names(observations) == "t"] <- ax$sfx
  curves <- do.call(rbind, curves[!vapply(curves, is.null, TRUE)])
  if (!is.null(curves)) names(curves)[names(curves) == "t"] <- ax$sfx
  params <- do.call(rbind, params[!vapply(params, is.null, TRUE)])
  names(metrics)[1] <- opt$plot_col
  names(observations)[1] <- opt$plot_col
  if (!is.null(curves)) names(curves)[1] <- opt$plot_col
  if (!is.null(params)) names(params)[1] <- opt$plot_col
  list(metrics = metrics, observations = observations, curves = curves, params = params,
       noise_sd = sigma)
}

# Dates on which most plots sit above or below their curves point to a
# whole-flight problem (radiometric calibration, illumination), not to plots.
date_diagnostics <- function(observations, sigma) {
  o <- observations[is.finite(observations$robust_residual), , drop = FALSE]
  if (!nrow(o)) return(NULL)
  dd <- split(o$robust_residual, o$date)
  res <- data.frame(date = names(dd), n_plots = vapply(dd, length, 1L),
                    median_residual = vapply(dd, stats::median, 0),
                    stringsAsFactors = FALSE)
  res$flag <- res$n_plots >= 3 & is.finite(sigma) &
    abs(res$median_residual) > DATE_SHIFT_SD * sigma
  rownames(res) <- NULL
  res
}

# Candidate prediction features and the last axis value each depends on, for
# claw-yield-prediction --feature-time. QC columns and the observation-window
# bounds are left out: the window identifies the flight schedule, i.e. the
# environment, and would leak environment identity into a model.
feature_time_table <- function(metrics, opt) {
  ax <- axis_info(opt$axis)
  num <- names(metrics)[vapply(metrics, is.numeric, TRUE)]
  num <- setdiff(num, c("n_dates", "rmse", "n_downweighted", "n_low_outliers",
                        "n_high_outliers", paste0(c("t_first_", "t_last_"), ax$sfx),
                        split_list(opt$keep_cols)))
  avail <- max(metrics[[paste0("t_last_", ax$sfx)]], na.rm = TRUE)
  data.frame(feature = num, available_t = avail, axis = ax$sfx,
             note = "computed from observations up to available_t",
             stringsAsFactors = FALSE)
}

# --- figures ----------------------------------------------------------------

plot_examples <- function(res, opt, out_dir) {
  ax <- axis_info(opt$axis)
  met <- res$metrics
  pc <- opt$plot_col
  show <- utils::head(met[[pc]][met$status == "fitted"], 11)
  other <- met[[pc]][met$status != "fitted"]
  if (length(other)) show <- c(show, other[1])
  if (!length(show)) return(invisible(NULL))
  nr <- ceiling(length(show) / 4)
  open_png(out_dir, "curves_examples.png", width = 10, height = 2.6 * nr + 0.6)
  op <- graphics::par(mfrow = c(nr, min(4, length(show))), mar = c(3.5, 3.5, 2, 0.5),
                      mgp = c(2.2, 0.7, 0), oma = c(0, 0, 1.5, 0))
  ob <- res$observations
  yl <- range(ob$observed, na.rm = TRUE)
  for (p in show) {
    o <- ob[ob[[pc]] == p, , drop = FALSE]
    mrow <- met[met[[pc]] == p, , drop = FALSE]
    pch <- ifelse(o$weight < 0.5, 1, 19)
    graphics::plot(o[[ax$sfx]], o$observed, pch = pch, ylim = yl, xlab = ax$label,
                   ylab = opt$vi, main = sprintf("%s (%s)", p, mrow$status), cex.main = 0.9)
    if (!is.null(res$curves)) {
      cc <- res$curves[res$curves[[pc]] == p, , drop = FALSE]
      if (nrow(cc)) graphics::lines(cc[[ax$sfx]], cc$fitted, col = "#1b7837", lwd = 2)
    }
    v <- c(mrow[[paste0("greenup_", ax$sfx)]], mrow[[paste0("peak_", ax$sfx)]],
           mrow[[paste0("senescence_", ax$sfx)]])
    graphics::abline(v = v[is.finite(v)], lty = c(2, 1, 2)[is.finite(v)], col = "grey40")
    hi <- o$outlier != ""
    if (any(hi)) graphics::points(o[[ax$sfx]][hi], o$observed[hi], pch = 4, col = "#d6604d", cex = 1.5)
  }
  graphics::mtext(sprintf("%s: %s fit, weighting %s (open circles weight < 0.5; red x outliers; dashed green-up/senescence, solid peak)",
                          opt$vi, opt$method, opt$weighting), outer = TRUE, cex = 0.8)
  graphics::par(op)
  grDevices::dev.off()
}

plot_all_curves <- function(res, opt, out_dir) {
  if (is.null(res$curves)) return(invisible(NULL))
  ax <- axis_info(opt$axis)
  cc <- res$curves
  pc <- opt$plot_col
  open_png(out_dir, "curves_all.png")
  graphics::plot(range(cc[[ax$sfx]], na.rm = TRUE), range(cc$fitted, na.rm = TRUE), type = "n",
                 xlab = ax$label, ylab = paste("fitted", opt$vi),
                 main = sprintf("Fitted %s curves, %d plots (observed windows only)",
                                opt$vi, length(unique(cc[[pc]]))))
  for (p in unique(cc[[pc]])) {
    s <- cc[cc[[pc]] == p, , drop = FALSE]
    graphics::lines(s[[ax$sfx]], s$fitted, col = grDevices::adjustcolor("#1b7837", 0.4))
  }
  grDevices::dev.off()
}

plot_distributions <- function(res, opt, out_dir) {
  ax <- axis_info(opt$axis)
  met <- res$metrics[res$metrics$status == "fitted", , drop = FALSE]
  cols <- c(paste0(c("greenup_", "peak_", "senescence_"), ax$sfx), paste0("auc_vi_", ax$len))
  cols <- cols[vapply(cols, function(k) sum(is.finite(met[[k]])) >= 2, TRUE)]
  if (!length(cols)) return(invisible(NULL))
  open_png(out_dir, "metric_distributions.png", width = 3 * length(cols), height = 3.4)
  op <- graphics::par(mfrow = c(1, length(cols)), mar = c(4, 4, 2, 0.5))
  for (k in cols) graphics::hist(met[[k]], main = k, xlab = k, col = "#a6dba0", border = "white")
  graphics::par(op)
  grDevices::dev.off()
}

# --- report -------------------------------------------------------------------

fmt <- function(x, d = 3) ifelse(is.finite(x), as.character(signif(x, d)), "NA")

summary_rows <- function(met, cols) {
  vapply(cols, function(k) {
    v <- met[[k]][is.finite(met[[k]])]
    if (!length(v)) return(sprintf("| %s | 0 | NA | NA | NA |", k))
    q <- stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE)
    sprintf("| %s | %d | %s | %s | %s |", k, length(v), fmt(q[2], 4), fmt(q[1], 4), fmt(q[3], 4))
  }, "")
}

write_report <- function(res, prep, opt, out_dir, ddiag) {
  ax <- axis_info(opt$axis)
  met <- res$metrics
  pc <- opt$plot_col
  warn <- prep$notes
  st <- table(factor(met$status, levels = c("fitted", "sparse", "limb_missing", "fit_failed")))
  if (st[["sparse"]])
    warn <- c(warn, sprintf("%d plot(s) have too few dates for a %s curve and got simple metrics only: %s.",
                            st[["sparse"]], opt$method,
                            paste(utils::head(met[[pc]][met$status == "sparse"], 10), collapse = ", ")))
  if (st[["limb_missing"]])
    warn <- c(warn, sprintf("%d plot(s) do not sample both the rise and the decline, so the double logistic was not fitted: %s.",
                            st[["limb_missing"]],
                            paste(utils::head(met[[pc]][met$status == "limb_missing"], 10), collapse = ", ")))
  if (st[["fit_failed"]])
    warn <- c(warn, sprintf("%d plot(s): the curve fit failed: %s.", st[["fit_failed"]],
                            paste(utils::head(met[[pc]][met$status == "fit_failed"], 10), collapse = ", ")))
  nlo <- sum(met$n_low_outliers, na.rm = TRUE)
  nhi <- sum(met$n_high_outliers, na.rm = TRUE)
  if (nlo + nhi > 0)
    warn <- c(warn, sprintf("%d low and %d high outlier observation(s): more than %g noise SD (%s) from a bisquare reference fit; see tables/observations.csv.%s%s",
                            nlo, nhi, OUTLIER_Z, fmt(res$noise_sd),
                            if (opt$weighting == "none") " Check them; for cloud or shadow drops use --weighting upper, for isolated errors either way --weighting robust." else "",
                            if (opt$weighting == "upper" && nhi > 0) " Upper-envelope weighting keeps full weight on high values, so the high outliers pull their curves up: check them, or use --weighting robust." else ""))
  if (!is.null(ddiag) && any(ddiag$flag))
    warn <- c(warn, sprintf("Date(s) where most plots sit off their curves (median residual over plots beyond %g noise SD): %s. Suspect a whole-flight problem (radiometric calibration, illumination) rather than the plots; with loess or sg, smoother bias near the window edges can also cause this.",
                            DATE_SHIFT_SD, paste(ddiag$date[ddiag$flag], collapse = ", ")))
  fit <- met[met$status == "fitted", , drop = FALSE]
  tf <- fit[[paste0("t_first_", ax$sfx)]]
  tl <- fit[[paste0("t_last_", ax$sfx)]]
  if (nrow(fit) > 1) {
    tol <- if (opt$axis == "gdd") 70 else 7
    if (diff(range(tf)) > tol || diff(range(tl)) > tol)
      warn <- c(warn, sprintf("Observed windows differ between plots (first %s %s to %s, last %s to %s). AUC and threshold dates depend on the window, so compare plots with similar windows only.",
                              ax$sfx, fmt(min(tf)), fmt(max(tf)), fmt(min(tl)), fmt(max(tl))))
  }
  nb <- sum(grepl("parameter\\(s\\) at a bound", met$notes))
  if (nb) warn <- c(warn, sprintf("%d plot(s) have curve parameters at a bound (see the notes column): treat those curves with caution.", nb))
  nsen <- sum(is.na(fit[[paste0("senescence_", ax$sfx)]]))
  if (nrow(fit) && nsen)
    warn <- c(warn, sprintf("Senescence was not reached within the observed window for %d of %d fitted plot(s); no value is extrapolated.", nsen, nrow(fit)))
  if (opt$axis == "doy")
    warn <- c(warn, "Axis is calendar day of year: plots planted on different dates are not comparable. Prefer --axis dap or gdd when planting dates are known.")
  if (opt$axis == "gdd")
    warn <- c(warn, sprintf("GDD axis from '%s' column '%s': confirm its base temperature and units (deg C or deg F days) with the user before comparing with other studies.",
                            basename(opt$gdd), opt$gdd_col))

  cols <- c(paste0("peak_", ax$sfx), "peak_value", paste0("greenup_", ax$sfx),
            paste0("senescence_", ax$sfx), paste0("season_length_", ax$len),
            paste0("senescence_rate_", ax$per), paste0("auc_vi_", ax$len), "rmse")
  lines <- c(sprintf("# %s report", SKILL), "",
             if (length(warn)) c("## Warnings", "", paste("-", warn), "") else NULL,
             "## Summary", "",
             sprintf("- Input: `%s`", opt$input),
             sprintf("- VI: `%s`; axis: %s; method: %s; weighting: %s", opt$vi, ax$label, opt$method, opt$weighting),
             sprintf("- Plots: %d (fitted %d, sparse %d, limb missing %d, fit failed %d); observations: %d",
                     nrow(met), st[["fitted"]], st[["sparse"]], st[["limb_missing"]], st[["fit_failed"]],
                     nrow(res$observations)),
             sprintf("- Fewest dates for a curve: %d; threshold: %g of amplitude; min limb fraction: %g",
                     min_dates_for(opt$method, opt$min_dates), opt$threshold, opt$min_limb_frac),
             sprintf("- Noise SD (pooled robust estimate, corrected for fitted degrees of freedom): %s %s units; outliers are flagged beyond %g noise SD from a bisquare reference fit",
                     fmt(res$noise_sd), opt$vi, OUTLIER_Z),
             if (!is.na(opt$cutoff)) sprintf("- Observations used up to %s %g (--cutoff).", ax$sfx, opt$cutoff),
             "", "### Metrics across fitted plots", "",
             "| metric | n | median | lower quartile | upper quartile |", "|---|---|---|---|---|",
             summary_rows(fit, cols), "",
             "### Definitions", "",
             sprintf("- **peak**: maximum of the fitted curve inside the observed window; near-peak window: where the curve is within %g%% of the amplitude of the peak (a wide window means the peak date itself is poorly defined).", 100 * NEAR_PEAK_FRAC),
             sprintf("- **green-up / senescence**: where the curve crosses base + %g x amplitude, using the minimum before the peak (green-up) or after it (senescence) as base. Reported only if that side changes by at least %g of the larger side; otherwise NA.", opt$threshold, opt$min_limb_frac),
             "- **rates**: steepest rise / decline of the fitted curve (VI units per axis unit); NA when the steepest point is at the window edge.",
             "- **AUC**: integral of the fitted curve over the observed window (trapezoid rule); `auc_obs` uses the raw observations.",
             "- No curve is evaluated outside each plot's first-to-last observation.", "",
             "## Figures", "",
             "![Example curves](figures/curves_examples.png)", "",
             if (file.exists(file.path(out_dir, "figures", "curves_all.png"))) c("![All curves](figures/curves_all.png)", "") else NULL,
             if (file.exists(file.path(out_dir, "figures", "metric_distributions.png"))) c("![Metric distributions](figures/metric_distributions.png)", "") else NULL,
             "## Limitations", "",
             "- Metrics are defined relative to each plot's observed window: series that start after emergence or stop before maturity give window-dependent green-up, senescence and AUC.",
             "- The double logistic assumes one growing cycle; regrowth, lodging or replanting break it.",
             "- Peak dates of flat-topped curves are poorly defined; prefer the near-peak window or the threshold dates.",
             "- Upper-envelope weighting trusts high values: an upward outlier (e.g. a mis-registered plot) is flagged but not down-weighted. Use --weighting robust for errors in either direction.", "",
             "## References", "",
             "- Beck PSA et al. (2006) Remote Sensing of Environment 100:321-334. https://doi.org/10.1016/j.rse.2005.10.021",
             "- Elmore AJ et al. (2012) Global Change Biology 18:656-674. https://doi.org/10.1111/j.1365-2486.2011.02521.x",
             "- Cleveland WS (1979) Journal of the American Statistical Association 74:829-836. https://doi.org/10.1080/01621459.1979.10481038",
             "- Jonsson P, Eklundh L (2002) IEEE Transactions on Geoscience and Remote Sensing 40:1824-1832. https://doi.org/10.1109/TGRS.2002.802519",
             "- Chen J et al. (2004) Remote Sensing of Environment 91:332-344. https://doi.org/10.1016/j.rse.2004.03.014",
             "- Savitzky A, Golay MJE (1964) Analytical Chemistry 36:1627-1639. https://doi.org/10.1021/ac60214a047",
             "- Elzhov TV et al. minpack.lm R package. https://doi.org/10.32614/CRAN.package.minpack.lm", "",
             sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(lines, file.path(out_dir, "report.md"))
  warn
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (!opt$method %in% METHODS)
    stop("--method must be one of: ", paste(METHODS, collapse = ", "), call. = FALSE)
  if (!opt$weighting %in% WEIGHTINGS)
    stop("--weighting must be one of: ", paste(WEIGHTINGS, collapse = ", "), call. = FALSE)
  if (opt$threshold <= 0 || opt$threshold >= 1) stop("--threshold must be between 0 and 1.", call. = FALSE)
  if (opt$min_limb_frac < 0 || opt$min_limb_frac > 1) stop("--min-limb-frac must be between 0 and 1.", call. = FALSE)
  if (opt$demo) {
    opt$input <- file.path(SKILL_DIR, "demo", "vi_long.csv")
    if (!nzchar(opt$keep_cols)) opt$keep_cols <- "env,genotype"
    if (opt$axis == "gdd" && !nzchar(opt$gdd)) {
      opt$gdd <- file.path(SKILL_DIR, "demo", "gdd_daily.csv")
      if (!nzchar(opt$env_col)) opt$env_col <- "env"
    }
  }
  if (!nzchar(opt$input)) stop("Provide --input FILE.csv, or --demo. See --help.", call. = FALSE)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

  dat <- read_csv_checked(opt$input)
  gdd_tab <- if (nzchar(opt$gdd)) read_csv_checked(opt$gdd, "GDD") else NULL
  prep <- prepare_series(dat, opt, gdd_tab)
  res <- run_timeseries(prep, opt)
  ax <- axis_info(opt$axis)
  ddiag <- date_diagnostics(res$observations, res$noise_sd)

  save_table(res$metrics, opt$out, "metrics.csv")
  save_table(res$observations, opt$out, "observations.csv")
  if (!is.null(res$curves)) save_table(res$curves, opt$out, "fitted_curves.csv")
  if (!is.null(res$params)) save_table(res$params, opt$out, "dl_parameters.csv")
  if (!is.null(ddiag)) save_table(ddiag, opt$out, "date_diagnostics.csv")
  save_table(feature_time_table(res$metrics, opt), opt$out, "feature_time.csv")

  plot_examples(res, opt, opt$out)
  plot_all_curves(res, opt, opt$out)
  plot_distributions(res, opt, opt$out)
  warn <- write_report(res, prep, opt, opt$out, ddiag)
  write_repro_bundle(opt$out, argv, inputs = c(opt$input, opt$gdd))
  cat(sprintf("Done: %d plots (%d fitted), %d warning(s). Report: %s\n", nrow(res$metrics),
              sum(res$metrics$status == "fitted"), length(warn),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(res)
}

if (sys.nframe() == 0L) main()
