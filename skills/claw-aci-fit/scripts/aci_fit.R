#!/usr/bin/env Rscript
# claw-aci-fit: fit C3 (Farquhar-von Caemmerer-Berry) or C4 (von Caemmerer)
# CO2-response (A-Ci) curves with PhotoGEA, report Vcmax, J, Jmax, TPU and RL
# at 25 C and at leaf temperature with confidence intervals, and run curve QC
# (Ci coverage, points per limitation, induction on return to reference CO2,
# stomatal closure at high CO2, residuals).
#
# Input: the tidy CSV written by claw-licor-6800 (or any CSV with the same
# column names and units).
#
# Usage:
#   Rscript scripts/aci_fit.R --demo
#   Rscript scripts/aci_fit.R --input licor_tidy.csv --curve-cols genotype,plot --cores 4
#   Rscript scripts/aci_fit.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT license.

SKILL   <- "claw-aci-fit"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/aci_fit.R [--demo] [--input FILE.csv] [--curve-cols COLS]
                             [--pathway C3|C4] [--temp-params bernacchi|sharkey|flat]
                             [--alpha-old fit|NUMBER] [--ca-atmospheric PPM]
                             [--optimizer deoptim|nmkb] [--generations N] [--seed N]
                             [--cores N] [--exclude-unstable] [--include-flagged]
                             [--keep-recovery-points] [--oxygen-pct PCT] [--pressure-kpa KPA]
                             [--induction-tol F] [--low-ci PPM] [--high-ci PPM]
                             [--min-points N] [--closure-ratio F] [--rmse-frac F] [--out DIR]

  --demo                  fit the bundled synthetic curves (demo/aci_demo_c3.csv, or
                          demo/aci_demo_c4.csv with --pathway C4)
  --input                 tidy CSV from claw-licor-6800 (one row per logged point)
  --curve-cols            column(s) that together identify one curve, comma-separated
  --pathway               C3 (default) or C4 (maize, sorghum, other NADP-ME C4 crops)
  --temp-params           C3 temperature responses: bernacchi (default), sharkey or flat
  --alpha-old             C3 TPU parameter alpha_old: fit (PhotoGEA default) or a fixed number
  --ca-atmospheric        growth CO2 (umol/mol) for the operating point; default: not computed
  --optimizer             deoptim (default, PhotoGEA's recommendation) or nmkb (fast, can
                          stop in local minima)
  --generations           DEoptim generations per curve (default 200, PhotoGEA default)
  --seed                  random seed; curve i (sorted IDs) uses seed + i - 1 (default 1234)
  --cores                 curves fitted in parallel (default 1)
  --exclude-unstable      also exclude points logged before stability (qc_unstable)
  --include-flagged       keep points with qc_hard_fail = TRUE in the fit (not advised)
  --keep-recovery-points  fit the repeated reference-CO2 points too (default: excluded)
  --oxygen-pct            O2 (%) when the CSV has no oxygen_pct column (e.g. 21)
  --pressure-kpa          total chamber pressure (kPa) when the CSV has no pa_kpa column
  --induction-tol         flag |A_return / A_initial - 1| above this (default 0.10)
  --low-ci, --high-ci     Ci thresholds for coverage checks (defaults C3 200/600, C4 100/400)
  --min-points            minimum points per limitation and per Ci range (default 3)
  --closure-ratio         flag gsw(highest CO2) / gsw(initial) below this (default 0.5)
  --rmse-frac             flag RMSE above this fraction of the curve's max A (default 0.05)
  --out                   output folder (default: agriclaw_output/claw-aci-fit)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, input = "", curve_cols = "", pathway = "C3",
                 temp_params = "bernacchi", alpha_old = "fit", ca_atmospheric = NA_real_,
                 optimizer = "deoptim", generations = 200, seed = 1234, cores = 1,
                 exclude_unstable = FALSE, include_flagged = FALSE,
                 keep_recovery_points = FALSE, oxygen_pct = NA_real_, pressure_kpa = NA_real_,
                 induction_tol = 0.10, low_ci = NA_real_, high_ci = NA_real_, min_points = 3,
                 closure_ratio = 0.5, rmse_frac = 0.05, out = "")

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

U_A <- "micromol m^(-2) s^(-1)"
U_C <- "micromol mol^(-1)"

split_cols <- function(s) {
  if (!nzchar(s)) return(character())
  trimws(strsplit(s, ",", fixed = TRUE)[[1L]])
}

read_tidy <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "", "txt"))
    stop("'", basename(path), "' looks like a raw LI-6800 log. Run claw-licor-6800 first and ",
         "pass its tables/licor_tidy.csv here.", call. = FALSE)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

# Check columns and units; fail loudly with the fix.
validate_input <- function(df, curve_cols, opt) {
  need <- c("a_umol_m2_s", "ci_umol_mol", "ca_umol_mol", "tleaf_c")
  miss <- setdiff(need, names(df))
  if (length(miss))
    stop("Missing column(s): ", paste(miss, collapse = ", "), ". The input must be the tidy CSV ",
         "from claw-licor-6800 (A in umol m-2 s-1, Ci and Ca in umol/mol, leaf temperature in C).",
         call. = FALSE)
  if (!length(curve_cols))
    stop("--curve-cols is required: which column(s) identify one A-Ci curve (e.g. genotype,plot)? ",
         "Candidate columns: ", paste(setdiff(names(df), c(need, grep("_umol|_mol_|_kpa|_pct|_c$|^qc_|_s$", names(df), value = TRUE))), collapse = ", "),
         call. = FALSE)
  miss <- setdiff(curve_cols, names(df))
  if (length(miss)) stop("--curve-cols column(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
  for (cn in need) if (!is.numeric(df[[cn]])) stop("Column ", cn, " must be numeric.", call. = FALSE)
  if (!"oxygen_pct" %in% names(df) && is.na(opt$oxygen_pct))
    stop("No oxygen_pct column. If the curves were measured in normal air, pass --oxygen-pct 21.",
         call. = FALSE)
  if ((!"pa_kpa" %in% names(df) || all(is.na(df$pa_kpa))) && is.na(opt$pressure_kpa))
    stop("No pa_kpa values (atmospheric pressure). Pass --pressure-kpa with the chamber pressure in kPa.",
         call. = FALSE)
  a <- df$a_umol_m2_s
  ci <- df$ci_umol_mol
  if (stats::median(ci, na.rm = TRUE) < 60)
    stop("Median Ci is ", signif(stats::median(ci, na.rm = TRUE), 3), ": Ci looks like Pa, not umol/mol. ",
         "Convert (umol/mol = Pa / total pressure in kPa * 1000).", call. = FALSE)
  if (max(a, na.rm = TRUE) > 100)
    stop("A exceeds 100 umol m-2 s-1: check the units of a_umol_m2_s.", call. = FALSE)
  tl <- df$tleaf_c
  if (any(tl > 60 | tl < -5, na.rm = TRUE))
    stop("Leaf temperatures outside -5..60 C: tleaf_c must be in degrees C.", call. = FALSE)
  invisible(TRUE)
}

curve_key <- function(df, curve_cols) {
  do.call(paste, c(lapply(curve_cols, function(cn) as.character(df[[cn]])), sep = " - "))
}

# Measurement order within a curve: time, then elapsed time, then obs, then
# input order.
order_rows <- function(df) {
  idx <- seq_len(nrow(df))
  keys <- list(df$curve_id)
  for (cn in c("time_unix_s", "elapsed_s", "obs"))
    if (cn %in% names(df) && !all(is.na(df[[cn]]))) { keys <- c(keys, list(df[[cn]])); break }
  df[do.call(order, c(keys, list(idx))), , drop = FALSE]
}

# The CO2 variable that was stepped: setpoint if logged, else reference CO2,
# else Ca. Returns the values and a relative/absolute tolerance for "same level".
driver_of <- function(d) {
  if ("co2_r_sp_umol_mol" %in% names(d) && !all(is.na(d$co2_r_sp_umol_mol)))
    return(list(x = d$co2_r_sp_umol_mol, abs_tol = 2, rel_tol = 0, name = "co2_r_sp_umol_mol"))
  if ("co2_r_umol_mol" %in% names(d) && !all(is.na(d$co2_r_umol_mol)))
    return(list(x = d$co2_r_umol_mol, abs_tol = 0, rel_tol = 0.03, name = "co2_r_umol_mol"))
  list(x = d$ca_umol_mol, abs_tol = 0, rel_tol = 0.05, name = "ca_umol_mol")
}

# Mark the points measured after the curve returned to its starting (reference)
# CO2 level, and compute the induction ratio A_return / A_initial from the last
# point of the initial block and the last point of the return block.
mark_recovery <- function(d) {
  n <- nrow(d)
  out <- list(recovery = rep(FALSE, n), ratio = NA_real_, gsw_ratio = NA_real_,
              a_initial = NA_real_, a_return = NA_real_, driver = NA_character_)
  if (n < 3L) return(out)
  dr <- driver_of(d)
  out$driver <- dr$name
  x <- dr$x
  ref <- x[1L]
  if (is.na(ref)) return(out)
  same <- !is.na(x) & abs(x - ref) <= pmax(dr$abs_tol, dr$rel_tol * abs(ref))
  depart <- which(!same)
  if (!length(depart)) return(out)
  first_dep <- depart[1L]
  later <- seq_len(n) > first_dep
  rec <- same & later
  out$recovery <- rec
  if (!any(rec)) return(out)
  # last point of the first return block (the one logged just before CO2 moved on)
  rec_idx <- which(rec)
  block_end <- rec_idx[1L]
  while ((block_end + 1L) %in% rec_idx) block_end <- block_end + 1L
  init_end <- first_dep - 1L
  out$a_initial <- d$a_umol_m2_s[init_end]
  out$a_return <- d$a_umol_m2_s[block_end]
  if (is.finite(out$a_initial) && out$a_initial > 0) out$ratio <- out$a_return / out$a_initial
  if ("gsw_mol_m2_s" %in% names(d)) {
    g0 <- d$gsw_mol_m2_s[init_end]
    if (is.finite(g0) && g0 > 0) out$gsw_ratio <- d$gsw_mol_m2_s[block_end] / g0
  }
  out
}

prepare_points <- function(df, opt, curve_cols) {
  df$curve_id <- curve_key(df, curve_cols)
  df <- order_rows(df)
  df$recovery_point <- FALSE
  df$excluded_reason <- ""
  ind <- list()
  for (id in unique(df$curve_id)) {
    rows <- which(df$curve_id == id)
    m <- mark_recovery(df[rows, , drop = FALSE])
    df$recovery_point[rows] <- m$recovery
    ind[[id]] <- data.frame(curve_id = id, a_initial_umol_m2_s = m$a_initial,
                            a_return_umol_m2_s = m$a_return, induction_ratio = m$ratio,
                            gsw_return_ratio = m$gsw_ratio, n_recovery_points = sum(m$recovery),
                            recovery_driver = m$driver, stringsAsFactors = FALSE)
  }
  reason <- character(nrow(df))
  add <- function(cond, why) {
    cond <- cond %in% TRUE
    reason[cond] <<- ifelse(nzchar(reason[cond]), paste(reason[cond], why, sep = ";"), why)
  }
  if (!opt$include_flagged && "qc_hard_fail" %in% names(df)) add(df$qc_hard_fail, "qc_hard_fail")
  if (opt$exclude_unstable && "qc_unstable" %in% names(df)) add(df$qc_unstable, "unstable")
  if (!opt$keep_recovery_points) add(df$recovery_point, "recovery_point")
  add(is.na(df$a_umol_m2_s) | is.na(df$ci_umol_mol) | is.na(df$tleaf_c), "missing_A_Ci_or_Tleaf")
  df$excluded_reason <- reason
  df$used_in_fit <- !nzchar(reason)
  list(points = df, induction = do.call(rbind, ind))
}

# Slim PhotoGEA exdf (few columns keeps the optimizer fast) with units. Rows
# are sorted by Ci within each curve, as PhotoGEA expects (its interpolated
# model curve runs from the first to the last row).
build_exdf <- function(df, opt) {
  df <- df[order(df$curve_id, df$ci_umol_mol), , drop = FALSE]
  n <- nrow(df)
  if ("pa_kpa" %in% names(df) && !all(is.na(df$pa_kpa))) {
    pa <- df$pa_kpa
    dp <- if ("delta_pcham_kpa" %in% names(df)) df$delta_pcham_kpa else rep(0, n)
  } else {
    pa <- rep(opt$pressure_kpa, n)
    dp <- rep(0, n)
  }
  if (any(is.na(pa) | is.na(dp)))
    stop("Missing pressure values (pa_kpa / delta_pcham_kpa) in some rows.", call. = FALSE)
  o2 <- if ("oxygen_pct" %in% names(df)) df$oxygen_pct else rep(opt$oxygen_pct, n)
  if (any(is.na(o2))) o2[is.na(o2)] <- opt$oxygen_pct
  if (any(is.na(o2)))
    stop("oxygen_pct is missing for some rows; pass --oxygen-pct to fill them explicitly.", call. = FALSE)
  md <- data.frame(curve_id = df$curve_id, A = df$a_umol_m2_s, Ci = df$ci_umol_mol,
                   Ca = df$ca_umol_mol, TleafCnd = df$tleaf_c, Pa = pa, DeltaPcham = dp,
                   Oxygen = o2, include_when_fitting = df$used_in_fit, stringsAsFactors = FALSE)
  units <- data.frame(A = U_A, Ci = U_C, Ca = U_C, TleafCnd = "degrees C", Pa = "kPa",
                      DeltaPcham = "kPa", Oxygen = "percent", stringsAsFactors = FALSE)
  # Curve means over the points used in the fit (needed for Jmax and reporting).
  use <- df$used_in_fit
  md$TleafCnd_avg <- stats::ave(ifelse(use, df$tleaf_c, NA), df$curve_id,
                                FUN = function(v) mean(v, na.rm = TRUE))
  units$TleafCnd_avg <- "degrees C"
  if ("qin_umol_m2_s" %in% names(df) && !all(is.na(df$qin_umol_m2_s))) {
    md$Qin_avg <- stats::ave(ifelse(use, df$qin_umol_m2_s, NA), df$curve_id,
                             FUN = function(v) mean(v, na.rm = TRUE))
    units$Qin_avg <- U_A
  }
  x <- PhotoGEA::exdf(md, units = units)
  x <- PhotoGEA::calculate_total_pressure(x)
  tp <- switch(opt$pathway,
               C3 = switch(opt$temp_params, bernacchi = PhotoGEA::c3_temperature_param_bernacchi,
                           sharkey = PhotoGEA::c3_temperature_param_sharkey,
                           flat = PhotoGEA::c3_temperature_param_flat),
               C4 = PhotoGEA::c4_temperature_param_vc)
  PhotoGEA::calculate_temperature_response(x, tp, "TleafCnd")
}

# Fit one curve; runs in the main process or a parallel worker, so it only
# uses PhotoGEA:: and base functions.
fit_one_curve <- function(job, pathway, optimizer, generations, ca_atm, alpha_old) {
  set.seed(job$seed)
  opt_fun <- if (optimizer == "nmkb") PhotoGEA::optimizer_nmkb(1e-7) else
    PhotoGEA::optimizer_deoptim(generations)
  ca <- if (is.na(ca_atm)) NA else ca_atm
  tryCatch({
    if (pathway == "C3") {
      fo <- if (identical(alpha_old, "fit")) list() else list(alpha_old = as.numeric(alpha_old))
      PhotoGEA::fit_c3_aci(job$data, Ca_atmospheric = ca, optim_fun = opt_fun, fit_options = fo)
    } else {
      PhotoGEA::fit_c4_aci(job$data, Ca_atmospheric = ca, optim_fun = opt_fun)
    }
  }, error = function(e) list(error = conditionMessage(e)))
}

run_fits <- function(x, opt) {
  ids <- sort(unique(x[, "curve_id"]))
  jobs <- lapply(seq_along(ids), function(i)
    list(id = ids[i], seed = opt$seed + i - 1L, data = x[x[, "curve_id"] == ids[i], , TRUE]))
  args <- list(pathway = opt$pathway, optimizer = opt$optimizer, generations = opt$generations,
               ca_atm = opt$ca_atmospheric, alpha_old = opt$alpha_old)
  cores <- max(1L, min(as.integer(opt$cores), length(jobs)))
  if (cores > 1L) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    res <- do.call(parallel::parLapply, c(list(cl = cl, X = jobs, fun = fit_one_curve), args))
  } else {
    res <- lapply(jobs, function(j) do.call(fit_one_curve, c(list(job = j), args)))
  }
  names(res) <- ids
  list(results = res, seeds = stats::setNames(vapply(jobs, `[[`, 0, "seed"), ids))
}

PARAMS <- list(C3 = c(Vcmax = "vcmax", J = "j", Tp = "tp", RL = "rl"),
               C4 = c(Vcmax = "vcmax", Vpmax = "vpmax", RL = "rl"))

# Rename PhotoGEA parameter columns to unit-labeled output names.
param_table <- function(pdf, pathway) {
  out <- data.frame(curve_id = pdf$curve_id, stringsAsFactors = FALSE)
  for (p in names(PARAMS[[pathway]])) {
    nm <- PARAMS[[pathway]][[p]]
    g <- function(cn) if (cn %in% names(pdf)) as.numeric(pdf[[cn]]) else rep(NA_real_, nrow(pdf))
    out[[paste0(nm, "_25_umol_m2_s")]] <- g(paste0(p, "_at_25"))
    out[[paste0(nm, "_25_lower_umol_m2_s")]] <- g(paste0(p, "_at_25_lower"))
    out[[paste0(nm, "_25_upper_umol_m2_s")]] <- g(paste0(p, "_at_25_upper"))
    out[[paste0(nm, "_tleaf_umol_m2_s")]] <- g(paste0(p, "_tl_avg"))
    tr <- paste0(p, "_trust")
    if (tr %in% names(pdf)) out[[paste0(nm, "_trust")]] <- as.character(pdf[[tr]])
  }
  if (pathway == "C3") {
    jm <- c(Jmax_at_25 = "jmax_25_umol_m2_s", Jmax_at_25_lower = "jmax_25_lower_umol_m2_s",
            Jmax_at_25_upper = "jmax_25_upper_umol_m2_s", Jmax_tl = "jmax_tleaf_umol_m2_s")
    for (cn in names(jm))
      out[[jm[[cn]]]] <- if (cn %in% names(pdf)) as.numeric(pdf[[cn]]) else NA_real_
    out$jmax_msg <- if ("Jmax_at_25_msg" %in% names(pdf)) as.character(pdf$Jmax_at_25_msg) else NA_character_
    out$alpha_old <- if ("alpha_old" %in% names(pdf)) as.numeric(pdf$alpha_old) else NA_real_
    out$alpha_old_trust <- if ("alpha_old_trust" %in% names(pdf)) as.character(pdf$alpha_old_trust) else NA_character_
    for (cn in c("n_Ac_limiting", "n_Aj_limiting", "n_Ap_limiting"))
      out[[tolower(cn)]] <- if (cn %in% names(pdf)) as.integer(pdf[[cn]]) else NA_integer_
    out$gamma_star_25_umol_mol <- if ("Gamma_star_at_25" %in% names(pdf)) pdf$Gamma_star_at_25 else NA
    out$kc_25_umol_mol <- if ("Kc_at_25" %in% names(pdf)) pdf$Kc_at_25 else NA
    out$ko_25_mmol_mol <- if ("Ko_at_25" %in% names(pdf)) pdf$Ko_at_25 else NA
  } else {
    out$gmc_25_mol_m2_s_bar <- if ("gmc_at_25" %in% names(pdf)) pdf$gmc_at_25 else NA
    out$vpr_umol_m2_s <- if ("Vpr" %in% names(pdf)) pdf$Vpr else NA
  }
  out$tleaf_mean_c <- if ("TleafCnd_avg" %in% names(pdf)) pdf$TleafCnd_avg else NA
  out$qin_mean_umol_m2_s <- if ("Qin_avg" %in% names(pdf)) pdf$Qin_avg else NA
  out$rmse_umol_m2_s <- if ("RMSE" %in% names(pdf)) pdf$RMSE else NA
  out$aic <- if ("AIC" %in% names(pdf)) pdf$AIC else NA
  out$dof <- if ("dof" %in% names(pdf)) pdf$dof else NA
  out$operating_ci_umol_mol <- if ("operating_Ci" %in% names(pdf)) pdf$operating_Ci else NA
  out$operating_an_model_umol_m2_s <- if ("operating_An_model" %in% names(pdf)) pdf$operating_An_model else NA
  out$convergence <- if ("convergence" %in% names(pdf)) pdf$convergence else NA
  out
}

collect_results <- function(fr, opt) {
  ok <- vapply(fr$results, function(r) is.null(r$error), logical(1))
  good <- fr$results[ok]
  params <- NULL
  fits <- NULL
  interp <- NULL
  if (length(good)) {
    pe <- PhotoGEA::consolidate(lapply(good, function(r) list(parameters = r$parameters)))$parameters
    # Jmax needs at least one usable J estimate (calculate_jmax stops if all are NA).
    if (opt$pathway == "C3" && "Qin_avg" %in% colnames(pe) && any(is.finite(pe[, "J_at_25"]))) {
      pe <- PhotoGEA::calculate_temperature_response(pe, PhotoGEA::jmax_temperature_param_bernacchi,
                                                     "TleafCnd_avg")
      pe <- PhotoGEA::calculate_jmax(pe)
    }
    params <- param_table(pe$main_data, opt$pathway)
    fits <- do.call(rbind, lapply(names(good), function(id) {
      f <- good[[id]]$fits$main_data
      data.frame(curve_id = id, ci_umol_mol = f$Ci, a_umol_m2_s = f$A, a_fit_umol_m2_s = f$A_fit,
                 residual_umol_m2_s = f$A_residuals, used_in_fit = f$include_when_fitting,
                 limiting = if (opt$pathway == "C3") limiting_process(f) else NA_character_,
                 stringsAsFactors = FALSE)
    }))
    interp <- lapply(good, function(r) r$fits_interpolated$main_data)
  }
  failed <- names(fr$results)[!ok]
  if (length(failed)) {
    fp <- data.frame(curve_id = failed, stringsAsFactors = FALSE)
    params <- if (is.null(params)) fp else rbind_fill(list(params, fp))
  }
  params$fit_error <- ""
  for (id in failed) params$fit_error[params$curve_id == id] <- fr$results[[id]]$error
  params$seed <- unname(fr$seeds[params$curve_id])
  params <- params[order(params$curve_id), , drop = FALSE]
  list(params = params, fits = fits, interp = interp)
}

# Limiting process = smallest carboxylation rate (Wc, Wj, Wp), as in PhotoGEA;
# comparing Ac/Aj/Ap instead would give the wrong answer below Gamma*.
limiting_process <- function(f) {
  if (!all(c("Wc", "Wj", "Wp") %in% names(f))) return(rep(NA_character_, nrow(f)))
  m <- cbind(Ac = f$Wc, Aj = f$Wj, Ap = f$Wp)
  m[!is.finite(m)] <- Inf
  colnames(m)[apply(m, 1L, which.min)]
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

# Curve-level QC on the points used in the fit. Thresholds are rules of thumb.
curve_qc <- function(points, params, induction, opt) {
  low <- if (is.na(opt$low_ci)) if (opt$pathway == "C3") 200 else 100 else opt$low_ci
  high <- if (is.na(opt$high_ci)) if (opt$pathway == "C3") 600 else 400 else opt$high_ci
  rows <- lapply(params$curve_id, function(id) {
    d <- points[points$curve_id == id, , drop = FALSE]
    u <- d[d$used_in_fit, , drop = FALSE]
    ind <- induction[induction$curve_id == id, , drop = FALSE]
    p <- params[params$curve_id == id, , drop = FALSE]
    flags <- character()
    n_low <- sum(u$ci_umol_mol < low, na.rm = TRUE)
    n_high <- sum(u$ci_umol_mol > high, na.rm = TRUE)
    if (n_low < opt$min_points) flags <- c(flags, "few_low_ci_points")
    if (n_high < opt$min_points) flags <- c(flags, "few_high_ci_points")
    r <- ind$induction_ratio
    if (length(r) && !is.na(r)) {
      if (r < 1 - opt$induction_tol) flags <- c(flags, "induction_drop")
      if (r > 1 + opt$induction_tol) flags <- c(flags, "induction_rise")
    }
    closure <- NA_real_
    if ("gsw_mol_m2_s" %in% names(d)) {
      dr <- driver_of(d)$x
      init <- d$gsw_mol_m2_s[1L]
      top <- d$gsw_mol_m2_s[which.max(dr)]
      if (is.finite(init) && init > 0 && length(top) && is.finite(top)) closure <- top / init
      if (!is.na(closure) && closure < opt$closure_ratio) flags <- c(flags, "stomatal_closure_high_co2")
    }
    max_a <- suppressWarnings(max(u$a_umol_m2_s, na.rm = TRUE))
    if (!is.na(p$rmse_umol_m2_s %||% NA) && is.finite(max_a) && max_a > 0 &&
        p$rmse_umol_m2_s > opt$rmse_frac * max_a) flags <- c(flags, "poor_fit_rmse")
    if (!is.na(p$rl_25_umol_m2_s %||% NA) && p$rl_25_umol_m2_s < 0) flags <- c(flags, "negative_rl")
    if (opt$pathway == "C3") {
      for (pr in list(c("n_ac_limiting", "vcmax"), c("n_aj_limiting", "j"), c("n_ap_limiting", "tp"))) {
        k <- p[[pr[1]]]
        if (!is.null(k) && !is.na(k) && k < opt$min_points) flags <- c(flags, paste0(pr[2], "_few_limiting_points"))
      }
    }
    if (nzchar(p$fit_error %||% "")) flags <- c(flags, "fit_failed")
    data.frame(curve_id = id, n_points = nrow(d), n_points_used = nrow(u),
               n_ci_below_low = n_low, n_ci_above_high = n_high,
               max_ci_umol_mol = if (nrow(u)) max(u$ci_umol_mol, na.rm = TRUE) else NA_real_,
               gsw_high_co2_ratio = closure, qc_flags = paste(flags, collapse = ";"),
               stringsAsFactors = FALSE)
  })
  qc <- do.call(rbind, rows)
  merge(qc, induction, by = "curve_id", all.x = TRUE, sort = FALSE)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# Columns that are constant within every curve (genotype, plot, instrument...).
curve_metadata <- function(points, curve_cols) {
  cand <- setdiff(names(points), c("curve_id", "used_in_fit", "excluded_reason", "recovery_point"))
  cand <- cand[!startsWith(cand, "qc_")]
  cand <- cand[!vapply(points[cand], function(v) is.numeric(v) || is.logical(v), logical(1)) |
                 cand %in% curve_cols]
  keep <- cand[vapply(cand, function(cn) all(tapply(points[[cn]], points$curve_id,
                                                    function(v) length(unique(v)) == 1L)), logical(1))]
  keep <- setdiff(keep, c("qc_flags", "date_local", "hhmmss", "user_remark", "file_type"))
  m <- points[!duplicated(points$curve_id), c("curve_id", keep), drop = FALSE]
  rownames(m) <- NULL
  m
}

plot_fits <- function(points, interp, params, opt, out_dir) {
  ids <- params$curve_id
  per <- 9L
  pages <- split(ids, ceiling(seq_along(ids) / per))
  files <- character()
  for (pg in seq_along(pages)) {
    nm <- sprintf("aci_fits_%02d.png", pg)
    open_png(out_dir, nm, width = 10, height = 9)
    graphics::par(mfrow = c(3, 3), mar = c(4, 4, 2.6, 0.8), cex = 0.75)
    for (id in pages[[pg]]) {
      d <- points[points$curve_id == id, , drop = FALSE]
      it <- interp[[id]]
      yl <- range(c(d$a_umol_m2_s, 0), na.rm = TRUE)
      yl[2] <- yl[2] * 1.08
      graphics::plot(d$ci_umol_mol, d$a_umol_m2_s, type = "n", ylim = yl,
                     xlab = "Ci (umol/mol)", ylab = "A (umol m-2 s-1)",
                     main = substr(id, 1L, 40L))
      graphics::abline(h = 0, col = "grey80")
      if (!is.null(it)) {
        lims <- if (opt$pathway == "C3") c("Ac", "Aj", "Ap") else c("Apc", "Ar", "Aj")
        cols <- c("#1f77b4", "#d95f02", "#1b9e77")
        for (k in seq_along(lims)) if (lims[k] %in% names(it) && any(is.finite(it[[lims[k]]])))
          graphics::lines(it$Ci, it[[lims[k]]], col = cols[k], lty = 2)
        graphics::lines(it$Ci, it$An, lwd = 2)
      }
      u <- d$used_in_fit
      graphics::points(d$ci_umol_mol[u], d$a_umol_m2_s[u], pch = 16)
      graphics::points(d$ci_umol_mol[!u], d$a_umol_m2_s[!u], pch = ifelse(d$recovery_point[!u], 2, 1),
                       col = "grey40")
      fl <- params$qc_flags[params$curve_id == id]
      if (length(fl) && nzchar(fl)) graphics::mtext(substr(fl, 1L, 60L), side = 3, line = 0.1,
                                                    cex = 0.55, col = "firebrick")
    }
    grDevices::dev.off()
    files <- c(files, nm)
  }
  files
}

plot_residuals <- function(fits, out_dir) {
  if (is.null(fits) || !nrow(fits)) return(invisible(NULL))
  open_png(out_dir, "residuals.png", width = 7, height = 4.5)
  on.exit(grDevices::dev.off())
  u <- fits[fits$used_in_fit %in% TRUE, , drop = FALSE]
  f <- factor(u$curve_id)
  graphics::plot(u$ci_umol_mol, u$residual_umol_m2_s, col = grDevices::hcl.colors(max(2L, nlevels(f)), "Dark 3")[as.integer(f)],
                 pch = 16, xlab = "Ci (umol/mol)", ylab = "Measured - fitted A (umol m-2 s-1)",
                 main = "Residuals of the points used in the fits")
  graphics::abline(h = 0, lty = 2)
}

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
  few_low_ci_points = "fewer than %d points below Ci %g: Vcmax (C3) or Vpmax (C4) is poorly constrained",
  few_high_ci_points = "fewer than %d points above Ci %g: J/TPU (C3) or the CO2-saturated plateau (C4) is poorly constrained",
  induction_drop = "A on return to the reference CO2 fell by more than the tolerance: the leaf lost activation or was not at steady state; this can masquerade as down-regulation",
  induction_rise = "A on return to the reference CO2 rose by more than the tolerance: the leaf was not fully induced when the curve started, so the first (Rubisco-limited) points underestimate capacity",
  stomatal_closure_high_co2 = "gsw at the highest CO2 fell below the closure ratio of its starting value: stomatal closure limited the Ci reached",
  poor_fit_rmse = "RMSE above the allowed fraction of the curve's maximum A",
  negative_rl = "negative RL estimate: RL is poorly constrained by this curve",
  vcmax_few_limiting_points = "Rubisco limits fewer than %d points",
  j_few_limiting_points = "RuBP regeneration limits fewer than %d points: J (and Jmax) unreliable or NA",
  tp_few_limiting_points = "TPU limits fewer than %d points: TPU unreliable or NA (common; needs high Ci)",
  fit_failed = "the fit failed (see fit_error)")

build_warnings <- function(params, opt, low, high) {
  w <- character()
  if (opt$pathway == "C3")
    w <- c(w, "All C3 parameters are apparent (Ci-based): the fits assume infinite mesophyll conductance (Cc = Ci), so Vcmax and J are lower than chloroplastic values and are comparable only with other Ci-based estimates made with the same temperature-response parameters.")
  else
    w <- c(w, "C4 fits use PhotoGEA's default mesophyll conductance (1 mol m-2 s-1 bar-1 at 25 C) and bundle-sheath conductance (0.003), with von Caemmerer (2021) temperature responses estimated for Setaria viridis; Vcmax, Vpmax and RL depend on these assumptions.")
  if (opt$optimizer == "nmkb")
    w <- c(w, "Optimizer nmkb was used: it is fast but can stop in local minima. Refit with the default (deoptim) before publishing.")
  for (fl in names(FLAG_TEXT)) {
    hit <- params$curve_id[grepl(paste0("(^|;)", fl, "($|;)"), params$qc_flags)]
    if (length(hit)) {
      txt <- FLAG_TEXT[[fl]]
      if (grepl("%d points below", txt)) txt <- sprintf(txt, opt$min_points, low)
      else if (grepl("%d points above", txt)) txt <- sprintf(txt, opt$min_points, high)
      else if (grepl("%d", txt)) txt <- sprintf(txt, opt$min_points)
      w <- c(w, sprintf("%s (%d curve(s): %s).", paste0(toupper(substr(txt, 1, 1)), substring(txt, 2)),
                        length(hit), paste(hit, collapse = "; ")))
    }
  }
  w
}

write_report <- function(res, opt, input, out_dir, figs) {
  p <- res$params
  low <- if (is.na(opt$low_ci)) if (opt$pathway == "C3") 200 else 100 else opt$low_ci
  high <- if (is.na(opt$high_ci)) if (opt$pathway == "C3") 600 else 400 else opt$high_ci
  w <- build_warnings(p, opt, low, high)
  pts <- res$points
  if (opt$pathway == "C3") {
    show <- p[, c("curve_id", "vcmax_25_umol_m2_s", "vcmax_25_lower_umol_m2_s", "vcmax_25_upper_umol_m2_s",
                  "j_25_umol_m2_s", "jmax_25_umol_m2_s", "tp_25_umol_m2_s", "rl_25_umol_m2_s",
                  "tleaf_mean_c", "induction_ratio", "qc_flags")]
    names(show) <- c("curve", "Vcmax25", "lower", "upper", "J25", "Jmax25", "TPU25", "RL25", "Tleaf", "A_return/A_init", "flags")
  } else {
    show <- p[, c("curve_id", "vcmax_25_umol_m2_s", "vpmax_25_umol_m2_s", "vpmax_25_lower_umol_m2_s",
                  "vpmax_25_upper_umol_m2_s", "rl_25_umol_m2_s", "tleaf_mean_c", "induction_ratio", "qc_flags")]
    names(show) <- c("curve", "Vcmax25", "Vpmax25", "lower", "upper", "RL25", "Tleaf", "A_return/A_init", "flags")
  }
  tp_txt <- switch(opt$temp_params,
    bernacchi = "PhotoGEA `c3_temperature_param_bernacchi`: Bernacchi et al. (2001) for Vcmax, RL, Kc, Ko and Gamma*, Bernacchi et al. (2003) for J, Bernacchi et al. (2002) for gm, Sharkey et al. (2007) for TPU",
    sharkey = "PhotoGEA `c3_temperature_param_sharkey` (Sharkey et al. 2007; chloroplastic Kc, Ko and Gamma* from Bernacchi et al. 2002)",
    flat = "PhotoGEA `c3_temperature_param_flat` (no temperature dependence; 25 C constants from Sharkey et al. 2007)")
  if (opt$pathway == "C4") tp_txt <- "PhotoGEA `c4_temperature_param_vc` (von Caemmerer 2021)"
  lines <- c(sprintf("# %s report", SKILL), "",
             "## Warnings", "", if (length(w)) paste("-", w) else "None.", "",
             "## Results", "",
             sprintf("Input: `%s`. %d curve(s), %d points, %d used in fits (%d recovery points and %d flagged points excluded).",
                     input, nrow(p), nrow(pts), sum(pts$used_in_fit),
                     sum(grepl("recovery_point", pts$excluded_reason)),
                     sum(grepl("qc_hard_fail|unstable|missing", pts$excluded_reason))), "",
             "Parameters in umol m-2 s-1, normalized to 25 C with the temperature responses below; leaf-temperature values and confidence limits are in `tables/aci_parameters.csv`. NA means PhotoGEA judged the estimate unreliable (the process never limited, or the upper confidence limit was infinite).", "",
             fmt_table(show), "",
             paste0("![A-Ci fits](figures/", figs[1L], ")"), "",
             "## Method", "",
             sprintf("Pathway %s. Fitted with PhotoGEA %s `%s` by maximum likelihood; optimizer %s%s; seed %d (curve i uses seed + i - 1). Confidence intervals are PhotoGEA relative-likelihood intervals (threshold 0.147, about 95%%). Temperature responses: %s.%s",
                     opt$pathway, as.character(utils::packageVersion("PhotoGEA")),
                     if (opt$pathway == "C3") "fit_c3_aci" else "fit_c4_aci", opt$optimizer,
                     if (opt$optimizer == "deoptim") sprintf(" (%d generations)", as.integer(opt$generations)) else "",
                     as.integer(opt$seed), tp_txt,
                     if (opt$pathway == "C3") sprintf(" alpha_old: %s. Jmax from J with `calculate_jmax` (Bernacchi et al. 2003 alpha_J and theta_J), using mean incident PPFD; Jmax depends on these assumed values.", opt$alpha_old) else ""),
             "",
             sprintf("Repeated points at the starting CO2 after the curve returned to it (recovery points) were %s. The induction ratio is A at the last recovery point divided by A at the last starting point; flagged outside 1 +/- %g.",
                     if (opt$keep_recovery_points) "kept in the fits" else "excluded from the fits", opt$induction_tol),
             "",
             sprintf("QC rules of thumb: at least %d points below Ci %g and above Ci %g; at least %d points per limitation; gsw at the highest CO2 >= %g x its starting value; RMSE <= %g x max A.",
                     opt$min_points, low, high, opt$min_points, opt$closure_ratio, opt$rmse_frac), "",
             "## Limitations", "",
             "- C3 parameters are apparent (infinite gm). Estimates of J need RuBP-regeneration-limited points and TPU needs TPU-limited points; many curves support neither.",
             "- Parameters normalized to 25 C inherit any error in the assumed temperature responses, which grows with distance from 25 C.",
             "- The FvCB model assumes steady state at every point; the induction ratio checks this only at the reference CO2.", "",
             "## References", "",
             "- Farquhar GD, von Caemmerer S, Berry JA (1980) A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species. Planta 149: 78-90. https://doi.org/10.1007/BF00386231",
             "- Sharkey TD, Bernacchi CJ, Farquhar GD, Singsaas EL (2007) Fitting photosynthetic carbon dioxide response curves for C3 leaves. Plant, Cell & Environment 30: 1035-1040. https://doi.org/10.1111/j.1365-3040.2007.01710.x",
             "- Bernacchi CJ, Singsaas EL, Pimentel C, Portis AR Jr, Long SP (2001) Improved temperature response functions for models of Rubisco-limited photosynthesis. Plant, Cell & Environment 24: 253-259. https://doi.org/10.1111/j.1365-3040.2001.00668.x",
             "- Bernacchi CJ, Pimentel C, Long SP (2003) In vivo temperature response functions of parameters required to model RuBP-limited photosynthesis. Plant, Cell & Environment 26: 1419-1430. https://doi.org/10.1046/j.0016-8025.2003.01050.x",
             "- Long SP, Bernacchi CJ (2003) Gas exchange measurements, what can they tell us about the underlying limitations to photosynthesis? Journal of Experimental Botany 54: 2393-2401. https://doi.org/10.1093/jxb/erg262",
             "- von Caemmerer S (2000) Biochemical Models of Leaf Photosynthesis. CSIRO Publishing. https://doi.org/10.1071/9780643103405",
             "- von Caemmerer S (2021) Updating the steady-state model of C4 photosynthesis. Journal of Experimental Botany 72: 6003-6017. https://doi.org/10.1093/jxb/erab266",
             "- Lochocki EB, Salesse-Smith CE, McGrath JM (2025) PhotoGEA: an R package for closer fitting of photosynthetic gas exchange data with non-Gaussian confidence interval estimation. Plant, Cell & Environment 48: 5104-5119. https://doi.org/10.1111/pce.15501",
             "", sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(lines, file.path(out_dir, "report.md"))
  invisible(w)
}

# Full analysis on a tidy data frame; returns tables (used by main and tests).
analyze_aci <- function(df, opt) {
  curve_cols <- split_cols(opt$curve_cols)
  validate_input(df, curve_cols, opt)
  prep <- prepare_points(df, opt, curve_cols)
  pts <- prep$points
  n_used <- tapply(pts$used_in_fit, pts$curve_id, sum)
  min_used <- if (opt$pathway == "C3") 6L else 5L
  too_few <- names(n_used)[n_used < min_used]
  if (length(too_few) == length(n_used))
    stop("No curve has at least ", min_used, " usable points; nothing to fit.", call. = FALSE)
  # Rows without A, Ci or leaf temperature cannot enter PhotoGEA; they stay in
  # point_use.csv with their exclusion reason.
  fit_pts <- pts[!pts$curve_id %in% too_few & !grepl("missing_", pts$excluded_reason), , drop = FALSE]
  x <- build_exdf(fit_pts, opt)
  fr <- run_fits(x, opt)
  cr <- collect_results(fr, opt)
  params <- cr$params
  if (length(too_few)) {
    params <- rbind_fill(list(params, data.frame(curve_id = too_few, fit_error = sprintf(
      "fewer than %d usable points", min_used), stringsAsFactors = FALSE)))
  }
  qc <- curve_qc(pts, params, prep$induction, opt)
  params <- merge(params, qc, by = "curve_id", all.x = TRUE, sort = FALSE)
  meta <- curve_metadata(pts, curve_cols)
  params <- merge(meta, params, by = "curve_id", all.y = TRUE, sort = FALSE)
  params <- params[order(params$curve_id), , drop = FALSE]
  rownames(params) <- NULL
  list(params = params, fits = cr$fits, interp = cr$interp, points = pts)
}

check_options <- function(opt) {
  opt$pathway <- toupper(opt$pathway)
  if (!opt$pathway %in% c("C3", "C4")) stop("--pathway must be C3 or C4.", call. = FALSE)
  opt$temp_params <- tolower(opt$temp_params)
  if (!opt$temp_params %in% c("bernacchi", "sharkey", "flat"))
    stop("--temp-params must be bernacchi, sharkey or flat.", call. = FALSE)
  if (!opt$optimizer %in% c("deoptim", "nmkb")) stop("--optimizer must be deoptim or nmkb.", call. = FALSE)
  if (!identical(opt$alpha_old, "fit") && is.na(suppressWarnings(as.numeric(opt$alpha_old))))
    stop("--alpha-old must be 'fit' or a number between 0 and 1.", call. = FALSE)
  if (opt$cores < 1) stop("--cores must be at least 1.", call. = FALSE)
  if (opt$generations < 10) stop("--generations must be at least 10.", call. = FALSE)
  opt
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- check_options(parse_args(argv))
  if (!requireNamespace("PhotoGEA", quietly = TRUE))
    stop("The R package PhotoGEA is required: install.packages('PhotoGEA').", call. = FALSE)
  if (opt$demo) {
    opt$input <- file.path(SKILL_DIR, "demo", if (opt$pathway == "C3") "aci_demo_c3.csv" else "aci_demo_c4.csv")
    if (!nzchar(opt$curve_cols)) opt$curve_cols <- "genotype,plot"
  }
  if (!nzchar(opt$input)) stop("Provide --input FILE.csv (claw-licor-6800 output), or --demo. See --help.", call. = FALSE)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

  df <- read_tidy(opt$input)
  res <- analyze_aci(df, opt)
  save_table(res$params, opt$out, "aci_parameters.csv")
  if (!is.null(res$fits)) save_table(res$fits, opt$out, "aci_fitted_points.csv")
  save_table(res$points[, intersect(c("curve_id", "source_file", "instrument_id", "obs", "ci_umol_mol",
                                      "a_umol_m2_s", "recovery_point", "used_in_fit", "excluded_reason"),
                                    names(res$points))], opt$out, "point_use.csv")
  figs <- plot_fits(res$points, res$interp, res$params, opt, opt$out)
  plot_residuals(res$fits, opt$out)
  w <- write_report(res, opt, opt$input, opt$out, figs)
  write_repro_bundle(opt$out, argv, inputs = opt$input)
  cat(sprintf("Done: %d curve(s), %d warning(s). Report: %s\n", nrow(res$params), length(w),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(res)
}

if (sys.nframe() == 0L) main()
