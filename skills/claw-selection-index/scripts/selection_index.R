#!/usr/bin/env Rscript
# claw-selection-index: response to selection (breeder's equation) and
# multi-trait selection indices (Smith-Hazel, Pesek-Baker desired gains,
# Kempthorne-Nordskog restricted index) from phenotypic (P) and genetic (G)
# covariance matrices.
#
# Usage:
#   Rscript scripts/selection_index.R --demo
#   Rscript scripts/selection_index.R --p-matrix P.csv --g-matrix G.csv \
#       --weights "yield_kg_ha=1,protein_pct=60" --proportion 0.05
#   Rscript scripts/selection_index.R --scenarios scenarios.csv
#   Rscript scripts/selection_index.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-selection-index"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/selection_index.R [--demo] [--p-matrix P.csv --g-matrix G.csv]
            [--g-matrix G.csv --e-matrix E.csv --ge-matrix GE.csv --n-env E --n-rep R]
            [--weights W] [--desired D] [--restrict TRAITS] [--proportion P]
            [--intensity I] [--n-candidates N] [--cycle-years L]
            [--scenarios FILE.csv] [--bend] [--bend-eps E] [--out DIR]

  --demo          run on the demo files in demo/ (all three indices + scenarios)
  --p-matrix      phenotypic covariance matrix CSV (first column = trait names,
                  header = the same trait names, trait units squared)
  --g-matrix      genetic (additive) covariance matrix CSV, same layout and order
  --e-matrix      instead of --p-matrix: residual (plot-level) covariance matrix from a
                  multi-environment analysis; P is built for entry means as
                  G + GE/n_env + E/(n_env n_rep)
  --ge-matrix     genotype x environment covariance matrix (needed if --n-env > 1)
  --n-env, --n-rep  environments and replicates per environment in the entry means
  --weights       economic weights per trait unit: 'trait=value,...' or a CSV with
                  columns trait,weight (Smith-Hazel index; also used by --restrict)
  --desired       desired gains per cycle in trait units: 'trait=value,...' or a
                  CSV with columns trait,desired_gain (Pesek-Baker index)
  --restrict      comma-separated traits whose genetic change is held at zero
                  (Kempthorne-Nordskog restricted index; needs --weights)
  --proportion    proportion of candidates selected, 0-1 (default 0.05)
  --intensity     selection intensity i; overrides --proportion (default: from p)
  --n-candidates  number of candidates; if > 0, i is the exact expected mean of the
                  best round(p*N) of N normal values (small-population correction)
  --cycle-years   years per breeding cycle, for gain per year (default: not used)
  --scenarios     CSV of breeder's-equation scenarios (see SKILL.md for columns)
  --bend          if G is not positive semi-definite, bend it (eigenvalue floor,
                  Jorjani et al. 2003) instead of stopping
  --bend-eps      eigenvalue floor used by --bend (default 1e-4)
  --out           output folder (default: agriclaw_output/claw-selection-index)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, p_matrix = "", g_matrix = "", e_matrix = "",
                 ge_matrix = "", n_env = 1, n_rep = 1, weights = "",
                 desired = "", restrict = "", proportion = 0.05, intensity = 0,
                 n_candidates = 0, cycle_years = 0, scenarios = "",
                 bend = FALSE, bend_eps = 1e-4, out = "")

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
# Selection intensity and the breeder's equation
# ---------------------------------------------------------------------------

# Truncation selection on a normal variable in an infinite population:
# i = phi(z) / p, with z the standardised truncation point (Falconer &
# Mackay 1996, ch. 11).
intensity_infinite <- function(p) {
  if (any(!is.finite(p) | p <= 0 | p >= 1))
    stop("Proportion selected must be strictly between 0 and 1.", call. = FALSE)
  stats::dnorm(stats::qnorm(p, lower.tail = FALSE)) / p
}

# Exact expected selection intensity when the best k of N candidates are kept
# (expected mean of the top k of N independent standard normal values; the
# quantity Burrows 1972 approximates). Uses
#   E[sum of top k] = N * integral x phi(x) P(Bin(N - 1, 1 - Phi(x)) <= k - 1) dx.
intensity_finite <- function(k, N) {
  k <- as.integer(round(k))
  N <- as.integer(round(N))
  if (is.na(N) || is.na(k) || N < 2L || k < 1L || k >= N)
    stop("Finite-population intensity needs 1 <= k < N and N >= 2 (got k = ", k,
         ", N = ", N, ").", call. = FALSE)
  f <- function(x) x * stats::dnorm(x) *
    stats::pbinom(k - 1L, N - 1L, stats::pnorm(x, lower.tail = FALSE))
  s <- stats::integrate(f, -12, 12, subdivisions = 2000L, rel.tol = 1e-10)$value
  N * s / k
}

# Selection intensity from the options: explicit i, or from p (and N).
selection_intensity <- function(proportion, n_candidates = 0, intensity = 0) {
  if (intensity > 0) return(list(i = intensity, k = NA_integer_,
                                 method = "given by --intensity"))
  if (n_candidates > 0) {
    k <- max(1L, as.integer(round(proportion * n_candidates)))
    return(list(i = intensity_finite(k, n_candidates), k = k,
                method = sprintf("exact expected intensity, best %d of %d", k,
                                 as.integer(n_candidates))))
  }
  list(i = intensity_infinite(proportion), k = NA_integer_,
       method = "truncation selection, infinite population: i = phi(z)/p")
}

# Response to selection per cycle and per year.
#   R = i h2 sigma_P (Falconer & Mackay 1996)
#   Delta G = i r sigma_A / L (Cobb et al. 2019)
breeders_equation <- function(i, h2 = NA, sigma_p = NA, accuracy = NA,
                              sigma_a = NA, cycle_years = NA) {
  chk <- function(x, lo, hi, what)
    if (!is.na(x) && (x <= lo || x > hi))
      stop(what, " must be in (", lo, ", ", hi, "], got ", x, ".", call. = FALSE)
  chk(h2, 0, 1, "h2")
  chk(accuracy, 0, 1, "accuracy")
  if (!is.na(sigma_p) && sigma_p <= 0) stop("sigma_p must be > 0.", call. = FALSE)
  if (!is.na(sigma_a) && sigma_a <= 0) stop("sigma_a must be > 0.", call. = FALSE)
  if (!is.na(cycle_years) && cycle_years <= 0)
    stop("cycle_years must be > 0.", call. = FALSE)
  r_pheno <- if (!is.na(h2) && !is.na(sigma_p)) i * h2 * sigma_p else NA_real_
  r_acc <- if (!is.na(accuracy) && !is.na(sigma_a)) i * accuracy * sigma_a else NA_real_
  if (is.na(r_pheno) && is.na(r_acc))
    stop("Each scenario needs h2 + sigma_p, or accuracy + sigma_a.", call. = FALSE)
  resp <- if (!is.na(r_acc)) r_acc else r_pheno
  data.frame(response_per_cycle_phenotypic = r_pheno,
             response_per_cycle_accuracy = r_acc,
             gain_per_year = if (is.na(cycle_years)) NA_real_ else resp / cycle_years)
}

# Phenotypic variance and heritability of entry means over n_env environments
# and n_rep replicates: var_P = var_G + var_GE/e + var_e/(r e), h2 = var_G/var_P
# (Holland et al. 2002). The same formula, applied to covariance matrices,
# gives the entry-mean P matrix.
entry_mean_h2 <- function(var_g, var_e, n_rep, var_ge = 0, n_env = 1) {
  if (any(c(var_g, var_e, var_ge) < 0) || n_rep < 1 || n_env < 1)
    stop("Variance components must be >= 0 and n_rep, n_env >= 1.", call. = FALSE)
  var_p <- var_g + var_ge / n_env + var_e / (n_env * n_rep)
  c(h2 = var_g / var_p, sigma_p = sqrt(var_p))
}

entry_mean_P <- function(G, E, GE = NULL, n_env = 1, n_rep = 1) {
  if (n_env < 1 || n_rep < 1) stop("--n-env and --n-rep must be at least 1.", call. = FALSE)
  P <- G + E / (n_env * n_rep)
  if (!is.null(GE)) P <- P + GE / n_env
  P
}

run_scenarios <- function(sc) {
  need <- c("scenario", "proportion_selected")
  miss <- setdiff(need, names(sc))
  if (length(miss))
    stop("Scenario file is missing column(s): ", paste(miss, collapse = ", "),
         ". Columns: scenario, proportion_selected, and h2 + sigma_p and/or ",
         "accuracy + sigma_a (or var_g, var_ge, var_e, n_env, n_rep); optional ",
         "n_candidates, cycle_years.", call. = FALSE)
  for (col in c("n_candidates", "h2", "sigma_p", "accuracy", "sigma_a", "cycle_years",
                "var_g", "var_ge", "var_e", "n_env", "n_rep"))
    if (!col %in% names(sc)) sc[[col]] <- NA_real_
  # Entry-mean h2 and sigma_P from variance components where h2/sigma_p are empty.
  for (r in seq_len(nrow(sc))) {
    if (is.na(sc$h2[r]) && is.na(sc$sigma_p[r]) && !is.na(sc$var_g[r])) {
      if (is.na(sc$var_e[r]) || is.na(sc$n_rep[r]))
        stop("Scenario '", sc$scenario[r], "': variance components need var_g, var_e and n_rep ",
             "(plus var_ge and n_env for several environments).", call. = FALSE)
      if (xor(is.na(sc$var_ge[r]), is.na(sc$n_env[r])))
        stop("Scenario '", sc$scenario[r], "': give var_ge and n_env together.", call. = FALSE)
      hm <- entry_mean_h2(sc$var_g[r], sc$var_e[r], sc$n_rep[r],
                          if (is.na(sc$var_ge[r])) 0 else sc$var_ge[r],
                          if (is.na(sc$n_env[r])) 1 else sc$n_env[r])
      sc$h2[r] <- hm[["h2"]]
      sc$sigma_p[r] <- hm[["sigma_p"]]
    }
  }
  rows <- lapply(seq_len(nrow(sc)), function(r) {
    s <- sc[r, ]
    nc <- if (is.na(s$n_candidates)) 0 else s$n_candidates
    si <- selection_intensity(s$proportion_selected, nc)
    be <- breeders_equation(si$i, s$h2, s$sigma_p, s$accuracy, s$sigma_a,
                            s$cycle_years)
    data.frame(scenario = s$scenario, proportion_selected = s$proportion_selected,
               n_candidates = s$n_candidates, n_selected = si$k,
               intensity = si$i, intensity_method = si$method,
               h2 = s$h2, sigma_p = s$sigma_p, accuracy = s$accuracy,
               sigma_a = s$sigma_a, cycle_years = s$cycle_years, be,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Matrices: reading and checks
# ---------------------------------------------------------------------------

read_matrix <- function(path, what) {
  if (!file.exists(path)) stop(what, " file not found: ", path, call. = FALSE)
  df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(df) < 2L || nrow(df) < 1L)
    stop(what, " must have a first column of trait names and one column per trait.",
         call. = FALSE)
  rn <- as.character(df[[1L]])
  M <- as.matrix(df[, -1L, drop = FALSE])
  if (!is.numeric(M))
    stop(what, " has non-numeric entries; it must be a numeric matrix.", call. = FALSE)
  if (nrow(M) != ncol(M))
    stop(what, " is ", nrow(M), " x ", ncol(M), "; it must be square.", call. = FALSE)
  if (anyNA(M)) stop(what, " has missing values.", call. = FALSE)
  cn <- colnames(M)
  if (anyDuplicated(rn) || !setequal(rn, cn))
    stop(what, ": row names (first column) and column names must be the same ",
         "unique trait names. Rows: ", paste(rn, collapse = ", "),
         "; columns: ", paste(cn, collapse = ", "), call. = FALSE)
  rownames(M) <- rn
  M[cn, cn, drop = FALSE]
}

symmetrise <- function(M, what) {
  asym <- max(abs(M - t(M)))
  if (asym > 1e-6 * max(abs(M)))
    stop(what, " is not symmetric (largest |M - t(M)| = ", signif(asym, 3),
         "). Check that it is a covariance matrix.", call. = FALSE)
  (M + t(M)) / 2
}

# Unweighted bending: eigenvalues below eps are set to eps and the matrix is
# rebuilt (Jorjani et al. 2003 with a unit weight matrix; Nilforooshan 2020).
bend_matrix <- function(M, eps = 1e-4) {
  e <- eigen(M, symmetric = TRUE)
  d <- pmax(e$values, eps)
  B <- e$vectors %*% diag(d, length(d)) %*% t(e$vectors)
  B <- (B + t(B)) / 2
  dimnames(B) <- dimnames(M)
  B
}

# Deviation statistics between a matrix and its bent version (upper triangle
# including the diagonal), as reported by the mbend package.
bend_stats <- function(M, B) {
  ut <- upper.tri(M, diag = TRUE)
  dv <- (B - M)[ut]
  data.frame(min_dev = min(dv), max_dev = max(dv), mean_dev = mean(dv),
             aad = mean(abs(dv)), correlation = stats::cor(B[ut], M[ut]),
             rmsd = sqrt(mean(dv^2)))
}

check_matrices <- function(P, G, bend = FALSE, bend_eps = 1e-4) {
  if (!identical(rownames(P), rownames(G))) {
    if (!setequal(rownames(P), rownames(G)))
      stop("P and G must contain the same traits. P: ",
           paste(rownames(P), collapse = ", "), "; G: ",
           paste(rownames(G), collapse = ", "), call. = FALSE)
    G <- G[rownames(P), rownames(P), drop = FALSE]
  }
  P <- symmetrise(P, "P")
  G <- symmetrise(G, "G")
  notes <- character()
  if (any(diag(P) <= 0)) stop("P has non-positive variances on its diagonal.", call. = FALSE)
  ev_p <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev_p) <= 1e-10 * max(ev_p))
    stop("P is not positive definite (smallest eigenvalue ", signif(min(ev_p), 3),
         "), so it cannot be inverted. Check for duplicated or linearly dependent ",
         "traits, or a P matrix assembled from different sources.", call. = FALSE)
  ev_g <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  tol <- 1e-8 * max(abs(ev_g))
  bent <- NULL
  if (min(ev_g) < -tol) {
    if (!bend)
      stop("G is not positive semi-definite (eigenvalues: ",
           paste(signif(ev_g, 4), collapse = ", "), "). This usually means sampling ",
           "error in the estimated genetic covariances, or G assembled from different ",
           "analyses. Re-estimate G jointly, or re-run with --bend to apply an ",
           "eigenvalue floor (Jorjani et al. 2003); bending changes G and is ",
           "reported.", call. = FALSE)
    G_orig <- G
    G <- bend_matrix(G, bend_eps)
    bent <- cbind(bend_stats(G_orig, G), eps = bend_eps)
    notes <- c(notes, sprintf(paste0(
      "G was not positive semi-definite (smallest eigenvalue %.4g); it was bent ",
      "with an eigenvalue floor of %g (Jorjani et al. 2003). Largest change in ",
      "any element: %.4g."), min(ev_g), bend_eps, max(abs(G - G_orig))))
    ev_g <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  }
  # Canonical heritabilities: eigenvalues of P^-1 G. All must lie in [0, 1]
  # if P and G describe the same population (P - G is then a covariance matrix).
  can_h2 <- Re(eigen(solve(P, G), only.values = TRUE)$values)
  h2 <- diag(G) / diag(P)
  if (any(diag(G) <= 0))
    notes <- c(notes, paste("Zero genetic variance for:",
                            paste(rownames(G)[diag(G) <= 0], collapse = ", ")))
  sd_g <- sqrt(pmax(diag(G), 0))
  rg <- G / outer(sd_g, sd_g)
  rg[!is.finite(rg)] <- NA_real_
  list(P = P, G = G, ev_p = ev_p, ev_g = ev_g, can_h2 = sort(can_h2, decreasing = TRUE),
       h2 = h2, rg = rg, bent = bent, notes = notes,
       ev_e = eigen(P - G, symmetric = TRUE, only.values = TRUE)$values)
}

# Parse 'trait=value,trait=value' or a two-column CSV (trait, value).
read_trait_values <- function(spec, value_col, traits, what) {
  if (!nzchar(spec)) return(NULL)
  if (file.exists(spec)) {
    df <- utils::read.csv(spec, check.names = FALSE, stringsAsFactors = FALSE)
    if (!all(c("trait", value_col) %in% names(df)))
      stop(what, " file needs columns 'trait' and '", value_col, "'.", call. = FALSE)
    v <- stats::setNames(as.numeric(df[[value_col]]), as.character(df$trait))
  } else {
    parts <- strsplit(strsplit(spec, ",", fixed = TRUE)[[1]], "=", fixed = TRUE)
    if (any(lengths(parts) != 2L))
      stop(what, " must look like 'trait=value,trait=value' or be a CSV file; got '",
           spec, "'.", call. = FALSE)
    v <- stats::setNames(suppressWarnings(as.numeric(vapply(parts, `[`, "", 2L))),
                         trimws(vapply(parts, `[`, "", 1L)))
  }
  if (anyNA(v)) stop(what, " contains non-numeric values.", call. = FALSE)
  extra <- setdiff(names(v), traits)
  if (length(extra))
    stop(what, " names trait(s) not in P and G: ", paste(extra, collapse = ", "),
         ". Traits in P and G: ", paste(traits, collapse = ", "), call. = FALSE)
  miss <- setdiff(traits, names(v))
  if (length(miss))
    stop(what, " has no value for: ", paste(miss, collapse = ", "),
         ". Give every trait in P and G a value (use 0 explicitly for a trait ",
         "with no economic value or no desired change).", call. = FALSE)
  v[traits]
}

# ---------------------------------------------------------------------------
# Selection indices
# ---------------------------------------------------------------------------

# Properties of any index I = b'x for selection intensity i.
#   sigma_I = sqrt(b'Pb); expected genetic gain per trait = i G b / sigma_I;
#   with weights w (H = w'g): r_IH = b'Gw / (sigma_I sigma_H), response in H =
#   i b'Gw / sigma_I (Hazel 1943; Ceron-Rojas & Crossa 2018, ch. 2).
index_metrics <- function(b, P, G, i, w = NULL) {
  b <- drop(b)
  s_i <- sqrt(drop(crossprod(b, P %*% b)))
  if (!is.finite(s_i) || s_i <= 0) stop("Index variance b'Pb is not positive.", call. = FALSE)
  gains <- i * drop(G %*% b) / s_i
  # Floating-point dust (e.g. 1e-16 for a restricted trait) is reported as 0.
  gains[abs(gains) < 1e-10 * max(abs(gains))] <- 0
  out <- list(b = stats::setNames(b, rownames(P)), sigma_I = s_i,
              gains = stats::setNames(gains, rownames(P)),
              sigma_H = NA_real_, r_IH = NA_real_, response_H = NA_real_)
  if (!is.null(w)) {
    s_h <- sqrt(drop(crossprod(w, G %*% w)))
    cov_hi <- drop(crossprod(b, G %*% w))
    out$sigma_H <- s_h
    out$r_IH <- cov_hi / (s_i * s_h)
    out$response_H <- i * cov_hi / s_i
  }
  out
}

# Smith (1936) - Hazel (1943) optimum index: b = P^-1 G w.
smith_hazel <- function(P, G, w) drop(solve(P, G %*% w))

# Pesek & Baker (1969) desired-gains index in the general form of Itoh &
# Yamada (1986): b = P^-1 G (G P^-1 G)^-1 d, which equals G^-1 d when G is
# invertible. Expected gains are proportional to d.
desired_gains_index <- function(P, G, d) {
  PiG <- solve(P, G)
  M <- G %*% PiG
  if (rcond(M) < 1e-12)
    stop("G is singular, so the desired-gains index is not defined. Check for ",
         "traits with zero genetic variance or perfectly correlated traits.", call. = FALSE)
  drop(PiG %*% solve(M, d))
}

# Kempthorne & Nordskog (1959) restricted index: b_R = [I - P^-1 C (C'P^-1 C)^-1 C'] P^-1 G w,
# with C the columns of G for the restricted traits; genetic change in those
# traits is zero.
restricted_index <- function(P, G, w, restrict) {
  C <- G[, restrict, drop = FALSE]
  PiC <- solve(P, C)
  CtPiC <- crossprod(C, PiC)
  if (rcond(CtPiC) < 1e-12)
    stop("The restricted traits have (near-)zero or collinear genetic variance; ",
         "the restricted index is not defined.", call. = FALSE)
  K <- diag(nrow(P)) - PiC %*% solve(CtPiC, t(C))
  drop(K %*% smith_hazel(P, G, w))
}

index_rows <- function(name, m, P, G, cycle_years) {
  data.frame(index = name, trait = names(m$gains), b = unname(m$b),
             expected_gain = unname(m$gains),
             gain_per_year = if (cycle_years > 0) unname(m$gains) / cycle_years else NA_real_,
             gain_in_sd_g = unname(m$gains) / sqrt(diag(G)),
             gain_in_sd_p = unname(m$gains) / sqrt(diag(P)),
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

fmt <- function(x, d = 4) ifelse(is.na(x), "", formatC(x, digits = d, format = "g"))

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
  if (opt$demo) {
    dd <- file.path(SKILL_DIR, "demo")
    opt$p_matrix <- file.path(dd, "P_matrix.csv")
    opt$g_matrix <- file.path(dd, "G_matrix.csv")
    opt$weights <- file.path(dd, "weights.csv")
    opt$desired <- file.path(dd, "desired_gains.csv")
    opt$restrict <- "height_cm"
    opt$scenarios <- file.path(dd, "scenarios.csv")
  }
  have_mats <- nzchar(opt$p_matrix) || nzchar(opt$g_matrix) || nzchar(opt$e_matrix)
  if (!have_mats && !nzchar(opt$scenarios))
    stop("Provide --p-matrix and --g-matrix (selection index), --scenarios ",
         "(breeder's equation), or --demo. See --help.", call. = FALSE)
  if (nzchar(opt$p_matrix) && nzchar(opt$e_matrix))
    stop("Give either --p-matrix, or --e-matrix (+ --ge-matrix) to build P, not both.", call. = FALSE)
  if (have_mats && !(nzchar(opt$g_matrix) && (nzchar(opt$p_matrix) || nzchar(opt$e_matrix))))
    stop("A selection index needs --g-matrix and either --p-matrix or --e-matrix.", call. = FALSE)
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

  warn <- character()
  body <- character()
  inputs <- c(opt$p_matrix, opt$g_matrix, opt$e_matrix, opt$ge_matrix, opt$scenarios,
              if (file.exists(opt$weights)) opt$weights,
              if (file.exists(opt$desired)) opt$desired)
  res <- list()

  # --- breeder's equation scenarios ---------------------------------------
  if (nzchar(opt$scenarios)) {
    if (!file.exists(opt$scenarios)) stop("Scenario file not found: ", opt$scenarios, call. = FALSE)
    sc <- utils::read.csv(opt$scenarios, check.names = FALSE, stringsAsFactors = FALSE)
    be <- run_scenarios(sc)
    save_table(be, opt$out, "breeders_equation.csv")
    res$scenarios <- be
    if (any(!is.na(be$n_candidates) & be$n_candidates > 0))
      warn <- c(warn, paste0("Some scenarios use a small candidate population; their ",
                             "intensity is the exact expected value for that N, which is ",
                             "lower than the infinite-population table value."))
    open_png(opt$out, "breeders_equation.png")
    y <- if (all(is.na(be$gain_per_year))) ifelse(is.na(be$response_per_cycle_accuracy),
           be$response_per_cycle_phenotypic, be$response_per_cycle_accuracy) else be$gain_per_year
    graphics::par(mar = c(5, 12, 3, 1))
    graphics::barplot(y, names.arg = be$scenario, horiz = TRUE, las = 1,
                      col = "steelblue", xlab = if (all(is.na(be$gain_per_year)))
                        "Response per cycle (trait units)" else "Genetic gain per year (trait units)",
                      main = "Breeder's equation scenarios")
    grDevices::dev.off()
    body <- c(body, "## Breeder's equation", "",
              "R = i h2 sigma_P per cycle (Falconer & Mackay 1996); Delta G = i r sigma_A / L per year (Cobb et al. 2019).", "",
              md_table(data.frame(scenario = be$scenario, p = fmt(be$proportion_selected),
                                  N = fmt(be$n_candidates), i = fmt(be$intensity),
                                  R_cycle_h2 = fmt(be$response_per_cycle_phenotypic),
                                  R_cycle_r = fmt(be$response_per_cycle_accuracy),
                                  gain_per_year = fmt(be$gain_per_year))),
              "", "![Scenarios](figures/breeders_equation.png)", "")
  }

  # --- selection indices ----------------------------------------------------
  if (have_mats) {
    G0 <- read_matrix(opt$g_matrix, "G matrix")
    if (nzchar(opt$p_matrix)) {
      P0 <- read_matrix(opt$p_matrix, "P matrix")
    } else {
      E0 <- read_matrix(opt$e_matrix, "E matrix")
      GE0 <- if (nzchar(opt$ge_matrix)) read_matrix(opt$ge_matrix, "GE matrix") else NULL
      if (opt$n_env > 1 && is.null(GE0))
        stop("--n-env > 1 needs --ge-matrix (genotype x environment covariances); without it ",
             "P would ignore G x E and overstate heritability.", call. = FALSE)
      for (M in list(E0, GE0)) if (!is.null(M) && !setequal(rownames(M), rownames(G0)))
        stop("G, E and GE matrices must contain the same traits.", call. = FALSE)
      E0 <- E0[rownames(G0), rownames(G0), drop = FALSE]
      if (!is.null(GE0)) GE0 <- GE0[rownames(G0), rownames(G0), drop = FALSE]
      P0 <- entry_mean_P(G0, E0, GE0, opt$n_env, opt$n_rep)
      save_table(data.frame(trait = rownames(P0), P0, check.names = FALSE), opt$out,
                 "P_entry_means.csv")
      warn <- c(warn, sprintf(paste0(
        "P was built for entry means over %g environment(s) x %g replicate(s): P = G + GE/%g + ",
        "E/%g (Holland et al. 2002). Check that this matches how candidates are evaluated."),
        opt$n_env, opt$n_rep, opt$n_env, opt$n_env * opt$n_rep))
    }
    ck <- check_matrices(P0, G0, opt$bend, opt$bend_eps)
    P <- ck$P
    G <- ck$G
    traits <- rownames(P)
    warn <- c(warn, ck$notes)
    if (any(ck$can_h2 > 1 + 1e-8))
      warn <- c(warn, sprintf(paste0(
        "P - G is not positive semi-definite: the largest canonical heritability ",
        "(eigenvalue of P^-1 G) is %.3f > 1. P and G are inconsistent - they may not ",
        "come from the same population or analysis. Index results are unreliable."),
        max(ck$can_h2)))
    if (any(ck$h2 > 1 + 1e-8))
      warn <- c(warn, paste("Heritability above 1 for:", paste(traits[ck$h2 > 1], collapse = ", ")))
    hi_rg <- which(!is.na(ck$rg) & abs(ck$rg) > 0.95 & upper.tri(ck$rg), arr.ind = TRUE)
    if (nrow(hi_rg))
      warn <- c(warn, paste0("Genetic correlation |r_G| > 0.95 between ",
                             paste(traits[hi_rg[, 1]], traits[hi_rg[, 2]], sep = " & ",
                                   collapse = "; "),
                             ": index weights will be unstable."))
    si <- selection_intensity(opt$proportion, opt$n_candidates, opt$intensity)
    i <- si$i
    w <- read_trait_values(opt$weights, "weight", traits, "Weights")
    d <- read_trait_values(opt$desired, "desired_gain", traits, "Desired gains")
    restrict <- if (nzchar(opt$restrict)) trimws(strsplit(opt$restrict, ",", fixed = TRUE)[[1]]) else character()
    if (length(restrict)) {
      bad <- setdiff(restrict, traits)
      if (length(bad)) stop("--restrict names trait(s) not in P and G: ",
                            paste(bad, collapse = ", "), call. = FALSE)
      if (is.null(w)) stop("--restrict needs --weights (the restricted index ",
                           "maximises the economic merit of the other traits).", call. = FALSE)
      if (length(restrict) >= length(traits))
        stop("Cannot restrict every trait.", call. = FALSE)
    }
    if (is.null(w) && is.null(d))
      stop("Give --weights (Smith-Hazel index) and/or --desired (desired-gains index).",
           call. = FALSE)

    gain_rows <- list()
    summ <- list()
    add <- function(name, m) {
      gain_rows[[length(gain_rows) + 1L]] <<- index_rows(name, m, P, G, opt$cycle_years)
      summ[[length(summ) + 1L]] <<- data.frame(
        index = name, sigma_I = m$sigma_I, r_IH = m$r_IH, response_H = m$response_H,
        response_H_per_year = if (opt$cycle_years > 0) m$response_H / opt$cycle_years else NA_real_,
        intensity = i, stringsAsFactors = FALSE)
      m
    }
    if (!is.null(w)) {
      res$smith_hazel <- add("smith_hazel", index_metrics(smith_hazel(P, G, w), P, G, i, w))
      res$base <- add("base_index_b_equals_w", index_metrics(w, P, G, i, w))
      if (length(restrict))
        res$restricted <- add(paste0("restricted_", paste(restrict, collapse = "+")),
                              index_metrics(restricted_index(P, G, w, restrict), P, G, i, w))
    }
    if (!is.null(d)) {
      if (all(d == 0)) stop("All desired gains are zero.", call. = FALSE)
      res$desired <- add("desired_gains", index_metrics(desired_gains_index(P, G, d), P, G, i, w))
      ratio <- res$desired$gains[d != 0] / d[d != 0]
      res$desired_cycles <- 1 / ratio[1]
    }
    gains <- do.call(rbind, gain_rows)
    summ <- do.call(rbind, summ)
    save_table(gains, opt$out, "index_gains.csv")
    save_table(summ, opt$out, "index_summary.csv")
    chk <- data.frame(trait = traits, var_p = diag(P), var_g = diag(G), h2 = ck$h2,
                      eigen_P = ck$ev_p, eigen_G = ck$ev_g, eigen_P_minus_G = ck$ev_e,
                      canonical_h2 = ck$can_h2)
    save_table(chk, opt$out, "input_checks.csv")
    rg_df <- data.frame(trait = traits, round(ck$rg, 4), check.names = FALSE)
    save_table(rg_df, opt$out, "genetic_correlations.csv")
    if (!is.null(ck$bent)) {
      save_table(ck$bent, opt$out, "bending_changes.csv")
      save_table(data.frame(trait = traits, G, check.names = FALSE), opt$out, "G_bent.csv")
    }
    res$checks <- ck
    res$intensity <- si

    # Figure: expected gains in genetic SD units, by index.
    mat <- do.call(cbind, lapply(split(gains$gain_in_sd_g, gains$index), function(x) x))
    rownames(mat) <- traits
    mat <- mat[, unique(gains$index), drop = FALSE]
    open_png(opt$out, "expected_gains.png", width = 8, height = 5)
    graphics::par(mar = c(8, 5, 3, 1))
    cols <- grDevices::hcl.colors(ncol(mat), "Dark 3")
    graphics::barplot(t(mat), beside = TRUE, col = cols, las = 2,
                      ylab = "Expected gain per cycle (genetic SD units)",
                      main = sprintf("Expected response per trait (i = %.3f)", i))
    graphics::abline(h = 0)
    graphics::legend("topright", legend = colnames(mat), fill = cols, bty = "n", cex = 0.8)
    grDevices::dev.off()

    body <- c(body, "## Inputs and checks", "",
              sprintf("Traits: %s.  ", paste(traits, collapse = ", ")),
              sprintf("Selection intensity i = %.4f (%s; proportion selected %s).",
                      i, si$method, fmt(opt$proportion)), "",
              md_table(data.frame(trait = traits, h2 = fmt(ck$h2, 3),
                                  sd_p = fmt(sqrt(diag(P))), sd_g = fmt(sqrt(diag(G))))),
              "",
              sprintf("Smallest eigenvalue of P: %.4g; of G: %.4g; of P - G: %.4g. Canonical heritabilities (eigenvalues of P^-1 G): %s.",
                      min(ck$ev_p), min(ck$ev_g), min(ck$ev_e),
                      paste(fmt(ck$can_h2, 3), collapse = ", ")), "",
              "## Selection indices", "",
              md_table(data.frame(index = summ$index, sigma_I = fmt(summ$sigma_I),
                                  r_IH = fmt(summ$r_IH, 3), response_H = fmt(summ$response_H))),
              "", "Expected genetic gain per cycle, trait units (b = index coefficients):", "",
              md_table(data.frame(index = gains$index, trait = gains$trait, b = fmt(gains$b),
                                  gain = fmt(gains$expected_gain),
                                  gain_sd_g = fmt(gains$gain_in_sd_g, 3))),
              "")
    if (!is.null(d))
      body <- c(body, sprintf(paste0(
        "Desired-gains index: expected gains are proportional to the desired gains d. ",
        "One cycle at this intensity gives %.3g x d, so reaching d takes about %.2g ",
        "cycle(s). This index ignores economic weights and does not maximise r_IH; ",
        "its r_IH and response_H above use the --weights merit, for comparison only."),
        1 / res$desired_cycles, res$desired_cycles), "")
    wd <- data.frame(trait = traits,
                     weight = if (is.null(w)) "" else fmt(w),
                     desired_gain = if (is.null(d)) "" else fmt(d),
                     restricted = ifelse(traits %in% restrict, "yes", ""))
    body <- c(body, "Economic weights, desired gains and restrictions used:", "",
              md_table(wd), "", "![Expected gains](figures/expected_gains.png)", "")
  }

  limits <- c("## Limitations", "",
              "- Expected responses assume P and G are known without error, multivariate normality, truncation selection on the index, and one generation of selection; realised responses are usually smaller because P and G are estimates.",
              "- Index results are only as good as P and G: both must describe the same population and the same unit of selection (for example, entry means over the same number of environments).",
              "- The desired-gains index fixes the direction of response, not its size, and does not maximise the correlation with economic merit.",
              "")
  refs <- c("## References", "",
            "- Falconer DS, Mackay TFC (1996) Introduction to Quantitative Genetics, 4th edn. Longman, Harlow. ISBN 978-0-582-24302-6.",
            "- Cobb JN et al. (2019) Theor Appl Genet 132:627-645. https://doi.org/10.1007/s00122-019-03317-0",
            "- Smith HF (1936) Ann Eugen 7:240-250. https://doi.org/10.1111/j.1469-1809.1936.tb02143.x",
            "- Hazel LN (1943) Genetics 28:476-490. https://doi.org/10.1093/genetics/28.6.476",
            "- Pesek J, Baker RJ (1969) Can J Plant Sci 49:803-804. https://doi.org/10.4141/cjps69-137",
            "- Kempthorne O, Nordskog AW (1959) Biometrics 15:10-19. https://doi.org/10.2307/2527598",
            "- Ceron-Rojas JJ, Crossa J (2018) Linear Selection Indices in Modern Plant Breeding. Springer. https://doi.org/10.1007/978-3-319-91223-3",
            "- Holland JB, Nyquist WE, Cervantes-Martinez CT (2002) Plant Breed Rev 22:9-112. https://doi.org/10.1002/9780470650202.ch2",
            "- Jorjani H, Klei L, Emanuelson U (2003) J Dairy Sci 86:677-679. https://doi.org/10.3168/jds.S0022-0302(03)73646-7")
  report <- c(sprintf("# %s report", SKILL), "",
              if (length(warn)) c("## Warnings", "", paste("-", warn), "") else NULL,
              body, limits, refs, "", sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(report, file.path(opt$out, "report.md"))
  write_repro_bundle(opt$out, argv, inputs = inputs)
  cat(sprintf("Done. Report: %s\n",
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(res)
}

if (sys.nframe() == 0L) main()
