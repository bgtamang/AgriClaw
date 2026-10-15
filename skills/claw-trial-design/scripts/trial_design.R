#!/usr/bin/env Rscript
# claw-trial-design: randomise a field trial with agricolae / FielDHub and lay it
# out on the field grid so that blocks are compact rectangles, whole plots are
# contiguous, strip-plot strips stay intact and every augmented block holds every
# check. Writes a field book, a field map, a design summary and design checks.
#
# Usage:
#   Rscript scripts/trial_design.R --demo
#   Rscript scripts/trial_design.R --design rcbd --entries 20 --reps 4 --seed 42
#   Rscript scripts/trial_design.R --design alpha --entry-list entries.csv --reps 3 --block-size 6
#   Rscript scripts/trial_design.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-trial-design"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/trial_design.R --design DESIGN [options]

Designs (package function used):
  crd        completely randomised (agricolae::design.crd)
  rcbd       randomised complete blocks (agricolae::design.rcbd)
  latin      Latin square(s) (agricolae::design.lsd)
  alpha      resolvable alpha / incomplete blocks (agricolae::design.alpha or
             FielDHub::alpha_lattice, whichever has the higher efficiency factor)
  rowcol     resolvable row-column (FielDHub::row_column)
  augmented  augmented RCBD, unreplicated tests + checks in every block (agricolae::design.dau)
  prep       partially replicated, spatially optimised (FielDHub::partially_replicated)
  diagonal   unreplicated with diagonal checks (FielDHub::diagonal_arrangement)
  split      split-plot in RCBD (agricolae::design.split)
  strip      strip-plot / split-block (agricolae::design.strip)
  factorial  factorial in RCBD (agricolae::design.ab)

Entries:
  --entries N        number of entries, named G01..GN
  --entry-list FILE  CSV with column 'entry'; optional 'is_check' (TRUE/FALSE) and,
                     for prep, 'reps' (copies per entry). Other columns are carried
                     into the field book.
  --checks N         add N check entries CHK1..CHKN (with --entries)
Design parameters:
  --reps N           replicates / complete blocks (squares for latin)
  --blocks N         number of blocks (augmented)
  --block-size K     plots per incomplete block (alpha; default: divisor of the
                     number of entries nearest its square root)
  --rep-rows N       rows in each replicate (rowcol; default: divisor nearest sqrt)
  --factor-a SPEC    'name:level1,level2,...' or 'name:3'. Whole-plot factor (split),
                     row strips (strip), first factor (factorial)
  --factor-b SPEC    sub-plot factor (split), column strips (strip), second factor
                     (factorial). For split and strip it defaults to the entries.
  --factor-c SPEC    third factor (factorial)
  --prep-fraction F  share of test entries planted twice (prep; default 0.25)
  --check-reps N     copies of each check (prep; default 4)
  --iterations N     optimisation iterations (rowcol; default 1000)
  --alpha-engine E   auto (default), agricolae or fieldhub
Field layout:
  --field-rows R     rows in the field grid (required for diagonal)
  --field-cols C     columns in the field grid
  --block-rows N     field rows spanned by each block (the rep for rcbd, split,
                     strip, factorial and rowcol; the incomplete block for alpha;
                     the block for augmented)
  --numbering S      plot numbering path: serpentine (default) or cartesian
  --plot-start N     first plot number (default 1)
General:
  --trial NAME       trial name used in plot_id (default TRIAL1)
  --seed N           random seed (default: drawn at random and recorded)
  --demo             alpha design for the 30 entries in demo/entries.csv
  --out DIR          output folder (default: agriclaw_output/claw-trial-design)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, design = "", entries = 0, entry_list = "", checks = 0,
                 reps = 0, blocks = 0, block_size = 0, rep_rows = 0,
                 factor_a = "", factor_b = "", factor_c = "",
                 prep_fraction = 0.25, check_reps = 4, iterations = 1000,
                 alpha_engine = "auto",
                 field_rows = 0, field_cols = 0, block_rows = 0,
                 numbering = "serpentine", plot_start = 1,
                 trial = "TRIAL1", seed = 0, out = "")

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
# Entries, factors and small utilities
# ---------------------------------------------------------------------------

DESIGNS <- c("crd", "rcbd", "latin", "alpha", "rowcol", "augmented", "prep",
             "diagonal", "split", "strip", "factorial")
ALIASES <- c(lsd = "latin", latinsquare = "latin", alphalattice = "alpha",
             lattice = "alpha", ibd = "alpha", rowcolumn = "rowcol",
             augmentedrcbd = "augmented", arcbd = "augmented", dau = "augmented",
             partiallyreplicated = "prep", diagonalarrangement = "diagonal",
             splitplot = "split", stripplot = "strip", splitblock = "strip",
             factorialrcbd = "factorial")
AGRICOLAE_KINDS <- "Super-Duper"   # agricolae's default RNG kind, recorded in the summary
RESERVED <- c("trial", "plot_id", "plot", "row", "col", "rep", "block", "design", "whole_plot",
              "entry_no", "entry", "is_check", "is_filler")

fail <- function(...) stop(paste0(...), call. = FALSE)

normalise_design <- function(x) {
  key <- gsub("[^a-z]", "", tolower(x))
  if (key %in% names(ALIASES)) key <- ALIASES[[key]]
  if (!key %in% DESIGNS)
    fail("Unknown --design '", x, "'. Choose one of: ", paste(DESIGNS, collapse = ", "), ".")
  key
}

parse_logical <- function(x, what) {
  if (is.logical(x)) return(!is.na(x) & x)
  v <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(v))
  out[v %in% c("TRUE", "T", "YES", "Y", "1", "CHECK")] <- TRUE
  out[v %in% c("FALSE", "F", "NO", "N", "0", "", "NA")] <- FALSE
  if (anyNA(out))
    fail("Column '", what, "' must be TRUE/FALSE (or yes/no, 1/0); found: ",
         paste(unique(x[is.na(out)]), collapse = ", "), ".")
  out
}

read_entry_list <- function(path) {
  if (!file.exists(path)) fail("Entry list not found: ", path)
  el <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                        strip.white = TRUE, na.strings = c("NA", ""))
  names(el) <- trimws(names(el))
  if (!"entry" %in% names(el))
    fail("The entry list needs a column named 'entry' (entry names). Columns found: ",
         paste(names(el), collapse = ", "), ".")
  el$entry <- trimws(as.character(el$entry))
  if (any(is.na(el$entry) | el$entry == ""))
    fail("The entry list has blank entry names in row(s) ",
         paste(which(is.na(el$entry) | el$entry == "") + 1L, collapse = ", "), " of the file.")
  dup <- unique(el$entry[duplicated(el$entry)])
  if (length(dup))
    fail("Entry names must be unique; duplicated: ", paste(utils::head(dup, 10), collapse = ", "), ".")
  el$is_check <- if ("is_check" %in% names(el)) parse_logical(el$is_check, "is_check") else FALSE
  if ("reps" %in% names(el)) {
    r <- suppressWarnings(as.numeric(el$reps))
    if (anyNA(r) || any(r < 1) || any(r != round(r)))
      fail("Column 'reps' must hold whole numbers >= 1 for every entry.")
    el$reps <- as.integer(r)
  }
  bad <- setdiff(intersect(names(el), RESERVED), c("entry", "is_check"))
  if (length(bad))
    fail("Entry-list column(s) ", paste(bad, collapse = ", "),
         " clash with field-book columns; rename them.")
  el
}

make_entries <- function(opt) {
  if (nzchar(opt$entry_list) && opt$entries > 0)
    fail("Give either --entries or --entry-list, not both.")
  if (nzchar(opt$entry_list)) {
    el <- read_entry_list(opt$entry_list)
    if (opt$checks > 0)
      fail("--checks generates check names; with --entry-list, flag checks in an is_check column instead.")
  } else {
    n <- opt$entries
    if (n < 1 || n != round(n)) fail("Give --entries N (a whole number) or --entry-list FILE.")
    test <- sprintf(paste0("G%0", nchar(n), "d"), seq_len(n))
    chk <- if (opt$checks > 0) paste0("CHK", seq_len(opt$checks)) else character()
    el <- data.frame(entry = c(chk, test), is_check = rep(c(TRUE, FALSE), c(length(chk), n)),
                     stringsAsFactors = FALSE)
  }
  el$entry_no <- seq_len(nrow(el))
  el
}

clean_name <- function(x) {
  y <- gsub("[^a-z0-9]+", "_", tolower(trimws(x)))
  gsub("^_|_$", "", y)
}

# 'irrigation:dry,wet' -> list(name = "irrigation", levels = c("dry", "wet"));
# 'nitrogen:3' -> levels nitrogen1..nitrogen3; 'dry,wet' -> name = default.
parse_factor <- function(spec, default_name) {
  spec <- trimws(spec)
  if (grepl(":", spec, fixed = TRUE)) {
    name <- clean_name(sub(":.*$", "", spec))
    rest <- sub("^[^:]*:", "", spec)
  } else {
    name <- default_name
    rest <- spec
  }
  if (!nzchar(name)) name <- default_name
  lev <- trimws(strsplit(rest, ",", fixed = TRUE)[[1]])
  lev <- lev[nzchar(lev)]
  if (length(lev) == 1L && grepl("^[0-9]+$", lev)) lev <- paste0(name, seq_len(as.integer(lev)))
  if (length(lev) < 2L) fail("Factor '", spec, "' needs at least 2 levels, e.g. 'irrigation:dry,wet'.")
  if (anyDuplicated(lev)) fail("Factor '", spec, "' has duplicated levels.")
  if (name %in% RESERVED) fail("Factor name '", name, "' clashes with a field-book column; rename it.")
  list(name = name, levels = lev)
}

# Run agricolae (which calls set.seed(seed, kinds) and so changes the session's
# RNG kind) and restore the previous RNG kind afterwards. Printed output is
# discarded; agricolae prints design parameters with cat().
with_agricolae <- function(expr) {
  old <- RNGkind()
  on.exit(RNGkind(old[1], old[2], old[3]), add = TRUE)
  res <- NULL
  utils::capture.output(res <- force(expr))
  res
}

# Run a FielDHub function, capturing what it prints (it prints its own error
# messages and returns NULL or an empty field book on bad input).
with_fieldhub <- function(expr) {
  res <- NULL
  printed <- utils::capture.output(
    res <- tryCatch(force(expr), error = function(e) structure(list(msg = conditionMessage(e)), class = "fh_error")))
  list(res = res, printed = printed)
}

divisors <- function(n) which(n %% seq_len(n) == 0)

# Divisor of n nearest sqrt(n), excluding 1 and n when possible.
near_sqrt_divisor <- function(n) {
  d <- divisors(n)
  d2 <- d[d > 1 & d < n]
  if (!length(d2)) return(NA_integer_)
  as.integer(d2[which.min(abs(d2 - sqrt(n)))])
}

# Efficiency factor of an equireplicate block design with equal block sizes:
# harmonic mean of the canonical efficiency factors (eigenvalues of C / r with
# C = r I - N N' / k). Returns 0 for a disconnected design.
efficiency_factor <- function(trt, blk) {
  N <- unclass(table(trt, blk))
  v <- nrow(N)
  r <- unique(rowSums(N))
  k <- unique(colSums(N))
  if (length(r) != 1L || length(k) != 1L)
    fail("efficiency_factor() needs equal replication and equal block sizes.")
  C <- r * diag(v) - (N %*% t(N)) / k
  ev <- sort(eigen(C / r, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)[seq_len(v - 1L)]
  if (any(ev < 1e-8)) return(0)
  (v - 1) / sum(1 / ev)
}

# Patterson & Williams' bound for resolvable designs, as printed by
# agricolae::design.alpha ("Efficiency factor").
alpha_upper_bound <- function(v, r, s) (v - 1) * (r - 1) / ((v - 1) * (r - 1) + r * (s - 1))

concurrence_range <- function(trt, blk) {
  N <- unclass(table(trt, blk))
  M <- N %*% t(N)
  off <- M[upper.tri(M)]
  c(min = min(off), max = max(off))
}

agricolae_alpha_ok <- function(r, k, s) {
  (r == 2 && k <= s) || (r == 3 && s %% 2 == 1 && k <= s) ||
    (r == 3 && s %% 2 == 0 && k < s) || (r == 4 && s %% 2 == 1 && s %% 3 != 0 && k <= s)
}

# ---------------------------------------------------------------------------
# Design generators. Each returns:
#   units  : one row per plot (entry, factor columns, rep, block, whole_plot,
#            idx_1..idx_L for free levels, fr/fc for a fixed inner grid, or
#            row/col when the package fixes the whole layout)
#   levels : nesting levels for the layout (outer to inner)
#   fixed  : TRUE when row/col come from the package (prep, diagonal)
#   info   : named list for the design summary
# ---------------------------------------------------------------------------

need <- function(ok, msg) if (!isTRUE(ok)) fail(msg)

gen_crd <- function(el, opt) {
  need(opt$reps >= 2, "CRD needs --reps of at least 2 (no error degrees of freedom otherwise).")
  d <- with_agricolae(agricolae::design.crd(el$entry, r = opt$reps, seed = opt$seed,
                                            kinds = AGRICOLAE_KINDS))
  bk <- d$book
  u <- data.frame(entry = as.character(bk[[3]]), rep = NA_integer_, block = NA_integer_,
                  idx_1 = seq_len(nrow(bk)), stringsAsFactors = FALSE)
  list(units = u, levels = list(list(name = "plot", n = nrow(u), weight = 0)),
       block_level = NA, info = list(package_function = "agricolae::design.crd"))
}

gen_rcbd <- function(el, opt) {
  need(opt$reps >= 2, "RCBD needs --reps of at least 2.")
  need(nrow(el) >= 2, "RCBD needs at least 2 entries.")
  d <- with_agricolae(agricolae::design.rcbd(el$entry, r = opt$reps, seed = opt$seed,
                                             kinds = AGRICOLAE_KINDS))
  bk <- d$book
  rep <- as.integer(as.character(bk$block))
  u <- data.frame(entry = as.character(bk[[3]]), rep = rep, block = NA_integer_,
                  idx_1 = rep, idx_2 = stats::ave(rep, rep, FUN = seq_along),
                  stringsAsFactors = FALSE)
  list(units = u,
       levels = list(list(name = "rep", n = opt$reps, weight = 3),
                     list(name = "plot", n = nrow(el), weight = 0)),
       block_level = 1, info = list(package_function = "agricolae::design.rcbd"))
}

gen_latin <- function(el, opt) {
  t <- nrow(el)
  need(t >= 3 && t <= 12, sprintf("A Latin square needs 3 to 12 treatments (got %d); use rcbd or rowcol otherwise.", t))
  nsq <- if (opt$reps >= 1) as.integer(opt$reps) else 1L
  us <- lapply(seq_len(nsq), function(q) {
    d <- with_agricolae(agricolae::design.lsd(el$entry, seed = opt$seed + q - 1,
                                              kinds = AGRICOLAE_KINDS))
    bk <- d$book
    data.frame(entry = as.character(bk[[4]]), rep = q, block = NA_integer_, idx_1 = q,
               fr = as.integer(as.character(bk$row)), fc = as.integer(as.character(bk$col)),
               stringsAsFactors = FALSE)
  })
  list(units = do.call(rbind, us),
       levels = list(list(name = "square", n = nsq, weight = 1),
                     list(name = "square_grid", n = t * t, weight = 0, gr = t, gc = t)),
       block_level = 1,
       info = list(package_function = "agricolae::design.lsd",
                   squares = nsq,
                   square_seeds = paste(opt$seed + seq_len(nsq) - 1, collapse = ", ")))
}

gen_alpha <- function(el, opt) {
  v <- nrow(el)
  r <- as.integer(opt$reps)
  need(r >= 2, "Alpha designs need --reps of at least 2 (resolvable replicates).")
  k <- as.integer(opt$block_size)
  auto_k <- FALSE
  if (k == 0) {
    k <- near_sqrt_divisor(v)
    auto_k <- TRUE
    if (is.na(k)) fail(sprintf("%d entries is prime: no block size divides it. Add or drop entries (e.g. checks) so the count factorises, or use rowcol/rcbd.", v))
  }
  need(k >= 2 && k < v, sprintf("--block-size must be between 2 and %d.", v - 1))
  if (v %% k != 0) {
    d <- divisors(v); d <- d[d > 1 & d < v]
    up <- ceiling(v / k) * k
    fail(sprintf(paste0("%d entries cannot be split into blocks of %d without dropping entries. ",
                        "Block sizes that divide %d: %s. Or add %d entries (e.g. checks) to reach %d = %d blocks of %d."),
                 v, k, v, if (length(d)) paste(d, collapse = ", ") else "none (prime)",
                 up - v, up, up / k, k))
  }
  s <- v / k
  engine <- tolower(opt$alpha_engine)
  need(engine %in% c("auto", "agricolae", "fieldhub"), "--alpha-engine must be auto, agricolae or fieldhub.")
  cands <- list()
  if (engine %in% c("auto", "agricolae")) {
    if (agricolae_alpha_ok(r, k, s)) {
      d <- with_agricolae(agricolae::design.alpha(el$entry, k = k, r = r, seed = opt$seed,
                                                  kinds = AGRICOLAE_KINDS))
      bk <- d$book
      blk <- as.integer(as.character(bk$block))
      rep <- as.integer(as.character(bk$replication))
      cands$agricolae <- data.frame(entry = as.character(bk[[4]]), rep = rep,
                                    block_in_rep = (blk - 1L) %% s + 1L,
                                    unit = as.integer(as.character(bk$cols)),
                                    stringsAsFactors = FALSE)
    } else if (engine == "agricolae") {
      fail(sprintf(paste0("agricolae::design.alpha has no alpha(0,1) generator for r = %d, k = %d, s = %d ",
                          "(it supports r = 2 with k <= s; r = 3 with s odd and k <= s, or s even and k < s; ",
                          "r = 4 with s odd, not a multiple of 3, and k <= s). Use --alpha-engine fieldhub or auto."), r, k, s))
    }
  }
  if (engine %in% c("auto", "fieldhub")) {
    fh <- with_fieldhub(FielDHub::alpha_lattice(t = v, k = k, r = r, seed = opt$seed,
                                                data = data.frame(ENTRY = seq_len(v), TREATMENT = el$entry)))
    if (inherits(fh$res, "fh_error") || is.null(fh$res$fieldBook)) {
      if (engine == "fieldhub" || !length(cands))
        fail("FielDHub::alpha_lattice could not build this design: ",
             if (inherits(fh$res, "fh_error")) fh$res$msg else paste(fh$printed, collapse = " "))
    } else {
      fb <- as.data.frame(fh$res$fieldBook)
      cands$fieldhub <- data.frame(entry = as.character(fb$TREATMENT), rep = as.integer(fb$REP),
                                   block_in_rep = as.integer(fb$IBLOCK), unit = as.integer(fb$UNIT),
                                   stringsAsFactors = FALSE)
    }
  }
  need(length(cands) > 0, "No alpha design could be generated with these parameters.")
  effs <- vapply(cands, function(cd) efficiency_factor(cd$entry, paste(cd$rep, cd$block_in_rep)), 0)
  pick <- names(effs)[which.max(effs + (names(effs) == "agricolae") * 1e-9)]
  cd <- cands[[pick]]
  cd <- cd[order(cd$rep, cd$block_in_rep, cd$unit), ]
  conc <- concurrence_range(cd$entry, paste(cd$rep, cd$block_in_rep))
  u <- data.frame(entry = cd$entry, rep = cd$rep, block = (cd$rep - 1L) * s + cd$block_in_rep,
                  idx_1 = cd$rep, idx_2 = cd$block_in_rep, idx_3 = cd$unit,
                  stringsAsFactors = FALSE)
  info <- list(package_function = if (pick == "agricolae") "agricolae::design.alpha" else "FielDHub::alpha_lattice",
               block_size = k, blocks_per_rep = s,
               block_size_chosen_automatically = auto_k,
               efficiency_factor = round(effs[[pick]], 4),
               efficiency_upper_bound = round(alpha_upper_bound(v, r, s), 4),
               concurrence_min = conc[["min"]], concurrence_max = conc[["max"]],
               concurrence_guarantee = if (pick == "agricolae") "alpha(0,1): every pair of entries shares at most one block" else
                 "none stated by FielDHub (A-optimised design); see concurrence_max")
  for (nm in names(effs)) info[[paste0("efficiency_factor_", nm)]] <- round(effs[[nm]], 4)
  if (pick == "fieldhub") {
    bm <- fh$res$blocksModel
    if (is.data.frame(bm) && "A-Efficiency" %in% names(bm))
      info$fieldhub_A_efficiency_blocks <- round(bm[["A-Efficiency"]][nrow(bm)], 4)
  }
  list(units = u,
       levels = list(list(name = "rep", n = r, weight = 2),
                     list(name = "block", n = s, weight = 3),
                     list(name = "plot", n = k, weight = 0)),
       block_level = 2, info = info)
}

gen_rowcol <- function(el, opt) {
  t <- nrow(el)
  r <- as.integer(opt$reps)
  need(r >= 2, "Row-column designs need --reps of at least 2.")
  nr <- as.integer(opt$rep_rows)
  auto <- FALSE
  if (nr == 0) {
    d <- divisors(t); d <- d[d > 1 & d <= sqrt(t)]
    if (!length(d)) fail(sprintf("%d entries cannot be arranged in a rows x columns replicate (prime). Add or drop entries.", t))
    nr <- as.integer(max(d))
    auto <- TRUE
  }
  need(nr >= 2 && nr < t && t %% nr == 0,
       sprintf("--rep-rows must divide the number of entries (%d) and lie between 2 and %d.", t, t - 1))
  # FielDHub 1.5.0 errors with "object 'best_efficiencies' not found" when no
  # swap improves its starting design (seed-dependent). Retry with seed + 1, ...
  # (deterministic, recorded in the summary).
  seed_used <- opt$seed
  for (attempt in 0:9) {
    seed_used <- opt$seed + attempt
    fh <- with_fieldhub(FielDHub::row_column(t = t, nrows = nr, r = r, seed = seed_used,
                                             iterations = as.integer(opt$iterations),
                                             data = data.frame(ENTRY = seq_len(t), TREATMENT = el$entry)))
    if (!(inherits(fh$res, "fh_error") && grepl("best_efficiencies", fh$res$msg, fixed = TRUE))) break
  }
  if (inherits(fh$res, "fh_error") || is.null(fh$res$fieldBook))
    fail("FielDHub::row_column could not build this design: ",
         if (inherits(fh$res, "fh_error")) fh$res$msg else paste(fh$printed, collapse = " "))
  fb <- as.data.frame(fh$res$fieldBook)
  nc <- t / nr
  u <- data.frame(entry = as.character(fb$TREATMENT), rep = as.integer(fb$REP), block = NA_integer_,
                  idx_1 = as.integer(fb$REP), fr = as.integer(fb$ROW), fc = as.integer(fb$COLUMN),
                  stringsAsFactors = FALSE)
  info <- list(package_function = "FielDHub::row_column", rows_per_rep = nr, cols_per_rep = nc,
               rows_per_rep_chosen_automatically = auto, iterations = as.integer(opt$iterations),
               fieldhub_seed_used = seed_used)
  bm <- fh$res$blocksModel
  if (is.list(bm) && !is.data.frame(bm)) bm <- bm[[1]]
  if (is.data.frame(bm) && "A-Efficiency" %in% names(bm)) {
    for (i in seq_len(nrow(bm)))
      info[[paste0("A_efficiency_", tolower(bm$Level[i]))]] <- round(bm[["A-Efficiency"]][i], 4)
  }
  list(units = u,
       levels = list(list(name = "rep", n = r, weight = 3),
                     list(name = "rep_grid", n = t, weight = 0, gr = nr, gc = nc)),
       block_level = 1, info = info)
}

gen_augmented <- function(el, opt) {
  b <- as.integer(if (opt$blocks > 0) opt$blocks else opt$reps)
  checks <- el$entry[el$is_check]
  new <- el$entry[!el$is_check]
  need(length(checks) >= 1, "Augmented designs need check entries: flag them in an is_check column or use --checks N.")
  need(b >= 2, "Augmented designs need --blocks of at least 2.")
  need(length(new) >= b, sprintf("Only %d new entries for %d blocks; use fewer blocks.", length(new), b))
  d <- with_agricolae(agricolae::design.dau(checks, new, r = b, seed = opt$seed,
                                            kinds = AGRICOLAE_KINDS))
  bk <- d$book
  blk <- as.integer(as.character(bk$block))
  u <- data.frame(entry = as.character(bk[[3]]), rep = NA_integer_, block = blk,
                  idx_1 = blk, idx_2 = stats::ave(blk, blk, FUN = seq_along),
                  stringsAsFactors = FALSE)
  sizes <- table(blk)
  m <- max(sizes)
  pad <- NULL
  if (any(sizes < m)) {
    pad <- do.call(rbind, lapply(names(sizes)[sizes < m], function(bl) {
      n0 <- sizes[[bl]]
      data.frame(entry = "FILLER", rep = NA_integer_, block = as.integer(bl),
                 idx_1 = as.integer(bl), idx_2 = seq(n0 + 1L, m), stringsAsFactors = FALSE)
    }))
    u <- rbind(u, pad)
  }
  c0 <- length(checks)
  info <- list(package_function = "agricolae::design.dau", blocks = b, checks = c0,
               new_entries = length(new),
               new_entries_per_block = paste(range(sizes - c0), collapse = "-"),
               filler_plots = if (is.null(pad)) 0L else nrow(pad),
               error_df_checks = (b - 1L) * (c0 - 1L))
  list(units = u,
       levels = list(list(name = "block", n = b, weight = 3),
                     list(name = "plot", n = m, weight = 0)),
       block_level = 1, info = info)
}

# Field dimensions for designs whose layout comes from FielDHub.
pick_dims <- function(total, opt) {
  R <- as.integer(opt$field_rows); C <- as.integer(opt$field_cols)
  if (R > 0 && C > 0) return(c(R, C))
  if (R > 0) { if (total %% R) fail(sprintf("%d plots do not fill %d rows.", total, R)); return(c(R, total / R)) }
  if (C > 0) { if (total %% C) fail(sprintf("%d plots do not fill %d columns.", total, C)); return(c(total / C, C)) }
  d <- divisors(total)
  d <- d[d <= sqrt(total)]
  if (max(d) == 1) fail(sprintf("%d plots is prime: no rectangular field fits. Change the number of entries, replicated entries or check copies.", total))
  c(max(d), total / max(d))
}

# Smallest |log(rows / cols)| over factorisations of n with 2 <= rows <= cols.
best_aspect <- function(n) {
  d <- divisors(n)
  d <- d[d >= 2 & d <= sqrt(n)]
  if (!length(d)) return(Inf)
  min(abs(log(d / (n / d))))
}

# Number of test entries to plant twice: fills the field when its size is
# given, otherwise stays near the requested fraction and picks the count whose
# plot total gives the most nearly square field.
prep_replicated_count <- function(n_test, base_total, opt) {
  R <- as.integer(opt$field_rows); C <- as.integer(opt$field_cols)
  target <- max(1L, as.integer(round(opt$prep_fraction * n_test)))
  if (R > 0 && C > 0) {
    n2 <- R * C - base_total
    if (n2 < 1 || n2 >= n_test)
      fail(sprintf(paste0("A %d x %d field holds %d plots; with %d unreplicated-plot slots for entries and checks, ",
                          "%d test entries would be planted twice (must be 1 to %d). Change the field size or check copies."),
                   R, C, R * C, base_total, n2, n_test - 1))
    return(n2)
  }
  span <- max(2L, as.integer(ceiling(0.1 * target)))
  cand <- unique(pmin(n_test - 1L, pmax(1L, target + seq(-span, span))))
  tot <- base_total + cand
  if (R > 0) cand <- cand[tot %% R == 0] else if (C > 0) cand <- cand[tot %% C == 0]
  if (!length(cand)) {
    dimv <- if (R > 0) R else C
    cand <- seq_len(n_test - 1L)[(base_total + seq_len(n_test - 1L)) %% dimv == 0]
    if (!length(cand)) fail(sprintf("No number of replicated entries fills a field with %d %s.", dimv, if (R > 0) "rows" else "columns"))
    return(cand[which.min(abs(cand - target))])
  }
  if (R > 0 || C > 0) return(cand[which.min(abs(cand - target))])
  asp <- vapply(base_total + cand, best_aspect, 0)
  if (all(!is.finite(asp))) fail("Could not find a plot total that forms a rectangular field; give --field-rows and --field-cols.")
  ord <- order(round(asp, 6), abs(cand - target))
  cand[ord[1]]
}

gen_prep <- function(el, opt) {
  test <- !el$is_check
  if ("reps" %in% names(el)) {
    reps <- el$reps
    src <- "reps column of the entry list"
  } else {
    need(opt$prep_fraction > 0 && opt$prep_fraction < 1, "--prep-fraction must be between 0 and 1.")
    need(sum(test) >= 2, "A p-rep design needs at least 2 test entries.")
    base_total <- sum(test) + sum(el$is_check) * as.integer(opt$check_reps)
    n2 <- prep_replicated_count(sum(test), base_total, opt)
    set.seed(opt$seed)
    idx <- which(test)
    pick <- idx[sample.int(length(idx), n2)]
    reps <- ifelse(el$is_check, as.integer(opt$check_reps), 1L)
    reps[pick] <- 2L
    src <- sprintf("%d of %d test entries (%.0f%%; requested %.0f%%) drawn at random with seed %d",
                   n2, sum(test), 100 * n2 / sum(test), 100 * opt$prep_fraction, opt$seed)
  }
  need(any(reps[test] >= 2), "A p-rep design needs some replicated test entries.")
  need(any(reps[test] == 1), "Every test entry is replicated: use rcbd or alpha instead of prep.")
  total <- sum(reps)
  dims <- pick_dims(total, opt)
  if (prod(dims) != total)
    fail(sprintf("The design needs %d plots but the field is %d x %d = %d. Feasible shapes: %s.",
                 total, dims[1], dims[2], prod(dims), feasible_shapes(total)))
  key <- sprintf("E%05d", seq_len(nrow(el)))
  fh <- with_fieldhub(FielDHub::partially_replicated(nrows = dims[1], ncols = dims[2], planter = "serpentine",
                                                     seed = opt$seed,
                                                     data = data.frame(ENTRY = seq_len(nrow(el)), NAME = key, REPS = reps)))
  fb <- if (!inherits(fh$res, "fh_error") && !is.null(fh$res$fieldBook)) as.data.frame(fh$res$fieldBook) else NULL
  if (is.null(fb) || nrow(fb) != total)
    fail("FielDHub::partially_replicated could not build this design. ",
         if (inherits(fh$res, "fh_error")) fh$res$msg else paste(trimws(fh$printed), collapse = " "))
  u <- data.frame(entry = el$entry[match(fb$TREATMENT, key)], rep = NA_integer_, block = NA_integer_,
                  row = as.integer(fb$ROW), col = as.integer(fb$COLUMN), stringsAsFactors = FALSE)
  need(!anyNA(u$entry), "Internal error: FielDHub returned unknown entry keys.")
  info <- list(package_function = "FielDHub::partially_replicated",
               replicated_test_entries = sum(reps[test] >= 2),
               unreplicated_test_entries = sum(reps[test] == 1),
               share_test_entries_replicated = round(mean(reps[test] >= 2), 3),
               checks = sum(el$is_check), check_copies = if (any(el$is_check)) paste(unique(reps[el$is_check]), collapse = ",") else "0",
               replication_source = src,
               fieldhub_min_distance_between_copies = if (!is.null(fh$res$infoDesign$min_distance)) round(fh$res$infoDesign$min_distance, 2) else NA)
  list(units = u, levels = NULL, fixed = TRUE, dims = dims, reps = setNames(reps, el$entry),
       block_level = NA, info = info)
}

gen_diagonal <- function(el, opt) {
  R <- as.integer(opt$field_rows); C <- as.integer(opt$field_cols)
  need(R > 0 && C > 0, "Diagonal designs need --field-rows and --field-cols (FielDHub fits the checks to the field).")
  chk <- el$entry[el$is_check]; lines <- el$entry[!el$is_check]
  need(length(chk) >= 1, "Diagonal designs need checks: flag them in is_check or use --checks N.")
  key <- c(sprintf("C%04d", seq_along(chk)), sprintf("L%05d", seq_along(lines)))
  fh <- with_fieldhub(FielDHub::diagonal_arrangement(nrows = R, ncols = C, lines = length(lines),
                                                     checks = length(chk), planter = "serpentine",
                                                     kindExpt = "SUDC", seed = opt$seed,
                                                     data = data.frame(ENTRY = seq_along(key), NAME = key)))
  fb <- if (!inherits(fh$res, "fh_error") && !is.null(fh$res$fieldBook)) as.data.frame(fh$res$fieldBook) else NULL
  if (is.null(fb) || !nrow(fb))
    fail(sprintf("FielDHub::diagonal_arrangement cannot fit %d lines and %d checks in a %d x %d field. ", length(lines), length(chk), R, C),
         paste(trimws(fh$printed[nzchar(trimws(fh$printed))]), collapse = " "),
         if (inherits(fh$res, "fh_error")) fh$res$msg else "")
  real <- c(chk, lines)
  ent <- real[match(fb$TREATMENT, key)]
  ent[is.na(ent) & grepl("fill", fb$TREATMENT, ignore.case = TRUE)] <- "FILLER"
  need(!anyNA(ent), "Internal error: FielDHub returned unknown entry keys.")
  u <- data.frame(entry = ent, rep = NA_integer_, block = NA_integer_,
                  row = as.integer(fb$ROW), col = as.integer(fb$COLUMN), stringsAsFactors = FALSE)
  info <- list(package_function = "FielDHub::diagonal_arrangement (SUDC)",
               lines = length(lines), checks = length(chk),
               check_plots = sum(ent %in% chk),
               percent_check_plots = round(100 * mean(ent %in% chk), 1),
               filler_plots = sum(ent == "FILLER"))
  list(units = u, levels = NULL, fixed = TRUE, dims = c(R, C), block_level = NA, info = info)
}

factor_or_entries <- function(spec, el, default_name) {
  if (nzchar(spec)) return(parse_factor(spec, default_name))
  list(name = "entry_factor", levels = el$entry, from_entries = TRUE)
}

gen_split <- function(el, opt) {
  need(nzchar(opt$factor_a), "Split-plot needs --factor-a (the whole-plot factor), e.g. 'irrigation:dry,wet'.")
  A <- parse_factor(opt$factor_a, "factor_a")
  B <- factor_or_entries(opt$factor_b, el, "factor_b")
  need(opt$reps >= 2, "Split-plot needs --reps of at least 2.")
  d <- with_agricolae(agricolae::design.split(A$levels, B$levels, r = opt$reps, design = "rcbd",
                                              seed = opt$seed, kinds = AGRICOLAE_KINDS))
  bk <- d$book
  rep <- as.integer(as.character(bk$block))
  wp <- as.integer(as.character(bk$plots))
  wp_id <- match(wp, unique(wp))
  wp_in_rep <- stats::ave(wp_id, rep, FUN = function(x) match(x, unique(x)))
  u <- data.frame(entry = paste(bk[[4]], bk[[5]], sep = "_"), rep = rep, block = NA_integer_,
                  whole_plot = wp_id, idx_1 = rep, idx_2 = wp_in_rep,
                  idx_3 = as.integer(as.character(bk$splots)), stringsAsFactors = FALSE)
  u[[A$name]] <- as.character(bk[[4]])
  u[[if (isTRUE(B$from_entries)) "entry_b" else B$name]] <- as.character(bk[[5]])
  list(units = u,
       levels = list(list(name = "rep", n = opt$reps, weight = 3),
                     list(name = "whole_plot", n = length(A$levels), weight = 0.5),
                     list(name = "sub_plot", n = length(B$levels), weight = 0)),
       block_level = 1, factors = list(A = A, B = B),
       info = list(package_function = "agricolae::design.split (RCBD whole plots)",
                   whole_plot_factor = A$name, sub_plot_factor = if (isTRUE(B$from_entries)) "entry" else B$name,
                   whole_plot_error_df = (opt$reps - 1) * (length(A$levels) - 1),
                   sub_plot_error_df = length(A$levels) * (opt$reps - 1) * (length(B$levels) - 1)))
}

gen_strip <- function(el, opt) {
  need(nzchar(opt$factor_a), "Strip-plot needs --factor-a (applied in row strips), e.g. 'tillage:conv,notill'.")
  A <- parse_factor(opt$factor_a, "factor_a")
  B <- factor_or_entries(opt$factor_b, el, "factor_b")
  need(opt$reps >= 2, "Strip-plot needs --reps of at least 2.")
  d <- with_agricolae(agricolae::design.strip(A$levels, B$levels, r = opt$reps, seed = opt$seed,
                                              kinds = AGRICOLAE_KINDS))
  bk <- d$book
  rep <- as.integer(as.character(bk$block))
  a <- as.character(bk[[3]]); b <- as.character(bk[[4]])
  # Row of each A level and column of each B level, in the order agricolae randomised them.
  fr <- fc <- integer(nrow(bk))
  for (q in unique(rep)) {
    i <- which(rep == q)
    a_ord <- unique(a[i]); b_ord <- unique(b[i])
    fr[i] <- match(a[i], a_ord); fc[i] <- match(b[i], b_ord)
  }
  u <- data.frame(entry = paste(a, b, sep = "_"), rep = rep, block = NA_integer_,
                  idx_1 = rep, fr = fr, fc = fc, stringsAsFactors = FALSE)
  u[[A$name]] <- a
  u[[if (isTRUE(B$from_entries)) "entry_b" else B$name]] <- b
  na <- length(A$levels); nb <- length(B$levels)
  list(units = u,
       levels = list(list(name = "rep", n = opt$reps, weight = 3),
                     list(name = "strip_grid", n = na * nb, weight = 0, gr = na, gc = nb)),
       block_level = 1, factors = list(A = A, B = B),
       info = list(package_function = "agricolae::design.strip",
                   row_strip_factor = A$name, column_strip_factor = if (isTRUE(B$from_entries)) "entry" else B$name,
                   row_strip_error_df = (opt$reps - 1) * (na - 1),
                   column_strip_error_df = (opt$reps - 1) * (nb - 1),
                   interaction_error_df = (opt$reps - 1) * (na - 1) * (nb - 1)))
}

gen_factorial <- function(el, opt) {
  specs <- c(opt$factor_a, opt$factor_b, opt$factor_c)
  specs <- specs[nzchar(specs)]
  need(length(specs) >= 2, "Factorial needs at least --factor-a and --factor-b.")
  need(opt$reps >= 2, "Factorial in RCBD needs --reps of at least 2.")
  fs <- Map(parse_factor, specs, c("factor_a", "factor_b", "factor_c")[seq_along(specs)])
  nm <- vapply(fs, `[[`, "", "name")
  need(!anyDuplicated(nm), "Factor names must differ.")
  d <- with_agricolae(agricolae::design.ab(vapply(fs, function(f) length(f$levels), 1L),
                                           r = opt$reps, design = "rcbd", seed = opt$seed,
                                           kinds = AGRICOLAE_KINDS))
  bk <- d$book
  rep <- as.integer(as.character(bk$block))
  labs <- lapply(seq_along(fs), function(j) fs[[j]]$levels[as.integer(as.character(bk[[2 + j]]))])
  u <- data.frame(entry = do.call(paste, c(labs, sep = "_")), rep = rep, block = NA_integer_,
                  idx_1 = rep, idx_2 = stats::ave(rep, rep, FUN = seq_along), stringsAsFactors = FALSE)
  for (j in seq_along(fs)) u[[fs[[j]]$name]] <- labs[[j]]
  ncomb <- prod(vapply(fs, function(f) length(f$levels), 1L))
  list(units = u,
       levels = list(list(name = "rep", n = opt$reps, weight = 3),
                     list(name = "plot", n = ncomb, weight = 0)),
       block_level = 1, factors = setNames(fs, LETTERS[seq_along(fs)]),
       info = list(package_function = "agricolae::design.ab (RCBD)", treatment_combinations = ncomb))
}

# ---------------------------------------------------------------------------
# Field layout: nested rectangles chosen by a small exhaustive search
# ---------------------------------------------------------------------------

divisor_pairs <- function(n) {
  d <- divisors(n)
  cbind(gr = d, gc = n %/% d)
}

# Every arrangement of the nested levels; one row per candidate.
layout_candidates <- function(levels) {
  opts <- lapply(levels, function(L) if (!is.null(L$gr)) cbind(gr = L$gr, gc = L$gc) else divisor_pairs(L$n))
  grid <- as.matrix(expand.grid(lapply(opts, function(o) seq_len(nrow(o)))))
  nl <- length(levels)
  gr <- matrix(0L, nrow(grid), nl); gc <- gr
  for (i in seq_len(nl)) {
    gr[, i] <- opts[[i]][grid[, i], "gr"]
    gc[, i] <- opts[[i]][grid[, i], "gc"]
  }
  list(gr = gr, gc = gc)
}

# Size (rows, cols) of one group at each level, plus the whole field.
group_sizes <- function(gr, gc) {
  nl <- length(gr)
  h <- vapply(seq_len(nl), function(i) prod(gr[seq_len(nl) > i]), 1)
  w <- vapply(seq_len(nl), function(i) prod(gc[seq_len(nl) > i]), 1)
  list(h = h, w = w, H = prod(gr), W = prod(gc))
}

layout_score <- function(gr, gc, levels, block_level) {
  g <- group_sizes(gr, gc)
  wts <- vapply(levels, function(L) L$weight, 1)
  s <- sum(wts * abs(log(g$h / g$w))) + abs(log(g$H / g$W))
  if (!is.na(block_level)) s <- s + 1e-3 * (g$h[block_level] > g$w[block_level])
  s + 1e-4 * (gc[1] > gr[1])
}

choose_layout <- function(levels, block_level, field_rows = 0, field_cols = 0, block_rows = 0) {
  cand <- layout_candidates(levels)
  n <- nrow(cand$gr)
  ok <- rep(TRUE, n); score <- rep(Inf, n); H <- W <- hb <- numeric(n)
  for (i in seq_len(n)) {
    g <- group_sizes(cand$gr[i, ], cand$gc[i, ])
    H[i] <- g$H; W[i] <- g$W
    hb[i] <- if (is.na(block_level)) NA else g$h[block_level]
    score[i] <- layout_score(cand$gr[i, ], cand$gc[i, ], levels, block_level)
  }
  if (block_rows > 0) {
    if (is.na(block_level)) fail("--block-rows does not apply to this design (it has no blocks).")
    ok <- ok & hb == block_rows
    if (!any(ok)) fail(sprintf("No layout gives blocks %d rows tall. Possible block heights: %s.",
                               block_rows, paste(sort(unique(hb)), collapse = ", ")))
  }
  feas <- ok
  if (field_rows > 0) ok <- ok & H == field_rows
  if (field_cols > 0) ok <- ok & W == field_cols
  if (!any(ok)) {
    shapes <- unique(data.frame(H = H[feas], W = W[feas]))
    shapes <- shapes[order(shapes$H), ]
    fail(sprintf("The design (%d plots) cannot be laid out as compact blocks in a field of %s x %s. Field shapes that work: %s.",
                 as.integer(H[1] * W[1]), if (field_rows > 0) field_rows else "any",
                 if (field_cols > 0) field_cols else "any",
                 paste(sprintf("%d x %d", shapes$H, shapes$W), collapse = ", ")))
  }
  i <- which(ok)[which.min(score[ok])]
  list(gr = cand$gr[i, ], gc = cand$gc[i, ])
}

apply_layout <- function(units, levels, choice) {
  nl <- length(levels)
  g <- group_sizes(choice$gr, choice$gc)
  row <- rep(1L, nrow(units)); col <- rep(1L, nrow(units))
  for (i in seq_len(nl)) {
    if (!is.null(levels[[i]]$gr)) {
      ri <- units$fr - 1L; ci <- units$fc - 1L
    } else {
      idx <- units[[paste0("idx_", i)]]
      ri <- (idx - 1L) %/% choice$gc[i]
      ci <- (idx - 1L) %% choice$gc[i]
    }
    row <- row + as.integer(ri * g$h[i])
    col <- col + as.integer(ci * g$w[i])
  }
  units$row <- row
  units$col <- col
  units
}

feasible_shapes <- function(total) {
  d <- divisors(total)
  paste(sprintf("%d x %d", d, total %/% d), collapse = ", ")
}

number_plots <- function(row, col, numbering = "serpentine", start = 1) {
  numbering <- tolower(numbering)
  if (!numbering %in% c("serpentine", "cartesian")) fail("--numbering must be serpentine or cartesian.")
  key_col <- if (numbering == "serpentine") ifelse(row %% 2 == 1, col, -col) else col
  ord <- order(row, key_col)
  plot <- integer(length(row))
  plot[ord] <- seq_along(row) + as.integer(start) - 1L
  plot
}

# ---------------------------------------------------------------------------
# Design-property checks (run on every design; a failure stops the run)
# ---------------------------------------------------------------------------

is_rectangle <- function(r, c) {
  length(r) == (diff(range(r)) + 1) * (diff(range(c)) + 1) && !anyDuplicated(paste(r, c))
}

check_design <- function(book, design, gen) {
  res <- list()
  add <- function(name, ok, detail = "") res[[length(res) + 1L]] <<- data.frame(check = name, passed = isTRUE(ok), detail = detail, stringsAsFactors = FALSE)
  pos <- paste(book$row, book$col)
  R <- max(book$row); C <- max(book$col)
  add("every plot has a unique field position", !anyDuplicated(pos))
  add("plots fill the field rectangle", nrow(book) == R * C && min(book$row) == 1 && min(book$col) == 1,
      sprintf("%d plots, %d x %d field", nrow(book), R, C))
  add("plot numbers are unique", !anyDuplicated(book$plot))
  real <- book[!book$is_filler, ]
  per_rep_once <- function(df) all(tapply(df$entry, df$rep, function(x) !anyDuplicated(x) && length(x) == length(unique(df$entry))))
  rect_by <- function(df, g) all(vapply(split(seq_len(nrow(df)), g), function(i) is_rectangle(df$row[i], df$col[i]), TRUE))
  if (design %in% c("rcbd", "factorial", "alpha", "rowcol", "split", "strip")) {
    add("every entry appears exactly once per replicate", per_rep_once(real))
    add("each replicate is a compact rectangle", rect_by(book, book$rep))
  }
  if (design == "latin")
    add("each square is a compact rectangle", rect_by(book, book$rep))
  if (design == "crd")
    add("every entry has the requested number of plots", length(unique(table(real$entry))) == 1L)
  if (design == "latin") {
    ok <- all(vapply(split(book, book$rep), function(s) all(tapply(s$entry, s$row, function(x) !anyDuplicated(x))) &&
                       all(tapply(s$entry, s$col, function(x) !anyDuplicated(x))) &&
                       length(unique(s$row)) == length(unique(s$entry)), TRUE))
    add("Latin square: each treatment once per row and once per column", ok)
  }
  if (design == "alpha") {
    k <- gen$info$block_size
    add("every incomplete block has the requested size", all(table(book$block) == k), sprintf("k = %d", k))
    add("incomplete blocks are nested in replicates", all(tapply(book$rep, book$block, function(x) length(unique(x)) == 1)))
    add("each incomplete block is a compact rectangle", rect_by(book, book$block))
    if (gen$info$package_function == "agricolae::design.alpha")
      add("alpha(0,1): no pair of entries shares more than one block", gen$info$concurrence_max <= 1,
          sprintf("max concurrence %d", gen$info$concurrence_max))
  }
  if (design == "rowcol") {
    ok <- all(vapply(split(book, book$rep), function(s) all(table(s$row) == length(unique(s$col))) &&
                       all(table(s$col) == length(unique(s$row))), TRUE))
    add("rows and columns of each replicate are complete", ok)
  }
  if (design == "augmented") {
    chk <- unique(book$entry[book$is_check])
    ok <- all(vapply(split(book$entry, book$block), function(x) all(chk %in% x) && sum(x %in% chk) == length(chk), TRUE))
    add("every block holds every check exactly once", ok, sprintf("%d checks", length(chk)))
    add("every new entry appears exactly once", all(table(real$entry[!real$is_check]) == 1))
    add("each block is a compact rectangle", rect_by(book, book$block))
  }
  if (design == "split") {
    A <- gen$factors$A$name
    wp_ok <- all(vapply(split(book, book$whole_plot), function(s) length(unique(s[[A]])) == 1 &&
                          !anyDuplicated(s$entry) && nrow(s) == length(gen$factors$B$levels), TRUE))
    add("each whole plot holds one whole-plot level and every sub-plot level once", wp_ok)
    add("each whole plot is contiguous (a rectangle)", rect_by(book, book$whole_plot))
    add("every whole-plot level appears once per replicate",
        all(tapply(book$whole_plot, book$rep, function(x) length(unique(x)) == length(gen$factors$A$levels))))
  }
  if (design == "strip") {
    A <- gen$factors$A$name
    B <- if (isTRUE(gen$factors$B$from_entries)) "entry_b" else gen$factors$B$name
    ok <- all(vapply(split(book, book$rep), function(s) all(tapply(s[[A]], s$row, function(x) length(unique(x)) == 1)) &&
                       all(tapply(s[[B]], s$col, function(x) length(unique(x)) == 1)), TRUE))
    add("strips intact: factor A constant along each row and factor B down each column of a replicate", ok)
  }
  if (design == "prep") {
    want <- gen$reps
    got <- table(book$entry)
    add("every entry has its requested number of copies", all(got[names(want)] == want))
  }
  if (design == "diagonal") {
    lines <- unique(book$entry[!book$is_check & !book$is_filler])
    add("every test line appears exactly once", all(table(book$entry[book$entry %in% lines]) == 1))
    add("checks are present", any(book$is_check))
  }
  do.call(rbind, res)
}

# Residual degrees of freedom of the fixed-effects skeleton ANOVA.
skeleton_error_df <- function(book, design) {
  b <- book[!book$is_filler, ]
  b$y <- stats::rnorm(nrow(b))
  f <- switch(design,
              crd = y ~ entry, prep = y ~ entry, diagonal = y ~ entry,
              rcbd = y ~ factor(rep) + entry, factorial = y ~ factor(rep) + entry,
              alpha = y ~ factor(rep) + factor(block) + entry,
              augmented = y ~ factor(block) + entry,
              latin = if (length(unique(b$rep)) > 1) y ~ factor(rep) / (factor(row) + factor(col)) + entry else y ~ factor(row) + factor(col) + entry,
              rowcol = y ~ factor(rep) / (factor(row) + factor(col)) + entry,
              NULL)
  if (is.null(f)) return(NA_integer_)
  as.integer(stats::lm(f, data = b)$df.residual)
}

# ---------------------------------------------------------------------------
# Build the design (no files written; tests call this directly)
# ---------------------------------------------------------------------------

build_design <- function(opt) {
  design <- normalise_design(opt$design)
  if (opt$seed == 0) {
    set.seed(NULL)
    opt$seed <- sample.int(999999L, 1L)
  }
  opt$seed <- as.integer(opt$seed)
  need(opt$seed >= 1, "--seed must be a positive whole number.")
  factor_only <- design == "factorial" || (design %in% c("split", "strip") && nzchar(opt$factor_b))
  if (factor_only && (nzchar(opt$entry_list) || opt$entries > 0))
    fail(if (design == "factorial") "Factorial designs take factors only (--factor-a/-b/-c), not --entries or --entry-list."
         else "Give either --factor-b or entries (--entries / --entry-list) for factor B, not both.")
  el <- if (factor_only) NULL else make_entries(opt)
  if (design %in% c("split", "strip") && !nzchar(opt$factor_b) && is.null(el))
    fail("Give --factor-b, or --entries/--entry-list to use the entries as the sub-plot factor.")
  gen <- switch(design, crd = gen_crd(el, opt), rcbd = gen_rcbd(el, opt), latin = gen_latin(el, opt),
                alpha = gen_alpha(el, opt), rowcol = gen_rowcol(el, opt), augmented = gen_augmented(el, opt),
                prep = gen_prep(el, opt), diagonal = gen_diagonal(el, opt), split = gen_split(el, opt),
                strip = gen_strip(el, opt), factorial = gen_factorial(el, opt))
  u <- gen$units
  if (isTRUE(gen$fixed)) {
    if (opt$block_rows > 0) fail("--block-rows does not apply to ", design, " designs (no blocks).")
    shape <- sprintf("%d rows x %d columns (layout from %s)", gen$dims[1], gen$dims[2], gen$info$package_function)
  } else {
    choice <- choose_layout(gen$levels, gen$block_level, opt$field_rows, opt$field_cols, opt$block_rows)
    u <- apply_layout(u, gen$levels, choice)
    g <- group_sizes(choice$gr, choice$gc)
    parts <- vapply(seq_along(gen$levels), function(i) {
      L <- gen$levels[[i]]
      if (length(gen$levels) == 1L)
        sprintf("%d plots arranged %d x %d", L$n, choice$gr[i], choice$gc[i])
      else if (i == length(gen$levels) && is.null(L$gr) && L$n > 1)
        sprintf("each %s: %d x %d plots", gen$levels[[max(1, i - 1)]]$name, choice$gr[i], choice$gc[i])
      else if (!is.null(L$gr)) sprintf("%s grid: %d x %d plots", gen$levels[[max(1, i - 1)]]$name, L$gr, L$gc)
      else sprintf("%d %s(s) arranged %d x %d", L$n, L$name, choice$gr[i], choice$gc[i])
    }, "")
    shape <- paste0(sprintf("%d rows x %d columns; ", as.integer(g$H), as.integer(g$W)), paste(parts, collapse = "; "))
  }
  # Assemble the field book
  if (is.null(el)) {
    entry_no <- match(u$entry, unique(u$entry[order(u$entry)]))
    is_check <- rep(FALSE, nrow(u))
    extra <- NULL
  } else if (design %in% c("split", "strip") && !nzchar(opt$factor_b)) {
    # entries are the sub-plot factor; combination label in 'entry'
    m <- match(u$entry_b, el$entry)
    entry_no <- as.integer(factor(u$entry, levels = unique(u$entry[order(m, u$entry)])))
    is_check <- el$is_check[m]
    extra <- NULL
  } else {
    m <- match(u$entry, el$entry)
    entry_no <- el$entry_no[m]
    is_check <- el$is_check[m]
    extra <- el[m, setdiff(names(el), c("entry", "is_check", "entry_no", "reps")), drop = FALSE]
  }
  is_filler <- u$entry == "FILLER"
  entry_no[is_filler] <- NA_integer_
  is_check[is_filler | is.na(is_check)] <- FALSE
  plot <- number_plots(u$row, u$col, opt$numbering, opt$plot_start)
  book <- data.frame(trial = opt$trial, plot_id = paste(opt$trial, plot, sep = "-"), plot = plot,
                     row = u$row, col = u$col, rep = u$rep, block = u$block, design = design,
                     stringsAsFactors = FALSE)
  if ("whole_plot" %in% names(u)) book$whole_plot <- u$whole_plot
  book$entry_no <- entry_no
  book$entry <- u$entry
  book$is_check <- is_check
  book$is_filler <- is_filler
  fcols <- setdiff(names(u), c("entry", "rep", "block", "whole_plot", "row", "col", "fr", "fc",
                               grep("^idx_", names(u), value = TRUE)))
  for (f in fcols) book[[f]] <- u[[f]]
  if (!is.null(extra) && ncol(extra)) book <- cbind(book, extra, stringsAsFactors = FALSE)
  if (grepl("[/\\\\]", opt$trial)) fail("--trial must not contain '/' or '\\' (Field Book rejects them in IDs).")
  book <- book[order(book$plot), ]
  rownames(book) <- NULL
  checks <- check_design(book, design, gen)
  if (!all(checks$passed))
    fail("Internal design check failed (please report): ",
         paste(checks$check[!checks$passed], collapse = "; "))
  set.seed(opt$seed)
  edf <- skeleton_error_df(book, design)
  list(book = book, design = design, gen = gen, checks = checks, shape = shape, seed = opt$seed,
       opt = opt, error_df = edf, entries = el)
}

# ---------------------------------------------------------------------------
# Warnings, summary, map and report
# ---------------------------------------------------------------------------

design_warnings <- function(res) {
  w <- character()
  i <- res$gen$info
  d <- res$design
  nent <- if (is.null(res$entries)) NA else nrow(res$entries)
  if (isTRUE(res$opt$seed_was_random))
    w <- c(w, sprintf("No --seed was given; seed %d was drawn at random. Keep it: it is the only way to regenerate this layout.", res$seed))
  if (d == "rcbd" && !is.na(nent) && nent > 30)
    w <- c(w, sprintf("RCBD with %d entries: blocks of %d plots are large and rarely uniform. An alpha (incomplete-block) or row-column design usually estimates entry differences more precisely.", nent, nent))
  if (d == "alpha") {
    if (isTRUE(i$block_size_chosen_automatically))
      w <- c(w, sprintf("--block-size not given: used k = %d (the divisor of %d nearest its square root). Confirm with the user.", i$block_size, nent))
    if (i$efficiency_factor < 0.9 * i$efficiency_upper_bound)
      w <- c(w, sprintf("Efficiency factor %.3f is well below the upper bound %.3f: blocks of %d are small for %d entries. A block size nearer sqrt(%d) = %.1f gives a more efficient design.",
                        i$efficiency_factor, i$efficiency_upper_bound, i$block_size, nent, nent, sqrt(nent)))
  }
  if (d == "rowcol" && isTRUE(i$rows_per_rep_chosen_automatically))
    w <- c(w, sprintf("--rep-rows not given: used %d rows x %d columns per replicate. Confirm with the user.", i$rows_per_rep, i$cols_per_rep))
  if (d == "rowcol" && i$fieldhub_seed_used != res$seed)
    w <- c(w, sprintf("FielDHub::row_column 1.5.0 fails for some seeds when its starting design cannot be improved ('best_efficiencies' not found); seed %d was used for FielDHub instead of %d. Re-running with --seed %d reproduces this layout.",
                      i$fieldhub_seed_used, res$seed, res$seed))
  if (d == "augmented") {
    if (i$error_df_checks < 10)
      w <- c(w, sprintf("Error df = (blocks - 1) x (checks - 1) = %d, below the rule-of-thumb minimum of about 10-12: tests of entry differences will have little power. Add checks or blocks.", i$error_df_checks))
    if (i$filler_plots > 0)
      w <- c(w, sprintf("New entries do not divide evenly over %d blocks: %d FILLER plot(s) pad the shorter blocks. Plant them with a border variety and do not analyse them.", i$blocks, i$filler_plots))
    if (i$checks < 2) w <- c(w, "Only one check: block effects are estimated from a single variety. Use at least 2 checks.")
  }
  if (d == "latin" && i$squares == 1 && !is.na(res$error_df) && res$error_df < 10)
    w <- c(w, sprintf("A single %dx%d Latin square has only %d error df. Consider --reps 2 or more squares.", nent, nent, res$error_df))
  if (d == "prep") {
    sh <- i$share_test_entries_replicated
    if (sh < 0.2 || sh > 0.3)
      w <- c(w, sprintf("%.0f%% of test entries are replicated; p-rep designs typically replicate about 20-30%% (25%% per FielDHub's reading of Cullis et al. 2006).", 100 * sh))
  }
  if (d == "diagonal" && i$filler_plots > 0)
    w <- c(w, sprintf("%d FILLER plot(s) were added by FielDHub to complete the field; plant them but do not analyse them.", i$filler_plots))
  if (d %in% c("split", "strip") && i[[grep("error_df", names(i))[1]]] < 4)
    w <- c(w, "Few error df for the whole-plot / row-strip factor: its effects are tested with little power. Add replicates if that factor matters.")
  if (!is.na(res$error_df) && res$error_df < 10 && !d %in% c("augmented", "latin", "split", "strip"))
    w <- c(w, sprintf("Only %d residual df in the skeleton ANOVA; consider more replicates.", res$error_df))
  w
}

design_summary <- function(res) {
  b <- res$book
  i <- res$gen$info
  base <- list(design = res$design,
               package_function = i$package_function,
               package_version = if (grepl("^FielDHub", i$package_function)) as.character(utils::packageVersion("FielDHub")) else as.character(utils::packageVersion("agricolae")),
               seed = res$seed,
               rng = if (grepl("^agricolae", i$package_function)) paste0(AGRICOLAE_KINDS, " (set by agricolae)") else "Mersenne-Twister (R default)",
               trial = res$opt$trial,
               entries = if (is.null(res$entries)) length(unique(b$entry)) else nrow(res$entries),
               checks = sum(!duplicated(b$entry[b$is_check])),
               plots = nrow(b),
               filler_plots = sum(b$is_filler),
               field_rows = max(b$row), field_cols = max(b$col),
               layout = res$shape,
               numbering = paste0(res$opt$numbering, ", plot ", min(b$plot), " at row 1 column 1"),
               residual_df_skeleton_anova = res$error_df)
  if (!all(is.na(b$rep))) base$replicates <- length(unique(stats::na.omit(b$rep)))
  extra <- i[setdiff(names(i), "package_function")]
  all <- c(base, extra)
  data.frame(parameter = names(all), value = vapply(all, function(x) paste(format(x), collapse = ", "), ""),
             stringsAsFactors = FALSE)
}

PALETTE <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
INK <- "#0b0b0b"; INK2 <- "#52514e"; SURFACE <- "#fcfcfb"

# Label colour with the higher WCAG contrast against the fill: white or ink.
text_on <- function(hex) {
  lum <- function(h) {
    rgb <- grDevices::col2rgb(h) / 255
    lin <- ifelse(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
    0.2126 * lin[1, ] + 0.7152 * lin[2, ] + 0.0722 * lin[3, ]
  }
  l <- lum(hex)
  ifelse((1.05) / (l + 0.05) > (l + 0.05) / (lum(INK) + 0.05), "#ffffff", INK)
}

outline_groups <- function(b, g, lwd, col = INK) {
  for (idx in split(seq_len(nrow(b)), g)) {
    graphics::rect(min(b$col[idx]) - 0.5, max(b$row[idx]) + 0.5, max(b$col[idx]) + 0.5,
                   min(b$row[idx]) - 0.5, border = col, lwd = lwd)
  }
}

draw_field_map <- function(res, out_dir) {
  b <- res$book
  d <- res$design
  R <- max(b$row); C <- max(b$col)
  w_in <- min(24, max(6.5, 0.5 * C + 1.6)); h_in <- min(24, max(4.5, 0.42 * R + 2.4))
  open_png(out_dir, "field_map.png", width = w_in, height = h_in)
  on.exit(grDevices::dev.off())
  graphics::par(mar = c(3.6, 3.2, 3.2, 0.8), bg = SURFACE, fg = INK2, col.axis = INK2)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0.5, C + 0.5), ylim = c(R + 0.5, 0.5), xaxs = "i", yaxs = "i")
  # Fill: factor levels for factorial-type designs, otherwise pale tints by replicate/block.
  fac <- NULL
  if (!is.null(res$gen$factors)) fac <- res$gen$factors$A$name
  if (d == "latin" && length(unique(b$entry)) <= length(PALETTE)) fac <- "entry"
  legend_lab <- character(); legend_col <- character()
  if (!is.null(fac) && length(unique(b[[fac]])) <= length(PALETTE)) {
    lv <- if (fac == "entry") sort(unique(b$entry)) else res$gen$factors$A$levels
    fill <- PALETTE[match(b[[fac]], lv)]
    legend_lab <- paste0(if (fac == "entry") "" else paste0(fac, " = "), lv); legend_col <- PALETTE[seq_along(lv)]
  } else {
    grp <- if (!all(is.na(b$block))) b$block else if (!all(is.na(b$rep))) b$rep else rep(1L, nrow(b))
    fill <- ifelse(grp %% 2 == 1, "#e3edfb", "#f4efe6")
  }
  fill[b$is_check] <- "#eb6834"
  fill[b$is_filler] <- "#d6d5d0"
  if (any(b$is_check)) { legend_lab <- c(legend_lab, "check"); legend_col <- c(legend_col, "#eb6834") }
  if (any(b$is_filler)) { legend_lab <- c(legend_lab, "filler"); legend_col <- c(legend_col, "#d6d5d0") }
  graphics::rect(b$col - 0.5, b$row + 0.5, b$col + 0.5, b$row - 0.5, col = fill, border = SURFACE, lwd = 1)
  # Outlines: whole plots / incomplete blocks thin, replicates / blocks thick.
  if ("whole_plot" %in% names(b)) outline_groups(b, b$whole_plot, 1.2, INK2)
  if (!all(is.na(b$block))) outline_groups(b, b$block, if (all(is.na(b$rep))) 2.6 else 1.4, INK)
  if (!all(is.na(b$rep))) outline_groups(b, b$rep, 3, INK)
  # Labels
  lab <- if (!is.null(res$gen$factors) && !d %in% c("factorial")) {
    bcol <- if (isTRUE(res$gen$factors$B$from_entries)) "entry_b" else res$gen$factors$B$name
    b[[bcol]]
  } else b$entry
  if (max(nchar(lab), na.rm = TRUE) > 7) lab <- ifelse(b$is_filler, "F", as.character(b$entry_no))
  lab[b$is_filler] <- "F"
  cell_w <- (w_in - 1.0) / C
  cell_h <- (h_in - 1.7) / R
  nch <- max(3, max(nchar(lab), na.rm = TRUE))
  cex <- max(0.3, min(0.9, 0.8 * cell_w / (nch * 0.075), cell_h / 0.3))
  if (C <= 60 && R <= 60)
    graphics::text(b$col, b$row + 0.08, lab, cex = cex, col = text_on(fill))
  if (C <= 30 && R <= 40)
    graphics::text(b$col - 0.45, b$row - 0.32, b$plot, cex = cex * 0.62, adj = c(0, 0.5), col = text_on(fill))
  at_c <- if (C <= 30) seq_len(C) else pretty(seq_len(C))
  at_r <- if (R <= 30) seq_len(R) else pretty(seq_len(R))
  graphics::axis(3, at = at_c[at_c >= 1 & at_c <= C], tick = FALSE, line = -0.6, cex.axis = 0.7)
  graphics::axis(2, at = at_r[at_r >= 1 & at_r <= R], tick = FALSE, line = -0.6, las = 1, cex.axis = 0.7)
  graphics::mtext("column", side = 3, line = 1.4, cex = 0.75, col = INK2)
  graphics::mtext("row", side = 2, line = 1.9, cex = 0.75, col = INK2)
  graphics::title(main = NULL)
  graphics::mtext(sprintf("%s design, %d plots (%d x %d), seed %d", res$design, nrow(b), R, C, res$seed),
                  side = 1, line = 0.6, adj = 0, cex = 0.8, col = INK, font = 2)
  outl <- c(if (!all(is.na(b$rep))) "thick outline = replicate" else NULL,
            if (!all(is.na(b$block))) (if (all(is.na(b$rep))) "thick outline = block" else "thin outline = incomplete block") else NULL,
            if ("whole_plot" %in% names(b)) "grey outline = whole plot" else NULL,
            sprintf("plot numbers %s from row 1, column 1 (top left)", res$opt$numbering))
  graphics::mtext(paste(outl, collapse = "; "), side = 1, line = 1.6, adj = 0, cex = 0.62, col = INK2)
  if (length(legend_lab)) {
    graphics::legend("bottomright", inset = c(0, -0.11 * 4.5 / h_in), legend = legend_lab, fill = legend_col,
                     border = NA, bty = "n", horiz = TRUE, cex = 0.62, xpd = TRUE, text.col = INK2)
  }
  invisible(TRUE)
}

md_table <- function(df) {
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

REFERENCES <- c(
  "- de Mendiburu F (2023). agricolae: Statistical Procedures for Agricultural Research. R package version 1.3-7. https://doi.org/10.32614/CRAN.package.agricolae",
  "- Murillo D, Gezan S, Heilman A, Walk T, Aparicio J, Horsley R (2021). FielDHub: A Shiny App for Design of Experiments in Life Sciences. Journal of Open Source Software 6(61):3122. https://doi.org/10.21105/joss.03122",
  "- Patterson HD, Williams ER (1976). A new class of resolvable incomplete block designs. Biometrika 63(1):83-92. https://doi.org/10.1093/biomet/63.1.83",
  "- Federer WT (1956). Augmented (or hoonuiaku) designs. Hawaiian Planters' Record 55(2):191-208. Also Cornell Biometrics Unit report BU-74-M, https://hdl.handle.net/1813/32841",
  "- Federer WT (1961). Augmented designs with one-way elimination of heterogeneity. Biometrics 17(3):447-473. https://doi.org/10.2307/2527837",
  "- Cullis BR, Smith AB, Coombes NE (2006). On the design of early generation variety trials with correlated data. Journal of Agricultural, Biological, and Environmental Statistics 11(4):381-393. https://doi.org/10.1198/108571106X154443",
  "- Clarke GPY, Stefanova KT (2011). Optimal design for early-generation plant-breeding trials with unreplicated or partially replicated test lines. Australian & New Zealand Journal of Statistics 53(4):461-480. https://doi.org/10.1111/j.1467-842X.2011.00642.x",
  "- Williams ER, John JA, Whitaker D (2006). Construction of resolvable spatial row-column designs. Biometrics 62(1):103-108. https://doi.org/10.1111/j.1541-0420.2005.00393.x")

write_report <- function(res, out_dir, warn, summ) {
  b <- res$book
  i <- res$gen$info
  key <- c("design", "package_function", "seed", "entries", "plots", "layout")
  eff <- NULL
  if (res$design == "alpha")
    eff <- c(sprintf("Efficiency factor E = %.3f (computed from the incidence matrix); upper bound for resolvable designs %.3f (Patterson & Williams 1976, as printed by agricolae).",
                     i$efficiency_factor, i$efficiency_upper_bound),
             sprintf("Pairwise concurrence: %d to %d. %s.", i$concurrence_min, i$concurrence_max, i$concurrence_guarantee),
             if (!is.null(i$efficiency_factor_agricolae) && !is.null(i$efficiency_factor_fieldhub))
               sprintf("Candidates compared: agricolae E = %.3f, FielDHub E = %.3f; the higher was kept.", i$efficiency_factor_agricolae, i$efficiency_factor_fieldhub) else NULL)
  if (res$design == "rowcol")
    eff <- sprintf("FielDHub A-efficiencies: %s.", paste(sprintf("%s %s", sub("A_efficiency_", "", grep("^A_efficiency_", names(i), value = TRUE)),
                                                                 unlist(i[grep("^A_efficiency_", names(i))])), collapse = ", "))
  if (res$design == "augmented")
    eff <- sprintf("Error df from checks = (blocks - 1) x (checks - 1) = (%d - 1) x (%d - 1) = %d.", i$blocks, i$checks, i$error_df_checks)
  lines <- c(sprintf("# %s report", SKILL), "",
             if (length(warn)) c("## Warnings", "", paste("-", warn), "") else NULL,
             "## Design", "",
             md_table(summ[summ$parameter %in% c(key, "replicates", "field_rows", "field_cols", "residual_df_skeleton_anova"), ]), "",
             if (length(eff)) c("## Efficiency", "", paste("-", eff), "") else NULL,
             "## Design checks", "",
             "Every check below was run on the final field book; the run stops if any fails.", "",
             md_table(transform(res$checks, passed = ifelse(passed, "yes", "NO"))), "",
             "## Field map", "", "![Field map](figures/field_map.png)", "",
             "Row 1 is drawn at the top and column 1 at the left. Orient the map to the real field before planting.", "",
             "## Files", "",
             "- `tables/field_book.csv`: one row per plot (plot_id, plot, row, col, rep, block, entry, is_check, is_filler, factor levels).",
             "- `tables/design_summary.csv`: parameters, package function, seed, efficiency.",
             "- `tables/design_checks.csv`: the design-property checks above.",
             "- Next step: `claw-fieldbook` turns the field book into a Field Book app import file.", "",
             "## References", "", REFERENCES, "",
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
    opt$entry_list <- file.path(SKILL_DIR, "demo", "entries.csv")
    if (!nzchar(opt$design)) opt$design <- "alpha"
    if (opt$reps == 0) opt$reps <- 3
    if (opt$block_size == 0) opt$block_size <- 5
    if (opt$seed == 0) opt$seed <- 2026
    if (opt$trial == DEFAULTS$trial) opt$trial <- "DEMO"
  }
  if (!nzchar(opt$design)) fail("Give --design (", paste(DESIGNS, collapse = ", "), ") or --demo. See --help.")
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  seed_given <- opt$seed > 0
  res <- build_design(opt)
  res$opt$seed_was_random <- !seed_given
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  warn <- design_warnings(res)
  summ <- design_summary(res)
  save_table(res$book, opt$out, "field_book.csv")
  save_table(summ, opt$out, "design_summary.csv")
  save_table(res$checks, opt$out, "design_checks.csv")
  draw_field_map(res, opt$out)
  write_report(res, opt$out, warn, summ)
  argv_rec <- if (seed_given) argv else c(argv, "--seed", as.character(res$seed))
  write_repro_bundle(opt$out, argv_rec, inputs = opt$entry_list)
  cat(sprintf("Done: %s design, %d plots (%d x %d), seed %d. Report: %s\n", res$design, nrow(res$book),
              max(res$book$row), max(res$book$col), res$seed,
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  if (length(warn)) cat(paste("Warning:", warn), sep = "\n")
  invisible(res)
}

if (sys.nframe() == 0L) main()
