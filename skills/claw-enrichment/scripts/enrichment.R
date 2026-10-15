#!/usr/bin/env Rscript
# claw-enrichment: over-representation analysis (hypergeometric test, BH) of GO,
# KEGG, MapMan, PANTHER or any custom gene sets with clusterProfiler::enricher,
# using annotation tables from crop genome projects (no OrgDb needed). The
# background (universe) is the set of genes that were actually tested, GO
# annotations can be propagated to ancestor terms over is_a/part_of (GO.db), and
# gene IDs are checked against the annotation before testing.
#
# Usage:
#   Rscript scripts/enrichment.R --demo
#   Rscript scripts/enrichment.R --de-table de_condition_drought_vs_control.csv \
#       --annotation Gmax_annotation_info.txt --gene-col locusName --term-col GO --propagate
#   Rscript scripts/enrichment.R --genes de_genes.txt --universe tested_genes.txt \
#       --annotation kegg_gene2pathway.csv --term-names kegg_names.csv
#   Rscript scripts/enrichment.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT license.

SKILL   <- "claw-enrichment"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/enrichment.R [--demo] (--de-table DE.csv | --genes FILE --universe FILE)
         --annotation FILE [options]

  --demo            run on the bundled synthetic DE table and GO annotation (demo/)
  --de-table        full DE results (e.g. claw-de-analysis tables/de_*.csv) with gene_id and
                    padj (and log2FoldChange for up/down). DE genes = padj < alpha; universe =
                    genes with non-missing padj (the genes actually tested)
  --genes           DE gene list (one ID per line, or a CSV with a gene_id column)
  --universe        all genes that were tested (same formats); required with --genes
  --annotation      gene-to-term table: two columns (gene, term) in long format, or one row per
                    gene with several terms separated by , ; | or spaces; or a GO GAF file
  --format          auto, table or gaf (default: auto)
  --gene-col        annotation column with gene IDs (name or number; default: first column,
                    or column 2 of a GAF)
  --term-col        annotation column with term IDs (default: second column)
  --term-names      optional CSV/TSV mapping term_id to a readable name (non-GO terms)
  --ontology        GO namespaces to test separately: any of BP,MF,CC (default: BP,MF,CC)
  --propagate       add ancestor GO terms over is_a and part_of (GO.db), for annotation files
                    that list only the most specific terms
  --direction       gene sets to test from --de-table: any of all,up,down (default: all,up,down)
  --alpha           FDR for calling DE genes from --de-table and for significant terms (default: 0.05)
  --lfc-min         also require |log2FoldChange| >= this for DE genes (default: 0)
  --min-size        smallest gene set tested, counted within the universe (default: 10)
  --max-size        largest gene set tested (default: 500)
  --min-overlap     stop if fewer than this fraction of universe or DE genes are in the
                    annotation (ID mismatch check; default: 0.2)
  --strip-isoform   annotation is per transcript/protein: drop .1, .1.p, _T001, _P001 suffixes
                    and merge terms per gene
  --out             output folder (default: agriclaw_output/claw-enrichment)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, de_table = "", genes = "", universe = "", annotation = "",
                 format = "auto", gene_col = "", term_col = "", term_names = "",
                 ontology = "BP,MF,CC", propagate = FALSE, direction = "all,up,down",
                 alpha = 0.05, lfc_min = 0, min_size = 10, max_size = 500,
                 min_overlap = 0.2, strip_isoform = FALSE, out = "")

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

fail <- function(...) stop(paste0(...), call. = FALSE)
split_list <- function(x, sep = ",") {
  if (!nzchar(x)) return(character())
  v <- trimws(strsplit(x, sep, fixed = TRUE)[[1L]])
  v[nzchar(v)]
}

read_delim_auto <- function(path, header = TRUE) {
  first <- readLines(path, n = 1L, warn = FALSE)
  sep <- if (grepl("\t", first, fixed = TRUE)) "\t" else ","
  utils::read.table(path, sep = sep, header = header, quote = "\"", comment.char = "",
                    colClasses = "character", check.names = FALSE, fill = TRUE)
}

# One ID per line, or a CSV/TSV with a gene_id column.
read_gene_list <- function(path, what) {
  if (!file.exists(path)) fail(what, " file not found: ", path)
  first <- readLines(path, n = 1L, warn = FALSE)
  if (grepl("[,\t]", first)) {
    d <- read_delim_auto(path)
    col <- if ("gene_id" %in% names(d)) "gene_id" else names(d)[[1L]]
    ids <- d[[col]]
  } else {
    ids <- readLines(path, warn = FALSE)
    if (length(ids) && tolower(ids[[1L]]) %in% c("gene_id", "gene", "id")) ids <- ids[-1L]
  }
  ids <- unique(trimws(ids))
  ids[nzchar(ids)]
}

# Gene sets from a full DE table: list(sets = named list, universe).
sets_from_de_table <- function(path, alpha, lfc_min, directions) {
  if (!file.exists(path)) fail("DE table not found: ", path)
  d <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  need <- c("gene_id", "padj")
  miss <- setdiff(need, names(d))
  if (length(miss)) fail("The DE table needs columns gene_id and padj (missing: ", paste(miss, collapse = ", "),
                         "). Use the full table from claw-de-analysis (tables/de_*.csv).")
  tested <- !is.na(d$padj)
  if (!any(tested)) fail("No gene in the DE table has a padj value, so there is no tested universe.")
  if (all(d$padj[tested] < alpha))
    fail("Every gene in the DE table is significant: this looks like a list of DE genes only. ",
         "Give the full results table (all tested genes), or --genes plus --universe.")
  has_lfc <- "log2FoldChange" %in% names(d)
  sig <- tested & d$padj < alpha
  if (lfc_min > 0) {
    if (!has_lfc) fail("--lfc-min needs a log2FoldChange column in the DE table.")
    sig <- sig & !is.na(d$log2FoldChange) & abs(d$log2FoldChange) >= lfc_min
  }
  sets <- list()
  for (dirn in directions) {
    if (dirn == "all") sets$all <- d$gene_id[sig]
    else {
      if (!has_lfc) next
      sets[[dirn]] <- d$gene_id[sig & if (dirn == "up") d$log2FoldChange > 0 else d$log2FoldChange < 0]
    }
  }
  list(sets = lapply(sets, unique), universe = unique(d$gene_id[tested]), n_rows = nrow(d),
       has_lfc = has_lfc)
}

# Every DE gene must be in the universe (the universe is the set of tested genes).
check_sets_in_universe <- function(sets, universe) {
  for (nm in names(sets)) {
    out_u <- setdiff(sets[[nm]], universe)
    if (length(out_u))
      fail(length(out_u), " gene(s) in the '", nm, "' list are not in the universe (e.g. ",
           paste(utils::head(out_u, 3), collapse = ", "), "). The universe must contain every tested gene, ",
           "including all DE genes.")
  }
  invisible(TRUE)
}

# --- Annotation -------------------------------------------------------------

col_pick <- function(d, spec, default_idx, what) {
  if (!nzchar(spec)) return(names(d)[[default_idx]])
  if (grepl("^[0-9]+$", spec)) {
    k <- as.integer(spec)
    if (k < 1 || k > ncol(d)) fail(what, " column number ", k, " is out of range (1..", ncol(d), ").")
    return(names(d)[[k]])
  }
  if (!spec %in% names(d)) fail(what, " column '", spec, "' not found. Columns: ", paste(names(d), collapse = ", "))
  spec
}

read_gaf <- function(path, gene_col = "") {
  con <- file(path, open = "r")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  lines <- lines[nzchar(lines) & !startsWith(lines, "!")]
  f <- strsplit(lines, "\t", fixed = TRUE)
  ok <- lengths(f) >= 15L
  if (!any(ok)) fail("No GAF annotation lines (17 tab-separated columns) in ", path)
  f <- f[ok]
  gcol <- if (nzchar(gene_col)) as.integer(gene_col) else 2L
  if (is.na(gcol) || gcol < 1 || gcol > 17) fail("For a GAF file, --gene-col is a column number (2 = DB Object ID, 3 = symbol).")
  qual <- vapply(f, `[`, "", 4L)
  d <- data.frame(gene = vapply(f, `[`, "", gcol), term = vapply(f, `[`, "", 5L), stringsAsFactors = FALSE)
  not <- grepl("NOT", qual, fixed = TRUE)
  list(t2g = unique(d[!not, c("term", "gene")]), n_not = sum(not))
}

# Returns list(t2g = data.frame(term, gene), t2n = data.frame(term, name) or NULL, notes)
read_annotation <- function(path, format = "auto", gene_col = "", term_col = "") {
  if (!file.exists(path)) fail("Annotation file not found: ", path)
  first <- readLines(path, n = 1L, warn = FALSE)
  if (format == "auto") format <- if (grepl("^!gaf-version", first) || grepl("\\.gaf(\\.gz)?$", path)) "gaf" else "table"
  if (format == "gaf") {
    g <- read_gaf(path, gene_col)
    return(list(t2g = g$t2g, t2n = NULL, notes = sprintf("GAF: %d NOT-qualified annotations excluded.", g$n_not)))
  }
  if (format != "table") fail("--format must be auto, table or gaf.")
  d <- read_delim_auto(path)
  if (ncol(d) < 2L) fail("The annotation needs at least two columns: gene ID and term ID(s).")
  gc <- col_pick(d, gene_col, 1L, "--gene-col")
  tc <- col_pick(d, term_col, if (gc == names(d)[[2L]]) 1L else 2L, "--term-col")
  terms <- split_terms(d[[tc]])
  n <- lengths(terms)
  t2g <- data.frame(term = trimws(unlist(terms, use.names = FALSE)), gene = rep(trimws(d[[gc]]), n),
                    stringsAsFactors = FALSE)
  t2g <- unique(t2g[nzchar(t2g$term) & nzchar(t2g$gene) & t2g$term != "NA", , drop = FALSE])
  if (!nrow(t2g)) fail("No gene-term pairs found in columns '", gc, "' and '", tc, "'.")
  t2n <- NULL
  name_col <- intersect(c("term_name", "name", "description", "Description"), names(d))
  if (length(name_col) && all(n <= 1)) {
    t2n <- unique(data.frame(term = d[[tc]], name = d[[name_col[1]]], stringsAsFactors = FALSE))
    t2n <- t2n[nzchar(t2n$term), , drop = FALSE]
  }
  list(t2g = t2g, t2n = t2n, notes = sprintf("Annotation table: gene column '%s', term column '%s'.", gc, tc))
}

# Split annotation cells on , ; | and also on whitespace between GO IDs. Term IDs that
# contain spaces (e.g. pathway names used as IDs) are kept whole.
split_terms <- function(x) {
  parts <- strsplit(ifelse(is.na(x), "", x), "[,;|]")
  lapply(parts, function(p) {
    p <- trimws(p)
    as.character(unlist(lapply(p, function(q) {
      tok <- strsplit(q, "[[:space:]]+")[[1L]]
      if (length(tok) > 1L && all(grepl("^GO:[0-9]{7}$", tok))) tok else q
    }), use.names = FALSE))
  })
}

read_term_names <- function(path) {
  if (!nzchar(path)) return(NULL)
  if (!file.exists(path)) fail("--term-names file not found: ", path)
  d <- read_delim_auto(path)
  if (ncol(d) < 2L) fail("--term-names needs two columns: term_id and name.")
  data.frame(term = d[[1L]], name = d[[2L]], stringsAsFactors = FALSE)
}

is_go <- function(terms) mean(grepl("^GO:[0-9]{7}$", terms)) > 0.9

# GO term table from GO.db, with secondary (alt) IDs mapped to primary IDs.
go_term_info <- function(terms) {
  terms <- unique(terms)
  info <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = terms,
                                                 columns = c("TERM", "ONTOLOGY"), keytype = "GOID"))
  info <- info[!is.na(info$ONTOLOGY), , drop = FALSE]
  missing <- setdiff(terms, info$GOID)
  remap <- stats::setNames(character(), character())
  if (length(missing)) {
    syn <- AnnotationDbi::toTable(GO.db::GOSYNONYM)
    syn <- unique(syn[syn$synonym %in% missing & grepl("^GO:", syn$synonym), c("synonym", "go_id")])
    syn <- syn[!duplicated(syn$synonym), , drop = FALSE]
    remap <- stats::setNames(syn$go_id, syn$synonym)
  }
  list(info = info, remap = remap, obsolete = setdiff(missing, names(remap)))
}

go_parent_edges <- function() {
  tabs <- list(BP = GO.db::GOBPPARENTS, MF = GO.db::GOMFPARENTS, CC = GO.db::GOCCPARENTS)
  do.call(rbind, lapply(names(tabs), function(o) {
    t <- AnnotationDbi::toTable(tabs[[o]])
    data.frame(child = t[[1L]], parent = t[[2L]], relation = t[[3L]], stringsAsFactors = FALSE)
  }))
}

PROPAGATE_RELATIONS <- c("isa", "is_a", "part of", "part_of")

# Add every ancestor over is_a/part_of edges (the GO true-path rule). Edges whose
# relation is regulates/positively_regulates/negatively_regulates are not followed.
propagate_annotation <- function(t2g, edges, relations = PROPAGATE_RELATIONS) {
  e <- edges[edges$relation %in% relations, , drop = FALSE]
  par <- split(e$parent, e$child)
  memo <- new.env(hash = TRUE, parent = emptyenv())
  anc <- function(t) {
    r <- get0(t, envir = memo, inherits = FALSE)
    if (!is.null(r)) return(r)
    p <- par[[t]]
    r <- if (is.null(p)) character() else unique(c(p, unlist(lapply(p, anc), use.names = FALSE)))
    assign(t, r, envir = memo)
    r
  }
  terms <- unique(t2g$term)
  A <- stats::setNames(lapply(terms, anc), terms)
  k <- lengths(A[t2g$term])
  add <- data.frame(term = unlist(A[t2g$term], use.names = FALSE), gene = rep(t2g$gene, k),
                    stringsAsFactors = FALSE)
  out <- unique(rbind(t2g, add))
  rownames(out) <- NULL
  out
}

# Fraction of annotated (gene, term) pairs whose is_a/part_of parents are also
# annotated to the gene. Near 1 = already propagated; low = most-specific terms only.
propagation_fraction <- function(t2g, edges, max_pairs = 20000L) {
  e <- edges[edges$relation %in% PROPAGATE_RELATIONS, , drop = FALSE]
  par <- split(e$parent, e$child)
  idx <- if (nrow(t2g) > max_pairs) unique(round(seq(1, nrow(t2g), length.out = max_pairs))) else seq_len(nrow(t2g))
  have <- paste(t2g$gene, t2g$term)
  ok <- vapply(idx, function(i) {
    p <- par[[t2g$term[i]]]
    if (is.null(p)) return(NA)
    all(paste(t2g$gene[i], p) %in% have)
  }, NA)
  mean(ok, na.rm = TRUE)
}

ISOFORM_RE <- "(\\.[0-9]+(\\.p)?|_[TP][0-9]+)$"

strip_isoform <- function(ids) sub(ISOFORM_RE, "", ids)

# Nomenclature hints for common crop gene ID systems (a heuristic, for messages only).
ID_STYLES <- list(
  "soybean Wm82.a1 (Glyma01g00010)" = "^Glyma[0-9]{2}g[0-9]{5}",
  "soybean Wm82.a2 or later (Glyma.01G000100)" = "^Glyma\\.[0-9]{2}G[0-9]{6}",
  "maize B73 RefGen_v3 (GRMZM2G...)" = "^(GRMZM[0-9]G[0-9]{6}|AC[0-9]+\\.[0-9]+_FG[0-9]+)",
  "maize B73 RefGen_v4 (Zm00001d...)" = "^Zm00001d[0-9]{6}",
  "maize B73 v5 (Zm00001eb...)" = "^Zm00001eb[0-9]{6}",
  "rice MSU7 (LOC_Os01g01010)" = "^LOC_Os[0-9]{2}g[0-9]{5}",
  "rice RAP-DB (Os01g0100100)" = "^Os[0-9]{2}[gt][0-9]{7}",
  "wheat IWGSC 01G-coded (TraesCS1A01G...)" = "^TraesCS[0-9U][ABDU]?01G",
  "wheat IWGSC 02G-coded (TraesCS1A02G...)" = "^TraesCS[0-9U][ABDU]?02G",
  "wheat IWGSC 03G-coded (TraesCS1A03G...)" = "^TraesCS[0-9U][ABDU]?03G",
  "Arabidopsis AGI (AT1G01010)" = "^AT[1-5CM]G[0-9]{5}",
  "Entrez Gene numeric ID" = "^[0-9]+$",
  "UniProt accession" = "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")

id_style <- function(ids) {
  ids <- utils::head(unique(ids[!is.na(ids)]), 2000)
  if (!length(ids)) return(NA_character_)
  hit <- vapply(ID_STYLES, function(p) mean(grepl(p, ids)), 0)
  if (max(hit) >= 0.5) names(ID_STYLES)[which.max(hit)] else NA_character_
}

diagnose_ids <- function(query, ann_genes) {
  base <- mean(query %in% ann_genes)
  hints <- character()
  f_iso <- mean(query %in% strip_isoform(ann_genes))
  if (f_iso > base + 0.2) hints <- c(hints, sprintf(paste0(
    "%.0f%% match after removing isoform suffixes (.1, .1.p, _T001) from the annotation IDs: the ",
    "annotation is per transcript or protein. Re-run with --strip-isoform."), 100 * f_iso))
  f_q <- mean(strip_isoform(query) %in% ann_genes)
  if (f_q > base + 0.2) hints <- c(hints, sprintf(paste0(
    "%.0f%% match after removing isoform suffixes from your gene list: the DE results are per ",
    "transcript. Summarize to genes first (claw-rnaseq-import), then re-run."), 100 * f_q))
  trim <- ann_genes
  for (k in 1:3) {
    trim <- sub("\\.[^.]*$", "", trim)
    f_t <- mean(query %in% trim)
    if (f_t > base + 0.2) {
      hints <- c(hints, sprintf(paste0("%.0f%% match after removing %d trailing '.'-separated field(s) ",
                                       "from the annotation IDs (e.g. an annotation-version suffix such as ",
                                       "'.Wm82.a4.v1'). Use the annotation column that holds plain gene IDs ",
                                       "(e.g. Phytozome locusName), or clean the IDs."), 100 * f_t, k))
      break
    }
  }
  f_case <- mean(tolower(query) %in% tolower(ann_genes))
  if (f_case > base + 0.2) hints <- c(hints, sprintf("%.0f%% match when case is ignored (e.g. Mercator/MapMan files are lower case).",
                                                      100 * f_case))
  f_pre <- mean(query %in% sub("^(gene|transcript)[:-]", "", ann_genes))
  if (f_pre > base + 0.2) hints <- c(hints, "Annotation IDs carry 'gene:'/'transcript:' prefixes; remove them.")
  sq <- id_style(query)
  sa <- id_style(ann_genes)
  if (!is.na(sq) && !is.na(sa) && sq != sa)
    hints <- c(hints, paste0("Your genes look like ", sq, " but the annotation looks like ", sa,
                             ". These are different gene-model releases or ID systems: convert with an ",
                             "official correspondence table (e.g. from SoyBase, MaizeGDB, Gramene or ",
                             "Phytozome), or get the annotation for the release you mapped to. String edits ",
                             "are not a valid conversion."))
  if (!length(hints)) hints <- "No simple fix found: the annotation probably belongs to another species, assembly or ID system."
  hints
}

# Over-representation test for one gene set. Returns the result table (terms with at
# least one DE gene) and the number of size-eligible terms used for the BH correction.
run_ora <- function(de, universe, t2g, t2n = NULL, min_size = 10, max_size = 500) {
  t2g <- unique(t2g[, c("term", "gene")])
  ext <- intersect(unique(t2g$gene), universe)
  sizes <- table(t2g$term[t2g$gene %in% ext])
  eligible <- names(sizes)[sizes >= min_size & sizes <= max_size]
  empty <- data.frame(term_id = character(), term_name = character(), n_de_in_term = integer(),
                      n_de_annotated = integer(), term_size = integer(), n_universe_annotated = integer(),
                      gene_ratio = character(), bg_ratio = character(), fold_enrichment = numeric(),
                      pvalue = numeric(), padj = numeric(), genes = character(), stringsAsFactors = FALSE)
  if (!length(eligible) || !length(intersect(de, ext)))
    return(list(table = empty, n_eligible = length(eligible), N = length(ext),
                n = length(intersect(de, ext))))
  muffle <- function(w) if (grepl("qvalue", conditionMessage(w))) invokeRestart("muffleWarning")
  e <- withCallingHandlers(suppressMessages(clusterProfiler::enricher(
    gene = de, universe = universe, TERM2GENE = t2g, TERM2NAME = if (is.null(t2n)) NA else t2n,
    pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = min_size,
    maxGSSize = max_size)), warning = muffle)
  r <- if (is.null(e)) NULL else as.data.frame(e)
  if (is.null(r) || !nrow(r))
    return(list(table = empty, n_eligible = length(eligible), N = length(ext), n = length(intersect(de, ext))))
  k <- as.integer(r$Count)
  n <- as.integer(sub("/.*", "", sub("^[0-9]+/", "", r$GeneRatio)))
  M <- as.integer(sub("/.*", "", r$BgRatio))
  N <- as.integer(sub(".*/", "", r$BgRatio))
  # BH over every size-eligible term, including terms with no DE gene (p = 1).
  m0 <- length(eligible) - nrow(r)
  padj <- stats::p.adjust(c(r$pvalue, rep(1, max(m0, 0))), method = "BH")[seq_len(nrow(r))]
  out <- data.frame(term_id = r$ID, term_name = r$Description, n_de_in_term = k, n_de_annotated = n,
                    term_size = M, n_universe_annotated = N, gene_ratio = r$GeneRatio, bg_ratio = r$BgRatio,
                    fold_enrichment = (k / n) / (M / N), pvalue = r$pvalue, padj = padj,
                    genes = r$geneID, stringsAsFactors = FALSE)
  out <- out[order(out$pvalue, -out$fold_enrichment), , drop = FALSE]
  rownames(out) <- NULL
  list(table = out, n_eligible = length(eligible), N = N[1], n = n[1])
}

# --- Figures ------------------------------------------------------------------

plot_dots <- function(tab, title, out, file, alpha, top = 15) {
  if (!nrow(tab)) return(invisible(NULL))
  d <- utils::head(tab, top)
  d <- d[rev(seq_len(nrow(d))), , drop = FALSE]
  lab <- ifelse(nchar(d$term_name) > 50, paste0(substr(d$term_name, 1, 47), "..."), d$term_name)
  lab <- ifelse(lab == d$term_id, lab, paste0(lab, " (", d$term_id, ")"))
  open_png(out, file, width = 10.5, height = 1.6 + 0.32 * nrow(d))
  op <- graphics::par(mar = c(4.5, 26, 3, 1))
  on.exit({graphics::par(op); grDevices::dev.off()})
  sig <- d$padj < alpha
  cex <- 0.8 + 2.2 * sqrt(d$n_de_in_term / max(d$n_de_in_term))
  graphics::plot(d$fold_enrichment, seq_len(nrow(d)), pch = 21, cex = cex, yaxt = "n", ylab = "",
                 bg = ifelse(sig, "#C0392B", "white"), col = ifelse(sig, "#7B241C", "#7F8C8D"),
                 xlab = "Fold enrichment", main = title,
                 xlim = c(0, max(d$fold_enrichment) * 1.08))
  graphics::axis(2, at = seq_len(nrow(d)), labels = lab, las = 1, cex.axis = 0.75)
  graphics::abline(v = 1, lty = 2, col = "grey50")
  graphics::mtext(sprintf("filled: padj < %g; size: number of DE genes in term", alpha), 3, 0.2, cex = 0.7)
}

# --- Main -------------------------------------------------------------------

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (opt$demo) {
    opt$de_table <- file.path(SKILL_DIR, "demo", "de_results.csv")
    opt$annotation <- file.path(SKILL_DIR, "demo", "go_annotation.tsv")
  }
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  if (!nzchar(opt$annotation)) fail("Give --annotation (gene-to-term table or GAF). See --help.")
  if (opt$min_size < 1 || opt$max_size < opt$min_size) fail("Need 1 <= --min-size <= --max-size.")
  warn <- character()
  notes <- character()

  # gene sets and universe
  if (nzchar(opt$de_table)) {
    if (nzchar(opt$genes) || nzchar(opt$universe)) fail("Give --de-table, or --genes with --universe, not both.")
    gs <- sets_from_de_table(opt$de_table, opt$alpha, opt$lfc_min, split_list(opt$direction))
    sets <- gs$sets
    universe <- gs$universe
    if (!gs$has_lfc && any(split_list(opt$direction) %in% c("up", "down")))
      warn <- c(warn, "The DE table has no log2FoldChange column: up/down sets were skipped.")
    notes <- c(notes, sprintf("DE genes: padj < %g%s from `%s`; universe: the %d genes with a padj (tested genes).",
                              opt$alpha, if (opt$lfc_min > 0) sprintf(" and |log2FC| >= %g", opt$lfc_min) else "",
                              basename(opt$de_table), length(universe)))
  } else {
    if (!nzchar(opt$genes)) fail("Give --de-table DE.csv, or --genes FILE with --universe FILE.")
    if (!nzchar(opt$universe))
      fail("--universe is required: list every gene that was tested (e.g. genes with a padj in the DE ",
           "results), not the whole genome. A genome background inflates enrichment for any category ",
           "of well-expressed genes (Timmons et al. 2015).")
    sets <- list(all = read_gene_list(opt$genes, "--genes"))
    universe <- read_gene_list(opt$universe, "--universe")
  }
  check_sets_in_universe(sets, universe)

  # annotation
  ann <- read_annotation(opt$annotation, opt$format, opt$gene_col, opt$term_col)
  t2g <- ann$t2g
  notes <- c(notes, ann$notes)
  if (opt$strip_isoform) {
    before <- length(unique(t2g$gene))
    t2g$gene <- strip_isoform(t2g$gene)
    t2g <- unique(t2g)
    notes <- c(notes, sprintf("--strip-isoform: %d annotation IDs merged into %d genes (union of terms).",
                              before, length(unique(t2g$gene))))
  }
  ann_genes <- unique(t2g$gene)
  cov <- data.frame(gene_set = c("universe", names(sets)),
                    n_genes = c(length(universe), lengths(sets)),
                    n_in_annotation = c(sum(universe %in% ann_genes),
                                        vapply(sets, function(s) sum(s %in% ann_genes), 1L)))
  cov$fraction <- round(cov$n_in_annotation / pmax(cov$n_genes, 1), 4)
  low <- cov$n_genes > 0 & cov$fraction < opt$min_overlap
  if (any(low)) {
    q <- if (low[1]) universe else sets[[cov$gene_set[which(low)[1]]]]
    fail(sprintf("Only %.1f%% of the %s genes are in the annotation (--min-overlap %g). ",
                 100 * cov$fraction[which(low)[1]], cov$gene_set[which(low)[1]], opt$min_overlap),
         "Gene IDs probably do not match.\n  Your IDs, e.g.: ", paste(utils::head(q, 3), collapse = ", "),
         "\n  Annotation IDs, e.g.: ", paste(utils::head(ann_genes, 3), collapse = ", "),
         "\n  Diagnosis:\n  - ", paste(diagnose_ids(q, ann_genes), collapse = "\n  - "))
  }
  if (cov$fraction[1] < 0.5)
    warn <- c(warn, sprintf(paste0("Only %.0f%% of tested genes have any annotation: results describe the ",
                                   "annotated minority, and terms for poorly annotated processes will be missed."),
                            100 * cov$fraction[1]))
  hint <- diagnose_ids(universe, ann_genes)
  if (cov$fraction[1] < 0.9 && !grepl("^No simple fix", hint[1]))
    warn <- c(warn, paste("ID check:", hint[1]))

  # GO handling
  go <- is_go(t2g$term)
  ontologies <- "custom"
  t2n <- read_term_names(opt$term_names)
  if (is.null(t2n)) t2n <- ann$t2n
  prop_frac <- NA_real_
  if (go) {
    gi <- go_term_info(t2g$term)
    if (length(gi$remap)) {
      hit <- t2g$term %in% names(gi$remap)
      t2g$term[hit] <- gi$remap[t2g$term[hit]]
      t2g <- unique(t2g)
      notes <- c(notes, sprintf("%d secondary GO IDs were mapped to their primary IDs (GO.db GOSYNONYM).",
                                length(gi$remap)))
    }
    if (length(gi$obsolete)) {
      n_pairs <- sum(t2g$term %in% gi$obsolete)
      t2g <- t2g[!t2g$term %in% gi$obsolete, , drop = FALSE]
      warn <- c(warn, sprintf(paste0("%d GO ID(s) in the annotation are not in GO.db %s (obsolete or newer ",
                                     "than this GO release; e.g. %s): %d gene-term pairs dropped."),
                              length(gi$obsolete), as.character(utils::packageVersion("GO.db")),
                              paste(utils::head(gi$obsolete, 3), collapse = ", "), n_pairs))
    }
    edges <- go_parent_edges()
    prop_frac <- propagation_fraction(t2g, edges)
    if (opt$propagate) {
      n0 <- nrow(t2g)
      t2g <- propagate_annotation(t2g, edges)
      notes <- c(notes, sprintf("Propagation over is_a/part_of (GO.db %s): %d -> %d gene-term pairs.",
                                as.character(utils::packageVersion("GO.db")), n0, nrow(t2g)))
    } else if (is.finite(prop_frac) && prop_frac < 0.5) {
      warn <- c(warn, sprintf(paste0("Only %.0f%% of annotations have their parent terms listed too: this file ",
                                     "holds mostly the most specific GO terms, so broader terms are under-counted. ",
                                     "Re-run with --propagate."), 100 * prop_frac))
    }
    info <- go_term_info(t2g$term)$info
    ont_of <- stats::setNames(info$ONTOLOGY, info$GOID)
    name_of <- stats::setNames(info$TERM, info$GOID)
    t2n <- data.frame(term = names(name_of), name = unname(name_of), stringsAsFactors = FALSE)
    ontologies <- intersect(toupper(split_list(opt$ontology)), c("BP", "MF", "CC"))
    if (!length(ontologies)) fail("--ontology must list BP, MF and/or CC.")
  }

  # run
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  runs <- list()
  for (ont in ontologies) {
    sub_t2g <- if (ont == "custom") t2g else t2g[ont_of[t2g$term] %in% ont, , drop = FALSE]
    for (nm in names(sets)) {
      lab <- if (ont == "custom") nm else paste(nm, ont, sep = "_")
      if (!length(sets[[nm]])) {
        warn <- c(warn, sprintf("Gene set '%s' is empty (no genes at padj < %g): not tested.", nm, opt$alpha))
        next
      }
      r <- run_ora(sets[[nm]], universe, sub_t2g, t2n, opt$min_size, opt$max_size)
      tab <- r$table
      tab$ontology <- if (ont == "custom") "custom" else ont
      tab$significant <- tab$padj < opt$alpha
      tab <- tab[, c("term_id", "term_name", "ontology", setdiff(names(tab), c("term_id", "term_name", "ontology")))]
      save_table(tab, opt$out, paste0("enrichment_", lab, ".csv"))
      plot_dots(tab, sprintf("%s genes, %s", nm, if (ont == "custom") "gene sets" else paste("GO", ont)),
                opt$out, paste0("dotplot_", lab, ".png"), opt$alpha)
      runs[[lab]] <- list(set = nm, ont = ont, table = tab, n_eligible = r$n_eligible, N = r$N, n = r$n,
                          n_input = length(sets[[nm]]))
      if (r$n < 10) warn <- c(warn, sprintf("%s: only %d annotated DE genes; enrichment has little power.", lab, r$n))
    }
  }
  if (!length(runs)) fail("No gene set could be tested.")
  for (nm in names(sets)) if (length(sets[[nm]]) > 0.5 * length(universe))
    warn <- c(warn, sprintf("'%s' contains %.0f%% of the universe: over-representation analysis is not informative when most genes are DE.",
                            nm, 100 * length(sets[[nm]]) / length(universe)))
  save_table(cov, opt$out, "annotation_coverage.csv")
  run_summary <- do.call(rbind, lapply(names(runs), function(l) {
    x <- runs[[l]]
    data.frame(analysis = l, gene_set = x$set, ontology = x$ont, n_genes = x$n_input,
               n_genes_annotated = x$n, n_universe_annotated = x$N, n_terms_tested = x$n_eligible,
               n_terms_significant = sum(x$table$significant), stringsAsFactors = FALSE)
  }))
  save_table(run_summary, opt$out, "run_summary.csv")

  rep <- c(sprintf("# %s report", SKILL), "",
           if (length(warn)) c("## Warnings", "", paste("-", warn), "") else NULL,
           "## Inputs", "", paste("-", notes),
           sprintf("- Annotation: `%s` (%d genes, %d gene-term pairs%s).", basename(opt$annotation),
                   length(unique(t2g$gene)), nrow(t2g), if (go) ", Gene Ontology" else ", custom gene sets"),
           if (go && is.finite(prop_frac)) sprintf("- Share of input annotations whose parent terms were already listed: %.0f%%.", 100 * prop_frac) else NULL,
           sprintf("- Gene-set size limits: %g to %g genes within the annotated universe.", opt$min_size, opt$max_size),
           "", "| gene set | genes | in annotation | fraction |", "|---|---|---|---|",
           sprintf("| %s | %d | %d | %.2f |", cov$gene_set, cov$n_genes, cov$n_in_annotation, cov$fraction), "",
           "## Results", "",
           sprintf("Hypergeometric test; BH-adjusted over all %s tested terms (padj < %g = significant).",
                   "size-eligible", opt$alpha), "",
           "| analysis | DE genes (annotated) | universe (annotated) | terms tested | significant |", "|---|---|---|---|---|",
           sprintf("| %s | %d (%d) | %d | %d | %d |", run_summary$analysis, run_summary$n_genes,
                   run_summary$n_genes_annotated, run_summary$n_universe_annotated, run_summary$n_terms_tested,
                   run_summary$n_terms_significant), "")
  for (l in names(runs)) {
    tab <- runs[[l]]$table
    top <- utils::head(tab[tab$significant, , drop = FALSE], 10)
    rep <- c(rep, sprintf("### %s", l), "")
    if (nrow(top)) {
      rep <- c(rep, "| term | name | DE in term / term size | fold enrichment | padj |", "|---|---|---|---|---|",
               sprintf("| %s | %s | %d / %d | %.2f | %.2e |", top$term_id, top$term_name, top$n_de_in_term,
                       top$term_size, top$fold_enrichment, top$padj), "",
               sprintf("![%s](figures/dotplot_%s.png)", l, l), "")
    } else rep <- c(rep, sprintf("No term at padj < %g.", opt$alpha), "")
  }
  rep <- c(rep, "## How to read the tables", "",
           "- `n_de_in_term` (k) of `n_de_annotated` (n) DE genes fall in the term; the term has `term_size` (M) of `n_universe_annotated` (N) annotated universe genes. p = P(X >= k), X ~ hypergeometric(M, N - M, n).",
           "- `fold_enrichment` = (k/n)/(M/N). `padj` is Benjamini-Hochberg over every term that passed the size limits, including terms without DE genes (p = 1), so it is slightly more conservative than clusterProfiler's default.",
           "- Only terms with at least one DE gene are listed.",
           "- Genes without any annotation are outside the test (N and n count annotated genes only).", "",
           "## Methods", "",
           paste0(sprintf("Over-representation was tested with clusterProfiler %s (enricher, hypergeometric test) ",
                          as.character(utils::packageVersion("clusterProfiler"))),
                  "against a background of the genes tested for differential expression. ",
                  if (go) sprintf("GO term names and namespaces came from GO.db %s%s. ",
                                  as.character(utils::packageVersion("GO.db")),
                                  if (opt$propagate) "; annotations were propagated to ancestor terms over is_a and part_of relations" else "") else "",
                  sprintf("Gene sets with %g to %g annotated universe genes were tested and p-values adjusted by Benjamini-Hochberg.",
                          opt$min_size, opt$max_size)),
           "", "## References", "",
           "- Wu T, Hu E, Xu S, et al. (2021). clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. The Innovation 2:100141. https://doi.org/10.1016/j.xinn.2021.100141",
           "- Ashburner M, Ball CA, Blake JA, et al. (2000). Gene Ontology: tool for the unification of biology. Nature Genetics 25:25-29. https://doi.org/10.1038/75556",
           "- The Gene Ontology Consortium (2023). The Gene Ontology knowledgebase in 2023. Genetics 224:iyad031. https://doi.org/10.1093/genetics/iyad031",
           "- Timmons JA, Szkop KJ, Gallagher IJ (2015). Multiple sources of bias confound functional enrichment analysis of global -omics data. Genome Biology 16:186. https://doi.org/10.1186/s13059-015-0761-7",
           "- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate. Journal of the Royal Statistical Society B 57:289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x",
           "", sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(rep, file.path(opt$out, "report.md"))
  write_repro_bundle(opt$out, argv, inputs = c(opt$de_table, opt$genes, opt$universe, opt$annotation, opt$term_names))
  cat(sprintf("Done: %s. Report: %s\n", paste(sprintf("%s %d", run_summary$analysis, run_summary$n_terms_significant),
                                               collapse = ", "),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(runs)
}

if (sys.nframe() == 0L) main()
