#!/usr/bin/env Rscript
# claw-rnaseq-import: import Salmon, kallisto or featureCounts quantifications
# into gene-level counts with tximport, build the transcript-to-gene table from
# a GFF3/GTF (or a user table), check that transcript IDs match, and run sample
# QC (mapping rate, library size, detected genes, correlations, PCA on a
# variance-stabilised transform, outlier flags).
#
# Usage:
#   Rscript scripts/rnaseq_import.R --demo
#   Rscript scripts/rnaseq_import.R --samples samples.csv --gff annotation.gff3
#   Rscript scripts/rnaseq_import.R --samples samples.csv --tx2gene tx2gene.csv
#   Rscript scripts/rnaseq_import.R --samples samples.csv --type featurecounts --counts counts.txt
#   Rscript scripts/rnaseq_import.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-rnaseq-import"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/rnaseq_import.R [--demo] --samples SAMPLES.csv
         (--gff ANNOTATION.gff3|.gtf | --tx2gene TX2GENE.csv | --type featurecounts --counts FILE)
         [options]

  --demo               run on the bundled synthetic Salmon demo (demo/)
  --samples            sample sheet CSV: columns 'sample' and 'path' (a quant.sf or
                       abundance.tsv file, or the folder holding it) plus design columns
                       such as condition and batch. Relative paths are read relative to
                       the sample sheet's folder. For featureCounts, 'path' is not
                       needed; an optional 'column' names the matrix column per sample.
  --type               auto, salmon, kallisto or featurecounts (default: auto)
  --gff                GFF3 or GTF annotation (plain or .gz) used to build tx2gene
  --tx2gene            CSV/TSV with transcript_id and gene_id columns (instead of --gff)
  --counts             featureCounts output matrix (with --type featurecounts)
  --fc-summary         featureCounts .summary file (optional; gives the assignment rate)
  --id-attr            GFF3 attribute holding transcript IDs: auto, ID, Name or
                       transcript_id (default: auto = the one matching the quant files)
  --strip-tx-version   remove a trailing '.N' version from transcript IDs on both sides.
                       Refused when '.N' is an isoform number (stripping creates duplicates)
  --min-id-match       minimum fraction of quantified transcripts found in tx2gene
                       (default: 0.95)
  --min-mapping-rate   flag samples below this mapping rate, in percent (default: 50)
  --outlier-mad        flag samples more than this many MADs below the median for
                       mapping rate, log library size, detected genes and median
                       correlation (default: 3)
  --detect-min-count   a gene counts as detected at this many estimated reads (default: 10)
  --transform          vst or rlog for the QC PCA and correlations (default: vst)
  --three-prime        3' tag-seq library (e.g. QuantSeq): counts have no length bias, so
                       no transcript-length offsets are used (tximport vignette)
  --color-by           sample-sheet column used to colour the PCA (default: first
                       design column)
  --out                output folder (default: agriclaw_output/claw-rnaseq-import)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, samples = "", type = "auto", gff = "", tx2gene = "",
                 counts = "", fc_summary = "", id_attr = "auto",
                 strip_tx_version = FALSE, min_id_match = 0.95,
                 min_mapping_rate = 50, outlier_mad = 3, detect_min_count = 10,
                 transform = "vst", three_prime = FALSE, color_by = "", out = "")

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

# Columns of the sample sheet that are bookkeeping, not design variables.
SHEET_SYSTEM_COLS <- c("sample", "path", "column")

read_sample_sheet <- function(path) {
  if (!nzchar(path)) fail("Provide --samples SAMPLES.csv (or --demo). See --help.")
  if (!file.exists(path)) fail("Sample sheet not found: ", path)
  s <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                       strip.white = TRUE)
  if (!"sample" %in% names(s))
    fail("The sample sheet needs a 'sample' column. Found: ", paste(names(s), collapse = ", "))
  s$sample <- as.character(s$sample)
  if (any(!nzchar(s$sample) | is.na(s$sample))) fail("The 'sample' column has empty entries.")
  if (anyDuplicated(s$sample))
    fail("Duplicated sample names: ", paste(unique(s$sample[duplicated(s$sample)]), collapse = ", "))
  bad <- s$sample[make.names(s$sample) != s$sample]
  if (length(bad))
    message("Note: sample names that are not syntactic R names are kept as given: ",
            paste(utils::head(bad, 5), collapse = ", "))
  s
}

is_absolute_path <- function(p) grepl("^(/|\\\\|[A-Za-z]:[/\\\\]|~)", p)

# Resolve each sample's quantification file and detect the quantifier.
resolve_quant_files <- function(sheet, base_dir, type = "auto") {
  if (!"path" %in% names(sheet))
    fail("The sample sheet needs a 'path' column pointing to each sample's quant.sf ",
         "(Salmon) or abundance.tsv (kallisto), or to the folder that holds it.")
  p <- as.character(sheet$path)
  p <- ifelse(is_absolute_path(p), p, file.path(base_dir, p))
  files <- character(length(p))
  types <- character(length(p))
  for (i in seq_along(p)) {
    f <- p[[i]]
    if (dir.exists(f)) {
      cand <- file.path(f, c("quant.sf", "quant.sf.gz", "abundance.tsv", "abundance.tsv.gz"))
      cand <- cand[file.exists(cand)]
      if (!length(cand))
        fail("No quant.sf or abundance.tsv inside folder '", f, "' (sample ", sheet$sample[[i]], ").")
      f <- cand[[1L]]
    }
    if (!file.exists(f)) fail("Quantification file not found for sample ", sheet$sample[[i]], ": ", f)
    b <- tolower(basename(f))
    types[[i]] <- if (grepl("^quant\\.sf", b)) "salmon" else if (grepl("^abundance\\.tsv", b))
      "kallisto" else if (grepl("\\.h5$", b)) "kallisto-h5" else "unknown"
    files[[i]] <- normalizePath(f, winslash = "/")
  }
  if (any(types == "kallisto-h5"))
    fail("kallisto .h5 files need the rhdf5 package; point 'path' at abundance.tsv instead.")
  if (type == "auto") {
    if (any(types == "unknown"))
      fail("Cannot tell the quantifier from file name(s): ",
           paste(basename(files[types == "unknown"]), collapse = ", "),
           ". Use --type salmon or --type kallisto.")
    if (length(unique(types)) > 1L)
      fail("The sample sheet mixes Salmon and kallisto files; quantify all samples the same way.")
    type <- types[[1L]]
  }
  names(files) <- sheet$sample
  list(files = files, type = type)
}

# Transcript IDs in the first column of a quant file.
read_quant_ids <- function(file) {
  x <- utils::read.delim(file, colClasses = "character", check.names = FALSE)[[1L]]
  x
}

# Mapping statistics written next to the quant file by Salmon or kallisto.
read_mapping_info <- function(file, type) {
  d <- dirname(file)
  out <- list(n_processed = NA_real_, n_mapped = NA_real_, mapping_rate_pct = NA_real_,
              library_type = NA_character_, source = NA_character_)
  if (type == "salmon") {
    j <- file.path(d, "aux_info", "meta_info.json")
    if (!file.exists(j)) return(out)
    m <- jsonlite::fromJSON(j)
    out$n_processed <- as.numeric(m$num_processed %||% NA)
    out$n_mapped <- as.numeric(m$num_mapped %||% NA)
    out$mapping_rate_pct <- as.numeric(m$percent_mapped %||% NA)
    if (is.na(out$mapping_rate_pct) && !is.na(out$n_mapped) && isTRUE(out$n_processed > 0))
      out$mapping_rate_pct <- 100 * out$n_mapped / out$n_processed
    lt <- m$library_types %||% NA_character_
    out$library_type <- paste(lt, collapse = ",")
    out$source <- "aux_info/meta_info.json"
  } else if (type == "kallisto") {
    j <- file.path(d, "run_info.json")
    if (!file.exists(j)) return(out)
    m <- jsonlite::fromJSON(j)
    out$n_processed <- as.numeric(m$n_processed %||% NA)
    out$n_mapped <- as.numeric(m$n_pseudoaligned %||% NA)
    out$mapping_rate_pct <- as.numeric(m$p_pseudoaligned %||% NA)
    out$source <- "run_info.json"
  }
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# --- Annotation: GFF3 / GTF to transcript-to-gene candidates -----------------

EXON_LIKE <- c("exon", "CDS", "five_prime_UTR", "three_prime_UTR", "UTR", "5UTR", "3UTR",
               "start_codon", "stop_codon", "intron", "polypeptide", "protein",
               "Selenocysteine")

strip_type_prefix <- function(x) sub("^(gene|transcript|mRNA|rna|cds|exon)[:-]", "", x)

attr_value <- function(attr, key, gtf = FALSE) {
  val <- rep(NA_character_, length(attr))
  if (gtf) {
    pat <- paste0("(^|;)\\s*", key, "\\s+\"([^\"]*)\"")
    has <- grepl(pat, attr)
    val[has] <- sub(paste0("^(.*;)?\\s*", key, "\\s+\"([^\"]*)\".*$"), "\\2", attr[has])
  } else {
    pat <- paste0("(^|;)\\s*", key, "=")
    has <- grepl(pat, attr)
    val[has] <- sub(paste0("^(.*;)?\\s*", key, "=([^;]*).*$"), "\\2", attr[has])
    enc <- !is.na(val) & grepl("%", val, fixed = TRUE)
    if (any(enc)) val[enc] <- vapply(val[enc], utils::URLdecode, "")
  }
  val
}

read_feature_lines <- function(path) {
  con <- file(path, open = "r")  # file() reads .gz transparently
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  fa <- which(startsWith(lines, "##FASTA"))
  if (length(fa)) lines <- lines[seq_len(fa[[1L]] - 1L)]
  lines[nzchar(lines) & !startsWith(lines, "#")]
}

# Returns list(format, candidates = named list of data.frames(transcript_id, gene_id)).
parse_annotation <- function(path) {
  if (!file.exists(path)) fail("Annotation file not found: ", path)
  lines <- read_feature_lines(path)
  if (!length(lines)) fail("No feature lines in annotation file: ", path)
  ncol9 <- lengths(regmatches(utils::head(lines, 200), gregexpr("\t", utils::head(lines, 200))))
  if (any(ncol9 != 8L))
    fail("Annotation file '", basename(path), "' is not tab-separated with 9 columns ",
         "(GFF3/GTF). Check the file format.")
  type_col <- sub("^[^\t]*\t[^\t]*\t([^\t]*)\t.*$", "\\1", lines)
  attr_col <- sub("^([^\t]*\t){8}", "", lines)
  is_gtf <- mean(grepl("^\\s*[A-Za-z_]+ \"", utils::head(attr_col, 200))) > 0.5
  if (is_gtf) {
    use <- type_col == "transcript"
    if (!any(use)) use <- type_col == "exon"
    tx <- attr_value(attr_col[use], "transcript_id", gtf = TRUE)
    gn <- attr_value(attr_col[use], "gene_id", gtf = TRUE)
    d <- unique(data.frame(transcript_id = tx, gene_id = gn, stringsAsFactors = FALSE))
    d <- d[!is.na(d$transcript_id) & !is.na(d$gene_id), , drop = FALSE]
    if (!nrow(d)) fail("The GTF has no transcript_id/gene_id attributes.")
    return(list(format = "GTF", candidates = list(transcript_id = d)))
  }
  keep <- !type_col %in% EXON_LIKE
  type_col <- type_col[keep]
  attr_col <- attr_col[keep]
  id <- attr_value(attr_col, "ID")
  parent <- attr_value(attr_col, "Parent")
  name <- attr_value(attr_col, "Name")
  txid <- attr_value(attr_col, "transcript_id")
  gid <- attr_value(attr_col, "gene_id")
  multi <- !is.na(parent) & grepl(",", parent, fixed = TRUE)
  parent <- sub(",.*$", "", parent)
  top <- is.na(parent) & !is.na(id)
  gene_rows <- which(top)
  tx_rows <- which(!is.na(parent) & !is.na(id) & parent %in% id[gene_rows])
  if (!length(tx_rows))
    fail("No transcript features found in the GFF3: expected features (mRNA, transcript, ",
         "ncRNA, ...) whose Parent is a gene. Check the file, or supply --tx2gene.")
  g <- match(parent[tx_rows], id[gene_rows])
  gi <- gene_rows[g]
  gene_clean <- strip_type_prefix(id[gi])
  mk <- function(tx, gene) {
    d <- unique(data.frame(transcript_id = tx, gene_id = gene, stringsAsFactors = FALSE))
    d[!is.na(d$transcript_id) & !is.na(d$gene_id), , drop = FALSE]
  }
  cands <- list(
    ID_clean = mk(strip_type_prefix(id[tx_rows]), gene_clean),
    transcript_id = mk(txid[tx_rows], ifelse(is.na(gid[gi]), gene_clean, gid[gi])),
    Name = mk(name[tx_rows], ifelse(is.na(name[gi]), gene_clean, name[gi])),
    ID = mk(id[tx_rows], id[gi]))
  cands <- cands[vapply(cands, nrow, 1L) > 0L]
  list(format = "GFF3", candidates = cands, n_multi_parent = sum(multi[tx_rows]),
       feature_types = sort(table(type_col[tx_rows]), decreasing = TRUE))
}

read_tx2gene_table <- function(path) {
  if (!file.exists(path)) fail("tx2gene file not found: ", path)
  first <- readLines(path, n = 1L, warn = FALSE)
  sep <- if (grepl("\t", first, fixed = TRUE)) "\t" else ","
  d <- utils::read.table(path, sep = sep, header = TRUE, quote = "\"", comment.char = "",
                         colClasses = "character", check.names = FALSE)
  if (ncol(d) < 2L) fail("The tx2gene table needs two columns: transcript_id and gene_id.")
  tx_col <- if ("transcript_id" %in% names(d)) "transcript_id" else names(d)[[1L]]
  gn_col <- if ("gene_id" %in% names(d)) "gene_id" else names(d)[[2L]]
  out <- unique(data.frame(transcript_id = d[[tx_col]], gene_id = d[[gn_col]],
                           stringsAsFactors = FALSE))
  out <- out[nzchar(out$transcript_id) & nzchar(out$gene_id), , drop = FALSE]
  dup <- unique(out$transcript_id[duplicated(out$transcript_id)])
  if (length(dup))
    fail("tx2gene maps ", length(dup), " transcript(s) to more than one gene, e.g. ",
         paste(utils::head(dup, 3), collapse = ", "), ". Each transcript must map to one gene.")
  out
}

# Pick the candidate whose transcript IDs match the quantified IDs best.
choose_tx2gene <- function(cands, quant_ids, id_attr = "auto") {
  frac <- vapply(cands, function(d) mean(quant_ids %in% d$transcript_id), 0)
  if (id_attr != "auto") {
    key <- if (id_attr == "ID") "ID" else id_attr
    if (!key %in% names(cands))
      fail("--id-attr ", id_attr, " is not available in this annotation. Available: ",
           paste(names(cands), collapse = ", "))
    return(list(tx2gene = cands[[key]], style = key, match = frac))
  }
  best <- names(cands)[which.max(frac)]
  list(tx2gene = cands[[best]], style = best, match = frac)
}

has_version_suffix <- function(x) grepl("\\.[0-9]+$", x)

strip_tx_version <- function(ids, label) {
  s <- sub("\\.[0-9]+$", "", ids)
  dup <- duplicated(s) & !duplicated(ids)
  if (any(dup)) {
    ex <- s[dup][1L]
    fail("--strip-tx-version would merge different ", label, " (e.g. ",
         paste(utils::head(ids[s == ex], 3), collapse = " and "), " both become '", ex,
         "'). Here the trailing '.N' is an isoform number, not a version suffix, so it must ",
         "not be removed. The quantification and the annotation come from different ",
         "annotation versions: rebuild the Salmon/kallisto index from the same annotation.")
  }
  s
}

# Nomenclature hints for common crop gene/transcript ID systems.
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
  "NCBI RefSeq transcript (XM_/NM_)" = "^[XN][MR]_[0-9]+",
  "NCBI LOC gene symbol (LOC100...)" = "^LOC[0-9]+$")

id_style <- function(ids) {
  ids <- utils::head(unique(ids[!is.na(ids)]), 2000)
  if (!length(ids)) return(NA_character_)
  hit <- vapply(ID_STYLES, function(p) mean(grepl(p, ids, ignore.case = FALSE)), 0)
  if (max(hit) >= 0.5) names(ID_STYLES)[which.max(hit)] else NA_character_
}

diagnose_id_mismatch <- function(quant_ids, tx_ids) {
  base <- mean(quant_ids %in% tx_ids)
  hints <- character()
  tryfix <- function(a, b) mean(a %in% b)
  sv <- function(x) sub("\\.[0-9]+$", "", x)
  checks <- list(
    list(f = tryfix(sv(quant_ids), tx_ids),
         msg = "the quantification IDs carry a trailing '.N' suffix that the annotation IDs lack"),
    list(f = tryfix(quant_ids, sv(tx_ids)),
         msg = "the annotation IDs carry a trailing '.N' suffix that the quantification IDs lack"),
    list(f = tryfix(sv(quant_ids), sv(tx_ids)),
         msg = "the trailing '.N' suffixes differ between quantification and annotation"))
  for (ck in checks) {
    if (ck$f > base + 0.2) {
      iso <- anyDuplicated(sv(quant_ids)) > 0L
      hints <- c(hints, paste0(
        "IDs match ", sprintf("%.0f%%", 100 * ck$f), " after removing a trailing '.N': ",
        ck$msg, ". ",
        if (iso) paste0("But in these IDs '.N' distinguishes isoforms (removing it merges ",
                        "transcripts), so this is an annotation-version mismatch: rebuild the ",
                        "index from the same annotation as the GFF/tx2gene.")
        else "If '.N' is an Ensembl-style version number, re-run with --strip-tx-version."))
      break
    }
  }
  fl <- tryfix(tolower(quant_ids), tolower(tx_ids))
  if (fl > base + 0.2) hints <- c(hints, sprintf(
    "IDs match %.0f%% when case is ignored: upper/lower case differs between files.", 100 * fl))
  fp <- tryfix(quant_ids, strip_type_prefix(tx_ids))
  if (fp > base + 0.2) hints <- c(hints, sprintf(paste0(
    "IDs match %.0f%% after removing 'transcript:'/'gene:' style prefixes from the ",
    "annotation IDs: use the GFF3 directly with --gff, or clean the tx2gene table."), 100 * fp))
  if (mean(grepl("|", utils::head(quant_ids, 500), fixed = TRUE)) > 0.5)
    hints <- c(hints, paste0("Quantification IDs contain '|' (GENCODE-style FASTA headers). ",
                             "Re-run Salmon with --gencode, or index a FASTA with plain IDs."))
  sq <- id_style(quant_ids)
  st <- id_style(tx_ids)
  if (!is.na(sq) && !is.na(st) && sq != st)
    hints <- c(hints, paste0("Quantification IDs look like ", sq, " but annotation IDs look like ",
                             st, ". These are different annotation versions; string edits ",
                             "cannot convert one to the other. Use the annotation that was ",
                             "used to build the index."))
  if (!length(hints))
    hints <- paste0("No simple fix found. The index and the annotation probably come from ",
                    "different sources or assembly versions (e.g. Ensembl Plants vs Phytozome).")
  hints
}

# --- Import -------------------------------------------------------------------

import_transcripts <- function(files, type) {
  ids1 <- read_quant_ids(files[[1L]])
  for (i in seq_along(files)[-1L]) {
    ids <- read_quant_ids(files[[i]])
    if (length(ids) != length(ids1) || any(ids != ids1))
      fail("Sample '", names(files)[[i]], "' was quantified against a different transcript ",
           "set or order than sample '", names(files)[[1L]], "' (", length(ids), " vs ",
           length(ids1), " transcripts). All samples must use the same Salmon/kallisto index.")
  }
  if (anyDuplicated(ids1)) fail("Duplicated transcript IDs in ", files[[1L]])
  suppressMessages(tximport::tximport(files, type = type, txOut = TRUE,
                                      countsFromAbundance = "no"))
}

# Summarise transcript-level estimates to genes. Returns the gene-level txi plus
# ID-matching diagnostics. Stops when too few transcripts match.
summarise_to_genes <- function(txi_tx, tx2gene, strip_version = FALSE, min_id_match = 0.95) {
  q_ids <- rownames(txi_tx$counts)
  if (strip_version) {
    q_ids <- strip_tx_version(q_ids, "quantified transcripts")
    tx2gene$transcript_id <- strip_tx_version(tx2gene$transcript_id, "annotation transcripts")
    for (m in c("abundance", "counts", "length")) rownames(txi_tx[[m]]) <- q_ids
  }
  matched <- q_ids %in% tx2gene$transcript_id
  frac <- mean(matched)
  reads_lost <- colSums(txi_tx$counts[!matched, , drop = FALSE]) /
    pmax(colSums(txi_tx$counts), 1e-12)
  if (frac < min_id_match) {
    hints <- diagnose_id_mismatch(q_ids, tx2gene$transcript_id)
    fail(sprintf("Only %.1f%% of quantified transcripts (%d of %d) are in tx2gene ",
                 100 * frac, sum(matched), length(q_ids)),
         sprintf("(--min-id-match %.2f).\n  Quantified IDs, e.g.: %s\n  Annotation IDs, e.g.: %s\n  Diagnosis:\n  - ",
                 min_id_match, paste(utils::head(q_ids, 3), collapse = ", "),
                 paste(utils::head(tx2gene$transcript_id, 3), collapse = ", ")),
         paste(hints, collapse = "\n  - "))
  }
  t2g <- tx2gene[tx2gene$transcript_id %in% q_ids, c("transcript_id", "gene_id")]
  txi <- suppressMessages(tximport::summarizeToGene(txi_tx, tx2gene = t2g,
                                                     countsFromAbundance = "no"))
  list(txi = txi, tx2gene = t2g, match_frac = frac, unmatched = q_ids[!matched],
       reads_lost_frac = reads_lost)
}

read_featurecounts <- function(path) {
  if (!file.exists(path)) fail("featureCounts file not found: ", path)
  d <- utils::read.delim(path, comment.char = "#", check.names = FALSE,
                         stringsAsFactors = FALSE)
  need <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  if (!all(need %in% names(d)[1:6]))
    fail("'", basename(path), "' does not look like featureCounts output (expected columns ",
         paste(need, collapse = ", "), " then one column per BAM).")
  m <- as.matrix(d[, -(1:6), drop = FALSE])
  if (!is.numeric(m) || any(m < 0, na.rm = TRUE) || any(is.na(m)))
    fail("featureCounts count columns must be non-negative numbers without NA.")
  if (any(m != round(m))) fail("featureCounts counts must be integers; got non-integer values.")
  rownames(m) <- d$Geneid
  storage.mode(m) <- "double"
  list(counts = m, length = stats::setNames(d$Length, d$Geneid))
}

clean_bam_name <- function(x) {
  x <- basename(gsub("\\\\", "/", x))
  x <- sub("\\.bam$", "", x, ignore.case = TRUE)
  x <- sub("(\\.sorted|\\.sortedByCoord|Aligned\\.sortedByCoord\\.out|Aligned\\.out)$", "", x)
  x
}

match_fc_columns <- function(sheet, cols) {
  if ("column" %in% names(sheet)) {
    want <- as.character(sheet$column)
    miss <- setdiff(want, cols)
    if (length(miss)) fail("Sample-sheet 'column' entries not found in the featureCounts matrix: ",
                           paste(miss, collapse = ", "))
    return(stats::setNames(want, sheet$sample))
  }
  cleaned <- clean_bam_name(cols)
  idx <- match(sheet$sample, cleaned)
  if (anyNA(idx) || anyDuplicated(cleaned))
    fail("Cannot match samples to featureCounts columns by name. Add a 'column' column to the ",
         "sample sheet with the exact matrix header for each sample. Matrix columns: ",
         paste(utils::head(cols, 6), collapse = ", "))
  stats::setNames(cols[idx], sheet$sample)
}

read_fc_summary <- function(path, colmap) {
  if (!nzchar(path)) return(NULL)
  if (!file.exists(path)) fail("featureCounts summary not found: ", path)
  d <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  tot <- colSums(d[, -1, drop = FALSE])
  ass <- unlist(d[d$Status == "Assigned", -1, drop = FALSE])
  out <- data.frame(sample = names(colmap), n_processed = tot[colmap], n_mapped = ass[colmap],
                    stringsAsFactors = FALSE)
  out$mapping_rate_pct <- 100 * out$n_mapped / out$n_processed
  out
}

# --- QC ---------------------------------------------------------------------

qc_transform <- function(counts, lengths = NULL, method = "vst") {
  keep <- rowSums(counts) > 0
  cd <- data.frame(sample = colnames(counts), row.names = colnames(counts))
  if (is.null(lengths)) {
    cnt <- round(counts[keep, , drop = FALSE])
    storage.mode(cnt) <- "integer"
    dds <- DESeq2::DESeqDataSetFromMatrix(cnt, cd, design = ~ 1)
  } else {
    txi <- list(counts = counts[keep, , drop = FALSE], length = lengths[keep, , drop = FALSE],
                abundance = NULL, countsFromAbundance = "no")
    dds <- suppressMessages(DESeq2::DESeqDataSetFromTximport(txi, cd, design = ~ 1))
  }
  dds <- suppressMessages(DESeq2::estimateSizeFactors(dds))
  if (method == "rlog") {
    tr <- suppressMessages(DESeq2::rlog(dds, blind = TRUE))
  } else {
    n_ok <- sum(rowMeans(DESeq2::counts(dds, normalized = TRUE)) > 5)
    tr <- if (n_ok >= 1000) suppressMessages(DESeq2::vst(dds, blind = TRUE)) else
      suppressMessages(DESeq2::varianceStabilizingTransformation(dds, blind = TRUE))
  }
  SummarizedExperiment::assay(tr)
}

# PCA as in DESeq2::plotPCA: top-variance genes, centred, not scaled.
qc_pca <- function(mat, ntop = 500) {
  rv <- apply(mat, 1, stats::var)
  sel <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  p <- stats::prcomp(t(mat[sel, , drop = FALSE]), center = TRUE, scale. = FALSE)
  ve <- p$sdev^2 / sum(p$sdev^2)
  list(scores = p$x, var_explained = ve, n_genes = length(sel))
}

# Lower-tail robust outlier rule: x < median - k * MAD.
mad_low <- function(x, k) {
  ok <- !is.na(x)
  out <- rep(FALSE, length(x))
  if (sum(ok) < 3L) return(out)
  m <- stats::median(x[ok])
  s <- stats::mad(x[ok])
  if (!is.finite(s) || s == 0) return(out)
  out[ok] <- x[ok] < m - k * s
  out
}

sample_qc <- function(counts, mapping, tmat, detect_min = 10, min_map = 50, k = 3) {
  s <- colnames(counts)
  cm <- stats::cor(tmat)
  med_cor <- vapply(seq_along(s), function(i) stats::median(cm[i, -i]), 0)
  qc <- data.frame(sample = s,
                   n_processed = mapping$n_processed[match(s, mapping$sample)],
                   n_mapped = mapping$n_mapped[match(s, mapping$sample)],
                   mapping_rate_pct = mapping$mapping_rate_pct[match(s, mapping$sample)],
                   assigned_reads = colSums(counts),
                   detected_genes = colSums(counts >= detect_min),
                   median_cor = med_cor, stringsAsFactors = FALSE)
  flags <- vector("list", length(s))
  add <- function(i, msg) for (j in which(i)) flags[[j]] <<- c(flags[[j]], msg)
  add(!is.na(qc$mapping_rate_pct) & qc$mapping_rate_pct < min_map,
      sprintf("mapping rate < %g%%", min_map))
  add(mad_low(qc$mapping_rate_pct, k), sprintf("mapping rate > %g MAD below median", k))
  add(mad_low(log10(pmax(qc$assigned_reads, 1)), k),
      sprintf("library size > %g MAD below median (log scale)", k))
  add(mad_low(qc$detected_genes, k), sprintf("detected genes > %g MAD below median", k))
  add(mad_low(qc$median_cor, k), sprintf("median correlation > %g MAD below median", k))
  qc$flags <- vapply(flags, function(f) if (is.null(f)) "" else paste(f, collapse = "; "), "")
  qc$flagged <- nzchar(qc$flags)
  list(qc = qc, cor = cm)
}

# --- Figures ------------------------------------------------------------------

plot_qc_overview <- function(qc, min_map, out) {
  open_png(out, "qc_overview.png", width = 9, height = 7)
  op <- graphics::par(mfrow = c(2, 2), mar = c(7, 4.5, 3, 1))
  on.exit({graphics::par(op); grDevices::dev.off()})
  col <- ifelse(qc$flagged, "#C0392B", "#7F8C8D")
  graphics::barplot(qc$assigned_reads / 1e6, names.arg = qc$sample, las = 2, col = col,
                    ylab = "Reads assigned to genes (millions)", main = "Library size")
  if (all(is.na(qc$mapping_rate_pct))) {
    graphics::plot.new(); graphics::title("Mapping rate: not available")
  } else {
    graphics::barplot(qc$mapping_rate_pct, names.arg = qc$sample, las = 2, col = col,
                      ylim = c(0, 100), ylab = "Mapping rate (%)", main = "Mapping rate")
    graphics::abline(h = min_map, lty = 2)
  }
  graphics::barplot(qc$detected_genes, names.arg = qc$sample, las = 2, col = col,
                    ylab = "Detected genes", main = "Detected genes")
  rng <- range(qc$median_cor)
  graphics::plot(seq_along(qc$sample), qc$median_cor, pch = 19, col = col, xaxt = "n",
                 xlab = "", ylab = "Median correlation with other samples",
                 main = "Sample agreement", ylim = rng + c(-0.05, 0.05) * diff(rng + c(0, 1e-6)))
  graphics::axis(1, at = seq_along(qc$sample), labels = qc$sample, las = 2)
}

plot_cor_heatmap <- function(cm, out) {
  hc <- stats::hclust(stats::as.dist(1 - cm), method = "average")
  o <- hc$order
  n <- ncol(cm)
  open_png(out, "sample_correlation.png", width = 7, height = 6.5)
  op <- graphics::par(mar = c(7, 7, 3, 5))
  on.exit({graphics::par(op); grDevices::dev.off()})
  pal <- grDevices::hcl.colors(50, "YlOrRd", rev = TRUE)
  graphics::image(seq_len(n), seq_len(n), cm[o, o], col = pal, axes = FALSE,
                  xlab = "", ylab = "", main = "Sample correlation (transformed counts)")
  graphics::axis(1, at = seq_len(n), labels = colnames(cm)[o], las = 2, cex.axis = 0.8)
  graphics::axis(2, at = seq_len(n), labels = colnames(cm)[o], las = 2, cex.axis = 0.8)
  if (n <= 30) graphics::text(rep(seq_len(n), n), rep(seq_len(n), each = n),
                              sprintf("%.2f", cm[o, o]), cex = 0.55)
  graphics::mtext(sprintf("range %.3f to %.3f", min(cm), max(cm[upper.tri(cm)])), 3, 0.2,
                  cex = 0.8)
}

# Legend position with the fewest points nearby.
emptiest_corner <- function(x, y) {
  rx <- (x - min(x)) / max(diff(range(x)), 1e-12)
  ry <- (y - min(y)) / max(diff(range(y)), 1e-12)
  n <- c(topright = sum(rx > 0.5 & ry > 0.5), topleft = sum(rx <= 0.5 & ry > 0.5),
         bottomright = sum(rx > 0.5 & ry <= 0.5), bottomleft = sum(rx <= 0.5 & ry <= 0.5))
  names(n)[which.min(n)]
}

plot_pca <- function(pca, sheet, color_by, qc, out) {
  sc <- pca$scores
  grp <- if (nzchar(color_by)) as.factor(sheet[[color_by]][match(rownames(sc), sheet$sample)])
         else factor(rep("all", nrow(sc)))
  pal <- grDevices::hcl.colors(max(3, nlevels(grp)), "Dark 3")[seq_len(nlevels(grp))]
  open_png(out, "pca.png", width = 7, height = 5.5)
  op <- graphics::par(mar = c(4.5, 4.5, 3, 1))
  on.exit({graphics::par(op); grDevices::dev.off()})
  pc2 <- if (ncol(sc) >= 2) sc[, 2] else rep(0, nrow(sc))
  graphics::plot(sc[, 1], pc2, pch = ifelse(qc$flagged[match(rownames(sc), qc$sample)], 17, 19),
                 col = pal[as.integer(grp)], cex = 1.4,
                 xlab = sprintf("PC1 (%.1f%%)", 100 * pca$var_explained[1]),
                 ylab = sprintf("PC2 (%.1f%%)", 100 * pca$var_explained[min(2, length(pca$var_explained))]),
                 main = sprintf("PCA of %d most variable genes", pca$n_genes))
  graphics::text(sc[, 1], pc2, rownames(sc), pos = 3, cex = 0.7)
  graphics::legend(emptiest_corner(sc[, 1], pc2), legend = c(levels(grp), "flagged sample"),
                   col = c(pal, "black"), pch = c(rep(19, nlevels(grp)), 17), bty = "n",
                   cex = 0.8, title = if (nzchar(color_by)) color_by else NULL)
}

# --- Main -------------------------------------------------------------------

design_columns <- function(sheet) setdiff(names(sheet), SHEET_SYSTEM_COLS)

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  if (opt$demo) {
    opt$samples <- file.path(SKILL_DIR, "demo", "samples.csv")
    opt$gff <- file.path(SKILL_DIR, "demo", "annotation.gff3")
  }
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  opt$type <- tolower(opt$type)
  if (!opt$type %in% c("auto", "salmon", "kallisto", "featurecounts"))
    fail("--type must be auto, salmon, kallisto or featurecounts.")
  if (!opt$transform %in% c("vst", "rlog")) fail("--transform must be vst or rlog.")
  if (!opt$id_attr %in% c("auto", "ID", "Name", "transcript_id", "ID_clean"))
    fail("--id-attr must be auto, ID, Name or transcript_id.")
  if (opt$min_id_match <= 0 || opt$min_id_match > 1) fail("--min-id-match must be in (0, 1].")
  sheet <- read_sample_sheet(opt$samples)
  base_dir <- dirname(normalizePath(opt$samples, winslash = "/"))
  if (opt$type != "featurecounts" && nzchar(opt$counts))
    fail("--counts is for featureCounts matrices; add --type featurecounts.")
  if (opt$type == "featurecounts" && !nzchar(opt$counts))
    fail("--type featurecounts needs --counts FEATURECOUNTS.txt.")

  warn <- character()
  notes <- character()
  inputs <- c(opt$samples, opt$gff, opt$tx2gene, opt$counts, opt$fc_summary)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  lengths <- NULL
  abundance <- NULL
  id_info <- NULL

  if (opt$type == "featurecounts") {
    fc <- read_featurecounts(opt$counts)
    colmap <- match_fc_columns(sheet, colnames(fc$counts))
    counts <- fc$counts[, colmap, drop = FALSE]
    colnames(counts) <- names(colmap)
    mapping <- read_fc_summary(opt$fc_summary, colmap)
    if (is.null(mapping)) {
      mapping <- data.frame(sample = sheet$sample, n_processed = NA_real_, n_mapped = NA_real_,
                            mapping_rate_pct = NA_real_)
      warn <- c(warn, "No --fc-summary given: the assignment rate per sample is not checked.")
    } else {
      notes <- c(notes, "For featureCounts, 'mapping_rate_pct' is the % of read pairs/reads assigned to genes (from the .summary file).")
    }
    quant_type <- "featureCounts"
    if (nzchar(opt$gff) || nzchar(opt$tx2gene))
      notes <- c(notes, "featureCounts output is already gene-level; --gff/--tx2gene were not used.")
  } else {
    rq <- resolve_quant_files(sheet, base_dir, opt$type)
    files <- rq$files
    quant_type <- rq$type
    inputs <- c(inputs, files)
    if (!nzchar(opt$gff) && !nzchar(opt$tx2gene))
      fail("Give --gff ANNOTATION (GFF3/GTF) or --tx2gene TABLE to map transcripts to genes. ",
           "Use the same annotation that was used to build the ", quant_type, " index.")
    if (nzchar(opt$gff) && nzchar(opt$tx2gene)) fail("Give either --gff or --tx2gene, not both.")
    txi_tx <- import_transcripts(files, quant_type)
    q_ids <- rownames(txi_tx$counts)
    if (nzchar(opt$gff)) {
      ann <- parse_annotation(opt$gff)
      ch <- choose_tx2gene(ann$candidates, if (opt$strip_tx_version)
        sub("\\.[0-9]+$", "", q_ids) else q_ids, opt$id_attr)
      if (opt$strip_tx_version && opt$id_attr == "auto") {
        # Choose the attribute on stripped IDs, as matching is done after stripping.
        fr <- vapply(ann$candidates, function(d)
          mean(sub("\\.[0-9]+$", "", q_ids) %in% sub("\\.[0-9]+$", "", d$transcript_id)), 0)
        ch$style <- names(ann$candidates)[which.max(fr)]
        ch$tx2gene <- ann$candidates[[ch$style]]
        ch$match <- fr
      }
      tx2gene <- ch$tx2gene
      id_info <- list(source = sprintf("%s (%s)", basename(opt$gff), ann$format),
                      style = ch$style, match = ch$match,
                      n_multi_parent = ann$n_multi_parent %||% 0)
      if (isTRUE(id_info$n_multi_parent > 0))
        warn <- c(warn, sprintf("%d transcript(s) in the GFF3 list several Parent genes; the first parent was used.",
                                id_info$n_multi_parent))
    } else {
      tx2gene <- read_tx2gene_table(opt$tx2gene)
      id_info <- list(source = basename(opt$tx2gene), style = "user table",
                      match = c(`user table` = mean(q_ids %in% tx2gene$transcript_id)))
    }
    sg <- summarise_to_genes(txi_tx, tx2gene, opt$strip_tx_version, opt$min_id_match)
    txi <- sg$txi
    counts <- txi$counts
    lengths <- txi$length
    abundance <- txi$abundance
    id_info$match_frac <- sg$match_frac
    id_info$n_unmatched <- length(sg$unmatched)
    id_info$n_quantified <- length(q_ids)
    if (length(sg$unmatched)) {
      save_table(data.frame(transcript_id = sg$unmatched), opt$out, "unmatched_transcripts.csv")
      warn <- c(warn, sprintf(paste0(
        "%d quantified transcript(s) (%.2f%%) are not in tx2gene and were dropped; they held ",
        "up to %.2f%% of a sample's reads. See tables/unmatched_transcripts.csv."),
        length(sg$unmatched), 100 * (1 - sg$match_frac), 100 * max(sg$reads_lost_frac)))
    }
    if (opt$strip_tx_version)
      notes <- c(notes, "Trailing '.N' version suffixes were removed from transcript IDs (--strip-tx-version).")
    mapping <- do.call(rbind, lapply(names(files), function(s) {
      m <- read_mapping_info(files[[s]], quant_type)
      data.frame(sample = s, n_processed = m$n_processed, n_mapped = m$n_mapped,
                 mapping_rate_pct = m$mapping_rate_pct, library_type = m$library_type,
                 stringsAsFactors = FALSE)
    }))
    miss <- mapping$sample[is.na(mapping$mapping_rate_pct)]
    if (length(miss))
      warn <- c(warn, sprintf("Mapping rate unavailable for %d sample(s) (%s): %s not found next to the quant file.",
                              length(miss), paste(utils::head(miss, 5), collapse = ", "),
                              if (quant_type == "salmon") "aux_info/meta_info.json" else "run_info.json"))
    lt <- unique(stats::na.omit(mapping$library_type))
    if (length(lt) > 1L)
      warn <- c(warn, paste0("Salmon inferred different library types across samples (",
                             paste(lt, collapse = ", "), "): check that all libraries were prepared the same way."))
    saveRDS(txi, file.path(opt$out, "txi.rds"))
    save_table(data.frame(gene_id = rownames(lengths), lengths, check.names = FALSE),
               opt$out, "gene_lengths.csv")
    save_table(data.frame(gene_id = rownames(abundance), abundance, check.names = FALSE),
               opt$out, "gene_tpm.csv")
    save_table(tx2gene, opt$out, "tx2gene.csv")
  }
  if (ncol(counts) < 2L) fail("Need at least 2 samples for QC.")

  # QC
  tmat <- qc_transform(counts, if (opt$three_prime) NULL else lengths, opt$transform)
  if (opt$three_prime && quant_type != "featureCounts")
    notes <- c(notes, "3' tag-seq (--three-prime): transcript-length offsets are not used for QC and must not be used for DE.")
  sq <- sample_qc(counts, mapping, tmat, opt$detect_min_count, opt$min_mapping_rate,
                  opt$outlier_mad)
  qc <- sq$qc
  pca <- qc_pca(tmat)
  color_by <- opt$color_by
  dcols <- design_columns(sheet)
  if (!nzchar(color_by) && length(dcols)) color_by <- dcols[[1L]]
  if (nzchar(color_by) && !color_by %in% names(sheet))
    fail("--color-by column '", color_by, "' is not in the sample sheet.")

  qc_out <- qc
  qc_out$median_cor <- round(qc_out$median_cor, 4)
  qc_out$PC1 <- round(pca$scores[qc$sample, 1], 3)
  if (ncol(pca$scores) >= 2) qc_out$PC2 <- round(pca$scores[qc$sample, 2], 3)
  save_table(qc_out, opt$out, "sample_qc.csv")
  save_table(data.frame(gene_id = rownames(counts), counts, check.names = FALSE),
             opt$out, "gene_counts.csv")
  save_table(data.frame(sample = rownames(sq$cor), round(sq$cor, 5), check.names = FALSE),
             opt$out, "sample_correlation.csv")
  samples_out <- sheet
  samples_out$qc_flagged <- qc$flagged[match(sheet$sample, qc$sample)]
  samples_out$qc_flags <- qc$flags[match(sheet$sample, qc$sample)]
  save_table(samples_out, opt$out, "samples.csv")

  plot_qc_overview(qc, opt$min_mapping_rate, opt$out)
  plot_cor_heatmap(sq$cor, opt$out)
  plot_pca(pca, sheet, color_by, qc, opt$out)

  # Warnings
  if (any(qc$flagged))
    warn <- c(warn, sprintf("Sample %s flagged: %s.", qc$sample[qc$flagged], qc$flags[qc$flagged]))
  if (nrow(qc) < 4)
    warn <- c(warn, "Fewer than 4 samples: the MAD-based outlier rules have little power.")

  # Report
  fmt_reads <- function(x) ifelse(is.na(x), "NA", sprintf("%.2f", x / 1e6))
  rep <- c(sprintf("# %s report", SKILL), "",
           if (length(warn)) c("## Warnings", "", paste("-", warn), "") else NULL,
           "## Summary", "",
           sprintf("- Input: %d samples, %s quantification (`%s`).", nrow(qc), quant_type,
                   basename(opt$samples)),
           sprintf("- Genes: %d (%d detected with >= %g reads in at least one sample).",
                   nrow(counts), sum(apply(counts >= opt$detect_min_count, 1, any)),
                   opt$detect_min_count))
  if (!is.null(id_info)) {
    mtxt <- paste(sprintf("%s %.0f%%", names(id_info$match), 100 * id_info$match), collapse = ", ")
    rep <- c(rep,
             sprintf("- Transcript-to-gene map: %s; transcript IDs taken from `%s`.",
                     id_info$source, id_info$style),
             sprintf("- ID match of quantified transcripts by attribute: %s.", mtxt),
             sprintf("- Transcripts used: %.2f%% of %d (unmatched: %d).",
                     100 * id_info$match_frac, id_info$n_quantified, id_info$n_unmatched))
  }
  if (length(notes)) rep <- c(rep, paste("-", notes))
  rep <- c(rep, "", "## Sample QC", "",
           sprintf("Flag rules: mapping rate < %g%%; or more than %g MADs below the median for mapping rate, log10 library size, detected genes (>= %g reads) or median between-sample correlation. Flags are prompts to inspect the sample, not automatic removals.",
                   opt$min_mapping_rate, opt$outlier_mad, opt$detect_min_count), "",
           "| sample | mapping rate (%) | assigned reads (M) | detected genes | median cor | flags |",
           "|---|---|---|---|---|---|",
           sprintf("| %s | %s | %s | %d | %.3f | %s |", qc$sample,
                   ifelse(is.na(qc$mapping_rate_pct), "NA", sprintf("%.1f", qc$mapping_rate_pct)),
                   fmt_reads(qc$assigned_reads), qc$detected_genes, qc$median_cor,
                   ifelse(qc$flagged, qc$flags, "")),
           "",
           sprintf("PCA on %s-transformed counts (blind to the design), %d most variable genes: PC1 %.1f%%, PC2 %.1f%% of variance.",
                   opt$transform, pca$n_genes, 100 * pca$var_explained[1],
                   100 * pca$var_explained[min(2, length(pca$var_explained))]),
           "", "![QC overview](figures/qc_overview.png)", "",
           "![PCA](figures/pca.png)", "", "![Sample correlation](figures/sample_correlation.png)", "",
           "## Next step", "")
  if (quant_type == "featureCounts") {
    rep <- c(rep, "Run claw-de-analysis with `--counts tables/gene_counts.csv --samples tables/samples.csv`.")
  } else {
    rep <- c(rep, if (opt$three_prime)
      "Run claw-de-analysis with `--txi txi.rds --samples tables/samples.csv --no-length-offset`: 3' tag-seq counts have no length bias, so the estimated counts are used without length offsets (tximport vignette). Never use `gene_tpm.csv` for differential expression."
      else "Run claw-de-analysis with `--txi txi.rds --samples tables/samples.csv`. It uses the estimated counts plus the average-transcript-length offsets (DESeqDataSetFromTximport). Never use `gene_tpm.csv` for differential expression.")
  }
  rep <- c(rep, "", "## Methods", "",
           paste0(if (quant_type == "featureCounts") "Gene-level counts were read from featureCounts output. " else
             sprintf("Transcript-level estimates from %s were imported and summarised to genes with tximport (countsFromAbundance = \"no\"; gene length = abundance-weighted mean of transcript effective lengths). ", quant_type),
             sprintf("For QC only, counts were transformed with DESeq2 %s (blind = TRUE); PCA used the %d most variable genes (as DESeq2::plotPCA).",
                     opt$transform, pca$n_genes)),
           "", "## References", "",
           "- Soneson C, Love MI, Robinson MD (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research 4:1521. https://doi.org/10.12688/f1000research.7563.1",
           "- Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods 14:417-419. https://doi.org/10.1038/nmeth.4197",
           "- Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15:550. https://doi.org/10.1186/s13059-014-0550-8",
           "", sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(rep, file.path(opt$out, "report.md"))
  write_repro_bundle(opt$out, argv, inputs = inputs)
  cat(sprintf("Done: %d samples, %d genes, %d flagged. Report: %s\n", nrow(qc), nrow(counts),
              sum(qc$flagged), normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(list(counts = counts, qc = qc, txi = if (quant_type == "featureCounts") NULL else txi))
}

if (sys.nframe() == 0L) main()
