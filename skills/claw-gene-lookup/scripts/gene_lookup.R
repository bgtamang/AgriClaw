#!/usr/bin/env Rscript
# claw-gene-lookup: list the genes inside GWAS or QTL intervals, with gene IDs,
# coordinates, descriptions and the distance to the peak SNP. Genes come from
# the Ensembl REST API (rest.ensembl.org, which serves Ensembl Plants genomes)
# or from a local GFF3 file for any assembly. The genome assembly of the user's
# coordinates must be stated, and the skill refuses to mix assemblies.
#
# Usage:
#   Rscript scripts/gene_lookup.R --demo
#   Rscript scripts/gene_lookup.R --intervals qtl_regions.csv --ensembl-species glycine_max \
#           --assembly Glycine_max_v2.1
#   Rscript scripts/gene_lookup.R --snps hits.csv --window-bp 250000 --gff genes.gff3 \
#           --assembly Wm82.a4.v1
#   Rscript scripts/gene_lookup.R --help
#
# Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.

SKILL   <- "claw-gene-lookup"
VERSION <- "1.0.0"

USAGE <- "Usage: Rscript scripts/gene_lookup.R [--demo] (--intervals FILE | --snps FILE)
                                  (--gff FILE | --ensembl-species NAME) --assembly NAME [options]

Regions (one of)
  --intervals FILE     CSV with chr, start, end (1-based, inclusive), optional name/region_id,
                       lead_snp, peak_pos or lead_pos_bp (qtl_regions.csv from claw-gwas-crop
                       works as is); or a .bed file (0-based start, converted)
  --snps FILE          CSV with snp, chr, pos (peak SNPs); windows are pos +/- --window-bp
  --window-bp N        half-width of SNP windows in bp
  --ld-table FILE      per-chromosome half-widths (chr, ld_distance_bp) from claw-ld-decay
  --coords C           1-based (default) or bed: coordinate convention of a CSV --intervals file

Genes (one of)
  --gff FILE           GFF3 annotation (.gff3/.gff, .gz ok) for the same assembly
  --ensembl-species S  Ensembl species name, e.g. glycine_max, zea_mays, triticum_aestivum,
                       oryza_sativa, sorghum_bicolor, arabidopsis_thaliana (needs internet)

Required
  --assembly NAME      assembly of YOUR coordinates (e.g. Glycine_max_v2.1, Wm82.a4.v1,
                       Zm-B73-REFERENCE-NAM-5.0). Checked against the gene source.

Options
  --feature-type T     GFF3 feature type(s) to report, comma-separated (default: gene)
  --biotype B          Ensembl only: keep one gene biotype, e.g. protein_coding
  --chr-map FILE       CSV with user_chr, target_chr to override chromosome-name matching
  --ensembl-server URL default https://rest.ensembl.org
  --out DIR            output folder (default: agriclaw_output/claw-gene-lookup)
"

# Option defaults. The type of each default decides how its value is parsed:
# logical = flag, numeric = number, character = text.
DEFAULTS <- list(demo = FALSE, intervals = "", snps = "", window_bp = NA_real_, ld_table = "",
                 coords = "1-based", gff = "", ensembl_species = "", assembly = "",
                 feature_type = "gene", biotype = "", chr_map = "",
                 ensembl_server = "https://rest.ensembl.org", out = "")

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

COL_BLUE <- "#2a78d6"; COL_ORANGE <- "#eb6834"; COL_INK <- "#0b0b0b"
COL_INK2 <- "#52514e"; COL_MUTED <- "#898781"

ENSEMBL_MAX_REGION <- 5e6   # overlap/region endpoint limit (5 Mb per request)

norm_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))

pick_col <- function(nm, aliases) {
  i <- match(TRUE, tolower(nm) %in% aliases)
  if (is.na(i)) NA_character_ else nm[i]
}

read_csv_checked <- function(path, what) {
  if (!nzchar(path) || !file.exists(path)) stop(what, " file not found: ", path, call. = FALSE)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
}

as_bp <- function(x, what) {
  v <- suppressWarnings(as.numeric(x))
  if (anyNA(v) || any(v < 0) || any(v != round(v)))
    stop(what, " must be whole, non-negative base-pair positions; check the file.", call. = FALSE)
  v
}

# --- Regions --------------------------------------------------------------------

# Intervals from CSV (1-based inclusive unless coords = "bed") or BED (0-based,
# half-open: start + 1, end unchanged).
read_intervals <- function(path, coords = "1-based") {
  if (!coords %in% c("1-based", "bed")) stop("--coords must be 1-based or bed.", call. = FALSE)
  is_bed <- grepl("\\.bed(\\.gz)?$", tolower(path))
  if (is_bed) {
    if (!file.exists(path)) stop("Intervals file not found: ", path, call. = FALSE)
    con <- gzfile(path, "rt")
    l <- readLines(con, warn = FALSE)
    close(con)
    l <- l[nzchar(l) & !grepl("^(#|track|browser)", l)]
    f <- strsplit(l, "\t", fixed = TRUE)
    if (any(lengths(f) < 3L)) stop("BED lines need at least chrom, chromStart, chromEnd.", call. = FALSE)
    iv <- data.frame(interval_id = vapply(seq_along(f), function(k) if (length(f[[k]]) >= 4L) f[[k]][4] else sprintf("I%02d", k), ""),
                     chr = vapply(f, `[`, "", 1L), start = as_bp(vapply(f, `[`, "", 2L), "BED start") + 1,
                     end = as_bp(vapply(f, `[`, "", 3L), "BED end"), stringsAsFactors = FALSE)
    iv$peak_pos <- NA_real_; iv$lead_snp <- NA_character_; iv$trait <- NA_character_
    attr(iv, "note") <- "BED input: 0-based starts were converted to 1-based (start + 1)."
  } else {
    df <- read_csv_checked(path, "Intervals")
    nm <- names(df)
    cc <- pick_col(nm, c("chr", "chrom", "chromosome", "seqid", "seqname"))
    cs <- pick_col(nm, c("start", "start_bp", "from", "begin"))
    ce <- pick_col(nm, c("end", "end_bp", "to", "stop"))
    if (anyNA(c(cc, cs, ce))) stop("The intervals file needs chr, start and end columns (start_bp/end_bp also work). Found: ",
                                   paste(nm, collapse = ", "), call. = FALSE)
    cid <- pick_col(nm, c("region_id", "interval_id", "name", "id"))
    cp <- pick_col(nm, c("peak_pos", "peak_pos_bp", "lead_pos_bp", "lead_pos", "peak"))
    cl <- pick_col(nm, c("lead_snp", "snp", "marker"))
    ct <- pick_col(nm, c("trait"))
    iv <- data.frame(interval_id = if (!is.na(cid)) df[[cid]] else sprintf("I%02d", seq_len(nrow(df))),
                     chr = df[[cc]], start = as_bp(df[[cs]], "Interval start"), end = as_bp(df[[ce]], "Interval end"),
                     peak_pos = if (!is.na(cp)) suppressWarnings(as.numeric(df[[cp]])) else NA_real_,
                     lead_snp = if (!is.na(cl)) df[[cl]] else NA_character_,
                     trait = if (!is.na(ct)) df[[ct]] else NA_character_, stringsAsFactors = FALSE)
    if (!is.na(ct) && anyDuplicated(iv$interval_id)) iv$interval_id <- paste(iv$trait, iv$interval_id, sep = ":")
    if (coords == "bed") {
      iv$start <- iv$start + 1
      attr(iv, "note") <- "--coords bed: 0-based starts were converted to 1-based (start + 1)."
    } else attr(iv, "note") <- "Intervals read as 1-based, inclusive coordinates (GFF3/VCF convention)."
  }
  if (!nrow(iv)) stop("The intervals file has no rows.", call. = FALSE)
  if (any(iv$start < 1)) stop("Interval starts must be >= 1 in 1-based coordinates. Was this a 0-based BED file? Use --coords bed.", call. = FALSE)
  if (any(iv$end < iv$start)) stop("Some intervals end before they start: ",
                                   paste(utils::head(iv$interval_id[iv$end < iv$start], 3), collapse = ", "), call. = FALSE)
  iv
}

# Peak SNPs to windows of +/- half-width (scalar or per chromosome).
snp_windows <- function(path, window_bp = NA_real_, ld_table = "") {
  df <- read_csv_checked(path, "SNP")
  nm <- names(df)
  cs <- pick_col(nm, c("snp", "marker", "id", "lead_snp", "rs#"))
  cc <- pick_col(nm, c("chr", "chrom", "chromosome"))
  cp <- pick_col(nm, c("pos", "position", "bp", "pos_bp", "lead_pos_bp"))
  if (anyNA(c(cc, cp))) stop("The SNP file needs chr and pos columns. Found: ", paste(nm, collapse = ", "), call. = FALSE)
  pos <- as_bp(df[[cp]], "SNP position")
  chr <- df[[cc]]
  if (nzchar(ld_table)) {
    lt <- read_csv_checked(ld_table, "LD table")
    if (!all(c("chr", "ld_distance_bp") %in% names(lt)))
      stop("--ld-table needs columns chr and ld_distance_bp (tables/ld_window.csv from claw-ld-decay).", call. = FALSE)
    dflt <- suppressWarnings(as.numeric(lt$ld_distance_bp[tolower(lt$chr) %in% c("all", "genome", "genome-wide")]))
    w <- suppressWarnings(as.numeric(lt$ld_distance_bp[match(chr, lt$chr)]))
    if (length(dflt)) w[is.na(w)] <- dflt[1]   # NA rows (not estimable) fall back to the 'all' row
    if (anyNA(w)) stop("--ld-table has no usable distance for chromosome(s) ", paste(unique(chr[is.na(w)]), collapse = ", "),
                       " (missing or NA, and no usable 'all' row; see the claw-ld-decay report). ",
                       "Ask the user for a window and pass --window-bp instead.", call. = FALSE)
    basis <- sprintf("per-chromosome LD distances from %s", ld_table)
  } else {
    if (is.na(window_bp) || window_bp < 0)
      stop("SNP mode needs a window half-width: --window-bp N, or --ld-table from claw-ld-decay. ",
           "Ask the user or run claw-ld-decay; do not guess a window.", call. = FALSE)
    w <- rep(window_bp, length(pos))
    basis <- sprintf("+/- %.0f bp given with --window-bp", window_bp)
  }
  iv <- data.frame(interval_id = if (!is.na(cs)) df[[cs]] else sprintf("SNP%02d", seq_along(pos)), chr = chr,
                   start = pmax(1, pos - w), end = pos + w, peak_pos = pos,
                   lead_snp = if (!is.na(cs)) df[[cs]] else NA_character_, trait = NA_character_,
                   stringsAsFactors = FALSE)
  attr(iv, "note") <- paste("Windows: peak SNP", basis, "(1-based, inclusive).")
  iv
}

# --- Chromosome names -------------------------------------------------------------

# Canonical key: lower case, without chromosome/species prefixes and leading
# zeros, so Gm01, Chr01, chr1 and 1 all become "1", and 1A stays "1a".
chrom_key <- function(x) {
  k <- tolower(trimws(as.character(x)))
  k <- sub("^(chromosome|chrom|chr|ch)[_ .-]?(?=[0-9])", "", k, perl = TRUE)
  k <- sub("^(gm|zm|os|sb|at|ta|sl)[_ .-]?(?=[0-9])", "", k, perl = TRUE)
  sub("^0+(?=[0-9])", "", k, perl = TRUE)
}

map_chromosomes <- function(user_chr, target_chr, chr_map = NULL) {
  u <- unique(as.character(user_chr))
  target_chr <- unique(as.character(target_chr))
  out <- data.frame(user_chr = u, target_chr = NA_character_, method = NA_character_, stringsAsFactors = FALSE)
  if (!is.null(chr_map)) {
    i <- match(u, chr_map$user_chr)
    out$target_chr <- chr_map$target_chr[i]
    out$method[!is.na(i)] <- "--chr-map"
  }
  exact <- is.na(out$target_chr) & u %in% target_chr
  out$target_chr[exact] <- u[exact]; out$method[exact] <- "exact"
  tk <- chrom_key(target_chr)
  amb <- character(); miss <- character()
  for (k in which(is.na(out$target_chr))) {
    hit <- target_chr[tk == chrom_key(u[k])]
    if (length(hit) == 1L) { out$target_chr[k] <- hit; out$method[k] <- "name normalised" }
    else if (length(hit) > 1L) amb <- c(amb, sprintf("%s -> %s", u[k], paste(hit, collapse = " or ")))
    else miss <- c(miss, u[k])
  }
  if (length(amb)) stop("Ambiguous chromosome names: ", paste(amb, collapse = "; "),
                        ". Give an explicit mapping with --chr-map (columns user_chr, target_chr).", call. = FALSE)
  if (length(miss)) stop("Chromosome(s) not found in the gene source: ", paste(miss, collapse = ", "),
                         ". The source uses names like ", paste(utils::head(target_chr, 5), collapse = ", "),
                         ". Check the assembly, or give --chr-map (columns user_chr, target_chr).", call. = FALSE)
  if (!is.null(chr_map)) {
    bad <- setdiff(out$target_chr, target_chr)
    if (length(bad)) stop("--chr-map points to chromosome(s) absent from the gene source: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  out
}

# --- Assembly checks ---------------------------------------------------------------

# Maize gene-model prefixes tie IDs to an assembly: Zm00001d = B73 RefGen_v4,
# Zm00001eb = Zm-B73-REFERENCE-NAM-5.0 (MaizeGDB; Ensembl archive).
maize_id_check <- function(gene_ids, assembly) {
  if (!length(gene_ids)) return(NA_character_)
  f_v5 <- mean(grepl("^Zm00001eb", gene_ids)); f_v4 <- mean(grepl("^Zm00001d", gene_ids))
  a <- norm_name(assembly)
  says_v4 <- grepl("v4|refgenv4|gramene40|agpv4", a)
  says_v5 <- grepl("v5|nam50|refgenv5", a)
  if (f_v5 > 0.5 && says_v4)
    stop("The gene IDs are maize B73 v5 models (Zm00001eb...), but --assembly says v4 ('", assembly,
         "'). v4 coordinates cannot be used with v5 genes: lift the coordinates over or use a v4 annotation.", call. = FALSE)
  if (f_v4 > 0.5 && says_v5)
    stop("The gene IDs are maize B73 RefGen_v4 models (Zm00001d...), but --assembly says v5 ('", assembly,
         "'). Use the annotation that matches the coordinates.", call. = FALSE)
  if (f_v5 > 0.5) return("gene IDs are maize B73 v5 models (Zm00001eb)")
  if (f_v4 > 0.5) return("gene IDs are maize B73 RefGen_v4 models (Zm00001d)")
  NA_character_
}

check_ensembl_assembly <- function(assembly, info) {
  accepted <- c(info$assembly_name, info$assembly_accession, info$default_coord_system_version)
  accepted <- accepted[!is.na(accepted) & nzchar(accepted)]
  if (!norm_name(assembly) %in% norm_name(accepted))
    stop("Assembly mismatch: your coordinates are on '", assembly, "', but Ensembl serves '", info$assembly_name,
         "' (", info$assembly_accession, ") for this species. Coordinates from different assemblies cannot be ",
         "mixed. If your positions really are on ", info$assembly_name, ", re-run with --assembly ",
         info$assembly_name, "; otherwise use --gff with the annotation of your assembly, or lift the ",
         "coordinates over first.", call. = FALSE)
  invisible(TRUE)
}

check_gff_assembly <- function(assembly, pragmas) {
  keys <- c("genome-build", "genome-version", "assembly", "annot-version", "genome-build-accession")
  vals <- pragmas[names(pragmas) %in% keys]
  if (!length(vals)) return(list(ok = NA, stated = NA_character_))
  if (norm_name(assembly) %in% norm_name(vals)) return(list(ok = TRUE, stated = paste(vals, collapse = "; ")))
  stop("Assembly mismatch: your coordinates are on '", assembly, "', but the GFF3 header says ",
       paste(sprintf("%s '%s'", names(vals), vals), collapse = " and "), ". Use the annotation for your ",
       "assembly, or confirm with the user and re-run with that exact name in --assembly.", call. = FALSE)
}

# --- GFF3 -----------------------------------------------------------------------------

url_decode <- function(x) {
  has <- !is.na(x) & grepl("%[0-9A-Fa-f]{2}", x)
  x[has] <- vapply(x[has], utils::URLdecode, "")
  x
}

parse_attributes <- function(attr_str) {
  parts <- strsplit(attr_str, ";", fixed = TRUE)
  get <- function(keys) vapply(parts, function(p) {
    kv <- strsplit(p, "=", fixed = TRUE)
    k <- vapply(kv, `[`, "", 1L)
    v <- vapply(kv, function(z) if (length(z) > 1L) paste(z[-1], collapse = "=") else NA_character_, "")
    for (key in keys) { i <- match(key, trimws(k)); if (!is.na(i)) return(v[i]) }
    NA_character_
  }, "")
  list(id = url_decode(get("ID")), name = url_decode(get(c("Name", "gene_name"))),
       description = url_decode(get(c("description", "Note", "product", "function"))),
       biotype = url_decode(get(c("biotype", "gene_biotype"))))
}

read_gff3 <- function(path, feature_types = "gene") {
  if (!file.exists(path)) stop("GFF3 file not found: ", path, call. = FALSE)
  con <- gzfile(path, "rt")
  lines <- readLines(con, warn = FALSE)
  close(con)
  fa <- match("##FASTA", lines)
  if (!is.na(fa)) lines <- lines[seq_len(fa - 1L)]
  hdr <- lines[startsWith(lines, "#")]
  pr <- regmatches(hdr, regexec("^#[#!]\\s*([A-Za-z-]+)\\s+(.+?)\\s*$", hdr))
  pr <- pr[lengths(pr) == 3L]
  pragmas <- stats::setNames(vapply(pr, `[`, "", 3L), tolower(vapply(pr, `[`, "", 2L)))
  sr <- regmatches(hdr, regexec("^##sequence-region\\s+(\\S+)\\s+([0-9]+)\\s+([0-9]+)", hdr))
  sr <- sr[lengths(sr) == 4L]
  seq_len_tab <- if (length(sr)) data.frame(name = vapply(sr, `[`, "", 2L), length = as.numeric(vapply(sr, `[`, "", 4L)),
                                            stringsAsFactors = FALSE) else NULL
  body <- lines[nzchar(lines) & !startsWith(lines, "#")]
  pat <- paste0("^[^\t]*\t[^\t]*\t(", paste(gsub("([.+])", "\\\\\\1", feature_types), collapse = "|"), ")\t")
  body <- body[grepl(pat, body)]
  if (!length(body)) stop("No '", paste(feature_types, collapse = "/"), "' features in the GFF3 file. Check --feature-type ",
                          "(column 3 of the file).", call. = FALSE)
  f <- strsplit(body, "\t", fixed = TRUE)
  if (any(lengths(f) != 9L)) stop("Malformed GFF3: every feature line needs 9 tab-separated columns.", call. = FALSE)
  m <- matrix(unlist(f), ncol = 9L, byrow = TRUE)
  a <- parse_attributes(m[, 9])
  id <- sub("^gene:", "", a$id)
  genes <- data.frame(gene_id = ifelse(is.na(id), paste0(m[, 1], ":", m[, 4], "-", m[, 5]), id),
                      gene_name = a$name, chr = m[, 1], start = as.numeric(m[, 4]), end = as.numeric(m[, 5]),
                      strand = m[, 7], biotype = ifelse(is.na(a$biotype), m[, 3], a$biotype),
                      description = a$description, stringsAsFactors = FALSE)
  if (anyNA(genes$start) || anyNA(genes$end) || any(genes$end < genes$start))
    stop("Malformed GFF3: feature start/end must be integers with start <= end.", call. = FALSE)
  list(genes = genes, pragmas = pragmas, seq_lengths = seq_len_tab)
}

# --- Ensembl REST -----------------------------------------------------------------------

http_get_text <- function(url) {
  req <- httr2::request(url)
  req <- httr2::req_headers(req, Accept = "application/json")
  req <- httr2::req_user_agent(req, "AgriClaw claw-gene-lookup (https://github.com/bgtamang/AgriClaw)")
  req <- httr2::req_timeout(req, 60)
  req <- httr2::req_throttle(req, capacity = 10, fill_time_s = 1)    # Ensembl allows 15 requests/s
  req <- httr2::req_retry(req, max_tries = 5,
                          is_transient = function(resp) httr2::resp_status(resp) %in% c(429, 500, 502, 503, 504))
  resp <- tryCatch(httr2::req_perform(req), error = function(e) {
    body <- tryCatch(httr2::resp_body_string(e$resp), error = function(e2) "")
    stop("Ensembl request failed: ", url, "\n", conditionMessage(e),
         if (nzchar(body)) paste0("\nServer said: ", substr(body, 1, 300)) else "", call. = FALSE)
  })
  httr2::resp_body_string(resp)
}

ensembl_get <- function(path, server, fetcher, raw_dir = NULL, raw_name = NULL) {
  url <- paste0(sub("/+$", "", server), path)
  txt <- fetcher(url)
  if (!is.null(raw_dir)) {
    dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(txt, file.path(raw_dir, raw_name), useBytes = TRUE)
  }
  jsonlite::fromJSON(txt, simplifyVector = TRUE)
}

ensembl_assembly_info <- function(species, server, fetcher, raw_dir = NULL) {
  info <- ensembl_get(sprintf("/info/assembly/%s?content-type=application/json", species), server, fetcher,
                      raw_dir, sprintf("info_assembly_%s.json", species))
  if (!is.null(info$error)) stop("Ensembl: ", info$error, call. = FALSE)
  tl <- info$top_level_region
  info$regions <- data.frame(name = as.character(tl$name), length = as.numeric(tl$length),
                             coord_system = as.character(tl$coord_system), stringsAsFactors = FALSE)
  info
}

# Split [start, end] into pieces of at most `size` bp (the endpoint limit).
region_chunks <- function(start, end, size = ENSEMBL_MAX_REGION) {
  s <- seq(start, end, by = size)
  data.frame(start = s, end = pmin(s + size - 1, end))
}

parse_ensembl_genes <- function(x) {
  if (is.null(x) || !length(x) || (is.data.frame(x) && !nrow(x)))
    return(data.frame(gene_id = character(), gene_name = character(), chr = character(), start = numeric(),
                      end = numeric(), strand = character(), biotype = character(), description = character(),
                      assembly_name = character(), stringsAsFactors = FALSE))
  col <- function(k) if (k %in% names(x)) as.character(x[[k]]) else rep(NA_character_, nrow(x))
  data.frame(gene_id = col("id"), gene_name = col("external_name"), chr = col("seq_region_name"),
             start = as.numeric(x$start), end = as.numeric(x$end),
             strand = ifelse(as.numeric(x$strand) > 0, "+", "-"), biotype = col("biotype"),
             description = sub("\\s*\\[Source:.*\\]$", "", col("description")),
             assembly_name = col("assembly_name"), stringsAsFactors = FALSE)
}

ensembl_genes <- function(species, chr, start, end, server, fetcher, biotype = "", raw_dir = NULL) {
  ch <- region_chunks(start, end)
  out <- list()
  for (k in seq_len(nrow(ch))) {
    path <- sprintf("/overlap/region/%s/%s:%.0f-%.0f?feature=gene;content-type=application/json%s",
                    species, chr, ch$start[k], ch$end[k], if (nzchar(biotype)) paste0(";biotype=", biotype) else "")
    x <- ensembl_get(path, server, fetcher, raw_dir,
                     sprintf("overlap_%s_%s_%.0f-%.0f.json", species, chr, ch$start[k], ch$end[k]))
    out[[k]] <- parse_ensembl_genes(x)
  }
  g <- do.call(rbind, out)
  g[!duplicated(g$gene_id), , drop = FALSE]
}

# --- Overlap and distances ------------------------------------------------------------------

# Genes overlapping [start, end] (both 1-based, inclusive): gene_start <= end
# and gene_end >= start. Distance to the peak is 0 inside the gene, otherwise
# the gap to the nearest gene end.
annotate_overlaps <- function(iv_row, genes) {
  g <- genes[genes$chr == iv_row$target_chr & genes$start <= iv_row$end & genes$end >= iv_row$start, , drop = FALSE]
  if (!nrow(g)) return(NULL)
  pk <- iv_row$peak_pos
  dist <- if (is.na(pk)) NA_real_ else ifelse(pk >= g$start & pk <= g$end, 0, pmin(abs(g$start - pk), abs(g$end - pk)))
  out <- data.frame(interval_id = iv_row$interval_id, trait = iv_row$trait, interval_chr = iv_row$chr,
                    interval_start_bp = iv_row$start, interval_end_bp = iv_row$end, lead_snp = iv_row$lead_snp,
                    peak_pos_bp = pk, gene_id = g$gene_id, gene_name = g$gene_name, chr = g$chr,
                    gene_start_bp = g$start, gene_end_bp = g$end, strand = g$strand, biotype = g$biotype,
                    description = g$description, distance_to_peak_bp = dist,
                    peak_in_gene = if (is.na(pk)) NA else dist == 0, stringsAsFactors = FALSE)
  out[order(out$distance_to_peak_bp, out$gene_start_bp), , drop = FALSE]
}

empty_result <- function() {
  data.frame(interval_id = character(), trait = character(), interval_chr = character(),
             interval_start_bp = numeric(), interval_end_bp = numeric(), lead_snp = character(),
             peak_pos_bp = numeric(), gene_id = character(), gene_name = character(), chr = character(),
             gene_start_bp = numeric(), gene_end_bp = numeric(), strand = character(), biotype = character(),
             description = character(), distance_to_peak_bp = numeric(), peak_in_gene = logical(),
             stringsAsFactors = FALSE)
}

lookup_genes <- function(iv, genes) {
  rows <- lapply(seq_len(nrow(iv)), function(k) annotate_overlaps(iv[k, ], genes))
  rows <- rows[!vapply(rows, is.null, NA)]
  if (!length(rows)) return(empty_result())
  do.call(rbind, rows)
}

# --- Plot ----------------------------------------------------------------------------------

plot_tracks <- function(iv, res, out, max_panels = 12) {
  k <- min(nrow(iv), max_panels)
  if (k == 0L) return(invisible(NULL))
  nc <- min(3L, k); nr <- ceiling(k / nc)
  open_png(out, "gene_tracks.png", width = 4 * nc, height = 2.4 * nr + 0.4)
  op <- graphics::par(mfrow = c(nr, nc), mar = c(4, 1, 2.5, 1), col.axis = COL_INK2, fg = COL_INK2)
  for (i in seq_len(k)) {
    r <- iv[i, ]
    g <- if (is.null(res) || !nrow(res)) res else res[res$interval_id == r$interval_id, , drop = FALSE]
    u <- if (r$end - r$start >= 2e6) list(div = 1e6, lab = "Mb") else list(div = 1e3, lab = "kb")
    graphics::plot(NA, xlim = c(r$start, r$end) / u$div, ylim = c(0, 3), yaxt = "n", bty = "n",
                   xlab = sprintf("Chromosome %s (%s)", r$target_chr, u$lab), ylab = "")
    graphics::title(main = sprintf("%s: %d gene(s)", r$interval_id, if (is.null(g)) 0L else nrow(g)),
                    cex.main = 0.9, col.main = COL_INK, adj = 0)
    if (!is.null(g) && nrow(g)) {
      y <- ifelse(g$strand == "-", 1, 2)
      graphics::segments(pmax(g$gene_start_bp, r$start) / u$div, y, pmin(g$gene_end_bp, r$end) / u$div, y,
                         lwd = 6, col = COL_BLUE, lend = 1)
    }
    graphics::mtext(c("- strand", "+ strand"), side = 2, at = c(1, 2), las = 1, line = -0.5, cex = 0.6, col = COL_INK2, adj = 0)
    if (!is.na(r$peak_pos)) graphics::abline(v = r$peak_pos / u$div, col = COL_ORANGE, lwd = 2, lty = 2)
  }
  graphics::par(op)
  grDevices::dev.off()
}

# --- Main ------------------------------------------------------------------------------------

run_lookup <- function(opt, fetcher = http_get_text) {
  if (!nzchar(opt$assembly))
    stop("State the genome assembly of your coordinates with --assembly (for example Glycine_max_v2.1, ",
         "Wm82.a4.v1, Zm-B73-REFERENCE-NAM-5.0, IRGSP-1.0). Ask the user; never guess it.", call. = FALSE)
  if (nzchar(opt$gff) == nzchar(opt$ensembl_species))
    stop("Give exactly one gene source: --gff FILE or --ensembl-species NAME.", call. = FALSE)
  if (nzchar(opt$intervals) == nzchar(opt$snps))
    stop("Give exactly one region source: --intervals FILE or --snps FILE.", call. = FALSE)
  iv <- if (nzchar(opt$intervals)) read_intervals(opt$intervals, opt$coords) else
    snp_windows(opt$snps, opt$window_bp, opt$ld_table)
  notes <- attr(iv, "note")
  warn <- character()
  chr_map <- if (nzchar(opt$chr_map)) {
    cm <- read_csv_checked(opt$chr_map, "--chr-map")
    if (!all(c("user_chr", "target_chr") %in% names(cm))) stop("--chr-map needs columns user_chr and target_chr.", call. = FALSE)
    cm
  } else NULL
  raw_dir <- file.path(opt$out, "raw_responses")

  if (nzchar(opt$gff)) {
    gff <- read_gff3(opt$gff, trimws(strsplit(opt$feature_type, ",", fixed = TRUE)[[1]]))
    genes <- gff$genes
    chk <- check_gff_assembly(opt$assembly, gff$pragmas)
    if (is.na(chk$ok)) warn <- c(warn, sprintf(
      "The GFF3 header does not name its genome build, so the assembly could not be checked. The user stated '%s': confirm the annotation file is for exactly that assembly.",
      opt$assembly))
    hint <- maize_id_check(genes$gene_id, opt$assembly)
    if (!is.na(hint)) notes <- c(notes, paste("Gene ID check:", hint))
    source_desc <- sprintf("GFF3 file %s (%s features%s)", basename(opt$gff), paste(unique(genes$biotype), collapse = "/"),
                           if (isTRUE(chk$ok)) sprintf("; header build: %s", chk$stated) else "")
    targets <- unique(genes$chr)
    seq_len_tab <- gff$seq_lengths
  } else {
    info <- ensembl_assembly_info(opt$ensembl_species, opt$ensembl_server, fetcher, raw_dir)
    check_ensembl_assembly(opt$assembly, info)
    targets <- info$regions$name
    seq_len_tab <- info$regions
    source_desc <- sprintf("Ensembl REST %s, species %s, assembly %s (%s)", opt$ensembl_server, opt$ensembl_species,
                           info$assembly_name, info$assembly_accession)
  }
  cm <- map_chromosomes(iv$chr, targets, chr_map)
  iv$target_chr <- cm$target_chr[match(iv$chr, cm$user_chr)]
  renamed <- cm[cm$method == "name normalised", , drop = FALSE]
  if (nrow(renamed)) notes <- c(notes, sprintf("Chromosome names matched by number: %s.",
                                               paste(sprintf("%s -> %s", renamed$user_chr, renamed$target_chr), collapse = ", ")))
  if (!is.null(seq_len_tab)) {
    L <- seq_len_tab$length[match(iv$target_chr, seq_len_tab$name)]
    probe <- ifelse(is.na(iv$peak_pos), iv$start, iv$peak_pos)
    beyond <- !is.na(L) & probe > L
    if (any(beyond)) stop("Position(s) beyond the end of the chromosome in this assembly (e.g. ", iv$interval_id[beyond][1],
                          " at ", format(probe[beyond][1], scientific = FALSE), " bp; chromosome ", iv$target_chr[beyond][1],
                          " is ", format(L[beyond][1], scientific = FALSE), " bp). The coordinates are probably from a ",
                          "different assembly.", call. = FALSE)
    clip <- !is.na(L) & iv$end > L
    if (any(clip)) {
      notes <- c(notes, sprintf("%d interval(s) extended past the chromosome end were clipped to its length.", sum(clip)))
      iv$end[clip] <- L[clip]
    }
  }
  if (nzchar(opt$ensembl_species)) {
    genes <- do.call(rbind, lapply(seq_len(nrow(iv)), function(k)
      ensembl_genes(opt$ensembl_species, iv$target_chr[k], iv$start[k], iv$end[k], opt$ensembl_server, fetcher,
                    opt$biotype, raw_dir)))
    genes <- genes[!duplicated(paste(genes$gene_id, genes$chr)), , drop = FALSE]
    wrong <- !is.na(genes$assembly_name) & genes$assembly_name != info$assembly_name
    if (any(wrong)) stop("Ensembl returned genes from assembly ", genes$assembly_name[wrong][1], ", not ",
                         info$assembly_name, ". Stopping rather than mixing assemblies.", call. = FALSE)
    hint <- maize_id_check(genes$gene_id, opt$assembly)
    if (!is.na(hint)) notes <- c(notes, paste("Gene ID check:", hint))
  }
  res <- lookup_genes(iv, genes)
  iv$n_genes <- vapply(iv$interval_id, function(id) sum(res$interval_id == id), 0L)
  iv$length_bp <- iv$end - iv$start + 1
  list(intervals = iv, genes = res, chr_map = cm, notes = notes, warn = warn, source = source_desc)
}

main <- function(argv = commandArgs(trailingOnly = TRUE), fetcher = http_get_text) {
  if (any(argv %in% c("-h", "--help"))) {
    cat(USAGE)
    return(invisible(NULL))
  }
  opt <- parse_args(argv)
  options(scipen = 100)
  if (opt$demo) {
    opt$intervals <- file.path(SKILL_DIR, "demo", "demo_qtl_regions.csv")
    opt$gff <- file.path(SKILL_DIR, "demo", "demo_annotation.gff3")
    if (!nzchar(opt$assembly)) opt$assembly <- "AgriClaw-demo-1.0"
  }
  if (!nzchar(opt$out)) opt$out <- file.path("agriclaw_output", SKILL)
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  r <- run_lookup(opt, fetcher)
  iv <- r$intervals; res <- r$genes
  save_table(res, opt$out, "genes_in_intervals.csv")
  save_table(iv[, c("interval_id", "trait", "chr", "target_chr", "start", "end", "length_bp", "peak_pos", "lead_snp", "n_genes")],
             opt$out, "intervals.csv")
  save_table(r$chr_map, opt$out, "chromosome_map.csv")
  plot_tracks(iv, res, opt$out)

  warn <- r$warn
  empty <- iv$interval_id[iv$n_genes == 0]
  if (length(empty)) warn <- c(warn, sprintf("%d interval(s) contain no genes: %s. Check the window size and the chromosome names.",
                                             length(empty), paste(utils::head(empty, 5), collapse = ", ")))
  big <- iv$interval_id[iv$n_genes > 300]
  if (length(big)) warn <- c(warn, sprintf("%d interval(s) contain more than 300 genes (%s): the window is too wide to suggest candidates; consider a narrower LD-based window.",
                                           length(big), paste(utils::head(big, 5), collapse = ", ")))
  top <- if (nrow(res)) {
    near <- res[!is.na(res$distance_to_peak_bp), , drop = FALSE]
    near <- near[order(near$interval_id, near$distance_to_peak_bp), , drop = FALSE]
    near <- near[!duplicated(near$interval_id), , drop = FALSE]
    if (nrow(near)) c("| interval | nearest gene | distance to peak (bp) | description |", "|---|---|---|---|",
                      sprintf("| %s | %s | %.0f | %s |", near$interval_id, near$gene_id, near$distance_to_peak_bp,
                              ifelse(is.na(near$description) | !nzchar(near$description), "-", substr(near$description, 1, 80))))
  } else NULL
  report <- c(
    sprintf("# %s report", SKILL), "",
    if (length(warn)) c("## Warnings", "", paste("-", warn), "") else NULL,
    "## Summary", "",
    sprintf("- Assembly of the user's coordinates: **%s**", opt$assembly),
    sprintf("- Gene source: %s", r$source),
    sprintf("- Intervals: %d; genes found: %d (%d unique)", nrow(iv), nrow(res), length(unique(res$gene_id))),
    paste("-", r$notes), "",
    "| interval | chr | start (bp) | end (bp) | genes |", "|---|---|---|---|---|",
    sprintf("| %s | %s | %.0f | %.0f | %d |", iv$interval_id, iv$target_chr, iv$start, iv$end, iv$n_genes), "",
    if (!is.null(top)) c("Nearest gene to each peak:", "", top, "") else NULL,
    "![Gene tracks](figures/gene_tracks.png)", "",
    "## Interpretation", "",
    "- A gene in an interval is a positional candidate only. Physical distance to the peak SNP is not evidence of causality: LD can extend across the whole interval, and the causal variant may be regulatory or in an unannotated gene.",
    "- Gene descriptions come from the annotation source and may be computational predictions.", "",
    "## References", "",
    "- Yates A et al. (2015) The Ensembl REST API: Ensembl data for any language. Bioinformatics 31:143-145. https://doi.org/10.1093/bioinformatics/btu613",
    "- Yates AD et al. (2022) Ensembl Genomes 2022: an expanding genome resource for non-vertebrates. Nucleic Acids Res 50:D996-D1003. https://doi.org/10.1093/nar/gkab1007",
    "- GFF3 specification: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md", "",
    sprintf("*%s %s*", SKILL, VERSION))
  write_utf8(report, file.path(opt$out, "report.md"))
  write_repro_bundle(opt$out, argv, inputs = c(opt$intervals, opt$snps, opt$gff, opt$ld_table, opt$chr_map))
  cat(sprintf("Done: %d interval(s), %d gene hit(s). Report: %s\n", nrow(iv), nrow(res),
              normalizePath(file.path(opt$out, "report.md"), winslash = "/")))
  invisible(r)
}

if (sys.nframe() == 0L) main()
