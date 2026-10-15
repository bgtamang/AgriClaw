# Tests for claw-gene-lookup. Run from anywhere: Rscript tests/test_gene_lookup.R
# Offline known answers: a hand-made GFF3 with genes on every interval boundary,
# BED versus 1-based coordinates, chromosome-name matching, assembly checks,
# and saved Ensembl REST responses parsed through a mock fetcher. A live
# Ensembl test runs only when rest.ensembl.org answers, and is skipped otherwise.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "gene_lookup.R"))  # loads functions; main() does not run
FIX <- file.path(SKILL_TEST_DIR, "tests", "fixtures")

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "gene_lookup.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
tmp <- file.path(tempdir(), "claw-gene-lookup-test")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE)
opts <- function(...) { o <- DEFAULTS; o$out <- file.path(tmp, "run"); m <- list(...); o[names(m)] <- m; o }
write_csv <- function(df, name) { p <- file.path(tmp, name); utils::write.csv(df, p, row.names = FALSE); p }
gff <- file.path(FIX, "mini.gff3")

# 1. GFF3 parsing ---------------------------------------------------------------------------
g <- read_gff3(gff)
check(setequal(g$genes$gene_id, c("geneA", "geneB", "geneC", "geneD", "geneE", "geneF", "geneG")),
      "gene features read; mRNA/exon/pseudogene skipped; line after ##FASTA ignored; 'gene:' prefix removed")
check(g$genes$description[g$genes$gene_id == "geneC"] == "Zinc finger, C2H2; nuclear", "URL-encoded attributes decoded")
check(g$genes$biotype[g$genes$gene_id == "geneF"] == "lncRNA" && g$genes$biotype[g$genes$gene_id == "geneA"] == "gene",
      "biotype attribute used when present, feature type otherwise")
check(g$pragmas[["genome-build"]] == "AgriClaw-mini-1.0" && identical(g$seq_lengths$length, c(1e5, 5e4)),
      "genome-build pragma and sequence-region lengths parsed")
g2 <- read_gff3(gff, c("gene", "pseudogene"))
check("pseudo1" %in% g2$genes$gene_id, "--feature-type gene,pseudogene includes pseudogenes")

# 2. Overlap boundaries and distances (1-based, inclusive) ---------------------------------------
f_iv <- write_csv(data.frame(name = "q1", chr = "Gm01", start = 1000, end = 2000, peak_pos = 1550), "iv.csv")
r <- run_lookup(opts(intervals = f_iv, gff = gff, assembly = "AgriClaw-mini-1.0"))
check(setequal(r$genes$gene_id, c("geneB", "geneC", "geneD", "geneF")),
      "interval 1000-2000 keeps genes touching either end and the spanning gene; drops 999 and 2001 neighbours and other chromosomes")
dist <- stats::setNames(r$genes$distance_to_peak_bp, r$genes$gene_id)
check(identical(unname(dist[c("geneB", "geneC", "geneD", "geneF")]), c(550, 0, 450, 0)),
      "distance to peak: 0 inside the gene, otherwise the gap to the nearest gene end")
check(identical(r$genes$gene_id[1:2], c("geneC", "geneF")) || identical(r$genes$gene_id[1:2], c("geneF", "geneC")),
      "genes are sorted by distance to the peak")

# 3. BED and coordinate conventions ---------------------------------------------------------------
f_bed <- file.path(tmp, "iv.bed")
writeLines(c("track name=test", "Gm01\t999\t2000\tq1"), f_bed)
rb <- run_lookup(opts(intervals = f_bed, gff = gff, assembly = "AgriClaw-mini-1.0"))
check(rb$intervals$start == 1000 && rb$intervals$end == 2000 && setequal(rb$genes$gene_id, c("geneB", "geneC", "geneD", "geneF")),
      "BED 999-2000 (0-based, half-open) is 1000-2000 in 1-based coordinates")
f_iv0 <- write_csv(data.frame(chr = "Gm01", start = 999, end = 2000), "iv0.csv")
r0 <- run_lookup(opts(intervals = f_iv0, gff = gff, assembly = "AgriClaw-mini-1.0", coords = "bed"))
check(r0$intervals$start == 1000, "--coords bed shifts CSV starts by +1")
r1 <- run_lookup(opts(intervals = f_iv0, gff = gff, assembly = "AgriClaw-mini-1.0"))
check("geneA" %in% r1$genes$gene_id, "without --coords bed the same CSV start 999 reaches the gene ending at 999")

# 4. SNP windows --------------------------------------------------------------------------------------
f_snp <- write_csv(data.frame(snp = "S1_1550", chr = "Gm01", pos = 1550), "snp.csv")
rs <- run_lookup(opts(snps = f_snp, window_bp = 450, gff = gff, assembly = "AgriClaw-mini-1.0"))
check(rs$intervals$start == 1100 && rs$intervals$end == 2000 && setequal(rs$genes$gene_id, c("geneC", "geneD", "geneF")),
      "SNP at 1550 +/- 450 bp gives 1100-2000 and genes C, D, F")
f_ld <- write_csv(data.frame(chr = c("all", "Gm01"), ld_distance_bp = c(10, 600)), "ld.csv")
rl <- run_lookup(opts(snps = f_snp, ld_table = f_ld, gff = gff, assembly = "AgriClaw-mini-1.0"))
check(rl$intervals$start == 950 && rl$intervals$end == 2150, "SNP windows from a claw-ld-decay ld_window.csv")
err <- tryCatch(run_lookup(opts(snps = f_snp, gff = gff, assembly = "AgriClaw-mini-1.0")), error = conditionMessage)
check(is.character(err) && grepl("window", err), "SNP mode without a window stops and asks for one")

# 5. Chromosome names -----------------------------------------------------------------------------------
check(identical(chrom_key(c("Gm01", "Chr01", "chr1", "1", "01", "Chromosome_1", "chr1A", "1A", "Mt")),
                c("1", "1", "1", "1", "1", "1", "1a", "1a", "mt")), "chromosome keys normalise prefixes and zeros")
for (nm in c("1", "Chr01", "chr1", "Gm1")) {
  f <- write_csv(data.frame(chr = nm, start = 1000, end = 2000), "nm.csv")
  rn <- run_lookup(opts(intervals = f, gff = gff, assembly = "AgriClaw-mini-1.0"))
  check(rn$intervals$target_chr == "Gm01" && nrow(rn$genes) == 4L, paste0("chromosome '", nm, "' matched to Gm01"))
}
err <- tryCatch(map_chromosomes("1", c("Chr1", "1A", "chr01")), error = conditionMessage)
check(is.character(err) && grepl("Ambiguous", err) && grepl("--chr-map", err), "ambiguous names stop and point to --chr-map")
err <- tryCatch(map_chromosomes("scaffold_9", c("Gm01", "Gm02")), error = conditionMessage)
check(is.character(err) && grepl("not found", err), "unknown chromosome stops with the source's naming")
check(map_chromosomes("1A", c("1A", "1B", "1D"))$target_chr == "1A" &&
        map_chromosomes("chr2D", c("1A", "2B", "2D"))$target_chr == "2D", "wheat sub-genome names (1A, 2D) kept apart")

# 6. Assembly checks -------------------------------------------------------------------------------------
check(nrow(run_lookup(opts(intervals = f_iv, gff = gff, assembly = "agriclaw_mini_1.0"))$genes) == 4L,
      "assembly names compared without case or punctuation")
err <- tryCatch(run_lookup(opts(intervals = f_iv, gff = gff, assembly = "Wm82.a4.v1")), error = conditionMessage)
check(is.character(err) && grepl("Assembly mismatch", err) && grepl("AgriClaw-mini-1.0", err),
      "a different assembly is refused, naming both assemblies")
err <- tryCatch(run_lookup(opts(intervals = f_iv, gff = gff, assembly = "")), error = conditionMessage)
check(is.character(err) && grepl("--assembly", err), "missing --assembly stops")
f_far <- write_csv(data.frame(chr = "Gm01", start = 140000, end = 160000, peak_pos = 150000), "far.csv")
err <- tryCatch(run_lookup(opts(intervals = f_far, gff = gff, assembly = "AgriClaw-mini-1.0")), error = conditionMessage)
check(is.character(err) && grepl("beyond the end", err), "a position past the chromosome end is refused (wrong assembly)")
f_clip <- write_csv(data.frame(chr = "Gm01", start = 99000, end = 120000, peak_pos = 99500), "clip.csv")
rc <- run_lookup(opts(intervals = f_clip, gff = gff, assembly = "AgriClaw-mini-1.0"))
check(rc$intervals$end == 1e5 && any(grepl("clipped", rc$notes)), "a window running past the chromosome end is clipped")
nopr <- file.path(tmp, "nopragma.gff3")
writeLines(readLines(gff)[-2], nopr)
rw <- run_lookup(opts(intervals = f_iv, gff = nopr, assembly = "AgriClaw-mini-1.0"))
check(length(rw$warn) == 1L && grepl("could not be checked", rw$warn), "GFF3 without a build pragma: warning, not silence")
zm <- file.path(tmp, "maize.gff3")
writeLines(c("##gff-version 3", "1\tx\tgene\t100\t900\t.\t+\t.\tID=gene:Zm00001eb000010",
             "1\tx\tgene\t1000\t1900\t.\t+\t.\tID=gene:Zm00001eb000020"), zm)
f_zm <- write_csv(data.frame(chr = "chr1", start = 1, end = 5000), "zm.csv")
err <- tryCatch(run_lookup(opts(intervals = f_zm, gff = zm, assembly = "B73_RefGen_v4")), error = conditionMessage)
check(is.character(err) && grepl("Zm00001eb", err), "maize v5 gene IDs with v4 coordinates are refused")

# 7. Ensembl REST, offline with saved responses ----------------------------------------------------------------
info_txt <- paste(readLines(file.path(FIX, "ensembl_info_assembly_glycine_max.json")), collapse = "\n")
ov_txt <- paste(readLines(file.path(FIX, "ensembl_overlap_glycine_max_1_1000000-1060000.json")), collapse = "\n")
urls <- character()
mock <- function(url) {
  urls <<- c(urls, url)
  if (grepl("/info/assembly/glycine_max", url, fixed = TRUE)) return(info_txt)
  if (grepl("/overlap/region/glycine_max/1:1000000-1060000", url, fixed = TRUE)) return(ov_txt)
  if (grepl("/overlap/region/", url, fixed = TRUE)) return("[]")
  stop("unexpected URL in offline test: ", url)
}
pg <- parse_ensembl_genes(jsonlite::fromJSON(ov_txt))
check(nrow(pg) == 9L && all(grepl("^GLYMA_01G", pg$gene_id)), "saved overlap response parses to 9 soybean genes")
k <- match("GLYMA_01G010400", pg$gene_id)
check(pg$start[k] == 1009200 && pg$end[k] == 1013846 && pg$strand[k] == "+" && pg$description[k] == "DEAD-box RNA helicase",
      "gene coordinates, strand and description (source tag removed) parsed")
f_soy <- write_csv(data.frame(chr = "Gm01", start = 1000000, end = 1060000, peak_pos = 1011000), "soy.csv")
re <- run_lookup(opts(intervals = f_soy, ensembl_species = "glycine_max", assembly = "Glycine_max_v2.1"), fetcher = mock)
check(nrow(re$genes) == 9L && re$intervals$target_chr == "1", "offline Ensembl lookup: Gm01 mapped to '1', 9 genes")
check(re$genes$gene_id[1] == "GLYMA_01G010400" && re$genes$distance_to_peak_bp[1] == 0, "peak inside GLYMA_01G010400 gives distance 0")
check(file.exists(file.path(tmp, "run", "raw_responses", "info_assembly_glycine_max.json")), "raw Ensembl responses are kept in --out")
err <- tryCatch(run_lookup(opts(intervals = f_soy, ensembl_species = "glycine_max", assembly = "Wm82.a4.v1"), fetcher = mock),
                error = conditionMessage)
check(is.character(err) && grepl("Glycine_max_v2.1", err) && grepl("GCA_000004515.4", err),
      "Ensembl assembly mismatch is refused, naming the served assembly and accession")
check(nrow(run_lookup(opts(intervals = f_soy, ensembl_species = "glycine_max", assembly = "GCA_000004515.4"), fetcher = mock)$genes) == 9L,
      "the GenBank accession is accepted as the assembly name")
f_out <- write_csv(data.frame(chr = "1", start = 59000000, end = 61000000, peak_pos = 6e7), "soy_far.csv")
err <- tryCatch(run_lookup(opts(intervals = f_out, ensembl_species = "glycine_max", assembly = "Glycine_max_v2.1"), fetcher = mock),
                error = conditionMessage)
check(is.character(err) && grepl("beyond the end", err), "position past the end of soybean chromosome 1 (56.8 Mb) is refused")
ch <- region_chunks(1, 12e6)
check(identical(ch$start, c(1, 5000001, 10000001)) && identical(ch$end, c(5e6, 1e7, 12e6)), "regions are split into <= 5 Mb requests")
urls <- character()
f_long <- write_csv(data.frame(chr = "1", start = 1, end = 12e6), "soy_long.csv")
invisible(run_lookup(opts(intervals = f_long, ensembl_species = "glycine_max", assembly = "Glycine_max_v2.1"), fetcher = mock))
ov_urls <- grep("/overlap/region/", urls, value = TRUE)
check(length(ov_urls) == 3L && grepl("1:1-5000000", ov_urls[1]) && grepl("1:10000001-12000000", ov_urls[3]),
      "a 12 Mb interval is fetched as three overlap requests")

# 8. End to end on the demo data -----------------------------------------------------------------------------
out <- file.path(tmp, "demo")
res <- run_cli(c("--demo", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "tables/genes_in_intervals.csv", "tables/intervals.csv", "tables/chromosome_map.csv",
            "figures/gene_tracks.png", "reproducibility/commands.sh", "reproducibility/sessionInfo.txt")) {
  check(file.exists(file.path(out, f)), paste("demo writes", f))
}
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
gi <- utils::read.csv(file.path(out, "tables", "genes_in_intervals.csv"), stringsAsFactors = FALSE)
ann <- utils::read.delim(file.path(SKILL_TEST_DIR, "demo", "demo_annotation.gff3"), header = FALSE, comment.char = "#",
                         stringsAsFactors = FALSE)
ann <- ann[ann$V3 == "gene", ]
qr <- utils::read.csv(file.path(SKILL_TEST_DIR, "demo", "demo_qtl_regions.csv"), stringsAsFactors = FALSE)
expected <- unlist(lapply(seq_len(nrow(qr)), function(i) {
  a <- ann[ann$V1 == sprintf("Chr%02d", as.integer(qr$chr[i])) & ann$V4 <= qr$end_bp[i] & ann$V5 >= qr$start_bp[i], ]
  paste(qr$trait[i], qr$region_id[i], sub("^ID=([^;]+);.*$", "\\1", a$V9))
}))
check(setequal(paste(gi$trait, sub("^.*:", "", gi$interval_id), gi$gene_id), expected),
      "demo genes equal an independent brute-force overlap of the GFF3")
bad <- run_cli(c("--demo", "--assembly", "IRGSP-1.0", "--out", shQuote(file.path(tmp, "bad"))))
check(bad$status != 0L && any(grepl("Assembly mismatch", bad$output)), "demo with the wrong assembly exits with an error")

# 9. Live Ensembl test (skipped when offline) ------------------------------------------------------------------
# Bounded live fetcher (20 s timeout, 2 tries) so a slow server cannot stall the suite.
live_fetch <- function(url) {
  req <- httr2::req_headers(httr2::request(url), Accept = "application/json")
  req <- httr2::req_retry(httr2::req_timeout(req, 20), max_tries = 2,
                          is_transient = function(resp) httr2::resp_status(resp) %in% c(429, 500, 502, 503, 504))
  tryCatch(httr2::resp_body_string(httr2::req_perform(req)),
           error = function(e) stop("Ensembl request failed: ", conditionMessage(e), call. = FALSE))
}
live <- tryCatch(grepl("ping", live_fetch("https://rest.ensembl.org/info/ping?content-type=application/json")),
                 error = function(e) FALSE)
if (isTRUE(live)) {
  lv <- tryCatch({
    inf <- ensembl_assembly_info("glycine_max", "https://rest.ensembl.org", live_fetch)
    list(info = inf, res = run_lookup(opts(intervals = f_soy, ensembl_species = "glycine_max", assembly = inf$assembly_name,
                                           out = file.path(tmp, "live")), fetcher = live_fetch))
  }, error = function(e) e)
  if (inherits(lv, "error") && grepl("Ensembl request failed", conditionMessage(lv))) {
    cat("SKIP: live Ensembl test (server error):", conditionMessage(lv), "\n")
  } else {
    if (inherits(lv, "error")) stop("FAIL: live Ensembl lookup: ", conditionMessage(lv), call. = FALSE)
    lg <- lv$res$genes
    check(nrow(lg) > 0L && all(lg$gene_start_bp <= 1060000 & lg$gene_end_bp >= 1000000) && all(lg$chr == "1"),
          sprintf("live Ensembl: %d genes, all overlapping 1:1000000-1060000 on %s", nrow(lg), lv$info$assembly_name))
  }
} else {
  cat("SKIP: live Ensembl test (rest.ensembl.org not reachable)\n")
}

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
