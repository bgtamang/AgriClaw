# Tests for claw-enrichment. Run from anywhere: Rscript tests/test_enrichment.R
# Known answers: hand-computed hypergeometric p-values and BH adjustment, planted
# enrichment recovered and null sets not, the effect of the universe on p-values,
# GO propagation on a hand-built DAG and on GO.db, ID-mismatch diagnosis, annotation
# parsing, and an end-to-end command-line run on the demo.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "enrichment.R"))  # loads functions only

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
expect_error <- function(expr, pattern, what) {
  msg <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  check(!is.na(msg) && grepl(pattern, msg), paste0(what, " (got: ", msg, ")"))
}
run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "enrichment.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
tmp <- file.path(tempdir(), "claw-enrichment-test")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE)

# 1. Hypergeometric p-values and BH adjustment by hand ---------------------------
g <- sprintf("G%03d", 1:200)
t2g <- rbind(data.frame(term = "T1", gene = g[1:20]), data.frame(term = "T2", gene = g[21:50]),
             data.frame(term = "T3", gene = g[51:60]), data.frame(term = "T4", gene = g[61:65]),
             data.frame(term = "T5", gene = g))
de <- c(g[1:10], g[51:52], g[100:107])
r <- run_ora(de, universe = g, t2g = t2g, min_size = 10, max_size = 150)
tab <- r$table
N <- 200; n <- 20
p1 <- stats::phyper(10 - 1, 20, N - 20, n, lower.tail = FALSE)
p3 <- stats::phyper(2 - 1, 10, N - 10, n, lower.tail = FALSE)
check(r$n_eligible == 3L, "size limits: T4 (5 genes) and T5 (200 genes) excluded, 3 terms eligible")
check(setequal(tab$term_id, c("T1", "T3")), "only terms with DE genes are listed (T2 has none)")
check(abs(tab$pvalue[tab$term_id == "T1"] - p1) < 1e-12 && abs(tab$pvalue[tab$term_id == "T3"] - p3) < 1e-12,
      "p-values equal phyper(k - 1, M, N - M, n, lower.tail = FALSE) by hand")
bh <- stats::p.adjust(c(p1, p3, 1), "BH")
check(abs(tab$padj[tab$term_id == "T1"] - bh[1]) < 1e-12 && abs(tab$padj[tab$term_id == "T3"] - bh[2]) < 1e-12,
      "BH adjustment counts every eligible term, including T2 with p = 1")
check(abs(tab$fold_enrichment[tab$term_id == "T1"] - (10 / 20) / (20 / 200)) < 1e-12, "fold enrichment (k/n)/(M/N)")
check(all(tab$n_universe_annotated == N) && all(tab$n_de_annotated == n), "N and n reported")

# 2. Planted enrichment is recovered; null gene sets are not ------------------------
set.seed(3)
genome <- sprintf("Gx%05d", 1:4000)
uni <- sample(genome, 3000)
terms <- lapply(1:50, function(i) sample(genome, if (i <= 3) 80 else sample(15:150, 1)))  # S01-S03 planted
names(terms) <- sprintf("S%02d", 1:50)
t2g2 <- data.frame(term = rep(names(terms), lengths(terms)), gene = unlist(terms))
t2g2 <- rbind(t2g2, data.frame(term = "ALL", gene = genome))      # every gene annotated
planted <- c("S01", "S02", "S03")
de2 <- unique(unlist(lapply(planted, function(p) {
  m <- intersect(terms[[p]], uni)
  sample(m, round(0.4 * length(m)))
})))
de2 <- unique(c(de2, sample(setdiff(uni, unlist(terms[planted])), 250 - length(de2))))
r2 <- run_ora(de2, uni, t2g2, min_size = 10, max_size = 500)$table
check(all(r2$padj[match(planted, r2$term_id)] < 1e-3), "planted terms recovered at padj < 0.001")
false_pos <- sum(r2$padj < 0.05 & !r2$term_id %in% planted)
cat(sprintf("  planted terms max padj %.1e; null terms at padj < 0.05: %d of 47\n",
            max(r2$padj[match(planted, r2$term_id)]), false_pos))
check(false_pos <= 1L, sprintf("null terms: %d of 47 significant at padj < 0.05 (allowed <= 1)", false_pos))

# 3. The universe matters: genome background inflates expression-linked terms -------
set.seed(4)
genome3 <- sprintf("Gy%05d", 1:4000)
expressed <- genome3[1:1000]
hk <- c(sample(expressed, 140), sample(setdiff(genome3, expressed), 10))   # well-expressed term
other <- lapply(1:30, function(i) sample(genome3, 80))
t2g3 <- rbind(data.frame(term = "HK", gene = hk),
              data.frame(term = rep(sprintf("O%02d", 1:30), 80), gene = unlist(other)),
              data.frame(term = "ALL", gene = genome3))
de3 <- sample(expressed, 150)                                      # unrelated to HK by construction
k <- sum(de3 %in% hk)
p_genome_hand <- stats::phyper(k - 1, 150, 4000 - 150, 150, lower.tail = FALSE)
p_tested_hand <- stats::phyper(k - 1, 140, 1000 - 140, 150, lower.tail = FALSE)
pg <- run_ora(de3, genome3, t2g3, max_size = 5000)$table
pt <- run_ora(de3, expressed, t2g3, max_size = 5000)$table
p_g <- pg$pvalue[pg$term_id == "HK"]
p_t <- pt$pvalue[pt$term_id == "HK"]
check(abs(p_g - p_genome_hand) < 1e-12 && abs(p_t - p_tested_hand) < 1e-12,
      "universe p-values equal the hand hypergeometric values")
check(p_g < 1e-4 && p_t > 0.05, sprintf(paste0("genome background makes a well-expressed term look enriched ",
      "(p = %.1e) while the tested-gene universe does not (p = %.2f)"), p_g, p_t))
expect_error(check_sets_in_universe(list(all = c(expressed[1:5], genome3[3999])), expressed),
             "not in the universe", "DE genes outside the universe are refused")

# 4. GO propagation on a hand-built DAG -------------------------------------------
edges <- data.frame(child = c("T4", "T4", "T2", "T3", "T5", "T5", "T6"),
                    parent = c("T2", "T3", "T1", "T1", "T2", "T6", "T1"),
                    relation = c("isa", "part of", "isa", "isa", "regulates", "isa", "isa"))
ann <- data.frame(term = c("T4", "T5", "T1"), gene = c("g1", "g2", "g3"))
pr <- propagate_annotation(ann, edges)
got <- function(gene) sort(pr$term[pr$gene == gene])
check(identical(got("g1"), c("T1", "T2", "T3", "T4")), "g1 (T4) gains is_a and part_of ancestors T2, T3, T1")
check(identical(got("g2"), c("T1", "T5", "T6")), "g2 (T5) gains T6, T1 but not T2 (regulates is not followed)")
check(identical(got("g3"), "T1"), "a root annotation is unchanged")
check(!anyDuplicated(paste(pr$term, pr$gene)), "no duplicated gene-term pairs after propagation")
check(propagation_fraction(pr, edges) == 1 && propagation_fraction(ann, edges) == 0,
      "propagation_fraction is 1 for a propagated file and 0 for most-specific terms only")

# 5. Propagation with GO.db follows is_a/part_of only -----------------------------
go_edges <- go_parent_edges()
pg2 <- propagate_annotation(data.frame(term = c("GO:0010468", "GO:0042631"), gene = c("r1", "w1")), go_edges)
anc_r <- pg2$term[pg2$gene == "r1"]
anc_w <- pg2$term[pg2$gene == "w1"]
check("GO:0008150" %in% anc_r && "GO:0008150" %in% anc_w, "GO.db propagation reaches the BP root GO:0008150")
check(!"GO:0010467" %in% anc_r, "'regulation of gene expression' does not propagate to 'gene expression'")
check("GO:0010467" %in% AnnotationDbi::get("GO:0010468", GO.db::GOBPANCESTOR),
      "(GO.db's GOBPANCESTOR would have added it: the reason for custom propagation)")
check(all(c("GO:0009414", "GO:0009415", "GO:0006950") %in% anc_w),
      "cellular response to water deprivation gains response to water deprivation, response to water, response to stress")
gi <- go_term_info(c("GO:0036422", "GO:0009414", "GO:0018298"))
check(identical(unname(gi$remap["GO:0036422"]), "GO:0000010") && identical(gi$obsolete, "GO:0018298"),
      "secondary GO IDs mapped to primary (GOSYNONYM); obsolete IDs reported")

# 6. Gene ID mismatch diagnosis ------------------------------------------------------
a2 <- sprintf("Glyma.%02dG%06d", rep(1:5, each = 40), rep(seq(100, 4000, 100), 5))
a1 <- sprintf("Glyma%02dg%05d", rep(1:5, each = 40), rep(seq(10, 400, 10), 5))
check(any(grepl("Wm82.a1", diagnose_ids(a2, a1))), "soybean Wm82.a2 genes vs Wm82.a1 annotation diagnosed")
check(any(grepl("--strip-isoform", diagnose_ids(a2, paste0(a2, ".1")))), "transcript-level annotation diagnosed")
check(any(grepl("trailing", diagnose_ids(a2, paste0(a2, ".Wm82.a4.v1")))), "annotation-version suffix diagnosed")
check(identical(strip_isoform(c("Glyma.01G000100.1", "Glyma.01G000100.1.p", "Zm00001eb000010_T001",
                                "Zm00001eb000010_P002", "AT1G01010.2")),
                c("Glyma.01G000100", "Glyma.01G000100", "Zm00001eb000010", "Zm00001eb000010", "AT1G01010")),
      "isoform suffixes stripped for soybean, maize and Arabidopsis IDs")
de_tab <- file.path(tmp, "de_a2.csv")
utils::write.csv(data.frame(gene_id = a2, log2FoldChange = rep(c(2, -2, 0.1, 0.2), 50),
                            padj = rep(c(0.001, 0.002, 0.5, 0.9), 50)), de_tab, row.names = FALSE)
ann_a1 <- file.path(tmp, "ann_a1.csv")
utils::write.csv(data.frame(gene_id = a1, GO = "GO:0009414"), ann_a1, row.names = FALSE)
expect_error(main(c("--de-table", de_tab, "--annotation", ann_a1, "--out", file.path(tmp, "o1"))),
             "Wm82.a1", "a1-vs-a2 mismatch stops the run with the nomenclature diagnosis")

# 7. Annotation parsing ---------------------------------------------------------------
wide <- file.path(tmp, "wide.tsv")
writeLines(c("locusName\tGO", "g1\tGO:0000001,GO:0000002;GO:0000003|GO:0000004 GO:0000005", "g2\t", "g3\tGO:0000001"), wide)
w <- read_annotation(wide)$t2g
check(nrow(w) == 6L && sum(w$gene == "g1") == 5L && !"g2" %in% w$gene,
      "wide annotation: several separators split, empty cells skipped")
st <- split_terms(c("GO:0000001,GO:0000002", "photosynthesis - antenna proteins|ko00196", "GO:0000004 GO:0000005", NA))
check(identical(st, list(c("GO:0000001", "GO:0000002"), c("photosynthesis - antenna proteins", "ko00196"),
                         c("GO:0000004", "GO:0000005"), character(0))),
      "term cells split on , ; | and between GO IDs, but term names with spaces stay whole")
gaf <- file.path(tmp, "a.gaf")
writeLines(c("!gaf-version: 2.2",
             paste("TAIR", "AT1G01010", "NAC001", "involved_in", "GO:0009414", "TAIR:1", "IEA", "", "P", "", "", "protein", "taxon:3702", "20200101", "TAIR", "", "", sep = "\t"),
             paste("TAIR", "AT1G01010", "NAC001", "NOT|involved_in", "GO:0009737", "TAIR:1", "IMP", "", "P", "", "", "protein", "taxon:3702", "20200101", "TAIR", "", "", sep = "\t")), gaf)
gg <- read_annotation(gaf)
check(nrow(gg$t2g) == 1L && gg$t2g$term == "GO:0009414" && gg$t2g$gene == "AT1G01010",
      "GAF: NOT-qualified annotations are excluded, DB Object ID used as gene")

# 8. End to end on the demo -------------------------------------------------------------
out <- file.path(tmp, "demo_out")
res <- run_cli(c("--demo", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(utils::tail(res$output, 5), collapse = "\n")))
for (f in c("report.md", "tables/enrichment_up_BP.csv", "tables/enrichment_down_MF.csv",
            "tables/enrichment_all_CC.csv", "tables/annotation_coverage.csv", "tables/run_summary.csv",
            "figures/dotplot_up_BP.png", "reproducibility/commands.sh", "reproducibility/sessionInfo.txt"))
  check(file.exists(file.path(out, f)), paste("demo writes", f))
planted_up <- c("GO:0009414", "GO:0042631", "GO:0009737", "GO:0009269")
planted_down <- c("GO:0009765", "GO:0009768", "GO:0016168", "GO:0009535", "GO:0009523")
up <- utils::read.csv(file.path(out, "tables", "enrichment_up_BP.csv"))
check(setequal(up$term_id[up$significant], planted_up), "demo: exactly the 4 planted drought terms are significant for up genes")
all_sig <- unlist(lapply(Sys.glob(file.path(out, "tables", "enrichment_*.csv")), function(f) {
  d <- utils::read.csv(f); d$term_id[d$significant] }))
check(all(all_sig %in% c(planted_up, planted_down)), "demo: no null term is significant in any analysis")
rep_txt <- readLines(file.path(out, "report.md"))
check(any(grepl("--propagate", rep_txt)) && any(grepl("GO:0018298", rep_txt)),
      "demo report warns about unpropagated annotation and the obsolete GO ID")
check(run_cli(c("--demo", "--color", "red"))$status != 0L, "unknown option exits with an error")
noun <- run_cli(c("--genes", shQuote(de_tab), "--annotation", shQuote(ann_a1), "--out", shQuote(file.path(tmp, "o2"))))
check(noun$status != 0L && any(grepl("--universe is required", noun$output)), "a gene list without --universe is refused")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
