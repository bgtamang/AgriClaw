# Generate the claw-gene-lookup demo data (synthetic, seed 21):
#   demo_annotation.gff3  a GFF3 annotation for an invented 3-chromosome genome
#                         ("AgriClaw-demo-1.0"; chromosomes Chr01-Chr03), with
#                         gene, mRNA and exon features and URL-encoded notes;
#   demo_qtl_regions.csv  QTL regions in the format written by claw-gwas-crop
#                         (chromosomes named 1-3, as in many marker files).
#
#   Rscript demo/make_demo.R

here <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/"))
})
set.seed(21)
chr_len <- c(10e6, 8e6, 12e6)
desc <- c("MYB transcription factor", "Cytochrome P450, family 71", "Auxin response factor",
          "Pentatricopeptide repeat protein", "Leucine-rich repeat receptor kinase",
          "Zinc finger protein, C2H2 type", "Hypothetical protein", "Sugar transporter; SWEET family",
          "Heat shock protein 70", "NAC domain protein")
enc <- function(x) gsub(";", "%3B", gsub(",", "%2C", gsub("=", "%3D", x, fixed = TRUE), fixed = TRUE), fixed = TRUE)
lines <- c("##gff-version 3", "#!genome-build AgriClaw-demo-1.0",
           sprintf("##sequence-region Chr%02d 1 %.0f", 1:3, chr_len))
for (ch in 1:3) {
  starts <- sort(sample.int(chr_len[ch] - 10000, 40))
  for (g in seq_along(starts)) {
    s <- starts[g]; e <- s + sample(800:8000, 1); str <- sample(c("+", "-"), 1)
    id <- sprintf("AcDemo%02dg%06d", ch, g * 100)
    lines <- c(lines,
               sprintf("Chr%02d\tdemo\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s;Note=%s", ch, s, e, str, id, id, enc(sample(desc, 1))),
               sprintf("Chr%02d\tdemo\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s.1;Parent=%s", ch, s, e, str, id, id),
               sprintf("Chr%02d\tdemo\texon\t%d\t%d\t.\t%s\t.\tParent=%s.1", ch, s, e, str, id))
  }
}
writeLines(lines, file.path(here, "demo_annotation.gff3"))
regions <- data.frame(trait = c("height_cm", "height_cm", "protein_pct"), region_id = c("R01", "R02", "R01"),
                      chr = c("1", "3", "2"), start_bp = c(2000000, 6500000, 1), end_bp = c(2600000, 7300000, 450000),
                      lead_snp = c("S1_2310440", "S3_6912003", "S2_150020"), lead_pos_bp = c(2310440, 6912003, 150020),
                      lead_p = c(3.1e-9, 8.4e-7, 2.2e-6), n_significant_snps = c(4L, 2L, 1L), ld_distance_bp = 300000)
utils::write.csv(regions, file.path(here, "demo_qtl_regions.csv"), row.names = FALSE)
cat("Wrote demo_annotation.gff3 and demo_qtl_regions.csv to", here, "\n")
