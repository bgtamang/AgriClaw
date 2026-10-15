# Tests for claw-licor-6800. Run from anywhere: Rscript tests/test_licor_6800.R
# Known answers: PhotoGEA's example LI-6800 files (row counts, serial numbers,
# setpoint sequence, values identical to PhotoGEA's reader), an "unevaluated
# Excel" copy of a real example file, hand-built rows with known QC flags, and
# designs with known confounding status. Ends with an end-to-end --demo run.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "licor_6800.R"))  # loads functions; main() does not run

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}

run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "licor_6800.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}

ex <- function(f) system.file("extdata", f, package = "PhotoGEA")
tmp <- file.path(tempdir(), "claw-licor-6800-test")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE)

# 1. PhotoGEA example files: structure and values ----------------------------
r1 <- read_licor_files(ex("c3_aci_1.xlsx"))
d1 <- r1$data
check(nrow(d1) == 48L, "c3_aci_1.xlsx has 48 logged points")
check(identical(r1$files$console_serial, "68C-831539") && identical(r1$files$head_serial, "68H-891539"),
      "serial numbers read from the file header")
pg <- PhotoGEA::read_gasex_file(ex("c3_aci_1.xlsx"))
check(isTRUE(all.equal(d1$a_umol_m2_s, pg[, "A"])) && isTRUE(all.equal(d1$ci_umol_mol, pg[, "Ci"])) &&
      isTRUE(all.equal(d1$gsw_mol_m2_s, pg[, "gsw"])) && isTRUE(all.equal(d1$tleaf_c, pg[, "TleafCnd"])),
      "A, Ci, gsw and leaf temperature identical to PhotoGEA's reader")
check(all(c("instrument", "species", "plot") %in% names(d1)), "user constants kept as columns")
key <- paste(d1$species, d1$plot)
check(length(unique(key)) == 3L && all(table(key) == 16L), "3 curves of 16 points")
check(identical(d1$co2_r_sp_umol_mol[1:16],
                c(400, 300, 200, 150, 100, 75, 50, 20, 400, 400, 600, 800, 1000, 1200, 1500, 1800)),
      "CO2 setpoint sequence of the first curve")

r2 <- read_licor_files(ex("c3_aci_2.xlsx"))
check(nrow(r2$data) == 48L && is.na(r2$files$console_serial), "c3_aci_2.xlsx: 48 rows, Remarks sheet deleted so no serial")
d2 <- assign_instrument(r2$data)
check(all(d2$instrument_id == "unknown:c3_aci_2.xlsx"), "missing serial gives instrument_id 'unknown:<file>'")
d1i <- assign_instrument(d1)
check(all(d1i$instrument_id == "68C-831539/68H-891539"), "instrument_id = console/head serial")
d1c <- assign_instrument(d1, "instrument")
check(all(d1c$instrument_id == "ripe1"), "--instrument-col overrides serial numbers")

rp <- read_licor_files(c(ex("plaintext_licor_file"), ex("plaintext_licor_file_v2")))
check(nrow(rp$data) == 192L, "plaintext files: 96 + 96 rows (v2 was closed and reopened mid-log)")
p1 <- rp$data[rp$data$source_file == "plaintext_licor_file", ]
p2 <- rp$data[rp$data$source_file == "plaintext_licor_file_v2", ]
check(isTRUE(all.equal(p1$a_umol_m2_s, p2$a_umol_m2_s)), "both plaintext versions give the same A values")
check(identical(unique(p1$head_serial), "68H-891539"), "plaintext header serial read")

both <- process_licor(c(ex("c3_aci_1.xlsx"), ex("c3_aci_2.xlsx")), c("instrument", "species", "plot"),
                      "species", instrument_col = "instrument")
check(nrow(both$data) == 96L && nrow(both$curves) == 6L, "two files combine to 96 rows and 6 curves")
check(all(both$curves$n_points == 16L) && !any(both$curves$spans_files), "each curve has 16 points in one file")
check(identical(both$confounding$summary$status, "crossed"), "demo design: species crossed with instrument")

# Duplicate observations and .xlsx / plaintext twins
dupres <- process_licor(c(ex("c3_aci_1.xlsx"), ex("plaintext_licor_file")))$data
isdup <- dupres$qc_duplicate %in% TRUE
check(sum(isdup) == 48L && all(dupres$source_file[isdup] == "plaintext_licor_file"),
      "48 rows of plaintext_licor_file repeat c3_aci_1.xlsx (same serial and timestamp): later copies flagged")
check(all(dupres$qc_hard_fail[isdup]) && !any(dupres$qc_duplicate[dupres$source_file == "c3_aci_1.xlsx"]),
      "duplicates are hard failures; the first copy is kept")
check(all(is.na(process_licor(ex("c3_aci_2.xlsx"))$data$qc_duplicate)),
      "duplicate check is not made when the instrument is unknown")
tw <- drop_xlsx_twins(c("d/log1.xlsx", "d/log1", "d/log2.xlsx"))
check(identical(as.character(tw), c("d/log1", "d/log2.xlsx")) && identical(attr(tw, "skipped"), "d/log1.xlsx"),
      "the .xlsx twin of a plaintext log is skipped")
logdir <- file.path(tmp, "logs")
dir.create(logdir)
file.copy(ex("c3_aci_1.xlsx"), file.path(logdir, "log1.xlsx"))
file.copy(ex("plaintext_licor_file"), file.path(logdir, "log1"))
writeLines("a,b", file.path(logdir, "notes.csv"))
writeLines("not a log", file.path(logdir, "readme"))
picked <- drop_xlsx_twins(list_input_files("", logdir))
check(identical(basename(as.character(picked)), "log1"),
      "--input-dir keeps LI-6800 logs only (no CSV or other text) and prefers the plaintext twin")

# 2. Unevaluated Excel formulas ------------------------------------------------
# Copy a real example file and strip the cached values of every formula cell,
# which is what an LI-6800 .xlsx looks like before it is opened in Excel.
make_unevaluated_copy <- function(src, dest, mode) {
  d <- tempfile("xlsx_")
  dir.create(d)
  utils::unzip(src, exdir = d)
  sheet <- file.path(d, "xl", "worksheets", "sheet1.xml")
  xml <- readChar(sheet, file.info(sheet)$size, useBytes = TRUE)
  repl <- if (mode == "zero") "\\1<v>0</v>" else "\\1"
  xml2 <- gsub("(<f[^>]*/>|<f[^>]*>[^<]*</f>)<v>[^<]*</v>", repl, xml, perl = TRUE, useBytes = TRUE)
  stopifnot(!identical(xml, xml2))
  con <- file(sheet, open = "wb")
  writeChar(xml2, con, eos = NULL, useBytes = TRUE)
  close(con)
  old <- setwd(d)
  on.exit(setwd(old))
  zip::zip(dest, files = list.files(".", recursive = TRUE, all.files = TRUE), mode = "mirror")
  dest
}
for (mode in c("zero", "empty")) {
  f <- make_unevaluated_copy(ex("c3_aci_1.xlsx"), file.path(tmp, paste0("uncalc_", mode, ".xlsx")), mode)
  raw <- PhotoGEA::read_gasex_file(f, check_for_zero = c())
  a <- raw[, "A"]
  check(all(is.na(a) | a == 0) && any(raw[, "CO2_r"] > 0),
        paste0("precondition (", mode, "): A unreadable but raw CO2 present"))
  err <- tryCatch(read_one_licor(f), error = function(e) conditionMessage(e))
  check(is.character(err) && grepl("Calculate Now", err) && grepl("plaintext", err) &&
        grepl("A, gsw, Ci", err), paste0("unevaluated Excel (", mode, ") stops with the explanation"))
}

# 3. QC flags on hand-built rows -------------------------------------------------
bad <- data.frame(
  a_umol_m2_s     = c(20,   20,   20,   -1,    20,   20,   20,   20,   NA),
  ci_umol_mol     = c(250,  -5,   420,  420,   250,  250,  250,  250,  250),
  ca_umol_mol     = c(400,  400,  400,  400,   400,  400,  400,  400,  400),
  gsw_mol_m2_s    = c(0.3,  0.3,  0.3,  0.3,  -0.01, 0.3,  0.3,  0.3,  0.3),
  tleaf_c         = c(25,   25,   25,   25,    25,   52,   25,   25,   25),
  stable_n        = c(2,    2,    2,    2,     2,    2,    1,    2,    2),
  stable_total    = c(2,    2,    2,    2,     2,    2,    2,    2,    2),
  co2_r_umol_mol  = c(410,  410,  410,  410,   410,  410,  410,  400.4, 410),
  co2_s_umol_mol  = c(400,  400,  400,  400,   400,  400,  400,  400,  400))
q <- add_qc_flags(bad, tleaf_min = 5, tleaf_max = 45, min_delta_co2 = 1)
check(identical(q$qc_ci_negative, c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)), "Ci < 0 flagged")
check(identical(q$qc_ci_gt_ca, c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)),
      "Ci > Ca flagged only when A > 0 (row 4 has A < 0, legitimate)")
check(identical(q$qc_gsw_negative, c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)), "gsw < 0 flagged")
check(identical(q$qc_tleaf_out_of_range, c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)), "leaf temperature 52 C flagged")
check(identical(q$qc_unstable, c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)), "Stable < Total flagged")
check(identical(q$qc_small_delta_co2, c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)), "|dCO2| = 0.4 < 1 flagged")
check(identical(q$qc_missing_values, c(rep(FALSE, 8), TRUE)), "missing A flagged")
check(identical(q$qc_hard_fail, c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE)),
      "hard failures = missing, Ci<0, Ci>Ca with A>0, gsw<0, Tleaf out of range")
check(identical(q$qc_flags[c(1, 2, 7)], c("", "ci_negative", "unstable")), "readable flag text")
q2 <- add_qc_flags(bad[, c("a_umol_m2_s", "ci_umol_mol")])
check(all(is.na(q2$qc_gsw_negative)) && all(is.na(q2$qc_unstable)), "checks without their columns are NA, not FALSE")

# 4. Confounding of design factors with instrument --------------------------------
mk <- function(geno, inst) data.frame(curve = paste0("c", seq_along(geno)), genotype = geno,
                                      instrument_id = inst, stringsAsFactors = FALSE)
st <- function(d) confounding_check(d, "genotype", "curve")$summary
s_comp <- st(mk(c("A", "A", "A", "B", "B", "B"), c("M1", "M1", "M1", "M2", "M2", "M2")))
check(s_comp$status == "completely confounded" && s_comp$warn && isTRUE(all.equal(s_comp$cramers_v, 1)),
      "genotype A only on M1, B only on M2: completely confounded, V = 1")
s_part <- st(mk(c("A", "A", "B", "B", "C", "C"), c("M1", "M2", "M1", "M2", "M2", "M2")))
check(s_part$status == "partially confounded" && s_part$warn, "C never on M1: partially confounded")
s_bal <- st(mk(c("A", "A", "B", "B"), c("M1", "M2", "M1", "M2")))
check(s_bal$status == "crossed" && !s_bal$warn && isTRUE(all.equal(s_bal$cramers_v, 0)),
      "balanced design: crossed, V = 0, no warning")
s_one <- st(mk(c("A", "B"), c("M1", "M1")))
check(s_one$status == "single instrument" && !s_one$warn, "one instrument: nothing to confound")
check(grepl("blocking factor", s_part$message), "partial confounding suggests instrument as a blocking factor")

# 5. End to end on the demo data ------------------------------------------------
out <- file.path(tmp, "demo")
res <- run_cli(c("--demo", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "tables/licor_tidy.csv", "tables/file_summary.csv", "tables/qc_flag_counts.csv",
            "tables/curve_summary.csv", "tables/confounding_summary.csv", "figures/qc_overview.png",
            "reproducibility/commands.sh", "reproducibility/sessionInfo.txt")) {
  check(file.exists(file.path(out, f)), paste("demo writes", f))
}
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
tidy <- utils::read.csv(file.path(out, "tables", "licor_tidy.csv"), check.names = FALSE)
check(nrow(tidy) == 96L && all(c("a_umol_m2_s", "ci_umol_mol", "instrument_id", "qc_hard_fail") %in% names(tidy)),
      "tidy CSV has 96 rows and the contract columns")
rep_txt <- readLines(file.path(out, "report.md"))
check(any(grepl("Oxygen is not 21%", rep_txt)), "report warns that c3_aci_2.xlsx logged 22% oxygen")

bad_run <- run_cli(c("--demo", "--color", "red"))
check(bad_run$status != 0L, "unknown option exits with an error")
csv_run <- run_cli(c("--input", shQuote(file.path(out, "tables", "licor_tidy.csv")), "--out", shQuote(file.path(tmp, "x"))))
check(csv_run$status != 0L && any(grepl("CSV", csv_run$output)), "a CSV input is rejected with a pointer to the next skill")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
