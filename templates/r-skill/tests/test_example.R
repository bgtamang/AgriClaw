# Tests for claw-example. Run from anywhere: Rscript tests/test_example.R
# Template for AgriClaw R tests: known-answer checks, a bad-input check, and an
# end-to-end run of the command line on the demo data.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "example.R"))  # loads functions; main() does not run

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}

run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "example.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}

# 1. Known answers, computed by hand
d <- data.frame(genotype = c("A", "A", "A", "B", "B"),
                yield_kg_ha = c(1, 2, 3, 10, 20))
r <- summarise_trait(d, "yield_kg_ha", "genotype", min_n = 3)
check(isTRUE(all.equal(r$mean, c(2, 15))), "group means")
check(isTRUE(all.equal(r$sd, c(1, sqrt(50)))), "group SDs")
check(isTRUE(all.equal(r$cv_pct, c(50, 100 * sqrt(50) / 15))), "CV %")
check(identical(r$low_n, c(FALSE, TRUE)), "low-n flag")

# 2. Bad input fails loudly and names the problem
err <- tryCatch(summarise_trait(d, "biomass_g", "genotype"),
                error = function(e) conditionMessage(e))
check(is.character(err) && grepl("biomass_g", err), "missing column is named in the error")

# 3. End to end on the demo data
out <- file.path(tempdir(), "claw-example-test")
unlink(out, recursive = TRUE)
res <- run_cli(c("--demo", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "tables/summary.csv", "figures/means.png",
            "reproducibility/commands.sh", "reproducibility/sessionInfo.txt")) {
  check(file.exists(file.path(out, f)), paste("demo writes", f))
}
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L,
      "demo writes a checksum file")
s <- utils::read.csv(file.path(out, "tables", "summary.csv"))
check(nrow(s) == 4L && all(s$n == 3L), "demo summary has 4 groups of 3")

# 4. Unknown options are rejected
bad <- run_cli(c("--demo", "--colour", "red"))
check(bad$status != 0L, "unknown option exits with an error")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
