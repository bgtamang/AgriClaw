# Run every R skill test (skills/*/tests/test_*.R), each in a fresh R process.
#
# Usage:
#   Rscript tests/run_r_tests.R                  # all skills
#   Rscript tests/run_r_tests.R claw-aci-fit     # skills whose name contains a pattern
#   Rscript tests/run_r_tests.R --templates      # also test templates/r-skill
#
# Exits with status 1 if any test file fails.

args <- commandArgs(trailingOnly = TRUE)
this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(this), ".."), winslash = "/")

files <- Sys.glob(file.path(root, "skills", "*", "tests", "test_*.R"))
if ("--templates" %in% args)
  files <- c(files, Sys.glob(file.path(root, "templates", "*", "tests", "test_*.R")))
patterns <- args[!startsWith(args, "--")]
if (length(patterns)) {
  keep <- Reduce(`|`, lapply(patterns, function(p) grepl(p, files, fixed = TRUE)))
  files <- files[keep]
}
if (!length(files)) {
  cat("No R test files found.\n")
  quit(status = 0)
}

rscript <- file.path(R.home("bin"), "Rscript")
res <- data.frame(test = character(), ok = logical(), seconds = numeric(),
                  summary = character(), stringsAsFactors = FALSE)
for (f in files) {
  label <- sub(paste0(root, "/"), "", normalizePath(f, winslash = "/"), fixed = TRUE)
  t0 <- Sys.time()
  out <- suppressWarnings(system2(rscript, shQuote(f), stdout = TRUE, stderr = TRUE))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  status <- attr(out, "status")
  ok <- is.null(status) || status == 0L
  pass_line <- grep("^PASS:", out, value = TRUE)
  summ <- if (ok && length(pass_line)) pass_line[length(pass_line)] else
    paste(utils::tail(out[nzchar(out)], 6), collapse = "\n    ")
  res[nrow(res) + 1L, ] <- list(label, ok, round(secs, 1), summ)
  cat(sprintf("%-4s %-60s %6.1fs  %s\n", if (ok) "ok" else "FAIL", label, secs,
              if (ok) summ else ""))
  if (!ok) cat("    ", summ, "\n", sep = "")
}
cat(sprintf("\n%d R test file(s): %d passed, %d failed.\n",
            nrow(res), sum(res$ok), sum(!res$ok)))
quit(status = if (all(res$ok)) 0L else 1L)
