# Tests for claw-trial-design. Run from anywhere: Rscript tests/test_trial_design.R
# Design-property checks for every design (written independently of the
# script's own internal checks), known-answer efficiency factors, failure modes
# that the old Python version got wrong, and an end-to-end --demo run.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "trial_design.R"))  # loads functions; main() does not run

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
expect_error <- function(expr, pattern, what) {
  msg <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  check(!is.na(msg) && grepl(pattern, msg), paste(what, "| got:", msg))
}

run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "trial_design.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}

O <- function(...) utils::modifyList(DEFAULTS, list(...))
is_rect <- function(r, c) length(r) == (diff(range(r)) + 1) * (diff(range(c)) + 1) && !anyDuplicated(paste(r, c))
all_rect <- function(b, g) all(vapply(split(seq_len(nrow(b)), g), function(i) is_rect(b$row[i], b$col[i]), TRUE))
once_per_rep <- function(b, n) all(tapply(b$entry, b$rep, function(x) length(x) == n && !anyDuplicated(x)))
grid_ok <- function(b) !anyDuplicated(paste(b$row, b$col)) && nrow(b) == max(b$row) * max(b$col)
max_concurrence <- function(trt, blk) { N <- unclass(table(trt, blk)); M <- tcrossprod(N); max(M[upper.tri(M)]) }

# Entry list with awkward names, checks and an extra column
tmp <- tempfile("claw_td_")
dir.create(tmp)
el_path <- file.path(tmp, "entries.csv")
el <- data.frame(entry = c("Check A", "Check, B", sprintf("line-%02d/x", 1:22)),
                 is_check = c("yes", "yes", rep("no", 22)), origin = c("ref", "ref", rep("cross", 22)))
utils::write.csv(el, el_path, row.names = FALSE)

# 1. RCBD: complete, compact blocks; names and extra columns propagate; seeds reproduce
r1 <- build_design(O(design = "rcbd", entry_list = el_path, reps = 4, seed = 11))
b <- r1$book
check(nrow(b) == 96 && grid_ok(b), "rcbd fills its field grid")
check(once_per_rep(b, 24), "rcbd: every entry exactly once per replicate")
check(all_rect(b, b$rep), "rcbd: every replicate is a compact rectangle")
check(setequal(unique(b$entry), el$entry), "rcbd: entry names (with spaces, commas, slashes) propagate unchanged")
check(all(b$is_check[b$entry %in% c("Check A", "Check, B")]) && !any(b$is_check[!b$entry %in% c("Check A", "Check, B")]),
      "rcbd: check flags propagate")
check(all(b$origin[b$entry == "Check A"] == "ref"), "rcbd: extra entry-list columns propagate")
r1b <- build_design(O(design = "rcbd", entry_list = el_path, reps = 4, seed = 11))
check(identical(r1$book, r1b$book), "rcbd: same seed gives an identical field book")
r1c <- build_design(O(design = "rcbd", entry_list = el_path, reps = 4, seed = 12))
check(!identical(r1$book$entry, r1c$book$entry), "rcbd: a different seed gives a different randomisation")
check(r1$error_df == (4 - 1) * (24 - 1), "rcbd: skeleton error df = (r - 1)(t - 1)")

# 2. Field dimensions and block shape are honoured; impossible fields fail with options
r2 <- build_design(O(design = "rcbd", entries = 12, reps = 4, seed = 2, field_rows = 4, field_cols = 12))
check(max(r2$book$row) == 4 && max(r2$book$col) == 12 && all_rect(r2$book, r2$book$rep),
      "rcbd: requested 4 x 12 field honoured with rectangular replicates")
r2b <- build_design(O(design = "rcbd", entries = 12, reps = 3, seed = 2, block_rows = 1))
check(all(tapply(r2b$book$row, r2b$book$rep, function(x) length(unique(x)) == 1)), "--block-rows 1 puts each replicate in one row")
expect_error(build_design(O(design = "rcbd", entries = 12, reps = 4, seed = 2, field_rows = 5, field_cols = 10)),
             "Field shapes that work: .*4 x 12", "infeasible field size fails and lists feasible shapes")

# 3. Serpentine and cartesian numbering
sp <- r2$book[order(r2$book$plot), ]
check(all(abs(diff(sp$row)) + abs(diff(sp$col)) == 1), "serpentine: consecutive plots are neighbours")
r3 <- build_design(O(design = "rcbd", entries = 12, reps = 4, seed = 2, field_rows = 4, field_cols = 12, numbering = "cartesian", plot_start = 101))
cp <- r3$book[order(r3$book$plot), ]
check(identical(cp$plot, 101:148) && all(diff(cp$row * 100 + cp$col) > 0), "cartesian: row-major numbering from --plot-start")
check(identical(r3$book[order(r3$book$row, r3$book$col), "entry"], r2$book[order(r2$book$row, r2$book$col), "entry"]),
      "numbering changes labels only, never which entry sits where")

# 4. Latin square
for (q in list(c(5, 1), c(4, 2))) {
  r4 <- build_design(O(design = "latin", entries = q[1], reps = q[2], seed = 3))
  b <- r4$book
  ok <- all(vapply(split(b, b$rep), function(s) all(table(s$entry, s$row) == 1) && all(table(s$entry, s$col) == 1), TRUE))
  check(ok && all_rect(b, b$rep), sprintf("latin %dx%d (%d square[s]): each treatment once per row and once per column", q[1], q[1], q[2]))
}

# 5. Efficiency factor: known answers for square lattices
sl <- with_agricolae(agricolae::design.alpha(1:25, k = 5, r = 2, seed = 5))$book
check(abs(efficiency_factor(sl[[4]], sl$block) - 6 / 8) < 1e-10, "simple lattice (k = 5): E = (k+1)/(k+3) = 0.75")
tl <- with_agricolae(agricolae::design.alpha(1:25, k = 5, r = 3, seed = 5))$book
check(abs(efficiency_factor(tl[[4]], tl$block) - 12 / 15) < 1e-10, "triple lattice (k = 5): E = 2(k+1)/(2k+5) = 0.8")
check(abs(alpha_upper_bound(25, 2, 5) - 0.75) < 1e-12, "upper bound equals the simple-lattice E when k = s")
check(efficiency_factor(c(1, 2, 3, 4), c(1, 1, 2, 2)) == 0, "a disconnected design has E = 0")

# 6. Alpha designs
ra <- build_design(O(design = "alpha", entries = 24, reps = 3, block_size = 4, seed = 4, alpha_engine = "agricolae"))
b <- ra$book
check(ra$gen$info$package_function == "agricolae::design.alpha", "alpha: agricolae engine used when requested")
check(once_per_rep(b, 24), "alpha: resolvable, every entry once per replicate")
check(all(table(b$block) == 4) && length(unique(b$block)) == 18, "alpha: 18 blocks of size 4")
check(all(tapply(b$rep, b$block, function(x) length(unique(x)) == 1)), "alpha: blocks nested in replicates")
check(all_rect(b, b$block) && all_rect(b, b$rep), "alpha: blocks and replicates are compact rectangles")
check(max_concurrence(b$entry, b$block) <= 1, "alpha(0,1): no pair of entries meets in more than one block")
check(ra$gen$info$efficiency_factor <= ra$gen$info$efficiency_upper_bound + 1e-9, "alpha: E does not exceed the upper bound")
check(ra$error_df == 72 - 1 - 2 - 15 - 23, "alpha: skeleton error df = N - 1 - (r-1) - r(s-1) - (v-1)")
rf <- build_design(O(design = "alpha", entries = 24, reps = 5, block_size = 4, seed = 4))
b <- rf$book
check(rf$gen$info$package_function == "FielDHub::alpha_lattice", "alpha: FielDHub used where agricolae has no generator (r = 5)")
check(once_per_rep(b, 24) && all(table(b$block) == 4) && all_rect(b, b$block), "alpha (FielDHub): resolvable with compact blocks of 4")
check(rf$gen$info$efficiency_factor <= rf$gen$info$efficiency_upper_bound + 1e-9, "alpha (FielDHub): E within the upper bound")
expect_error(build_design(O(design = "alpha", entries = 22, reps = 2, block_size = 4, seed = 1)),
             "cannot be split into blocks of 4 without dropping entries.*2, 11", "alpha: entries never dropped silently when k does not divide v")
expect_error(build_design(O(design = "alpha", entries = 23, reps = 2, seed = 1)), "prime", "alpha: prime entry count is refused with advice")
rauto <- build_design(O(design = "alpha", entries = 36, reps = 2, seed = 8))
check(rauto$gen$info$block_size == 6, "alpha: default block size is the divisor nearest sqrt(v)")

# 7. Augmented RCBD (agricolae::design.dau)
rg <- build_design(O(design = "augmented", entries = 23, checks = 4, blocks = 5, seed = 6))
b <- rg$book
chk <- paste0("CHK", 1:4)
check(all(vapply(split(b$entry, b$block), function(x) all(table(factor(x[x %in% chk], levels = chk)) == 1), TRUE)),
      "augmented: every block holds every check exactly once")
check(all(table(b$entry[grepl("^G", b$entry)]) == 1) && sum(grepl("^G", b$entry)) == 23, "augmented: every new entry exactly once, none dropped")
check(all_rect(b, b$block), "augmented: blocks are compact rectangles")
check(sum(b$is_filler) == 2 && all(b$entry[b$is_filler] == "FILLER"), "augmented: uneven blocks padded with flagged FILLER plots")
set.seed(1)
bb <- b[!b$is_filler, ]; bb$y <- rnorm(nrow(bb))
check(stats::lm(y ~ factor(block) + entry, data = bb)$df.residual == (5 - 1) * (4 - 1) && rg$gen$info$error_df_checks == 12,
      "augmented: error df = (b - 1)(c - 1), confirmed by lm")
expect_error(build_design(O(design = "augmented", entries = 20, blocks = 4, seed = 1)), "need check entries", "augmented without checks is refused")

# 8. Split-plot: whole plots contiguous and complete; strata df from aov
rs <- build_design(O(design = "split", factor_a = "irrigation:dry,wet,full", factor_b = "variety:V1,V2,V3,V4", reps = 3, seed = 9))
b <- rs$book
check(all(vapply(split(b, b$whole_plot), function(s) length(unique(s$irrigation)) == 1 && setequal(s$variety, paste0("V", 1:4)) && nrow(s) == 4, TRUE)),
      "split: each whole plot = one irrigation level with all four varieties")
check(all_rect(b, b$whole_plot) && all_rect(b, b$rep), "split: whole plots and replicates are contiguous rectangles")
check(all(table(b$rep, b$irrigation) == 4), "split: each irrigation level once per replicate (4 sub-plots)")
set.seed(2); b$y <- rnorm(nrow(b))
sm <- summary(stats::aov(y ~ irrigation * variety + Error(factor(rep) / irrigation), data = b))
dfs <- vapply(sm, function(s) { s <- s[[1]]; s["Residuals", "Df"] }, 0)
check(dfs[["Error: factor(rep):irrigation"]] == rs$gen$info$whole_plot_error_df && dfs[["Error: Within"]] == rs$gen$info$sub_plot_error_df,
      "split: whole-plot (4) and sub-plot (18) error df match aov strata")

# 9. Strip-plot: strips stay intact (the old serpentine bug)
rt <- build_design(O(design = "strip", factor_a = "tillage:conv,notill,ridge", factor_b = "N:0,60,120,180", reps = 3, seed = 10))
b <- rt$book
check(all(vapply(split(b, b$rep), function(s) all(tapply(s$tillage, s$row, function(x) length(unique(x)) == 1)) &&
                   all(tapply(s$n, s$col, function(x) length(unique(x)) == 1)), TRUE)),
      "strip: tillage constant along each row and N constant down each column within replicates")
check(once_per_rep(b, 12) && all_rect(b, b$rep), "strip: all 12 combinations once per compact replicate")

# 10. Factorial in RCBD
rx <- build_design(O(design = "factorial", factor_a = "N:0,60,120", factor_b = "P:0,40", factor_c = "K:0,1", reps = 2, seed = 11))
b <- rx$book
check(once_per_rep(b, 12) && all(table(b$n, b$p, b$k) == 2), "factorial: every N x P x K combination once per replicate")

# 11. CRD
rc <- build_design(O(design = "crd", entries = 6, reps = 4, seed = 1))
check(all(table(rc$book$entry) == 4) && grid_ok(rc$book), "crd: each entry 4 times; field filled")

# 12. Row-column (FielDHub::row_column)
rr <- build_design(O(design = "rowcol", entries = 20, reps = 2, rep_rows = 4, iterations = 50, seed = 5))
b <- rr$book
check(once_per_rep(b, 20) && all_rect(b, b$rep), "rowcol: resolvable with compact replicates")
check(all(vapply(split(b, b$rep), function(s) all(table(s$row) == 5) && all(table(s$col) == 4), TRUE)), "rowcol: 4 complete rows x 5 complete columns per replicate")

# 13. P-rep (FielDHub::partially_replicated)
rp <- build_design(O(design = "prep", entries = 60, checks = 2, check_reps = 4, seed = 7))
b <- rp$book
tb <- table(b$entry)
ntwo <- sum(tb[grepl("^G", names(tb))] == 2)
check(all(tb[c("CHK1", "CHK2")] == 4) && all(tb[grepl("^G", names(tb))] %in% 1:2), "prep: checks x4, test entries x1 or x2")
check(ntwo >= 12 && ntwo <= 18 && ntwo == rp$gen$info$replicated_test_entries, "prep: about 25% of test entries replicated (13-17 of 60)")
check(grid_ok(b) && nrow(b) == 60 + ntwo + 8, "prep: plots fill the field exactly")
rp2 <- build_design(O(design = "prep", entries = 100, checks = 3, check_reps = 4, field_rows = 12, field_cols = 11, seed = 3))
check(max(rp2$book$row) == 12 && max(rp2$book$col) == 11 && sum(table(rp2$book$entry) == 2) == 132 - 112,
      "prep: given field size is filled by adjusting the number of replicated entries")

# 14. Diagonal (FielDHub::diagonal_arrangement)
rd <- build_design(O(design = "diagonal", entries = 100, checks = 4, field_rows = 10, field_cols = 12, seed = 8))
b <- rd$book
check(all(table(b$entry[grepl("^G", b$entry)]) == 1) && sum(grepl("^G", b$entry)) == 100, "diagonal: every line exactly once")
check(sum(b$is_check) + 100 + sum(b$is_filler) == 120 && all(paste0("CHK", 1:4) %in% b$entry), "diagonal: checks fill the remaining plots")
expect_error(build_design(O(design = "diagonal", entries = 100, checks = 4, seed = 8)), "field-rows", "diagonal without field size is refused")

# 15. Input validation
dup <- file.path(tmp, "dup.csv")
utils::write.csv(data.frame(entry = c("A", "B", "A")), dup, row.names = FALSE)
expect_error(build_design(O(design = "rcbd", entry_list = dup, reps = 2, seed = 1)), "unique.*A", "duplicated entry names are refused")
expect_error(build_design(O(design = "rcbd", entries = 10, reps = 1, seed = 1)), "reps", "RCBD with one replicate is refused")
expect_error(build_design(O(design = "magic", entries = 10, seed = 1)), "Unknown --design", "unknown design is refused")

# 16. End to end on the demo data
out <- file.path(tmp, "demo_out")
res <- run_cli(c("--demo", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "tables/field_book.csv", "tables/design_summary.csv", "tables/design_checks.csv",
            "figures/field_map.png", "reproducibility/commands.sh", "reproducibility/sessionInfo.txt")) {
  check(file.exists(file.path(out, f)), paste("demo writes", f))
}
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
fb <- utils::read.csv(file.path(out, "tables", "field_book.csv"), stringsAsFactors = FALSE)
summ <- utils::read.csv(file.path(out, "tables", "design_summary.csv"), stringsAsFactors = FALSE)
check(nrow(fb) == 90 && all(table(fb$entry) == 3) && "pedigree" %in% names(fb), "demo: 30 entries x 3 replicates with pedigree carried over")
check(summ$value[summ$parameter == "seed"] == "2026", "demo: seed recorded in the design summary")
check(all(utils::read.csv(file.path(out, "tables", "design_checks.csv"))$passed), "demo: all design checks pass")
nos <- run_cli(c("--design", "rcbd", "--entries", "6", "--reps", "2", "--out", shQuote(file.path(tmp, "noseed"))))
cmd <- readLines(file.path(tmp, "noseed", "reproducibility", "commands.sh"))
check(nos$status == 0L && any(grepl("--seed", cmd)), "a random seed is drawn and written into commands.sh")
bad <- run_cli(c("--design", "rcbd", "--colour", "red"))
check(bad$status != 0L, "unknown option exits with an error")

unlink(tmp, recursive = TRUE)
cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
