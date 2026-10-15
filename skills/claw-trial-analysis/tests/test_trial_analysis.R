# Tests for claw-trial-analysis. Run from anywhere: Rscript tests/test_trial_analysis.R
# Known answers: the Cullis H2 implementation against the closed form for
# balanced RCBD; simulated fields (100 genotypes x 2 replicates, known genetic
# and residual variances, with and without a smooth spatial trend) where BLUEs
# must track the true genetic values, better under the spatial model, and H2
# must land near the truth on average; a planted outlier; input failures; a
# smoke test on statgenSTA's dropsRaw; and an end-to-end --demo run.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "trial_analysis.R"))  # loads functions; main() does not run

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
  script <- file.path(SKILL_TEST_DIR, "scripts", "trial_analysis.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
O <- function(...) utils::modifyList(DEFAULTS, list(...))
quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# Simulated field: 100 genotypes x 2 replicates (RCBD), each replicate a 10 x 10
# block, stacked to a 20-row x 10-column field. sigma2_g = sigma2_e = 1, so the
# true entry-mean H2 = 1 / (1 + 1/2) = 0.667.
simulate_field <- function(seed, amp = 1.5, ng = 100) {
  set.seed(seed)
  g <- stats::setNames(stats::rnorm(ng), sprintf("G%03d", seq_len(ng)))
  d <- do.call(rbind, lapply(1:2, function(k) data.frame(rep = k, entry = sample(names(g)),
                                                         row = rep(1:10, each = 10) + 10 * (k - 1), col = rep(1:10, 10),
                                                         stringsAsFactors = FALSE)))
  trend <- amp * (sin(d$row / 3.5) + cos(d$col / 2.5) + 0.05 * d$row)
  d$y <- 10 + g[d$entry] + trend + stats::rnorm(nrow(d))
  d$trial <- paste0("S", seed)
  d$plot_id <- paste0("P", seq_len(nrow(d)))
  list(d = d, g = g)
}
TRUE_H2 <- 1 / (1 + 1 / 2)

# 1. Cullis H2 against the closed form for balanced RCBD
s0 <- simulate_field(101, amp = 0)
m0 <- quiet(lme4::lmer(y ~ factor(rep) + (1 | entry), data = s0$d))
vc <- as.data.frame(lme4::VarCorr(m0))
closed <- vc$vcov[vc$grp == "entry"] / (vc$vcov[vc$grp == "entry"] + vc$vcov[vc$grp == "Residual"] / 2)
check(abs(cullis_h2(m0, "entry") - closed) < 1e-8, sprintf("Cullis H2 = sigma2_g / (sigma2_g + sigma2_e / r) for balanced RCBD (%.4f)", closed))
s0b <- s0$d; s0b$y[1:3] <- NA
m0b <- quiet(lme4::lmer(y ~ factor(rep) + (1 | entry), data = s0b))
h_unbal <- cullis_h2(m0b, "entry")
check(h_unbal > 0 && h_unbal < 1, "Cullis H2 stays in (0, 1) with missing plots")

# 2. Simulated fields with a smooth trend: spatial model recovers genotypes better.
# Seed 6 also carries a planted +6 SD outlier in plot P57.
seeds <- 1:6
outlier_run <- NULL
res_t <- lapply(seeds, function(sd) {
  s <- simulate_field(sd)
  if (sd == 6) s$d$y[s$d$plot_id == "P57"] <- s$d$y[s$d$plot_id == "P57"] + 6
  r <- quiet(analyze(s$d, O(traits = "y", design = "rcbd", engine = "both")))
  if (sd == 6) outlier_run <<- r
  bs <- r$blues[r$blues$model == "SpATS", ]; bl <- r$blues[r$blues$model == "lme4", ]
  sm <- r$summary
  c(cor_s = stats::cor(bs$blue, s$g[bs$genotype]), cor_l = stats::cor(bl$blue, s$g[bl$genotype]),
    h2_s = sm$h2[sm$model == "SpATS"], h2_l = sm$h2[sm$model == "lme4"], ve_s = sm$var_residual[sm$model == "SpATS"])
})
rt <- do.call(rbind, res_t)
cm <- colMeans(rt)
cat(sprintf("  trend: mean cor SpATS %.3f, lme4 %.3f; mean H2 SpATS %.3f, lme4 %.3f (truth %.3f); mean residual var SpATS %.3f (truth 1)
",
            cm["cor_s"], cm["cor_l"], cm["h2_s"], cm["h2_l"], TRUE_H2, cm["ve_s"]))
check(cm["cor_s"] >= 0.75, "trend: spatial BLUEs correlate with true genetic values (mean r >= 0.75 over 6 seeds)")
check(cm["cor_s"] - cm["cor_l"] >= 0.05 && sum(rt[, "cor_s"] > rt[, "cor_l"]) >= 5,
      "trend: spatial BLUEs beat the RCBD model (mean gain >= 0.05, better in >= 5 of 6 seeds)")
check(abs(cm["h2_s"] - TRUE_H2) <= 0.1, "trend: mean SpATS H2 within 0.1 of the true 0.667")
check(cm["h2_s"] > cm["h2_l"], "trend: ignoring the trend lowers H2 (lme4 below SpATS on average)")
check(abs(cm["ve_s"] - 1) <= 0.25, "trend: SpATS residual variance within 0.25 of the true 1")

# 3. Same seeds without a trend: the two models agree and both recover H2
res_n <- do.call(rbind, lapply(seeds, function(sd) {
  s <- simulate_field(sd, amp = 0)
  r <- quiet(analyze(s$d, O(traits = "y", design = "rcbd", engine = "both")))
  sm <- r$summary
  c(h2_s = sm$h2[sm$model == "SpATS"], h2_l = sm$h2[sm$model == "lme4"])
}))
cn <- colMeans(res_n)
cat(sprintf("  no trend: mean H2 SpATS %.3f, lme4 %.3f (truth %.3f)
", cn["h2_s"], cn["h2_l"], TRUE_H2))
check(abs(cn["h2_l"] - TRUE_H2) <= 0.1, "no trend: mean lme4 H2 within 0.1 of the true 0.667")
check(abs(cn["h2_s"] - cn["h2_l"]) <= 0.05, "no trend: spatial and non-spatial H2 agree within 0.05")

# 4. The planted outlier is found and traced to its plot; outputs are consistent
r <- outlier_run
o <- r$outliers
check(!is.null(o) && any(o$plot_id == "P57" & o$model == "SpATS"), "planted +6 SD outlier flagged by outlierSTA on the spatial model")
check(all(c("trial", "genotype", "y", "se_y", "wt_y") %in% names(r$met)) && nrow(r$met) == 100,
      "met_input: one row per genotype with BLUE, se and weight")
check(isTRUE(all.equal(r$met$wt_y, 1 / r$met$se_y^2)), "met_input weights = 1 / se^2")
check(r$info[[1]]$blues_from == "SpATS", "BLUEs exported from the spatial model by default")
check(sum(r$blups$model == "SpATS") == 100 && sum(r$blups$model == "lme4") == 100, "BLUPs for every genotype from both models")
check(!is.null(r$comparison) && nrow(r$comparison) == 1 && r$comparison$cor_blues > 0.5 && r$comparison$more_precise == "spatial (SpATS)",
      sprintf("comparison table written; spatial model more precise under a trend (BLUE r = %.2f)", r$comparison$cor_blues))
s <- simulate_field(7)

# 5. Input checks fail loudly
bad <- s$d; bad$row[5] <- NA
expect_error(analyze(bad, O(traits = "y", design = "rcbd")), "Coordinates are incomplete", "incomplete coordinates are refused for spatial models")
rl <- quiet(analyze(bad, O(traits = "y", design = "rcbd", engine = "lme4")))
check(identical(rl$info[[1]]$engines, "lme4") && any(grepl("Coordinates are incomplete", rl$warnings)) && nrow(rl$met) == 100,
      "lme4 block model runs with incomplete coordinates, which are ignored with a warning")
dup <- s$d; dup$row[2] <- dup$row[1]; dup$col[2] <- dup$col[1]
expect_error(analyze(dup, O(traits = "y", design = "rcbd")), "share a row/column position", "duplicated coordinates are refused")
unrep <- s$d[s$d$rep == 1, c("entry", "y", "trial", "plot_id")]
expect_error(analyze(unrep, O(traits = "y")), "No replicate, block or coordinate", "data with no design structure are refused")
unrep$block <- rep(1:10, each = 10)
expect_error(analyze(unrep, O(traits = "y")), "unreplicated.*need replicated checks or a spatial model",
             "unreplicated trial without checks or coordinates is refused")
txt <- s$d; txt$y <- as.character(txt$y); txt$y[3] <- "lodged"
expect_error(analyze(txt, O(traits = "y")), "not numeric", "non-numeric trait is refused")
expect_error(analyze(s$d, O(traits = "biomass_g")), "biomass_g", "missing trait is named in the error")

# 6. Augmented-style trials: unreplicated lines + 3 checks in each of 6 blocks
set.seed(5)
chk <- c("C1", "C2", "C3")
aug <- do.call(rbind, lapply(1:6, function(b) data.frame(block = b, entry = sample(c(chk, sprintf("L%02d", (b - 1) * 8 + 1:8))),
                                                         row = b, col = 1:11, stringsAsFactors = FALSE)))
aug$is_check <- aug$entry %in% chk
geff <- stats::setNames(stats::rnorm(length(unique(aug$entry)), 0, 1), unique(aug$entry))
aug$y <- 20 + geff[aug$entry] + rep(stats::rnorm(6, 0, 0.5), each = 11) + stats::rnorm(nrow(aug), 0, 0.5)
ra <- quiet(analyze(aug[, c("block", "entry", "is_check", "y")], O(traits = "y")))
check(ra$info[[1]]$design == "ibd" && identical(ra$info[[1]]$engines, "lme4") && !length(ra$info[[1]]$checkid_engines),
      "augmented, no coordinates: blocks model (ibd) with lme4; checks stay random genotypes in the BLUP model")
check(nrow(ra$met) == 51 && all(chk %in% ra$met$genotype), "augmented, no coordinates: BLUEs for all 51 entries")
rs <- quiet(analyze(aug, O(traits = "y", engine = "SpATS")))
check(identical(rs$info[[1]]$checkid_engines, "SpATS") && !any(chk %in% rs$blups$genotype) && sum(rs$blups$model == "SpATS") == 48,
      "augmented with coordinates: SpATS fits checks as fixed checkId effects; BLUPs for the 48 test lines only")
check(rs$summary$status == "ok" && rs$summary$h2 > 0, "augmented with coordinates: spatial model fits and gives H2 > 0")

# 7. Smoke test on statgenSTA's dropsRaw (one DROPS trial, spatial model only)
data("dropsRaw", package = "statgenSTA", envir = environment())
dr <- dropsRaw[dropsRaw$Experiment == "Gai12W", ]
rd <- quiet(analyze(dr, O(traits = "grain.yield", genotype = "Variety_ID", trial = "Experiment", rep = "Replicate",
                          block = "block", row = "Row", col = "Column", engine = "SpATS")))
sm <- rd$summary
check(rd$info[[1]]$design == "res.ibd", "dropsRaw: blocks nested in replicates detected (res.ibd)")
check(nrow(rd$met) == length(unique(dr$Variety_ID[!is.na(dr$grain.yield)])), "dropsRaw: one BLUE per genotype with data")
check(sm$h2 > 0 && sm$h2 < 1 && sm$cv_pct > 0 && sm$var_genetic > 0, sprintf("dropsRaw: sane H2 (%.2f), CV (%.1f%%), genetic variance", sm$h2, sm$cv_pct))
check(abs(mean(rd$met$grain.yield) - mean(dr$grain.yield, na.rm = TRUE)) < 0.5, "dropsRaw: BLUE mean close to the raw mean (t/ha)")

# 8. End to end on the demo data
tmp <- tempfile("claw_ta_")
res <- run_cli(c("--demo", "--traits", "yield_kg_ha", "--out", shQuote(tmp)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "tables/blues.csv", "tables/blups.csv", "tables/model_summary.csv", "tables/model_comparison.csv",
            "tables/outliers.csv", "tables/met_input.csv", "figures/spatial_DEMO-AI_yield_kg_ha.png",
            "reproducibility/commands.sh", "reproducibility/sessionInfo.txt")) {
  check(file.exists(file.path(tmp, f)), paste("demo writes", f))
}
dout <- utils::read.csv(file.path(tmp, "tables", "outliers.csv"), stringsAsFactors = FALSE)
check(any(dout$plot_id == "DEMO-AI-57"), "demo: the planted outlier plot DEMO-AI-57 is reported")
tv <- utils::read.csv(file.path(SKILL_TEST_DIR, "demo", "true_genetic_values.csv"), stringsAsFactors = FALSE)
met <- utils::read.csv(file.path(tmp, "tables", "met_input.csv"), stringsAsFactors = FALSE)
check(stats::cor(met$yield_kg_ha, tv$true_yield_effect[match(met$genotype, tv$entry)]) > 0.8, "demo: BLUEs track the true genetic values (r > 0.8)")
bad <- run_cli(c("--input", shQuote(file.path(SKILL_TEST_DIR, "demo", "trial_clean.csv"))))
check(bad$status != 0L, "missing --traits exits with an error")
unlink(tmp, recursive = TRUE)

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
