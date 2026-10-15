# Tests for claw-met-gxe. Run from anywhere: Rscript tests/test_met_gxe.R
# Known answers from simulated METs (Finlay-Wilkinson slopes, variance
# components with and without stage-1 weights, nesting and location x year
# models, a rank-1 AMMI interaction, two mega-environments), Shukla's
# stability variance against agricolae, missing-cell handling, a smoke test
# on statgenGxE's dropsPheno, bad inputs, and an end-to-end --demo run.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "met_gxe.R"))  # loads functions; main() does not run

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
  script <- file.path(SKILL_TEST_DIR, "scripts", "met_gxe.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
O <- function(...) utils::modifyList(DEFAULTS, list(...))
grid_met <- function(ng, ne) expand.grid(genotype = sprintf("G%03d", seq_len(ng)), env = sprintf("E%02d", seq_len(ne)),
                                         stringsAsFactors = FALSE)
gid <- function(n) sprintf("G%03d", seq_len(n))

# 1. Finlay-Wilkinson recovers known sensitivities (y = mu + g + b * e + eps)
for (s in 1:2) {
  set.seed(s)
  d <- grid_met(80, 8); gi <- match(d$genotype, gid(80)); ej <- match(d$env, sprintf("E%02d", 1:8))
  g <- rnorm(80); b <- rnorm(80, 1, 0.25); e <- seq(-3, 3, length.out = 8)
  d$y <- 10 + g[gi] + b[gi] * e[ej] + rnorm(nrow(d), 0, 0.3); d$se2 <- NA
  est <- run_fw(make_td(d, "y"), "y", character())$estimates
  r <- cor(est$Sens, b[match(est$Genotype, gid(80))])
  check(r > 0.9, sprintf("FW seed %d: cor(estimated, true sensitivity) = %.3f > 0.9", s, r))
  check(abs(mean(est$Sens) - 1) < 1e-3, "FW: sensitivities average 1 for a complete table")
}

# 2. Variance components and heritability, unweighted (means: residual = GxE + error)
vcr <- t(vapply(1:3, function(s) {
  set.seed(100 + s)
  d <- grid_met(100, 8); gi <- match(d$genotype, gid(100)); ej <- match(d$env, sprintf("E%02d", 1:8))
  d$y <- 5 + rnorm(8, 0, 2)[ej] + rnorm(100)[gi] + rnorm(nrow(d), 0, sqrt(0.5)) + rnorm(nrow(d), 0, sqrt(0.25)); d$se2 <- NA
  v <- run_varcomp(make_td(d, "y"), d, "y")
  c(s2g = v$components$variance[v$components$component == "genotype"],
    res = v$components$variance[v$components$component == "residual"], h2 = v$h2)
}, numeric(3)))
m <- colMeans(vcr)
check(abs(m[["s2g"]] - 1) < 0.2, sprintf("gxeVarComp: mean genotype variance %.3f within 0.2 of 1", m[["s2g"]]))
check(abs(m[["res"]] - 0.75) < 0.08, sprintf("gxeVarComp: mean residual %.3f within 0.08 of 0.75 (GxE 0.5 + error 0.25)", m[["res"]]))
check(abs(m[["h2"]] - 1 / (1 + 0.75 / 8)) < 0.03, sprintf("herit: mean H2 %.3f within 0.03 of the entry-mean truth 0.914", m[["h2"]]))

# 3. Weighted stage 2 with known stage-1 variances (nlme) separates GxE from error
se2 <- c(0.05, 0.05, 0.1, 0.1, 0.3, 0.3, 1, 2)
wr <- t(vapply(1:3, function(s) {
  set.seed(200 + s)
  d <- grid_met(80, 8); gi <- match(d$genotype, gid(80)); ej <- match(d$env, sprintf("E%02d", 1:8))
  g <- rnorm(80)
  d$y <- 5 + rnorm(8, 0, 2)[ej] + g[gi] + rnorm(nrow(d), 0, sqrt(0.5)) + rnorm(nrow(d), 0, sqrt(se2[ej])); d$se2 <- se2[ej]
  w <- run_weighted_stage2(d); v <- run_varcomp(make_td(d, "y"), d, "y")
  c(s2g = w$sigma2_g, s2ge = w$sigma2_ge,
    accw = cor(w$blup$blup_weighted, g[match(w$blup$genotype, gid(80))]),
    accu = cor(v$blup$predictedValue, g[match(as.character(v$blup$genotype), gid(80))]))
}, numeric(4)))
m <- colMeans(wr)
check(abs(m[["s2ge"]] - 0.5) < 0.1, sprintf("weighted stage 2: mean GxE variance %.3f within 0.1 of 0.5", m[["s2ge"]]))
check(abs(m[["s2g"]] - 1) < 0.2, sprintf("weighted stage 2: mean genotype variance %.3f within 0.2 of 1", m[["s2g"]]))
check(m[["accw"]] >= m[["accu"]], sprintf("weighted BLUPs at least as accurate as unweighted (%.3f vs %.3f) with unequal trial precision", m[["accw"]], m[["accu"]]))

# 4. Nested trials (scenario) and location x year models
nr <- t(vapply(1:3, function(s) {
  set.seed(300 + s)
  d <- grid_met(80, 9); gi <- match(d$genotype, gid(80)); ej <- match(d$env, sprintf("E%02d", 1:9)); sc <- rep(1:3, each = 3)
  gs <- matrix(rnorm(80 * 3, 0, sqrt(0.4)), 80, 3)
  d$y <- 5 + rnorm(9, 0, 2)[ej] + rnorm(80)[gi] + gs[cbind(gi, sc[ej])] + rnorm(nrow(d), 0, sqrt(0.5)); d$se2 <- NA
  d$nesting <- paste0("S", sc[ej])
  v <- run_varcomp(make_td(d, "y"), d, "y")
  setNames(v$components$variance, v$components$component)[c("genotype", "genotype:nesting", "residual")]
}, numeric(3)))
m <- colMeans(nr)
check(abs(m[[1]] - 1) < 0.25 && abs(m[[2]] - 0.4) < 0.15 && abs(m[[3]] - 0.5) < 0.08,
      sprintf("nesting model recovers G 1, G x scenario 0.4, residual 0.5 (got %.2f, %.2f, %.2f)", m[[1]], m[[2]], m[[3]]))
lr <- t(vapply(1:3, function(s) {
  set.seed(400 + s)
  tr <- expand.grid(l = 1:4, yy = 1:3)
  d <- grid_met(80, 12); gi <- match(d$genotype, gid(80)); ej <- match(d$env, sprintf("E%02d", 1:12))
  gl <- matrix(rnorm(80 * 4, 0, sqrt(0.4)), 80, 4); gy <- matrix(rnorm(80 * 3, 0, sqrt(0.3)), 80, 3)
  d$y <- 5 + rnorm(12, 0, 2)[ej] + rnorm(80)[gi] + gl[cbind(gi, tr$l[ej])] + gy[cbind(gi, tr$yy[ej])] + rnorm(nrow(d), 0, sqrt(0.5)); d$se2 <- NA
  d$location <- paste0("L", tr$l[ej]); d$year <- paste0("Y", tr$yy[ej])
  v <- run_varcomp(make_td(d, "y"), d, "y")
  setNames(v$components$variance, v$components$component)[c("genotype", "genotype:loc", "genotype:year", "residual")]
}, numeric(4)))
m <- colMeans(lr)
check(abs(m[[1]] - 1) < 0.3 && abs(m[[2]] - 0.4) < 0.12 && abs(m[[3]] - 0.3) < 0.12 && abs(m[[4]] - 0.5) < 0.06,
      sprintf("location x year model recovers G 1, GL 0.4, GY 0.3, residual 0.5 (got %.2f, %.2f, %.2f, %.2f)", m[[1]], m[[2]], m[[3]], m[[4]]))

# 5. Stability measures: Shukla against agricolae, formula properties
stb <- utils::read.csv(system.file("external/dataStb.csv", package = "agricolae"))
ag <- suppressWarnings(agricolae::stability.par(stb, rep = 4, MSerror = 1.8))$statistics
ds <- data.frame(genotype = rep(rownames(stb), ncol(stb)), env = rep(names(stb), each = nrow(stb)),
                 y = unlist(stb), se2 = NA, stringsAsFactors = FALSE)
st <- run_stability(make_td(ds, "y"), ds, "y", "max", completeness(ds))
i <- match(rownames(ag), st$genotype)
check(max(abs(4 * st$wricke_ecovalence[i] - ag$Ecovalence)) < 1e-6, "Wricke's ecovalence (statgenGxE) x rep = agricolae's ecovalence")
check(max(abs(4 * st$shukla_variance[i] - ag[["Sigma-square"]])) < 1e-5, "Shukla's stability variance x rep = agricolae::stability.par Sigma-square")
y <- as.matrix(stb); ge <- sweep(sweep(y, 1, rowMeans(y)), 2, colMeans(y)) + mean(y)
check(abs(mean(st$shukla_variance) - sum(ge^2) / ((nrow(y) - 1) * (ncol(y) - 1))) < 1e-8, "mean Shukla variance = G x E mean square")
check(max(abs(st$static_env_variance[i] - apply(y, 1, stats::var))) < 1e-8, "statgenGxE 'static' = each genotype's variance across environments")
check(max(abs(st$superiority[i] - rowSums(sweep(y, 2, apply(y, 2, max))^2) / (2 * ncol(y)))) < 1e-8, "superiority = sum (y - best)^2 / (2 q) (Lin & Binns)")

# 6. AMMI recovers a rank-1 interaction; missing cells imputed or excluded; too many missing -> skipped
set.seed(7)
d <- grid_met(60, 8); gi <- match(d$genotype, gid(60)); ej <- match(d$env, sprintf("E%02d", 1:8))
u <- rnorm(60); v <- rnorm(8)
d$y <- 10 + rnorm(60)[gi] + rnorm(8, 0, 2)[ej] + 1.5 * u[gi] * v[ej] + rnorm(nrow(d), 0, 0.2); d$se2 <- NA
td <- make_td(d, "y")
am <- ammi_like(td, "y", 2, FALSE, character())
check(am$importance["Proportion of Variance", "PC1"] > 0.8, "AMMI: PC1 carries most of a rank-1 interaction")
check(abs(cor(am$genoScores[gid(60), "PC1"], u)) > 0.95, "AMMI: PC1 genotype scores track the true interaction scores")
gge <- ammi_like(td, "y", 2, TRUE, character())
check(gge$importance["Proportion of Variance", "PC1"] > 0 && nrow(scores_table(gge, 2)) == 68, "GGE: scores for 60 genotypes + 8 environments")
dm <- d[-sample(nrow(d), 24), ]                      # 5% of cells missing
cm <- completeness(dm)
check(cm$missing == 24 && cm$cells == 480, "completeness counts missing genotype x environment cells")
ami <- ammi_like(make_td(dm, "y"), "y", 2, FALSE, character())
check(nrow(ami$fitted) == 480 && !anyNA(ami$fitted$fittedValue), "AMMI: missing cells are imputed, fitted values for all 480 cells")
exc <- setdiff(rownames(cm$table), cm$complete_genotypes)
ame <- ammi_like(make_td(dm, "y"), "y", 2, FALSE, exc)
check(nrow(ame$genoScores) == length(cm$complete_genotypes), "--complete-only: incomplete genotypes excluded from AMMI")
dbig <- d[-sample(nrow(d), 180), ]                   # 37.5% missing
check(inherits(ammi_like(make_td(dbig, "y"), "y", 2, FALSE, character()), "skipped"), "AMMI with > 30% missing cells is skipped, not crashed")

# 7. Mega-environments: two groups of environments with different winners
set.seed(9)
d <- grid_met(30, 8); gi <- match(d$genotype, gid(30)); ej <- match(d$env, sprintf("E%02d", 1:8))
grp <- rep(c(1, 2), each = 4)
spec <- ifelse(seq_len(30) == 1, 1, ifelse(seq_len(30) == 2, 2, 0))   # G001 wins group 1, G002 group 2
d$y <- 10 + rnorm(8)[ej] + rnorm(30, 0, 0.3)[gi] + 4 * (spec[gi] == grp[ej]) + rnorm(nrow(d), 0, 0.2); d$se2 <- NA
me <- run_megaenv(make_td(d, "y"), "y", "max", character(), d)
tab <- me$summTab
win <- setNames(as.character(tab$Winning_genotype), as.character(tab$Trial))
check(all(win[sprintf("E%02d", 1:4)] == "G001") && all(win[sprintf("E%02d", 5:8)] == "G002") && length(unique(tab$Mega_factor)) == 2,
      "mega-environments: the two simulated groups and their winners are recovered")

# 8. Smoke test on statgenGxE's dropsPheno (50 genotypes x 10 experiments)
data("dropsPheno", package = "statgenGxE", envir = environment())
sel <- sort(unique(dropsPheno$Variety_ID))[1:50]
dp <- dropsPheno[dropsPheno$Variety_ID %in% sel, c("Experiment", "Variety_ID", "grain.yield")]
names(dp) <- c("trial", "genotype", "grain_yield_t_ha")
prep <- prepare_met(dp, O(trait = "grain_yield_t_ha"))
res <- analyze_met(prep$dat, O(trait = "grain_yield_t_ha"), prep$weight_source)
check(res$n_env == 10 && res$n_geno == 50, "dropsPheno: 10 environments x 50 genotypes")
check(res$varcomp$h2 > 0 && res$varcomp$h2 < 1 && all(res$varcomp$components$variance > 0), "dropsPheno: positive variances, 0 < H2 < 1")
check(cor(res$genotypes$blup, res$genotypes$raw_mean) > 0.98, "dropsPheno: MET BLUPs agree with genotype means for a complete table")
check(abs(mean(res$fw$estimates$Sens) - 1) < 1e-3 && !inherits(res$ammi, "skipped") && !inherits(res$gge, "skipped"),
      "dropsPheno: FW, AMMI and GGE run")
check(res$gge$importance["Proportion of Variance", "PC1"] > res$ammi$importance["Proportion of Variance", "PC1"],
      "dropsPheno: GGE PC1 (G + GxE) takes a larger share than AMMI PC1 (GxE only)")
check(nrow(res$stability) == 50 && !inherits(res$megaenv, "skipped"), "dropsPheno: stability for 50 genotypes and mega-environments")

# 9. Bad inputs fail loudly
base <- data.frame(trial = rep(c("A", "B", "C"), each = 4), genotype = rep(paste0("g", 1:4), 3), yld = rnorm(12))
expect_error(prepare_met(rbind(base, base[1, ]), O(trait = "yld")), "more than once", "duplicated genotype x environment refused")
expect_error(prepare_met(base, O(trait = "yield")), "Column 'yield' not found", "missing trait column named in the error")
expect_error(prepare_met(base[base$trial == "A", ], O(trait = "yld")), "at least 2", "a single environment is refused")
bad <- base; bad$yld <- as.character(bad$yld); bad$yld[2] <- "lodged"
expect_error(prepare_met(bad, O(trait = "yld")), "non-numeric", "non-numeric trait refused")
ly <- cbind(base, loc = "L1", yr = "2025")
expect_error(prepare_met(ly, O(trait = "yld", location = "loc", year = "yr")), "location-year", "location x year with shared location-years refused")
check(abs(shukla_variance(c(1, 2, 3, 4), 5)[1] - (4 * 1 / (2 * 4) - 10 / (3 * 2 * 4))) < 1e-12, "Shukla formula: hand-computed value")

# 10. End to end on the demo; two input files stack to the same result
tmp <- tempfile("claw_met_"); dir.create(tmp)
out <- file.path(tmp, "demo")
r <- run_cli(c("--demo", "--out", shQuote(out)))
check(r$status == 0L, paste("demo run exits cleanly:", paste(utils::tail(r$output, 5), collapse = "\n")))
for (f in c("report.md", "tables/genotype_blups.csv", "tables/variance_components.csv", "tables/heritability.csv",
            "tables/ammi_anova.csv", "tables/pc_importance.csv", "tables/finlay_wilkinson.csv", "tables/stability.csv",
            "tables/mega_environments.csv", "figures/gge_biplot.png", "figures/ammi2_biplot.png",
            "figures/finlay_wilkinson.png", "figures/stability_vs_mean.png", "reproducibility/commands.sh", "reproducibility/sessionInfo.txt"))
  check(file.exists(file.path(out, f)), paste("demo writes", f))
gb <- utils::read.csv(file.path(out, "tables", "genotype_blups.csv"), stringsAsFactors = FALSE)
truth <- utils::read.csv(file.path(SKILL_TEST_DIR, "demo", "met_truth.csv"), stringsAsFactors = FALSE)
check(nrow(gb) == 40 && "blup_weighted" %in% names(gb), "demo: 40 genotypes with weighted BLUPs (SEs used)")
rg <- cor(gb$blup_weighted, truth$g_kg_ha[match(gb$genotype, truth$genotype)])
check(rg > 0.8, sprintf("demo: weighted BLUPs track the simulated genotype effects (r = %.2f > 0.8)", rg))
demo <- utils::read.csv(file.path(SKILL_TEST_DIR, "demo", "met_input.csv"), stringsAsFactors = FALSE)
p1 <- file.path(tmp, "y2025.csv"); p2 <- file.path(tmp, "y2026.csv")
utils::write.csv(demo[demo$year == 2025, ], p1, row.names = FALSE)
utils::write.csv(demo[demo$year == 2026, ], p2, row.names = FALSE)
rd <- read_inputs(paste(p1, p2, sep = ","))
pp <- prepare_met(rd$data, O(trait = "yield_kg_ha"))
check(nrow(pp$dat) == nrow(demo) && pp$weight_source == "wt_yield_kg_ha", "two input files stack into one MET with weights")
nw <- prepare_met(rd$data, O(trait = "yield_kg_ha", no_weights = TRUE))
check(nw$weight_source == "none" && all(is.na(nw$dat$se2)), "--no-weights ignores the weight column")
bad <- run_cli(c("--input", shQuote(p1)))
check(bad$status != 0L && any(grepl("--trait", bad$output)), "missing --trait exits with a helpful error")

unlink(tmp, recursive = TRUE)
cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
