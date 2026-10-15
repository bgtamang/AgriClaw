# Generate the claw-met-gxe demo data (run from anywhere: Rscript demo/make_demo.R).
# Writes demo/met_input.csv: per-trial genotype BLUEs in the layout written by
# claw-trial-analysis (trial, genotype, <trait>, se_<trait>, wt_<trait>), here
# for 40 soybean genotypes in 8 trials (4 locations x 2 years), plus location
# and year columns. Simulated with known structure (seed 2026):
#   y_ij = mu + e_j * b_i + g_i + u_i * v_j + ge_ij + stage-1 error
# b_i = Finlay-Wilkinson sensitivity, u_i * v_j = a crossover interaction that
# splits the trials into two groups, ge_ij = residual G x E, and the stage-1
# error has the trial's standard error. Three genotype x trial cells are missing.
# The true values are written to demo/met_truth.csv.

demo_dir <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else "demo"
})
set.seed(2026)
geno <- sprintf("SB%02d", 1:40)
loc <- c("AMES", "NASHUA", "LEWIS", "CRAWFORDSVILLE")
yr <- c(2025, 2026)
trials <- expand.grid(location = loc, year = yr, stringsAsFactors = FALSE)
trials$trial <- paste(trials$location, trials$year, sep = "-")
ne <- nrow(trials)
e <- c(600, -200, 300, -700, 900, -400, 100, -600)          # environment effects (kg/ha)
v <- c(1, 1, -1, -1, 1, 1, -1, -1)                          # crossover: locations 1-2 vs 3-4
se <- c(180, 220, 260, 300, 160, 240, 280, 380)             # stage-1 SE of a BLUE in each trial
g <- stats::rnorm(40, 0, 300)
b <- stats::rnorm(40, 1, 0.25)
u <- stats::rnorm(40, 0, 180)
d <- expand.grid(genotype = geno, j = seq_len(ne), stringsAsFactors = FALSE)
i <- match(d$genotype, geno)
d$trial <- trials$trial[d$j]
d$location <- trials$location[d$j]
d$year <- trials$year[d$j]
se_ij <- se[d$j] * stats::runif(nrow(d), 0.9, 1.1)
d$yield_kg_ha <- round(3800 + e[d$j] * b[i] + g[i] + u[i] * v[d$j] + stats::rnorm(nrow(d), 0, 120) +
                         stats::rnorm(nrow(d), 0, se_ij))
d$se_yield_kg_ha <- round(se_ij, 1)
d$wt_yield_kg_ha <- signif(1 / d$se_yield_kg_ha^2, 6)
drop <- c(which(d$genotype == "SB07" & d$trial == "LEWIS-2025"),
          which(d$genotype == "SB19" & d$trial == "AMES-2026"),
          which(d$genotype == "SB33" & d$trial == "NASHUA-2025"))
d <- d[-drop, c("trial", "location", "year", "genotype", "yield_kg_ha", "se_yield_kg_ha", "wt_yield_kg_ha")]
utils::write.csv(d, file.path(demo_dir, "met_input.csv"), row.names = FALSE)
utils::write.csv(data.frame(genotype = geno, g_kg_ha = round(g, 1), sensitivity_b = round(b, 3), crossover_u = round(u, 1)),
                 file.path(demo_dir, "met_truth.csv"), row.names = FALSE)
cat("Wrote", nrow(d), "rows.\n")
