# Generate the claw-trial-analysis demo data (run from anywhere: Rscript demo/make_demo.R).
# A synthetic resolvable incomplete-block trial in the layout written by
# claw-fieldbook (clean_data.csv): 100 entries (98 lines + 2 checks), 2
# replicates side by side in a 10-row x 20-column field, blocks of 10 plots
# (one field row within a replicate). Yield has known genetic and residual
# variances plus a smooth spatial trend; one plot has a planted outlier and one
# plot is missing. Seed 2026.

demo_dir <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else "demo"
})
set.seed(2026)
entries <- c("CHECK-A", "CHECK-B", sprintf("L%03d", 1:98))
d <- do.call(rbind, lapply(1:2, function(r) {
  data.frame(rep = r, entry = sample(entries), row = rep(1:10, each = 10),
             col = rep(1:10, 10) + 10 * (r - 1), stringsAsFactors = FALSE)
}))
d$block <- (d$rep - 1) * 10 + d$row
d <- d[order(d$row, d$col), ]
d$plot <- seq_len(nrow(d))
d$trial <- "DEMO-AI"
d$plot_id <- paste(d$trial, d$plot, sep = "-")
d$design <- "alpha"
d$entry_no <- match(d$entry, entries)
d$is_check <- d$entry %in% c("CHECK-A", "CHECK-B")
d$is_filler <- FALSE

# True values: sigma2_g = 400^2, sigma2_e = 300^2 (yield); entry-mean H2 with 2 reps = 0.78
g_y <- stats::setNames(stats::rnorm(length(entries), 0, 400), entries)
g_h <- stats::setNames(stats::rnorm(length(entries), 0, 6), entries)
trend <- 500 * sin(d$row / 3) + 300 * cos(d$col / 4)
d$yield_kg_ha <- round(4000 + g_y[d$entry] + trend + stats::rnorm(nrow(d), 0, 300))
d$plant_height_cm <- round(90 + g_h[d$entry] + trend / 100 + stats::rnorm(nrow(d), 0, 4), 1)
d$yield_kg_ha[d$plot == 57] <- d$yield_kg_ha[d$plot == 57] + 1800   # planted outlier
d$yield_kg_ha[d$plot == 133] <- NA                                    # missing plot
d$qc_issues <- ""

cols <- c("trial", "plot_id", "plot", "row", "col", "rep", "block", "design", "entry_no", "entry",
          "is_check", "is_filler", "yield_kg_ha", "plant_height_cm", "qc_issues")
utils::write.csv(d[, cols], file.path(demo_dir, "trial_clean.csv"), row.names = FALSE, na = "")
utils::write.csv(data.frame(entry = entries, true_yield_effect = round(g_y, 1), true_height_effect = round(g_h, 2)),
                 file.path(demo_dir, "true_genetic_values.csv"), row.names = FALSE)
cat("Wrote", nrow(d), "plots.\n")
