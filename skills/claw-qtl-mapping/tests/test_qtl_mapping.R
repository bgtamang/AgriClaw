# Tests for claw-qtl-mapping. Run from anywhere: Rscript tests/test_qtl_mapping.R
# Known answers: the peaks, intervals and permutation thresholds documented in
# the R/qtl2 user guide for the bundled iron data; simulated RIL populations
# with a known QTL (position, additive effect); a simulated 4-way MAGIC-type RIL
# population with known founder allele effects; the algebraic identity between
# Haley-Knott LOD and R^2; QC design checks (duplicates, distortion, codes, units).

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "qtl_mapping.R"))   # loads functions only
source(file.path(SKILL_TEST_DIR, "demo", "make_demo.R"))         # RIL simulators

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
near <- function(x, y, tol) isTRUE(all(abs(x - y) <= tol))

run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "qtl_mapping.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
dd <- file.path(SKILL_TEST_DIR, "demo")

# 0. The RIL simulator used below and in the demo: one crossover per Morgan per
#    meiosis (Haldane), and in RILs by selfing the recombination fraction between
#    loci r cM apart expands to R = 2r/(1+2r) (Haldane & Waddington 1931), so the
#    genotype correlation is 1 - 2R.
set.seed(17)
pos0 <- seq(0, 100, by = 1)
sw <- replicate(3000, sum(diff(meiosis(rep(1L, 101), rep(2L, 101), pos0, 100)) != 0))
check(near(mean(sw), 1, 0.05), sprintf("simulator: %.3f crossovers per meiosis per Morgan", mean(sw)))
lab0 <- simulate_ril_labels(1500, list(pos0), 100, generations = 10L)[[1]]
for (dcm in c(5, 20)) {
  r <- 0.5 * (1 - exp(-2 * dcm / 100))
  R <- 2 * r / (1 + 2 * r)
  rho <- stats::cor(lab0[, 31], lab0[, 31 + dcm], use = "complete.obs")
  check(near(rho, 1 - 2 * R, 0.05), sprintf("simulator: RIL genotype correlation at %g cM %.3f vs %.3f", dcm, rho, 1 - 2 * R))
}

# 1. R/qtl2 user guide, iron data (F2 intercross with an X chromosome):
#    insert_pseudomarkers(step = 1), calc_genoprob(error_prob = 0.002),
#    scan1 with Xcovar, then find_peaks(threshold = 4, drop = 1.5).
iron <- read_cross_control(system.file("extdata", "iron.zip", package = "qtl2"))$cross
sc_i <- run_scan(iron, c("liver", "spleen"), step = 1, error_prob = 0.002)
pk <- get_peaks(sc_i, threshold = 4, interval = "lod", drop = 1.5)
check(nrow(pk) == 5L, "iron: five peaks above LOD 4 (user guide)")
check(identical(as.character(pk$chr), c("2", "7", "16", "8", "9")), "iron: peak chromosomes")
check(near(pk$pos, c(56.8, 50.1, 28.6, 13.6, 56.6), 1e-8), "iron: peak positions (cM)")
check(near(pk$lod, c(4.9576, 4.0508, 7.6816, 4.3029, 12.0632), 5e-5), "iron: peak LOD scores")
check(near(pk$ci_lo, c(48.1, 13.1, 6.6, 0.0, 53.6), 1e-8) &&
        near(pk$ci_hi, c(73.2, 53.6, 40.4, 32.7, 61.2), 1e-8),
      "iron: 1.5-LOD support intervals (user guide)")
pk2 <- get_peaks(sc_i, threshold = 4, interval = "lod", drop = 1.5, peakdrop = 1.8)
p7 <- pk2[pk2$chr == "7", ]
check(nrow(pk2) == 6L && near(p7$pos, c(25.1, 50.1), 1e-8) && near(p7$lod, c(4.0400, 4.0508), 5e-5) &&
        near(p7$ci_lo, c(13.1, 31.7), 1e-8) && near(p7$ci_hi, c(28.4, 53.6), 1e-8),
      "iron: two chr 7 peaks with peakdrop = 1.8 (user guide)")
pb <- get_peaks(sc_i, threshold = 4, interval = "bayes", prob = 0.95)
p9 <- pb[pb$chr == "9", ]
check(near(c(p9$ci_lo, p9$pos, p9$ci_hi), c(53.6, 56.6, 61.2), 1e-8),
      "iron: 95% Bayes credible interval on chr 9 for spleen (user guide)")

# 2. Effects and variance explained at a peak.
e2 <- peak_effects(sc_i, "liver", "2", 56.8)
check(near(attr(e2, "lod"), 4.957564, 1e-5), "fit1 at the peak reproduces the scan LOD")
gm <- stats::setNames(e2$estimate[e2$type == "genotype_mean"], e2$effect[e2$type == "genotype_mean"])
a2 <- e2$estimate[e2$effect == "a"]
check(near(a2, (gm[["BB"]] - gm[["SS"]]) / 2, 1e-6), "additive effect a = (BB - AA)/2")
check(near(e2$estimate[e2$effect == "d"], gm[["SB"]] - (gm[["SS"]] + gm[["BB"]]) / 2, 1e-6),
      "dominance effect d = AB - (AA + BB)/2")
pr16 <- qtl2::pull_genoprobpos(sc_i$probs, sc_i$map, "16", 28.6)
y16 <- iron$pheno[, "liver"]
r2 <- summary(stats::lm(y16 ~ pr16[, -1]))$r.squared
check(near(pve_from_lod(7.681569, length(y16)), r2, 1e-6),
      "PVE = 1 - 10^(-2 LOD/n) equals the regression R^2 at the peak (Haley-Knott)")

# 2b. Additive covariates: sex as a factor covariate (dummy-coded by the script)
#     must give the same LODs as the numeric sex covariate of the ?scan1 example.
sc_sex <- run_scan(iron, "liver", step = 1, error_prob = 0.002,
                   addcovar = covariate_matrix(iron, "sex"))
sexnum <- stats::setNames(match(iron$covar$sex, c("f", "m")), rownames(iron$covar))
ref <- qtl2::scan1(sc_sex$probs, iron$pheno[, "liver", drop = FALSE], addcovar = sexnum,
                   Xcovar = qtl2::get_x_covar(iron))
check(max(abs(unclass(sc_sex$out) - unclass(ref))) < 1e-8,
      "factor covariate (sex) reproduces qtl2::scan1 with a numeric covariate")
err <- tryCatch(covariate_matrix(iron, "batch"), error = conditionMessage)
check(grepl("not found", err), "unknown covariate name is rejected")

# 3. Permutation thresholds. The guide reports 1000 permutations stratified by
#    the X covariates: 20% threshold 2.63 and 5% threshold 3.46 for liver. With
#    100 permutations the Monte Carlo error is larger, hence tolerances 0.3 / 0.45.
sc_liver <- sc_i
sc_liver$pheno <- sc_i$pheno[, "liver", drop = FALSE]
pm <- run_perm(sc_liver, n_perm = 100, seed = 7, alpha = c(0.2, 0.05), x_specific = FALSE)
check(near(pm$thrA["0.2", "liver"], 2.63, 0.3), sprintf("iron liver 20%% threshold %.2f vs 2.63", pm$thrA["0.2", "liver"]))
check(near(pm$thrA["0.05", "liver"], 3.46, 0.45), sprintf("iron liver 5%% threshold %.2f vs 3.46", pm$thrA["0.05", "liver"]))
check(near(pm$thrA["0.05", "liver"], stats::quantile(pm$operm$liver[, 1], 0.95), 1e-12),
      "threshold equals the 95% quantile of the permutation maxima")
cat(sprintf("iron liver thresholds (100 permutations): 20%% %.2f (guide 2.63), 5%% %.2f (guide 3.46)\n",
            pm$thrA["0.2", "liver"], pm$thrA["0.05", "liver"]))

# 4. Simulated RIL populations (selfing, F8; 150 lines; 2 chromosomes; markers
#    every 2 cM; QTL on chr 2 at 47.3 cM between markers, additive effect 1,
#    PVE 0.20). Over 20 populations the 1.5-LOD support interval must contain the
#    true position in at least 18 (nominal coverage about 95%; P(>= 18 of 20) =
#    0.92 at 95%), and the 95% Bayes interval in at least 17.
sim_ril_cross <- function(seed, n = 150, qpos = 47.3, a = 1, pve = 0.2) {
  set.seed(seed)
  mpos <- list(seq(0, 100, by = 2), seq(0, 100, by = 2))
  pos_all <- mpos
  pos_all[[2]] <- sort(c(mpos[[2]], qpos))
  lab <- simulate_ril_labels(n, pos_all, c(100, 100), generations = 8L)
  ids <- sprintf("L%03d", seq_len(n))
  geno <- do.call(cbind, lapply(1:2, function(c) {
    g <- lab[[c]][, pos_all[[c]] %in% mpos[[c]], drop = FALSE]
    matrix(ifelse(is.na(g), "-", ifelse(g == 1L, "A", "B")), n, ncol(g),
           dimnames = list(ids, sprintf("c%d_%03d", c, mpos[[c]])))
  }))
  gm <- data.frame(marker = colnames(geno), chr = rep(1:2, lengths(mpos)), pos = unlist(mpos))
  q <- lab[[2]][, match(qpos, pos_all[[2]])]
  x <- ifelse(is.na(q), 0, ifelse(q == 2L, 1, -1))
  y <- a * x + stats::rnorm(n, 0, sqrt(a^2 * (1 - pve) / pve))
  d <- tempfile("sim")
  dir.create(d)
  utils::write.csv(data.frame(id = ids, geno, check.names = FALSE), file.path(d, "g.csv"),
                   row.names = FALSE, quote = FALSE)
  utils::write.csv(gm, file.path(d, "m.csv"), row.names = FALSE, quote = FALSE)
  utils::write.csv(data.frame(id = ids, y = y), file.path(d, "p.csv"), row.names = FALSE)
  build_cross_from_csv(file.path(d, "g.csv"), file.path(d, "m.csv"), file.path(d, "p.csv"),
                       "riself", work_dir = file.path(d, "q"))$cross
}
sims <- t(vapply(1:20, function(s) {
  cr <- sim_ril_cross(2000 + s)
  sc <- run_scan(cr, "y", step = 1)
  pl <- get_peaks(sc, threshold = 0, interval = "lod", drop = 1.5)
  pl <- pl[pl$chr == "2", ]
  pb <- get_peaks(sc, threshold = 0, interval = "bayes", prob = 0.95)
  pb <- pb[pb$chr == "2", ]
  e <- peak_effects(sc, "y", "2", pl$pos)
  c(lod = pl$lod, pos = pl$pos, cover_lod = pl$ci_lo <= 47.3 && 47.3 <= pl$ci_hi,
    cover_bayes = pb$ci_lo <= 47.3 && 47.3 <= pb$ci_hi, a = e$estimate[e$effect == "a"])
}, numeric(5)))
check(sum(sims[, "cover_lod"]) >= 18, sprintf("1.5-LOD interval coverage %d/20", sum(sims[, "cover_lod"])))
check(sum(sims[, "cover_bayes"]) >= 17, sprintf("95%% Bayes interval coverage %d/20", sum(sims[, "cover_bayes"])))
check(all(sims[, "lod"] > 3), "QTL detected (LOD > 3) in every simulated population")
check(mean(abs(sims[, "pos"] - 47.3)) < 5, "mean position error under 5 cM")
check(near(mean(sims[, "a"]), 1, 0.1), sprintf("mean additive effect %.3f within 0.1 of the true 1", mean(sims[, "a"])))
cat(sprintf("RIL simulation: coverage LOD %d/20, Bayes %d/20; mean |pos error| %.2f cM; mean a %.3f\n",
            sum(sims[, "cover_lod"]), sum(sims[, "cover_bayes"]), mean(abs(sims[, "pos"] - 47.3)),
            mean(sims[, "a"])))

# 5. Simulated 4-way RIL by selfing (qtl2 cross type riself4) from founders
#    A-D with random funnels, SNP genotypes, and a QTL on chr 1 at 31.5 cM with
#    founder allele effects (+1.5, 0, 0, -1.5). Control-file route.
make_riself4 <- function(seed = 44, n = 200, qpos = 31.5) {
  set.seed(seed)
  mpos <- list(seq(0, 100, by = 2), seq(0, 100, by = 2))
  pos_all <- mpos
  pos_all[[1]] <- sort(c(mpos[[1]], qpos))
  funnel <- replicate(n, sample(4L), simplify = FALSE)
  lab <- simulate_ril_labels(n, pos_all, c(100, 100), generations = 8L, n_founders = 4L,
                             funnel = funnel)
  ids <- sprintf("M%03d", seq_len(n))
  fg <- list()
  geno <- list()
  for (c in 1:2) {
    nm <- sprintf("s%d_%03d", c, mpos[[c]])
    pat <- t(replicate(length(nm), {
      repeat { p <- sample(0:1, 4, replace = TRUE); if (any(p == 0) && any(p == 1)) break }
      p
    }))
    fg[[c]] <- matrix(ifelse(t(pat) == 1L, "B", "A"), 4, length(nm), dimnames = list(LETTERS[1:4], nm))
    g <- lab[[c]][, pos_all[[c]] %in% mpos[[c]], drop = FALSE]
    geno[[c]] <- matrix("-", n, length(nm), dimnames = list(ids, nm))
    ok <- !is.na(g)
    geno[[c]][ok] <- ifelse(pat[cbind(col(g)[ok], g[ok])] == 1L, "B", "A")
  }
  eff <- c(1.5, 0, 0, -1.5)
  q <- lab[[1]][, match(qpos, pos_all[[1]])]
  y <- ifelse(is.na(q), 0, eff[q]) + stats::rnorm(n, 0, 1.5)
  d <- tempfile("riself4")
  dir.create(d)
  utils::write.csv(data.frame(id = ids, do.call(cbind, geno), check.names = FALSE),
                   file.path(d, "geno.csv"), row.names = FALSE, quote = FALSE)
  utils::write.csv(data.frame(id = LETTERS[1:4], do.call(cbind, fg), check.names = FALSE),
                   file.path(d, "fgeno.csv"), row.names = FALSE, quote = FALSE)
  utils::write.csv(data.frame(marker = unlist(lapply(fg, colnames)), chr = rep(1:2, lengths(mpos)),
                              pos = unlist(mpos)), file.path(d, "gmap.csv"), row.names = FALSE, quote = FALSE)
  utils::write.csv(data.frame(id = ids, y = round(y, 3)), file.path(d, "pheno.csv"), row.names = FALSE)
  ci <- do.call(rbind, funnel)
  utils::write.csv(data.frame(id = ids, a = ci[, 1], b = ci[, 2], c = ci[, 3], d = ci[, 4]),
                   file.path(d, "crossinfo.csv"), row.names = FALSE, quote = FALSE)
  ctrl <- file.path(d, "magic4.yaml")
  qtl2::write_control_file(ctrl, crosstype = "riself4", geno_file = "geno.csv",
                           founder_geno_file = "fgeno.csv", gmap_file = "gmap.csv",
                           pheno_file = "pheno.csv", crossinfo_file = "crossinfo.csv",
                           geno_codes = c(A = 1L, B = 3L), alleles = LETTERS[1:4],
                           na.strings = c("-", "NA"), overwrite = TRUE)
  list(control = ctrl, eff = eff, qpos = qpos)
}
m4 <- make_riself4()
cr4 <- read_cross_control(m4$control)$cross
sc4 <- run_scan(cr4, "y", step = 1)
check(sc4$multiparent && dim(sc4$probs[["1"]])[2] == 4L, "riself4: scan uses 4 founder allele probabilities")
pk4 <- get_peaks(sc4, threshold = 3, interval = "lod", drop = 1.5)
p4 <- pk4[pk4$chr == "1", ]
check(nrow(pk4) == 1L && nrow(p4) == 1L && p4$ci_lo <= m4$qpos && m4$qpos <= p4$ci_hi,
      sprintf("riself4: one QTL, on chr 1, interval %.1f-%.1f contains 31.5", p4$ci_lo, p4$ci_hi))
e4 <- peak_effects(sc4, "y", "1", p4$pos)
ae <- stats::setNames(e4$estimate[e4$effect %in% LETTERS[1:4]], e4$effect[e4$effect %in% LETTERS[1:4]])
check(stats::cor(ae[LETTERS[1:4]], m4$eff) > 0.95 && ae[["A"]] > 0.8 && ae[["D"]] < -0.8,
      sprintf("riself4: founder allele effects recovered (A %.2f, B %.2f, C %.2f, D %.2f)",
              ae[["A"]], ae[["B"]], ae[["C"]], ae[["D"]]))
cat(sprintf("riself4: peak chr 1 at %.1f cM (interval %.1f-%.1f, LOD %.1f); allele effects A %.2f B %.2f C %.2f D %.2f (true +1.5 0 0 -1.5, centred)\n",
            p4$pos, p4$ci_lo, p4$ci_hi, p4$lod, ae[["A"]], ae[["B"]], ae[["C"]], ae[["D"]]))
sc4l <- run_scan(cr4, "y", step = 1, kinship = "loco")
pl4 <- get_peaks(sc4l, threshold = 3, interval = "lod", drop = 1.5)
check(nrow(pl4) == 1L && abs(pl4$pos - p4$pos) <= 5, "LOCO mixed-model scan finds the same QTL")

# 6. QC and input checks.
qcd <- qc_cross(iron)
check(all(c("marker", "missing_rate", "seg_dist_p", "distorted") %in% names(qcd$markers)),
      "marker QC table")
check(all(is.na(qcd$markers$seg_dist_p[qcd$markers$chr == "X"])),
      "no segregation test on the X chromosome")
gd <- utils::read.csv(file.path(dd, "geno.csv"), colClasses = "character", check.names = FALSE)
gd[2, -1] <- gd[1, -1]                                # line 2 duplicates line 1
gd[, "m1_05"] <- ifelse(seq_len(nrow(gd)) <= 140, "A", "B")   # 140 A : 20 B
tq <- tempfile("qc")
dir.create(tq)
utils::write.csv(gd, file.path(tq, "g.csv"), row.names = FALSE, quote = FALSE)
crq <- build_cross_from_csv(file.path(tq, "g.csv"), file.path(dd, "gmap.csv"), file.path(dd, "pheno.csv"),
                            "riself", work_dir = file.path(tq, "w"))$cross
qq <- qc_cross(crq)
check(!is.null(qq$duplicates) && any(qq$duplicates$ind1 == "RIL001" & qq$duplicates$ind2 == "RIL002"),
      "duplicated individuals are flagged")
check(qq$markers$distorted[qq$markers$marker == "m1_05"] && sum(qq$markers$distorted) == 1L,
      "a 7:1 marker is flagged as distorted and normal markers are not")
gh <- utils::read.csv(file.path(dd, "geno.csv"), colClasses = "character", check.names = FALSE)
gh[1, 2] <- "H"
utils::write.csv(gh, file.path(tq, "gh.csv"), row.names = FALSE, quote = FALSE)
err <- tryCatch(build_cross_from_csv(file.path(tq, "gh.csv"), file.path(dd, "gmap.csv"),
                                     file.path(dd, "pheno.csv"), "riself", work_dir = tempfile()),
                error = conditionMessage)
check(grepl("Heterozygous", err), "heterozygous calls in a riself cross are rejected with advice")
pm_bp <- utils::read.csv(file.path(dd, "pmap.csv"))
pm_bp$pos <- round(pm_bp$pos * 1e6)
utils::write.csv(pm_bp, file.path(tq, "map_bp.csv"), row.names = FALSE, quote = FALSE)
err <- tryCatch(build_cross_from_csv(file.path(dd, "geno.csv"), file.path(tq, "map_bp.csv"),
                                     file.path(dd, "pheno.csv"), "riself", work_dir = tempfile()),
                error = conditionMessage)
check(grepl("base pairs or Mb", err), "base-pair positions given as cM are detected")
err <- tryCatch(build_cross_from_csv(file.path(dd, "geno.csv"), file.path(tq, "map_bp.csv"),
                                     file.path(dd, "pheno.csv"), "riself", map_units = "bp",
                                     work_dir = tempfile()), error = conditionMessage)
check(grepl("genetic map", err) && grepl("cm-per-mb", err),
      "physical map without --cm-per-mb asks for a genetic map")
conv <- build_cross_from_csv(file.path(dd, "geno.csv"), file.path(tq, "map_bp.csv"),
                             file.path(dd, "pheno.csv"), "riself", map_units = "bp",
                             cm_per_mb = 2, work_dir = tempfile())
check(any(grepl("genome-average", conv$notes)) &&
        near(max(conv$cross$gmap[["1"]]), 2 * max(pm_bp$pos[pm_bp$chr == 1]) / 1e6, 1e-6),
      "--cm-per-mb converts positions and warns that it is an approximation")
err <- tryCatch(build_cross_from_csv(file.path(dd, "geno.csv"), file.path(dd, "gmap.csv"),
                                     file.path(dd, "pheno.csv"), "riself4", work_dir = tempfile()),
                error = conditionMessage)
check(grepl("--control", err), "multiparent cross types are directed to the control-file route")

# 7. End to end: demo (CSV route); every simulated QTL lies in its interval.
out <- file.path(tempdir(), "claw-qtl-mapping-test")
unlink(out, recursive = TRUE)
res <- run_cli(c("--demo", "--n-perm", "200", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "tables/peaks.csv", "tables/peak_effects.csv", "tables/thresholds.csv",
            "tables/lod_scores.csv", "tables/marker_qc.csv", "tables/individual_qc.csv",
            "figures/lod_days_to_flowering.png", "figures/lod_plant_height_cm.png",
            "qtl2_input/cross.yaml", "reproducibility/commands.sh", "reproducibility/sessionInfo.txt")) {
  check(file.exists(file.path(out, f)), paste("demo writes", f))
}
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
pkd <- utils::read.csv(file.path(out, "tables", "peaks.csv"), colClasses = c(chr = "character"))
tru <- utils::read.csv(file.path(dd, "truth.csv"))
hit <- vapply(seq_len(nrow(tru)), function(k) {
  p <- pkd[pkd$trait == tru$trait[k] & pkd$chr == as.character(tru$chr[k]), ]
  nrow(p) == 1L && p$ci_lo_cM <= tru$pos_cM[k] && tru$pos_cM[k] <= p$ci_hi_cM &&
    abs(p$a_effect - tru$additive_effect[k]) < 3 * p$a_se
}, TRUE)
check(all(hit) && nrow(pkd) == 3L, "demo: the 3 simulated QTL are found, inside their intervals, effects within 3 SE")
check(all(c("ci_lo_Mb", "ci_hi_Mb") %in% names(pkd)) && all(pkd$ci_lo_Mb <= pkd$ci_hi_Mb),
      "intervals are converted to Mb with the physical map")
thr <- utils::read.csv(file.path(out, "tables", "thresholds.csv"))
check(all(thr$threshold_autosomes[thr$alpha == 0.05] > thr$threshold_autosomes[thr$alpha == 0.1]),
      "5% thresholds exceed 10% thresholds")
rep_txt <- readLines(file.path(out, "report.md"))
check(any(grepl("Only 200 permutations", rep_txt)) && any(grepl("Single-QTL genome scans", rep_txt)),
      "report warns about few permutations and states the single-QTL limitation")

# 8. End to end: multiparent control file (riself4) through the command line.
out4 <- file.path(tempdir(), "claw-qtl-mapping-test-magic")
unlink(out4, recursive = TRUE)
res4 <- run_cli(c("--control", shQuote(m4$control), "--n-perm", "50", "--out", shQuote(out4)))
check(res4$status == 0L, paste("multiparent run exits cleanly:", paste(res4$output, collapse = "\n")))
pe4 <- utils::read.csv(file.path(out4, "tables", "peak_effects.csv"))
check(all(LETTERS[1:4] %in% pe4$effect[pe4$type == "allele_effect"]),
      "multiparent run reports founder allele effects")

# 9. Seed reproducibility of permutation thresholds.
crd <- build_cross_from_csv(file.path(dd, "geno.csv"), file.path(dd, "gmap.csv"), file.path(dd, "pheno.csv"),
                            "riself", work_dir = tempfile())$cross
scd <- run_scan(crd, "plant_height_cm", step = 0)
check(identical(run_perm(scd, 50, seed = 3)$table, run_perm(scd, 50, seed = 3)$table),
      "same seed gives the same thresholds")

# 10. Unknown options are rejected.
bad <- run_cli(c("--demo", "--colour", "red"))
check(bad$status != 0L, "unknown option exits with an error")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
