# Tests for claw-ld-decay. Run from anywhere: Rscript tests/test_ld_decay.R
# Known answers: r2 equals cor()^2; the Hill-Weir formula against hand
# arithmetic; recovery of a known rho from data generated on the Hill-Weir
# curve; exact pair enumeration and uniform pair sampling; the analytical
# background r2 of independent markers; and forward-simulated chromosomes with
# known recombination rates, where decay distance must shrink as the rate grows.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "ld_decay.R"))   # loads functions; main() does not run
source(file.path(SKILL_TEST_DIR, "demo", "make_demo.R"))      # simulate_ld_population()

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "ld_decay.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
tmp <- file.path(tempdir(), "claw-ld-decay-test")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE)

# 1. r2 equals cor()^2 (pairwise-complete observations) -------------------------------------
set.seed(1)
X <- matrix(sample(c(0, 1, 2), 40 * 6, TRUE), 40, 6)
X[sample(length(X), 20)] <- NA
I <- c(1, 1, 2, 3, 4); J <- c(2, 3, 5, 6, 6)
rr <- r2_pairs(X, I, J)
ref <- mapply(function(i, j) stats::cor(X[, i], X[, j], use = "complete.obs")^2, I, J)
nref <- mapply(function(i, j) sum(!is.na(X[, i]) & !is.na(X[, j])), I, J)
check(max(abs(rr$r2 - ref)) < 1e-12, "r2 equals cor(use = 'complete.obs')^2")
check(all(rr$n_obs == nref), "n_obs counts samples called at both markers")
rr2 <- r2_pairs(X, I, J, chunk = 2L)
check(identical(rr$r2, rr2$r2), "r2 does not depend on the chunk size")
Xm <- cbind(X[, 1], 2)
check(is.na(r2_pairs(Xm, 1, 2)$r2), "r2 with a monomorphic marker is NA, not 0")

# 2. Hill-Weir expectation against hand arithmetic ---------------------------------------------
check(abs(hw_expect(0, 1, 1e12) - 10 / 22) < 1e-10, "E(r2) at C = 0 and large n is 10/22")
check(abs(hw_expect(0, 1, 100) - (10 / 22) * (1 + 36 / 2200)) < 1e-12, "E(r2) at C = 0, n = 100")
check(abs(hw_expect(5, 1, 50) - (15 / 112) * (1 + 776 / 5600)) < 1e-12, "E(r2) at C = 5, n = 50")
check(abs(hw_expect(2e6, 2.5e-6, 50) - hw_expect(5, 1, 50)) < 1e-12, "C = rho x distance")
check(all(diff(hw_expect(seq(0, 1e4, by = 10), 1, 100)) < 0), "E(r2) decreases with distance")

# 3. Fit recovers a known rho; distances solve E(r2) = target ------------------------------------
set.seed(2)
rho_true <- 2e-6; n <- 120
d <- stats::runif(5000, 1e3, 5e6)
r2_sim <- hw_expect(d, rho_true, n) + stats::rnorm(5000, 0, 0.05)
fit <- fit_hw(d, r2_sim, n)
check(abs(fit$rho / rho_true - 1) < 0.05 && !fit$at_boundary, sprintf("fitted rho %.3g within 5%% of %.3g", fit$rho, rho_true))
d02 <- hw_distance(0.2, rho_true, n)
check(abs(hw_expect(d02, rho_true, n) - 0.2) < 1e-8, "distance at r2 = 0.2 solves the Hill-Weir equation")
check(is.na(hw_distance(0.5 / n, rho_true, n)), "a target below the curve's floor (about 1/n) gives NA")
check(hw_distance(0.9, rho_true, n) == 0, "a target above the curve's maximum gives 0")

# 4. Binned curve and model-free distance (hand-worked) ---------------------------------------------
bins <- data.frame(mean_dist_bp = c(1, 2, 3, 4), mean_r2 = c(0.6, 0.4, 0.3, 0.1))
check(abs(binned_distance(bins, 0.2) - 3.5) < 1e-12, "binned distance interpolates between bins (3.5)")
check(is.na(binned_distance(bins, 0.05)) && is.na(binned_distance(bins, 0.7)),
      "binned distance is NA when the curve never crosses or starts below the target")
b <- bin_r2(c(rep(10, 30), rep(1000, 30)), c(rep(0.8, 30), rep(0.2, 30)), n_bins = 4, min_n = 20)
check(nrow(b) == 2L && all(abs(b$mean_r2 - c(0.8, 0.2)) < 1e-12), "bin means on a two-cluster example")

# 4b. LD window choice: usable Hill-Weir, extrapolation refused, fallback, not estimable --------------------
sw <- data.frame(chr = c("1", "2", "3", "4", "all"), dist_r2_threshold_hw_bp = c(2e6, 5e3, NA, 4e4, 3e3),
                 hw_extrapolated = c(FALSE, TRUE, FALSE, FALSE, TRUE), hw_fit_at_boundary = FALSE,
                 dist_r2_threshold_binned_bp = c(2.2e6, 8e4, 1.5e6, NA, NA), r2_max_binned = c(0.6, 0.5, 0.4, 0.05, 0.04))
w <- ld_window(sw, "hw", 0.2)
check(identical(w$ld_distance_bp, c(2e6, 8e4, 1.5e6, NA, NA)),
      "LD window: Hill-Weir when usable; binned when Hill-Weir is extrapolated or missing; NA when the data never reach r2 = 0.2")
check(grepl("fallback", w$basis[2]) && grepl("not estimable", w$basis[4]) && grepl("max 0.050", w$basis[4]),
      "LD window basis says when a fallback was used and why a window is not estimable")
check(identical(ld_window(sw, "binned", 0.2)$ld_distance_bp, c(2.2e6, 8e4, 1.5e6, NA, NA)), "--window-basis binned prefers binned values")

# 5. Pair sampling: exact enumeration, distance limit, uniform subsampling -----------------------------
pos <- c(500, 100, 300, 200, 900, 400)
sp <- sample_pairs(pos, max_dist = 250, max_pairs = 1000)
bf <- which(outer(pos, pos, function(a, b) b - a > 0 & b - a <= 250), arr.ind = TRUE)
key <- function(i, j) sort(paste(pmin(i, j), pmax(i, j)))
check(identical(key(sp$i, sp$j), key(bf[, 1], bf[, 2])), "all pairs within max-dist are enumerated exactly")
set.seed(3)
pos20 <- sort(sample.int(1e6, 20)); total <- choose(20, 2)
inc <- matrix(0, 20, 20)
for (r in 1:300) { s <- sample_pairs(pos20, Inf, 50); inc[cbind(s$i, s$j)] <- inc[cbind(s$i, s$j)] + 1 }
freq <- inc[upper.tri(inc)] / 300
check(all(abs(freq - 50 / total) < 0.1), "subsampled pairs are uniform: every pair included with p = 50/190")
s1 <- sample_pairs(pos20, 2e5, 30)
check(nrow(s1) <= 30 && !anyDuplicated(paste(s1$i, s1$j)) && all(abs(pos20[s1$j] - pos20[s1$i]) <= 2e5),
      "subsampled pairs are unique and respect max-dist")

# 6. Background r2 of independent markers matches the half-normal expectation ------------------------
set.seed(4)
n_l <- 200
Xi <- matrix(2 * stats::rbinom(n_l * 400, 1, 0.5), n_l, 400)
bg <- background_r2(Xi, rep(1:4, each = 100), n_pairs = 4000)
expected_bg <- (sqrt(2 / (pi * n_l)) + stats::qnorm(0.95) * sqrt((1 - 2 / pi) / n_l))^2
check(abs(bg$r2_p95 / expected_bg - 1) < 0.15,
      sprintf("background r2 %.4f within 15%% of the analytical %.4f for unlinked markers", bg$r2_p95, expected_bg))

# 7. Forward-simulated chromosomes with known recombination rates ------------------------------------
rates <- c(0.5, 1, 2)
G_gen <- 20
theta_half <- 1 - 0.5^(1 / (2 * G_gen))                        # (1 - theta)^(2G) = 1/2
theory_half <- (-0.5 * log(1 - 2 * theta_half)) / (rates * 1e-8)  # Haldane map, bp
opt <- DEFAULTS
opt$max_pairs <- 10000
for (seed in c(4, 9)) {
  sim <- simulate_ld_population(rates_cM_per_Mb = rates, generations = G_gen, seed = seed)
  r <- run_ld(t(sim$G), sim$map, opt)
  pc <- r$per_chr[match(c("1", "2", "3"), r$per_chr$chr), ]
  check(all(diff(pc$dist_r2_threshold_hw_bp) < 0),
        sprintf("seed %d: Hill-Weir r2 = 0.2 distance shrinks as recombination rises (%s Mb)", seed,
                paste(round(pc$dist_r2_threshold_hw_bp / 1e6, 2), collapse = " > ")))
  check(all(diff(pc$dist_r2_threshold_binned_bp) < 0),
        sprintf("seed %d: binned r2 = 0.2 distance shrinks as recombination rises", seed))
  check(all(diff(pc$dist_half_max_binned_bp) < 0), sprintf("seed %d: binned half-max distance shrinks as recombination rises", seed))
  ratio <- pc$dist_r2_threshold_binned_bp[1] / pc$dist_r2_threshold_binned_bp[3]
  check(ratio > 2 && ratio < 8, sprintf("seed %d: 4-fold rate difference gives a %.1f-fold distance difference", seed, ratio))
  scale <- mean(pc$dist_half_max_binned_bp / theory_half)
  check(scale > 0.6 && scale < 1.4,
        sprintf("seed %d: half-max distances average %.2f x the closed-form (1 - theta)^(2G) = 1/2 value", seed, scale))
  check(r$background$r2_p95 < 0.05, sprintf("seed %d: no structure, so background r2 is small (%.3f)", seed, r$background$r2_p95))
}

# 8. End to end on the demo data ------------------------------------------------------------------------
out <- file.path(tmp, "demo")
res <- run_cli(c("--demo", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "tables/ld_decay_summary.csv", "tables/ld_binned.csv", "tables/ld_window.csv",
            "figures/ld_decay_genome.png", "figures/ld_decay_by_chr.png",
            "reproducibility/commands.sh", "reproducibility/sessionInfo.txt")) {
  check(file.exists(file.path(out, f)), paste("demo writes", f))
}
check(length(Sys.glob(file.path(out, "reproducibility", "checksums.*"))) == 1L, "demo writes a checksum file")
win <- utils::read.csv(file.path(out, "tables", "ld_window.csv"), stringsAsFactors = FALSE)
truth <- utils::read.csv(file.path(SKILL_TEST_DIR, "demo", "demo_truth.csv"), stringsAsFactors = FALSE)
wd <- win$ld_distance_bp[match(as.character(truth$chr[order(truth$rate_cM_per_Mb)]), win$chr)]
check(identical(names(win), c("chr", "ld_distance_bp", "basis")) && "all" %in% win$chr,
      "ld_window.csv has chr, ld_distance_bp, basis and a genome-wide 'all' row (claw-gwas-crop format)")
check(all(diff(wd) < 0), "demo: LD window shrinks from the 0.5 to the 2 cM/Mb chromosome")
res2 <- run_cli(c("--demo", "--max-dist", "2000000", "--window-basis", "binned", "--max-pairs", "5000",
                   "--n-unlinked", "1000", "--out", shQuote(file.path(tmp, "demo2"))))
check(res2$status == 0L, "demo with --max-dist and --window-basis binned runs")
s2 <- utils::read.csv(file.path(tmp, "demo2", "tables", "ld_decay_summary.csv"), stringsAsFactors = FALSE)
check(all(s2$max_pair_dist_bp <= 2e6), "--max-dist limits pair distances")
bad <- run_cli(c("--demo", "--r2-threshold", "1.5"))
check(bad$status != 0L, "an r2 threshold outside (0, 1) is rejected")

cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
