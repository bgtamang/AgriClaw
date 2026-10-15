# Generate the claw-gwas-crop demo data: a synthetic structured panel of inbred
# lines with local LD, and two traits:
#   height_cm   3 QTNs (about 10% of the phenotypic variance each) + polygenic
#               background + subpopulation effect + noise;
#   protein_pct no marker effects at all: subpopulation effect + noise, so any
#               hit is a false positive and lambda should be close to 1.
#
#   Rscript demo/make_demo.R
#
# Sourcing this file only defines the simulation functions; files are written
# only when it is run with Rscript.

# Inbred lines in K subpopulations. Founder haplotypes per subpopulation are
# drawn from Balding-Nichols allele frequencies (Fst); each line is a mosaic of
# its subpopulation's founders with Poisson breakpoints, which creates LD that
# decays with distance.
simulate_panel <- function(n_per_pop = c(100, 100, 100), n_chr = 5, m_per_chr = 400,
                           chr_len = 50e6, fst = 0.15, n_founders = 12,
                           breaks_per_chr = 8, seed = 1) {
  set.seed(seed)
  K <- length(n_per_pop); n <- sum(n_per_pop); m <- n_chr * m_per_chr
  pop <- rep(seq_len(K), n_per_pop)
  chr <- rep(seq_len(n_chr), each = m_per_chr)
  pos <- unlist(lapply(seq_len(n_chr), function(i) sort(sample.int(chr_len, m_per_chr))))
  p0 <- stats::runif(m, 0.05, 0.95)
  a <- (1 - fst) / fst
  pk <- sapply(seq_len(K), function(k) stats::rbeta(m, p0 * a, (1 - p0) * a))
  G <- matrix(0L, n, m)
  for (k in seq_len(K)) {
    fh <- matrix(stats::rbinom(n_founders * m, 1, rep(pk[, k], each = n_founders)), n_founders, m)
    for (i in which(pop == k)) {
      for (ch in seq_len(n_chr)) {
        idx <- which(chr == ch)
        nb <- stats::rpois(1, breaks_per_chr)
        bp <- sort(stats::runif(nb, 0, chr_len))
        seg <- findInterval(pos[idx], bp) + 1L
        f <- sample.int(n_founders, nb + 1L, replace = TRUE)
        G[i, idx] <- 2L * fh[cbind(f[seg], idx)]
      }
    }
  }
  dimnames(G) <- list(sprintf("L%03d", seq_len(n)), sprintf("S%d_%d", chr, pos))
  af <- colMeans(G) / 2
  keep <- pmin(af, 1 - af) >= 0.05          # drop monomorphic and rare markers
  list(G = G[, keep], map = data.frame(marker = colnames(G)[keep], chr = as.character(chr[keep]),
                                       pos = pos[keep], stringsAsFactors = FALSE), pop = pop)
}

# Traits with known architecture (variance shares of the standardized trait).
simulate_traits <- function(panel, seed = 1, qtn_share = 0.10, poly_share = 0.20,
                            pop_share = 0.15, null_pop_share = 0.30, qtn_chr = c("1", "3", "5")) {
  set.seed(seed + 1000)
  G <- panel$G; map <- panel$map; pop <- panel$pop; n <- nrow(G)
  z <- function(v) (v - mean(v)) / stats::sd(v)
  af <- colMeans(G) / 2
  pk <- sapply(sort(unique(pop)), function(k) colMeans(G[pop == k, , drop = FALSE]) / 2)
  dif <- apply(pk, 1, function(x) diff(range(x)))
  ok <- which(pmin(af, 1 - af) > 0.3 & dif < 0.25)   # common, weakly differentiated QTNs
  qtn <- vapply(qtn_chr, function(ch) { cand <- ok[map$chr[ok] == ch]; cand[sample.int(length(cand), 1)] }, 1L)
  poly <- sample(setdiff(seq_len(ncol(G)), qtn), 200)
  gpoly <- z(drop(scale(G[, poly]) %*% stats::rnorm(200)))
  popeff <- z(c(1, 0, -1)[pop])
  qsum <- Reduce(`+`, lapply(qtn, function(j) z(G[, j])))
  y_qtn <- sqrt(qtn_share) * qsum + sqrt(poly_share) * gpoly + sqrt(pop_share) * popeff +
    sqrt(1 - length(qtn) * qtn_share - poly_share - pop_share) * z(stats::rnorm(n))
  y_null <- sqrt(null_pop_share) * popeff + sqrt(1 - null_pop_share) * z(stats::rnorm(n))
  list(y_qtn = y_qtn, y_null = y_null, qtn = map$marker[qtn],
       qtn_effect_sign = sign(vapply(qtn, function(j) stats::cor(G[, j], y_qtn), 0)))
}

if (sys.nframe() == 0L) {
  here <- local({
    f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    dirname(normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/"))
  })
  panel <- simulate_panel(n_per_pop = c(84, 83, 83), m_per_chr = 300, seed = 2)
  tr <- simulate_traits(panel, seed = 2)
  G <- t(panel$G)
  set.seed(3)
  G[sample(length(G), round(0.01 * length(G)))] <- NA     # 1% missing calls
  bases <- c("A", "C", "G", "T")
  ref <- sample(bases, nrow(G), TRUE)
  alt <- vapply(ref, function(b) sample(setdiff(bases, b), 1), "")
  geno <- data.frame(marker = panel$map$marker, chr = panel$map$chr, pos = panel$map$pos,
                     ref = ref, alt = unname(alt), G, check.names = FALSE)
  con <- gzfile(file.path(here, "demo_geno_012.csv.gz"), "w")
  utils::write.csv(geno, con, row.names = FALSE, na = "")
  close(con)
  pheno <- data.frame(genotype = rownames(panel$G),
                      height_cm = round(85 + 10 * tr$y_qtn, 2),
                      protein_pct = round(40 + 1.5 * tr$y_null, 3))
  utils::write.csv(pheno, file.path(here, "demo_pheno_blues.csv"), row.names = FALSE)
  utils::write.csv(data.frame(trait = "height_cm", qtn = tr$qtn,
                              chr = panel$map$chr[match(tr$qtn, panel$map$marker)],
                              pos = panel$map$pos[match(tr$qtn, panel$map$marker)]),
                   file.path(here, "demo_truth.csv"), row.names = FALSE)
  cat("Wrote demo_geno_012.csv.gz, demo_pheno_blues.csv and demo_truth.csv to", here, "\n")
}
