# Generate the claw-ld-decay demo data: inbred (doubled-haploid) lines sampled
# from a forward-simulated random-mating population founded by a few inbred
# founders. Each chromosome has its own known recombination rate, so LD must
# decay faster (shorter distance) on chromosomes with more recombination.
#
#   Rscript demo/make_demo.R
#
# Sourcing this file only defines simulate_ld_population(); files are written
# only when it is run with Rscript.

# Haploid forward simulation: every generation, each of pop_size gametes is a
# recombinant of two random parental gametes, with Poisson(L) crossovers placed
# uniformly (Haldane map, no interference). Gametes are tracked as founder IDs
# at each marker; alleles come from a random non-trivial split of the founders.
simulate_ld_population <- function(rates_cM_per_Mb = c(0.5, 1, 2), chr_len_bp = 10e6,
                                   m_per_chr = 300, n_founders = 3, pop_size = 400,
                                   generations = 20, n_lines = 200, seed = 1) {
  set.seed(seed)
  geno <- list(); map <- list()
  for (k in seq_along(rates_cM_per_Mb)) {
    pos <- sort(sample.int(chr_len_bp, m_per_chr))
    L <- rates_cM_per_Mb[k] * chr_len_bp / 1e6 / 100      # Morgans per meiosis
    H <- matrix(rep(sample.int(n_founders, pop_size, replace = TRUE), m_per_chr), pop_size, m_per_chr)
    for (g in seq_len(generations)) {
      H2 <- H
      for (r in seq_len(pop_size)) {
        par <- sample.int(pop_size, 2L)
        nx <- stats::rpois(1L, L)
        if (nx == 0L) { H2[r, ] <- H[par[1], ]; next }
        side <- findInterval(pos, sort(stats::runif(nx, 0, chr_len_bp))) %% 2L
        H2[r, ] <- ifelse(side == 0L, H[par[1], ], H[par[2], ])
      }
      H <- H2
    }
    lines <- H[sample.int(pop_size, n_lines), , drop = FALSE]
    FA <- matrix(0L, n_founders, m_per_chr)
    for (j in seq_len(m_per_chr)) {
      repeat { a <- stats::rbinom(n_founders, 1L, 0.5); if (any(a == 0L) && any(a == 1L)) break }
      FA[, j] <- a
    }
    Gk <- matrix(2L * FA[cbind(as.vector(lines), rep(seq_len(m_per_chr), each = n_lines))], n_lines, m_per_chr)
    geno[[k]] <- Gk
    map[[k]] <- data.frame(marker = sprintf("S%d_%d", k, pos), chr = as.character(k), pos = pos,
                           rate_cM_per_Mb = rates_cM_per_Mb[k], stringsAsFactors = FALSE)
  }
  G <- do.call(cbind, geno)
  mp <- do.call(rbind, map)
  dimnames(G) <- list(sprintf("DH%03d", seq_len(n_lines)), mp$marker)
  list(G = G, map = mp)
}

if (sys.nframe() == 0L) {
  here <- local({
    f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    dirname(normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/"))
  })
  sim <- simulate_ld_population(rates_cM_per_Mb = c(0.5, 1, 2), n_lines = 150, seed = 7)
  G <- t(sim$G)
  set.seed(8)
  G[sample(length(G), round(0.01 * length(G)))] <- NA
  tab <- data.frame(marker = sim$map$marker, chr = sim$map$chr, pos = sim$map$pos, G, check.names = FALSE)
  con <- gzfile(file.path(here, "demo_ld_012.csv.gz"), "w")
  utils::write.csv(tab, con, row.names = FALSE, na = "")
  close(con)
  utils::write.csv(unique(sim$map[, c("chr", "rate_cM_per_Mb")]), file.path(here, "demo_truth.csv"), row.names = FALSE)
  cat("Wrote demo_ld_012.csv.gz and demo_truth.csv to", here, "\n")
}
