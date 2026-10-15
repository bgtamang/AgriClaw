#!/usr/bin/env Rscript
# Generate the claw-qtl-mapping demo data (synthetic, fixed seed), and provide
# the forward simulators used by the tests.
#
#   Rscript demo/make_demo.R
#
# A recombinant inbred line (RIL) population by selfing is simulated from two
# inbred parents (A, B) by single-seed descent to F8 with an explicit model of
# meiosis (Haldane map, no crossover interference). Residual heterozygous calls
# are written as missing ("-"), as genotyping pipelines do for RILs, and 1% of
# calls are set missing at random.
#
# Files written next to this script (CSV, comma-separated):
#   geno.csv   id + 150 markers coded A / B / - (missing)
#   gmap.csv   marker, chr, pos (cM)
#   pmap.csv   marker, chr, pos (Mb; invented physical positions)
#   pheno.csv  id, days_to_flowering, plant_height_cm
#   truth.csv  the simulated QTL (chromosome, position, additive effect)

# Gamete from homologues h1, h2 (allele labels at positions pos, cM) on a
# chromosome of length L cM: Poisson number of crossovers at uniform positions.
meiosis <- function(h1, h2, pos, L) {
  n_xo <- stats::rpois(1L, L / 100)
  first <- sample.int(2L, 1L)
  if (n_xo == 0L) return(if (first == 1L) h1 else h2)
  k <- findInterval(pos, sort(stats::runif(n_xo, 0, L)))
  if (first == 1L) ifelse(k %% 2L == 0L, h1, h2) else ifelse(k %% 2L == 0L, h2, h1)
}

# Founder-label haplotypes of RILs by selfing. founders_per_line: list (one per
# line) of integer founder labels in funnel order (2 for a biparental cross,
# 4 for a 4-way cross [(1 x 2) x (3 x 4)]). Returns, for each chromosome, an
# n_lines x n_positions matrix of founder labels; NA = still heterozygous.
simulate_ril_labels <- function(n_lines, pos_list, chr_len, generations = 8L,
                                n_founders = 2L, funnel = NULL) {
  lapply(seq_along(pos_list), function(c) {
    pos <- pos_list[[c]]
    L <- chr_len[[c]]
    t(vapply(seq_len(n_lines), function(l) {
      f <- if (is.null(funnel)) seq_len(n_founders) else funnel[[l]]
      hom <- function(a) rep(a, length(pos))
      if (n_founders == 2L) {
        ind <- rbind(hom(f[1]), hom(f[2]))                      # F1
      } else {
        g1 <- meiosis(hom(f[1]), hom(f[2]), pos, L)             # gamete of F1 (1 x 2)
        g2 <- meiosis(hom(f[3]), hom(f[4]), pos, L)             # gamete of F1 (3 x 4)
        ind <- rbind(g1, g2)                                    # 4-way hybrid
      }
      for (g in seq_len(generations - 1L))                      # selfing
        ind <- rbind(meiosis(ind[1, ], ind[2, ], pos, L), meiosis(ind[1, ], ind[2, ], pos, L))
      ifelse(ind[1, ] == ind[2, ], ind[1, ], NA_integer_)
    }, integer(length(pos))))
  })
}

make_demo <- function(out_dir = dirname(this_file()), seed = 20261015L) {
  set.seed(seed)
  n <- 160L
  n_chr <- 5L
  chr_len <- rep(100, n_chr)
  n_mark <- 30L
  # Markers about every 3.4 cM with jitter, like a SNP array (no large gaps).
  gpos <- lapply(seq_len(n_chr), function(c) {
    p <- seq(0, 100, length.out = n_mark)
    p[2:(n_mark - 1L)] <- p[2:(n_mark - 1L)] + stats::runif(n_mark - 2L, -1.2, 1.2)
    round(p, 2)
  })
  # QTL positions on a fine grid together with the markers.
  qtl <- data.frame(trait = c("days_to_flowering", "days_to_flowering", "plant_height_cm"),
                    chr = c(1L, 3L, 2L), pos_cM = c(45.3, 70.8, 21.6),
                    additive_effect = c(2.0, -1.4, 3.0))
  pos_all <- lapply(seq_len(n_chr), function(c) sort(c(gpos[[c]], qtl$pos_cM[qtl$chr == c])))
  lab <- simulate_ril_labels(n, pos_all, chr_len, generations = 8L)
  ids <- sprintf("RIL%03d", seq_len(n))
  geno <- list()
  gmap <- list()
  for (c in seq_len(n_chr)) {
    is_m <- pos_all[[c]] %in% gpos[[c]] & !duplicated(pos_all[[c]])
    g <- lab[[c]][, is_m, drop = FALSE]
    nm <- sprintf("m%d_%02d", c, seq_len(ncol(g)))
    geno[[c]] <- matrix(ifelse(is.na(g), "-", ifelse(g == 1L, "A", "B")), n, ncol(g),
                        dimnames = list(ids, nm))
    gmap[[c]] <- data.frame(marker = nm, chr = c, pos = pos_all[[c]][is_m])
  }
  G <- do.call(cbind, geno)
  G[matrix(stats::runif(length(G)) < 0.01, nrow(G))] <- "-"
  gm <- do.call(rbind, gmap)
  # Genotype (+1 for B, -1 for A) at each QTL; residual heterozygotes count 0.
  qg <- vapply(seq_len(nrow(qtl)), function(k) {
    c <- qtl$chr[k]
    j <- match(qtl$pos_cM[k], pos_all[[c]])
    x <- lab[[c]][, j]
    ifelse(is.na(x), 0, ifelse(x == 2L, 1, -1))
  }, numeric(n))
  dtf <- 55 + qg[, 1] * qtl$additive_effect[1] + qg[, 2] * qtl$additive_effect[2] +
    stats::rnorm(n, 0, 2.2)
  ht <- 80 + qg[, 3] * qtl$additive_effect[3] + stats::rnorm(n, 0, 8)
  dtf[sample(n, 3)] <- NA
  utils::write.csv(data.frame(id = ids, G, check.names = FALSE),
                   file.path(out_dir, "geno.csv"), row.names = FALSE, quote = FALSE)
  utils::write.csv(gm, file.path(out_dir, "gmap.csv"), row.names = FALSE, quote = FALSE)
  # Invented physical map: about 0.45 Mb per cM with local variation.
  pm <- gm
  pm$pos <- round(unlist(lapply(split(gm$pos, gm$chr), function(p)
    cumsum(c(0.2, diff(p) * stats::runif(length(p) - 1L, 0.2, 0.7))))), 3)
  utils::write.csv(pm, file.path(out_dir, "pmap.csv"), row.names = FALSE, quote = FALSE)
  utils::write.csv(data.frame(id = ids, days_to_flowering = round(dtf, 1),
                              plant_height_cm = round(ht, 1)),
                   file.path(out_dir, "pheno.csv"), row.names = FALSE, quote = FALSE, na = "NA")
  utils::write.csv(qtl, file.path(out_dir, "truth.csv"), row.names = FALSE, quote = FALSE)
  invisible(list(geno = G, gmap = gm, qtl = qtl))
}

this_file <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) return(normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/"))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(of, winslash = "/"))
  }
  file.path(getwd(), "make_demo.R")
}

if (sys.nframe() == 0L) {
  make_demo()
  cat("Demo files written to", dirname(this_file()), "\n")
}
