#!/usr/bin/env Rscript
# Generate the claw-genomic-prediction demo data (synthetic, fixed seed).
#
#   Rscript demo/make_demo.R
#
# Eight inbred founders with random SNP alleles are crossed in seven
# biparental families. Recombinant inbred lines are derived by single-seed
# descent to F6 with a forward simulation of meiosis (Haldane map, no
# interference). Genetic values come from 60 additive QTL among the SNPs.
# Six families (240 lines) have phenotypes; the seventh (30 lines) is an
# unphenotyped set of selection candidates.
#
# Files written next to this script:
#   genotypes.csv          id + 500 SNPs coded 0/1/2 (count of allele B)
#   phenotypes.csv         id, family, yield_kg_ha (entry-mean BLUEs), h2 = 0.5
#   phenotypes_multienv.csv id, family, env, yield_kg_ha (BLUEs per environment)
#   truth.csv              id, family, true genetic value (for checking only)

# One meiosis on one chromosome: gamete from haplotypes h1, h2 at marker
# positions pos (cM) on a chromosome of length L cM (Haldane, no interference).
meiosis <- function(h1, h2, pos, L) {
  n_xo <- stats::rpois(1L, L / 100)
  first <- sample.int(2L, 1L)
  if (n_xo == 0L) return(if (first == 1L) h1 else h2)
  k <- findInterval(pos, sort(stats::runif(n_xo, 0, L)))
  if (first == 1L) ifelse(k %% 2L == 0L, h1, h2) else ifelse(k %% 2L == 0L, h2, h1)
}

# Recombinant inbred lines by single-seed descent from the cross p1 x p2.
# p1, p2: lists (one element per chromosome) of 0/1 allele vectors.
# Returns a list of lines, each a list over chromosomes of 2-row matrices.
make_rils <- function(p1, p2, n_lines, map, chr_len, generations = 6L) {
  lapply(seq_len(n_lines), function(l) {
    lapply(seq_along(map), function(c) {
      pos <- map[[c]]
      ind <- rbind(p1[[c]], p2[[c]])                  # F1
      for (g in seq_len(generations - 1L))            # F2 ... F_generations
        ind <- rbind(meiosis(ind[1, ], ind[2, ], pos, chr_len),
                     meiosis(ind[1, ], ind[2, ], pos, chr_len))
      ind
    })
  })
}

make_demo <- function(out_dir = dirname(this_file()), seed = 20261015L) {
  set.seed(seed)
  n_chr <- 5L
  chr_len <- 100
  n_mark <- 100L                                   # per chromosome
  map <- replicate(n_chr, sort(stats::runif(n_mark, 0, chr_len)), simplify = FALSE)
  n_founder <- 8L
  freq <- lapply(seq_len(n_chr), function(c) stats::runif(n_mark, 0.2, 0.8))
  founders <- lapply(seq_len(n_founder), function(f)
    lapply(seq_len(n_chr), function(c) as.integer(stats::runif(n_mark) < freq[[c]])))
  crosses <- rbind(c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(1, 3), c(2, 5), c(4, 7))
  n_per <- c(40, 40, 40, 40, 40, 40, 30)
  lines <- list()
  fam <- character()
  for (k in seq_len(nrow(crosses))) {
    rils <- make_rils(founders[[crosses[k, 1]]], founders[[crosses[k, 2]]],
                      n_per[k], map, chr_len)
    lines <- c(lines, rils)
    fam <- c(fam, rep(sprintf("F%d", k), n_per[k]))
  }
  # Genotype matrix: count of allele "1" over the two homologues (0/1/2).
  X <- t(vapply(lines, function(l) unlist(lapply(l, colSums)), numeric(n_chr * n_mark)))
  ids <- sprintf("RIL%03d", seq_len(nrow(X)))
  mk <- sprintf("snp_%d_%03d", rep(seq_len(n_chr), each = n_mark),
                rep(seq_len(n_mark), n_chr))
  dimnames(X) <- list(ids, mk)
  # Additive genetic values from 60 QTL among the SNPs.
  qtl <- sort(sample(ncol(X), 60L))
  eff <- stats::rnorm(60L)
  g <- drop(scale(X[, qtl], scale = FALSE) %*% eff)
  g <- 250 * (g - mean(g)) / stats::sd(g)          # genetic SD 250 kg/ha
  # Entry-mean BLUEs with h2 = 0.5 (error SD 250 kg/ha).
  y <- 3500 + g + stats::rnorm(length(g), 0, 250)
  phen <- fam != "F7"                               # family 7 = unphenotyped candidates
  geno_df <- data.frame(id = ids, X, check.names = FALSE)
  utils::write.csv(geno_df, file.path(out_dir, "genotypes.csv"), row.names = FALSE)
  utils::write.csv(data.frame(id = ids[phen], family = fam[phen],
                              yield_kg_ha = round(y[phen], 1)),
                   file.path(out_dir, "phenotypes.csv"), row.names = FALSE)
  # Three environments: main effect + G x E (genetic correlation 0.8 between
  # environments) + error; BLUE per line and environment.
  envs <- c("E1", "E2", "E3")
  env_mean <- c(E1 = 3000, E2 = 3500, E3 = 4100)
  me <- do.call(rbind, lapply(envs, function(e) {
    ge <- stats::rnorm(length(g), 0, 125)
    data.frame(id = ids, family = fam, env = e,
               yield_kg_ha = round(env_mean[[e]] + g + ge + stats::rnorm(length(g), 0, 300), 1))
  }))
  me <- me[me$family != "F7", ]
  utils::write.csv(me, file.path(out_dir, "phenotypes_multienv.csv"), row.names = FALSE)
  utils::write.csv(data.frame(id = ids, family = fam, true_genetic_value = round(g, 3)),
                   file.path(out_dir, "truth.csv"), row.names = FALSE)
  invisible(list(X = X, g = g, y = y, family = fam))
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
