# Generate the claw-marker-qc demo data: a synthetic inbred panel with planted
# problems, written as demo_panel.vcf.gz, plus demo_truth.csv listing every
# planted problem so the tests can check that each one is detected exactly.
#
#   Rscript demo/make_demo.R
#
# Sourcing this file only defines simulate_qc_panel(); files are written only
# when it is run with Rscript.

simulate_qc_panel <- function(seed = 42, n = 120, n_chr = 3, m_per_chr = 200,
                              miss_rate = 0.02, err_rate = 0.003) {
  set.seed(seed)
  m <- n_chr * m_per_chr
  chr <- rep(seq_len(n_chr), each = m_per_chr)
  pos <- unlist(lapply(seq_len(n_chr), function(i) sort(sample.int(40e6, m_per_chr))))
  marker <- sprintf("S%d_%d", chr, pos)
  bases <- c("A", "C", "G", "T")
  ref <- sample(bases, m, replace = TRUE)
  alt <- vapply(ref, function(b) sample(setdiff(bases, b), 1), "")
  samples <- sprintf("G%03d", seq_len(n))
  # Two subpopulations; ALT frequencies kept in [0.2, 0.8] so that no marker is
  # rare by chance (only the planted ones are).
  pop <- rep(1:2, length.out = n)
  p0 <- stats::runif(m, 0.3, 0.7)
  pk <- cbind(pmin(pmax(p0 + stats::rnorm(m, 0, 0.1), 0.2), 0.8),
              pmin(pmax(p0 + stats::rnorm(m, 0, 0.1), 0.2), 0.8))
  G <- matrix(0L, m, n, dimnames = list(marker, samples))
  for (j in seq_len(n)) G[, j] <- 2L * stats::rbinom(m, 1, pk[, pop[j]])
  truth <- data.frame(problem = character(), id = character(), stringsAsFactors = FALSE)
  add <- function(problem, id) truth <<- rbind(truth, data.frame(problem = problem, id = id))

  # Planted sample problems (indices chosen without overlap).
  idx <- sample.int(n)
  dup_src <- idx[1:3]; dup_cp <- idx[4:6]
  contam <- idx[7:9]; contam_b <- idx[10:12]
  lowcall <- idx[13:16]
  for (k in 1:3) {            # contamination: a mixture of two lines is heterozygous where they differ
    a <- G[, contam[k]]; b <- G[, contam_b[k]]
    G[, contam[k]] <- ifelse(a == b, a, 1L)
    add("sample_heterozygosity", samples[contam[k]])
  }
  for (k in 1:3) {            # duplicates: copy with 0.5% discordant calls
    x <- G[, dup_src[k]]
    flip <- sample.int(m, round(0.005 * m))
    x[flip] <- 2L - x[flip]
    G[, dup_cp[k]] <- x
    pr <- sort(samples[c(dup_src[k], dup_cp[k])])
    add("duplicate_pair", paste(pr, collapse = "~"))
  }
  # Planted marker problems.
  midx <- sample.int(m)
  mono <- midx[1:10]; rare <- midx[11:20]; highmiss <- midx[21:40]
  G[mono, ] <- 0L
  for (i in rare) { G[i, ] <- 0L; G[i, sample.int(n, 3)] <- 2L }
  for (i in mono) add("marker_monomorphic", marker[i])
  for (i in rare) add("marker_low_maf", marker[i])
  # Genotyping error: a few heterozygous calls and allele flips.
  err <- matrix(stats::runif(m * n) < err_rate, m, n)
  err[c(mono, rare), ] <- FALSE
  G[err] <- 1L
  # Missing data: background, high-missing markers, low-call-rate samples.
  miss <- matrix(stats::runif(m * n) < miss_rate, m, n)
  miss[highmiss, ] <- matrix(stats::runif(length(highmiss) * n) < 0.4, length(highmiss), n)
  miss[, lowcall] <- matrix(stats::runif(m * length(lowcall)) < 0.5, m, length(lowcall))
  miss[rare, ] <- FALSE
  G[miss] <- NA
  for (i in highmiss) add("marker_call_rate", marker[i])
  for (j in lowcall) add("sample_call_rate", samples[j])
  map <- data.frame(marker = marker, chr = as.character(chr), pos = pos, ref = ref, alt = unname(alt),
                    stringsAsFactors = FALSE)
  list(geno = G, map = map, truth = truth)
}

if (sys.nframe() == 0L) {
  here <- local({
    f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    dirname(normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/"))
  })
  source(file.path(dirname(here), "scripts", "marker_qc.R"))  # for write_vcf()
  sim <- simulate_qc_panel()
  write_vcf(sim$geno, sim$map, file.path(here, "demo_panel.vcf.gz"),
            source = "claw-marker-qc demo/make_demo.R (synthetic, seed 42)")
  utils::write.csv(sim$truth, file.path(here, "demo_truth.csv"), row.names = FALSE)
  cat("Wrote demo_panel.vcf.gz and demo_truth.csv to", here, "\n")
}
