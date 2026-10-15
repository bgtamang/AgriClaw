#!/usr/bin/env Rscript
# Generate the claw-selection-index demo inputs (synthetic, illustrative only).
#
#   Rscript demo/make_demo.R
#
# The phenotypic (P) and genetic (G) covariance matrices are built from
# invented standard deviations, heritabilities and correlations for four
# soybean-like traits. They are NOT estimates from any real population; they
# only exercise the script. The generator checks that P is positive definite,
# G is positive semi-definite and P - G is positive semi-definite, so every
# implied heritability (including for linear combinations) lies in [0, 1].

make_demo <- function(out_dir = dirname(this_file())) {
  traits <- c("yield_kg_ha", "protein_pct", "oil_pct", "height_cm")
  sd_p <- c(300, 1.2, 0.8, 8)          # phenotypic SDs, trait units
  h2 <- c(0.45, 0.70, 0.65, 0.60)      # heritabilities (entry-mean basis)
  r_g <- matrix(c( 1.00, -0.40,  0.20,  0.30,
                  -0.40,  1.00, -0.60, -0.10,
                   0.20, -0.60,  1.00,  0.00,
                   0.30, -0.10,  0.00,  1.00), 4, 4, byrow = TRUE)
  r_p <- matrix(c( 1.00, -0.25,  0.15,  0.25,
                  -0.25,  1.00, -0.45, -0.05,
                   0.15, -0.45,  1.00,  0.00,
                   0.25, -0.05,  0.00,  1.00), 4, 4, byrow = TRUE)
  sd_g <- sqrt(h2) * sd_p
  P <- diag(sd_p) %*% r_p %*% diag(sd_p)
  G <- diag(sd_g) %*% r_g %*% diag(sd_g)
  dimnames(P) <- dimnames(G) <- list(traits, traits)
  stopifnot(min(eigen(P, symmetric = TRUE)$values) > 0,
            min(eigen(G, symmetric = TRUE)$values) > 0,
            min(eigen(P - G, symmetric = TRUE)$values) > 0)

  write_matrix <- function(M, file) {
    df <- data.frame(trait = rownames(M), round(M, 6), check.names = FALSE)
    utils::write.csv(df, file.path(out_dir, file), row.names = FALSE)
  }
  write_matrix(P, "P_matrix.csv")
  write_matrix(G, "G_matrix.csv")

  # Relative economic weights per trait unit (invented): value of +1 kg/ha
  # yield, +1 % protein, +1 % oil; height has no economic value.
  utils::write.csv(data.frame(trait = traits, weight = c(1, 60, 40, 0)),
                   file.path(out_dir, "weights.csv"), row.names = FALSE)
  # Desired gains per cycle, trait units (Pesek & Baker index).
  utils::write.csv(data.frame(trait = traits, desired_gain = c(60, 0.3, 0, 0)),
                   file.path(out_dir, "desired_gains.csv"), row.names = FALSE)
  # Breeder's-equation scenarios (invented): phenotypic selection on a
  # 4-year cycle (given h2, a small population, or variance components from a
  # multi-environment trial with 4 environments x 2 replicates) versus genomic
  # selection on a 2-year cycle.
  sc <- data.frame(
    scenario = c("phenotypic_4yr", "phenotypic_4yr_small", "phenotypic_4yr_varcomp", "genomic_2yr"),
    proportion_selected = c(0.10, 0.10, 0.10, 0.05),
    n_candidates = c(NA, 30, NA, NA),
    h2 = c(0.45, 0.45, NA, NA),
    sigma_p = c(300, 300, NA, NA),
    var_g = c(NA, NA, 40000, NA),
    var_ge = c(NA, NA, 30000, NA),
    var_e = c(NA, NA, 160000, NA),
    n_env = c(NA, NA, 4, NA),
    n_rep = c(NA, NA, 2, NA),
    accuracy = c(NA, NA, NA, 0.50),
    sigma_a = c(NA, NA, NA, 201.2),
    cycle_years = c(4, 4, 4, 2))
  utils::write.csv(sc, file.path(out_dir, "scenarios.csv"), row.names = FALSE,
                   na = "")
  invisible(list(P = P, G = G))
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
