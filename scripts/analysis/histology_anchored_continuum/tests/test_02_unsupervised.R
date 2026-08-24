#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_path <- file.path(
  Sys.getenv("MASLD_PROJECT_ROOT", getwd()),
  "scripts/analysis/histology_anchored_continuum/02_build_unsupervised_axes.R"
)
source(script_path)

set.seed(HAC_SEED)

# The released row-sum criterion is strict: a sum equal to n samples is removed.
counts <- matrix(rpois(1005L * 6L, lambda = 5), nrow = 1005L, ncol = 6L)
rownames(counts) <- sprintf("ENSG%011d", seq_len(nrow(counts)))
colnames(counts) <- paste0("S", seq_len(ncol(counts)))
counts[nrow(counts), ] <- 1L
normalised <- normalise_cohort(counts, "SYNTHETIC")
stopifnot(nrow(normalised) == 1004L)
stopifnot(!rownames(counts)[nrow(counts)] %in% rownames(normalised))
stopifnot(identical(colnames(normalised), colnames(counts)))

# Named frozen vectors must retain their names after numeric coercion.
named <- setNames(c(1, 2, 3), c("ENSG1.1", "ENSG2.2", "ENSG3.3"))
coerced <- as_frozen_vector(named, "synthetic")
stopifnot(identical(names(coerced), c("ENSG1", "ENSG2", "ENSG3")))

# A reversed frozen loading must reverse the arbitrary PCA sign and yield a
# positive post-orientation cosine without consulting any phenotype.
expression <- matrix(rnorm(145L * 12L), nrow = 145L, ncol = 12L)
rownames(expression) <- sprintf("ENSG%011d", seq_len(nrow(expression)))
colnames(expression) <- paste0("P", seq_len(ncol(expression)))
reference_fit <- prcomp(t(expression), center = TRUE, scale. = FALSE, rank. = 2L)
reference_loading <- reference_fit$rotation[, 1L]
names(reference_loading) <- rownames(reference_fit$rotation)
aligned <- fit_oriented_pc1(
  expression, -reference_loading, "signature_pc1", "SYNTHETIC"
)
stopifnot(aligned$orientation_sign == -1)
stopifnot(aligned$oriented_cosine > 1 - 1e-12)
stopifnot(cor(aligned$score, reference_fit$x[, 1L]) < -1 + 1e-12)

# Percentile ranks are bounded, deterministic, and average exact ties.
pct <- percentile_rank(c(1, 1, 3, 4))
stopifnot(isTRUE(all.equal(pct, c(1 / 6, 1 / 6, 2 / 3, 1))))

cat("test_02_unsupervised: PASS\n")
