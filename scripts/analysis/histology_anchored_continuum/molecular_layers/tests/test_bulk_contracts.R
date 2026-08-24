#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

args <- commandArgs(trailingOnly = FALSE)
test_file <- sub("^--file=", "", args[grepl("^--file=", args)])
stopifnot(length(test_file) == 1L)
script_dir <- dirname(dirname(normalizePath(test_file)))
source(file.path(script_dir, "lib_molecular_layers.R"))
source(file.path(script_dir, "10_bulk_lib.R"))

assert_equal <- function(x, y, tolerance = 1e-10, message = "values differ") {
  if (is.numeric(x) && is.numeric(y)) {
    stopifnot(length(x) == length(y), all(abs(x - y) < tolerance))
  } else {
    stopifnot(identical(x, y))
  }
  invisible(TRUE)
}

# Complete-family BH keeps untestable rows explicit and uses declared n.
p <- c(0.001, NA_real_, 0.02)
assert_equal(
  ml_complete_bh(p, 10L)[c(1L, 3L)],
  p.adjust(p[c(1L, 3L)], method = "BH", n = 10L)
)
stopifnot(is.na(ml_complete_bh(p, 10L)[2L]))

# Fixed and random meta-analysis agree when there is no heterogeneity.
fixed <- ml_fixed_meta(c(1, 1), c(0.2, 0.2))
random <- ml_random_meta(c(1, 1), c(0.2, 0.2))
assert_equal(fixed$beta, 1)
assert_equal(random$beta, 1)
assert_equal(random$tau2, 0)

# The fixed projection LOO score exactly removes the target contribution.
test_genes <- paste0("g", seq_len(101L))
expression <- matrix(
  seq_len(202L), nrow = 101L,
  dimnames = list(test_genes, c("s1", "s2"))
)
projection <- data.table(
  gene_id_base = test_genes,
  discovery_center = rep(0, 101L),
  discovery_loading_oriented = seq_len(101L),
  used_for_fixed_projection = TRUE
)
dummy <- data.table(
  axis_id = "signature_pc1", gene_id_base = test_genes,
  frozen_loading = rep(1, 101L)
)
score <- bulk_loo_axis("fixed_projection", "g2", expression, dummy, projection)
use <- setdiff(test_genes, "g2")
expected <- as.numeric(crossprod(seq_len(101L)[match(use, test_genes)], expression[use, ]))
assert_equal(score, expected)
projection[, discovery_loading_oriented := -discovery_loading_oriented]
flipped_score <- bulk_loo_axis("fixed_projection", "g2", expression, dummy, projection)
assert_equal(flipped_score, -score)

# Weighted coefficient recovers an exact linear signal.
design <- cbind(intercept = 1, axis_z = seq(-2, 2, length.out = 10))
y <- 3 + 2 * design[, "axis_z"] + c(-0.1, 0.1, 0, 0.05, -0.05, 0, 0.1, -0.1, 0.05, -0.05)
fit <- bulk_weighted_coefficient(y, design, rep(1, 10), "axis_z")
stopifnot(fit$estimable, abs(fit$beta - 2) < 0.05, fit$p_value < 1e-6)

# Nested timing test has two additional df for ns(x, df=3) versus linear, and
# returns a null F statistic when the larger model does not reduce SSE.
nested <- bulk_nested_f_test(
  linear_sse = c(100, 90), nonlinear_sse = c(80, 95),
  df_difference = 2L, nonlinear_residual_df = c(90, 90)
)
stopifnot(nested$f_value[[1L]] > 0, nested$p_value[[1L]] < 0.001)
assert_equal(nested$f_value[[2L]], 0)
assert_equal(nested$p_value[[2L]], 1)

# Shape concordance is computed for every row and is invariant to positive
# scaling within a cohort.
shape_a <- rbind(seq_len(101L), rev(seq_len(101L)))
shape_b <- rbind(2 * seq_len(101L), rev(seq_len(101L)))
assert_equal(bulk_row_spearman(shape_a, shape_b), c(1, 1))

# Membership requires both axes, all four cohort-axis directions, and excludes
# every signature component from ordinary inference.
meta <- data.table(
  gene_id_versioned = rep(c("v1", "v2"), each = 2),
  gene_id_base = rep(c("g1", "g2"), each = 2),
  gene_name = rep(c("A", "B"), each = 2),
  axis_component = rep(c(FALSE, TRUE), each = 2),
  axis_id = rep(c("signature_pc1", "fixed_projection"), 2),
  beta = 1,
  bh_q_value = 0.01,
  direction_concordant = TRUE
)
cohort <- CJ(
  gene_id_versioned = c("v1", "v2"),
  dataset = c("GSE162694", "GSE213621"),
  axis_id = c("signature_pc1", "fixed_projection")
)
cohort[, `:=`(
  gene_id_base = ifelse(gene_id_versioned == "v1", "g1", "g2"),
  beta = 1
)]
membership <- bulk_continuum_membership(
  meta, cohort, c("signature_pc1", "fixed_projection"), 23370L, "g2"
)
stopifnot(membership[gene_id_base == "g1", continuum_associated])
stopifnot(!membership[gene_id_base == "g2", continuum_associated])
stopifnot(membership[gene_id_base == "g2", membership_class] ==
            "axis_constituent_not_independently_testable")

# Frozen real-input contracts: 145 published genes, 139 observed genes, six
# explicit missing genes, and disjoint source-overlap/evaluation participants.
contract <- ml_read_contract()
source_root <- ml_source_root(contract)
signature <- fread(file.path(source_root, "reproduction", "signature_genes.tsv"))
stopifnot(
  nrow(signature) == 145L,
  sum(as.logical(signature$in_resource)) == 139L,
  sum(!as.logical(signature$in_resource)) == 6L
)
manifest <- fread(ml_resolve(contract$sample_manifest))
discovery_ids <- manifest[dataset %in% contract$source_overlap_cohorts, sample_id]
evaluation_ids <- manifest[dataset %in% contract$evaluation_cohorts, sample_id]
stopifnot(length(intersect(discovery_ids, evaluation_ids)) == 0L)

cat("BULK_CONTRACT_TESTS_PASS\n")
