#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
source(file.path(
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))),
  "..", "lncrna_biotype_association_common.R"
))

set.seed(20260811)
n <- 1000L
fixture <- data.table(
  gene_id_versioned = sprintf("ENSG%011d.1", seq_len(n)),
  gene_type = rep(c("protein_coding", "lncRNA"), each = n / 2L),
  AveExpr = rnorm(n),
  expression_variability = rexp(n),
  gene_length_bp = sample(500:100000, n, replace = TRUE),
  cohort_detection_count = sample(1:5, n, replace = TRUE)
)
linear <- -3 + 0.6 * (fixture$gene_type == "lncRNA") + 0.3 * fixture$AveExpr
fixture[, treat_positive := rbinom(.N, 1L, plogis(linear)) == 1L]
result <- fit_biotype_association(fixture)
stopifnot(
  nrow(result$result) == 1L,
  result$result$comparison == "lncRNA_vs_protein_coding",
  is.finite(result$result$odds_ratio),
  result$result$odds_ratio > 1,
  nrow(result$model_data) == n
)
cat("PASS: lncRNA biotype association fixture\n")
