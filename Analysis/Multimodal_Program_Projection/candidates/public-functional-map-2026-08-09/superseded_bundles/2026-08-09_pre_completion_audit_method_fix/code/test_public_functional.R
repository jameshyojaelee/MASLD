#!/usr/bin/env Rscript

script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])
source(file.path(dirname(normalizePath(script_path)), "functional_analysis_common.R"))

set.seed(20260809)
genes <- paste0("G", seq_len(10))
expression <- matrix(rnorm(10 * 8), nrow = 10, dimnames = list(genes, paste0("S", seq_len(8))))
membership <- data.frame(
  program_uid = rep("fixture_program", 10), mapped_symbol = genes,
  original_l1_weight = seq_len(10) / sum(seq_len(10)), stringsAsFactors = FALSE
)

weighted <- score_programs(expression, membership, "weighted")
z <- t(scale(t(expression)))
expected_weighted <- as.numeric(crossprod(membership$original_l1_weight, z))
stopifnot(max(abs(weighted$scores[1, ] - expected_weighted)) < 1e-12)

equal <- score_programs(expression, membership, "equal")
expected_equal <- colMeans(z)
stopifnot(max(abs(equal$scores[1, ] - expected_equal)) < 1e-12)

leave_top <- score_programs(expression, membership, "leave_top")
retained <- membership[-which.max(membership$original_l1_weight), ]
expected_leave <- as.numeric(crossprod(retained$original_l1_weight / sum(retained$original_l1_weight), z[retained$mapped_symbol, ]))
stopifnot(max(abs(leave_top$scores[1, ] - expected_leave)) < 1e-12)
stopifnot(weighted$audit$testable, weighted$audit$n_mapped_genes == 10L, abs(weighted$audit$retained_l1_weight - 1) < 1e-12)

metadata <- data.frame(group = rep(c("control", "treated"), each = 4), donor = rep(paste0("D", 1:4), 2), stringsAsFactors = FALSE)
effects <- matrix(rnorm(25 * 8), nrow = 25)
effects[, metadata$group == "treated"] <- effects[, metadata$group == "treated"] + 0.5
rownames(effects) <- paste0("P", seq_len(25))
fit <- fit_group_contrast(effects, metadata, "group", "treated", "control", "donor")
stopifnot(nrow(fit) == 25L, all(is.finite(fit$p)), median(fit$estimate) > 0)

frozen <- load_frozen_objects()
stopifnot(nrow(frozen$registry) == 117L, nrow(frozen$external) == 2L, nrow(frozen$membership) == 7093L)
cat("PUBLIC_FUNCTIONAL_FIXTURES_PASS\tweighted;equal;leave_top;paired_design;frozen_universes\n")
