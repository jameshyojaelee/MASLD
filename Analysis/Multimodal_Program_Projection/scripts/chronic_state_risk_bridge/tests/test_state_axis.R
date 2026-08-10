#!/usr/bin/env Rscript

script_path <- normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]), mustWork = TRUE)
source(file.path(dirname(dirname(script_path)), "state_axis_common.R"))

genes <- paste0("G", 1:16)
expression <- rbind(
  matrix(rep(c(0, 2, 4, 6), each = 8), nrow = 8),
  matrix(rep(c(6, 4, 2, 0), each = 8), nrow = 8)
)
rownames(expression) <- genes
colnames(expression) <- c("B1", "B2", "F1", "F2")
membership <- rbind(
  data.frame(program_uid = "P1", mapped_symbol = genes[1:8], original_l1_weight = rep(1 / 8, 8)),
  data.frame(program_uid = "P2", mapped_symbol = genes[9:16], original_l1_weight = rep(1 / 8, 8))
)
axis <- data.frame(
  program_uid = c("P1", "P2"), cell_type = c("A", "B"),
  primary_loading = c(0.75, -0.25), sign_loading = c(1, -1),
  lineage_balanced_loading = c(0.5, -0.5), stringsAsFactors = FALSE
)

programs <- program_scores(expression, c("B1", "B2"), membership)
stopifnot(all(programs$audit$testable))
state <- state_score(programs$scores, axis, "primary", membership,
                     minimum_programs = 2L, minimum_axis_mass = 1)
expected_z <- c(-sqrt(0.5), sqrt(0.5), 3 * sqrt(0.5), 5 * sqrt(0.5))
expected <- 0.75 * expected_z - 0.25 * (-expected_z)
stopifnot(max(abs(state$scores - expected)) < 1e-12)

permuted <- expression[, c("F1", "F2", "B1", "B2")]
permuted_programs <- program_scores(permuted, c("B1", "B2"), membership)
stopifnot(max(abs(permuted_programs$scores[, colnames(expression)] - programs$scores)) < 1e-12)
cat("STATE_AXIS_FIXTURE_PASS\n")
