#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: check_locked_substrate.R RUN_ROOT ARM", call. = FALSE)
run_root <- normalizePath(args[[1]], mustWork = TRUE)
arm <- args[[2]]
if (!arm %in% c("R0", "F_locked")) stop("ARM must be R0 or F_locked", call. = FALSE)
if (!file.exists(file.path(run_root, ".bg001_candidate_root"))) stop("Missing candidate sentinel", call. = FALSE)
source(file.path(
  run_root, "source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation/scientific_gate_helpers.R"
))

arm_root <- file.path(run_root, "arms", arm)
locked_qc <- fread(file.path(run_root, "frozen_sets/locked_sample_qc_report.csv"))
candidate_qc <- fread(file.path(arm_root, "qc/sample_qc_report.csv"))
locked_meta <- as.data.table(readRDS(file.path(run_root, "frozen_sets/locked_meta_matched.rds")))
candidate_meta <- as.data.table(readRDS(file.path(arm_root, "results/integration/meta_matched.rds")))
result <- data.table(
  arm = arm,
  full_substrate_samples = nrow(candidate_qc),
  locked_samples = nrow(locked_qc),
  qc_equal_to_frozen = bg001_table_equal_by_key(
    locked_qc, candidate_qc, "sample_id",
    ignored = c("canonical_total_counts", "candidate_total_counts")
  ),
  meta_equal_to_frozen = bg001_table_equal_by_key(locked_meta, candidate_meta, "sample_id")
)
result[, passed := qc_equal_to_frozen & meta_equal_to_frozen]
output <- file.path(arm_root, "provenance/locked_substrate_check.tsv")
if (file.exists(output)) stop("Refusing existing locked-substrate check", call. = FALSE)
fwrite(result, output, sep = "\t")
if (!result$passed) stop(arm, " does not equal the frozen complete QC/meta substrate", call. = FALSE)
cat("PASS locked full-substrate QC/meta equality: ", arm, "\n", sep = "")
