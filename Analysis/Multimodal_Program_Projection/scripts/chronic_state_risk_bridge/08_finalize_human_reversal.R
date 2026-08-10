#!/usr/bin/env Rscript

script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])
source(file.path(dirname(normalizePath(script_path)), "state_axis_common.R"))
require_plan43_seal()

pre <- read_tsv(file.path(candidate_root, "human_reversal_effects_pre_exact_nas.tsv"))
changes <- read_tsv(file.path(candidate_root, "human_reversal/participant_state_changes.tsv"))
pairs_active <- jsonlite::read_json(file.path(candidate_root, "source_gates/HUMAN_GATE_ACTIVE.json"))
pairs <- read_tsv(file.path(candidate_root, pairs_active$active_pair_manifest))

primary <- changes[changes$dataset_id == "GSE48452" & changes$score_scheme == "weighted__primary" & changes$primary_eligible == "true", ]
input <- data.frame(delta_state_score = primary$delta_state_score, delta_nas = as.numeric(primary$delta_nas))
input_path <- file.path(candidate_root, "human_reversal/gse48452_exact_nas_input.tsv")
write_tsv(input, input_path)
binary <- Sys.getenv("PLAN43_EXACT_NAS_BINARY")
if (!nzchar(binary) || !file.exists(binary)) stop("PLAN43_EXACT_NAS_BINARY is missing")
output_path <- file.path(candidate_root, "human_reversal/gse48452_exact_nas_permutation.tsv")
status <- system2(binary, c(input_path, output_path))
if (status != 0L) stop("Exact NAS permutation failed with status ", status)
exact <- read_tsv(output_path)
target <- pre$dataset_id == "GSE48452" & pre$score_scheme == "weighted__primary"
if (sum(target) != 1L || nrow(exact) != 1L || exact$n_permutations != 10810800) stop("Exact NAS output drift")
if (abs(pre$estimate[target] - exact$estimate) > 1e-12) stop("R/C++ Spearman mismatch")
pre$p[target] <- exact$p

# Re-derive expected directions here as an independent guard against a stale
# pre-exact file. Responder-minus-nonresponder delta-S should be negative;
# delta-S versus delta-NAS should be positive.
pre$expected_direction <- ifelse(pre$dataset_id %in% c("GSE106737", "GSE83452"), "negative", "positive")
pre$direction_agrees <- ifelse(pre$expected_direction == "negative", pre$estimate < 0, pre$estimate > 0)

primary_rows <- pre[pre$score_scheme == "weighted__primary", ]
primary_rows$q <- p.adjust(primary_rows$p, "BH")
pre$q <- NA_real_
pre$q[match(paste(primary_rows$dataset_id, primary_rows$contrast_id), paste(pre$dataset_id, pre$contrast_id))] <- primary_rows$q
write_tsv(pre, file.path(candidate_root, "human_reversal_effects.tsv"))

# GSE106737 and GSE83452 are fingerprint-confirmed overlapping accessions. The
# outcome-unseen GSE83452 analysis represents their one shared cohort in meta-analysis.
meta_rows <- primary_rows[primary_rows$dataset_id %in% c("GSE83452", "GSE48452"), ]
orientation <- ifelse(meta_rows$dataset_id == "GSE83452", -1, 1)
meta_rows$aligned_estimate <- orientation * meta_rows$estimate
z <- qnorm(meta_rows$p / 2, lower.tail = FALSE) * sign(meta_rows$aligned_estimate)
meta_z <- sum(z) / sqrt(length(z)); meta_p <- 2 * pnorm(-abs(meta_z))
meta <- data.frame(meta_analysis_id = "effective_independent_cohorts",
  included_datasets = paste(meta_rows$dataset_id, collapse = ";"), n_effective_cohorts = nrow(meta_rows),
  signed_stouffer_z = meta_z, p = meta_p, expected_direction = "aligned_reversal",
  direction_agrees = all(meta_rows$aligned_estimate > 0), three_cohort_preregistered_gate_evaluable = FALSE,
  reason = "GSE106737 and GSE83452 share 78 fingerprint-confirmed arrays and are one cohort contribution",
  stringsAsFactors = FALSE)
write_tsv(meta, file.path(candidate_root, "human_reversal/human_reversal_meta_analysis.tsv"))

# Correct the same sign convention in the already-computed leave-participant
# sensitivities without changing their estimates.
sensitivity_path <- file.path(candidate_root, "human_reversal/human_reversal_sensitivity.tsv")
sensitivity <- read_tsv(sensitivity_path)
is_leave <- sensitivity$sensitivity == "leave_one_participant_out"
expected_negative <- sensitivity$dataset_id %in% c("GSE106737", "GSE83452")
sensitivity$direction_agrees[is_leave] <- ifelse(expected_negative[is_leave],
  sensitivity$estimate[is_leave] < 0, sensitivity$estimate[is_leave] > 0)
write_tsv(sensitivity, sensitivity_path)
cat("HUMAN_REVERSAL_FINALIZED\n")
