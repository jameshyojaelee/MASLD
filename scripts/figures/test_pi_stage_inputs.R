#!/usr/bin/env Rscript
# Synthetic input checks only; no disease outcomes are fitted.
suppressPackageStartupMessages({library(data.table); library(edgeR)})
root <- normalizePath(Sys.getenv("MASLD_PROJECT_ROOT", "."))
script <- file.path(root, "scripts/figures/build_pi_stage_extensions.R")
fixture <- tempfile("stage-qualified-fixture-")
dir.create(fixture)
stage <- file.path(fixture, "stage")
dir.create(stage)
ids <- paste0("sample", 1:10)
metadata <- data.table(sample_id = ids, dataset = rep(paste0("cohort", 1:5), each = 2),
  analysis_unit_id = paste0("source_donor", 1:10), unit_evidence = "fixture_source_key",
  inferred_sex = rep(c("F", "M"), 5))
dge <- DGEList(matrix(10L, nrow = 2, ncol = 10,
                     dimnames = list(c("ENSG1.1", "ENSG2.1"), ids)))
saveRDS(dge, file.path(stage, "dge.rds"))
saveRDS(metadata, file.path(stage, "metadata.rds"))
annotation <- data.table(gene_id = rownames(dge), ensembl_base = c("ENSG1", "ENSG2"),
  gene_name = c("fixture1", "fixture2"), gene_biotype = "protein_coding", chromosome = "chr1")
fwrite(annotation, file.path(stage, "annotation.tsv"), sep = "\t")
fwrite(data.table(n_samples = 10L, n_genes = 2L, n_cohorts = 5L, identity_qualified = TRUE,
  dge_path = file.path(stage, "dge.rds"), meta_path = file.path(stage, "metadata.rds"),
  gene_metadata_path = file.path(stage, "annotation.tsv")),
  file.path(stage, "model_input_manifest.tsv"), sep = "\t")
writeLines('{"status":"VALIDATED"}', file.path(stage, "VALIDATED.json"))
writeLines('{"n_samples":10,"n_genes":2}', file.path(stage, "expected_counts.json"))
roster <- file.path(fixture, "roster.tsv")
scores <- file.path(fixture, "scores.tsv")
fwrite(metadata, roster, sep = "\t")
fwrite(data.table(sample_id = ids, P1 = 1:10, P2 = 1:10, P3 = 1:10,
                 P4 = 1:10, P5 = 1:10, P6 = 1:10), scores, sep = "\t")
check_preflight <- function(label, success) {
  log <- file.path(fixture, paste0(label, ".log"))
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
    stdout = log, stderr = log, env = c(
      paste0("MASLD_PROJECT_ROOT=", root), paste0("STAGE_RELEASE_ROOT=", stage),
      paste0("FIGURE_CANDIDATE_ROOT=", file.path(fixture, label)),
      paste0("STAGE_CONTEXT_ROSTER=", roster), paste0("STAGE_FIXED_NMF_SCORES=", scores),
      "STAGE_EXTENSION_PREFLIGHT_ONLY=1"))
  stopifnot(if (success) status == 0L else status != 0L)
}
check_preflight("dynamic_ten_source_units", TRUE)
changed <- copy(metadata); changed[1, analysis_unit_id := "wrong_donor"]
fwrite(changed, roster, sep = "\t")
check_preflight("reject_roster_identity_disagreement", FALSE)
fwrite(metadata, roster, sep = "\t")
fwrite(data.table(sample_id = ids[-1], P1 = 1:9), scores, sep = "\t")
check_preflight("reject_incomplete_projection_roster", FALSE)
# A whole axis that failed the fixed reconstruction criteria is NA and is accepted;
# a partly missing axis, an infinite score and an all-NA table are rejected.
fwrite(data.table(sample_id = ids, P1 = 1:10, P2 = 1:10, P3 = 1:10, P4 = 1:10, P5 = 1:10, P6 = NA_real_), scores, sep = "\t")
check_preflight("accept_whole_axis_unavailable", TRUE)
fwrite(data.table(sample_id = ids, P1 = 1:10, P2 = 1:10, P3 = 1:10, P4 = 1:10, P5 = 1:10, P6 = c(NA, 2:10)), scores, sep = "\t")
check_preflight("reject_partly_missing_axis", FALSE)
fwrite(data.table(sample_id = ids, P1 = 1:10, P2 = 1:10, P3 = 1:10, P4 = 1:10, P5 = 1:10, P6 = c(Inf, 2:10)), scores, sep = "\t")
check_preflight("reject_infinite_score", FALSE)
fwrite(data.table(sample_id = ids, P1 = NA_real_, P2 = NA_real_, P3 = NA_real_, P4 = NA_real_, P5 = NA_real_, P6 = NA_real_), scores, sep = "\t")
check_preflight("reject_all_axes_unavailable", FALSE)
fwrite(data.table(sample_id = ids, P1 = 1:10, P2 = 1:10, P3 = 1:10, P4 = 1:10, P5 = 1:10, P6 = 1:10), scores, sep = "\t")

# Load only the contrast function, excluding all script input/output execution.
expressions <- parse(script)
for (expression in expressions) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) && is.symbol(expression[[2L]]) &&
      as.character(expression[[2L]]) %in% c("fit_contrast", "fit_cohort", "program_group_summary")) eval(expression)
}
unit_qualified <- TRUE
dge_all <- dge
cohort_design_audit <- list()
samples <- copy(metadata)[, fib_group := "F0"]
failed <- fit_contrast(samples, "fib_group", "F0", "F1", "F0_to_F1", "fibrosis_adjacent")
stopifnot(nrow(failed$result) == 0L, !failed$audit$estimable,
          failed$audit$bh_family_size == 0L, all(failed$samples$unit_evidence == "fixture_source_key"))
eligibility <- copy(metadata)[, group_binary := "Disease"]
failed_cohort <- fit_cohort("cohort1")
stopifnot(nrow(failed_cohort) == 0L, !cohort_design_audit$cohort1$estimable,
          cohort_design_audit$cohort1$n_control == 0L)
programs <- paste0("P", 1:6)
program_labels <- data.table(program_code = programs, biological_label = programs)
nmf_long <- CJ(sample_id = ids, program = programs)[, `:=`(score = NA_real_, centered_score = NaN, program_label = program)]
failed_program <- program_group_summary(samples, "fib_group", c("F0", "F1"), "fibrosis")
stopifnot(nrow(failed_program) == 12L, all(!failed_program$estimable),
          all(is.na(failed_program$centered_activity)), all(failed_program[group == "F1", n_samples] == 0L),
          all(failed_program$failure_reason == "fixed_nmf_axis_unavailable_failed_reconstruction"))
nmf_long[, score := rep(1:10, 6)]
nonvarying <- program_group_summary(samples, "fib_group", c("F0", "F1"), "fibrosis")
stopifnot(all(nonvarying$failure_reason == "program_score_has_no_finite_standardized_variation"))
cat("PASS: dynamic qualified roster, identity mismatch, projection mismatch, missing stage comparator, missing cohort controls; no models fitted\n")

system_script <- file.path(root, "scripts/analysis/histology_anchored_continuum/molecular_layers/70_infer_molecular_systems.R")
for (expression in parse(system_script)) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) && is.symbol(expression[[2L]]) &&
      as.character(expression[[2L]]) %in% c("system_context_units", "system_score_lineage")) eval(expression)
}
ml_assert <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)
expect_failure <- function(expression) stopifnot(inherits(tryCatch({force(expression); NULL}, error = identity), "error"))
stopifnot(nrow(system_context_units(metadata)) == 10L)
expect_failure(system_context_units(rbind(metadata, metadata[1])))
expect_failure(system_context_units(metadata[, unit_evidence := NULL]))
lineage_root <- file.path(fixture, "scores_lineage")
dir.create(file.path(lineage_root, "provenance"), recursive = TRUE)
fwrite(data.table(path = "fixture_expression", sha256 = "fixture_current_hash"),
       file.path(lineage_root, "provenance/input_manifest.tsv"), sep = "\t")
system_score_lineage(lineage_root, "fixture_current_hash")
expect_failure(system_score_lineage(lineage_root, "fixture_old_hash"))
cat("PASS: system unit uniqueness/evidence and corrected expression lineage; no system outcomes fitted\n")
