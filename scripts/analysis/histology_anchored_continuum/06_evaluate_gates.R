#!/usr/bin/env Rscript

# Key decision: evaluate only the frozen promotion gates and leave all promotion
# actions outside this candidate workflow.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "lib_continuum.R"))

pre <- read_prespec()
out <- out_root()

co_primary <- as.character(pre$axes$co_primary)
independent <- as.character(pre$cohorts$independent_primary)
all_resource_cohorts <- c(
  as.character(pre$cohorts$source_overlap),
  independent
)
focal_map <- unlist(pre$programs$focal_programs, use.names = TRUE)
focal_ids <- unname(focal_map)

input_paths <- c(
  prespecification = file.path(script_dir, "00_prespecification.json"),
  prepared_input_manifest = file.path(out, "inputs", "input_manifest.tsv"),
  protected_scope_snapshot = file.path(out, "inputs", "protected_scope_snapshot.tsv"),
  reproduction_input_checksums = file.path(out, "reproduction", "input_checksums.tsv"),
  reproduction_summary = file.path(out, "reproduction", "reproduction_summary.tsv"),
  unsupervised_input_checksums = file.path(out, "unsupervised", "input_checksums.tsv"),
  signature_coverage = file.path(out, "unsupervised", "coverage.tsv"),
  projection_input_checksums = file.path(out, "projection", "input_checksums.tsv"),
  co_primary_agreement = file.path(out, "programs", "co_primary_agreement.tsv"),
  histology_anchors = file.path(out, "programs", "axis_histology_anchors.tsv"),
  program_meta_analysis = file.path(out, "programs", "program_meta_analysis.tsv"),
  program_testability = file.path(out, "programs", "program_testability.tsv"),
  paired_tests = file.path(out, "paired", "paired_tests.tsv")
)
absent <- names(input_paths)[!file.exists(input_paths)]
assert_true(!length(absent), paste0(
  "Gate inputs are absent: ", paste(absent, collapse = ", ")
))
input_paths <- vapply(input_paths, normalizePath, character(1L), mustWork = TRUE)

strict_logical <- function(x, label) {
  if (is.logical(x)) {
    ans <- x
  } else {
    value <- toupper(trimws(as.character(x)))
    ans <- rep(NA, length(value))
    ans[value == "TRUE"] <- TRUE
    ans[value == "FALSE"] <- FALSE
  }
  assert_true(!anyNA(ans), paste0(label, " contains a non-logical value"))
  ans
}

assert_columns <- function(tab, required, label) {
  missing <- setdiff(required, names(tab))
  assert_true(!length(missing), paste0(
    label, " schema drift; missing: ", paste(missing, collapse = ", ")
  ))
}

single_row <- function(tab, label) {
  assert_true(nrow(tab) == 1L, paste0(label, " must contain exactly one row"))
  tab[1L]
}

gate_row <- function(gate_id, gate_family, gate_level, required_for_overall,
                     axis_id = NA_character_, dataset = NA_character_,
                     program_uid = NA_character_, observed_primary = NA_real_,
                     observed_secondary = NA_real_, threshold = NA_real_,
                     comparison, pass, evidence_file, detail) {
  data.table(
    gate_id = gate_id,
    gate_family = gate_family,
    gate_level = gate_level,
    required_for_overall = required_for_overall,
    axis_id = axis_id,
    dataset = dataset,
    program_uid = program_uid,
    observed_primary = as.numeric(observed_primary),
    observed_secondary = as.numeric(observed_secondary),
    threshold = as.numeric(threshold),
    comparison = comparison,
    pass = as.logical(pass),
    evidence_file = evidence_file,
    detail = detail
  )
}

reproduction <- fread(input_paths[["reproduction_summary"]], na.strings = c("", "NA"))
assert_columns(reproduction, c("metric", "value", "pass"), "Reproduction summary")
released <- single_row(
  reproduction[metric == "released_pc1_spearman"],
  "Released-PC1 reproduction metric"
)
released_rho <- suppressWarnings(as.numeric(released$value))
assert_true(is.finite(released_rho), "Released-PC1 Spearman rho is not numeric")
reproduction_threshold <- as.numeric(
  pre$promotion_gates$reproduction_abs_spearman_minimum
)

coverage <- fread(input_paths[["signature_coverage"]], na.strings = c("", "NA"))
assert_columns(
  coverage,
  c("dataset", "n_signature_total", "n_signature_used", "coverage_fraction"),
  "Signature coverage"
)
assert_true(!anyDuplicated(coverage$dataset), "Signature coverage duplicates a cohort")
assert_true(all(coverage$n_signature_total == pre$signature$published_gene_count),
            "Signature coverage denominator drift")
assert_true(setequal(coverage$dataset, all_resource_cohorts), paste0(
  "Signature coverage cohort set drift: ",
  paste(sort(setdiff(all_resource_cohorts, coverage$dataset)), collapse = ", ")
))
coverage <- coverage[match(all_resource_cohorts, dataset)]
coverage_threshold <- as.numeric(pre$signature$minimum_coverage_fraction)

agreement <- fread(input_paths[["co_primary_agreement"]], na.strings = c("", "NA"))
assert_columns(
  agreement,
  c("dataset", "axis_x", "axis_y", "n_donors", "spearman_rho"),
  "Co-primary agreement"
)
agreement <- agreement[
  dataset %in% independent & axis_x == co_primary[[1L]] & axis_y == co_primary[[2L]]
]
assert_true(nrow(agreement) == length(independent) &&
              !anyDuplicated(agreement$dataset) &&
              setequal(agreement$dataset, independent),
            "Independent-cohort co-primary agreement is incomplete")
agreement <- agreement[match(independent, dataset)]
agreement_threshold <- as.numeric(
  pre$promotion_gates$co_primary_agreement_spearman_minimum
)

anchors <- fread(input_paths[["histology_anchors"]], na.strings = c("", "NA"))
assert_columns(
  anchors,
  c("axis_id", "dataset", "anchor", "n_donors", "spearman_rho",
    "bootstrap_ci_low", "bootstrap_ci_high", "bootstrap_replicates"),
  "Histology anchors"
)
anchors <- anchors[
  axis_id %in% co_primary & dataset %in% independent & anchor == "fibrosis_stage"
]
assert_true(nrow(anchors) == length(co_primary) * length(independent) &&
              !anyDuplicated(anchors[, .(axis_id, dataset)]),
            "Co-primary fibrosis anchors are incomplete")
assert_true(all(anchors$bootstrap_replicates == pre$resampling$bootstrap_replicates),
            "Anchor bootstrap replicate count drift")
anchors[, `:=`(
  axis_order__ = match(axis_id, co_primary),
  dataset_order__ = match(dataset, independent)
)]
setorder(anchors, axis_order__, dataset_order__)
anchors[, c("axis_order__", "dataset_order__") := NULL]

program_meta_all <- fread(
  input_paths[["program_meta_analysis"]], na.strings = c("", "NA")
)
assert_columns(
  program_meta_all,
  c("axis_id", "program_uid", "estimable", "beta", "p_value", "q_value",
    "n_cohorts", "direction_concordant", "module_name", "cohort_betas"),
  "Program meta-analysis"
)
co_primary_meta <- program_meta_all[axis_id %in% co_primary]
assert_true(
  nrow(co_primary_meta) == length(co_primary) * pre$programs$family_size &&
    all(co_primary_meta[, uniqueN(program_uid), by = axis_id]$V1 ==
          pre$programs$family_size) &&
    !anyDuplicated(co_primary_meta[, .(axis_id, program_uid)]),
  "Co-primary program meta-analysis does not contain two complete 117-test families"
)
co_primary_meta[, q_recomputed__ := p_adjust_complete_family(
  p_value, pre$programs$family_size
), by = axis_id]
q_match <- (is.na(co_primary_meta$q_value) & is.na(co_primary_meta$q_recomputed__)) |
  (is.finite(co_primary_meta$q_value) & is.finite(co_primary_meta$q_recomputed__) &
     abs(co_primary_meta$q_value - co_primary_meta$q_recomputed__) < 1e-12)
assert_true(all(q_match), "Program q-values do not reproduce BH adjustment with n=117")
co_primary_meta[, q_recomputed__ := NULL]
program_meta <- co_primary_meta[
  axis_id %in% co_primary & program_uid %in% focal_ids
]
assert_true(nrow(program_meta) == length(co_primary) * length(focal_ids) &&
              !anyDuplicated(program_meta[, .(axis_id, program_uid)]),
            "Focal co-primary program meta-analysis is incomplete")
program_meta[, `:=`(
  estimable = strict_logical(estimable, "Program estimable"),
  direction_concordant = strict_logical(
    direction_concordant, "Program direction concordance"
  )
)]

testability <- fread(
  input_paths[["program_testability"]], na.strings = c("", "NA")
)
assert_columns(
  testability,
  c("dataset", "program_uid", "n_observed_genes", "retained_l1_fraction",
    "observed_l1_fraction", "testable", "n_excluded_signature_genes",
    "excluded_l1_weight"),
  "Program testability"
)
testability <- testability[
  dataset %in% independent & program_uid %in% focal_ids
]
assert_true(nrow(testability) == length(independent) * length(focal_ids) &&
              !anyDuplicated(testability[, .(dataset, program_uid)]),
            "Independent-cohort focal testability is incomplete")
testability[, testable := strict_logical(testable, "Program testability")]

program_q_threshold <- as.numeric(pre$promotion_gates$program_meta_bh_maximum)
focal_rows <- lapply(seq_along(focal_map), function(i) {
  program_key <- names(focal_map)[[i]]
  program_id <- unname(focal_map[[i]])
  a <- single_row(
    program_meta[program_uid == program_id & axis_id == co_primary[[1L]]],
    paste0(program_key, "/", co_primary[[1L]])
  )
  b <- single_row(
    program_meta[program_uid == program_id & axis_id == co_primary[[2L]]],
    paste0(program_key, "/", co_primary[[2L]])
  )
  cov <- testability[program_uid == program_id]
  sign_concordant <- is.finite(a$beta) && is.finite(b$beta) &&
    sign(a$beta) != 0 && sign(a$beta) == sign(b$beta)
  a_pass <- a$estimable && a$n_cohorts == length(independent) &&
    a$direction_concordant && is.finite(a$q_value) && a$q_value < program_q_threshold
  b_pass <- b$estimable && b$n_cohorts == length(independent) &&
    b$direction_concordant && is.finite(b$q_value) && b$q_value < program_q_threshold
  all_testable <- nrow(cov) == length(independent) && all(cov$testable)
  program_pass <- a_pass && b_pass && sign_concordant && all_testable
  failure <- c(
    if (!all_testable) "not_testable_in_both_independent_cohorts",
    if (!a_pass) paste0(co_primary[[1L]], "_meta_gate_failed"),
    if (!b_pass) paste0(co_primary[[2L]], "_meta_gate_failed"),
    if (!sign_concordant) "cross_axis_meta_beta_sign_discordant"
  )
  data.table(
    focal_program = program_key,
    program_uid = program_id,
    module_name = unique(na.omit(c(a$module_name, b$module_name)))[1L],
    all_independent_cohorts_testable = all_testable,
    minimum_observed_genes = min(cov$n_observed_genes),
    minimum_retained_l1_fraction = min(cov$retained_l1_fraction),
    minimum_observed_l1_fraction = min(cov$observed_l1_fraction),
    n_excluded_signature_genes = unique(cov$n_excluded_signature_genes)[1L],
    excluded_l1_weight = unique(cov$excluded_l1_weight)[1L],
    signature_pc1_beta = a$beta,
    signature_pc1_q_value = a$q_value,
    signature_pc1_direction_concordant = a$direction_concordant,
    signature_pc1_cohort_betas = a$cohort_betas,
    fixed_projection_beta = b$beta,
    fixed_projection_q_value = b$q_value,
    fixed_projection_direction_concordant = b$direction_concordant,
    fixed_projection_cohort_betas = b$cohort_betas,
    cross_axis_meta_beta_sign_concordant = sign_concordant,
    program_gate_pass = program_pass,
    failure_reason = if (length(failure)) paste(failure, collapse = ";") else NA_character_
  )
})
focal_gate <- rbindlist(focal_rows, fill = TRUE)
assert_true(nrow(focal_gate) == 2L, "Focal program gate must contain two programs")

paired <- fread(input_paths[["paired_tests"]], na.strings = c("", "NA"))
assert_columns(
  paired,
  c("axis_id", "endpoint", "n_donors", "spearman_rho",
    "permutation_p_value", "permutation_replicates", "holm_p_value", "alternative"),
  "Paired tests"
)
paired <- paired[axis_id %in% co_primary & endpoint == "delta_fibrosis"]
assert_true(nrow(paired) == length(co_primary) && !anyDuplicated(paired$axis_id) &&
              setequal(paired$axis_id, co_primary),
            "Paired delta-fibrosis gate rows are incomplete")
assert_true(all(paired$alternative == "greater") &&
              all(paired$permutation_replicates == pre$resampling$permutation_replicates),
            "Paired test direction or permutation count drift")
assert_true(all(paired$n_donors == 54L),
            "Paired promotion rows must use all 54 QC-valid donors")
paired[, holm_recomputed__ := p.adjust(permutation_p_value, method = "holm")]
assert_true(all(is.finite(paired$holm_p_value)) &&
              all(abs(paired$holm_p_value - paired$holm_recomputed__) < 1e-12),
            "Paired Holm p-values do not reproduce the two-score family")
paired[, holm_recomputed__ := NULL]
paired <- paired[match(co_primary, axis_id)]
paired_threshold <- as.numeric(pre$promotion_gates$paired_holm_maximum)

atomic <- list()
atomic[[length(atomic) + 1L]] <- gate_row(
  "released_pc1_abs_spearman", "exact_reproduction", "atomic", FALSE,
  observed_primary = abs(released_rho), threshold = reproduction_threshold,
  comparison = ">=", pass = abs(released_rho) >= reproduction_threshold,
  evidence_file = "reproduction/reproduction_summary.tsv",
  detail = paste0("signed Spearman rho=", signif(released_rho, 8))
)
for (i in seq_len(nrow(coverage))) {
  row <- coverage[i]
  atomic[[length(atomic) + 1L]] <- gate_row(
    paste0("signature_coverage_", row$dataset), "signature_coverage", "atomic", FALSE,
    dataset = row$dataset, observed_primary = row$coverage_fraction,
    observed_secondary = row$n_signature_used, threshold = coverage_threshold,
    comparison = ">=", pass = is.finite(row$coverage_fraction) &&
      row$coverage_fraction >= coverage_threshold,
    evidence_file = "unsupervised/coverage.tsv",
    detail = paste0(row$n_signature_used, "/", row$n_signature_total, " genes")
  )
}
for (i in seq_len(nrow(agreement))) {
  row <- agreement[i]
  atomic[[length(atomic) + 1L]] <- gate_row(
    paste0("co_primary_agreement_", row$dataset), "co_primary_agreement",
    "atomic", FALSE, axis_id = paste(co_primary, collapse = "_vs_"),
    dataset = row$dataset, observed_primary = row$spearman_rho,
    observed_secondary = row$n_donors, threshold = agreement_threshold,
    comparison = ">=", pass = is.finite(row$spearman_rho) &&
      row$spearman_rho >= agreement_threshold,
    evidence_file = "programs/co_primary_agreement.tsv",
    detail = "Spearman agreement between the two co-primary scores"
  )
}
for (i in seq_len(nrow(anchors))) {
  row <- anchors[i]
  atomic[[length(atomic) + 1L]] <- gate_row(
    paste0("fibrosis_anchor_", row$axis_id, "_", row$dataset),
    "fibrosis_anchor", "atomic", FALSE, axis_id = row$axis_id,
    dataset = row$dataset, observed_primary = row$spearman_rho,
    observed_secondary = row$bootstrap_ci_low, threshold = 0,
    comparison = "rho>0_and_bootstrap_ci_low>0",
    pass = is.finite(row$spearman_rho) && row$spearman_rho > 0 &&
      is.finite(row$bootstrap_ci_low) && row$bootstrap_ci_low > 0,
    evidence_file = "programs/axis_histology_anchors.tsv",
    detail = paste0(
      "bootstrap 95% CI [", signif(row$bootstrap_ci_low, 6), ", ",
      signif(row$bootstrap_ci_high, 6), "]; n=", row$n_donors
    )
  )
}
for (i in seq_len(nrow(focal_gate))) {
  row <- focal_gate[i]
  atomic[[length(atomic) + 1L]] <- gate_row(
    paste0("focal_program_", row$focal_program), "focal_program", "alternative",
    FALSE, axis_id = paste(co_primary, collapse = "+"),
    program_uid = row$program_uid,
    observed_primary = max(row$signature_pc1_q_value, row$fixed_projection_q_value),
    observed_secondary = min(abs(row$signature_pc1_beta), abs(row$fixed_projection_beta)),
    threshold = program_q_threshold,
    comparison = "both_q<0.05_and_each_cohort_direction_concordant_and_cross_axis_sign_concordant",
    pass = row$program_gate_pass,
    evidence_file = "programs/program_meta_analysis.tsv",
    detail = ifelse(is.na(row$failure_reason), "all focal-program conditions pass",
                    row$failure_reason)
  )
}
for (i in seq_len(nrow(paired))) {
  row <- paired[i]
  atomic[[length(atomic) + 1L]] <- gate_row(
    paste0("paired_delta_fibrosis_", row$axis_id), "paired_replication", "atomic",
    FALSE, axis_id = row$axis_id, dataset = pre$cohorts$paired_replication,
    observed_primary = row$spearman_rho, observed_secondary = row$holm_p_value,
    threshold = paired_threshold, comparison = "rho>0_and_Holm_p<0.05",
    pass = is.finite(row$spearman_rho) && row$spearman_rho > 0 &&
      is.finite(row$holm_p_value) && row$holm_p_value < paired_threshold,
    evidence_file = "paired/paired_tests.tsv",
    detail = paste0(
      "one-sided permutation p=", signif(row$permutation_p_value, 6),
      "; Holm p=", signif(row$holm_p_value, 6), "; n=", row$n_donors
    )
  )
}
atomic <- rbindlist(atomic, fill = TRUE)

family_spec <- list(
  exact_reproduction = list(expected = 1L, mode = "all"),
  signature_coverage = list(expected = length(all_resource_cohorts), mode = "all"),
  co_primary_agreement = list(expected = length(independent), mode = "all"),
  fibrosis_anchor = list(expected = length(co_primary) * length(independent), mode = "all"),
  focal_program = list(expected = 1L, mode = "any"),
  paired_replication = list(expected = length(co_primary), mode = "all")
)
family_rows <- lapply(names(family_spec), function(family) {
  spec <- family_spec[[family]]
  block <- atomic[gate_family == family]
  n_pass <- sum(block$pass)
  family_pass <- if (spec$mode == "any") n_pass >= spec$expected else
    nrow(block) == spec$expected && n_pass == spec$expected
  gate_row(
    paste0("family_", family), family, "family", TRUE,
    observed_primary = n_pass, observed_secondary = nrow(block),
    threshold = spec$expected,
    comparison = if (spec$mode == "any") ">=" else "all_expected_components_pass",
    pass = family_pass,
    evidence_file = paste(unique(block$evidence_file), collapse = ";"),
    detail = paste0(n_pass, "/", nrow(block), " component or alternative rows pass")
  )
})
family_gates <- rbindlist(family_rows, fill = TRUE)
overall_pass <- all(family_gates$pass)
overall <- gate_row(
  "overall_promotion_gate", "overall", "overall", FALSE,
  observed_primary = sum(family_gates$pass), observed_secondary = nrow(family_gates),
  threshold = nrow(family_gates), comparison = "all_gate_families_pass",
  pass = overall_pass, evidence_file = "decision/gate_decisions.tsv",
  detail = paste0(sum(family_gates$pass), "/", nrow(family_gates),
                  " prespecified gate families pass")
)
gate_decisions <- rbindlist(list(atomic, family_gates, overall), fill = TRUE)

current_placement <- "candidate_only_not_promoted"
next_action <- if (overall_pass) {
  paste(
    "All prespecified gates pass, but promotion remains prohibited until explicit approval",
    "and scripts/manuscript/validate_resource_scope.py succeeds."
  )
} else {
  paste(
    "At least one prespecified gate fails. Retain this analysis in Figure S3 and",
    "leave Figure 3 unchanged."
  )
}
recommended_placement <- if (overall_pass) {
  "Figure_3_only_after_explicit_approval_and_scope_validation"
} else {
  "Figure_S3_only"
}
disposition <- if (overall_pass) {
  "eligible_pending_explicit_approval"
} else {
  "supplement_only_gate_failure"
}
figure_route <- if (overall_pass) "figure3_proposal" else "figure_s3_only"

decision_tsv <- data.table(
  workstream_id = pre$workstream_id,
  candidate_root = out,
  release_state = pre$release_state,
  overall_pass = overall_pass,
  disposition = disposition,
  figure_route = figure_route,
  explicit_approval_required = TRUE,
  promotion_status = current_placement,
  all_prespecified_gates_pass = overall_pass,
  figure_3_eligible_pending_explicit_approval = overall_pass,
  automatically_promoted = FALSE,
  figure_3_modified = FALSE,
  current_placement = current_placement,
  recommended_placement = recommended_placement,
  n_gate_families = nrow(family_gates),
  n_gate_families_pass = sum(family_gates$pass),
  composition_adjustment_state =
    "deferred_pending_synchronized_five_cohort_composition_release",
  required_next_action = next_action,
  scientific_claim = pre$scientific_claim
)
decision_json <- as.list(decision_tsv[1L])
decision_json$passing_focal_programs <- unname(
  focal_gate[program_gate_pass == TRUE, program_uid]
)
decision_json$failed_gate_families <- family_gates[pass == FALSE, gate_family]
decision_json$competitor_release <- pre$competitor_release
decision_json$prior_slingshot_benchmark_role <- paste(
  "preserved alternative sensitivity; it is not the released operational PC1",
  "and cannot support a histology-wins conclusion"
)

input_provenance <- data.table(
  artifact_role = names(input_paths),
  path = unname(input_paths),
  size_bytes = as.numeric(file.size(unname(input_paths))),
  sha256 = vapply(unname(input_paths), sha256_file, character(1L))
)

gate_dir <- file.path(out, "decision")
ensure_new_dir(gate_dir)
write_tsv_once(gate_decisions, file.path(gate_dir, "gate_decisions.tsv"))
write_tsv_once(focal_gate, file.path(gate_dir, "focal_program_gate.tsv"))
write_tsv_once(decision_tsv, file.path(gate_dir, "promotion_decision.tsv"))
write_json_once(decision_json, file.path(gate_dir, "promotion_decision.json"))
write_tsv_once(input_provenance, file.path(gate_dir, "input_checksums.tsv"))
write_session_info(file.path(gate_dir, "sessionInfo.txt"))

output_paths <- c(
  gate_decisions = file.path(gate_dir, "gate_decisions.tsv"),
  focal_program_gate = file.path(gate_dir, "focal_program_gate.tsv"),
  promotion_decision_tsv = file.path(gate_dir, "promotion_decision.tsv"),
  promotion_decision_json = file.path(gate_dir, "promotion_decision.json"),
  gate_session_info = file.path(gate_dir, "sessionInfo.txt")
)
output_provenance <- data.table(
  artifact_role = names(output_paths),
  path = unname(output_paths),
  size_bytes = as.numeric(file.size(unname(output_paths))),
  sha256 = vapply(unname(output_paths), sha256_file, character(1L))
)
write_tsv_once(output_provenance, file.path(gate_dir, "output_checksums.tsv"))

message(
  "GATE_EVALUATION_COMPLETE: all_prespecified_gates_pass=", overall_pass,
  "; promotion_status=", current_placement
)
