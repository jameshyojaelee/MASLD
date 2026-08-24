#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "lib_molecular_layers.R"))

contract <- ml_read_contract()
candidate <- ml_out_root()
decision_out <- file.path(candidate, "decision")
ml_assert(!dir.exists(decision_out), paste0("Refusing to overwrite decision namespace: ", decision_out))
ml_ensure_dir(decision_out)

required_namespaces <- c("inputs", "bulk", "pathway_tf/pathways", "pathway_tf/tf",
                         "programs/nmf", "programs/hotspot", "paired")
namespace_state <- data.table(
  namespace = required_namespaces,
  path = file.path(candidate, required_namespaces)
)
namespace_state[, exists := dir.exists(path)]
ml_assert(all(namespace_state$exists), "A required parallel workstream namespace is absent")

axis_support <- function(data, id_columns, q_column, beta_column,
                         direction_column = NULL) {
  axes <- as.character(contract$co_primary_axes)
  data[, {
    present <- all(axes %in% axis_id)
    q <- get(q_column)
    beta <- get(beta_column)
    direction <- if (!is.null(direction_column)) get(direction_column) else rep(TRUE, .N)
    list(
      both_axes_present = present,
      both_axes_q_below_005 = present && all(is.finite(q) & q < 0.05),
      both_axes_direction_concordant = present && all(direction) &&
        uniqueN(sign(beta[is.finite(beta) & beta != 0])) == 1L,
      mean_meta_beta = mean(beta, na.rm = TRUE),
      maximum_q_value = if (all(is.finite(q))) max(q) else NA_real_
    )
  }, by = id_columns]
}

message("Integrating transcript membership")
transcript <- fread(file.path(candidate, "bulk", "continuum_membership.tsv.gz"))
ml_assert(nrow(transcript) == contract$gene_family_size, "Transcript membership family drift")
transcript_meta <- fread(file.path(candidate, "bulk", "continuum_gene_meta.tsv.gz"))
transcript_effect <- transcript_meta[, .(
  mean_meta_beta = mean(beta, na.rm = TRUE),
  maximum_q_value = if (all(is.finite(bh_q_value))) max(bh_q_value) else NA_real_
), by = gene_id_base]
transcript_registry <- merge(transcript, transcript_effect, by = "gene_id_base", all.x = TRUE)
transcript_loo_meta <- fread(file.path(candidate, "bulk", "signature_loo_gene_meta.tsv"))
transcript_loo_effect <- transcript_loo_meta[, .(
  loo_mean_meta_beta = mean(beta, na.rm = TRUE),
  loo_maximum_q_value = if (all(is.finite(bh_q_value))) max(bh_q_value) else NA_real_
), by = gene_id_base]
transcript_registry <- merge(
  transcript_registry, transcript_loo_effect, by = "gene_id_base", all.x = TRUE
)
transcript_registry[axis_component == TRUE, `:=`(
  mean_meta_beta = loo_mean_meta_beta,
  maximum_q_value = loo_maximum_q_value
)]
transcript_registry <- transcript_registry[, .(
  molecular_layer = "transcript",
  feature_id = gene_id_base,
  feature_label = gene_name,
  feature_key = paste0("transcript::", gene_id_base),
  testable = !axis_component,
  continuum_associated,
  membership_class,
  mean_meta_beta,
  maximum_q_value,
  family_size = contract$gene_family_size,
  evidence_role = fifelse(
    axis_component, "axis_constituent_leave_one_out_reported_separately",
    "stage_sex_adjusted_continuum_association"
  )
)]

message("Integrating six frozen pathway collections")
pathway_registry_rows <- list()
pathway_checks <- list()
for (collection in names(contract$pathway_collections)) {
  family_size <- as.integer(contract$pathway_collections[[collection]])
  collection_root <- file.path(candidate, "pathway_tf", "pathways", collection)
  meta <- fread(file.path(collection_root, "meta_analysis.tsv"))
  coverage <- fread(file.path(collection_root, "testability.tsv"))
  ml_assert(uniqueN(meta$set_id) == family_size && nrow(meta) == family_size * 2L,
            paste0(collection, " pathway meta family drift"))
  support <- axis_support(meta, "set_id", "fixed_q_value", "fixed_beta",
                          "direction_concordant")
  coverage_support <- coverage[dataset %in% contract$evaluation_cohorts, .(
    testable = all(testable),
    minimum_retained_fraction = min(observed_retained_fraction, na.rm = TRUE),
    minimum_observed_genes = min(n_observed_genes, na.rm = TRUE)
  ), by = set_id]
  support <- merge(support, coverage_support, by = "set_id", all.x = TRUE)
  support[, continuum_associated := testable & both_axes_q_below_005 &
            both_axes_direction_concordant]
  support[, `:=`(
    molecular_layer = paste0("pathway_", collection),
    feature_id = set_id,
    feature_label = set_id,
    feature_key = paste0("pathway_", collection, "::", set_id),
    membership_class = fifelse(
      !testable, "untestable", fifelse(continuum_associated,
        "continuum_associated", "not_continuum_associated")
    ),
    family_size = family_size,
    evidence_role = if (collection == "hallmark") {
      "confirmatory_pathway_family"
    } else {
      "exploratory_complete_pathway_family"
    }
  )]
  pathway_registry_rows[[collection]] <- support[, .(
    molecular_layer, feature_id, feature_label, feature_key, testable,
    continuum_associated, membership_class, mean_meta_beta, maximum_q_value,
    family_size, evidence_role
  )]
  pathway_checks[[collection]] <- data.table(
    check_id = paste0("pathway_family_", collection),
    pass = uniqueN(meta$set_id) == family_size && nrow(meta) == family_size * 2L,
    observed = uniqueN(meta$set_id), expected = family_size
  )
}
pathway_registry <- rbindlist(pathway_registry_rows)

message("Integrating TF, NMF, and Hotspot layers")
tf_meta <- fread(file.path(candidate, "pathway_tf", "tf", "regulon_meta_analysis.tsv"))
ml_assert(uniqueN(tf_meta$tf) == contract$tf_family_size &&
            nrow(tf_meta) == contract$tf_family_size * 2L, "TF family drift")
tf_support <- axis_support(tf_meta, "tf", "fixed_q_value", "fixed_beta",
                           "direction_concordant")
tf_coverage <- fread(file.path(candidate, "pathway_tf", "tf", "regulon_testability.tsv"))
tf_coverage <- tf_coverage[dataset %in% contract$evaluation_cohorts,
                           .(testable = all(testable)), by = tf]
tf_support <- merge(tf_support, tf_coverage, by = "tf", all.x = TRUE)
tf_support[, continuum_associated := testable & both_axes_q_below_005 &
             both_axes_direction_concordant]
tf_registry <- tf_support[, .(
  molecular_layer = "tf_regulon",
  feature_id = tf,
  feature_label = tf,
  feature_key = paste0("tf_regulon::", tf),
  testable,
  continuum_associated,
  membership_class = fifelse(!testable, "untestable", fifelse(
    continuum_associated, "continuum_associated", "not_continuum_associated"
  )),
  mean_meta_beta,
  maximum_q_value,
  family_size = contract$tf_family_size,
  evidence_role = "signed_regulon_score_not_direct_tf_activity"
)]

nmf_membership <- fread(file.path(candidate, "programs", "nmf",
                                  "nmf_continuum_membership.tsv"))
nmf_meta <- fread(file.path(candidate, "programs", "nmf", "nmf_meta_analysis.tsv"))
ml_assert(nrow(nmf_membership) == contract$nmf_family_size, "NMF family drift")
nmf_effect <- nmf_meta[
  score_variant == "signature_excluded_fixed_w" & inference_eligible == TRUE,
  .(
    mean_meta_beta = mean(fixed_beta, na.rm = TRUE),
    maximum_q_value = if (all(is.finite(fixed_q_value))) max(fixed_q_value) else NA_real_
  ),
  by = outcome_id
]
nmf_registry <- merge(nmf_membership, nmf_effect, by = "outcome_id", all.x = TRUE)
nmf_registry <- nmf_registry[, .(
  molecular_layer = "nmf_fixed_basis",
  feature_id = outcome_id,
  feature_label = outcome_id,
  feature_key = paste0("nmf_fixed_basis::", outcome_id),
  testable = fixed_w_gate_pass,
  continuum_associated,
  membership_class = evidence_label,
  mean_meta_beta,
  maximum_q_value,
  family_size = contract$nmf_family_size,
  evidence_role = "fixed_w_signature_excluded_projection"
)]

hotspot_membership <- fread(file.path(candidate, "programs", "hotspot",
                                      "hotspot_continuum_membership.tsv"))
hotspot_meta <- fread(file.path(candidate, "programs", "hotspot",
                                "hotspot_meta_analysis.tsv"))
hotspot_testability <- fread(file.path(candidate, "programs", "hotspot",
                                       "hotspot_testability.tsv"))
ml_assert(nrow(hotspot_membership) == contract$program_family_size, "Hotspot family drift")
hotspot_effect <- hotspot_meta[, .(
  mean_meta_beta = mean(fixed_beta, na.rm = TRUE),
  maximum_q_value = if (all(is.finite(fixed_q_value))) max(fixed_q_value) else NA_real_
), by = program_uid]
hotspot_test <- hotspot_testability[, .(testable = all(testable)), by = program_uid]
hotspot_registry <- Reduce(
  function(x, y) merge(x, y, by = "program_uid", all.x = TRUE),
  list(hotspot_membership, hotspot_effect, hotspot_test)
)
hotspot_registry <- hotspot_registry[, .(
  molecular_layer = "hotspot_program",
  feature_id = program_uid,
  feature_label = module_name,
  feature_key = paste0("hotspot_program::", program_uid),
  testable,
  continuum_associated,
  membership_class = evidence_label,
  mean_meta_beta,
  maximum_q_value,
  family_size = contract$program_family_size,
  evidence_role = "frozen_signature_excluded_weighted_program"
)]

registry <- rbindlist(list(
  transcript_registry, pathway_registry, tf_registry, nmf_registry, hotspot_registry
), use.names = TRUE, fill = TRUE)

message("Adding paired evidence as a non-vetoing upgrade")
paired_gene <- fread(file.path(candidate, "paired", "paired_gene_results.tsv"))
paired_hotspot <- fread(file.path(candidate, "paired", "paired_program_results.tsv"))
paired_nmf <- fread(file.path(candidate, "paired", "paired_nmf_results.tsv"))
paired_pathway <- fread(file.path(candidate, "paired", "paired_pathway_results.tsv.gz"))
paired_gene_key_map <- transcript_registry[
  feature_label %in% unique(paired_gene$feature_id),
  .(feature_id = feature_label, feature_key)
]
ml_assert(nrow(paired_gene_key_map) == uniqueN(paired_gene$feature_id) &&
            !anyDuplicated(paired_gene_key_map$feature_id),
          "Fixed paired manuscript genes do not map one-to-one to the transcript registry")
paired_gene <- merge(
  paired_gene, paired_gene_key_map,
  by = "feature_id", all.x = TRUE, sort = FALSE
)
ml_assert(!anyNA(paired_gene$feature_key),
          "A fixed paired manuscript gene does not map to the GENCODE transcript registry")
paired_hotspot[, feature_key := paste0("hotspot_program::", feature_id)]
paired_nmf[, feature_key := paste0("nmf_native::", feature_id)]
paired_pathway[, feature_key := paste0("pathway_", collection, "::", feature_id)]
paired_all <- rbindlist(list(
  paired_gene[, .(feature_key, axis_id, paired_beta = beta,
                  paired_adjusted_p = adjusted_p_value)],
  paired_hotspot[, .(feature_key, axis_id, paired_beta = beta,
                     paired_adjusted_p = adjusted_p_value)],
  paired_nmf[, .(feature_key, axis_id, paired_beta = beta,
                 paired_adjusted_p = adjusted_p_value)],
  paired_pathway[, .(feature_key, axis_id, paired_beta = beta,
                     paired_adjusted_p = adjusted_p_value)]
), use.names = TRUE, fill = TRUE)
paired_upgrade <- paired_all[, .(
  paired_axes_present = all(contract$co_primary_axes %in% axis_id),
  paired_both_axes_adjusted_p_below_005 = all(
    is.finite(paired_adjusted_p) & paired_adjusted_p < 0.05
  ),
  paired_axis_direction_concordant = uniqueN(
    sign(paired_beta[is.finite(paired_beta) & paired_beta != 0])
  ) == 1L,
  mean_paired_beta = mean(paired_beta, na.rm = TRUE)
), by = feature_key]
registry <- merge(registry, paired_upgrade, by = "feature_key", all.x = TRUE)
registry[, paired_direction_matches_cross_sectional :=
           is.finite(mean_paired_beta) & is.finite(mean_meta_beta) &
           sign(mean_paired_beta) == sign(mean_meta_beta)]
registry[, within_person_supported :=
           paired_axes_present %in% TRUE & paired_both_axes_adjusted_p_below_005 %in% TRUE &
           paired_axis_direction_concordant %in% TRUE &
           paired_direction_matches_cross_sectional %in% TRUE]
registry[, paired_evidence_role := fifelse(
  within_person_supported, "within_person_supported",
  fifelse(is.na(paired_axes_present), "not_prespecified_for_paired_testing",
          "paired_tested_not_supported")
)]
ml_write_tsv_once(registry, file.path(decision_out, "continuum_membership_registry.tsv.gz"))

focal_gate <- fread(file.path(candidate, "programs", "hotspot",
                              "hotspot_focal_figure4f_gate.tsv"))
figure4f_gate <- nrow(focal_gate) == length(contract$focal_programs) &&
  all(focal_gate$figure4f_cross_sectional_gate_pass)
paired_manifest <- fread(file.path(candidate, "paired", "paired_participant_manifest.tsv"))

checks <- rbindlist(list(
  data.table(
    check_id = c(
      "required_namespaces", "transcript_family_23370", "tf_family_268",
      "nmf_family_10", "hotspot_family_117", "paired_54_participants",
      "figure4f_two_focal_programs", "composition_sensitivity"
    ),
    pass = c(
      all(namespace_state$exists), nrow(transcript) == contract$gene_family_size,
      uniqueN(tf_meta$tf) == contract$tf_family_size,
      nrow(nmf_membership) == contract$nmf_family_size,
      nrow(hotspot_membership) == contract$program_family_size,
      uniqueN(paired_manifest$donor_id) == 54L,
      figure4f_gate, NA
    ),
    blocking = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
    state = c(
      rep("tested", 7L), "blocked_upstream_composition"
    )
  ),
  rbindlist(pathway_checks)[, .(
    check_id, pass, blocking = TRUE, state = "tested"
  )]
), use.names = TRUE, fill = TRUE)
checks[, candidate_valid := fifelse(is.na(pass), !blocking, pass)]
ml_assert(all(checks$candidate_valid), "A blocking molecular-layer integration gate failed")
ml_write_tsv_once(checks, file.path(decision_out, "integration_gates.tsv"))
ml_write_tsv_once(focal_gate, file.path(decision_out, "figure4f_focal_gate.tsv"))

layer_summary <- registry[, .(
  family_size = unique(family_size),
  n_features = .N,
  n_testable = sum(testable %in% TRUE),
  n_continuum_associated = sum(continuum_associated %in% TRUE),
  n_within_person_supported = sum(within_person_supported %in% TRUE)
), by = molecular_layer]
ml_write_tsv_once(layer_summary, file.path(decision_out, "layer_summary.tsv"))

decision <- list(
  workstream_id = contract$workstream_id,
  release_state = "candidate_only",
  pipeline_valid = TRUE,
  figure4f_cross_sectional_gate_pass = figure4f_gate,
  current_figure3_modified = FALSE,
  continuum_membership_is_separate_from_canonical_degs = TRUE,
  paired_evidence_is_non_vetoing_upgrade = TRUE,
  composition_state = "blocked_upstream_composition_nonblocking",
  automatic_promotion = FALSE,
  explicit_approval_required = TRUE,
  claim = contract$claim
)
ml_write_json_once(decision, file.path(decision_out, "candidate_decision.json"))
ml_write_session_info(file.path(decision_out, "sessionInfo.txt"))
message("MOLECULAR_LAYER_INTEGRATION_COMPLETE: ", decision_out)
