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
out <- ml_out_root(must_exist = TRUE)
input_out <- file.path(out, "inputs")
ml_assert(!dir.exists(input_out), paste0("Refusing to overwrite input freeze: ", input_out))
ml_ensure_dir(input_out)

input_keys <- c(
  "dge_rds", "canonical_deg_results", "sample_manifest", "stage_results",
  "stage_extension_results", "gene_annotation", "program_registry",
  "program_membership", "nmf_loadings", "nmf_cache", "paired_crosswalk",
  "paired_counts", "paired_metadata", "current_fig4e"
)
paths <- setNames(lapply(input_keys, function(key) ml_resolve(contract[[key]])), input_keys)
source_root <- ml_source_root(contract)
source_files <- c(
  source_axes_unsupervised = file.path(source_root, "unsupervised", "participant_scores.tsv"),
  source_axes_projection = file.path(source_root, "projection", "participant_scores.tsv"),
  source_cohort_loadings = file.path(source_root, "unsupervised", "loadings.tsv"),
  source_projection_loadings = file.path(source_root, "projection", "fixed_projection_loadings.tsv"),
  source_normalized_GSE162694 = file.path(
    source_root, "unsupervised", "normalized_expression_by_cohort", "GSE162694.rds"
  ),
  source_normalized_GSE213621 = file.path(
    source_root, "unsupervised", "normalized_expression_by_cohort", "GSE213621.rds"
  ),
  source_discovery_signature_model = file.path(
    source_root, "reproduction", "discovery_signature_model.rds"
  ),
  source_signature = file.path(source_root, "reproduction", "signature_genes.tsv"),
  source_signature_mapping = file.path(source_root, "programs", "signature_gencode_v49_mapping.tsv"),
  source_program_scores = file.path(source_root, "programs", "program_scores.tsv.gz"),
  source_program_testability = file.path(source_root, "programs", "program_testability.tsv"),
  source_program_models = file.path(source_root, "programs", "program_fibrosis_models.tsv"),
  source_program_meta = file.path(source_root, "programs", "program_meta_analysis.tsv"),
  source_paired_scores = file.path(source_root, "paired", "participant_scores.tsv")
)
ml_assert(all(file.exists(source_files)), "The frozen HAC source candidate is incomplete")

pathway_dir <- ml_resolve(contract$pathway_dir)
gmt_files <- file.path(pathway_dir, paste0(names(contract$pathway_collections), ".gmt"))
names(gmt_files) <- paste0("gmt_", names(contract$pathway_collections))
ml_assert(all(file.exists(gmt_files)), "A frozen pathway GMT is missing")

all_paths <- c(unlist(paths), source_files, gmt_files)
manifest <- data.table(
  input_id = names(all_paths),
  path = normalizePath(unname(all_paths), mustWork = TRUE)
)
manifest[, `:=`(
  size_bytes = file.info(path)$size,
  sha256 = vapply(path, ml_sha256, character(1)),
  frozen_before_outcome_analysis = TRUE
)]
ml_write_tsv_once(manifest, file.path(input_out, "input_manifest.tsv"))

code_files <- list.files(
  script_dir, recursive = TRUE, full.names = TRUE,
  pattern = "(\\.R|\\.py|\\.sbatch|\\.sh|\\.json|README\\.md)$"
)
code_files <- sort(normalizePath(code_files, mustWork = TRUE))
code_manifest <- data.table(
  relative_path = sub(paste0("^", script_dir, "/?"), "", code_files),
  path = code_files,
  size_bytes = file.info(code_files)$size,
  sha256 = vapply(code_files, ml_sha256, character(1)),
  frozen_before_outcome_analysis = TRUE
)
ml_write_tsv_once(code_manifest, file.path(input_out, "analysis_code_manifest.tsv"))

contract_copy <- contract
contract_copy$resolved_source_hac_candidate <- source_root
contract_copy$input_manifest_sha256 <- ml_sha256(file.path(input_out, "input_manifest.tsv"))
contract_copy$analysis_code_manifest_sha256 <- ml_sha256(
  file.path(input_out, "analysis_code_manifest.tsv")
)
ml_write_json_once(contract_copy, file.path(input_out, "frozen_contract.json"))

required_packages <- c(
  "data.table", "jsonlite", "edgeR", "limma", "splines", "ggplot2",
  "patchwork", "png"
)
optional_packages <- c("fgsea", "decoupleR", "dorothea", "nnls", "NMF")
package_state <- data.table(
  package = c(required_packages, optional_packages),
  required = c(rep(TRUE, length(required_packages)), rep(FALSE, length(optional_packages)))
)
package_state[, available := vapply(package, requireNamespace, logical(1), quietly = TRUE)]
package_state[, version := vapply(package, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}, character(1))]
ml_write_tsv_once(package_state, file.path(input_out, "package_state.tsv"))
ml_assert(all(package_state[required == TRUE, available]), "A required R package is unavailable")

mapping <- fread(source_files[["source_signature_mapping"]])
signature <- fread(source_files[["source_signature"]])
ml_assert(nrow(signature) == contract$signature_published_size, "Published signature size drift")
observed_column <- if ("in_resource" %in% names(signature)) "in_resource" else NULL
if (!is.null(observed_column)) {
  ml_assert(sum(as.logical(signature[[observed_column]])) == contract$signature_observed_size,
            "Observed signature size drift")
}
ml_assert(nrow(mapping) == contract$signature_published_size, "Signature mapping size drift")
ml_write_tsv_once(mapping, file.path(input_out, "signature_mapping.tsv"))

axes <- ml_load_axes(contract)
axis_census <- axes[, .(
  n_participants = uniqueN(sample_id),
  n_finite = sum(is.finite(axis_raw))
), by = .(dataset, axis_id)]
ml_write_tsv_once(axis_census, file.path(input_out, "axis_census.tsv"))

roles <- data.table(
  dataset = c(
    contract$source_overlap_cohorts,
    contract$evaluation_cohorts,
    contract$paired_cohort
  ),
  continuum_role = c(
    rep("source_overlap_discovery", length(contract$source_overlap_cohorts)),
    rep("competitor_independent_evaluation", length(contract$evaluation_cohorts)),
    "within_person_replication"
  ),
  overlaps_canonical_resource_discovery = c(
    rep(TRUE, length(contract$source_overlap_cohorts) + length(contract$evaluation_cohorts)),
    FALSE
  ),
  eligible_as_independent_validation_of_canonical_degs = FALSE,
  eligible_as_independent_validation_of_source_overlap_signature = c(
    rep(FALSE, length(contract$source_overlap_cohorts)),
    rep(TRUE, length(contract$evaluation_cohorts)),
    FALSE
  ),
  eligible_as_within_person_replication = c(
    rep(FALSE, length(contract$source_overlap_cohorts) + length(contract$evaluation_cohorts)),
    TRUE
  )
)
ml_write_tsv_once(roles, file.path(input_out, "source_overlap_matrix.tsv"))

composition_state <- data.table(
  analysis = "composition_adjustment",
  state = "blocked_upstream_composition",
  promotion_blocking = FALSE,
  reason = paste(
    "The synchronized corrected composition release is not an input to this candidate;",
    "no legacy composition table may be mixed into the analysis."
  )
)
ml_write_tsv_once(composition_state, file.path(input_out, "composition_state.tsv"))
ml_write_session_info(file.path(input_out, "sessionInfo.txt"))
message("MOLECULAR_LAYER_INPUT_FREEZE_COMPLETE: ", input_out)
