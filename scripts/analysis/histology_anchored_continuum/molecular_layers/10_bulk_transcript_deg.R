#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

bulk_args <- commandArgs(trailingOnly = FALSE)
bulk_file <- sub("^--file=", "", bulk_args[grepl("^--file=", bulk_args)])
stopifnot(length(bulk_file) == 1L)
bulk_script_dir <- dirname(normalizePath(bulk_file))
source(file.path(bulk_script_dir, "lib_molecular_layers.R"))
source(file.path(bulk_script_dir, "10_bulk_lib.R"))

contract <- ml_read_contract()
set.seed(contract$seed)
final_out <- file.path(ml_out_root(), "bulk")
ml_assert(!dir.exists(final_out), paste0("Refusing to overwrite bulk namespace: ", final_out))
temporary_out <- file.path(
  ml_out_root(), paste0(".bulk.tmp.", Sys.getenv("SLURM_JOB_ID", Sys.getpid()), ".", Sys.getpid())
)
ml_assert(!dir.exists(temporary_out), paste0("Temporary bulk namespace exists: ", temporary_out))
ml_ensure_dir(temporary_out)

dge_path <- ml_resolve(contract$dge_rds)
manifest_path <- ml_resolve(contract$sample_manifest)
annotation_path <- ml_resolve(contract$gene_annotation)
canonical_path <- ml_resolve(contract$canonical_deg_results)
stage_path <- ml_resolve(contract$stage_results)
source_root <- ml_source_root(contract)
signature_path <- file.path(source_root, "reproduction", "signature_genes.tsv")
axes_unsupervised_path <- file.path(source_root, "unsupervised", "participant_scores.tsv")
axes_projection_path <- file.path(source_root, "projection", "participant_scores.tsv")
loading_path <- file.path(source_root, "unsupervised", "loadings.tsv")
projection_loading_path <- file.path(source_root, "projection", "fixed_projection_loadings.tsv")
normalized_dir <- file.path(source_root, "unsupervised", "normalized_expression_by_cohort")
input_paths <- c(
  dge = dge_path, participant_manifest = manifest_path, gene_annotation = annotation_path,
  canonical_coefficients = canonical_path, stage_coefficients = stage_path,
  signature = signature_path, axes_unsupervised = axes_unsupervised_path,
  axes_projection = axes_projection_path, cohort_pc1_loadings = loading_path,
  fixed_projection_loadings = projection_loading_path,
  normalized_GSE162694 = file.path(normalized_dir, "GSE162694.rds"),
  normalized_GSE213621 = file.path(normalized_dir, "GSE213621.rds")
)
ml_assert(all(file.exists(input_paths)), "A required bulk-layer input is absent")
input_manifest <- data.table(
  input_id = names(input_paths),
  path = normalizePath(unname(input_paths), mustWork = TRUE)
)
input_manifest[, `:=`(
  size_bytes = file.info(path)$size,
  sha256 = vapply(path, ml_sha256, character(1)),
  frozen_before_bulk_outcomes = TRUE
)]
ml_write_tsv_once(input_manifest, file.path(temporary_out, "input_manifest.tsv"))

message("Loading the selected 23,370-gene F_five substrate")
dge <- readRDS(dge_path)
ml_assert(inherits(dge, "DGEList"), "The selected bulk substrate is not a DGEList")
ml_assert(identical(dim(dge), c(contract$gene_family_size, 844L)),
          "F_five gene/participant census drift")
gene_base <- ml_base_gene_id(rownames(dge))
ml_assert(!anyDuplicated(gene_base), "Version stripping creates duplicate gene IDs")

manifest <- fread(manifest_path, na.strings = c("", "NA"))
bulk_required_columns(
  manifest,
  c("sample_id", "analysis_unit_id", "dataset", "group_binary", "inferred_sex",
    "fibrosis_stage", "nas_score"),
  "five-cohort participant manifest"
)
ml_assert(nrow(manifest) == 844L && !anyDuplicated(manifest$sample_id) &&
            !anyDuplicated(manifest$analysis_unit_id),
          "Five-cohort manifest is not one row per participant")
ml_assert(setequal(manifest$sample_id, colnames(dge)),
          "F_five DGE and participant manifest sample sets differ")
manifest <- manifest[match(colnames(dge), sample_id)]
ml_assert(identical(manifest$sample_id, colnames(dge)), "Manifest ordering failed")

annotation <- fread(annotation_path)
bulk_required_columns(annotation, c("gene_id", "gene_name"), "GENCODE v49 annotation")
annotation[, gene_id_base := ml_base_gene_id(gene_id)]
annotation <- unique(annotation[, .(gene_id_base, gene_name)])
ml_assert(annotation[, uniqueN(gene_name), by = gene_id_base][V1 > 1L, .N] == 0L,
          "One Ensembl base ID maps to multiple GENCODE symbols")
gene_annotation <- annotation[match(gene_base, gene_id_base)]
gene_annotation[, gene_id_base := gene_base]
gene_annotation[is.na(gene_name) | gene_name == "", gene_name := gene_id_base]

signature <- fread(signature_path)
bulk_required_columns(signature, c("gene_id_base", "gene_symbol", "in_resource"),
                      "frozen signature")
signature[, gene_id_base := ml_base_gene_id(gene_id_base)]
ml_assert(nrow(signature) == contract$signature_published_size &&
            sum(as.logical(signature$in_resource)) == contract$signature_observed_size,
          "Frozen 145/139 signature census drift")
signature_ids <- signature[as.logical(in_resource), gene_id_base]
# The frozen axes use the 139 signature genes flagged in_resource; all must be in
# the bulk universe. A refit with a larger expression-filtered universe (the
# corrected-control fit adds DPEP1) can observe a signature gene outside the frozen
# axis; it stays out of the axis and, because the axis does not use it, is tested
# as an ordinary gene like any other.
ml_assert(all(signature_ids %in% gene_base),
          "Frozen signature genes are missing from the F_five gene universe")
signature_outside_axis <- setdiff(intersect(signature$gene_id_base, gene_base), signature_ids)
if (length(signature_outside_axis)) {
  message("Signature genes observed but outside the frozen axis: ",
          paste(signature[gene_id_base %in% signature_outside_axis, gene_symbol], collapse = ", "))
}
signature[, `:=`(
  observed_in_F_five = gene_id_base %in% gene_base,
  ordinary_continuum_inference_eligible = !gene_id_base %in% signature_ids,
  evidence_role = fifelse(
    gene_id_base %in% signature_ids,
    "axis_component_leave_one_out_only",
    fifelse(gene_id_base %in% gene_base,
            "signature_gene_outside_frozen_axis",
            "unavailable_in_F_five")
  )
)]
ml_write_tsv_once(signature, file.path(temporary_out, "signature_gene_ledger.tsv"))

axes <- ml_load_axes(contract)
evaluation <- as.character(contract$evaluation_cohorts)
source_overlap <- as.character(contract$source_overlap_cohorts)
co_primary <- as.character(contract$co_primary_axes)
ml_assert(setequal(unique(axes$dataset), c(source_overlap, evaluation)),
          "Frozen continuum axes do not cover exactly the five Resource cohorts")
ml_assert(!anyDuplicated(axes[, .(sample_id, axis_id)]), "Duplicate continuum score rows")

message("Fitting cohort-specific stage/sex-adjusted continuum transcript models")
continuum_rows <- list()
voom_cache <- list()
timing_rows <- list()
timing_prediction_cache <- list()
for (cohort in evaluation) {
  for (axis_name in co_primary) {
    model_meta <- merge(
      manifest[dataset == cohort],
      axes[dataset == cohort & axis_id == axis_name, .(sample_id, axis_raw)],
      by = "sample_id", all = FALSE, sort = FALSE
    )
    model_meta <- model_meta[
      is.finite(axis_raw) & !is.na(fibrosis_stage) & !is.na(inferred_sex)
    ]
    model_meta[, `:=`(
      axis_z = ml_standardize(axis_raw),
      fibrosis_factor = factor(fibrosis_stage),
      sex_factor = factor(inferred_sex)
    )]
    ml_assert(nrow(model_meta) >= 20L && uniqueN(model_meta$fibrosis_factor) >= 2L &&
                uniqueN(model_meta$sex_factor) >= 2L,
              paste0(cohort, " has insufficient complete stage/sex continuum data"))
    design <- model.matrix(~ fibrosis_factor + sex_factor + axis_z, data = model_meta)
    rownames(design) <- model_meta$sample_id
    fitted <- bulk_voom_fit(
      dge, model_meta$sample_id, design, "axis_z", gene_annotation,
      cohort, "continuum_stage_sex_adjusted", axis_name
    )
    result <- fitted$table
    result[, `:=`(
      axis_component = gene_id_base %in% signature_ids,
      inference_eligible = !gene_id_base %in% signature_ids,
      p_value = fifelse(gene_id_base %in% signature_ids, NA_real_, raw_p_value)
    )]
    result[, bh_q_value := ml_complete_bh(p_value, contract$gene_family_size)]
    result[, `:=`(
      bh_family_size = contract$gene_family_size,
      model = "expression ~ continuum_z + factor(fibrosis_stage) + inferred_sex",
      replication_unit = "participant",
      overlaps_continuum_score_construction = TRUE,
      overlaps_canonical_resource_discovery = TRUE,
      overlaps_source_overlap_discovery = FALSE,
      overlaps_paired_replication = FALSE,
      eligible_as_independent_validation_of_canonical_degs = FALSE,
      evaluation_role = "competitor_independent_continuum_evaluation"
    )]
    continuum_rows[[paste(cohort, axis_name)]] <- result
    voom_cache[[paste(cohort, axis_name)]] <- list(
      voom = fitted$voom,
      fit = fitted$fit,
      metadata = model_meta,
      design = design
    )
    timing <- bulk_spline_timing(
      fitted$voom, fitted$fit, design, model_meta$axis_z,
      gene_annotation, cohort, axis_name, signature_ids,
      contract$gene_family_size
    )
    timing$table[, `:=`(
      replication_unit = "participant",
      overlaps_continuum_score_construction = TRUE,
      overlaps_canonical_resource_discovery = TRUE,
      overlaps_source_overlap_discovery = FALSE,
      overlaps_paired_replication = FALSE,
      eligible_as_independent_validation_of_canonical_degs = FALSE
    )]
    timing_rows[[paste(cohort, axis_name)]] <- timing$table
    timing_prediction_cache[[paste(cohort, axis_name)]] <- timing$predictions
  }
}
continuum_cohort <- rbindlist(continuum_rows, use.names = TRUE)
ml_assert(nrow(continuum_cohort) == contract$gene_family_size * length(evaluation) *
            length(co_primary), "Continuum cohort ledger is incomplete")
continuum_meta <- bulk_meta_table(
  continuum_cohort, contract$gene_family_size, signature_ids
)
continuum_meta[, `:=`(
  evaluation_cohorts = paste(evaluation, collapse = ";"),
  overlaps_continuum_score_construction = TRUE,
  overlaps_canonical_resource_discovery = TRUE,
  overlaps_source_overlap_discovery = FALSE,
  overlaps_paired_replication = FALSE,
  eligible_as_independent_validation_of_canonical_degs = FALSE
)]
membership <- bulk_continuum_membership(
  continuum_meta, continuum_cohort, co_primary,
  contract$gene_family_size, signature_ids
)
ml_write_tsv_once(continuum_cohort, file.path(temporary_out, "continuum_gene_cohort.tsv.gz"))
ml_write_tsv_once(continuum_meta, file.path(temporary_out, "continuum_gene_meta.tsv.gz"))
ml_write_tsv_once(membership, file.path(temporary_out, "continuum_membership.tsv.gz"))

timing_cohort <- rbindlist(timing_rows, use.names = TRUE)
ml_assert(nrow(timing_cohort) == contract$gene_family_size * length(evaluation) *
            length(co_primary), "Nonlinear timing ledger is incomplete")
timing_shape_rows <- lapply(co_primary, function(axis_name) {
  first <- timing_prediction_cache[[paste(evaluation[[1L]], axis_name)]]
  second <- timing_prediction_cache[[paste(evaluation[[2L]], axis_name)]]
  ml_assert(identical(rownames(first), rownames(second)),
            paste0(axis_name, " spline gene orders differ between cohorts"))
  q_values <- dcast(
    timing_cohort[axis_id == axis_name],
    gene_id_versioned + gene_id_base + gene_name + axis_component ~ dataset,
    value.var = "bh_q_value"
  )
  ml_assert(all(evaluation %in% names(q_values)),
            paste0(axis_name, " timing table lacks an evaluation cohort"))
  shape_values <- data.table(
    gene_id_versioned = rownames(first),
    shape_spearman = bulk_row_spearman(first, second)
  )
  q_values <- merge(
    q_values, shape_values, by = "gene_id_versioned", all.x = TRUE, sort = FALSE
  )
  q_values[, `:=`(
    axis_id = axis_name,
    prediction_grid = "101_within_cohort_axis_percentiles",
    shape_threshold = 0.70
  )]
  q_values[, nonlinear_timing_supported := !axis_component &
    get(evaluation[[1L]]) < 0.05 & get(evaluation[[2L]]) < 0.05 &
    shape_spearman >= 0.70]
  q_values[, evidence_role := fifelse(
    axis_component,
    "axis_constituent_not_independently_testable",
    "cross_cohort_nonlinear_shape_concordance"
  )]
  setnames(q_values, evaluation, paste0("bh_q_", evaluation))
  q_values
})
timing_shape <- rbindlist(timing_shape_rows, use.names = TRUE)
ml_write_tsv_once(
  timing_cohort,
  file.path(temporary_out, "continuum_gene_nonlinear_timing_cohort.tsv.gz")
)
ml_write_tsv_once(
  timing_shape,
  file.path(temporary_out, "continuum_gene_nonlinear_shape_concordance.tsv.gz")
)

message("Fitting the prespecified GSE162694 NAS-replacement sensitivity")
nas_rows <- list()
nas_cohort <- "GSE162694"
for (axis_name in co_primary) {
  nas_meta <- merge(
    manifest[dataset == nas_cohort],
    axes[dataset == nas_cohort & axis_id == axis_name, .(sample_id, axis_raw)],
    by = "sample_id", all = FALSE, sort = FALSE
  )
  nas_meta <- nas_meta[is.finite(axis_raw) & !is.na(nas_score) & !is.na(inferred_sex)]
  nas_meta[, `:=`(
    axis_z = ml_standardize(axis_raw),
    nas_factor = factor(nas_score),
    sex_factor = factor(inferred_sex)
  )]
  ml_assert(nrow(nas_meta) >= 20L && uniqueN(nas_meta$nas_factor) >= 2L &&
              uniqueN(nas_meta$sex_factor) >= 2L,
            "GSE162694 has insufficient complete NAS/sex continuum data")
  nas_design <- model.matrix(~ nas_factor + sex_factor + axis_z, data = nas_meta)
  rownames(nas_design) <- nas_meta$sample_id
  fitted <- bulk_voom_fit(
    dge, nas_meta$sample_id, nas_design, "axis_z", gene_annotation,
    nas_cohort, "continuum_nas_sex_adjusted_sensitivity", axis_name
  )
  result <- fitted$table
  result[, `:=`(
    axis_component = gene_id_base %in% signature_ids,
    inference_eligible = !gene_id_base %in% signature_ids,
    p_value = fifelse(gene_id_base %in% signature_ids, NA_real_, raw_p_value)
  )]
  result[, bh_q_value := ml_complete_bh(p_value, contract$gene_family_size)]
  result[, `:=`(
    bh_family_size = contract$gene_family_size,
    model = "expression ~ continuum_z + factor(NAS) + inferred_sex",
    sensitivity_role = "NAS_replaces_fibrosis_GSE162694_only",
    replication_unit = "participant",
    overlaps_continuum_score_construction = TRUE,
    overlaps_canonical_resource_discovery = TRUE,
    overlaps_source_overlap_discovery = FALSE,
    overlaps_paired_replication = FALSE,
    eligible_as_independent_validation_of_canonical_degs = FALSE,
    evidence_role = fifelse(
      axis_component,
      "axis_constituent_not_independently_testable",
      "continuum_association_NAS_sex_adjusted_sensitivity"
    )
  )]
  nas_rows[[axis_name]] <- result
}
nas_sensitivity <- rbindlist(nas_rows, use.names = TRUE)
ml_assert(nrow(nas_sensitivity) == contract$gene_family_size * length(co_primary),
          "NAS-replacement sensitivity ledger is incomplete")
ml_write_tsv_once(
  nas_sensitivity,
  file.path(temporary_out, "continuum_gene_nas_sensitivity.tsv.gz")
)

message("Building target-specific leave-one-gene-out axes for 139 signature genes")
cohort_loadings <- fread(loading_path)
projection_loadings <- fread(projection_loading_path)
loo_rows <- list()
for (cohort in evaluation) {
  expression <- readRDS(file.path(normalized_dir, paste0(cohort, ".rds")))
  ml_assert(is.matrix(expression) && setequal(colnames(expression),
                                               manifest[dataset == cohort, sample_id]),
            paste0(cohort, " normalized expression schema drift"))
  cohort_signature <- intersect(signature_ids, rownames(expression))
  ml_assert(length(cohort_signature) == contract$signature_observed_size,
            paste0(cohort, " lacks an observed signature gene"))
  loading_subset <- cohort_loadings[dataset == cohort]
  for (axis_name in co_primary) {
    cache <- voom_cache[[paste(cohort, axis_name)]]
    model_meta <- cache$metadata
    for (target in cohort_signature) {
      axis <- bulk_loo_axis(
        axis_name, target, expression, loading_subset, projection_loadings
      )
      axis_z <- ml_standardize(axis[model_meta$sample_id])
      loo_design <- model.matrix(
        ~ fibrosis_factor + sex_factor + axis_z,
        data = model_meta
      )
      rownames(loo_design) <- model_meta$sample_id
      target_row <- match(target, ml_base_gene_id(rownames(cache$voom$E)))
      ml_assert(!is.na(target_row), paste0("LOO target missing from voom object: ", target))
      coefficient <- bulk_weighted_coefficient(
        cache$voom$E[target_row, ], loo_design,
        cache$voom$weights[target_row, ], "axis_z"
      )
      loo_rows[[length(loo_rows) + 1L]] <- cbind(
        data.table(
          gene_id_versioned = rownames(cache$voom$E)[target_row],
          gene_id_base = target,
          gene_name = gene_annotation$gene_name[match(target, gene_annotation$gene_id_base)],
          dataset = cohort,
          axis_id = axis_name,
          analysis = "axis_component_leave_one_out",
          axis_target_excluded = TRUE,
          model = "voom_expression ~ leave_one_gene_out_continuum_z + factor(fibrosis_stage) + inferred_sex",
          evidence_role = "axis_component_leave_one_out_not_independent_validation",
          replication_unit = "participant",
          overlaps_continuum_score_construction = TRUE,
          overlaps_canonical_resource_discovery = TRUE,
          overlaps_source_overlap_discovery = FALSE,
          overlaps_paired_replication = FALSE,
          eligible_as_independent_validation_of_canonical_degs = FALSE
        ),
        coefficient
      )
    }
  }
}
loo_cohort <- rbindlist(loo_rows, use.names = TRUE, fill = TRUE)
ml_assert(nrow(loo_cohort) == contract$signature_observed_size * length(evaluation) *
            length(co_primary), "LOO cohort ledger is incomplete")
loo_meta <- loo_cohort[, {
  fixed <- ml_fixed_meta(beta, se)
  list(
    estimable = fixed$estimable, beta = fixed$beta, se = fixed$se,
    z_value = fixed$z, p_value = fixed$p_value,
    ci_low = fixed$ci_low, ci_high = fixed$ci_high,
    n_cohorts = fixed$n_cohorts,
    direction_concordant = uniqueN(sign(beta[is.finite(beta) & beta != 0])) == 1L
  )
}, by = .(gene_id_versioned, gene_id_base, gene_name, axis_id)]
loo_meta[, bh_q_value := ml_complete_bh(p_value, contract$signature_observed_size),
         by = axis_id]
loo_meta[, `:=`(
  bh_family_size = contract$signature_observed_size,
  membership_class = fifelse(
    direction_concordant & bh_q_value < 0.05,
    "axis_component_leave_one_out_supported",
    "axis_component_leave_one_out_not_supported"
  ),
  evidence_role = "axis_component_leave_one_out_not_independent_validation"
)]
ml_write_tsv_once(loo_cohort, file.path(temporary_out, "signature_loo_gene_cohort.tsv.gz"))
ml_write_tsv_once(loo_meta, file.path(temporary_out, "signature_loo_gene_meta.tsv"))

message("Refitting source-overlap disease-control discovery and held-out evaluation")
discovery_meta <- manifest[dataset %in% source_overlap &
                             !is.na(group_binary) & !is.na(inferred_sex)]
discovery_meta[, `:=`(
  dataset_factor = factor(dataset),
  sex_factor = factor(inferred_sex),
  group_factor = factor(group_binary, levels = c("Control", "Disease"))
)]
discovery_design <- model.matrix(
  ~ dataset_factor + sex_factor + group_factor, data = discovery_meta
)
rownames(discovery_design) <- discovery_meta$sample_id
discovery_fit <- bulk_voom_fit(
  dge, discovery_meta$sample_id, discovery_design, "group_factorDisease",
  gene_annotation, "source_overlap_pooled", "heldout_disease_discovery"
)$table
discovery_fit[, `:=`(
  p_value = raw_p_value,
  bh_q_value = ml_complete_bh(raw_p_value, contract$gene_family_size),
  bh_family_size = contract$gene_family_size,
  abs_log2fc_threshold = 0.5,
  heldout_signature_member = ml_complete_bh(raw_p_value, contract$gene_family_size) < 0.05 &
    abs(beta) > 0.5,
  discovery_cohorts = paste(source_overlap, collapse = ";"),
  evaluation_cohorts = paste(evaluation, collapse = ";"),
  evaluation_sample_overlap = FALSE,
  replication_unit = "participant",
  overlaps_continuum_score_construction = TRUE,
  overlaps_canonical_resource_discovery = TRUE,
  overlaps_source_overlap_discovery = TRUE,
  overlaps_paired_replication = FALSE,
  eligible_as_independent_validation_of_canonical_degs = FALSE,
  eligible_as_independent_validation_of_source_overlap_signature = FALSE
)]
ml_assert(sum(discovery_fit$heldout_signature_member) > 0L,
          "Source-overlap held-out signature is empty")

heldout_rows <- list()
heldout_voom <- list()
for (cohort in evaluation) {
  cohort_meta <- manifest[dataset == cohort & !is.na(group_binary) & !is.na(inferred_sex)]
  cohort_meta[, `:=`(
    sex_factor = factor(inferred_sex),
    group_factor = factor(group_binary, levels = c("Control", "Disease"))
  )]
  design <- model.matrix(~ sex_factor + group_factor, data = cohort_meta)
  rownames(design) <- cohort_meta$sample_id
  fitted <- bulk_voom_fit(
    dge, cohort_meta$sample_id, design, "group_factorDisease", gene_annotation,
    cohort, "heldout_disease_validation"
  )
  heldout_rows[[cohort]] <- fitted$table
  heldout_rows[[cohort]][, `:=`(
    overlaps_continuum_score_construction = TRUE,
    overlaps_canonical_resource_discovery = TRUE,
    overlaps_source_overlap_discovery = FALSE,
    overlaps_paired_replication = FALSE,
    eligible_as_independent_validation_of_canonical_degs = FALSE,
    eligible_as_independent_validation_of_source_overlap_signature = TRUE
  )]
  heldout_voom[[cohort]] <- list(voom = fitted$voom, metadata = cohort_meta)
}
heldout_cohort <- rbindlist(heldout_rows)
heldout_meta <- heldout_cohort[, {
  fixed <- ml_fixed_meta(beta, se)
  list(
    beta = fixed$beta, se = fixed$se, z_value = fixed$z,
    p_value = fixed$p_value, ci_low = fixed$ci_low, ci_high = fixed$ci_high,
    n_cohorts = fixed$n_cohorts,
    direction_concordant = uniqueN(sign(beta[is.finite(beta) & beta != 0])) == 1L
  )
}, by = .(gene_id_versioned, gene_id_base, gene_name)]
heldout_meta[, bh_q_value := ml_complete_bh(p_value, contract$gene_family_size)]
heldout_meta[, bh_family_size := contract$gene_family_size]

selected <- discovery_fit[heldout_signature_member == TRUE]
heldout_score_rows <- list()
heldout_score_models <- list()
for (cohort in evaluation) {
  cache <- heldout_voom[[cohort]]
  indices <- match(selected$gene_id_versioned, rownames(cache$voom$E))
  ml_assert(!anyNA(indices), "A discovery signature gene is absent from held-out expression")
  values <- cache$voom$E[indices, , drop = FALSE]
  z <- t(apply(values, 1L, ml_standardize))
  direction <- sign(selected$beta)
  score <- as.numeric(crossprod(direction / length(direction), z))
  score_table <- data.table(
    sample_id = colnames(values), dataset = cohort,
    heldout_disease_signature_score = score,
    n_signature_genes = length(direction),
    signature_frozen_before_evaluation = TRUE,
    evaluation_sample_overlap_with_discovery = FALSE,
    overlaps_continuum_score_construction = TRUE,
    overlaps_canonical_resource_discovery = TRUE,
    overlaps_source_overlap_discovery = FALSE,
    overlaps_paired_replication = FALSE,
    eligible_as_independent_validation_of_canonical_degs = FALSE,
    eligible_as_independent_validation_of_source_overlap_signature = TRUE
  )
  heldout_score_rows[[cohort]] <- score_table
  score_meta <- merge(cache$metadata, score_table, by = c("sample_id", "dataset"))
  score_fit <- lm(
    heldout_disease_signature_score ~ group_factor + sex_factor,
    data = score_meta
  )
  coefficient <- coef(summary(score_fit))["group_factorDisease", ]
  heldout_score_models[[cohort]] <- data.table(
    dataset = cohort, beta = coefficient[["Estimate"]],
    se = coefficient[["Std. Error"]], p_value = coefficient[["Pr(>|t|)"]],
    ci_low = coefficient[["Estimate"]] - qt(0.975, df.residual(score_fit)) *
      coefficient[["Std. Error"]],
    ci_high = coefficient[["Estimate"]] + qt(0.975, df.residual(score_fit)) *
      coefficient[["Std. Error"]],
    n_participants = nrow(score_meta), model = "heldout_score ~ disease_control + inferred_sex"
  )
}
heldout_scores <- rbindlist(heldout_score_rows)
score_models <- rbindlist(heldout_score_models)
score_meta <- ml_fixed_meta(score_models$beta, score_models$se)
score_models <- rbindlist(list(
  score_models[, analysis_level := "cohort"],
  data.table(
    dataset = "fixed_effect_meta", beta = score_meta$beta, se = score_meta$se,
    p_value = score_meta$p_value, ci_low = score_meta$ci_low,
    ci_high = score_meta$ci_high, n_participants = sum(score_models$n_participants),
    model = "fixed_effect_inverse_variance", analysis_level = "meta"
  )
), use.names = TRUE, fill = TRUE)

ml_write_tsv_once(discovery_fit, file.path(temporary_out, "heldout_discovery_gene_results.tsv.gz"))
ml_write_tsv_once(heldout_cohort, file.path(temporary_out, "heldout_validation_gene_cohort.tsv.gz"))
ml_write_tsv_once(heldout_meta, file.path(temporary_out, "heldout_validation_gene_meta.tsv.gz"))
ml_write_tsv_once(heldout_scores, file.path(temporary_out, "heldout_participant_scores.tsv.gz"))
ml_write_tsv_once(score_models, file.path(temporary_out, "heldout_score_validation.tsv"))

message("Writing canonical DEG and stage same-substrate concordance ledgers")
canonical <- fread(canonical_path)
bulk_required_columns(canonical, c("gene", "symbol", "logFC", "padj"),
                      "canonical F_five coefficient table")
canonical[, gene_id_base := ml_base_gene_id(gene)]
canonical[, canonical_deg := padj < 0.05 & abs(logFC) > 0.5]
ml_assert(nrow(canonical) == contract$gene_family_size &&
            sum(canonical$canonical_deg) == contract$canonical_deg_size,
          "Canonical 23,370/1,347 DEG contract drift")

canonical_join <- merge(
  continuum_meta,
  canonical[, .(
    gene_id_base, canonical_logFC = logFC, canonical_padj = padj,
    canonical_deg, canonical_symbol = symbol
  )],
  by = "gene_id_base", all.x = TRUE, sort = FALSE
)
canonical_join[, direction_agrees :=
                 sign(beta) == sign(canonical_logFC) & beta != 0 & canonical_logFC != 0]
canonical_join <- merge(
  canonical_join,
  membership[, .(gene_id_base, continuum_associated, membership_class)],
  by = "gene_id_base", all.x = TRUE, sort = FALSE
)
canonical_summary <- canonical_join[, .(
  n_genes = .N,
  spearman_all_genes = bulk_rank_cor(beta, canonical_logFC),
  n_canonical_degs = sum(canonical_deg),
  canonical_deg_direction_agreement_fraction = mean(direction_agrees[canonical_deg]),
  n_canonical_degs_continuum_associated = sum(canonical_deg & continuum_associated),
  n_continuum_associated = sum(continuum_associated)
), by = axis_id]
canonical_summary[, interpretation :=
  "same_substrate_concordance_not_independent_validation"]

stage <- fread(stage_path)
bulk_required_columns(stage, c("gene_id_base", "logFC", "axis", "contrast"),
                      "current fibrosis-stage transcript results")
stage_summary <- merge(
  continuum_meta[inference_eligible == TRUE, .(gene_id_base, axis_id, continuum_beta = beta)],
  stage[, .(gene_id_base = ml_base_gene_id(gene_id_base), stage_axis = axis,
            stage_contrast = contrast, stage_logFC = logFC)],
  by = "gene_id_base", allow.cartesian = TRUE
)[, .(
  n_genes = sum(is.finite(continuum_beta) & is.finite(stage_logFC)),
  spearman = bulk_rank_cor(continuum_beta, stage_logFC)
), by = .(axis_id, stage_axis, stage_contrast)]
stage_summary[, interpretation :=
  "same_substrate_stage_concordance_not_independent_validation"]

heldout_concordance <- merge(
  discovery_fit[, .(gene_id_base, discovery_beta = beta, discovery_member = heldout_signature_member)],
  heldout_meta[, .(gene_id_base, heldout_beta = beta, heldout_q = bh_q_value,
                   heldout_direction_concordant = direction_concordant)],
  by = "gene_id_base"
)
heldout_summary <- heldout_concordance[, .(
  n_genes = .N,
  spearman_all_genes = bulk_rank_cor(discovery_beta, heldout_beta),
  n_discovery_signature_genes = sum(discovery_member),
  discovery_signature_direction_agreement_fraction = mean(
    sign(discovery_beta[discovery_member]) == sign(heldout_beta[discovery_member])
  ),
  n_discovery_signature_heldout_bh_supported = sum(
    discovery_member & heldout_direction_concordant & heldout_q < 0.05
  )
)]
heldout_summary[, interpretation :=
  "heldout_evaluation_of_source_overlap_disease_signature"]

ml_write_tsv_once(canonical_join, file.path(temporary_out, "canonical_deg_continuum.tsv.gz"))
ml_write_tsv_once(canonical_summary, file.path(temporary_out, "canonical_deg_concordance.tsv"))
ml_write_tsv_once(stage_summary, file.path(temporary_out, "stage_continuum_concordance.tsv"))
ml_write_tsv_once(heldout_concordance, file.path(temporary_out, "heldout_gene_concordance.tsv.gz"))
ml_write_tsv_once(heldout_summary, file.path(temporary_out, "heldout_concordance_summary.tsv"))

run_summary <- data.table(
  metric = c(
    "gene_family_size", "signature_published", "signature_observed",
    "canonical_deg_count", "source_overlap_discovery_participants",
    "evaluation_participants", "continuum_associated_count",
    "nonlinear_timing_supported_axis_gene_pairs",
    "nas_sensitivity_gene_axis_rows"
  ),
  value = c(
    contract$gene_family_size, nrow(signature), length(signature_ids),
    sum(canonical$canonical_deg), nrow(discovery_meta),
    manifest[dataset %in% evaluation, .N], sum(membership$continuum_associated),
    sum(timing_shape$nonlinear_timing_supported, na.rm = TRUE),
    nrow(nas_sensitivity)
  )
)
ml_write_tsv_once(run_summary, file.path(temporary_out, "run_summary.tsv"))
ml_write_session_info(file.path(temporary_out, "sessionInfo.txt"))

ml_write_json_once(
  list(
    status = "BULK_TRANSCRIPT_LAYER_COMPLETE",
    release_state = "candidate_only",
    output_namespace = "bulk",
    gene_family_size = contract$gene_family_size,
    signature_observed_size = contract$signature_observed_size,
    canonical_deg_size = contract$canonical_deg_size,
    nonlinear_spline_df = 3L,
    nonlinear_shape_grid_percentiles = 101L,
    nas_replacement_cohort = "GSE162694",
    seed = contract$seed,
    claim_language = "continuum-associated transcripts; not continuum DEGs"
  ),
  file.path(temporary_out, "COMPLETE.json")
)
artifacts <- list.files(temporary_out, full.names = TRUE, recursive = TRUE)
artifacts <- artifacts[!file.info(artifacts)$isdir]
artifact_manifest <- data.table(
  artifact = substring(artifacts, nchar(temporary_out) + 2L),
  size_bytes = file.info(artifacts)$size,
  sha256 = vapply(artifacts, ml_sha256, character(1))
)
ml_write_tsv_once(artifact_manifest, file.path(temporary_out, "output_manifest.tsv"))
ml_assert(file.rename(temporary_out, final_out), "Atomic bulk namespace publication failed")
message("HAC_MOLECULAR_BULK_COMPLETE: ", final_out)
