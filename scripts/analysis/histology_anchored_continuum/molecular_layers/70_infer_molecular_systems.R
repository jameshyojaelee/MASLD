#!/usr/bin/env Rscript
# KEY MESSAGE: Donor-level inference determines which frozen membership-defined molecular systems remain associated with the external continuum after recorded fibrosis stage and sex adjustment.

suppressPackageStartupMessages({
  library(data.table)
  library(splines)
})

options(digits = 17, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    paste(
      "Usage: 70_infer_molecular_systems.R",
      "MOLECULAR_ANALYSIS_CANDIDATE PROMOTED_ATLAS_SOURCE OUTPUT_CANDIDATE"
    ),
    call. = FALSE
  )
}

analysis_candidate <- normalizePath(args[[1L]], mustWork = TRUE)
atlas_source <- normalizePath(args[[2L]], mustWork = TRUE)
output_candidate <- args[[3L]]
if (file.exists(output_candidate)) {
  stop("Refusing to overwrite output candidate: ", output_candidate, call. = FALSE)
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "lib_molecular_layers.R"))
contract <- ml_read_contract()

table_dir <- file.path(output_candidate, "tables")
provenance_dir <- file.path(output_candidate, "provenance")
validation_dir <- file.path(output_candidate, "validation")
for (path in c(table_dir, provenance_dir, validation_dir)) ml_ensure_dir(path)

system_family_size <- 43L
evaluation <- as.character(contract$evaluation_cohorts)
axes_primary <- as.character(contract$co_primary_axes)
aggregation_methods <- c("collection_balanced", "feature_balanced_median")
collection_order <- c(
  "hallmark", "kegg", "reactome", "go_bp", "go_mf", "go_cc", "hotspot"
)

community_path <- file.path(atlas_source, "geometry", "community_registry.tsv")
node_path <- file.path(atlas_source, "overlays", "all_node_stage_continuum_overlays.tsv.gz")
atlas_checksums_path <- file.path(atlas_source, "promoted_source_checksums.tsv")
legacy_window_path <- file.path(
  atlas_source, "overlays", "community_balanced_fixed_windows.tsv.gz"
)
axis_paths <- c(
  file.path(ml_source_root(contract), "unsupervised", "participant_scores.tsv"),
  file.path(ml_source_root(contract), "projection", "participant_scores.tsv")
)
pathway_paths <- file.path(
  analysis_candidate, "pathway_tf", "pathways", collection_order[1:6],
  "donor_scores.tsv.gz"
)
hotspot_path <- file.path(
  analysis_candidate, "programs", "hotspot", "hotspot_participant_scores.tsv.gz"
)
required_inputs <- c(
  community_path, node_path, atlas_checksums_path, legacy_window_path, axis_paths,
  pathway_paths, hotspot_path
)
ml_assert(all(file.exists(required_inputs)), "A molecular-system inference input is missing")

atlas_checksums <- fread(atlas_checksums_path)
for (path in c(community_path, node_path, legacy_window_path)) {
  target_basename <- basename(path)
  expected <- atlas_checksums[
    basename(get("path")) == target_basename, unique(sha256)
  ]
  ml_assert(length(expected) == 1L,
            paste0("Promoted atlas checksum is absent or ambiguous: ", path))
  ml_assert(identical(ml_sha256(path), expected),
            paste0("Promoted atlas source hash drift: ", path))
}

input_manifest <- data.table(
  path = normalizePath(required_inputs, mustWork = TRUE),
  input_role = c(
    "frozen_system_registry", "frozen_system_node_membership",
    "promoted_atlas_checksum_registry", "legacy_feature_first_windows",
    "signature_pc1_scores",
    "fixed_projection_scores", rep("pathway_donor_scores", 6L),
    "hotspot_donor_scores_and_metadata"
  )
)
input_manifest[, sha256 := vapply(path, ml_sha256, character(1))]
input_manifest[, hash_frozen_before_outcome_read := TRUE]
ml_write_tsv_once(input_manifest, file.path(provenance_dir, "input_manifest.tsv"))

communities <- fread(community_path, na.strings = c("", "NA"))
nodes <- fread(node_path, na.strings = c("", "NA"))[
  component_role == "main_system_component",
  .(feature_id, collection, community_id, system_display)
]
ml_assert(nrow(communities) == system_family_size &&
            uniqueN(communities$community_id) == system_family_size,
          "Frozen system registry is not the complete 43-system family")
ml_assert(nrow(nodes) == 5439L && uniqueN(nodes$feature_id) == 5439L,
          "Frozen main-component node family drift")
ml_assert(setequal(unique(nodes$collection), collection_order),
          "Frozen system collections drift")

hotspot <- fread(hotspot_path, na.strings = c("", "NA"))
required_hotspot <- c(
  "program_uid", "sample_id", "analysis_unit_id", "dataset", "inferred_sex",
  "fibrosis_stage", "nas_score", "outcome_z"
)
ml_assert(all(required_hotspot %in% names(hotspot)),
          "Hotspot donor-score schema drift")
metadata <- unique(hotspot[, .(
  sample_id, analysis_unit_id, dataset, inferred_sex, fibrosis_stage, nas_score
)], by = c("sample_id", "dataset"))
ml_assert(nrow(metadata) == 844L && !anyDuplicated(metadata$analysis_unit_id),
          "System inference requires the frozen 844-participant roster")
ml_assert(metadata[dataset %in% evaluation, uniqueN(sample_id)] == 501L,
          "Evaluation-cohort participant census drift")

score_rows <- vector("list", length(collection_order))
for (index in seq_len(6L)) {
  collection_name <- collection_order[[index]]
  values <- fread(pathway_paths[[index]], na.strings = c("", "NA"))
  ml_assert(all(c("set_id", "sample_id", "pathway_score", "dataset") %in% names(values)),
            paste0("Pathway donor-score schema drift: ", collection_name))
  values <- values[
    set_id %in% nodes[collection == collection_name, feature_id],
    .(
      feature_id = set_id, collection = collection_name,
      sample_id, dataset, feature_score = as.numeric(pathway_score)
    )
  ]
  score_rows[[index]] <- values
}
score_rows[[7L]] <- hotspot[
  program_uid %in% nodes[collection == "hotspot", feature_id],
  .(
    feature_id = program_uid, collection = "hotspot",
    sample_id, dataset, feature_score = as.numeric(outcome_z)
  )
]
feature_scores <- rbindlist(score_rows, use.names = TRUE)
rm(score_rows)
feature_scores <- merge(
  feature_scores, nodes,
  by = c("feature_id", "collection"), all = FALSE, sort = FALSE
)
ml_assert(uniqueN(feature_scores$feature_id) == 5439L,
          "At least one frozen system node lacks donor scores")
ml_assert(!anyDuplicated(feature_scores[, .(feature_id, sample_id, dataset)]),
          "Feature donor scores are duplicated")
ml_assert(all(is.finite(feature_scores$feature_score)),
          "Feature donor scores contain non-finite values")

feature_score_scale_audit <- feature_scores[, .(
  mean_score = mean(feature_score), sd_score = sd(feature_score),
  n_participants = .N
), by = .(feature_id, collection, dataset)]
ml_assert(all(feature_score_scale_audit$n_participants > 1L) &&
            all(is.finite(feature_score_scale_audit$sd_score)) &&
            all(feature_score_scale_audit$sd_score > 0),
          "A feature lacks donor-level score variation")

collection_components <- feature_scores[, .(
  collection_median = median(feature_score),
  n_features_in_collection = .N
), by = .(
  sample_id, dataset, community_id, system_display, collection
)]
collection_balanced <- collection_components[, .(
  system_score_raw = mean(collection_median),
  n_collections = .N,
  n_features = sum(n_features_in_collection)
), by = .(sample_id, dataset, community_id, system_display)]
collection_balanced[, aggregation_method := "collection_balanced"]

feature_balanced <- feature_scores[, .(
  system_score_raw = median(feature_score),
  n_collections = uniqueN(collection),
  n_features = .N
), by = .(sample_id, dataset, community_id, system_display)]
feature_balanced[, aggregation_method := "feature_balanced_median"]

system_scores <- rbindlist(
  list(collection_balanced, feature_balanced), use.names = TRUE
)
system_scores[, system_score_z := ml_standardize(system_score_raw),
              by = .(dataset, aggregation_method, community_id)]
system_scores <- merge(
  system_scores, metadata,
  by = c("sample_id", "dataset"), all.x = TRUE, sort = FALSE
)
setorder(system_scores, aggregation_method, community_id, dataset, sample_id)
expected_score_rows <- 844L * system_family_size * length(aggregation_methods)
ml_assert(nrow(system_scores) == expected_score_rows &&
            !anyDuplicated(system_scores[, .(
              aggregation_method, community_id, sample_id, dataset
            )]),
          "Participant-level system-score family is incomplete or duplicated")
ml_assert(all(is.finite(system_scores$system_score_z)) &&
            all(!is.na(system_scores$inferred_sex)),
          "Participant-level system scores or sex metadata are incomplete")
ml_write_tsv_once(
  system_scores,
  file.path(table_dir, "participant_system_scores.tsv.gz")
)
ml_write_tsv_once(
  collection_components,
  file.path(table_dir, "participant_system_collection_components.tsv.gz")
)
ml_write_tsv_once(
  feature_score_scale_audit,
  file.path(table_dir, "feature_score_scale_audit.tsv.gz")
)

# The promoted heatmaps summarize each feature within a window before taking
# collection medians. Inference instead requires one pre-outcome system score per
# participant. Quantify the non-commutativity of these aggregation orders so an
# inferential label is never silently attached to a different plotted estimand.
all_axes <- ml_load_axes(contract)
fixed_projection <- all_axes[axis_id == "fixed_projection"]
ml_assert(nrow(fixed_projection) == 844L &&
            !anyDuplicated(fixed_projection[, .(sample_id, dataset)]),
          "Fixed-projection roster is incomplete for window reconciliation")
window_score_rows <- vector("list", length(contract$window_centers))
for (window_index in seq_along(contract$window_centers)) {
  center <- as.numeric(contract$window_centers[[window_index]])
  lower <- center - as.numeric(contract$window_width) / 2
  upper <- center + as.numeric(contract$window_width) / 2
  selected_axes <- if (window_index == length(contract$window_centers)) {
    fixed_projection[axis_percentile >= lower & axis_percentile <= upper]
  } else {
    fixed_projection[axis_percentile >= lower & axis_percentile < upper]
  }
  donor_rows <- merge(
    system_scores[aggregation_method == "collection_balanced"],
    selected_axes[, .(sample_id, dataset)],
    by = c("sample_id", "dataset"), all = FALSE, sort = FALSE
  )
  window_score_rows[[window_index]] <- donor_rows[, .(
    donor_first_window_score = mean(system_score_raw),
    donor_first_se = sd(system_score_raw) / sqrt(.N),
    n_participants_donor_first = .N
  ), by = .(community_id, system_display, dataset)][, `:=`(
    window_id = window_index, center = center, lower = lower, upper = upper
  )]
}
donor_first_windows <- rbindlist(window_score_rows, use.names = TRUE)
legacy_windows <- fread(legacy_window_path, na.strings = c("", "NA"))[, .(
  community_id, system_display, dataset, window_id,
  feature_first_window_score = collection_balanced_score,
  n_participants_feature_first_min = n_participants_min,
  n_participants_feature_first_max = n_participants_max
)]
window_reconciliation <- merge(
  donor_first_windows, legacy_windows,
  by = c("community_id", "system_display", "dataset", "window_id"),
  all = TRUE, sort = FALSE
)
window_reconciliation[, `:=`(
  score_difference_donor_minus_feature_first =
    donor_first_window_score - feature_first_window_score,
  aggregation_level = "cohort_window"
)]
pooled_reconciliation <- window_reconciliation[, .(
  donor_first_window_score = mean(donor_first_window_score),
  feature_first_window_score = mean(feature_first_window_score),
  score_difference_donor_minus_feature_first =
    mean(donor_first_window_score) - mean(feature_first_window_score),
  n_cohorts = .N
), by = .(community_id, system_display, window_id, center, lower, upper)]
pooled_reconciliation[, `:=`(
  dataset = "equal_weight_five_cohort_pool",
  aggregation_level = "pooled_window"
)]
window_reconciliation_summary <- rbindlist(list(
  window_reconciliation[, .(
    comparison_scope = "all_cohort_windows",
    n_tiles = .N,
    pearson_r = cor(donor_first_window_score, feature_first_window_score),
    spearman_rho = cor(
      donor_first_window_score, feature_first_window_score, method = "spearman"
    ),
    maximum_absolute_difference = max(abs(
      score_difference_donor_minus_feature_first
    )),
    median_absolute_difference = median(abs(
      score_difference_donor_minus_feature_first
    ))
  )],
  pooled_reconciliation[, .(
    comparison_scope = "equal_weight_five_cohort_pool",
    n_tiles = .N,
    pearson_r = cor(donor_first_window_score, feature_first_window_score),
    spearman_rho = cor(
      donor_first_window_score, feature_first_window_score, method = "spearman"
    ),
    maximum_absolute_difference = max(abs(
      score_difference_donor_minus_feature_first
    )),
    median_absolute_difference = median(abs(
      score_difference_donor_minus_feature_first
    ))
  )]
), use.names = TRUE)
ml_assert(nrow(window_reconciliation) == 43L * 5L * 9L &&
            nrow(pooled_reconciliation) == 43L * 9L,
          "Window-estimand reconciliation family is incomplete")
ml_write_tsv_once(
  window_reconciliation,
  file.path(table_dir, "window_estimand_reconciliation.tsv.gz")
)
ml_write_tsv_once(
  pooled_reconciliation,
  file.path(table_dir, "pooled_window_estimand_reconciliation.tsv")
)
ml_write_tsv_once(
  window_reconciliation_summary,
  file.path(table_dir, "window_estimand_reconciliation_summary.tsv")
)

axes <- all_axes[
  dataset %in% evaluation & axis_id %in% axes_primary
]
ml_assert(nrow(axes) == 501L * length(axes_primary) &&
            !anyDuplicated(axes[, .(sample_id, dataset, axis_id)]),
          "Evaluation-cohort co-primary axis family drift")

fit_continuum <- function(data) {
  d <- copy(data)[
    is.finite(system_score_z) & is.finite(axis_raw) &
      !is.na(fibrosis_stage) & !is.na(inferred_sex)
  ]
  empty <- data.table(
    estimable = FALSE, failure_reason = "insufficient_complete_design",
    n_participants = nrow(d), beta = NA_real_, se = NA_real_,
    ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
    nonlinear_p_value = NA_real_
  )
  if (nrow(d) < 20L || uniqueN(d$fibrosis_stage) < 2L ||
      uniqueN(d$inferred_sex) < 2L) return(empty)
  d[, `:=`(
    axis_z = ml_standardize(axis_raw),
    stage_factor = factor(fibrosis_stage),
    sex_factor = factor(inferred_sex)
  )]
  linear <- tryCatch(
    lm(system_score_z ~ axis_z + stage_factor + sex_factor, data = d),
    error = function(e) e
  )
  if (inherits(linear, "error") ||
      !"axis_z" %in% rownames(coef(summary(linear)))) {
    empty[, failure_reason := if (inherits(linear, "error")) {
      conditionMessage(linear)
    } else {
      "axis_not_estimable"
    }]
    return(empty)
  }
  coefficient <- coef(summary(linear))["axis_z", ]
  nonlinear <- tryCatch(
    lm(
      system_score_z ~ ns(axis_z, df = 3) + stage_factor + sex_factor,
      data = d
    ),
    error = function(e) NULL
  )
  nonlinear_p <- NA_real_
  if (!is.null(nonlinear)) {
    comparison <- tryCatch(anova(linear, nonlinear), error = function(e) NULL)
    if (!is.null(comparison) && nrow(comparison) == 2L) {
      nonlinear_p <- comparison$`Pr(>F)`[[2L]]
    }
  }
  data.table(
    estimable = TRUE, failure_reason = NA_character_, n_participants = nrow(d),
    beta = unname(coefficient[["Estimate"]]),
    se = unname(coefficient[["Std. Error"]]),
    ci_low = unname(coefficient[["Estimate"]] -
                      qt(0.975, df.residual(linear)) * coefficient[["Std. Error"]]),
    ci_high = unname(coefficient[["Estimate"]] +
                       qt(0.975, df.residual(linear)) * coefficient[["Std. Error"]]),
    p_value = unname(coefficient[["Pr(>|t|)"]]),
    nonlinear_p_value = nonlinear_p
  )
}

continuum_model_rows <- list()
for (aggregation in aggregation_methods) {
  for (axis_name in axes_primary) {
    for (cohort in evaluation) {
      joined <- merge(
        system_scores[
          aggregation_method == aggregation & dataset == cohort
        ],
        axes[axis_id == axis_name & dataset == cohort,
             .(sample_id, dataset, axis_raw)],
        by = c("sample_id", "dataset"), all = FALSE, sort = FALSE
      )
      estimates <- joined[, fit_continuum(.SD), by = .(community_id, system_display)]
      estimates[, `:=`(
        aggregation_method = aggregation, axis_id = axis_name, dataset = cohort,
        model = "system_score_z ~ continuum_z + factor(fibrosis_stage) + inferred_sex",
        biological_unit = "participant"
      )]
      continuum_model_rows[[length(continuum_model_rows) + 1L]] <- estimates
    }
  }
}
continuum_models <- rbindlist(continuum_model_rows, use.names = TRUE, fill = TRUE)
continuum_models[, `:=`(
  q_value = ml_complete_bh(p_value, system_family_size),
  nonlinear_q_value = ml_complete_bh(nonlinear_p_value, system_family_size),
  bh_family_size = system_family_size
), by = .(aggregation_method, axis_id, dataset)]
ml_assert(nrow(continuum_models) ==
            length(aggregation_methods) * length(axes_primary) *
            length(evaluation) * system_family_size,
          "Continuum system-model family size drift")
ml_assert(all(continuum_models$estimable),
          "At least one primary system continuum model is not estimable")
ml_write_tsv_once(
  continuum_models,
  file.path(table_dir, "system_continuum_cohort_models.tsv")
)

continuum_meta <- continuum_models[, {
  fixed <- ml_fixed_meta(beta, se)
  random <- ml_random_meta(beta, se)
  .(
    estimable = fixed$estimable,
    fixed_beta = fixed$beta, fixed_se = fixed$se, fixed_z = fixed$z,
    fixed_p_value = fixed$p_value, fixed_ci_low = fixed$ci_low,
    fixed_ci_high = fixed$ci_high,
    random_beta = random$beta, random_se = random$se,
    random_p_value = random$p_value, random_ci_low = random$ci_low,
    random_ci_high = random$ci_high, random_tau2 = random$tau2,
    n_cohorts = fixed$n_cohorts,
    direction_concordant = .N == length(evaluation) &&
      (all(beta > 0) || all(beta < 0)),
    cohort_betas = paste(dataset, signif(beta, 7), sep = "=", collapse = ";")
  )
}, by = .(aggregation_method, axis_id, community_id, system_display)]
continuum_meta[, `:=`(
  fixed_q_value = ml_complete_bh(fixed_p_value, system_family_size),
  random_q_value = ml_complete_bh(random_p_value, system_family_size),
  bh_family_size = system_family_size
), by = .(aggregation_method, axis_id)]
ml_write_tsv_once(
  continuum_meta,
  file.path(table_dir, "system_continuum_meta_analysis.tsv")
)

fit_stage <- function(data) {
  d <- copy(data)[
    is.finite(system_score_z) & !is.na(fibrosis_stage) & !is.na(inferred_sex)
  ]
  empty <- data.table(
    estimable = FALSE, failure_reason = "insufficient_complete_design",
    n_participants = nrow(d), beta = NA_real_, se = NA_real_,
    ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_
  )
  if (nrow(d) < 20L || uniqueN(d$fibrosis_stage) < 2L ||
      uniqueN(d$inferred_sex) < 2L) return(empty)
  d[, `:=`(
    stage_z = ml_standardize(as.numeric(fibrosis_stage)),
    sex_factor = factor(inferred_sex)
  )]
  fit <- tryCatch(
    lm(system_score_z ~ stage_z + sex_factor, data = d),
    error = function(e) e
  )
  if (inherits(fit, "error") ||
      !"stage_z" %in% rownames(coef(summary(fit)))) {
    empty[, failure_reason := if (inherits(fit, "error")) {
      conditionMessage(fit)
    } else {
      "stage_not_estimable"
    }]
    return(empty)
  }
  coefficient <- coef(summary(fit))["stage_z", ]
  data.table(
    estimable = TRUE, failure_reason = NA_character_, n_participants = nrow(d),
    beta = unname(coefficient[["Estimate"]]),
    se = unname(coefficient[["Std. Error"]]),
    ci_low = unname(coefficient[["Estimate"]] -
                      qt(0.975, df.residual(fit)) * coefficient[["Std. Error"]]),
    ci_high = unname(coefficient[["Estimate"]] +
                       qt(0.975, df.residual(fit)) * coefficient[["Std. Error"]]),
    p_value = unname(coefficient[["Pr(>|t|)"]])
  )
}

stage_models <- system_scores[
  aggregation_method %in% aggregation_methods & dataset %in% evaluation,
  fit_stage(.SD),
  by = .(aggregation_method, dataset, community_id, system_display)
]
stage_models[, `:=`(
  q_value = ml_complete_bh(p_value, system_family_size),
  bh_family_size = system_family_size,
  model = "system_score_z ~ fibrosis_stage_z + inferred_sex",
  biological_unit = "participant"
), by = .(aggregation_method, dataset)]
ml_assert(nrow(stage_models) ==
            length(aggregation_methods) * length(evaluation) * system_family_size &&
            all(stage_models$estimable),
          "Fibrosis-stage system-model family size or estimability drift")
ml_write_tsv_once(
  stage_models,
  file.path(table_dir, "system_fibrosis_stage_cohort_models.tsv")
)

stage_meta <- stage_models[, {
  fixed <- ml_fixed_meta(beta, se)
  random <- ml_random_meta(beta, se)
  .(
    estimable = fixed$estimable,
    fixed_beta = fixed$beta, fixed_se = fixed$se, fixed_z = fixed$z,
    fixed_p_value = fixed$p_value, fixed_ci_low = fixed$ci_low,
    fixed_ci_high = fixed$ci_high,
    random_beta = random$beta, random_se = random$se,
    random_p_value = random$p_value, random_tau2 = random$tau2,
    n_cohorts = fixed$n_cohorts,
    direction_concordant = .N == length(evaluation) &&
      (all(beta > 0) || all(beta < 0)),
    cohort_betas = paste(dataset, signif(beta, 7), sep = "=", collapse = ";")
  )
}, by = .(aggregation_method, community_id, system_display)]
stage_meta[, `:=`(
  fixed_q_value = ml_complete_bh(fixed_p_value, system_family_size),
  random_q_value = ml_complete_bh(random_p_value, system_family_size),
  bh_family_size = system_family_size
), by = aggregation_method]
ml_write_tsv_once(
  stage_meta,
  file.path(table_dir, "system_fibrosis_stage_meta_analysis.tsv")
)

permutation_family <- function(score_data, axis_data, replicates, seed) {
  system_ids <- sort(unique(score_data$community_id))
  wide <- dcast(
    score_data,
    sample_id + dataset + fibrosis_stage + inferred_sex ~ community_id,
    value.var = "system_score_z"
  )
  joined <- merge(
    wide, axis_data[, .(sample_id, dataset, axis_raw)],
    by = c("sample_id", "dataset"), all = FALSE, sort = FALSE
  )
  complete <- is.finite(joined$axis_raw) & !is.na(joined$fibrosis_stage) &
    !is.na(joined$inferred_sex) & complete.cases(joined[, ..system_ids])
  joined <- joined[complete]
  x <- ml_standardize(joined$axis_raw)
  y <- as.matrix(joined[, ..system_ids])
  nuisance <- model.matrix(
    ~ factor(fibrosis_stage) + factor(inferred_sex), data = joined
  )
  nuisance_qr <- qr(nuisance)
  x_residual <- as.numeric(qr.resid(nuisance_qr, x))
  y_residual <- qr.resid(nuisance_qr, y)
  x_residual <- x_residual - mean(x_residual)
  y_residual <- sweep(y_residual, 2L, colMeans(y_residual), "-")
  denominator <- sqrt(sum(x_residual^2) * colSums(y_residual^2))
  observed <- as.numeric(crossprod(x_residual, y_residual)) / denominator
  strata <- interaction(
    factor(joined$fibrosis_stage), factor(joined$inferred_sex), drop = TRUE
  )
  groups <- split(seq_along(x_residual), strata)
  exceed <- integer(length(system_ids))
  max_exceed <- integer(length(system_ids))
  finite_replicates <- 0L
  set.seed(seed)
  for (replicate_index in seq_len(replicates)) {
    permuted <- x_residual
    for (group in groups) {
      permuted[group] <- sample(x_residual[group], length(group), replace = FALSE)
    }
    null <- as.numeric(crossprod(permuted, y_residual)) / denominator
    if (all(is.finite(null))) {
      finite_replicates <- finite_replicates + 1L
      exceed <- exceed + as.integer(abs(null) >= abs(observed))
      maximum <- max(abs(null))
      max_exceed <- max_exceed + as.integer(maximum >= abs(observed))
    }
  }
  data.table(
    community_id = system_ids,
    n_participants = nrow(joined),
    residual_pearson = observed,
    permutation_p_value = (1 + exceed) / (1 + finite_replicates),
    maxT_fwer_p_value = (1 + max_exceed) / (1 + finite_replicates),
    permutation_replicates = replicates,
    n_finite_permutations = finite_replicates,
    permutation_strata = "fibrosis_stage_by_inferred_sex",
    two_sided = TRUE
  )
}

requested_replicates <- as.integer(contract$permutation_replicates)
if (identical(Sys.getenv("HAC_ML_TEST_MODE", unset = "0"), "1")) {
  requested_replicates <- as.integer(Sys.getenv(
    "HAC_ML_PERMUTATIONS", unset = as.character(requested_replicates)
  ))
}
ml_assert(requested_replicates == 10000L ||
            identical(Sys.getenv("HAC_ML_TEST_MODE", unset = "0"), "1"),
          "Production system permutations must use 10,000 replicates")

permutation_rows <- list()
for (axis_index in seq_along(axes_primary)) {
  axis_name <- axes_primary[[axis_index]]
  for (cohort_index in seq_along(evaluation)) {
    cohort <- evaluation[[cohort_index]]
    result <- permutation_family(
      system_scores[
        aggregation_method == "collection_balanced" & dataset == cohort
      ],
      axes[axis_id == axis_name & dataset == cohort],
      requested_replicates,
      contract$seed + axis_index * 100000L + cohort_index * 10000L
    )
    result[, `:=`(
      axis_id = axis_name, dataset = cohort,
      aggregation_method = "collection_balanced"
    )]
    permutation_rows[[length(permutation_rows) + 1L]] <- result
  }
}
permutations <- rbindlist(permutation_rows, use.names = TRUE)
permutations[, permutation_q_value := ml_complete_bh(
  permutation_p_value, system_family_size
), by = .(axis_id, dataset)]
permutations <- merge(
  permutations,
  communities[, .(community_id, system_display)],
  by = "community_id", all.x = TRUE
)
ml_assert(nrow(permutations) ==
            length(axes_primary) * length(evaluation) * system_family_size &&
            all(permutations$n_finite_permutations == requested_replicates),
          "System permutation family size or replicate count drift")
ml_write_tsv_once(
  permutations,
  file.path(table_dir, "system_within_stage_permutations.tsv")
)

feature_balanced_membership <- continuum_meta[
  aggregation_method == "feature_balanced_median",
  .(
    feature_balanced_both_axes_q = .N == length(axes_primary) &&
      all(is.finite(fixed_q_value) & fixed_q_value < 0.05),
    feature_balanced_all_four_directions = .N == length(axes_primary) &&
      all(direction_concordant) && length(unique(sign(fixed_beta))) == 1L,
    feature_balanced_maximum_meta_q = max(fixed_q_value)
  ),
  by = .(community_id, system_display)
]
feature_balanced_membership[, feature_balanced_supported :=
  feature_balanced_both_axes_q & feature_balanced_all_four_directions]

random_effects_membership <- continuum_meta[
  aggregation_method == "collection_balanced",
  .(
    random_effects_both_axes_q = .N == length(axes_primary) &&
      all(is.finite(random_q_value) & random_q_value < 0.05),
    random_effects_maximum_meta_q = max(random_q_value)
  ),
  by = .(community_id, system_display)
]

permutation_membership <- permutations[, .(
  permutation_all_four_directions = .N ==
    length(axes_primary) * length(evaluation) &&
    length(unique(sign(residual_pearson))) == 1L,
  permutation_all_four_bh_q = .N ==
    length(axes_primary) * length(evaluation) &&
    all(is.finite(permutation_q_value) & permutation_q_value < 0.05),
  permutation_all_four_maxT = .N ==
    length(axes_primary) * length(evaluation) &&
    all(is.finite(maxT_fwer_p_value) & maxT_fwer_p_value < 0.05),
  maximum_permutation_q = max(permutation_q_value),
  maximum_maxT_fwer_p = max(maxT_fwer_p_value)
), by = .(community_id, system_display)]
permutation_membership[, `:=`(
  permutation_bh_supported_all_four =
    permutation_all_four_directions & permutation_all_four_bh_q,
  permutation_maxT_supported_all_four =
    permutation_all_four_directions & permutation_all_four_maxT
)]

continuum_membership <- continuum_meta[
  aggregation_method == "collection_balanced",
  .(
    both_axes_present = setequal(axis_id, axes_primary),
    both_axes_meta_q = .N == length(axes_primary) &&
      all(is.finite(fixed_q_value) & fixed_q_value < 0.05),
    all_four_cohort_axis_directions = .N == length(axes_primary) &&
      all(direction_concordant) && length(unique(sign(fixed_beta))) == 1L,
    maximum_meta_q = max(fixed_q_value),
    mean_meta_beta = mean(fixed_beta),
    axis_fixed_betas = paste(axis_id, signif(fixed_beta, 7), sep = "=", collapse = ";"),
    axis_fixed_q_values = paste(axis_id, signif(fixed_q_value, 7), sep = "=", collapse = ";")
  ),
  by = .(community_id, system_display)
]
continuum_membership[, continuum_associated_system :=
  both_axes_present & both_axes_meta_q & all_four_cohort_axis_directions]
continuum_membership[, evidence_label := fifelse(
  continuum_associated_system,
  "continuum_associated_system",
  "not_continuum_associated_system"
)]
stage_membership <- stage_meta[
  aggregation_method == "collection_balanced",
  .(
    stage_fixed_beta = fixed_beta,
    stage_fixed_q_value = fixed_q_value,
    stage_direction_concordant = direction_concordant,
    fibrosis_stage_associated_system = is.finite(fixed_q_value) &&
      fixed_q_value < 0.05 && direction_concordant
  ),
  by = .(community_id, system_display)
]
membership <- merge(
  continuum_membership, stage_membership,
  by = c("community_id", "system_display"), all = TRUE
)
membership <- merge(
  membership, feature_balanced_membership,
  by = c("community_id", "system_display"), all = TRUE
)
membership <- merge(
  membership, random_effects_membership,
  by = c("community_id", "system_display"), all = TRUE
)
membership <- merge(
  membership, permutation_membership,
  by = c("community_id", "system_display"), all = TRUE
)
membership <- merge(
  communities[, .(
    community_id, n_features, n_pathways, n_hotspot_programs,
    medoid_feature_id, medoid_feature_label, medoid_collection,
    representative_hotspot_programs, representative_hotspot_cell_types
  )],
  membership, by = "community_id", all.x = TRUE
)
membership[, cross_layer_class := fcase(
  continuum_associated_system & fibrosis_stage_associated_system,
  "stage_and_continuum_supported",
  continuum_associated_system,
  "continuum_supported_only",
  fibrosis_stage_associated_system,
  "stage_supported_only",
  default = "neither_supported"
)]
setorder(membership, -continuum_associated_system, -mean_meta_beta, community_id)
ml_assert(nrow(membership) == system_family_size,
          "System inference membership registry is not the complete 43-system family")
ml_write_tsv_once(
  membership,
  file.path(table_dir, "system_inference_membership.tsv")
)

summary <- data.table(
  metric = c(
    "participants_all_five_cohorts", "participants_primary_evaluation_cohorts",
    "molecular_system_family_size", "continuum_axes", "primary_cohorts",
    "collection_balanced_continuum_supported_systems",
    "collection_balanced_stage_supported_systems",
    "stage_and_continuum_supported_systems",
    "feature_balanced_continuum_supported_systems",
    "permutation_bh_supported_all_four_systems",
    "permutation_maxT_supported_all_four_systems",
    "permutation_replicates_per_axis_cohort",
    "primary_multiple_testing_family",
    "window_means_inferentially_tested",
    "pooled_window_donor_first_vs_feature_first_spearman"
  ),
  value = c(
    "844", "501", as.character(system_family_size),
    paste(axes_primary, collapse = ";"), paste(evaluation, collapse = ";"),
    as.character(sum(membership$continuum_associated_system)),
    as.character(sum(membership$fibrosis_stage_associated_system)),
    as.character(sum(membership$cross_layer_class == "stage_and_continuum_supported")),
    as.character(sum(membership$feature_balanced_supported)),
    as.character(sum(membership$permutation_bh_supported_all_four)),
    as.character(sum(membership$permutation_maxT_supported_all_four)),
    as.character(requested_replicates),
    "BH n=43 separately per aggregation method, axis, and model",
    "FALSE",
    as.character(window_reconciliation_summary[
      comparison_scope == "equal_weight_five_cohort_pool", spearman_rho
    ])
  ),
  interpretation_scope = c(
    "descriptive_score_roster", "primary_inference", "complete_family",
    "co_primary", "competitor_independent", "primary_inference",
    "secondary_stage_alignment", "cross_layer_summary",
    "aggregation_sensitivity", "permutation_sensitivity",
    "permutation_fwer_sensitivity", "sensitivity",
    "complete_family", "windows_remain_descriptive", "estimand_reconciliation"
  )
)
ml_write_tsv_once(summary, file.path(table_dir, "system_inference_summary.tsv"))

bh_recomputed <- continuum_meta[, .(
  maximum_absolute_difference = max(abs(
    fixed_q_value - ml_complete_bh(fixed_p_value, system_family_size)
  ))
), by = .(aggregation_method, axis_id)]
bh_pass <- all(bh_recomputed$maximum_absolute_difference < 1e-15)
participant_count_pass <- system_scores[, uniqueN(sample_id)] == 844L
complete_family_pass <-
  nrow(continuum_meta) == length(aggregation_methods) *
    length(axes_primary) * system_family_size &&
  nrow(stage_meta) == length(aggregation_methods) * system_family_size
permutation_pass <-
  nrow(permutations) == length(axes_primary) * length(evaluation) *
    system_family_size &&
  all(permutations$permutation_p_value >= 1 / (requested_replicates + 1)) &&
  all(permutations$maxT_fwer_p_value >= 1 / (requested_replicates + 1))

audit <- data.table(
  check = c(
    "input_hashes_frozen_before_outcomes", "complete_frozen_system_family",
    "one_score_per_participant_system", "primary_evaluation_roster",
    "complete_model_families", "exact_bh_family_size_43",
    "within_stage_sex_permutations", "both_axis_direction_gate",
    "window_estimands_reconciled", "window_means_not_tested",
    "deterministic_seed", "candidate_only"
  ),
  passed = c(
    all(input_manifest$hash_frozen_before_outcome_read),
    nrow(communities) == system_family_size && nrow(nodes) == 5439L,
    participant_count_pass && nrow(system_scores) == expected_score_rows,
    metadata[dataset %in% evaluation, uniqueN(sample_id)] == 501L,
    complete_family_pass, bh_pass, permutation_pass,
    all(continuum_membership$both_axes_present),
    nrow(window_reconciliation) == 43L * 5L * 9L &&
      nrow(pooled_reconciliation) == 43L * 9L,
    TRUE,
    identical(as.integer(contract$seed), 20260817L),
    !file.exists(file.path(output_candidate, "PROMOTED.ok"))
  ),
  detail = c(
    paste(basename(required_inputs), collapse = ";"),
    "43 systems; 5,439 main-component feature nodes",
    "844 participants x 43 systems x 2 aggregation methods",
    paste(evaluation, collapse = ";"),
    "43 systems x 2 axes x 2 cohorts x 2 aggregations; plus stage models",
    paste(bh_recomputed$aggregation_method, bh_recomputed$axis_id,
          format(bh_recomputed$maximum_absolute_difference, scientific = TRUE),
          sep = ":", collapse = ";"),
    paste0(requested_replicates, " permutations per axis-cohort; BH and maxT retained"),
    "Both co-primary meta q<0.05 and all four cohort-axis slopes share direction",
    paste(
      "Donor-first versus feature-first pooled-window Spearman =",
      signif(window_reconciliation_summary[
        comparison_scope == "equal_weight_five_cohort_pool", spearman_rho
      ], 6)
    ),
    "Nine overlapping fixed windows retain visualization_only status",
    as.character(contract$seed),
    normalizePath(output_candidate, mustWork = TRUE)
  )
)
ml_assert(all(audit$passed), "Molecular-system inference audit failed")
ml_write_tsv_once(audit, file.path(validation_dir, "audit.tsv"))
ml_write_session_info(file.path(provenance_dir, "sessionInfo.txt"))

command_record <- data.table(
  field = c(
    "script", "analysis_candidate", "atlas_source", "output_candidate",
    "seed", "permutation_replicates"
  ),
  value = c(
    normalizePath(script_file, mustWork = TRUE), analysis_candidate, atlas_source,
    normalizePath(output_candidate, mustWork = TRUE), as.character(contract$seed),
    as.character(requested_replicates)
  )
)
ml_write_tsv_once(command_record, file.path(provenance_dir, "command.tsv"))

output_files <- c(
  file.path(table_dir, "participant_system_scores.tsv.gz"),
  file.path(table_dir, "participant_system_collection_components.tsv.gz"),
  file.path(table_dir, "feature_score_scale_audit.tsv.gz"),
  file.path(table_dir, "window_estimand_reconciliation.tsv.gz"),
  file.path(table_dir, "pooled_window_estimand_reconciliation.tsv"),
  file.path(table_dir, "window_estimand_reconciliation_summary.tsv"),
  file.path(table_dir, "system_continuum_cohort_models.tsv"),
  file.path(table_dir, "system_continuum_meta_analysis.tsv"),
  file.path(table_dir, "system_fibrosis_stage_cohort_models.tsv"),
  file.path(table_dir, "system_fibrosis_stage_meta_analysis.tsv"),
  file.path(table_dir, "system_within_stage_permutations.tsv"),
  file.path(table_dir, "system_inference_membership.tsv"),
  file.path(table_dir, "system_inference_summary.tsv"),
  file.path(validation_dir, "audit.tsv"),
  file.path(provenance_dir, "input_manifest.tsv"),
  file.path(provenance_dir, "sessionInfo.txt"),
  file.path(provenance_dir, "command.tsv")
)
output_checksums <- data.table(
  path = normalizePath(output_files, mustWork = TRUE),
  sha256 = vapply(output_files, ml_sha256, character(1))
)
ml_write_tsv_once(output_checksums, file.path(provenance_dir, "output_checksums.tsv"))
writeLines("candidate_ready_for_statistical_review",
           file.path(output_candidate, "READY_FOR_REVIEW.ok"))
