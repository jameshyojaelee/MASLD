#!/usr/bin/env Rscript
# KEY MESSAGE: Molecular-system continuum associations are stress tested against a jointly powered cross-cohort permutation, size-matched competitive nulls, system redundancy, and the measurement coarseness of recorded fibrosis stage.

suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
  library(igraph)
  library(parallel)
})

options(digits = 17, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    paste(
      "Usage: 73_stress_test_molecular_systems.R",
      "MOLECULAR_ANALYSIS_CANDIDATE PROMOTED_ATLAS_SOURCE",
      "SYSTEM_INFERENCE_CANDIDATE OUTPUT_CANDIDATE"
    ),
    call. = FALSE
  )
}

analysis_candidate <- normalizePath(args[[1L]], mustWork = TRUE)
atlas_source <- normalizePath(args[[2L]], mustWork = TRUE)
inference_candidate <- normalizePath(args[[3L]], mustWork = TRUE)
output_candidate <- args[[4L]]
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
main_component_nodes <- 5439L
evaluation <- as.character(contract$evaluation_cohorts)
axes_primary <- as.character(contract$co_primary_axes)
collection_order <- ml_node_collection_order()

test_mode <- identical(Sys.getenv("HAC_ML_TEST_MODE", unset = "0"), "1")
permutation_replicates <- as.integer(contract$permutation_replicates)
competitive_draws <- 5000L
if (test_mode) {
  permutation_replicates <- as.integer(Sys.getenv("HAC_ML_PERMUTATIONS", unset = "200"))
  competitive_draws <- as.integer(Sys.getenv("HAC_ML_COMPETITIVE_DRAWS", unset = "50"))
}
ml_assert(
  test_mode || (permutation_replicates == 10000L && competitive_draws == 5000L),
  "Production stress tests must use 10,000 permutations and 5,000 competitive draws"
)
worker_count <- max(1L, as.integer(Sys.getenv("HAC_ML_WORKERS", unset = "1")))

# --------------------------------------------------------------------------
# Frozen, checksummed inputs
# --------------------------------------------------------------------------

community_path <- file.path(atlas_source, "geometry", "community_registry.tsv")
node_path <- file.path(atlas_source, "overlays", "all_node_stage_continuum_overlays.tsv.gz")
edge_path <- file.path(atlas_source, "geometry", "mutual_30nn_edges.tsv.gz")
atlas_checksums_path <- file.path(atlas_source, "promoted_source_checksums.tsv")
frozen_scores_path <- file.path(
  inference_candidate, "tables", "participant_system_scores.tsv.gz"
)
frozen_meta_path <- file.path(
  inference_candidate, "tables", "system_continuum_meta_analysis.tsv"
)
frozen_stage_meta_path <- file.path(
  inference_candidate, "tables", "system_fibrosis_stage_meta_analysis.tsv"
)
frozen_membership_path <- file.path(
  inference_candidate, "tables", "system_inference_membership.tsv"
)
frozen_permutation_path <- file.path(
  inference_candidate, "tables", "system_within_stage_permutations.tsv"
)
inference_checksums_path <- file.path(
  inference_candidate, "provenance", "output_checksums.tsv"
)
axis_paths <- c(
  file.path(ml_source_root(contract), "unsupervised", "participant_scores.tsv"),
  file.path(ml_source_root(contract), "projection", "participant_scores.tsv")
)
hotspot_path <- file.path(
  analysis_candidate, "programs", "hotspot", "hotspot_participant_scores.tsv.gz"
)
pathway_paths <- file.path(
  analysis_candidate, "pathway_tf", "pathways", collection_order[1:6],
  "donor_scores.tsv.gz"
)

required_inputs <- c(
  community_path, node_path, edge_path, atlas_checksums_path,
  frozen_scores_path, frozen_meta_path, frozen_stage_meta_path,
  frozen_membership_path, frozen_permutation_path, inference_checksums_path,
  axis_paths, pathway_paths, hotspot_path
)
ml_assert(all(file.exists(required_inputs)), "A stress-test input is missing")

atlas_checksums <- fread(atlas_checksums_path)
for (path in c(community_path, node_path)) {
  target_basename <- basename(path)
  expected <- atlas_checksums[
    basename(get("path")) == target_basename, unique(sha256)
  ]
  ml_assert(length(expected) == 1L,
            paste0("Promoted atlas checksum is absent or ambiguous: ", path))
  ml_assert(identical(ml_sha256(path), expected),
            paste0("Promoted atlas source hash drift: ", path))
}
inference_checksums <- fread(inference_checksums_path)
for (path in c(frozen_scores_path, frozen_meta_path, frozen_stage_meta_path,
               frozen_membership_path, frozen_permutation_path)) {
  target_basename <- basename(path)
  expected <- inference_checksums[
    basename(get("path")) == target_basename, unique(sha256)
  ]
  ml_assert(length(expected) == 1L,
            paste0("Frozen inference checksum is absent or ambiguous: ", path))
  ml_assert(identical(ml_sha256(path), expected),
            paste0("Frozen inference source hash drift: ", path))
}

input_manifest <- data.table(
  path = normalizePath(required_inputs, mustWork = TRUE),
  input_role = c(
    "frozen_system_registry", "frozen_system_node_membership",
    "frozen_mutual_30nn_graph", "promoted_atlas_checksum_registry",
    "frozen_participant_system_scores", "frozen_continuum_meta_analysis",
    "frozen_stage_meta_analysis", "frozen_system_membership",
    "frozen_within_stage_permutations", "frozen_inference_checksum_registry",
    "signature_pc1_scores", "fixed_projection_scores",
    rep("pathway_donor_scores", 6L), "hotspot_donor_scores_and_metadata"
  )
)
input_manifest[, sha256 := vapply(path, ml_sha256, character(1))]
ml_write_tsv_once(input_manifest, file.path(provenance_dir, "input_manifest.tsv"))

# --------------------------------------------------------------------------
# Frozen substrate
# --------------------------------------------------------------------------

communities <- fread(community_path, na.strings = c("", "NA"))
nodes <- fread(node_path, na.strings = c("", "NA"))[
  component_role == "main_system_component",
  .(feature_id, collection, community_id, system_display)
]
ml_assert(nrow(communities) == system_family_size,
          "Frozen system registry is not the complete 43-system family")
ml_assert(nrow(nodes) == main_component_nodes &&
            uniqueN(nodes$feature_id) == main_component_nodes,
          "Frozen main-component node family drift")
ml_assert(setequal(unique(nodes$collection), collection_order),
          "Frozen system collections drift")

hotspot <- fread(hotspot_path, na.strings = c("", "NA"))
metadata <- unique(hotspot[, .(
  sample_id, analysis_unit_id, dataset, inferred_sex, fibrosis_stage
)], by = c("sample_id", "dataset"))
ml_assert(nrow(metadata) == 844L && !anyDuplicated(metadata$analysis_unit_id),
          "Stress tests require the frozen 844-participant roster")

axes <- ml_load_axes(contract)[
  dataset %in% evaluation & axis_id %in% axes_primary
]
ml_assert(nrow(axes) == 501L * length(axes_primary),
          "Evaluation-cohort co-primary axis family drift")

feature_scores <- ml_load_node_donor_scores(analysis_candidate, nodes, hotspot)
ml_assert(uniqueN(feature_scores$feature_id) == main_component_nodes,
          "At least one frozen system node lacks donor scores")

# --------------------------------------------------------------------------
# Equivalence against the frozen 70 candidate
# --------------------------------------------------------------------------

rebuilt <- ml_collection_balanced_scores(feature_scores, "community_id")
frozen_scores <- fread(frozen_scores_path, na.strings = c("", "NA"))[
  aggregation_method == "collection_balanced",
  .(sample_id, dataset, community_id, frozen_raw = system_score_raw,
    frozen_z = system_score_z)
]
equivalence <- merge(
  rebuilt, frozen_scores,
  by = c("sample_id", "dataset", "community_id"), all = TRUE, sort = TRUE
)
ml_assert(nrow(equivalence) == 844L * system_family_size &&
            !anyNA(equivalence$system_score_raw) && !anyNA(equivalence$frozen_raw),
          "Rebuilt system scores do not align with the frozen 70 candidate")
maximum_score_difference <- max(abs(equivalence$system_score_raw - equivalence$frozen_raw))
ml_assert(maximum_score_difference < 1e-9,
          paste0("Rebuilt system scores diverge from the frozen candidate: ",
                 maximum_score_difference))

# --------------------------------------------------------------------------
# Donor-by-node score matrix and cohort analysis contexts
# --------------------------------------------------------------------------

donor_key <- unique(feature_scores[, .(sample_id, dataset)])
setorder(donor_key, dataset, sample_id)
donor_key[, donor_index := seq_len(.N)]
node_key <- unique(nodes[, .(feature_id, collection, community_id)])
setorder(node_key, collection, feature_id)
node_key[, node_index := seq_len(.N)]

wide <- merge(feature_scores, donor_key, by = c("sample_id", "dataset"), sort = FALSE)
wide <- merge(wide, node_key[, .(feature_id, collection, node_index)],
              by = c("feature_id", "collection"), sort = FALSE)
score_matrix <- matrix(
  NA_real_, nrow = nrow(node_key), ncol = nrow(donor_key),
  dimnames = list(node_key$feature_id, NULL)
)
score_matrix[cbind(wide$node_index, wide$donor_index)] <- wide$feature_score
ml_assert(!anyNA(score_matrix), "Node-by-donor score matrix is incomplete")
rm(wide)

collection_code <- as.integer(factor(node_key$collection, levels = collection_order))
community_code <- match(node_key$community_id, communities$community_id)
system_node_sets <- split(node_key$node_index, node_key$community_id)
system_node_sets <- system_node_sets[communities$community_id]
ml_assert(identical(lengths(system_node_sets, use.names = FALSE),
                    as.integer(communities$n_features)),
          "System node counts disagree with the frozen registry")

# Collection-balanced score for an arbitrary node set, restricted to donor columns.
set_scores <- function(node_indices, donor_columns) {
  block <- score_matrix[node_indices, donor_columns, drop = FALSE]
  groups <- split(seq_along(node_indices), collection_code[node_indices])
  medians <- vapply(
    groups,
    function(rows) colMedians(block[rows, , drop = FALSE]),
    numeric(length(donor_columns))
  )
  if (is.null(dim(medians))) medians <- matrix(medians, nrow = length(donor_columns))
  rowMeans(medians)
}

donor_meta <- merge(donor_key, metadata, by = c("sample_id", "dataset"), sort = FALSE)
setorder(donor_meta, donor_index)

build_context <- function(axis_name, cohort) {
  axis_values <- axes[axis_id == axis_name & dataset == cohort,
                      .(sample_id, dataset, axis_raw)]
  cohort_donors <- donor_key[dataset == cohort]
  setorder(cohort_donors, donor_index)
  frame <- merge(donor_meta[dataset == cohort], axis_values,
                 by = c("sample_id", "dataset"), sort = FALSE)
  frame <- frame[is.finite(axis_raw) & !is.na(fibrosis_stage) & !is.na(inferred_sex)]
  setorder(frame, donor_index)
  nuisance <- model.matrix(
    ~ factor(fibrosis_stage) + factor(inferred_sex), data = frame
  )
  nuisance_qr <- qr(nuisance)
  x_residual <- as.numeric(qr.resid(nuisance_qr, ml_standardize(frame$axis_raw)))
  x_residual <- x_residual - mean(x_residual)
  strata <- interaction(
    factor(frame$fibrosis_stage), factor(frame$inferred_sex), drop = TRUE
  )
  list(
    axis_id = axis_name, dataset = cohort,
    # Scores are built over every donor in the cohort and standardized there,
    # matching 70_infer_molecular_systems.R, then restricted to the
    # stage-complete analysis rows. Standardizing on the subset instead would
    # leave each cohort's z unchanged but rescale the inverse-variance meta.
    donor_columns = cohort_donors$donor_index,
    analysis_positions = match(frame$donor_index, cohort_donors$donor_index),
    n_cohort_donors = nrow(cohort_donors),
    n_participants = nrow(frame),
    nuisance_qr = nuisance_qr,
    nuisance_rank = ncol(nuisance),
    x_residual = x_residual,
    sxx = sum(x_residual^2),
    residual_df = nrow(frame) - ncol(nuisance) - 1L,
    fisher_weight = nrow(frame) - (ncol(nuisance) - 1L) - 3L,
    strata_groups = split(seq_len(nrow(frame)), strata)
  )
}

contexts <- list()
for (axis_name in axes_primary) {
  for (cohort in evaluation) {
    contexts[[paste(axis_name, cohort, sep = "|")]] <- build_context(axis_name, cohort)
  }
}
ml_assert(all(vapply(contexts, function(x) x$n_participants, integer(1)) > 20L),
          "A cohort analysis context is too small to fit")

# Standardize the whole-cohort score vector, then keep the analysis rows.
context_outcome <- function(context, y_cohort) {
  ml_standardize(y_cohort)[context$analysis_positions]
}

# Exact lm(y ~ x + nuisance) coefficient and standard error via residualization.
fit_fast <- function(context, y_cohort) {
  y <- context_outcome(context, y_cohort)
  y_residual <- as.numeric(qr.resid(context$nuisance_qr, y))
  y_residual <- y_residual - mean(y_residual)
  beta <- sum(context$x_residual * y_residual) / context$sxx
  rss <- sum(y_residual^2) - beta^2 * context$sxx
  se <- sqrt(rss / context$residual_df / context$sxx)
  c(beta = beta, se = se)
}

fixed_meta_z <- function(beta, se) {
  ok <- is.finite(beta) & is.finite(se) & se > 0
  if (sum(ok) < 2L) return(NA_real_)
  weight <- 1 / se[ok]^2
  estimate <- sum(weight * beta[ok]) / sum(weight)
  estimate / sqrt(1 / sum(weight))
}

# --------------------------------------------------------------------------
# Fast-path calibration against the frozen parametric fits
# --------------------------------------------------------------------------

observed_fits <- rbindlist(lapply(names(contexts), function(key) {
  context <- contexts[[key]]
  fits <- vapply(
    system_node_sets,
    function(ix) fit_fast(context, set_scores(ix, context$donor_columns)),
    numeric(2L)
  )
  data.table(
    community_id = communities$community_id,
    axis_id = context$axis_id, dataset = context$dataset,
    fast_beta = fits["beta", ], fast_se = fits["se", ]
  )
}))
frozen_meta <- fread(frozen_meta_path, na.strings = c("", "NA"))[
  aggregation_method == "collection_balanced",
  .(community_id, axis_id, frozen_z = fixed_z, frozen_beta = fixed_beta)
]
observed_meta <- observed_fits[, .(
  fast_z = fixed_meta_z(fast_beta, fast_se)
), by = .(community_id, axis_id)]
calibration <- merge(observed_meta, frozen_meta, by = c("community_id", "axis_id"))
ml_assert(nrow(calibration) == system_family_size * length(axes_primary),
          "Fast-path calibration family drift")
maximum_z_difference <- max(abs(calibration$fast_z - calibration$frozen_z))
ml_assert(maximum_z_difference < 1e-6,
          paste0("Fast-path meta z diverges from the frozen fit: ",
                 maximum_z_difference))
ml_write_tsv_once(calibration, file.path(validation_dir, "fast_path_calibration.tsv"))

# --------------------------------------------------------------------------
# 1. Joint cross-cohort meta-permutation
# --------------------------------------------------------------------------

joint_permutation <- function(axis_name, replicates, seed) {
  cohort_contexts <- contexts[paste(axis_name, evaluation, sep = "|")]
  prepared <- lapply(cohort_contexts, function(context) {
    y <- vapply(
      system_node_sets,
      function(ix) context_outcome(context, set_scores(ix, context$donor_columns)),
      numeric(context$n_participants)
    )
    y_residual <- qr.resid(context$nuisance_qr, y)
    y_residual <- sweep(y_residual, 2L, colMeans(y_residual), "-")
    list(
      context = context,
      y_residual = y_residual,
      syy = colSums(y_residual^2)
    )
  })
  combine <- function(x_list) {
    weight_total <- 0
    accumulated <- numeric(system_family_size)
    for (index in seq_along(prepared)) {
      item <- prepared[[index]]
      r <- as.numeric(crossprod(x_list[[index]], item$y_residual)) /
        sqrt(sum(x_list[[index]]^2) * item$syy)
      r <- pmin(pmax(r, -0.999999999), 0.999999999)
      weight <- item$context$fisher_weight
      accumulated <- accumulated + weight * atanh(r)
      weight_total <- weight_total + weight
    }
    accumulated / sqrt(weight_total)
  }
  observed_x <- lapply(prepared, function(item) item$context$x_residual)
  observed <- combine(observed_x)
  exceed <- integer(system_family_size)
  max_exceed <- integer(system_family_size)
  set.seed(seed)
  for (replicate_index in seq_len(replicates)) {
    permuted <- lapply(prepared, function(item) {
      values <- item$context$x_residual
      for (group in item$context$strata_groups) {
        values[group] <- sample(values[group], length(group), replace = FALSE)
      }
      values
    })
    null <- combine(permuted)
    exceed <- exceed + as.integer(abs(null) >= abs(observed))
    max_exceed <- max_exceed + as.integer(max(abs(null)) >= abs(observed))
  }
  data.table(
    community_id = communities$community_id,
    axis_id = axis_name,
    joint_fisher_z = observed,
    joint_permutation_p_value = (1 + exceed) / (1 + replicates),
    joint_maxT_fwer_p_value = (1 + max_exceed) / (1 + replicates),
    permutation_replicates = replicates,
    n_participants = sum(vapply(prepared, function(x) x$context$n_participants, integer(1))),
    permutation_strata = "dataset_by_fibrosis_stage_by_inferred_sex",
    two_sided = TRUE
  )
}

joint_rows <- list()
for (axis_index in seq_along(axes_primary)) {
  joint_rows[[axis_index]] <- joint_permutation(
    axes_primary[[axis_index]], permutation_replicates,
    contract$seed + axis_index * 700000L
  )
}
joint <- rbindlist(joint_rows)
joint[, joint_permutation_q_value := ml_complete_bh(
  joint_permutation_p_value, system_family_size
), by = axis_id]
joint <- merge(joint, communities[, .(community_id, system_display)],
               by = "community_id", all.x = TRUE)
setorder(joint, axis_id, community_id)
ml_assert(nrow(joint) == system_family_size * length(axes_primary),
          "Joint permutation family size drift")
ml_write_tsv_once(joint, file.path(table_dir, "system_joint_permutations.tsv"))

joint_membership <- joint[, .(
  joint_both_axes_directions = .N == length(axes_primary) &&
    length(unique(sign(joint_fisher_z))) == 1L,
  joint_both_axes_bh = .N == length(axes_primary) &&
    all(is.finite(joint_permutation_q_value) & joint_permutation_q_value < 0.05),
  joint_both_axes_maxT = .N == length(axes_primary) &&
    all(is.finite(joint_maxT_fwer_p_value) & joint_maxT_fwer_p_value < 0.05),
  maximum_joint_permutation_q = max(joint_permutation_q_value),
  maximum_joint_maxT_fwer_p = max(joint_maxT_fwer_p_value)
), by = .(community_id, system_display)]
joint_membership[, `:=`(
  joint_bh_supported = joint_both_axes_directions & joint_both_axes_bh,
  joint_maxT_supported = joint_both_axes_directions & joint_both_axes_maxT
)]

# --------------------------------------------------------------------------
# 2. Competitive nulls
# --------------------------------------------------------------------------

edges <- fread(edge_path, na.strings = c("", "NA"))
edges <- edges[from %in% node_key$feature_id & to %in% node_key$feature_id]
graph <- graph_from_data_frame(
  edges[, .(from, to)], directed = FALSE, vertices = node_key$feature_id
)
graph <- simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)
ml_assert(vcount(graph) == main_component_nodes,
          "Frozen graph vertex count drift")
ml_assert(count_components(graph) == 1L,
          "Frozen main component is not connected")
adjacency <- lapply(as_adj_list(graph, mode = "all"), as.integer)

draw_connected_set <- function(target_size) {
  ml_draw_connected_subgraph(adjacency, target_size)
}

draw_uniform_set <- function(target_size) {
  sample.int(main_component_nodes, target_size, replace = FALSE)
}

statistic_for_set <- function(node_indices) {
  vapply(axes_primary, function(axis_name) {
    beta <- numeric(length(evaluation))
    se <- numeric(length(evaluation))
    for (index in seq_along(evaluation)) {
      context <- contexts[[paste(axis_name, evaluation[[index]], sep = "|")]]
      fit <- fit_fast(context, set_scores(node_indices, context$donor_columns))
      beta[index] <- fit[["beta"]]
      se[index] <- fit[["se"]]
    }
    abs(fixed_meta_z(beta, se))
  }, numeric(1L))
}

competitive_for_system <- function(system_index) {
  set.seed(contract$seed + 900000L + system_index)
  node_indices <- system_node_sets[[system_index]]
  target_size <- length(node_indices)
  observed <- statistic_for_set(node_indices)
  results <- list()
  draw_rows <- list()
  for (null_name in c("connected_subgraph", "size_matched_uniform")) {
    exceed <- integer(length(axes_primary))
    finite_draws <- 0L
    composition <- numeric(length(collection_order))
    degree_total <- 0
    # Retained so the supplementary panel can show the null distribution
    # itself rather than only the exceedance count.
    null_statistics <- matrix(
      NA_real_, nrow = competitive_draws, ncol = length(axes_primary)
    )
    for (draw_index in seq_len(competitive_draws)) {
      drawn <- if (identical(null_name, "connected_subgraph")) {
        draw_connected_set(target_size)
      } else {
        draw_uniform_set(target_size)
      }
      if (is.null(drawn)) next
      statistic <- statistic_for_set(drawn)
      if (!all(is.finite(statistic))) next
      finite_draws <- finite_draws + 1L
      null_statistics[finite_draws, ] <- statistic
      exceed <- exceed + as.integer(statistic >= observed)
      composition <- composition +
        tabulate(collection_code[drawn], nbins = length(collection_order))
      degree_total <- degree_total + mean(lengths(adjacency[drawn]))
    }
    results[[null_name]] <- data.table(
      community_id = communities$community_id[[system_index]],
      null_model = null_name,
      axis_id = axes_primary,
      n_features = target_size,
      observed_absolute_meta_z = observed,
      competitive_p_value = (1 + exceed) / (1 + finite_draws),
      finite_draws = finite_draws,
      requested_draws = competitive_draws,
      null_median_absolute_meta_z = apply(
        null_statistics[seq_len(finite_draws), , drop = FALSE], 2L, median
      ),
      null_p95_absolute_meta_z = apply(
        null_statistics[seq_len(finite_draws), , drop = FALSE], 2L,
        quantile, probs = 0.95
      ),
      mean_null_hotspot_nodes = composition[[7L]] / max(1L, finite_draws),
      mean_null_degree = degree_total / max(1L, finite_draws),
      observed_mean_degree = mean(lengths(adjacency[node_indices]))
    )
    draw_rows[[null_name]] <- data.table(
      community_id = communities$community_id[[system_index]],
      null_model = null_name,
      axis_id = rep(axes_primary, each = finite_draws),
      draw_index = rep(seq_len(finite_draws), times = length(axes_primary)),
      null_absolute_meta_z = as.numeric(
        null_statistics[seq_len(finite_draws), , drop = FALSE]
      )
    )
  }
  list(summary = rbindlist(results), draws = rbindlist(draw_rows))
}

competitive_parts <- mclapply(
  seq_len(system_family_size), competitive_for_system, mc.cores = worker_count
)
ml_assert(all(vapply(competitive_parts, is.list, logical(1))),
          "A competitive-null worker failed to return results")
competitive <- rbindlist(lapply(competitive_parts, `[[`, "summary"))
competitive_draw_statistics <- rbindlist(
  lapply(competitive_parts, `[[`, "draws")
)
rm(competitive_parts)
ml_assert(nrow(competitive) ==
            system_family_size * 2L * length(axes_primary),
          "Competitive null family size drift")
ml_assert(all(competitive$finite_draws > 0L),
          "A competitive null produced no usable draws")
competitive[, competitive_q_value := ml_complete_bh(
  competitive_p_value, system_family_size
), by = .(null_model, axis_id)]
competitive <- merge(competitive, communities[, .(community_id, system_display)],
                     by = "community_id", all.x = TRUE)
setorder(competitive, null_model, axis_id, community_id)
ml_write_tsv_once(competitive, file.path(table_dir, "system_competitive_null.tsv"))

setorder(competitive_draw_statistics, null_model, axis_id, community_id, draw_index)
ml_assert(nrow(competitive_draw_statistics) == sum(competitive$finite_draws),
          "Retained competitive draw statistics do not match the usable draw counts")
ml_write_tsv_once(
  competitive_draw_statistics,
  file.path(table_dir, "system_competitive_null_draws.tsv.gz")
)

competitive_membership <- competitive[
  null_model == "connected_subgraph", .(
    connected_null_both_axes = .N == length(axes_primary) &&
      all(is.finite(competitive_q_value) & competitive_q_value < 0.05),
    maximum_connected_null_q = max(competitive_q_value)
  ), by = .(community_id, system_display)
]
uniform_membership <- competitive[
  null_model == "size_matched_uniform", .(
    uniform_null_both_axes = .N == length(axes_primary) &&
      all(is.finite(competitive_q_value) & competitive_q_value < 0.05),
    maximum_uniform_null_q = max(competitive_q_value)
  ), by = .(community_id, system_display)
]

# --------------------------------------------------------------------------
# 3. Redundancy across the 43 systems
# --------------------------------------------------------------------------

redundancy_rows <- list()
score_frame <- merge(rebuilt, donor_meta[, .(sample_id, dataset)],
                     by = c("sample_id", "dataset"), sort = FALSE)
for (cohort in c(evaluation, "pooled_centred")) {
  subset_frame <- if (identical(cohort, "pooled_centred")) {
    copy(score_frame[dataset %in% evaluation])
  } else {
    copy(score_frame[dataset == cohort])
  }
  subset_frame[, system_score_c := ml_standardize(system_score_raw),
               by = .(dataset, community_id)]
  matrix_form <- dcast(subset_frame, sample_id + dataset ~ community_id,
                       value.var = "system_score_c")
  values <- as.matrix(matrix_form[, communities$community_id, with = FALSE])
  correlation <- cor(values, use = "complete.obs")
  eigenvalues <- eigen(correlation, symmetric = TRUE, only.values = TRUE)$values
  redundancy_rows[[cohort]] <- data.table(
    scope = cohort,
    n_participants = nrow(values),
    n_systems = system_family_size,
    effective_independent_systems = ml_effective_tests(correlation),
    pc1_variance_explained = eigenvalues[[1L]] / sum(eigenvalues),
    pc1_to_pc3_variance_explained = sum(eigenvalues[1:3]) / sum(eigenvalues),
    mean_absolute_offdiagonal_correlation =
      mean(abs(correlation[upper.tri(correlation)])),
    interpretation_scope = "descriptive_redundancy_no_inferential_p_value"
  )
}
redundancy <- rbindlist(redundancy_rows)
ml_write_tsv_once(redundancy, file.path(table_dir, "system_redundancy.tsv"))

effect_vectors <- dcast(
  joint[, .(community_id, axis_id, joint_fisher_z)],
  community_id ~ axis_id, value.var = "joint_fisher_z"
)
effect_summary <- data.table(
  scope = "continuum_effect_vector",
  n_systems = system_family_size,
  positive_systems = sum(effect_vectors[[axes_primary[[1L]]]] > 0),
  axis_agreement_spearman = cor(
    effect_vectors[[axes_primary[[1L]]]], effect_vectors[[axes_primary[[2L]]]],
    method = "spearman"
  ),
  interpretation_scope = "descriptive_redundancy_no_inferential_p_value"
)
ml_write_tsv_once(effect_summary, file.path(table_dir, "system_effect_vector_summary.tsv"))

# --------------------------------------------------------------------------
# 4. Measurement-matched stage versus continuum comparison
# --------------------------------------------------------------------------

stage_context <- function(cohort) {
  frame <- donor_meta[dataset == cohort & !is.na(fibrosis_stage) & !is.na(inferred_sex)]
  setorder(frame, donor_index)
  frame
}

attenuation_rows <- list()
for (axis_name in axes_primary) {
  for (cohort in evaluation) {
    context <- contexts[[paste(axis_name, cohort, sep = "|")]]
    frame <- stage_context(cohort)
    axis_values <- axes[axis_id == axis_name & dataset == cohort,
                        .(sample_id, dataset, axis_raw)]
    frame <- merge(frame, axis_values, by = c("sample_id", "dataset"), sort = FALSE)
    frame <- frame[is.finite(axis_raw)]
    setorder(frame, donor_index)
    ml_assert(identical(frame$donor_index,
                        donor_key[dataset == cohort][context$analysis_positions,
                                                     donor_index]),
              "Attenuation frame does not match the fitted analysis rows")
    stage_levels <- sort(unique(frame$fibrosis_stage))
    stage_counts <- as.numeric(table(factor(frame$fibrosis_stage, levels = stage_levels)))
    breaks <- unique(c(0, cumsum(stage_counts) / sum(stage_counts)))
    binned <- cut(
      rank(frame$axis_raw, ties.method = "first") / nrow(frame),
      breaks = c(0, breaks[-1]), labels = FALSE, include.lowest = TRUE
    )
    sex_factor <- factor(frame$inferred_sex)
    stage_z <- ml_standardize(as.numeric(frame$fibrosis_stage))
    continuous_z <- ml_standardize(frame$axis_raw)
    binned_z <- ml_standardize(as.numeric(binned))
    for (system_index in seq_len(system_family_size)) {
      y <- context_outcome(
        context, set_scores(system_node_sets[[system_index]], context$donor_columns)
      )
      stage_fit <- coef(summary(lm(y ~ stage_z + sex_factor)))["stage_z", ]
      continuous_fit <- coef(summary(lm(y ~ continuous_z + sex_factor)))["continuous_z", ]
      binned_fit <- coef(summary(lm(y ~ binned_z + sex_factor)))["binned_z", ]
      attenuation_rows[[length(attenuation_rows) + 1L]] <- data.table(
        community_id = communities$community_id[[system_index]],
        axis_id = axis_name, dataset = cohort,
        n_participants = nrow(frame),
        n_stage_levels = length(stage_levels),
        marginal_stage_beta = stage_fit[[1L]],
        marginal_continuum_beta = continuous_fit[[1L]],
        measurement_matched_continuum_beta = binned_fit[[1L]]
      )
    }
  }
}
attenuation <- rbindlist(attenuation_rows)
attenuation[, `:=`(
  raw_continuum_to_stage_ratio = marginal_continuum_beta / marginal_stage_beta,
  matched_continuum_to_stage_ratio = measurement_matched_continuum_beta / marginal_stage_beta
)]
attenuation <- merge(attenuation, communities[, .(community_id, system_display)],
                     by = "community_id", all.x = TRUE)
setorder(attenuation, axis_id, dataset, community_id)
ml_assert(nrow(attenuation) ==
            system_family_size * length(axes_primary) * length(evaluation),
          "Attenuation family size drift")
ml_write_tsv_once(attenuation, file.path(table_dir, "stage_continuum_attenuation.tsv"))

same_sign <- attenuation[sign(marginal_stage_beta) == sign(marginal_continuum_beta)]
reliability_grid <- c(0.6, 0.7, 0.8, 0.9)
attenuation_summary <- rbindlist(lapply(reliability_grid, function(reliability) {
  data.table(
    assumed_stage_reliability = reliability,
    median_raw_ratio = median(same_sign$raw_continuum_to_stage_ratio),
    median_measurement_matched_ratio = median(same_sign$matched_continuum_to_stage_ratio),
    median_reliability_corrected_ratio =
      median(same_sign$matched_continuum_to_stage_ratio) * reliability,
    n_same_sign_rows = nrow(same_sign),
    interpretation_scope = "estimate_assumed_reliability_not_measured"
  )
}))
ml_write_tsv_once(
  attenuation_summary, file.path(table_dir, "stage_continuum_attenuation_summary.tsv")
)

# --------------------------------------------------------------------------
# Combined membership and summary
# --------------------------------------------------------------------------

frozen_membership <- fread(frozen_membership_path, na.strings = c("", "NA"))[, .(
  community_id, system_display, n_features, medoid_feature_label,
  representative_hotspot_programs, cross_layer_class,
  continuum_associated_system, fibrosis_stage_associated_system,
  permutation_bh_supported_all_four, permutation_maxT_supported_all_four
)]
membership <- Reduce(
  function(x, y) merge(x, y, by = c("community_id", "system_display"), all = TRUE),
  list(frozen_membership, joint_membership, competitive_membership, uniform_membership)
)
membership[, stress_test_class := fifelse(
  joint_maxT_supported & connected_null_both_axes,
  "joint_strict_and_competitive",
  fifelse(joint_maxT_supported, "joint_strict_only",
          fifelse(connected_null_both_axes, "competitive_only", "unsupported"))
)]
ml_assert(nrow(membership) == system_family_size, "Stress-test membership drift")
setorder(membership, community_id)
ml_write_tsv_once(membership, file.path(table_dir, "system_stress_membership.tsv"))

strict_comparison <- membership[, .N, by = .(
  frozen_four_test_maxT = permutation_maxT_supported_all_four,
  joint_maxT_supported, cross_layer_class
)]
setorder(strict_comparison, -N)
ml_write_tsv_once(strict_comparison, file.path(table_dir, "strict_criterion_comparison.tsv"))

summary_table <- data.table(
  metric = c(
    "molecular_system_family_size",
    "stage_adjusted_participants_GSE162694",
    "stage_adjusted_participants_GSE213621",
    "stage_adjusted_participants_total",
    "frozen_four_test_maxT_supported_systems",
    "joint_permutation_bh_supported_systems",
    "joint_permutation_maxT_supported_systems",
    "joint_maxT_supported_continuum_only_systems",
    "connected_subgraph_null_supported_systems",
    "size_matched_uniform_null_supported_systems",
    "joint_strict_and_competitive_systems",
    "effective_independent_systems_pooled",
    "pooled_pc1_variance_explained",
    "median_raw_continuum_to_stage_beta_ratio",
    "median_measurement_matched_continuum_to_stage_beta_ratio",
    "competitive_draws_per_system_per_null",
    "permutation_replicates_per_axis"
  ),
  value = c(
    as.character(system_family_size),
    as.character(contexts[[paste(axes_primary[[1L]], evaluation[[1L]], sep = "|")]]$n_participants),
    as.character(contexts[[paste(axes_primary[[1L]], evaluation[[2L]], sep = "|")]]$n_participants),
    as.character(sum(vapply(
      contexts[paste(axes_primary[[1L]], evaluation, sep = "|")],
      function(x) x$n_participants, integer(1)
    ))),
    as.character(sum(membership$permutation_maxT_supported_all_four)),
    as.character(sum(membership$joint_bh_supported)),
    as.character(sum(membership$joint_maxT_supported)),
    as.character(sum(membership$joint_maxT_supported &
                       membership$cross_layer_class == "continuum_supported_only")),
    as.character(sum(membership$connected_null_both_axes)),
    as.character(sum(membership$uniform_null_both_axes)),
    as.character(sum(membership$stress_test_class == "joint_strict_and_competitive")),
    as.character(redundancy[scope == "pooled_centred", effective_independent_systems]),
    as.character(redundancy[scope == "pooled_centred", pc1_variance_explained]),
    as.character(attenuation_summary$median_raw_ratio[[1L]]),
    as.character(attenuation_summary$median_measurement_matched_ratio[[1L]]),
    as.character(competitive_draws),
    as.character(permutation_replicates)
  ),
  interpretation_scope = c(
    "complete_family", rep("stage_adjusted_complete_case", 3L),
    "frozen_four_test_conjunction", "joint_cross_cohort_permutation",
    "joint_cross_cohort_fwer", "joint_cross_cohort_fwer",
    "competitive_conservative_reference", "competitive_permissive_reference",
    "conjunction_of_strict_and_competitive",
    rep("descriptive_redundancy_no_inferential_p_value", 2L),
    rep("measurement_matched_comparison", 2L),
    rep("sensitivity", 2L)
  )
)
ml_write_tsv_once(summary_table, file.path(table_dir, "system_stress_summary.tsv"))

# --------------------------------------------------------------------------
# Validation, provenance
# --------------------------------------------------------------------------

audit <- data.table(
  check = c(
    "rebuilt_scores_match_frozen_candidate",
    "fast_path_matches_frozen_parametric_fits",
    "joint_permutation_family_complete",
    "competitive_null_family_complete",
    "competitive_draws_all_usable",
    "competitive_null_distributions_retained",
    "attenuation_family_complete",
    "redundancy_scopes_complete",
    "no_outcome_used_to_define_geometry"
  ),
  detail = c(
    paste0("maximum absolute score difference ", maximum_score_difference),
    paste0("maximum absolute meta z difference ", maximum_z_difference),
    paste0(nrow(joint), " rows across ", length(axes_primary), " axes"),
    paste0(nrow(competitive), " rows across two null models"),
    paste0("minimum usable draws ", min(competitive$finite_draws)),
    paste0(nrow(competitive_draw_statistics), " retained null statistics"),
    paste0(nrow(attenuation), " rows"),
    paste0(nrow(redundancy), " scopes"),
    "geometry, community IDs, and node membership read from frozen checksummed sources only"
  ),
  passed = c(
    maximum_score_difference < 1e-9,
    maximum_z_difference < 1e-6,
    nrow(joint) == system_family_size * length(axes_primary),
    nrow(competitive) == system_family_size * 2L * length(axes_primary),
    min(competitive$finite_draws) > 0L,
    nrow(competitive_draw_statistics) == sum(competitive$finite_draws),
    nrow(attenuation) == system_family_size * length(axes_primary) * length(evaluation),
    nrow(redundancy) == length(evaluation) + 1L,
    TRUE
  )
)
ml_write_tsv_once(audit, file.path(validation_dir, "audit.tsv"))
ml_assert(all(audit$passed), "A stress-test validation check failed")

command_record <- data.table(
  field = c("script", "analysis_candidate", "atlas_source", "inference_candidate",
            "output_candidate", "seed", "permutation_replicates",
            "competitive_draws", "workers"),
  value = c(script_file, analysis_candidate, atlas_source, inference_candidate,
            normalizePath(output_candidate), as.character(contract$seed),
            as.character(permutation_replicates), as.character(competitive_draws),
            as.character(worker_count))
)
ml_write_tsv_once(command_record, file.path(provenance_dir, "command.tsv"))

output_files <- list.files(
  c(table_dir, validation_dir), full.names = TRUE, recursive = TRUE
)
ml_write_tsv_once(
  data.table(
    path = output_files,
    sha256 = vapply(output_files, ml_sha256, character(1))
  ),
  file.path(provenance_dir, "output_checksums.tsv")
)
ml_write_session_info(file.path(provenance_dir, "sessionInfo.txt"))

cat("Stress tests complete:", normalizePath(output_candidate), "\n")
