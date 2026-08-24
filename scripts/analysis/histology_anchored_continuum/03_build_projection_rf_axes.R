#!/usr/bin/env Rscript
# Build the two presentation/benchmark axes that depend on the frozen discovery
# fit: an out-of-sample projection on the published 145-gene direction and a
# discovery-only random-forest fibrosis benchmark. The cohort-specific 145-gene
# PC1 is supplied by 02; its percentile and the projection percentile define the
# consensus. No Resource-cohort histology is read until every score is frozen.

suppressPackageStartupMessages({
  library(data.table)
  library(ranger)
})

options(digits = 17, scipen = 999)

SEED <- 20260817L
N_TREES <- 1000L
MTRY <- 11L
MIN_NODE_SIZE <- 5L
SIGNATURE_SIZE <- 145L
MIN_SIGNATURE_COVERAGE <- 0.90
EXPECTED_RESOURCE_SIGNATURE <- 139L
TARGET_COHORTS <- c("GSE126848", "GSE130970", "GSE135251",
                    "GSE162694", "GSE213621")
INDEPENDENT_COHORTS <- c("GSE162694", "GSE213621")

fail <- function(...) stop(paste0(...), call. = FALSE)

assert_true <- function(value, message) {
  if (!isTRUE(value)) fail(message)
  invisible(TRUE)
}

project_root <- Sys.getenv("MASLD_PROJECT_ROOT", "")
out_root <- Sys.getenv("HAC_OUT_ROOT", "")
assert_true(nzchar(project_root), "MASLD_PROJECT_ROOT is unset")
assert_true(nzchar(out_root), "HAC_OUT_ROOT is unset")
assert_true(dir.exists(out_root), paste0("HAC_OUT_ROOT does not exist: ", out_root))

output_dir <- file.path(out_root, "projection")
assert_true(!file.exists(output_dir),
            paste0("Refusing to overwrite an existing output namespace: ", output_dir))

write_tsv_once <- function(x, path) {
  assert_true(!file.exists(path), paste0("Refusing to overwrite output: ", path))
  fwrite(x, path, sep = "\t", na = "NA")
  invisible(path)
}

write_lines_once <- function(x, path) {
  assert_true(!file.exists(path), paste0("Refusing to overwrite output: ", path))
  writeLines(x, path)
  invisible(path)
}

save_rds_once <- function(x, path) {
  assert_true(!file.exists(path), paste0("Refusing to overwrite output: ", path))
  saveRDS(x, path)
  invisible(path)
}

sha256_file <- function(path) {
  assert_true(file.exists(path), paste0("Cannot hash missing input: ", path))
  out <- system2("sha256sum", shQuote(path), stdout = TRUE)
  assert_true(length(out) == 1L, paste0("sha256sum failed for ", path))
  sub("[[:space:]].*$", "", out)
}

rank_percentile <- function(x) {
  assert_true(all(is.finite(x)), "Cannot rank a score containing non-finite values")
  if (length(x) == 1L) return(0.5)
  (rank(x, ties.method = "average") - 1) / (length(x) - 1)
}

base_gene_ids <- function(x) sub("\\..*$", "", as.character(x))

validate_expression <- function(x, label) {
  assert_true(is.matrix(x) && is.numeric(x), paste0(label, " is not a numeric matrix"))
  assert_true(!is.null(rownames(x)) && !is.null(colnames(x)),
              paste0(label, " lacks gene or participant names"))
  genes <- base_gene_ids(rownames(x))
  assert_true(!anyDuplicated(genes), paste0(label, " has duplicate Ensembl base IDs"))
  assert_true(!anyDuplicated(colnames(x)), paste0(label, " has duplicate participant IDs"))
  assert_true(all(is.finite(x)), paste0(label, " contains non-finite expression values"))
  rownames(x) <- genes
  x
}

required_columns <- function(x, columns, label) {
  absent <- setdiff(columns, names(x))
  assert_true(!length(absent),
              paste0(label, " lacks required columns: ", paste(absent, collapse = ", ")))
}

evaluation_role <- function(dataset) {
  fifelse(dataset %chin% INDEPENDENT_COHORTS, "independent_validation",
          "source_overlap_replication")
}

reproduction_dir <- file.path(out_root, "reproduction")
unsupervised_dir <- file.path(out_root, "unsupervised")
model_path <- file.path(reproduction_dir, "discovery_signature_model.rds")
discovery_expression_path <- file.path(reproduction_dir, "discovery_normalized_expression.rds")
discovery_metadata_path <- file.path(reproduction_dir, "discovery_metadata.tsv")
signature_path <- file.path(reproduction_dir, "signature_genes.tsv")
unsupervised_scores_path <- file.path(unsupervised_dir, "participant_scores.tsv")
matrix_dir <- file.path(unsupervised_dir, "normalized_expression_by_cohort")
manifest_path <- Sys.getenv(
  "HAC_MANIFEST_PATH",
  Sys.getenv(
    "HAC_RESOURCE_MANIFEST",
    file.path(project_root,
              "figures/candidates/pi-figure-redesign-2026-08-13-v3/analysis/",
              "stage_extensions/five_cohort_sample_manifest.tsv")
  )
)

input_paths <- c(
  discovery_signature_model = model_path,
  discovery_normalized_expression = discovery_expression_path,
  discovery_metadata = discovery_metadata_path,
  signature_genes = signature_path,
  unsupervised_participant_scores = unsupervised_scores_path,
  resource_manifest = manifest_path,
  setNames(file.path(matrix_dir, paste0(TARGET_COHORTS, ".rds")),
           paste0("normalized_expression_", TARGET_COHORTS))
)
missing_inputs <- input_paths[!file.exists(input_paths)]
assert_true(!length(missing_inputs),
            paste0("Required inputs are absent: ", paste(missing_inputs, collapse = ", ")))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
assert_true(dir.exists(output_dir), paste0("Could not create output namespace: ", output_dir))

input_manifest <- data.table(
  input_role = names(input_paths),
  path = unname(input_paths),
  size_bytes = as.numeric(file.size(input_paths)),
  sha256 = vapply(input_paths, sha256_file, character(1))
)
write_tsv_once(input_manifest, file.path(output_dir, "input_checksums.tsv"))

# --------------------------------------------------------- frozen discovery --
model <- readRDS(model_path)
model_fields <- c("signature_gene_ids_all", "common_gene_ids", "discovery_center",
                  "discovery_loading_oriented")
assert_true(all(model_fields %in% names(model)),
            paste0("Discovery model lacks fields: ",
                   paste(setdiff(model_fields, names(model)), collapse = ", ")))

signature_all <- unique(base_gene_ids(model$signature_gene_ids_all))
common_model_genes <- unique(base_gene_ids(model$common_gene_ids))
assert_true(length(signature_all) == SIGNATURE_SIZE,
            sprintf("Expected %d unique signature genes; found %d",
                    SIGNATURE_SIZE, length(signature_all)))
assert_true(length(common_model_genes) == EXPECTED_RESOURCE_SIGNATURE,
            sprintf("Expected %d signature genes in the Resource; found %d",
                    EXPECTED_RESOURCE_SIGNATURE, length(common_model_genes)))
assert_true(all(common_model_genes %in% signature_all),
            "The discovery model's common genes are not a subset of the 145-gene signature")

center <- as.numeric(model$discovery_center)
loading <- as.numeric(model$discovery_loading_oriented)
names(center) <- base_gene_ids(names(model$discovery_center))
names(loading) <- base_gene_ids(names(model$discovery_loading_oriented))
assert_true(!anyDuplicated(names(center)) && !anyDuplicated(names(loading)),
            "Discovery center or loading has duplicate Ensembl base IDs")
assert_true(all(common_model_genes %in% names(center)) &&
              all(common_model_genes %in% names(loading)),
            "A common signature gene lacks a frozen center or oriented loading")
assert_true(all(is.finite(center[common_model_genes])) &&
              all(is.finite(loading[common_model_genes])),
            "Frozen centers or loadings contain non-finite values")
assert_true(sqrt(sum(loading[common_model_genes]^2)) > 0,
            "Frozen discovery loading has zero norm")

signature <- fread(signature_path)
required_columns(signature, c("gene_id_base", "in_resource"), "signature_genes.tsv")
signature[, gene_id_base := base_gene_ids(gene_id_base)]
assert_true(nrow(signature) == SIGNATURE_SIZE && !anyDuplicated(signature$gene_id_base),
            "signature_genes.tsv must contain exactly 145 unique genes")
assert_true(setequal(signature$gene_id_base, signature_all),
            "signature_genes.tsv and discovery model disagree on signature membership")
signature[, in_resource := as.logical(in_resource)]
assert_true(sum(signature$in_resource) == EXPECTED_RESOURCE_SIGNATURE,
            "signature_genes.tsv does not record exactly 139 Resource genes")
assert_true(sum(!signature$in_resource) == SIGNATURE_SIZE - EXPECTED_RESOURCE_SIGNATURE,
            "signature_genes.tsv does not record the expected six missing genes")

discovery_expression <- validate_expression(readRDS(discovery_expression_path),
                                            "discovery_normalized_expression.rds")
discovery_metadata <- fread(discovery_metadata_path)
required_columns(discovery_metadata,
                 c("sample_id", "dataset", "sex", "NAS", "fibrosis"),
                 "discovery_metadata.tsv")
assert_true(nrow(discovery_metadata) == 135L,
            sprintf("RF training requires exactly 135 discovery donors; found %d",
                    nrow(discovery_metadata)))
assert_true(!anyDuplicated(discovery_metadata$sample_id),
            "Discovery metadata has duplicate participants")
assert_true(setequal(discovery_metadata$sample_id, colnames(discovery_expression)),
            "Discovery expression and metadata participant sets differ")
assert_true(all(is.finite(discovery_metadata$fibrosis)) &&
              all(discovery_metadata$fibrosis %in% 0:4),
            "Discovery fibrosis must be complete numeric F0-F4")
discovery_expression <- discovery_expression[, discovery_metadata$sample_id, drop = FALSE]

# ---------------------------------------------------------- Resource scores --
resource_matrices <- setNames(lapply(TARGET_COHORTS, function(dataset) {
  validate_expression(readRDS(file.path(matrix_dir, paste0(dataset, ".rds"))),
                      paste0("normalized expression for ", dataset))
}), TARGET_COHORTS)

training_overlap <- rbindlist(lapply(TARGET_COHORTS, function(dataset) {
  overlap <- intersect(discovery_metadata$sample_id,
                       colnames(resource_matrices[[dataset]]))
  data.table(
    dataset = dataset,
    evaluation_role = evaluation_role(dataset),
    n_target_participants = ncol(resource_matrices[[dataset]]),
    n_exact_sample_id_overlap_with_discovery = length(overlap),
    overlapping_sample_ids = if (length(overlap)) paste(overlap, collapse = ";") else NA_character_
  )
}))
assert_true(all(training_overlap[
  evaluation_role == "independent_validation",
  n_exact_sample_id_overlap_with_discovery
] == 0L), "An independent validation participant occurs in the RF discovery training set")

unsupervised_scores <- fread(unsupervised_scores_path)
required_columns(unsupervised_scores,
                 c("sample_id", "dataset", "axis_id", "axis_raw", "axis_percentile"),
                 "unsupervised/participant_scores.tsv")
cohort_pc1 <- unsupervised_scores[axis_id == "signature_pc1"]
assert_true(nrow(cohort_pc1) > 0L,
            "No signature_pc1 rows exist in unsupervised participant scores")
assert_true(!anyDuplicated(cohort_pc1[, .(sample_id, dataset)]),
            "signature_pc1 has duplicate participant rows")
assert_true(setequal(unique(cohort_pc1$dataset), TARGET_COHORTS),
            "signature_pc1 does not cover exactly the five Resource cohorts")

coverage <- rbindlist(lapply(TARGET_COHORTS, function(dataset) {
  mat <- resource_matrices[[dataset]]
  present <- intersect(signature_all, rownames(mat))
  data.table(
    dataset = dataset,
    evaluation_role = evaluation_role(dataset),
    n_participants = ncol(mat),
    n_signature_total = SIGNATURE_SIZE,
    n_signature_present = length(present),
    coverage_fraction = length(present) / SIGNATURE_SIZE,
    n_frozen_loading_present = sum(common_model_genes %in% rownames(mat))
  )
}))
assert_true(all(coverage$coverage_fraction >= MIN_SIGNATURE_COVERAGE),
            "At least one Resource cohort has less than 90% signature coverage")
assert_true(all(coverage$n_frozen_loading_present >= ceiling(MIN_SIGNATURE_COVERAGE * SIGNATURE_SIZE)),
            "At least one Resource cohort has too few frozen-loading genes")

fixed_scores <- rbindlist(lapply(TARGET_COHORTS, function(dataset) {
  mat <- resource_matrices[[dataset]]
  genes <- common_model_genes[common_model_genes %in% rownames(mat)]
  projected <- as.numeric(crossprod(
    loading[genes],
    sweep(mat[genes, , drop = FALSE], 1L, center[genes], FUN = "-")
  ))
  data.table(
    sample_id = colnames(mat),
    dataset = dataset,
    fixed_projection_raw = projected,
    fixed_projection_percentile = rank_percentile(projected),
    n_projection_genes = length(genes),
    projection_coverage_fraction = length(genes) / SIGNATURE_SIZE
  )
}))

for (cohort_id in TARGET_COHORTS) {
  matrix_samples <- colnames(resource_matrices[[cohort_id]])
  pc1_samples <- cohort_pc1[cohort_pc1$dataset == cohort_id, sample_id]
  projection_samples <- fixed_scores[fixed_scores$dataset == cohort_id, sample_id]
  assert_true(setequal(matrix_samples, pc1_samples),
              paste0(cohort_id, ": matrix and cohort-PC1 participant sets differ"))
  assert_true(setequal(matrix_samples, projection_samples),
              paste0(cohort_id, ": matrix and projection participant sets differ"))
}

pc1_for_join <- cohort_pc1[, .(
  sample_id, dataset,
  signature_pc1_raw = axis_raw,
  signature_pc1_percentile = axis_percentile
)]
score_registry <- merge(pc1_for_join, fixed_scores,
                        by = c("sample_id", "dataset"), all = TRUE, sort = FALSE)
assert_true(nrow(score_registry) == sum(vapply(resource_matrices, ncol, integer(1))),
            "Score join changed the Resource participant count")
assert_true(all(complete.cases(score_registry)),
            "Cohort-PC1/projection score registry contains missing values")
score_registry[, consensus_mean_percentile :=
                 (signature_pc1_percentile + fixed_projection_percentile) / 2]

# ------------------------------------------------ discovery-only RF benchmark --
# The RF predictor family is frozen from discovery alone. A target cohort cannot
# change which variables enter the fitted forest. If a frozen gene failed a
# target cohort's row-sum filter, its centered target value is set to zero (the
# frozen discovery mean) and the event remains explicit in coverage.tsv.
rf_genes <- common_model_genes
assert_true(all(rf_genes %in% rownames(discovery_expression)) &&
              all(rf_genes %in% names(center)),
            "A discovery-frozen RF predictor lacks expression or a frozen center")
assert_true(MTRY <= length(rf_genes), "RF mtry exceeds the available predictor count")

x_train <- t(sweep(discovery_expression[rf_genes, , drop = FALSE],
                   1L, center[rf_genes], FUN = "-"))
colnames(x_train) <- rf_genes
assert_true(all(is.finite(x_train)), "RF training matrix contains non-finite values")

set.seed(SEED)
rf_fit <- ranger(
  x = x_train,
  y = as.numeric(discovery_metadata$fibrosis),
  num.trees = N_TREES,
  mtry = MTRY,
  min.node.size = MIN_NODE_SIZE,
  seed = SEED,
  num.threads = 1L,
  write.forest = TRUE,
  importance = "none",
  verbose = FALSE
)

rf_scores <- rbindlist(lapply(TARGET_COHORTS, function(dataset) {
  mat <- resource_matrices[[dataset]]
  present <- intersect(rf_genes, rownames(mat))
  x_target <- matrix(
    0, nrow = ncol(mat), ncol = length(rf_genes),
    dimnames = list(colnames(mat), rf_genes)
  )
  x_target[, present] <- t(sweep(mat[present, , drop = FALSE],
                                 1L, center[present], FUN = "-"))
  prediction <- as.numeric(predict(rf_fit, data = x_target,
                                   num.threads = 1L)$predictions)
  assert_true(length(prediction) == ncol(mat) && all(is.finite(prediction)),
              paste0("RF prediction failed for ", dataset))
  data.table(
    sample_id = colnames(mat),
    dataset = dataset,
    rf_fibrosis_prediction = prediction,
    rf_fibrosis_percentile = rank_percentile(prediction),
    n_rf_genes_observed = length(present),
    n_rf_genes_imputed_to_discovery_center = length(rf_genes) - length(present)
  )
}))
rf_target_coverage <- unique(rf_scores[, .(
  dataset, n_rf_genes_observed, n_rf_genes_imputed_to_discovery_center
)])
assert_true(nrow(rf_target_coverage) == length(TARGET_COHORTS),
            "RF target-predictor coverage audit is incomplete")
coverage <- merge(coverage, rf_target_coverage, by = "dataset", all.x = TRUE,
                  sort = FALSE)
assert_true(!anyNA(coverage$n_rf_genes_observed),
            "RF target-predictor coverage join failed")

score_registry <- merge(score_registry, rf_scores,
                        by = c("sample_id", "dataset"), all = TRUE, sort = FALSE)
assert_true(all(complete.cases(score_registry)), "RF score join introduced missing values")
score_registry[, evaluation_role := evaluation_role(dataset)]

# Annotation is deliberately merged only after all scores are final. It is never
# passed to PCA projection, consensus construction, or Resource prediction.
resource_metadata <- fread(manifest_path)
required_columns(resource_metadata,
                 c("sample_id", "dataset", "inferred_sex", "fibrosis_stage", "nas_score"),
                 "five-cohort Resource manifest")
assert_true(!anyDuplicated(resource_metadata$sample_id),
            "Resource manifest has duplicate biological participants")
assert_true(setequal(resource_metadata$dataset, TARGET_COHORTS),
            "Resource manifest does not contain exactly the five frozen cohorts")
annotation <- resource_metadata[, .(sample_id, dataset,
                                    inferred_sex, fibrosis_stage, nas_score)]
assert_true(setequal(
  paste(annotation$dataset, annotation$sample_id, sep = "::"),
  paste(score_registry$dataset, score_registry$sample_id, sep = "::")
), "Resource manifest and frozen score registry participant sets differ")
score_registry <- merge(score_registry, annotation,
                        by = c("sample_id", "dataset"), all.x = TRUE, sort = FALSE)
assert_true(!anyNA(score_registry$inferred_sex),
            "Resource sex annotation is missing after score construction")
setorder(score_registry, dataset, sample_id)

long_scores <- rbindlist(list(
  score_registry[, .(
    sample_id, dataset, evaluation_role,
    axis_id = "fixed_projection",
    axis_role = "co_primary",
    axis_raw = fixed_projection_raw,
    axis_percentile = fixed_projection_percentile,
    n_genes = n_projection_genes,
    coverage_fraction = projection_coverage_fraction
  )],
  score_registry[, .(
    sample_id, dataset, evaluation_role,
    axis_id = "consensus_rank",
    axis_role = "presentation_only",
    axis_raw = consensus_mean_percentile,
    axis_percentile = consensus_mean_percentile,
    n_genes = n_projection_genes,
    coverage_fraction = projection_coverage_fraction
  )],
  score_registry[, .(
    sample_id, dataset, evaluation_role,
    axis_id = "rf_fibrosis",
    axis_role = "histology_anchored_secondary",
    axis_raw = rf_fibrosis_prediction,
    axis_percentile = rf_fibrosis_percentile,
    n_genes = length(rf_genes),
    coverage_fraction = length(rf_genes) / SIGNATURE_SIZE
  )]
), use.names = TRUE)
long_scores <- merge(
  long_scores,
  annotation,
  by = c("sample_id", "dataset"), all.x = TRUE, sort = FALSE
)
long_scores[, `:=`(
  target_histology_used_for_score = FALSE,
  seed = fifelse(axis_id == "rf_fibrosis", SEED, NA_integer_)
)]
setorder(long_scores, axis_id, dataset, sample_id)

axis_vectors <- c(
  signature_pc1 = "signature_pc1_percentile",
  fixed_projection = "fixed_projection_percentile",
  consensus_rank = "consensus_mean_percentile",
  rf_fibrosis = "rf_fibrosis_percentile"
)
agreement_pairs <- combn(names(axis_vectors), 2L, simplify = FALSE)
score_agreement <- rbindlist(lapply(TARGET_COHORTS, function(cohort_id) {
  d <- score_registry[score_registry$dataset == cohort_id]
  rbindlist(lapply(agreement_pairs, function(pair) {
    x_column <- axis_vectors[[pair[[1]]]]
    y_column <- axis_vectors[[pair[[2]]]]
    data.table(
      dataset = cohort_id,
      evaluation_role = evaluation_role(cohort_id),
      axis_x = pair[[1]],
      axis_y = pair[[2]],
      n = nrow(d),
      spearman_rho = cor(d[[x_column]], d[[y_column]], method = "spearman"),
      pearson_r = cor(d[[x_column]], d[[y_column]], method = "pearson")
    )
  }))
}))
co_primary_agreement <- score_agreement[
  axis_x == "signature_pc1" & axis_y == "fixed_projection",
  .(dataset, evaluation_role, n, spearman_rho, pearson_r,
    threshold = 0.70,
    pass = is.finite(spearman_rho) & spearman_rho >= 0.70)
]

# -------------------------------------------------------------- provenance --
training_samples <- discovery_metadata[, .(
  sample_id,
  dataset,
  fibrosis,
  in_training = TRUE,
  target_resource_histology_used = FALSE
)]
setorder(training_samples, dataset, sample_id)

training_manifest <- data.table(
  model_id = "rf_fibrosis",
  engine = "ranger",
  ranger_version = as.character(packageVersion("ranger")),
  response = "numeric_fibrosis_F0_F4",
  training_scope = "published_discovery_135_only",
  n_training_donors = nrow(discovery_metadata),
  n_training_datasets = uniqueN(discovery_metadata$dataset),
  n_predictors = length(rf_genes),
  predictor_freeze_scope = "published_discovery_signature_only",
  missing_target_predictor_rule = "centered_value_zero_equals_discovery_center",
  num_trees = N_TREES,
  mtry = MTRY,
  min_node_size = MIN_NODE_SIZE,
  seed = SEED,
  num_threads = 1L,
  target_histology_used_for_preprocessing = FALSE,
  target_histology_used_for_training = FALSE,
  n_independent_sample_id_overlap = sum(training_overlap[
    evaluation_role == "independent_validation",
    n_exact_sample_id_overlap_with_discovery
  ]),
  oob_prediction_mse = rf_fit$prediction.error
)

projection_loadings <- merge(
  signature,
  data.table(
    gene_id_base = common_model_genes,
    discovery_center = center[common_model_genes],
    discovery_loading_oriented = loading[common_model_genes]
  ),
  by = "gene_id_base", all.x = TRUE, sort = FALSE
)
projection_loadings[, used_for_fixed_projection :=
                      gene_id_base %in% common_model_genes]
projection_loadings[, used_for_rf := gene_id_base %in% rf_genes]
setorder(projection_loadings, gene_id_base)

write_tsv_once(score_registry, file.path(output_dir, "score_registry.tsv"))
write_tsv_once(long_scores, file.path(output_dir, "participant_scores.tsv"))
write_tsv_once(coverage, file.path(output_dir, "coverage.tsv"))
write_tsv_once(score_agreement, file.path(output_dir, "score_agreement.tsv"))
write_tsv_once(co_primary_agreement,
               file.path(output_dir, "co_primary_agreement.tsv"))
write_tsv_once(training_manifest,
               file.path(output_dir, "rf_training_manifest.tsv"))
write_tsv_once(training_samples,
               file.path(output_dir, "rf_training_samples.tsv"))
write_tsv_once(training_overlap,
               file.path(output_dir, "rf_training_overlap_audit.tsv"))
write_tsv_once(data.table(gene_id_base = rf_genes),
               file.path(output_dir, "rf_training_genes.tsv"))
write_tsv_once(projection_loadings,
               file.path(output_dir, "fixed_projection_loadings.tsv"))

model_bundle <- list(
  model_id = "rf_fibrosis",
  ranger_fit = rf_fit,
  training_gene_ids = rf_genes,
  discovery_center = center[rf_genes],
  parameters = list(
    num.trees = N_TREES,
    mtry = MTRY,
    min.node.size = MIN_NODE_SIZE,
    seed = SEED,
    num.threads = 1L
  ),
  training_scope = "published_discovery_135_only",
  target_histology_used = FALSE
)
save_rds_once(model_bundle, file.path(output_dir, "rf_model.rds"))

run_parameters <- data.table(
  parameter = c("seed", "num_trees", "mtry", "min_node_size",
                "signature_size", "minimum_signature_coverage",
                "expected_resource_signature", "n_rf_genes",
                "target_histology_used"),
  value = c(as.character(SEED), as.character(N_TREES), as.character(MTRY),
            as.character(MIN_NODE_SIZE), as.character(SIGNATURE_SIZE),
            as.character(MIN_SIGNATURE_COVERAGE),
            as.character(EXPECTED_RESOURCE_SIGNATURE),
            as.character(length(rf_genes)), "FALSE")
)
write_tsv_once(run_parameters, file.path(output_dir, "run_parameters.tsv"))
write_lines_once(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

output_files <- list.files(output_dir, full.names = TRUE)
output_files <- output_files[basename(output_files) != "output_checksums.tsv"]
output_manifest <- data.table(
  relative_path = basename(output_files),
  size_bytes = as.numeric(file.size(output_files)),
  sha256 = vapply(output_files, sha256_file, character(1))
)
setorder(output_manifest, relative_path)
write_tsv_once(output_manifest, file.path(output_dir, "output_checksums.tsv"))

cat(sprintf(
  "PROJECTION_RF_COMPLETE: %d Resource participants; %d RF genes; output=%s\n",
  nrow(score_registry), length(rf_genes), output_dir
))
