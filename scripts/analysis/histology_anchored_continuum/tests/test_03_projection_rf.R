#!/usr/bin/env Rscript

# Small end-to-end contract test for the frozen projection, discovery-only RF,
# and percentile consensus. The fixture preserves the production cardinalities
# that 03 treats as provenance invariants (145/139 genes and 135 training donors).

suppressPackageStartupMessages(library(data.table))

project_root <- Sys.getenv("MASLD_PROJECT_ROOT", getwd())
script_path <- file.path(
  project_root,
  "scripts/analysis/histology_anchored_continuum/03_build_projection_rf_axes.R"
)
stopifnot(file.exists(script_path))

TARGET_COHORTS <- c("GSE126848", "GSE130970", "GSE135251",
                    "GSE162694", "GSE213621")

make_fixture <- function(root) {
  dir.create(file.path(root, "reproduction"), recursive = TRUE)
  dir.create(file.path(root, "unsupervised", "normalized_expression_by_cohort"),
             recursive = TRUE)

  set.seed(20260817L)
  signature <- sprintf("ENSG%011d", seq_len(145L))
  common <- signature[seq_len(139L)]
  loading <- setNames(rnorm(length(common)), common)
  loading <- loading / sqrt(sum(loading^2))
  center <- setNames(runif(length(common), 5, 10), common)
  saveRDS(list(
    signature_gene_ids_all = signature,
    common_gene_ids = common,
    discovery_center = center,
    discovery_loading_oriented = loading
  ), file.path(root, "reproduction", "discovery_signature_model.rds"))

  n_discovery <- 135L
  discovery_samples <- sprintf("D%03d", seq_len(n_discovery))
  fibrosis <- rep(0:4, length.out = n_discovery)
  discovery_expression <- matrix(
    rnorm(length(common) * n_discovery, sd = 0.5),
    nrow = length(common), ncol = n_discovery,
    dimnames = list(common, discovery_samples)
  ) + center + outer(loading, fibrosis * 2)
  saveRDS(discovery_expression,
          file.path(root, "reproduction", "discovery_normalized_expression.rds"))
  fwrite(data.table(
    sample_id = discovery_samples,
    dataset = rep(c("UCAM", "SANYAL"), length.out = n_discovery),
    sex = rep(c("F", "M"), length.out = n_discovery),
    NAS = rep(0:8, length.out = n_discovery),
    fibrosis = fibrosis
  ), file.path(root, "reproduction", "discovery_metadata.tsv"), sep = "\t")
  fwrite(data.table(
    gene_id_base = signature,
    gene_symbol = paste0("G", seq_along(signature)),
    in_resource = signature %in% common
  ), file.path(root, "reproduction", "signature_genes.tsv"), sep = "\t")

  score_rows <- list()
  metadata_rows <- list()
  for (i in seq_along(TARGET_COHORTS)) {
    dataset <- TARGET_COHORTS[[i]]
    samples <- paste0(dataset, "_", sprintf("%02d", seq_len(12L)))
    latent <- seq(-2, 2, length.out = length(samples)) + rnorm(length(samples), sd = 0.1)
    expression <- matrix(
      rnorm(length(common) * length(samples), sd = 0.6),
      nrow = length(common), ncol = length(samples),
      dimnames = list(common, samples)
    ) + center + outer(loading, latent * 3)
    saveRDS(expression, file.path(
      root, "unsupervised", "normalized_expression_by_cohort", paste0(dataset, ".rds")
    ))
    raw <- as.numeric(crossprod(loading, sweep(expression, 1L, center, "-"))) +
      rnorm(length(samples), sd = 0.1)
    score_rows[[i]] <- data.table(
      sample_id = samples,
      dataset = dataset,
      axis_id = "signature_pc1",
      axis_raw = raw,
      axis_percentile = (rank(raw) - 1) / (length(raw) - 1)
    )
    metadata_rows[[i]] <- data.table(
      sample_id = samples,
      dataset = dataset,
      inferred_sex = rep(c("F", "M"), length.out = length(samples)),
      fibrosis_stage = rep(0:4, length.out = length(samples)),
      nas_score = rep(0:8, length.out = length(samples))
    )
  }
  fwrite(rbindlist(score_rows),
         file.path(root, "unsupervised", "participant_scores.tsv"), sep = "\t")
  manifest_path <- file.path(root, "manifest.tsv")
  fwrite(rbindlist(metadata_rows), manifest_path, sep = "\t")
  manifest_path
}

run_fixture <- function(root, manifest_path) {
  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script_path)),
    env = c(
      paste0("MASLD_PROJECT_ROOT=", shQuote(project_root)),
      paste0("HAC_OUT_ROOT=", shQuote(root)),
      paste0("HAC_MANIFEST_PATH=", shQuote(manifest_path))
    ),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop(paste(output, collapse = "\n"), call. = FALSE)
  }
  invisible(output)
}

root_one <- tempfile("hac03_test_one_")
root_two <- tempfile("hac03_test_two_")
root_labels <- tempfile("hac03_test_labels_")
root_overlap <- tempfile("hac03_test_overlap_")
dir.create(root_one)
dir.create(root_two)
dir.create(root_labels)
dir.create(root_overlap)
on.exit(unlink(c(root_one, root_two, root_labels, root_overlap),
               recursive = TRUE, force = TRUE), add = TRUE)
manifest_one <- make_fixture(root_one)
manifest_two <- make_fixture(root_two)
manifest_labels <- make_fixture(root_labels)
run_fixture(root_one, manifest_one)
run_fixture(root_two, manifest_two)

# Target-cohort histology, NAS, and sex are annotation only. Permuting them must
# not change any PCA, projection, consensus, or RF score.
target_labels <- fread(manifest_labels)
set.seed(9001L)
target_labels[, `:=`(
  inferred_sex = sample(inferred_sex),
  fibrosis_stage = sample(fibrosis_stage),
  nas_score = sample(nas_score)
), by = dataset]
fwrite(target_labels, manifest_labels, sep = "\t")
run_fixture(root_labels, manifest_labels)

scores_one <- fread(file.path(root_one, "projection", "participant_scores.tsv"))
scores_two <- fread(file.path(root_two, "projection", "participant_scores.tsv"))
registry <- fread(file.path(root_one, "projection", "score_registry.tsv"))
agreement <- fread(file.path(root_one, "projection", "score_agreement.tsv"))
training <- fread(file.path(root_one, "projection", "rf_training_manifest.tsv"))
training_samples <- fread(file.path(root_one, "projection", "rf_training_samples.tsv"))
overlap <- fread(file.path(root_one, "projection", "rf_training_overlap_audit.tsv"))
coverage <- fread(file.path(root_one, "projection", "coverage.tsv"))
registry_labels <- fread(file.path(root_labels, "projection", "score_registry.tsv"))

stopifnot(nrow(scores_one) == 5L * 12L * 3L)
stopifnot(setequal(unique(scores_one$axis_id),
                   c("fixed_projection", "consensus_rank", "rf_fibrosis")))
stopifnot(all(scores_one$axis_percentile >= 0 & scores_one$axis_percentile <= 1))
stopifnot(isTRUE(all.equal(
  registry$consensus_mean_percentile,
  (registry$signature_pc1_percentile + registry$fixed_projection_percentile) / 2,
  tolerance = 1e-14
)))

# Re-derive one fixed projection from the frozen center and loading rather than
# trusting the script's own intermediate output.
model <- readRDS(file.path(root_one, "reproduction", "discovery_signature_model.rds"))
first_dataset <- TARGET_COHORTS[[1L]]
expression <- readRDS(file.path(
  root_one, "unsupervised", "normalized_expression_by_cohort",
  paste0(first_dataset, ".rds")
))
expected_projection <- as.numeric(crossprod(
  model$discovery_loading_oriented,
  sweep(expression, 1L, model$discovery_center, "-")
))
observed_projection <- registry[dataset == first_dataset][
  match(colnames(expression), sample_id), fixed_projection_raw
]
stopifnot(isTRUE(all.equal(observed_projection, expected_projection, tolerance = 1e-12)))

stopifnot(all(c("dataset", "axis_x", "axis_y", "spearman_rho", "n") %in%
                names(agreement)))
stopifnot(training$n_training_donors == 135L)
stopifnot(training$num_trees == 1000L && training$mtry == 11L &&
            training$min_node_size == 5L && training$seed == 20260817L)
stopifnot(!training$target_histology_used_for_preprocessing)
stopifnot(!training$target_histology_used_for_training)
stopifnot(nrow(training_samples) == 135L && all(training_samples$in_training))
stopifnot(all(overlap[evaluation_role == "independent_validation",
                      n_exact_sample_id_overlap_with_discovery] == 0L))
stopifnot(all(coverage$n_signature_present == 139L))

score_columns <- c(
  "sample_id", "dataset", "signature_pc1_raw", "signature_pc1_percentile",
  "fixed_projection_raw", "fixed_projection_percentile",
  "consensus_mean_percentile", "rf_fibrosis_prediction",
  "rf_fibrosis_percentile"
)
registry_scores <- registry[, ..score_columns]
label_permuted_scores <- registry_labels[, ..score_columns]
setorder(registry_scores, dataset, sample_id)
setorder(label_permuted_scores, dataset, sample_id)
stopifnot(identical(registry_scores, label_permuted_scores))

# One thread, a fixed seed, and identical inputs must produce identical scores.
setcolorder(scores_two, names(scores_one))
stopifnot(identical(scores_one, scores_two))

# Write-once protection must reject reuse of a completed namespace.
rerun <- system2(
  file.path(R.home("bin"), "Rscript"),
  c("--vanilla", shQuote(script_path)),
  env = c(
    paste0("MASLD_PROJECT_ROOT=", shQuote(project_root)),
    paste0("HAC_OUT_ROOT=", shQuote(root_one)),
    paste0("HAC_MANIFEST_PATH=", shQuote(manifest_one))
  ),
  stdout = TRUE, stderr = TRUE
)
stopifnot(!is.null(attr(rerun, "status")) && attr(rerun, "status") != 0L)

# The integration must fail before model fitting if an independent-validation
# participant ID occurs in the 135-donor discovery training set.
manifest_overlap <- make_fixture(root_overlap)
overlap_dataset <- "GSE162694"
matrix_path <- file.path(
  root_overlap, "unsupervised", "normalized_expression_by_cohort",
  paste0(overlap_dataset, ".rds")
)
overlap_matrix <- readRDS(matrix_path)
old_sample_id <- colnames(overlap_matrix)[[1L]]
colnames(overlap_matrix)[[1L]] <- "D001"
saveRDS(overlap_matrix, matrix_path)

overlap_scores <- fread(file.path(
  root_overlap, "unsupervised", "participant_scores.tsv"
))
overlap_scores[dataset == overlap_dataset & sample_id == old_sample_id,
               sample_id := "D001"]
fwrite(overlap_scores, file.path(
  root_overlap, "unsupervised", "participant_scores.tsv"
), sep = "\t")
overlap_manifest <- fread(manifest_overlap)
overlap_manifest[dataset == overlap_dataset & sample_id == old_sample_id,
                 sample_id := "D001"]
fwrite(overlap_manifest, manifest_overlap, sep = "\t")

overlap_run <- system2(
  file.path(R.home("bin"), "Rscript"),
  c("--vanilla", shQuote(script_path)),
  env = c(
    paste0("MASLD_PROJECT_ROOT=", shQuote(project_root)),
    paste0("HAC_OUT_ROOT=", shQuote(root_overlap)),
    paste0("HAC_MANIFEST_PATH=", shQuote(manifest_overlap))
  ),
  stdout = TRUE, stderr = TRUE
)
stopifnot(!is.null(attr(overlap_run, "status")) &&
            attr(overlap_run, "status") != 0L)
stopifnot(any(grepl(
  "An independent validation participant occurs in the RF discovery training set",
  overlap_run, fixed = TRUE
)))

cat("test_03_projection_rf: PASS\n")
