#!/usr/bin/env Rscript
# Build label-blind, cohort-specific molecular-continuum axes for the five-
# cohort Resource substrate.
#
# This intentionally does not reproduce the discovery-cohort ComBat step. Each
# Resource cohort is processed independently, so neither another cohort nor a
# target histology label can affect a participant's expression transform. The
# expression filter and transform follow the released Kamzolas code:
#
#   rowSums(fragment_counts) > number_of_participants
#   log2(fragment_counts + 1)
#   quantile normalisation
#
# Two axes are emitted. `signature_pc1` is co-primary; it uses only the
# frozen 145-gene signature genes that were mapped by the reproduction stage.
# `full_transcriptome_pc1` is an unsupervised control. PCA signs are arbitrary,
# so both signs are fixed solely by cosine alignment to the frozen discovery
# loadings. Histology is never read by this script.

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

options(digits = 17, scipen = 999)

HAC_SEED <- 20260817L
EXPECTED_COHORTS <- c(
  "GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"
)
EXPECTED_COHORT_SIZES <- c(
  GSE126848 = 53L, GSE130970 = 76L, GSE135251 = 214L,
  GSE162694 = 140L, GSE213621 = 361L
)
SIGNATURE_TOTAL_EXPECTED <- 145L
SIGNATURE_RESOURCE_EXPECTED <- 139L
MIN_SIGNATURE_COVERAGE <- 0.90

fail <- function(...) stop(paste0(...), call. = FALSE)

assert_true <- function(value, ...) {
  if (!isTRUE(value)) fail(...)
  invisible(TRUE)
}

required_env <- function(name) {
  value <- Sys.getenv(name, "")
  assert_true(nzchar(value), name, " is unset")
  value
}

strip_ensembl_version <- function(x) sub("\\..*$", "", as.character(x))

sha256_file <- function(path) {
  assert_true(file.exists(path), "Cannot hash absent input: ", path)
  output <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  assert_true(is.null(status) || identical(status, 0L),
              "sha256sum failed for ", path, ": ", paste(output, collapse = " "))
  sub("[[:space:]].*$", "", output[[1L]])
}

write_tsv_once <- function(value, path) {
  assert_true(!file.exists(path), "Refusing to overwrite output: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(value, path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

save_rds_once <- function(value, path) {
  assert_true(!file.exists(path), "Refusing to overwrite output: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(value, path, version = 3L)
  invisible(path)
}

percentile_rank <- function(x) {
  assert_true(all(is.finite(x)), "Cannot rank a non-finite axis")
  if (length(x) == 1L) return(0.5)
  (rank(x, ties.method = "average") - 1) / (length(x) - 1)
}

as_frozen_vector <- function(value, field, fallback_names = NULL) {
  original_names <- names(value)
  value <- as.numeric(value)
  # as.numeric() drops names, so capture them before coercion.
  if (!is.null(original_names)) names(value) <- strip_ensembl_version(original_names)
  if (is.null(names(value)) && !is.null(fallback_names)) {
    assert_true(length(value) == length(fallback_names),
                field, " has no names and does not match its declared gene list")
    names(value) <- strip_ensembl_version(fallback_names)
  }
  assert_true(!is.null(names(value)) && all(nzchar(names(value))),
              field, " must be a gene-named numeric vector")
  assert_true(!anyDuplicated(names(value)), field, " has duplicated Ensembl base IDs")
  assert_true(all(is.finite(value)), field, " contains non-finite values")
  value
}

cosine_similarity <- function(x, y) {
  assert_true(length(x) == length(y) && length(x) >= 2L,
              "Cosine inputs have incompatible dimensions")
  denom <- sqrt(sum(x * x)) * sqrt(sum(y * y))
  assert_true(is.finite(denom) && denom > 0, "Cannot align against a zero loading vector")
  sum(x * y) / denom
}

fit_oriented_pc1 <- function(expression, frozen_loading, axis_id, cohort_id) {
  assert_true(is.matrix(expression), "PCA expression input must be a matrix")
  assert_true(nrow(expression) >= 2L && ncol(expression) >= 2L,
              cohort_id, " / ", axis_id, " has insufficient dimensions")
  assert_true(!anyNA(expression) && all(is.finite(expression)),
              cohort_id, " / ", axis_id, " expression contains non-finite values")
  assert_true(!anyDuplicated(rownames(expression)),
              cohort_id, " / ", axis_id, " has duplicated gene IDs")
  assert_true(!anyDuplicated(colnames(expression)),
              cohort_id, " / ", axis_id, " has duplicated sample IDs")

  # rank.=2 avoids materialising the full right-singular-vector matrix while
  # retaining the first two singular values for a degeneracy diagnostic.
  fit <- prcomp(t(expression), center = TRUE, scale. = FALSE, rank. = 2L)
  assert_true(ncol(fit$x) >= 1L && ncol(fit$rotation) >= 1L,
              cohort_id, " / ", axis_id, " PCA did not return PC1")
  raw_score <- fit$x[, 1L]
  raw_loading <- fit$rotation[, 1L]
  names(raw_loading) <- rownames(fit$rotation)

  shared <- intersect(names(raw_loading), names(frozen_loading))
  min_alignment_genes <- if (identical(axis_id, "signature_pc1")) {
    ceiling(MIN_SIGNATURE_COVERAGE * SIGNATURE_TOTAL_EXPECTED)
  } else {
    1000L
  }
  assert_true(length(shared) >= min_alignment_genes,
              cohort_id, " / ", axis_id, " has only ", length(shared),
              " genes for discovery-loading alignment; requires ", min_alignment_genes)

  raw_cosine <- cosine_similarity(raw_loading[shared], frozen_loading[shared])
  assert_true(is.finite(raw_cosine) && abs(raw_cosine) > sqrt(.Machine$double.eps),
              cohort_id, " / ", axis_id,
              " PC1 is orthogonal to the frozen discovery loading; orientation is undefined")
  orientation_sign <- if (raw_cosine < 0) -1 else 1
  oriented_score <- as.numeric(raw_score) * orientation_sign
  names(oriented_score) <- rownames(fit$x)
  oriented_loading <- raw_loading * orientation_sign
  oriented_cosine <- cosine_similarity(oriented_loading[shared], frozen_loading[shared])
  assert_true(oriented_cosine > 0, cohort_id, " / ", axis_id,
              " failed to align positively with the discovery loading")

  variance_fraction <- fit$sdev^2 / sum(fit$sdev^2)
  assert_true(length(variance_fraction) >= 1L && is.finite(variance_fraction[[1L]]) &&
                variance_fraction[[1L]] > 0,
              cohort_id, " / ", axis_id, " has zero or invalid PC1 variance")

  list(
    score = oriented_score,
    raw_loading = raw_loading,
    oriented_loading = oriented_loading,
    raw_cosine = raw_cosine,
    oriented_cosine = oriented_cosine,
    orientation_sign = orientation_sign,
    n_alignment_genes = length(shared),
    pc1_variance_pct = 100 * variance_fraction[[1L]],
    pc1_pc2_sdev_ratio = if (length(fit$sdev) >= 2L && fit$sdev[[2L]] > 0) {
      fit$sdev[[1L]] / fit$sdev[[2L]]
    } else {
      Inf
    }
  )
}

normalise_cohort <- function(counts, cohort_id) {
  assert_true(is.matrix(counts), cohort_id, " counts are not a matrix")
  assert_true(ncol(counts) >= 2L, cohort_id, " has fewer than two participants")
  assert_true(!anyNA(counts) && all(is.finite(counts)),
              cohort_id, " counts contain non-finite values")
  assert_true(all(counts >= 0) && all(counts == floor(counts)),
              cohort_id, " input is not a raw non-negative integer fragment matrix")

  # Exact released-code threshold: keep genes whose cohort-wide sum exceeds the
  # number of samples (mean raw fragment count > 1).
  keep <- rowSums(counts) > ncol(counts)
  assert_true(sum(keep) >= 1000L, cohort_id,
              " retains fewer than 1,000 genes after the paper's row-sum filter")
  log_expression <- log2(counts[keep, , drop = FALSE] + 1)

  # limma's pure-R implementation is numerically equivalent to
  # preprocessCore::normalize.quantiles for this complete matrix and does not
  # attempt to spawn threads outside a SLURM cgroup.
  normalised <- normalizeQuantiles(log_expression)
  dimnames(normalised) <- dimnames(log_expression)
  assert_true(identical(dim(normalised), dim(log_expression)),
              cohort_id, " quantile normalisation changed matrix dimensions")
  assert_true(all(is.finite(normalised)),
              cohort_id, " quantile normalisation returned non-finite values")
  normalised
}

main <- function() {
  out_root <- required_env("HAC_OUT_ROOT")
  dge_path <- required_env("HAC_DGE_PATH")
  manifest_path <- required_env("HAC_MANIFEST_PATH")
  reproduction_dir <- file.path(out_root, "reproduction")
  model_path <- file.path(reproduction_dir, "discovery_signature_model.rds")
  signature_path <- file.path(reproduction_dir, "signature_genes.tsv")
  output_dir <- file.path(out_root, "unsupervised")
  matrix_dir <- file.path(output_dir, "normalized_expression_by_cohort")

  for (path in c(dge_path, manifest_path, model_path, signature_path)) {
    assert_true(file.exists(path), "Required input is absent: ", path)
  }
  assert_true(!file.exists(output_dir),
              "Refusing to reuse an existing unsupervised output directory: ", output_dir)
  dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
  assert_true(dir.exists(matrix_dir), "Could not create output directory: ", matrix_dir)
  set.seed(HAC_SEED)

  # Select only identifiers. Histology columns in the shared manifest never
  # enter this process, including for axis orientation.
  meta <- fread(manifest_path, select = c("sample_id", "dataset"))
  required_meta <- c("sample_id", "dataset")
  assert_true(all(required_meta %in% names(meta)),
              "Manifest lacks required columns: ",
              paste(setdiff(required_meta, names(meta)), collapse = ", "))
  assert_true(nrow(meta) > 0L && !anyNA(meta$sample_id) && !anyNA(meta$dataset),
              "Manifest has missing participant or cohort identifiers")
  assert_true(!anyDuplicated(meta$sample_id), "Manifest sample_id is not unique")
  assert_true(setequal(unique(meta$dataset), EXPECTED_COHORTS),
              "Manifest cohorts differ from the frozen five-cohort Resource set")
  observed_sizes <- table(meta$dataset)
  assert_true(nrow(meta) == sum(EXPECTED_COHORT_SIZES) &&
                all(observed_sizes[names(EXPECTED_COHORT_SIZES)] == EXPECTED_COHORT_SIZES),
              "Manifest participant counts differ from the frozen 844-participant substrate")

  dge <- readRDS(dge_path)
  counts <- if (inherits(dge, "DGEList")) dge$counts else if (is.matrix(dge)) dge else NULL
  assert_true(is.matrix(counts), "HAC_DGE_PATH must contain a DGEList or count matrix")
  assert_true(!is.null(rownames(counts)) && !is.null(colnames(counts)),
              "Count matrix requires gene and sample dimnames")
  assert_true(setequal(meta$sample_id, colnames(counts)),
              "DGE samples do not exactly match the manifest participants")
  counts <- counts[, meta$sample_id, drop = FALSE]
  assert_true(identical(colnames(counts), meta$sample_id),
              "Count columns did not reorder to manifest order")
  base_ids <- strip_ensembl_version(rownames(counts))
  assert_true(all(grepl("^ENSG[0-9]+$", base_ids)),
              "DGE row names are not Ensembl gene IDs")
  assert_true(!anyDuplicated(base_ids),
              "Stripping Ensembl versions produced duplicate gene IDs")
  rownames(counts) <- base_ids

  model <- readRDS(model_path)
  required_model <- c(
    "signature_gene_ids_all", "common_gene_ids", "discovery_center",
    "discovery_loading_oriented", "discovery_full_center",
    "discovery_full_loading_oriented"
  )
  assert_true(is.list(model) && all(required_model %in% names(model)),
              "Discovery model lacks required fields: ",
              paste(setdiff(required_model, names(model)), collapse = ", "))

  signature <- fread(signature_path)
  assert_true(all(c("gene_id_base", "gene_symbol", "in_resource") %in% names(signature)),
              "signature_genes.tsv has an unexpected schema")
  signature[, gene_id_base := strip_ensembl_version(gene_id_base)]
  assert_true(!anyDuplicated(signature$gene_id_base),
              "signature_genes.tsv has duplicate Ensembl base IDs")

  signature_all <- unique(strip_ensembl_version(model$signature_gene_ids_all))
  common_genes <- unique(strip_ensembl_version(model$common_gene_ids))
  assert_true(length(signature_all) == SIGNATURE_TOTAL_EXPECTED,
              "Frozen signature contains ", length(signature_all), " rather than 145 genes")
  assert_true(length(common_genes) == SIGNATURE_RESOURCE_EXPECTED,
              "Frozen Resource intersection contains ", length(common_genes),
              " rather than 139 genes")
  assert_true(setequal(signature_all, signature$gene_id_base),
              "Signature table and discovery model contain different 145-gene sets")
  assert_true(all(common_genes %in% signature_all),
              "common_gene_ids is not a subset of the frozen signature")
  in_resource <- signature$gene_id_base[as.logical(signature$in_resource)]
  assert_true(setequal(common_genes, in_resource),
              "signature_genes.tsv in_resource flags disagree with common_gene_ids")
  assert_true(setequal(common_genes, intersect(signature_all, rownames(counts))),
              "Frozen common_gene_ids disagrees with the current DGE substrate")

  signature_loading <- as_frozen_vector(
    model$discovery_loading_oriented,
    "discovery_loading_oriented",
    fallback_names = common_genes
  )
  full_loading <- as_frozen_vector(
    model$discovery_full_loading_oriented,
    "discovery_full_loading_oriented"
  )
  assert_true(all(common_genes %in% names(signature_loading)),
              "Frozen signature loading omits common signature genes")
  assert_true(length(intersect(rownames(counts), names(full_loading))) >= 1000L,
              "Frozen full-transcriptome loading overlaps fewer than 1,000 Resource genes")

  input_manifest <- data.table(
    input_role = c("raw_fragment_dge", "participant_manifest",
                   "frozen_discovery_model", "frozen_signature_table"),
    path = c(dge_path, manifest_path, model_path, signature_path)
  )
  input_manifest[, `:=`(
    size_bytes = file.info(path)$size,
    sha256 = vapply(path, sha256_file, character(1L))
  )]
  write_tsv_once(input_manifest, file.path(output_dir, "input_checksums.tsv"))

  score_rows <- list()
  loading_rows <- list()
  coverage_rows <- list()
  gene_coverage_rows <- list()
  alignment_rows <- list()
  row_index <- 0L

  for (cohort_id in EXPECTED_COHORTS) {
    cohort_samples <- meta[dataset == cohort_id, sample_id]
    assert_true(length(cohort_samples) >= 2L, cohort_id, " has insufficient participants")
    cohort_counts <- counts[, cohort_samples, drop = FALSE]
    expression <- normalise_cohort(cohort_counts, cohort_id)
    matrix_path <- file.path(matrix_dir, paste0(cohort_id, ".rds"))
    save_rds_once(expression, matrix_path)

    in_raw <- signature_all %in% rownames(cohort_counts)
    passed_filter <- signature_all %in% rownames(expression)
    used_signature <- common_genes[common_genes %in% rownames(expression)]
    coverage_fraction <- length(used_signature) / length(signature_all)
    coverage_pass <- coverage_fraction >= MIN_SIGNATURE_COVERAGE
    assert_true(coverage_pass, cohort_id, " retains ", length(used_signature),
                "/145 signature genes after filtering (",
                sprintf("%.1f%%", 100 * coverage_fraction), "); requires at least 90%")

    coverage_rows[[length(coverage_rows) + 1L]] <- data.table(
      dataset = cohort_id,
      n_participants = length(cohort_samples),
      n_genes_raw = nrow(cohort_counts),
      n_genes_after_row_sum_filter = nrow(expression),
      n_signature_total = length(signature_all),
      n_signature_in_raw_matrix = sum(in_raw),
      n_signature_after_row_sum_filter = sum(passed_filter),
      n_signature_used = length(used_signature),
      coverage_fraction = coverage_fraction,
      minimum_coverage_fraction = MIN_SIGNATURE_COVERAGE,
      coverage_pass = coverage_pass
    )
    gene_coverage_rows[[length(gene_coverage_rows) + 1L]] <- data.table(
      dataset = cohort_id,
      gene_id_base = signature_all,
      gene_symbol = signature$gene_symbol[match(signature_all, signature$gene_id_base)],
      in_raw_matrix = in_raw,
      passed_row_sum_filter = passed_filter,
      used_in_signature_pc1 = signature_all %in% used_signature,
      missing_reason = fifelse(!in_raw, "absent_from_resource_dge",
                        fifelse(!passed_filter, "failed_cohort_row_sum_filter", "used"))
    )

    axis_inputs <- list(
      signature_pc1 = list(
        expression = expression[used_signature, , drop = FALSE],
        frozen_loading = signature_loading,
        n_signature_present = length(used_signature),
        coverage_fraction = coverage_fraction
      ),
      full_transcriptome_pc1 = list(
        expression = expression,
        frozen_loading = full_loading,
        n_signature_present = length(used_signature),
        coverage_fraction = coverage_fraction
      )
    )

    for (axis_id in names(axis_inputs)) {
      axis_input <- axis_inputs[[axis_id]]
      fit <- fit_oriented_pc1(
        axis_input$expression, axis_input$frozen_loading, axis_id, cohort_id
      )
      row_index <- row_index + 1L
      score_rows[[row_index]] <- data.table(
        sample_id = names(fit$score),
        dataset = cohort_id,
        axis_id = axis_id,
        axis_raw = as.numeric(fit$score),
        axis_percentile = percentile_rank(fit$score),
        pc1_variance_pct = fit$pc1_variance_pct,
        n_genes = nrow(axis_input$expression),
        n_signature_present = axis_input$n_signature_present,
        n_signature_total = length(signature_all),
        coverage_fraction = axis_input$coverage_fraction,
        orientation_cosine = fit$oriented_cosine,
        orientation_sign = fit$orientation_sign
      )

      gene_ids <- names(fit$raw_loading)
      loading_rows[[row_index]] <- data.table(
        dataset = cohort_id,
        axis_id = axis_id,
        gene_id_base = gene_ids,
        loading = as.numeric(fit$raw_loading),
        oriented_loading = as.numeric(fit$oriented_loading),
        frozen_loading = as.numeric(axis_input$frozen_loading[gene_ids])
      )
      alignment_rows[[row_index]] <- data.table(
        dataset = cohort_id,
        axis_id = axis_id,
        n_participants = length(cohort_samples),
        n_axis_genes = nrow(axis_input$expression),
        n_alignment_genes = fit$n_alignment_genes,
        raw_orientation_cosine = fit$raw_cosine,
        orientation_sign = fit$orientation_sign,
        oriented_cosine = fit$oriented_cosine,
        pc1_variance_pct = fit$pc1_variance_pct,
        pc1_pc2_sdev_ratio = fit$pc1_pc2_sdev_ratio,
        orientation_source = "frozen_discovery_loading",
        histology_used = FALSE
      )
    }
  }

  scores <- rbindlist(score_rows, use.names = TRUE)
  loadings <- rbindlist(loading_rows, use.names = TRUE)
  coverage <- rbindlist(coverage_rows, use.names = TRUE)
  gene_coverage <- rbindlist(gene_coverage_rows, use.names = TRUE)
  alignment <- rbindlist(alignment_rows, use.names = TRUE)

  assert_true(nrow(scores) == 2L * nrow(meta),
              "Participant score output is not complete for both axes")
  assert_true(!anyDuplicated(scores[, .(sample_id, axis_id)]),
              "Participant score output has duplicate sample-axis rows")
  assert_true(all(scores$axis_percentile >= 0 & scores$axis_percentile <= 1),
              "Within-cohort percentiles must lie between zero and one")
  assert_true(all(scores[, .N, by = .(dataset, axis_id)]$N ==
                    scores[, uniqueN(sample_id), by = .(dataset, axis_id)]$V1),
              "A participant was duplicated within an axis/cohort")
  assert_true(nrow(coverage) == length(EXPECTED_COHORTS) && all(coverage$coverage_pass),
              "Coverage audit is incomplete or contains a failure")
  assert_true(nrow(alignment) == 2L * length(EXPECTED_COHORTS) &&
                all(alignment$oriented_cosine > 0) &&
                !any(alignment$histology_used),
              "Alignment audit failed")

  setorder(scores, dataset, axis_id, sample_id)
  setorder(loadings, dataset, axis_id, gene_id_base)
  setorder(coverage, dataset)
  setorder(gene_coverage, dataset, gene_id_base)
  setorder(alignment, dataset, axis_id)
  write_tsv_once(scores, file.path(output_dir, "participant_scores.tsv"))
  write_tsv_once(loadings, file.path(output_dir, "loadings.tsv"))
  write_tsv_once(coverage, file.path(output_dir, "coverage.tsv"))
  write_tsv_once(gene_coverage, file.path(output_dir, "signature_gene_coverage.tsv"))
  write_tsv_once(alignment, file.path(output_dir, "alignment_diagnostics.tsv"))
  write_tsv_once(data.table(
    parameter = c(
      "seed", "row_sum_filter", "transform", "quantile_normalization_scope",
      "cross_cohort_combat", "histology_protection", "signature_total",
      "signature_resource_observed", "minimum_signature_coverage",
      "signature_axis_role", "full_transcriptome_axis_role"
    ),
    value = c(
      as.character(HAC_SEED), "rowSums(counts) > n_samples", "log2(counts + 1)",
      "within_cohort", "none", "none", as.character(SIGNATURE_TOTAL_EXPECTED),
      as.character(SIGNATURE_RESOURCE_EXPECTED), as.character(MIN_SIGNATURE_COVERAGE),
      "co_primary", "control"
    )
  ), file.path(output_dir, "parameters.tsv"))
  session_path <- file.path(output_dir, "sessionInfo.txt")
  assert_true(!file.exists(session_path), "Refusing to overwrite output: ", session_path)
  writeLines(capture.output(sessionInfo()), session_path, useBytes = TRUE)

  cat("UNSUPERVISED_AXES_COMPLETE\t", output_dir, "\n", sep = "")
}

if (sys.nframe() == 0L) main()
