#!/usr/bin/env Rscript
# Native-count diagnostic of the proposed normalization and projection variants.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(NMF); library(nnls); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || file.exists(args[2L])) stop("Usage: prepared_protocol.rds new_output_directory")
protocol <- readRDS(args[1L])
state <- protocol$historical_preprocessing
source("scripts/analysis/histology_anchored_continuum/molecular_layers/fixed_nmf_count_transport.R")
cache <- readRDS(state$source_paths$cache)
original <- readRDS(state$source_paths$dge)
metadata <- fread(state$source_paths$metadata)
metadata <- metadata[match(colnames(original), sample_id)]
stopifnot(identical(metadata$sample_id, colnames(original)))
sizes <- setNames(original$samples$lib.size, colnames(original))
factors <- fixed_reference_tmm(original$counts[protocol$TMM_reference_gene_ids, , drop = FALSE], sizes,
  protocol$TMM_reference_counts, protocol$TMM_reference_library_size, protocol$TMM_reference_norm_factor)
retained_ids <- protocol$basis_feature_mask$gene_id[protocol$basis_feature_mask$retained]
expression <- fixed_nmf_transport_expression(original$counts[retained_ids, , drop = FALSE], sizes,
  factors, state, metadata$dataset)
contract <- fromJSON("scripts/analysis/histology_anchored_continuum/molecular_layers/00_contract.json")
signature_path <- file.path(contract$source_hac_candidate, "programs/signature_gencode_v49_mapping.tsv")
signature <- fread(signature_path)
column <- if ("gencode_v49_gene_symbol" %in% names(signature)) "gencode_v49_gene_symbol" else "released_gene_symbol"
signature_symbols <- unique(toupper(signature[[column]]))
annotation <- fread(contract$gene_annotation)
symbols <- toupper(annotation$gene_name[match(rownames(cache$mat_nn), annotation$gene_id)])
exclude_signature <- !is.na(symbols) & symbols %in% signature_symbols
sample_index <- match(colnames(cache$mat_nn), colnames(expression))
stopifnot(!anyNA(sample_index), !anyNA(metadata$dataset))
history_cohorts <- metadata$dataset[sample_index]
nnls_matrix <- function(w, x) vapply(seq_len(ncol(x)), function(i) nnls::nnls(w, x[, i])$x, numeric(ncol(w)))
checks <- list()
score_rows <- list()
for (k in c(4L, 6L)) {
  full_basis <- NMF::basis(cache$nmf_results[[as.character(k)]])
  for (variant in c("full_stage", "signature_excluded_continuum")) {
    original_keep <- if (variant == "full_stage") rep(TRUE, nrow(full_basis)) else !exclude_signature
    revised_keep <- original_keep & protocol$basis_feature_mask$retained
    target_ids <- rownames(full_basis)[revised_keep]
    w <- full_basis[revised_keep, , drop = FALSE]
    reference <- if (variant == "full_stage") NMF::coef(cache$nmf_results[[as.character(k)]]) else
      nnls_matrix(full_basis[original_keep, , drop = FALSE], cache$mat_nn[original_keep, , drop = FALSE])
    h <- nnls_matrix(w, expression[target_ids, , drop = FALSE])
    for (j in seq_len(k)) {
      global <- cor(h[j, sample_index], reference[j, ], method = "spearman")
      per_cohort <- vapply(sort(unique(history_cohorts)), function(cohort) {
        selected <- history_cohorts == cohort
        cor(h[j, sample_index[selected]], reference[j, selected], method = "spearman")
      }, numeric(1))
      mass <- sum(w[, j]) / sum(full_basis[original_keep, j])
      checks[[length(checks)+1L]] <- data.table(k = k, program = paste0("P", j), variant = variant,
        retained_rows = nrow(w), retained_basis_mass = mass, retained_basis_rank = qr(w)$rank,
        historical_global_spearman = global, historical_min_cohort_spearman = min(per_cohort),
        transport_reconstruction_eligible = global >= .99 && all(per_cohort >= .9) && mass >= .8 && qr(w)$rank == k,
        all_native_scores_finite_nonnegative = all(is.finite(h[j, ])) && all(h[j, ] >= 0))
      score_rows[[length(score_rows)+1L]] <- data.table(sample_id = colnames(expression), dataset = metadata$dataset,
        k = k, program = paste0("P", j), variant = variant, native_count_diagnostic_loading = h[j, ],
        diagnostic_not_corrected_projection = TRUE)
    }
  }
}
checks <- rbindlist(checks)
dir.create(args[2L], recursive = TRUE)
fwrite(checks, file.path(args[2L], "native_normalization_transport_checks.tsv"), sep = "\t")
fwrite(rbindlist(score_rows), file.path(args[2L], "native_diagnostic_projection_coefficients.tsv.gz"), sep = "\t")
fwrite(data.table(sample_id = names(factors), TMM_factor = unname(factors), historical_library_size = sizes),
  file.path(args[2L], "native_fixed_reference_normalization.tsv"), sep = "\t")
report <- list(status = "native_count_transport_diagnostic_only", corrected_scores_produced = FALSE, bases_refitted = FALSE,
  normalization_reference = protocol$TMM_reference_sample_id, original_feature_definition_changed = TRUE,
  native_libraries = ncol(original), original_disease_training_libraries = ncol(cache$mat_nn),
  other_native_libraries_without_original_disease_TMM = ncol(original) - ncol(cache$mat_nn),
  full_stage_axes_meeting_criteria = sum(checks[variant == "full_stage", transport_reconstruction_eligible]),
  continuum_axes_meeting_criteria = sum(checks[variant == "signature_excluded_continuum", transport_reconstruction_eligible]),
  no_outcome_model_fitted = TRUE, signature_mapping_path = normalizePath(signature_path),
  limitation = "Native-count reconstruction cannot certify corrected-count transport. All revised axes retain their original fixed weights; the mask and normalization change projection definition and require scientific review.")
writeLines(toJSON(report, auto_unbox = TRUE, pretty = TRUE), file.path(args[2L], "native_transport_review.json"))
writeLines(capture.output(sessionInfo()), file.path(args[2L], "sessionInfo.txt"))
print(checks)
cat(toJSON(report, auto_unbox = TRUE), "\n")
