#!/usr/bin/env Rscript
# Prepare a concrete revised projection protocol; never emit corrected scores.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(NMF); library(nnls); library(jsonlite)})
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 2L || file.exists(argv[2L])) stop("Usage: annotation_directory new_output_directory")
root <- normalizePath(".")
candidate <- file.path(root, "figures/candidates/resource-corrections-20261005-v1")
state_path <- file.path(candidate, "bulk_execution_sources/fixed_nmf_preprocessing_review/frozen_preprocessing_state.rds")
state <- readRDS(state_path)
sha <- function(p) strsplit(system2("sha256sum", shQuote(p), stdout = TRUE), "[[:space:]]+")[[1]][[1]]
stopifnot(all(vapply(unlist(state$source_paths), sha, character(1)) == unlist(state$source_sha256)))
helper <- file.path(root, "scripts/analysis/histology_anchored_continuum/molecular_layers/fixed_nmf_count_transport.R")
source(helper)
annotation_path <- file.path(argv[1L], "annotation_correspondence.tsv")
summary <- fromJSON(file.path(argv[1L], "resolution_summary.json"))
stopifnot(sha(annotation_path) == summary$annotation_correspondence_sha256)
annotation <- as.data.frame(fread(annotation_path, colClasses = "character"))
cache <- readRDS(state$source_paths$cache)
stopifnot(identical(rownames(cache$mat_nn), state$exact_gene_ids))
mask <- fixed_nmf_counting_mask(state$exact_gene_ids, annotation)
original <- readRDS(state$source_paths$dge)
normalization_mask <- fixed_nmf_counting_mask(rownames(original), annotation)
normalization_ids <- normalization_mask$gene_id[normalization_mask$retained]
roster_path <- file.path(candidate, "identity_resource_qualified/bulk_execution_inputs/primary_bulk_metadata.tsv")
roster <- fread(roster_path)
reference_candidates <- roster[group_binary == "Disease" & sample_id %in% state$historical_sample_ids]
reference_candidates <- reference_candidates[order(-full_called, sample_id)]
stopifnot(nrow(reference_candidates) > 0)
reference_id <- reference_candidates$sample_id[1L]
reference_counts <- original$counts[normalization_ids, reference_id]
names(reference_counts) <- normalization_ids
reference_library_size <- state$historical_library_sizes[[reference_id]]
reference_norm_factor <- state$historical_TMM_sample_factors[[reference_id]]
metadata <- fread(state$source_paths$metadata)
cohorts <- metadata$dataset[match(colnames(cache$mat_nn), metadata$sample_id)]
stopifnot(!anyNA(cohorts))
protocol_basis <- list()
checks <- list()
for (k in c(4L, 6L)) {
  fit <- cache$nmf_results[[as.character(k)]]
  basis <- NMF::basis(fit)
  native <- NMF::coef(fit)
  stopifnot(identical(rownames(basis), mask$gene_id), identical(colnames(native), colnames(cache$mat_nn)))
  retained <- basis[mask$retained, , drop = FALSE]
  fraction <- colSums(retained) / colSums(basis)
  coefficients <- vapply(seq_len(ncol(cache$mat_nn)), function(i) {
    nnls::nnls(retained, cache$mat_nn[mask$retained, i])$x
  }, numeric(k))
  stopifnot(all(is.finite(coefficients)), all(coefficients >= 0))
  for (j in seq_len(k)) {
    global_rho <- cor(coefficients[j, ], native[j, ], method = "spearman")
    per_cohort <- vapply(sort(unique(cohorts)), function(cohort) {
      selected <- cohorts == cohort
      if (sum(selected) < 3L) return(NA_real_)
      cor(coefficients[j, selected], native[j, selected], method = "spearman")
    }, numeric(1))
    checks[[length(checks) + 1L]] <- data.table(k = k, program = paste0("P", j),
      retained_rows = nrow(retained), retained_basis_mass = fraction[j],
      retained_basis_rank = qr(retained)$rank, historical_global_spearman = global_rho,
      historical_min_cohort_spearman = min(per_cohort, na.rm = TRUE),
      reconstruction_eligible = is.finite(global_rho) && global_rho >= .99 &&
        fraction[j] >= .8 && qr(retained)$rank == k &&
        all(is.finite(per_cohort) & per_cohort >= .9))
  }
  protocol_basis[[as.character(k)]] <- retained
}
checks <- rbindlist(checks)
protocol <- list(state = "revised_fixed_W_transport_prepared_for_review_not_original_axis_definition",
  bases_refitted = FALSE, corrected_scores_produced = FALSE,
  basis_feature_mask = mask, basis = protocol_basis, historical_preprocessing = state,
  TMM_reference_sample_id = reference_id, TMM_reference_gene_ids = normalization_ids,
  TMM_reference_counts = reference_counts, TMM_reference_library_size = reference_library_size,
  TMM_reference_norm_factor = reference_norm_factor,
  reference_selection = "highest full_called among primary disease units in original training; lexical sample ID tie break",
  signature_excluded_continuum = "intersect unchanged-object mask with the existing fixed signature exclusion; recompute that variant's checks separately",
  correction_requirements = "complete annotation, certified full corrected count matrix, exact versions and eligible roster; no zero substitution",
  normalization = "edgeR pairwise TMM against serialized historical reference on unchanged one-to-one objects, ratio anchored to original reference factor; full corrected assigned library totals",
  preprocessing = "voom log2((count+0.5)/(effective library size+1)*1e6); subtract fixed historical cohort shifts and minima; clip below zero",
  adoption = "original feature definition changes; review required before scientific adoption even when reconstruction criteria pass")
dir.create(argv[2L], recursive = TRUE)
saveRDS(protocol, file.path(argv[2L], "revised_fixed_nmf_projection_protocol.rds"))
fwrite(as.data.table(mask), file.path(argv[2L], "basis_counting_object_mask.tsv"), sep = "\t")
fwrite(checks, file.path(argv[2L], "historical_masked_projection_checks.tsv"), sep = "\t")
fwrite(as.data.table(normalization_mask), file.path(argv[2L], "normalization_counting_object_mask.tsv"), sep = "\t")
report <- list(status = protocol$state, bases_refitted = FALSE, corrected_scores_produced = FALSE,
  full_stage_basis_rows = nrow(mask), retained_rows = sum(mask$retained), masked_rows = sum(!mask$retained),
  normalization_rows = length(normalization_ids), TMM_reference_sample_id = reference_id,
  n_axes_meeting_historical_reconstruction_criteria = sum(checks$reconstruction_eligible),
  n_axes = nrow(checks), annotation_unresolved = summary$n_unresolved_annotations,
  original_axes_replaced = FALSE, protocol_for_review = TRUE,
  source_sha256 = list(annotation = sha(annotation_path), preprocessing = sha(state_path), helper = sha(helper),
    primary_roster = sha(roster_path)),
  limitation = "Historical masked reconstruction is not corrected-count transport validation. Continuum signature-excluded reconstruction and controls must be checked separately; no outcome fit or corrected projection has run.")
writeLines(toJSON(report, auto_unbox = TRUE, pretty = TRUE), file.path(argv[2L], "projection_preparation.json"))
writeLines(capture.output(sessionInfo()), file.path(argv[2L], "sessionInfo.txt"))
print(checks)
cat(toJSON(report, auto_unbox = TRUE), "\n")
