#!/usr/bin/env Rscript
# Corrected-count fixed-W NMF scores under the reviewed r2 transport protocol.
# Weights are the original k4/k6 bases restricted to the protocol's retained rows;
# nothing is refitted and no outcome is read. Axes that failed the fixed
# reconstruction criteria on historical counts keep their rows and carry NA.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(NMF); library(nnls); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L || file.exists(args[5L])) {
  stop("Usage: protocol_r2.rds assembled_count_directory unit_metadata.tsv model_input_manifest.tsv new_output_directory")
}
set.seed(20260923)
options(warn = 1)
root <- normalizePath(".")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "fixed_nmf_count_transport.R"))
sha <- function(p) strsplit(system2("sha256sum", shQuote(p), stdout = TRUE), "[[:space:]]+")[[1L]][[1L]]
as_count_matrix <- function(x) {
  if (inherits(x, "DGEList")) x <- x$counts
  if (inherits(x, "Matrix")) x <- as.matrix(x)
  stopifnot(is.matrix(x), is.numeric(x))
  x
}

# ---- protocol: re-derive the axis decisions instead of trusting the stored vector
protocol_path <- normalizePath(args[1L], mustWork = TRUE)
protocol <- readRDS(protocol_path)
state <- protocol$historical_preprocessing
checks <- as.data.frame(protocol$historical_normalization_transport_checks)
stopifnot(identical(protocol$state, "revised_fixed_W_transport_prepared_for_review_not_original_axis_definition"),
  isFALSE(protocol$bases_refitted), nrow(checks) == 20L, is.character(protocol$unavailable_score_policy))
rederived <- with(checks, fixed_nmf_axis_eligible(historical_global_spearman,
  historical_min_cohort_spearman, retained_basis_mass, retained_basis_rank, k))
stopifnot(identical(rederived, checks$transport_reconstruction_eligible))
variants <- c("full_stage", "signature_excluded_continuum")
unavailable <- setNames(lapply(variants, function(v) fixed_nmf_unavailable_axes(checks, v)), variants)
stopifnot(identical(unavailable$full_stage, protocol$stage_unavailable_axes),
          identical(unavailable$signature_excluded_continuum, protocol$continuum_unavailable_axes))

# ---- certified corrected counts and the unit roster
count_dir <- normalizePath(args[2L], mustWork = TRUE)
counts_path <- file.path(count_dir, "counts_raw.rds")
completion <- fromJSON(file.path(count_dir, "count_completion.json"))
stopifnot(isTRUE(completion$counts_complete), isTRUE(completion$annotation_resolution_complete),
  completion$n_unresolved_annotations == 0L, identical(completion$count_arm, "primary"),
  identical(sha(counts_path), completion$counts_sha256))
counts <- as_count_matrix(readRDS(counts_path))
stopifnot(ncol(counts) == completion$n_samples, nrow(counts) == completion$n_genes,
  all(is.finite(counts)), all(counts >= 0), !anyDuplicated(colnames(counts)), !anyDuplicated(rownames(counts)))
roster_path <- normalizePath(args[3L], mustWork = TRUE)
roster <- fread(roster_path)
unit_columns <- c("sample_id", "dataset", "analysis_unit_id", "unit_evidence")
stopifnot(all(unit_columns %in% names(roster)), !anyNA(roster[, ..unit_columns]),
  !anyDuplicated(roster$sample_id), !anyDuplicated(roster$analysis_unit_id),
  setequal(roster$sample_id, colnames(counts)))
roster <- roster[match(colnames(counts), sample_id), ..unit_columns]
stopifnot(all(roster$dataset %in% colnames(state$cohort_batch_subtraction)))

# ---- the refit's model counts must equal the certified counts on every shared gene
model_manifest_path <- normalizePath(args[4L], mustWork = TRUE)
model <- fread(model_manifest_path)
stopifnot(nrow(model) == 1L, isTRUE(as.logical(model$identity_qualified)), model$n_samples == ncol(counts))
model_dge_path <- normalizePath(model$dge_path, mustWork = TRUE)
model_dge <- readRDS(model_dge_path)
stopifnot(setequal(colnames(model_dge), colnames(counts)), all(rownames(model_dge) %in% rownames(counts)))
model_counts <- as_count_matrix(model_dge)[, colnames(counts), drop = FALSE]
stopifnot(all(model_counts == counts[rownames(model_counts), , drop = FALSE]))
rm(model_dge, model_counts)

# ---- normalization against the serialized historical reference; full assigned totals
sizes <- colSums(counts)
stopifnot(all(sizes > 0))
reference_ids <- protocol$TMM_reference_gene_ids
retained_ids <- protocol$basis_feature_mask$gene_id[protocol$basis_feature_mask$retained]
stopifnot(all(reference_ids %in% rownames(counts)), all(retained_ids %in% rownames(counts)))
factors <- fixed_reference_tmm(counts[reference_ids, , drop = FALSE], sizes, protocol$TMM_reference_counts,
  protocol$TMM_reference_library_size, protocol$TMM_reference_norm_factor)
expression <- fixed_nmf_transport_expression(counts[retained_ids, , drop = FALSE], sizes, factors, state, roster$dataset)

# ---- fixed bases, signature exclusion for the continuum variant
cache <- readRDS(state$source_paths$cache)
contract <- fromJSON(file.path(script_dir, "00_contract.json"))
signature_path <- file.path(contract$source_hac_candidate, "programs/signature_gencode_v49_mapping.tsv")
signature <- fread(signature_path)
signature_column <- if ("gencode_v49_gene_symbol" %in% names(signature)) "gencode_v49_gene_symbol" else "released_gene_symbol"
annotation <- fread(contract$gene_annotation)
symbols <- toupper(annotation$gene_name[match(rownames(cache$mat_nn), annotation$gene_id)])
exclude_signature <- !is.na(symbols) & symbols %in% unique(toupper(signature[[signature_column]]))
nnls_matrix <- function(w, x) vapply(seq_len(ncol(x)), function(i) nnls::nnls(w, x[, i])$x, numeric(ncol(w)))
overlap_historical <- intersect(colnames(cache$mat_nn), roster$sample_id)

native_path <- file.path(protocol$transport_check_path, "native_diagnostic_projection_coefficients.tsv.gz")
native <- as.data.frame(fread(native_path))

wide <- list(); long <- list(); agreement <- list()
for (k in c(4L, 6L)) {
  fit <- cache$nmf_results[[as.character(k)]]
  full_basis <- NMF::basis(fit)
  stopifnot(identical(rownames(full_basis), protocol$basis_feature_mask$gene_id))
  for (variant in variants) {
    original_keep <- if (variant == "full_stage") rep(TRUE, nrow(full_basis)) else !exclude_signature
    revised_keep <- original_keep & protocol$basis_feature_mask$retained
    w <- full_basis[revised_keep, , drop = FALSE]
    if (variant == "full_stage") stopifnot(identical(w, protocol$basis[[as.character(k)]]))
    stopifnot(nrow(w) == checks$retained_rows[checks$k == k & checks$variant == variant][1L], qr(w)$rank == k)
    h <- fixed_nmf_nnls_scores(w, expression[rownames(w), , drop = FALSE])
    scores <- fixed_nmf_mask_unavailable(h, k, unavailable[[variant]])
    reference <- if (variant == "full_stage") NMF::coef(fit) else
      nnls_matrix(full_basis[original_keep, , drop = FALSE], cache$mat_nn[original_keep, , drop = FALSE])
    if (variant != "full_stage") colnames(reference) <- colnames(cache$mat_nn)
    for (j in seq_len(k)) {
      program <- paste0("P", j)
      axis <- paste0("k", k, "-", program)
      available <- !axis %in% unavailable[[variant]]
      row <- checks[checks$k == k & checks$program == program & checks$variant == variant, ]
      reason <- if (available) NA_character_ else paste(c(
        if (row$historical_global_spearman < .99) sprintf("historical global Spearman %.4f < 0.99", row$historical_global_spearman),
        if (row$historical_min_cohort_spearman < .9) sprintf("minimum cohort Spearman %.4f < 0.90", row$historical_min_cohort_spearman),
        if (row$retained_basis_mass < .8) sprintf("retained basis mass %.4f < 0.80", row$retained_basis_mass),
        if (row$retained_basis_rank != k) sprintf("retained basis rank %d < %d", row$retained_basis_rank, k)), collapse = "; ")
      long[[length(long) + 1L]] <- data.table(roster, k = k, program = program, variant = variant,
        score = scores[j, ], axis_available = available, unavailable_reason = reason)
      if (available) {
        ids <- overlap_historical
        global <- cor(scores[j, ids], reference[j, ids], method = "spearman")
        by_cohort <- vapply(sort(unique(roster$dataset)), function(d) {
          s <- intersect(ids, roster$sample_id[roster$dataset == d])
          if (length(s) < 3L) return(NA_real_)
          cor(scores[j, s], reference[j, s], method = "spearman")
        }, numeric(1))
        nat <- native[native$k == k & native$program == program & native$variant == variant, ]
        shared <- intersect(ids, nat$sample_id)
        native_rho <- cor(scores[j, shared], nat$native_count_diagnostic_loading[match(shared, nat$sample_id)], method = "spearman")
        agreement[[length(agreement) + 1L]] <- data.table(k = k, program = program, variant = variant,
          n_overlap_historical_libraries = length(ids), spearman_vs_historical_H = global,
          min_cohort_spearman_vs_historical_H = suppressWarnings(min(by_cohort, na.rm = TRUE)),
          n_overlap_native_diagnostic = length(shared), spearman_vs_native_count_diagnostic = native_rho,
          fraction_zero_scores = mean(scores[j, ] == 0))
      }
    }
    wide[[paste(variant, k, sep = "_k")]] <- data.table(sample_id = colnames(scores), t(scores))
  }
}
long <- rbindlist(long)
agreement <- rbindlist(agreement)

dir.create(args[5L], recursive = TRUE)
out <- normalizePath(args[5L])
for (name in names(wide)) {
  stopifnot(!anyDuplicated(wide[[name]]$sample_id))
  fwrite(wide[[name]], file.path(out, paste0("fixed_nmf_scores_", name, ".tsv")), sep = "\t", na = "NA")
}
fwrite(long, file.path(out, "fixed_nmf_scores_long.tsv.gz"), sep = "\t", na = "NA")
availability <- unique(long[, .(k, program, variant, axis_available, unavailable_reason)])
availability <- merge(availability, as.data.table(checks)[, .(k, program, variant, retained_rows,
  historical_global_spearman, historical_min_cohort_spearman, retained_basis_mass, retained_basis_rank)],
  by = c("k", "program", "variant"))
fwrite(availability, file.path(out, "axis_availability.tsv"), sep = "\t", na = "NA")
fwrite(agreement, file.path(out, "corrected_vs_historical_diagnostics.tsv"), sep = "\t")
fwrite(data.table(sample_id = colnames(counts), dataset = roster$dataset, assigned_library_size = sizes,
  TMM_factor = unname(factors[colnames(counts)]), effective_library_size = sizes * unname(factors[colnames(counts)])),
  file.path(out, "corrected_normalization.tsv"), sep = "\t")
summary_by_cohort <- long[axis_available == TRUE, .(n = .N, mean = mean(score), sd = sd(score),
  minimum = min(score), maximum = max(score), n_zero = sum(score == 0)), by = .(k, program, variant, dataset)]
fwrite(summary_by_cohort, file.path(out, "score_summary_by_cohort.tsv"), sep = "\t")

inputs <- c(protocol = protocol_path, corrected_counts = counts_path,
  count_certificate = file.path(count_dir, "count_completion.json"), unit_metadata = roster_path,
  refit_model_manifest = model_manifest_path, refit_model_dge = model_dge_path,
  historical_nmf_cache = state$source_paths$cache, signature_mapping = normalizePath(signature_path),
  gene_annotation = normalizePath(contract$gene_annotation), native_count_diagnostic = native_path,
  helper = file.path(script_dir, "fixed_nmf_count_transport.R"), scorer = normalizePath(script_file),
  contract = file.path(script_dir, "00_contract.json"))
dir.create(file.path(out, "provenance"))
fwrite(data.table(role = names(inputs), path = unname(inputs), sha256 = vapply(inputs, sha, character(1))),
  file.path(out, "provenance", "input_manifest.tsv"), sep = "\t")
report <- list(status = "corrected_fixed_nmf_scores_from_reviewed_sensitivity_protocol", bases_refitted = FALSE,
  fixed_criteria_relaxed = FALSE, outcomes_read = FALSE, original_axes_replaced = FALSE, adoption = "not adopted",
  n_units = nrow(roster), n_cohorts = uniqueN(roster$dataset), cohorts = as.list(table(roster$dataset)),
  annotation_unresolved = completion$n_unresolved_annotations, n_held_untestable_annotations = completion$n_held_untestable_annotations,
  stage_unavailable_axes = unavailable$full_stage, continuum_unavailable_axes = unavailable$signature_excluded_continuum,
  stage_rows = sum(protocol$basis_feature_mask$retained),
  continuum_rows = checks$retained_rows[checks$variant == "signature_excluded_continuum"][1L],
  unavailable_score_policy = protocol$unavailable_score_policy, seed = 20260923L,
  scores_finite_where_available = all(is.finite(long[axis_available == TRUE, score])),
  scores_na_only_on_unavailable_axes = all(is.na(long$score) == !long$axis_available))
stopifnot(report$scores_finite_where_available, report$scores_na_only_on_unavailable_axes)
writeLines(toJSON(report, auto_unbox = TRUE, pretty = TRUE), file.path(out, "scoring_report.json"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(availability[order(variant, k, program)])
print(agreement[order(variant, k, program)])
cat(toJSON(report, auto_unbox = TRUE), "\n")
