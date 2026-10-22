# Outcome-independent preparation for a revised fixed-W projection.
# These projections are a sensitivity definition, not the original 5,000-row axes.

fixed_nmf_counting_mask <- function(basis_ids, correspondence) {
  stopifnot(!anyDuplicated(basis_ids), !anyDuplicated(correspondence$gene_id))
  idx <- match(basis_ids, correspondence$gene_id)
  changed_targets <- unique(correspondence$counting_gene_id[
    correspondence$disposition == "equivalent" &
      correspondence$gene_id != correspondence$counting_gene_id])
  d <- correspondence[idx, , drop = FALSE]
  reason <- rep("unchanged_one_to_one_counting_object", length(basis_ids))
  reason[!is.na(idx) & d$counting_gene_id %in% changed_targets] <- "counting_target_receives_merged_ALT_features"
  reason[!is.na(idx) & d$counting_gene_id != basis_ids] <- "source_feature_changes_counting_object"
  reason[!is.na(idx) & !d$disposition %in% c("primary_gene", "distinct_locus", "alternative_only")] <- "excluded_unresolved_or_nonunique_source_object"
  reason[is.na(idx)] <- "exact_basis_gene_version_absent"
  data.frame(gene_id = basis_ids, retained = reason == "unchanged_one_to_one_counting_object",
             reason = reason, stringsAsFactors = FALSE)
}

fixed_reference_tmm <- function(counts, library_sizes, reference_counts,
                                reference_library_size, reference_norm_factor) {
  stopifnot(is.matrix(counts), identical(rownames(counts), names(reference_counts)),
    all(is.finite(counts)), all(counts >= 0), all(is.finite(reference_counts)),
    all(reference_counts >= 0), length(library_sizes) == ncol(counts),
    all(is.finite(library_sizes)), all(library_sizes > 0),
    is.finite(reference_library_size), reference_library_size > 0,
    is.finite(reference_norm_factor), reference_norm_factor > 0)
  # Pairwise factors are centered geometrically by edgeR. Their ratio cancels
  # that centering, so a sample's normalization cannot depend on other samples.
  factors <- vapply(seq_len(ncol(counts)), function(i) {
    pair <- cbind(reference_counts, counts[, i])
    nf <- edgeR::calcNormFactors(pair, lib.size = c(reference_library_size, library_sizes[i]),
                               method = "TMM", refColumn = 1L)
    nf[2L] / nf[1L] * reference_norm_factor
  }, numeric(1))
  stopifnot(all(is.finite(factors)), all(factors > 0))
  setNames(factors, colnames(counts))
}

fixed_nmf_transport_expression <- function(counts, library_sizes, factors, state, cohorts) {
  stopifnot(identical(names(library_sizes), colnames(counts)),
    identical(names(factors), colnames(counts)), length(cohorts) == ncol(counts),
    all(rownames(counts) %in% names(state$row_minimum)),
    all(cohorts %in% colnames(state$cohort_batch_subtraction)))
  # Same 0.5-count and one-library-count offsets as installed limma::voom.
  expression <- log2(sweep(counts + .5, 2L, library_sizes * factors + 1, "/") * 1e6)
  expression <- expression - state$cohort_batch_subtraction[rownames(counts), cohorts, drop = FALSE]
  expression <- sweep(expression, 1L, state$row_minimum[rownames(counts)], "-")
  expression[expression < 0] <- 0
  stopifnot(all(is.finite(expression)))
  expression
}

# Reconstruction criteria fixed before any corrected count was scored. They are
# never relaxed; an axis that fails keeps its row and carries NA scores.
fixed_nmf_axis_eligible <- function(global_spearman, min_cohort_spearman, basis_mass, basis_rank, k) {
  global_spearman >= .99 & min_cohort_spearman >= .9 & basis_mass >= .8 & basis_rank == k
}

# Axis names are "k<k>-P<j>", matching protocol$stage_unavailable_axes.
fixed_nmf_unavailable_axes <- function(checks, variant_name) {
  checks <- as.data.frame(checks)
  selected <- checks[checks$variant == variant_name & !checks$transport_reconstruction_eligible, , drop = FALSE]
  if (!nrow(selected)) return(character())
  paste0("k", selected$k, "-", selected$program)
}

fixed_nmf_nnls_scores <- function(w, x) {
  stopifnot(requireNamespace("nnls", quietly = TRUE), nrow(w) == nrow(x),
            identical(rownames(w), rownames(x)), all(is.finite(w)), all(is.finite(x)))
  h <- vapply(seq_len(ncol(x)), function(i) nnls::nnls(w, x[, i])$x, numeric(ncol(w)))
  dimnames(h) <- list(paste0("P", seq_len(ncol(w))), colnames(x))
  stopifnot(all(is.finite(h)), all(h >= 0))
  h
}

# scores: programs x samples for one k. Unavailable axes become NA for every
# sample; no zero or outcome-derived replacement is ever written.
fixed_nmf_mask_unavailable <- function(scores, k, unavailable) {
  stopifnot(identical(rownames(scores), paste0("P", seq_len(nrow(scores)))), nrow(scores) == k)
  scores[paste0("k", k, "-", rownames(scores)) %in% unavailable, ] <- NA_real_
  scores
}
