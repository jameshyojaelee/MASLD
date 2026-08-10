# Pure helpers shared by the BG-001 scientific comparison and its mutation tests.

bg001_table_equal_by_key <- function(reference, candidate, key, ignored = character()) {
  reference <- data.table::as.data.table(data.table::copy(reference))
  candidate <- data.table::as.data.table(data.table::copy(candidate))
  fields <- setdiff(union(names(reference), names(candidate)), ignored)
  if (!key %in% fields || !all(fields %in% names(reference)) || !all(fields %in% names(candidate)) ||
      anyDuplicated(reference[[key]]) || anyDuplicated(candidate[[key]]) ||
      !setequal(reference[[key]], candidate[[key]])) {
    return(FALSE)
  }
  data.table::setorderv(reference, key)
  data.table::setorderv(candidate, key)
  # Normalize only container/row-name attributes introduced by data.table.
  # Column order, storage types, missingness, and every serialized scalar value
  # remain exact; no floating tolerance is permitted for a locked substrate.
  left <- as.data.frame(reference[, ..fields], stringsAsFactors = FALSE)
  right <- as.data.frame(candidate[, ..fields], stringsAsFactors = FALSE)
  rownames(left) <- NULL
  rownames(right) <- NULL
  identical(left, right)
}

bg001_missing_model_sex <- function(metadata, sample_ids) {
  metadata <- data.table::as.data.table(metadata)
  matched <- metadata$inferred_sex[match(sample_ids, metadata$sample_id)]
  sum(is.na(matched) | !nzchar(trimws(as.character(matched))))
}

bg001_nonfinite_counts <- function(table, fields) {
  missing <- setdiff(fields, names(table))
  if (length(missing)) stop("Missing required model fields: ", paste(missing, collapse = ", "), call. = FALSE)
  vapply(fields, function(field) {
    values <- suppressWarnings(as.numeric(table[[field]]))
    sum(!is.finite(values))
  }, integer(1))
}

bg001_expression_gate_summary <- function(model, selector, median_bound, p95_bound) {
  rows <- model[selector]
  delta <- abs(rows$delta_logFC[is.finite(rows$delta_logFC)])
  if (!length(delta)) {
    return(list(
      n = 0L,
      median_abs_delta_logFC = NA_real_,
      p95_abs_delta_logFC = NA_real_,
      median_pass = FALSE,
      p95_pass = FALSE,
      overall_pass = FALSE
    ))
  }
  median_value <- median(delta)
  p95_value <- as.numeric(quantile(delta, 0.95, names = FALSE, type = 8))
  median_pass <- is.finite(median_value) && median_value <= median_bound
  p95_pass <- is.finite(p95_value) && p95_value <= p95_bound
  list(
    n = length(delta),
    median_abs_delta_logFC = median_value,
    p95_abs_delta_logFC = p95_value,
    median_pass = median_pass,
    p95_pass = p95_pass,
    overall_pass = median_pass && p95_pass
  )
}

bg001_protected_direction_reversals <- function(model) {
  required <- c("protected", "logFC_old", "logFC_new")
  missing <- setdiff(required, names(model))
  if (length(missing)) {
    stop(
      "Protected-direction check lacks fields: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  eligible <- !is.na(model$protected) & model$protected &
    is.finite(model$logFC_old) & is.finite(model$logFC_new)
  as.integer(sum(sign(model$logFC_old[eligible]) != sign(model$logFC_new[eligible])))
}

bg001_gate_failure_specification <- function(gate) {
  specifications <- list(
    G1 = list(
      job = "qc_sex_reconstruction_and_refit",
      action = "Inspect QC/sex reconstruction and refit the implicated sample/cohort substrate."
    ),
    G2 = list(
      job = "dependency_full_raw_dge_rebuild",
      action = "Rebuild every active raw-count/DGE consumer; global coefficient correlation cannot waive gene-universe failure."
    ),
    G3 = list(
      job = "normalization_weight_localization",
      action = "Localize by cohort/group and inspect RLE factors, weights, and voom curves before targeted refits."
    ),
    G4 = list(
      job = "limma_leave_one_cohort_out",
      action = "Localize coefficient instability and run leave-one-cohort-out fits for implicated cohorts."
    ),
    G5 = list(
      job = "inspect_treat_threshold_and_dependencies",
      action = "Inspect threshold proximity and rebuild all dependency-selected TREAT consumers."
    )
  )
  if (length(gate) != 1L || !gate %in% names(specifications)) {
    stop("Unknown scientific-gate failure: ", paste(gate, collapse = ","), call. = FALSE)
  }
  specifications[[gate]]
}
