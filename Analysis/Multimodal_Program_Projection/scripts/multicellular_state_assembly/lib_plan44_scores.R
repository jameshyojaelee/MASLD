suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

read_frozen_plan44 <- function(candidate_root) {
  membership <- fread(file.path(
    candidate_root,
    "../chronic-state-risk-bridge-2026-08-09/frozen_inputs/",
    "program_membership__program_membership_v2.tsv"
  ))
  geometry <- fread(file.path(candidate_root, "frozen_lineage_geometry.tsv"))
  stopifnot(nrow(geometry) == 117L, uniqueN(geometry$program_uid) == 117L)
  stopifnot(uniqueN(membership$program_uid) == 117L)
  list(membership = membership, geometry = geometry)
}

aggregate_counts_by_symbol <- function(count_matrix, symbols) {
  stopifnot(nrow(count_matrix) == length(symbols))
  keep <- !is.na(symbols) & nzchar(symbols)
  count_matrix <- count_matrix[keep, , drop = FALSE]
  symbols <- symbols[keep]
  aggregated <- rowsum(count_matrix, group = symbols, reorder = FALSE, na.rm = TRUE)
  storage.mode(aggregated) <- "numeric"
  aggregated
}

tmm_logcpm <- function(count_matrix) {
  stopifnot(all(is.finite(count_matrix)), all(count_matrix >= 0))
  dge <- DGEList(counts = count_matrix)
  dge <- calcNormFactors(dge, method = "TMM")
  cpm(dge, log = TRUE, prior.count = 2)
}

deseq_sizefactor_log <- function(count_matrix) {
  stopifnot(all(is.finite(count_matrix)), all(count_matrix >= 0))
  rounded <- round(count_matrix)
  size_factors <- DESeq2::estimateSizeFactorsForMatrix(rounded)
  log2(sweep(rounded, 2L, size_factors, "/") + 1)
}

standardize_by_reference <- function(expression, strata, reference) {
  stopifnot(ncol(expression) == length(strata), length(strata) == length(reference))
  z <- matrix(
    NA_real_,
    nrow = nrow(expression),
    ncol = ncol(expression),
    dimnames = dimnames(expression)
  )
  for (level in unique(strata)) {
    target_columns <- which(strata == level)
    reference_columns <- which(strata == level & reference)
    if (length(reference_columns) < 2L) {
      stop(sprintf("fewer than two reference samples in stratum %s", level))
    }
    center <- rowMeans(expression[, reference_columns, drop = FALSE])
    scale <- apply(expression[, reference_columns, drop = FALSE], 1L, sd)
    usable <- is.finite(center) & is.finite(scale) & scale > 1e-8
    z[usable, target_columns] <- sweep(
      sweep(expression[usable, target_columns, drop = FALSE], 1L, center[usable], "-"),
      1L,
      scale[usable],
      "/"
    )
  }
  z
}

standardize_global_reference <- function(expression, reference) {
  stopifnot(ncol(expression) == length(reference), sum(reference) >= 3L)
  center <- rowMeans(expression[, reference, drop = FALSE])
  scale <- apply(expression[, reference, drop = FALSE], 1L, sd)
  usable <- is.finite(center) & is.finite(scale) & scale > 1e-8
  z <- matrix(
    NA_real_,
    nrow = nrow(expression),
    ncol = ncol(expression),
    dimnames = dimnames(expression)
  )
  z[usable, ] <- sweep(
    sweep(expression[usable, , drop = FALSE], 1L, center[usable], "-"),
    1L,
    scale[usable],
    "/"
  )
  z
}

score_frozen_geometry <- function(z, membership, geometry, mode = "primary") {
  allowed_modes <- c(
    "primary",
    "equal_gene",
    "leave_top_gene",
    "sign_loading",
    "leave_top_program",
    "shared_gene_removed"
  )
  stopifnot(mode %in% allowed_modes)
  available_symbols <- rownames(z)[rowSums(is.finite(z)) == ncol(z)]
  work <- copy(membership)
  work <- work[mapped_symbol %in% available_symbols]

  if (mode == "shared_gene_removed") {
    primary_members <- unique(
      membership[
        cell_type %in% c("hepatocytes", "fibroblasts", "macrophages"),
        .(cell_type, mapped_symbol)
      ]
    )
    shared <- primary_members[, .(n_lineages = uniqueN(cell_type)), by = mapped_symbol][
      n_lineages > 1L,
      mapped_symbol
    ]
    work <- work[!mapped_symbol %in% shared]
  }
  if (mode == "leave_top_gene") {
    work[, remove_gene := mapped_symbol[which.max(original_l1_weight)], by = program_uid]
    work <- work[mapped_symbol != remove_gene]
    work[, remove_gene := NULL]
  }

  program_rows <- vector("list", nrow(geometry))
  program_scores <- matrix(
    NA_real_,
    nrow = nrow(geometry),
    ncol = ncol(z),
    dimnames = list(geometry$program_uid, colnames(z))
  )
  for (index in seq_len(nrow(geometry))) {
    uid <- geometry$program_uid[[index]]
    members <- work[program_uid == uid]
    retained_weight <- sum(members$original_l1_weight)
    n_genes <- uniqueN(members$mapped_symbol)
    testable <- n_genes >= 8L && retained_weight >= 0.20
    program_rows[[index]] <- data.table(
      program_uid = uid,
      cell_type = geometry$cell_type[[index]],
      module = geometry$module[[index]],
      mode = mode,
      n_mapped_genes = n_genes,
      retained_original_l1_weight = retained_weight,
      testable = testable
    )
    if (!testable) next
    members <- members[!duplicated(mapped_symbol)]
    weights <- if (mode == "equal_gene") {
      rep(1, nrow(members))
    } else {
      members$original_l1_weight
    }
    values <- z[members$mapped_symbol, , drop = FALSE]
    program_scores[index, ] <- colSums(values * weights) / sum(weights)
  }
  program_testability <- rbindlist(program_rows)

  geometry_work <- copy(geometry)
  geometry_work[, lineage_loading := as.numeric(within_lineage_loading)]
  if (mode == "sign_loading") {
    geometry_work[, lineage_loading := sign(as.numeric(stage_beta)) / .N, by = cell_type]
  }
  if (mode == "leave_top_program") {
    geometry_work[
      , remove_program := program_uid[which.max(abs(lineage_loading))],
      by = cell_type
    ]
    geometry_work <- geometry_work[program_uid != remove_program]
  }

  lineage_scores <- matrix(
    NA_real_,
    nrow = uniqueN(geometry$cell_type),
    ncol = ncol(z),
    dimnames = list(sort(unique(geometry$cell_type)), colnames(z))
  )
  lineage_rows <- list()
  for (lineage in sort(unique(geometry$cell_type))) {
    all_lineage <- geometry[cell_type == lineage]
    selected <- geometry_work[cell_type == lineage]
    tested_uids <- program_testability[testable == TRUE, program_uid]
    selected <- selected[program_uid %in% tested_uids]
    retained_loading <- sum(abs(selected$lineage_loading))
    loading_denominator <- sum(
      abs(geometry_work[cell_type == lineage, lineage_loading])
    )
    n_testable <- nrow(selected)
    n_total <- nrow(all_lineage)
    coverage <- if (loading_denominator > 0) retained_loading / loading_denominator else 0
    testable <- n_testable >= ceiling(0.75 * n_total) && coverage >= 0.80
    lineage_rows[[lineage]] <- data.table(
      cell_type = lineage,
      mode = mode,
      n_programs_total = n_total,
      n_programs_testable = n_testable,
      program_fraction = n_testable / n_total,
      retained_absolute_loading = coverage,
      testable = testable
    )
    if (!testable) next
    scores <- program_scores[selected$program_uid, , drop = FALSE]
    lineage_scores[lineage, ] <- colSums(scores * selected$lineage_loading) /
      sum(abs(selected$lineage_loading))
  }
  list(
    program_scores = program_scores,
    program_testability = program_testability,
    lineage_scores = lineage_scores,
    lineage_testability = rbindlist(lineage_rows)
  )
}

matrix_to_long <- function(matrix_value, row_name, column_name, value_name) {
  output <- as.data.table(as.table(matrix_value))
  setnames(output, c("V1", "V2", "N"), c(row_name, column_name, value_name))
  output
}
