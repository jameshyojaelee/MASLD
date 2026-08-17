#!/usr/bin/env Rscript
# Per-dataset decontX correction for every frozen program lineage.

suppressMessages({
  library(celda)
  library(jsonlite)
  library(Matrix)
})

seed <- 20260815L
set.seed(seed)

project_root <- Sys.getenv("MASLD_PROJECT_ROOT")
candidate_root <- Sys.getenv("CAND_ROOT")
dataset <- Sys.getenv("DECONTX_ONLY_DATASET")
stopifnot(nzchar(project_root), nzchar(candidate_root), nzchar(dataset))

upstream_root <- file.path(
  project_root,
  "Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z"
)
upstream_work <- file.path(upstream_root, "work")
candidate_work <- file.path(candidate_root, "work")
result_root <- file.path(candidate_root, "results", "decontx")
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)

manifest <- fromJSON(file.path(candidate_work, "manifest.json"), simplifyVector = FALSE)
stopifnot(dataset %in% names(manifest$datasets))
genes_all <- readLines(file.path(upstream_work, "genes.txt"))
score_genes <- readLines(file.path(candidate_work, "score_genes.txt"))
score_index <- match(score_genes, genes_all)
stopifnot(!anyNA(score_index))

prefix <- file.path(result_root, dataset)
principal_outputs <- c(
  paste0(prefix, "__run_log.tsv"),
  paste0(prefix, "__ambient_by_gene_lineage.tsv.gz"),
  paste0(prefix, "__contamination_per_program_cell.tsv.gz")
)
candidate_ds <- file.path(candidate_work, dataset)
binary_outputs <- c(
  file.path(candidate_ds, "dec_data.bin"),
  file.path(candidate_ds, "dec_indices.bin"),
  file.path(candidate_ds, "dec_indptr.bin"),
  file.path(candidate_ds, "raw_data.bin"),
  file.path(candidate_ds, "raw_indices.bin"),
  file.path(candidate_ds, "raw_indptr.bin"),
  file.path(candidate_ds, "dec_cells.csv.gz"),
  file.path(candidate_ds, "dec_dims.json")
)
if (any(file.exists(c(principal_outputs, binary_outputs)))) {
  stop(sprintf("Refusing to overwrite an existing %s decontX artifact", dataset))
}

read_matrix <- function(dsdir, info, ngenes) {
  matrices <- list()
  for (chunk in info$chunks) {
    chunk_id <- chunk$chunk
    nnz <- chunk$nnz
    ncells <- chunk$ncells
    i <- readBin(
      file.path(dsdir, sprintf("chunk%d_indices.bin", chunk_id)),
      what = integer(), n = nnz, size = 4, endian = "little", signed = TRUE
    )
    p <- readBin(
      file.path(dsdir, sprintf("chunk%d_indptr.bin", chunk_id)),
      what = integer(), n = ncells + 1L, size = 4, endian = "little", signed = TRUE
    )
    x <- as.double(readBin(
      file.path(dsdir, sprintf("chunk%d_data.bin", chunk_id)),
      what = integer(), n = nnz, size = 4, endian = "little", signed = TRUE
    ))
    stopifnot(length(i) == nnz, length(p) == ncells + 1L, length(x) == nnz)
    matrices[[length(matrices) + 1L]] <- new(
      "dgCMatrix", i = i, p = p, x = x,
      Dim = c(as.integer(ngenes), as.integer(ncells))
    )
  }
  if (length(matrices) == 1L) matrices[[1L]] else do.call(cbind, matrices)
}

group_sums <- function(matrix, group) {
  levels <- sort(unique(group))
  indicator <- sparseMatrix(
    i = match(group, levels), j = seq_along(group), x = 1,
    dims = c(length(levels), length(group))
  )
  result <- as.matrix(matrix %*% t(indicator))
  colnames(result) <- levels
  result
}

upstream_ds <- file.path(upstream_work, dataset)
info <- fromJSON(file.path(upstream_ds, "dims.json"), simplifyVector = FALSE)
meta <- read.csv(
  file.path(candidate_ds, "meta.csv.gz"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
cat(sprintf("[decontX] %s: reconstructing %d cells\n", dataset, info$ncells))
raw <- read_matrix(upstream_ds, info, length(genes_all))
stopifnot(ncol(raw) == nrow(meta), nrow(raw) == length(genes_all))
stopifnot(abs(sum(raw@x) - as.numeric(info$total_count_sum)) <=
            1e-3 * max(1, abs(as.numeric(info$total_count_sum))))

minimum_cluster_cells <- 10L
cluster_counts <- table(meta$cell_type)
small_clusters <- names(cluster_counts)[cluster_counts < minimum_cluster_cells]
cluster_used <- meta$cell_type
cluster_used[cluster_used %in% small_clusters] <- "__pooled_small__"
z <- as.integer(factor(cluster_used))

minimum_batch_cells <- 100L
batch_counts <- table(meta$sample)
small_batches <- names(batch_counts)[batch_counts < minimum_batch_cells]
batch_used <- meta$sample
batch_used[batch_used %in% small_batches] <- "__pooled_small_samples__"

set.seed(seed)
start_time <- Sys.time()
fit <- tryCatch(
  celda::decontX(
    x = raw,
    z = z,
    batch = batch_used,
    seed = seed,
    verbose = TRUE
  ),
  error = function(error) {
    cat(sprintf("[decontX] %s failed: %s\n", dataset, conditionMessage(error)))
    NULL
  }
)
elapsed_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
if (is.null(fit)) {
  corrected <- raw
  contamination <- rep(NA_real_, ncol(raw))
  correction_status <- "uncorrected_passthrough"
  failure_reason <- "decontX_error"
} else {
  corrected <- fit$decontXcounts
  contamination <- as.numeric(fit$contamination)
  correction_status <- "corrected"
  failure_reason <- ""
}

native_labels <- c(
  "Cholangiocytes", "Fibroblasts", "Hepatocytes", "Macrophages", "T cells"
)
report_selection <- meta$cell_type %in% native_labels
report_group <- meta$cell_type[report_selection]
raw_report <- raw[score_index, report_selection, drop = FALSE]
corrected_report <- corrected[score_index, report_selection, drop = FALSE]
raw_sums <- group_sums(raw_report, report_group)
corrected_sums <- group_sums(corrected_report, report_group)
ambient_rows <- list()
for (column in seq_len(ncol(raw_sums))) {
  raw_sum <- raw_sums[, column]
  corrected_sum <- corrected_sums[, column]
  ambient_rows[[length(ambient_rows) + 1L]] <- data.frame(
    dataset = dataset,
    correction_status = correction_status,
    failure_reason = failure_reason,
    cell_type = colnames(raw_sums)[column],
    gene = score_genes,
    n_cells = sum(report_group == colnames(raw_sums)[column]),
    sum_raw = raw_sum,
    sum_corrected = corrected_sum,
    ambient_fraction = ifelse(raw_sum > 0, 1 - corrected_sum / raw_sum, NA_real_),
    stringsAsFactors = FALSE
  )
}
ambient_table <- do.call(rbind, ambient_rows)
ambient_connection <- gzfile(paste0(prefix, "__ambient_by_gene_lineage.tsv.gz"), "w")
write.table(
  ambient_table, ambient_connection, sep = "\t", quote = FALSE, row.names = FALSE
)
close(ambient_connection)

program_cells <- which(nzchar(meta$program_lineage))
stopifnot(length(program_cells) > 0)
cell_table <- meta[program_cells, c(
  "cell_id", "sample", "cell_type", "program_lineage"
)]
cell_table$dataset <- dataset
cell_table$contamination <- contamination[program_cells]
cell_connection <- gzfile(
  paste0(prefix, "__contamination_per_program_cell.tsv.gz"), "w"
)
write.table(cell_table, cell_connection, sep = "\t", quote = FALSE, row.names = FALSE)
close(cell_connection)

corrected_score <- as(corrected[score_index, program_cells, drop = FALSE], "dgCMatrix")
raw_score <- as(raw[score_index, program_cells, drop = FALSE], "dgCMatrix")
writeBin(as.integer(corrected_score@i), file.path(candidate_ds, "dec_indices.bin"),
         size = 4, endian = "little")
writeBin(as.integer(corrected_score@p), file.path(candidate_ds, "dec_indptr.bin"),
         size = 4, endian = "little")
writeBin(as.double(corrected_score@x), file.path(candidate_ds, "dec_data.bin"),
         size = 8, endian = "little")
writeBin(as.integer(raw_score@i), file.path(candidate_ds, "raw_indices.bin"),
         size = 4, endian = "little")
writeBin(as.integer(raw_score@p), file.path(candidate_ds, "raw_indptr.bin"),
         size = 4, endian = "little")
writeBin(as.double(raw_score@x), file.path(candidate_ds, "raw_data.bin"),
         size = 8, endian = "little")
cell_index_connection <- gzfile(file.path(candidate_ds, "dec_cells.csv.gz"), "w")
write.csv(
  cell_table[, c("cell_id", "sample", "dataset", "program_lineage")],
  cell_index_connection,
  row.names = FALSE
)
close(cell_index_connection)

dims <- list(
  dataset = dataset,
  correction_status = correction_status,
  failure_reason = failure_reason,
  n_genes_score = length(score_genes),
  n_cells = length(program_cells),
  corrected_nnz = length(corrected_score@x),
  raw_nnz = length(raw_score@x),
  corrected_sum = sum(corrected_score@x),
  raw_sum = sum(raw_score@x),
  lineage_cells = as.list(table(meta$program_lineage[program_cells]))
)
write(
  toJSON(dims, auto_unbox = TRUE, digits = 16, pretty = TRUE),
  file.path(candidate_ds, "dec_dims.json")
)

finite_contamination <- contamination[is.finite(contamination)]
summary_or_na <- function(function_name, values, ...) {
  if (!length(values)) return(NA_real_)
  do.call(function_name, c(list(values), list(...)))
}
run_log <- data.frame(
  dataset = dataset,
  correction_status = correction_status,
  failure_reason = failure_reason,
  n_cells = ncol(raw),
  n_program_universe_cells = length(program_cells),
  n_program_lineages = length(unique(meta$program_lineage[program_cells])),
  n_clusters = length(unique(z)),
  pooled_small_types = paste(small_clusters, collapse = ";"),
  n_ambient_batches = length(unique(batch_used)),
  n_pooled_small_samples = length(small_batches),
  contamination_median = summary_or_na("median", finite_contamination),
  contamination_mean = summary_or_na("mean", finite_contamination),
  contamination_q10 = summary_or_na("quantile", finite_contamination, probs = 0.10),
  contamination_q90 = summary_or_na("quantile", finite_contamination, probs = 0.90),
  decontx_minutes = elapsed_minutes,
  seed = seed,
  stringsAsFactors = FALSE
)
write.table(run_log, paste0(prefix, "__run_log.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
writeLines(capture.output(sessionInfo()), paste0(prefix, "__sessionInfo.txt"))
cat(sprintf(
  "[decontX] %s complete: status=%s; score matrix=%d genes x %d cells\n",
  dataset, correction_status, length(score_genes), length(program_cells)
))
