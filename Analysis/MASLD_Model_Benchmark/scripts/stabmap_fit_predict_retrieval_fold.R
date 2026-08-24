#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(StabMap)
})

fail <- function(message) stop(message, call. = FALSE)

read_ids <- function(path, expected_fields) {
  value <- read.delim(
    path, header = TRUE, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE
  )
  if (!identical(colnames(value), expected_fields)) {
    fail(paste("identifier schema differs:", path))
  }
  value
}

read_counts <- function(path, feature_ids, cell_ids) {
  value <- as(readMM(path), "dgCMatrix")
  if (!identical(dim(value), c(length(feature_ids), length(cell_ids)))) {
    fail(paste("count matrix axes differ:", path))
  }
  if (length(value@x) == 0L || any(!is.finite(value@x)) || any(value@x < 0)) {
    fail(paste("count matrix values differ:", path))
  }
  rownames(value) <- feature_ids
  colnames(value) <- cell_ids
  value
}

log_cpm <- function(counts) {
  depth <- Matrix::colSums(counts)
  if (any(!is.finite(depth)) || any(depth <= 0)) fail("RNA depth differs")
  value <- counts %*% Diagonal(x = 10000 / depth)
  value@x <- log1p(value@x)
  dimnames(value) <- dimnames(counts)
  value
}

fit_atac_tfidf <- function(bridge_counts, query_counts) {
  bridge_binary <- bridge_counts
  query_binary <- query_counts
  bridge_binary@x[] <- 1
  query_binary@x[] <- 1
  bridge_depth <- Matrix::colSums(bridge_binary)
  query_depth <- Matrix::colSums(query_binary)
  if (any(bridge_depth <= 0) || any(query_depth <= 0)) {
    fail("ATAC detected-feature depth differs")
  }
  detection <- Matrix::rowSums(bridge_binary)
  idf <- log1p(ncol(bridge_binary) / pmax(detection, 1))
  bridge <- Diagonal(x = idf) %*% bridge_binary %*%
    Diagonal(x = 10000 / bridge_depth)
  query <- Diagonal(x = idf) %*% query_binary %*%
    Diagonal(x = 10000 / query_depth)
  dimnames(bridge) <- dimnames(bridge_counts)
  dimnames(query) <- dimnames(query_counts)
  list(
    bridge = bridge,
    query = query,
    idf = idf
  )
}

prefix_rows <- function(value, prefix) {
  rownames(value) <- paste0(prefix, rownames(value))
  value
}

prefix_columns <- function(value, prefix) {
  if (is.null(colnames(value)) || length(colnames(value)) != ncol(value)) {
    fail(paste("column names missing before prefix", prefix))
  }
  colnames(value) <- paste0(prefix, colnames(value))
  value
}

select_training_features <- function(value, number) {
  variance <- MatrixGenerics::rowVars(value)
  names(variance) <- rownames(value)
  eligible <- which(is.finite(variance) & variance > 0)
  ordered <- eligible[order(-variance[eligible], names(variance)[eligible])]
  if (length(ordered) < number) fail("training-only feature selection is too small")
  names(variance)[ordered[seq_len(number)]]
}

reference_rebase <- function(embedding, reference_ids) {
  if (!all(reference_ids %in% rownames(embedding))) {
    fail("reference embedding identities differ")
  }
  center <- colMeans(embedding[reference_ids, , drop = FALSE])
  sweep(embedding, 2L, center, "-")
}

run_stabmap <- function(reference, bridge, query_rna, query_atac, seed) {
  set.seed(seed)
  assays <- list(
    reference_rna = reference,
    bridge_multiome = bridge,
    query_rna = query_rna,
    query_atac = query_atac
  )
  assays <- lapply(assays, as.matrix)
  value <- StabMap::stabMap(
    assays,
    labels_list = NULL,
    reference_list = c("reference_rna"),
    ncomponentsReference = 16L,
    ncomponentsSubset = 16L,
    suppressMessages = TRUE,
    projectAll = FALSE,
    restrictFeatures = FALSE,
    maxFeatures = 1000L,
    plot = FALSE,
    scale.center = FALSE,
    scale.scale = FALSE,
    BPPARAM = BiocParallel::SerialParam(),
    verbose = FALSE
  )
  if (ncol(value) != 16L || any(!is.finite(value))) {
    fail("StabMap embedding differs")
  }
  reference_rebase(value, colnames(reference))
}

write_embedding <- function(path, identifiers, prefix, embedding) {
  internal_ids <- paste0(prefix, identifiers)
  if (!all(internal_ids %in% rownames(embedding))) {
    fail("query embedding identities differ")
  }
  value <- embedding[internal_ids, , drop = FALSE]
  if (any(!is.finite(value)) || any(sqrt(rowSums(value ^ 2)) <= 0)) {
    fail("query embedding values differ")
  }
  rownames(value) <- identifiers
  output <- data.frame(
    query_id = identifiers, value, check.names = FALSE, stringsAsFactors = FALSE
  )
  write.table(
    output, path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE
  )
}

run_fold <- function(input, output, fold, seed) {
  genes <- read_ids(file.path(input, "genes.tsv"), c("feature_index", "gene_id"))
  peaks <- read_ids(file.path(input, "peaks.tsv"), c("feature_index", "peak_id"))
  training <- read_ids(
    file.path(input, "training_rows.tsv"),
    c("training_id", "rna_state", "atac_state")
  )
  query_rna_rows <- read_ids(
    file.path(input, "query_rna_rows.tsv"), c("rna_query_id", "rna_state")
  )
  query_atac_rows <- read_ids(
    file.path(input, "query_atac_rows.tsv"), c("atac_query_id", "atac_state")
  )
  if (any(training$rna_state != "observed") ||
      any(training$atac_state != "observed") ||
      any(query_rna_rows$rna_state != "observed") ||
      any(query_atac_rows$atac_state != "observed")) {
    fail("observed-state contract differs")
  }
  training_ids <- training$training_id
  rna_ids <- query_rna_rows$rna_query_id
  atac_ids <- query_atac_rows$atac_query_id
  if (length(intersect(rna_ids, atac_ids)) != 0L) {
    fail("held RNA and ATAC identities are not blinded")
  }
  training_rna <- read_counts(
    file.path(input, "training_rna.mtx"), genes$gene_id, training_ids
  )
  training_atac <- read_counts(
    file.path(input, "training_atac.mtx"), peaks$peak_id, training_ids
  )
  query_rna_counts <- read_counts(
    file.path(input, "query_rna.mtx"), genes$gene_id, rna_ids
  )
  query_atac_counts <- read_counts(
    file.path(input, "query_atac.mtx"), peaks$peak_id, atac_ids
  )

  ordered <- order(training_ids)
  reference_indices <- ordered[seq.int(1L, length(ordered), by = 2L)]
  bridge_indices <- ordered[seq.int(2L, length(ordered), by = 2L)]
  if (length(reference_indices) < 100L || length(bridge_indices) < 100L) {
    fail("training mosaic split is too small")
  }
  reference_rna <- prefix_columns(
    prefix_rows(log_cpm(training_rna[, reference_indices, drop = FALSE]), "rna:"),
    "reference::"
  )
  bridge_rna <- prefix_rows(
    log_cpm(training_rna[, bridge_indices, drop = FALSE]), "rna:"
  )
  query_rna <- prefix_columns(
    prefix_rows(log_cpm(query_rna_counts), "rna:"), "rna_query::"
  )
  atac <- fit_atac_tfidf(
    training_atac[, bridge_indices, drop = FALSE], query_atac_counts
  )
  bridge_atac <- prefix_rows(atac$bridge, "atac:")
  query_atac <- prefix_columns(prefix_rows(atac$query, "atac:"), "atac_query::")
  rna_features <- select_training_features(reference_rna, 1000L)
  atac_features <- select_training_features(bridge_atac, 1000L)
  reference_rna <- reference_rna[rna_features, , drop = FALSE]
  bridge_rna <- bridge_rna[rna_features, , drop = FALSE]
  query_rna <- query_rna[rna_features, , drop = FALSE]
  bridge_atac <- bridge_atac[atac_features, , drop = FALSE]
  query_atac <- query_atac[atac_features, , drop = FALSE]
  bridge <- prefix_columns(rbind(bridge_rna, bridge_atac), "bridge::")

  dir.create(output, recursive = TRUE, showWarnings = FALSE, mode = "0750")
  state_path <- file.path(output, "training_state.rds")
  saveRDS(
    list(
      reference_rna = reference_rna,
      bridge = bridge,
      atac_idf = atac$idf,
      rna_features = rna_features,
      atac_features = atac_features,
      reference_indices = reference_indices,
      bridge_indices = bridge_indices
    ),
    state_path, version = 3
  )
  state <- readRDS(state_path)
  if (!identical(state$reference_rna, reference_rna) ||
      !identical(state$bridge, bridge) || !identical(state$atac_idf, atac$idf)) {
    fail("StabMap training state resume differs")
  }

  embedding <- run_stabmap(
    state$reference_rna, state$bridge, query_rna, query_atac, seed
  )
  reverse_rna <- rev(seq_len(ncol(query_rna)))
  reverse_atac <- rev(seq_len(ncol(query_atac)))
  repeated <- run_stabmap(
    state$reference_rna,
    state$bridge,
    query_rna[, reverse_rna, drop = FALSE],
    query_atac[, reverse_atac, drop = FALSE],
    seed
  )
  common <- intersect(rownames(embedding), rownames(repeated))
  difference <- max(abs(embedding[common, , drop = FALSE] -
                        repeated[common, , drop = FALSE]))
  if (!is.finite(difference) || difference > 1e-8) {
    fail("StabMap query-order invariance differs")
  }
  write_embedding(
    file.path(output, "query_rna_embeddings.tsv"), rna_ids, "rna_query::", embedding
  )
  write_embedding(
    file.path(output, "query_atac_embeddings.tsv"), atac_ids, "atac_query::", embedding
  )
  receipt <- data.frame(
    key = c(
      "status", "model_id", "exact_package", "dataset_id", "fold",
      "seed",
      "training_nuclei", "reference_rna_nuclei", "bridge_multiome_nuclei",
      "query_rna_nuclei", "query_atac_nuclei", "genes", "peaks",
      "selected_rna_features", "selected_atac_features",
      "embedding_dimensions", "training_state_resume", "query_order_invariant",
      "max_query_order_difference", "hidden_pair_map_read",
      "retrieval_metrics_calculated", "outcomes_read", "test_outcomes_read"
    ),
    value = c(
      "pass", "stabmap", "StabMap 1.6.0", "gse296875", as.character(fold),
      as.character(seed),
      length(training_ids), length(reference_indices), length(bridge_indices),
      length(rna_ids), length(atac_ids), nrow(genes), nrow(peaks), "1000", "1000",
      "16", "true",
      "true", format(difference, scientific = TRUE, digits = 17), "false", "false",
      "false", "false"
    ),
    stringsAsFactors = FALSE
  )
  write.table(
    receipt, file.path(output, "receipt.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE
  )
  cat(paste(receipt$key, receipt$value, sep = "\t"), sep = "\n")
}

self_test <- function() {
  counts <- Matrix(c(1, 0, 2, 3, 1, 0, 0, 2, 1, 1, 0, 1), nrow = 3, sparse = TRUE)
  normalized <- log_cpm(counts)
  stopifnot(identical(dim(normalized), dim(counts)))
  tfidf <- fit_atac_tfidf(counts, counts)
  stopifnot(identical(dim(tfidf$bridge), dim(counts)))
  cat("status\tpass\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--self-test")) {
  self_test()
  quit(save = "no", status = 0L)
}
if (!(length(args) %in% c(3L, 4L))) {
  fail("usage: stabmap_fit_predict_retrieval_fold.R INPUT OUTPUT FOLD [BASE_SEED]")
}
input <- normalizePath(args[[1L]], mustWork = TRUE)
output <- args[[2L]]
fold <- as.integer(args[[3L]])
if (is.na(fold) || fold < 0L || fold > 4L) fail("fold differs")
base_seed <- if (length(args) == 4L) as.integer(args[[4L]]) else 4817L
if (is.na(base_seed) || base_seed < 0L) fail("base seed differs")
if (file.exists(output)) fail("output exists")
run_fold(input, output, fold, base_seed + fold * 101L)
