#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

fail <- function(message) {
  stop(message, call. = FALSE)
}

read_ids <- function(path, expected_fields) {
  value <- read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
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

log_normalize <- function(counts) {
  depth <- Matrix::colSums(counts)
  if (any(!is.finite(depth)) || any(depth <= 0)) {
    fail("RNA depth differs")
  }
  normalized <- counts %*% Diagonal(x = 10000 / depth)
  normalized@x <- log1p(normalized@x)
  normalized
}

fit_rna_projection <- function(training, query, dimensions = 30L) {
  training_norm <- log_normalize(training)
  query_norm <- log_normalize(query)
  training_dense <- t(as.matrix(training_norm))
  query_dense <- t(as.matrix(query_norm))
  center <- colMeans(training_dense)
  scale <- apply(training_dense, 2L, sd)
  scale[!is.finite(scale) | scale <= 0] <- 1
  training_scaled <- sweep(training_dense, 2L, center, "-")
  training_scaled <- sweep(training_scaled, 2L, scale, "/")
  query_scaled <- sweep(query_dense, 2L, center, "-")
  query_scaled <- sweep(query_scaled, 2L, scale, "/")
  fit <- irlba::prcomp_irlba(
    training_scaled,
    n = dimensions,
    center = FALSE,
    scale. = FALSE
  )
  training_embedding <- fit$x[, seq_len(dimensions), drop = FALSE]
  query_embedding <- query_scaled %*% fit$rotation[, seq_len(dimensions), drop = FALSE]
  rownames(training_embedding) <- colnames(training)
  rownames(query_embedding) <- colnames(query)
  colnames(training_embedding) <- paste0("PC_", seq_len(dimensions))
  colnames(query_embedding) <- paste0("PC_", seq_len(dimensions))
  list(training = training_embedding, query = query_embedding)
}

fit_atac_projection <- function(training, query, dimensions = 30L) {
  training_depth <- Matrix::colSums(training)
  query_depth <- Matrix::colSums(query)
  if (any(training_depth <= 0) || any(query_depth <= 0)) {
    fail("ATAC depth differs")
  }
  detection <- Matrix::rowSums(training > 0)
  inverse <- log1p(ncol(training) / pmax(detection, 1))
  training_tfidf <- Diagonal(x = inverse) %*% training %*%
    Diagonal(x = 10000 / training_depth)
  query_tfidf <- Diagonal(x = inverse) %*% query %*%
    Diagonal(x = 10000 / query_depth)
  fit <- irlba::irlba(t(training_tfidf), nv = dimensions, nu = dimensions)
  training_embedding <- fit$u[, seq_len(dimensions), drop = FALSE] %*%
    diag(fit$d[seq_len(dimensions)], nrow = dimensions)
  query_embedding <- t(query_tfidf) %*% fit$v[, seq_len(dimensions), drop = FALSE]
  rownames(training_embedding) <- colnames(training)
  rownames(query_embedding) <- colnames(query)
  colnames(training_embedding) <- paste0("LSI_", seq_len(dimensions))
  colnames(query_embedding) <- paste0("LSI_", seq_len(dimensions))
  list(training = training_embedding, query = query_embedding)
}

make_reduction <- function(embedding, assay, key, projected = FALSE) {
  embedding <- as.matrix(embedding)
  misc <- if (projected) list(ref.dims = seq_len(ncol(embedding))) else list()
  CreateDimReducObject(
    embeddings = embedding, assay = assay, key = key, misc = misc
  )
}

map_queries <- function(bridge, query_rna_counts, query_atac_counts,
                        query_rna_embedding, query_atac_embedding) {
  query_rna <- CreateSeuratObject(counts = query_rna_counts, assay = "RNA")
  query_rna[["pca"]] <- make_reduction(
    query_rna_embedding, "RNA", "PC_", projected = TRUE
  )
  query_atac <- CreateSeuratObject(counts = query_atac_counts, assay = "ATAC")
  query_atac[["lsi"]] <- make_reduction(
    query_atac_embedding, "ATAC", "LSI_", projected = TRUE
  )
  bridge_cells_representation <- getFromNamespace(
    "BridgeCellsRepresentation", "Seurat"
  )
  mapped <- bridge_cells_representation(
    object.list = list(query_rna, query_atac),
    bridge.object = bridge,
    object.reduction = c("pca", "lsi"),
    bridge.reduction = c("pca", "lsi"),
    laplacian.reduction = "lap",
    laplacian.dims = 1:16,
    bridge.assay.name = "Bridge",
    verbose = FALSE
  )
  extract <- function(object) {
    value <- t(as.matrix(GetAssayData(object, assay = "Bridge", layer = "data")))
    if (ncol(value) != 16L || any(!is.finite(value))) {
      fail("Seurat WNN bridge query embedding differs")
    }
    norms <- sqrt(rowSums(value ^ 2))
    if (any(norms <= 0)) {
      fail("Seurat WNN bridge query embedding is zero")
    }
    value
  }
  list(rna = extract(mapped[[1L]]), atac = extract(mapped[[2L]]))
}

write_embedding <- function(path, identifiers, embedding) {
  if (!identical(rownames(embedding), identifiers)) {
    fail("query embedding identity order differs")
  }
  value <- data.frame(
    query_id = identifiers,
    embedding,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write.table(
    value,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

run_fold <- function(input, output, fold, seed) {
  gene_rows <- read_ids(file.path(input, "genes.tsv"), c("feature_index", "gene_id"))
  peak_rows <- read_ids(file.path(input, "peaks.tsv"), c("feature_index", "peak_id"))
  training_rows <- read_ids(
    file.path(input, "training_rows.tsv"),
    c("training_id", "rna_state", "atac_state")
  )
  rna_rows <- read_ids(
    file.path(input, "query_rna_rows.tsv"), c("rna_query_id", "rna_state")
  )
  atac_rows <- read_ids(
    file.path(input, "query_atac_rows.tsv"), c("atac_query_id", "atac_state")
  )
  training_ids <- training_rows$training_id
  rna_ids <- rna_rows$rna_query_id
  atac_ids <- atac_rows$atac_query_id
  if (length(intersect(rna_ids, atac_ids)) != 0L) {
    fail("held query IDs are not disjoint")
  }
  training_rna <- read_counts(
    file.path(input, "training_rna.mtx"), gene_rows$gene_id, training_ids
  )
  training_atac <- read_counts(
    file.path(input, "training_atac.mtx"), peak_rows$peak_id, training_ids
  )
  query_rna <- read_counts(
    file.path(input, "query_rna.mtx"), gene_rows$gene_id, rna_ids
  )
  query_atac <- read_counts(
    file.path(input, "query_atac.mtx"), peak_rows$peak_id, atac_ids
  )
  set.seed(seed)
  rna <- fit_rna_projection(training_rna, query_rna, dimensions = 30L)
  set.seed(seed + 1L)
  atac <- fit_atac_projection(training_atac, query_atac, dimensions = 30L)
  bridge <- CreateSeuratObject(counts = training_rna, assay = "RNA")
  bridge[["ATAC"]] <- CreateAssayObject(counts = training_atac)
  bridge[["pca"]] <- make_reduction(rna$training, "RNA", "PC_")
  bridge[["lsi"]] <- make_reduction(atac$training, "ATAC", "LSI_")
  set.seed(seed + 2L)
  bridge <- FindMultiModalNeighbors(
    bridge,
    reduction.list = list("pca", "lsi"),
    dims.list = list(1:30, 1:30),
    k.nn = min(20L, ncol(bridge) - 1L),
    knn.range = min(100L, ncol(bridge) - 1L),
    knn.graph.name = "wknn",
    snn.graph.name = "wsnn",
    weighted.nn.name = "weighted.nn",
    verbose = FALSE
  )
  set.seed(seed + 3L)
  bridge <- RunGraphLaplacian(
    bridge,
    graph = "wsnn",
    reduction.name = "lap",
    reduction.key = "lap_",
    verbose = FALSE
  )
  embeddings <- map_queries(
    bridge, query_rna, query_atac, rna$query, atac$query
  )
  dir.create(output, recursive = TRUE, showWarnings = FALSE, mode = "0750")
  state_path <- file.path(output, "training_bridge.rds")
  saveRDS(bridge, state_path, version = 3)
  resumed <- readRDS(state_path)
  if (
    !identical(Embeddings(bridge, "pca"), Embeddings(resumed, "pca")) ||
    !identical(Embeddings(bridge, "lsi"), Embeddings(resumed, "lsi")) ||
    !identical(Embeddings(bridge, "lap"), Embeddings(resumed, "lap"))
  ) {
    fail("Seurat WNN bridge strict state resume differs")
  }
  repeated <- map_queries(
    resumed, query_rna, query_atac, rna$query, atac$query
  )
  if (!identical(embeddings$rna, repeated$rna) ||
      !identical(embeddings$atac, repeated$atac)) {
    fail("Seurat WNN bridge repeated query mapping differs")
  }
  write_embedding(file.path(output, "query_rna_embeddings.tsv"), rna_ids, embeddings$rna)
  write_embedding(file.path(output, "query_atac_embeddings.tsv"), atac_ids, embeddings$atac)
  list(
    fold = fold,
    training_nuclei = length(training_ids),
    query_rna_nuclei = length(rna_ids),
    query_atac_nuclei = length(atac_ids),
    genes = length(gene_rows$gene_id),
    peaks = length(peak_rows$peak_id),
    pca_dimensions = 30L,
    lsi_dimensions = 30L,
    bridge_laplacian_dimensions = 16L,
    wnn_k = min(20L, length(training_ids) - 1L),
    strict_state_resume = TRUE,
    repeated_query_mapping_bit_identical = TRUE,
    hidden_pair_map_read = FALSE,
    retrieval_metrics_calculated = FALSE
  )
}

self_test <- function() {
  counts <- Matrix(
    c(1, 0, 2, 3, 1, 0, 0, 2, 1, 1, 0, 1),
    nrow = 3,
    sparse = TRUE
  )
  normalized <- log_normalize(counts)
  stopifnot(identical(dim(normalized), dim(counts)))
  stopifnot(all(normalized@x >= 0), all(is.finite(normalized@x)))
  cat("status\tpass\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--self-test")) {
  self_test()
  quit(save = "no", status = 0L)
}
if (length(args) != 2L) {
  fail("usage: seurat_wnn_bridge_fit_predict_smoke.R INPUT OUTPUT")
}
input_root <- normalizePath(args[[1L]], mustWork = TRUE)
output_root <- args[[2L]]
if (file.exists(output_root)) {
  fail("output exists")
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE, mode = "0750")
seed <- 1709L
folds <- list()
for (fold in 0:4) {
  fold_id <- as.character(fold)
  folds[[fold_id]] <- run_fold(
    file.path(input_root, paste0("fold_", fold)),
    file.path(output_root, paste0("fold_", fold), "seurat_wnn_bridge"),
    fold,
    seed + fold * 101L
  )
}
receipt <- list(
  schema_version = "masld-bench-seurat-wnn-bridge-fit-predict-v1",
  status = "pass",
  model_id = "seurat_wnn_bridge",
  exact_package = "Seurat 5.5.1",
  dataset_id = "gse296875",
  pairing_topology = "same_nucleus",
  outer_unit = "donor",
  folds = folds,
  training_pairing_available = TRUE,
  query_modality_ids_disjoint = TRUE,
  query_atac_row_order_permuted = TRUE,
  hidden_pair_map_read = FALSE,
  retrieval_metrics_calculated = FALSE,
  cells_used_as_biological_replicates = FALSE,
  outcomes_read = FALSE,
  test_outcomes_read = FALSE,
  rna_to_atac_prediction_executed = FALSE,
  sealed_rna_conditioned_atac_eligible = FALSE,
  champion_claim_allowed = FALSE
)
jsonlite::write_json(
  receipt,
  file.path(output_root, "receipt.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  null = "null",
  digits = NA
)
cat(jsonlite::toJSON(receipt, auto_unbox = TRUE, digits = NA), "\n")
