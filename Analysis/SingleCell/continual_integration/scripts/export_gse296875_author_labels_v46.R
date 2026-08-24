#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(SeuratObject))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(paste(
    "usage: export_gse296875_author_labels_v46.R",
    "SOURCE_RDS LABEL_COLUMN EXPECTED_CELLS EXPECTED_LABELS EXPECTED_ROUTING_SHA256 ROUTING_TOKEN OUTPUT_TSV"
  ))
}

source <- normalizePath(args[[1L]], mustWork = TRUE)
label_column <- args[[2L]]
expected_cells <- as.integer(args[[3L]])
expected_labels <- as.integer(args[[4L]])
expected_routing_sha256 <- args[[5L]]
routing_token <- normalizePath(args[[6L]], mustWork = TRUE)
output <- normalizePath(args[[7L]], mustWork = FALSE)
if (file.exists(output)) {
  stop("refusing to overwrite GSE296875 author-label export")
}
if (!dir.exists(dirname(output))) {
  stop("author-label output directory does not exist")
}

token <- readLines(routing_token, warn = FALSE)
if (length(token) != 1L || token[[1L]] != expected_routing_sha256) {
  stop("GSE296875 routing-completion token differs")
}

object <- readRDS(source)
if (!inherits(object, "Seurat") || !(label_column %in% colnames(object@meta.data))) {
  stop("locked GSE296875 author-label column is absent")
}
cell_id <- as.character(rownames(object@meta.data))
author_label <- as.character(object@meta.data[[label_column]])
if (
  length(cell_id) != expected_cells ||
  length(unique(cell_id)) != expected_cells ||
  length(unique(author_label)) != expected_labels ||
  anyNA(author_label) || any(author_label == "")
) {
  stop("GSE296875 author-label census differs")
}

connection <- file(output, open = "wx")
on.exit(close(connection), add = TRUE)
write.table(
  data.frame(cell_id = cell_id, author_label = author_label),
  file = connection,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
