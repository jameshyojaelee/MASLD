#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(SeuratObject))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: inspect_gse296875_structure_v44.R SOURCE_RDS OUTPUT_TSV")
}

source <- normalizePath(args[[1L]], mustWork = TRUE)
output <- normalizePath(args[[2L]], mustWork = FALSE)
if (file.exists(output)) {
  stop("refusing to overwrite structural inspection output")
}
if (!dir.exists(dirname(output))) {
  stop("structural inspection output directory does not exist")
}

object <- readRDS(source)
if (!inherits(object, "Seurat")) {
  stop("GSE296875 processed object is not a Seurat object")
}

rows <- list(
  data.frame(kind = "object", name = "class", value = paste(class(object), collapse = ";")),
  data.frame(kind = "object", name = "cells", value = as.character(ncol(object))),
  data.frame(kind = "object", name = "features_active_assay", value = as.character(nrow(object)))
)

for (index in seq_along(colnames(object@meta.data))) {
  rows[[length(rows) + 1L]] <- data.frame(
    kind = "metadata_column",
    name = as.character(index),
    value = colnames(object@meta.data)[[index]]
  )
}
for (assay_name in names(object@assays)) {
  assay <- object@assays[[assay_name]]
  rows[[length(rows) + 1L]] <- data.frame(
    kind = "assay",
    name = assay_name,
    value = paste(class(assay), collapse = ";")
  )
  rows[[length(rows) + 1L]] <- data.frame(
    kind = "assay_dimensions",
    name = assay_name,
    value = paste(nrow(assay), ncol(assay), sep = "x")
  )
  for (layer_name in Layers(assay)) {
    rows[[length(rows) + 1L]] <- data.frame(
      kind = "assay_layer",
      name = assay_name,
      value = layer_name
    )
  }
}
for (reduction_name in names(object@reductions)) {
  reduction <- object@reductions[[reduction_name]]
  rows[[length(rows) + 1L]] <- data.frame(
    kind = "reduction",
    name = reduction_name,
    value = paste(nrow(reduction), ncol(reduction), sep = "x")
  )
}

result <- do.call(rbind, rows)
connection <- file(output, open = "wx")
on.exit(close(connection), add = TRUE)
write.table(
  result,
  file = connection,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
