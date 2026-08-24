#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(SeuratObject))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: inspect_gse296875_metadata_structure_v44a.R SOURCE_RDS OUTPUT_TSV")
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
object_class <- as.character(attr(object, "class"))
if (!isS4(object) || length(object_class) != 1L || object_class[[1L]] != "Seurat") {
  stop("GSE296875 processed object is not a Seurat S4 object")
}

metadata <- slot(object, "meta.data")
active_ident <- slot(object, "active.ident")
if (nrow(metadata) != length(active_ident)) {
  stop("GSE296875 top-level cell dimensions differ")
}
rows <- list(
  data.frame(kind = "object", name = "class", value = "Seurat"),
  data.frame(kind = "object", name = "cells", value = as.character(nrow(metadata)))
)
for (index in seq_along(colnames(metadata))) {
  rows[[length(rows) + 1L]] <- data.frame(
    kind = "metadata_column",
    name = as.character(index),
    value = colnames(metadata)[[index]]
  )
}
for (assay_name in names(slot(object, "assays"))) {
  rows[[length(rows) + 1L]] <- data.frame(
    kind = "assay_name",
    name = assay_name,
    value = assay_name
  )
}
for (reduction_name in names(slot(object, "reductions"))) {
  rows[[length(rows) + 1L]] <- data.frame(
    kind = "reduction_name",
    name = reduction_name,
    value = reduction_name
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
