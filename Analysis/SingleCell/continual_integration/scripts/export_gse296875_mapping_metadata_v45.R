#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(SeuratObject))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop(paste(
    "usage: export_gse296875_mapping_metadata_v45.R",
    "SOURCE_RDS DONOR_COLUMN WELL_COLUMN BARCODE_COLUMN EXPECTED_CELLS EXPECTED_DONORS EXPECTED_WELLS OUTPUT_TSV SCRATCH_CACHE_RDS"
  ))
}

source <- normalizePath(args[[1L]], mustWork = TRUE)
donor_column <- args[[2L]]
well_column <- args[[3L]]
barcode_column <- args[[4L]]
expected_cells <- as.integer(args[[5L]])
expected_donors <- as.integer(args[[6L]])
expected_wells <- as.integer(args[[7L]])
output <- normalizePath(args[[8L]], mustWork = FALSE)
scratch_cache <- normalizePath(args[[9L]], mustWork = FALSE)
if (file.exists(output)) {
  stop("refusing to overwrite GSE296875 mapping metadata")
}
if (!dir.exists(dirname(output))) {
  stop("mapping-metadata output directory does not exist")
}
if (file.exists(scratch_cache) || !dir.exists(dirname(scratch_cache))) {
  stop("scratch cache exists or its parent directory is absent")
}

object <- readRDS(source)
if (!inherits(object, "Seurat")) {
  stop("GSE296875 processed object is not a Seurat object")
}
if (!all(c(donor_column, well_column, barcode_column) %in% colnames(object@meta.data))) {
  stop("locked GSE296875 mapping metadata columns are absent")
}

cell_id <- as.character(rownames(object@meta.data))
donor_id <- as.character(object@meta.data[[donor_column]])
well_id <- as.character(object@meta.data[[well_column]])
raw_barcode <- as.character(object@meta.data[[barcode_column]])
if (
  length(cell_id) != expected_cells ||
  length(unique(cell_id)) != expected_cells ||
  length(unique(donor_id)) != expected_donors ||
  length(unique(well_id)) != expected_wells ||
  anyNA(donor_id) || anyNA(well_id) ||
  anyNA(raw_barcode) || any(donor_id == "") || any(well_id == "") ||
  any(raw_barcode == "") || anyDuplicated(paste(well_id, raw_barcode, sep = "\r"))
) {
  stop("GSE296875 mapping metadata census or identifiers differ")
}

connection <- file(output, open = "wx")
on.exit(close(connection), add = TRUE)
write.table(
  data.frame(
    cell_id = cell_id,
    donor_id = donor_id,
    well_id = well_id,
    raw_barcode = raw_barcode
  ),
  file = connection,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
saveRDS(object, file = scratch_cache, compress = FALSE)
