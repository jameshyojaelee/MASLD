#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(paste(
    "usage: export_gse296875_atac_smoke.R",
    "SOURCE_RDS SELECTION_TSV OUTPUT_DIRECTORY"
  ))
}

source_rds <- normalizePath(args[[1L]], mustWork = TRUE)
selection_path <- normalizePath(args[[2L]], mustWork = TRUE)
output <- normalizePath(args[[3L]], mustWork = FALSE)
if (dir.exists(output) || file.exists(output)) {
  stop("refusing to overwrite GSE296875 ATAC smoke output")
}
if (!dir.exists(dirname(output))) {
  stop("GSE296875 ATAC smoke output parent does not exist")
}
dir.create(output, mode = "0750")

selection <- read.delim(
  selection_path,
  header = TRUE,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required <- c("cell_id", "donor_id", "well_id", "broad_label")
if (!all(required %in% colnames(selection))) {
  stop("GSE296875 selection lacks required columns")
}
if (
  nrow(selection) != 1000L ||
  anyNA(selection$cell_id) ||
  anyDuplicated(selection$cell_id)
) {
  stop("GSE296875 selection does not contain 1,000 unique nuclei")
}

object <- readRDS(source_rds)
if (!inherits(object, "Seurat")) {
  stop("GSE296875 processed object is not a Seurat object")
}
assays <- slot(object, "assays")
if (!("peaks" %in% names(assays))) {
  stop("GSE296875 processed object lacks the peaks assay")
}
peak_assay <- assays[["peaks"]]
peak_storage <- attributes(peak_assay)
storage_names <- names(peak_storage)
if (is.null(storage_names)) {
  stop("GSE296875 peaks assay has no serialized storage attributes")
}
if ("counts" %in% storage_names) {
  counts <- peak_storage[["counts"]]
} else if ("layers" %in% storage_names) {
  layers <- peak_storage[["layers"]]
  if (!("counts" %in% names(layers))) {
    stop("GSE296875 peaks assay lacks a counts layer")
  }
  counts <- layers[["counts"]]
} else {
  stop("GSE296875 peaks assay has no recognized count storage")
}
if (!inherits(counts, "sparseMatrix")) {
  stop("GSE296875 ATAC counts are not a sparse matrix")
}
if (is.null(rownames(counts)) || is.null(colnames(counts))) {
  stop("GSE296875 ATAC counts lack feature or nucleus identifiers")
}
missing_cells <- setdiff(selection$cell_id, colnames(counts))
if (length(missing_cells) > 0L) {
  stop("selected GSE296875 nuclei are absent from the ATAC assay")
}

parse_peak <- function(value) {
  match <- regexec("^(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y))[:-]([0-9]+)-([0-9]+)$", value)
  pieces <- regmatches(value, match)
  valid <- lengths(pieces) == 4L
  chromosome <- rep(NA_character_, length(value))
  start <- rep(NA_real_, length(value))
  end <- rep(NA_real_, length(value))
  chromosome[valid] <- vapply(pieces[valid], `[[`, character(1L), 2L)
  start[valid] <- as.numeric(vapply(pieces[valid], `[[`, character(1L), 3L))
  end[valid] <- as.numeric(vapply(pieces[valid], `[[`, character(1L), 4L))
  data.frame(
    peak_id = value,
    chromosome = chromosome,
    source_start_1based_closed = start,
    source_end_1based_closed = end,
    valid = valid & is.finite(start) & is.finite(end) & start >= 1 & end >= start,
    stringsAsFactors = FALSE
  )
}

peaks <- parse_peak(rownames(counts))
if (!any(peaks$valid)) {
  stop("GSE296875 processed peak identifiers do not match the primary-contig contract")
}
peaks <- peaks[
  peaks$valid,
  c(
    "peak_id",
    "chromosome",
    "source_start_1based_closed",
    "source_end_1based_closed"
  )
]
if (anyDuplicated(peaks$peak_id)) {
  stop("GSE296875 processed primary peak identifiers are not unique")
}
selected_counts <- counts[peaks$peak_id, selection$cell_id, drop = FALSE]
if (
  nrow(selected_counts) != nrow(peaks) ||
  ncol(selected_counts) != nrow(selection) ||
  any(selected_counts@x < 0) ||
  any(selected_counts@x != floor(selected_counts@x)) ||
  any(Matrix::colSums(selected_counts) <= 0)
) {
  stop("GSE296875 selected ATAC matrix violates count or census constraints")
}

Matrix::writeMM(selected_counts, file.path(output, "atac_counts.mtx"))
write.table(
  peaks,
  file = file.path(output, "atac_features.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
write.table(
  data.frame(cell_id = colnames(selected_counts)),
  file = file.path(output, "atac_cells.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  con = file.path(output, "R_sessionInfo.txt"),
  useBytes = TRUE
)
