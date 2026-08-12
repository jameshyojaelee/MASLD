#!/usr/bin/env Rscript
# Re-derive the Resource candidate's raw and fitted bulk censuses without fitting a model.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: preflight_bulk_census.R F_LEGACY_RAW F_LEGACY_DGE F_FIVE_DGE")
}
for (path in args) {
  if (!file.exists(path) || nzchar(Sys.readlink(path))) {
    stop("Missing or symlinked census input: ", path)
  }
}

raw <- readRDS(args[[1L]])
legacy <- readRDS(args[[2L]])
five <- readRDS(args[[3L]])

object_dim <- function(object, label) {
  value <- dim(object)
  if (is.null(value) && is.list(object) && !is.null(object$counts)) {
    value <- dim(object$counts)
  }
  if (length(value) != 2L || anyNA(value)) {
    stop("Cannot derive two-dimensional census for ", label)
  }
  as.integer(value)
}

cohort_count <- function(object, label) {
  if (is.null(object$samples) || !"dataset" %in% colnames(object$samples)) {
    stop("Missing dataset field in ", label)
  }
  value <- unique(as.character(object$samples$dataset))
  value <- value[!is.na(value) & nzchar(value)]
  length(value)
}

raw_dim <- object_dim(raw, "F_legacy raw")
legacy_dim <- object_dim(legacy, "F_legacy DGE")
five_dim <- object_dim(five, "F_five DGE")
metrics <- c(
  raw_n_genes = raw_dim[[1L]],
  raw_n_samples = raw_dim[[2L]],
  legacy_n_genes = legacy_dim[[1L]],
  legacy_n_samples = legacy_dim[[2L]],
  legacy_n_cohorts = cohort_count(legacy, "F_legacy DGE"),
  five_n_genes = five_dim[[1L]],
  five_n_samples = five_dim[[2L]],
  five_n_cohorts = cohort_count(five, "F_five DGE")
)

write.table(
  data.frame(metric = names(metrics), value = unname(metrics)),
  file = stdout(), sep = "\t", quote = FALSE, row.names = FALSE,
  col.names = TRUE, eol = "\n"
)
