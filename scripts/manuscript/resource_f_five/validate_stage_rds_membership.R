#!/usr/bin/env Rscript
# Emit an independent sample/sex/count-membership audit for stage validation.

suppressPackageStartupMessages(library(data.table))

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) {
  stop(
    "Usage: validate_stage_rds_membership.R <raw.rds> <dge.rds> <meta.rds>",
    call. = FALSE
  )
}

raw <- readRDS(arguments[[1L]])
dge <- readRDS(arguments[[2L]])
meta <- as.data.table(readRDS(arguments[[3L]]))
required <- c("sample_id", "inferred_sex")
if (!all(required %in% names(meta)) || anyDuplicated(meta$sample_id)) {
  stop("Independent stage RDS metadata schema drift", call. = FALSE)
}
if (!identical(dim(raw), c(86369L, 1281L)) ||
    !identical(dim(dge), c(24196L, 1257L)) ||
    anyDuplicated(colnames(raw)) ||
    anyDuplicated(colnames(dge))) {
  stop("Independent stage RDS dimension/key drift", call. = FALSE)
}
audit <- meta[, .(
  sample_id = as.character(sample_id),
  inferred_sex = as.character(inferred_sex),
  raw_present = sample_id %in% colnames(raw),
  corrected_dge_present = sample_id %in% colnames(dge)
)]
if (nrow(audit) != 1281L ||
    anyNA(audit$inferred_sex) ||
    any(!audit$inferred_sex %in% c("F", "M")) ||
    sum(audit$raw_present) != 1281L ||
    sum(audit$corrected_dge_present) != 1257L) {
  stop("Independent stage RDS membership census drift", call. = FALSE)
}
write.table(
  as.data.frame(audit),
  file = stdout(),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  na = "NA",
  eol = "\n"
)
