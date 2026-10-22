#!/usr/bin/env Rscript
# H30: write the 774-unit roster with stage and covariates as plain text. Metadata only; no read counts.
# Usage: 01_export_roster_stage.R META_MATCHED_RDS OUT_TSV
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 01_export_roster_stage.R META_MATCHED_RDS OUT_TSV")
if (file.exists(args[[2]])) stop("Refusing to overwrite ", args[[2]])
meta <- as.data.frame(readRDS(args[[1]]))
keep <- c("sample_id", "analysis_unit_id", "dataset", "group_binary", "fibrosis_stage", "nas_score",
          "inferred_sex", "age", "diagnosis_harmonized")
stopifnot(all(keep %in% colnames(meta)), nrow(meta) == 774L, !anyDuplicated(meta$sample_id),
          !anyDuplicated(meta$analysis_unit_id))
write.table(meta[, keep], args[[2]], sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", nrow(meta), "units\n")
