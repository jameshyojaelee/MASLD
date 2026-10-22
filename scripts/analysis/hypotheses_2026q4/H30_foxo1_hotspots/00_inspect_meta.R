#!/usr/bin/env Rscript
# H30 step 0: list the columns and stage counts of the 774-unit metadata. Metadata only; no read counts.
args <- commandArgs(trailingOnly = TRUE)
meta <- readRDS(args[[1]])
cat("class:", class(meta), "\n")
meta <- as.data.frame(meta)
cat("dim:", dim(meta), "\n")
print(colnames(meta))
for (v in intersect(colnames(meta), c("dataset", "group_binary", "fibrosis", "fibrosis_stage", "stage", "NAS", "nas", "group", "inferred_sex"))) {
  cat("\n==", v, "\n"); print(table(meta[[v]], useNA = "ifany"))
}
