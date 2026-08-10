#!/usr/bin/env Rscript

script_path <- normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "../../../.."), mustWork = TRUE)
candidate <- file.path(root, "Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09")
if (!all(file.exists(file.path(candidate, c("SEALED.json", "SEAL_VALIDATED.json"))))) stop("validated seal required")

read_matrix <- function(path) {
  x <- read.delim(gzfile(path), sep = "\t", comment.char = "!", quote = "\"",
                  check.names = FALSE, stringsAsFactors = FALSE)
  values <- as.matrix(x[, -1, drop = FALSE]); storage.mode(values) <- "double"
  rownames(values) <- as.character(x[[1]])
  values
}

a <- read_matrix(file.path(root, "Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/sources/GSE106737/GSE106737_series_matrix.txt.gz"))
b <- read_matrix(file.path(candidate, "sources/GSE83452/GSE83452_series_matrix.txt.gz"))
common <- sort(intersect(rownames(a), rownames(b)))
if (length(common) < 5000L) stop("insufficient common probes for sample fingerprint")
probes <- common[seq_len(5000L)]
a <- a[probes, , drop = FALSE]; b <- b[probes, , drop = FALSE]
correlations <- cor(a, b, use = "pairwise.complete.obs", method = "pearson")
rows <- vector("list", ncol(a))
for (i in seq_len(ncol(a))) {
  best <- which.max(correlations[i, ])
  reciprocal <- which.max(correlations[, best]) == i
  rows[[i]] <- data.frame(
    gse106737_sample = colnames(a)[i], gse83452_best_sample = colnames(b)[best],
    pearson_r = correlations[i, best], reciprocal_best = reciprocal,
    duplicate_threshold_pass = reciprocal && correlations[i, best] >= 0.9999,
    n_fixed_probes = length(probes), probe_selection = "first_5000_lexicographic_common_probe_ids",
    stringsAsFactors = FALSE
  )
}
output <- do.call(rbind, rows)
path <- file.path(candidate, "source_gates/sample_fingerprint_matches.tsv")
dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
write.table(output, path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("FINGERPRINT_COMPLETE\t", sum(output$duplicate_threshold_pass), "/", nrow(output), "\n", sep = "")
