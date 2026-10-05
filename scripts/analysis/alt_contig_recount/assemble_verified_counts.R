#!/usr/bin/env Rscript
# Combine independently counted molecules only after all annotation relationships are resolved.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 4L) stop("Usage: annotation_summary.json primary_sample_ids.tsv sample_count_root output_root")
annotation_path <- normalizePath(argv[[1]], mustWork = TRUE)
roster_path <- normalizePath(argv[[2]], mustWork = TRUE)
count_root <- normalizePath(argv[[3]], mustWork = TRUE)
out_root <- argv[[4]]
if (file.exists(out_root)) stop("Refusing to overwrite output root")
sha <- function(path) strsplit(system2("sha256sum", shQuote(path), stdout = TRUE), "[[:space:]]+")[[1]][[1]]
annotation <- fromJSON(annotation_path)
if (!isTRUE(annotation$annotation_resolution_complete) || annotation$n_unresolved_annotations != 0L) stop("Annotation correspondence unresolved")
roster <- fread(roster_path, colClasses = "character")
unit_fields <- c("sample_id", "analysis_unit_id", "unit_evidence")
supported_unit_evidence <- c("source_labels_plus_stable_filtered_genetic_identity_cluster_not_verified_person",
                             "authoritative_source_participant_pairing_not_screened_in_A1")
if (!all(unit_fields %in% names(roster)) || !nrow(roster) ||
    anyNA(roster[, ..unit_fields]) ||
    any(vapply(roster[, ..unit_fields], function(x) any(!nzchar(x)), logical(1))) ||
    anyDuplicated(roster$sample_id) || anyDuplicated(roster$analysis_unit_id) ||
    !all(roster$unit_evidence %in% supported_unit_evidence)) {
  stop("Primary roster must use one library per supported source-linked or genetic-identity unit")
}
annotation_sha <- sha(annotation_path)
genes <- sort(unique(fread(file.path(dirname(annotation_path), "annotation_correspondence.tsv"))$counting_gene_id))
counts <- matrix(0L, nrow = length(genes), ncol = nrow(roster), dimnames = list(genes, roster$sample_id))
for (i in seq_len(nrow(roster))) {
  sample <- roster$sample_id[[i]]
  path <- file.path(count_root, sample, "gene_counts.tsv")
  certificate <- fromJSON(file.path(count_root, sample, "count_completion.json"))
  if (certificate$sample_id != sample || certificate$annotation_summary_sha256 != annotation_sha ||
      certificate$counts_tsv_sha256 != sha(path) || certificate$maximum_contribution_per_molecule != 1L ||
      isTRUE(certificate$fractional_renormalization)) stop("Invalid molecule certificate: ", sample)
  x <- fread(path)
  if (!identical(names(x), c("gene_id", sample)) || anyDuplicated(x$gene_id) ||
      !identical(sort(x$gene_id), genes)) stop("Gene universe differs: ", sample)
  values <- x[[sample]][match(genes, x$gene_id)]
  if (anyNA(values) || any(values < 0) || any(values != floor(values)) || any(values > .Machine$integer.max)) stop("Invalid integer fragment counts")
  counts[, i] <- as.integer(values)
}
dir.create(out_root, recursive = TRUE)
matrix_path <- file.path(out_root, "counts_raw.rds")
saveRDS(counts, matrix_path)
certificate <- list(annotation_resolution_complete = TRUE, n_unresolved_annotations = 0L, counts_complete = TRUE,
                    counts_sha256 = sha(matrix_path), annotation_summary = annotation_path,
                    annotation_summary_sha256 = annotation_sha, primary_sample_roster = roster_path,
                    primary_sample_roster_sha256 = sha(roster_path), n_genes = nrow(counts), n_samples = ncol(counts),
                    unit = "one library per source-linked participant or conservative genetic identity cluster",
                    unit_evidence_counts = as.list(table(roster$unit_evidence)),
                    unit_definitions = sort(unique(roster$unit_evidence)),
                    counting = "integer molecules; complete within-gene alignment sets only",
                    candidate_not_adopted = TRUE)
writeLines(toJSON(certificate, auto_unbox = TRUE, pretty = TRUE), file.path(out_root, "count_completion.json"))
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
cat(toJSON(certificate, auto_unbox = TRUE), "\n")
