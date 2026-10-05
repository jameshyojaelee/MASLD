#!/usr/bin/env Rscript
# Exact feature-unit comparison only; no projection or gene-name matching.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: frozen_preprocessing_state.rds annotation_correspondence.tsv new_output_directory")
if (file.exists(args[[3L]])) stop("Refusing to overwrite review output")
sha <- function(path) strsplit(system2("sha256sum", shQuote(path), stdout = TRUE), "[[:space:]]+")[[1L]][[1L]]
state <- readRDS(args[[1L]])
ids <- state$exact_gene_ids
stopifnot(length(ids) == 5000L, !anyDuplicated(ids))
correspondence <- fread(args[[2L]], colClasses = "character")
stopifnot(!anyDuplicated(correspondence$gene_id))
match_index <- match(ids, correspondence$gene_id)
review <- correspondence[match_index]
review[, frozen_basis_gene_id := ids]
review[, exact_version_present := !is.na(match_index)]
review[, frozen_base_id := sub("[.][0-9]+$", "", frozen_basis_gene_id)]
available_bases <- sub("[.][0-9]+$", "", correspondence$gene_id)
review[, missing_exact_version_but_base_present := !exact_version_present & frozen_base_id %in% available_bases]
review[, gene_object_changed := exact_version_present & counting_gene_id != frozen_basis_gene_id]
target_disposition <- setNames(correspondence$disposition, correspondence$gene_id)
review[, ALT_merged_into_primary := gene_object_changed & target_disposition[counting_gene_id] == "primary_gene"]
review[, unresolved_source := exact_version_present & disposition == "unresolved"]
review[, target_id_known := exact_version_present & !is.na(counting_gene_id) & nzchar(counting_gene_id)]
target_counts <- review[target_id_known == TRUE, .N, by = counting_gene_id]
duplicates <- target_counts[N > 1L]
review[, basis_rows_per_counting_target := target_counts$N[match(counting_gene_id, target_counts$counting_gene_id)]]
review[, duplicate_counting_target_among_basis_rows := !is.na(basis_rows_per_counting_target) & basis_rows_per_counting_target > 1L]
report <- list(status = "exact_feature_units_review_only_no_projection",
  frozen_basis_rows = length(ids), exact_version_missing = sum(!review$exact_version_present),
  missing_exact_version_but_same_base_present = sum(review$missing_exact_version_but_base_present),
  changed_gene_objects = sum(review$gene_object_changed, na.rm = TRUE),
  ALT_representations_merged_into_primary = sum(review$ALT_merged_into_primary, na.rm = TRUE),
  duplicate_counting_target_groups = nrow(duplicates),
  basis_rows_in_duplicate_target_groups = sum(duplicates$N),
  unresolved_source_rows = sum(review$unresolved_source),
  excluded_or_unassigned_target_rows = sum(!review$target_id_known),
  source_dispositions = as.list(table(review$disposition, useNA = "ifany")),
  source_sha256 = list(frozen_preprocessing_state = sha(args[[1L]]), annotation_correspondence = sha(args[[2L]])),
  symbol_matching_used = FALSE, missing_gene_zero_substitution = FALSE,
  interpretation = "Any changed or duplicated counting object prevents treating the old basis rows as unchanged expression features. Exact-version presence alone does not certify the eventual corrected DGE retains all basis genes. Normalization and control-library transport remain unresolved.")
dir.create(args[[3L]], recursive = TRUE)
fwrite(review, file.path(args[[3L]], "basis_5000_exact_counting_objects.tsv"), sep = "\t")
fwrite(duplicates, file.path(args[[3L]], "duplicate_counting_target_groups.tsv"), sep = "\t")
writeLines(toJSON(report, auto_unbox = TRUE, pretty = TRUE), file.path(args[[3L]], "feature_unit_review.json"))
cat(toJSON(report, auto_unbox = TRUE), "\n")
