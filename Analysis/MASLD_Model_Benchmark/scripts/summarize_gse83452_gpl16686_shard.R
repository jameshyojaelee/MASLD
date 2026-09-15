#!/usr/bin/env Rscript

# Label-blind single-array summarization of one GSE83452 shard.
#
# Every array is read and summarized on its own: oligo::read.celfiles is called
# with exactly one path, so the median polish inside oligo::rma has a single
# column and nothing to borrow from. That is the whole point, and it is what
# model-check-091 proved on these CEL files. A batched oligo::rma moves 97
# percent of features by up to 3.85 log2, which is why the batch is held at one
# array by construction rather than by convention.
#
# Sharding is safe for the same reason. Each shard is an independent job, and
# the sentinel arrays that appear in every shard let the merge confirm that
# independent jobs on different nodes produce bit-identical summaries.
#
# Matrices are written at 17 significant digits and read back with a bitwise
# assertion. utils::write.table serialises doubles at 15 significant digits,
# which silently rounds roughly 91 percent of values on a log2 scale; that was
# measured on GSE49541 and frozen as model-check-093.
#
# Reads no label, no GEO series matrix; runs no all-sample RMA, no across-array
# normalization; fits no model.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: summarize_gse83452_gpl16686_shard.R CEL_DIR SHARD_TXT CROSSWALK_TSV SHARD_ID OUTPUT_DIR")
}
cel_dir <- args[[1L]]
shard_path <- args[[2L]]
crosswalk_path <- args[[3L]]
shard_id <- args[[4L]]
output_dir <- args[[5L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260825L)
suppressPackageStartupMessages({
  library(oligo)
  library(pd.hugene.2.0.st)
  library(jsonlite)
})

EXPECTED <- c(
  "oligo" = "1.74.0", "affyio" = "1.80.0", "affxparser" = "1.82.0",
  "pd.hugene.2.0.st" = "3.14.1"
)
observed <- vapply(
  names(EXPECTED), function(package) as.character(packageVersion(package)), character(1)
)
if (!identical(unname(observed), unname(EXPECTED))) {
  stop(sprintf(
    "GPL16686 summarization package versions differ: %s",
    paste(sprintf("%s=%s", names(observed), observed), collapse = ",")
  ))
}

sha256_of <- function(path) substr(system2("sha256sum", shQuote(path), stdout = TRUE), 1L, 64L)

digest_doubles <- function(values) {
  path <- tempfile(fileext = ".f8")
  on.exit(unlink(path), add = TRUE)
  connection <- file(path, "wb")
  writeBin(as.double(values), connection, size = 8L, endian = "little")
  close(connection)
  sha256_of(path)
}

## ---------------------------------------------------- frozen crosswalk axis
crosswalk <- utils::read.delim(
  crosswalk_path, sep = "\t", quote = "", check.names = FALSE,
  colClasses = "character", na.strings = character(0)
)
platform_features <- crosswalk$platform_feature_id
if (length(platform_features) != 53617L || anyDuplicated(platform_features) > 0L) {
  stop("frozen GPL16686 crosswalk axis is not 53617 unique features")
}
eligible <- crosswalk[crosswalk$eligible_for_gene_matrix == "true", , drop = FALSE]
gene_ids <- sort(unique(eligible$ensembl_gene_id))
feature_row <- match(eligible$platform_feature_id, platform_features)
gene_index <- match(eligible$ensembl_gene_id, gene_ids)
gene_groups <- split(feature_row, gene_index)
names(gene_groups) <- gene_ids[as.integer(names(gene_groups))]
gene_groups <- gene_groups[gene_ids]

accessions <- sort(unique(readLines(shard_path)))
accessions <- accessions[nzchar(accessions)]
cel_path <- function(accession) file.path(cel_dir, paste0(accession, ".CEL"))
if (!all(file.exists(cel_path(accessions)))) stop("a staged GSE83452 CEL is absent")

## ---------------------------------------- one array per call, never a batch
platform_matrix <- matrix(
  NA_real_, nrow = length(platform_features), ncol = length(accessions),
  dimnames = list(platform_features, accessions)
)
gene_matrix <- matrix(
  NA_real_, nrow = length(gene_ids), ncol = length(accessions),
  dimnames = list(gene_ids, accessions)
)
per_array <- vector("list", length(accessions))

for (index in seq_along(accessions)) {
  accession <- accessions[[index]]
  feature_set <- oligo::read.celfiles(cel_path(accession), verbose = FALSE)
  if (ncol(Biobase::exprs(feature_set)) != 1L) {
    stop(sprintf("%s was not read as a single array", accession))
  }
  raw <- as.numeric(Biobase::exprs(feature_set))
  raw_quantiles <- unname(stats::quantile(raw, probs = c(0, 0.25, 0.5, 0.75, 1)))

  summary <- oligo::rma(feature_set, background = TRUE, normalize = FALSE, target = "core")
  values <- Biobase::exprs(summary)
  if (ncol(values) != 1L) stop(sprintf("%s summary is not a single column", accession))
  if (!setequal(rownames(values), platform_features)) {
    stop(sprintf("%s core axis differs from the frozen crosswalk axis", accession))
  }
  column <- as.numeric(values[platform_features, 1L])
  if (!all(is.finite(column))) stop(sprintf("%s summary carries a nonfinite value", accession))
  platform_matrix[, index] <- column
  gene_matrix[, index] <- vapply(
    gene_groups, function(rows) stats::median(column[rows]), numeric(1)
  )

  summary_quantiles <- unname(stats::quantile(column, probs = c(0, 0.25, 0.5, 0.75, 1)))
  per_array[[index]] <- list(
    sample_accession = accession,
    shard = shard_id,
    arrays_read_in_this_call = 1L,
    array_digest = digest_doubles(column),
    raw_probes = length(raw),
    raw_nonpositive = sum(raw <= 0),
    raw_min = raw_quantiles[[1L]], raw_q25 = raw_quantiles[[2L]],
    raw_median = raw_quantiles[[3L]], raw_q75 = raw_quantiles[[4L]],
    raw_max = raw_quantiles[[5L]],
    summary_features = length(column),
    summary_finite = sum(is.finite(column)),
    summary_min = summary_quantiles[[1L]], summary_q25 = summary_quantiles[[2L]],
    summary_median = summary_quantiles[[3L]], summary_q75 = summary_quantiles[[4L]],
    summary_max = summary_quantiles[[5L]],
    descriptive_flag_only = TRUE,
    excluded = FALSE
  )
  rm(feature_set, summary, values, raw)
  invisible(gc(verbose = FALSE))
}

if (any(!is.finite(platform_matrix))) stop("shard platform matrix carries a nonfinite value")
if (any(!is.finite(gene_matrix))) stop("shard gene matrix carries a nonfinite value")

## --------------------------- 17 significant digits, then a bitwise round trip
write_matrix <- function(values, id_column, path) {
  text <- matrix(sprintf("%.17g", values), nrow = nrow(values), ncol = ncol(values))
  frame <- data.frame(
    id = rownames(values), text, check.names = FALSE, stringsAsFactors = FALSE
  )
  names(frame) <- c(id_column, colnames(values))
  utils::write.table(frame, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  restored <- as.matrix(utils::read.delim(
    path, sep = "\t", quote = "", check.names = FALSE, row.names = 1L
  ))
  if (!identical(dim(restored), dim(values)) ||
      !identical(rownames(restored), rownames(values)) ||
      !identical(colnames(restored), colnames(values)) ||
      !identical(as.numeric(restored), as.numeric(values))) {
    stop(sprintf("written matrix does not round-trip bitwise: %s", path))
  }
  invisible(TRUE)
}
platform_path <- file.path(output_dir, sprintf("shard_%s_platform_feature_matrix.tsv", shard_id))
gene_path <- file.path(output_dir, sprintf("shard_%s_gene_matrix.tsv", shard_id))
write_matrix(platform_matrix, "platform_feature_id", platform_path)
write_matrix(gene_matrix, "ensembl_gene_id", gene_path)

utils::write.table(
  do.call(rbind, lapply(per_array, function(r) as.data.frame(r, stringsAsFactors = FALSE))),
  file.path(output_dir, sprintf("shard_%s_qc_metrics.tsv", shard_id)),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

receipt <- list(
  schema_version = "masld-bench-gse83452-gpl16686-shard-summary-v1",
  status = "pass_label_blind_single_array_shard",
  series = "GSE83452",
  platform_id = "GPL16686",
  cohort_family_id = "antwerp_inserm_shared",
  shard = shard_id,
  method = "single_array_core_summary_normalize_false",
  arrays_read_per_summarization_call = 1L,
  arrays_summarized = length(accessions),
  platform_features = nrow(platform_matrix),
  gene_features = nrow(gene_matrix),
  package_versions = as.list(observed),
  crosswalk_sha256 = sha256_of(crosswalk_path),
  platform_matrix_sha256 = sha256_of(platform_path),
  gene_matrix_sha256 = sha256_of(gene_path),
  matrices_round_trip_bitwise = TRUE,
  per_array = per_array,
  label_blind = TRUE,
  samples_excluded = 0L,
  all_sample_RMA_run = FALSE,
  across_array_quantile_normalization_run = FALSE,
  GEO_series_matrix_read = FALSE,
  labels_read = FALSE,
  model_training_activated = FALSE,
  sealed_outcomes_read = FALSE
)
writeLines(
  jsonlite::toJSON(receipt, auto_unbox = TRUE, pretty = TRUE, digits = NA),
  file.path(output_dir, sprintf("shard_%s_receipt.json", shard_id))
)
writeLines(
  capture.output(print(utils::sessionInfo())),
  file.path(output_dir, sprintf("shard_%s_R_sessionInfo.txt", shard_id))
)
cat(sprintf("shard=%s arrays=%d features=%d genes=%d\n",
            shard_id, length(accessions), nrow(platform_matrix), nrow(gene_matrix)))
