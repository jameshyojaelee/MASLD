#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: merge_gse213621.R RUN_ROOT BAM_MANIFEST", call. = FALSE)
}
run_root <- normalizePath(args[[1]], mustWork = TRUE)
manifest_path <- normalizePath(args[[2]], mustWork = TRUE)
sentinel <- file.path(run_root, ".bg001_candidate_root")
if (!file.exists(sentinel)) stop("Missing BG-001 candidate sentinel", call. = FALSE)

chunk_root <- file.path(run_root, "counts/GSE213621/chunks")
chunk_dirs <- file.path(chunk_root, paste0("chunk", 0:3))
observed_chunk_dirs <- sort(basename(list.dirs(chunk_root, recursive = FALSE, full.names = TRUE)))
if (!identical(observed_chunk_dirs, paste0("chunk", 0:3))) {
  stop("GSE213621 requires exactly chunk0, chunk1, chunk2, and chunk3", call. = FALSE)
}
required <- unlist(lapply(chunk_dirs, function(x) {
  file.path(x, c(
    "gene_counts.txt", "gene_counts.txt.summary", "featureCounts.log", "environment.txt",
    "bam_checksums.tsv", "validation.json", "provenance.sha256", "COMPLETE"
  ))
}))
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing required chunk artifacts: ", paste(missing, collapse = ", "), call. = FALSE)

manifest <- fread(manifest_path)
expected_samples <- manifest[dataset == "GSE213621" & expected_status == "included", sample_id]
if (length(expected_samples) != 367L || anyDuplicated(expected_samples)) {
  stop("GSE213621 manifest must contain 367 unique included samples", call. = FALSE)
}
sha256_file <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("sha256sum failed for ", path, call. = FALSE)
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
verify_provenance <- function(chunk_dir) {
  provenance_path <- file.path(chunk_dir, "provenance.sha256")
  entries <- readLines(provenance_path, warn = FALSE)
  expected_names <- c(
    "bam_checksums.tsv", "environment.txt", "gene_counts.txt", "gene_counts.txt.summary",
    "featureCounts.log", "validation.json"
  )
  observed_names <- character()
  for (entry in entries) {
    match <- regexec("^([0-9a-f]{64})[[:space:]]+(.+)$", entry)
    fields <- regmatches(entry, match)[[1L]]
    if (length(fields) != 3L) stop("Malformed chunk provenance: ", provenance_path, call. = FALSE)
    artifact <- normalizePath(fields[[3L]], mustWork = TRUE)
    if (dirname(artifact) != normalizePath(chunk_dir, mustWork = TRUE)) {
      stop("Chunk provenance path escapes its chunk: ", artifact, call. = FALSE)
    }
    observed_names <- c(observed_names, basename(artifact))
    if (sha256_file(artifact) != fields[[2L]]) stop("Chunk provenance hash drift: ", artifact, call. = FALSE)
  }
  if (!setequal(observed_names, expected_names) || anyDuplicated(observed_names)) {
    stop("Chunk provenance does not cover the exact required artifacts: ", chunk_dir, call. = FALSE)
  }
}
for (chunk_index in 0:3) {
  chunk_dir <- chunk_dirs[[chunk_index + 1L]]
  verify_provenance(chunk_dir)
  payload <- jsonlite::read_json(file.path(chunk_dir, "validation.json"), simplifyVector = TRUE)
  start <- chunk_index * 92L + 1L
  end <- min(start + 91L, length(expected_samples))
  expected_chunk <- expected_samples[start:end]
  if (!identical(payload$status, "PASS") || payload$n_genes != 86369L ||
      payload$n_samples != length(expected_chunk) || !identical(payload$samples, expected_chunk)) {
    stop("Chunk ", chunk_index, " validation payload does not match the frozen manifest", call. = FALSE)
  }
  observed_hashes <- c(
    counts_sha256 = sha256_file(file.path(chunk_dir, "gene_counts.txt")),
    summary_sha256 = sha256_file(file.path(chunk_dir, "gene_counts.txt.summary")),
    log_sha256 = sha256_file(file.path(chunk_dir, "featureCounts.log")),
    environment_sha256 = sha256_file(file.path(chunk_dir, "environment.txt"))
  )
  if (!all(vapply(names(observed_hashes), function(field) {
    identical(payload[[field]], unname(observed_hashes[[field]]))
  }, logical(1)))) stop("Chunk validation hashes are stale for chunk ", chunk_index, call. = FALSE)
}

ann <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
tables <- lapply(chunk_dirs, function(d) fread(file.path(d, "gene_counts.txt"), skip = 1L, check.names = FALSE))
for (i in seq_along(tables)) {
  x <- tables[[i]]
  if (!identical(names(x)[seq_len(6L)], ann)) stop("Chunk ", i - 1L, " annotation schema mismatch", call. = FALSE)
  if (nrow(x) != 86369L || anyDuplicated(x$Geneid)) stop("Chunk ", i - 1L, " gene contract failure", call. = FALSE)
  if (anyNA(x[, ..ann])) stop("Chunk ", i - 1L, " contains missing annotations", call. = FALSE)
}

key_string <- function(x) do.call(paste, c(x[, ..ann], sep = "\034"))
reference_keys <- key_string(tables[[1]])
if (anyDuplicated(reference_keys)) stop("Reference annotation keys are not unique", call. = FALSE)

ordered <- lapply(seq_along(tables), function(i) {
  x <- tables[[i]]
  keys <- key_string(x)
  if (!setequal(keys, reference_keys)) stop("Chunk ", i - 1L, " annotation keys differ", call. = FALSE)
  x[match(reference_keys, keys)]
})

sample_names <- unlist(lapply(ordered, function(x) basename(names(x)[-(1:6)])), use.names = FALSE)
sample_names <- sub("[.]Aligned[.]sortedByCoord[.]out[.]bam$", "", sample_names)
if (length(sample_names) != 367L || anyDuplicated(sample_names) || !identical(sample_names, expected_samples)) {
  stop("Merged count sample order does not exactly match the 367-sample manifest", call. = FALSE)
}

merged <- copy(ordered[[1]][, ..ann])
for (x in ordered) merged <- cbind(merged, x[, -(1:6), with = FALSE])
if (ncol(merged) != 373L) stop("Merged count matrix has wrong column count", call. = FALSE)
merged_column <- 7L
for (x in ordered) {
  for (source_column in 7:ncol(x)) {
    if (!identical(merged[[merged_column]], x[[source_column]])) {
      stop("Keyed merge changed source sample column ", names(x)[[source_column]], call. = FALSE)
    }
    merged_column <- merged_column + 1L
  }
}

summaries <- lapply(chunk_dirs, function(d) fread(file.path(d, "gene_counts.txt.summary"), check.names = FALSE))
statuses <- summaries[[1]][[1]]
if (anyDuplicated(statuses)) stop("Duplicate summary status", call. = FALSE)
for (i in seq_along(summaries)) {
  x <- summaries[[i]]
  if (!setequal(x[[1]], statuses)) stop("Chunk ", i - 1L, " summary statuses differ", call. = FALSE)
  summaries[[i]] <- x[match(statuses, x[[1]])]
}
merged_summary <- copy(summaries[[1]][, 1, with = FALSE])
for (x in summaries) merged_summary <- cbind(merged_summary, x[, -1, with = FALSE])
summary_samples <- basename(names(merged_summary)[-1])
summary_samples <- sub("[.]Aligned[.]sortedByCoord[.]out[.]bam$", "", summary_samples)
if (!identical(summary_samples, sample_names)) stop("Merged count and summary sample orders differ", call. = FALSE)

out_dir <- file.path(run_root, "counts/GSE213621/.merge_staging")
if (file.exists(out_dir)) stop("Refusing existing GSE213621 merge staging directory", call. = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_counts <- file.path(out_dir, "gene_counts.txt")
out_summary <- paste0(out_counts, ".summary")
if (file.exists(out_counts) || file.exists(out_summary)) stop("Refusing to overwrite staged GSE213621 output", call. = FALSE)
tmp_counts <- paste0(out_counts, ".tmp")
tmp_summary <- paste0(out_summary, ".tmp")
comment <- paste0("# Program:featureCounts v2.1.1; BG-001 keyed merge; Command:\"featureCounts\" \"-p\" \"--countReadPairs\" \"-B\" \"-s\" \"2\" \"-a\" \"/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz\"; chunks=4")
writeLines(comment, tmp_counts)
fwrite(merged, tmp_counts, sep = "\t", append = TRUE, col.names = TRUE, quote = FALSE)
fwrite(merged_summary, tmp_summary, sep = "\t", quote = FALSE)
if (!file.rename(tmp_counts, out_counts)) stop("Atomic staged-count rename failed", call. = FALSE)
if (!file.rename(tmp_summary, out_summary)) stop("Atomic staged-summary rename failed", call. = FALSE)

logs <- unlist(lapply(chunk_dirs, function(d) readLines(file.path(d, "featureCounts.log"), warn = FALSE)))
writeLines(logs, file.path(out_dir, "featureCounts.log"))
environments <- unlist(lapply(seq_along(chunk_dirs), function(i) c(
  paste0("chunk_environment_begin\t", i - 1L),
  readLines(file.path(chunk_dirs[[i]], "environment.txt"), warn = FALSE),
  paste0("chunk_environment_end\t", i - 1L)
)))
writeLines(environments, file.path(out_dir, "environment.txt"))
cat(sprintf("Staged %d genes x %d samples at %s for full validation before publication\n",
            nrow(merged), ncol(merged) - 6L, out_counts))
