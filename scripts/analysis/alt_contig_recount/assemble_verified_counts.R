#!/usr/bin/env Rscript
# Combine read/fragment counts only after all annotation relationships are resolved.
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
grouped_path <- normalizePath(annotation$grouped_counting_annotation, mustWork = TRUE)
correspondence_path <- file.path(dirname(annotation_path), "annotation_correspondence.tsv")
grouped_sha <- sha(grouped_path)
correspondence_sha <- sha(correspondence_path)
if (!identical(grouped_sha, annotation$grouped_counting_annotation_sha256) ||
    !identical(correspondence_sha, annotation$annotation_correspondence_sha256)) stop("Used annotation/correspondence hash differs")
genes <- sort(unique(fread(correspondence_path)$counting_gene_id))
counts <- matrix(0L, nrow = length(genes), ncol = nrow(roster), dimnames = list(genes, roster$sample_id))
sample_evidence <- vector("list", nrow(roster))
for (i in seq_len(nrow(roster))) {
  sample <- roster$sample_id[[i]]
  path <- file.path(count_root, sample, "gene_counts.tsv")
  certificate <- fromJSON(file.path(count_root, sample, "count_completion.json"))
  if (length(certificate$count_reader_source_sha256) != 1L ||
      !grepl("^[a-f0-9]{64}$", certificate$count_reader_source_sha256)) stop("Count reader source hash missing: ", sample)
  if (certificate$sample_id != sample || certificate$annotation_summary_sha256 != annotation_sha ||
      certificate$counts_tsv_sha256 != sha(path) || certificate$maximum_contribution_per_molecule != 1L ||
      isTRUE(certificate$fractional_renormalization)) stop("Invalid molecule certificate: ", sample)
  evidence_path <- file.path(count_root, sample, "assignment_execution.json")
  if (!identical(certificate$grouped_counting_annotation_sha256, grouped_sha) ||
      !identical(certificate$annotation_correspondence_sha256, correspondence_sha) ||
      !identical(certificate$assignment_execution_sha256, sha(evidence_path))) stop("Counting execution lineage differs: ", sample)
  execution <- fromJSON(evidence_path)
  argv_fc <- execution$featurecounts_argv
  allowed_flags <- c("-T", "-s", "-M", "-R", "-a", "-o", "--Rpath", "-p", "--countReadPairs", "-B")
  if (!length(argv_fc) || basename(argv_fc[[1L]]) != "featureCounts" ||
      any(!argv_fc[startsWith(argv_fc, "-")] %in% allowed_flags) ||
      !all(c("python", "featureCounts", "samtools") %in% names(execution$versions)) ||
      any(!nzchar(unlist(execution$versions[c("python", "featureCounts", "samtools")])))) {
    stop("Unknown assignment options or executable versions missing: ", sample)
  }
  option <- function(name) {
    where <- which(argv_fc == name)
    if (length(where) != 1L || where == length(argv_fc)) stop("Missing/repeated execution flag: ", name)
    argv_fc[[where + 1L]]
  }
  value_options <- c("-T", "-s", "-R", "-a", "-o", "--Rpath")
  switches <- c("-M", "-p", "--countReadPairs", "-B")
  position <- 2L
  while (position < length(argv_fc)) {
    token <- argv_fc[[position]]
    if (token %in% value_options) {
      option(token)
      position <- position + 2L
    } else if (token %in% switches) {
      position <- position + 1L
    } else stop("Unrecognized featureCounts eligibility argument: ", token)
  }
  if (length(certificate$layout) != 1L || !certificate$layout %in% c("SINGLE", "PAIRED")) {
    stop("Invalid library layout: ", sample)
  }
  paired <- identical(certificate$layout, "PAIRED")
  if (!identical(execution$status, "complete") || execution$sample_id != sample ||
      !identical(execution$layout, certificate$layout) || execution$featurecounts_returncode != 0L ||
      execution$samtools_sort_returncode != 0L || !identical(option("-a"), grouped_path) ||
      !identical(option("-s"), "2") || !identical(option("-R"), "BAM") ||
      sum(argv_fc == "-M") != 1L ||
      any(vapply(c("-p", "--countReadPairs", "-B"), function(flag) sum(argv_fc == flag) != as.integer(paired), logical(1))) ||
      !identical(execution$annotation_summary_sha256, annotation_sha) ||
      !identical(execution$grouped_counting_annotation_sha256, grouped_sha) ||
      !identical(execution$annotation_correspondence_sha256, correspondence_sha) ||
      !identical(execution$source_bam_before, execution$source_bam_after) ||
      !identical(certificate$source_bam, execution$source_bam_before) ||
      !identical(tail(argv_fc, 1L), execution$source_bam_before$path)) stop("Actual assignment command/source evidence differs: ", sample)
  sample_evidence[[i]] <- data.table(sample_id = sample,
    count_certificate_sha256 = sha(file.path(count_root, sample, "count_completion.json")),
    count_reader_source_sha256 = certificate$count_reader_source_sha256,
    assignment_execution_path = normalizePath(evidence_path),
    assignment_execution_sha256 = sha(evidence_path), source_bam = execution$source_bam_before$path,
    source_bam_size_bytes = execution$source_bam_before$size_bytes,
    source_bam_mtime_ns = execution$source_bam_before$mtime_ns,
    layout = certificate$layout, featurecounts_argv = toJSON(argv_fc))
  x <- fread(path)
  if (!identical(names(x), c("gene_id", sample)) || anyDuplicated(x$gene_id) ||
      !identical(sort(x$gene_id), genes)) stop("Gene universe differs: ", sample)
  values <- x[[sample]][match(genes, x$gene_id)]
  if (anyNA(values) || any(values < 0) || any(values != floor(values)) || any(values > .Machine$integer.max)) stop("Invalid integer fragment counts")
  counts[, i] <- as.integer(values)
}
if (uniqueN(rbindlist(sample_evidence)$count_reader_source_sha256) != 1L) stop("Count matrix would mix different molecule-reading rules")
dir.create(out_root, recursive = TRUE)
matrix_path <- file.path(out_root, "counts_raw.rds")
saveRDS(counts, matrix_path)
sample_manifest <- file.path(out_root, "sample_count_evidence.tsv")
fwrite(rbindlist(sample_evidence), sample_manifest, sep = "\t")
certificate <- list(annotation_resolution_complete = TRUE, n_unresolved_annotations = 0L, counts_complete = TRUE,
                    counts_sha256 = sha(matrix_path), annotation_summary = annotation_path,
                    grouped_counting_annotation = grouped_path, grouped_counting_annotation_sha256 = grouped_sha,
                    annotation_correspondence_sha256 = correspondence_sha,
                    sample_count_evidence = normalizePath(sample_manifest), sample_count_evidence_sha256 = sha(sample_manifest),
                    count_reader_source_sha256 = unique(rbindlist(sample_evidence)$count_reader_source_sha256),
                    annotation_summary_sha256 = annotation_sha, primary_sample_roster = roster_path,
                    primary_sample_roster_sha256 = sha(roster_path), n_genes = nrow(counts), n_samples = ncol(counts),
                    unit = "one library per source-linked participant or conservative genetic identity cluster",
                    unit_evidence_counts = as.list(table(roster$unit_evidence)),
                    unit_definitions = sort(unique(roster$unit_evidence)),
                    counting = "integer read/fragment counts; complete within-gene alignment sets only; not UMI-deduplicated molecules",
                    candidate_not_adopted = TRUE)
writeLines(toJSON(certificate, auto_unbox = TRUE, pretty = TRUE), file.path(out_root, "count_completion.json"))
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
cat(toJSON(certificate, auto_unbox = TRUE), "\n")
