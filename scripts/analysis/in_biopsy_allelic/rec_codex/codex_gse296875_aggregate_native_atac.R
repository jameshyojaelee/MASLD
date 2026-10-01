#!/usr/bin/env Rscript
# Native author count aggregation only. No analytical access to meta.data.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !grepl("^[0-9]+$", Sys.getenv("SLURM_JOB_ID")) ||
    !nzchar(Sys.getenv("SLURM_JOB_NODELIST"))) stop("requires one owned output path and a compute allocation")
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
rec <- file.path(root, "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z")
output <- normalizePath(args[[1L]], mustWork = TRUE)
expected_output <- file.path(rec, paste0("codex_gse296875_aligned_native_pseudobulk_", Sys.getenv("SLURM_JOB_ID")))
if (!identical(output, normalizePath(expected_output, mustWork = TRUE))) stop("unexpected output directory")
source <- file.path(root, "Analysis/SingleCell/continual_integration/candidates/external-gse296875-download-v44-2026-08-18T000000Z/GSE296875_processed.rds")
expected_rds_sha <- "6af74b7dbdd3839841f90cb7f3e47a3184b4dd3dd02242f6f3bdbd14533cfcab"
expected_geometry_sha <- "50a52d097a5232edcf136dd53d6549ca2d1aa5dce9a57cf055ff10b7b2e0f3e3"
expected_total <- 341091266
cap <- 2 * 1024^3
phase <- "installed_reader_evidence"
require_true <- function(value, message) if (!isTRUE(value)) stop(message, call. = FALSE)
sha <- function(path) {
  x <- system2("sha256sum", c("--", shQuote(path)), stdout = TRUE, stderr = TRUE)
  require_true(is.null(attr(x, "status")) && length(x) == 1L && grepl("^[0-9a-f]{64} ", x), "SHA256 failed")
  substr(x, 1L, 64L)
}
fingerprint <- function(path) {
  x <- file.info(path)
  require_true(nrow(x) == 1L && !is.na(x$size) && !x$isdir, "missing regular input")
  list(bytes = unname(x$size), mtime = as.numeric(x$mtime), ctime = as.numeric(x$ctime))
}
output_bytes <- function() {
  paths <- list.files(output, all.files = TRUE, recursive = TRUE, full.names = TRUE, no.. = TRUE)
  sum(file.info(paths)$size)
}
check_cap <- function() require_true(output_bytes() <= cap, "inclusive output exceeds 2 GiB")
save_json <- function(value, name) {
  final <- file.path(output, name); temporary <- paste0(final, ".partial")
  require_true(!file.exists(final) && !file.exists(temporary), "refusing metadata overwrite")
  jsonlite::write_json(value, temporary, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  check_cap()
  require_true(file.rename(temporary, final), "atomic metadata rename failed")
}
write_npy_header <- function(connection) {
  # NPY v1.0: little-endian int32, [unit, peak], Fortran layout (unit fastest).
  dictionary <- "{'descr': '<i4', 'fortran_order': True, 'shape': (273, 306706), }"
  padding <- (64L - ((10L + nchar(dictionary, type = "bytes") + 1L) %% 64L)) %% 64L
  header <- paste0(dictionary, strrep(" ", padding), "\n")
  size <- nchar(header, type = "bytes")
  require_true(size < 65536L, "NPY v1 header length exceeds uint16")
  writeBin(c(as.raw(147L), charToRaw("NUMPY"), as.raw(c(1L, 0L))), connection)
  writeBin(as.raw(c(size %% 256L, size %/% 256L)), connection)
  writeBin(charToRaw(header), connection)
  10L + size
}
binaries <- list(); geometry <- NULL
tryCatch({
  reader_files <- character()
  for (package in c("Matrix", "SeuratObject", "jsonlite")) {
    location <- find.package(package, quiet = TRUE)
    require_true(length(location) == 1L, paste("installed reader unavailable:", package))
    reader_files <- c(reader_files, file.path(location, "DESCRIPTION"),
      list.files(file.path(location, "R"), pattern = "\\.(rdb|rdx)$", full.names = TRUE),
      list.files(file.path(location, "libs"), pattern = "\\.so$", full.names = TRUE))
  }
  reader_hashes <- setNames(vapply(reader_files, sha, character(1L)), reader_files)
  suppressPackageStartupMessages({ library(Matrix); library(SeuratObject); library(jsonlite) })
  phase <- "canonical_author_RDS_compute_hash"
  before <- fingerprint(source)
  require_true(before$bytes == 8534583098 && sha(source) == expected_rds_sha, "canonical author RDS differs")
  unit_path <- file.path(output, "unit_identity.tsv")
  cell_path <- file.path(output, "cell_identity.tsv")
  unit_sha <- sha(unit_path); cell_sha <- sha(cell_path)
  units <- read.delim(unit_path, colClasses = "character", check.names = FALSE, quote = "")
  cells <- read.delim(cell_path, colClasses = "character", check.names = FALSE, quote = "")
  require_true(identical(names(units), c("unit_index", "donor_id", "well_id", "native_source_celltype", "nuclei", "sampling_state", "halfA_nuclei", "halfB_nuclei", "halfA_sampling_state", "halfB_sampling_state", "halfA_missing", "halfB_missing", "split_diagnostic_available")), "unit identity schema differs")
  require_true(identical(names(cells), c("well_id", "raw_barcode", "cell_id", "donor_id", "source_label", "half_sort_sha256", "nucleus_half")), "cell identity schema differs")
  require_true(nrow(units) == 273L && nrow(cells) == 68398L &&
    identical(as.integer(units$unit_index), 0:272), "identity dimensions or unit order differ")
  phase <- "whole_RDS_deserialization_peaks_access_only"
  message("ATAC: readRDS deserializes entire object; analytical access only peaks counts and dimnames")
  object <- readRDS(source)
  require_true(inherits(object, "Seurat"), "author object is not Seurat")
  assays <- slot(object, "assays")
  require_true("peaks" %in% names(assays), "author object lacks native peaks")
  storage <- attributes(assays[["peaks"]])
  if ("counts" %in% names(storage)) {
    counts <- storage[["counts"]]
    layout <- "peaks serialized counts attribute"
  } else if ("layers" %in% names(storage) && "counts" %in% names(storage[["layers"]])) {
    counts <- storage[["layers"]][["counts"]]
    layout <- "peaks serialized layers/counts attribute"
  } else stop("unsupported native peaks count storage")
  require_true(inherits(counts, "dgCMatrix") && identical(dim(counts), c(306706L, 68398L)), "native sparse count layout differs")
  require_true(identical(colnames(counts), cells$cell_id) && !anyDuplicated(colnames(counts)) &&
    !anyNA(rownames(counts)) && !anyDuplicated(rownames(counts)), "full native axes do not match exact identity partition")
  unit_key <- paste(units$donor_id, units$native_source_celltype, sep = "\034")
  assignment <- match(paste(cells$donor_id, cells$source_label, sep = "\034"), unit_key)
  require_true(!anyNA(assignment) && !anyDuplicated(unit_key) &&
    identical(tabulate(assignment, nbins = 273L), as.integer(units$nuclei)), "native stratum assignment differs")
  require_true(all(cells$nucleus_half %in% c("A", "B")), "invalid nucleus-half identity")
  half_assignment <- assignment + ifelse(cells$nucleus_half == "A", 273L, 546L)
  require_true(identical(tabulate(half_assignment, nbins = 819L)[274:546], as.integer(units$halfA_nuclei)) &&
    identical(tabulate(half_assignment, nbins = 819L)[547:819], as.integer(units$halfB_nuclei)) &&
    all(as.integer(units$halfA_nuclei) + as.integer(units$halfB_nuclei) == as.integer(units$nuclei)), "half identity assignment differs")
  group <- sparseMatrix(i = rep(seq_along(assignment), 2L), j = c(assignment, half_assignment), x = 1,
    dims = c(68398L, 819L))
  # Geometry-only serialized attributes; never inspect GRanges mcols/data.
  range_evidence <- list(serialized_ranges_present = "ranges" %in% names(storage),
    geometry_readable = FALSE, coordinate_metadata_only = TRUE,
    independently_confirmed_original_fragment_assignment_convention = FALSE)
  range_chr <- range_start <- range_end <- NULL
  if ("ranges" %in% names(storage)) {
    native_ranges <- storage[["ranges"]]
    range_evidence$serialized_class <- class(native_ranges)
    fields <- attributes(native_ranges)
    if (inherits(native_ranges, "GRanges") && all(c("seqnames", "ranges", "seqinfo") %in% names(fields))) {
      sequence <- attributes(fields[["seqnames"]])
      intervals <- attributes(fields[["ranges"]])
      seqinfo <- attributes(fields[["seqinfo"]])
      if (all(c("values", "lengths") %in% names(sequence)) &&
          all(c("start", "width") %in% names(intervals))) {
        range_chr <- rep(as.character(sequence[["values"]]), times = sequence[["lengths"]])
        range_start <- as.numeric(intervals[["start"]])
        range_end <- range_start + as.numeric(intervals[["width"]]) - 1
        require_true(length(range_chr) == length(range_start) && length(range_start) == 306706L &&
          all(is.finite(range_start) & is.finite(range_end) & range_start >= 1 & range_end >= range_start), "serialized GRanges coordinates invalid")
        range_evidence$geometry_readable <- TRUE
        range_evidence$interval_class <- class(fields[["ranges"]])
        range_evidence$author_object_coordinate_semantics <- "GRanges/IRanges integer start and width; closed end=start+width-1"
        range_evidence$genome_assembly_seqinfo_tags <- if ("genome" %in% names(seqinfo)) unique(as.character(seqinfo[["genome"]])) else character()
        range_evidence$seqinfo_class <- class(fields[["seqinfo"]])
      }
    }
  }
  partitions <- c("full", "halfA", "halfB")
  binary_paths <- setNames(file.path(output, paste0("ATAC_", partitions, "_counts_int32.npy")), partitions)
  geometry_path <- file.path(output, "native_peak_geometry.tsv")
  require_true(!any(file.exists(binary_paths)) && !any(file.exists(paste0(binary_paths, ".partial"))) &&
    !file.exists(geometry_path) && !file.exists(paste0(geometry_path, ".partial")), "refusing ATAC data overwrite")
  require_true(output_bytes() + 3 * 306706 * 273 * 4 + 4096 + 25 * 1024^2 < cap, "ATAC outputs would exceed saved-size cap")
  header_bytes <- integer(3L)
  for (index in seq_along(partitions)) {
    binaries[[partitions[[index]]]] <- file(paste0(binary_paths[[index]], ".partial"), "wb")
    header_bytes[[index]] <- write_npy_header(binaries[[partitions[[index]]]])
  }
  geometry <- file(paste0(geometry_path, ".partial"), "wb")
  fields <- c("peak_id", "chromosome", "source_start_1based_closed", "source_end_1based_closed", "bed_start_0based", "bed_end_half_open")
  writeLines(paste(fields, collapse = "\t"), geometry, sep = "\r\n", useBytes = TRUE)
  source_unit_totals <- aggregate_unit_totals <- numeric(819L)
  source_total <- 0; maximum_aggregate <- 0
  phase <- "all_native_peak_blocks_and_exact_integer_conservation"
  for (first in seq.int(1L, 306706L, by = 8192L)) {
    last <- min(306706L, first + 8192L - 1L)
    ids <- rownames(counts)[first:last]
    pieces <- regmatches(ids, regexec("^(chr(?:[1-9]|1[0-9]|2[0-2]))[:-]([0-9]+)-([0-9]+)$", ids))
    require_true(all(lengths(pieces) == 4L), "native autosomal peak identifier geometry differs; no peaks filtered")
    chromosomes <- vapply(pieces, `[[`, character(1L), 2L)
    start <- as.numeric(vapply(pieces, `[[`, character(1L), 3L))
    end <- as.numeric(vapply(pieces, `[[`, character(1L), 4L))
    require_true(all(is.finite(start) & is.finite(end) & start >= 1 & end >= start), "native peak bounds invalid")
    if (isTRUE(range_evidence$geometry_readable)) {
      require_true(identical(chromosomes, range_chr[first:last]) &&
        identical(start, range_start[first:last]) && identical(end, range_end[first:last]), "serialized GRanges and unchanged peak identifier geometry differ")
    }
    writeLines(paste(ids, chromosomes, sprintf("%.0f", start), sprintf("%.0f", end),
      sprintf("%.0f", start - 1), sprintf("%.0f", end), sep = "\t"), geometry, sep = "\r\n", useBytes = TRUE)
    block <- counts[first:last, , drop = FALSE]
    if (length(block@x)) for (entry in seq.int(1, length(block@x), by = 1000000)) {
      values <- block@x[entry:min(length(block@x), entry + 999999)]
      require_true(all(is.finite(values) & values >= 0 & values == floor(values)), "native entries violate integer count validity")
    }
    aggregated <- as.matrix(block %*% group)
    require_true(all(is.finite(aggregated)) && all(aggregated >= 0) &&
      all(aggregated == floor(aggregated)) && max(aggregated) <= .Machine$integer.max,
      "native aggregate cannot be represented exactly as int32")
    source_feature_totals <- Matrix::rowSums(block)
    require_true(all(aggregated[,1:273,drop=FALSE] == aggregated[,274:546,drop=FALSE] + aggregated[,547:819,drop=FALSE]), "ATAC full=A+B failed")
    require_true(identical(as.numeric(rowSums(aggregated[,1:273,drop=FALSE])), as.numeric(source_feature_totals)),
      "ATAC per-native-peak conservation failed")
    source_unit_totals <- source_unit_totals + as.numeric(crossprod(group, Matrix::colSums(block)))
    aggregate_unit_totals <- aggregate_unit_totals + colSums(aggregated)
    source_total <- source_total + sum(source_feature_totals)
    maximum_aggregate <- max(maximum_aggregate, max(aggregated))
    for (index in seq_along(partitions)) {
      columns <- (index - 1L) * 273L + seq_len(273L)
      writeBin(as.integer(t(aggregated[,columns,drop=FALSE])), binaries[[partitions[[index]]]], size = 4L, endian = "little")
    }
    check_cap()
    message("ATAC native peaks aggregated: ", last, "/306706")
  }
  for (connection in binaries) close(connection)
  binaries <- list()
  close(geometry); geometry <- NULL
  require_true(source_total == expected_total && sum(aggregate_unit_totals[1:273]) == expected_total &&
    identical(source_unit_totals, aggregate_unit_totals), "ATAC grand/group/raw denominator conservation failed")
  require_true(all(aggregate_unit_totals[1:273][units$sampling_state == "unsampled"] == 0) &&
    all(aggregate_unit_totals[274:546][units$halfA_sampling_state == "unsampled"] == 0) &&
    all(aggregate_unit_totals[547:819][units$halfB_sampling_state == "unsampled"] == 0), "unsampled native ATAC full/half stratum is not empty")
  require_true(sha(paste0(geometry_path, ".partial")) == expected_geometry_sha, "full native geometry differs from frozen fixture metadata")
  for (index in seq_along(partitions)) {
    require_true(fingerprint(paste0(binary_paths[[index]], ".partial"))$bytes == 306706 * 273 * 4 + header_bytes[[index]], "ATAC NPY output length differs")
    require_true(file.rename(paste0(binary_paths[[index]], ".partial"), binary_paths[[index]]), "atomic ATAC NPY rename failed")
  }
  require_true(file.rename(paste0(geometry_path, ".partial"), geometry_path), "atomic geometry rename failed")
  phase <- "final_source_reader_and_metadata_checks"
  require_true(identical(fingerprint(source), before) && sha(unit_path) == unit_sha && sha(cell_path) == cell_sha,
    "ATAC source or assignment changed during aggregation")
  require_true(identical(unname(vapply(reader_files, sha, character(1L))), unname(reader_hashes)), "installed R reader source changed")
  range_evidence$all_ranges_match_unchanged_peak_identifiers <- if (isTRUE(range_evidence$geometry_readable)) TRUE else NA
  save_json(range_evidence, "ATAC_serialized_ranges_geometry.json")
  writeLines(capture.output(sessionInfo()), file.path(output, "R_sessionInfo.txt"), useBytes = TRUE)
  save_json(list(completed = TRUE, source_path = source, full_RDS_computed_sha256 = expected_rds_sha,
    source_stat_before = before, source_stat_after = fingerprint(source), native_sparse_layout = layout,
    native_peaks = 306706L, same_nuclei = 68398L, strata = 273L,
    total_native_peak_events = source_total, unit_totals_full = unname(aggregate_unit_totals[1:273]),
    unit_totals_halfA = unname(aggregate_unit_totals[274:546]), unit_totals_halfB = unname(aggregate_unit_totals[547:819]),
    maximum_aggregate = maximum_aggregate, dtype = "little-endian int32",
    cast_bound = "all aggregates <= .Machine$integer.max and grand total 341091266 < 2^31",
    NPY_shape = c(273L,306706L), NPY_fortran_order = TRUE, NPY_header_bytes = unname(header_bytes),
    per_feature_and_group_conservation = TRUE, full_equals_halfA_plus_halfB = TRUE,
    native_geometry_sha256 = expected_geometry_sha,
    native_coordinate_convention = "existing exporter GRanges_1_based_closed, start-1/end only for identity metadata",
    independently_confirmed_source_coordinate_origin = FALSE,
    count_unit = "native author peak-count events", unique_fragment_or_molecule_units_verified = FALSE,
    whole_RDS_deserialized = TRUE, analytical_meta_data_accessed = FALSE,
    installed_reader_hashes = as.list(reader_hashes), no_filtering_normalization_fits_or_target_selection = TRUE), "ATAC_receipt.json")
}, error = function(error) {
  for (connection in binaries) try(close(connection), silent = TRUE)
  if (!is.null(geometry)) try(close(geometry), silent = TRUE)
  if (requireNamespace("jsonlite", quietly = TRUE)) try(save_json(list(completed = FALSE, phase = phase,
    error = conditionMessage(error), call = paste(deparse(conditionCall(error)), collapse = " "),
    partial_files_preserved = TRUE), "ATAC_failure.json"), silent = TRUE)
  message("ATAC FAILED phase=", phase, ": ", conditionMessage(error))
  quit(save = "no", status = 1L)
})
