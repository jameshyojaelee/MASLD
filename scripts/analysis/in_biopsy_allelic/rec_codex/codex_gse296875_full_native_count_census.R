#!/usr/bin/env Rscript
# Compute-only, exposed-development native author peak-count census.
# readRDS deserializes the whole object; analytical access excludes meta.data.
suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
  library(jsonlite)
})

root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
rec <- file.path(root, "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z")
source <- file.path(root, "Analysis/SingleCell/continual_integration/candidates/external-gse296875-download-v44-2026-08-18T000000Z/GSE296875_processed.rds")
lock_path <- file.path(root, "Analysis/SingleCell/continual_integration/reference/gse296875_source_lock_v44.json")
fixture_axis <- file.path(rec, "codex_gse296875_metadata_census_21999159/native_custom_merged_atac_geometry.tsv")
helper <- file.path(root, "scripts/analysis/in_biopsy_allelic/rec_codex/codex_gse296875_full_native_count_census.R")
launcher <- file.path(root, "scripts/analysis/in_biopsy_allelic/rec_codex/run_codex_gse296875_full_native_count_census.sbatch")
expected_lock_sha <- "6be686b4187309f038097cf73d3f2a8dd34e72d666bc92733a53e99a1ccbb620"
expected_rds_sha <- "6af74b7dbdd3839841f90cb7f3e47a3184b4dd3dd02242f6f3bdbd14533cfcab"
expected_geometry_sha <- "50a52d097a5232edcf136dd53d6549ca2d1aa5dce9a57cf055ff10b7b2e0f3e3"
cap <- 2 * 1024^2
chunk_size <- 1000000L
axis_chunk_size <- 8192L
job <- Sys.getenv("SLURM_JOB_ID")
if (!grepl("^[0-9]+$", job) || !nzchar(Sys.getenv("SLURM_JOB_NODELIST"))) {
  stop("native count census requires a SLURM compute allocation")
}
if (length(commandArgs(trailingOnly = TRUE))) stop("this fixed-source census accepts no overrides")
output <- file.path(rec, paste0("codex_gse296875_full_native_count_census_", job))
if (file.exists(output) || !dir.create(output, mode = "0750")) {
  stop("refusing an existing or unavailable output directory")
}
phase <- "source_evidence"
require_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}
sha <- function(path) {
  result <- system2("sha256sum", c("--", shQuote(path)), stdout = TRUE, stderr = TRUE)
  require_true(is.null(attr(result, "status")) && length(result) == 1L &&
    grepl("^[0-9a-f]{64} ", result), paste("SHA256 failed:", path))
  substr(result, 1L, 64L)
}
fingerprint <- function(path) {
  x <- file.info(path)
  require_true(nrow(x) == 1L && !is.na(x$size) && !x$isdir, paste("missing file:", path))
  list(bytes = unname(x$size), mtime = as.numeric(x$mtime), ctime = as.numeric(x$ctime))
}
output_bytes <- function() {
  paths <- list.files(output, all.files = TRUE, recursive = TRUE, full.names = TRUE, no.. = TRUE)
  if (!length(paths)) return(0)
  sum(file.info(paths)$size)
}
save_json <- function(value, name) {
  target <- file.path(output, name)
  temporary <- paste0(target, ".partial")
  require_true(!file.exists(target) && !file.exists(temporary), "refusing metadata overwrite")
  write_json(value, temporary, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  require_true(output_bytes() < cap, "inclusive saved output reaches the 2 MiB cap")
  require_true(file.rename(temporary, target), "atomic metadata rename failed")
}
open_digest <- function(name) {
  path <- file.path(output, name)
  require_true(!file.exists(path), "refusing digest overwrite")
  list(path = path, connection = pipe(paste("sha256sum >", shQuote(path)), "w"))
}
finish_digest <- function(stream) {
  status <- close(stream$connection)
  require_true(is.null(status) || identical(status, 0L), "axis SHA256 subprocess failed")
  value <- readLines(stream$path, warn = FALSE)
  require_true(length(value) == 1L && grepl("^[0-9a-f]{64} ", value), "axis SHA256 output malformed")
  substr(value, 1L, 64L)
}

tryCatch({
  before <- fingerprint(source)
  require_true(before$bytes == 8534583098, "full author RDS byte size differs")
  require_true(sha(lock_path) == expected_lock_sha, "source lock SHA256 differs")
  lock <- read_json(lock_path, simplifyVector = TRUE)
  require_true(lock$processed_source$bytes == before$bytes &&
    lock$processed_source$sha256 == expected_rds_sha, "processed source lock differs")
  fixture_before <- fingerprint(fixture_axis)
  require_true(fixture_before$bytes == 20635809 && sha(fixture_axis) == expected_geometry_sha,
    "completed fixture identity geometry differs")
  source_hashes <- list(helper = sha(helper), launcher = sha(launcher), source_lock = sha(lock_path))
  phase <- "full_RDS_byte_hash_on_compute"
  computed_rds_sha <- sha(source)
  require_true(computed_rds_sha == expected_rds_sha, "full author RDS computed SHA256 differs")
  manifest <- list(
    schema_version = "codex-gse296875-full-native-peak-count-census-v1",
    job_id = job, hostname = Sys.info()[["nodename"]], started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    argv = commandArgs(), resources = list(cpus = 1, memory_GiB = 64, hours = 90),
    scientific_role = "project-exposed public development; count validity only",
    source = list(path = source, stat_before = before, official_url = lock$processed_source$url,
      recorded_sha256 = expected_rds_sha, computed_sha256 = computed_rds_sha,
      full_RDS_byte_hash_computed_on_compute = TRUE),
    fixture_identity_metadata = list(path = fixture_axis, sha256 = expected_geometry_sha,
      bytes = fixture_before$bytes, molecular_payload_read = FALSE),
    source_hashes = source_hashes, output_cap_bytes_exclusive = cap,
    count_unit = "native author peak-count events",
    unique_fragment_or_molecule_units_verified = FALSE,
    independently_confirmed_source_coordinate_origin = FALSE,
    deserialization_limit = paste("readRDS deserializes the entire author object, including stored metadata.",
      "Analytical accesses are restricted to peaks assay count storage and dimnames; meta.data is never accessed."),
    exclusions = c("clinical values", "fitting", "target ranking", "axis selection or filtering",
      "protected sources", "receiving outcomes", "fragment or UMI equivalence"))
  save_json(manifest, "manifest.json")
  phase <- "readRDS_and_native_sparse_accessor"
  message("loading canonical full author RDS; no clinical metadata access")
  object <- readRDS(source)
  require_true(inherits(object, "Seurat"), "processed object is not Seurat")
  assays <- slot(object, "assays")
  require_true("peaks" %in% names(assays), "author object lacks peaks assay")
  storage <- attributes(assays[["peaks"]])
  if ("counts" %in% names(storage)) {
    counts <- storage[["counts"]]
    layout <- "peaks serialized counts attribute"
  } else if ("layers" %in% names(storage) && "counts" %in% names(storage[["layers"]])) {
    counts <- storage[["layers"]][["counts"]]
    layout <- "peaks serialized layers/counts attribute"
  } else stop("native peaks count storage unsupported")
  require_true(inherits(counts, "dgCMatrix"), "bounded native census requires dgCMatrix layout")
  nr <- nrow(counts); nc <- ncol(counts); nnz <- length(counts@x)
  peaks <- rownames(counts); cells <- colnames(counts)
  axis <- list(rows = nr, columns = nc, stored_entries = nnz, matrix_class = class(counts),
    storage_layout = layout, duplicate_peaks = sum(duplicated(peaks)),
    duplicate_cells = sum(duplicated(cells)),
    invalid_cell_ids = sum(is.na(cells) | !grepl("^well[1-8]_[ACGT]+-1$", cells)))
  save_json(axis, "native_axis_layout.json")
  require_true(nr == 306706L && nc == 68398L && length(peaks) == nr && length(cells) == nc,
    "full native dimensions or dimnames differ from the frozen identity census")
  require_true(!anyNA(peaks) && !anyNA(cells) && axis$duplicate_peaks == 0 &&
    axis$duplicate_cells == 0 && axis$invalid_cell_ids == 0, "duplicate or invalid native axes")
  require_true(length(counts@i) == nnz && length(counts@p) == nc + 1L &&
    counts@p[1L] == 0 && tail(counts@p, 1L) == nnz && all(diff(counts@p) >= 0),
    "native sparse pointers differ from dgCMatrix layout")
  phase <- "all_peak_geometry_and_native_axis_hash"
  geometry <- open_digest("native_geometry.sha256")
  fields <- c("peak_id", "chromosome", "source_start_1based_closed", "source_end_1based_closed",
    "bed_start_0based", "bed_end_half_open")
  writeLines(paste(fields, collapse = "\t"), geometry$connection, sep = "\r\n", useBytes = TRUE)
  contigs <- setNames(integer(22L), paste0("chr", 1:22))
  invalid_geometry <- 0L
  for (first in seq.int(1L, nr, by = axis_chunk_size)) {
    last <- min(nr, first + axis_chunk_size - 1L)
    ids <- peaks[first:last]
    pieces <- regmatches(ids, regexec("^(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y))[:-]([0-9]+)-([0-9]+)$", ids))
    valid <- lengths(pieces) == 4L
    chromosome <- rep(NA_character_, length(ids)); start <- end <- rep(NA_real_, length(ids))
    chromosome[valid] <- vapply(pieces[valid], `[[`, character(1L), 2L)
    start[valid] <- as.numeric(vapply(pieces[valid], `[[`, character(1L), 3L))
    end[valid] <- as.numeric(vapply(pieces[valid], `[[`, character(1L), 4L))
    valid <- valid & is.finite(start) & is.finite(end) & start >= 1 & end >= start &
      chromosome %in% names(contigs)
    invalid_geometry <- invalid_geometry + sum(!valid)
    if (all(valid)) {
      tab <- table(chromosome)
      contigs[names(tab)] <- contigs[names(tab)] + as.integer(tab)
      rows <- paste(ids, chromosome, sprintf("%.0f", start), sprintf("%.0f", end),
        sprintf("%.0f", start - 1), sprintf("%.0f", end), sep = "\t")
      writeLines(rows, geometry$connection, sep = "\r\n", useBytes = TRUE)
    }
  }
  geometry_sha <- finish_digest(geometry)
  save_json(list(invalid_geometry = invalid_geometry, contigs = as.list(contigs),
    parsed_coordinates_follow_existing_exporter_convention = "GRanges_1_based_closed",
    independently_confirmed_source_coordinate_origin = FALSE,
    bed_conversion_for_identity_digest_only = "start-1,end; no overlaps or count transformation",
    native_geometry_sha256 = geometry_sha, fixture_geometry_sha256 = expected_geometry_sha,
    complete_geometry_hash = invalid_geometry == 0L,
    hash_codec = "six-field TSV with header and CRLF for each row, native order, UTF8",
    exact_fixture_axis_agreement = geometry_sha == expected_geometry_sha), "geometry_check.json")
  require_true(invalid_geometry == 0L && geometry_sha == expected_geometry_sha,
    "full native geometry does not exactly match the frozen fixture identity axis; no peaks filtered")
  barcode <- open_digest("native_barcode_axis.sha256")
  for (first in seq.int(1L, nc, by = axis_chunk_size)) {
    writeLines(cells[first:min(nc, first + axis_chunk_size - 1L)], barcode$connection, useBytes = TRUE)
  }
  barcode_sha <- finish_digest(barcode)
  phase <- "bounded_all_native_sparse_entry_checks"
  invalid <- c(nonfinite = 0, negative = 0, noninteger = 0, invalid_row_index = 0)
  explicit_zeros <- 0; positive_entries <- 0; maximum_entry <- 0
  if (nnz) for (first in seq.int(1, nnz, by = chunk_size)) {
    last <- min(nnz, first + chunk_size - 1)
    values <- counts@x[first:last]
    indices <- counts@i[first:last]
    finite <- is.finite(values)
    invalid["nonfinite"] <- invalid["nonfinite"] + sum(!finite)
    invalid["negative"] <- invalid["negative"] + sum(values[finite] < 0)
    invalid["noninteger"] <- invalid["noninteger"] + sum(values[finite] != floor(values[finite]))
    invalid["invalid_row_index"] <- invalid["invalid_row_index"] + sum(indices < 0L | indices >= nr)
    explicit_zeros <- explicit_zeros + sum(values[finite] == 0)
    positive_entries <- positive_entries + sum(values[finite] > 0)
    if (any(finite)) maximum_entry <- max(maximum_entry, values[finite])
  }
  entries <- list(invalid = as.list(invalid), explicit_stored_zeros = explicit_zeros,
    positive_entries = positive_entries, maximum_entry = maximum_entry, chunk_entries = chunk_size)
  save_json(entries, "entry_validity.json")
  require_true(all(invalid == 0), "native stored entries fail finite/nonnegative/integer/index checks")
  phase <- "raw_native_column_totals_summary"
  totals <- Matrix::colSums(counts)
  require_true(length(totals) == nc && all(is.finite(totals)) && all(totals >= 0),
    "raw native column totals are invalid")
  raw_totals <- list(n_columns = length(totals), zero_columns = sum(totals == 0),
    positive_columns = sum(totals > 0), min = min(totals), max = max(totals),
    mean = mean(totals), total_peak_count_events = sum(totals),
    quantiles = as.list(quantile(totals, c(0, .01, .25, .5, .75, .99, 1), names = TRUE)),
    all_totals_integer = all(totals == floor(totals)),
    totals_above_exact_double_integer_range = sum(totals > 2^53),
    global_total_within_exact_double_integer_range = sum(totals) <= 2^53,
    zero_columns_preserved = TRUE, count_unit = "native author peak-count events")
  phase <- "final_source_and_output_checks"
  require_true(identical(fingerprint(source), before), "author RDS changed during census")
  require_true(identical(fingerprint(fixture_axis), fixture_before), "fixture identity metadata changed")
  require_true(sha(helper) == source_hashes$helper && sha(launcher) == source_hashes$launcher &&
    sha(lock_path) == source_hashes$source_lock, "executed source or lock changed during census")
  writeLines(capture.output(sessionInfo()), file.path(output, "R_sessionInfo.txt"), useBytes = TRUE)
  require_true(output_bytes() < cap, "inclusive output reaches the 2 MiB cap")
  summary <- list(schema_version = manifest$schema_version, native_count_validity_census_completed = TRUE,
    scientific_role = manifest$scientific_role, finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    axis = axis, entry_checks = entries, raw_column_totals = raw_totals,
    native_barcode_axis_sha256 = barcode_sha, barcode_hash_codec = "UTF8 cell_id plus LF, native order",
    native_geometry_sha256 = geometry_sha, exact_fixture_axis_agreement = TRUE,
    count_unit = manifest$count_unit, unique_fragment_or_molecule_units_verified = FALSE,
    parsed_coordinates_follow_existing_exporter_convention = "GRanges_1_based_closed",
    independently_confirmed_source_coordinate_origin = FALSE,
    full_RDS_computed_sha256 = computed_rds_sha, source_stat_after = fingerprint(source),
    analytical_meta_data_accessed = FALSE, whole_RDS_deserialized = TRUE,
    clinical_joins_performed = FALSE, counts_or_axes_exported = FALSE, model_or_target_scoring_performed = FALSE,
    output_cap_bytes_exclusive = cap, saved_bytes_before_final_summary = output_bytes(),
    final_inclusive_cap_checked_before_atomic_completion = TRUE)
  save_json(summary, "summary.json")
  message("native count validity census complete; saved bytes=", output_bytes())
}, error = function(error) {
  failure <- list(native_count_validity_census_completed = FALSE, phase = phase,
    message = conditionMessage(error), call = paste(deparse(conditionCall(error)), collapse = " "),
    failed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), partial_files_preserved = TRUE)
  try(save_json(failure, "failure.json"), silent = TRUE)
  message("FAILED phase=", phase, ": ", conditionMessage(error))
  quit(save = "no", status = 1L)
})
