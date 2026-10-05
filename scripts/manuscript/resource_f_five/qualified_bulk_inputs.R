# Inputs for the count- and identity-qualified Resource bulk candidate.
# A completed annotation map alone does not certify a corrected count matrix.

read_qualified_bulk_inputs <- function(counts_path, roster_path, certificate_path,
                                       original_sample_ids) {
  paths <- c(counts_path, roster_path, certificate_path)
  if (any(!nzchar(paths)) || any(!file.exists(paths))) {
    stop("Qualified bulk analysis requires counts, unit roster and count certificate", call. = FALSE)
  }
  certificate <- jsonlite::fromJSON(certificate_path)
  if (!isTRUE(certificate$annotation_resolution_complete) ||
      !identical(as.integer(certificate$n_unresolved_annotations), 0L) ||
      !isTRUE(certificate$counts_complete)) {
    stop("Annotation correspondence or corrected fragment counting is incomplete", call. = FALSE)
  }
  sha256 <- function(path) {
    strsplit(system2("sha256sum", shQuote(normalizePath(path)),
                    stdout = TRUE)[[1L]], "[[:space:]]+")[[1L]][[1L]]
  }
  digest <- sha256(counts_path)
  if (!identical(certificate$counts_sha256, digest)) {
    stop("Corrected counts do not match the completed count certificate", call. = FALSE)
  }
  annotation_path <- certificate$annotation_summary
  if (!is.character(annotation_path) || length(annotation_path) != 1L ||
      !nzchar(annotation_path) || !file.exists(annotation_path) ||
      !identical(certificate$annotation_summary_sha256, sha256(annotation_path)) ||
      !identical(certificate$primary_sample_roster_sha256, sha256(roster_path))) {
    stop("Count certificate must match the completed annotation evidence and selected unit roster", call. = FALSE)
  }
  annotation <- jsonlite::fromJSON(annotation_path)
  if (!isTRUE(annotation$annotation_resolution_complete) ||
      !identical(as.integer(annotation$n_unresolved_annotations), 0L)) {
    stop("Source annotation correspondence is incomplete", call. = FALSE)
  }
  roster <- data.table::fread(roster_path, colClasses = "character")
  fields <- c("sample_id", "analysis_unit_id", "unit_evidence")
  if (!all(fields %in% names(roster)) || !nrow(roster) ||
      anyNA(roster[, ..fields]) ||
      any(vapply(roster[, ..fields], function(x) any(!nzchar(x)), logical(1))) ||
      anyDuplicated(roster$sample_id) || anyDuplicated(roster$analysis_unit_id) ||
      !all(roster$sample_id %in% original_sample_ids)) {
    stop("Primary roster must contain supported, unique units from the original model roster", call. = FALSE)
  }
  counts <- readRDS(counts_path)
  if (!is.matrix(counts) || !is.numeric(counts) || !nrow(counts) || !ncol(counts) ||
      is.null(rownames(counts)) || is.null(colnames(counts)) ||
      anyNA(rownames(counts)) || anyNA(colnames(counts)) ||
      any(!nzchar(rownames(counts))) || any(!nzchar(colnames(counts))) ||
      anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)) ||
      !all(roster$sample_id %in% colnames(counts))) {
    stop("Corrected counts must have unique gene and sample identifiers covering the primary roster", call. = FALSE)
  }
  # Bound temporary allocations while checking a potentially large matrix.
  for (start in seq.int(1L, nrow(counts), by = 1024L)) {
    block <- counts[start:min(start + 1023L, nrow(counts)), , drop = FALSE]
    if (any(!is.finite(block)) || any(block < 0) || any(block != floor(block))) {
      stop("Corrected counts must be nonnegative integer fragment contributions", call. = FALSE)
    }
  }
  list(counts = counts, roster = roster, certificate = certificate)
}

assert_bulk_unit_alignment <- function(meta, sample_ids) {
  if (!"analysis_unit_id" %in% names(meta)) return(invisible(FALSE))
  if (!all(c("sample_id", "unit_evidence") %in% names(meta)) ||
      anyDuplicated(meta$sample_id)) {
    stop("A qualified model requires unique sample rows and unit evidence", call. = FALSE)
  }
  units <- meta$analysis_unit_id[match(sample_ids, meta$sample_id)]
  evidence <- meta$unit_evidence[match(sample_ids, meta$sample_id)]
  if (anyNA(units) || any(!nzchar(units)) || anyDuplicated(units) ||
      anyNA(evidence) || any(!nzchar(evidence))) {
    stop("A model sample lacks a supported unique analysis unit", call. = FALSE)
  }
  invisible(TRUE)
}
