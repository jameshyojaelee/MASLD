#!/usr/bin/env Rscript
# Inspect deposited donor ATAC measurements; no fitting or genotype access.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L, nzchar(Sys.getenv("SLURM_JOB_ID")))
asset <- args[[1]]
peaks_file <- args[[2]]
out <- args[[3]]
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
env <- new.env(parent = baseenv())
objects <- load(asset, envir = env)
if (!all(c("peaks.dds", "peaks.vst") %in% objects)) {
  write_json(list(status = "expected_objects_missing", objects = objects),
             file.path(out, "inspection.json"), pretty = TRUE, auto_unbox = TRUE)
  stop("Deposit lacks expected peaks.dds / peaks.vst objects")
}
peaks <- read.delim(gzfile(peaks_file), check.names = FALSE, comment.char = "",
                   stringsAsFactors = FALSE)
stopifnot(all(c("#chr", "start", "end", "peakID") %in% names(peaks)),
          !anyDuplicated(peaks$peakID), all(peaks$end > peaks$start))
matrix_stats <- function(x) {
  if (is.null(x)) return(list(available = FALSE))
  z <- list(available = TRUE, class = class(x), dimensions = dim(x),
            elements = length(x), missing = 0, nonfinite = 0, negative = 0,
            noninteger = 0, minimum = Inf, maximum = -Inf)
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  for (i in seq.int(1L, nrow(x), by = 10000L)) {
    block <- as.matrix(x[i:min(i + 9999L, nrow(x)), , drop = FALSE])
    z$missing <- z$missing + sum(is.na(block))
    z$nonfinite <- z$nonfinite + sum(!is.finite(block))
    finite <- block[is.finite(block)]
    z$negative <- z$negative + sum(finite < 0)
    z$noninteger <- z$noninteger + sum(finite != floor(finite))
    if (length(finite)) {
      z$minimum <- min(z$minimum, min(finite))
      z$maximum <- max(z$maximum, max(finite))
    }
  }
  if (!is.finite(z$minimum)) z$minimum <- NA_real_
  if (!is.finite(z$maximum)) z$maximum <- NA_real_
  z
}
axis_info <- function(x) {
  list(n = length(x), present = !is.null(x), missing_or_blank = sum(is.na(x) | x == ""),
       duplicate_instances = sum(duplicated(x)),
       ordered_R_serialized_sha256 = digest(x, algo = "sha256"))
}
field_info <- function(x) lapply(x, function(y) list(class = class(y), length = length(y)))
numeric_summary <- function(x) {
  finite <- x[is.finite(x)]
  list(n = length(x), finite = length(finite),
       quantiles = if (length(finite)) as.list(quantile(finite, c(0, .25, .5, .75, 1))) else NULL)
}
inspect <- function(x) {
  if (!is(x, "SummarizedExperiment")) {
    return(list(class = class(x), supported_inspection = FALSE, dimensions = dim(x)))
  }
  ids <- rownames(x)
  m <- match(ids, peaks$peakID)
  good <- !is.na(m)
  rr <- if (is(x, "RangedSummarizedExperiment")) rowRanges(x) else NULL
  coord <- list(available = is(rr, "GRanges"), exact_BED_to_GRanges_identity = NA,
                mismatches = NA_integer_)
  if (is(rr, "GRanges") && length(rr) == nrow(x)) {
    equal <- as.character(seqnames(rr))[good] == peaks[["#chr"]][m[good]] &
      start(rr)[good] - 1L == peaks$start[m[good]] & end(rr)[good] == peaks$end[m[good]]
    coord$exact_BED_to_GRanges_identity <- all(equal) && all(good)
    coord$mismatches <- sum(!equal)
  }
  list(class = class(x), supported_inspection = TRUE, dimensions = dim(x),
       row_identity = axis_info(ids), donor_column_identity = axis_info(colnames(x)),
       peak_ID_matches = sum(good), peak_ID_unmatched = sum(!good),
       matched_peak_chromosomes = as.list(table(peaks[["#chr"]][m[good]])),
       matched_peak_widths_bp = numeric_summary(peaks$end[m[good]] - peaks$start[m[good]]),
       coordinates = coord, assay_names = assayNames(x),
       assay_summaries = setNames(lapply(assayNames(x), function(a) matrix_stats(assay(x, a))), assayNames(x)),
       row_annotation_fields = field_info(as.list(rowData(x))),
       donor_annotation_fields = field_info(as.list(colData(x))),
       metadata_fields = field_info(metadata(x)))
}
dds <- env[["peaks.dds"]]
vst <- env[["peaks.vst"]]
summaries <- list(peaks.dds = inspect(dds), peaks.vst = inspect(vst))
is_dds <- is(dds, "DESeqDataSet")
norm <- if (is_dds) list(size_factors = matrix_stats(sizeFactors(dds)),
                        normalization_factors = matrix_stats(normalizationFactors(dds)),
                        design = paste(deparse(design(dds)), collapse = " ")) else list(available = FALSE)
compatible_axes <- identical(rownames(dds), rownames(vst)) && identical(colnames(dds), colnames(vst))
raw_ok <- is_dds && "counts" %in% assayNames(dds) &&
  summaries$peaks.dds$assay_summaries$counts$nonfinite == 0 &&
  summaries$peaks.dds$assay_summaries$counts$negative == 0 &&
  summaries$peaks.dds$assay_summaries$counts$noninteger == 0
libraries <- list(available = FALSE)
if (is_dds && "counts" %in% assayNames(dds)) {
  libraries <- list(available = TRUE,
    full_deposited_peak_matrix_donor_count_totals = numeric_summary(colSums(assay(dds, "counts"))),
    semantics = "Totals sum every deposited autosomal peak row. They are not verified total sequenced reads, all mapped reads, or off-peak-inclusive library sizes.")
  cd <- as.list(colData(dds))
  named <- grep("lib|depth|read|fragment|count|sizeFactor|offset|gc", names(cd), ignore.case = TRUE, value = TRUE)
  libraries$possible_library_or_offset_fields <- setNames(lapply(named, function(k) {
    if (is.numeric(cd[[k]])) numeric_summary(cd[[k]]) else list(class = class(cd[[k]]), values_not_exported = TRUE)
  }), named)
}
ids_ok <- summaries$peaks.dds$supported_inspection &&
  summaries$peaks.dds$row_identity$present && summaries$peaks.dds$donor_column_identity$present &&
  summaries$peaks.dds$row_identity$duplicate_instances == 0 &&
  summaries$peaks.dds$donor_column_identity$duplicate_instances == 0 &&
  summaries$peaks.dds$row_identity$missing_or_blank == 0 &&
  summaries$peaks.dds$donor_column_identity$missing_or_blank == 0 &&
  summaries$peaks.dds$peak_ID_unmatched == 0
summary <- list(status = "aggregate_inspection_complete_no_profile_accuracy_run", job_id = Sys.getenv("SLURM_JOB_ID"),
  loaded_object_classes = setNames(lapply(objects, function(n) class(env[[n]])), objects),
  objects = summaries, normalization = norm, dds_vst_axes_identical = compatible_axes,
  library_metadata = libraries,
  measurement_identity_candidate = isTRUE(raw_ok && ids_ok && compatible_axes),
  original_source_peak_rows = nrow(peaks), reference_build = "GRCh38_from_source_deposit_and_Methods",
  inferential_unit = "Source describes 138 donors after merging libraries; inspect actual column count/identity above. No donor crosswalk was independently verified by this inspection.",
  native_endpoint = "WASP-filtered read counts overlapping source consensus peaks; not caQTL beta",
  source_declared_processing = "featureCounts; EDASeq GC correction; DESeq2 library-size normalization and variance stabilization",
  GC_correction_contents = "Inspect explicit normalization and annotation fields above; source declaration alone does not identify which stored object/offset implements GC correction",
  preprocessing_limit = "Deposited normalization/VST used the full original cohort. A future donor-held fitted comparison must learn permitted normalization/calibration within training donors. Do not fit to source-wide transformed outcomes and claim fold-local preprocessing.",
  source_family = "same Currin2025 donors and peaks as development caQTL labels; not independent external validation",
  exposure_limit = "Foundation-checkpoint source exposure unresolved; inherited reference-locus exposure and this source's label development are separate limitations",
  genotype_crosswalk = "Not verified; donor labels in this object cannot be assumed to match pre-swap GSE26105 sample titles",
  values_exported = FALSE, genotype_values_accessed = FALSE, protected_outcomes_accessed = FALSE,
  inference_or_fitting_run = FALSE)
write_json(summary, file.path(out, "inspection.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
cat("Aggregate inspection complete; measurement identity candidate:", summary$measurement_identity_candidate, "\n")
