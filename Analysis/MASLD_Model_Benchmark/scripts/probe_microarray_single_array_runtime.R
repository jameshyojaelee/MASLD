#!/usr/bin/env Rscript

# Runtime capability probe for the frozen microarray Bioconductor runtime.
#
# Exercises exactly three capabilities on ONE array per platform:
#   1. Calvin CEL header parse.
#   2. Genuine single-array summarization (fRMA with the exact frozen GPL570
#      vectors; SCAN and a normalize-free oligo core summary for GPL16686).
#   3. probe -> Entrez -> Ensembl resolution through the platform .db package
#      and org.Hs.eg.db.
#
# It reads one array per platform, never runs an all-sample RMA or an
# across-array quantile normalization, never reads a GEO series matrix, never
# reads a label, never fits a model, and never exports an expression value.
# Only shapes, finite counts, and distribution summaries leave this script.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: probe_microarray_single_array_runtime.R GPL570.CEL GPL16686.CEL OUTPUT_DIR")
}
fixtures <- c(GPL570 = args[[1L]], GPL16686 = args[[2L]])
output_dir <- args[[3L]]
if (!all(file.exists(fixtures))) stop("CEL fixture is absent")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260825L)
suppressPackageStartupMessages(library(jsonlite))

PACKAGES <- c(
  "affy", "affyio", "affxparser", "oligo", "frma", "SCAN.UPC",
  "hgu133plus2cdf", "hgu133plus2.db", "hgu133plus2frmavecs",
  "pd.hugene.2.0.st", "pd.hg.u133.plus.2",
  "hugene20sttranscriptcluster.db", "AnnotationDbi", "org.Hs.eg.db"
)

installed <- rownames(installed.packages())
package_versions <- as.list(vapply(
  PACKAGES,
  function(package) {
    if (package %in% installed) as.character(packageVersion(package)) else "not_installed"
  },
  character(1)
))

attempt <- function(expression) {
  tryCatch(
    list(status = "pass", error = "not_applicable", value = expression),
    error = function(condition) {
      list(status = "failed", error = conditionMessage(condition), value = NULL)
    }
  )
}

# Only shape and distribution summaries are retained; no per-feature expression
# vector leaves this function.
summarize_matrix <- function(values) {
  values <- as.matrix(values)
  finite <- is.finite(values)
  quantiles <- if (any(finite)) {
    unname(stats::quantile(values[finite], probs = c(0, 0.25, 0.5, 0.75, 1)))
  } else {
    rep(NA_real_, 5L)
  }
  identifiers <- rownames(values)
  list(
    features = nrow(values),
    arrays = ncol(values),
    finite_values = sum(finite),
    nonfinite_values = sum(!finite),
    min = quantiles[[1L]],
    q25 = quantiles[[2L]],
    median = quantiles[[3L]],
    q75 = quantiles[[4L]],
    max = quantiles[[5L]],
    feature_ids_present = !is.null(identifiers),
    duplicate_feature_ids = !is.null(identifiers) && anyDuplicated(identifiers) > 0L,
    expression_values_exported = FALSE
  )
}

## ------------------------------------------------------- capability 1: header
header_probe <- lapply(names(fixtures), function(platform) {
  result <- attempt(affyio::read.celfile.header(fixtures[[platform]], info = "full"))
  record <- list(
    platform = platform,
    cel_bytes = as.numeric(file.size(fixtures[[platform]])),
    status = result$status,
    error = result$error,
    cdf_name = "not_available",
    grid_columns = NA_integer_,
    grid_rows = NA_integer_,
    header_fields = "not_available"
  )
  if (identical(result$status, "pass")) {
    header <- result$value
    dimensions <- header[["CEL dimensions"]]
    record$cdf_name <- as.character(header[["cdfName"]])
    record$grid_columns <- as.integer(dimensions[[1L]])
    record$grid_rows <- as.integer(dimensions[[2L]])
    record$header_fields <- paste(names(header), collapse = ";")
  }
  record
})
names(header_probe) <- names(fixtures)

## ----------------------------------------- capability 2: single-array summary
summarizer_probe <- list()

# GPL570: fRMA with the exact hgu133plus2frmavecs frozen vectors. fRMA is a
# single-array method by construction: the normalization target and the probe
# effects come from the frozen vector package, never from the other 71 arrays.
gpl570_frma <- attempt({
  batch <- affy::ReadAffy(filenames = unname(fixtures[["GPL570"]]))
  if (length(Biobase::sampleNames(batch)) != 1L) stop("GPL570 fixture is not a single array")
  loader <- new.env(parent = emptyenv())
  utils::data("hgu133plus2frmavecs", package = "hgu133plus2frmavecs", envir = loader)
  vectors <- loader[["hgu133plus2frmavecs"]]
  summary <- frma::frma(batch, input.vecs = vectors, verbose = FALSE)
  Biobase::exprs(summary)
})
summarizer_probe$GPL570_frma_hgu133plus2frmavecs <- c(
  list(
    method = "frma_with_exact_hgu133plus2frmavecs",
    arrays_read = 1L,
    across_array_normalization_used = FALSE,
    status = gpl570_frma$status,
    error = gpl570_frma$error
  ),
  if (identical(gpl570_frma$status, "pass")) summarize_matrix(gpl570_frma$value) else list()
)

# GPL570: SCAN, the requirements' second named candidate. Recorded whether or not
# it executes; fRMA above is the primary GPL570 path.
gpl570_scan <- attempt({
  summary <- SCAN.UPC::SCAN(unname(fixtures[["GPL570"]]), outFilePath = NA, verbose = FALSE)
  Biobase::exprs(summary)
})
summarizer_probe$GPL570_scan <- c(
  list(
    method = "SCAN_single_array",
    arrays_read = 1L,
    across_array_normalization_used = FALSE,
    status = gpl570_scan$status,
    error = gpl570_scan$error
  ),
  if (identical(gpl570_scan$status, "pass")) summarize_matrix(gpl570_scan$value) else list()
)

# GPL16686: SCAN with pd.hugene.2.0.st, the requirements' named candidate.
gpl16686_scan <- attempt({
  summary <- SCAN.UPC::SCAN(unname(fixtures[["GPL16686"]]), outFilePath = NA, verbose = FALSE)
  Biobase::exprs(summary)
})
summarizer_probe$GPL16686_scan <- c(
  list(
    method = "SCAN_with_pd.hugene.2.0.st",
    arrays_read = 1L,
    across_array_normalization_used = FALSE,
    status = gpl16686_scan$status,
    error = gpl16686_scan$error
  ),
  if (identical(gpl16686_scan$status, "pass")) summarize_matrix(gpl16686_scan$value) else list()
)

# GPL16686: normalize-free oligo transcript-cluster summary. RMA background and
# median-polish summarization run on this array alone; normalize = FALSE keeps
# every cross-array quantile step out of the path.
gpl16686_core <- attempt({
  raw <- oligo::read.celfiles(unname(fixtures[["GPL16686"]]), verbose = FALSE)
  if (length(Biobase::sampleNames(raw)) != 1L) stop("GPL16686 fixture is not a single array")
  summary <- oligo::rma(raw, target = "core", background = TRUE, normalize = FALSE)
  Biobase::exprs(summary)
})
summarizer_probe$GPL16686_oligo_core_single_array <- c(
  list(
    method = "single_array_core_summary_normalize_false",
    arrays_read = 1L,
    across_array_normalization_used = FALSE,
    status = gpl16686_core$status,
    error = gpl16686_core$error
  ),
  if (identical(gpl16686_core$status, "pass")) summarize_matrix(gpl16686_core$value) else list()
)

## ------------------------------- capability 3: probe -> Entrez -> Ensembl
crosswalk_probe <- function(platform, db_name) {
  attempt({
    database <- getExportedValue(db_name, db_name)
    probes <- sort(AnnotationDbi::keys(database, keytype = "PROBEID"))
    sample_probes <- utils::head(probes, 200L)
    entrez <- AnnotationDbi::select(
      database, keys = sample_probes, keytype = "PROBEID", columns = "ENTREZID"
    )
    entrez <- entrez[!is.na(entrez$ENTREZID), , drop = FALSE]
    probe_counts <- table(entrez$PROBEID)
    one_to_one_probes <- names(probe_counts)[probe_counts == 1L]
    identifiers <- sort(unique(entrez$ENTREZID))
    ensembl <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = identifiers, keytype = "ENTREZID", columns = "ENSEMBL"
    )
    ensembl <- ensembl[!is.na(ensembl$ENSEMBL), , drop = FALSE]
    entrez_counts <- table(ensembl$ENTREZID)
    list(
      platform = platform,
      annotation_package = db_name,
      annotation_package_version = as.character(packageVersion(db_name)),
      org_hs_eg_db_version = as.character(packageVersion("org.Hs.eg.db")),
      total_probes_in_annotation = length(probes),
      probes_queried = length(sample_probes),
      probes_with_entrez = length(probe_counts),
      probes_with_exactly_one_entrez = length(one_to_one_probes),
      distinct_entrez_queried = length(identifiers),
      entrez_with_ensembl = length(entrez_counts),
      entrez_with_exactly_one_ensembl = sum(entrez_counts == 1L),
      gencode_v49_intersection_performed = FALSE,
      crosswalk_frozen_for_transfer = FALSE
    )
  })
}
mapping_probe <- list(
  GPL570 = crosswalk_probe("GPL570", "hgu133plus2.db"),
  GPL16686 = crosswalk_probe("GPL16686", "hugene20sttranscriptcluster.db")
)
mapping_probe <- lapply(mapping_probe, function(result) {
  record <- list(status = result$status, error = result$error)
  if (identical(result$status, "pass")) c(record, result$value) else record
})

## --------------------------------------------------------------------- verdict
passed <- function(node) identical(node$status, "pass")
capabilities <- list(
  calvin_header_parse_both_platforms =
    passed(header_probe$GPL570) && passed(header_probe$GPL16686),
  GPL570_single_array_summarization =
    passed(summarizer_probe$GPL570_frma_hgu133plus2frmavecs) ||
      passed(summarizer_probe$GPL570_scan),
  GPL16686_single_array_summarization =
    passed(summarizer_probe$GPL16686_scan) ||
      passed(summarizer_probe$GPL16686_oligo_core_single_array),
  probe_to_entrez_to_ensembl_both_platforms =
    passed(mapping_probe$GPL570) && passed(mapping_probe$GPL16686)
)

session <- utils::capture.output(utils::sessionInfo())
writeLines(session, file.path(output_dir, "sessionInfo.txt"))

report <- list(
  schema_version = "masld-bench-microarray-single-array-capability-v1",
  status = if (all(unlist(capabilities))) "pass_single_array_runtime_capability" else "failed",
  R_version = R.version.string,
  R_platform = R.version$platform,
  R_home = R.home(),
  library_paths = paste(.libPaths(), collapse = ";"),
  package_versions = package_versions,
  calvin_header = header_probe,
  single_array_summarizers = summarizer_probe,
  probe_to_gene_mapping = mapping_probe,
  capabilities = capabilities,
  arrays_read_per_platform = 1L,
  all_sample_RMA_run = FALSE,
  across_array_quantile_normalization_run = FALSE,
  GEO_series_matrix_read = FALSE,
  expression_values_exported = FALSE,
  raw_probe_grid_rank_as_gene_expression_used = FALSE,
  gencode_v49_crosswalk_frozen = FALSE,
  labels_read = FALSE,
  model_training_activated = FALSE,
  sealed_outcomes_read = FALSE
)
writeLines(
  jsonlite::toJSON(report, auto_unbox = TRUE, null = "null", na = "string", pretty = TRUE),
  file.path(output_dir, "capability_probe.json")
)
cat("CAPABILITY_PROBE\t", report$status, "\n", sep = "")
if (!identical(report$status, "pass_single_array_runtime_capability")) {
  quit(status = 1L)
}
