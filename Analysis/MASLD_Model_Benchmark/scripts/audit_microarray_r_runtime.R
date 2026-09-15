#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: audit_microarray_r_runtime.R RUNTIME_ID GPL570.CEL GPL16686.CEL OUTPUT.tsv")
}

runtime_id <- args[[1L]]
fixture_paths <- c(GPL570 = args[[2L]], GPL16686 = args[[3L]])
output_path <- args[[4L]]
if (!all(file.exists(fixture_paths))) stop("CEL reader fixture is absent")
`%||%` <- function(left, right) if (is.null(left)) right else left

packages <- c(
  "affy", "affyio", "affxparser", "oligo", "frma", "SCAN.UPC",
  "hgu133plus2cdf", "hgu133plus2.db", "hgu133plus2frmavecs",
  "pd.hugene.2.0.st", "hugene20sttranscriptcluster.db",
  "AnnotationDbi", "org.Hs.eg.db"
)
installed <- installed.packages()
inventory <- lapply(packages, function(package) {
  list(
    package = package,
    installed = package %in% rownames(installed),
    version = if (package %in% rownames(installed)) unname(installed[package, "Version"]) else "not_installed"
  )
})
names(inventory) <- packages

reader <- list(
  package = "affyio",
  function_name = "read.celfile.header",
  function_present = FALSE,
  fixtures = list(
    GPL570 = list(status = "not_run", error = "reader_unavailable"),
    GPL16686 = list(status = "not_run", error = "reader_unavailable")
  )
)
if ("affyio" %in% rownames(installed)) {
  namespace <- asNamespace("affyio")
  reader$function_present <- exists("read.celfile.header", envir = namespace, inherits = FALSE)
  if (isTRUE(reader$function_present)) {
    header_reader <- get("read.celfile.header", envir = namespace, inherits = FALSE)
    reader$formals <- names(formals(header_reader))
    for (platform in names(fixture_paths)) {
      result <- tryCatch(
        {
          value <- header_reader(fixture_paths[[platform]])
          list(
            status = "pass",
            error = "not_applicable",
            class = paste(class(value), collapse = ";"),
            names = paste(names(value), collapse = ";"),
            dimensions = paste(dim(value), collapse = "x")
          )
        },
        error = function(error) list(status = "failed", error = conditionMessage(error))
      )
      reader$fixtures[[platform]] <- result
    }
  }
}

has <- function(package) package %in% rownames(installed)
capabilities <- list(
  Calvin_header_reader_both_platforms = isTRUE(reader$function_present) &&
    all(vapply(reader$fixtures, function(value) identical(value$status, "pass"), logical(1))),
  GPL570_exact_frma = all(vapply(c("frma", "hgu133plus2frmavecs"), has, logical(1))),
  GPL570_SCAN = all(vapply(c("SCAN.UPC", "hgu133plus2cdf"), has, logical(1))),
  GPL570_affy_summary = all(vapply(c("affy", "hgu133plus2cdf"), has, logical(1))),
  GPL16686_SCAN = all(vapply(c("SCAN.UPC", "oligo", "pd.hugene.2.0.st"), has, logical(1))),
  GPL16686_oligo_core_summary = all(vapply(c("oligo", "pd.hugene.2.0.st"), has, logical(1))),
  GPL16686_transcript_cluster_annotation = has("hugene20sttranscriptcluster.db")
)
capabilities$GPL570_expression_rank_baseline <- capabilities$GPL570_exact_frma ||
  capabilities$GPL570_SCAN || capabilities$GPL570_affy_summary
capabilities$GPL16686_expression_rank_baseline <- capabilities$GPL16686_SCAN ||
  capabilities$GPL16686_oligo_core_summary

status <- if (isTRUE(capabilities$Calvin_header_reader_both_platforms)) {
  "pass_Calvin_header_reader_biological_summarization_blocked"
} else {
  "blocked_no_executable_Calvin_reader_or_biological_summarizer"
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
records <- character()
add_record <- function(key, value) {
  if (length(key) != 1L || length(value) != 1L || is.na(value) || !nzchar(as.character(value))) {
    stop("runtime audit record must contain one nonempty key and one explicit value")
  }
  if (key %in% names(records)) stop("runtime audit key collision")
  records[[key]] <<- as.character(value)
}
add_record("schema_version", "masld-bench-microarray-R-runtime-audit-v1")
add_record("status", status)
add_record("runtime_id", runtime_id)
add_record("R_version", R.version.string)
add_record("R_home", R.home())
add_record("library_paths", paste(.libPaths(), collapse = ";"))
add_record("reader.package", reader$package)
add_record("reader.function_name", reader$function_name)
add_record("reader.function_present", tolower(as.character(reader$function_present)))
add_record("reader.formals", if (length(reader$formals %||% character())) paste(reader$formals, collapse = ";") else "not_available")
for (package in packages) {
  add_record(paste0("package.", package, ".installed"), tolower(as.character(inventory[[package]]$installed)))
}
for (package in packages) {
  add_record(paste0("package.", package, ".version"), inventory[[package]]$version)
}
for (platform in c("GPL570", "GPL16686")) {
  add_record(paste0("reader.", platform, ".status"), reader$fixtures[[platform]]$status)
}
for (platform in c("GPL570", "GPL16686")) {
  add_record(paste0("reader.", platform, ".error"), reader$fixtures[[platform]]$error)
}
for (capability in names(capabilities)) {
  add_record(paste0("capability.", capability), tolower(as.character(capabilities[[capability]])))
}
for (key in c(
  "raw_probe_grid_rank_as_gene_expression_allowed",
  "ordinary_all_sample_RMA_run", "quantile_normalization_run",
  "CEL_expression_values_exported", "model_training_activated", "labels_read"
)) {
  add_record(key, "false")
}
rows <- data.frame(key = names(records), value = unname(records), stringsAsFactors = FALSE)
if (nrow(rows) != 55L) stop("runtime audit fixed schema must contain exactly 55 records")
write.table(rows, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("RUNTIME_AUDIT\t", runtime_id, "\t", status, "\n", sep = "")
