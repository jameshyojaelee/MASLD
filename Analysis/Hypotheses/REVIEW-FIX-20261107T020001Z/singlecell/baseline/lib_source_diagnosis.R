# Source-qualified RNA labels. Historical release tables remain read-only.
source_diagnosis_root <- function(base) Sys.getenv(
  "SINGLECELL_DIAGNOSIS_ROOT",
  file.path(base, "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/metadata")
)

load_source_diagnosis <- function(base, level = c("donor", "library")) {
  level <- match.arg(level)
  path <- file.path(source_diagnosis_root(base), paste0(level, "_metadata.tsv"))
  corrected <- file.path(source_diagnosis_root(base), "library_metadata_corrected.tsv")
  if (level == "library" && file.exists(corrected)) path <- corrected
  measured <- file.path(source_diagnosis_root(base), "library_metadata_source_F.tsv")
  if (level == "library" && file.exists(measured)) path <- measured
  if (!file.exists(path)) stop("Missing source-qualified RNA metadata: ", path)
  x <- data.table::fread(path)
  key <- if (level == "donor") "donor" else "sample"
  if (anyDuplicated(x[[key]])) stop("Ambiguous source-qualified RNA key: ", key)
  x
}
