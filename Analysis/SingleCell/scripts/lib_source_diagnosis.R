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

# Count only lineage cells in the eligible library columns actually present in
# this pseudobulk. The upstream export omits library-lineages with <10 cells.
source_lineage_cell_counts <- function(libraries, cell_type, contributing_runs, donor_map) {
  canonical <- function(x) gsub("[^A-Za-z0-9]", "_", x)
  key <- names(libraries)[canonical(names(libraries)) == canonical(paste0("n_", cell_type))]
  if (length(key) != 1L) stop("Missing or ambiguous lineage cell census: ", cell_type)
  if (anyDuplicated(contributing_runs)) stop("Duplicate pseudobulk library columns")
  j <- match(contributing_runs, libraries$sample)
  if (anyNA(j)) stop("Pseudobulk library is missing from the source cell census")
  donor <- ifelse(contributing_runs %in% names(donor_map), donor_map[contributing_runs], contributing_runs)
  if (!identical(as.character(donor), as.character(libraries$donor[j]))) {
    stop("Source census and pseudobulk donor map disagree")
  }
  n <- as.numeric(libraries[[key]][j])
  if (any(!is.finite(n) | n < 0 | n != floor(n))) stop("Invalid lineage cell counts")
  data.table::data.table(sample = donor, n_cells = n)[,
    .(n_cells = sum(n_cells), n_libraries = .N), by = sample]
}
