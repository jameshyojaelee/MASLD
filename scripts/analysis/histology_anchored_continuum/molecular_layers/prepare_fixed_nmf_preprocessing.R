#!/usr/bin/env Rscript
# Recover the original frozen NMF preprocessing; this does not produce corrected scores.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(jsonlite)})
root <- normalizePath(Sys.getenv("MASLD_PROJECT_ROOT", "."))
output <- Sys.getenv("FIXED_NMF_PREPROCESSING_OUTPUT", "")
if (!nzchar(output) || file.exists(output)) stop("A new FIXED_NMF_PREPROCESSING_OUTPUT is required")
sha <- function(path) strsplit(system2("sha256sum", shQuote(path), stdout = TRUE), "[[:space:]]+")[[1L]][[1L]]
paths <- c(
  cache = file.path(root, "RNA-seq/results/subtypes/nmf_results_cache_clean.rds"),
  dge = file.path(root, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds"),
  metadata = file.path(root, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"),
  qc = file.path(root, "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv"),
  preprocessing_script = file.path(root, "RNA-seq/95_nmf_clean_ksweep.R"),
  projection_script = file.path(root, "scripts/analysis/histology_anchored_continuum/molecular_layers/30_nmf_continuum.R")
)
hashes <- vapply(paths, sha, character(1))
cache <- readRDS(paths[["cache"]])
stopifnot(identical(names(cache), c("nmf_results", "metrics", "mat_nn", "mat")),
          identical(dimnames(cache$mat), dimnames(cache$mat_nn)))
metadata <- fread(paths[["metadata"]])
qc <- fread(paths[["qc"]])
metadata <- merge(metadata, qc[, .(sample_id, pass_technical)], by = "sample_id", sort = FALSE)
eligible <- metadata[group_binary == "Disease" & pass_technical == TRUE, sample_id]
stopifnot(setequal(eligible, colnames(cache$mat)), !anyDuplicated(metadata$sample_id))
dge <- readRDS(paths[["dge"]])
stopifnot(all(rownames(cache$mat) %in% rownames(dge)), all(colnames(cache$mat) %in% colnames(dge)))
dge <- dge[, colnames(cache$mat)]
dge <- calcNormFactors(dge, method = "TMM")
v <- voom(dge, design = NULL, plot = FALSE)
logcpm <- v$E[rownames(cache$mat), , drop = FALSE]
meta <- metadata[match(colnames(cache$mat), sample_id)]
reconstructed <- removeBatchEffect(logcpm, batch = meta$dataset)
matrix_difference <- max(abs(reconstructed - cache$mat))
stopifnot(is.finite(matrix_difference), matrix_difference < 1e-10)
minimum <- apply(cache$mat, 1, min)
shifted <- sweep(cache$mat, 1, minimum, "-")
shifted[shifted < 0] <- 0
nonnegative_difference <- max(abs(shifted - cache$mat_nn))
stopifnot(nonnegative_difference < 1e-12)
shift <- logcpm - cache$mat
cohorts <- sort(unique(meta$dataset))
batch_shift <- vapply(cohorts, function(cohort) rowMeans(shift[, meta$dataset == cohort, drop = FALSE]),
                     numeric(nrow(shift)))
rownames(batch_shift) <- rownames(cache$mat)
cohort_constancy <- max(abs(shift - batch_shift[, match(meta$dataset, cohorts), drop = FALSE]))
stopifnot(cohort_constancy < 1e-10)

# Use the existing NNLS function; do not refit bases or replace H with new scores.
previous_test_mode <- Sys.getenv("HAC_ML_TEST_MODE", unset = "")
Sys.setenv(HAC_ML_TEST_MODE = "1")
source(paths[["projection_script"]])
if (nzchar(previous_test_mode)) Sys.setenv(HAC_ML_TEST_MODE = previous_test_mode) else Sys.unsetenv("HAC_ML_TEST_MODE")
stopifnot(requireNamespace("NMF", quietly = TRUE))
projection_checks <- list()
sample_indices <- c(1L, as.integer(ncol(cache$mat_nn) / 2), ncol(cache$mat_nn))
for (k in c(4L, 6L)) {
  basis <- NMF::basis(cache$nmf_results[[as.character(k)]])
  stopifnot(identical(rownames(basis), rownames(cache$mat_nn)))
  colnames(basis) <- paste0("P", seq_len(k))
  coefficient <- nmf_nnls_coefficients(basis, cache$mat_nn[, sample_indices, drop = FALSE])
  stopifnot(identical(dim(coefficient), c(k, length(sample_indices))), all(is.finite(coefficient)), all(coefficient >= 0))
  projection_checks[[as.character(k)]] <- list(k = k, historical_samples = length(sample_indices),
    basis_refitted = FALSE, algorithm = "existing nmf_nnls_coefficients / nnls::nnls",
    finite_nonnegative = TRUE)
}
state <- list(schema_version = 1L, state = "historical_preprocessing_reproduced_not_corrected_projection",
  exact_gene_ids = rownames(cache$mat), historical_sample_ids = colnames(cache$mat),
  row_minimum = minimum, cohort_batch_subtraction = batch_shift,
  historical_TMM_sample_factors = setNames(dge$samples$norm.factors, colnames(dge)),
  historical_library_sizes = setNames(dge$samples$lib.size, colnames(dge)),
  source_paths = as.list(paths), source_sha256 = as.list(hashes))
dir.create(output, recursive = TRUE)
state_path <- file.path(output, "frozen_preprocessing_state.rds")
saveRDS(state, state_path)
fwrite(data.table(source_id = names(paths), path = unname(paths), sha256 = unname(hashes)),
       file.path(output, "source_manifest.tsv"), sep = "\t")
report <- list(status = "historical_preprocessing_reproduced_corrected_projection_unresolved",
  n_frozen_genes = nrow(cache$mat), n_historical_disease_libraries = ncol(cache$mat),
  maximum_original_matrix_difference = matrix_difference,
  maximum_nonnegative_transform_difference = nonnegative_difference,
  maximum_within_cohort_batch_shift_difference = cohort_constancy,
  recovered_frozen_parameters = c("row minima", "gene-by-cohort batch subtraction", "historical disease-library TMM factors and library sizes"),
  frozen_parameter_file_sha256 = sha(state_path), projection_checks = projection_checks,
  corrected_scores_produced = FALSE, bases_refitted = FALSE,
  unresolved_choices = c(
    "Original TMM population is 1,104 disease libraries. Qualified context includes controls without frozen training TMM factors; corrected-expression normalization population/reference must be specified before scoring.",
    "Corrected counting changes gene/locus memberships and library totals. Frozen basis requires exact gene IDs and all 5,000 genes; missing genes cannot be substituted with zero or folded by symbol.",
    "Full-W stage axes and signature-excluded continuum projections have different applicability. Neither may silently replace the other."))
writeLines(toJSON(report, auto_unbox = TRUE, pretty = TRUE), file.path(output, "preprocessing_review.json"))
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
cat(toJSON(report, auto_unbox = TRUE), "\n")
