#!/usr/bin/env Rscript
# Build a genome-wide, donor-collapsed baseline lineage-expression reference.
# This replaces the run-weighted legacy cell-origin table for Plan 45 routing.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

project_root <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
candidate_id <- Sys.getenv(
  "PLAN45_CANDIDATE_ID",
  unset = "source-independent-risk-state-relay-lineage-observability-2026-08-10"
)
if (!grepl("^source-independent-risk-state-relay-[A-Za-z0-9._-]+$", candidate_id)) {
  stop("Invalid PLAN45_CANDIDATE_ID")
}
candidate_root <- file.path(
  project_root,
  "Analysis/Multimodal_Program_Projection/candidates",
  candidate_id
)
if (dir.exists(candidate_root)) stop("Refusing to overwrite candidate: ", candidate_root)

pb_dir <- file.path(project_root, "Analysis/SingleCell/results_gpu_v2/pseudobulk")
collapse_library <- file.path(project_root, "Analysis/SingleCell/scripts/lib_donor_collapse.R")
if (!dir.exists(pb_dir) || !file.exists(collapse_library)) stop("Lineage source inputs absent")
source(collapse_library)

pb_files <- sort(list.files(pb_dir, pattern = "_pseudobulk\\.csv$", full.names = TRUE))
expected_cell_types <- sort(c(
  "B_cells", "Basophils", "Cholangiocytes", "Circulating_NK_NKT",
  "Endothelial_cells", "Fibroblasts", "Hepatocytes", "Macrophages",
  "Mono+mono_derived_cells", "Neutrophils", "Plasma_cells", "Resident_NK",
  "T_cells", "cDC1s", "cDC2s", "pDCs"
))
observed_cell_types <- sort(sub("_pseudobulk\\.csv$", "", basename(pb_files)))
if (!identical(observed_cell_types, expected_cell_types)) {
  stop(
    "Cell-type pseudobulk inventory drift; missing=",
    paste(setdiff(expected_cell_types, observed_cell_types), collapse = ";"),
    " unexpected=",
    paste(setdiff(observed_cell_types, expected_cell_types), collapse = ";")
  )
}

sha256 <- function(path) {
  output <- system2("sha256sum", args = path, stdout = TRUE, stderr = TRUE)
  if (length(output) != 1L) stop("sha256sum failed for ", path)
  strsplit(output, "[[:space:]]+")[[1L]][[1L]]
}

atomic_fwrite <- function(x, path) {
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "")
  if (!file.rename(tmp, path)) stop("Atomic rename failed for ", path)
}

safe_cell_type <- function(path) {
  sub("_pseudobulk\\.csv$", "", basename(path))
}

srr_to_donor <- build_srr_to_donor_map(project_root)
if (length(srr_to_donor) == 0L) stop("Donor-collapse map is empty")

long_parts <- vector("list", length(pb_files))
source_audit <- vector("list", length(pb_files))
for (index in seq_along(pb_files)) {
  path <- pb_files[[index]]
  cell_type <- safe_cell_type(path)
  message("[lineage] ", cell_type, " <- ", basename(path))
  counts_dt <- fread(path, check.names = FALSE)
  if (ncol(counts_dt) < 3L) stop("Pseudobulk matrix has fewer than two samples: ", path)
  genes <- as.character(counts_dt[[1L]])
  counts <- as.matrix(counts_dt[, -1L, with = FALSE])
  storage.mode(counts) <- "numeric"
  rownames(counts) <- genes
  if (anyDuplicated(genes)) {
    counts <- rowsum(counts, group = genes, reorder = FALSE)
    genes <- rownames(counts)
  }
  run_ids <- colnames(counts)
  donor_ids <- ifelse(
    run_ids %in% names(srr_to_donor),
    srr_to_donor[run_ids],
    run_ids
  )
  collapsed <- collapse_counts_to_donor(counts, project_root)
  if (anyDuplicated(colnames(collapsed))) stop("Donor collapse left duplicate columns: ", path)
  library_size <- colSums(collapsed)
  usable <- is.finite(library_size) & library_size > 0
  if (sum(usable) < 5L) stop("Fewer than five usable donors for ", cell_type)
  cpm <- sweep(collapsed[, usable, drop = FALSE], 2L, library_size[usable], "/") * 1e6
  long_parts[[index]] <- data.table(
    gene_symbol = rownames(cpm),
    cell_type = cell_type,
    n_donors = ncol(cpm),
    mean_cpm = rowMeans(cpm),
    median_cpm = apply(cpm, 1L, median),
    donor_detection_fraction_cpm1 = rowMeans(cpm >= 1),
    donor_detection_fraction_nonzero = rowMeans(cpm > 0)
  )
  source_audit[[index]] <- data.table(
    cell_type = cell_type,
    source_path = sub(paste0("^", project_root, "/"), "", path),
    n_gene_rows = nrow(cpm),
    n_run_columns = length(run_ids),
    n_donor_columns = ncol(collapsed),
    n_usable_donors = ncol(cpm),
    n_run_ids_in_authoritative_pairing = sum(run_ids %in% names(srr_to_donor)),
    n_collapsed_columns = length(run_ids) - length(unique(donor_ids)),
    source_sha256 = sha256(path)
  )
  rm(counts_dt, counts, collapsed, cpm)
  gc(verbose = FALSE)
}

long <- rbindlist(long_parts, use.names = TRUE)
audit <- rbindlist(source_audit, use.names = TRUE)
if (uniqueN(long$cell_type) != 16L) stop("Lineage reference lost a cell type")

wide <- dcast(long, gene_symbol ~ cell_type, value.var = "mean_cpm", fill = 0)
genes <- wide$gene_symbol
matrix_values <- as.matrix(wide[, -1L, with = FALSE])
cell_types <- colnames(matrix_values)
row_total <- rowSums(matrix_values)
row_max <- apply(matrix_values, 1L, max)
dominant_index <- max.col(matrix_values, ties.method = "first")
dominant_celltype <- cell_types[dominant_index]
second_cpm <- apply(matrix_values, 1L, function(values) {
  if (length(values) < 2L) return(NA_real_)
  sort(values, partial = length(values) - 1L)[length(values) - 1L]
})
tau <- rowSums(1 - matrix_values / ifelse(row_max > 0, row_max, NA_real_)) /
  (ncol(matrix_values) - 1L)
summary <- data.table(
  gene_symbol = genes,
  n_celltypes = ncol(matrix_values),
  total_mean_cpm = row_total,
  dominant_celltype = dominant_celltype,
  dominant_mean_cpm = row_max,
  second_mean_cpm = second_cpm,
  dominant_to_second_ratio = fifelse(
    second_cpm > 0,
    row_max / second_cpm,
    fifelse(row_max > 0, Inf, NA_real_)
  ),
  dominant_fraction = fifelse(row_total > 0, row_max / row_total, NA_real_),
  specificity_tau = tau
)
setkey(long, gene_symbol, cell_type)
summary[, dominant_detection_fraction_cpm1 := long[
  summary,
  on = .(gene_symbol, cell_type = dominant_celltype),
  donor_detection_fraction_cpm1
]]
summary[, lineage_expression_gate :=
  dominant_mean_cpm >= 1 & dominant_detection_fraction_cpm1 >= 0.25]
summary[, lineage_specificity_gate := fifelse(
  lineage_expression_gate & !is.na(dominant_to_second_ratio) &
    dominant_to_second_ratio >= 2,
  TRUE,
  FALSE
)]
summary[, interpretation := fifelse(
  lineage_specificity_gate,
  "candidate_source_lineage_requires_accessibility_validation",
  fifelse(
    lineage_expression_gate,
    "multilineage_expression_no_unique_source",
    "lineage_expression_insufficient"
  )
)]

dir.create(candidate_root, recursive = TRUE)
long_path <- file.path(candidate_root, "donor_lineage_observability.tsv")
summary_path <- file.path(candidate_root, "donor_lineage_gene_summary.tsv")
audit_path <- file.path(candidate_root, "lineage_source_audit.tsv")
atomic_fwrite(long, long_path)
atomic_fwrite(summary, summary_path)
atomic_fwrite(audit, audit_path)

pairing_files <- file.path(project_root, .DONOR_PAIRING_FILES)
manifest <- rbindlist(list(
  data.table(
    source_id = paste0("pseudobulk_", audit$cell_type),
    source_path = audit$source_path,
    size_bytes = file.info(pb_files)$size,
    sha256 = audit$source_sha256,
    role = "run-level raw-count pseudobulk; donor collapsed in this candidate"
  ),
  data.table(
    source_id = paste0("donor_pairing_", names(.DONOR_PAIRING_FILES)),
    source_path = .DONOR_PAIRING_FILES,
    size_bytes = file.info(pairing_files)$size,
    sha256 = vapply(pairing_files, sha256, character(1L)),
    role = "authoritative run-to-biological-donor mapping"
  ),
  data.table(
    source_id = "donor_collapse_implementation",
    source_path = sub(paste0("^", project_root, "/"), "", collapse_library),
    size_bytes = file.info(collapse_library)$size,
    sha256 = sha256(collapse_library),
    role = "shared donor-collapse implementation"
  )
), use.names = TRUE)
manifest_path <- file.path(candidate_root, "lineage_source_manifest.tsv")
atomic_fwrite(manifest, manifest_path)

seal <- list(
  status = "donor_collapsed_lineage_observability_reference",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  candidate_id = candidate_id,
  n_celltypes = uniqueN(long$cell_type),
  n_genes = uniqueN(long$gene_symbol),
  n_long_rows = nrow(long),
  lineage_gate_definition = paste(
    "dominant mean CPM >=1 and donor detection fraction CPM>=1 >=0.25;",
    "specificity additionally requires dominant/second mean CPM >=2"
  ),
  disease_or_genetic_outcomes_inspected = FALSE,
  experimental_targets_frozen = FALSE,
  known_boundary = paste(
    "baseline snRNA expression can route a lineage but cannot establish",
    "regulatory-element accessibility or cell-autonomous action"
  ),
  output_sha256 = list(
    donor_lineage_observability = sha256(long_path),
    donor_lineage_gene_summary = sha256(summary_path),
    lineage_source_audit = sha256(audit_path),
    lineage_source_manifest = sha256(manifest_path)
  )
)
seal_path <- file.path(candidate_root, "LINEAGE_REFERENCE_SEALED.json")
tmp_seal <- tempfile(pattern = ".LINEAGE_REFERENCE_SEALED.", tmpdir = candidate_root)
writeLines(toJSON(seal, auto_unbox = TRUE, pretty = TRUE), tmp_seal)
if (!file.rename(tmp_seal, seal_path)) stop("Failed to seal lineage reference")
message(
  "[lineage] sealed ", uniqueN(long$gene_symbol), " genes x ",
  uniqueN(long$cell_type), " cell types; targets not frozen"
)
