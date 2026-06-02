#!/usr/bin/env Rscript
# ============================================================================
# run_ucell_F_transitions.R
# ----------------------------------------------------------------------------
# UCell sanity-check cross-validation for scDRS F-transition scores.
# Andreatta & Carmona, CSBJ 2021, PMID 34285779.
#
# - Loads .gs file written by run_scdrs_F_transitions.py.
# - Loads integrated_atlas.h5ad (raw counts) via zellkonverter.
# - Restricts to preparation_method in {unsorted, nuclei}.
# - For each F transition, builds signed signature (up genes positive, down
#   genes negative — UCell handles signed via leading "-" prefix).
# - Computes UCell score per cell with ScoreSignatures_UCell (chunked by
#   cell-type to bound memory; w_neg=1, maxRank=1500, ncores=8).
# - Aggregates median per cell type.
# - Compares against scDRS group-level enrichment (median z per cell-type)
#   from scdrs_cell_scores.parquet — Spearman rho across (cell_type ×
#   transition) pairs.
#
# Outputs:
#   ucell_F_transitions_scores.csv         (per-cell UCell scores, long form)
#   ucell_celltype_summary.csv             (median per cell-type × transition)
#   scdrs_vs_ucell_concordance.csv         (paired comparison + Spearman)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Matrix)
  library(UCell)
})

t0 <- Sys.time()
log_msg <- function(...) cat(sprintf("[ucell %s] %s\n",
                                     format(Sys.time(), "%H:%M:%S"),
                                     paste0(..., collapse = "")))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DEG_DIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/disease_signatures")
ATLAS_H5 <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad")
UMAP_CSV <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/atlas_umap_for_fig2.csv.gz")

GS_FILE          <- file.path(DEG_DIR, "F_transitions_scdrs.gs")
SCDRS_CELL_SCORE <- file.path(DEG_DIR, "scdrs_cell_scores.parquet")
SCDRS_CT_ENRICH  <- file.path(DEG_DIR, "scdrs_celltype_enrichment.csv")

OUT_CELL <- file.path(DEG_DIR, "ucell_F_transitions_scores.csv")
OUT_CT   <- file.path(DEG_DIR, "ucell_celltype_summary.csv")
OUT_CONC <- file.path(DEG_DIR, "scdrs_vs_ucell_concordance.csv")

PREP_KEEP <- c("unsorted", "nuclei")
TRANSITIONS <- c("F1_vs_F0", "F2_vs_F1", "F3_vs_F2", "F4_vs_F3")
NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
MAXRANK <- 1500

log_msg("ncores=", NCORES, "  maxRank=", MAXRANK)

# ----------------------------------------------------------------------------
# Parse .gs to build signed UCell signatures
# ----------------------------------------------------------------------------
log_msg("Loading .gs from ", GS_FILE)
gs_lines <- readLines(GS_FILE)
gs_lines <- gs_lines[nchar(gs_lines) > 0 & !startsWith(gs_lines, "TRAIT")]

parse_gs <- function(line) {
  parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
  trait <- parts[1]
  payload <- parts[2]
  pairs <- strsplit(payload, ",", fixed = TRUE)[[1]]
  m <- do.call(rbind, strsplit(pairs, ":", fixed = TRUE))
  list(trait = trait, genes = m[, 1], weights = as.numeric(m[, 2]))
}

gs_list <- lapply(gs_lines, parse_gs)
names(gs_list) <- sapply(gs_list, `[[`, "trait")
log_msg("Parsed ", length(gs_list), " .gs traits")

# Build UCell signed signatures: list per transition (signed mode = up + "-down")
ucell_sigs <- list()
for (tr in TRANSITIONS) {
  up_key <- paste0(tr, "__up")
  dn_key <- paste0(tr, "__down")
  up_g <- if (up_key %in% names(gs_list)) gs_list[[up_key]]$genes else character(0)
  dn_g <- if (dn_key %in% names(gs_list)) gs_list[[dn_key]]$genes else character(0)
  if (length(up_g) == 0 && length(dn_g) == 0) next
  ucell_sigs[[tr]] <- c(up_g, paste0("-", dn_g))
}
log_msg("UCell signatures: ", paste(names(ucell_sigs), sapply(ucell_sigs, length),
                                    sep = "=", collapse = ", "))

# ----------------------------------------------------------------------------
# Load atlas; use raw counts (assayName='counts' from .raw via zellkonverter)
# ----------------------------------------------------------------------------
log_msg("Loading h5ad via zellkonverter (use_hdf5=FALSE)")
sce <- readH5AD(ATLAS_H5, use_hdf5 = FALSE, reader = "python", verbose = FALSE)
log_msg("SCE: ", ncol(sce), " cells x ", nrow(sce), " genes;  assays=",
        paste(assayNames(sce), collapse = ","))

# .raw becomes assayName 'X' or 'raw_X' depending on zellkonverter mode; pick the
# layer that holds counts (integers). If multiple, prefer 'counts' or 'raw_X'.
layer_to_use <- NULL
for (an in assayNames(sce)) {
  m <- assay(sce, an)
  if (is(m, "dgCMatrix") || is(m, "dgRMatrix") || is(m, "dgTMatrix")) {
    sample_vals <- m@x[seq_len(min(1000, length(m@x)))]
    if (all(sample_vals == round(sample_vals))) {
      layer_to_use <- an
      break
    }
  }
}
if (is.null(layer_to_use)) layer_to_use <- assayNames(sce)[1]
log_msg("Using assay layer for UCell: ", layer_to_use)
counts_mat <- assay(sce, layer_to_use)

# Set rownames to var_names (already done by zellkonverter)
rownames(sce) <- rownames(counts_mat)

# ----------------------------------------------------------------------------
# Apply prep filter via UMAP CSV (atlas obs lacks preparation_method)
# ----------------------------------------------------------------------------
log_msg("Loading UMAP meta and applying positional alignment ",
        "(atlas_umap_for_fig2.csv.gz has 8 cols, no cell_id; same row order ",
        "as integrated_atlas.h5ad obs)")
umap_dt <- fread(UMAP_CSV)
log_msg("  UMAP meta cols: ", paste(names(umap_dt), collapse = ", "))
if (nrow(umap_dt) != ncol(sce)) {
  stop("UMAP rows (", nrow(umap_dt), ") != SCE cells (", ncol(sce), ")")
}
# Defensive consistency check
ct_match <- mean(as.character(umap_dt$cell_type) ==
                 as.character(colData(sce)$cell_type))
log_msg("  positional cell_type match: ", sprintf("%.4f", ct_match))
if (ct_match < 0.99) stop("UMAP CSV does not align positionally with atlas")
colData(sce)$preparation_method <- as.character(umap_dt$preparation_method)
if ("disease_stage_coarse" %in% names(umap_dt)) {
  colData(sce)$disease_stage_coarse <- as.character(umap_dt$disease_stage_coarse)
}

keep <- colData(sce)$preparation_method %in% PREP_KEEP
log_msg("Prep filter: keep ", sum(keep, na.rm = TRUE), " / ", ncol(sce))
sce <- sce[, which(keep)]
counts_mat <- assay(sce, layer_to_use)

# ----------------------------------------------------------------------------
# UCell scoring per cell-type chunk (memory-bounded)
# ----------------------------------------------------------------------------
ct_vec <- as.character(colData(sce)$cell_type)
ct_levels <- unique(ct_vec)
log_msg("Cell types: ", paste(ct_levels, collapse = ", "))

# Score each cell-type chunk separately to avoid building a huge dense
# rank matrix in memory.
ucell_chunks <- list()
for (ct in ct_levels) {
  idx <- which(ct_vec == ct)
  if (length(idx) == 0) next
  log_msg("  UCell ", ct, " n=", length(idx))
  m <- counts_mat[, idx, drop = FALSE]
  scores <- ScoreSignatures_UCell(
    matrix = m,
    features = ucell_sigs,
    maxRank = MAXRANK,
    w_neg = 1,
    chunk.size = 5000,
    ncores = NCORES,
    name = ""
  )
  # scores: rows = cells, cols = signature names (no '_UCell' suffix because name="")
  scores_dt <- as.data.table(scores, keep.rownames = "cell_id")
  scores_dt[, cell_type := ct]
  ucell_chunks[[ct]] <- scores_dt
  rm(scores, m); gc(verbose = FALSE)
}

ucell_dt <- rbindlist(ucell_chunks, use.names = TRUE, fill = TRUE)
log_msg("UCell scored cells: ", nrow(ucell_dt))

# Long form
val_cols <- setdiff(names(ucell_dt), c("cell_id", "cell_type"))
ucell_long <- melt(ucell_dt, id.vars = c("cell_id", "cell_type"),
                   measure.vars = val_cols,
                   variable.name = "transition", value.name = "ucell_score")
ucell_long[, transition := as.character(transition)]
fwrite(ucell_long, OUT_CELL)
log_msg("Wrote per-cell UCell: ", OUT_CELL, " (", nrow(ucell_long), " rows)")

# Per cell-type median
ucell_ct <- ucell_long[, .(
  ucell_median = median(ucell_score, na.rm = TRUE),
  ucell_mean   = mean(ucell_score, na.rm = TRUE),
  n_cells      = .N
), by = .(cell_type, transition)]
fwrite(ucell_ct, OUT_CT)
log_msg("Wrote per-celltype UCell: ", OUT_CT)

# ----------------------------------------------------------------------------
# Compare against scDRS per cell-type enrichment
# ----------------------------------------------------------------------------
log_msg("Loading scDRS celltype enrichment for comparison")
if (!file.exists(SCDRS_CT_ENRICH)) {
  log_msg("ERROR: ", SCDRS_CT_ENRICH, " not found; run scDRS first")
  quit(status = 0)
}
scdrs_ct <- fread(SCDRS_CT_ENRICH)
# We compare UCell-signed (transition-level) vs scDRS-signed (mode=='signed')
scdrs_signed <- scdrs_ct[mode == "signed"]
# Build NES proxy: -log10(assoc_mcp + 1e-4) signed by sign of (n_fdr_0.05 - n_cells/2 * (n_fdr_0.05/n_cells))
# Simpler: use 'assoc_zscore' if present, else -log10(assoc_mcp).
nes_col <- if ("assoc_zscore" %in% names(scdrs_signed)) "assoc_zscore" else "assoc_mcp"
log_msg("scDRS NES proxy column: ", nes_col)

if (nes_col == "assoc_mcp") {
  scdrs_signed[, scdrs_nes := -log10(pmax(assoc_mcp, 1e-4))]
} else {
  scdrs_signed[, scdrs_nes := get(nes_col)]
}

merged <- merge(
  scdrs_signed[, .(cell_type, transition, scdrs_nes,
                   assoc_mcp = if ("assoc_mcp" %in% names(scdrs_signed)) assoc_mcp else NA_real_)],
  ucell_ct[, .(cell_type, transition, ucell_median, ucell_mean, n_cells)],
  by = c("cell_type", "transition"),
  all = TRUE
)
fwrite(merged, OUT_CONC)

# Spearman rho across all (ct x transition) pairs
ok <- !is.na(merged$scdrs_nes) & !is.na(merged$ucell_median)
rho_all <- if (sum(ok) >= 5) cor(merged$scdrs_nes[ok], merged$ucell_median[ok], method = "spearman") else NA_real_
log_msg("Spearman rho (scDRS NES vs UCell median, all pairs): ",
        sprintf("%.3f", rho_all), "  n=", sum(ok))

# Per-transition rho
for (tr in TRANSITIONS) {
  sub <- merged[transition == tr & !is.na(scdrs_nes) & !is.na(ucell_median)]
  rho <- if (nrow(sub) >= 5) cor(sub$scdrs_nes, sub$ucell_median, method = "spearman") else NA_real_
  log_msg("  ", tr, ": rho=", sprintf("%.3f", rho), "  n=", nrow(sub))
}

log_msg("DONE in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
