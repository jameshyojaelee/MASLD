#!/usr/bin/env Rscript
# 89_toast_celltype_de.R
#
# Analysis A1 v2 (partial) — TOAST-based cell-type-specific DE from bulk.
#
# TOAST (Tools for the Analysis of mixed Samples in omicS data) performs
# reference-based CT-specific DE from bulk RNA-seq using scRNA-derived
# cell-type proportions.
#
# Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md
#
# NOTE: CARseq, bMIND, speckle, MOFA2, CellChat, nichenetr, multinichenetr
# remain unavailable (conda env install + GitHub rate limit failures).
# TOAST partially fills the A1 v2 goal — per-cell-type DE from bulk data.
#
# Env: rnaseq (has TOAST)
# Outputs: RNA-seq/results/celltype_attribution/

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(TOAST)
  library(edgeR)
  library(limma)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
DECONV <- file.path(BASE, "Analysis/Deconvolution/results")
OUTDIR <- file.path(BASE, "RNA-seq/results/celltype_attribution")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("[1] Loading merged DGE + metadata...")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
ycfg <- yaml::read_yaml(file.path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design", "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
keep <- dge$samples$dataset %in% mega_cohorts
dge <- dge[, keep]

message("[2] Loading MuSiC proportions for ALL cell types...")
datasets <- unique(dge$samples$dataset)
prop_list <- lapply(datasets, function(ds) {
  f <- file.path(DECONV, ds, paste0(ds, "_music_prop_weighted.tsv"))
  if (!file.exists(f)) return(NULL)
  df <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE)
  df$sample_id <- rownames(df)
  setDT(df)
  df
})
prop <- rbindlist(prop_list, fill = TRUE, use.names = TRUE)
ct_cols <- setdiff(names(prop), "sample_id")
message(sprintf("  Cell types: %s", paste(ct_cols, collapse = ", ")))

# Match to samples in dge
prop_mat <- as.matrix(prop[match(colnames(dge), sample_id), ..ct_cols])
# Clean: impute missing proportions with column median
for (j in seq_len(ncol(prop_mat))) {
  na_idx <- is.na(prop_mat[, j])
  if (any(na_idx)) prop_mat[na_idx, j] <- median(prop_mat[, j], na.rm = TRUE)
}
# Renormalize rows to sum to 1
prop_mat <- prop_mat / rowSums(prop_mat)

# Group
grp <- factor(dge$samples$group_binary, levels = c("Control","Disease"))
dataset <- factor(dge$samples$dataset)

# Keep only cell types with mean proportion > 0.01 (otherwise unstable)
ct_means <- colMeans(prop_mat, na.rm = TRUE)
ct_keep <- names(ct_means)[ct_means > 0.01]
message(sprintf("  Cell types kept (mean prop > 1%%): %s",
                paste(ct_keep, collapse = ", ")))
prop_mat <- prop_mat[, ct_keep, drop = FALSE]
# Sanitize column names for TOAST (no spaces/special chars)
ct_sanitize_map <- setNames(ct_keep, make.names(ct_keep))
colnames(prop_mat) <- names(ct_sanitize_map)
ct_keep_sanitized <- colnames(prop_mat)

message("[3] Filter genes + voom transform...")
keep_genes <- filterByExpr(dge, group = grp)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
v <- voom(dge, design = model.matrix(~ grp + dataset), plot = FALSE)
Y <- v$E  # log2 CPM
message(sprintf("  Genes: %d x Samples: %d", nrow(Y), ncol(Y)))

message("[4] TOAST cell-type-specific DE (disease vs control)...")
# Build design matrix (sample-level) — TOAST csDE needs Y (Gene x Sample),
# Prop (Sample x CellType), and design_data
# TOAST requires design_info to NOT contain variables that are collinear with
# the cell-type proportions. Dropping dataset causes singular solve; use group
# only. Cell-type proportions sum to 1 — TOAST handles this internally.
design_info <- data.frame(group = grp)
rownames(design_info) <- colnames(Y)

# Drop smallest proportion cell type to avoid perfect collinearity (proportions
# sum to 1). This is a standard precaution.
keep_n <- ncol(prop_mat) - 1
ct_order <- order(colMeans(prop_mat, na.rm = TRUE), decreasing = TRUE)
prop_mat <- prop_mat[, ct_order[1:keep_n], drop = FALSE]
# Renormalize rows
prop_mat <- prop_mat / rowSums(prop_mat)
ct_keep_sanitized <- colnames(prop_mat)
message(sprintf("  TOAST using cell types: %s",
                paste(ct_keep_sanitized, collapse = ", ")))

design_out <- makeDesign(design_info, prop_mat)

# Fit full model
fit <- fitModel(design_out, Y)

# For each cell type, test cell-type-specific disease effect
results_all <- rbindlist(lapply(ct_keep_sanitized, function(ct) {
  tryCatch({
    # TOAST syntax: coef = "group" (the group variable name)
    res <- csTest(fit, coef = "group", cell_type = ct,
                  contrast_matrix = NULL, verbose = FALSE)
    dt <- as.data.table(res, keep.rownames = "gene")
    dt[, cell_type := ct]
    dt[, cell_type_original := ct_sanitize_map[[ct]]]
    dt
  }, error = function(e) {
    message(sprintf("  Failed for %s: %s", ct, e$message))
    NULL
  })
}), fill = TRUE)

fwrite(results_all, file.path(OUTDIR, "toast_celltype_de_all.csv"))
message(sprintf("  TOAST results rows: %d (across %d cell types)",
                nrow(results_all), uniqueN(results_all$cell_type)))

# Per-cell-type DEG counts
if ("fdr" %in% names(results_all)) {
  deg_counts <- results_all[, .(n_sig_padj05 = sum(fdr < 0.05, na.rm = TRUE),
                                 n_sig_padj10 = sum(fdr < 0.10, na.rm = TRUE)),
                             by = cell_type][order(-n_sig_padj05)]
} else if ("adj.P.Val" %in% names(results_all)) {
  deg_counts <- results_all[, .(n_sig_padj05 = sum(adj.P.Val < 0.05, na.rm = TRUE),
                                 n_sig_padj10 = sum(adj.P.Val < 0.10, na.rm = TRUE)),
                             by = cell_type][order(-n_sig_padj05)]
} else {
  # Fallback to raw pvalue
  pval_col <- intersect(c("P.Value","p.value","pvalue","pval"), names(results_all))[1]
  if (!is.na(pval_col)) {
    results_all[, bh_fdr := p.adjust(get(pval_col), method = "BH"), by = cell_type]
    deg_counts <- results_all[, .(n_sig_padj05 = sum(bh_fdr < 0.05, na.rm = TRUE),
                                   n_sig_padj10 = sum(bh_fdr < 0.10, na.rm = TRUE)),
                               by = cell_type][order(-n_sig_padj05)]
  } else deg_counts <- data.table()
}
fwrite(deg_counts, file.path(OUTDIR, "toast_celltype_deg_counts.csv"))

summary_lines <- c(
  sprintf("TOAST cell-type-specific DE from bulk"),
  sprintf("Samples: %d; cell types tested: %d", ncol(Y), length(ct_keep)),
  sprintf("Genes: %d", nrow(Y)),
  "",
  "=== Cell types tested ===",
  capture.output(print(data.table(cell_type = ct_keep,
                                   mean_prop = ct_means[ct_keep]))),
  "",
  "=== DEG counts per cell type (TOAST) ===",
  capture.output(print(deg_counts)))
writeLines(summary_lines, file.path(OUTDIR, "toast_celltype_de_summary.txt"))
writeLines(summary_lines)

message("Done. Outputs in: ", OUTDIR)
