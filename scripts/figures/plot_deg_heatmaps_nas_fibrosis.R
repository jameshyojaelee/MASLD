#!/usr/bin/env Rscript
# plot_deg_heatmaps_nas_fibrosis.R
# Generate heatmaps of logFC for DEGs across NAS scores and fibrosis stages.
# Handles All DEGs and Stage-Unique DEGs.
# Includes control group (NAS0/F0) with zero logFC, Viridis palette, and K-Means clustering.
# Uses anno_mark to label 10 important MASLD genes.

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(viridis)
})

cat("Starting DEG Heatmap generation with updated aesthetics...\n")

BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- file.path(FIG2_DIR, "panels")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PADJ_THRESH <- 0.1
LFC_THRESH <- 0.5

nas_dt <- load_nas_score_progression()
fib_dt <- load_fibrosis_stage_dream()

# Helper: use gene symbols as row names (deduplicate by keeping highest variance)
use_symbols <- function(mat_dt, dt_with_symbols) {
  # Map gene (Ensembl) → symbol
  sym_map <- unique(dt_with_symbols[, .(gene, symbol)])
  mat_dt <- merge(mat_dt, sym_map, by = "gene", all.x = TRUE)
  mat_dt[is.na(symbol), symbol := gene]
  mat_dt[, gene := NULL]

  # Deduplicate symbols: keep row with highest variance across columns
  num_cols <- setdiff(names(mat_dt), "symbol")
  mat_dt[, .row_var := apply(.SD, 1, var, na.rm = TRUE), .SDcols = num_cols]
  mat_dt <- mat_dt[mat_dt[, .I[which.max(.row_var)], by = symbol]$V1]
  mat_dt[, .row_var := NULL]
  mat_dt
}

# 1. Process NAS
nas_sig <- nas_dt[padj < PADJ_THRESH & abs(logFC) > LFC_THRESH]
nas_all_degs <- unique(nas_sig$gene)

nas_freq <- nas_sig[, .N, by = gene]
nas_unique_degs <- nas_freq[N == 1, gene]
# Map unique DEG Ensembl IDs to symbols for later filtering
nas_unique_sym <- unique(nas_dt[gene %in% nas_unique_degs, symbol])

nas_mat_dt <- dcast(nas_dt[gene %in% nas_all_degs], gene ~ as.character(nas_level), value.var = "logFC")
nas_mat_dt <- use_symbols(nas_mat_dt, nas_dt)
nas_mat <- as.matrix(nas_mat_dt[, -"symbol", with = FALSE])
rownames(nas_mat) <- nas_mat_dt$symbol
nas_mat[is.na(nas_mat)] <- 0
nas_mat <- nas_mat[, setdiff(colnames(nas_mat), "0"), drop = FALSE]
nas_cols <- as.character(sort(as.integer(colnames(nas_mat))))
nas_mat <- nas_mat[, nas_cols, drop = FALSE]
colnames(nas_mat) <- paste0("NAS", colnames(nas_mat))

# 2. Process Fibrosis
fib_sig <- fib_dt[padj < PADJ_THRESH & abs(logFC) > LFC_THRESH]
fib_all_degs <- unique(fib_sig$gene)

fib_freq <- fib_sig[, .N, by = gene]
fib_unique_degs <- fib_freq[N == 1, gene]
fib_unique_sym <- unique(fib_dt[gene %in% fib_unique_degs, symbol])

fib_mat_dt <- dcast(fib_dt[gene %in% fib_all_degs], gene ~ as.character(fib_stage), value.var = "logFC")
fib_mat_dt <- use_symbols(fib_mat_dt, fib_dt)
fib_mat <- as.matrix(fib_mat_dt[, -"symbol", with = FALSE])
rownames(fib_mat) <- fib_mat_dt$symbol
fib_mat[is.na(fib_mat)] <- 0
fib_mat <- fib_mat[, setdiff(colnames(fib_mat), "0"), drop = FALSE]
fib_cols <- as.character(sort(as.integer(colnames(fib_mat))))
fib_mat <- fib_mat[, fib_cols, drop = FALSE]
colnames(fib_mat) <- paste0("F", colnames(fib_mat))

# Magenta-white-cyan diverging palette (matches publication_theme.R)
col_fun <- colorRamp2(
  c(-2.5, 0, 2.5),
  c(masld_colors$down, "white", masld_colors$up)
)

metric_colors <- c("NAS" = "#AD1457", "Fibrosis" = "#1565C0")

# Stage annotation: NAS uses magenta gradient, fibrosis uses blue gradient
# Use full labels (NAS0, F0) to avoid duplicate keys
stage_colors <- c(
  "NAS0" = "#FCE4EC", "NAS1" = "#F8BBD0", "NAS2" = "#F48FB1", "NAS3" = "#EC407A",
  "NAS4" = "#D81B60", "NAS5" = "#C2185B", "NAS6" = "#AD1457", "NAS7" = "#880E4F",
  "F0" = "#E3F2FD", "F1" = "#90CAF9", "F2" = "#42A5F5", "F3" = "#1565C0", "F4" = "#0D47A1"
)

# Known MASLD markers (20 genes from literature)
target_markers <- c(
  # Genetic risk / steatosis
  "PNPLA3", "TM6SF2", "HSD17B13", "MBOAT7", "GCKR",

  # Inflammation / NASH progression
  "CXCL10", "CCL2", "TREM2", "SPP1", "IL1B",
  # Fibrosis / ECM remodeling
  "COL1A1", "TIMP1", "TGFB1", "ACTA2", "LUM",
  # Lipid metabolism / stress
  "FASN", "SCD", "PPARA", "CYP7A1",
  # Clinical / staging biomarker
  "AKR1B10"
)

# Build set of significant target markers per context
nas_sig_markers <- unique(nas_dt[padj < PADJ_THRESH & abs(logFC) > LFC_THRESH & symbol %in% target_markers, symbol])
fib_sig_markers <- unique(fib_dt[padj < PADJ_THRESH & abs(logFC) > LFC_THRESH & symbol %in% target_markers, symbol])
all_sig_markers <- union(nas_sig_markers, fib_sig_markers)

cat("Significant target markers (NAS):", paste(sort(nas_sig_markers), collapse=", "), "\n")
cat("Significant target markers (Fib):", paste(sort(fib_sig_markers), collapse=", "), "\n")

get_anno_mark <- function(mat, sig_markers, n_labels = 20) {
  # Only label known genes that are actually significant DEGs
  found <- intersect(sig_markers, rownames(mat))
  if (length(found) < n_labels && nrow(mat) >= n_labels) {
    # Pad with top variable genes (skip any remaining Ensembl IDs)
    rv <- apply(mat, 1, var)
    top_v <- names(sort(rv, decreasing = TRUE))
    top_v <- top_v[!grepl("^ENSG", top_v)]
    top_c <- setdiff(top_v, found)
    found <- c(found, head(top_c, n_labels - length(found)))
  } else if (length(found) > n_labels) {
    found <- head(found, n_labels)
  }

  if (length(found) == 0) return(NULL)

  idx <- match(found, rownames(mat))
  ha <- rowAnnotation(
    genes = anno_mark(at = idx, labels = found,
                      labels_gp = gpar(fontsize = 8, fontface = "italic"),
                      padding = unit(1, "mm"))
  )
  return(ha)
}

# Classify each gene by its peak activation stage
get_peak_stage <- function(mat) {
  peak_col <- apply(abs(mat), 1, which.max)
  peak_stage <- colnames(mat)[peak_col]
  # Direction at peak
  peak_dir <- ifelse(mat[cbind(seq_len(nrow(mat)), peak_col)] > 0, "Up", "Down")
  # Ordered factor: columns in order, Up before Down within each
  all_stages <- colnames(mat)
  stage_dir <- paste0(peak_stage, " (", peak_dir, ")")
  # Create ordered levels
  lvls <- as.vector(outer(all_stages, c("Up", "Down"), function(s, d) paste0(s, " (", d, ")")))
  # Keep only levels that exist
  lvls <- lvls[lvls %in% unique(stage_dir)]
  factor(stage_dir, levels = lvls)
}

draw_heatmap <- function(mat, out_file, title, sig_markers, show_row_names = FALSE) {
  if (nrow(mat) == 0) return()

  mat_cap <- pmax(pmin(mat, 2.5), -2.5)

  right_ha <- get_anno_mark(mat, sig_markers)
  row_split <- get_peak_stage(mat)

  pdf(out_file, width = ncol(mat)*0.5 + 4, height = 7)
  ht <- Heatmap(mat_cap,
                name = "logFC\n(capped)",
                col = col_fun,
                cluster_columns = FALSE,
                show_row_names = show_row_names,
                show_row_dend = FALSE,
                right_annotation = right_ha,
                row_split = row_split,
                cluster_row_slices = FALSE,
                row_title_gp = gpar(fontsize = 7),
                row_title_rot = 0,
                row_gap = unit(0.5, "mm"),
                border = TRUE,
                column_title = title,
                use_raster = TRUE,
                raster_quality = 4)
  draw(ht)
  dev.off()
  cat("Saved:", out_file, "\n")
}

draw_heatmap(nas_mat, file.path(OUT_DIR, "nas_all_degs_heatmap.pdf"), "NAS (All DEGs)", nas_sig_markers)
draw_heatmap(nas_mat[rownames(nas_mat) %in% nas_unique_sym, , drop=FALSE],
             file.path(OUT_DIR, "nas_unique_degs_heatmap.pdf"), "NAS (Unique DEGs)", nas_sig_markers)

draw_heatmap(fib_mat, file.path(OUT_DIR, "fib_all_degs_heatmap.pdf"), "Fibrosis (All DEGs)", fib_sig_markers)
draw_heatmap(fib_mat[rownames(fib_mat) %in% fib_unique_sym, , drop=FALSE],
             file.path(OUT_DIR, "fib_unique_degs_heatmap.pdf"), "Fibrosis (Unique DEGs)", fib_sig_markers)

# Combined
combined_genes <- union(nas_all_degs, fib_all_degs)

nas_mat_full <- dcast(nas_dt[gene %in% combined_genes], gene ~ as.character(nas_level), value.var = "logFC")
if ("0" %in% names(nas_mat_full)) nas_mat_full[, `0` := NULL]
num_cols_nas <- setdiff(names(nas_mat_full), "gene")
setnames(nas_mat_full, num_cols_nas, paste0("NAS", num_cols_nas))

fib_mat_full <- dcast(fib_dt[gene %in% combined_genes], gene ~ as.character(fib_stage), value.var = "logFC")
if ("0" %in% names(fib_mat_full)) fib_mat_full[, `0` := NULL]
num_cols_fib <- setdiff(names(fib_mat_full), "gene")
setnames(fib_mat_full, num_cols_fib, paste0("F", num_cols_fib))

comb_dt <- merge(nas_mat_full, fib_mat_full, by = "gene", all = TRUE)
# Convert to symbols
comb_sym <- unique(rbind(nas_dt[, .(gene, symbol)], fib_dt[, .(gene, symbol)]))
comb_dt <- merge(comb_dt, comb_sym, by = "gene", all.x = TRUE)
comb_dt[is.na(symbol), symbol := gene]
# Deduplicate symbols by variance
num_comb_cols <- setdiff(names(comb_dt), c("gene", "symbol"))
comb_dt[, .row_var := apply(.SD, 1, var, na.rm = TRUE), .SDcols = num_comb_cols]
comb_dt <- comb_dt[comb_dt[, .I[which.max(.row_var)], by = symbol]$V1]
comb_dt[, c("gene", ".row_var") := NULL]

comb_mat <- as.matrix(comb_dt[, -"symbol", with = FALSE])
rownames(comb_mat) <- comb_dt$symbol
comb_mat[is.na(comb_mat)] <- 0

nas_full_cols <- grep("^NAS", colnames(comb_mat), value = TRUE)
nas_full_cols <- nas_full_cols[order(as.integer(gsub("NAS", "", nas_full_cols)))]
fib_full_cols <- grep("^F", colnames(comb_mat), value = TRUE)
fib_full_cols <- fib_full_cols[order(as.integer(gsub("F", "", fib_full_cols)))]

comb_mat_ordered <- cbind(
  comb_mat[, nas_full_cols, drop = FALSE],
  comb_mat[, fib_full_cols, drop = FALSE]
)

metric_vec <- c(rep("NAS", length(nas_full_cols)), rep("Fibrosis", length(fib_full_cols)))
stage_vec <- c(nas_full_cols, fib_full_cols)  # Use full labels: NAS0, F0, etc.

ha_top <- HeatmapAnnotation(
  Group = metric_vec,
  Stage = stage_vec,
  col = list(
    Group = metric_colors,
    Stage = stage_colors
  ),
  border = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 8),
  simple_anno_size = unit(3, "mm")
)

split_vec <- factor(metric_vec, levels=c("NAS", "Fibrosis"))

draw_combined <- function(mat, out_file, title) {
  if (nrow(mat) == 0) return()
  mat_cap <- pmax(pmin(mat, 2.5), -2.5)

  right_ha <- get_anno_mark(mat, all_sig_markers)
  row_split <- get_peak_stage(mat)

  pdf(out_file, width = ncol(mat)*0.4 + 5, height = 8)
  ht <- Heatmap(mat_cap,
                name = "logFC\n(capped)",
                col = col_fun,
                top_annotation = ha_top,
                cluster_columns = FALSE,
                column_split = split_vec,
                column_title = NULL,
                show_row_names = FALSE,
                show_column_names = FALSE,
                show_row_dend = FALSE,
                right_annotation = right_ha,
                row_split = row_split,
                cluster_row_slices = FALSE,
                row_title_gp = gpar(fontsize = 7),
                row_title_rot = 0,
                row_gap = unit(0.5, "mm"),
                border = TRUE,
                use_raster = TRUE,
                raster_quality = 4)
  draw(ht)
  dev.off()
  cat("Saved:", out_file, "\n")
}

draw_combined(comb_mat_ordered, file.path(OUT_DIR, "combined_nas_fib_all_degs_heatmap.pdf"), "NAS & Fibrosis (All DEGs)")

comb_unique_sym <- union(nas_unique_sym, fib_unique_sym)
draw_combined(comb_mat_ordered[rownames(comb_mat_ordered) %in% comb_unique_sym, , drop=FALSE],
              file.path(OUT_DIR, "combined_nas_fib_unique_degs_heatmap.pdf"), "NAS & Fibrosis (Unique DEGs)")

# --- Unified heatmap: NAS + Fibrosis columns together, no column_split ---
# Genes classified by BOTH NAS and fibrosis DEG status in row annotation

# Row annotation: which metric(s) each gene is a DEG in
nas_sym_set <- unique(nas_dt[gene %in% nas_all_degs, symbol])
fib_sym_set <- unique(fib_dt[gene %in% fib_all_degs, symbol])

deg_source <- ifelse(
  rownames(comb_mat_ordered) %in% nas_sym_set & rownames(comb_mat_ordered) %in% fib_sym_set,
  "Both",
  ifelse(rownames(comb_mat_ordered) %in% nas_sym_set, "NAS only", "Fibrosis only")
)

draw_unified <- function(mat, deg_class, out_file, title) {
  if (nrow(mat) == 0) return()
  mat_cap <- pmax(pmin(mat, 2.5), -2.5)

  right_ha <- get_anno_mark(mat, all_sig_markers)

  # Row split by peak activation stage (clear blocks)
  row_split <- get_peak_stage(mat)

  # Column annotation: metric + stage (no column_split)
  col_metric <- ifelse(grepl("^NAS", colnames(mat)), "NAS", "Fibrosis")
  col_stage <- colnames(mat)  # Full labels: NAS0, F0, etc.

  ha_top_unified <- HeatmapAnnotation(
    Metric = col_metric,
    Stage = col_stage,
    col = list(Metric = metric_colors, Stage = stage_colors),
    border = TRUE,
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 8),
    simple_anno_size = unit(3, "mm")
  )

  # Left annotation: DEG source
  deg_colors <- c("Both" = "#7B1FA2", "NAS only" = "#AD1457", "Fibrosis only" = "#1565C0")
  ha_left <- rowAnnotation(
    Source = deg_class,
    col = list(Source = deg_colors),
    border = TRUE,
    show_annotation_name = TRUE,
    annotation_name_side = "top",
    simple_anno_size = unit(3, "mm")
  )

  pdf(out_file, width = ncol(mat)*0.4 + 5.5, height = 8)
  ht <- Heatmap(mat_cap,
                name = "logFC\n(capped)",
                col = col_fun,
                top_annotation = ha_top_unified,
                left_annotation = ha_left,
                cluster_columns = FALSE,
                show_row_names = FALSE,
                show_column_names = TRUE,
                column_names_rot = 45,
                column_names_gp = gpar(fontsize = 9),
                show_row_dend = FALSE,
                right_annotation = right_ha,
                row_split = row_split,
                cluster_row_slices = FALSE,
                row_title_gp = gpar(fontsize = 7),
                row_title_rot = 0,
                row_gap = unit(0.5, "mm"),
                border = TRUE,
                column_title = title,
                use_raster = TRUE,
                raster_quality = 4)
  draw(ht)
  dev.off()
  cat("Saved:", out_file, "\n")
}

draw_unified(comb_mat_ordered, deg_source,
             file.path(OUT_DIR, "unified_nas_fib_all_degs_heatmap.pdf"),
             "NAS & Fibrosis DEGs (Unified)")

# Unified unique DEGs
unified_unique_idx <- rownames(comb_mat_ordered) %in% comb_unique_sym
if (sum(unified_unique_idx) > 0) {
  draw_unified(comb_mat_ordered[unified_unique_idx, , drop=FALSE],
               deg_source[unified_unique_idx],
               file.path(OUT_DIR, "unified_nas_fib_unique_degs_heatmap.pdf"),
               "NAS & Fibrosis Unique DEGs (Unified)")
}

cat("Heatmap generation complete.\n")
