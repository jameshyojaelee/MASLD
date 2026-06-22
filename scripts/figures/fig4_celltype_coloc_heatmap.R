#!/usr/bin/env Rscript
# =============================================================================
# Fig 4 Panel: Cell-type COLOC Resolution Heatmap
# =============================================================================
# Shows cell-type-resolved MR evidence from sc-eQTL TWAS across 4 liver cell
# types. Rows = genes significant (p<0.05) in >=1 cell type; columns = cell
# types. Fill = signed -log10(p) from IVW/Wald ratio MR, capped +-10.
# Annotations: bulk dream logFC (left), COLOC PP.H4 + drug target (right),
# heritability fold enrichment (top).
# =============================================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

cat("=== Cell-type COLOC Resolution Heatmap ===\n")
cat("Timestamp:", format(Sys.time()), "\n\n")

# ---- Output directory -------------------------------------------------------
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
dir.create(file.path(FIG3_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)
out_dir <- file.path(FIG3_DIR, "panels")

# =============================================================================
# 1. Load sc-eQTL TWAS MR results
# =============================================================================
cat("Loading sc-eQTL TWAS results...\n")
mr_raw <- fread(file.path(BASE,
  "RNA-seq/results/causal_inference/sceqtl_twas/sceqtl_twas_all_results.csv"))
cat("  Raw rows:", nrow(mr_raw), "\n")

# Keep IVW or Wald ratio (primary MR estimates)
mr <- mr_raw[method %in% c("Inverse variance weighted", "Wald ratio")]
cat("  IVW/Wald rows:", nrow(mr), "\n")

# Best p-value per gene x cell_type (across GWAS and condition)
# Take the entry with the smallest p-value
mr_best <- mr[, {
  idx <- which.min(pval)
  list(best_p = pval[idx], best_b = b[idx], best_se = se[idx],
       best_gwas = gwas[idx], best_fdr = fdr[idx])
}, by = .(exposure, cell_type)]
cat("  Gene x cell_type pairs:", nrow(mr_best), "\n")

# =============================================================================
# 2. Select genes: p < 0.05 in >= 1 cell type
# =============================================================================
sig_genes <- mr_best[best_p < 0.05, .(n_sig_ct = .N), by = exposure]
cat("  Genes sig (p<0.05) in >=1 cell type:", nrow(sig_genes), "\n")

# With ~2971 genes, we need to filter for a readable heatmap (50-200 target)
# Strategy: use FDR < 0.1 in >= 1 cell type as primary filter
fdr_genes <- mr_best[best_fdr < 0.1, .(n_sig_ct = .N), by = exposure]
cat("  Genes sig (FDR<0.1) in >=1 cell type:", nrow(fdr_genes), "\n")

# If FDR<0.1 gives 50-200, use that; otherwise adjust
if (nrow(fdr_genes) >= 50 && nrow(fdr_genes) <= 200) {
  selected_genes <- fdr_genes$exposure
  cat("  Using FDR<0.1 filter:", length(selected_genes), "genes\n")
} else if (nrow(fdr_genes) > 200) {
  # Tighten: require FDR < 0.1 in >= 2 cell types, or use stricter FDR
  fdr2_genes <- mr_best[best_fdr < 0.1, .(n_sig_ct = .N), by = exposure][n_sig_ct >= 2]
  if (nrow(fdr2_genes) >= 50) {
    selected_genes <- fdr2_genes$exposure
    cat("  Using FDR<0.1 in >=2 cell types:", length(selected_genes), "genes\n")
  } else {
    # Use FDR < 0.05
    fdr05_genes <- mr_best[best_fdr < 0.05, .(n_sig_ct = .N), by = exposure]
    if (nrow(fdr05_genes) >= 50 && nrow(fdr05_genes) <= 200) {
      selected_genes <- fdr05_genes$exposure
      cat("  Using FDR<0.05 filter:", length(selected_genes), "genes\n")
    } else {
      # Take top 150 by minimum p-value across cell types
      top_genes <- mr_best[, .(min_p = min(best_p)), by = exposure][order(min_p)][1:150]
      selected_genes <- top_genes$exposure
      cat("  Using top 150 by min p-value\n")
    }
  }
} else {
  # Too few with FDR<0.1 -- relax to p < 0.01
  p001_genes <- mr_best[best_p < 0.01, .(n_sig_ct = .N), by = exposure]
  if (nrow(p001_genes) >= 50 && nrow(p001_genes) <= 200) {
    selected_genes <- p001_genes$exposure
    cat("  Using p<0.01 filter:", length(selected_genes), "genes\n")
  } else if (nrow(p001_genes) > 200) {
    top_genes <- mr_best[, .(min_p = min(best_p)), by = exposure][order(min_p)][1:150]
    selected_genes <- top_genes$exposure
    cat("  Using top 150 by min p-value\n")
  } else {
    # Use p < 0.05 but cap at 200
    p05_genes <- mr_best[best_p < 0.05, .(min_p = min(best_p)), by = exposure][order(min_p)]
    selected_genes <- head(p05_genes$exposure, 150)
    cat("  Using top 150 from p<0.05 genes\n")
  }
}

cat("  Final gene count:", length(selected_genes), "\n")

# =============================================================================
# 3. Build heatmap matrix: signed -log10(p), capped at +-10
# =============================================================================
mat_data <- mr_best[exposure %in% selected_genes]

# Compute signed -log10(p)
mat_data[, signed_logp := -log10(pmax(best_p, 1e-20)) * sign(best_b)]

# Cap at +-10
mat_data[, signed_logp := pmin(pmax(signed_logp, -10), 10)]

# Pivot to wide matrix
mat_wide <- dcast(mat_data, exposure ~ cell_type, value.var = "signed_logp", fill = 0)
gene_names <- mat_wide$exposure
mat <- as.matrix(mat_wide[, -1, with = FALSE])
rownames(mat) <- gene_names

# Order columns by total number of significant genes (descending)
col_order <- mr_best[best_p < 0.05, .N, by = cell_type][order(-N)]$cell_type
col_order <- col_order[col_order %in% colnames(mat)]
# Add any missing columns
remaining <- setdiff(colnames(mat), col_order)
col_order <- c(col_order, remaining)
mat <- mat[, col_order, drop = FALSE]

# Pretty column names
col_labels <- gsub("_", " ", colnames(mat))
col_labels <- gsub("\\b(\\w)", "\\U\\1", col_labels, perl = TRUE)

cat("  Matrix dimensions:", nrow(mat), "x", ncol(mat), "\n")

# =============================================================================
# 4. Load annotation data
# =============================================================================

# ---- 4a. Bulk dream logFC (ENSEMBL ID mapping via COLOC gene-level) ---------
cat("Loading COLOC gene-level for ID mapping...\n")
coloc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc <- coloc[gene != "" & !is.na(gene)]
coloc[, ensembl_clean := sub("\\.\\d+$", "", ensembl)]

cat("Loading dream results...\n")
dream <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"))
dream[, ensembl_clean := sub("\\.\\d+$", "", gene)]

# Map dream ENSEMBL to gene symbols via COLOC table
ens2sym <- coloc[, .(gene, ensembl_clean)]
ens2sym <- ens2sym[!duplicated(ensembl_clean)]

dream_mapped <- merge(dream, ens2sym, by = "ensembl_clean", suffixes = c("_dream", "_sym"))
dream_mapped <- dream_mapped[gene_sym %in% selected_genes]

# If gene_sym column not created properly, try alternative
if (!"gene_sym" %in% names(dream_mapped)) {
  # The merge creates gene_dream (ENSEMBL) and gene (symbol)
  setnames(dream_mapped, "gene", "gene_sym", skip_absent = TRUE)
}

# Get logFC per gene (some genes may have multiple ENSEMBL IDs - take the one with smallest padj)
dream_lfc <- dream_mapped[, {
  idx <- which.min(padj)
  list(logFC = logFC[idx], padj = padj[idx])
}, by = gene_sym]
setnames(dream_lfc, "gene_sym", "gene")

cat("  Dream logFC mapped for", nrow(dream_lfc), "of", length(selected_genes), "genes\n")

# ---- 4b. COLOC PP.H4 from gene_level_coloc.csv -----------------------------
coloc_pp4 <- coloc[gene %in% selected_genes, .(gene, coloc_best_pp4, coloc_n_gwas_h4_05)]
cat("  COLOC PP.H4 available for", nrow(coloc_pp4), "genes\n")

# ---- 4c. Drug target status (DGIdb + OpenTargets + LINCS + curated) ----------
atlas_drug <- fread(file.path(BASE,
  "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
  select = c("human_symbol", "dgidb_druggable", "opentargets_drug", "progression_drug_reversal_n"))
# Broad druggability: DGIdb OR OpenTargets OR has LINCS reversal OR curated MASLD targets
curated_targets <- c("THRB", "GLP1R", "FGFR1", "KLB", "GCGR", "DGAT2", "HSD17B13",
                     "NR1H4", "PPARA", "PPARD", "PPARG", "SCD", "FASN", "ACC1",
                     "FXR", "PNPLA3", "MBOAT7", "MARC1", "GPAM")
drug_genes <- unique(c(
  atlas_drug[dgidb_druggable == TRUE]$human_symbol,
  atlas_drug[!is.na(opentargets_drug) & opentargets_drug != ""]$human_symbol,
  atlas_drug[!is.na(progression_drug_reversal_n) & progression_drug_reversal_n > 0]$human_symbol,
  curated_targets
))
drug_in_heatmap <- intersect(drug_genes, selected_genes)
cat("  Druggable genes in heatmap:", length(drug_in_heatmap),
    ifelse(length(drug_in_heatmap) > 0,
           paste0(" (", paste(head(drug_in_heatmap, 15), collapse = ", "),
                  ifelse(length(drug_in_heatmap) > 15, ", ...", ""), ")"),
           "(none)"), "\n")

# ---- 4d. INTACT scores (cell-type) ------------------------------------------
cat("Loading INTACT scores...\n")
intact <- fread(file.path(BASE, "RNA-seq/results/gwas_rna_integration/intact_scores.csv"))
intact_ann <- intact[gene %in% selected_genes, .(gene, multi_intact_score, intact_n_celltypes_05)]
cat("  INTACT scores for", nrow(intact_ann), "genes\n")

# ---- 4e. Cell-type heritability enrichment -----------------------------------
cat("Loading cell-type heritability...\n")
herit <- fread(file.path(BASE,
  "RNA-seq/results/gwas_rna_integration/celltype_heritability_results.csv"))

# Map heritability cell-type names to sc-eQTL cell-type names
ct_name_map <- c(
  "Hepatocytes"             = "hepatocyte",
  "Cholangiocytes"          = "cholangiocyte",
  "Endothelial_cells"       = "endothelial_cell",
  "Fibroblasts"             = "stellate_cell"  # fibroblasts ~ stellate cells in liver
)
herit[, ct_mapped := ct_name_map[cell_type]]
herit_mapped <- herit[!is.na(ct_mapped) & ct_mapped %in% colnames(mat)]
cat("  Heritability enrichment mapped for", nrow(herit_mapped), "cell types\n")

# =============================================================================
# 5. Build annotations
# =============================================================================

# ---- Left annotation: bulk dream logFC --------------------------------------
dream_vec <- rep(NA_real_, nrow(mat))
names(dream_vec) <- rownames(mat)
m <- match(dream_lfc$gene, rownames(mat))
dream_vec[m[!is.na(m)]] <- dream_lfc$logFC[!is.na(m)]

# Dream significance
dream_sig_vec <- rep(NA_real_, nrow(mat))
names(dream_sig_vec) <- rownames(mat)
m2 <- match(dream_lfc$gene, rownames(mat))
dream_sig_vec[m2[!is.na(m2)]] <- dream_lfc$padj[!is.na(m2)]

lfc_max <- max(abs(dream_vec), na.rm = TRUE)
lfc_cap <- min(lfc_max, 3)

ha_left <- rowAnnotation(
  `Integrated\nlogFC` = anno_barplot(
    dream_vec,
    gp = gpar(fill = ifelse(dream_vec > 0, masld_colors$up, masld_colors$down)),
    width = unit(1.5, "cm"),
    axis_param = list(gp = gpar(fontsize = 5)),
    ylim = c(-lfc_cap, lfc_cap)
  ),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 6),
  annotation_name_rot = 90
)

# ---- Right annotations: COLOC PP.H4 + drug target + INTACT ------------------
pp4_vec <- rep(0, nrow(mat))
names(pp4_vec) <- rownames(mat)
m3 <- match(coloc_pp4$gene, rownames(mat))
pp4_vec[m3[!is.na(m3)]] <- coloc_pp4$coloc_best_pp4[!is.na(m3)]

# Drug target indicator
drug_vec <- ifelse(rownames(mat) %in% drug_genes, "Yes", "No")

# INTACT score
intact_vec <- rep(NA_real_, nrow(mat))
names(intact_vec) <- rownames(mat)
m4 <- match(intact_ann$gene, rownames(mat))
intact_vec[m4[!is.na(m4)]] <- intact_ann$multi_intact_score[!is.na(m4)]

ha_right <- rowAnnotation(
  `COLOC PP.H4` = anno_barplot(
    pp4_vec,
    gp = gpar(fill = "#00695C"),
    width = unit(1.2, "cm"),
    axis_param = list(gp = gpar(fontsize = 5)),
    ylim = c(0, 1)
  ),
  `INTACT` = anno_barplot(
    intact_vec,
    gp = gpar(fill = "#7B1FA2"),
    width = unit(1.0, "cm"),
    axis_param = list(gp = gpar(fontsize = 5)),
    ylim = c(0, 1)
  ),
  `Drug\nTarget` = anno_simple(
    drug_vec,
    col = c("Yes" = "#C9265E", "No" = "white"),
    width = unit(0.3, "cm"),
    border = TRUE
  ),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 6)
)

# ---- Top annotation: cell-type heritability fold enrichment ------------------
herit_vals <- rep(NA_real_, ncol(mat))
names(herit_vals) <- colnames(mat)
for (i in seq_len(nrow(herit_mapped))) {
  ct <- herit_mapped$ct_mapped[i]
  if (ct %in% names(herit_vals)) {
    herit_vals[ct] <- herit_mapped$fold_enrichment_coloc[i]
  }
}

# Significance stars for heritability
herit_stars <- rep("", ncol(mat))
names(herit_stars) <- colnames(mat)
for (i in seq_len(nrow(herit_mapped))) {
  ct <- herit_mapped$ct_mapped[i]
  if (ct %in% names(herit_stars)) {
    fdr_val <- herit_mapped$fdr_coloc[i]
    if (!is.na(fdr_val)) {
      if (fdr_val < 0.001) herit_stars[ct] <- "***"
      else if (fdr_val < 0.01) herit_stars[ct] <- "**"
      else if (fdr_val < 0.05) herit_stars[ct] <- "*"
    }
  }
}

# Cell-type colors matching ct_palette (map lowercase to display names)
ct_color_map <- c(
  hepatocyte       = "#0D47A1",
  cholangiocyte    = "#1565C0",
  endothelial_cell = "#2E7D32",
  stellate_cell    = "#F57F17"
)

ha_top <- HeatmapAnnotation(
  `GWAS enrichment\n(fold vs background)` = anno_barplot(
    herit_vals,
    gp = gpar(fill = ct_color_map[colnames(mat)]),
    height = unit(1.5, "cm"),
    axis_param = list(gp = gpar(fontsize = 5)),
    baseline = 1  # reference line at fold = 1
  ),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 5),
  annotation_name_side = "left"
)

# =============================================================================
# 6. Color scale: diverging blue-white-red for signed -log10(p)
# =============================================================================
col_fun <- colorRamp2(
  breaks = c(-10, -5, -2, 0, 2, 5, 10),
  colors = c("#0D47A1", "#1565C0", "#90CAF9", "white", "#F4A674", "#C9265E", "#880E4F")
)

# =============================================================================
# 7. Draw heatmap
# =============================================================================
cat("Drawing heatmap...\n")

# Determine whether to show row names based on gene count
show_row_names <- length(selected_genes) <= 100
row_fontsize <- if (length(selected_genes) <= 60) 5 else if (length(selected_genes) <= 100) 4 else 3

ht <- Heatmap(
  mat,
  name = "Signed\n-log10(p)",
  col = col_fun,

  # Clustering
  clustering_method_rows    = "ward.D2",
  clustering_distance_rows  = "euclidean",
  clustering_method_columns = "ward.D2",
  clustering_distance_columns = "euclidean",

  # Row settings
  show_row_names = show_row_names,
  row_names_gp   = gpar(fontsize = row_fontsize, fontface = ifelse(
    rownames(mat) %in% drug_genes, "bold", "plain"
  )),
  row_names_side = "left",
  row_dend_width = unit(1.5, "cm"),


  # Column settings
  column_labels  = col_labels,
  column_names_gp = gpar(fontsize = 5),
  column_names_rot = 45,
  show_column_dend = TRUE,
  column_dend_height = unit(0.8, "cm"),

  # Cell settings
  rect_gp = gpar(col = "grey90", lwd = 0.25),

  # Annotations
  left_annotation  = ha_left,
  right_annotation = ha_right,
  top_annotation   = ha_top,

  # Legend
  heatmap_legend_param = list(
    title = "scTWAS\nsign(beta) x\n-log10(p)",
    title_gp = gpar(fontsize = 5, fontface = "bold"),
    labels_gp = gpar(fontsize = 5),
    legend_height = unit(3, "cm"),
    at = c(-10, -5, 0, 5, 10),
    labels = c("-10\n(down in\ndisease)", "-5", "0", "5", "10\n(up in\ndisease)")
  ),

  # Row split by k-means for visual grouping if many genes
  row_km = if (length(selected_genes) > 80) 4 else if (length(selected_genes) > 40) 3 else 2,
  row_km_repeats = 50,
  row_title_gp = gpar(fontsize = 7, fontface = "bold"),
  row_gap = unit(1, "mm"),

  # Width/height
  width = unit(4, "cm"),
  border = TRUE
)

# =============================================================================
# 8. Save figure
# =============================================================================
out_file <- file.path(out_dir, "panel_celltype_coloc_heatmap.pdf")
cat("Saving to:", out_file, "\n")

# Calculate height based on gene count
fig_height <- max(6, min(14, length(selected_genes) * 0.08 + 3))
fig_width  <- 8  # inches

# RETIRED 2026-06-17: dead/illegible panel, cut in the Fig 2 review and archived.
# This script misroutes to FIG3_DIR (=fig2_genetics). Archived copy:
# figures/main/fig2_genetics/panels/_archive/panel_celltype_coloc_heatmap.pdf
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
if (FALSE) {
pdf_device(out_file, width = fig_width, height = fig_height)
draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     padding = unit(c(5, 5, 5, 5), "mm"),
     merge_legend = TRUE)
dev.off()
cat("  Saved:", out_file, "\n")
}

# =============================================================================
# 9. Summary statistics
# =============================================================================
cat("\n=== Summary ===\n")
cat("Genes in heatmap:", nrow(mat), "\n")
cat("Cell types:", ncol(mat), "\n")
cat("Signed -log10(p) range:", round(range(mat), 2), "\n")

# Gene with strongest signal per cell type
for (ct in colnames(mat)) {
  idx <- which.max(abs(mat[, ct]))
  cat(sprintf("  %s: top gene = %s (%.2f)\n", ct, rownames(mat)[idx], mat[idx, ct]))
}

# Count genes with COLOC support
cat("\nGenes with COLOC PP.H4 > 0.5:", sum(pp4_vec > 0.5), "\n")
cat("Genes with COLOC PP.H4 > 0.8:", sum(pp4_vec > 0.8), "\n")
cat("Genes with INTACT > 0.5:", sum(intact_vec > 0.5, na.rm = TRUE), "\n")

# Drug targets present
cat("Drug targets in heatmap:", paste(drug_in_heatmap, collapse = ", "),
    ifelse(length(drug_in_heatmap) == 0, "(none)", ""), "\n")

# Dream logFC coverage
cat("Genes with dream logFC:", sum(!is.na(dream_vec)), "/", length(dream_vec), "\n")

cat("\nDone.\n")
