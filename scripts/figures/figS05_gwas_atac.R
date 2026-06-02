#!/usr/bin/env Rscript
# ==========================================================================
# figS05_gwas_atac.R — Supplementary Figure 5: GWAS-ATAC Regulatory Variants
#
# Generates panels: DA volcanos, peak annotations, chromVAR heatmap,
# mouse bulk ATAC panels, SCENIC+ correlation, cross-species conservation.
#
# Output: figures/supplementary/figS05_epigenomic_spatial/panel_*.pdf
# ==========================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

# --- Paths ----------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts", "figures", "publication_theme.R"))
source(file.path(BASE, "scripts", "figures", "load_figure_data.R"))

ATAC_HM  <- file.path(BASE, "Analysis", "ATAC", "Human_Multiome")
ATAC_MB  <- file.path(BASE, "Analysis", "ATAC", "Mouse_Bulk")
ATAC_INT <- file.path(BASE, "Analysis", "ATAC", "Integration")
DREAM    <- file.path(BASE, "RNA-seq", "Human", "Patient_Cohorts",
                      "analysis", "integration", "results", "integration")
OUT      <- FIGS05_DIR
PANELS   <- file.path(FIGS05_DIR, "panels")
dir.create(PANELS, showWarnings = FALSE, recursive = TRUE)

save_panel <- function(p, name, w = fig_half_width, h = 3) {
  path <- file.path(PANELS, name)
  save_fig(p, path, width = w, height = h)
  message("  Saved ", path)
}

cat("==============================================================\n")
cat("Generating supplementary ATAC panels (R)\n")
cat(paste0("Output: ", OUT, "\n"))
cat("==============================================================\n")

# ===================================================================
# Panel 9: scATAC DA volcano plot (hepatocyte-only, corrected)
# ===================================================================
message("\nPanel 09: scATAC DA volcano (hepatocyte, corrected)...")
# Prefer corrected results; fall back to original
da_file_corrected <- file.path(ATAC_HM, "results", "l8_annotated_corrected",
                               "scatac_da_gene_annotated.csv")
da_file_original  <- file.path(ATAC_HM, "results", "l8_annotated",
                               "scatac_da_gene_annotated.csv")
da_file <- if (file.exists(da_file_corrected)) da_file_corrected else da_file_original
da_source <- if (file.exists(da_file_corrected)) "corrected" else "original"
message("  Using ", da_source, " DA results: ", da_file)

if (file.exists(da_file)) {
  da <- fread(da_file)

  # Subset to hepatocytes only
  if ("cell_type" %in% names(da)) {
    da <- da[cell_type == "Hepatocyte"]
  }

  da[, neg_log10_padj := -log10(pmin(padj, 1))]
  da[, sig := fifelse(padj < 0.05 & abs(logFC) > 0.25,
                       fifelse(logFC > 0, "Up", "Down"), "NS")]

  n_up   <- sum(da$sig == "Up", na.rm = TRUE)
  n_down <- sum(da$sig == "Down", na.rm = TRUE)

  # Winsorize logFC for display
  lfc_cap <- quantile(abs(da$logFC), 0.995, na.rm = TRUE)
  lfc_cap <- max(lfc_cap, 1.0)
  da[, logFC_plot := pmin(pmax(logFC, -lfc_cap), lfc_cap)]

  # Label top 5 genes per direction — deduplicate: keep best peak per gene
  da[, rank_score := abs(logFC) * neg_log10_padj]
  da_best <- da[sig != "NS" & !is.na(gene_symbol) & gene_symbol != ""]
  da_best <- da_best[!grepl("^ENSG", gene_symbol)]
  da_best <- da_best[da_best[, .I[which.max(rank_score)], by = gene_symbol]$V1]
  top_down <- da_best[sig == "Down"][order(-rank_score)][1:min(5, .N)]
  top_up   <- da_best[sig == "Up"][order(-rank_score)][1:min(5, .N)]
  top_genes <- rbind(top_down, top_up)
  top_genes[, logFC_plot := pmin(pmax(logFC, -lfc_cap), lfc_cap)]

  # Symmetric x-axis
  xlim_val <- lfc_cap + 0.3

  subtitle_text <- paste0(
    formatC(nrow(da), big.mark = ","), " peaks, ",
    n_up, " opened, ", n_down, " closed"
  )

  p9 <- ggplot(da, aes(x = logFC_plot, y = neg_log10_padj, color = sig)) +
    geom_point(size = 0.3, alpha = 0.4, shape = 16) +
    scale_color_manual(values = c(Down = masld_colors$down,
                                  NS   = masld_colors$ns,
                                  Up   = masld_colors$up),
                       name = NULL) +
    geom_text_repel(data = top_genes,
                    aes(label = gene_symbol),
                    size = 1.8, max.overlaps = 25,
                    segment.size = 0.2, min.segment.length = 0) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed",
               linewidth = 0.3, color = "gray60") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               linewidth = 0.3, color = "gray60") +
    coord_cartesian(xlim = c(-xlim_val, xlim_val)) +
    labs(x = expression(log[2]~"fold change (MASLD / Normal)"),
         y = expression(-log[10]~"(adjusted p-value)"),
         title = "Hepatocyte Differential Accessibility",
         subtitle = subtitle_text) +
    annotate("text", x = xlim_val * 0.7,
             y = max(da$neg_log10_padj, na.rm = TRUE) * 0.95,
             label = paste0(n_up, " opened"), color = masld_colors$up,
             size = 2.2, hjust = 0.5) +
    annotate("text", x = -xlim_val * 0.7,
             y = max(da$neg_log10_padj, na.rm = TRUE) * 0.95,
             label = paste0(n_down, " closed"), color = masld_colors$down,
             size = 2.2, hjust = 0.5) +
    theme_masld() +
    theme(legend.position = "none")

  save_panel(p9, "panel_09_scatac_da_volcano.pdf", w = fig_half_width, h = 3.5)
} else {
  message("  WARNING: ", da_file, " not found, skipping panel 09")
}

# ===================================================================
# Panel 10: Peak genomic annotation distribution (scATAC)
# ===================================================================
message("\nPanel 10: Peak genomic annotation distribution...")
annot_file <- file.path(ATAC_MB, "results", "annotated_peaks_primary.csv")

# Use all scATAC annotated peaks (all cell types)
da_annot <- if (file.exists(da_file)) fread(da_file) else NULL

if (!is.null(da_annot) && "distance_to_tss" %in% names(da_annot)) {
  # Classify peaks by distance to TSS
  da_annot[, annotation := fcase(
    abs(distance_to_tss) <= 1000,  "Promoter (<1kb)",
    abs(distance_to_tss) <= 3000,  "Promoter (1-3kb)",
    abs(distance_to_tss) <= 10000, "Proximal (3-10kb)",
    abs(distance_to_tss) <= 100000, "Distal (10-100kb)",
    default = "Intergenic (>100kb)"
  )]

  annot_counts <- da_annot[, .N, by = annotation]
  annot_counts[, pct := N / sum(N) * 100]
  annot_counts[, annotation := factor(annotation,
    levels = c("Promoter (<1kb)", "Promoter (1-3kb)", "Proximal (3-10kb)",
               "Distal (10-100kb)", "Intergenic (>100kb)"))]

  annot_colors <- c("Promoter (<1kb)" = "#0D47A1", "Promoter (1-3kb)" = "#1565C0",
                     "Proximal (3-10kb)" = "#42A5F5", "Distal (10-100kb)" = "#E91E63",
                     "Intergenic (>100kb)" = "#BDBDBD")

  p10 <- ggplot(annot_counts, aes(x = annotation, y = pct, fill = annotation)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = annot_colors) +
    geom_text(aes(label = paste0(round(pct, 1), "%")),
              vjust = -0.3, size = 2) +
    labs(x = NULL, y = "Percentage of DA peaks",
         title = "scATAC DA peak genomic distribution") +
    theme_masld() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1, size = 5))

  save_panel(p10, "panel_10_peak_annotation.pdf", w = fig_half_width, h = 3)
} else {
  message("  WARNING: No annotated scATAC peaks available, skipping panel 10")
}

# ===================================================================
# Panel 11: chromVAR TF motif heatmap (top TFs across cell types)
# ===================================================================
message("\nPanel 11: chromVAR TF heatmap...")
# Prefer label-transfer-corrected chromVAR v2; fall back to original
cv_file <- file.path(ATAC_HM, "results", "chromvar_v2", "chromvar_tf_activity.csv")
if (!file.exists(cv_file))
  cv_file <- file.path(ATAC_HM, "results", "chromvar", "chromvar_tf_activity.csv")

if (file.exists(cv_file)) {
  cv <- fread(cv_file)

  # Get top 30 TFs by significance in hepatocytes (handle both naming conventions)
  cv_hep <- cv[cell_type %in% c("Hepatocyte", "Hepatocytes")][order(padj)][1:min(30, .N)]
  hep_ct_name <- cv_hep$cell_type[1]  # "Hepatocyte" or "Hepatocytes"
  top_tfs <- cv_hep$tf_name

  # Subset to these TFs across all cell types
  cv_sub <- cv[tf_name %in% top_tfs]
  cv_sub[, logFC_capped := pmin(pmax(logFC_deviation, -1.5), 1.5)]

  # Create matrix (mean per TF × cell type if duplicates)
  mat <- dcast(cv_sub, tf_name ~ cell_type, value.var = "logFC_capped",
               fun.aggregate = mean, fill = 0)
  tf_names <- mat$tf_name
  mat_vals <- as.matrix(mat[, -1, with = FALSE])
  rownames(mat_vals) <- tf_names

  # Order TFs by hepatocyte deviation
  hep_order <- order(mat_vals[, hep_ct_name], decreasing = TRUE)
  mat_vals <- mat_vals[hep_order, , drop = FALSE]

  # Plot as tile
  plot_df <- as.data.table(reshape2::melt(mat_vals,
                           varnames = c("TF", "Cell_Type"),
                           value.name = "Deviation"))

  p11 <- ggplot(plot_df, aes(x = Cell_Type, y = TF, fill = Deviation)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = masld_colors$down, mid = "white",
                         high = masld_colors$up, midpoint = 0,
                         name = expression(Delta~"deviation"),
                         limits = c(-1.5, 1.5)) +
    labs(x = NULL, y = NULL,
         title = "chromVAR TF motif deviations (MASLD vs Normal)") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
          axis.text.y = element_text(size = 4.5))

  save_panel(p11, "panel_11_chromvar_heatmap.pdf",
             w = fig_half_width + 0.5, h = 5)
} else {
  message("  WARNING: ", cv_file, " not found, skipping panel 11")
}

# ===================================================================
# Panel 12: Mouse bulk ATAC peak annotation distribution
# ===================================================================
message("\nPanel 12: Mouse peak annotation distribution...")

if (file.exists(annot_file)) {
  mouse_peaks <- fread(annot_file)

  # Simplify ChIPseeker annotation categories
  mouse_peaks[, annot_simple := fcase(
    grepl("^Promoter", annotation),          "Promoter",
    grepl("Intron",    annotation),          "Intron",
    grepl("Exon",      annotation),          "Exon",
    grepl("UTR",       annotation),          "UTR",
    grepl("Downstream", annotation),         "Downstream",
    grepl("Distal Intergenic", annotation),  "Intergenic",
    default = "Other"
  )]

  annot_counts_m <- mouse_peaks[, .N, by = annot_simple]
  annot_counts_m[, pct := N / sum(N) * 100]
  annot_counts_m[, annot_simple := factor(annot_simple,
    levels = c("Promoter", "Exon", "UTR", "Intron", "Downstream",
               "Intergenic", "Other"))]

  annot_colors_m <- c(Promoter = "#0D47A1", Exon = "#1565C0",
                       UTR = "#42A5F5", Intron = "#7B1FA2",
                       Downstream = "#E91E63", Intergenic = "#BDBDBD",
                       Other = "#E0E0E0")

  p12 <- ggplot(annot_counts_m, aes(x = annot_simple, y = pct,
                                     fill = annot_simple)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = annot_colors_m) +
    geom_text(aes(label = paste0(round(pct, 1), "%")),
              vjust = -0.3, size = 2) +
    labs(x = NULL, y = "Percentage of peaks",
         title = paste0("Mouse bulk ATAC peak annotations (",
                        formatC(nrow(mouse_peaks), big.mark = ","),
                        " peaks)")) +
    theme_masld() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1))

  save_panel(p12, "panel_12_mouse_peak_annotation.pdf",
             w = fig_half_width, h = 3)
} else {
  message("  WARNING: ", annot_file, " not found, skipping panel 12")
}

# ===================================================================
# Panel 13: Mouse DA volcano (HFD vs Control)
# ===================================================================
message("\nPanel 13: Mouse DA volcano...")
da_mouse_file <- file.path(ATAC_MB, "results", "da_HFD_vs_Control.csv")

if (file.exists(da_mouse_file)) {
  da_m <- fread(da_mouse_file)
  da_m[, neg_log10_padj := -log10(pmin(padj, 1))]
  da_m[, sig := fifelse(padj < 0.05 & abs(logFC) > 0.5,
                         fifelse(logFC > 0, "Opened", "Closed"), "NS")]

  n_open   <- sum(da_m$sig == "Opened", na.rm = TRUE)
  n_closed <- sum(da_m$sig == "Closed", na.rm = TRUE)

  # Merge with annotations to get gene names for labeling
  if (file.exists(annot_file)) {
    annot_m <- fread(annot_file)
    da_m <- merge(da_m, annot_m[, .(peak_id, SYMBOL)],
                  by = "peak_id", all.x = TRUE)
  } else {
    da_m[, SYMBOL := NA_character_]
  }

  # Top peaks for labeling
  da_m[, rank_score := abs(logFC) * neg_log10_padj]
  top_peaks <- da_m[sig != "NS" & !is.na(SYMBOL) & SYMBOL != ""][
    order(-rank_score)][1:min(15, .N)]

  p13 <- ggplot(da_m, aes(x = logFC, y = neg_log10_padj, color = sig)) +
    geom_point(size = 0.3, alpha = 0.5, shape = 16) +
    scale_color_manual(values = c(Closed = masld_colors$down,
                                  NS     = masld_colors$ns,
                                  Opened = masld_colors$up)) +
    geom_text_repel(data = top_peaks,
                    aes(label = SYMBOL), size = 1.8,
                    max.overlaps = 20, segment.size = 0.2,
                    min.segment.length = 0) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
               linewidth = 0.3, color = "gray60") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               linewidth = 0.3, color = "gray60") +
    labs(x = expression(log[2]~"fold change (HFD / Control)"),
         y = expression(-log[10]~"(adjusted p-value)"),
         title = "Mouse bulk ATAC: HFD vs Control") +
    annotate("text", x = max(da_m$logFC, na.rm = TRUE) * 0.6,
             y = max(da_m$neg_log10_padj, na.rm = TRUE) * 0.95,
             label = paste0(n_open, " opened"), color = masld_colors$up,
             size = 2.2, hjust = 1) +
    annotate("text", x = min(da_m$logFC, na.rm = TRUE) * 0.6,
             y = max(da_m$neg_log10_padj, na.rm = TRUE) * 0.95,
             label = paste0(n_closed, " closed"), color = masld_colors$down,
             size = 2.2, hjust = 0) +
    theme_masld() +
    theme(legend.position = "none")

  save_panel(p13, "panel_13_mouse_da_volcano.pdf", w = fig_half_width, h = 3)
} else {
  message("  WARNING: ", da_mouse_file, " not found, skipping panel 13")
}

# ===================================================================
# Panel 14: Mouse promoter accessibility heatmap (top DA genes)
# ===================================================================
message("\nPanel 14: Mouse promoter accessibility heatmap...")
prom_file <- file.path(ATAC_MB, "results", "promoter_accessibility.csv")

if (file.exists(prom_file)) {
  prom <- fread(prom_file)

  # Get all DA files for diet models
  da_files <- list.files(file.path(ATAC_MB, "results"),
                         pattern = "^da_.*_vs_Control\\.csv$",
                         full.names = TRUE)

  if (length(da_files) > 0) {
    all_da <- rbindlist(lapply(da_files, function(f) {
      dt <- fread(f)
      # Keep only shared columns across contrasts
      dt[, .(peak_id, logFC, padj, contrast)]
    }))

    # Merge promoter peaks with DA results
    # peak_coordinate uses "chr:start-end" but DA peak_id uses "chr_start_end"
    if ("peak_coordinate" %in% names(prom)) {
      prom[, peak_id := gsub("[:-]", "_", peak_coordinate)]
    }
    if ("peak_id" %in% names(prom) && "peak_id" %in% names(all_da)) {
      da_prom <- merge(all_da, prom[, .(peak_id, gene_symbol)],
                       by = "peak_id", all.x = FALSE)
    } else {
      da_prom <- data.table()
    }

    if (nrow(da_prom) > 0 && "gene_symbol" %in% names(da_prom)) {
      # Get top 40 genes by max absolute logFC across any contrast
      gene_max <- da_prom[, .(max_lfc = max(abs(logFC), na.rm = TRUE)),
                          by = gene_symbol][order(-max_lfc)][1:min(40, .N)]

      da_top <- da_prom[gene_symbol %in% gene_max$gene_symbol]
      da_top[, logFC_capped := pmin(pmax(logFC, -3), 3)]

      mat_df <- dcast(da_top, gene_symbol ~ contrast,
                      value.var = "logFC_capped",
                      fun.aggregate = mean, fill = 0)

      if (ncol(mat_df) > 1) {
        gene_names <- mat_df$gene_symbol
        mat <- as.matrix(mat_df[, -1, with = FALSE])
        rownames(mat) <- gene_names

        # Cluster rows
        if (nrow(mat) > 2) {
          hc <- hclust(dist(mat))
          mat <- mat[hc$order, , drop = FALSE]
        }

        plot_df <- as.data.table(reshape2::melt(mat,
          varnames = c("Gene", "Contrast"), value.name = "logFC"))

        p14 <- ggplot(plot_df, aes(x = Contrast, y = Gene, fill = logFC)) +
          geom_tile(color = "white", linewidth = 0.2) +
          scale_fill_gradient2(low = masld_colors$down, mid = "white",
                               high = masld_colors$up, midpoint = 0,
                               name = expression(log[2]~FC),
                               limits = c(-3, 3)) +
          labs(x = NULL, y = NULL,
               title = "Mouse promoter accessibility (top DA genes)") +
          theme_masld() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
                axis.text.y = element_text(size = 4))

        save_panel(p14, "panel_14_mouse_promoter_heatmap.pdf",
                   w = fig_half_width, h = 5)
      }
    }
  }
} else {
  message("  WARNING: ", prom_file, " not found, skipping panel 14")
}

# ===================================================================
# Panel 15: SCENIC+ regulon activity vs bulk RNA LFC
# ===================================================================
message("\nPanel 15: SCENIC+ regulon-expression correlation...")
reg_file <- file.path(ATAC_HM, "scenic_plus", "disease_regulons.csv")
dream_file <- file.path(DREAM, "dream_results.csv")

if (file.exists(reg_file) && file.exists(dream_file)) {
  reg <- fread(reg_file)
  dream <- fread(dream_file)

  # Map Ensembl IDs to symbols using load_figure_data.R helper
  dream <- add_symbols(dream, gene_col = "gene")

  # Get unique regulon TFs with activity diff
  reg_unique <- unique(reg[, .(tf_name, regulon_activity_diff)])

  # Merge with dream logFC by TF symbol
  merged <- merge(reg_unique, dream[!is.na(symbol), .(symbol, logFC)],
                  by.x = "tf_name", by.y = "symbol", all.x = TRUE)
  merged <- merged[!is.na(logFC) & !is.na(regulon_activity_diff)]

  if (nrow(merged) > 3) {
    cor_test <- cor.test(merged$logFC, merged$regulon_activity_diff,
                         method = "spearman")

    p15 <- ggplot(merged, aes(x = logFC, y = regulon_activity_diff)) +
      geom_point(color = masld_colors$up, size = 2, alpha = 0.7, shape = 16) +
      geom_smooth(method = "lm", se = TRUE, color = "gray40",
                  linewidth = 0.5, linetype = "dashed") +
      geom_text_repel(aes(label = tf_name), size = 2, max.overlaps = 15,
                      segment.size = 0.2) +
      labs(x = expression("Bulk RNA integrated"~log[2]~FC),
           y = expression(Delta~"regulon activity (MASLD - Normal)"),
           title = "SCENIC+ regulon activity vs transcriptomic change") +
      annotate("text", x = min(merged$logFC) * 0.8,
               y = max(merged$regulon_activity_diff) * 0.95,
               label = paste0("rho = ", round(cor_test$estimate, 3),
                              "\np = ", format.pval(cor_test$p.value, digits = 2)),
               hjust = 0, size = 2.2, color = "gray30") +
      theme_masld()

    save_panel(p15, "panel_15_scenic_rna_correlation.pdf",
               w = fig_half_width, h = fig_half_width * 0.8)
  } else {
    message("  WARNING: Too few matched regulons for correlation, skipping")
  }
} else {
  message("  WARNING: Missing regulon or dream file, skipping panel 15")
}

# ===================================================================
# Panel 16: Cross-species promoter conservation
# ===================================================================
message("\nPanel 16: Cross-species promoter conservation...")
l8_file <- file.path(ATAC_INT, "results", "l8_atac_columns.csv")

if (file.exists(l8_file)) {
  l8 <- fread(l8_file)

  # Compare mouse vs human promoter accessibility and DA
  l8_both <- l8[!is.na(mouse_da_logFC) & !is.na(hepatocyte_da_logFC)]

  if (nrow(l8_both) > 10) {
    # Classify conservation
    l8_both[, conservation := fcase(
      cross_species_promoter_conserved == TRUE, "Conserved",
      mouse_promoter_accessible == TRUE & human_promoter_accessible == TRUE,
        "Both accessible",
      mouse_promoter_accessible == TRUE, "Mouse only",
      human_promoter_accessible == TRUE, "Human only",
      default = "Neither"
    )]

    # Correlation of DA logFC
    cor_val <- cor(l8_both$mouse_da_logFC, l8_both$hepatocyte_da_logFC,
                   use = "complete.obs", method = "spearman")

    cons_colors <- c("Conserved" = "#00695C", "Both accessible" = "#42A5F5",
                      "Mouse only" = "#F48FB1", "Human only" = "#7B1FA2",
                      "Neither" = "#BDBDBD")

    p16 <- ggplot(l8_both, aes(x = mouse_da_logFC, y = hepatocyte_da_logFC,
                                color = conservation)) +
      geom_point(size = 0.5, alpha = 0.5, shape = 16) +
      scale_color_manual(values = cons_colors, name = "Promoter status") +
      geom_smooth(method = "lm", se = FALSE, color = "gray40",
                  linewidth = 0.5, linetype = "dashed") +
      geom_hline(yintercept = 0, linewidth = 0.2, color = "gray70") +
      geom_vline(xintercept = 0, linewidth = 0.2, color = "gray70") +
      labs(x = expression("Mouse DA"~log[2]~FC~"(HFD / Control)"),
           y = expression("Human scATAC DA"~log[2]~FC~"(MASLD / Normal)"),
           title = "Cross-species chromatin accessibility") +
      annotate("text", x = min(l8_both$mouse_da_logFC, na.rm = TRUE) * 0.5,
               y = max(l8_both$hepatocyte_da_logFC, na.rm = TRUE) * 0.9,
               label = paste0("rho = ", round(cor_val, 3)),
               size = 2.5, color = "gray30") +
      theme_masld() +
      theme(legend.position = c(0.85, 0.2),
            legend.key.size = unit(0.25, "cm"))

    save_panel(p16, "panel_16_cross_species_conservation.pdf",
               w = fig_half_width, h = fig_half_width * 0.8)
  } else {
    message("  WARNING: Too few genes with both mouse+human DA, skipping")
  }
} else {
  message("  WARNING: ", l8_file, " not found, skipping panel 16")
}

cat("\n==============================================================\n")
cat("R panels complete.\n")
cat("==============================================================\n")
