#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_PDF <- file.path(FIG4_SC_DIR, "panels", "supplementary", "pseudobulk_bulk.pdf")
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

scvi_dir <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2")
pb_de_dir   <- file.path(scvi_dir, "pseudobulk_de")
hep_de_file <- file.path(pb_de_dir, "Hepatocytes_de.csv")

dream <- load_dream_results()

if (file.exists(hep_de_file) && !is.null(dream)) {
  hep_de <- fread(hep_de_file)
  hep_de[, ensembl_clean := sub("\\..*", "", gene)]

  dream_slim <- copy(dream)
  bulk_padj_col <- intersect(c("bulk_padj", "padj"), names(dream_slim))[1]
  bulk_lfc_col  <- intersect(c("bulk_logFC", "logFC"), names(dream_slim))[1]
  dream_slim[, ensembl_clean := sub("\\..*", "", gene)]
  dream_slim <- dream_slim[, .(ensembl_clean,
                                bulk_lfc = get(bulk_lfc_col),
                                bulk_padj = get(bulk_padj_col))]

  merged_f <- merge(
    hep_de[, .(ensembl_clean, hep_lfc = logFC, hep_padj = padj)],
    dream_slim,
    by = "ensembl_clean"
  )

  # Calculate Spearman correlation
  valid <- merged_f[is.finite(hep_lfc) & is.finite(bulk_lfc)]
  cor_s <- cor(valid$hep_lfc, valid$bulk_lfc, method = "spearman", use = "complete.obs")
  
  cor_label <- sprintf("Spearman rho = %.3f\n(n=%s genes)", 
                       cor_s, format(nrow(valid), big.mark = ","))

  lim_f_hep <- max(abs(valid$hep_lfc), na.rm = TRUE)
  lim_f_dream <- max(abs(valid$bulk_lfc), na.rm = TRUE)

  # Trim extreme 1% for visualization limits to prevent squishing
  lim_f_hep <- quantile(abs(valid$hep_lfc), 0.99, na.rm=TRUE)
  lim_f_dream <- quantile(abs(valid$bulk_lfc), 0.99, na.rm=TRUE)

  P_CUTOFF_SC <- 0.1
  P_CUTOFF_BULK <- 0.1

  # Color by concordance
  merged_f[, quad_class := fcase(
      (!is.na(hep_padj) & hep_padj < P_CUTOFF_SC & !is.na(bulk_padj) & bulk_padj < P_CUTOFF_BULK), "Both DE",
      (!is.na(hep_padj) & hep_padj < P_CUTOFF_SC), "Pseudobulk only",
      (!is.na(bulk_padj) & bulk_padj < P_CUTOFF_BULK), "Bulk only",
      default = "Non-sig"
  )]

  sig_colors <- c("Both DE" = "#C2185B", 
                  "Pseudobulk only" = "#42A5F5",
                  "Bulk only" = "#E1BEE7",
                  "Non-sig" = "gray85")
                  
  merged_f[, quad_class := factor(quad_class, levels=c("Non-sig", "Bulk only", "Pseudobulk only", "Both DE"))]
  merged_f <- merged_f[order(quad_class)]

  p_f <- ggplot(merged_f, aes(x = hep_lfc, y = bulk_lfc, color = quad_class)) +
    # Draw reference axes
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "gray30") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "gray30") +
    geom_point(size = 0.8, alpha = 0.9, shape = 16, stroke = 0) +
    scale_color_manual(values = sig_colors, name = NULL) +
    xlim(c(-lim_f_hep, lim_f_hep)) + ylim(c(-lim_f_dream, lim_f_dream)) +
    annotate("label", x = -lim_f_hep * 0.9, y = lim_f_dream * 0.9, label = cor_label,
             hjust = 0, vjust=1, size = 2.5, fill = "white", label.padding = unit(0.2, "lines")) +
    labs(x = expression("scRNA-seq (Pseudobulk) log"[2]*"FC"),
         y = expression("Bulk RNA-seq log"[2]*"FC"),
         title = sprintf("Hepatocyte vs Bulk LFC Concordance\n(padj < %s)", P_CUTOFF_SC)) +
    theme_masld() +
    theme(legend.key.size = unit(0.2, "cm"),
          legend.text = element_text(size = 6),
          legend.position = "bottom",
          aspect.ratio = 1) +
    guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1), nrow=2))
    
    save_fig(p_f, OUT_PDF, width = 3.5, height = 4.0)
    cat("Successfully saved updated plot to", OUT_PDF, "\n")
} else {
    cat("Error loading data\n")
}
