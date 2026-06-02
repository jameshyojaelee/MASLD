#!/usr/bin/env Rscript
# fig_composition_features.R — 4-panel figure: deconvolution panels + feature importance
# Output: figures/supplementary/figS10_prediction/fig_composition_features.pdf + individual panels
# Layout: 2x2 grid (a-d)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
})

# --- Paths ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

RESULTS_DECONV <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/deconvolution_panel")
RESULTS_PROG   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/prognosis")
META_PATH      <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
OUT_DIR        <- FIGS10_DIR
dir.create(file.path(OUT_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

# ===========================================================================
# Panel (a): Deconvolution panel performance — grouped barplot
# ===========================================================================
pa <- tryCatch({
  perf <- read_csv(file.path(RESULTS_DECONV, "deconv_panel_performance.csv"), show_col_types = FALSE)
  stel_bin <- read_csv(file.path(RESULTS_DECONV, "deconv_panel_stellate_binary.csv"), show_col_types = FALSE)

  # Filter to the 3 panel sizes of interest (10, 25, 50)
  perf_sub <- perf %>%
    filter(panel_size %in% c(10, 25, 50)) %>%
    mutate(
      panel_size = factor(panel_size, levels = c(10, 25, 50)),
      cell_type  = factor(cell_type, levels = c("Hepatocyte", "Macrophage", "Endothelial", "Stellate"))
    )

  # Cell-type colors — colorblind-friendly, from publication theme
  ct_bar_colors <- c(
    "10" = "#42A5F5",  # Light blue
    "25" = "#1565C0",  # Deep blue
    "50" = "#0D47A1"   # Navy
  )

  # Stellate binary AUROC annotation
  stel_auroc <- stel_bin %>% filter(metric == "auroc_stellate_binary") %>% pull(value)

  p <- ggplot(perf_sub, aes(x = cell_type, y = spearman, fill = panel_size)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
    scale_fill_manual(values = ct_bar_colors, name = "Panel size") +
    scale_y_continuous(breaks = seq(-0.4, 0.8, 0.2), limits = c(-0.4, 0.75)) +
    labs(x = NULL, y = "Spearman rho") +
    # Annotate stellate binary AUROC
    annotate("text", x = 4.35, y = 0.65, size = 2,
             label = paste0("Stellate binary\nAUROC = ", round(stel_auroc, 3)),
             hjust = 1, color = "#880E4F", fontface = "italic") +
    theme_masld() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
      legend.position = c(0.15, 0.85),
      legend.background = element_blank()
    )
  p
}, error = function(e) {
  message("[Panel a] Error: ", e$message)
  placeholder("Panel a: Deconv performance\n(data not available)")
})


# ===========================================================================
# Panel (b): Predicted vs actual stellate — scatter by fibrosis stage
# ===========================================================================
pb <- tryCatch({
  pred <- read_csv(file.path(RESULTS_DECONV, "deconv_panel_predictions.csv"), show_col_types = FALSE)
  meta <- read_csv(META_PATH, show_col_types = FALSE) %>%
    select(sample_id, fibrosis_stage)

  # Filter to Stellate cell type (this is the 25-gene panel predictions)
  stel_pred <- pred %>%
    filter(cell_type == "Stellate") %>%
    left_join(meta, by = "sample_id")

  # Compute Spearman rho
  rho_val <- cor(stel_pred$predicted, stel_pred$actual, method = "spearman", use = "complete.obs")

  # Fibrosis stage as ordered factor
  stel_pred <- stel_pred %>%
    mutate(fibrosis_stage = factor(fibrosis_stage, levels = c("F0", "F1", "F2", "F3", "F4")))

  p <- ggplot(stel_pred, aes(x = actual, y = predicted, color = fibrosis_stage)) +
    rasterize_layer(geom_point(size = 0.6, alpha = 0.6)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = 0.04, linetype = "dotted", color = "#880E4F", linewidth = 0.4) +
    scale_color_manual(values = fibrosis_stage_colors, name = "Fibrosis", na.value = "grey70") +
    annotate("text", x = 0.02, y = max(stel_pred$predicted, na.rm = TRUE) * 0.95,
             label = paste0("rho = ", round(rho_val, 3)),
             hjust = 0, size = 2.2, fontface = "italic") +
    annotate("text", x = max(stel_pred$actual, na.rm = TRUE) * 0.6, y = 0.048,
             label = "4% threshold", size = 1.8, color = "#880E4F", fontface = "italic") +
    labs(x = "Actual stellate fraction", y = "Predicted stellate fraction") +
    theme_masld() +
    theme(
      legend.position = c(0.15, 0.80),
      legend.background = element_blank(),
      legend.key.size = unit(0.25, "cm")
    )
  p
}, error = function(e) {
  message("[Panel b] Error: ", e$message)
  placeholder("Panel b: Stellate scatter\n(data not available)")
})


# ===========================================================================
# Panel (c): PRS feature importance — horizontal bar chart, top 20
# ===========================================================================
pc <- tryCatch({
  feat_imp <- read_csv(file.path(RESULTS_PROG, "prs_feature_importance.csv"), show_col_types = FALSE)

  # Use M5_full config (the integrated model)
  feat_m5 <- feat_imp %>%
    filter(config == "M5_full", abs_mean_coefficient > 0) %>%
    arrange(desc(abs_mean_coefficient)) %>%
    slice_head(n = 20)

  # Assign tier based on feature prefix
  feat_m5 <- feat_m5 %>%
    mutate(
      tier = case_when(
        grepl("^div_", feature)        ~ "T1 Divergence",
        grepl("^ct_", feature)         ~ "T2 Cell-type",
        grepl("^bulkformer_|^nasvae_", feature) ~ "T3 Embedding",
        grepl("^trans_", feature)      ~ "T5 Transition",
        grepl("^clin_|^hcc_", feature) ~ "T7 Clinical",
        TRUE                           ~ "Other"
      ),
      # Clean up feature name for display
      display_name = gsub("^div_|^ct_|^bulkformer_|^nasvae_|^trans_|^clin_|^hcc_", "", feature),
      # Stability label
      stability_label = paste0(n_folds_nonzero, "/", n_folds_total)
    )

  # Tier colors — colorblind-friendly
  tier_colors_prs <- c(
    "T1 Divergence" = "#00695C",   # Teal
    "T2 Cell-type"  = "#F57F17",   # Orange
    "T3 Embedding"  = "#7B1FA2",   # Purple
    "T5 Transition" = "#2E7D32",   # Green
    "T7 Clinical"   = "#757575",   # Grey
    "Other"         = "#BDBDBD"
  )

  # Order by coefficient
  feat_m5 <- feat_m5 %>%
    mutate(display_name = factor(display_name, levels = rev(display_name)))

  p <- ggplot(feat_m5, aes(x = abs_mean_coefficient, y = display_name, fill = tier)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = stability_label), hjust = -0.1, size = 1.8, color = "grey30") +
    scale_fill_manual(values = tier_colors_prs, name = "Feature tier") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Mean |coefficient|", y = NULL) +
    theme_masld() +
    theme(
      legend.position = c(0.72, 0.25),
      legend.background = element_blank(),
      legend.key.size = unit(0.25, "cm"),
      axis.text.y = element_text(size = 5.5)
    )
  p
}, error = function(e) {
  message("[Panel c] Error: ", e$message)
  placeholder("Panel c: PRS features\n(data not available)")
})


# ===========================================================================
# Panel (d): 5-gene prognosis panel — violin plots by S1 vs S2
# ===========================================================================
pd <- tryCatch({
  feat_mat <- read_csv(file.path(RESULTS_PROG, "prognosis_feature_matrix.csv"), show_col_types = FALSE)
  labels   <- read_csv(file.path(RESULTS_PROG, "prognosis_pseudo_labels.csv"), show_col_types = FALSE)

  # 5-gene panel
  panel_genes <- c("CCDC80", "SH3YL1", "BDKRB2", "SHISA9", "LINC02099")
  panel_cols  <- paste0("div_", panel_genes)

  # Extract divergence scores for these genes
  gene_data <- feat_mat %>%
    select(sample_id, all_of(panel_cols)) %>%
    left_join(labels %>% select(sample_id, s2_binary), by = "sample_id") %>%
    filter(!is.na(s2_binary))

  # Pivot to long
  gene_long <- gene_data %>%
    pivot_longer(cols = all_of(panel_cols), names_to = "gene", values_to = "divergence") %>%
    mutate(
      gene = gsub("^div_", "", gene),
      gene = factor(gene, levels = panel_genes),
      subtype = ifelse(s2_binary == 1, "S2 (progressor)", "S1 (non-progressor)")
    )

  # Gene annotations (brief biological function)
  gene_annot <- data.frame(
    gene = factor(panel_genes, levels = panel_genes),
    func = c("ECM remodeling", "PI3K signaling", "Bradykinin receptor", "Synapse adhesion", "lncRNA"),
    stringsAsFactors = FALSE
  )

  # Stability: all 5 genes are 6/6 from M3_divergence config
  feat_imp <- read_csv(file.path(RESULTS_PROG, "prs_feature_importance.csv"), show_col_types = FALSE)
  stab_info <- feat_imp %>%
    filter(config == "M3_divergence", feature %in% panel_cols) %>%
    mutate(gene = gsub("^div_", "", feature)) %>%
    select(gene, n_folds_nonzero, n_folds_total)

  # Subtype colors
  subtype_colors <- c(
    "S1 (non-progressor)" = "#1565C0",  # Blue
    "S2 (progressor)"     = "#C2185B"   # Magenta
  )

  p <- ggplot(gene_long, aes(x = gene, y = divergence, fill = subtype)) +
    geom_violin(scale = "width", width = 0.7, linewidth = 0.2, alpha = 0.8,
                position = position_dodge(width = 0.75)) +
    geom_boxplot(width = 0.12, outlier.size = 0.3, linewidth = 0.2, alpha = 0.5,
                 position = position_dodge(width = 0.75)) +
    scale_fill_manual(values = subtype_colors, name = "Subtype") +
    # Annotate stability and function below
    geom_text(data = gene_annot, aes(x = gene, y = -Inf, label = func),
              inherit.aes = FALSE, vjust = 1.5, size = 1.6, color = "grey40",
              fontface = "italic") +
    labs(x = NULL, y = "S1-S2 divergence score") +
    theme_masld() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = 6, face = "italic"),
      legend.position = c(0.82, 0.90),
      legend.background = element_blank(),
      legend.key.size = unit(0.25, "cm"),
      plot.margin = margin(3, 3, 12, 3)  # extra bottom for annotation
    ) +
    coord_cartesian(clip = "off")
  p
}, error = function(e) {
  message("[Panel d] Error: ", e$message)
  placeholder("Panel d: 5-gene panel\n(data not available)")
})


# ===========================================================================
# Assemble composite figure
# ===========================================================================
fig <- (pa | pb) / (pc | pd) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

# Save composite
save_fig(fig, file.path(OUT_DIR, "fig_composition_features.pdf"),
         width = 180 / 25.4, height = 140 / 25.4)
message("Saved: ", file.path(OUT_DIR, "fig_composition_features.pdf"))

# Save individual panels
save_fig(pa, file.path(OUT_DIR, "panels", "panel_a_deconv_performance.pdf"),
         width = 88 / 25.4, height = 70 / 25.4)
save_fig(pb, file.path(OUT_DIR, "panels", "panel_b_stellate_scatter.pdf"),
         width = 88 / 25.4, height = 70 / 25.4)
save_fig(pc, file.path(OUT_DIR, "panels", "panel_c_prs_features.pdf"),
         width = 88 / 25.4, height = 70 / 25.4)
save_fig(pd, file.path(OUT_DIR, "panels", "panel_d_5gene_panel.pdf"),
         width = 88 / 25.4, height = 70 / 25.4)
message("Saved individual panels to: ", file.path(OUT_DIR, "panels"))

message("Done.")
