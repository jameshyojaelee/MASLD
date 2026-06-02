#!/usr/bin/env Rscript
# 04_generate_plots.R
# Generate diagnostic and result visualizations for combined MCD analysis

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

root <- normalizePath(".")
metadata_path <- file.path(root, "metadata", "samples.tsv")
vst_path <- file.path(root, "analysis_mcd_vs_control", "vst_counts.csv")
results_path <- file.path(root, "analysis_mcd_vs_control", "deseq2_mcd_vs_control.tsv")
plot_dir <- file.path(root, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

message("=== Generating Plots ===")

# --- Load data ---
samples <- read_tsv(metadata_path, col_types = cols()) %>%
  mutate(
    diet = factor(diet, levels = c("Control", "MCD")),
    batch = factor(batch)
  )

vst_data <- read_csv(vst_path, col_types = cols()) %>%
  column_to_rownames("gene") %>%
  as.matrix()

results <- read_tsv(results_path, col_types = cols())

# Align samples
common_samples <- intersect(samples$sample_id, colnames(vst_data))
samples <- samples %>% filter(sample_id %in% common_samples)
vst_data <- vst_data[, samples$sample_id]

message("Loaded ", nrow(samples), " samples, ", nrow(vst_data), " genes")

# --- 1. PCA Plot ---
message("Creating PCA plot...")
pca <- prcomp(t(vst_data), scale. = FALSE)
pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("sample_id") %>%
  left_join(samples, by = "sample_id")

var_explained <- round(100 * summary(pca)$importance[2, 1:2], 1)

# Use dataset_label for shape instead of batch
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = diet, shape = dataset_label)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("Control" = "#2E86AB", "MCD" = "#E94F37")) +
  # Map shapes: In-House (16=circle), External (17=triangle)
  scale_shape_manual(values = c("In-House MCD" = 16, "External MCD" = 17)) +
  labs(
    title = "PCA: Combined MCD Datasets (Batch-Corrected)",
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)"),
    color = "Diet",
    shape = "Dataset Source"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16)
  )

ggsave(file.path(plot_dir, "pca_batch_diet.pdf"), pca_plot, width = 10, height = 7)
message("  Saved: pca_batch_diet.pdf")

# ... [Volcano and MA plots remain same] ...

# --- 4. Sample Distance Heatmap ---
message("Creating sample heatmap...")
sample_dists <- dist(t(vst_data))
sample_dist_matrix <- as.matrix(sample_dists)

# Use dataset_label in annotation
annotation_df <- samples %>%
  select(sample_id, diet, dataset_label) %>%
  column_to_rownames("sample_id")

ann_colors <- list(
  diet = c("Control" = "#2E86AB", "MCD" = "#E94F37"),
  dataset_label = c("In-House MCD" = "#8E44AD", "External MCD" = "#F39C12")
)

pdf(file.path(plot_dir, "sample_heatmap.pdf"), width = 12, height = 10)
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  annotation_row = annotation_df,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  main = "Sample Distance Heatmap (VST)",
  fontsize = 8
)
dev.off()
message("  Saved: sample_heatmap.pdf")

message("\n=== All plots generated ===")
