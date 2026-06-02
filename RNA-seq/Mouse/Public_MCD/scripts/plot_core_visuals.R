#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

root <- normalizePath(".")
analysis_dir <- file.path(root, "analysis_mcd_vs_control")
res_path <- file.path(analysis_dir, "deseq2_mcd_vs_control.tsv")
norm_path <- file.path(analysis_dir, "normalized_counts.csv")
meta_path <- file.path(root, "metadata", "samples.tsv")
plots_dir <- file.path(analysis_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(res_path)) {
  stop("Missing DESeq2 results at ", res_path)
}
if (!file.exists(norm_path)) {
  stop("Missing normalized counts at ", norm_path)
}
if (!file.exists(meta_path)) {
  stop("Missing metadata at ", meta_path)
}

res <- read_tsv(res_path, col_types = cols())
norm_counts <- read_csv(norm_path, col_types = cols())
meta <- read_tsv(meta_path, col_types = cols())

label_col <- if ("gene_name" %in% colnames(res)) {
  ifelse(is.na(res$gene_name) | res$gene_name == "", res$gene, res$gene_name)
} else {
  res$gene
}

res <- res %>%
  mutate(
    label = label_col,
    padj = ifelse(is.na(padj), 1, padj),
    log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
    sig = padj < 0.1 & abs(log2FoldChange) > 1,
    neglog10 = -log10(padj),
    significance = case_when(
      padj < 0.1 & log2FoldChange > 1 ~ "Up (padj<0.1, log2FC>1)",
      padj < 0.1 & log2FoldChange < -1 ~ "Down (padj<0.1, log2FC<-1)",
      padj < 0.1 ~ "Significant (padj<0.1)",
      TRUE ~ "Not significant"
    )
  )

# Configure PDF device
pdf.options(useDingbats = FALSE)
options(bitmapType = "cairo")

save_plot <- function(plot, filename, width = 7, height = 6) {
  ggsave(filename, plot = plot, width = width, height = height, units = "in", dpi = 300)
}

# Volcano plot (match in-house week-pooled combined styling)
top_labels <- res %>%
  filter(abs(log2FoldChange) > 2) %>%
  arrange(padj) %>%
  slice_head(n = 15)

volcano <- res %>%
  ggplot(aes(x = log2FoldChange, y = neglog10)) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1.0) +
  geom_point(
    data = top_labels,
    aes(x = log2FoldChange, y = neglog10),
    color = "black", size = 2, shape = 21, fill = "yellow", stroke = 0.5
  ) +
  ggrepel::geom_text_repel(
    data = top_labels,
    aes(label = label),
    size = 2,
    fontface = "bold",
    box.padding = 0.5,
    min.segment.length = 0
  ) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "gray30", linewidth = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray30", linewidth = 0.5) +
  scale_color_manual(
    values = c(
      "Up (padj<0.1, log2FC>1)" = "#C23B75",
      "Down (padj<0.1, log2FC<-1)" = "#F2A45E",
      "Significant (padj<0.1)" = "#E07BB6",
      "Not significant" = "gray70"
    ),
    breaks = c("Up (padj<0.1, log2FC>1)", "Down (padj<0.1, log2FC<-1)", "Significant (padj<0.1)", "Not significant")
  ) +
  theme_bw(base_size = 7, base_family = "Helvetica") +
  labs(
    title = "Volcano Plot", x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "Significance"
  ) +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

save_plot(volcano, file.path(plots_dir, "volcano.pdf"), width = 7, height = 6)

# MA plot
if ("baseMean" %in% colnames(res)) {
  ma <- res %>%
    ggplot(aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = sig), alpha = 0.7, size = 1.2) +
    scale_x_log10() +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#2C7FB8")) +
    theme_minimal(base_size = 12) +
    labs(
      title = "MA: MCD vs Control",
      x = "mean of normalized counts",
      y = "log2 fold change",
      color = "Significant"
    )
  save_plot(ma, file.path(plots_dir, "ma.pdf"))
}

# PCA
count_mat <- norm_counts %>%
  column_to_rownames("gene") %>%
  as.matrix()
log_counts <- log2(count_mat + 1)
sample_ids <- colnames(log_counts)
gene_vars <- apply(log_counts, 1, var, na.rm = TRUE)
log_counts_pca <- log_counts[gene_vars > 0, , drop = FALSE]

if (nrow(log_counts_pca) >= 2) {
  pca <- prcomp(t(log_counts_pca), scale. = TRUE)
  pca_df <- as.data.frame(pca$x[, 1:2]) %>%
    rownames_to_column("sample_id") %>%
    left_join(meta, by = c("sample_id" = "sample_id"))
  percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = diet, shape = genotype)) +
    geom_point(size = 3, alpha = 0.9) +
    theme_minimal(base_size = 12) +
    labs(
      title = "PCA (log2 normalized counts)",
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)")
    )
  save_plot(pca_plot, file.path(plots_dir, "pca.pdf"))
} else {
  message("Skipping PCA: insufficient variable genes after filtering.")
}

# Sample distance heatmap
vars <- apply(log_counts, 1, var, na.rm = TRUE)
top_var_genes <- names(sort(vars, decreasing = TRUE))[seq_len(min(500, length(vars)))]
dist_mat <- dist(t(log_counts[top_var_genes, , drop = FALSE]))
annotation <- meta %>%
  select(any_of(c("sample_id", "diet", "genotype"))) %>%
  column_to_rownames("sample_id")
  pdf(file.path(plots_dir, "sample_distance_heatmap.pdf"), width = 8, height = 7)
pheatmap(
  as.matrix(dist_mat),
  annotation_col = annotation,
  annotation_row = annotation,
  main = "Sample distance (top variable genes)"
)
dev.off()

# Top DE heatmap
top_de <- res %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 50) %>%
  pull(gene)
top_mat <- log_counts[top_de, , drop = FALSE]
pdf(file.path(plots_dir, "top_de_heatmap.pdf"), width = 8, height = 8)
pheatmap(
  top_mat,
  scale = "row",
  annotation_col = annotation,
  main = "Top DE genes (padj)"
)
dev.off()
