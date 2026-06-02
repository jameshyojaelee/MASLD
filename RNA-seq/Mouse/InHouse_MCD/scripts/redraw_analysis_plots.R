#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
})

root <- normalizePath("./")

# Unified palette (peach & magenta)
diet_colors <- c("Control" = "#F2A45E", "MCD" = "#C23B75")
sex_shapes  <- c("Female" = 16, "Male" = 17)

# Configure PDF device
pdf.options(useDingbats = FALSE)
options(bitmapType = "cairo")

theme_clean <- function() {
  theme_bw(base_size = 7, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#efefef"),
      plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 6),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7)
    )
}

plot_ma <- function(res_tbl, outfile_prefix) {
  df <- res_tbl %>% mutate(signif = padj < 0.1 & !is.na(padj))
  p <- ggplot(df, aes(baseMean, log2FoldChange, colour = signif)) +
    geom_point(alpha = 0.5, size = 2) +
    scale_x_log10() +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#C23B75")) +
    labs(title = "MA Plot", x = "Mean Expression", y = "log2 Fold Change", colour = "padj < 0.1") +
    theme_clean()
  ggsave(paste0(outfile_prefix, "_MA.pdf"), p, width = 10, height = 8)
}

plot_volcano <- function(res_tbl, outfile_prefix) {
  df <- res_tbl %>%
    mutate(neglog10 = -log10(padj),
           significance = case_when(
             padj < 0.1 & log2FoldChange > 1  ~ "Up",
             padj < 0.1 & log2FoldChange < -1 ~ "Down",
             padj < 0.1                        ~ "Sig",
             TRUE                                ~ "NS"))
  cols <- c("Up" = "#C23B75", "Down" = "#F2A45E", "Sig" = "#E07BB6", "NS" = "grey70")
  p <- ggplot(df, aes(log2FoldChange, neglog10, colour = significance)) +
    geom_point(alpha = 0.6, size = 1.8) +
    geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
    scale_color_manual(values = cols, breaks = c("Up", "Down", "Sig", "NS")) +
    labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10(padj)") +
    theme_clean() + theme(legend.position = "bottom")
  ggsave(paste0(outfile_prefix, "_volcano.pdf"), p, width = 6, height = 5)
}

plot_pca <- function(dds_path, sample_info, outfile) {
  dds <- readRDS(dds_path)
  vsd <- vst(dds, blind = TRUE)
  pca <- plotPCA(vsd, intgroup = c("diet", "sex", "week"), returnData = TRUE)
  percentVar <- round(100 * attr(pca, "percentVar"), 1)
  p <- ggplot(pca, aes(PC1, PC2, colour = diet, shape = sex, label = week)) +
    geom_point(size = 5, alpha = 0.85) +
    geom_text_repel(size = 4, show.legend = FALSE) +
    scale_color_manual(values = diet_colors) +
    scale_shape_manual(values = sex_shapes) +
    labs(x = paste0("PC1 (", percentVar[1], "%)"), y = paste0("PC2 (", percentVar[2], "%)"), title = "PCA Plot") +
    theme_clean()
  ggsave(outfile, p, width = 10, height = 7)
}

samples <- read_tsv(file.path(root, "metadata", "samples.tsv"), show_col_types = FALSE)

analyses <- list(
  list(dir = file.path(root, "results/female_pooled"),  type = "wald"),
  list(dir = file.path(root, "results/male_pooled"),    type = "wald"),
  list(dir = file.path(root, "results/combined_pooled"), type = "wald")
)

for (a in analyses) {
  d <- a$dir
  message("Redrawing plots in ", d)
  dds_path <- file.path(d, "dds.rds")
  if (file.exists(dds_path)) {
    plot_pca(dds_path, samples, file.path(d, "PCA.pdf"))
  }
  if (a$type == "wald") {
    res_path <- file.path(d, "deseq2_results.tsv")
    if (file.exists(res_path)) {
      res_tbl <- read_tsv(res_path, show_col_types = FALSE)
      plot_ma(res_tbl, file.path(d, "deseq2"))
      plot_volcano(res_tbl, file.path(d, "deseq2"))
    }
  }
}

message("Done. Updated PCA/MA/Volcano with new colors.")
