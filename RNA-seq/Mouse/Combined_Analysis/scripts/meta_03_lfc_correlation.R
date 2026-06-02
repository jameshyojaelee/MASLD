#!/usr/bin/env Rscript
# meta_03_lfc_correlation.R
# Compare log2FC between combined and individual datasets

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})

root <- normalizePath(".")
out_dir <- file.path(root, "meta_analysis")
plot_dir <- file.path(root, "plots")

message("=== LFC Correlation Analysis ===\n")

# Load DEG lists
all_degs <- readRDS(file.path(out_dir, "all_degs_list.rds"))

# Prepare combined dataset for joining
combined <- all_degs$Combined %>%
  select(gene_id, combined_lfc = log2FoldChange, combined_padj = padj, 
         combined_baseMean = baseMean, gene_name) %>%
  filter(!is.na(combined_lfc))

# Function to create scatter plot
create_lfc_scatter <- function(df, x_col, y_col, x_label, y_label, 
                                out_file, highlight_diff = 2) {
  # Filter to genes in both
  df <- df %>% filter(!is.na(!!sym(x_col)), !is.na(!!sym(y_col)))
  
  # Compute correlation
  cor_val <- cor(df[[x_col]], df[[y_col]], method = "spearman", use = "complete.obs")
  
  # Identify genes with large differences
  df <- df %>%
    mutate(
      diff = abs(!!sym(x_col) - !!sym(y_col)),
      highlight = diff > highlight_diff & combined_padj < 0.1 & combined_lfc > 0
    )
  
  top_diff <- df %>%
    filter(highlight) %>%
    arrange(desc(diff)) %>%
    head(10)
  
  p <- ggplot(df, aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_point(aes(color = combined_padj < 0.1 & combined_lfc > 0), 
               alpha = 0.4, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#E94F37"),
                       labels = c("FALSE" = "Not sig up", "TRUE" = "Sig up in Combined"),
                       name = "") +
    geom_text_repel(data = top_diff, aes(label = gene_name), 
                    size = 3, max.overlaps = 15) +
    labs(
      title = paste0("LFC Correlation: ", x_label, " vs ", y_label),
      subtitle = sprintf("Spearman ρ = %.3f (n = %d genes)", cor_val, nrow(df)),
      x = paste0("log2FC (", x_label, ")"),
      y = paste0("log2FC (", y_label, ")")
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom") +
    coord_fixed(ratio = 1) +
    xlim(c(-10, 15)) + ylim(c(-10, 15))
  
  ggsave(out_file, p, width = 8, height = 8)
  message("Saved: ", out_file)
  
  return(tibble(comparison = paste0(x_label, " vs ", y_label), 
                spearman_rho = cor_val, n_genes = nrow(df)))
}

# Compare combined vs each individual
correlations <- list()

# vs InHouse
if ("InHouse" %in% names(all_degs)) {
  inhouse <- all_degs$InHouse %>%
    select(gene_id, inhouse_lfc = log2FoldChange)
  
  merged <- combined %>% inner_join(inhouse, by = "gene_id")
  correlations[[1]] <- create_lfc_scatter(
    merged, "inhouse_lfc", "combined_lfc",
    "InHouse", "Combined",
    file.path(plot_dir, "lfc_inhouse_vs_combined.pdf")
  )
}

# vs GSE156918
if ("GSE156918" %in% names(all_degs)) {
  gse156918 <- all_degs$GSE156918 %>%
    select(gene_id, gse156918_lfc = log2FoldChange)
  
  merged <- combined %>% inner_join(gse156918, by = "gene_id")
  correlations[[2]] <- create_lfc_scatter(
    merged, "gse156918_lfc", "combined_lfc",
    "GSE156918", "Combined",
    file.path(plot_dir, "lfc_gse156918_vs_combined.pdf")
  )
}

# vs GSE205974
if ("GSE205974" %in% names(all_degs)) {
  gse205974 <- all_degs$GSE205974 %>%
    select(gene_id, gse205974_lfc = log2FoldChange)
  
  merged <- combined %>% inner_join(gse205974, by = "gene_id")
  correlations[[3]] <- create_lfc_scatter(
    merged, "gse205974_lfc", "combined_lfc",
    "GSE205974", "Combined",
    file.path(plot_dir, "lfc_gse205974_vs_combined.pdf")
  )
}

# Save correlation summary
cor_summary <- bind_rows(correlations)
write_csv(cor_summary, file.path(out_dir, "lfc_correlations.csv"))
message("\n=== LFC Correlation Summary ===")
print(cor_summary)
