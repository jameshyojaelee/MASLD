#!/usr/bin/env Rscript
# Cross-Dataset Comparative Analysis: Diet Mouse Models vs MCD Models
# Generates DEG count bar plots, Jaccard heatmaps, LFC correlations, UpSet plots

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(UpSetR)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(RColorBrewer)
  library(gridExtra)
  library(grid)
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================

ROOT_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
# Update output directory to use the correct path structure
OUTPUT_DIR <- file.path(ROOT_DIR, "Mouse/Public_Diet_Models/analysis/comparative")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Publication settings
pdf.options(useDingbats = FALSE)

# Colors (Diet = Blues/Greens, MCD = Reds/Pinks)
# Keys must match the NEW LABELS used in data loading
DIET_COLORS <- c(
  "20wk FPC" = "#4BAEEF",
  "7wk CDAHFD" = "#528FC5",
  "52wk HFD" = "#30D796",
  "7wk HFD" = "#40B499",
  "LIDPAD" = "#8BC163",
  "19wk HFD" = "#7CB9E8"
)

MCD_COLORS <- c(
  "MCD In-House" = "#E14B9D",
  "MCD (Yue)" = "#C23B75",
  "MCD (Paquette)" = "#D358C7"
)

ALL_COLORS <- c(DIET_COLORS, MCD_COLORS)

# Thresholds
P_CUTOFF <- 0.1
LFC_CUTOFF <- 0.5  # |log2FC| > 0.5

# ==============================================================================
# DATA LOADING
# ==============================================================================

cat("========================================\n")
cat("Cross-Dataset Comparative Analysis\n")
cat("========================================\n\n")

# Diet Mouse Results
# Keys are the CLEAN LABELS to be used in plots
diet_files <- list(
  "20wk FPC" = file.path(ROOT_DIR, "Mouse/Public_Diet_Models/analysis/deseq2/GSE162876_20wks_FPC_vs_LFD_results.csv"),
  "7wk CDAHFD" = file.path(ROOT_DIR, "Mouse/Public_Diet_Models/analysis/deseq2/GSE162876_7wks_CDAHFD_vs_LFD_results.csv"),
  "52wk HFD" = file.path(ROOT_DIR, "Mouse/Public_Diet_Models/analysis/deseq2/GSE274914_52w_HFD_vs_LFD_results.csv"),
  "7wk HFD" = file.path(ROOT_DIR, "Mouse/Public_Diet_Models/analysis/deseq2/GSE274914_7w_HFD_vs_LFD_results.csv"),
  "LIDPAD" = file.path(ROOT_DIR, "Mouse/Public_Diet_Models/analysis/deseq2/GSE159911_LIDPAD_vs_Control_results.csv"),
  "19wk HFD" = file.path(ROOT_DIR, "Mouse/Public_Diet_Models/analysis/deseq2/GSE224069_HFD_vs_Chow_results.csv")
)

# MCD Results
# Keys are the CLEAN LABELS
mcd_files <- list(
  "MCD In-House" = file.path(ROOT_DIR, "Mouse/InHouse_MCD/results/combined_pooled/deseq2_results.tsv"),
  "MCD (Yue)" = file.path(ROOT_DIR, "Mouse/Public_MCD/GSE205974/analysis_mcd_vs_control/deseq2_mcd_vs_control.tsv"),
  "MCD (Paquette)" = file.path(ROOT_DIR, "Mouse/Public_MCD/GSE156918/analysis_mcd_vs_control/deseq2_mcd_vs_control.tsv")
)

# Load all datasets
load_deseq_results <- function(file_path, name) {
  if (!file.exists(file_path)) {
    cat("  Warning: File not found:", file_path, "\n")
    return(NULL)
  }
  
  # Detect file type
  if (grepl("\\.csv$", file_path)) {
    df <- read.csv(file_path, stringsAsFactors = FALSE)
  } else {
    df <- read.delim(file_path, stringsAsFactors = FALSE)
  }
  
  # Standardize columns
  if (!"gene_id" %in% colnames(df)) {
    if ("gene" %in% colnames(df)) {
      df$gene_id <- df$gene  # MCD files use 'gene' column
    } else if ("X" %in% colnames(df)) {
      df$gene_id <- df$X
    } else if (!is.null(rownames(df)) && !identical(rownames(df), as.character(1:nrow(df)))) {
      df$gene_id <- rownames(df)
    } else {
      cat("  Warning: No gene ID column found in", name, "\n")
      return(NULL)
    }
  }
  
  # Clean gene_id (remove version)
  df$gene_id_clean <- sub("\\..*", "", df$gene_id)
  df$dataset <- name
  
  return(df)
}

cat("Loading Diet Mouse Results...\n")
diet_data <- lapply(names(diet_files), function(n) {
  cat("  Loading:", n, "\n")
  load_deseq_results(diet_files[[n]], n)
})
names(diet_data) <- names(diet_files)
diet_data <- diet_data[!sapply(diet_data, is.null)]

cat("\nLoading MCD Results...\n")
mcd_data <- lapply(names(mcd_files), function(n) {
  cat("  Loading:", n, "\n")
  load_deseq_results(mcd_files[[n]], n)
})
names(mcd_data) <- names(mcd_files)
mcd_data <- mcd_data[!sapply(mcd_data, is.null)]

all_data <- c(diet_data, mcd_data)

cat("\nLoaded", length(all_data), "datasets\n\n")

# ==============================================================================
# 1. DEG COUNT COMPARISON
# ==============================================================================

cat("Generating DEG count comparison...\n")

count_degs <- function(df) {
  df <- df %>%
    filter(!is.na(padj) & !is.na(log2FoldChange))
  
  n_up <- sum(df$padj < P_CUTOFF & df$log2FoldChange > LFC_CUTOFF, na.rm = TRUE)
  n_down <- sum(df$padj < P_CUTOFF & df$log2FoldChange < -LFC_CUTOFF, na.rm = TRUE)
  
  return(data.frame(up = n_up, down = n_down))
}

deg_counts <- do.call(rbind, lapply(names(all_data), function(n) {
  counts <- count_degs(all_data[[n]])
  data.frame(
    dataset = n,
    direction = c("Up", "Down"),
    count = c(counts$up, counts$down),
    type = ifelse(grepl("MCD", n), "MCD", "Diet")
  )
}))

# Order datasets
deg_counts$dataset <- factor(deg_counts$dataset, 
                              levels = c(names(diet_files), names(mcd_files)))

# Bar plot
p_deg <- ggplot(deg_counts, aes(x = dataset, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(
    values = c(
      "Up" = "#C23B75",
      "Down" = "#F2A45E"
    ),
    labels = c("Upregulated", "Downregulated")
  ) +
  theme_minimal(base_family = "Helvetica", base_size = 7) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.title = element_text(size = 8),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "DEG Counts: Diet Models vs MCD Models",
    subtitle = sprintf("|LFC| > %.1f, padj < %.2f", LFC_CUTOFF, P_CUTOFF),
    x = NULL,
    y = "Number of DEGs"
  )

cairo_pdf(file.path(OUTPUT_DIR, "deg_counts_comparison.pdf"), width = 5, height = 4, family = "Helvetica")
print(p_deg)
dev.off()
cat("  Created: deg_counts_comparison.pdf\n")

# Save counts table
deg_summary <- deg_counts %>%
  pivot_wider(names_from = direction, values_from = count) %>%
  mutate(total = Up + Down)
write.csv(deg_summary, file.path(OUTPUT_DIR, "deg_counts_summary.csv"), row.names = FALSE)

# ==============================================================================
# 2. JACCARD SIMILARITY HEATMAP
# ==============================================================================

cat("Calculating Jaccard similarity...\n")

# Get DEG sets
get_deg_set <- function(df, direction = "up") {
  df <- df %>% filter(!is.na(padj) & !is.na(log2FoldChange))
  
  if (direction == "up") {
    df %>% filter(padj < P_CUTOFF & log2FoldChange > LFC_CUTOFF) %>% pull(gene_id_clean)
  } else {
    df %>% filter(padj < P_CUTOFF & log2FoldChange < -LFC_CUTOFF) %>% pull(gene_id_clean)
  }
}

# Calculate Jaccard for all pairs
jaccard <- function(a, b) {
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

datasets <- names(all_data)
n <- length(datasets)

# Upregulated genes
up_sets <- lapply(all_data, get_deg_set, direction = "up")
jaccard_up <- matrix(0, n, n, dimnames = list(datasets, datasets))
for (i in 1:n) {
  for (j in 1:n) {
    jaccard_up[i, j] <- jaccard(up_sets[[i]], up_sets[[j]])
  }
}

# Downregulated genes
down_sets <- lapply(all_data, get_deg_set, direction = "down")
jaccard_down <- matrix(0, n, n, dimnames = list(datasets, datasets))
for (i in 1:n) {
  for (j in 1:n) {
    jaccard_down[i, j] <- jaccard(down_sets[[i]], down_sets[[j]])
  }
}

# Row annotations
row_annot <- data.frame(
  Type = ifelse(grepl("MCD", datasets), "MCD", "Diet"),
  row.names = datasets
)

# Heatmaps
cairo_pdf(file.path(OUTPUT_DIR, "jaccard_heatmap_upregulated.pdf"), width = 6, height = 5, family = "Helvetica")
pheatmap(jaccard_up, 
         main = "Jaccard Similarity: Upregulated DEGs",
         color = colorRampPalette(c("white", "#C23B75"))(50),
         annotation_row = row_annot,
         annotation_col = row_annot,
         fontsize = 7,
         fontsize_row = 6,
         fontsize_col = 6,
         display_numbers = TRUE,
         number_format = "%.2f",
         number_color = "black",
         fontsize_number = 5)
dev.off()
cat("  Created: jaccard_heatmap_upregulated.pdf\n")

cairo_pdf(file.path(OUTPUT_DIR, "jaccard_heatmap_downregulated.pdf"), width = 6, height = 5, family = "Helvetica")
pheatmap(jaccard_down, 
         main = "Jaccard Similarity: Downregulated DEGs",
         color = colorRampPalette(c("white", "#F2A45E"))(50),
         annotation_row = row_annot,
         annotation_col = row_annot,
         fontsize = 7,
         fontsize_row = 6,
         fontsize_col = 6,
         display_numbers = TRUE,
         number_format = "%.2f",
         number_color = "black",
         fontsize_number = 5)
dev.off()
cat("  Created: jaccard_heatmap_downregulated.pdf\n")

# Save matrices
write.csv(jaccard_up, file.path(OUTPUT_DIR, "jaccard_matrix_upregulated.csv"))
write.csv(jaccard_down, file.path(OUTPUT_DIR, "jaccard_matrix_downregulated.csv"))

# ==============================================================================
# 3. LFC CORRELATION WITH MCD POOLED
# ==============================================================================

cat("Calculating LFC correlations with MCD Pooled...\n")

mcd_pooled <- all_data[["MCD In-House"]]
if (!is.null(mcd_pooled)) {
  
  mcd_lfc <- mcd_pooled %>%
    dplyr::select(gene_id_clean, log2FoldChange) %>%
    rename(lfc_mcd = log2FoldChange) %>%
    filter(!is.na(lfc_mcd))
  
  correlations <- data.frame(
    dataset = character(),
    pearson = numeric(),
    spearman = numeric(),
    n_genes = integer(),
    stringsAsFactors = FALSE
  )
  
  for (name in names(diet_data)) {
    df <- diet_data[[name]]
    diet_lfc <- df %>%
      dplyr::select(gene_id_clean, log2FoldChange) %>%
      rename(lfc_diet = log2FoldChange) %>%
      filter(!is.na(lfc_diet))
    
    merged <- inner_join(mcd_lfc, diet_lfc, by = "gene_id_clean")
    
    if (nrow(merged) > 100) {
      # Filter outliers for visualization
      merged <- merged %>%
        filter(abs(lfc_mcd) < 10 & abs(lfc_diet) < 10)
      
      cor_p <- cor(merged$lfc_mcd, merged$lfc_diet, method = "pearson")
      cor_s <- cor(merged$lfc_mcd, merged$lfc_diet, method = "spearman")
      
      correlations <- rbind(correlations, data.frame(
        dataset = name,
        pearson = cor_p,
        spearman = cor_s,
        n_genes = nrow(merged)
      ))
      
      # Scatter plot
      p <- ggplot(merged, aes(x = lfc_mcd, y = lfc_diet)) +
        geom_point(alpha = 0.3, size = 0.5, color = ALL_COLORS[[name]]) +
        geom_smooth(method = "lm", color = "#C23B75", se = FALSE, linewidth = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        theme_minimal(base_family = "Helvetica", base_size = 7) +
        labs(
          title = paste(name, "vs MCD In-House"),
          subtitle = sprintf("r = %.3f, ρ = %.3f (n = %d genes)", cor_p, cor_s, nrow(merged)),
          x = "MCD In-House LFC",
          y = paste(name, "LFC")
        ) +
        theme(
          plot.title = element_text(size = 8, face = "bold"),
          plot.subtitle = element_text(size = 6)
        )
      
      # Clean filename for saving (remove spaces/parens)
      clean_name <- gsub("[ ()]", "_", name)
      clean_name <- gsub("__", "_", clean_name)
      cairo_pdf(file.path(OUTPUT_DIR, paste0("lfc_correlation_", clean_name, ".pdf")), 
                width = 4, height = 4, family = "Helvetica")
      print(p)
      dev.off()
    }
  }
  
  write.csv(correlations, file.path(OUTPUT_DIR, "lfc_correlations_summary.csv"), row.names = FALSE)
  cat("  Created LFC correlation plots and summary\n")
}

# ==============================================================================
# 4. UPSET PLOT
# ==============================================================================

cat("Generating UpSet plot...\n")

# Combine up + down for UpSet
all_deg_sets <- lapply(all_data, function(df) {
  df %>%
    filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
    filter(padj < P_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF) %>%
    pull(gene_id_clean) %>%
    unique()
})

# Filter to non-empty sets
all_deg_sets <- all_deg_sets[sapply(all_deg_sets, length) > 0]

if (length(all_deg_sets) >= 2) {
  cairo_pdf(file.path(OUTPUT_DIR, "upset_plot.pdf"), width = 8, height = 5, family = "Helvetica")
  print(upset(fromList(all_deg_sets), 
              nsets = length(all_deg_sets),
              order.by = "freq",
              decreasing = TRUE,
              mainbar.y.label = "Intersection Size",
              sets.x.label = "DEGs per Dataset",
              text.scale = c(1.3, 1.3, 1.0, 1.0, 1.3, 1.0),
              point.size = 3,
              line.size = 1))
  dev.off()
  cat("  Created: upset_plot.pdf\n")
}

# ==============================================================================
# 5. PCA OF LFC PROFILES
# ==============================================================================

cat("Generating PCA of LFC profiles...\n")

# Build LFC matrix
all_lfc <- lapply(names(all_data), function(name) {
  df <- all_data[[name]]
  df %>%
    dplyr::select(gene_id_clean, log2FoldChange) %>%
    rename(!!name := log2FoldChange)
})

lfc_matrix <- Reduce(function(x, y) full_join(x, y, by = "gene_id_clean"), all_lfc)
rownames(lfc_matrix) <- lfc_matrix$gene_id_clean
lfc_matrix <- lfc_matrix[, -1]

# Filter genes with data in most datasets and high variance
lfc_matrix <- lfc_matrix[rowSums(!is.na(lfc_matrix)) >= ncol(lfc_matrix) * 0.7, ]
lfc_matrix[is.na(lfc_matrix)] <- 0

# Filter top variable genes
vars <- apply(lfc_matrix, 1, var)
lfc_matrix <- lfc_matrix[order(vars, decreasing = TRUE)[1:min(5000, nrow(lfc_matrix))], ]

# Transpose for PCA (columns = genes, rows = datasets)
lfc_t <- t(as.matrix(lfc_matrix))

# PCA
pca <- prcomp(lfc_t, scale. = TRUE, center = TRUE)
pca_df <- as.data.frame(pca$x[, 1:2])
pca_df$dataset <- rownames(pca_df)
pca_df$type <- ifelse(grepl("MCD", pca_df$dataset), "MCD", "Diet")

var_explained <- summary(pca)$importance[2, 1:2] * 100

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = type, label = dataset)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, hjust = 0.5, size = 2.5) +
  scale_color_manual(values = c("Diet" = "#4BAEEF", "MCD" = "#C23B75")) +
  theme_minimal(base_family = "Helvetica", base_size = 7) +
  labs(
    title = "PCA of LFC Profiles",
    subtitle = "Based on top 5000 variable genes",
    x = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
    y = sprintf("PC2 (%.1f%% variance)", var_explained[2]),
    color = "Model Type"
  ) +
  theme(
    plot.title = element_text(size = 8, face = "bold"),
    legend.position = "bottom"
  )

cairo_pdf(file.path(OUTPUT_DIR, "pca_lfc_profiles.pdf"), width = 5, height = 4, family = "Helvetica")
print(p_pca)
dev.off()
cat("  Created: pca_lfc_profiles.pdf\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n========================================\n")
cat("Analysis Complete!\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("========================================\n")
