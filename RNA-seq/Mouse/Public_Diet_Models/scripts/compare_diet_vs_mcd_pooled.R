#!/usr/bin/env Rscript
# Cross-Dataset Comparative Analysis: Diet Mouse Models vs MCD Models (Pooled Only)
# Compares diet models to MCD_Pooled datasets only (no individual weeks)
# Includes: in-house MCD pooled, external MCD (GSE156918+GSE205974) pooled

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(UpSetR)
  library(RColorBrewer)
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================

ROOT_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
OUTPUT_DIR <- file.path(ROOT_DIR, "other_diet_mouse_RNAseq/analysis/comparative")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Publication settings
pdf.options(useDingbats = FALSE)

# Colors - dataset labels include diet type
DIET_COLORS <- c(
  "GSE162876_FPC_20wk" = "#4BAEEF",
  "GSE162876_CDAHFD_7wk" = "#528FC5",
  "GSE274914_HFD_52w" = "#30D796",
  "GSE274914_HFD_7w" = "#40B499",
  "GSE159911_LIDPAD" = "#8BC163",
  "GSE224069_HFD" = "#7CB9E8"
)

MCD_COLORS <- c(
  "MCD_InHouse" = "#C23B75",
  "MCD_External" = "#E14B9D"
)

ALL_COLORS <- c(DIET_COLORS, MCD_COLORS)

# Thresholds (relaxed for HFD datasets with smaller effect sizes)
P_CUTOFF <- 0.1   # padj threshold
LFC_CUTOFF <- 0.5 # log2 fold change threshold

# ==============================================================================
# DATA LOADING
# ==============================================================================

cat("========================================\n")
cat("Cross-Dataset Comparative Analysis\n")
cat("(Diet Models vs MCD Pooled Only)\n")
cat("========================================\n\n")

# Diet Mouse Results
# Excluded: GSE263273 (pharmacological interventions), GSE225616 (siRNA vs Vehicle within GAN diet)
# Labels include diet type for clarity
diet_files <- list(
  "GSE162876_FPC_20wk" = file.path(ROOT_DIR, "other_diet_mouse_RNAseq/analysis/deseq2/GSE162876_20wks_FPC_vs_LFD_results.csv"),
  "GSE162876_CDAHFD_7wk" = file.path(ROOT_DIR, "other_diet_mouse_RNAseq/analysis/deseq2/GSE162876_7wks_CDAHFD_vs_LFD_results.csv"),
  "GSE274914_HFD_52w" = file.path(ROOT_DIR, "other_diet_mouse_RNAseq/analysis/deseq2/GSE274914_52w_HFD_vs_LFD_results.csv"),
  "GSE274914_HFD_7w" = file.path(ROOT_DIR, "other_diet_mouse_RNAseq/analysis/deseq2/GSE274914_7w_HFD_vs_LFD_results.csv"),
  "GSE159911_LIDPAD" = file.path(ROOT_DIR, "other_diet_mouse_RNAseq/analysis/deseq2/GSE159911_LIDPAD_vs_Control_results.csv"),
  "GSE224069_HFD" = file.path(ROOT_DIR, "other_diet_mouse_RNAseq/analysis/deseq2/GSE224069_HFD_vs_Chow_results.csv")
)

# MCD Results (Pooled Only)
mcd_files <- list(
  "MCD_InHouse" = file.path(ROOT_DIR, "in-house_MCD_RNAseq/results/combined_pooled/deseq2_results.tsv"),
  "MCD_External" = file.path(ROOT_DIR, "other_MCD_RNAseq/analysis_mcd_vs_control/deseq2_mcd_vs_control.tsv")
)

# Load datasets
load_deseq_results <- function(file_path, name) {
  if (!file.exists(file_path)) {
    cat("  Warning: File not found:", file_path, "\n")
    return(NULL)
  }
  
  if (grepl("\\.csv$", file_path)) {
    df <- read.csv(file_path, stringsAsFactors = FALSE)
  } else {
    df <- read.delim(file_path, stringsAsFactors = FALSE)
  }
  
  # Standardize gene_id column
  if (!"gene_id" %in% colnames(df)) {
    if ("gene" %in% colnames(df)) {
      df$gene_id <- df$gene
    } else if ("X" %in% colnames(df)) {
      df$gene_id <- df$X
    } else if (!is.null(rownames(df)) && !identical(rownames(df), as.character(1:nrow(df)))) {
      df$gene_id <- rownames(df)
    } else {
      cat("  Warning: No gene ID column found in", name, "\n")
      return(NULL)
    }
  }
  
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

cat("\nLoading MCD Pooled Results...\n")
mcd_data <- lapply(names(mcd_files), function(n) {
  cat("  Loading:", n, "\n")
  load_deseq_results(mcd_files[[n]], n)
})
names(mcd_data) <- names(mcd_files)
mcd_data <- mcd_data[!sapply(mcd_data, is.null)]

all_data <- c(diet_data, mcd_data)

cat("\nLoaded", length(all_data), "datasets\n")
cat("  Diet models:", length(diet_data), "\n")
cat("  MCD pooled:", length(mcd_data), "\n\n")

# ==============================================================================
# 1. DEG COUNT COMPARISON
# ==============================================================================

cat("Generating DEG count comparison...\n")

count_degs <- function(df) {
  df <- df %>% filter(!is.na(padj) & !is.na(log2FoldChange))
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

deg_counts$dataset <- factor(deg_counts$dataset, levels = c(names(diet_files), names(mcd_files)))

# All datasets get same Up/Down colors (MCD is also a diet model)
p_deg <- ggplot(deg_counts, aes(x = dataset, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(
    values = c("Up" = "#C23B75", "Down" = "#F2A45E"),
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
    title = "DEG Counts: Diet and MCD Mouse Models",
    subtitle = sprintf("|LFC| > %.1f, padj < %.2f", LFC_CUTOFF, P_CUTOFF),
    x = NULL, y = "Number of DEGs"
  )

cairo_pdf(file.path(OUTPUT_DIR, "deg_counts_comparison.pdf"), width = 5, height = 4, family = "Helvetica")
print(p_deg)
dev.off()
cat("  Created: deg_counts_comparison.pdf\n")

# Save summary
deg_summary <- deg_counts %>%
  pivot_wider(names_from = direction, values_from = count) %>%
  mutate(total = Up + Down)
write.csv(deg_summary, file.path(OUTPUT_DIR, "deg_counts_summary.csv"), row.names = FALSE)

# ==============================================================================
# 2. JACCARD SIMILARITY HEATMAP
# ==============================================================================

cat("Calculating Jaccard similarity...\n")

get_deg_set <- function(df, direction = "up") {
  df <- df %>% filter(!is.na(padj) & !is.na(log2FoldChange))
  if (direction == "up") {
    df %>% filter(padj < P_CUTOFF & log2FoldChange > LFC_CUTOFF) %>% pull(gene_id_clean)
  } else {
    df %>% filter(padj < P_CUTOFF & log2FoldChange < -LFC_CUTOFF) %>% pull(gene_id_clean)
  }
}

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

row_annot <- data.frame(
  Type = ifelse(grepl("MCD", datasets), "MCD", "Diet"),
  row.names = datasets
)

cairo_pdf(file.path(OUTPUT_DIR, "jaccard_heatmap_upregulated.pdf"), width = 6, height = 5, family = "Helvetica")
pheatmap(jaccard_up, 
         main = "Jaccard Similarity: Upregulated DEGs",
         color = colorRampPalette(c("white", "#C23B75"))(50),
         annotation_row = row_annot,
         annotation_col = row_annot,
         fontsize = 7, fontsize_row = 6, fontsize_col = 6,
         display_numbers = TRUE, number_format = "%.2f",
         number_color = "black", fontsize_number = 5)
dev.off()
cat("  Created: jaccard_heatmap_upregulated.pdf\n")

cairo_pdf(file.path(OUTPUT_DIR, "jaccard_heatmap_downregulated.pdf"), width = 6, height = 5, family = "Helvetica")
pheatmap(jaccard_down, 
         main = "Jaccard Similarity: Downregulated DEGs",
         color = colorRampPalette(c("white", "#F2A45E"))(50),
         annotation_row = row_annot,
         annotation_col = row_annot,
         fontsize = 7, fontsize_row = 6, fontsize_col = 6,
         display_numbers = TRUE, number_format = "%.2f",
         number_color = "black", fontsize_number = 5)
dev.off()
cat("  Created: jaccard_heatmap_downregulated.pdf\n")

write.csv(jaccard_up, file.path(OUTPUT_DIR, "jaccard_matrix_upregulated.csv"))
write.csv(jaccard_down, file.path(OUTPUT_DIR, "jaccard_matrix_downregulated.csv"))

# ==============================================================================
# 3. LFC CORRELATION WITH MCD POOLED
# ==============================================================================

cat("Calculating LFC correlations with MCD Pooled datasets...\n")

for (mcd_name in names(mcd_data)) {
  mcd_df <- mcd_data[[mcd_name]]
  
  mcd_lfc <- mcd_df %>%
    dplyr::select(gene_id_clean, log2FoldChange) %>%
    rename(lfc_mcd = log2FoldChange) %>%
    filter(!is.na(lfc_mcd))
  
  for (diet_name in names(diet_data)) {
    df <- diet_data[[diet_name]]
    diet_lfc <- df %>%
      dplyr::select(gene_id_clean, log2FoldChange) %>%
      rename(lfc_diet = log2FoldChange) %>%
      filter(!is.na(lfc_diet))
    
    merged <- inner_join(mcd_lfc, diet_lfc, by = "gene_id_clean")
    
    if (nrow(merged) > 100) {
      merged <- merged %>% filter(abs(lfc_mcd) < 10 & abs(lfc_diet) < 10)
      
      cor_p <- cor(merged$lfc_mcd, merged$lfc_diet, method = "pearson")
      cor_s <- cor(merged$lfc_mcd, merged$lfc_diet, method = "spearman")
      
      p <- ggplot(merged, aes(x = lfc_mcd, y = lfc_diet)) +
        geom_point(alpha = 0.3, size = 0.5, color = DIET_COLORS[[diet_name]]) +
        geom_smooth(method = "lm", color = "#C23B75", se = FALSE, linewidth = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        theme_minimal(base_family = "Helvetica", base_size = 7) +
        labs(
          title = paste(diet_name, "vs", mcd_name),
          subtitle = sprintf("r = %.3f, ρ = %.3f (n = %d genes)", cor_p, cor_s, nrow(merged)),
          x = paste(mcd_name, "LFC"),
          y = paste(diet_name, "LFC")
        )
      
      out_file <- paste0("lfc_correlation_", diet_name, "_vs_", mcd_name, ".pdf")
      cairo_pdf(file.path(OUTPUT_DIR, out_file), width = 4, height = 4, family = "Helvetica")
      print(p)
      dev.off()
    }
  }
}
cat("  Created LFC correlation plots\n")

# ==============================================================================
# 4. UPSET PLOT
# ==============================================================================

cat("Generating UpSet plot...\n")

all_deg_sets <- lapply(all_data, function(df) {
  df %>%
    filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
    filter(padj < P_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF) %>%
    pull(gene_id_clean) %>% unique()
})
all_deg_sets <- all_deg_sets[sapply(all_deg_sets, length) > 0]

if (length(all_deg_sets) >= 2) {
  cairo_pdf(file.path(OUTPUT_DIR, "upset_plot.pdf"), width = 8, height = 5, family = "Helvetica")
  print(upset(fromList(all_deg_sets), 
              nsets = length(all_deg_sets),
              order.by = "freq", decreasing = TRUE,
              mainbar.y.label = "Intersection Size",
              sets.x.label = "DEGs per Dataset",
              text.scale = c(1.3, 1.3, 1.0, 1.0, 1.3, 1.0),
              point.size = 3, line.size = 1))
  dev.off()
  cat("  Created: upset_plot.pdf\n")
}

# ==============================================================================
# 5. PCA OF LFC PROFILES
# ==============================================================================

cat("Generating PCA of LFC profiles...\n")

all_lfc <- lapply(names(all_data), function(name) {
  df <- all_data[[name]]
  df %>%
    dplyr::select(gene_id_clean, log2FoldChange) %>%
    rename(!!name := log2FoldChange)
})

lfc_matrix <- Reduce(function(x, y) full_join(x, y, by = "gene_id_clean"), all_lfc)
rownames(lfc_matrix) <- lfc_matrix$gene_id_clean
lfc_matrix <- lfc_matrix[, -1]
lfc_matrix <- lfc_matrix[rowSums(!is.na(lfc_matrix)) >= ncol(lfc_matrix) * 0.7, ]
lfc_matrix[is.na(lfc_matrix)] <- 0

vars <- apply(lfc_matrix, 1, var)
lfc_matrix <- lfc_matrix[order(vars, decreasing = TRUE)[1:min(5000, nrow(lfc_matrix))], ]

lfc_t <- t(as.matrix(lfc_matrix))
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
  )

cairo_pdf(file.path(OUTPUT_DIR, "pca_lfc_profiles.pdf"), width = 5, height = 4, family = "Helvetica")
print(p_pca)
dev.off()
cat("  Created: pca_lfc_profiles.pdf\n")

cat("\n========================================\n")
cat("Analysis Complete!\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("========================================\n")
