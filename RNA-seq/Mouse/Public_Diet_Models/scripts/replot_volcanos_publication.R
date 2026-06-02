#!/usr/bin/env Rscript
# Publication-Quality Volcano Plots for Diet Mouse RNA-seq
# Uses gene symbols, Helvetica fonts, Sanjana Lab color scheme

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(stringr)
  library(AnnotationDbi)
  library(org.Mm.eg.db)  # Mouse annotations
})

# ==============================================================================
# PUBLICATION THEME CONFIGURATION (from publication_theme_guidelines.md)
# ==============================================================================

# Font Configuration for Illustrator (Type 42)
pdf.options(useDingbats = FALSE)

# Color Scheme (MCD Style)
COLOR_UP <- "#C23B75"     # Magenta for upregulated
COLOR_DOWN <- "#F2A45E"   # Orange for downregulated  
COLOR_SIG <- "#E07BB6"    # Pink for significant but |LFC| < 1
COLOR_NS <- "gray70"      # Gray for not significant

# Thresholds
P_CUTOFF <- 0.05
LFC_CUTOFF <- 1.0
TOP_N_LABELS <- 3

# Publication Theme
theme_publication <- function() {
  theme_minimal(base_family = "Helvetica", base_size = 7) +
    theme(
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 6),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7),
      legend.position = "bottom",
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ==============================================================================
# CONFIGURATION
# ==============================================================================

ROOT_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq"
DESEQ2_DIR <- file.path(ROOT_DIR, "analysis", "deseq2")

# Dataset labels for nice plot titles
DATASET_LABELS <- list(
  "GSE162876_20wks_FPC_vs_LFD" = "GSE162876: 20wk FPC vs LFD",
  "GSE162876_7wks_CDAHFD_vs_LFD" = "GSE162876: 7wk CDAHFD vs LFD",
  "GSE274914_52w_HFD_vs_LFD" = "GSE274914: 52w HFD vs LFD (Female)",
  "GSE274914_7w_HFD_vs_LFD" = "GSE274914: 7w HFD vs LFD",
  "GSE159911_LIDPAD_vs_Control" = "GSE159911: LIDPAD vs Control",
  "GSE225616_GAN_siRNA_vs_Vehicle" = "GSE225616: GAN siRNA vs Vehicle",
  "GSE263273_OCA_vs_Vehicle" = "GSE263273: OCA vs Vehicle",
  "GSE263273_INT787_vs_Vehicle" = "GSE263273: INT787 vs Vehicle"
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Map Ensembl IDs to Gene Symbols
map_symbols <- function(df) {
  # Identify gene_id column
  if (!"gene_id" %in% colnames(df)) {
    df$gene_id <- rownames(df)
  }
  
  # Clean Ensembl IDs (remove version suffix)
  ids <- df$gene_id
  ids_clean <- sub("\\..*", "", ids)
  
  cat("  Mapping", length(ids), "Ensembl IDs to gene symbols...\n")
  
  tryCatch({
    syms <- mapIds(org.Mm.eg.db, 
                   keys = ids_clean, 
                   column = "SYMBOL", 
                   keytype = "ENSEMBL", 
                   multiVals = "first")
    df$gene_symbol <- syms
    # Fill failed mappings with original ID (stripped version)
    df$gene_symbol[is.na(df$gene_symbol)] <- ids_clean[is.na(df$gene_symbol)]
    
    n_mapped <- sum(!is.na(syms))
    cat("  Successfully mapped", n_mapped, "of", length(ids), "genes\n")
  }, error = function(e) {
    cat("  Warning: Annotation mapping failed:", e$message, "\n")
    df$gene_symbol <- ids_clean
  })
  
  return(df)
}

# Generate Publication-Quality Volcano Plot
generate_volcano <- function(results_file, output_prefix, title) {
  cat("Processing:", results_file, "\n")
  
  # Load data
  df <- read.csv(results_file, stringsAsFactors = FALSE)
  
  # Required columns
  if (!all(c("log2FoldChange", "padj") %in% colnames(df))) {
    cat("  Skipping: Required columns missing\n")
    return(NULL)
  }
  
  # Map gene symbols
  df <- map_symbols(df)
  
  # Prepare for plotting
  df <- df %>%
    mutate(
      padj = ifelse(is.na(padj), 1, padj),
      neglog10 = -log10(padj),
      log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
      significance = case_when(
        padj <= P_CUTOFF & log2FoldChange >= LFC_CUTOFF ~ "Up",
        padj <= P_CUTOFF & log2FoldChange <= -LFC_CUTOFF ~ "Down",
        padj <= P_CUTOFF ~ "Significant",
        TRUE ~ "Not significant"
      )
    )
  
  # Cap extreme -log10 values for visualization
  max_neglog <- quantile(df$neglog10[is.finite(df$neglog10)], 0.995, na.rm = TRUE)
  df$neglog10 <- pmin(df$neglog10, max_neglog * 1.2)
  
  # Select top genes for labeling
  top_labels <- df %>%
    filter(significance %in% c("Up", "Down")) %>%
    arrange(padj) %>%
    slice_head(n = TOP_N_LABELS)
  
  # Count DEGs
  n_up <- sum(df$significance == "Up", na.rm = TRUE)
  n_down <- sum(df$significance == "Down", na.rm = TRUE)
  subtitle <- sprintf("Up: %d | Down: %d (|LFC| > %.1f, padj < %.2f)", 
                      n_up, n_down, LFC_CUTOFF, P_CUTOFF)
  
  # Create plot
  p <- ggplot(df, aes(x = log2FoldChange, y = neglog10)) +
    geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
    
    # Highlight top genes
    geom_point(data = top_labels, 
               aes(x = log2FoldChange, y = neglog10),
               color = "black", size = 2.5, shape = 21, fill = "#FFE066", stroke = 0.8) +
    
    # Gene labels
    geom_text_repel(
      data = top_labels, 
      aes(label = gene_symbol),
      size = 2.5, 
      fontface = "bold",
      box.padding = 0.4,
      point.padding = 0.3,
      segment.size = 0.3,
      segment.color = "gray40",
      max.overlaps = 25,
      force = 1.5
    ) +
    
    # Threshold lines
    geom_hline(yintercept = -log10(P_CUTOFF), linetype = "dashed", 
               color = "gray50", linewidth = 0.4) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", 
               color = "gray50", linewidth = 0.4) +
    
    # Colors
    scale_color_manual(
      values = c(
        "Up" = COLOR_UP,
        "Down" = COLOR_DOWN,
        "Significant" = COLOR_SIG,
        "Not significant" = COLOR_NS
      ),
      breaks = c("Up", "Down", "Significant", "Not significant")
    ) +
    
    # Theme
    theme_publication() +
    
    labs(
      title = title,
      subtitle = subtitle,
      x = expression(log[2]~"Fold Change"),
      y = expression(-log[10]~"(adjusted p-value)"),
      color = "Significance"
    ) +
    
    guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1)))
  
  # Save
  out_file <- paste0(output_prefix, "_volcano.pdf")
  cairo_pdf(out_file, width = 6, height = 5, family = "Helvetica")
  print(p)
  dev.off()
  
  cat("  Created:", out_file, "\n")
  
  return(list(n_up = n_up, n_down = n_down))
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

cat("========================================\n")
cat("Publication Volcano Plot Generator\n")
cat("========================================\n\n")

# Find all results files
results_files <- list.files(DESEQ2_DIR, pattern = "_results\\.csv$", full.names = TRUE)

cat("Found", length(results_files), "results files\n\n")

# Process each
summary_stats <- data.frame(
  dataset = character(),
  n_up = integer(),
  n_down = integer(),
  stringsAsFactors = FALSE
)

for (f in results_files) {
  # Extract dataset name
  basename <- sub("_results\\.csv$", "", basename(f))
  
  # Get nice title
  title <- DATASET_LABELS[[basename]]
  if (is.null(title)) title <- basename
  
  output_prefix <- file.path(DESEQ2_DIR, basename)
  
  stats <- generate_volcano(f, output_prefix, title)
  
  if (!is.null(stats)) {
    summary_stats <- rbind(summary_stats, data.frame(
      dataset = basename,
      n_up = stats$n_up,
      n_down = stats$n_down
    ))
  }
  
  cat("\n")
}

# Save summary
summary_file <- file.path(DESEQ2_DIR, "deg_summary.csv")
write.csv(summary_stats, summary_file, row.names = FALSE)
cat("Saved DEG summary to:", summary_file, "\n")

cat("\n========================================\n")
cat("Done!\n")
cat("========================================\n")
