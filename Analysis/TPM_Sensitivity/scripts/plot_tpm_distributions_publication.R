#!/usr/bin/env Rscript
#' Publication-quality TPM Distribution Plots
#'
#' Uses cat_palette from publication_color_themes.R and generates informative,
#' aesthetically pleasing distribution visualizations in PDF format.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggridges)
  library(patchwork)
  library(scales)
})

# =============================================================================
# Configuration
# =============================================================================

# Hardcoded ROOT for SLURM compatibility
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

# Source color palette
source(file.path(ROOT, "publication_color_themes.R"))

# Output directory
OUTPUT_DIR <- file.path(ROOT, "RNA-seq", "TPM_analysis", "filtered")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Dataset configuration
datasets <- list(
  list(
    label = "MCD Week pooled (combined)",
    display = "Cas13 mouse MCD",
    file = "mcd_week_pooled_combined.tsv.gz",
    species = "Mouse"
  ),
  list(
    label = "GSE156918 (external MCD)",
    display = "Mouse MCD (Paquette)",
    file = "other_mcd_gse156918.tsv.gz",
    species = "Mouse"
  ),
  list(
    label = "GSE205974 (external MCD)",
    display = "Mouse MCD (Yue)",
    file = "other_mcd_gse205974.tsv.gz",
    species = "Mouse"
  ),
  list(
    label = "GSE130970 NAS high",
    display = "Human (Hoang)",
    file = "gse130970_nas_high.csv.gz",
    species = "Human"
  ),
  list(
    label = "GSE135251 NAS high",
    display = "Human (Govaere)",
    file = "gse135251_nas_high.csv.gz",
    species = "Human"
  )
)

DATA_DIR <- file.path(ROOT, "streamlit_deg_explorer", "data")

# Theme settings
theme_publication <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = "sans"),
      plot.title = element_text(size = base_size + 4, face = "bold", hjust = 0.5, margin = margin(b = 15)),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, color = "grey40", margin = margin(b = 10)),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 2),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.margin = margin(20, 20, 20, 20),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

load_tpm_data <- function() {
  all_data <- list()
  
  for (ds in datasets) {
    path <- file.path(DATA_DIR, ds$file)
    if (!file.exists(path)) {
      message(sprintf("Skipping %s: file not found", ds$label))
      next
    }
    
    sep <- if (grepl("\\.tsv", ds$file)) "\t" else ","
    df <- read.csv(gzfile(path), sep = sep, stringsAsFactors = FALSE)
    
    if (!"tpm_mean" %in% colnames(df)) {
      message(sprintf("Skipping %s: no tpm_mean column", ds$label))
      next
    }
    
    # FILTERING: padj < 0.1, Log2FC > 0.8, and remove zero expression
    if ("padj" %in% colnames(df) && "log2FoldChange" %in% colnames(df)) {
        df <- df %>% filter(padj < 0.1, log2FoldChange >= 0.58, tpm_mean > 0)
    } else {
        message(sprintf("Warning: %s missing padj/log2FoldChange, filtering only for tpm > 0", ds$label))
        df <- df %>% filter(tpm_mean > 0)
    }

    tpm_vals <- df$tpm_mean[!is.na(df$tpm_mean)]
    
    all_data[[ds$label]] <- data.frame(
      dataset = ds$display,
      species = ds$species,
      tpm = tpm_vals,
      log2_tpm = log2(tpm_vals + 1),
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows(all_data)
}

compute_summary_stats <- function(data) {
  data %>%
    group_by(dataset, species) %>%
    summarize(
      n_genes = n(),
      median_tpm = median(tpm, na.rm = TRUE),
      mean_tpm = mean(tpm, na.rm = TRUE),
      q25_tpm = quantile(tpm, 0.25, na.rm = TRUE),
      q75_tpm = quantile(tpm, 0.75, na.rm = TRUE),
      pct_zero = mean(tpm == 0, na.rm = TRUE) * 100,
      pct_low = mean(tpm > 0 & tpm < 1, na.rm = TRUE) * 100,
      pct_med = mean(tpm >= 1 & tpm < 10, na.rm = TRUE) * 100,
      pct_high = mean(tpm >= 10, na.rm = TRUE) * 100,
      .groups = "drop"
    )
}

# =============================================================================
# Plot Functions
# =============================================================================

# Color setup using cat_palette
get_dataset_colors <- function(data) {
  unique_datasets <- unique(data$dataset)
  n <- length(unique_datasets)
  colors <- cat_palette[1:n]
  names(colors) <- unique_datasets
  colors
}

get_species_colors <- function() {
  c("Mouse" = cat_palette[2], "Human" = cat_palette[1])
}

# ... (skip to plot_ridgeline)

# 1. Ridgeline Plot with density peaks
plot_ridgeline <- function(data, colors) {
  # Order by species then by name
  data$dataset <- factor(data$dataset, levels = rev(unique(data$dataset)))
  
  # Cap at 5 for visualization
  data$log2_tpm <- pmin(data$log2_tpm, 5)

  ggplot(data, aes(x = log2_tpm, y = dataset, fill = dataset, color = dataset)) +
    geom_density_ridges(
      alpha = 0.85,
      scale = 2.5,
      rel_min_height = 0.01,
      quantile_lines = TRUE,
      quantiles = 2,
      color = "white",
      linewidth = 0.8
    ) +
    geom_vline(xintercept = 5, linetype = "dashed", color = "grey50") +
    annotate("text", x = 4.8, y = length(unique(data$dataset)) + 0.5, 
             label = "Capped at 5", hjust = 1, vjust = 0, size = 4, fontface = "italic", color = "grey40") +
    scale_fill_manual(values = colors) +
    scale_x_continuous(
      name = expression(bold(log[2](TPM + 1))),
      limits = c(0, 5.2),
      breaks = seq(0, 5, 1),
      expand = c(0, 0)
    ) +
    labs(
      title = "TPM Distribution (Filtered: padj < 0.1, LFC >= 0.58)",
      subtitle = "Ridgeline density plots (Capped at log2(TPM+1) = 5)",
      y = NULL
    ) +
    theme_publication() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 12, face = "bold"),
      panel.grid.major.y = element_blank()
    )
}

# 2. Violin + Boxplot Hybrid
plot_violin_box <- function(data, colors) {
  ggplot(data, aes(x = dataset, y = log2_tpm, fill = dataset)) +
    geom_violin(
      alpha = 0.7,
      color = "white",
      linewidth = 0.5,
      trim = FALSE,
      scale = "width"
    ) +
    geom_boxplot(
      width = 0.15,
      fill = "white",
      alpha = 0.9,
      color = "grey30",
      outlier.shape = NA
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 18,
      size = 3,
      color = "grey20"
    ) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(
      name = expression(bold(log[2](TPM + 1))),
      limits = c(0, 12),
      breaks = seq(0, 12, 2)
    ) +
    labs(
      title = "Expression Distribution by Dataset",
      subtitle = "Violin plots with embedded boxplots and median markers",
      x = NULL
    ) +
    theme_publication() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 25, hjust = 1, size = 11, face = "bold")
    )
}

# 3. Species Comparison (combined density)
plot_species_density <- function(data) {
  species_colors <- get_species_colors()
  
  ggplot(data, aes(x = log2_tpm, fill = species, color = species)) +
    geom_density(alpha = 0.5, linewidth = 1.2) +
    scale_fill_manual(values = species_colors, name = "Species") +
    scale_color_manual(values = species_colors, name = "Species") +
    scale_x_continuous(
      name = expression(bold(log[2](TPM + 1))),
      limits = c(0, 12),
      breaks = seq(0, 12, 2)
    ) +
    labs(
      title = "Expression Distribution by Species",
      subtitle = "Aggregated density comparison: Mouse MCD vs Human MASLD",
      y = "Density"
    ) +
    theme_publication() +
    theme(
      legend.position = c(0.85, 0.85),
      legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3)
    )
}

# 4. Cumulative Distribution (ECDF)
plot_ecdf <- function(data, colors) {
  ggplot(data, aes(x = log2_tpm, color = dataset)) +
    stat_ecdf(linewidth = 1.2, alpha = 0.9) +
    scale_color_manual(values = colors, name = "Dataset") +
    scale_x_continuous(
      name = expression(bold(log[2](TPM + 1))),
      limits = c(0, 12),
      breaks = seq(0, 12, 2)
    ) +
    scale_y_continuous(
      name = "Cumulative Proportion",
      labels = percent_format(),
      breaks = seq(0, 1, 0.2)
    ) +
    labs(
      title = "Cumulative Expression Distribution",
      subtitle = "Empirical cumulative distribution functions by dataset"
    ) +
    theme_publication() +
    theme(
      legend.position = c(0.85, 0.35),
      legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3)
    )
}

# 5. Expression Bins Stacked Bar
plot_expression_bins <- function(summary_stats, colors) {
  bins_df <- summary_stats %>%
    select(dataset, pct_zero, pct_low, pct_med, pct_high) %>%
    pivot_longer(
      cols = starts_with("pct_"),
      names_to = "bin",
      values_to = "percentage"
    ) %>%
    mutate(
      bin = factor(
        bin,
        levels = c("pct_zero", "pct_low", "pct_med", "pct_high"),
        labels = c("Zero (TPM = 0)", "Low (0 < TPM < 1)", "Medium (1 ≤ TPM < 10)", "High (TPM ≥ 10)")
      )
    )
  
  bin_colors <- c(
    "Zero (TPM = 0)" = cat_palette[18],
    "Low (0 < TPM < 1)" = cat_palette[5],
    "Medium (1 ≤ TPM < 10)" = cat_palette[1],
    "High (TPM ≥ 10)" = cat_palette[2]
  )
  
  ggplot(bins_df, aes(x = dataset, y = percentage, fill = bin)) +
    geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = bin_colors, name = "Expression Level") +
    scale_y_continuous(
      name = "Percentage of Genes",
      labels = percent_format(scale = 1),
      expand = c(0, 0)
    ) +
    labs(
      title = "Gene Expression Level Distribution",
      subtitle = "Proportion of genes in each expression category",
      x = NULL
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1, size = 11, face = "bold"),
      legend.position = "right"
    )
}

# 6. Heatmap-style Summary
plot_summary_heatmap <- function(summary_stats) {
  # Create a metric matrix
  metrics_df <- summary_stats %>%
    select(dataset, median_tpm, mean_tpm, pct_zero, pct_high) %>%
    pivot_longer(
      cols = -dataset,
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(
        metric,
        levels = c("median_tpm", "mean_tpm", "pct_zero", "pct_high"),
        labels = c("Median TPM", "Mean TPM", "% Zero Expression", "% High Expression")
      )
    )
  
  # Scale values within each metric for coloring
  metrics_df <- metrics_df %>%
    group_by(metric) %>%
    mutate(scaled = (value - min(value)) / (max(value) - min(value) + 1e-9)) %>%
    ungroup()
  
  ggplot(metrics_df, aes(x = metric, y = dataset, fill = scaled)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(
      aes(label = ifelse(grepl("%", metric), sprintf("%.1f%%", value), sprintf("%.2f", value))),
      size = 4,
      fontface = "bold",
      color = ifelse(metrics_df$scaled > 0.5, "white", "grey20")
    ) +
    scale_fill_gradientn(
      colors = c(cat_palette[1], cat_palette[3], cat_palette[2]),
      guide = "none"
    ) +
    labs(
      title = "Expression Summary Metrics",
      subtitle = "Key statistics across datasets (color-scaled within metric)",
      x = NULL,
      y = NULL
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = 11, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      panel.grid = element_blank()
    )
}

# 7. Multi-panel Summary Figure
create_summary_figure <- function(data, summary_stats, colors) {
  p1 <- plot_ridgeline(data, colors)
  p2 <- plot_violin_box(data, colors)
  p3 <- plot_species_density(data)
  p4 <- plot_expression_bins(summary_stats, colors)
  
  combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = "TPM Expression Distribution Analysis",
      subtitle = "Comprehensive overview of gene expression levels across MASLD/MCD datasets",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey40")
      )
    ) +
    plot_layout(heights = c(1, 1))
  
  combined
}

# =============================================================================
# Main Execution
# =============================================================================

main <- function() {
  message("Loading TPM data...")
  data <- load_tpm_data()
  
  if (nrow(data) == 0) {
    stop("No data loaded. Check file paths.")
  }
  
  message(sprintf("Loaded %d observations from %d datasets", nrow(data), length(unique(data$dataset))))
  
  # Compute summary statistics
  summary_stats <- compute_summary_stats(data)
  
  # Get colors
  colors <- get_dataset_colors(data)
  
  # Save summary stats
  stats_path <- file.path(OUTPUT_DIR, "tpm_distribution_summary_stats.tsv")
  write.table(summary_stats, stats_path, sep = "\t", row.names = FALSE, quote = FALSE)
  message(sprintf("Saved summary statistics to %s", stats_path))
  
  # Generate individual plots
  message("Generating plots...")
  
  # 1. Ridgeline
  pdf(file.path(OUTPUT_DIR, "tpm_distribution_ridgeline.pdf"), width = 10, height = 7)
  print(plot_ridgeline(data, colors))
  dev.off()
  message("  - Ridgeline plot saved")
  
  # 2. Violin + Box
  pdf(file.path(OUTPUT_DIR, "tpm_distribution_violin.pdf"), width = 10, height = 7)
  print(plot_violin_box(data, colors))
  dev.off()
  message("  - Violin plot saved")
  
  # 3. Species Density
  pdf(file.path(OUTPUT_DIR, "tpm_distribution_species_comparison.pdf"), width = 10, height = 6)
  print(plot_species_density(data))
  dev.off()
  message("  - Species comparison saved")
  
  # 4. ECDF
  pdf(file.path(OUTPUT_DIR, "tpm_distribution_ecdf.pdf"), width = 10, height = 6)
  print(plot_ecdf(data, colors))
  dev.off()
  message("  - ECDF plot saved")
  
  # 5. Expression Bins
  pdf(file.path(OUTPUT_DIR, "tpm_distribution_expression_bins.pdf"), width = 10, height = 7)
  print(plot_expression_bins(summary_stats, colors))
  dev.off()
  message("  - Expression bins plot saved")
  
  # 6. Summary Heatmap
  pdf(file.path(OUTPUT_DIR, "tpm_distribution_summary_heatmap.pdf"), width = 10, height = 6)
  print(plot_summary_heatmap(summary_stats))
  dev.off()
  message("  - Summary heatmap saved")
  
  # 7. Multi-panel Summary
  pdf(file.path(OUTPUT_DIR, "tpm_distribution_comprehensive_summary.pdf"), width = 16, height = 12)
  print(create_summary_figure(data, summary_stats, colors))
  dev.off()
  message("  - Comprehensive summary figure saved")
  
  message("\nAll plots saved to: ", OUTPUT_DIR)
  message("Done!")
}

# Run if called directly
if (sys.nframe() == 0) {
  main()
}
