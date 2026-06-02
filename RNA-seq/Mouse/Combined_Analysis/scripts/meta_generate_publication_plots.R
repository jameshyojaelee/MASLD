#!/usr/bin/env Rscript
# meta_generate_publication_plots.R
# Generate publication-quality plots for meta-analysis report
# Following Sanjana Lab publication guidelines

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
})

# Configure PDF device for Illustrator compatibility (Type 42 fonts)
pdf.options(useDingbats = FALSE)

# Source color palettes
color_file <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/publication_color_themes.R"
if (file.exists(color_file)) source(color_file)

# Define colors from lab palette
magenta <- "#e14b9d"
pink <- "#e35070"
purple <- "#d358c7"
blue <- "#4baeef"
green <- "#30d796"
orange <- "#e1b172"
gray <- "#656D77"

# Dataset-specific colors
dataset_colors <- c(
  "Combined" = magenta,
  "InHouse" = blue,
  "GSE156918" = green,
  "GSE205974" = orange
)

# Publication theme (following guidelines: Helvetica, 7-8pt fonts)
theme_publication <- theme_minimal(base_family = "Helvetica", base_size = 7) +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 7, color = "gray40"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(10, 10, 10, 10)
  )

root <- normalizePath(".")
out_dir <- file.path(root, "meta_analysis")
plot_dir <- file.path(root, "plots", "publication")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

message("=== Generating Publication-Quality Plots ===\n")

# Load data
gene_sets <- readRDS(file.path(out_dir, "upregulated_gene_sets.rds"))
all_degs <- readRDS(file.path(out_dir, "all_degs_list.rds"))

# ============================================================
# PLOT 1: DEG Counts Comparison Bar Chart
# ============================================================
message("Creating DEG counts comparison...")

deg_counts <- read_csv(file.path(out_dir, "deg_counts_summary.csv"), col_types = cols())

# Reshape for plotting
deg_long <- deg_counts %>%
  pivot_longer(cols = c(lenient, stringent), names_to = "threshold", values_to = "count") %>%
  mutate(
    threshold = factor(threshold, levels = c("lenient", "stringent"),
                       labels = c("padj<0.1, LFC>0", "padj<0.1, LFC>=0.58")),
    dataset = factor(dataset, levels = c("Combined", "InHouse", "GSE156918", "GSE205974"))
  )

p1 <- ggplot(deg_long, aes(x = dataset, y = count, fill = dataset)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = format(count, big.mark = ",")), 
            vjust = -0.3, size = 2.5, family = "Helvetica") +
  facet_wrap(~threshold, scales = "free_y") +
  scale_fill_manual(values = dataset_colors, guide = "none") +
  labs(
    title = "Upregulated DEGs by Dataset",
    subtitle = "MCD vs Control (padj < 0.1)",
    x = NULL, y = "Number of DEGs"
  ) +
  theme_publication +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_dir, "deg_counts_comparison.pdf"), p1, 
       width = 5, height = 3.5, device = cairo_pdf)
message("  Saved: deg_counts_comparison.pdf")

# ============================================================
# PLOT 2: Jaccard Similarity Heatmap
# ============================================================
message("Creating Jaccard similarity heatmap...")

jaccard_mat <- read.csv(file.path(out_dir, "jaccard_matrix_lenient.csv"), row.names = 1)

# Convert to long format
jaccard_long <- jaccard_mat %>%
  rownames_to_column("dataset1") %>%
  pivot_longer(-dataset1, names_to = "dataset2", values_to = "jaccard") %>%
  mutate(
    dataset1 = factor(dataset1, levels = c("Combined", "InHouse", "GSE156918", "GSE205974")),
    dataset2 = factor(dataset2, levels = c("Combined", "InHouse", "GSE156918", "GSE205974"))
  )

p2 <- ggplot(jaccard_long, aes(x = dataset1, y = dataset2, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", jaccard)), size = 3, family = "Helvetica") +
  scale_fill_gradient2(low = "white", mid = pink, high = magenta, 
                       midpoint = 0.5, limits = c(0, 1),
                       name = "Jaccard Index") +
  labs(
    title = "Overlap Similarity Between Datasets",
    subtitle = "Jaccard Index of upregulated DEGs (lenient threshold)",
    x = NULL, y = NULL
  ) +
  theme_publication +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  coord_fixed()

ggsave(file.path(plot_dir, "jaccard_heatmap.pdf"), p2, 
       width = 4.5, height = 4, device = cairo_pdf)
message("  Saved: jaccard_heatmap.pdf")

# ============================================================
# PLOT 3: Library Coverage Bar Chart
# ============================================================
message("Creating library coverage chart...")

library_cov <- read_csv(file.path(out_dir, "library_coverage.csv"), col_types = cols()) %>%
  mutate(dataset = factor(dataset, levels = c("Combined", "InHouse", "GSE156918", "GSE205974")))

p3 <- ggplot(library_cov, aes(x = dataset, y = lenient_pct, fill = dataset)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = paste0(lenient_pct, "%")), vjust = -0.3, size = 2.5, family = "Helvetica") +
  scale_fill_manual(values = dataset_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 75), expand = c(0, 0)) +
  labs(
    title = "Library Gene Detection by Dataset",
    subtitle = "% of final_core_degs.csv genes detected as upregulated",
    x = NULL, y = "Library Coverage (%)"
  ) +
  theme_publication +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_dir, "library_coverage.pdf"), p3, 
       width = 4, height = 3.5, device = cairo_pdf)
message("  Saved: library_coverage.pdf")

# ============================================================
# PLOT 4: LFC Correlation Scatter (Combined vs InHouse)
# ============================================================
message("Creating LFC correlation plots...")

combined <- all_degs$Combined %>%
  select(gene_id, combined_lfc = log2FoldChange, combined_padj = padj, gene_name)

inhouse <- all_degs$InHouse %>%
  select(gene_id, inhouse_lfc = log2FoldChange)

merged <- combined %>% 
  inner_join(inhouse, by = "gene_id") %>%
  filter(!is.na(combined_lfc), !is.na(inhouse_lfc))

cor_val <- cor(merged$inhouse_lfc, merged$combined_lfc, method = "spearman")

# Identify genes with large changes
merged <- merged %>%
  mutate(
    sig_up = combined_padj < 0.1 & combined_lfc > 0,
    highlight = abs(inhouse_lfc - combined_lfc) > 2 & sig_up
  )

top_genes <- merged %>% filter(highlight) %>% arrange(desc(abs(combined_lfc))) %>% head(8)

p4 <- ggplot(merged, aes(x = inhouse_lfc, y = combined_lfc)) +
  geom_point(aes(color = sig_up), alpha = 0.3, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_text_repel(data = top_genes, aes(label = gene_name), 
                  size = 2, family = "Helvetica", max.overlaps = 10,
                  box.padding = 0.3, segment.size = 0.2) +
  scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = magenta),
                     labels = c("Not sig", "Sig up in Combined"), name = NULL) +
  labs(
    title = "LFC Correlation: InHouse vs Combined",
    subtitle = sprintf("Spearman rho = %.3f", cor_val),
    x = "log2FC (InHouse)", y = "log2FC (Combined)"
  ) +
  theme_publication +
  theme(legend.position = c(0.15, 0.85)) +
  coord_fixed(ratio = 1, xlim = c(-8, 12), ylim = c(-8, 12))

ggsave(file.path(plot_dir, "lfc_correlation_inhouse.pdf"), p4, 
       width = 4, height = 4, device = cairo_pdf)
message("  Saved: lfc_correlation_inhouse.pdf")

# ============================================================
# PLOT 5: Gained vs Lost Genes Summary
# ============================================================
message("Creating gained/lost genes summary...")

overlap_summary <- read_csv(file.path(out_dir, "overlap_summary.csv"), col_types = cols())

# Extract gained and lost
gained_lost <- overlap_summary %>%
  filter(str_detect(comparison, "Gained|Lost")) %>%
  mutate(
    type = ifelse(str_detect(comparison, "Gained"), "Gained", "Lost"),
    type = factor(type, levels = c("Gained", "Lost"))
  )

p5 <- ggplot(gained_lost, aes(x = type, y = lenient, fill = type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = format(lenient, big.mark = ",")), vjust = -0.3, size = 3, family = "Helvetica") +
  scale_fill_manual(values = c("Gained" = green, "Lost" = pink), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Genes Gained & Lost in Combined Analysis",
    subtitle = "Compared to union of individual datasets",
    x = NULL, y = "Number of Genes"
  ) +
  theme_publication

ggsave(file.path(plot_dir, "gained_lost_genes.pdf"), p5, 
       width = 3.5, height = 3.5, device = cairo_pdf)
message("  Saved: gained_lost_genes.pdf")

# ============================================================
# PLOT 6: Cross-Species Summary
# ============================================================
message("Creating cross-species summary...")

cross_species <- read_csv(file.path(out_dir, "cross_species_overlap.csv"), col_types = cols()) %>%
  mutate(
    dataset = str_replace(comparison, "Combined vs ", ""),
    dataset = factor(dataset, levels = c("GSE130970_NAS", "GSE135251_NAS", "GSE130970_Fibrosis"))
  )

p6 <- ggplot(cross_species, aes(x = dataset, y = overlap_lenient, fill = dataset)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = overlap_lenient), vjust = -0.3, size = 2.5, family = "Helvetica") +
  scale_fill_manual(values = c(
    "GSE130970_NAS" = purple,
    "GSE135251_NAS" = blue,
    "GSE130970_Fibrosis" = orange
  ), guide = "none") +
  labs(
    title = "Mouse-Human Concordant Genes",
    subtitle = "Combined MCD DEGs overlapping with human patient DEGs",
    x = NULL, y = "Concordant Genes"
  ) +
  theme_publication +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_dir, "cross_species_overlap.pdf"), p6, 
       width = 4, height = 3.5, device = cairo_pdf)
message("  Saved: cross_species_overlap.pdf")

# ============================================================
# PLOT 7: Combined Summary Figure (Multi-panel)
# ============================================================
message("Creating combined summary figure...")

# Use patchwork if available, otherwise skip
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  
  combined_fig <- (p1 | p3) / (p2 | p5) +
    plot_annotation(
      title = "Combined MCD RNA-seq Meta-Analysis Summary",
      theme = theme(plot.title = element_text(size = 10, face = "bold", family = "Helvetica"))
    )
  
  ggsave(file.path(plot_dir, "summary_figure.pdf"), combined_fig, 
         width = 8, height = 7, device = cairo_pdf)
  message("  Saved: summary_figure.pdf")
} else {
  message("  Skipping summary figure (patchwork not installed)")
}

message("\n=== All publication plots saved to: ", plot_dir, " ===")
