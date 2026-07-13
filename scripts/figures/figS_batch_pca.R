#!/usr/bin/env Rscript
# =============================================================================
# figS_batch_pca.R
# Batch-effect PCA visualization: PC1 vs PC2 colored by cohort and disease status
# Output: figures/supplementary/figS_batch_pca.pdf (two-panel figure)
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

dge_path <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")
meta_path <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
out_path <- file.path(BASE, "figures/supplementary/figS_batch_pca.pdf")

# --- Publication theme (minimal) ---
theme_pub <- theme_bw(base_size = 6) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    strip.background = element_blank(),
    strip.text = element_text(face = "plain"),
    legend.position = "right",
    legend.background = element_blank()
  )

cat("=== Batch-effect PCA Visualization ===\n")
cat("Loading DGE object...\n")

dge <- readRDS(dge_path)
cat(sprintf("  DGE: %d genes x %d samples\n", nrow(dge), ncol(dge)))

cat("Loading metadata...\n")
meta <- fread(meta_path)

# --- Normalize ---
cat("Computing logCPM (TMM-normalized)...\n")
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

# --- Top 5000 most variable genes ---
cat("Selecting top 5000 most variable genes...\n")
gene_vars <- apply(logcpm, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(5000, length(gene_vars))]
logcpm_top <- logcpm[top_genes, ]
cat(sprintf("  Using %d genes for PCA\n", length(top_genes)))

# --- PCA ---
cat("Running PCA...\n")
pca <- prcomp(t(logcpm_top), center = TRUE, scale. = TRUE)
var_explained <- summary(pca)$importance[2, 1:5] * 100

pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)

# --- Merge metadata ---
# Identify the dataset/cohort column
if ("dataset" %in% names(dge$samples)) {
  pca_df$dataset <- dge$samples$dataset
} else if ("group" %in% names(dge$samples)) {
  # Try to get from unified metadata
  pca_df$dataset <- NA
}

if ("group_binary" %in% names(dge$samples)) {
  pca_df$group_binary <- dge$samples$group_binary
} else {
  pca_df$group_binary <- NA
}

# If dataset not in DGE samples, try from metadata
if (all(is.na(pca_df$dataset))) {
  cat("  Dataset not in DGE samples, merging from unified metadata...\n")
  meta_sub <- meta[, .(sample_id, dataset)]
  meta_sub <- unique(meta_sub)
  pca_df <- merge(pca_df, meta_sub, by = "sample_id", all.x = TRUE, suffixes = c("", ".meta"))
  if ("dataset.meta" %in% names(pca_df)) {
    pca_df$dataset <- pca_df$dataset.meta
    pca_df$dataset.meta <- NULL
  }
}

if (all(is.na(pca_df$group_binary))) {
  cat("  group_binary not in DGE samples, merging from unified metadata...\n")
  meta_sub2 <- meta[, .(sample_id, group_binary)]
  meta_sub2 <- unique(meta_sub2)
  pca_df <- merge(pca_df, meta_sub2, by = "sample_id", all.x = TRUE, suffixes = c("", ".meta2"))
  if ("group_binary.meta2" %in% names(pca_df)) {
    pca_df$group_binary <- pca_df$group_binary.meta2
    pca_df$group_binary.meta2 <- NULL
  }
}

cat(sprintf("  Datasets: %s\n", paste(sort(unique(pca_df$dataset)), collapse = ", ")))
cat(sprintf("  Groups: %s\n", paste(sort(unique(pca_df$group_binary)), collapse = ", ")))
cat(sprintf("  PC1: %.1f%% variance | PC2: %.1f%% variance\n", var_explained[1], var_explained[2]))

# --- Color palettes ---
n_datasets <- length(unique(na.omit(pca_df$dataset)))
# Qualitative palette for cohorts
cohort_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
  "#D55E00", "#CC79A7", "#999999", "#000000"
)[seq_len(n_datasets)]
names(cohort_colors) <- sort(unique(na.omit(pca_df$dataset)))

# Disease status: Control = gray, Disease = colored
group_colors <- c("Control" = "#9E9E9E", "Disease" = "#D32F2F")

# --- Panel A: PCA colored by cohort ---
p_cohort <- ggplot(pca_df, aes(x = PC1, y = PC2, color = dataset)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_manual(values = cohort_colors, name = "Cohort") +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_pub +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

# --- Panel B: PCA colored by disease status ---
p_group <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group_binary)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_manual(values = group_colors, name = "Status") +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_pub +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

# --- Combine and save ---
p_combined <- p_cohort + p_group + plot_layout(ncol = 2, widths = c(1, 1))

message("[caption] Panel A: PCA by cohort. Panel B: PCA by disease status.")

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
ggsave(out_path, p_combined, width = 7.09, height = 3.04, device = cairo_pdf)
cat(sprintf("\nSaved: %s\n", out_path))

# --- Also save PCA coordinates for reference ---
coord_path <- file.path(BASE, "figures/supplementary/figS_batch_pca_coords.csv")
fwrite(pca_df, coord_path)
cat(sprintf("Saved PCA coordinates: %s\n", coord_path))

cat("\n=== Done ===\n")
