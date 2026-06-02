#!/usr/bin/env Rscript
# M06_mouse_batch_correction_umap.R
# ---------------------------------------------------------------------------
# ComBat batch correction → PCA + UMAP visualization of all 9 mouse datasets.
# Input:  results/merged_dge.rds, results/meta_matched.rds
# Output: results/pca_corrected.pdf, umap_dataset.pdf, umap_diet.pdf,
#         umap_condition.pdf, corrected_logcpm.rds, umap_coordinates.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(sva)
  library(uwot)
  library(ggplot2)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration"
RDIR <- file.path(BASE, "results")

# --- Load merged DGE and metadata ---
dge  <- readRDS(file.path(RDIR, "merged_dge.rds"))
meta <- readRDS(file.path(RDIR, "meta_matched.rds"))
cat("Loaded DGE:", ncol(dge), "samples x", nrow(dge), "genes\n")

# Align metadata to DGE column order
meta <- meta[match(colnames(dge), meta$sample_id)]
stopifnot(all(meta$sample_id == colnames(dge)))

# Log-CPM
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

# --- ComBat batch correction ---
cat("\nRunning ComBat batch correction (batch = dataset)...\n")
batch <- meta$dataset
mod   <- model.matrix(~ group_binary, data = meta)

corrected <- ComBat(dat = logcpm, batch = batch, mod = mod)
cat("ComBat complete.\n")

# --- PCA on corrected expression ---
cat("\nComputing PCA...\n")
vars <- apply(corrected, 1, var)
top_genes <- names(sort(vars, decreasing = TRUE))[1:min(2000, length(vars))]
pca <- prcomp(t(corrected[top_genes, ]), scale. = TRUE, center = TRUE)
pve <- summary(pca)$importance[2, 1:2] * 100

pca_df <- data.table(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  dataset = meta$dataset,
  diet = meta$diet_model,
  group = meta$group_binary
)

# Dynamic color palettes — scale to however many datasets/diets are present
# Uses a colorRampPalette based on the publication magenta/blue family
.make_named_palette <- function(labels, seed_colors = c("#C2185B", "#AD1457",
    "#880E4F", "#7B1FA2", "#6A1B9A", "#1A237E", "#1565C0", "#42A5F5", "#26C6DA")) {
  labels <- sort(unique(labels))
  n <- length(labels)
  if (n <= length(seed_colors)) {
    pal <- seed_colors[seq_len(n)]
  } else {
    pal <- colorRampPalette(seed_colors)(n)
  }
  setNames(pal, labels)
}

datasets_present <- sort(unique(meta$dataset))
diets_present    <- sort(unique(meta$diet_model))

dataset_colors <- .make_named_palette(datasets_present)
diet_colors    <- .make_named_palette(diets_present,
  seed_colors = c("#0D47A1", "#C2185B", "#7B1FA2", "#E91E63", "#42A5F5",
                  "#00BCD4", "#FF9800", "#4CAF50", "#795548"))

theme_pub <- theme_minimal(base_size = 12) +
  theme(text = element_text(family = "sans"))

# PCA plot — color by dataset
pdf(file.path(RDIR, "pca_corrected.pdf"), width = 10, height = 7)
p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = dataset, shape = group)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = dataset_colors) +
  scale_shape_manual(values = c("Control" = 16, "Disease" = 17)) +
  labs(
    title = sprintf("PCA after ComBat Batch Correction (%d Samples, %d Datasets)",
                    ncol(dge), length(datasets_present)),
    x = sprintf("PC1 (%.1f%%)", pve[1]),
    y = sprintf("PC2 (%.1f%%)", pve[2])
  ) +
  theme_pub
print(p1)

# PCA — color by condition
p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Control" = "#1565C0", "Disease" = "#C2185B")) +
  labs(
    title = "PCA - Disease vs Control (after ComBat)",
    x = sprintf("PC1 (%.1f%%)", pve[1]),
    y = sprintf("PC2 (%.1f%%)", pve[2])
  ) +
  theme_pub
print(p2)
dev.off()
cat("Saved: pca_corrected.pdf\n")

# --- UMAP ---
cat("\nComputing UMAP...\n")
set.seed(42)
umap_out <- umap(t(corrected[top_genes, ]),
                  n_neighbors = 30, min_dist = 0.3,
                  n_components = 2, metric = "cosine")

umap_df <- data.table(
  UMAP1 = umap_out[, 1],
  UMAP2 = umap_out[, 2],
  dataset = meta$dataset,
  diet = meta$diet_model,
  group = meta$group_binary
)

# UMAP by dataset
pdf(file.path(RDIR, "umap_dataset.pdf"), width = 9, height = 7)
p3 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = dataset)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = dataset_colors) +
  labs(
    title = "UMAP - Mouse MASLD Datasets (ComBat-corrected)",
    subtitle = sprintf("%d samples, %d genes (top 2000 variable)", ncol(dge), length(top_genes))
  ) +
  theme_pub
print(p3)
dev.off()
cat("Saved: umap_dataset.pdf\n")

# UMAP by diet
pdf(file.path(RDIR, "umap_diet.pdf"), width = 9, height = 7)
p4 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = diet)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = diet_colors) +
  labs(
    title = "UMAP - Mouse MASLD by Diet Model (ComBat-corrected)",
    subtitle = sprintf("%d samples, %d diet models", ncol(dge), length(diets_present))
  ) +
  theme_pub
print(p4)
dev.off()
cat("Saved: umap_diet.pdf\n")

# UMAP by condition
pdf(file.path(RDIR, "umap_condition.pdf"), width = 9, height = 7)
p5 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = group)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Control" = "#1565C0", "Disease" = "#C2185B")) +
  labs(
    title = "UMAP - Disease vs Control (ComBat-corrected)",
    subtitle = sprintf("%d samples, %d diet models", ncol(dge), length(diets_present))
  ) +
  theme_pub
print(p5)
dev.off()
cat("Saved: umap_condition.pdf\n")

# Save corrected expression and UMAP coordinates
saveRDS(corrected, file.path(RDIR, "corrected_logcpm.rds"))
fwrite(umap_df, file.path(RDIR, "umap_coordinates.csv"))

cat("\nBatch correction and UMAP complete.\n")
