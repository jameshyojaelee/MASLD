#!/usr/bin/env Rscript
# 09_batch_correction_umap.R
# ---------------------------------------------------------------------------
# Harmony batch correction → PCA + UMAP visualization of all 6 cohorts.
# Uses dataset as the batch variable for Harmony correction.
# Input:  results/integration/merged_dge.rds
# Output: results/integration/pca_corrected.pdf, umap_dataset.pdf, umap_condition.pdf
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(harmony)
  library(uwot)
  library(ggplot2)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/integration")

# --- Load merged DGE (all QC-passing samples, 6 cohorts) ---
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
cat("Loaded DGE:", ncol(dge), "samples x", nrow(dge), "genes\n")

# Log-CPM
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

# --- PCA (pre-correction) ---
cat("\nComputing PCA on top 2000 variable genes...\n")
vars <- apply(logcpm, 1, var)
top_genes <- names(sort(vars, decreasing = TRUE))[1:min(2000, length(vars))]
pca <- prcomp(t(logcpm[top_genes, ]), scale. = TRUE, center = TRUE)

# --- Harmony batch correction in PCA space ---
# Uses dataset as the batch variable.
# (PRJNA512027 — with its L0/S0 library-prep confound — was permanently removed
# from the pipeline 2026-05-15.)
cat("\nRunning Harmony batch correction...\n")
batch <- as.character(dge$samples$dataset)
cat("Batch levels:", paste(sort(unique(batch)), collapse = ", "), "\n")

harmony_out <- HarmonyMatrix(
  pca$x[, 1:30],            # top 30 PCs
  meta_data = data.frame(batch = batch),
  vars_use = "batch",
  do_pca = FALSE,
  max.iter.harmony = 20,
  verbose = TRUE
)
cat("Harmony complete.\n")

# --- PCA on corrected embeddings ---
pve <- summary(pca)$importance[2, 1:2] * 100

pca_df <- data.table(
  PC1 = harmony_out[, 1],
  PC2 = harmony_out[, 2],
  dataset = dge$samples$dataset,
  condition = dge$samples$condition,
  group = dge$samples$group_binary
)

# PCA plots
pdf(file.path(RDIR, "pca_corrected.pdf"), width = 10, height = 7)

# Color by dataset
p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = dataset, shape = group)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c(
    "GSE126848"   = "#AD1457",
    "GSE130970"   = "#7B1FA2",
    "GSE135251"   = "#1565C0",
    "GSE167523"   = "#00695C",
    "GSE213621"   = "#C2185B"
  )) +
  scale_shape_manual(values = c("Control" = 16, "Disease" = 17)) +
  labs(
    title = sprintf("PCA after Harmony Batch Correction (%d cohorts)", length(unique(batch))),
    subtitle = sprintf("%d samples, Harmony on top 30 PCs", ncol(dge)),
    x = "Harmony PC1", y = "Harmony PC2"
  ) +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = "sans"))
print(p1)

# Color by condition
p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, shape = dataset)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Control" = "#1565C0", "Disease" = "#C2185B")) +
  labs(
    title = "PCA - Disease vs Control (after Harmony)",
    x = "Harmony PC1", y = "Harmony PC2"
  ) +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = "sans"))
print(p2)
dev.off()
cat("Saved: pca_corrected.pdf\n")

# --- UMAP on Harmony-corrected PCs ---
cat("\nComputing UMAP on Harmony-corrected PCs...\n")
set.seed(42)
umap_out <- umap(harmony_out,
                 n_neighbors = 30, min_dist = 0.3,
                 n_components = 2, metric = "cosine")

umap_df <- data.table(
  UMAP1 = umap_out[, 1],
  UMAP2 = umap_out[, 2],
  dataset = dge$samples$dataset,
  condition = dge$samples$condition,
  group = dge$samples$group_binary
)

# UMAP by dataset
pdf(file.path(RDIR, "umap_dataset.pdf"), width = 9, height = 7)
p3 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = dataset)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c(
    "GSE126848"   = "#AD1457",
    "GSE130970"   = "#7B1FA2",
    "GSE135251"   = "#1565C0",
    "GSE162694"   = "#0277BD",
    "GSE167523"   = "#00695C",
    "GSE174478"   = "#2E7D32",
    "GSE193066"   = "#F57F17",
    "GSE213621"   = "#E64A19",
    "GSE240729"   = "#BF360C"
  )) +
  labs(
    title = sprintf("UMAP - Human MASLD Cohorts (Harmony-corrected, %d cohorts)", length(unique(batch))),
    subtitle = sprintf("%d samples, Harmony on top 30 PCs (2000 var genes)", ncol(dge))
  ) +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = "sans"))
print(p3)
dev.off()
cat("Saved: umap_dataset.pdf\n")

# UMAP by condition
pdf(file.path(RDIR, "umap_condition.pdf"), width = 9, height = 7)
p4 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = group)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Control" = "#1565C0", "Disease" = "#C2185B")) +
  labs(
    title = "UMAP - Disease vs Control (Harmony-corrected)",
    subtitle = sprintf("%d samples, Harmony on top 30 PCs (2000 var genes)", ncol(dge))
  ) +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = "sans"))
print(p4)
dev.off()
cat("Saved: umap_condition.pdf\n")

# Save UMAP coordinates for downstream use
fwrite(umap_df, file.path(RDIR, "umap_coordinates.csv"))
cat("\nHarmony batch correction and UMAP complete.\n")
