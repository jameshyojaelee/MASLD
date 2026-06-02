#!/usr/bin/env Rscript
# run_advanced_enrichment.R - Single-Sample Enrichment (GSVA) & TF Activity (DoRothEA)
#
# This script performs "single-sample" (per-dataset) enrichment analysis to reveal
# heterogeneity across the meta-analyzed cohorts.
#
# Methods:
# 1. GSVA (Gene Set Variation Analysis) using Hallmark gene sets
# 2. DoRothEA (via decoupleR) for Transcription Factor activity inference
#
# Usage: Rscript scripts/run_advanced_enrichment.R

# Critical: Exclude user library
.libPaths(.libPaths()[!grepl("jameslee", .libPaths())])

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(GSVA)
    library(decoupleR)
    library(OmnipathR)
    library(ComplexHeatmap)
    library(circlize)
    library(msigdbr)
})

# ============================================================================
# Configuration
# ============================================================================
BASE_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
STREAMLIT_DIR <- file.path(BASE_DIR, "streamlit_deg_explorer")
PATHWAY_DIR <- file.path(BASE_DIR, "Analysis/downstream_analysis/pathway_analysis")
OUTPUT_DIR <- file.path(PATHWAY_DIR, "results/advanced")
PLOT_DIR <- file.path(PATHWAY_DIR, "plots/advanced")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load publication theme coloring
source(file.path(PATHWAY_DIR, "scripts/publication_theme.R"))

cat("=== Advanced Enrichment Analysis (GSVA & DoRothEA) ===\n\n")

# ============================================================================
# 1. Prepare Expression Matrix (LFC per Dataset)
# ============================================================================
cat("[1/5] Loading and preparing LFC matrix...\n")

master_file <- file.path(STREAMLIT_DIR, "master_ortholog_matrix.csv.gz")
master_df <- fread(master_file) %>% as_tibble()

# Extract LFC columns
lfc_cols <- grep("_lfc$", names(master_df), value = TRUE)
cat(sprintf("  Found %d datasets\n", length(lfc_cols)))

# Create matrix: Rows = Symbol, Cols = Dataset
# Handle duplicates by taking max magnitude
mat_df <- master_df %>%
    select(Symbol, all_of(lfc_cols)) %>%
    filter(!is.na(Symbol) & Symbol != "") %>%
    group_by(Symbol) %>%
    # Taking the row with max absolute mean LFC across datasets
    mutate(mean_abs = rowMeans(abs(across(all_of(lfc_cols))), na.rm = TRUE)) %>%
    slice_max(mean_abs, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-mean_abs)

# Convert to matrix
mat <- as.matrix(mat_df %>% select(-Symbol))
rownames(mat) <- mat_df$Symbol

# Clean dataset names
colnames(mat) <- gsub("_lfc$", "", colnames(mat))
colnames(mat) <- gsub("human_", "H_", colnames(mat)) # Shorten prefix
colnames(mat) <- gsub("mouse_", "M_", colnames(mat))

# NOTE: We do NOT impute NA values. GSVA handles NAs natively.
# Setting NA to 0 would incorrectly treat "not measured" as "no change",
# biasing results toward the null for datasets with more missing data.
# GSVA will ignore genes with NA values when computing enrichment scores.
na_count <- sum(is.na(mat))
cat(sprintf("  NA values retained: %d (%.1f%% of matrix)\n", 
            na_count, 100 * na_count / length(mat)))

cat(sprintf("  Matrix dimensions: %d genes x %d datasets\n", nrow(mat), ncol(mat)))

# ============================================================================
# 2. GSVA Analysis (Hallmark)
# ============================================================================
cat("[2/5] Running GSVA (Hallmark)...\n")

# Get Hallmark sets
hallmark_gs <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_list <- split(hallmark_gs$gene_symbol, hallmark_gs$gs_name)

# Run GSVA
# method="ssgsea" or "gsva". "gsva" is generally better for differential comparisons
gsva_res <- gsva(mat, hallmark_list, method = "gsva", verbose = FALSE)

# Save result
saveRDS(gsva_res, file.path(OUTPUT_DIR, "gsva_hallmark_res.rds"))
write.csv(gsva_res, file.path(OUTPUT_DIR, "gsva_hallmark_scores.csv"))

cat("  GSVA complete. Generating heatmap...\n")

# Plot Heatmap
pdf(file.path(PLOT_DIR, "gsva_hallmark_heatmap.pdf"), width = 12, height = 10)

# Clean pathway names
rownames(gsva_res) <- gsub("HALLMARK_", "", rownames(gsva_res))
rownames(gsva_res) <- gsub("_", " ", rownames(gsva_res))

# Colors: Blue (Down) - White - Magenta (Up)
col_fun <- colorRamp2(c(-0.5, 0, 0.5), c(palette1[6], "white", MAIN_COLOR))

Heatmap(gsva_res, 
        name = "GSVA Score", 
        col = col_fun,
        row_title = "Hallmark Pathways",
        column_title = "Datasets (Human & Mouse)",
        column_title_gp = gpar(fontface = "bold"),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2")

dev.off()
cat("  Saved: gsva_hallmark_heatmap.pdf\n")

# ============================================================================
# 3. DoRothEA TF Activity (decoupleR)
# ============================================================================
cat("[3/5] Inferring TF Activity (DoRothEA)...\n")

# Load Network (Human)
net <- get_dorothea(organism = "human", levels = c("A", "B", "C"))

# Run decoupleR (Multivariate Linear Model - mlm - is fast and robust)
# We treat the LFC matrix as "expression" signatures
tf_acts <- run_mlm(mat, net, .source = "source", .target = "target", .mor = "mor", minsize = 5)

# Extract scores (source x condition)
tf_mat <- tf_acts %>%
    pivot_wider(id_cols = source, names_from = condition, values_from = score) %>%
    column_to_rownames("source") %>%
    as.matrix()

# Save result
saveRDS(tf_mat, file.path(OUTPUT_DIR, "dorothea_tf_activity.rds"))
write.csv(tf_mat, file.path(OUTPUT_DIR, "dorothea_tf_scores.csv"))

cat(sprintf("  Inferred activity for %d TFs. Generating heatmap...\n", nrow(tf_mat)))

# ============================================================================
# 4. Top TF Visualization
# ============================================================================
cat("[4/5] Plotting Top TF Activities...\n")

# Filter for most variable TFs
top_n <- 25
tf_var <- apply(tf_mat, 1, var)
top_tfs <- names(sort(tf_var, decreasing = TRUE))[1:top_n]
tf_subset <- tf_mat[top_tfs, ]

pdf(file.path(PLOT_DIR, "dorothea_top25_tf_heatmap.pdf"), width = 12, height = 8)

# Colors: Stronger scale for TF activity (t-values from mlm)
# Usually ranges -3 to 3 or more
limit <- max(abs(tf_subset)) * 0.8
col_fun_tf <- colorRamp2(c(-limit, 0, limit), c("navy", "white", "firebrick"))

Heatmap(tf_subset,
        name = "TF Activity",
        col = col_fun_tf,
        row_title = "Top Variable Transcription Factors",
        column_title = "Datasets",
        rect_gp = gpar(col = "white", lwd = 1), # Grid lines
        row_names_gp = gpar(fontsize = 10, fontface = "bold"),
        column_names_gp = gpar(fontsize = 8),
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "pearson")

dev.off()
cat("  Saved: dorothea_top25_tf_heatmap.pdf\n")

# ============================================================================
# 5. Summary
# ============================================================================
cat("[5/5] Advanced Analysis Complete.\n")
cat(sprintf("Results in: %s\n", OUTPUT_DIR))
cat(sprintf("Plots in:   %s\n", PLOT_DIR))
