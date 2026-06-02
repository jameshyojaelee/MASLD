#!/usr/bin/env Rscript
# 23_pseudotime_trajectory.R
# ---------------------------------------------------------------------------
# MASLD Disease Progression Trajectory via Slingshot
# Objective: Model continuous gene activation waves from Healthy -> NASH
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(slingshot)
  library(tradeSeq)
  library(ggplot2)
})

cat("=== Phase 10: Pseudotime Trajectory (Slingshot) ===\n\n")

WD  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
RES <- file.path(WD, "Human/Patient_Cohorts/analysis/integration/results")
OUT <- file.path(RES, "pseudotime_trajectory")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Load Data
# ============================================================
cat("Loading PCA Coordinates and Metadata...\n")
# pca_res <- readRDS(file.path(RES, "harmony_pca_coordinates.rds"))
# meta <- readRDS(file.path(RES, "meta_matched.rds"))
# Ensure proper ordering
# meta <- meta[match(rownames(pca_res), rownames(meta)), ]

# Filter to disease progression relevant sample types (e.g. Healthy -> NAFL -> NASH)
# Also consider focusing on the BayesPrism deconvolved matrix later
# subset_idx <- meta$condition %in% c("Healthy", "NAFL", "NASH")
# pca_sub <- pca_res[subset_idx, ]
# meta_sub <- meta[subset_idx, ]

# ============================================================
# 2. Fit Slingshot Trajectory
# ============================================================
cat("Fitting Slingshot Lineages...\n")
# sds <- slingshot(pca_sub, clusterLabels = meta_sub$condition, start.clus = "Healthy")
#
# Pseudotime extraction
# pt <- slingPseudotime(sds)
#
# Save object
# saveRDS(sds, file.path(OUT, "slingshot_disease_trajectory.rds"))

# ============================================================
# 3. Model Dynamic Gene Expression (tradeSeq)
# ============================================================
cat("Fitting continuous expression models (tradeSeq)...\n")
# Note: tradeSeq requires count data.
# load counts matrix for subset
# counts_sub <- ...
# sce <- fitGAM(counts = counts_sub, pseudotime = pt, cellWeights = slingCurveWeights(sds))
#
# Association Test (finding genes that change significantly over pseudotime)
# assocRes <- associationTest(sce)
#
# Save results
# saveRDS(assocRes, file.path(OUT, "tradeseq_association_results.rds"))

cat("\nPipeline configured. Awaiting live execution post-GSE213621 integration.\n")
