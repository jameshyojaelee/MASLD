#!/usr/bin/env Rscript
# 24_tf_grn_analysis.R
# ---------------------------------------------------------------------------
# Inferring Transcription Factor Activity using decoupleR / Viper
# Objective: Identify upstream Master Regulators driving consensus DEGs
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(decoupleR)
  library(progeny)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

cat("=== Phase 11: Transcription Factor Activity (decoupleR + viper) ===\n\n")

WD  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
RES <- file.path(WD, "Human/Patient_Cohorts/analysis/integration/results")
OUT <- file.path(RES, "tf_activity")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Load Data
# ============================================================
cat("Loading Disease Signatures and Regulons...\n")
# meta_stats <- fread(file.path(RES, "disease_signatures/unified_disease_signatures.csv"))

# Build input matrix: rows = genes, cols = contrast mapping (e.g., LFC)
# We use the meta-analysis LFC as our signature
# lfc_vec <- setNames(meta_stats$meta_logFC, meta_stats$symbol)
# lfc_mat <- as.matrix(lfc_vec)
# colnames(lfc_mat) <- "Meta_NASH_vs_Healthy"

# ============================================================
# 2. Extract DoRothEA or CollecTRI network
# ============================================================
cat("Querying CollecTRI Prior Knowledge Network...\n")
# net <- get_collectri(organism='human', split_complexes=FALSE)

# ============================================================
# 3. Model Transcription Factor Activity
# ============================================================
cat("Running univariate linear model (ulm) vs VIPER on LFC...\n")
#
# Run decoupleR methods
# acts <- run_wmean(mat=lfc_mat, net=net, .source="source", .target="target", .mor="weight", times=100)
# acts_ulm <- run_ulm(mat=lfc_mat, net=net, .source="source", .target="target", .mor="weight")
# acts_viper <- run_viper(mat=lfc_mat, net=net, .source="source", .target="target", .mor="weight")

# ============================================================
# 4. Integrate and Save
# ============================================================
cat("Finding Master Regulators...\n")

# Identify TFs whose activity is highly upregulated (driving NASH)
# or highly downregulated
# top_tfs <- head(acts_ulm[order(acts_ulm$score, decreasing = TRUE), ])

# saveRDS(acts_ulm, file=file.path(OUT, "TF_activity_ulm.rds"))
# fwrite(as.data.table(acts_ulm), file=file.path(OUT, "TF_activity_ulm.csv"))

# Run PROGENy for Pathway Level Summaries
# progeny_net <- get_progeny(organism = 'human', top = 100)
# progeny_acts <- run_wmean(mat=lfc_mat, net=progeny_net, .source="source", .target="target", .mor="weight")

cat("\nPipeline configured. Awaiting live execution post-GSE213621 integration.\n")
