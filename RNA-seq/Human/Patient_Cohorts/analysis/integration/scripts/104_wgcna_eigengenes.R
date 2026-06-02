#!/usr/bin/env Rscript
# 104_wgcna_eigengenes.R
# ---------------------------------------------------------------------------
# Compute WGCNA module eigengenes from the top 5,000 most variable genes.
#
# Input:  results/integration/merged_dge.rds (34,453 genes x 1,444 samples)
# Output: results/staging_classifier/wgcna_eigengenes.rds       (samples x modules matrix)
#         results/staging_classifier/wgcna_module_assignments.csv (gene → module mapping)
#         results/staging_classifier/wgcna_soft_threshold.csv     (power selection table)
#
# Usage: Rscript 104_wgcna_eigengenes.R
# SLURM: cpu, 8 CPUs, 64G, 48h
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(WGCNA)
})

# WGCNA settings
allowWGCNAThreads(nThreads = 8)
options(stringsAsFactors = FALSE)
set.seed(42)

BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR   <- file.path(INT, "results/integration")
OUTDIR <- file.path(INT, "results/staging_classifier")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 104: WGCNA Module Eigengenes ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ============================================================
# STEP 1: Load expression data and compute TMM logCPM
# ============================================================
cat("--- Step 1: Load DGE and compute logCPM ---\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
cat("DGE loaded:", nrow(dge), "genes x", ncol(dge), "samples\n")

dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
cat("logCPM matrix:", nrow(logcpm), "x", ncol(logcpm), "\n\n")

# ============================================================
# STEP 2: Select top 5,000 most variable genes
# ============================================================
cat("--- Step 2: Select top 5,000 most variable genes ---\n")
gene_vars <- apply(logcpm, 1, var)
top_idx <- order(gene_vars, decreasing = TRUE)[1:5000]
logcpm_top <- logcpm[top_idx, ]
cat("Selected top 5,000 genes by variance\n")
cat("Variance range of selected genes:", round(range(gene_vars[top_idx]), 3), "\n")
cat("Variance range of excluded genes:", round(range(gene_vars[-top_idx]), 3), "\n\n")

# ============================================================
# STEP 3: Transpose for WGCNA (samples as rows, genes as columns)
# ============================================================
cat("--- Step 3: Transpose matrix ---\n")
datExpr <- t(logcpm_top)
cat("datExpr:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")

# Check for problematic genes/samples
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  cat("Removing problematic genes or samples...\n")
  if (sum(!gsg$goodGenes) > 0) {
    cat("  Removed genes:", sum(!gsg$goodGenes), "\n")
  }
  if (sum(!gsg$goodSamples) > 0) {
    cat("  Removed samples:", sum(!gsg$goodSamples), "\n")
  }
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
cat("After QC:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n\n")

# ============================================================
# STEP 4: Pick soft-thresholding power
# ============================================================
cat("--- Step 4: Pick soft-thresholding power ---\n")
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 3,
                          networkType = "signed")

sft_table <- as.data.table(sft$fitIndices)
setnames(sft_table, c("Power", "SFT.R.sq", "slope", "truncated.R.sq",
                        "mean.k.", "median.k.", "max.k."))
fwrite(sft_table, file.path(OUTDIR, "wgcna_soft_threshold.csv"))
cat("\nSoft threshold table:\n")
print(sft_table[, .(Power, SFT.R.sq = round(SFT.R.sq, 3),
                      slope = round(slope, 2), mean.k. = round(mean.k., 1))])

# Select power: MUST have negative slope (scale-free criterion) AND R^2 > 0.8
# The standard pickSoftThreshold often picks power=1 for multi-cohort data
# because R^2 is high but the SLOPE IS WRONG (positive = not scale-free)
chosen_power <- sft$powerEstimate

# Override if slope is positive (wrong scale-free direction) or power < 6
# For multi-cohort bulk RNA-seq, power 10-14 is typical for signed networks
neg_slope <- sft_table[slope < 0 & SFT.R.sq >= 0.25]
if (!is.na(chosen_power) && chosen_power < 6) {
  cat("\nWARNING: Auto-selected power", chosen_power, "is too low (slope likely positive).\n")
  cat("Multi-cohort data requires higher power to suppress batch-driven correlations.\n")
  chosen_power <- NA  # Force manual selection
}
if (is.na(chosen_power) || is.null(chosen_power)) {
  if (nrow(neg_slope) > 0) {
    # Pick first power with negative slope and R^2 >= 0.25
    chosen_power <- neg_slope$Power[1]
    cat("Selected power", chosen_power, "based on negative slope criterion\n")
  } else {
    # No negative slope found — use power=12 as a sensible default for signed networks
    chosen_power <- 12
    cat("No power achieved negative slope. Using default power=12 for signed network.\n")
    cat("(This is expected for multi-cohort data with batch effects.)\n")
  }
}
cat("\nFinal chosen power:", chosen_power, "\n")
cat("\n")

# ============================================================
# STEP 5: Build network and detect modules
# ============================================================
cat("--- Step 5: blockwiseModules ---\n")
cat("Power:", chosen_power, "\n")
cat("This may take 15-30 minutes for", ncol(datExpr), "genes...\n")
t0 <- Sys.time()

net <- blockwiseModules(
  datExpr,
  power = chosen_power,
  TOMType = "signed",
  networkType = "signed",
  minModuleSize = 50,       # Slightly larger to avoid noise modules
  deepSplit = 2,            # More aggressive splitting (0-4 scale, 2 is default WGCNA recommended)
  reassignThreshold = 0,
  mergeCutHeight = 0.30,    # Slightly more permissive merging (was 0.25)
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  maxBlockSize = 5000,
  verbose = 3,
  nThreads = 8
)

t1 <- Sys.time()
cat("\nblockwiseModules completed in", round(difftime(t1, t0, units = "mins"), 1), "minutes\n")

# Module summary
module_colors <- labels2colors(net$colors)
n_modules <- length(unique(net$colors)) - 1  # exclude module 0 (unassigned)
cat("Modules detected:", n_modules, "(excluding grey/unassigned)\n")
cat("Genes per module:\n")
print(table(net$colors))
cat("\n")

# ============================================================
# STEP 6: Extract module eigengenes
# ============================================================
cat("--- Step 6: Extract module eigengenes ---\n")

MEs <- net$MEs
cat("Eigengene matrix:", nrow(MEs), "samples x", ncol(MEs), "modules\n")

# Clean column names: ME0 -> ME_grey, ME1 -> ME_turquoise, etc.
# Keep numeric labels for clarity
colnames(MEs) <- paste0("ME", gsub("^ME", "", colnames(MEs)))

cat("Module eigengenes:\n")
cat("  Columns:", paste(head(colnames(MEs), 10), collapse = ", "),
    if (ncol(MEs) > 10) "..." else "", "\n")
cat("  Value range:", round(range(MEs), 3), "\n\n")

# ============================================================
# STEP 7: Module membership (kME) for each gene
# ============================================================
cat("--- Step 7: Compute module membership (kME) ---\n")

# kME = correlation between each gene's expression and each eigengene
# Only compute for non-grey modules
non_grey_MEs <- MEs[, colnames(MEs) != "ME0", drop = FALSE]
kME <- cor(datExpr, non_grey_MEs, use = "pairwise.complete.obs")

# Build gene-module assignment table
module_assignments <- data.table(
  gene = colnames(datExpr),
  module_numeric = net$colors,
  module_color = module_colors
)

# Add kME for the assigned module
module_assignments[, kME_own := {
  me_col <- paste0("ME", module_numeric)
  sapply(seq_len(.N), function(i) {
    col <- me_col[i]
    if (col %in% colnames(kME)) kME[gene[i], col] else NA_real_
  })
}]

# Sort by module then kME
setorder(module_assignments, module_numeric, -kME_own)
cat("Module membership table:", nrow(module_assignments), "genes\n")
cat("Top hub genes per module:\n")
for (m in sort(unique(module_assignments$module_numeric))) {
  if (m == 0) next
  top3 <- module_assignments[module_numeric == m][1:min(3, sum(module_assignments$module_numeric == m))]
  cat(sprintf("  Module %d (%s): %s\n", m,
              unique(top3$module_color),
              paste(top3$gene, collapse = ", ")))
}
cat("\n")

# ============================================================
# STEP 8: Save outputs
# ============================================================
cat("--- Step 8: Save outputs ---\n")

# Save eigengene matrix (samples x modules)
ME_matrix <- as.matrix(MEs)
saveRDS(ME_matrix, file.path(OUTDIR, "wgcna_eigengenes.rds"))
cat("Saved:", file.path(OUTDIR, "wgcna_eigengenes.rds"), "\n")

# Save module assignments
fwrite(module_assignments, file.path(OUTDIR, "wgcna_module_assignments.csv"))
cat("Saved:", file.path(OUTDIR, "wgcna_module_assignments.csv"), "\n")

# Soft threshold already saved in step 4
cat("Saved:", file.path(OUTDIR, "wgcna_soft_threshold.csv"), "\n")

cat("\n=== Summary ===\n")
cat("Input:", ncol(datExpr), "genes x", nrow(datExpr), "samples\n")
cat("Soft threshold power:", chosen_power, "\n")
cat("Modules:", n_modules, "(+1 unassigned)\n")
cat("Eigengene matrix:", nrow(ME_matrix), "samples x", ncol(ME_matrix), "modules\n")
cat("Unassigned genes (module 0):", sum(net$colors == 0), "\n")
cat("\nFinished:", as.character(Sys.time()), "\n")
