# WGCNA Pipeline - Step 02: Soft Thresholding Power Selection
# This script analyzes the scale-free topology fit to help select the soft-thresholding power.

library(WGCNA)
library(tidyverse)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  DATASET_NAME <- "Patient" # Default
  cat("No arguments provided. Using default: Dataset =", DATASET_NAME, "\n")
  cat("Usage: Rscript 02_soft_thresholding.R <Dataset_Name>\n")
} else {
  DATASET_NAME <- args[1]
}

# ==============================================================================
# CONFIGURATION
# ==============================================================================
INPUT_DATA <- file.path("../results/01_preprocessing", paste0(DATASET_NAME, ".RData"))
OUTPUT_DIR <- "../results/02_soft_thresholding"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(paste("Processing Dataset:", DATASET_NAME, "\n"))

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("Loading preprocessed data...\n")
load(INPUT_DATA)
# Creates 'datExpr' and 'metadata_aligned'

# ==============================================================================
# ANALYSIS
# ==============================================================================
# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

cat("Picking soft thresholding power (this may take a while)...\n")
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# ==============================================================================
# PLOTTING
# ==============================================================================
cat("Generating plots...\n")
pdf(file = file.path(OUTPUT_DIR, paste0("soft_thresholding_power_", DATASET_NAME, ".pdf")), width = 9, height = 5)
par(mfrow = c(1, 2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
# Determine an approximate threshold line (usually 0.8 or 0.9)
abline(h = 0.90, col = "red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")

dev.off()

cat("Analysis complete. Check plots in:", file.path(OUTPUT_DIR, "soft_thresholding_power.pdf"), "\n")
cat("Remember the chosen power for the next step (Step 03).\n")
