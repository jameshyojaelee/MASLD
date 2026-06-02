# WGCNA Pipeline - Step 01: Data Preprocessing
# This script loads expression data, filters outliers, and prepares it for WGCNA.

library(WGCNA)
library(tidyverse)
library(DESeq2)

options(stringsAsFactors = FALSE)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
# Path to variance stabilized data (rlog or vst transformed)
# Expects a matrix or dataframe where rows are GENES and columns are SAMPLES
DATA_FILE <- "path/to/your/vst_expression_data.csv" 

# Path to metadata file
# Expects CSV where rows are SAMPLES and columns are traits/conditions
METADATA_FILE <- "path/to/your/metadata.csv"

# Output directory for processed data
OUTPUT_DIR <- "../results/01_preprocessing"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# DATA LOADING
# ==============================================================================
cat("Loading data...\n")

# Load Expression Data
if (grepl(".rds$", DATA_FILE, ignore.case = TRUE)) {
  exp_data <- readRDS(DATA_FILE)
} else {
  exp_data <- read.csv(DATA_FILE, row.names = 1, check.names = FALSE)
}

# Load Metadata
if (grepl(".rds$", METADATA_FILE, ignore.case = TRUE)) {
  metadata <- readRDS(METADATA_FILE)
} else {
  metadata <- read.csv(METADATA_FILE, row.names = 1)
}

# ==============================================================================
# DATA CLEANING AND ALIGNMENT
# ==============================================================================
cat("Checking data integrity...\n")

# Transpose expression data: WGCNA expects Rows = SAMPLES, Columns = GENES
datExpr0 <- as.data.frame(t(exp_data))

# Check for genes and samples with too many missing values
gsg <- goodSamplesGenes(datExpr0, verbose = 3)

if (!gsg$allOK) {
  # Optionally remove genes/samples that don't pass the check
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Align expression data with metadata
common_samples <- intersect(rownames(datExpr0), rownames(metadata))
datExpr <- datExpr0[common_samples, ]
metadata_aligned <- metadata[common_samples, ]

cat(paste("Number of samples retained:", nrow(datExpr), "\n"))
cat(paste("Number of genes retained:", ncol(datExpr), "\n"))

# ==============================================================================
# SAMPLE CLUSTERING (OUTLIER DETECTION)
# ==============================================================================
cat("Clustering samples to detect outliers...\n")

sampleTree <- hclust(dist(datExpr), method = "average")

pdf(file = file.path(OUTPUT_DIR, "sample_clustering_tree.pdf"), width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# Check for outliers: Plot a line if you want to set a cut height
# abline(h = 15, col = "red") 
dev.off()

# Save cleaned data for next steps
save(datExpr, metadata_aligned, file = file.path(OUTPUT_DIR, "WGCNA_data_input.RData"))

cat("Preprocessing complete. Data saved to:", file.path(OUTPUT_DIR, "WGCNA_data_input.RData"), "\n")
