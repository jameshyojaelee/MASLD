# WGCNA Pipeline - Step 03: Network Construction and Module Detection
# This script constructs the co-expression network and identifies modules.

library(WGCNA)
library(tidyverse)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  # Defaults for testing if run manually without args
  DATASET_NAME <- "Patient"
  POWER <- 6
  cat("No arguments provided. Using defaults: Dataset =", DATASET_NAME, ", Power =", POWER, "\n")
  cat("Usage: Rscript 03_network_construction.R <Dataset_Name> <Power>\n")
} else {
  DATASET_NAME <- args[1]
  POWER <- as.numeric(args[2])
}

# ==============================================================================
# CONFIGURATION
# ==============================================================================
INPUT_DATA <- file.path("../results/01_preprocessing", paste0(DATASET_NAME, ".RData"))
OUTPUT_DIR <- "../results/03_network_construction"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(paste("Processing Dataset:", DATASET_NAME, "\n"))


# Network construction parameters
MIN_MODULE_SIZE <- 30
MERGE_CUT_HEIGHT <- 0.25 # Merge modules with > 0.75 correlation

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("Loading preprocessed data...\n")
load(INPUT_DATA)

# Handle different variable names from preprocessing steps
if (exists("datExpr_patient")) {
  datExpr <- datExpr_patient
} else if (exists("datExpr_other")) {
  datExpr <- datExpr_other
} else if (exists("datExpr_inhouse")) {
  datExpr <- datExpr_inhouse
}

if (!exists("datExpr")) {
  stop("Could not find 'datExpr' variable. Checked: datExpr_patient, datExpr_other, datExpr_inhouse.")
}

# ==============================================================================
# NETWORK CONSTRUCTION
# ==============================================================================
cat(paste("Constructing network with power =", POWER, "...\n"))

# blockwiseModules is efficient for large datasets and handles memory automatically
net <- blockwiseModules(datExpr, power = POWER,
                        TOMType = "unsigned", minModuleSize = MIN_MODULE_SIZE,
                        reassignThreshold = 0, mergeCutHeight = MERGE_CUT_HEIGHT,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = file.path(OUTPUT_DIR, "TOM"),
                        verbose = 3)

# ==============================================================================
# VISUALIZATION
# ==============================================================================
cat("Generating dendrogram plots...\n")

# Convert labels to colors for plotting
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

pdf(file = file.path(OUTPUT_DIR, "module_dendrogram.pdf"), width = 12, height = 9)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# ==============================================================================
# SAVE RESULTS
# ==============================================================================
save(MEs, moduleLabels, moduleColors, geneTree, net, 
     file = file.path(OUTPUT_DIR, paste0("WGCNA_network_", DATASET_NAME, ".RData")))

cat("Network construction complete. Data saved to:", file.path(OUTPUT_DIR, paste0("WGCNA_network_", DATASET_NAME, ".RData")), "\n")
