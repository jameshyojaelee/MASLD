# WGCNA Pipeline - Step 04: Module-Trait Relationship
# This script correlates module eigengenes with external traits.

library(WGCNA)
library(tidyverse)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  DATASET_NAME <- "Patient" # Default
  cat("No arguments provided. Using default: Dataset =", DATASET_NAME, "\n")
  cat("Usage: Rscript 04_module_trait_relationship.R <Dataset_Name>\n")
} else {
  DATASET_NAME <- args[1]
}

# ==============================================================================
# CONFIGURATION
# ==============================================================================
INPUT_DATA_PREP <- file.path("../results/01_preprocessing", paste0(DATASET_NAME, ".RData"))
INPUT_DATA_NET <- file.path("../results/03_network_construction", paste0("WGCNA_network_", DATASET_NAME, ".RData"))
OUTPUT_DIR <- "../results/04_module_trait"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(paste("Processing Dataset:", DATASET_NAME, "\n"))

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("Loading data...\n")
load(INPUT_DATA_PREP) # datExpr, metadata_aligned
load(INPUT_DATA_NET)  # MEs, moduleLabels, moduleColors, net

# Handle dataset-specific variables
if (exists("metadata_patient")) {
  datExpr <- datExpr_patient
  metadata <- metadata_patient
} else if (exists("metadata_other")) {
  datExpr <- datExpr_other
  metadata <- metadata_other
} else if (exists("metadata_inhouse")) {
  datExpr <- datExpr_inhouse
  metadata <- metadata_inhouse
}

# Ensure metadata consists of numeric traits for correlation
traits <- metadata %>%
  mutate_if(is.character, as.factor) %>%
  mutate_if(is.factor, as.numeric) # robust conversion for correlation

# ==============================================================================
# CORRELATION ANALYSIS
# ==============================================================================
cat("Calculating module-trait correlations...\n")

# Order MEs for visualization
MEs_ordered <- orderMEs(MEs)

moduleTraitCor <- cor(MEs_ordered, traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# ==============================================================================
# VISUALIZATION
# ==============================================================================
cat("Generating heatmap...\n")

pdf(file = file.path(OUTPUT_DIR, paste0("module_trait_relationships_", DATASET_NAME, ".pdf")), width = 10, height = 10)
# Will display correlations and their p-values
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs_ordered),
               ySymbols = names(MEs_ordered),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# ==============================================================================
# GENE SIGNIFICANCE & MODULE MEMBERSHIP
# ==============================================================================
# Calculate Module Membership (MM) and Gene Significance (GS)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs_ordered, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))

# Calculate Gene Significance for EACH trait
geneTraitSignificance <- as.data.frame(cor(datExpr, traits, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(datExpr)))

# Save comprehensive results
save(geneModuleMembership, MMPvalue, geneTraitSignificance, GSPvalue, 
     file = file.path(OUTPUT_DIR, paste0("WGCNA_gene_significance_", DATASET_NAME, ".RData")))

cat("Module-trait analysis complete. Plot saved to:", file.path(OUTPUT_DIR, paste0("module_trait_relationships_", DATASET_NAME, ".pdf")), "\n")
