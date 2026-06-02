#!/usr/bin/env Rscript
# Export single-cell reference to CIBERSORTx-compatible format
# Input: SingleCellExperiment or H5AD file
# Output: Tab-separated matrix (genes x cells) for CIBERSORTx upload

# Force conda library path to prevent user library contamination
conda_lib <- Sys.getenv("R_LIBS_SITE")
if (nzchar(conda_lib) && dir.exists(conda_lib)) {
  .libPaths(conda_lib)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 11_export_cibersortx.R <input_sce.rds> <output_prefix> <label_key>")
}

input_sce <- args[[1]]
output_prefix <- args[[2]]
label_key <- args[[3]]

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
})

cat("Loading SCE:", input_sce, "\n")
sce <- readRDS(input_sce)

# Get counts matrix
counts <- assay(sce, "counts")
if (is.null(counts)) {
  counts <- assay(sce, 1)  
}

# Get cell type labels
labels <- colData(sce)[[label_key]]
if (is.null(labels)) {
  # Try alternative keys
  for (alt_key in c("celltype", "consensus_label", "celltypist_label")) {
    if (alt_key %in% colnames(colData(sce))) {
      labels <- colData(sce)[[alt_key]]
      cat("Using alternative label key:", alt_key, "\n")
      break
    }
  }
}

if (is.null(labels)) {
  stop("No valid label key found in SCE")
}

# Create sample names with cell type suffix for CIBERSORTx
cell_names <- paste0(labels, "_", seq_len(ncol(sce)))

# Export 1: Single-cell reference matrix (genes x cells, tab-separated)
# CIBERSORTx format: Gene\tCell1\tCell2\t...
sc_matrix_file <- paste0(output_prefix, "_sc_matrix.txt")
cat("Exporting single-cell matrix to:", sc_matrix_file, "\n")

# Convert sparse to dense if needed (CIBERSORTx requires dense)
counts_dense <- as.matrix(counts)

# Add gene names as first column
out_df <- data.frame(Gene = rownames(counts_dense), counts_dense, check.names = FALSE)
colnames(out_df)[-1] <- cell_names

write.table(out_df, sc_matrix_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Export 2: Cell type labels file
labels_file <- paste0(output_prefix, "_labels.txt")
cat("Exporting labels to:", labels_file, "\n")
labels_df <- data.frame(Cell = cell_names, CellType = labels)
write.table(labels_df, labels_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Export 3: Phenotype classes (for signature matrix generation)
phenoclasses_file <- paste0(output_prefix, "_phenoclasses.txt")
cat("Exporting phenotype classes to:", phenoclasses_file, "\n")
# CIBERSORTx phenotype file has 2 rows: sample names and class labels
pheno_df <- rbind(cell_names, as.character(labels))
write.table(pheno_df, phenoclasses_file, sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

# Summary statistics
cat("\n=== Summary ===\n")
cat("Total cells:", ncol(sce), "\n")
cat("Total genes:", nrow(sce), "\n")
cat("Cell types:\n")
print(table(labels))

cat("\nOutput files:\n")
cat(" - Single-cell matrix:", sc_matrix_file, "\n")
cat(" - Cell labels:", labels_file, "\n")
cat(" - Phenotype classes:", phenoclasses_file, "\n")
cat("\nUpload", sc_matrix_file, "and", phenoclasses_file, "to CIBERSORTx\n")
cat("(https://cibersortx.stanford.edu/) to generate signature matrix.\n")
