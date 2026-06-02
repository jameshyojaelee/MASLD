#!/usr/bin/env Rscript
# BayesPrism / InstaPrism Deconvolution
# Uses InstaPrism for faster deconvolution with BayesPrism-equivalent results
# 
# References:
# - BayesPrism: https://github.com/Danko-Lab/BayesPrism
# - InstaPrism: https://github.com/humengying0907/InstaPrism (faster implementation)

# Force conda library path to prevent user library contamination
conda_lib <- Sys.getenv("R_LIBS_SITE")
if (nzchar(conda_lib) && dir.exists(conda_lib)) {
  .libPaths(conda_lib)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript 12_run_bayesprism.R <ref_sce.rds> <bulk_counts.tsv> <output_dir> <dataset_name> [label_key]")
}

ref_sce_path <- args[[1]]
bulk_counts_path <- args[[2]]
output_dir <- args[[3]]
dataset_name <- args[[4]]
label_key <- if (length(args) >= 5) args[[5]] else "celltype"

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
})

# Check if InstaPrism is available, otherwise use BayesPrism
use_instaprism <- requireNamespace("InstaPrism", quietly = TRUE)
if (use_instaprism) {
  library(InstaPrism)
  cat("Using InstaPrism (fast mode)\n")
} else if (requireNamespace("BayesPrism", quietly = TRUE)) {
  library(BayesPrism)
  cat("Using BayesPrism\n")
} else {
  stop("Neither InstaPrism nor BayesPrism is installed. Please install one:\n",
       "  devtools::install_github('humengying0907/InstaPrism')  # Fast\n",
       "  devtools::install_github('Danko-Lab/BayesPrism')       # Original")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading reference SCE:", ref_sce_path, "\n")
sce <- readRDS(ref_sce_path)

# Extract reference counts and labels
ref_counts <- assay(sce, "counts")
if (is.null(ref_counts)) ref_counts <- assay(sce, 1)

cell_types <- colData(sce)[[label_key]]
if (is.null(cell_types)) {
  for (alt in c("celltype", "consensus_label", "celltypist_label")) {
    if (alt %in% colnames(colData(sce))) {
      cell_types <- colData(sce)[[alt]]
      cat("Using label key:", alt, "\n")
      break
    }
  }
}

# Load bulk counts
cat("Loading bulk counts:", bulk_counts_path, "\n")
bulk <- read.delim(bulk_counts_path, row.names = 1, check.names = FALSE)
bulk <- as.matrix(bulk)

# Handle bulk counts matrix orientation
# InstaPrism expects genes x samples
# BayesPrism expects samples x genes
if (use_instaprism) {
  # For InstaPrism: ensure genes x samples
  if (nrow(bulk) < ncol(bulk) && ncol(bulk) > 100) {
    # Assuming fewer samples than genes, so if rows < cols, it might be samples x genes
    cat("Transposing bulk matrix for InstaPrism (samples as rows -> genes as rows)\n")
    bulk <- t(bulk)
  }
} else {
  # For BayesPrism: ensure samples x genes
  if (nrow(bulk) > ncol(bulk) && nrow(bulk) > 100) {
    cat("Transposing bulk matrix for BayesPrism (genes as rows -> samples as rows)\n")
    bulk <- t(bulk)
  }
}

# Handle reference matrix orientation
# InstaPrism expects genes x cells (standard SCE counts)
# BayesPrism expects cells x genes (needs transpose)

if (use_instaprism) {
  ref_counts_final <- as.matrix(ref_counts)
  cat("Using untransposed reference for InstaPrism: genes x cells\n")
} else {
  ref_counts_final <- t(as.matrix(ref_counts))
  cat("Transposing reference for BayesPrism: cells x genes\n")
}

cat("Reference dimensions:", nrow(ref_counts_final), "x", ncol(ref_counts_final), "\n")
cat("Bulk dimensions:", nrow(bulk), "samples x", ncol(bulk), "genes\n")
cat("Cell types:\n")
print(table(cell_types))

# Run deconvolution
if (use_instaprism) {
  # InstaPrism workflow
  cat("\nPreparing InstaPrism reference...\n")
  
  # InstaPrism refPrepare - note: uses sc_Expr not sc.eset
  ref_obj <- refPrepare(
    sc_Expr = ref_counts_final,  # genes x cells matrix
    cell.type.labels = as.character(cell_types),
    cell.state.labels = as.character(cell_types)  # Use cell type as state
  )
  
  cat("Running InstaPrism deconvolution...\n")
  result <- InstaPrism(
    bulk_Expr = bulk,  # note: bulk_Expr not bulk.eset
    refPhi_cs = ref_obj  # note: refPhi_cs not ref
  )
  
  # Extract proportions - cell type level
  proportions <- t(result@Post.ini.ct@theta)
  
} else {
  # BayesPrism workflow
  cat("\nCreating BayesPrism object...\n")
  
  prism <- new.prism(
    reference = ref_counts_final,
    mixture = bulk,
    cell.type.labels = as.character(cell_types),
    cell.state.labels = as.character(cell_types),
    key = NULL,
    input.type = "count.matrix"
  )
  
  cat("Running BayesPrism deconvolution (this may take a while)...\n")
  result <- run.prism(prism, n.cores = parallel::detectCores() - 1)
  
  # Extract proportions
  proportions <- get.fraction(result, which = "final", state.or.type = "type")
}

# Save results
prop_file <- file.path(output_dir, paste0(dataset_name, "_bayesprism_proportions.tsv"))
cat("Saving proportions to:", prop_file, "\n")
write.table(proportions, prop_file, sep = "\t", quote = FALSE, row.names = TRUE)

# Save full result object
rds_file <- file.path(output_dir, paste0(dataset_name, "_bayesprism_result.rds"))
cat("Saving result object to:", rds_file, "\n")
saveRDS(result, rds_file)

cat("\n=== BayesPrism Deconvolution Complete ===\n")
cat("Proportions file:", prop_file, "\n")
cat("Result object:", rds_file, "\n")
