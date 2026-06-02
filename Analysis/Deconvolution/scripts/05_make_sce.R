args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 05_make_sce.R <input_h5ad> <output_rds> <label_key>")
}

input_h5ad <- args[[1]]
output_rds <- args[[2]]
label_key  <- args[[3]]

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

cat("Reading h5ad via rhdf5:", input_h5ad, "\n")

# ---------- Read sparse matrix X (CSR format in h5ad) ----------
x_data    <- h5read(input_h5ad, "/X/data")
x_indices <- h5read(input_h5ad, "/X/indices")
x_indptr  <- h5read(input_h5ad, "/X/indptr")

# h5ad stores cells × genes (CSR). R's dgCMatrix is CSC (genes × cells).
# Build CSR then transpose to get genes × cells dgCMatrix.
n_cells <- length(x_indptr) - 1L
n_genes <- as.integer(max(x_indices) + 1L)  # 0-based indices

cat("  Sparse matrix:", n_cells, "cells x", n_genes, "genes,",
    length(x_data), "non-zero entries\n")

# h5ad CSR (cells×genes) has the same structure as CSC (genes×cells):
#   CSR row_ptr → CSC col_ptr, CSR col_idx → CSC row_idx
# So we can construct genes×cells dgCMatrix directly without transposing.
counts <- new("dgCMatrix",
              i = as.integer(x_indices),
              p = as.integer(x_indptr),
              x = as.numeric(x_data),
              Dim = c(n_genes, n_cells))

rm(x_data, x_indices, x_indptr)
gc()

# ---------- Read obs (cell metadata) ----------
read_categorical <- function(h5file, path) {
  cats  <- h5read(h5file, paste0(path, "/categories"))
  codes <- h5read(h5file, paste0(path, "/codes"))
  cats[codes + 1L]  # 0-based → 1-based
}

obs_index <- h5read(input_h5ad, "/obs/_index")

# Discover obs columns
obs_items <- h5ls(input_h5ad, recursive = TRUE)
obs_items <- obs_items[obs_items$group == "/obs" & obs_items$name != "_index", ]

obs_df <- data.frame(row.names = obs_index)
for (i in seq_len(nrow(obs_items))) {
  col_name <- obs_items$name[i]
  col_path <- paste0("/obs/", col_name)
  if (obs_items$otype[i] == "H5I_GROUP") {
    # Categorical
    obs_df[[col_name]] <- read_categorical(input_h5ad, col_path)
  } else {
    # Direct dataset
    obs_df[[col_name]] <- h5read(input_h5ad, col_path)
  }
}

cat("  obs columns:", paste(colnames(obs_df), collapse = ", "), "\n")

# ---------- Read var (gene metadata) ----------
var_index <- h5read(input_h5ad, "/var/_index")
gene_ids  <- tryCatch(h5read(input_h5ad, "/var/gene_ids"), error = function(e) NULL)

var_df <- data.frame(row.names = var_index)
if (!is.null(gene_ids)) var_df$gene_ids <- gene_ids

# Set dimension names
rownames(counts) <- var_index
colnames(counts) <- obs_index

# ---------- Validate label_key ----------
if (!label_key %in% colnames(obs_df)) {
  stop(paste0("Label key '", label_key, "' not found in obs. Available: ",
              paste(colnames(obs_df), collapse = ", ")))
}

# ---------- Standardize MuSiC-required fields ----------
obs_df$celltype <- as.character(obs_df[[label_key]])

if ("sample" %in% colnames(obs_df)) {
  obs_df$sampleID <- as.character(obs_df$sample)
} else if ("run_id" %in% colnames(obs_df)) {
  obs_df$sampleID <- as.character(obs_df$run_id)
} else {
  obs_df$sampleID <- as.character(seq_len(nrow(obs_df)))
}

cat("  Cell types:", length(unique(obs_df$celltype)), "\n")
cat("  Samples:", length(unique(obs_df$sampleID)), "\n")

# ---------- Build SCE ----------
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = DataFrame(obs_df),
  rowData = DataFrame(var_df)
)

cat("SCE:", nrow(sce), "genes x", ncol(sce), "cells\n")

saveRDS(sce, output_rds)
cat("Saved:", output_rds, "\n")
