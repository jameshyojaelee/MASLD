#!/usr/bin/env Rscript
# Pseudobulk DE analysis per cell type using limma-voom
#
# Input:  results_gpu_v2/pseudobulk/{CellType}_pseudobulk.csv (genes x samples, raw counts)
# Output: results_gpu_v2/pseudobulk_de/{CellType}_de.csv (limma-voom DE results)
#
# Filters:
#   - Excludes sorted samples (CD45+, DCs, monos+macs) from Liver_Atlas
#   - Uses condition_harmonized (Healthy vs MASLD) for the contrast
#   - Includes dataset as a covariate to control for batch effects

suppressPackageStartupMessages({
  library(rhdf5)
  library(limma)
  library(edgeR)
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
scvi_dir <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2")
h5ad_file <- file.path(scvi_dir, "integrated_atlas.h5ad")
pb_dir <- file.path(scvi_dir, "pseudobulk")
out_dir <- file.path(scvi_dir, "pseudobulk_de")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# 1. Extract sample -> condition mapping from h5ad
# --------------------------------------------------------------------------
message("Extracting sample metadata from h5ad...")
stopifnot(file.exists(h5ad_file))

items <- h5ls(h5ad_file)

read_categorical <- function(h5file, key) {
  group_path <- paste0("/obs/", key)
  obs_groups <- subset(items, group == "/obs" & otype == "H5I_GROUP")$name
  if (key %in% obs_groups) {
    codes <- h5read(h5file, paste0(group_path, "/codes"))
    cats  <- h5read(h5file, paste0(group_path, "/categories"))
    return(cats[codes + 1L])
  }
  vals <- h5read(h5file, group_path)
  if (is.list(vals) && "codes" %in% names(vals)) {
    return(vals$categories[vals$codes + 1L])
  }
  as.vector(vals)
}

sample_vec  <- read_categorical(h5ad_file, "sample")
dataset_vec <- read_categorical(h5ad_file, "dataset")

# Read condition_harmonized (added by add_sample_metadata.py)
cond_vec <- tryCatch(
  read_categorical(h5ad_file, "condition_harmonized"),
  error = function(e) {
    message("  condition_harmonized not found, falling back to condition")
    read_categorical(h5ad_file, "condition")
  }
)

# Read preparation_method (added by add_sample_metadata.py)
prep_vec <- tryCatch(
  read_categorical(h5ad_file, "preparation_method"),
  error = function(e) {
    message("  preparation_method not found, no filtering applied")
    rep("unsorted", length(sample_vec))
  }
)

# Build cell-level metadata
cell_meta <- data.table(
  sample = sample_vec,
  dataset = dataset_vec,
  condition = cond_vec,
  preparation_method = prep_vec
)

# Collapse to unique sample-level metadata (majority vote for condition)
sample_meta <- cell_meta[, .(
  dataset   = dataset[1],
  condition = names(sort(table(condition), decreasing = TRUE))[1],
  preparation_method = preparation_method[1]
), by = sample]

message("  ", nrow(sample_meta), " unique samples across ",
        length(unique(sample_meta$condition)), " conditions: ",
        paste(unique(sample_meta$condition), collapse = ", "))

# Filter out sorted samples for DE analysis
de_preps <- c("nuclei", "cd45_negative", "unsorted")
n_before <- nrow(sample_meta)
sample_meta <- sample_meta[preparation_method %in% de_preps]
n_after <- nrow(sample_meta)
message("  Filtered: ", n_before, " -> ", n_after, " samples ",
        "(dropped ", n_before - n_after, " sorted samples)")

# Keep only Healthy vs MASLD for binary contrast
sample_meta <- sample_meta[condition %in% c("Healthy", "MASLD")]
message("  After condition filter (Healthy/MASLD): ", nrow(sample_meta), " samples")
message("    Condition breakdown: ",
        paste(names(table(sample_meta$condition)), table(sample_meta$condition),
              sep = "=", collapse = ", "))
message("    Dataset breakdown: ",
        paste(names(table(sample_meta$dataset)), table(sample_meta$dataset),
              sep = "=", collapse = ", "))

# --------------------------------------------------------------------------
# 2. Run limma-voom DE per cell type
# --------------------------------------------------------------------------
pb_files <- list.files(pb_dir, pattern = "_pseudobulk\\.csv$", full.names = TRUE)
message("Found ", length(pb_files), " pseudobulk files")

for (pb_file in pb_files) {
  ct_name <- gsub("_pseudobulk\\.csv$", "", basename(pb_file))
  message("\n--- Processing: ", ct_name, " ---")

  # Load count matrix (genes x samples)
  counts_dt <- fread(pb_file)
  gene_names <- counts_dt[[1]]
  counts_mat <- as.matrix(counts_dt[, -1, with = FALSE])
  rownames(counts_mat) <- gene_names
  storage.mode(counts_mat) <- "numeric"

  # Match columns (samples) to DE-eligible metadata
  col_samples <- colnames(counts_mat)
  matched_idx <- match(col_samples, sample_meta$sample)
  keep <- !is.na(matched_idx)
  if (sum(keep) < 3) {
    message("  Skipping: fewer than 3 DE-eligible samples")
    next
  }
  counts_mat <- counts_mat[, keep, drop = FALSE]
  matched <- sample_meta[matched_idx[keep]]

  # Filter samples with fewer than 50 cells contributing to the pseudobulk.
  # Cell counts per sample are derived from cell_meta (one row per cell); the
  # pseudobulk CSV columns are named by sample ID, so we count how many cells
  # from cell_meta match each sample column for this cell type.  The pseudobulk
  # files are named {CellType}_pseudobulk.csv and the h5ad obs must contain a
  # 'cell_type' or equivalent field; here we conservatively count from the
  # colnames present in this cell-type's count matrix (i.e. only samples that
  # contributed ≥1 count already appear as columns).
  cell_counts <- table(cell_meta$sample)
  sample_cell_n <- as.integer(cell_counts[colnames(counts_mat)])
  sample_cell_n[is.na(sample_cell_n)] <- 0L
  keep_mincells <- sample_cell_n >= 50L
  n_dropped_mincells <- sum(!keep_mincells)
  if (n_dropped_mincells > 0) {
    message("  Min-cells filter: dropping ", n_dropped_mincells,
            " sample(s) with <50 cells (kept ", sum(keep_mincells), ")")
    counts_mat <- counts_mat[, keep_mincells, drop = FALSE]
    matched    <- matched[keep_mincells]
  }
  if (ncol(counts_mat) < 3) {
    message("  Skipping: fewer than 3 samples after min-cells filter")
    next
  }

  # Need at least 2 samples per condition
  cond_tab <- table(matched$condition)
  valid_conds <- names(cond_tab[cond_tab >= 2])
  if (length(valid_conds) < 2) {
    message("  Skipping: need >=2 conditions with >=2 samples each")
    next
  }

  # Filter to valid conditions
  keep_cond <- matched$condition %in% valid_conds
  counts_mat <- counts_mat[, keep_cond, drop = FALSE]
  matched <- matched[keep_cond]

  message("  Samples: ", ncol(counts_mat),
          " (", paste(names(table(matched$condition)), table(matched$condition),
                      sep = "=", collapse = ", "), ")")

  # Filter low-expression genes: >=5 counts in >=3 samples
  keep_genes <- rowSums(counts_mat >= 5) >= 3
  counts_mat <- counts_mat[keep_genes, , drop = FALSE]
  message("  Genes after filtering: ", nrow(counts_mat))

  if (nrow(counts_mat) < 100) {
    message("  Skipping: too few genes (", nrow(counts_mat), ")")
    next
  }

  # Set up factors — MASLD vs Healthy (Healthy as reference)
  condition <- factor(matched$condition, levels = c("Healthy", "MASLD"))
  dataset <- factor(matched$dataset)

  # Build design matrix: include dataset as covariate if >1 dataset
  n_datasets <- length(unique(matched$dataset))
  if (n_datasets > 1) {
    design <- model.matrix(~ dataset + condition)
    coef_name <- "conditionMASLD"
    message("  Design: ~ dataset + condition (", n_datasets, " datasets)")
  } else {
    design <- model.matrix(~ condition)
    coef_name <- "conditionMASLD"
    message("  Design: ~ condition (single dataset)")
  }

  tryCatch({
    dge <- DGEList(counts = counts_mat)
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)

    # Extract MASLD coefficient
    coef_idx <- which(colnames(design) == coef_name)
    if (length(coef_idx) == 0) {
      message("  WARNING: coefficient '", coef_name, "' not found in design")
      message("  Available: ", paste(colnames(design), collapse = ", "))
      next
    }

    res <- topTable(fit, coef = coef_idx, number = Inf, sort.by = "none")
    res$gene <- rownames(res)

    setnames(setDT(res),
             old = c("logFC", "adj.P.Val", "P.Value", "t", "B"),
             new = c("logFC", "padj", "pvalue", "t_stat", "B"),
             skip_absent = TRUE)

    res$cell_type <- ct_name
    res$contrast <- "MASLD_vs_Healthy"

    out_file <- file.path(out_dir, paste0(ct_name, "_de.csv"))
    fwrite(res, out_file)
    message("  Saved: ", out_file, " (", nrow(res), " genes, ",
            sum(res$padj < 0.05, na.rm = TRUE), " sig at padj<0.05)")
  }, error = function(e) {
    message("  ERROR in limma-voom: ", e$message)
  })
}

message("\nPseudobulk DE complete. Results in: ", out_dir)
