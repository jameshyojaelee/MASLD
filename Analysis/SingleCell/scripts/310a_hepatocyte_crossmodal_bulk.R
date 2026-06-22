#!/usr/bin/env Rscript
# 310a: Cross-modal — Bulk RNA-seq pseudobulk DE per hepatocyte subtype
#
# For each hepatocyte subtype, build pseudobulk count matrices, run
# limma-voom DE (MASLD vs Healthy), and compute Spearman correlation
# with the dream mega-analysis t-statistics.
#
# Input:
#   hepatocyte_subtypes/hepatocyte_subtype_metadata.csv  (from 309)
#   hepatocyte_subtypes/hepatocyte_atlas.h5ad  (raw counts via rhdf5)
#   canonical_deg_results.csv  (from integration pipeline; LVQW canonical)
#
# Output (to hepatocyte_subtypes/crossmodal/bulk/):
#   {subtype}_de.csv          — limma-voom DE results per subtype
#   subtype_bulk_correlation.csv  — Spearman rho per subtype vs dream
#
# Environment: rnaseq (CPU)

suppressPackageStartupMessages({
  library(rhdf5)
  library(limma)
  library(edgeR)
  library(data.table)
  library(Matrix)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SC_DIR   <- file.path(BASE, "Analysis/SingleCell")
RESULTS  <- file.path(SC_DIR, "results_gpu_v2")
SUB_DIR  <- file.path(RESULTS, "hepatocyte_subtypes")
OUT_DIR  <- file.path(SUB_DIR, "crossmodal", "bulk")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Dream results for correlation
DREAM_PATH <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")

# ---------------------------------------------------------------------------
# 1. Load subtype metadata
# ---------------------------------------------------------------------------
message("Loading subtype metadata...")
meta_path <- file.path(SUB_DIR, "hepatocyte_subtype_metadata.csv")
stopifnot(file.exists(meta_path))
meta <- fread(meta_path)
# Script 309 uses pandas to_csv(index=True), so the first column is the cell
# barcode index (unnamed -> V1 in fread)
if ("V1" %in% names(meta)) {
  setnames(meta, "V1", "cell_barcode")
} else if (!("cell_barcode" %in% names(meta))) {
  # First column is the barcode regardless of its name
  setnames(meta, names(meta)[1], "cell_barcode")
}
message("  ", nrow(meta), " cells, ", length(unique(meta$hepatocyte_subtype)), " subtypes")

# ---------------------------------------------------------------------------
# 2. Load raw counts from h5ad (hepatocyte atlas)
# ---------------------------------------------------------------------------
message("Loading raw counts from hepatocyte_atlas.h5ad via rhdf5...")
h5_path <- file.path(SUB_DIR, "hepatocyte_atlas.h5ad")
stopifnot(file.exists(h5_path))

items <- h5ls(h5_path)

# Helper: read a sparse matrix from an h5ad group path, handling CSR/CSC
read_sparse_h5 <- function(h5file, group_path, n_obs, n_var) {
  data_vec <- h5read(h5file, paste0(group_path, "/data"))
  indices  <- h5read(h5file, paste0(group_path, "/indices"))
  indptr   <- h5read(h5file, paste0(group_path, "/indptr"))

  # Detect encoding type: CSR vs CSC
  enc_type <- tryCatch(
    h5readAttributes(h5file, group_path)[["encoding-type"]],
    error = function(e) NULL
  )

  if (!is.null(enc_type) && grepl("csr", enc_type, ignore.case = TRUE)) {
    # CSR: indptr is row pointer (length n_obs+1), indices are column indices
    # sparseMatrix with j= creates dgRMatrix-like via i/j/x triplet then stores as dgCMatrix
    mat <- sparseMatrix(
      j = as.integer(indices) + 1L,
      p = as.integer(indptr),
      x = as.numeric(data_vec),
      dims = c(n_obs, n_var),
      repr = "C"  # Return as CSC (dgCMatrix) — internally transposes the CSR
    )
  } else if (!is.null(enc_type) && grepl("csc", enc_type, ignore.case = TRUE)) {
    # CSC: indptr is column pointer, indices are row indices
    mat <- sparseMatrix(
      i = as.integer(indices) + 1L,
      p = as.integer(indptr),
      x = as.numeric(data_vec),
      dims = c(n_obs, n_var)
    )
  } else {
    # Default: assume CSR (most common h5ad format)
    # Infer from indptr length: if length = n_obs+1, it's row-oriented (CSR)
    if (length(indptr) == n_obs + 1L) {
      mat <- sparseMatrix(
        j = as.integer(indices) + 1L,
        p = as.integer(indptr),
        x = as.numeric(data_vec),
        dims = c(n_obs, n_var),
        repr = "C"
      )
    } else {
      # Column-oriented (CSC)
      mat <- sparseMatrix(
        i = as.integer(indices) + 1L,
        p = as.integer(indptr),
        x = as.numeric(data_vec),
        dims = c(n_obs, n_var)
      )
    }
  }
  return(mat)
}

# Gene names and cell barcodes
gene_names <- as.character(h5read(h5_path, "var/_index"))
obs_names  <- as.character(h5read(h5_path, "obs/_index"))
n_obs <- length(obs_names)
n_var <- length(gene_names)

# Try layers/counts first (stored by Script 308), fall back to raw/X, then X
layer_names <- items[items$group == "/layers", "name"]
raw_names   <- items[items$group == "/raw", "name"]

if ("counts" %in% layer_names) {
  message("  Reading from layers/counts...")
  layer_path <- "layers/counts"
  test_read <- h5read(h5_path, layer_path)
  if (is.list(test_read)) {
    # Sparse format
    counts_mat <- read_sparse_h5(h5_path, layer_path, n_obs, n_var)
    # h5ad stores obs x var; we want genes x cells
    counts_mat <- t(counts_mat)
  } else {
    # Dense: h5ad stores obs x var
    counts_mat <- t(test_read)
  }
  rm(test_read); gc()
  message("  Read from layers/counts")
} else if ("X" %in% raw_names) {
  message("  Reading from raw/X...")
  test_read <- h5read(h5_path, "raw/X")
  if (is.list(test_read)) {
    counts_mat <- read_sparse_h5(h5_path, "raw/X", n_obs, n_var)
    counts_mat <- t(counts_mat)
  } else {
    counts_mat <- t(test_read)
  }
  rm(test_read); gc()
  message("  Read from raw/X")
} else {
  message("  Reading from X...")
  test_read <- h5read(h5_path, "X")
  if (is.list(test_read)) {
    counts_mat <- read_sparse_h5(h5_path, "X", n_obs, n_var)
    counts_mat <- t(counts_mat)
  } else {
    counts_mat <- t(test_read)
  }
  rm(test_read); gc()
  message("  Read from X")
}

# Now counts_mat is genes x cells
rownames(counts_mat) <- gene_names
colnames(counts_mat) <- obs_names

message("  Count matrix: ", nrow(counts_mat), " genes x ", ncol(counts_mat), " cells")

# Ensure counts_mat is a regular matrix for rowSums performance
if (is(counts_mat, "sparseMatrix")) {
  message("  Note: keeping sparse format for memory efficiency; converting subsets on the fly")
  is_sparse <- TRUE
} else {
  is_sparse <- FALSE
}

# ---------------------------------------------------------------------------
# 3. Load dream results
# ---------------------------------------------------------------------------
message("Loading dream results...")
stopifnot(file.exists(DREAM_PATH))
dream <- fread(DREAM_PATH)
# dream has columns: gene, logFC, AveExpr, t, P.Value, padj, ...
# Use 'symbol' column (gene symbols) — 'gene' has Ensembl IDs which
# won't match the scRNA var_names (gene symbols like TSPAN6)
dream_t <- dream[, .(gene = symbol, dream_t = t, bulk_logFC = logFC, bulk_padj = padj)]
setkey(dream_t, gene)
message("  ", nrow(dream_t), " genes with dream t-statistics")

# ---------------------------------------------------------------------------
# 4. Pseudobulk DE per subtype
# ---------------------------------------------------------------------------
subtypes <- as.character(sort(unique(meta$hepatocyte_subtype)))
message("Running pseudobulk DE for ", length(subtypes), " subtypes...")

corr_rows <- list()

for (st in subtypes) {
  message("\n--- Subtype: ", st, " ---")

  # Cells in this subtype
  cell_meta_st <- meta[hepatocyte_subtype == st]
  cells <- cell_meta_st$cell_barcode
  message("  ", length(cells), " cells")

  # Build pseudobulk: sum counts per sample
  samples <- unique(cell_meta_st$sample)

  pb_list <- list()
  sample_info <- list()
  for (samp in samples) {
    samp_cells <- cell_meta_st[sample == samp]
    if (nrow(samp_cells) < 10) next  # min 10 cells per sample

    samp_idx <- match(samp_cells$cell_barcode, obs_names)
    samp_idx <- samp_idx[!is.na(samp_idx)]
    if (length(samp_idx) < 10) next

    if (is_sparse) {
      pb_list[[samp]] <- Matrix::rowSums(counts_mat[, samp_idx, drop = FALSE])
    } else {
      pb_list[[samp]] <- rowSums(counts_mat[, samp_idx, drop = FALSE])
    }
    sample_info[[samp]] <- data.table(
      sample = samp,
      dataset = samp_cells$dataset[1],
      condition = samp_cells$condition[1],
      n_cells = nrow(samp_cells)
    )
  }

  if (length(pb_list) < 3) {
    message("  Skipping: fewer than 3 samples with >= 10 cells")
    next
  }

  pb_mat <- do.call(cbind, pb_list)
  sample_dt <- rbindlist(sample_info)

  # Map condition to binary: NAFLD, NASH, Cirrhotic -> MASLD
  sample_dt[, condition_binary := fifelse(
    condition %in% c("MASLD", "NAFLD", "NASH", "Cirrhotic"),
    "MASLD", condition)]
  sample_dt <- sample_dt[condition_binary %in% c("Healthy", "MASLD")]
  pb_mat <- pb_mat[, sample_dt$sample, drop = FALSE]

  cond_tab <- table(sample_dt$condition_binary)
  if (length(cond_tab) < 2 || any(cond_tab < 2)) {
    message("  Skipping: need >=2 per condition (",
            paste(names(cond_tab), cond_tab, sep = "=", collapse = ", "), ")")
    next
  }

  message("  Samples: ", ncol(pb_mat), " (",
          paste(names(cond_tab), cond_tab, sep = "=", collapse = ", "), ")")

  # Filter low-expression genes: >= 5 counts in >= 3 samples
  keep_genes <- rowSums(pb_mat >= 5) >= 3
  pb_mat <- pb_mat[keep_genes, , drop = FALSE]
  message("  Genes after filtering: ", nrow(pb_mat))

  if (nrow(pb_mat) < 100) {
    message("  Skipping: too few genes (", nrow(pb_mat), ")")
    next
  }

  # Design matrix
  condition <- factor(sample_dt$condition_binary, levels = c("Healthy", "MASLD"))
  dataset <- factor(sample_dt$dataset)
  n_datasets <- length(unique(sample_dt$dataset))

  if (n_datasets > 1) {
    design <- model.matrix(~ dataset + condition)
    message("  Design: ~ dataset + condition (", n_datasets, " datasets)")
  } else {
    design <- model.matrix(~ condition)
    message("  Design: ~ condition (single dataset)")
  }

  tryCatch({
    dge <- DGEList(counts = pb_mat)
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)

    coef_idx <- which(colnames(design) == "conditionMASLD")
    if (length(coef_idx) == 0) {
      message("  WARNING: coefficient 'conditionMASLD' not found in design")
      message("  Available: ", paste(colnames(design), collapse = ", "))
      next
    }

    res <- topTable(fit, coef = coef_idx, number = Inf, sort.by = "none")
    res$gene <- rownames(res)
    setnames(setDT(res),
             old = c("logFC", "adj.P.Val", "P.Value", "t", "B"),
             new = c("logFC", "padj", "pvalue", "t_stat", "B"),
             skip_absent = TRUE)
    res$subtype <- st
    res$contrast <- "MASLD_vs_Healthy"

    fwrite(res, file.path(OUT_DIR, paste0(st, "_de.csv")))
    message("  Saved: ", nrow(res), " genes, ",
            sum(res$padj < 0.05, na.rm = TRUE), " sig at padj<0.05")

    # Correlation with dream t-statistics
    shared <- intersect(res$gene, dream_t$gene)
    if (length(shared) > 100) {
      res_sub <- res[gene %in% shared]
      dream_sub <- dream_t[shared]
      # Align by gene name
      m_res   <- match(shared, res_sub$gene)
      m_dream <- match(shared, dream_sub$gene)
      rho <- cor(res_sub[m_res, t_stat],
                 dream_sub[m_dream, dream_t],
                 method = "spearman", use = "complete.obs")

      # Direction concordance among co-significant genes (padj < 0.05 in both)
      res_sig_genes <- res_sub[padj < 0.05, gene]
      dream_sig_genes <- dream_t[bulk_padj < 0.05, gene]
      cosig_genes <- intersect(res_sig_genes, dream_sig_genes)
      n_cosig <- length(cosig_genes)

      if (n_cosig > 0) {
        cosig_res   <- res_sub[match(cosig_genes, gene)]
        cosig_dream <- dream_t[match(cosig_genes, gene)]
        dir_conc <- mean(sign(cosig_res$logFC) == sign(cosig_dream$bulk_logFC),
                         na.rm = TRUE)
      } else {
        dir_conc <- NA_real_
      }

      corr_rows[[st]] <- data.table(
        subtype = st,
        n_shared = length(shared),
        spearman_rho = rho,
        n_cosig = n_cosig,
        direction_concordance = dir_conc
      )
      message("  Spearman rho vs dream: ", round(rho, 3),
              " (", length(shared), " shared genes)")
      if (!is.na(dir_conc)) {
        message("  Direction concordance: ", round(dir_conc * 100, 1),
                "% among ", n_cosig, " co-significant genes")
      }
    } else {
      message("  Skipping correlation: only ", length(shared), " shared genes")
    }
  }, error = function(e) {
    message("  ERROR in limma-voom: ", e$message)
  })
}

# ---------------------------------------------------------------------------
# 5. Save correlation summary
# ---------------------------------------------------------------------------
if (length(corr_rows) > 0) {
  corr_dt <- rbindlist(corr_rows)
  corr_file <- file.path(OUT_DIR, "subtype_bulk_correlation.csv")
  fwrite(corr_dt, corr_file)
  message("\nBulk correlation summary:")
  print(corr_dt)
  message("\nSaved: ", corr_file)
} else {
  message("\nNo subtypes produced valid DE results")
}

message("\n310a complete. Results in: ", OUT_DIR)
