#!/usr/bin/env Rscript
#' 304: Monocle 3 Trajectory + tradeSeq DE.
#'
#' For each of 5 core cell types:
#'   1. Load preprocessed data (exported from h5ad by helper)
#'   2. Build Monocle 3 CDS, learn graph, order cells
#'   3. Run tradeSeq: startVsEndTest, patternTest, conditionTest
#'
#' Inputs:
#'   - Cell-type subsets from 300: pseudotime/{celltype}_subset.h5ad
#'     (read via anndata R bridge or pre-exported RDS)
#'
#' Outputs (to results_gpu_v2/pseudotime/):
#'   - monocle3_pseudotime_{celltype}.csv
#'   - monocle3_branch_points_{celltype}.csv
#'   - tradeseq_startend_{celltype}.csv
#'   - tradeseq_condition_{celltype}.csv
#'
#' Usage:
#'   sbatch run_pseudotime_pipeline.sh

suppressPackageStartupMessages({
  library(monocle3)
  library(tradeSeq)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(BiocParallel)
})

# Register parallel backend for tradeSeq (uses SLURM-allocated CPUs)
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
register(MulticoreParam(workers = n_cores))
cat("Registered BiocParallel backend:", n_cores, "cores\n")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SC_DIR <- file.path(BASE, "Analysis/SingleCell")
RESULTS <- file.path(SC_DIR, "results_gpu_v2")
PT_DIR <- file.path(RESULTS, "pseudotime")

CELL_TYPES <- c("Hepatocytes", "Macrophages", "Fibroblasts",
                "Endothelial_cells", "Cholangiocytes")

# Configure reticulate to use rapids_singlecell Python (has scanpy/anndata)
rapids_python <- "/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python3"
if (file.exists(rapids_python)) {
  Sys.setenv(RETICULATE_PYTHON = rapids_python)
}

cat("==========================================================\n")
cat("304: Monocle 3 + tradeSeq\n")
cat("==========================================================\n")

# ---------------------------------------------------------------------------
# Helper: load h5ad via reticulate + metadata CSV from Script 300
# ---------------------------------------------------------------------------
load_subset_as_cds <- function(ct_name) {
  cat("Loading h5ad for:", ct_name, "\n")

  h5ad_path <- file.path(PT_DIR, paste0(ct_name, "_subset.h5ad"))
  meta_path <- file.path(PT_DIR, paste0(ct_name, "_metadata.csv"))

  if (!file.exists(h5ad_path)) {
    stop("h5ad not found: ", h5ad_path, ". Run Script 300 first.")
  }

  # Use reticulate to read h5ad
  library(reticulate)
  sc <- import("scanpy")
  sp <- import("scipy.sparse")

  cat("Reading h5ad via reticulate...\n")
  adata <- sc$read_h5ad(h5ad_path)
  np <- import("numpy")

  # Subset to HVGs in Python BEFORE converting to R (avoids 2^31 sparse limit)
  hvg_mask <- py_to_r(adata$var[["highly_variable"]])
  if (is.null(hvg_mask) || all(is.na(hvg_mask))) {
    cat("No HVG column found — using top 5000 variable genes\n")
    X_full <- adata$X
    if (sp$issparse(X_full)) {
      gene_var <- py_to_r(np$array(X_full$power(2L)$mean(axis = 0L))$flatten()) -
                  py_to_r(np$array(X_full$mean(axis = 0L))$flatten())^2
    } else {
      gene_var <- py_to_r(np$var(X_full, axis = 0L))
    }
    top_idx <- order(gene_var, decreasing = TRUE)[1:min(5000, length(gene_var))]
    adata_sub <- adata[, as.integer(top_idx - 1L)]  # Python 0-indexed
  } else {
    hvg_idx <- which(as.logical(hvg_mask))
    cat("Subsetting to", length(hvg_idx), "HVGs in Python...\n")
    adata_sub <- adata[, as.integer(hvg_idx - 1L)]
  }

  # Now extract the smaller matrix (cells x HVGs)
  X <- adata_sub$X
  if (sp$issparse(X)) {
    X_coo <- sp$coo_matrix(X)
    X_r <- sparseMatrix(
      i = as.integer(py_to_r(X_coo$row)) + 1L,
      j = as.integer(py_to_r(X_coo$col)) + 1L,
      x = as.numeric(py_to_r(X_coo$data)),
      dims = as.integer(py_to_r(adata_sub$shape))
    )
  } else {
    X_r <- as(as.matrix(py_to_r(X)), "dgCMatrix")
  }

  var_names <- as.character(py_to_r(adata_sub$var_names$tolist()))
  cell_names <- as.character(py_to_r(adata$obs_names$tolist()))

  # Free Python objects
  rm(adata_sub, X, X_coo)
  gc()

  # Transpose to genes x cells (Monocle format)
  X_r <- t(X_r)
  rownames(X_r) <- var_names
  colnames(X_r) <- cell_names

  cat("Matrix shape:", nrow(X_r), "genes x", ncol(X_r), "cells\n")

  # Cell metadata from pre-exported CSV (more reliable than reticulate DataFrame)
  if (file.exists(meta_path)) {
    cell_metadata <- fread(meta_path, data.table = FALSE)
    rownames(cell_metadata) <- cell_metadata[[1]]
  } else {
    cell_metadata <- as.data.frame(py_to_r(adata$obs))
    rownames(cell_metadata) <- cell_names
  }

  # Gene metadata
  gene_metadata <- data.frame(
    gene_short_name = var_names,
    row.names = var_names
  )

  # Extract UMAP
  umap <- NULL
  if ("UMAP_1" %in% colnames(cell_metadata) && "UMAP_2" %in% colnames(cell_metadata)) {
    umap <- as.matrix(cell_metadata[, c("UMAP_1", "UMAP_2")])
    rownames(umap) <- cell_names
  } else if (!is.null(adata$obsm) && "X_umap" %in% names(adata$obsm)) {
    umap <- as.matrix(py_to_r(adata$obsm[["X_umap"]]))
    rownames(umap) <- cell_names
  }

  # Build CDS
  cat("Building CDS object...\n")
  cds <- new_cell_data_set(
    expression_data = X_r,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )

  # Transfer UMAP
  if (!is.null(umap)) {
    reducedDims(cds)[["UMAP"]] <- umap
  }

  return(cds)
}


# ---------------------------------------------------------------------------
# Process each cell type
# ---------------------------------------------------------------------------
for (ct in CELL_TYPES) {
  cat("\n==================================================\n")
  cat("Processing:", ct, "\n")
  cat("==================================================\n")

  # Check for h5ad subset from Script 300
  h5ad_path <- file.path(PT_DIR, paste0(ct, "_subset.h5ad"))
  if (!file.exists(h5ad_path)) {
    cat("h5ad subset not found for", ct, "— skipping\n")
    cat("  Expected:", h5ad_path, "\n")
    next
  }

  # Load as CDS from MatrixMarket + CSV
  cds <- NULL
  tryCatch({
    cds <- load_subset_as_cds(ct)
    cat("CDS shape:", dim(cds), "\n")
  }, error = function(e) {
    cat("Failed to load data:", conditionMessage(e), "\n")
  })

  if (is.null(cds)) next

  n_cells <- ncol(cds)
  cat("Working with", n_cells, "cells\n")

  # ===================================================================
  # Monocle 3 trajectory
  # ===================================================================
  cat("--- Monocle 3: Learn graph ---\n")

  # Preprocess if PCA not transferred
  if (is.null(reducedDims(cds)[["PCA"]])) {
    cds <- preprocess_cds(cds, num_dim = 50)
  }

  # Reduce dimensions if UMAP not transferred
  if (is.null(reducedDims(cds)[["UMAP"]])) {
    cds <- reduce_dimension(cds)
  }

  # Cluster cells (needed for UMAP graph)
  cds <- cluster_cells(cds, resolution = 1e-3)

  # Learn graph
  cds <- learn_graph(cds, use_partition = FALSE)

  # Order cells: root = Healthy condition cells
  condition <- colData(cds)$condition
  healthy_cells <- colnames(cds)[condition == "Healthy"]

  if (length(healthy_cells) > 0) {
    # Find the graph node closest to Healthy cells
    cds <- order_cells(cds, root_cells = healthy_cells[1:min(10, length(healthy_cells))])
  } else {
    # Interactive fallback — use first cell
    cat("WARNING: No Healthy cells found — using first cell as root\n")
    cds <- order_cells(cds, root_cells = colnames(cds)[1])
  }

  # Extract pseudotime
  pt <- pseudotime(cds)
  cat("Monocle 3 pseudotime: median =", median(pt, na.rm = TRUE),
      ", range = [", min(pt, na.rm = TRUE), ",", max(pt, na.rm = TRUE), "]\n")

  # Save pseudotime
  pt_df <- data.frame(
    cell = names(pt),
    monocle3_pseudotime = as.numeric(pt),
    condition = as.character(colData(cds)$condition),
    stringsAsFactors = FALSE
  )
  out_path <- file.path(PT_DIR, paste0("monocle3_pseudotime_", ct, ".csv"))
  fwrite(pt_df, out_path)
  cat("Saved Monocle 3 pseudotime to", out_path, "\n")

  # Branch points
  graph <- principal_graph(cds)[["UMAP"]]
  if (!is.null(graph)) {
    bp <- which(igraph::degree(graph) > 2)
    if (length(bp) > 0) {
      bp_coords <- reducedDims(cds)[["UMAP"]]  # approximate
      bp_df <- data.frame(
        branch_point = seq_along(bp),
        node = bp,
        stringsAsFactors = FALSE
      )
      out_path <- file.path(PT_DIR, paste0("monocle3_branch_points_", ct, ".csv"))
      fwrite(bp_df, out_path)
      cat("Found", length(bp), "branch points\n")
    } else {
      cat("No branch points found (linear trajectory)\n")
    }
  }

  # ===================================================================
  # tradeSeq DE analysis
  # ===================================================================
  cat("\n--- tradeSeq: Differential expression along trajectory ---\n")

  # Fit GAMs (use top 5000 variable genes for tractability)
  counts_mat <- counts(cds)
  gene_var <- apply(counts_mat, 1, var)
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(5000, nrow(counts_mat))]
  counts_sub <- counts_mat[top_genes, ]

  # Get pseudotime and cell weights from Monocle 3
  pt_numeric <- as.numeric(pseudotime(cds))
  valid <- is.finite(pt_numeric) & pt_numeric < Inf

  if (sum(valid) < 100) {
    cat("Too few valid pseudotime cells (", sum(valid), ") — skipping tradeSeq\n")
    rm(cds)
    gc()
    next
  }

  cat("Fitting GAMs on", length(top_genes), "genes,", sum(valid), "cells...\n")

  # tradeSeq expects a pseudotime matrix (cells x lineages)
  # For single lineage, it's a single column
  pt_mat <- matrix(pt_numeric[valid], ncol = 1)
  cellWeights <- matrix(1, nrow = sum(valid), ncol = 1)
  counts_valid <- counts_sub[, valid]

  tryCatch({
    sce <- fitGAM(
      counts = as.matrix(counts_valid),
      pseudotime = pt_mat,
      cellWeights = cellWeights,
      nknots = 6,
      verbose = TRUE,
      parallel = TRUE
    )

    # startVsEndTest: genes changing from start to end
    cat("Running startVsEndTest...\n")
    se_res <- startVsEndTest(sce)
    se_res$gene <- rownames(se_res)
    se_res <- se_res[order(se_res$pvalue), ]
    out_path <- file.path(PT_DIR, paste0("tradeseq_startend_", ct, ".csv"))
    fwrite(as.data.frame(se_res), out_path)
    cat("startVsEndTest: ", sum(se_res$pvalue < 0.05, na.rm = TRUE),
        " significant genes (p<0.05)\n")

    # conditionTest: genes DE between conditions at matched pseudotime
    # This requires conditions as a factor in the GAM
    cat("Running conditionTest...\n")
    cond_valid <- colData(cds)$condition[valid]
    cond_factor <- factor(cond_valid)

    if (nlevels(cond_factor) >= 2) {
      tryCatch({
        # Re-fit with conditions
        sce_cond <- fitGAM(
          counts = as.matrix(counts_valid),
          pseudotime = pt_mat,
          cellWeights = cellWeights,
          nknots = 6,
          conditions = cond_factor,
          verbose = TRUE,
          parallel = TRUE
        )

        cond_res <- conditionTest(sce_cond)
        cond_res$gene <- rownames(cond_res)
        cond_res <- cond_res[order(cond_res$pvalue), ]
        out_path <- file.path(PT_DIR, paste0("tradeseq_condition_", ct, ".csv"))
        fwrite(as.data.frame(cond_res), out_path)
        cat("conditionTest: ", sum(cond_res$pvalue < 0.05, na.rm = TRUE),
            " significant genes (p<0.05)\n")
      }, error = function(e) {
        cat("conditionTest failed:", conditionMessage(e), "\n")
      })
    }

  }, error = function(e) {
    cat("tradeSeq failed:", conditionMessage(e), "\n")
  })

  # Compare with Palantir pseudotime if available
  palantir_path <- file.path(PT_DIR, paste0("palantir_pseudotime_", ct, ".csv"))
  if (file.exists(palantir_path)) {
    pr_df <- fread(palantir_path)
    common <- intersect(pt_df$cell, pr_df$V1)
    if (length(common) > 30) {
      m3_pt <- pt_df$monocle3_pseudotime[match(common, pt_df$cell)]
      pr_pt <- pr_df$palantir_pseudotime[match(common, pr_df$V1)]
      valid_both <- is.finite(m3_pt) & is.finite(pr_pt)
      if (sum(valid_both) > 30) {
        rho <- cor(m3_pt[valid_both], pr_pt[valid_both], method = "spearman")
        cat("Monocle 3 vs Palantir correlation: rho =", round(rho, 3), "\n")
      }
    }
  }

  rm(cds)
  gc()
}

cat("\n=== 304: Monocle 3 + tradeSeq COMPLETE ===\n")
