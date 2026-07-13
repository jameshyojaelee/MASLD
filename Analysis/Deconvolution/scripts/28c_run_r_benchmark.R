#!/usr/bin/env Rscript
# 28c_run_r_benchmark.R  --  Pseudobulk benchmark, step 3/4 (MuSiC + InstaPrism).
#
# Reads the 100k atlas h5ad via rhdf5 (mirrors 05_make_sce.R), builds a genes x cells
# counts matrix + celltype/sampleID vectors, then for each fold builds a reference from
# OUT-of-fold donors and deconvolves the IN-fold donor pseudobulk (RAW COUNTS) with BOTH
# MuSiC and InstaPrism. Estimates reindexed to the 16 canonical cell types.
#
# MuSiC / InstaPrism consume RAW COUNTS (negative-binomial-aware); Rectangle (28b) uses CPM.
# Fold membership + pseudobulk come from 28a, so all methods share identical inputs.
#
# Outputs: $BENCH_OUT_DIR/music_estimates.tsv, $BENCH_OUT_DIR/instaprism_estimates.tsv
#
# Run in the `music_deconv` env via its bin/R (the env's Rscript shim is broken):
#   RNA-seq/.mamba/music_deconv/bin/R --vanilla --no-echo -f 28c_run_r_benchmark.R

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(MuSiC)
  library(InstaPrism)
})

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
H5AD <- file.path(ROOT, "Analysis/Deconvolution/reference/reference_human_scalesc_100k.h5ad")
OUT_DIR <- Sys.getenv("BENCH_OUT_DIR",
                      file.path(ROOT, "Analysis/Deconvolution/rectangle_comparison/pseudobulk"))

canonical <- readLines(file.path(OUT_DIR, "canonical_celltypes.txt"))
canonical <- canonical[nzchar(canonical)]

# ---------------------------------------------------------------- read atlas via rhdf5
cat("[28c] reading", H5AD, "\n")
x_data    <- h5read(H5AD, "/X/data")
x_indices <- h5read(H5AD, "/X/indices")
x_indptr  <- h5read(H5AD, "/X/indptr")
n_cells <- length(x_indptr) - 1L
n_genes <- as.integer(max(x_indices) + 1L)
# h5ad CSR (cells x genes) shares structure with CSC (genes x cells) -> build directly
counts <- new("dgCMatrix",
              i = as.integer(x_indices),
              p = as.integer(x_indptr),
              x = as.numeric(x_data),
              Dim = c(n_genes, n_cells))
rm(x_data, x_indices, x_indptr); gc()

read_categorical <- function(path) {
  cats  <- h5read(H5AD, paste0(path, "/categories"))
  codes <- h5read(H5AD, paste0(path, "/codes"))
  cats[codes + 1L]
}
genes     <- as.character(h5read(H5AD, "/var/_index"))
cell_type <- as.character(read_categorical("/obs/cell_type"))
donor     <- as.character(read_categorical("/obs/sample"))
rownames(counts) <- genes
colnames(counts) <- paste0("cell", seq_len(n_cells))
cat("[28c] counts:", nrow(counts), "genes x", ncol(counts), "cells\n")

# ---------------------------------------------------------------- benchmark inputs
folds <- read.delim(file.path(OUT_DIR, "donor_folds.tsv"), colClasses = c("character", "integer", "integer"))
pb <- read.delim(file.path(OUT_DIR, "pseudobulk_counts.tsv"), row.names = 1, check.names = FALSE)
pb <- as.matrix(pb)                                    # genes x donors (raw counts)
storage.mode(pb) <- "numeric"
fold_ids <- sort(unique(folds$fold))

# reindex an estimate matrix (samples x types) to samples x canonical (0 fill)
reindex_canon <- function(mat) {
  out <- matrix(0.0, nrow = nrow(mat), ncol = length(canonical),
                dimnames = list(rownames(mat), canonical))
  shared <- intersect(colnames(mat), canonical)
  out[, shared] <- mat[, shared, drop = FALSE]
  out
}

music_list <- list()
insta_list <- list()

for (fold in fold_ids) {
  ref_donors  <- folds$donor[folds$fold != fold]
  test_donors <- folds$donor[folds$fold == fold]
  ref_cells   <- which(donor %in% ref_donors)
  cat(sprintf("[28c] fold %d: ref cells=%d (donors=%d, types=%d); test donors=%d\n",
              fold, length(ref_cells), length(ref_donors),
              length(unique(cell_type[ref_cells])), length(test_donors)))

  ref_counts <- counts[, ref_cells, drop = FALSE]      # genes x cells
  ref_ct     <- cell_type[ref_cells]
  bulk       <- pb[, test_donors, drop = FALSE]         # genes x samples

  # ---- MuSiC ----
  sce <- SingleCellExperiment(
    assays  = list(counts = ref_counts),
    colData = DataFrame(celltype = ref_ct,
                        sampleID = donor[ref_cells],
                        row.names = colnames(ref_counts))
  )
  common <- intersect(rownames(bulk), rownames(sce))
  mres <- MuSiC::music_prop(
    bulk.mtx = bulk[common, , drop = FALSE],
    sc.sce   = sce[common, , drop = FALSE],
    clusters = "celltype",
    samples  = "sampleID",
    verbose  = FALSE
  )
  music_list[[as.character(fold)]] <- reindex_canon(as.matrix(mres$Est.prop.weighted))
  cat(sprintf("[28c] fold %d: MuSiC done\n", fold))

  # ---- InstaPrism ----
  ref_obj <- refPrepare(
    sc_Expr          = as.matrix(ref_counts),           # genes x cells
    cell.type.labels = as.character(ref_ct),
    cell.state.labels = as.character(ref_ct)
  )
  ires <- InstaPrism(bulk_Expr = bulk, refPhi_cs = ref_obj)   # bulk genes x samples
  iprop <- t(ires@Post.ini.ct@theta)                    # samples x types
  insta_list[[as.character(fold)]] <- reindex_canon(iprop)
  cat(sprintf("[28c] fold %d: InstaPrism done\n", fold))
  rm(ref_counts, sce, ref_obj, ires); gc()
}

write_est <- function(lst, fname) {
  mat <- do.call(rbind, lst)
  mat <- mat[folds$donor, , drop = FALSE]               # original donor order
  df <- data.frame(donor = rownames(mat), mat, check.names = FALSE)
  write.table(df, file.path(OUT_DIR, fname), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("[28c] wrote", fname, nrow(mat), "x", ncol(mat), "\n")
}

write_est(music_list, "music_estimates.tsv")
write_est(insta_list, "instaprism_estimates.tsv")
cat("[28c] DONE\n")
