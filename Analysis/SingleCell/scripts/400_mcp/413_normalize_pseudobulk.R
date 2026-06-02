#!/usr/bin/env Rscript
# 413_normalize_pseudobulk.R
# TMM normalize + log2(CPM+1) per cell type independently (DIALOGUE expects per-CT variance).
# Apply same confounder strip as cNMF (already done upstream but re-verify).
# Intersect gene universe across cell types (genes with CPM>1 in >=5% cells in >=2 cell types).
# Write per-cell-type log-CPM matrices + final gene universe table.

suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
})

root <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
in_dir <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/inputs/dialogue_pseudobulk")
out_dir <- in_dir
stopifnot(dir.exists(in_dir))

celltypes <- c("hepatocytes", "endothelial_cells", "fibroblasts", "macrophages", "cholangiocytes")

# --- Pass 1: detection filter per cell type ---
detected <- list()
raw_cts  <- list()
for (ct in celltypes) {
  f <- file.path(in_dir, sprintf("%s_counts.tsv.gz", ct))
  cat(sprintf("[413] Load %s\n", basename(f)))
  m <- data.table::fread(f, data.table = FALSE)
  rn <- m[[1]]; m[[1]] <- NULL
  rownames(m) <- rn
  raw_cts[[ct]] <- as.matrix(m)
  # CPM filter: >=1 CPM in >=5% of samples
  cpm_m <- edgeR::cpm(raw_cts[[ct]])
  frac_detected <- rowMeans(cpm_m >= 1)
  detected[[ct]] <- names(frac_detected[frac_detected >= 0.05])
}

# --- Gene universe: detected in >=2 cell types ---
all_genes <- unique(unlist(detected))
det_mat <- sapply(detected, function(g) as.integer(all_genes %in% g))
rownames(det_mat) <- all_genes
gene_universe <- rownames(det_mat)[rowSums(det_mat) >= 2]
cat(sprintf("[413] Gene universe (detected in >=2 cell types): %d\n", length(gene_universe)))

write.table(
  data.frame(gene_name = gene_universe),
  file.path(out_dir, "gene_universe.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# --- Pass 2: TMM + log2(CPM+1) restricted to gene universe ---
for (ct in celltypes) {
  m <- raw_cts[[ct]]
  m <- m[rownames(m) %in% gene_universe, , drop = FALSE]
  dge <- edgeR::DGEList(counts = m)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1) # log2(CPM + prior)
  out <- file.path(out_dir, sprintf("%s_logcpm.tsv.gz", ct))
  cat(sprintf("[413] %s -> %s (%d x %d)\n", ct, basename(out), nrow(logcpm), ncol(logcpm)))
  dt <- data.table::as.data.table(logcpm, keep.rownames = "gene_name")
  data.table::fwrite(dt, out, sep = "\t", compress = "gzip")
}

cat("[413] DONE.\n")
