#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(NMF); library(data.table)})
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
nmf_obj <- readRDS(file.path(BASE, "RNA-seq/results/subtypes/nmf_results_cache_clean.rds"))
nmf_k6 <- nmf_obj$nmf_results[["6"]]
W <- basis(nmf_k6)
gencode <- fread(cmd = paste("zcat", file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")))
gencode[, ens_clean := sub("\\..*", "", gene_id)]
ens2sym <- setNames(gencode$gene_name, gencode$ens_clean)
W_sym <- ens2sym[sub("\\..*", "", rownames(W))]
keep <- !is.na(W_sym)
W <- W[keep, ]; rownames(W) <- W_sym[keep]

hs_mg <- fread(file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/hotspot_modules/hepatocytes/module_genes.tsv"))
hep24_top50 <- hs_mg[module == 24][order(-weight)][1:50, gene]

cat("=== Rebuilt NMF col 1 top-15 (symbols) ===\n")
ord1 <- order(W[, 1], decreasing = TRUE)
print(data.table(gene = rownames(W)[ord1[1:15]], weight = round(W[ord1[1:15], 1], 2)))

for (j in 1:6) {
  ord <- order(W[, j], decreasing = TRUE)
  top50 <- rownames(W)[ord[1:50]]
  inter <- intersect(hep24_top50, top50)
  cat(sprintf("\n=== NMF col %d - hep__24 overlap: %d genes ===\n", j, length(inter)))
  if (length(inter) > 0) cat("  ", paste(inter, collapse = ", "), "\n")
}
