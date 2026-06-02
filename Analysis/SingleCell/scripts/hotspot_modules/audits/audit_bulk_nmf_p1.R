#!/usr/bin/env Rscript
# Audit: does Hotspot hep__24 (HKDC1 stress) actually match bulk NMF k=6 P1?
# all_modules.tsv claims best_match_panel=bulk_nmf_k6, best_match_program=P1
# with jaccard=0.041 and cosine=1.0 (now known to be intersection artifact).
suppressPackageStartupMessages({ library(NMF); library(data.table) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT <- file.path(BASE, "Analysis/SingleCell/scripts/hotspot_modules/audits")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Load bulk NMF k=6 cache + program labels
nmf_obj <- readRDS(file.path(BASE, "RNA-seq/results/subtypes/nmf_results_cache_clean.rds"))
# Cache is a list of NMF fits indexed by k (3..8). k=6 is canonical.
k_target <- 6
nmf_k6 <- nmf_obj$nmf_results[[as.character(k_target)]]
W <- basis(nmf_k6)
H <- coef(nmf_k6)

# Convert Ensembl -> symbols using gencode v49 metadata
gencode <- fread(cmd = paste("zcat", file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")))
# Try canonical column names
gene_id_col <- intersect(c("gene_id", "ensembl_id"), names(gencode))[1]
symbol_col  <- intersect(c("gene_name", "symbol", "human_symbol"), names(gencode))[1]
cat(sprintf("gencode cols: id=%s symbol=%s\n", gene_id_col, symbol_col))
gencode[, ens_clean := sub("\\..*", "", get(gene_id_col))]
ens2sym <- setNames(gencode[[symbol_col]], gencode$ens_clean)

W_ens <- sub("\\..*", "", rownames(W))
W_sym <- ens2sym[W_ens]
n_mapped <- sum(!is.na(W_sym))
cat(sprintf("Ensembl -> symbol: %d / %d mapped (%.1f%%)\n",
            n_mapped, length(W_ens), 100 * n_mapped / length(W_ens)))
# Keep only mapped rows
keep <- !is.na(W_sym)
W <- W[keep, , drop = FALSE]
rownames(W) <- W_sym[keep]
cat("Bulk NMF basis matrix W:", nrow(W), "genes x", ncol(W), "programs\n")

# Try to find human-readable program labels
labels_path <- file.path(BASE, "RNA-seq/results/subtypes/program_labels_protonly.csv")
labels <- if (file.exists(labels_path)) {
  fread(labels_path)
} else {
  data.table(program = paste0("P", seq_len(ncol(W))),
             label = paste0("P", seq_len(ncol(W))))
}
cat("Program labels available:", paste(names(labels), collapse=", "), "\n")
print(labels)

# Top-50 genes per NMF program
ng_per_prog <- 50
nmf_top <- data.table()
for (j in seq_len(ncol(W))) {
  ord <- order(W[, j], decreasing = TRUE)
  top <- head(rownames(W)[ord], ng_per_prog)
  nmf_top <- rbind(nmf_top,
    data.table(program = paste0("P", j), rank = seq_along(top),
               gene = top, weight = W[ord[seq_len(ng_per_prog)], j]))
}
fwrite(nmf_top, file.path(OUT, "bulk_nmf_k6_top50.tsv"), sep = "\t")
cat("\nWrote bulk_nmf_k6_top50.tsv\n")

# Where does HKDC1 sit in each program?
cat("\n=== HKDC1 placement across 6 bulk NMF programs ===\n")
hkdc1_rows <- which(rownames(W) == "HKDC1")
if (length(hkdc1_rows) == 0) {
  cat("HKDC1 not in W rownames. Top of rownames:\n")
  print(head(rownames(W), 20))
} else {
  for (j in seq_len(ncol(W))) {
    rank_in_p <- rank(-W[, j])[hkdc1_rows]
    cat(sprintf("  P%d: HKDC1 weight=%.4f, rank=%d/%d\n",
                j, W[hkdc1_rows, j], rank_in_p, nrow(W)))
  }
}

# Hep__24 top-50 from Hotspot
hs_mg <- fread(file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/hotspot_modules/hepatocytes/module_genes.tsv"))
hep24_top50 <- hs_mg[module == 24][order(-weight)][1:50, gene]
cat("\nHotspot hep__24 top-50 (head):", paste(head(hep24_top50, 10), collapse=", "), "...\n")

# Compute clean Jaccard@50 + hypergeometric per program
cat("\n=== hep__24 top-50 vs each NMF program top-50 ===\n")
N_universe <- length(unique(c(rownames(W), hs_mg$gene)))
cat(sprintf("Gene universe = %d (union of NMF rows + Hotspot genes)\n", N_universe))
match_tbl <- data.table()
for (j in seq_len(ncol(W))) {
  p_top50 <- nmf_top[program == paste0("P", j), gene]
  inter <- length(intersect(hep24_top50, p_top50))
  union_ <- length(union(hep24_top50, p_top50))
  jacc <- inter / union_
  # hypergeometric
  p_hyper <- phyper(inter - 1, length(p_top50),
                    N_universe - length(p_top50), length(hep24_top50),
                    lower.tail = FALSE)
  match_tbl <- rbind(match_tbl, data.table(
    nmf_program = paste0("P", j), intersect_count = inter,
    jaccard_50 = round(jacc, 4),
    hyper_p = signif(p_hyper, 3)))
}
print(match_tbl[order(-jaccard_50)])
fwrite(match_tbl, file.path(OUT, "hep24_vs_bulk_nmf_match.tsv"), sep = "\t")
cat("\nBest match by Jaccard:", match_tbl[order(-jaccard_50)][1, nmf_program], "\n")
