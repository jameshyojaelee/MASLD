#!/usr/bin/env Rscript
# Check DIAMOND DEG replication across existing 5 diet models
suppressPackageStartupMessages({library(data.table)})
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

# Load dream pooled
dream <- fread(file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/meta_analysis/dream_pooled_results.csv"))
dream[, gene_id_bare := gsub("[.][0-9]+$", "", gene)]

# Our MASH DE
de <- fread(file.path(BASE, "RNA-seq/Mouse/Western_Diet_Datasets/GSE220575/results/de_mash_vs_control_annotated.csv"))

# Per-diet results
per_diet_dir <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")
diets <- c("MCD", "HFD", "CDAHFD", "FPC", "LIDPAD")
per_diet_up <- list()
for (d in diets) {
  f <- file.path(per_diet_dir, paste0(d, "_de_results.csv"))
  if (file.exists(f)) {
    dd <- fread(f)
    if ("gene" %in% names(dd)) dd[, gene_id_bare := gsub("[.][0-9]+$", "", gene)]
    else if ("gene_id" %in% names(dd)) dd[, gene_id_bare := gsub("[.][0-9]+$", "", gene_id)]
    per_diet_up[[d]] <- dd[adj.P.Val < 0.05 & logFC > 0, gene_id_bare]
  }
}

# Our sig up genes
our_up <- de[adj.P.Val < 0.05 & logFC > 0, gene_id]

# Count replication
n_diets_rep <- sapply(our_up, function(g) {
  sum(sapply(per_diet_up, function(d_genes) g %in% d_genes))
})

cat("=== DIAMOND MASH-up DEG Replication Across Existing 5 Diets ===\n")
cat("Total DIAMOND MASH-up DEGs:", length(our_up), "\n\n")
tab <- table(factor(n_diets_rep, levels = 0:5))
for (n in 0:5) {
  ct <- tab[as.character(n)]
  cat(sprintf("  %d diets: %d genes (%0.1f%%)\n", n, ct, ct / length(our_up) * 100))
}

# 4+ replication
rep4plus <- our_up[n_diets_rep >= 4]
cat("\nGenes up in DIAMOND + 4+ existing diets:", length(rep4plus), "\n")
de_rep4 <- de[gene_id %in% rep4plus & adj.P.Val < 0.05 & logFC > 0][order(-logFC)]
if (nrow(de_rep4) > 0) {
  cat("\nTop 20 replicated in 4+ diets:\n")
  print(de_rep4[1:min(20, nrow(de_rep4)),
    .(symbol = mouse_symbol_gtf, logFC = round(logFC, 2), padj = signif(adj.P.Val, 3))])
}

# Novel: DIAMOND-only
novel <- our_up[n_diets_rep == 0]
de_novel <- de[gene_id %in% novel & adj.P.Val < 0.05 & logFC > 0][order(-logFC)]
cat("\n\nNovel DIAMOND-specific up-DEGs:", length(novel), "\n")
if (nrow(de_novel) > 0) {
  cat("Top 20 novel:\n")
  print(de_novel[1:min(20, nrow(de_novel)),
    .(symbol = mouse_symbol_gtf, logFC = round(logFC, 2),
      padj = signif(adj.P.Val, 3), biotype = mouse_biotype)])
}

# 5/5 + DIAMOND = 6-model concordant
rep5 <- our_up[n_diets_rep == 5]
de_rep5 <- de[gene_id %in% rep5 & adj.P.Val < 0.05 & logFC > 0][order(-logFC)]
cat("\n\n=== 6-Model Concordant (all 5 existing + DIAMOND) ===\n")
cat("Count:", length(rep5), "\n")
if (nrow(de_rep5) > 0) {
  cat("Top 30:\n")
  print(de_rep5[1:min(30, nrow(de_rep5)),
    .(symbol = mouse_symbol_gtf, logFC = round(logFC, 2), padj = signif(adj.P.Val, 3))])
}
