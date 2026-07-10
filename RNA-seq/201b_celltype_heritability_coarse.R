#!/usr/bin/env Rscript
# 201b_celltype_heritability_coarse.R
# COARSE 7-lineage version of the per-cell-type GWAS-COLOC marker enrichment
# (the x-axis of the Fig 3E genetics<->expression disconnect panel).
# Replicates 201_celltype_heritability.R's COLOC enrichment EXACTLY, but on the
# 7 coarse-lineage pseudobulk DE (pseudobulk_de_coarse/) instead of the 11 fine.
#
# Method (identical to 201): per lineage, take the top-10% most disease-specific
# genes (specificity = |logFC| x -log10(padj)); test whether their canonical
# gene-level COLOC PP4 is enriched vs all other genes (one-sided Wilcoxon);
# fold = mean_in / mean_out. Canonical coloc = gene_level_coloc.csv.

suppressPackageStartupMessages({ library(data.table) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

# --- gene-level COLOC PP4 (canonical) ---
coloc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc <- coloc[gene != "" & !is.na(gene)]
gwas_genes <- coloc[, .(gene, coloc_pp4 = coloc_best_pp4)]

# --- coarse-lineage pseudobulk DE ---
sc_dir   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de_coarse")
sc_files <- list.files(sc_dir, pattern = "_de\\.csv$", full.names = TRUE)
sc_de <- rbindlist(lapply(sc_files, function(f) {
  dt <- fread(f); dt[, .(gene, cell_type, logFC, padj)]
}), fill = TRUE)
cell_types <- unique(sc_de$cell_type)
cat("Lineages:", paste(cell_types, collapse = ", "), "\n")

# --- per-lineage top-10% disease-specific gene sets (as in 201) ---
ct_gene_sets <- lapply(cell_types, function(ct) {
  d <- sc_de[cell_type == ct & !is.na(padj)]
  d[, spec := abs(logFC) * pmin(-log10(padj + 1e-300), 300)]
  d[order(-spec)][seq_len(min(nrow(d), round(nrow(d) * 0.1)))]$gene
})
names(ct_gene_sets) <- cell_types

# --- COLOC enrichment per lineage (Wilcoxon greater; fold in/out) ---
res <- rbindlist(lapply(cell_types, function(ct) {
  in_set <- gwas_genes$gene %in% ct_gene_sets[[ct]]
  if (sum(in_set) < 5) return(NULL)
  wt <- wilcox.test(gwas_genes$coloc_pp4[in_set], gwas_genes$coloc_pp4[!in_set],
                    alternative = "greater")
  mi <- mean(gwas_genes$coloc_pp4[in_set], na.rm = TRUE)
  mo <- mean(gwas_genes$coloc_pp4[!in_set], na.rm = TRUE)
  data.table(cell_type = ct, n_genes_in_set = sum(in_set),
             mean_coloc_in = mi, mean_coloc_out = mo,
             fold_enrichment_coloc = mi / max(mo, 1e-10),
             wilcox_p_coloc = wt$p.value)
}), fill = TRUE)
res[, fdr_coloc := p.adjust(wilcox_p_coloc, method = "BH")]
setorder(res, -fold_enrichment_coloc)

out <- file.path(BASE, "RNA-seq/results/gwas_rna_integration/celltype_heritability_coarse_results.csv")
fwrite(res, out)
cat("\nSaved:", out, "\n")
print(res[, .(cell_type, fold_enrichment_coloc = round(fold_enrichment_coloc, 3),
              wilcox_p_coloc = signif(wilcox_p_coloc, 3), fdr_coloc = signif(fdr_coloc, 3))])
