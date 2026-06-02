#!/usr/bin/env Rscript
# 201_celltype_heritability.R — Cell-Type Heritability Enrichment
#
# Tests which scRNA cell types are enriched for MASLD GWAS signal,
# using gene-level COLOC/TWAS scores stratified by cell-type specificity.
#
# Approach (equivalent to MAGMA.Celltyping concept):
#   1. Gene-level GWAS scores: COLOC PP.H4, TWAS z-scores (already available)
#   2. Cell-type specificity: tau index from scRNA pseudobulk DE
#   3. Enrichment: linear model + competitive rank-sum test
#
# This avoids the MAGMA binary dependency while testing the same hypothesis:
# "Are genes specifically expressed in cell type X enriched for GWAS signal?"
#
# Inputs:
#   - COLOC gene-level: GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv
#   - COLOC per-GWAS: GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv
#   - TWAS results: RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv
#   - scRNA pseudobulk DE: Analysis/SingleCell/results_gpu_v2/{CellType}_de.csv
#   - Dream DEGs: RNA-seq/Human/.../results/integration/dream_results.csv
#
# Outputs:
#   - RNA-seq/results/gwas_rna_integration/celltype_heritability_results.csv
#   - RNA-seq/results/gwas_rna_integration/celltype_specificity_scores.csv
#
# SLURM: cpu partition, 8 CPUs, 64GB, 48h

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

outdir <- file.path(BASE, "RNA-seq/results/gwas_rna_integration")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# 1. Load gene-level GWAS scores
# ===========================================================================
cat("Loading gene-level GWAS scores...\n")

# COLOC PP.H4 per gene (best across GWAS)
coloc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc <- coloc[gene != "" & !is.na(gene)]
cat("  COLOC:", nrow(coloc), "genes\n")

# TWAS z-scores
twas <- fread(file.path(BASE,
  "RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv"))
# Best TWAS z per gene across GWAS
twas_best <- twas[, .(
  twas_z = zscore[which.min(pvalue)],
  twas_p = min(pvalue)
), by = .(gene = gene_name)]
cat("  TWAS:", nrow(twas_best), "genes\n")

# Merge into unified gene-level GWAS score table
gwas_genes <- merge(
  coloc[, .(gene, coloc_pp4 = coloc_best_pp4, coloc_n_gwas = coloc_n_gwas_h4_05)],
  twas_best,
  by = "gene", all = TRUE
)

# Composite GWAS score (normalized rank average of COLOC + TWAS)
gwas_genes[, coloc_rank := frank(coloc_pp4, na.last = "keep") / sum(!is.na(coloc_pp4))]
gwas_genes[, twas_rank := frank(abs(twas_z), na.last = "keep") / sum(!is.na(twas_z))]
gwas_genes[, gwas_score := rowMeans(cbind(coloc_rank, twas_rank), na.rm = TRUE)]

cat("  Unified gene-level scores:", nrow(gwas_genes), "genes\n")

# ===========================================================================
# 2. Compute cell-type specificity from scRNA pseudobulk DE
# ===========================================================================
cat("\nComputing cell-type specificity from scRNA...\n")

sc_dir <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")
sc_files <- list.files(sc_dir, pattern = "_de\\.csv$", full.names = TRUE)
# D4 HOLD closure 2026-04-22: exclude allcell pseudobulk (different schema, uses
# `lfc`/`ave` columns rather than the per-cell-type `logFC`/`AveExpr`). Without
# this filter the glob pulls in allcell_pseudobulk_de.csv and breaks the
# downstream logFC / AveExpr select.
sc_files <- sc_files[!grepl("^allcell_", basename(sc_files))]
cat("  Found", length(sc_files), "cell-type DE files\n")

sc_list <- lapply(sc_files, function(f) {
  dt <- fread(f)
  if (!"cell_type" %in% names(dt)) {
    ct <- gsub("_de\\.csv$", "", basename(f))
    dt[, cell_type := ct]
  }
  return(dt[, .(gene, cell_type, logFC, padj, AveExpr)])
})

sc_de <- rbindlist(sc_list, fill = TRUE)
cell_types <- unique(sc_de$cell_type)
cat("  Cell types:", paste(cell_types, collapse = ", "), "\n")
cat("  Total DE entries:", nrow(sc_de), "\n")

# Compute tau specificity index
# tau = proportion of cell types where the gene is NOT DE
# High tau = cell-type-specific expression changes
sc_sig <- sc_de[, .(is_sig = padj < 0.05 & abs(logFC) > 0.25), by = .(gene, cell_type)]
tau_dt <- sc_sig[, .(
  n_ct_sig = sum(is_sig, na.rm = TRUE),
  n_ct_tested = .N,
  tau = 1 - sum(is_sig, na.rm = TRUE) / .N
), by = gene]

# Cell-type specificity: for each gene, which cell type has the strongest effect?
sc_best <- sc_de[, .(
  best_ct = cell_type[which.max(abs(logFC))],
  best_ct_logFC = logFC[which.max(abs(logFC))],
  best_ct_padj = padj[which.max(abs(logFC))]
), by = gene]

specificity_dt <- merge(tau_dt, sc_best, by = "gene")
cat("  Genes with specificity scores:", nrow(specificity_dt), "\n")

# Create per-cell-type gene sets (top 10% most specific genes)
ct_gene_sets <- list()
for (ct in cell_types) {
  ct_de <- sc_de[cell_type == ct & !is.na(padj)]
  # Rank by |logFC| × -log10(padj)
  ct_de[, specificity_score := abs(logFC) * pmin(-log10(padj + 1e-300), 300)]
  top_genes <- ct_de[order(-specificity_score)][1:min(nrow(ct_de), round(nrow(ct_de) * 0.1))]
  ct_gene_sets[[ct]] <- top_genes$gene
}

cat("  Gene set sizes (top 10% per cell type):\n")
for (ct in names(ct_gene_sets)) {
  cat("    ", ct, ":", length(ct_gene_sets[[ct]]), "\n")
}

# ===========================================================================
# 3. Enrichment test: cell-type-specific genes × GWAS scores
# ===========================================================================
cat("\n=== Cell-type heritability enrichment ===\n")

enrichment_results <- lapply(names(ct_gene_sets), function(ct) {
  ct_genes <- ct_gene_sets[[ct]]

  # Method 1: Wilcoxon rank-sum test
  # Compare GWAS scores of cell-type-specific genes vs all others
  in_set <- gwas_genes$gene %in% ct_genes
  if (sum(in_set) < 5) return(NULL)

  # Test on COLOC PP.H4
  wt_coloc <- wilcox.test(
    gwas_genes$coloc_pp4[in_set],
    gwas_genes$coloc_pp4[!in_set],
    alternative = "greater"
  )

  # Test on |TWAS z|
  wt_twas <- wilcox.test(
    abs(gwas_genes$twas_z[in_set & !is.na(gwas_genes$twas_z)]),
    abs(gwas_genes$twas_z[!in_set & !is.na(gwas_genes$twas_z)]),
    alternative = "greater"
  )

  # Test on composite GWAS score
  wt_composite <- wilcox.test(
    gwas_genes$gwas_score[in_set & !is.na(gwas_genes$gwas_score)],
    gwas_genes$gwas_score[!in_set & !is.na(gwas_genes$gwas_score)],
    alternative = "greater"
  )

  # Method 2: Linear model (COLOC ~ specificity_score)
  ct_de_sub <- sc_de[cell_type == ct & gene %in% gwas_genes$gene]
  ct_merged <- merge(
    gwas_genes[, .(gene, coloc_pp4, gwas_score)],
    ct_de_sub[, .(gene, ct_logFC = logFC, ct_padj = padj)],
    by = "gene"
  )

  lm_res <- tryCatch({
    fit <- lm(coloc_pp4 ~ abs(ct_logFC), data = ct_merged)
    coef(summary(fit))["abs(ct_logFC)", ]
  }, error = function(e) c(Estimate = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA))

  # Effect sizes
  mean_coloc_in <- mean(gwas_genes$coloc_pp4[in_set], na.rm = TRUE)
  mean_coloc_out <- mean(gwas_genes$coloc_pp4[!in_set], na.rm = TRUE)
  fold_enrichment <- mean_coloc_in / max(mean_coloc_out, 1e-10)

  # Proportion of COLOC genes
  n_coloc_in <- sum(gwas_genes$coloc_pp4[in_set] > 0.5, na.rm = TRUE)
  n_coloc_out <- sum(gwas_genes$coloc_pp4[!in_set] > 0.5, na.rm = TRUE)
  fisher_test <- fisher.test(matrix(c(
    n_coloc_in, sum(in_set) - n_coloc_in,
    n_coloc_out, sum(!in_set) - n_coloc_out
  ), nrow = 2))

  data.table(
    cell_type = ct,
    n_genes_in_set = sum(in_set),
    n_coloc_in_set = n_coloc_in,
    mean_coloc_in = mean_coloc_in,
    mean_coloc_out = mean_coloc_out,
    fold_enrichment_coloc = fold_enrichment,
    wilcox_p_coloc = wt_coloc$p.value,
    wilcox_p_twas = wt_twas$p.value,
    wilcox_p_composite = wt_composite$p.value,
    lm_beta = lm_res["Estimate"],
    lm_p = lm_res["Pr(>|t|)"],
    fisher_or = fisher_test$estimate,
    fisher_p = fisher_test$p.value
  )
})

enrichment_dt <- rbindlist(enrichment_results, fill = TRUE)

# FDR correction
enrichment_dt[, fdr_coloc := p.adjust(wilcox_p_coloc, method = "BH")]
enrichment_dt[, fdr_twas := p.adjust(wilcox_p_twas, method = "BH")]
enrichment_dt[, fdr_composite := p.adjust(wilcox_p_composite, method = "BH")]

# D4 HOLD closure 2026-04-22: bonferroni family-wise MTC across 11 celltypes × 24 regulons
# Rationale: "HNF4A regulon enriched in hepatocytes" claim needs family-wise MTC
# across all tested cell types. Bonferroni is the strictest control; significant
# calls here are robust to multiple comparisons. BH (above) retained for power.
enrichment_dt[, wilcox_p_coloc_bonf     := p.adjust(wilcox_p_coloc,     method = "bonferroni")]
enrichment_dt[, wilcox_p_twas_bonf      := p.adjust(wilcox_p_twas,      method = "bonferroni")]
enrichment_dt[, wilcox_p_composite_bonf := p.adjust(wilcox_p_composite, method = "bonferroni")]

# Sort by enrichment
setorder(enrichment_dt, wilcox_p_composite)

cat("\n=== Results ===\n")
print(enrichment_dt[, .(cell_type, fold_enrichment_coloc,
  wilcox_p_coloc, wilcox_p_twas, fdr_composite)])

# ===========================================================================
# 4. Stratified enrichment by DEG quintiles
# ===========================================================================
cat("\n=== Stratified enrichment by DEG significance ===\n")

dream <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results.csv"))

# Dream uses ENSEMBL IDs — map to symbols using COLOC data
coloc_map <- coloc[gene != "" & ensembl != "", .(ensembl, symbol = gene)]
coloc_map[, ensembl_base := sub("\\.\\d+$", "", ensembl)]
dream[, ensembl_base := sub("\\.\\d+$", "", gene)]
dream <- merge(dream, coloc_map[, .(ensembl_base, symbol)], by = "ensembl_base", all.x = TRUE)
dream[!is.na(symbol), gene_symbol := symbol]
dream[is.na(symbol), gene_symbol := gene]

# Create quintiles by |t-statistic|
dream[, t_quintile := cut(abs(t), quantile(abs(t), probs = 0:5/5, na.rm = TRUE),
                           include.lowest = TRUE, labels = paste0("Q", 1:5))]

# For each quintile, what's the mean COLOC PP.H4?
dream_gwas <- merge(dream[, .(gene = gene_symbol, t_quintile, logFC, padj)],
                    gwas_genes[, .(gene, coloc_pp4)], by = "gene")

strat_results <- dream_gwas[, .(
  n_genes = .N,
  mean_coloc = mean(coloc_pp4, na.rm = TRUE),
  n_coloc_05 = sum(coloc_pp4 > 0.5, na.rm = TRUE),
  prop_coloc_05 = sum(coloc_pp4 > 0.5, na.rm = TRUE) / .N
), by = t_quintile][order(t_quintile)]

cat("DEG quintile enrichment for COLOC:\n")
print(strat_results)

# Trend test
dream_gwas_valid <- dream_gwas[!is.na(t_quintile) & !is.na(coloc_pp4)]
if (nrow(dream_gwas_valid) > 10) {
  dream_gwas_valid[, quintile_num := as.integer(factor(t_quintile))]
  trend_test <- cor.test(dream_gwas_valid$quintile_num, dream_gwas_valid$coloc_pp4,
                         method = "spearman")
  cat("Trend test (DEG quintile vs COLOC):", "rho =", round(trend_test$estimate, 4),
      ", p =", format.pval(trend_test$p.value), "\n")
} else {
  cat("Trend test: insufficient data (", nrow(dream_gwas_valid), "rows)\n")
}

# ===========================================================================
# 5. Save results
# ===========================================================================
cat("\n=== Saving results ===\n")

fwrite(enrichment_dt, file.path(outdir, "celltype_heritability_results.csv"))
cat("  Saved celltype_heritability_results.csv:", nrow(enrichment_dt), "rows\n")

fwrite(specificity_dt, file.path(outdir, "celltype_specificity_scores.csv"))
cat("  Saved celltype_specificity_scores.csv:", nrow(specificity_dt), "rows\n")

fwrite(strat_results, file.path(outdir, "deg_quintile_coloc_enrichment.csv"))
cat("  Saved deg_quintile_coloc_enrichment.csv:", nrow(strat_results), "rows\n")

# Save gene sets for downstream use
saveRDS(ct_gene_sets, file.path(outdir, "celltype_gene_sets.rds"))
cat("  Saved celltype_gene_sets.rds\n")

cat("\nDone.\n")
