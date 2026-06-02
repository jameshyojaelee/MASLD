#!/usr/bin/env Rscript
# 205_proportion_coloc.R — Cell-Type Proportion × COLOC Gene Expression
#
# Tests whether expression of COLOC-confirmed genes correlates with
# deconvolved cell-type proportions, linking genetic causal evidence
# to cellular composition changes.
#
# Inputs:
#   - MuSiC deconvolution proportions: Analysis/Deconvolution/results/{dataset}/
#   - COLOC gene-level: GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv
#   - sc-eQTL cell-type assignments: RNA-seq/results/causal_inference/sceqtl/
#   - Bulk expression: merged_dge.rds (corrected logCPM)
#   - Sample metadata: unified_metadata.csv
#
# Outputs:
#   - RNA-seq/results/gwas_rna_integration/proportion_coloc_correlation.csv
#   - RNA-seq/results/gwas_rna_integration/proportion_coloc_summary.csv
#
# SLURM: cpu partition, 4 CPUs, 32GB, 48h

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INTEGRATION <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")

outdir <- file.path(BASE, "RNA-seq/results/gwas_rna_integration")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# 1. Load MuSiC deconvolution proportions (all 10 datasets)
# ===========================================================================
cat("Loading MuSiC deconvolution proportions...\n")
music_dir <- file.path(BASE, "Analysis/Deconvolution/results")
datasets <- list.dirs(music_dir, full.names = FALSE, recursive = FALSE)
datasets <- datasets[!datasets %in% c("summary_plots")]

music_list <- lapply(datasets, function(ds) {
  prop_file <- file.path(music_dir, ds, paste0(ds, "_music_prop_weighted.tsv"))
  if (!file.exists(prop_file)) {
    cat("  WARNING: No MuSiC proportions for", ds, "\n")
    return(NULL)
  }
  # Use read.table which handles row names in TSVs properly
  mat <- read.table(prop_file, header = TRUE, sep = "\t",
                    check.names = FALSE, row.names = 1)
  df <- as.data.table(mat, keep.rownames = "sample_id")
  df[, dataset := ds]
  return(df)
})

music_props <- rbindlist(music_list, fill = TRUE, use.names = TRUE)
cat("  Loaded proportions for", nrow(music_props), "samples across",
    uniqueN(music_props$dataset), "datasets\n")

# Key cell types to test (column names may have spaces)
ct_cols <- intersect(names(music_props),
  c("Hepatocytes", "Macrophages", "Fibroblasts", "Endothelial cells",
    "T cells", "B cells", "Cholangiocytes", "Circulating NK/NKT",
    "Resident NK", "Mono+mono derived cells", "Plasma cells",
    "Cholangiocytes", "Basophils", "Neutrophils"))
# If empty, try all numeric columns (except sample_id, dataset)
if (length(ct_cols) == 0) {
  ct_cols <- setdiff(names(music_props)[sapply(music_props, is.numeric)],
                     c("sample_id", "dataset"))
}
cat("  Cell types available:", paste(ct_cols, collapse = ", "), "\n")

# ===========================================================================
# 2. Load bulk expression (corrected logCPM)
# ===========================================================================
cat("Loading bulk expression matrix...\n")
dge <- readRDS(file.path(INTEGRATION, "results/integration/merged_dge.rds"))
logcpm <- edgeR::cpm(dge, log = TRUE)
cat("  Expression matrix:", nrow(logcpm), "genes x", ncol(logcpm), "samples\n")

# Load metadata
meta <- fread(file.path(INTEGRATION, "metadata/unified_metadata.csv"))
cat("  Metadata:", nrow(meta), "samples\n")

# ===========================================================================
# 3. Load COLOC + sc-eQTL cell-type assignments
# ===========================================================================
cat("Loading COLOC results...\n")
coloc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc <- coloc[gene != "" & !is.na(gene)]

# sc-eQTL cell-type COLOC
sceqtl_file <- file.path(BASE,
  "RNA-seq/results/causal_inference/sceqtl/sceqtl_multi_celltype_summary.csv")
if (file.exists(sceqtl_file)) {
  sceqtl <- fread(sceqtl_file)
  cat("  sc-eQTL cell-type assignments:", nrow(sceqtl), "genes\n")
} else {
  cat("  WARNING: sc-eQTL summary not found, using bulk COLOC only\n")
  sceqtl <- data.table(gene = character(), best_cell_type = character())
}

# COLOC genes with PP.H4 > 0.5
coloc_genes <- coloc[coloc_best_pp4 > 0.5, .(gene, coloc_pp4 = coloc_best_pp4)]
cat("  COLOC genes (PP.H4 > 0.5):", nrow(coloc_genes), "\n")

# ===========================================================================
# 4. Match samples across expression + proportions + metadata
# ===========================================================================
cat("\nMatching samples...\n")

# Samples in both expression and proportions
common_samples <- intersect(colnames(logcpm), music_props$sample_id)
cat("  Samples in both expression and proportions:", length(common_samples), "\n")

# Subset and align
expr_mat <- logcpm[, common_samples]
prop_dt <- music_props[sample_id %in% common_samples]
prop_dt <- prop_dt[match(common_samples, sample_id)]

# Add disease status from metadata
meta_matched <- meta[match(common_samples, meta$sample_id)]
disease_status <- meta_matched$group_binary  # Control vs Disease

cat("  Disease status: Control =", sum(disease_status == "Control", na.rm = TRUE),
    ", Disease =", sum(disease_status == "Disease", na.rm = TRUE), "\n")

# ===========================================================================
# 5. Correlation analysis: COLOC gene expression × cell-type proportion
# ===========================================================================
cat("\n=== Running correlation analysis ===\n")

# Function to test one gene × one cell type
test_gene_ct <- function(gene_symbol, ct, expr_vec, prop_vec, disease_vec) {
  # Remove NAs
  valid <- !is.na(expr_vec) & !is.na(prop_vec)
  if (sum(valid) < 30) return(NULL)

  e <- expr_vec[valid]
  p <- prop_vec[valid]
  d <- disease_vec[valid]

  # Simple correlation
  cor_test <- cor.test(e, p, method = "spearman")

  # Partial correlation adjusting for disease status
  d_num <- as.numeric(factor(d))
  resid_e <- residuals(lm(e ~ d_num))
  resid_p <- residuals(lm(p ~ d_num))
  partial_cor <- cor.test(resid_e, resid_p, method = "spearman")

  data.table(
    gene = gene_symbol,
    cell_type = ct,
    rho = cor_test$estimate,
    pvalue = cor_test$p.value,
    partial_rho = partial_cor$estimate,
    partial_pvalue = partial_cor$p.value,
    n_samples = sum(valid)
  )
}

# Map gene symbols to rownames (may be ENSEMBL IDs)
# Check if rownames are gene symbols or ENSEMBL
sample_row <- rownames(expr_mat)[1:5]
is_ensembl <- any(grepl("^ENSG", sample_row))
cat("  Expression matrix row format:", ifelse(is_ensembl, "ENSEMBL", "symbol"), "\n")

if (is_ensembl) {
  # Need to map gene symbols to ENSEMBL IDs
  # Use COLOC ensembl column; expression rownames may have version suffix (ENSG.XX)
  gene_map <- coloc[gene != "" & ensembl != "", .(gene, ensembl)]
  coloc_genes_mapped <- merge(coloc_genes, gene_map, by = "gene")
  # Try exact match first
  coloc_genes_mapped <- coloc_genes_mapped[ensembl %in% rownames(expr_mat)]
  if (nrow(coloc_genes_mapped) == 0) {
    # Try matching with version suffix stripped or added
    expr_base <- sub("\\.\\d+$", "", rownames(expr_mat))
    coloc_base <- sub("\\.\\d+$", "", gene_map$ensembl)
    # Create map: base ENSEMBL → full rowname
    base_to_full <- setNames(rownames(expr_mat), expr_base)
    gene_map[, ensembl_base := sub("\\.\\d+$", "", ensembl)]
    gene_map[, expr_id := base_to_full[ensembl_base]]
    gene_map <- gene_map[!is.na(expr_id)]
    coloc_genes_mapped <- merge(coloc_genes, gene_map[, .(gene, ensembl = expr_id)], by = "gene")
  }
  cat("  Mapped COLOC genes to expression:", nrow(coloc_genes_mapped), "\n")
  gene_id_col <- "ensembl"
} else {
  coloc_genes_mapped <- coloc_genes[gene %in% rownames(expr_mat)]
  coloc_genes_mapped[, ensembl := gene]
  cat("  COLOC genes in expression:", nrow(coloc_genes_mapped), "\n")
  gene_id_col <- "gene"
}

# Run for all COLOC genes × all cell types
results_list <- list()
for (i in seq_len(nrow(coloc_genes_mapped))) {
  g <- coloc_genes_mapped$gene[i]
  g_id <- coloc_genes_mapped[[gene_id_col]][i]

  # Get expression vector - handle potential partial matches in ENSEMBL IDs
  expr_rows <- grep(paste0("^", g_id), rownames(expr_mat), value = TRUE)
  if (length(expr_rows) == 0) next
  expr_vec <- expr_mat[expr_rows[1], ]

  for (ct in ct_cols) {
    prop_vec <- prop_dt[[ct]]
    if (is.null(prop_vec) || all(is.na(prop_vec))) next

    res <- test_gene_ct(g, ct, expr_vec, prop_vec, disease_status)
    if (!is.null(res)) {
      res[, coloc_pp4 := coloc_genes_mapped$coloc_pp4[i]]
      results_list[[length(results_list) + 1]] <- res
    }
  }
}

results_dt <- rbindlist(results_list, fill = TRUE)
cat("\nTotal gene × cell-type tests:", nrow(results_dt), "\n")

if (nrow(results_dt) == 0) {
  cat("WARNING: No gene × cell-type tests completed. Check sample ID matching.\n")
  # Save empty results
  fwrite(data.table(gene = character(), cell_type = character()),
         file.path(outdir, "proportion_coloc_correlation.csv"))
  fwrite(data.table(cell_type = character(), note = "No results"),
         file.path(outdir, "proportion_coloc_summary.csv"))
  cat("Done (no results).\n")
  quit(save = "no", status = 0)
}

# FDR correction
results_dt[, fdr := p.adjust(pvalue, method = "BH")]
results_dt[, partial_fdr := p.adjust(partial_pvalue, method = "BH")]

# Add sc-eQTL cell-type assignment
if (nrow(sceqtl) > 0) {
  results_dt <- merge(results_dt,
    sceqtl[, .(gene, sceqtl_best_celltype = best_cell_type)],
    by = "gene", all.x = TRUE)
  # Flag if correlation cell type matches sc-eQTL cell type
  results_dt[, celltype_concordant := (cell_type == sceqtl_best_celltype)]
}

# ===========================================================================
# 6. Compare COLOC vs non-COLOC genes (permutation test)
# ===========================================================================
cat("\n=== Permutation test: COLOC vs non-COLOC genes ===\n")

# For hepatocytes: compare mean |rho| of COLOC genes vs random gene sets
hep_results <- results_dt[cell_type == "Hepatocytes"]
mean_rho_coloc <- mean(abs(hep_results$rho), na.rm = TRUE)

# Random gene sets (1000 permutations)
set.seed(42)
n_coloc <- nrow(coloc_genes_mapped)
all_genes <- rownames(expr_mat)
n_perm <- 1000

random_rhos <- sapply(seq_len(n_perm), function(iter) {
  random_genes <- sample(all_genes, min(n_coloc, length(all_genes)))
  rhos <- sapply(random_genes, function(g_id) {
    expr_vec <- expr_mat[g_id, ]
    hep_prop <- prop_dt[["Hepatocytes"]]
    if (is.null(hep_prop) || all(is.na(hep_prop))) return(NA)
    valid <- !is.na(expr_vec) & !is.na(hep_prop)
    if (sum(valid) < 30) return(NA)
    cor(expr_vec[valid], hep_prop[valid], method = "spearman")
  })
  mean(abs(rhos), na.rm = TRUE)
})

perm_pvalue <- mean(random_rhos >= mean_rho_coloc)
cat("  Mean |rho| COLOC genes:", round(mean_rho_coloc, 4), "\n")
cat("  Mean |rho| random genes:", round(mean(random_rhos), 4), "\n")
cat("  Permutation p-value:", perm_pvalue, "\n")

# ===========================================================================
# 7. Summary by cell type
# ===========================================================================
cat("\n=== Summary by cell type ===\n")
summary_dt <- results_dt[, .(
  n_genes = .N,
  n_sig_005 = sum(fdr < 0.05, na.rm = TRUE),
  n_sig_01 = sum(fdr < 0.1, na.rm = TRUE),
  mean_rho = mean(rho, na.rm = TRUE),
  mean_abs_rho = mean(abs(rho), na.rm = TRUE),
  mean_partial_rho = mean(partial_rho, na.rm = TRUE)
), by = cell_type]

print(summary_dt[order(-mean_abs_rho)])

# ===========================================================================
# 8. Save results
# ===========================================================================
cat("\n=== Saving results ===\n")

fwrite(results_dt, file.path(outdir, "proportion_coloc_correlation.csv"))
cat("  Saved proportion_coloc_correlation.csv:", nrow(results_dt), "rows\n")

# Summary with permutation result
summary_dt[, perm_pvalue_vs_random := perm_pvalue]
fwrite(summary_dt, file.path(outdir, "proportion_coloc_summary.csv"))
cat("  Saved proportion_coloc_summary.csv:", nrow(summary_dt), "rows\n")

cat("\nDone.\n")
