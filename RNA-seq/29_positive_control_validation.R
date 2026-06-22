#!/usr/bin/env Rscript
# 29_positive_control_validation.R
# Cross-reference 65 positive control genes against all evidence layers

library(data.table)
library(biomaRt)

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUTDIR <- file.path(BASE, "RNA-seq/results/validation")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# --- Load positive controls ---
pc <- fread(file.path(BASE, "results/library/positive_control.csv"))
cat("Positive control genes:", nrow(pc), "\n")
# First column is "Gene symbol"
pc_symbols <- pc[["Gene symbol"]]

# --- Load evidence layers ---
consensus <- fread(file.path(RDIR, "consensus_degs.csv"))
concordance <- fread(file.path(BASE, "Analysis/Cross_Species_Concordance/results/concordance_atlas_unified.csv"))
# mr_file removed 2026-04-22 — MR ditched from paper.

# --- Map Ensembl to symbols ---
cat("Mapping Ensembl IDs to HGNC symbols...\n")
# Use local GENCODE v49 annotation cache (avoids biomaRt network dependency)
ann_cache <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation/human_ensg_to_symbol.tsv")
if (file.exists(ann_cache)) {
  ann <- fread(ann_cache)
  consensus[, ensembl_clean := sub("\\..*", "", gene)]
  if ("symbol" %in% names(consensus)) {
    setnames(consensus, "symbol", "hgnc_symbol")
  } else {
    consensus <- merge(consensus, ann[, .(gene_base, symbol)],
                       by.x = "ensembl_clean", by.y = "gene_base", all.x = TRUE)
    setnames(consensus, "symbol", "hgnc_symbol")
  }
  cat("  Mapped", sum(!is.na(consensus$hgnc_symbol)), "/", nrow(consensus), "genes\n")
} else {
  # Fallback: try biomaRt
  mart <- tryCatch(useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"), error = function(e) NULL)
  if (!is.null(mart)) {
    ens_ids <- sub("\\..*", "", consensus$gene)
    map_df <- as.data.table(getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                                   filters = "ensembl_gene_id", values = ens_ids, mart = mart))
    consensus[, ensembl_clean := sub("\\..*", "", gene)]
    consensus <- merge(consensus, map_df, by.x = "ensembl_clean", by.y = "ensembl_gene_id", all.x = TRUE)
  } else {
    cat("  Local annotation cache and biomaRt both unavailable. Symbol mapping skipped.\n")
    consensus[, hgnc_symbol := NA_character_]
  }
}

# --- Build validation table ---
results <- data.table(gene = pc_symbols)

# Consensus DEGs
if ("hgnc_symbol" %in% names(consensus)) {
  # Select available columns (meta_logFC/meta_padj may not exist post-threshold change)
  avail_cols <- intersect(names(consensus), c("bulk_logFC", "bulk_padj", "bulk_sig", "bulk_dir"))
  cons_match <- consensus[hgnc_symbol %in% pc_symbols, c("hgnc_symbol", avail_cols), with = FALSE]
  setnames(cons_match, "hgnc_symbol", "gene")
  cons_match <- cons_match[!duplicated(gene)]
  results <- merge(results, cons_match, by = "gene", all.x = TRUE)
} else {
  results[, c("bulk_logFC", "bulk_padj", "bulk_sig") := .(NA, NA, NA)]
}

# Concordance atlas (column is n_concordant, not best_n_concordant)
conc_match <- concordance[human_symbol %in% pc_symbols, .(gene = human_symbol, concordance_category = primary_category, n_concordant)]
conc_match <- conc_match[!duplicated(gene)]
results <- merge(results, conc_match, by = "gene", all.x = TRUE)

# TWAS/MR merge REMOVED 2026-04-22 — MR ditched from paper.
# TWAS-only signal is surfaced via the multi-evidence atlas (Script 27a).

# Classify: GWAS-type vs expression-type
gwas_genes <- c("PNPLA3", "TM6SF2", "GCKR", "HSD17B13", "MBOAT7", "FTO", "SERPINA1")
results[, control_type := fifelse(gene %in% gwas_genes, "GWAS_variant", "Expression_driven")]

# Summary
cat("\nValidation summary:\n")
cat("  In concordance atlas:", sum(!is.na(results$concordance_category)), "/", nrow(results), "\n")
# TWAS signal tally skipped here 2026-04-22 (MR-joined summary retired).
if ("twas_fdr" %in% names(results)) {
  cat("  With TWAS signal (fdr<0.1):", sum(!is.na(results$twas_fdr) & results$twas_fdr < 0.1), "/", nrow(results), "\n")
}
if ("bulk_sig" %in% names(results)) {
  cat("  In consensus DEGs (bulk_sig):", sum(results$bulk_sig == TRUE, na.rm = TRUE), "/", nrow(results), "\n")
}

fwrite(results, file.path(OUTDIR, "positive_control_validation.csv"))
cat("Saved:", file.path(OUTDIR, "positive_control_validation.csv"), "\n")
