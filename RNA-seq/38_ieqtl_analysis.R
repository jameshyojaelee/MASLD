#!/usr/bin/env Rscript
# 38_ieqtl_analysis.R
# ---------------------------------------------------------------------------
# Interaction eQTL (ieQTL) Analysis for Disease-Specific Gene Regulation
#
# Identifies genes where eQTL effects significantly differ between MASLD
# and control conditions using sc-eQTL data from Zenodo. These disease-
# interacting eQTLs reveal genes whose regulation CHANGES in disease --
# a stronger biological signal than standard eQTLs.
#
# Steps:
#   1. Load significant ieQTL annotations from the master annotations file
#   2. Load disease-interaction ieQTL files for hepatocyte + other cell types
#   3. Extract significant interaction effects (FDR < 0.05)
#   4. Cross-reference with dream DEGs, deconvolution, and drug targets
#   5. Build per-gene output table and per-cell-type summary
#
# Inputs:
#   - sceQTL_output/significant_sc_eQTLs_with_annotations.txt.gz
#   - sceQTL_output/ieQTLs/interaction.{celltype}_disease_group3.txt.gz
#   - dream_results.csv, deconv_attribution_scores.csv, mr_convergent_drug_targets.csv
#
# Outputs:
#   - results/causal_inference/sceqtl/ieqtl_disease_genes.csv
#   - results/causal_inference/sceqtl/ieqtl_summary.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

cat("=== Script 38: Interaction eQTL Analysis ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# Configuration
# ==============================================================================
BASE_DIR    <- Sys.getenv("MASLD_PROJECT_ROOT",
                 "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
EQTL_DIR   <- file.path(BASE_DIR,
                 "RNA-seq/References/eQTL/MASLD_sc_eQTL_2025/Zenodo/sceQTL_output")
IEQTL_DIR  <- file.path(EQTL_DIR, "ieQTLs")
RESULTS_DIR <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/sceqtl")

ANNOT_FILE <- file.path(EQTL_DIR,
                 "significant_sc_eQTLs_with_annotations.txt.gz")
DREAM_FILE <- file.path(BASE_DIR,
                 "RNA-seq/Human/Patient_Cohorts/analysis/integration/results",
                 "integration/dream_results.csv")
GENE_CACHE <- file.path(BASE_DIR,
                 "RNA-seq/Human/Patient_Cohorts/analysis/integration/results",
                 "gene_annotation/human_ensg_to_symbol.tsv")
DECONV_FILE <- file.path(BASE_DIR,
                 "RNA-seq/results/causal_inference/deconv_attribution_scores.csv")
DRUG_FILE  <- file.path(BASE_DIR,
                 "RNA-seq/results/drug_repurposing/mr_convergent_drug_targets.csv")

# Thresholds
IEQTL_FDR_THR <- 0.05
DREAM_PADJ_THR <- 0.10

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Cell types with disease-interaction files
CELL_TYPES <- c("hepatocyte", "cholangiocyte", "stellate_cell", "endothelial_cell")

# Known MASLD genes for highlighting
MASLD_GENES <- c("PNPLA3", "TM6SF2", "HSD17B13", "MBOAT7", "GCKR", "MARC1",
                 "MTARC2", "CIDEB", "CIDEC", "COL1A1", "ACTA2", "TGFB1",
                 "TNF", "IL6", "IL1B", "CCL2", "PPARA", "PPARG", "NR1H4",
                 "THRB", "FGF19", "FGF21", "GDF15", "ADIPOQ", "LEP")

# ==============================================================================
# 1. Load gene annotation (Ensembl -> symbol)
# ==============================================================================
cat("--- Step 1: Loading gene annotation ---\n")

if (file.exists(GENE_CACHE)) {
  ann <- fread(GENE_CACHE)
  ann_map <- ann[symbol != "" & !is.na(symbol) & !duplicated(gene_base),
                 .(ensembl_id = gene_base, symbol)]
  cat("  Gene annotation loaded:", nrow(ann_map), "mappings\n")
} else {
  stop("Gene annotation file not found: ", GENE_CACHE)
}

# ==============================================================================
# 2. Load significant ieQTL annotations from master file
# ==============================================================================
cat("\n--- Step 2: Loading significant ieQTL annotations ---\n")

annot <- fread(ANNOT_FILE)
cat("  Total annotated sc-eQTLs:", nrow(annot), "\n")

# Extract per-cell-type ieQTL gene lists from boolean columns
ieqtl_cols <- grep("^is_ieQTL_", names(annot), value = TRUE)
cat("  ieQTL annotation columns:", paste(ieqtl_cols, collapse = ", "), "\n")

annot_ieqtl_summary <- lapply(ieqtl_cols, function(col) {
  ct <- sub("^is_ieQTL_", "", col)
  genes <- unique(annot[get(col) == TRUE, gene])
  data.table(cell_type = ct, n_annotated_ieqtl_genes = length(genes),
             genes = list(genes))
})
annot_ieqtl_summary <- rbindlist(annot_ieqtl_summary)

for (i in seq_len(nrow(annot_ieqtl_summary))) {
  cat(sprintf("  %s: %d annotated ieQTL genes\n",
              annot_ieqtl_summary$cell_type[i],
              annot_ieqtl_summary$n_annotated_ieqtl_genes[i]))
}

# ==============================================================================
# 3. Load raw disease-interaction ieQTL files
# ==============================================================================
cat("\n--- Step 3: Loading raw disease ieQTL files ---\n")

# Scan for all disease-interaction files
disease_files <- list.files(IEQTL_DIR, pattern = "disease", full.names = TRUE)
cat("  Found", length(disease_files), "disease interaction files:\n")
for (f in disease_files) cat("    ", basename(f), "\n")

# Process each disease-interaction file
ieqtl_all <- list()

for (f in disease_files) {
  # Extract cell type from filename: interaction.{celltype}_disease_group3.txt.gz
  fname <- basename(f)
  ct <- sub("^interaction\\.", "", fname)
  ct <- sub("_disease_group3\\.txt\\.gz$", "", ct)
  cat(sprintf("\n  Processing: %s (cell type: %s)\n", fname, ct))

  dt <- fread(f)
  cat(sprintf("    Total rows: %s\n", format(nrow(dt), big.mark = ",")))

  # Filter to interaction terms (G:phenoNAFLD)
  dt_int <- dt[grepl(":pheno", term, fixed = FALSE)]
  cat(sprintf("    Interaction term rows: %s\n",
              format(nrow(dt_int), big.mark = ",")))

  if (nrow(dt_int) == 0) {
    cat("    WARNING: No interaction terms found. Skipping.\n")
    next
  }

  # Apply FDR correction on lrt_pval (likelihood ratio test)
  dt_int[, lrt_fdr := p.adjust(lrt_pval, method = "BH")]

  # Also apply FDR on pval_full (interaction term p-value from full model)
  dt_int[, interaction_fdr := p.adjust(pval_full, method = "BH")]

  cat(sprintf("    Significant (lrt_fdr < %g): %d SNP-gene pairs\n",
              IEQTL_FDR_THR, sum(dt_int$lrt_fdr < IEQTL_FDR_THR, na.rm = TRUE)))
  cat(sprintf("    Significant (interaction_fdr < %g): %d SNP-gene pairs\n",
              IEQTL_FDR_THR,
              sum(dt_int$interaction_fdr < IEQTL_FDR_THR, na.rm = TRUE)))

  # Filter to significant ieQTLs (use lrt_fdr as primary)
  dt_sig <- dt_int[lrt_fdr < IEQTL_FDR_THR]
  cat(sprintf("    Unique significant genes: %d\n",
              uniqueN(dt_sig$gene)))

  if (nrow(dt_sig) == 0) {
    cat("    No significant ieQTLs at FDR <", IEQTL_FDR_THR, ". Skipping.\n")
    next
  }

  # Per-gene summary: pick top SNP (lowest lrt_pval), count SNPs
  gene_summary <- dt_sig[, .(
    interaction_beta = Estimate_full[which.min(lrt_pval)],
    interaction_pval = pval_full[which.min(lrt_pval)],
    interaction_fdr  = interaction_fdr[which.min(lrt_pval)],
    lrt_pval         = min(lrt_pval, na.rm = TRUE),
    lrt_fdr          = min(lrt_fdr, na.rm = TRUE),
    n_snps_ieqtl     = .N,
    top_snp          = snp_chr_pos_ref_alt[which.min(lrt_pval)]
  ), by = gene]

  gene_summary[, cell_type := ct]
  gene_summary[, fdr_threshold_used := "FDR_0.05"]
  ieqtl_all[[ct]] <- gene_summary

  cat(sprintf("    Genes with ieQTLs: %d\n", nrow(gene_summary)))
}

if (length(ieqtl_all) == 0) {
  cat("\nWARNING: No significant ieQTLs found across any cell type.\n")
  cat("Attempting with relaxed threshold (FDR < 0.10)...\n")

  # Retry with relaxed threshold
  IEQTL_FDR_THR_RELAXED <- 0.10
  for (f in disease_files) {
    fname <- basename(f)
    ct <- sub("^interaction\\.", "", fname)
    ct <- sub("_disease_group3\\.txt\\.gz$", "", ct)

    dt <- fread(f)
    dt_int <- dt[grepl(":pheno", term, fixed = FALSE)]
    if (nrow(dt_int) == 0) next

    dt_int[, lrt_fdr := p.adjust(lrt_pval, method = "BH")]
    dt_int[, interaction_fdr := p.adjust(pval_full, method = "BH")]

    dt_sig <- dt_int[lrt_fdr < IEQTL_FDR_THR_RELAXED]
    cat(sprintf("  %s: %d significant at FDR < %g\n",
                ct, uniqueN(dt_sig$gene), IEQTL_FDR_THR_RELAXED))

    if (nrow(dt_sig) == 0) next

    gene_summary <- dt_sig[, .(
      interaction_beta = Estimate_full[which.min(lrt_pval)],
      interaction_pval = pval_full[which.min(lrt_pval)],
      interaction_fdr  = interaction_fdr[which.min(lrt_pval)],
      lrt_pval         = min(lrt_pval, na.rm = TRUE),
      lrt_fdr          = min(lrt_fdr, na.rm = TRUE),
      n_snps_ieqtl     = .N,
      top_snp          = snp_chr_pos_ref_alt[which.min(lrt_pval)]
    ), by = gene]

    gene_summary[, cell_type := ct]
    gene_summary[, fdr_threshold_used := "FDR_0.10"]
    ieqtl_all[[ct]] <- gene_summary
  }

  if (length(ieqtl_all) > 0) {
    cat("  Using relaxed FDR threshold of", IEQTL_FDR_THR_RELAXED, "\n")
    IEQTL_FDR_THR <- IEQTL_FDR_THR_RELAXED
  }
}

# Combine all cell types
if (length(ieqtl_all) == 0) {
  # Still no results -- produce nominal (uncorrected) output for reference
  cat("\nNo FDR-significant ieQTLs found. Producing nominal (p < 0.05) output.\n")

  for (f in disease_files) {
    fname <- basename(f)
    ct <- sub("^interaction\\.", "", fname)
    ct <- sub("_disease_group3\\.txt\\.gz$", "", ct)

    dt <- fread(f)
    dt_int <- dt[grepl(":pheno", term, fixed = FALSE)]
    if (nrow(dt_int) == 0) next

    dt_int[, lrt_fdr := p.adjust(lrt_pval, method = "BH")]
    dt_int[, interaction_fdr := p.adjust(pval_full, method = "BH")]

    # Use nominal p-value for selection, still report FDR
    dt_sig <- dt_int[pval_full < 0.05]
    cat(sprintf("  %s: %d genes at nominal p < 0.05\n",
                ct, uniqueN(dt_sig$gene)))

    if (nrow(dt_sig) == 0) next

    gene_summary <- dt_sig[, .(
      interaction_beta = Estimate_full[which.min(lrt_pval)],
      interaction_pval = pval_full[which.min(lrt_pval)],
      interaction_fdr  = interaction_fdr[which.min(lrt_pval)],
      lrt_pval         = min(lrt_pval, na.rm = TRUE),
      lrt_fdr          = min(lrt_fdr, na.rm = TRUE),
      n_snps_ieqtl     = .N,
      top_snp          = snp_chr_pos_ref_alt[which.min(lrt_pval)]
    ), by = gene]

    gene_summary[, cell_type := ct]
    gene_summary[, fdr_threshold_used := "nominal_0.05"]
    ieqtl_all[[ct]] <- gene_summary
  }

  cat("  NOTE: Results use nominal p < 0.05 (not FDR-corrected).\n")
}

ieqtl <- rbindlist(ieqtl_all, use.names = TRUE)
cat(sprintf("\n  Total ieQTL gene-celltype entries: %d\n", nrow(ieqtl)))
cat(sprintf("  Unique ieQTL genes: %d\n", uniqueN(ieqtl$gene)))
cat(sprintf("  Cell types with results: %s\n",
            paste(unique(ieqtl$cell_type), collapse = ", ")))

# ==============================================================================
# 4. Load dream DEGs
# ==============================================================================
cat("\n--- Step 4: Loading dream DEGs ---\n")

dream <- fread(DREAM_FILE)
dream[, ensembl_id := sub("\\.\\d+$", "", gene)]
dream <- merge(dream, ann_map, by = "ensembl_id", all.x = TRUE)
cat("  Dream results:", nrow(dream), "genes\n")
cat("  DEGs (padj <", DREAM_PADJ_THR, "):",
    sum(dream$padj < DREAM_PADJ_THR, na.rm = TRUE), "\n")

# Build lookup: symbol -> dream stats
dream_lookup <- dream[!is.na(symbol) & symbol != "",
                      .(dream_logFC = logFC, dream_padj = padj,
                        dream_tstat = t),
                      by = symbol]
# Keep one entry per symbol (lowest padj)
dream_lookup <- dream_lookup[order(dream_padj)]
dream_lookup <- dream_lookup[!duplicated(symbol)]

# ==============================================================================
# 5. Load consensus tier from multi-evidence atlas (if available)
# ==============================================================================
cat("\n--- Step 5: Loading consensus tier info ---\n")

atlas_file <- file.path(BASE_DIR,
                        "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
if (file.exists(atlas_file)) {
  atlas <- fread(atlas_file,
                 select = c("human_symbol", "human_consensus_tier"))
  setnames(atlas, "human_symbol", "symbol")
  atlas <- atlas[!is.na(symbol) & symbol != ""]
  atlas <- atlas[!duplicated(symbol)]
  cat("  Atlas loaded:", nrow(atlas), "genes with tier info\n")
} else {
  cat("  Multi-evidence atlas not found. Consensus tier will be NA.\n")
  atlas <- data.table(symbol = character(), human_consensus_tier = character())
}

# ==============================================================================
# 6. Load deconvolution attribution
# ==============================================================================
cat("\n--- Step 6: Loading deconvolution attribution ---\n")

if (file.exists(DECONV_FILE)) {
  deconv <- fread(DECONV_FILE)
  deconv[, ensembl_id := sub("\\.\\d+$", "", gene)]
  deconv <- merge(deconv[, .(ensembl_id, category)],
                  ann_map, by = "ensembl_id", all.x = TRUE)
  hep_intrinsic <- unique(deconv[category == "Hepatocyte_intrinsic", symbol])
  hep_intrinsic <- hep_intrinsic[!is.na(hep_intrinsic)]
  cat("  Hepatocyte-intrinsic genes:", length(hep_intrinsic), "\n")
} else {
  cat("  WARNING: Deconvolution file not found. Skipping.\n")
  hep_intrinsic <- character(0)
}

# ==============================================================================
# 7. Load drug targets
# ==============================================================================
cat("\n--- Step 7: Loading drug targets ---\n")

if (file.exists(DRUG_FILE)) {
  drugs <- fread(DRUG_FILE)
  drug_genes <- unique(drugs$gene)
  cat("  Multi-layer drug targets:", length(drug_genes), "\n")
} else {
  cat("  WARNING: Drug targets file not found. Skipping.\n")
  drug_genes <- character(0)
}

# ==============================================================================
# 8. Cross-reference and build output table
# ==============================================================================
cat("\n--- Step 8: Cross-referencing ieQTL genes ---\n")

# Merge dream info
ieqtl <- merge(ieqtl, dream_lookup, by.x = "gene", by.y = "symbol",
               all.x = TRUE)

# Merge consensus tier
ieqtl <- merge(ieqtl, atlas, by.x = "gene", by.y = "symbol",
               all.x = TRUE)
setnames(ieqtl, "human_consensus_tier", "consensus_tier",
         skip_absent = TRUE)

# Add flags
ieqtl[, is_deg := !is.na(dream_padj) & dream_padj < DREAM_PADJ_THR]
ieqtl[, is_hep_intrinsic := gene %in% hep_intrinsic]
ieqtl[, is_drug_target := gene %in% drug_genes]

# Reorder columns
out_cols <- c("gene", "cell_type", "fdr_threshold_used",
              "interaction_beta", "interaction_pval",
              "interaction_fdr", "lrt_pval", "lrt_fdr", "n_snps_ieqtl",
              "top_snp", "is_deg", "dream_logFC", "dream_padj",
              "consensus_tier", "is_hep_intrinsic", "is_drug_target")
# Only keep columns that exist
out_cols <- intersect(out_cols, names(ieqtl))
ieqtl_out <- ieqtl[order(lrt_pval), ..out_cols]

# Save
out_file <- file.path(RESULTS_DIR, "ieqtl_disease_genes.csv")
fwrite(ieqtl_out, out_file)
cat("  Saved:", out_file, "\n")
cat("  Rows:", nrow(ieqtl_out), "\n")

# ==============================================================================
# 9. Build summary table
# ==============================================================================
cat("\n--- Step 9: Building summary ---\n")

summary_dt <- ieqtl[, .(
  n_ieqtl_genes    = uniqueN(gene),
  n_overlap_deg    = uniqueN(gene[is_deg == TRUE]),
  n_hep_intrinsic  = uniqueN(gene[is_hep_intrinsic == TRUE]),
  n_drug_target    = uniqueN(gene[is_drug_target == TRUE]),
  top_genes        = paste(head(gene[order(lrt_pval)], 10), collapse = ",")
), by = cell_type]

# Add annotated ieQTL gene counts from the master annotations file
for (i in seq_len(nrow(annot_ieqtl_summary))) {
  ct <- annot_ieqtl_summary$cell_type[i]
  n  <- annot_ieqtl_summary$n_annotated_ieqtl_genes[i]
  if (ct %in% summary_dt$cell_type) {
    summary_dt[cell_type == ct, n_annotated_ieqtl := n]
  }
}

# Save
sum_file <- file.path(RESULTS_DIR, "ieqtl_summary.csv")
fwrite(summary_dt, sum_file)
cat("  Saved:", sum_file, "\n\n")

# ==============================================================================
# 10. Print summary
# ==============================================================================
cat("=== ieQTL Analysis Summary ===\n\n")
cat(sprintf("FDR threshold used: %g\n", IEQTL_FDR_THR))
cat(sprintf("Total unique ieQTL genes: %d\n\n", uniqueN(ieqtl$gene)))

for (ct in unique(summary_dt$cell_type)) {
  row <- summary_dt[cell_type == ct]
  cat(sprintf("--- %s ---\n", ct))
  cat(sprintf("  ieQTL genes (disease interaction): %d\n", row$n_ieqtl_genes))
  if ("n_annotated_ieqtl" %in% names(row)) {
    cat(sprintf("  Annotated ieQTL genes (master file): %d\n",
                row$n_annotated_ieqtl))
  }
  cat(sprintf("  Overlap with dream DEGs: %d (%.1f%%)\n",
              row$n_overlap_deg,
              100 * row$n_overlap_deg / max(row$n_ieqtl_genes, 1)))
  cat(sprintf("  Hepatocyte-intrinsic: %d\n", row$n_hep_intrinsic))
  cat(sprintf("  Drug targets: %d\n", row$n_drug_target))
  cat(sprintf("  Top genes: %s\n\n", row$top_genes))
}

# Highlight known MASLD genes among ieQTLs
masld_hits <- ieqtl[gene %in% MASLD_GENES]
if (nrow(masld_hits) > 0) {
  cat("=== Known MASLD genes with disease ieQTLs ===\n")
  for (i in seq_len(nrow(masld_hits))) {
    r <- masld_hits[i]
    cat(sprintf("  %s [%s]: interaction beta = %.3f, lrt_pval = %.2e, FDR = %.2e",
                r$gene, r$cell_type, r$interaction_beta,
                r$lrt_pval, r$lrt_fdr))
    if (r$is_deg) {
      cat(sprintf(", DEG logFC = %.2f", r$dream_logFC))
    }
    cat("\n")
  }
} else {
  cat("No known MASLD genes found among disease ieQTLs.\n")
}

cat("\n=== Script 38 Complete ===\n")
cat("End time:", format(Sys.time()), "\n")
