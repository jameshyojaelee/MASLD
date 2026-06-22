#!/usr/bin/env Rscript
# 206_gwas_rna_integration.R — Unified GWAS-RNA-scRNA Integration Table
#
# Merges all Phase A results (200-205) with the existing multi-evidence atlas
# into a single co-analysis table with convergence tier classification.
#
# Run AFTER all Phase A scripts have completed.
#
# Inputs (all from RNA-seq/results/gwas_rna_integration/):
#   - intact_scores.csv (Script 200)
#   - celltype_heritability_results.csv (Script 201)
#   - geneset_enrichment_results.csv (Script 202)
#   - mr_mediation_results.csv (Script 203)
#   - scdrs_celltype_enrichment.csv (Script 204)
#   - proportion_coloc_correlation.csv (Script 205)
#   - Multi-evidence atlas: RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
#
# Outputs:
#   - RNA-seq/results/gwas_rna_integration/gwas_rna_scrna_integrated.csv
#   - RNA-seq/results/gwas_rna_integration/convergence_summary.csv
#
# SLURM: cpu partition, 4 CPUs, 16GB, 48h

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

indir  <- file.path(BASE, "RNA-seq/results/gwas_rna_integration")
outdir <- indir

# ===========================================================================
# 1. Load multi-evidence atlas (base table)
# ===========================================================================
cat("Loading multi-evidence atlas...\n")
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
cat("  Atlas:", nrow(atlas), "genes x", ncol(atlas), "columns\n")
stopifnot(all(c("bulk_padj", "bulk_logFC") %in% names(atlas)))  # C2: human-bulk channel (was dream_*)

# Standardize gene column name
gene_col <- intersect(names(atlas), c("human_symbol", "gene", "gene_symbol"))[1]
if (gene_col != "gene") {
  setnames(atlas, gene_col, "gene")
}

# ===========================================================================
# 2. Load INTACT scores (Script 200)
# ===========================================================================
cat("Loading INTACT scores...\n")
intact_file <- file.path(indir, "intact_scores.csv")
if (file.exists(intact_file)) {
  intact <- fread(intact_file)
  cat("  INTACT:", nrow(intact), "genes\n")
  # Select available columns (FDR removed in v4 — scores are posteriors, not p-values)
  intact_cols <- intersect(names(intact),
    c("gene", "intact_score_bulk", "intact_score_ct", "intact_best_celltype",
      "intact_n_celltypes_05", "twas_z_best", "coloc_pp4"))
  atlas <- merge(atlas, intact[, ..intact_cols], by = "gene", all.x = TRUE)
} else {
  cat("  WARNING: INTACT scores not found\n")
}

# ===========================================================================
# 3. Load cell-type heritability (Script 201) — gene-level not per-gene,
#    but cell-type specificity scores are per-gene
# ===========================================================================
cat("Loading cell-type specificity scores...\n")
spec_file <- file.path(indir, "celltype_specificity_scores.csv")
if (file.exists(spec_file)) {
  spec <- fread(spec_file)
  cat("  Specificity:", nrow(spec), "genes\n")
  atlas <- merge(atlas,
    spec[, .(gene, sc_tau = tau, sc_n_ct_sig = n_ct_sig, sc_best_ct = best_ct,
             sc_best_ct_logFC = best_ct_logFC)],
    by = "gene", all.x = TRUE)
} else {
  cat("  WARNING: Specificity scores not found\n")
}

# ===========================================================================
# 4. (Removed — MR mediation excluded from pipeline)
# ===========================================================================

# ===========================================================================
# 5. Load scDRS cell-type enrichment (Script 204) — summary per cell type,
#    plus within-cell-type disease vs control comparison.
#    Annotate each gene with scDRS info from its best cell type.
# ===========================================================================
cat("Loading scDRS results...\n")
scdrs_ct_file <- file.path(indir, "scdrs_celltype_enrichment.csv")
scdrs_wt_file <- file.path(indir, "scdrs_within_celltype.csv")

if (file.exists(scdrs_ct_file)) {
  scdrs_ct <- fread(scdrs_ct_file)
  cat("  scDRS cell-type enrichment:", nrow(scdrs_ct), "cell types\n")

  # CAVEAT (M3, mega-review A6.3/A6.7 remediation 2026-06-21):
  # scdrs_celltype_enrichment.csv carries PER-CELL scDRS statistics
  # (n_cells ~ 6.6e5; mc_pvalue at the 1/n_mc floor, fdr ~0.008) — these are
  # PSEUDOREPLICATED: every cell is treated as an independent observation, so
  # the p-values/FDR are inflated to the Monte-Carlo floor and carry no valid
  # per-cell-type significance. We therefore retire ALL per-cell p-value-derived
  # outputs (mc_pvalue / prop_sig / fdr / enriched flags / "best enriched" by
  # min-p) from the atlas. Only the DESCRIPTIVE mean scDRS score is kept, and it
  # is NOT used to gate the convergence tier. The CANONICAL cell-type
  # heritability / GWAS-enrichment arm is the DONOR-level MAGMA analysis
  # (Script 201, celltype_heritability_results.csv) — cite that for any
  # cell-type GWAS-enrichment significance claim, not these columns.

  # Build per-cell-type DESCRIPTIVE scDRS lookup (mean score only — no p-values)
  scdrs_lookup <- scdrs_ct[, .(cell_type, scdrs_mean_score = mean_scdrs)]

  # If genes have a best cell type (from specificity, section 3), annotate
  # with that cell type's DESCRIPTIVE mean scDRS score (no significance gate)
  if ("sc_best_ct" %in% names(atlas)) {
    atlas <- merge(atlas,
      scdrs_lookup[, .(sc_best_ct = cell_type, scdrs_ct_mean_score = scdrs_mean_score)],
      by = "sc_best_ct", all.x = TRUE)
    cat("  Genes annotated with descriptive scDRS mean score:",
        sum(!is.na(atlas$scdrs_ct_mean_score)), "\n")
  } else {
    cat("  NOTE: sc_best_ct not available; adding descriptive scDRS column only\n")
    atlas[, scdrs_ct_mean_score := NA_real_]
  }
} else {
  cat("  WARNING: scDRS cell-type enrichment not found\n")
}

# 5b. Load scDRS within-cell-type disease vs control (Script 204)
if (file.exists(scdrs_wt_file)) {
  scdrs_wt <- fread(scdrs_wt_file)
  cat("  scDRS within-cell-type:", nrow(scdrs_wt), "cell types\n")

  # Annotate genes with disease vs control scDRS shift in their best cell type
  if ("sc_best_ct" %in% names(atlas)) {
    wt_lookup <- scdrs_wt[, .(sc_best_ct = cell_type,
                               scdrs_disease_mean = mean_disease,
                               scdrs_control_mean = mean_control,
                               scdrs_disease_p = mannwhitney_p)]
    atlas <- merge(atlas, wt_lookup, by = "sc_best_ct", all.x = TRUE)
    atlas[, scdrs_disease_shift := scdrs_disease_mean - scdrs_control_mean]
    cat("  Genes with disease scDRS shift:", sum(!is.na(atlas$scdrs_disease_shift)), "\n")
  }
} else {
  cat("  WARNING: scDRS within-cell-type results not found\n")
}

# ===========================================================================
# 6. Load proportion × COLOC (Script 205)
# ===========================================================================
cat("Loading proportion × COLOC results...\n")
prop_file <- file.path(indir, "proportion_coloc_correlation.csv")
if (file.exists(prop_file)) {
  prop <- fread(prop_file)
  # Best correlation per gene (strongest cell-type association)
  prop_best <- prop[, .(
    prop_coloc_rho = rho[which.max(abs(rho))],
    prop_coloc_p = pvalue[which.max(abs(rho))],
    prop_coloc_ct = cell_type[which.max(abs(rho))]
  ), by = gene]
  cat("  Proportion COLOC:", nrow(prop_best), "genes\n")
  atlas <- merge(atlas, prop_best, by = "gene", all.x = TRUE)
} else {
  cat("  WARNING: Proportion COLOC results not found\n")
}

# ===========================================================================
# 7. Derive convergence tier
# ===========================================================================
cat("\n=== Deriving convergence tiers ===\n")

# Check which columns exist
has_coloc <- "coloc_best_pp4" %in% names(atlas) || "broadaway_coloc_pp4" %in% names(atlas)
has_dream <- "bulk_padj" %in% names(atlas)
has_intact <- "intact_score_bulk" %in% names(atlas)
has_mediation <- "mediation_proportion" %in% names(atlas)

# Get the COLOC column name
coloc_col <- intersect(names(atlas), c("coloc_best_pp4", "broadaway_coloc_pp4"))[1]

# Count evidence sources per gene
atlas[, n_evidence_gwas_rna := 0L]

# 1. Bulk DEG (limma-voom-qw C2)
if (has_dream) {
  # NOTE: Exploratory DEG flag for GWAS overlap annotation. Primary threshold: padj<0.05 + |logFC|>0.3.
  atlas[, is_deg := bulk_padj < 0.1 & !is.na(bulk_padj)]
  atlas[is_deg == TRUE, n_evidence_gwas_rna := n_evidence_gwas_rna + 1L]
}

# 2. COLOC
if (!is.null(coloc_col) && coloc_col %in% names(atlas)) {
  atlas[, is_coloc := get(coloc_col) > 0.5 & !is.na(get(coloc_col))]
  atlas[is_coloc == TRUE, n_evidence_gwas_rna := n_evidence_gwas_rna + 1L]
}

# 3. INTACT
if (has_intact) {
  atlas[, is_intact := intact_score_bulk > 0.5 & !is.na(intact_score_bulk)]
  atlas[is_intact == TRUE, n_evidence_gwas_rna := n_evidence_gwas_rna + 1L]
}

# 4. Cell-type resolved (sc-eQTL or scDRS)
if ("intact_score_ct" %in% names(atlas)) {
  atlas[, is_ct_resolved := intact_score_ct > 0.5 & !is.na(intact_score_ct)]
  atlas[is_ct_resolved == TRUE, n_evidence_gwas_rna := n_evidence_gwas_rna + 1L]
}

# 5. (MR removed from pipeline)

# Convergence tier
atlas[, gwas_rna_convergence_tier := fcase(
  n_evidence_gwas_rna >= 4, "Strong",
  n_evidence_gwas_rna >= 3, "Moderate",
  n_evidence_gwas_rna >= 2, "Suggestive",
  n_evidence_gwas_rna >= 1, "Single_source",
  default = "No_evidence"
)]

cat("Convergence tier distribution:\n")
print(atlas[, .N, by = gwas_rna_convergence_tier][order(-N)])

# ===========================================================================
# 8. Summary statistics
# ===========================================================================
cat("\n=== Integration summary ===\n")
cat("Total genes:", nrow(atlas), "\n")
cat("Total columns:", ncol(atlas), "\n")

new_cols <- c("intact_score_bulk", "intact_score_ct",
              "intact_best_celltype", "sc_tau", "sc_best_ct",
              "prop_coloc_rho", "gwas_rna_convergence_tier",
              # scdrs_ct_mean_score is DESCRIPTIVE only (M3 remediation); the
              # per-cell p-value-derived cols (mc_pvalue/prop_sig/fdr/enriched/
              # n_sig_celltypes/best_celltype) were retired — see section 5.
              "scdrs_ct_mean_score",
              "scdrs_disease_shift", "scdrs_disease_p")
for (col in new_cols) {
  if (col %in% names(atlas)) {
    n_non_na <- sum(!is.na(atlas[[col]]))
    cat(sprintf("  %-30s: %d non-NA values\n", col, n_non_na))
  }
}

# Top convergent genes
cat("\nTop 20 convergent genes (Strong tier):\n")
strong <- atlas[gwas_rna_convergence_tier == "Strong"]
if (nrow(strong) > 0) {
  # Sort by INTACT score if available
  if ("intact_score_bulk" %in% names(strong)) {
    setorder(strong, -intact_score_bulk)
  }
  print(strong[1:min(20, nrow(strong)),
    .(gene, bulk_logFC, get(coloc_col), intact_score_bulk,
      mediation_proportion, n_evidence_gwas_rna)])
}

# ===========================================================================
# 9. Save
# ===========================================================================
cat("\n=== Saving integrated table ===\n")
fwrite(atlas, file.path(outdir, "gwas_rna_scrna_integrated.csv"))
cat("  Saved gwas_rna_scrna_integrated.csv:", nrow(atlas), "genes x", ncol(atlas), "cols\n")

# Convergence summary
convergence_summary <- atlas[n_evidence_gwas_rna >= 1, .(
  gene,
  gwas_rna_convergence_tier,
  n_evidence_gwas_rna
)]
if ("bulk_logFC" %in% names(atlas)) convergence_summary <- merge(convergence_summary, atlas[, .(gene, bulk_logFC, bulk_padj)], by = "gene")
if (!is.null(coloc_col) && coloc_col %in% names(atlas)) {
  convergence_summary <- merge(convergence_summary, atlas[, .(gene, coloc_pp4 = get(coloc_col))], by = "gene")
}
if ("intact_score_bulk" %in% names(atlas)) {
  convergence_summary <- merge(convergence_summary, atlas[, .(gene, intact_score_bulk)], by = "gene")
}

setorder(convergence_summary, -n_evidence_gwas_rna)
fwrite(convergence_summary, file.path(outdir, "convergence_summary.csv"))
cat("  Saved convergence_summary.csv:", nrow(convergence_summary), "rows\n")

cat("\nDone.\n")
