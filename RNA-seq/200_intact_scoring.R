#!/usr/bin/env Rscript
# 200_intact_scoring.R — Multi-INTACT composite scoring (bulk TWAS + scTWAS × COLOC)
#
# Uses multi_intact() (Okamoto et al., AJHG 2024 extension) to jointly
# integrate bulk TWAS + scTWAS z-scores with their respective COLOC PP.H4
# into a single causal probability per gene. The EM algorithm prevents
# BF saturation that occurs with single-product intact() at extreme z-scores.
#
# Product 1: Bulk TWAS (S-PrediXcan) + ABF COLOC PP.H4 (Broadaway eQTL)
# Product 2: scTWAS (sc-eQTL MR) + ABF COLOC PP.H4 (same, for the sc gene)
#
# Inputs:
#   - Bulk TWAS: RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv
#   - scTWAS:    RNA-seq/results/causal_inference/sceqtl_twas/sceqtl_twas_all_results.csv
#   - COLOC:     GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv
#
# Outputs:
#   - RNA-seq/results/gwas_rna_integration/intact_scores.csv
#   - RNA-seq/results/gwas_rna_integration/intact_bulk_per_gwas.csv
#   - RNA-seq/results/gwas_rna_integration/intact_celltype_per_ct.csv
#
# SLURM: cpu partition, 4 CPUs, 32GB, 48h

suppressPackageStartupMessages({
  library(data.table)
  library(INTACT)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

outdir <- file.path(BASE, "RNA-seq/results/gwas_rna_integration")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# 1. Load bulk TWAS
# ===========================================================================
cat("Loading bulk TWAS results...\n")
twas <- fread(file.path(BASE, "RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv"))
# Best z per gene (across GWAS)
twas_best <- twas[, .(bulk_z = zscore[which.min(pvalue)],
                       bulk_p = min(pvalue)),
                  by = .(gene = gene_name)]
cat("  Bulk TWAS:", nrow(twas_best), "genes\n")

# ===========================================================================
# 2. Load scTWAS (cell-type)
# ===========================================================================
cat("Loading scTWAS results...\n")
sceqtl <- fread(file.path(BASE,
  "RNA-seq/results/causal_inference/sceqtl_twas/sceqtl_twas_all_results.csv"))
sceqtl[, twas_z := b / se]
sceqtl <- sceqtl[!is.na(twas_z) & is.finite(twas_z)]
# Best z per gene (across cell types and GWAS)
sceqtl_best <- sceqtl[, .(sc_z = twas_z[which.min(pval)],
                           sc_p = min(pval),
                           best_ct = cell_type[which.min(pval)]),
                       by = .(gene = exposure)]
cat("  scTWAS:", nrow(sceqtl_best), "genes\n")

# ===========================================================================
# 3. Load COLOC PP.H4
# ===========================================================================
cat("Loading COLOC results...\n")
coloc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc <- coloc[gene != "" & !is.na(gene)]
cat("  COLOC:", nrow(coloc), "genes\n")

# ===========================================================================
# 4. Build multi_intact input: genes with BOTH bulk + scTWAS
# ===========================================================================
cat("\n=== Building Multi-INTACT input ===\n")

# Merge all three sources
mi_dt <- merge(twas_best, sceqtl_best, by = "gene")
mi_dt <- merge(mi_dt, coloc[, .(gene, GLCP = coloc_best_pp4)], by = "gene")
mi_dt <- mi_dt[!is.na(bulk_z) & !is.na(sc_z) & !is.na(GLCP)]

cat("  Genes with bulk TWAS + scTWAS + COLOC:", nrow(mi_dt), "\n")

# multi_intact requires:
#   gene, GLCP_1, GLCP_2, z_1, z_2, chisq
# DISCLOSURE (review A10#1): GLCP_2 reuses the bulk COLOC PP.H4 because no sc-eQTL COLOC
# table exists for the scTWAS product. With glcp_aggreg="max" this is non-inflating
# (max(x,x)=x); the per-gene score is driven by the product z-magnitudes (GPRP_1/GPRP_2).
# Treat the sc product as a sensitivity arm, not an independent COLOC line of evidence.
# chisq = z_1^2 + z_2^2 (approximation assuming independence)
mi_input <- data.frame(
  gene = mi_dt$gene,
  GLCP_1 = mi_dt$GLCP,     # COLOC for bulk TWAS product
  GLCP_2 = mi_dt$GLCP,     # scTWAS product — reuses bulk COLOC (no sc-eQTL COLOC; see disclosure)
  z_1 = mi_dt$bulk_z,      # bulk TWAS z
  z_2 = mi_dt$sc_z,        # scTWAS z
  chisq = mi_dt$bulk_z^2 + mi_dt$sc_z^2  # multivariate Wald
)

cat("  Input dimensions:", nrow(mi_input), "genes x", ncol(mi_input), "columns\n")
cat("  |z_1| range:", round(range(abs(mi_input$z_1)), 2), "\n")
cat("  |z_2| range:", round(range(abs(mi_input$z_2)), 2), "\n")
cat("  GLCP range:", round(range(mi_input$GLCP_1), 4), "\n")

# ===========================================================================
# 5. Run Multi-INTACT with EM algorithm
# ===========================================================================
cat("\n=== Running Multi-INTACT (EM algorithm) ===\n")

mi_result <- multi_intact(
  df = mi_input,
  prior_fun = linear,
  em_algorithm = TRUE,
  glcp_aggreg = "max",
  return_model_posteriors = TRUE
)

# Extract results from multi_intact return list
# Returns: list with [[1]] = data.frame (gene, GPPC, posterior_0/1/2/12, GPRP_1, GPRP_2)
#   GPPC = Gene Product Pair Causal probability (both products relevant)
#   GPRP_1 = Gene Product Relevance Probability for product 1 (bulk TWAS)
#   GPRP_2 = Gene Product Relevance Probability for product 2 (scTWAS)
mi_df <- mi_result[[1]]
cat("  EM converged:", mi_result[[3]], "\n")
cat("  Result columns:", paste(names(mi_df), collapse = ", "), "\n")
cat("  Prior parameter estimates:", paste(round(mi_result[[2]], 4), collapse = ", "), "\n")

# Match back to mi_dt by gene
mi_df <- as.data.table(mi_df)
mi_dt <- merge(mi_dt, mi_df[, .(gene, multi_intact_gppc = GPPC,
                                  multi_intact_gprp1 = GPRP_1,
                                  multi_intact_gprp2 = GPRP_2,
                                  posterior_0, posterior_12)],
               by = "gene", all.x = TRUE)

# Primary score: max of GPRP_1 and GPRP_2 (either product being causal)
# This is more informative than GPPC (which requires BOTH products)
mi_dt[, multi_intact_score := pmax(multi_intact_gprp1, multi_intact_gprp2, na.rm = TRUE)]

cat("\n=== Multi-INTACT distribution ===\n")
cat("  Exactly 1.0:", sum(mi_dt$multi_intact_score == 1, na.rm = TRUE), "\n")
cat("  > 0.9:", sum(mi_dt$multi_intact_score > 0.9, na.rm = TRUE), "\n")
cat("  > 0.5:", sum(mi_dt$multi_intact_score > 0.5, na.rm = TRUE), "\n")
cat("  > 0.1:", sum(mi_dt$multi_intact_score > 0.1, na.rm = TRUE), "\n")
cat("  = 0:", sum(mi_dt$multi_intact_score == 0, na.rm = TRUE), "\n")
print(quantile(mi_dt$multi_intact_score, probs = seq(0, 1, 0.1), na.rm = TRUE))

# ===========================================================================
# 6. Also run single-product INTACT per cell type (for cell-type resolution)
# ===========================================================================
cat("\n=== Running single INTACT per cell type (for cell-type assignment) ===\n")

cell_types <- unique(sceqtl$cell_type)
intact_ct_list <- list()

for (ct in cell_types) {
  ct_sub <- sceqtl[cell_type == ct]
  merged <- merge(
    ct_sub[, .(gene = exposure, twas_z, twas_p = pval, gwas)],
    coloc[, .(gene, GLCP = coloc_best_pp4)],
    by = "gene"
  )
  merged <- merged[!is.na(twas_z) & !is.na(GLCP)]

  if (nrow(merged) < 10) {
    cat("  Skipping", ct, "— only", nrow(merged), "genes\n")
    next
  }

  cat("  INTACT for", ct, ":", nrow(merged), "genes")

  scores <- intact(GLCP_vec = merged$GLCP, z_vec = merged$twas_z, prior_fun = linear)
  merged[, intact_score := scores]
  merged[, intact_rank := frank(-intact_score)]
  merged[, cell_type := ct]

  cat(" -> >0.5:", sum(scores > 0.5), "\n")
  intact_ct_list[[ct]] <- merged
}

intact_ct <- rbindlist(intact_ct_list, fill = TRUE)

# Best cell-type per gene
intact_ct_best <- intact_ct[, .(
  intact_score_ct = max(intact_score, na.rm = TRUE),
  intact_best_celltype = cell_type[which.max(intact_score)],
  intact_n_celltypes_05 = sum(intact_score > 0.5)
), by = gene]

# ===========================================================================
# 7. Merge Multi-INTACT + cell-type INTACT
# ===========================================================================
cat("\n=== Merging results ===\n")

intact_merged <- merge(
  mi_dt[, .(gene, multi_intact_score, bulk_z, sc_z, best_ct, coloc_pp4 = GLCP)],
  intact_ct_best,
  by = "gene", all = TRUE
)

# Provenance flag (review A10#5): distinguish genuine multi-product INTACT (bulk+sc,
# EM null-mixture) from single-product cell-type INTACT. Do NOT overwrite
# multi_intact_score with the weaker single-product intact_score_ct — that silently
# promoted 74 scTWAS-only genes into the Tier-1 "Genetic_validated" tier via the bulk
# gate. intact_score_ct remains available as a supplementary annotation but is never
# alone-sufficient for the 46d Tier-1 INTACT gate.
intact_merged[, intact_score_source := fcase(
  !is.na(multi_intact_score), "multi",
  !is.na(intact_score_ct),    "ct_single",
  default = NA_character_)]

# Legacy column name for downstream compatibility — multi-product score ONLY
# (NA for bulk-absent genes; ct-only genes carry their score in intact_score_ct).
intact_merged[, intact_score_bulk := multi_intact_score]

cat("  Total genes:", nrow(intact_merged), "\n")
cat("  Multi-INTACT > 0.5:", sum(intact_merged$multi_intact_score > 0.5, na.rm = TRUE), "\n")
cat("  Multi-INTACT > 0.9:", sum(intact_merged$multi_intact_score > 0.9, na.rm = TRUE), "\n")
cat("  Multi-INTACT = 1.0:", sum(intact_merged$multi_intact_score == 1, na.rm = TRUE), "\n")
cat("  Cell-type INTACT > 0.5:", sum(intact_merged$intact_score_ct > 0.5, na.rm = TRUE), "\n")

# ===========================================================================
# 8. Validation
# ===========================================================================
cat("\n=== Validation: Known MASLD genes ===\n")
known <- c("PNPLA3", "TM6SF2", "HSD17B13", "MBOAT7", "GCKR",
           "MARC1", "THRB", "GPAM", "MTARC1")

for (g in known) {
  row <- intact_merged[gene == g]
  if (nrow(row) > 0) {
    cat(sprintf("  %-12s  multi_INTACT=%.3f  ct_INTACT=%.3f  best_ct=%s  COLOC=%.3f\n",
                g,
                ifelse(is.na(row$multi_intact_score[1]), 0, row$multi_intact_score[1]),
                ifelse(is.na(row$intact_score_ct[1]), 0, row$intact_score_ct[1]),
                ifelse(is.na(row$intact_best_celltype[1]), "NA", row$intact_best_celltype[1]),
                ifelse(is.na(row$coloc_pp4[1]), 0, row$coloc_pp4[1])))
  } else {
    cat(sprintf("  %-12s  NOT FOUND\n", g))
  }
}

# ===========================================================================
# 9. Save
# ===========================================================================
cat("\n=== Saving ===\n")

fwrite(intact_merged, file.path(outdir, "intact_scores.csv"))
cat("  intact_scores.csv:", nrow(intact_merged), "rows\n")

fwrite(intact_ct, file.path(outdir, "intact_celltype_per_ct.csv"))
cat("  intact_celltype_per_ct.csv:", nrow(intact_ct), "rows\n")

# Save per-GWAS bulk for reference
twas_all <- twas[, .(gene = gene_name, twas_z = zscore, twas_p = pvalue, gwas)]
twas_coloc <- merge(twas_all, coloc[, .(gene, GLCP = coloc_best_pp4)], by = "gene")
twas_coloc <- twas_coloc[!is.na(twas_z) & !is.na(GLCP)]
bulk_scores <- twas_coloc[, {
  s <- intact(GLCP_vec = GLCP, z_vec = twas_z, prior_fun = linear)
  .(twas_z, twas_p, GLCP, intact_score = s, gwas)
}, by = gene]
fwrite(bulk_scores, file.path(outdir, "intact_bulk_per_gwas.csv"))
cat("  intact_bulk_per_gwas.csv:", nrow(bulk_scores), "rows\n")

cat("\nDone.\n")
