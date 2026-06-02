#!/usr/bin/env Rscript
# 211_progression_driver_genetics.R — Annotate progression drivers with genetic evidence
#
# For each top progression driver, annotate with:
#   - COLOC PP.H4 (best across 24 GWAS, and per-phenotype class)
#   - TWAS z-score (best across multi-GWAS)
#   - GWAS-ATAC motif disruption (if any)
#
# "Genetically validated progression drivers" = high progression tau AND COLOC PP.H4 > 0.9
# Progression causal score = PP4 x tau x |transition_logFC|
# Drug target classification: onset vs progression based on GWAS phenotype colocalization
#
# Inputs:
#   - Progression drivers: .../results/progression/driver_summary.csv (top 350)
#   - Full driver scores: .../results/progression/driver_scores.csv (all genes x transitions)
#   - Transition tau: .../results/progression/transition_tau_index.csv
#   - Transition dream results: .../results/progression/transition_fib_dream_results.csv
#   - COLOC gene-level: GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv
#   - COLOC full: GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv
#   - TWAS: RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv
#   - GWAS-ATAC motif: GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv
#   - GWAS-ATAC gene-level: GWAS/finemapping/results/gwas_atac/gene_level_gwas_atac.csv
#   - Drug validation: RNA-seq/results/drug_repurposing/clinical_drug_validation_table.csv
#
# Outputs:
#   - RNA-seq/results/stratified_causal/progression_driver_genetics.csv
#   - RNA-seq/results/stratified_causal/drug_target_progression_classification.csv
#
# SLURM: cpu partition, 4 CPUs, 32GB, 48h

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

outdir <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

prog_base <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression")

# ===========================================================================
# 1. Load progression driver data
# ===========================================================================
cat("=== Loading progression driver data ===\n")

# Top 350 drivers (ranked by composite driver_score)
drivers <- fread(file.path(prog_base, "driver_summary.csv"))
cat("  Driver summary:", nrow(drivers), "top drivers\n")

# Full driver scores (all genes x transitions) — use for broader annotation
# This is large (~238K rows), load selectively
driver_scores <- fread(file.path(prog_base, "driver_scores.csv"))
cat("  Full driver scores:", nrow(driver_scores), "gene-transition pairs\n")

# Transition tau (stage-specificity index per gene)
transition_tau <- fread(file.path(prog_base, "transition_tau_index.csv"))
cat("  Transition tau:", nrow(transition_tau), "genes\n")

# Transition fibrosis dream results (per-transition logFC)
trans_fib <- fread(file.path(prog_base, "transition_fib_dream_results.csv"))
cat("  Transition fibrosis dream:", nrow(trans_fib), "rows,",
    uniqueN(trans_fib$transition), "transitions\n")

# ===========================================================================
# 2. Load COLOC data
# ===========================================================================
cat("\n=== Loading COLOC data ===\n")

# Gene-level summary (best PP4 across all GWAS)
coloc_gene <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
# Remove the row with empty gene name if present
coloc_gene <- coloc_gene[gene != "" & !is.na(gene)]
# T0.4 (2026-04-22): prefer SuSiE PP4 over legacy ABF.
if ("coloc_best_susie_pp4" %in% names(coloc_gene)) {
  coloc_gene[, coloc_best_abf_pp4 := coloc_best_pp4]
  coloc_gene[, coloc_best_pp4 := fifelse(!is.na(coloc_best_susie_pp4),
                                         coloc_best_susie_pp4, coloc_best_pp4)]
  cat("  T0.4: coloc_best_pp4 now sourced from SuSiE with ABF fallback\n")
}
cat("  COLOC gene-level:", nrow(coloc_gene), "genes;",
    sum(coloc_gene$coloc_best_pp4 > 0.9), "with PP.H4 > 0.9\n")

# Full per-GWAS COLOC (for phenotype-specific PP4)
cat("  Loading full per-GWAS COLOC (436K rows)...\n")
coloc_full <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
cat("  Full COLOC:", nrow(coloc_full), "gene-GWAS pairs\n")

# Classify GWAS into phenotype categories for onset vs progression
# ARCHIVED 2026-04-09: Whitfield 2023 (36653562) — provenance unverified, Cirrhosis/HCC duplicate Ghouse cases
# Removed: 2023_36653562_NAFLD_EUR, 2023_36653562_NASH_EUR, 2023_36653562_Cirrhosis_EUR,
#           2023_36653562_HCC_EUR, 2023_36653562_Obesity_EUR
# ARCHIVED 2026-04-09: Pazoki_PDFF removed — duplicate of 2022_36402844_PDFF_EUR (same GCST90267352)
gwas_class <- data.table(
  gwas_name = c(
    # Disease onset GWAS (NAFLD diagnosis)
    "2019_31311600_NAFLD_EUR", "2020_32298765_NAFLD_EUR",
    "2021_34841290_NAFLD_EUR",
    "2023_36280732_NAFLD_deCode_EUR", "2023_36280732_NAFLD_Intermountain_EUR",
    "2023_36280732_NAFLD_UKBB_EUR",
    "FinnGen_NAFLD", "FinnGen_NASH",  # RESTORED 2026-04-09: verified FinnGen R12 provenance
    # Steatosis GWAS (PDFF — early/metabolic)
    "2021_34128465_PDFF_EUR", "2021_34957434_PDFF_EUR",
    "2022_36402844_PDFF_EUR",
    # Progression GWAS (cirrhosis, HCC)
    "Ghouse_Cirrhosis",
    "Ghouse_HCC",
    "FinnGen_HCC",  # RESTORED 2026-04-09
    # Liver enzyme GWAS (indirect/continuous)
    "UKBB_ALT", "UKBB_AST", "UKBB_GGT"
  ),
  gwas_category = c(
    rep("onset", 8),
    rep("steatosis", 3),
    rep("progression", 3),
    rep("enzyme", 3)
  )
)

# Compute per-gene per-category best PP4
coloc_by_cat <- merge(coloc_full[, .(gene, gwas_name, PP.H4.abf)],
                      gwas_class, by = "gwas_name", all.x = TRUE)
# Fill any unclassified GWAS as "other"
coloc_by_cat[is.na(gwas_category), gwas_category := "other"]

cat("  GWAS category distribution:\n")
print(coloc_by_cat[, .N, by = gwas_category][order(-N)])

# Best PP4 per gene per category
pp4_by_cat <- coloc_by_cat[, .(
  pp4_best = max(PP.H4.abf, na.rm = TRUE),
  best_gwas = gwas_name[which.max(PP.H4.abf)],
  n_gwas_h4_05 = sum(PP.H4.abf > 0.5, na.rm = TRUE)
), by = .(gene, gwas_category)]

# Pivot to wide format
pp4_wide <- dcast(pp4_by_cat, gene ~ gwas_category,
                  value.var = c("pp4_best", "n_gwas_h4_05"),
                  fill = 0)
cat("  Per-category PP4 computed for", nrow(pp4_wide), "genes\n")

# ===========================================================================
# 3. Load TWAS data
# ===========================================================================
cat("\n=== Loading TWAS data ===\n")
twas <- fread(file.path(BASE,
  "RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv"))

# Best TWAS per gene (minimum p-value across GWAS)
twas_best <- twas[, .(
  twas_best_z     = zscore[which.min(pvalue)],
  twas_best_p     = min(pvalue, na.rm = TRUE),
  twas_best_gwas  = gwas[which.min(pvalue)],
  twas_n_sig_gwas = sum(pvalue < 0.05, na.rm = TRUE)
), by = .(gene_symbol = gene_name)]
cat("  TWAS best per gene:", nrow(twas_best), "genes;",
    sum(twas_best$twas_best_p < 0.05), "nominally significant\n")

# ===========================================================================
# 4. Load GWAS-ATAC motif disruption data
# ===========================================================================
cat("\n=== Loading GWAS-ATAC motif disruption data ===\n")

motif <- fread(file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv"))
cat("  Motif disruptions:", nrow(motif), "variant-TF pairs\n")

# Gene-level GWAS-ATAC summary
gwas_atac_gene <- fread(file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/gene_level_gwas_atac.csv"))
cat("  GWAS-ATAC gene-level:", nrow(gwas_atac_gene), "genes\n")

# Summarize motif disruption per gene (using gene-level file)
# The gene_level_gwas_atac.csv has: assigned_gene, gwas_variant_in_peak,
#   gwas_motif_disrupted, gwas_atac_regulatory_score, etc.
# Some assigned_gene values may be ensembl IDs — we use them as-is and merge by both
gwas_atac_gene[, gene_symbol := assigned_gene]

# Also summarize key disease-regulon TF disruptions per gene from raw motif data
# Count how many disease-regulon TFs are disrupted near each gene
disease_tf_disrupted <- motif[motif_in_disease_regulon == TRUE,
  .(n_disease_tf_disrupted = uniqueN(tf_name),
    disease_tfs_disrupted = paste(unique(tf_name), collapse = ";")),
  by = .(SNP_id)]
cat("  Variants disrupting disease regulon TFs:", nrow(disease_tf_disrupted), "\n")

# ===========================================================================
# 5. Build gene-symbol lookup from COLOC
# ===========================================================================
# COLOC uses 'gene' as symbol; driver_summary uses 'gene_symbol'
# transition_tau uses ensembl IDs ('gene')
# We need a mapping: driver_summary$gene_symbol -> coloc$gene

# The driver_summary already has gene_symbol. COLOC gene_level has gene (= symbol).
# transition_tau has ensembl IDs only — we can get symbols from driver_scores
# which has gene_symbol and matches transition_tau on ensembl via transition_programs

# For the main merge, use gene_symbol from drivers -> gene in coloc
cat("\n=== Merging progression drivers with genetic evidence ===\n")

# ===========================================================================
# 6. Annotate all driver_summary genes with COLOC
# ===========================================================================

# Get best overall PP4 from gene-level COLOC
drivers_ann <- merge(drivers, coloc_gene[, .(gene, coloc_best_pp4, coloc_best_gwas,
                                              coloc_n_gwas_h4_05, coloc_n_gwas_h4_08)],
                     by.x = "gene_symbol", by.y = "gene", all.x = TRUE)

# Merge per-category PP4
# Select only the columns we need from pp4_wide
pp4_cols <- names(pp4_wide)[names(pp4_wide) != "gene"]
drivers_ann <- merge(drivers_ann, pp4_wide,
                     by.x = "gene_symbol", by.y = "gene", all.x = TRUE)

# Fill NAs with 0 for PP4 columns
pp4_num_cols <- grep("^pp4_best_|^n_gwas_h4_05_", names(drivers_ann), value = TRUE)
for (col in pp4_num_cols) {
  set(drivers_ann, which(is.na(drivers_ann[[col]])), col, 0)
}

cat("  After COLOC merge:", nrow(drivers_ann), "drivers\n")
cat("    With PP.H4 > 0.9:", sum(drivers_ann$coloc_best_pp4 > 0.9, na.rm = TRUE), "\n")

# ===========================================================================
# 7. Annotate with TWAS
# ===========================================================================
drivers_ann <- merge(drivers_ann, twas_best,
                     by = "gene_symbol", all.x = TRUE)
cat("  After TWAS merge:", nrow(drivers_ann), "drivers\n")
cat("    With TWAS p < 0.05:", sum(drivers_ann$twas_best_p < 0.05, na.rm = TRUE), "\n")

# ===========================================================================
# 8. Annotate with GWAS-ATAC
# ===========================================================================
drivers_ann <- merge(drivers_ann,
  gwas_atac_gene[, .(gene_symbol,
                      gwas_variant_in_peak,
                      gwas_variant_cell_types,
                      gwas_n_variants,
                      gwas_max_pip_in_peak,
                      gwas_motif_disrupted,
                      gwas_atac_regulatory_score)],
  by = "gene_symbol", all.x = TRUE)

cat("  After GWAS-ATAC merge:", nrow(drivers_ann), "drivers\n")
cat("    With GWAS variant in peak:",
    sum(drivers_ann$gwas_variant_in_peak == TRUE, na.rm = TRUE), "\n")
cat("    With motif disruption:",
    sum(!is.na(drivers_ann$gwas_motif_disrupted) & drivers_ann$gwas_motif_disrupted != "",
        na.rm = TRUE), "\n")

# ===========================================================================
# 9. Compute progression_causal_score
# ===========================================================================
cat("\n=== Computing progression_causal_score ===\n")

# progression_causal_score = PP4 x tau x |transition_logFC|
# Use best PP4 (coloc_best_pp4), tau from driver_summary, trans_logFC from driver_summary
drivers_ann[, progression_causal_score := fifelse(
  is.na(coloc_best_pp4) | is.na(tau) | is.na(trans_logFC),
  0,
  coloc_best_pp4 * tau * abs(trans_logFC)
)]

cat("  Score distribution:\n")
cat("    Min:", min(drivers_ann$progression_causal_score, na.rm = TRUE), "\n")
cat("    Median:", median(drivers_ann$progression_causal_score, na.rm = TRUE), "\n")
cat("    Max:", max(drivers_ann$progression_causal_score, na.rm = TRUE), "\n")
cat("    > 0.1:", sum(drivers_ann$progression_causal_score > 0.1, na.rm = TRUE), "\n")

# ===========================================================================
# 10. Identify genetically validated progression drivers
# ===========================================================================
cat("\n=== Identifying genetically validated progression drivers ===\n")

# Genetically validated = high tau (> 0.5) AND COLOC PP.H4 > 0.9
drivers_ann[, genetically_validated := (
  !is.na(coloc_best_pp4) & coloc_best_pp4 > 0.9 &
  !is.na(tau) & tau > 0.5
)]

n_validated <- sum(drivers_ann$genetically_validated, na.rm = TRUE)
cat("  Genetically validated progression drivers:", n_validated,
    "out of", nrow(drivers_ann), "\n")

if (n_validated > 0) {
  top_validated <- drivers_ann[genetically_validated == TRUE][
    order(-progression_causal_score)][1:min(20, n_validated)]
  cat("\n  Top genetically validated drivers:\n")
  print(top_validated[, .(gene_symbol, transition, driver_rank, tau,
                          coloc_best_pp4, coloc_best_gwas,
                          twas_best_z, progression_causal_score)])
}

# ===========================================================================
# 11. Classify as onset vs progression target
# ===========================================================================
cat("\n=== Classifying onset vs progression targets ===\n")

# For each driver, check if COLOC is stronger in onset GWAS vs progression GWAS
drivers_ann[, onset_pp4 := pmax(pp4_best_onset, pp4_best_steatosis, na.rm = TRUE)]
drivers_ann[, progression_pp4 := pp4_best_progression]

# Classification logic
drivers_ann[, progression_class := fifelse(
  is.na(coloc_best_pp4) | coloc_best_pp4 < 0.1, "not_causal",
  fifelse(
    onset_pp4 > 0.9 & progression_pp4 > 0.9, "pan_stage",
    fifelse(
      onset_pp4 > 0.9, "onset",
      fifelse(
        progression_pp4 > 0.9, "progression",
        fifelse(
          pp4_best_enzyme > 0.9, "enzyme_associated",
          "low_confidence"
        )
      )
    )
  )
)]

cat("  Progression classification:\n")
print(drivers_ann[, .N, by = progression_class][order(-N)])

# ===========================================================================
# 12. Save progression_driver_genetics.csv
# ===========================================================================
cat("\n=== Saving progression_driver_genetics.csv ===\n")

# Order by progression_causal_score descending
setorder(drivers_ann, -progression_causal_score)

fwrite(drivers_ann,
       file.path(outdir, "progression_driver_genetics.csv"))
cat("  Written:", nrow(drivers_ann), "rows to",
    file.path(outdir, "progression_driver_genetics.csv"), "\n")

# ===========================================================================
# 13. Drug target progression classification
# ===========================================================================
cat("\n=== Drug target progression classification ===\n")

drug_table <- fread(file.path(BASE,
  "RNA-seq/results/drug_repurposing/clinical_drug_validation_table.csv"))
cat("  Drug validation table:", nrow(drug_table), "entries\n")

# Get unique drug targets
drug_targets <- unique(drug_table[, .(drug, target_gene, stage, moa)])
cat("  Unique drug-target pairs:", nrow(drug_targets), "\n")

# For each drug target gene, get per-phenotype-class COLOC PP4
drug_genetics <- merge(drug_targets, pp4_wide,
                       by.x = "target_gene", by.y = "gene", all.x = TRUE)

# Merge overall COLOC
drug_genetics <- merge(drug_genetics,
  coloc_gene[, .(gene, coloc_best_pp4, coloc_best_gwas,
                  coloc_n_gwas_h4_05, coloc_n_gwas_h4_08)],
  by.x = "target_gene", by.y = "gene", all.x = TRUE)

# Merge TWAS
drug_genetics <- merge(drug_genetics, twas_best,
                       by.x = "target_gene", by.y = "gene_symbol", all.x = TRUE)

# Merge GWAS-ATAC
drug_genetics <- merge(drug_genetics,
  gwas_atac_gene[, .(gene_symbol, gwas_variant_in_peak,
                      gwas_motif_disrupted, gwas_atac_regulatory_score)],
  by.x = "target_gene", by.y = "gene_symbol", all.x = TRUE)

# Fill NAs in PP4 columns
pp4_drug_cols <- grep("^pp4_best_|^n_gwas_h4_05_", names(drug_genetics), value = TRUE)
for (col in pp4_drug_cols) {
  set(drug_genetics, which(is.na(drug_genetics[[col]])), col, 0)
}

# Classify onset vs progression for drug targets
drug_genetics[, onset_pp4 := pmax(pp4_best_onset, pp4_best_steatosis, na.rm = TRUE)]
drug_genetics[, progression_pp4 := pp4_best_progression]

drug_genetics[, target_class := fifelse(
  is.na(coloc_best_pp4) | coloc_best_pp4 < 0.1, "not_genetically_causal",
  fifelse(
    onset_pp4 > 0.9 & progression_pp4 > 0.9, "pan_stage_target",
    fifelse(
      onset_pp4 > 0.9, "onset_target",
      fifelse(
        progression_pp4 > 0.9, "progression_target",
        fifelse(
          pp4_best_enzyme > 0.9, "enzyme_associated_target",
          "low_confidence_causal"
        )
      )
    )
  )
)]

# Also check per-transition driver score for drug targets
# Merge driver info for drug target genes (best transition)
drug_driver_info <- driver_scores[gene_symbol %in% drug_targets$target_gene,
  .(peak_driver_score = max(driver_score, na.rm = TRUE),
    peak_driver_transition = transition[which.max(driver_score)],
    peak_driver_rank = driver_rank[which.max(driver_score)]),
  by = .(gene_symbol)]

drug_genetics <- merge(drug_genetics, drug_driver_info,
                       by.x = "target_gene", by.y = "gene_symbol", all.x = TRUE)

# Compute progression_causal_score for drug targets
# Use tau from transition_tau if available
# First get tau for drug target genes — from driver_scores (already has tau)
drug_tau <- driver_scores[gene_symbol %in% drug_targets$target_gene,
  .(tau = max(tau, na.rm = TRUE),
    peak_logFC = trans_logFC[which.max(driver_score)]),
  by = .(gene_symbol)]
# Handle -Inf from empty max
drug_tau[is.infinite(tau), tau := NA_real_]

drug_genetics <- merge(drug_genetics, drug_tau,
                       by.x = "target_gene", by.y = "gene_symbol", all.x = TRUE)

drug_genetics[, progression_causal_score := fifelse(
  is.na(coloc_best_pp4) | is.na(tau) | is.na(peak_logFC), 0,
  coloc_best_pp4 * tau * abs(peak_logFC)
)]

cat("\n  Drug target classification:\n")
print(drug_genetics[, .(drug, target_gene, stage, target_class, coloc_best_pp4,
                         onset_pp4, progression_pp4, pp4_best_enzyme,
                         progression_causal_score)])

# Save
setorder(drug_genetics, -progression_causal_score)
fwrite(drug_genetics,
       file.path(outdir, "drug_target_progression_classification.csv"))
cat("\n  Written:", nrow(drug_genetics), "rows to",
    file.path(outdir, "drug_target_progression_classification.csv"), "\n")

# ===========================================================================
# 14. Summary statistics
# ===========================================================================
cat("\n============================================\n")
cat("=== SUMMARY ===\n")
cat("============================================\n")
cat("Progression drivers annotated:", nrow(drivers_ann), "\n")
cat("  With COLOC PP.H4 > 0.9:", sum(drivers_ann$coloc_best_pp4 > 0.9, na.rm = TRUE), "\n")
cat("  With COLOC PP.H4 > 0.8:", sum(drivers_ann$coloc_best_pp4 > 0.8, na.rm = TRUE), "\n")
cat("  With TWAS p < 0.05:", sum(drivers_ann$twas_best_p < 0.05, na.rm = TRUE), "\n")
cat("  With GWAS-ATAC variant:", sum(drivers_ann$gwas_variant_in_peak == TRUE, na.rm = TRUE), "\n")
cat("  Genetically validated (tau>0.5 & PP4>0.9):", n_validated, "\n")
cat("\nProgression classification of drivers:\n")
print(drivers_ann[, .N, by = progression_class][order(-N)])
cat("\nDrug target classifications:\n")
print(drug_genetics[, .(drug, target_gene, target_class, coloc_best_pp4,
                         progression_causal_score)][order(-progression_causal_score)])
cat("\nOutputs:\n")
cat("  ", file.path(outdir, "progression_driver_genetics.csv"), "\n")
cat("  ", file.path(outdir, "drug_target_progression_classification.csv"), "\n")
cat("\nDone.\n")
