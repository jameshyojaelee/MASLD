#!/usr/bin/env Rscript
# 213_spatial_coloc.R — Zonation-stratified COLOC enrichment + spatial SVG × causal overlap
#
# Part of Module C (Spatial Genetic Architecture) in the Stratified Causal Pipeline.
#
# Analyses:
#   1. COLOC gene enrichment among periportal vs pericentral genes (Fisher's exact)
#   2. Spatial SVGs × COLOC overlap: Are spatially variable genes also genetically causal?
#   3. GWAS-ATAC × zonation: Regulatory variants in hepatocyte peaks × zonation-specific genes
#   4. Zone-specific causal scores: PP.H4 weighted by zonation strength
#
# Inputs:
#   - GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv          (19,030 genes)
#   - RNA-seq/results/zonation/deg_zonation_classification.csv           (14,080 genes; 24 PP, 21 PC)
#   - Analysis/Spatial/results/svg/differential_svgs.csv                 (3,000 genes; 150+ SVGs)
#   - Analysis/Spatial/results/svg/svgs_Healthy.csv                      (per-condition SVGs)
#   - Analysis/Spatial/results/svg/svgs_Steatotic.csv                    (per-condition SVGs)
#   - GWAS/finemapping/results/gwas_atac/gwas_atac_variant_annotation.csv (3,176 variant-peak overlaps)
#
# Outputs:
#   - RNA-seq/results/stratified_causal/spatial_coloc_enrichment.csv
#   - RNA-seq/results/stratified_causal/zone_causal_scores.csv
#
# SLURM: cpu partition, 4 CPUs, 32GB, 48h

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

outdir <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# 1. Load COLOC gene-level results
# ===========================================================================
cat("Loading COLOC gene-level results...\n")
coloc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
# Remove rows with empty gene name
coloc <- coloc[gene != "" & !is.na(gene)]
# T0.4 (2026-04-22): prefer SuSiE PP4 over legacy ABF.
if ("coloc_best_susie_pp4" %in% names(coloc)) {
  coloc[, coloc_best_abf_pp4 := coloc_best_pp4]
  coloc[, coloc_best_pp4 := fifelse(!is.na(coloc_best_susie_pp4),
                                    coloc_best_susie_pp4, coloc_best_pp4)]
  cat("  T0.4: coloc_best_pp4 now sourced from SuSiE with ABF fallback\n")
}
cat("  COLOC genes loaded:", nrow(coloc), "\n")
cat("  PP.H4 > 0.9:", sum(coloc$coloc_best_pp4 > 0.9), "\n")
cat("  PP.H4 > 0.8:", sum(coloc$coloc_best_pp4 > 0.8), "\n")

# ===========================================================================
# 2. Load zonation classification
# ===========================================================================
cat("\nLoading zonation classification...\n")
zonation <- fread(file.path(BASE, "RNA-seq/results/zonation/deg_zonation_classification.csv"))
cat("  Zonation genes:", nrow(zonation), "\n")
# Guard (2026-06-01, F-213-P0): the regenerated 2026-05-29
# deg_zonation_classification.csv no longer carries `spatial_zonation_class`
# (43_zonation_classification.R only joins it conditionally, and that join did
# not fire at the May-29 re-run). Downstream code (lines below + zone_direction
# fcase + final print) references this column, so re-running 213 against the
# current input previously crashed with "object spatial_zonation_class not
# found". Create it as all-NA when absent so the Visium-derived spatial arm
# degrades to empty (no spatial PP/PC genes) instead of erroring; the
# reference-marker zonation arm (zonation_class) is unaffected.
if (!("spatial_zonation_class" %in% names(zonation))) {
  warning("spatial_zonation_class missing from zonation input; ",
          "spatial-zonation arm disabled (all-NA). Re-run 43 to restore it.")
  zonation[, spatial_zonation_class := NA_character_]
}
cat("  Zonation classes:\n")
print(table(zonation$zonation_class))
cat("  Spatial zonation classes:\n")
print(table(zonation$spatial_zonation_class, useNA = "ifany"))

# ===========================================================================
# 3. Load spatial SVG results
# ===========================================================================
cat("\nLoading spatial SVG results...\n")
diff_svgs <- fread(file.path(BASE, "Analysis/Spatial/results/svg/differential_svgs.csv"))
# The first column is the gene name (unnamed in CSV); rename it
if (names(diff_svgs)[1] == "V1") setnames(diff_svgs, "V1", "gene")
svg_healthy <- fread(file.path(BASE, "Analysis/Spatial/results/svg/svgs_Healthy.csv"))
if (names(svg_healthy)[1] == "V1") setnames(svg_healthy, "V1", "gene")
svg_steatotic <- fread(file.path(BASE, "Analysis/Spatial/results/svg/svgs_Steatotic.csv"))
if (names(svg_steatotic)[1] == "V1") setnames(svg_steatotic, "V1", "gene")

# Define SVG gene sets
svg_any <- diff_svgs[svg_healthy == TRUE | svg_masld == TRUE, gene]
svg_disease_emergent <- diff_svgs[category == "disease_emergent_SVG", gene]
svg_disease_lost <- diff_svgs[category == "disease_lost_SVG", gene]
svg_stable <- diff_svgs[category == "stable" & (svg_healthy == TRUE | svg_masld == TRUE), gene]

cat("  Total genes in differential SVG table:", nrow(diff_svgs), "\n")
cat("  SVGs (any condition):", length(svg_any), "\n")
cat("  Disease-emergent SVGs:", length(svg_disease_emergent), "\n")
cat("  Disease-lost SVGs:", length(svg_disease_lost), "\n")
cat("  Stable SVGs:", length(svg_stable), "\n")

# ===========================================================================
# 4. Load GWAS-ATAC variant annotation (hepatocyte rows)
# ===========================================================================
cat("\nLoading GWAS-ATAC variant annotation...\n")
gwas_atac <- fread(file.path(BASE, "GWAS/finemapping/results/gwas_atac/gwas_atac_variant_annotation.csv"))
cat("  Total variant-peak overlaps:", nrow(gwas_atac), "\n")
cat("  Hepatocyte overlaps:", sum(gwas_atac$cell_type == "Hepatocytes"), "\n")

# Get genes linked to hepatocyte ATAC peaks
hep_atac <- gwas_atac[cell_type == "Hepatocytes"]
# Use nearest_gene as the gene annotation (most reliable)
hep_atac_genes <- unique(hep_atac$nearest_gene)
hep_atac_genes <- hep_atac_genes[hep_atac_genes != "" & !is.na(hep_atac_genes)]
cat("  Unique genes near hepatocyte ATAC peaks:", length(hep_atac_genes), "\n")

# ===========================================================================
# 5. Build universe (intersection of genes in COLOC + zonation)
# ===========================================================================
cat("\nBuilding gene universes...\n")
universe_coloc_zonation <- intersect(coloc$gene, zonation$human_symbol)
cat("  COLOC-zonation universe:", length(universe_coloc_zonation), "\n")

universe_coloc_svg <- intersect(coloc$gene, diff_svgs$gene)
cat("  COLOC-SVG universe:", length(universe_coloc_svg), "\n")

# COLOC genes at PP.H4 > 0.9
coloc_sig <- coloc[coloc_best_pp4 > 0.9, gene]
coloc_sig_08 <- coloc[coloc_best_pp4 > 0.8, gene]

# ===========================================================================
# 6. Fisher's exact enrichment: COLOC × zonation
# ===========================================================================
cat("\n=== COLOC × Zonation Enrichment (Fisher's exact) ===\n")

# Helper function for Fisher's exact test
fisher_enrichment <- function(coloc_genes, test_genes, universe, label) {
  coloc_in_u <- intersect(coloc_genes, universe)
  test_in_u <- intersect(test_genes, universe)

  a <- length(intersect(coloc_in_u, test_in_u))     # both
  b <- length(setdiff(coloc_in_u, test_in_u))        # COLOC only
  c <- length(setdiff(test_in_u, coloc_in_u))        # test only
  d <- length(universe) - a - b - c                   # neither

  mat <- matrix(c(a, b, c, d), nrow = 2)
  ft <- fisher.test(mat)

  overlap_genes <- intersect(coloc_in_u, test_in_u)

  data.table(
    test_label = label,
    n_coloc_in_universe = length(coloc_in_u),
    n_test_in_universe = length(test_in_u),
    n_overlap = a,
    odds_ratio = ft$estimate,
    pvalue = ft$p.value,
    ci_lower = ft$conf.int[1],
    ci_upper = ft$conf.int[2],
    overlap_genes = paste(sort(overlap_genes), collapse = ";")
  )
}

# Define zonation gene sets
pp_genes <- zonation[zonation_class == "Periportal", human_symbol]
pc_genes <- zonation[zonation_class == "Pericentral", human_symbol]
spatial_pp_genes <- zonation[spatial_zonation_class == "Periportal-enriched", human_symbol]
spatial_pc_genes <- zonation[spatial_zonation_class == "Pericentral-enriched", human_symbol]

enrichment_results <- rbindlist(list(
  # Zonation class (from reference markers)
  # Note: coloc_sig is filtered at PP4>0.9 (line 115), coloc_sig_08 at PP4>0.8 (line 116).
  # Labels relabeled 2026-04-22 (T2.11) to match the underlying filter thresholds.
  fisher_enrichment(coloc_sig, pp_genes, universe_coloc_zonation, "COLOC_PP4>0.9 × Periportal"),
  fisher_enrichment(coloc_sig, pc_genes, universe_coloc_zonation, "COLOC_PP4>0.9 × Pericentral"),
  fisher_enrichment(coloc_sig_08, pp_genes, universe_coloc_zonation, "COLOC_PP4>0.8 × Periportal"),
  fisher_enrichment(coloc_sig_08, pc_genes, universe_coloc_zonation, "COLOC_PP4>0.8 × Pericentral"),

  # Spatial zonation class (from Visium data)
  fisher_enrichment(coloc_sig, spatial_pp_genes, universe_coloc_zonation, "COLOC_PP4>0.9 × Spatial_Periportal"),
  fisher_enrichment(coloc_sig, spatial_pc_genes, universe_coloc_zonation, "COLOC_PP4>0.9 × Spatial_Pericentral"),
  fisher_enrichment(coloc_sig_08, spatial_pp_genes, universe_coloc_zonation, "COLOC_PP4>0.8 × Spatial_Periportal"),
  fisher_enrichment(coloc_sig_08, spatial_pc_genes, universe_coloc_zonation, "COLOC_PP4>0.8 × Spatial_Pericentral")
))

# BH correction across all tests
enrichment_results[, padj := p.adjust(pvalue, method = "BH")]

cat("\nEnrichment results:\n")
print(enrichment_results[, .(test_label, n_overlap, odds_ratio = round(odds_ratio, 2),
                              pvalue = signif(pvalue, 3), padj = signif(padj, 3))])

# ===========================================================================
# 7. Spatial SVGs × COLOC overlap
# ===========================================================================
cat("\n=== Spatial SVGs × COLOC Overlap ===\n")

# T2.11 fix (2026-04-22): labels relabeled "PP4>0.5" → "PP4>0.9" to match
# the underlying filter (coloc_sig is filtered at PP.H4 > 0.9 on line 115).
svg_enrichment <- rbindlist(list(
  fisher_enrichment(coloc_sig, svg_any, universe_coloc_svg, "COLOC_PP4>0.9 × SVG_any"),
  fisher_enrichment(coloc_sig, svg_disease_emergent, universe_coloc_svg, "COLOC_PP4>0.9 × SVG_disease_emergent"),
  fisher_enrichment(coloc_sig, svg_disease_lost, universe_coloc_svg, "COLOC_PP4>0.9 × SVG_disease_lost"),
  fisher_enrichment(coloc_sig, svg_stable, universe_coloc_svg, "COLOC_PP4>0.9 × SVG_stable"),
  fisher_enrichment(coloc_sig_08, svg_any, universe_coloc_svg, "COLOC_PP4>0.8 × SVG_any"),
  fisher_enrichment(coloc_sig_08, svg_disease_emergent, universe_coloc_svg, "COLOC_PP4>0.8 × SVG_disease_emergent")
))

svg_enrichment[, padj := p.adjust(pvalue, method = "BH")]

cat("\nSVG enrichment results:\n")
print(svg_enrichment[, .(test_label, n_overlap, odds_ratio = round(odds_ratio, 2),
                          pvalue = signif(pvalue, 3), padj = signif(padj, 3))])

# ===========================================================================
# 8. GWAS-ATAC × zonation: hepatocyte peaks in zonation-specific genes
# ===========================================================================
cat("\n=== GWAS-ATAC Hepatocyte Peaks × Zonation ===\n")

# Universe for this test: genes that appear in both zonation and GWAS-ATAC
universe_atac_zonation <- intersect(hep_atac_genes, zonation$human_symbol)
cat("  GWAS-ATAC × zonation universe:", length(universe_atac_zonation), "\n")

atac_zonation_enrichment <- rbindlist(list(
  fisher_enrichment(hep_atac_genes, pp_genes,
                    union(hep_atac_genes, zonation$human_symbol),
                    "Hep_ATAC_peaks × Periportal"),
  fisher_enrichment(hep_atac_genes, pc_genes,
                    union(hep_atac_genes, zonation$human_symbol),
                    "Hep_ATAC_peaks × Pericentral"),
  fisher_enrichment(hep_atac_genes, spatial_pp_genes,
                    union(hep_atac_genes, zonation$human_symbol),
                    "Hep_ATAC_peaks × Spatial_Periportal"),
  fisher_enrichment(hep_atac_genes, spatial_pc_genes,
                    union(hep_atac_genes, zonation$human_symbol),
                    "Hep_ATAC_peaks × Spatial_Pericentral")
))

atac_zonation_enrichment[, padj := p.adjust(pvalue, method = "BH")]

cat("\nGWAS-ATAC × zonation enrichment:\n")
print(atac_zonation_enrichment[, .(test_label, n_overlap, odds_ratio = round(odds_ratio, 2),
                                    pvalue = signif(pvalue, 3), padj = signif(padj, 3))])

# ===========================================================================
# 9. Combine all enrichment results
# ===========================================================================
all_enrichment <- rbindlist(list(
  enrichment_results[, module := "COLOC_zonation"],
  svg_enrichment[, module := "SVG_COLOC"],
  atac_zonation_enrichment[, module := "GWAS_ATAC_zonation"]
))

fwrite(all_enrichment, file.path(outdir, "spatial_coloc_enrichment.csv"))
cat("\nSaved spatial_coloc_enrichment.csv:", nrow(all_enrichment), "tests\n")

# ===========================================================================
# 10. Zone-specific causal scores
# ===========================================================================
cat("\n=== Computing zone-specific causal scores ===\n")

# Merge COLOC PP.H4 onto zonation table
zone_scores <- merge(zonation, coloc[, .(gene, coloc_best_pp4, coloc_best_gwas,
                                          coloc_n_gwas_h4_05, coloc_n_gwas_h4_08)],
                     by.x = "human_symbol", by.y = "gene", all.x = TRUE)

# Add SVG status
diff_svgs_summary <- diff_svgs[, .(gene, svg_healthy, svg_masld, delta_I, category)]
zone_scores <- merge(zone_scores, diff_svgs_summary,
                     by.x = "human_symbol", by.y = "gene", all.x = TRUE)

# Add GWAS-ATAC hepatocyte peak status
zone_scores[, in_hep_atac_peak := human_symbol %in% hep_atac_genes]

# Compute zone-causal score: PP.H4 weighted by zonation specificity
# Higher score = stronger genetic causality AND zonation specificity
zone_scores[, pp4 := fifelse(is.na(coloc_best_pp4), 0, coloc_best_pp4)]

# Periportal score: positive = periportal
zone_scores[, zone_direction := fcase(
  zonation_class == "Periportal", 1,
  zonation_class == "Pericentral", -1,
  spatial_zonation_class == "Periportal-enriched", 0.5,
  spatial_zonation_class == "Pericentral-enriched", -0.5,
  default = 0
)]

# Zone-causal composite: |zone_direction| * PP.H4
zone_scores[, zone_causal_score := abs(zone_direction) * pp4]

# Classify into zone-causal categories
zone_scores[, zone_coloc_class := fcase(
  pp4 > 0.9 & zone_direction > 0, "Periportal_causal",
  pp4 > 0.9 & zone_direction < 0, "Pericentral_causal",
  pp4 > 0.9 & zone_direction == 0, "Pan_lobular_causal",
  pp4 <= 0.5, "Not_causal"
)]

# Add SVG-causal flag
zone_scores[, svg_causal := pp4 > 0.9 & (svg_healthy == TRUE | svg_masld == TRUE)]

cat("\nZone-causal classification:\n")
print(table(zone_scores$zone_coloc_class))
cat("\nSVG-causal genes:", sum(zone_scores$svg_causal, na.rm = TRUE), "\n")

# Clean up temporary column
zone_scores[, pp4 := NULL]

# Sort by zone_causal_score descending
setorder(zone_scores, -zone_causal_score)

fwrite(zone_scores, file.path(outdir, "zone_causal_scores.csv"))
cat("Saved zone_causal_scores.csv:", nrow(zone_scores), "genes\n")

# ===========================================================================
# 11. Summary statistics
# ===========================================================================
cat("\n=== Summary ===\n")
cat("Total enrichment tests:", nrow(all_enrichment), "\n")
cat("Significant (padj < 0.05):", sum(all_enrichment$padj < 0.05), "\n")
cat("Zone-causal genes (PP.H4 > 0.9 + zonation):\n")
cat("  Periportal causal:", sum(zone_scores$zone_coloc_class == "Periportal_causal"), "\n")
cat("  Pericentral causal:", sum(zone_scores$zone_coloc_class == "Pericentral_causal"), "\n")
cat("  Pan-lobular causal:", sum(zone_scores$zone_coloc_class == "Pan_lobular_causal"), "\n")
cat("  SVG + causal:", sum(zone_scores$svg_causal, na.rm = TRUE), "\n")
cat("  In hepatocyte ATAC peak + causal:",
    sum(zone_scores$in_hep_atac_peak & zone_scores$zone_coloc_class != "Not_causal"), "\n")

# Top zone-causal genes
cat("\nTop 20 zone-causal genes:\n")
top_zone <- zone_scores[zone_coloc_class != "Not_causal"][order(-zone_causal_score)][1:min(20, .N)]
print(top_zone[, .(human_symbol, zonation_class, spatial_zonation_class,
                     coloc_best_pp4 = round(coloc_best_pp4, 3),
                     zone_coloc_class, zone_causal_score = round(zone_causal_score, 3),
                     in_hep_atac_peak, svg_causal)])

cat("\nScript 213 complete.\n")
