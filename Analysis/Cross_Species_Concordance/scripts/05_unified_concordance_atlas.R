#!/usr/bin/env Rscript
# 05_unified_concordance_atlas.R
# ---------------------------------------------------------------------------
# Integrate all concordance metrics into a unified per-gene atlas
# Generates translatability scores for library prioritization
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

cat("=== Phase 5: Unified Concordance Atlas ===\n\n")

WD  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance"
RES <- file.path(WD, "results")

BASE    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
H_INT   <- file.path(BASE, "Human/Patient_Cohorts/analysis/integration")
ANNOT   <- file.path(H_INT, "results/gene_annotation")
DS_DIR  <- file.path(H_INT, "results/disease_signatures")

ortho <- fread(file.path(ANNOT, "ortholog_mapping.tsv"))

# ============================================================
#  Load Phase 1: gene concordance
# ============================================================
cat("Loading Phase 1 results...\n")
gene_conc <- tryCatch(fread(file.path(RES, "gene_concordance_per_gene.csv")), error = function(e) NULL)
gene_matrix <- tryCatch(fread(file.path(RES, "gene_concordance_matrix_20x.csv")), error = function(e) NULL)

# ============================================================
#  Load Phase 2: pathway results
# ============================================================
cat("Loading Phase 2 results...\n")
fgsea_conc   <- tryCatch(fread(file.path(RES, "fgsea_pathway_concordance.csv")), error = function(e) NULL)
ssgsea_conc  <- tryCatch(fread(file.path(RES, "ssgsea_concordance.csv")), error = function(e) NULL)
ora_conc     <- tryCatch(fread(file.path(RES, "ora_term_concordance.csv")), error = function(e) NULL)

# ============================================================
#  Load Phase 3: network results
# ============================================================
cat("Loading Phase 3 results...\n")
wgcna_mods   <- tryCatch(fread(file.path(RES, "wgcna_module_assignments.csv")), error = function(e) NULL)
wgcna_pres   <- tryCatch(fread(file.path(RES, "wgcna_preservation_stats.csv")), error = function(e) NULL)
tf_human     <- tryCatch(fread(file.path(RES, "tf_activity_human.csv")), error = function(e) NULL)
tf_conc      <- tryCatch(fread(file.path(RES, "tf_concordance.csv")), error = function(e) NULL)

# ============================================================
#  Build unified per-gene table
# ============================================================
cat("\n=== Building Unified Atlas ===\n")

# Start with gene concordance as backbone
if (!is.null(gene_conc)) {
  atlas <- copy(gene_conc)
} else {
  cat("ERROR: gene_concordance_per_gene.csv not found, cannot build atlas\n")
  quit(status = 1)
}

# Add WGCNA module info
if (!is.null(wgcna_mods) && !is.null(wgcna_pres)) {
  atlas <- merge(atlas, wgcna_mods[, .(gene, wgcna_module = module)],
    by.x = "human_symbol", by.y = "gene", all.x = TRUE)

  # Get average preservation per module
  if (nrow(wgcna_pres) > 0) {
    mod_pres <- wgcna_pres[, .(
      mean_Zsummary = mean(Zsummary, na.rm = TRUE),
      n_preserved = sum(preservation == "Highly_Preserved")
    ), by = module]
    atlas <- merge(atlas, mod_pres[, .(module, mean_Zsummary, n_preserved)],
      by.x = "wgcna_module", by.y = "module", all.x = TRUE)
  }
  cat("  Added WGCNA module info\n")
}

# ============================================================
#  Compute Translatability Score (4-component, sensitivity-tested)
# ============================================================
cat("\n=== Computing Translatability Score ===\n")

# Components (all normalized 0-1):
# 1. Gene concordance: best_n_concordant / 5 (across all signatures)
# 2. Module preservation: mean_Zsummary / 20 (capped at 1)
# 3. Multi-signature agreement: n_signatures_concordant / 4
# 4. Category bonus: Conserved = 0.3, Moderate = 0.15 (from primary_category)

atlas[, score_gene := best_n_concordant / 5]

if ("mean_Zsummary" %in% names(atlas)) {
  atlas[, score_module := pmin(mean_Zsummary / 20, 1)]
  atlas[is.na(score_module), score_module := 0]
} else {
  atlas[, score_module := 0]
}

atlas[, score_multi_sig := n_signatures_concordant / 4]

atlas[, score_category := fcase(
  primary_category == "Conserved", 0.3,
  primary_category == "Moderate_Concordance", 0.15,
  default = 0
)]

# Default weights (documented and sensitivity-tested below)
W_GENE <- 0.40
W_MODULE <- 0.20
W_MULTI <- 0.25
W_CAT <- 0.15

atlas[, translatability_score := round(
  W_GENE * score_gene + W_MODULE * score_module +
  W_MULTI * score_multi_sig + W_CAT * score_category, 3)]

# -----------------------------------------------------------------------
# Cross-anchor confidence tier
# Combines NAFL-vs-NASH and disease_vs_ctrl anchors to classify genes by
# the robustness of their cross-species conservation.
#
# Dual_Conserved  — Conserved in BOTH anchors: the highest-confidence
#                   Cas13 targets. These genes are elevated in disease vs
#                   healthy AND increase during NAFL→NASH, and both signals
#                   are concordant across ≥3 mouse diet models.
#
# NASH_Conserved  — Conserved only in NAFL-vs-NASH: captures the
#                   NASH progression program well-modeled by mice, but the
#                   gene may not be broadly elevated vs healthy controls.
#
# MASLD_Conserved — Conserved only in disease_vs_ctrl: broadly elevated
#                   in MASLD vs healthy and conserved in mice, but does not
#                   specifically increase in the NAFL→NASH step (may be
#                   stably elevated throughout all disease stages).
#
# Other           — neither anchor classifies the gene as Conserved.
# -----------------------------------------------------------------------
has_dvc <- "dvc_category" %in% names(atlas)

if (has_dvc) {
  atlas[, cross_anchor_tier := fcase(
    primary_category == "Conserved" & dvc_category == "Conserved",  "Dual_Conserved",
    primary_category == "Conserved" & dvc_category != "Conserved",  "NASH_Conserved",
    primary_category != "Conserved" & dvc_category == "Conserved",  "MASLD_Conserved",
    default = "Other"
  )]

  # Dual_Conserved genes earn a small score bonus (capped at 1.0) to reflect
  # that both disease anchors agree on their cross-species conservation.
  atlas[cross_anchor_tier == "Dual_Conserved",
        translatability_score := pmin(round(translatability_score + 0.05, 3), 1.0)]

  cat("Cross-anchor confidence tiers:\n")
  print(atlas[, .N, by = cross_anchor_tier][order(-N)])
} else {
  cat("NOTE: dvc_category not found in gene_concordance_per_gene.csv —\n")
  cat("      re-run 01_corrected_gene_concordance.R to generate it.\n")
  atlas[, cross_anchor_tier := NA_character_]
}

# Re-classify tiers after Dual_Conserved bonus
atlas[, translatability_tier := fifelse(
  translatability_score >= 0.5, "High",
  fifelse(translatability_score >= 0.25, "Medium",
  fifelse(translatability_score > 0, "Low", "None"))
)]

cat("\nTranslatability tiers (after cross-anchor adjustment):\n")
print(atlas[, .N, by = translatability_tier][order(-N)])

# ============================================================
#  Sensitivity analysis: vary weights and report tier sizes
# ============================================================
cat("\n=== Translatability Score Sensitivity Analysis ===\n")

weight_sets <- data.table(
  name    = c("default", "gene_heavy", "module_heavy", "multi_heavy", "balanced", "no_category"),
  w_gene  = c(0.40,       0.60,         0.30,           0.30,          0.25,       0.45),
  w_module= c(0.20,       0.15,         0.40,           0.15,          0.25,       0.25),
  w_multi = c(0.25,       0.15,         0.15,           0.40,          0.25,       0.30),
  w_cat   = c(0.15,       0.10,         0.15,           0.15,          0.25,       0.00)
)

sensitivity <- weight_sets[, {
  s <- w_gene * atlas$score_gene + w_module * atlas$score_module +
       w_multi * atlas$score_multi_sig + w_cat * atlas$score_category
  .(n_high = sum(s >= 0.5), n_medium = sum(s >= 0.25 & s < 0.5),
    n_low = sum(s > 0 & s < 0.25), n_none = sum(s == 0))
}, by = .(name, w_gene, w_module, w_multi, w_cat)]

fwrite(sensitivity, file.path(RES, "translatability_sensitivity.csv"))
cat("Sensitivity analysis (6 weight sets):\n")
print(sensitivity[, .(name, n_high, n_medium, n_low)])

# ============================================================
#  Top translatable targets
# ============================================================
cat("\nTop 30 High-Confidence Translatable Targets:\n")
top_targets <- atlas[translatability_tier == "High"][order(-translatability_score)]
if (nrow(top_targets) > 0) {
  show_cols <- intersect(
    c("human_symbol", "cross_anchor_tier", "primary_category", "dvc_category",
      "n_concordant", "dvc_n_concordant", "n_signatures_concordant",
      "translatability_score", "diets_concordant"),
    names(atlas)
  )
  print(top_targets[1:min(30, nrow(top_targets)), ..show_cols])
}

# ============================================================
#  Save
# ============================================================
# Order columns — cross_anchor_tier and dvc_* columns placed after primary
front_cols <- c("human_symbol", "mouse_gene_id",
  "cross_anchor_tier",
  "primary_category", "dvc_category",
  "n_concordant", "n_discordant", "n_diets_sig", "mean_h_lfc",
  "diets_concordant", "diets_discordant",
  "dvc_n_concordant", "dvc_n_discordant", "dvc_diets_concordant",
  "dvc_h_significant", "dvc_mean_h_lfc",
  "best_n_concordant", "n_signatures_concordant", "best_category", "best_signature",
  "translatability_score", "translatability_tier")
other_cols <- setdiff(names(atlas), front_cols)
setcolorder(atlas, c(intersect(front_cols, names(atlas)), other_cols))

fwrite(atlas[order(-translatability_score)],
  file.path(RES, "concordance_atlas_unified.csv"))
cat("\nSaved: concordance_atlas_unified.csv\n")

# ============================================================
#  Update unified disease signatures
# ============================================================
cat("\n=== Updating Unified Disease Signatures ===\n")
uds_path <- file.path(DS_DIR, "unified_disease_signatures.csv")
if (file.exists(uds_path)) {
  uds <- fread(uds_path)

  # Strip old concordance columns
  drop_cols <- intersect(names(uds), c("concordance_category", "n_diets_concordant",
    "translatability_score", "translatability_tier"))
  if (length(drop_cols) > 0) uds[, (drop_cols) := NULL]

  # Map via symbol
  conc_info <- atlas[, .(human_symbol, concordance_category = primary_category,
    n_diets_concordant = n_concordant,
    cross_anchor_tier,
    dvc_category,
    translatability_score, translatability_tier)]
  conc_info <- conc_info[!duplicated(human_symbol)]

  uds <- merge(uds, conc_info, by.x = "symbol", by.y = "human_symbol", all.x = TRUE)
  fwrite(uds, uds_path)
  cat("Updated: unified_disease_signatures.csv\n")
}

# ============================================================
#  Summary statistics
# ============================================================
cat("\n=== ATLAS SUMMARY ===\n")
cat(sprintf("  Total genes in atlas: %d\n", nrow(atlas)))
cat(sprintf("  Genes with translatability score > 0: %d\n", sum(atlas$translatability_score > 0)))
cat(sprintf("  High-confidence targets (score >= 0.5): %d\n", sum(atlas$translatability_tier == "High")))
cat(sprintf("  Medium-confidence targets (0.25-0.5):   %d\n", sum(atlas$translatability_tier == "Medium")))

cat("\nCross-anchor confidence tiers:\n")
cat("  (Dual = conserved in BOTH NAFL-vs-NASH AND disease-vs-ctrl anchors)\n")
print(atlas[, .N, by = cross_anchor_tier][order(-N)])

cat("\nPrimary (NAFL-vs-NASH) category distribution:\n")
print(atlas[, .N, by = primary_category][order(-N)])

if ("dvc_category" %in% names(atlas)) {
  cat("\nDisease-vs-Control anchor category distribution:\n")
  print(atlas[, .N, by = dvc_category][order(-N)])
}

if (!is.null(gene_matrix)) {
  cat("\nBest model per human signature (by ρ):\n")
  best <- gene_matrix[, .SD[which.max(rho_all)], by = human_signature]
  print(best[, .(human_signature, diet, rho_all, concordance_pct)])
}

cat("\n=== Phase 5 complete ===\n")
cat("=== ALL PHASES DONE ===\n")
