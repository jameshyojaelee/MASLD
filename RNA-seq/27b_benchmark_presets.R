#!/usr/bin/env Rscript
# 27b_benchmark_presets.R
# Benchmark preset configurations against positive controls
#
# Reads the multi_evidence_atlas.csv and evaluates 5 preset filtering
# configurations by computing precision, recall, F1, and drug recovery.
# Grid-searches the "Therapeutic recovery" preset for optimal thresholds.
#
# Output:
#   - benchmark_results.csv (per-preset metrics)
#   - preset_configurations.json (threshold definitions)
#   - pr_curve_data.csv (for Figure 7d)

library(data.table)
library(jsonlite)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTDIR <- file.path(BASE, "RNA-seq/results/multi_evidence")

cat("=== Multi-Evidence Atlas Benchmarking ===\n\n")

# ================================================================
# Load atlas
# ================================================================
atlas_file <- file.path(OUTDIR, "multi_evidence_atlas.csv")
if (!file.exists(atlas_file)) {
  stop("Atlas file not found. Run 27a_assemble_evidence_atlas.R first.")
}
atlas <- fread(atlas_file)
cat("Atlas loaded:", nrow(atlas), "genes x", ncol(atlas), "columns\n")
stopifnot(all(c("bulk_padj","bulk_logFC") %in% names(atlas)))

# ================================================================
# Load positive controls
# ================================================================
pc_file <- file.path(BASE, "RNA-seq/results/validation/positive_control_validation.csv")
pc <- fread(pc_file)

# Expression-driven controls (DE expected)
expr_controls <- pc[control_type == "Expression_driven"]$gene
# Known drug targets (pharmacological mechanisms — NOT all expected to be DE)
# These are tracked separately as an informational metric
pharma_targets <- c("THRB", "NR1H4", "PPARG", "GLP1R", "ACACB", "ACACA")
# Convergent drug targets (high convergence score, most ARE DE) used as actual benchmark
drug_targets <- character(0)  # populated from convergent_drug_targets.csv below
# Genetic-risk controls (NOT expected to be DE — these should NOT be selected)
genetic_controls <- pc[control_type == "Genetic_risk"]$gene

# Combine expression controls for positive set
positive_in_atlas <- expr_controls[expr_controls %in% atlas$human_symbol]
cat("\nExpression controls in atlas:", length(positive_in_atlas), "/", length(expr_controls), "\n")
cat("Genetic-risk controls (expected absent):", paste(genetic_controls[genetic_controls %in% atlas$human_symbol], collapse=", "), "\n")

# Convergent drug targets (these ARE DE-expected — high convergence score)
drug_file <- file.path(BASE, "RNA-seq/results/drug_repurposing/convergent_drug_targets.csv")
if (file.exists(drug_file)) {
  conv_drugs <- fread(drug_file)
  # Top convergent targets with convergence_score >= 5
  drug_targets <- conv_drugs[convergence_score >= 5]$symbol
  drug_targets <- drug_targets[drug_targets %in% atlas$human_symbol]
  cat("Convergent drug targets (score>=5) in atlas:", length(drug_targets), "\n")
} else {
  drug_targets <- character(0)
}

# Pharmacological targets (NOT expected to be DE — track separately)
cat("Pharmacological mechanism targets:", paste(pharma_targets, collapse=", "), "\n")
pharma_in_atlas <- pharma_targets[pharma_targets %in% atlas$human_symbol]
cat("  In atlas:", length(pharma_in_atlas), "\n")
pharma_de <- atlas[human_symbol %in% pharma_in_atlas & !is.na(bulk_padj) & bulk_padj < 0.1]
cat("  With padj < 0.1:", nrow(pharma_de), "(", paste(pharma_de$human_symbol, collapse=", "), ")\n")
cat("  NOTE: Most pharma targets are NOT expected to pass DE filters\n")

# Full positive set: expression controls + convergent drug targets
full_positive <- unique(c(positive_in_atlas, drug_targets))
full_positive <- full_positive[full_positive %in% atlas$human_symbol]
cat("Full positive set in atlas:", length(full_positive), "\n")

# ================================================================
# Filtering function
# ================================================================
apply_preset <- function(atlas, config) {
  # Start with all genes
  selected <- rep(TRUE, nrow(atlas))

  # L1 filter (always applied as AND)
  if (!is.null(config$l1_padj)) {
    selected <- selected & !is.na(atlas$bulk_padj) & atlas$bulk_padj < config$l1_padj
  }
  if (!is.null(config$l1_lfc)) {
    selected <- selected & !is.na(atlas$bulk_logFC) & abs(atlas$bulk_logFC) > config$l1_lfc
  }

  # Additional layer filters — collect as booleans
  layer_passes <- list()

  # L2: Mouse — EXCLUDED from filtering (supplementary only)

  # L3
  if (!is.null(config$l3_categories)) {
    layer_passes$L3 <- !is.na(atlas$primary_category) & atlas$primary_category %in% config$l3_categories
  }
  if (!is.null(config$l3_min_concordant)) {
    l3_conc <- !is.na(atlas$n_concordant_diets) & atlas$n_concordant_diets >= config$l3_min_concordant
    if ("L3" %in% names(layer_passes)) {
      layer_passes$L3 <- layer_passes$L3 & l3_conc
    } else {
      layer_passes$L3 <- l3_conc
    }
  }

  # L4 (includes sc-eQTL, Broadaway, Zenodo, ieQTL + Phase 5 overhaul evidence)
  # (MR columns — mr_pval, mr_sig, mr_ivw_pval — removed 2026-04-22; MR ditched
  #  from paper. Preset `l4_mr_or_twas_sig` retained as an alias that now only
  #  considers TWAS + COLOC + INTACT signals.)
  if (isTRUE(config$l4_any_causal) || isTRUE(config$l4_mr_or_twas_sig)) {
    layer_passes$L4 <- (!is.na(atlas$twas_pval) & atlas$twas_pval < 0.05) |
                        (!is.na(atlas$coloc_pp4) & atlas$coloc_pp4 > 0.5) |
                        (if ("sceqtl_coloc_best_pp4" %in% names(atlas))
                          !is.na(atlas$sceqtl_coloc_best_pp4) & atlas$sceqtl_coloc_best_pp4 > 0.5
                         else rep(FALSE, nrow(atlas))) |
                        (if ("broadaway_coloc_pp4" %in% names(atlas))
                          !is.na(atlas$broadaway_coloc_pp4) & atlas$broadaway_coloc_pp4 > 0.5
                         else rep(FALSE, nrow(atlas))) |
                        (if ("zenodo_nafld_coloc" %in% names(atlas))
                          !is.na(atlas$zenodo_nafld_coloc) & atlas$zenodo_nafld_coloc == TRUE
                         else rep(FALSE, nrow(atlas))) |
                        (if ("ieqtl_disease_interaction" %in% names(atlas))
                          !is.na(atlas$ieqtl_disease_interaction) & atlas$ieqtl_disease_interaction == TRUE
                         else rep(FALSE, nrow(atlas))) |
                        # Phase 5 overhaul columns
                        (if ("otters_broadaway_pval" %in% names(atlas))
                          !is.na(atlas$otters_broadaway_pval) & atlas$otters_broadaway_pval < 0.05
                         else rep(FALSE, nrow(atlas))) |
                        (if ("coloc_susie_best_pp4" %in% names(atlas))
                          !is.na(atlas$coloc_susie_best_pp4) & atlas$coloc_susie_best_pp4 > 0.5
                         else rep(FALSE, nrow(atlas))) |
                        (if ("hyprcoloc_posterior" %in% names(atlas))
                          !is.na(atlas$hyprcoloc_posterior) & atlas$hyprcoloc_posterior > 0.5
                         else rep(FALSE, nrow(atlas)))
                        # cTWAS dropped (review A09#1) — branch removed
  }

  # L4 causal_robust (from Phase 5 overhaul)
  if (isTRUE(config$l4_causal_robust)) {
    if ("causal_robustness" %in% names(atlas)) {
      layer_passes$L4 <- !is.na(atlas$causal_robustness) & atlas$causal_robustness == "Robust"
    } else if ("causal_methods_sig" %in% names(atlas)) {
      layer_passes$L4 <- !is.na(atlas$causal_methods_sig) & atlas$causal_methods_sig >= 3
    }
  }

  # L5
  if (!is.null(config$l5_exclude_divergent) && config$l5_exclude_divergent) {
    layer_passes$L5 <- is.na(atlas$sex_class) | atlas$sex_class != "Sex_divergent"
  }
  if (!is.null(config$l5_classes)) {
    layer_passes$L5 <- !is.na(atlas$sex_class) & atlas$sex_class %in% config$l5_classes
  }

  # L6
  if (!is.null(config$l6_min_pathways)) {
    layer_passes$L6 <- !is.na(atlas$n_leading_edge_pathways) &
                        atlas$n_leading_edge_pathways >= config$l6_min_pathways
  }

  # L7
  if (isTRUE(config$l7_non_essential)) {
    layer_passes$L7 <- is.na(atlas$is_essential) | !atlas$is_essential
  }

  # Combine layer passes according to logic
  if (length(layer_passes) > 0) {
    if (!is.null(config$logic) && config$logic == "AND") {
      # All specified layers must pass
      for (lp in layer_passes) {
        selected <- selected & lp
      }
    } else if (!is.null(config$logic) && config$logic == "OR_with_min") {
      # OR across secondary layers, but must pass min_layers_required
      layer_matrix <- do.call(cbind, layer_passes)
      layer_sum <- rowSums(layer_matrix, na.rm = TRUE)
      min_req <- ifelse(is.null(config$min_layers_required), 1, config$min_layers_required)
      selected <- selected & (layer_sum >= min_req)
    } else {
      # Default: AND with L1, OR across secondary layers, AND L7 safety
      secondary <- layer_passes[!names(layer_passes) %in% "L7"]
      if (length(secondary) > 0) {
        any_secondary <- Reduce(`|`, secondary)
        selected <- selected & any_secondary
      }
      if ("L7" %in% names(layer_passes)) {
        selected <- selected & layer_passes$L7
      }
    }
  }

  # Min layers_active
  if (!is.null(config$min_layers_active)) {
    selected <- selected & atlas$layers_active >= config$min_layers_active
  }

  return(which(selected))
}

# ================================================================
# Evaluate a selection against positive controls
# ================================================================
evaluate_selection <- function(atlas, selected_idx, positive_set, known_drugs = NULL) {
  if (is.null(known_drugs)) known_drugs <- drug_targets
  selected_genes <- atlas$human_symbol[selected_idx]
  n_selected <- length(selected_genes)

  tp <- sum(selected_genes %in% positive_set)
  precision <- ifelse(n_selected > 0, tp / n_selected, 0)
  recall <- ifelse(length(positive_set) > 0, tp / length(positive_set), 0)
  f1 <- ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)

  # Known drug recovery
  drugs_recovered <- sum(known_drugs %in% selected_genes)
  drug_recovery <- if (length(known_drugs) == 0) 0 else drugs_recovered / length(known_drugs)

  # Genetic controls captured (should be 0 ideally)
  genetic_captured <- sum(genetic_controls %in% selected_genes)

  list(
    n_selected = n_selected,
    true_positives = tp,
    precision = precision,
    recall = recall,
    f1 = f1,
    drugs_recovered = drugs_recovered,
    drug_recovery_rate = drug_recovery,
    genetic_controls_captured = genetic_captured,
    selected_genes = selected_genes
  )
}

# ================================================================
# Define 5 preset configurations
# ================================================================
presets <- list(
  therapeutic_recovery = list(
    name = "Therapeutic recovery",
    description = "Optimized to recover known MASLD therapeutics; AND(L1), OR(L4-L6), AND(L7)",
    l1_padj = 0.1,  # Will be grid-searched
    l1_lfc = 0.5,   # Will be grid-searched
    l4_any_causal = TRUE,
    l6_min_pathways = 1,
    l7_non_essential = TRUE,
    logic = "default",
    min_layers_active = 2  # Will be grid-searched
  ),
  stringent_multi_evidence = list(
    name = "Stringent multi-evidence",
    description = "High-confidence targets with convergent evidence across multiple layers",
    l1_padj = 0.05,
    l1_lfc = 1.0,
    l3_categories = "Conserved",
    l5_exclude_divergent = TRUE,
    l6_min_pathways = 3,
    l7_non_essential = TRUE,
    logic = "OR_with_min",
    min_layers_required = 2
  ),
  broad_discovery = list(
    name = "Broad discovery",
    description = "Inclusive selection for screening libraries; L1 significance + any supporting layer",
    l1_padj = 0.1,
    l1_lfc = 0.5,
    min_layers_active = 2,
    logic = "AND"
  ),
  cross_species_conserved = list(
    name = "Cross-species conserved",
    description = "Genes significant in human with conserved cross-species directionality",
    l1_padj = 0.1,
    l3_categories = "Conserved",
    logic = "AND"
  ),
  causal_supported = list(
    name = "Causal-supported",
    description = "Genes with both transcriptomic DE and causal genetic evidence",
    l1_padj = 0.1,
    l4_mr_or_twas_sig = TRUE,
    logic = "AND"
  ),
  causal_robust = list(
    name = "Causal-robust",
    description = "Genes with 3+ concordant causal methods (TWAS/COLOC/HyPrColoc/cTWAS)",
    l1_padj = 0.1,
    l4_causal_robust = TRUE,
    logic = "AND"
  )
)

# ================================================================
# Grid search for therapeutic_recovery preset
# ================================================================
cat("\n=== Grid search: Therapeutic recovery preset ===\n")

padj_grid <- c(0.01, 0.05, 0.1)
lfc_grid <- c(0.3, 0.5, 0.8, 1.0)
min_layers_grid <- c(2, 3, 4)

grid_results <- data.table()

for (padj_thresh in padj_grid) {
  for (lfc_thresh in lfc_grid) {
    for (min_layers in min_layers_grid) {
      config <- presets$therapeutic_recovery
      config$l1_padj <- padj_thresh
      config$l1_lfc <- lfc_thresh
      config$min_layers_active <- min_layers

      selected_idx <- apply_preset(atlas, config)
      metrics <- evaluate_selection(atlas, selected_idx, full_positive)

      grid_results <- rbindlist(list(grid_results, data.table(
        l1_padj = padj_thresh,
        l1_lfc = lfc_thresh,
        min_layers = min_layers,
        n_selected = metrics$n_selected,
        precision = metrics$precision,
        recall = metrics$recall,
        f1 = metrics$f1,
        drug_recovery = metrics$drug_recovery_rate,
        drugs_recovered = metrics$drugs_recovered
      )))
    }
  }
}

# Find F1-optimal
best_idx <- which.max(grid_results$f1)
best_config <- grid_results[best_idx]
cat("\nBest therapeutic recovery configuration:\n")
print(best_config)

# Update therapeutic_recovery preset with optimal values
presets$therapeutic_recovery$l1_padj <- best_config$l1_padj
presets$therapeutic_recovery$l1_lfc <- best_config$l1_lfc
presets$therapeutic_recovery$min_layers_active <- best_config$min_layers

# ================================================================
# Evaluate all presets
# ================================================================
cat("\n=== Evaluating all presets ===\n")

benchmark_results <- data.table()

for (preset_name in names(presets)) {
  config <- presets[[preset_name]]
  selected_idx <- apply_preset(atlas, config)
  metrics <- evaluate_selection(atlas, selected_idx, full_positive)

  cat(sprintf("\n%s: %d genes selected, precision=%.3f, recall=%.3f, F1=%.3f, drugs=%d/%d\n",
              config$name, metrics$n_selected, metrics$precision, metrics$recall,
              metrics$f1, metrics$drugs_recovered, length(drug_targets)))

  # Show which drugs recovered
  drugs_found <- drug_targets[drug_targets %in% metrics$selected_genes]
  if (length(drugs_found) > 0) {
    cat("  Drugs recovered:", paste(drugs_found, collapse = ", "), "\n")
  }

  benchmark_results <- rbindlist(list(benchmark_results, data.table(
    preset = preset_name,
    preset_label = config$name,
    n_selected = metrics$n_selected,
    true_positives = metrics$true_positives,
    precision = round(metrics$precision, 4),
    recall = round(metrics$recall, 4),
    f1 = round(metrics$f1, 4),
    drugs_recovered = metrics$drugs_recovered,
    drug_recovery_rate = round(metrics$drug_recovery_rate, 4),
    genetic_controls_captured = metrics$genetic_controls_captured
  )))
}

# ================================================================
# PR curve data (varying stringency for Figure 7d)
# ================================================================
cat("\n=== Generating PR curve data ===\n")

pr_data <- data.table()

# Vary L1 padj threshold
for (padj_thresh in c(0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3)) {
  for (min_layers in c(1, 2, 3, 4, 5)) {
    config <- list(l1_padj = padj_thresh, min_layers_active = min_layers, logic = "AND")
    selected_idx <- apply_preset(atlas, config)
    metrics <- evaluate_selection(atlas, selected_idx, full_positive)
    pr_data <- rbindlist(list(pr_data, data.table(
      padj_threshold = padj_thresh,
      min_layers = min_layers,
      n_selected = metrics$n_selected,
      precision = metrics$precision,
      recall = metrics$recall,
      f1 = metrics$f1
    )))
  }
}

# ================================================================
# Output
# ================================================================

# Benchmark results
fwrite(benchmark_results, file.path(OUTDIR, "benchmark_results.csv"))
cat("\nSaved:", file.path(OUTDIR, "benchmark_results.csv"), "\n")

# Grid search results
fwrite(grid_results, file.path(OUTDIR, "grid_search_results.csv"))
cat("Saved:", file.path(OUTDIR, "grid_search_results.csv"), "\n")

# PR curve data
fwrite(pr_data, file.path(OUTDIR, "pr_curve_data.csv"))
cat("Saved:", file.path(OUTDIR, "pr_curve_data.csv"), "\n")

# Preset configurations as JSON
preset_json <- lapply(presets, function(p) {
  # Remove NULL values for clean JSON
  p[!sapply(p, is.null)]
})
write_json(preset_json, file.path(OUTDIR, "preset_configurations.json"),
           pretty = TRUE, auto_unbox = TRUE)
cat("Saved:", file.path(OUTDIR, "preset_configurations.json"), "\n")

# ================================================================
# Summary
# ================================================================
cat("\n=== Benchmark Summary ===\n")
print(benchmark_results[, .(preset_label, n_selected, precision, recall, f1, drug_recovery_rate)])

# Suggested default (therapeutic_recovery)
suggested <- benchmark_results[preset == "therapeutic_recovery"]
cat(sprintf("\nSuggested default (Therapeutic recovery):\n"))
cat(sprintf("  L1 padj < %s, |LFC| > %s, min_layers >= %d\n",
            presets$therapeutic_recovery$l1_padj,
            presets$therapeutic_recovery$l1_lfc,
            presets$therapeutic_recovery$min_layers_active))
cat(sprintf("  Selects %d genes, recovers %d/%d known drug targets\n",
            suggested$n_selected, suggested$drugs_recovered, length(drug_targets)))
cat(sprintf("  Precision=%.3f, Recall=%.3f, F1=%.3f\n",
            suggested$precision, suggested$recall, suggested$f1))

cat("\nDone.\n")
