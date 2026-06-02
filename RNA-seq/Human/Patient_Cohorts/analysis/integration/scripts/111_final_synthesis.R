#!/usr/bin/env Rscript
# =========================================================================
# 111_final_synthesis.R
# Comprehensive comparison of ALL foundation-model staging classifiers +
# 12-panel publication figure.
#
# Loads every result file produced by scripts 97-110, builds a unified
# comparison table, and generates a Nature-compatible supplementary figure.
#
# Input  (all from results/staging_classifier/):
#   v3_model_summary.csv, nas_vae_model_summary.csv,
#   concept_bottleneck_summary.csv, multitask_nas_summary.csv,
#   fibrosis_transfer_results.csv, stacking_ensemble_results.csv,
#   ablation_summary.csv, sex_auroc_comparison.csv,
#   plasma_classifier_results.csv, sensitivity_bootstrap_ci.csv,
#   lto_summary.csv, v3_vs_v2_comparison.csv,
#   concept_activations.csv, nas_embeddings_all_samples.csv,
#   modeling_metadata.csv, plasma_panel_curve.csv
#
# Output:
#   figures/supplementary/figS_staging_classifier/fig_model_comparison.pdf  (12 panels)
#   results/staging_classifier/foundation_model_comparison.csv
#   results/staging_classifier/foundation_model_summary.csv
#
# SLURM: cpu partition, 8 CPUs, 32G RAM, 48h
# Env:   micromamba activate rnaseq
#
# Usage:
#   sbatch --job-name=stg111_synthesis \
#          --partition=cpu --cpus-per-task=8 --mem=32G --time=48:00:00 \
#          --output=logs/111_synthesis_%j.out \
#          --error=logs/111_synthesis_%j.err \
#          --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && \
#                  micromamba activate rnaseq && \
#                  cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \
#                  Rscript 111_final_synthesis.R'"
# =========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

cat("=======================================================================\n")
cat("111: Final Synthesis — All-Model Comparison + Publication Figure\n")
cat("=======================================================================\n")

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/staging_classifier")
# Canonical figures root (project rule: PDF only, outputs under figures/).
FIGDIR <- file.path(BASE, "figures/supplementary/figS_staging_classifier")
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Publication theme
# -------------------------------------------------------------------------
theme_pub <- function(base_size = 8) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1,
                                       hjust = 0),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size - 1),
      legend.text      = element_text(size = base_size - 1),
      legend.key.size  = unit(0.3, "cm"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(4, 4, 4, 4)
    )
}

# Helper: create a display-friendly model label that stays unique.
# sub("\\(.*\\)", ...) can collapse different models to the same string
# (e.g. "Transfer (A)" and "Transfer (B)" both become "Transfer ").
# This function keeps the full parenthetical or abbreviates it, ensuring uniqueness.
make_model_short <- function(model_vec) {
  # Abbreviate parenthetical: keep first word only
  short <- sub("\\(([A-Za-z0-9_]+)[^)]*\\)", "(\\1)", model_vec)
  # If still duplicated, keep the original
  dups <- duplicated(short) | duplicated(short, fromLast = TRUE)
  short[dups] <- model_vec[dups]
  # Trim trailing whitespace
  trimws(short)
}

# Colour palette for models
model_cols <- c(
  "V3 Elastic Net"       = "#2166ac",
  "NAS-VAE + RF"         = "#4393c3",
  "Concept Bottleneck"   = "#92c5de",
  "Multi-task CORN"      = "#d6604d",
  "Transfer Learning"    = "#b2182b",
  "Stacking Ensemble"    = "#1b7837",
  "Plasma (Olink)"       = "#762a83",
  "NFASC+GDF15"          = "#c2a5cf"
)

# -------------------------------------------------------------------------
# Helper: safe fread
# -------------------------------------------------------------------------
safe_fread <- function(path, ...) {
  if (file.exists(path)) {
    tryCatch(fread(path, ...), error = function(e) {
      message("  Warning: could not read ", basename(path), ": ", e$message)
      NULL
    })
  } else {
    message("  Not found: ", basename(path))
    NULL
  }
}

# =========================================================================
# 1. Load all result files
# =========================================================================
cat("\n--- Loading result files ---\n")

v3_sum     <- safe_fread(file.path(RDIR, "v3_model_summary.csv"))
vae_sum    <- safe_fread(file.path(RDIR, "nas_vae_model_summary.csv"))
cbm_sum    <- safe_fread(file.path(RDIR, "concept_bottleneck_summary.csv"))
mtl_sum    <- safe_fread(file.path(RDIR, "multitask_nas_summary.csv"))
transfer   <- safe_fread(file.path(RDIR, "fibrosis_transfer_results.csv"))
stacking   <- safe_fread(file.path(RDIR, "clean_stacking_results.csv"))
ablation   <- safe_fread(file.path(RDIR, "ablation_summary.csv"))
sex_comp   <- safe_fread(file.path(RDIR, "sex_auroc_comparison.csv"))
plasma     <- safe_fread(file.path(RDIR, "plasma_classifier_results.csv"))
bootstrap  <- safe_fread(file.path(RDIR, "sensitivity_bootstrap_ci.csv"))
lto        <- safe_fread(file.path(RDIR, "lto_summary.csv"))
v3v2       <- safe_fread(file.path(RDIR, "v3_vs_v2_comparison.csv"))
concepts   <- safe_fread(file.path(RDIR, "concept_activations.csv"))
embed      <- safe_fread(file.path(RDIR, "nas_embeddings_all_samples.csv"))
meta       <- safe_fread(file.path(RDIR, "modeling_metadata.csv"))
panel_curve<- safe_fread(file.path(RDIR, "plasma_panel_curve.csv"))
stacking_cmp <- safe_fread(file.path(RDIR, "stacking_vs_single_comparison.csv"))

# =========================================================================
# 2. Build unified comparison table
# =========================================================================
cat("\n--- Building unified comparison table ---\n")

comparison_rows <- list()

# --- V3 Elastic Net ---
if (!is.null(v3_sum)) {
  for (i in seq_len(nrow(v3_sum))) {
    r <- v3_sum[i]
    comparison_rows[[length(comparison_rows) + 1]] <- data.table(
      model = "V3 Elastic Net",
      target = r$target,
      metric_auroc = r$mean_auroc,
      metric_auroc_sd = r$sd_auroc,
      metric_qwk = r$mean_qwk,
      metric_qwk_sd = r$sd_qwk,
      n_folds = r$n_folds,
      n_samples = r$n_samples
    )
  }
}

# --- NAS-VAE + RF (best model per target) ---
if (!is.null(vae_sum)) {
  best_vae <- vae_sum[, .SD[which.max(mean_auroc)], by = .(target)]
  for (i in seq_len(nrow(best_vae))) {
    r <- best_vae[i]
    comparison_rows[[length(comparison_rows) + 1]] <- data.table(
      model = paste0("NAS-VAE + ", r$model),
      target = r$target,
      metric_auroc = r$mean_auroc,
      metric_auroc_sd = r$std_auroc,
      metric_qwk = r$mean_qwk,
      metric_qwk_sd = NA_real_,
      n_folds = r$n_folds,
      n_samples = NA_integer_
    )
  }
}

# --- Concept Bottleneck ---
if (!is.null(cbm_sum)) {
  for (i in seq_len(nrow(cbm_sum))) {
    r <- cbm_sum[i]
    comparison_rows[[length(comparison_rows) + 1]] <- data.table(
      model = "Concept Bottleneck",
      target = r$task,
      metric_auroc = r$auroc_mean,
      metric_auroc_sd = r$auroc_std,
      metric_qwk = r$qwk_mean,
      metric_qwk_sd = r$qwk_std,
      n_folds = r$n_folds,
      n_samples = NA_integer_
    )
  }
}

# --- Multi-task CORN ---
if (!is.null(mtl_sum)) {
  for (i in seq_len(nrow(mtl_sum))) {
    r <- mtl_sum[i]
    comparison_rows[[length(comparison_rows) + 1]] <- data.table(
      model = "Multi-task CORN",
      target = r$task,
      metric_auroc = r$auroc_mean,
      metric_auroc_sd = r$auroc_std,
      metric_qwk = r$qwk_mean,
      metric_qwk_sd = r$qwk_std,
      n_folds = r$n_folds,
      n_samples = NA_integer_
    )
  }
}

# --- Transfer Learning ---
if (!is.null(transfer)) {
  # Aggregate by approach + target (take best approach per target)
  tr_agg <- transfer[, .(
    mean_auroc = mean(auroc, na.rm = TRUE),
    mean_qwk   = mean(qwk, na.rm = TRUE),
    n_folds    = .N
  ), by = .(approach, target)]

  best_tr <- tr_agg[, .SD[which.max(fifelse(is.na(mean_auroc), mean_qwk, mean_auroc))],
                     by = .(target)]
  for (i in seq_len(nrow(best_tr))) {
    r <- best_tr[i]
    comparison_rows[[length(comparison_rows) + 1]] <- data.table(
      model = paste0("Transfer (", r$approach, ")"),
      target = r$target,
      metric_auroc = r$mean_auroc,
      metric_auroc_sd = NA_real_,
      metric_qwk = r$mean_qwk,
      metric_qwk_sd = NA_real_,
      n_folds = r$n_folds,
      n_samples = NA_integer_
    )
  }
}

# --- Stacking Ensemble (clean, leakage-free summary) ---
if (!is.null(stacking) && nrow(stacking) > 0) {
  # clean_stacking_results.csv has columns:
  #   nas_4group_qwk, nas_ge5_auroc, n_samples_ordinal, n_samples_binary, n_folds, method, leakage_free
  if ("nas_ge5_auroc" %in% names(stacking)) {
    comparison_rows[[length(comparison_rows) + 1]] <- data.table(
      model = "Stacking Ensemble",
      target = "nas_ge5",
      metric_auroc = stacking$nas_ge5_auroc[1],
      metric_auroc_sd = NA_real_,
      metric_qwk = NA_real_,
      metric_qwk_sd = NA_real_,
      n_folds = stacking$n_folds[1],
      n_samples = stacking$n_samples_binary[1]
    )
  }
  if ("nas_4group_qwk" %in% names(stacking)) {
    comparison_rows[[length(comparison_rows) + 1]] <- data.table(
      model = "Stacking Ensemble",
      target = "nas_4group",
      metric_auroc = NA_real_,
      metric_auroc_sd = NA_real_,
      metric_qwk = stacking$nas_4group_qwk[1],
      metric_qwk_sd = NA_real_,
      n_folds = stacking$n_folds[1],
      n_samples = stacking$n_samples_ordinal[1]
    )
  }
}

# --- Plasma ---
if (!is.null(plasma)) {
  pl_5fold <- plasma[cv_method == "5fold_repeated"]
  for (i in seq_len(nrow(pl_5fold))) {
    r <- pl_5fold[i]
    comparison_rows[[length(comparison_rows) + 1]] <- data.table(
      model = paste0("Plasma (", r$experiment, ")"),
      target = "fib_ge3_plasma",
      metric_auroc = r$auroc,
      metric_auroc_sd = r$auroc_sd,
      metric_qwk = NA_real_,
      metric_qwk_sd = NA_real_,
      n_folds = NA_integer_,
      n_samples = NA_integer_
    )
  }
}

# Combine
if (length(comparison_rows) > 0) {
  comp_dt <- rbindlist(comparison_rows, fill = TRUE)
  cat(sprintf("Unified comparison table: %d rows x %d cols\n",
              nrow(comp_dt), ncol(comp_dt)))
} else {
  comp_dt <- data.table()
  cat("WARNING: No comparison rows generated\n")
}

# =========================================================================
# 3. Generate 12-panel publication figure
# =========================================================================
cat("\n--- Generating 12-panel figure ---\n")

panels <- list()

# ---- (a) F>=3 AUROC comparison (horizontal bars) ----
cat("  Panel (a): F>=3 AUROC comparison\n")
fib_comp <- comp_dt[grepl("fib_ge3|F>=3|disease", target, ignore.case = TRUE) &
                      !is.na(metric_auroc)]
if (nrow(fib_comp) > 0) {
  fib_comp[, model_short := make_model_short(model)]
  fib_comp <- fib_comp[order(metric_auroc)]
  fib_comp[, model_short := factor(model_short, levels = unique(model_short))]

  panels[["a"]] <- ggplot(fib_comp, aes(x = metric_auroc, y = model_short)) +
    geom_col(aes(fill = model_short), width = 0.7, show.legend = FALSE) +
    geom_errorbarh(aes(xmin = metric_auroc - fifelse(is.na(metric_auroc_sd), 0, metric_auroc_sd),
                       xmax = metric_auroc + fifelse(is.na(metric_auroc_sd), 0, metric_auroc_sd)),
                   height = 0.3, linewidth = 0.3) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(limits = c(0.5, 1.0), breaks = seq(0.5, 1.0, 0.1)) +
    scale_fill_manual(values = rep("#2166ac", 20)) +
    labs(x = "AUROC", y = NULL, title = "F>=3 Binary Classification") +
    theme_pub()
} else {
  panels[["a"]] <- ggplot() + theme_void() + ggtitle("(a) No F>=3 data")
}

# ---- (b) NAS QWK comparison ----
cat("  Panel (b): NAS QWK comparison\n")
nas_comp <- comp_dt[grepl("nas_4|nas_group4|NAS 4", target, ignore.case = TRUE) &
                      !is.na(metric_qwk)]
if (nrow(nas_comp) == 0) {
  # Try broader match
  nas_comp <- comp_dt[grepl("nas", target, ignore.case = TRUE) & !is.na(metric_qwk)]
}
if (nrow(nas_comp) > 0) {
  nas_comp <- nas_comp[order(metric_qwk)]
  nas_comp[, model_short := make_model_short(model)]
  nas_comp[, model_short := factor(model_short, levels = unique(model_short))]

  panels[["b"]] <- ggplot(nas_comp, aes(x = metric_qwk, y = model_short)) +
    geom_col(aes(fill = model_short), width = 0.7, show.legend = FALSE) +
    geom_errorbarh(aes(xmin = metric_qwk - fifelse(is.na(metric_qwk_sd), 0, metric_qwk_sd),
                       xmax = metric_qwk + fifelse(is.na(metric_qwk_sd), 0, metric_qwk_sd)),
                   height = 0.3, linewidth = 0.3) +
    scale_fill_manual(values = rep("#d6604d", 20)) +
    labs(x = "QWK", y = NULL, title = "NAS Ordinal (QWK)") +
    theme_pub()
} else {
  panels[["b"]] <- ggplot() + theme_void() + ggtitle("(b) No NAS QWK data")
}

# ---- (c) NAS 9-class QWK (VAE breakthrough) ----
cat("  Panel (c): NAS 9-class QWK\n")
nas9_comp <- comp_dt[grepl("nas_9|9.class", target, ignore.case = TRUE) &
                       !is.na(metric_qwk)]
if (nrow(nas9_comp) > 0) {
  nas9_comp <- nas9_comp[order(metric_qwk)]
  nas9_comp[, model_short := make_model_short(model)]
  nas9_comp[, model_short := factor(model_short, levels = unique(model_short))]
  nas9_comp[, highlight := grepl("VAE", model)]

  panels[["c"]] <- ggplot(nas9_comp, aes(x = metric_qwk, y = model_short)) +
    geom_col(aes(fill = highlight), width = 0.7, show.legend = FALSE) +
    scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = "#1b7837")) +
    labs(x = "QWK", y = NULL, title = "NAS 9-class Ordinal") +
    theme_pub()
} else {
  panels[["c"]] <- ggplot() + theme_void() + ggtitle("(c) No NAS 9-class data")
}

# ---- (d) Concept bottleneck activations heatmap ----
cat("  Panel (d): Concept activations heatmap\n")
if (!is.null(concepts) && !is.null(meta) && nrow(concepts) > 0) {
  # Merge with metadata for grouping
  concept_cols <- setdiff(names(concepts), "sample_id")
  # Limit to first 20 concepts if many
  if (length(concept_cols) > 20) concept_cols <- concept_cols[1:20]

  concept_m <- merge(concepts, meta[, .(sample_id, nas_group4, fib_stage)],
                     by = "sample_id")
  concept_m <- concept_m[nas_group4 >= 0]

  if (nrow(concept_m) > 0) {
    # Mean per NAS group
    concept_means <- concept_m[, lapply(.SD, mean, na.rm = TRUE),
                                by = .(nas_group4), .SDcols = concept_cols]
    # Melt for heatmap
    cm_long <- melt(concept_means, id.vars = "nas_group4",
                    variable.name = "concept", value.name = "activation")
    # Truncate concept names
    cm_long[, concept_short := substr(concept, 1, 20)]

    panels[["d"]] <- ggplot(cm_long, aes(x = factor(nas_group4), y = concept_short,
                                          fill = activation)) +
      geom_tile(colour = "white", linewidth = 0.2) +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                           midpoint = 0.5, name = "Act.") +
      labs(x = "NAS Group", y = NULL, title = "Concept Activations") +
      theme_pub() +
      theme(axis.text.y = element_text(size = 5))
  } else {
    panels[["d"]] <- ggplot() + theme_void() + ggtitle("(d) No concept data")
  }
} else {
  panels[["d"]] <- ggplot() + theme_void() + ggtitle("(d) Concepts unavailable")
}

# ---- (e) Sex-stratified AUROC ----
cat("  Panel (e): Sex-stratified AUROC\n")
if (!is.null(sex_comp) && nrow(sex_comp) > 0) {
  sex_long <- melt(sex_comp,
                   id.vars = c("strategy", "sex_subset"),
                   measure.vars = c("mean_auroc_female", "mean_auroc_male"),
                   variable.name = "sex", value.name = "auroc")
  sex_long <- sex_long[!is.na(auroc)]
  sex_long[, sex := fifelse(grepl("female", sex), "Female", "Male")]

  panels[["e"]] <- ggplot(sex_long, aes(x = strategy, y = auroc, fill = sex)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("Female" = "#c51b7d", "Male" = "#4d9221"),
                      name = "Sex") +
    labs(x = NULL, y = "AUROC", title = "Sex-Stratified Performance") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5))
} else {
  panels[["e"]] <- ggplot() + theme_void() + ggtitle("(e) No sex data")
}

# ---- (f) Modality ablation lollipop ----
cat("  Panel (f): Modality ablation\n")
if (!is.null(ablation) && nrow(ablation) > 0) {
  ablation[, combination_short := substr(combination, 1, 25)]
  ablation <- ablation[order(mean_auroc)]
  ablation[, combination_short := factor(combination_short,
                                          levels = combination_short)]

  panels[["f"]] <- ggplot(ablation, aes(x = mean_auroc, y = combination_short)) +
    geom_segment(aes(x = 0.5, xend = mean_auroc,
                     y = combination_short, yend = combination_short),
                 colour = "grey60", linewidth = 0.4) +
    geom_point(size = 2, colour = "#2166ac") +
    geom_errorbarh(aes(xmin = mean_auroc - fifelse(is.na(sd_auroc), 0, sd_auroc),
                       xmax = mean_auroc + fifelse(is.na(sd_auroc), 0, sd_auroc)),
                   height = 0.25, linewidth = 0.3) +
    scale_x_continuous(limits = c(0.5, 1.0)) +
    labs(x = "AUROC", y = NULL, title = "Feature Ablation") +
    theme_pub()
} else {
  panels[["f"]] <- ggplot() + theme_void() + ggtitle("(f) No ablation data")
}

# ---- (g) VAE UMAP coloured by NAS score ----
cat("  Panel (g): VAE UMAP by NAS score\n")
if (!is.null(embed) && !is.null(meta) && nrow(embed) > 0) {
  z_cols <- grep("^z", names(embed), value = TRUE)
  embed_m <- merge(embed, meta[, .(sample_id, nas_score)], by = "sample_id")
  embed_m <- embed_m[nas_score >= 0]

  if (nrow(embed_m) > 20) {
    # Quick PCA -> 2D (UMAP would need uwot)
    mat <- as.matrix(embed_m[, ..z_cols])
    pca <- prcomp(mat, center = TRUE, scale. = FALSE, rank. = 2)
    embed_m[, pc1 := pca$x[, 1]]
    embed_m[, pc2 := pca$x[, 2]]

    panels[["g"]] <- ggplot(embed_m, aes(x = pc1, y = pc2, colour = nas_score)) +
      geom_point(size = 0.4, alpha = 0.7) +
      scale_colour_viridis_c(option = "magma", name = "NAS") +
      labs(x = "PC 1", y = "PC 2", title = "VAE Embedding (NAS score)") +
      theme_pub() +
      theme(legend.position = c(0.85, 0.2))
  } else {
    panels[["g"]] <- ggplot() + theme_void() + ggtitle("(g) Too few NAS samples")
  }
} else {
  panels[["g"]] <- ggplot() + theme_void() + ggtitle("(g) No embeddings")
}

# ---- (h) Tissue-to-plasma funnel ----
cat("  Panel (h): Tissue-to-plasma funnel\n")
if (!is.null(plasma) && nrow(plasma) > 0) {
  plasma_5f <- plasma[cv_method == "5fold_repeated"]
  plasma_5f[, label := paste0(experiment, " (", n_proteins, ")")]
  plasma_5f <- plasma_5f[order(-auroc)]

  panels[["h"]] <- ggplot(plasma_5f, aes(x = reorder(label, auroc), y = auroc)) +
    geom_col(fill = "#762a83", width = 0.6) +
    geom_errorbar(aes(ymin = auroc - fifelse(is.na(auroc_sd), 0, auroc_sd),
                      ymax = auroc + fifelse(is.na(auroc_sd), 0, auroc_sd)),
                  width = 0.25, linewidth = 0.3) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    coord_flip() +
    scale_y_continuous(limits = c(0.5, 1.0)) +
    labs(x = NULL, y = "AUROC", title = "Tissue -> Plasma Funnel") +
    theme_pub() +
    theme(axis.text.y = element_text(size = 5))
} else {
  panels[["h"]] <- ggplot() + theme_void() + ggtitle("(h) No plasma data")
}

# ---- (i) Leave-two-out stress test ----
cat("  Panel (i): Leave-two-out stress test\n")
if (!is.null(lto) && nrow(lto) > 0) {
  lto_df <- data.table(
    metric = c("Ref AUROC", "Mean LTO", "Min LTO", "Max LTO"),
    value  = c(lto$ref_auroc_v3_loco, lto$mean_lto_auroc,
               lto$min_lto_auroc, lto$max_lto_auroc)
  )
  lto_df[, metric := factor(metric, levels = rev(metric))]

  panels[["i"]] <- ggplot(lto_df, aes(x = value, y = metric)) +
    geom_col(fill = "#4393c3", width = 0.6) +
    geom_vline(xintercept = lto$ref_auroc_v3_loco,
               linetype = "dashed", colour = "#b2182b") +
    scale_x_continuous(limits = c(0.5, 1.0)) +
    labs(x = "AUROC", y = NULL, title = "Leave-Two-Out Stress") +
    theme_pub()
} else {
  panels[["i"]] <- ggplot() + theme_void() + ggtitle("(i) No LTO data")
}

# ---- (j) Feature stability Jaccard ----
cat("  Panel (j): Feature stability (bootstrap CI)\n")
if (!is.null(bootstrap) && nrow(bootstrap) > 0) {
  # Show bootstrap CI as a single annotated panel
  bs <- bootstrap[1]
  bs_df <- data.table(
    x = 1,
    auroc = bs$observed_auroc,
    ci_lo = bs$ci_lower,
    ci_hi = bs$ci_upper,
    label = sprintf("AUROC = %.3f [%.3f, %.3f]",
                    bs$observed_auroc, bs$ci_lower, bs$ci_upper)
  )

  panels[["j"]] <- ggplot(bs_df, aes(x = x, y = auroc)) +
    geom_point(size = 3, colour = "#2166ac") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, linewidth = 0.5) +
    geom_text(aes(label = label), hjust = -0.2, size = 2.5) +
    scale_y_continuous(limits = c(0.6, 1.0)) +
    scale_x_continuous(limits = c(0.5, 2.0)) +
    labs(x = NULL, y = "AUROC", title = "Bootstrap 95% CI") +
    theme_pub() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
} else {
  panels[["j"]] <- ggplot() + theme_void() + ggtitle("(j) No bootstrap data")
}

# ---- (k) V2 vs V3 leakage ----
cat("  Panel (k): V2 vs V3 leakage comparison\n")
if (!is.null(v3v2) && nrow(v3v2) > 0) {
  v3v2_long <- melt(v3v2, id.vars = "target",
                    measure.vars = c("v3_auroc", "v2_auroc"),
                    variable.name = "version", value.name = "auroc")
  v3v2_long <- v3v2_long[!is.na(auroc)]
  v3v2_long[, version := fifelse(version == "v3_auroc", "V3 (proper)", "V2 (leaked)")]

  panels[["k"]] <- ggplot(v3v2_long, aes(x = target, y = auroc, fill = version)) +
    geom_col(position = position_dodge(width = 0.6), width = 0.5) +
    scale_fill_manual(values = c("V3 (proper)" = "#2166ac", "V2 (leaked)" = "#d6604d"),
                      name = "Version") +
    labs(x = NULL, y = "AUROC / QWK", title = "V2 (leaked) vs V3 (proper)") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 5))
} else {
  panels[["k"]] <- ggplot() + theme_void() + ggtitle("(k) No V2/V3 data")
}

# ---- (l) Panel size vs AUROC curve ----
cat("  Panel (l): Panel size vs AUROC\n")
if (!is.null(panel_curve) && nrow(panel_curve) > 0) {
  panels[["l"]] <- ggplot(panel_curve, aes(x = panel_size, y = auroc_5fold)) +
    geom_line(colour = "#762a83", linewidth = 0.7) +
    geom_point(colour = "#762a83", size = 2) +
    geom_errorbar(aes(ymin = auroc_5fold - fifelse(is.na(auroc_5fold_sd), 0, auroc_5fold_sd),
                      ymax = auroc_5fold + fifelse(is.na(auroc_5fold_sd), 0, auroc_5fold_sd)),
                  width = 0.05, linewidth = 0.3) +
    scale_x_log10() +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    labs(x = "Plasma Panel Size (proteins)", y = "AUROC",
         title = "Panel Size vs AUROC") +
    theme_pub()
} else {
  panels[["l"]] <- ggplot() + theme_void() + ggtitle("(l) No panel curve data")
}

# =========================================================================
# 4. Assemble and save
# =========================================================================
cat("\n--- Assembling figure ---\n")

# Compose layout:
#   a b c d
#   e f g h
#   i j k l
combined <- (
  panels[["a"]] + panels[["b"]] + panels[["c"]] + panels[["d"]] +
  panels[["e"]] + panels[["f"]] + panels[["g"]] + panels[["h"]] +
  panels[["i"]] + panels[["j"]] + panels[["k"]] + panels[["l"]]
) + plot_layout(ncol = 4, nrow = 3) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(face = "bold", size = 10))
  )

fig_path <- file.path(FIGDIR, "fig_model_comparison.pdf")
ggsave(fig_path, combined, width = 16, height = 12, units = "in", dpi = 300)
cat(sprintf("Saved figure: %s\n", fig_path))

# ---- Save tables ----
comp_path <- file.path(RDIR, "foundation_model_comparison.csv")
fwrite(comp_dt, comp_path)
cat(sprintf("Saved comparison table: %s (%d rows)\n", comp_path, nrow(comp_dt)))

# Summary: best model per target
if (nrow(comp_dt) > 0) {
  summary_dt <- comp_dt[, .(
    best_model_auroc = model[which.max(fifelse(is.na(metric_auroc), -Inf, metric_auroc))],
    best_auroc       = max(metric_auroc, na.rm = TRUE),
    best_model_qwk   = model[which.max(fifelse(is.na(metric_qwk), -Inf, metric_qwk))],
    best_qwk         = max(metric_qwk, na.rm = TRUE),
    n_models         = .N
  ), by = .(target)]

  summary_path <- file.path(RDIR, "foundation_model_summary.csv")
  fwrite(summary_dt, summary_path)
  cat(sprintf("Saved summary table: %s (%d targets)\n", summary_path, nrow(summary_dt)))

  cat("\n--- Best Model per Target ---\n")
  for (i in seq_len(nrow(summary_dt))) {
    r <- summary_dt[i]
    cat(sprintf("  %-25s  AUROC: %-30s (%.3f)  QWK: %-30s (%.3f)\n",
                r$target,
                r$best_model_auroc, fifelse(is.finite(r$best_auroc), r$best_auroc, NA_real_),
                r$best_model_qwk, fifelse(is.finite(r$best_qwk), r$best_qwk, NA_real_)))
  }
}

cat("\n=======================================================================\n")
cat("111_final_synthesis.R complete\n")
cat("=======================================================================\n")
