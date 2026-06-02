#!/usr/bin/env Rscript
# ============================================================================
# 155_prognosis_figures.R
# Prognosis Publication Figures — 8-panel composite figure for the prognosis
# (S1/S2 molecular subtype) analysis.
#
# Panels:
#   (a) PRS distribution by fibrosis stage — violin of P(S2) by F0-F4
#   (b) LOCO-CV AUROC barplot — 5 feature configs + XGBoost, error bars
#   (c) Risk stratification waterfall — F1-F2 sorted by P(S2)
#   (d) Clinical stratification — S2 proportion per risk quartile (F1-F2)
#   (e) Feature importance — top 20 stable features with coefficient bars
#   (f) HCC molecular score — boxplot by fibrosis stage
#   (g) Pseudo-survival curves — CPS quartile survival curves
#   (h) Panel performance — AUROC vs gene panel size
#
# Input:
#   results/prognosis/prs_all_predictions.csv
#   results/prognosis/prs_model_comparison.csv
#   results/prognosis/prs_loco_results.csv
#   results/prognosis/prs_feature_importance.csv
#   results/prognosis/prs_clinical_stratification.csv
#   results/prognosis/hcc_molecular_score.csv
#   results/prognosis/survival_curves_by_cps.csv
#   results/prognosis/panel_performance_curve.csv  (from Script 154)
#   results/staging_classifier/modeling_metadata.csv
#
# Output:
#   figures/prognosis/figS_prognosis.pdf          (composite)
#   figures/prognosis/panels/panel_{a-h}.pdf       (individual)
#
# SLURM: io partition, 4 CPUs, 8G RAM, 48h
# Env:   micromamba activate rnaseq
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# Publication theme
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

# Paths
INT <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
PROG_DIR <- file.path(INT, "results/prognosis")
STAGING_DIR <- file.path(INT, "results/staging_classifier")
FIG_DIR <- file.path(BASE, "figures/prognosis")
PANEL_DIR <- file.path(FIG_DIR, "panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 155_prognosis_figures.R ===\n")
cat("Date:", format(Sys.time()), "\n\n")

# ---------------------------------------------------------------------------
# Subtype colors (consistent with NMF subtype palette)
# ---------------------------------------------------------------------------
subtype_colors <- c(
  S1 = "#1565C0",   # Blue (stable)
  S2 = "#C2185B"    # Magenta (progressor)
)

# Config display names
config_labels <- c(
  M1_clinical   = "M1: Clinical",
  M2_celltype   = "M2: Cell-type",
  M3_divergence = "M3: Divergence",
  M4_embeddings = "M4: Embeddings",
  M5_full       = "M5: Full"
)

# ===================================================================
# Helper: safe read
# ===================================================================
safe_read <- function(path, desc = "") {
  if (!file.exists(path)) {
    warning(sprintf("File not found: %s (%s)", path, desc))
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE)
}

# ===================================================================
# Load data
# ===================================================================
cat("Loading data...\n")

predictions <- safe_read(file.path(PROG_DIR, "prs_all_predictions.csv"), "predictions")
model_comp  <- safe_read(file.path(PROG_DIR, "prs_model_comparison.csv"), "model comparison")
loco_res    <- safe_read(file.path(PROG_DIR, "prs_loco_results.csv"), "LOCO results")
feat_imp    <- safe_read(file.path(PROG_DIR, "prs_feature_importance.csv"), "feature importance")
clin_strat  <- safe_read(file.path(PROG_DIR, "prs_clinical_stratification.csv"), "stratification")
hcc_scores  <- safe_read(file.path(PROG_DIR, "hcc_molecular_score.csv"), "HCC scores")
surv_curves <- safe_read(file.path(PROG_DIR, "survival_curves_by_cps.csv"), "survival")
panel_perf  <- safe_read(file.path(PROG_DIR, "panel_performance_curve.csv"), "panel perf")
metadata    <- safe_read(file.path(STAGING_DIR, "modeling_metadata.csv"), "metadata")

# ===================================================================
# Panel (a): PRS distribution by fibrosis stage
# ===================================================================
cat("Panel (a): PRS by fibrosis stage...\n")
pa <- tryCatch({
  # Get M2_celltype elastic net predictions (best single model)
  pred_m2 <- predictions %>%
    filter(config == "M2_celltype", method == "elastic_net")

  # Merge with metadata for fibrosis stage
  pred_meta <- pred_m2 %>%
    inner_join(metadata %>% select(sample_id, fibrosis_stage), by = "sample_id") %>%
    filter(!is.na(fibrosis_stage), fibrosis_stage >= 0)

  pred_meta$fibrosis_label <- paste0("F", pred_meta$fibrosis_stage)
  pred_meta$fibrosis_label <- factor(pred_meta$fibrosis_label,
                                     levels = c("F0", "F1", "F2", "F3", "F4"))

  ggplot(pred_meta, aes(x = fibrosis_label, y = p_s2_predicted, fill = fibrosis_label)) +
    geom_violin(alpha = 0.7, scale = "width", linewidth = 0.3) +
    geom_boxplot(width = 0.15, outlier.size = 0.5, linewidth = 0.3, fill = "white",
                 alpha = 0.8) +
    scale_fill_manual(values = fibrosis_stage_colors, guide = "none") +
    labs(x = "Fibrosis stage", y = "P(S2 progressor)") +
    ggtitle("Progressor risk by fibrosis stage") +
    theme_masld()
}, error = function(e) {
  cat("  Panel (a) failed:", conditionMessage(e), "\n")
  placeholder("Panel a: PRS by fibrosis\n(data unavailable)")
})

# ===================================================================
# Panel (b): LOCO-CV AUROC barplot
# ===================================================================
cat("Panel (b): Model comparison barplot...\n")
pb <- tryCatch({
  comp <- model_comp %>%
    mutate(
      model_label = case_when(
        method == "xgboost" ~ paste0(config_labels[config], "\n(XGBoost)"),
        TRUE                ~ config_labels[config]
      ),
      model_label = factor(model_label,
        levels = c(config_labels, paste0(config_labels["M5_full"], "\n(XGBoost)")))
    ) %>%
    filter(!is.na(model_label))

  ggplot(comp, aes(x = model_label, y = mean_auroc, fill = config)) +
    geom_col(width = 0.7, alpha = 0.85, show.legend = FALSE) +
    geom_errorbar(aes(ymin = mean_auroc - sd_auroc, ymax = mean_auroc + sd_auroc),
                  width = 0.2, linewidth = 0.3) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_fill_manual(values = c(
      M1_clinical = "#BDBDBD", M2_celltype = "#C2185B",
      M3_divergence = "#7B1FA2", M4_embeddings = "#42A5F5",
      M5_full = "#0D47A1"
    )) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "LOCO-CV AUROC") +
    ggtitle("S1/S2 classification performance") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 5))
}, error = function(e) {
  cat("  Panel (b) failed:", conditionMessage(e), "\n")
  placeholder("Panel b: AUROC comparison\n(data unavailable)")
})

# ===================================================================
# Panel (c): Risk stratification waterfall (F1-F2)
# ===================================================================
cat("Panel (c): Waterfall plot...\n")
pc <- tryCatch({
  # M2_celltype predictions for F1-F2 patients
  pred_m2 <- predictions %>%
    filter(config == "M2_celltype", method == "elastic_net")

  pred_f12 <- pred_m2 %>%
    inner_join(metadata %>% select(sample_id, fibrosis_stage), by = "sample_id") %>%
    filter(fibrosis_stage %in% c(1, 2)) %>%
    arrange(p_s2_predicted) %>%
    mutate(
      rank = row_number(),
      subtype = ifelse(s2_binary_true == 1, "S2", "S1")
    )

  ggplot(pred_f12, aes(x = rank, y = p_s2_predicted, fill = subtype)) +
    geom_col(width = 1, linewidth = 0) +
    scale_fill_manual(values = subtype_colors, name = "True subtype") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray30", linewidth = 0.3) +
    labs(x = "Patients (F1-F2, ranked)", y = "P(S2)") +
    ggtitle("Risk stratification waterfall (F1-F2)") +
    theme_masld() +
    theme(legend.position = c(0.15, 0.85))
}, error = function(e) {
  cat("  Panel (c) failed:", conditionMessage(e), "\n")
  placeholder("Panel c: Waterfall\n(data unavailable)")
})

# ===================================================================
# Panel (d): Clinical stratification — S2 proportion per quartile
# ===================================================================
cat("Panel (d): Clinical stratification...\n")
pd_plot <- tryCatch({
  strat <- clin_strat %>%
    filter(grepl("^Q", quartile)) %>%
    mutate(
      quartile = factor(quartile, levels = c("Q1_low", "Q2", "Q3", "Q4_high")),
      quartile_label = c(Q1_low = "Q1\n(low risk)", Q2 = "Q2", Q3 = "Q3",
                         Q4_high = "Q4\n(high risk)")[as.character(quartile)]
    )

  ggplot(strat, aes(x = quartile_label, y = s2_fraction)) +
    geom_col(aes(fill = s2_fraction), width = 0.7, show.legend = FALSE) +
    scale_fill_gradient(low = "#1565C0", high = "#C2185B") +
    geom_text(aes(label = sprintf("%.0f%%", s2_fraction * 100)),
              vjust = -0.5, size = 2) +
    labs(x = "Risk quartile (F1-F2)", y = "Fraction S2 progressor") +
    scale_y_continuous(limits = c(0, 0.55), expand = expansion(mult = c(0, 0.05))) +
    ggtitle("S2 enrichment by risk quartile") +
    theme_masld()
}, error = function(e) {
  cat("  Panel (d) failed:", conditionMessage(e), "\n")
  placeholder("Panel d: Risk quartiles\n(data unavailable)")
})

# ===================================================================
# Panel (e): Feature importance — top 20 stable features
# ===================================================================
cat("Panel (e): Feature importance...\n")
pe <- tryCatch({
  # Use M5_full features for the most comprehensive view
  imp_full <- feat_imp %>%
    filter(config == "M5_full") %>%
    arrange(desc(n_folds_nonzero), desc(abs_mean_coefficient)) %>%
    head(20) %>%
    mutate(
      display_name = sub("^div_", "", feature),
      display_name = sub("^ct_", "CT: ", display_name),
      display_name = sub("^clin_", "Clin: ", display_name),
      display_name = sub("^nasvae_", "VAE: ", display_name),
      display_name = sub("^bulkformer_", "BF: ", display_name),
      display_name = sub("^trans_", "Trans: ", display_name),
      display_name = sub("^hcc_", "HCC: ", display_name),
      feature_type = case_when(
        grepl("^div_", feature) ~ "Divergence gene",
        grepl("^ct_", feature)  ~ "Cell-type fraction",
        grepl("^clin_", feature) ~ "Clinical",
        grepl("^nasvae_|^bulkformer_", feature) ~ "Embedding",
        grepl("^trans_", feature) ~ "Transition",
        grepl("^hcc_", feature) ~ "HCC score",
        TRUE ~ "Other"
      )
    )

  imp_full$display_name <- factor(imp_full$display_name,
                                  levels = rev(imp_full$display_name))

  ggplot(imp_full, aes(x = mean_coefficient, y = display_name, fill = feature_type)) +
    geom_col(alpha = 0.85, width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "gray30") +
    scale_fill_manual(values = c(
      "Divergence gene"  = "#C2185B",
      "Cell-type fraction" = "#7B1FA2",
      "Clinical"         = "#BDBDBD",
      "Embedding"        = "#42A5F5",
      "Transition"       = "#00695C",
      "HCC score"        = "#F57F17",
      "Other"            = "#78909C"
    ), name = "Feature type") +
    labs(x = "Mean coefficient", y = NULL) +
    ggtitle("Top 20 prognostic features (M5)") +
    theme_masld() +
    theme(legend.position = "right")
}, error = function(e) {
  cat("  Panel (e) failed:", conditionMessage(e), "\n")
  placeholder("Panel e: Feature importance\n(data unavailable)")
})

# ===================================================================
# Panel (f): HCC molecular score by fibrosis stage
# ===================================================================
cat("Panel (f): HCC molecular score...\n")
pf <- tryCatch({
  hcc <- hcc_scores %>%
    filter(!is.na(fibrosis_stage), fibrosis_stage >= 0) %>%
    mutate(
      fibrosis_label = paste0("F", fibrosis_stage),
      fibrosis_label = factor(fibrosis_label, levels = c("F0", "F1", "F2", "F3", "F4")),
      subtype = ifelse(s2_binary == 1, "S2", "S1")
    )

  ggplot(hcc, aes(x = fibrosis_label, y = hcc_molecular_score)) +
    geom_boxplot(aes(fill = fibrosis_label), outlier.size = 0.5, linewidth = 0.3,
                 alpha = 0.7, show.legend = FALSE) +
    geom_jitter(aes(color = subtype), width = 0.15, size = 0.3, alpha = 0.4) +
    scale_fill_manual(values = fibrosis_stage_colors) +
    scale_color_manual(values = subtype_colors, name = "Subtype") +
    labs(x = "Fibrosis stage", y = "HCC molecular score") +
    ggtitle("HCC risk signature by stage") +
    theme_masld() +
    theme(legend.position = c(0.15, 0.85))
}, error = function(e) {
  cat("  Panel (f) failed:", conditionMessage(e), "\n")
  placeholder("Panel f: HCC score\n(data unavailable)")
})

# ===================================================================
# Panel (g): Pseudo-survival curves by CPS quartile
# ===================================================================
cat("Panel (g): Pseudo-survival curves...\n")
pg <- tryCatch({
  surv <- surv_curves %>%
    mutate(
      cps_quartile = factor(cps_quartile,
                            levels = c("Q1_low", "Q2", "Q3", "Q4_high"))
    )

  cps_colors <- c(
    Q1_low  = "#1565C0",
    Q2      = "#42A5F5",
    Q3      = "#E91E63",
    Q4_high = "#880E4F"
  )

  ggplot(surv, aes(x = pseudotime_bin, y = survival_fraction,
                   color = cps_quartile, group = cps_quartile)) +
    geom_step(linewidth = 0.6) +
    scale_color_manual(values = cps_colors,
                       labels = c("Q1 (low risk)", "Q2", "Q3", "Q4 (high risk)"),
                       name = "CPS quartile") +
    labs(x = "Pseudotime (disease progression)",
         y = "Fraction remaining non-fibrotic") +
    scale_y_continuous(limits = c(0.7, 1.02), expand = expansion(mult = c(0, 0))) +
    ggtitle("Pseudo-survival by progression score") +
    theme_masld() +
    theme(legend.position = c(0.25, 0.25))
}, error = function(e) {
  cat("  Panel (g) failed:", conditionMessage(e), "\n")
  placeholder("Panel g: Survival curves\n(data unavailable)")
})

# ===================================================================
# Panel (h): Panel performance — AUROC vs gene panel size
# ===================================================================
cat("Panel (h): Panel performance curve...\n")
ph <- tryCatch({
  if (is.null(panel_perf)) {
    stop("panel_performance_curve.csv not found (Script 154 may not have run yet)")
  }

  pp <- panel_perf %>%
    filter(panel_size <= 100)

  elbow_row <- pp %>% filter(is_elbow == TRUE)

  # Reference line: full divergence model AUROC
  full_auroc <- pp %>% filter(panel_size == max(panel_size)) %>% pull(mean_auroc)
  threshold <- full_auroc - 0.05

  ggplot(pp, aes(x = panel_size, y = mean_auroc)) +
    geom_ribbon(aes(ymin = mean_auroc - sd_auroc, ymax = mean_auroc + sd_auroc),
                fill = "#C2185B", alpha = 0.15) +
    geom_line(color = "#C2185B", linewidth = 0.6) +
    geom_point(color = "#C2185B", size = 1.5) +
    {if (nrow(elbow_row) > 0) geom_point(data = elbow_row, color = "#0D47A1",
                                          size = 3, shape = 18)} +
    geom_hline(yintercept = full_auroc, linetype = "solid",
               color = "gray50", linewidth = 0.3) +
    geom_hline(yintercept = threshold, linetype = "dashed",
               color = "gray50", linewidth = 0.3) +
    annotate("text", x = max(pp$panel_size) * 0.7, y = full_auroc + 0.015,
             label = "Full model", size = 2, color = "gray40") +
    annotate("text", x = max(pp$panel_size) * 0.7, y = threshold - 0.015,
             label = "90% threshold", size = 2, color = "gray40") +
    scale_x_continuous(breaks = c(5, 10, 15, 25, 50, 100)) +
    labs(x = "Number of genes in panel", y = "LOCO-CV AUROC") +
    ggtitle("Minimal prognosis panel performance") +
    theme_masld()
}, error = function(e) {
  cat("  Panel (h) failed:", conditionMessage(e), "\n")
  placeholder("Panel h: Panel performance\n(awaiting Script 154)")
})

# ===================================================================
# Assemble composite figure
# ===================================================================
cat("\nAssembling composite figure...\n")

composite <- (pa | pb) / (pc | pd_plot) / (pe | pf) / (pg | ph) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(size = 8, face = "bold"))
  )

save_fig_tall(composite,
              file.path(FIG_DIR, "figS_prognosis.pdf"),
              width = fig_full_width,
              height = 10)
cat("  Saved figS_prognosis.pdf\n")

# ===================================================================
# Save individual panels
# ===================================================================
cat("Saving individual panels...\n")

panel_list <- list(
  panel_a = pa,
  panel_b = pb,
  panel_c = pc,
  panel_d = pd_plot,
  panel_e = pe,
  panel_f = pf,
  panel_g = pg,
  panel_h = ph
)

for (name in names(panel_list)) {
  out_path <- file.path(PANEL_DIR, paste0(name, ".pdf"))
  tryCatch({
    save_fig(panel_list[[name]], out_path, width = fig_half_width, height = 3)
    cat(sprintf("  Saved %s\n", out_path))
  }, error = function(e) {
    cat(sprintf("  Failed to save %s: %s\n", name, conditionMessage(e)))
  })
}

cat("\n=== 155_prognosis_figures.R COMPLETE ===\n")
