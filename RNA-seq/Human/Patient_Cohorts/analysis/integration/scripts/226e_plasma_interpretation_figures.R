#!/usr/bin/env Rscript
# QUARANTINED 2026-06-01: plots supervised Olink plasma-classifier results that are WITHDRAWN (unrecoverable subject-label join). Do not regenerate for the manuscript. See docs/reviews/2026-06-01-proteomics-extreme-review.md.
# 226e_plasma_interpretation_figures.R
# ---------------------------------------------------------------------------
# Publication figures for plasma proteomics interpretation (Script 226 series).
#
# Panels:
#   A  — Top-20 protein importance lollipop (T5 etiology binary)
#   B  — Per-class SHAP bar chart (T4 etiology 3-class, top 10 per class)
#   C  — Tissue-plasma bridge evidence heatmap (top-20 × 3 targets)
#   D  — Panel-size AUROC / F1 curve (T5 + T1, SHAP vs random)
#   E  — Cross-platform transfer waterfall
#   F  — Effect size concordance scatter (Olink vs DIA-MS)
#   G  — Pathway NES heatmap (Hallmark, top 15 by mean |NES|)
#
# Inputs  (INDIR = results/multiprogram/):
#   226a_protein_importance.csv
#   226a_panel_curves.csv
#   226b_tissue_bridge.csv
#   226b_enrichment.csv
#   226c_cross_platform.csv
#   226c_effect_concordance.csv
#   226d_etiology_pathways.csv
#   226d_etiology_protein_lists.csv
#
# Outputs (OUTDIR = figures/prediction/):
#   226e_panel_A.pdf … 226e_panel_G.pdf
#   226e_plasma_interpretation.pdf  (composite)
#
# Usage: Rscript 226e_plasma_interpretation_figures.R
# SLURM: --partition=cpu --cpus-per-task=4 --mem=16G --time=02:00:00
# Env:   micromamba activate rnaseq
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INDIR  <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram")
OUTDIR <- file.path(BASE, "figures/prediction")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Source publication theme (defines theme_masld, masld_colors, save_fig)
source(file.path(BASE, "scripts/figures/publication_theme.R"))

cat("=== 226e: Plasma Interpretation Figures ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ---------------------------------------------------------------------------
# Helper: save individual panel
# ---------------------------------------------------------------------------
save_panel <- function(p, name, width, height) {
  path <- file.path(OUTDIR, paste0("226e_panel_", name, ".pdf"))
  pdf_dev <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
  ggsave(path, p, width = width, height = height, device = pdf_dev)
  cat("  Saved:", path, "\n")
}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
cat("Loading input CSVs...\n")

importance  <- fread(file.path(INDIR, "226a_protein_importance.csv"))
panel_curves <- fread(file.path(INDIR, "226a_panel_curves.csv"))
bridge      <- fread(file.path(INDIR, "226b_tissue_bridge.csv"))
enrichment  <- fread(file.path(INDIR, "226b_enrichment.csv"))
cross_plat  <- fread(file.path(INDIR, "226c_cross_platform.csv"))
effect_conc <- fread(file.path(INDIR, "226c_effect_concordance.csv"))
pathways    <- fread(file.path(INDIR, "226d_etiology_pathways.csv"))
prot_lists  <- fread(file.path(INDIR, "226d_etiology_protein_lists.csv"))

cat("  importance rows:", nrow(importance), "\n")
cat("  panel_curves rows:", nrow(panel_curves), "\n")
cat("  bridge rows:", nrow(bridge), "\n")
cat("  cross_plat rows:", nrow(cross_plat), "\n")
cat("  effect_conc rows:", nrow(effect_conc), "\n")
cat("  pathways rows:", nrow(pathways), "\n")
cat("  prot_lists rows:", nrow(prot_lists), "\n\n")

# ---------------------------------------------------------------------------
# Panel A: Top-20 protein importance lollipop — T5 etiology binary
# ---------------------------------------------------------------------------
cat("Building Panel A: Top-20 lollipop (T5)...\n")

bridge_t5 <- bridge[target == "T5_etiology_binary"][order(shap_rank)]
# Keep rank <= 20; fall back to importance CSV if bridge is smaller
if (nrow(bridge_t5) == 0 || max(bridge_t5$shap_rank, na.rm = TRUE) < 20) {
  cat("  NOTE: using 226a_protein_importance.csv for T5 (bridge has <20 rows)\n")
  imp_t5 <- importance[target == "T5_etiology_binary"][order(rank)][seq_len(20)]
  # Construct minimal bridge-like table
  bridge_t5 <- imp_t5[, .(
    protein      = protein,
    shap_rank    = rank,
    fold_stability = fold_stability,
    mean_abs_shap  = mean_abs_shap,
    is_deg         = FALSE,
    f2_switch_class = "none"
  )]
} else {
  bridge_t5 <- bridge_t5[shap_rank <= 20]
}

# Coerce types
bridge_t5[, is_deg := as.logical(is_deg)]
if (!"f2_switch_class" %in% names(bridge_t5)) {
  bridge_t5[, f2_switch_class := "none"]
}
bridge_t5[is.na(f2_switch_class), f2_switch_class := "none"]
bridge_t5[, f2_switch_class := as.character(f2_switch_class)]

# Order proteins by importance (ascending so most important is at top in coord_flip)
# Use unique() to prevent duplicate factor levels when a protein appears in multiple targets
protein_order <- unique(bridge_t5[order(shap_rank), protein])
bridge_t5[, protein := factor(protein, levels = rev(protein_order))]

# Stability label
bridge_t5[, stab_label := paste0(round(fold_stability * 100), "%")]

# Map f2 class to shape
shape_map <- c("none" = 21, "switch" = 24, "onset" = 23, "late" = 23)
shape_vals <- shape_map[as.character(bridge_t5$f2_switch_class)]
shape_vals[is.na(shape_vals)] <- 21L

bridge_t5[, pt_shape := shape_vals]

# Colors: is_deg
col_deg    <- masld_colors$up        # "#C2185B"
col_nodeg  <- masld_colors$ns        # "#BDBDBD"
bridge_t5[, pt_color := ifelse(is_deg == TRUE, col_deg, col_nodeg)]

p_A <- ggplot(bridge_t5, aes(x = protein, y = mean_abs_shap)) +
  geom_segment(aes(xend = protein, y = 0, yend = mean_abs_shap),
               color = "gray70", linewidth = 0.4) +
  geom_point(aes(fill = is_deg, shape = f2_switch_class),
             size = 2, color = "gray30", stroke = 0.3) +
  geom_text(aes(y = max(mean_abs_shap, na.rm = TRUE) * 1.05,
                label = stab_label),
            hjust = 0, size = 2, color = "gray40") +
  scale_fill_manual(
    name   = "Tissue DEG",
    values = c("TRUE" = col_deg, "FALSE" = col_nodeg),
    labels = c("TRUE" = "Yes", "FALSE" = "No")
  ) +
  scale_shape_manual(
    name   = "F2 switch class",
    values = c("none" = 21, "switch" = 24, "onset" = 23, "late" = 22),
    labels = c("none" = "None", "switch" = "Switch", "onset" = "Onset", "late" = "Late")
  ) +
  coord_flip() +
  expand_limits(y = max(bridge_t5$mean_abs_shap, na.rm = TRUE) * 1.2) +
  labs(
    title = "Top-20 proteins: MASLD vs other etiology",
    x     = NULL,
    y     = "Mean |SHAP| (permutation importance)"
  ) +
  theme_masld(base_size = 9) +
  theme(
    legend.position = "right",
    legend.box      = "vertical"
  )

save_panel(p_A, "A", width = fig_half_width + 0.5, height = 4.5)

# ---------------------------------------------------------------------------
# Panel B: Per-class SHAP bar chart — T4 etiology 3-class, top 10 per class
# ---------------------------------------------------------------------------
cat("Building Panel B: Per-class SHAP bars (T4)...\n")

if (nrow(prot_lists) == 0) {
  cat("  WARNING: 226d_etiology_protein_lists.csv is empty — skipping Panel B\n")
  p_B <- placeholder("Panel B: No data")
} else {
  # Build long table: top 10 per class by |SHAP|, with signed SHAP value
  classes <- c("MASLD", "CVH", "ARLD")
  col_abs  <- c("MASLD" = "masld_abs_shap", "CVH" = "cvh_abs_shap", "ARLD" = "arld_abs_shap")
  col_sign <- c("MASLD" = "masld_shap",     "CVH" = "cvh_shap",     "ARLD" = "arld_shap")

  b_list <- lapply(classes, function(cl) {
    abs_col  <- col_abs[cl]
    sign_col <- col_sign[cl]
    if (!abs_col %in% names(prot_lists) || !sign_col %in% names(prot_lists)) {
      return(NULL)
    }
    dt <- copy(prot_lists)[order(-get(abs_col))][seq_len(min(10, .N))]
    data.table(
      protein        = dt$protein,
      etiology_class = cl,
      abs_shap       = dt[[abs_col]],
      signed_shap    = dt[[sign_col]]
    )
  })
  b_long <- rbindlist(b_list[!sapply(b_list, is.null)])

  # Order proteins within each facet by abs_shap descending
  b_long[, protein_f := reorder(protein, abs_shap)]

  p_B <- ggplot(b_long, aes(x = signed_shap, y = protein_f, fill = signed_shap)) +
    geom_col(width = 0.7) +
    scale_fill_gradient2(
      low      = masld_colors$down,
      mid      = "white",
      high     = masld_colors$up,
      midpoint = 0,
      name     = "Mean SHAP"
    ) +
    facet_wrap(~etiology_class, scales = "free_y", ncol = 3) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "black") +
    labs(
      title = "Etiology-specific protein programs (T4, top 10 per class)",
      x     = "Mean signed SHAP value",
      y     = NULL
    ) +
    theme_masld(base_size = 9) +
    theme(
      legend.position = "bottom",
      axis.text.y     = element_text(size = 7)
    )
}

save_panel(p_B, "B", width = fig_full_width, height = 4.5)

# ---------------------------------------------------------------------------
# Panel C: Tissue-plasma bridge evidence dot/tile plot
# ---------------------------------------------------------------------------
cat("Building Panel C: Tissue-bridge evidence heatmap...\n")

# Define evidence columns to show (handle missing columns gracefully)
ev_cols <- c(
  "is_deg"           = "Tissue DEG",
  "is_conserved" = "Conserved",
  "f2_switch_class"  = "F2 switch",
  "dgidb_druggable"  = "DGIdb druggable",
  "opentargets_drug" = "OpenTargets drug"
)
present_cols <- intersect(names(ev_cols), names(bridge))
# Drop columns where all values are NA
present_cols <- present_cols[sapply(present_cols, function(cc) {
  !all(is.na(bridge[[cc]]))
})]

if (length(present_cols) == 0) {
  cat("  WARNING: no evidence columns available — skipping Panel C\n")
  p_C <- placeholder("Panel C: No evidence columns")
} else {
  # Build long format for tile plot — top-20 per target
  c_list <- lapply(unique(bridge$target), function(tgt) {
    dt <- bridge[target == tgt][order(shap_rank)][seq_len(min(20, .N))]
    # For each evidence column, convert to a display value
    rows <- lapply(present_cols, function(cc) {
      val_raw <- as.character(dt[[cc]])
      # Numeric → binary threshold
      if (cc == "f2_switch_class") {
        # Keep as categorical string: none/switch/onset/late
        val_raw[is.na(val_raw)] <- "none"
        value_label <- val_raw
      } else {
        val_raw[is.na(val_raw)] <- "FALSE"
        val_raw[val_raw %in% c("True", "TRUE", "1")] <- "Yes"
        val_raw[val_raw %in% c("False", "FALSE", "0")] <- "No"
        value_label <- val_raw
      }
      data.table(
        protein     = dt$protein,
        shap_rank   = dt$shap_rank,
        target      = tgt,
        evidence    = ev_cols[cc],
        value_label = value_label
      )
    })
    rbindlist(rows)
  })
  c_long <- rbindlist(c_list)

  # Assign numeric score for tile color
  c_long[, score := fcase(
    value_label == "Yes"    , 1,
    value_label == "switch" , 0.8,
    value_label == "onset"  , 0.6,
    value_label == "late"   , 0.4,
    value_label == "No"     , 0,
    value_label == "none"   , 0,
    default = NA_real_
  )]

  # Target labels for facets
  target_labels <- c(
    "T1_binary"        = "Fibrosis F≥2",
    "T4_etiology_3"    = "Etiology 3-class",
    "T5_etiology_binary" = "MASLD vs other"
  )
  c_long[, target_label := target_labels[target]]
  c_long[is.na(target_label), target_label := target]

  # Order proteins by shap_rank within each target
  c_long[, protein_f := factor(protein,
    levels = rev(unique(c_long[order(target, shap_rank), protein])))]

  p_C <- ggplot(c_long,
    aes(x = evidence, y = protein_f, fill = score)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient(
      low      = "gray92",
      high     = masld_colors$up,
      na.value = "gray85",
      name     = "Evidence\npresence",
      limits   = c(0, 1)
    ) +
    facet_wrap(~target_label, scales = "free_y", ncol = 3) +
    labs(
      title = "Tissue-plasma bridge: evidence layers (top-20 per target)",
      x     = NULL,
      y     = NULL
    ) +
    theme_masld(base_size = 9) +
    theme(
      axis.text.x     = element_text(angle = 35, hjust = 1, size = 7),
      axis.text.y     = element_text(size = 6.5),
      legend.position = "right"
    )
}

save_panel(p_C, "C", width = fig_full_width, height = 6)

# ---------------------------------------------------------------------------
# Panel D: Panel-size AUROC curve — T5 + T1
# ---------------------------------------------------------------------------
cat("Building Panel D: Panel-size curves (T5 + T1)...\n")

# Hardcoded NFASC+GDF15 benchmarks (size=2 known biomarker combo)
benchmarks <- data.table(
  target  = c("T1_binary", "T5_etiology_binary"),
  bench   = c(0.734,        NA_real_)
)

targets_d <- c("T1_binary", "T5_etiology_binary")
pc_sub <- panel_curves[target %in% targets_d]

if (nrow(pc_sub) == 0) {
  cat("  WARNING: no panel curve data for T1/T5 — skipping Panel D\n")
  p_D <- placeholder("Panel D: No data")
} else {
  # SHAP-selected line: is_random == FALSE
  shap_line <- pc_sub[is_random == FALSE,
    .(metric_value = mean(metric_value, na.rm = TRUE)),
    by = .(target, panel_size, metric_name)]

  # Random ribbon: mean ± SD across draw_id × fold
  rand_ribbon <- pc_sub[is_random == TRUE,
    .(mean_m = mean(metric_value, na.rm = TRUE),
      sd_m   = sd(metric_value, na.rm = TRUE)),
    by = .(target, panel_size, metric_name)]

  # Metric label per target
  metric_labels <- c(
    "auroc"    = "AUROC",
    "macro_f1" = "Macro F1"
  )
  shap_line[,  metric_label := metric_labels[metric_name]]
  rand_ribbon[, metric_label := metric_labels[metric_name]]

  target_labels_d <- c(
    "T1_binary"          = "Fibrosis F≥2 (AUROC)",
    "T5_etiology_binary" = "MASLD vs other (AUROC)"
  )
  shap_line[,  target_label := target_labels_d[target]]
  rand_ribbon[, target_label := target_labels_d[target]]

  # Panel sizes for x axis
  size_breaks <- c(2, 5, 10, 20, 50, 100, 200, 1461)

  p_D <- ggplot() +
    # Random ribbon
    geom_ribbon(
      data    = rand_ribbon,
      aes(x = panel_size, ymin = mean_m - sd_m, ymax = mean_m + sd_m),
      fill    = "gray80", alpha = 0.5
    ) +
    # Random mean line (dashed)
    geom_line(
      data  = rand_ribbon,
      aes(x = panel_size, y = mean_m),
      color = "gray60", linetype = "dashed", linewidth = 0.5
    ) +
    # SHAP-selected line
    geom_line(
      data  = shap_line,
      aes(x = panel_size, y = metric_value),
      color = masld_colors$down, linewidth = 0.7
    ) +
    geom_point(
      data  = shap_line,
      aes(x = panel_size, y = metric_value),
      color = masld_colors$down, size = 1.8
    ) +
    # NFASC+GDF15 benchmark (T1 only)
    geom_hline(
      data    = benchmarks[target == "T1_binary" & !is.na(bench)][
        , target_label := target_labels_d["T1_binary"]],
      aes(yintercept = bench),
      color   = masld_colors$up, linetype = "dashed", linewidth = 0.5
    ) +
    scale_x_log10(
      breaks = size_breaks,
      labels = as.character(size_breaks)
    ) +
    facet_wrap(~target_label, scales = "free_y", ncol = 2) +
    labs(
      title = "SHAP-selected panel vs random (shaded = random mean ± SD)",
      x     = "Panel size (proteins)",
      y     = "Metric value"
    ) +
    theme_masld(base_size = 9) +
    theme(legend.position = "none")
}

save_panel(p_D, "D", width = fig_full_width, height = 3)

# ---------------------------------------------------------------------------
# Panel E: Cross-platform transfer waterfall
# ---------------------------------------------------------------------------
cat("Building Panel E: Cross-platform waterfall...\n")

# Exclude effect_concordance row; keep AUROC rows only
cp_plot <- cross_plat[!experiment %in% c("effect_concordance")]
if (nrow(cp_plot) == 0) {
  cat("  WARNING: cross_plat has no usable rows — skipping Panel E\n")
  p_E <- placeholder("Panel E: No data")
} else {
  cp_plot[, sd_auroc := suppressWarnings(as.numeric(sd_auroc))]
  cp_plot[, auroc    := suppressWarnings(as.numeric(auroc))]

  # Assign display labels and group
  label_map <- c(
    "full_olink"      = "Full Olink\n(N=1,461)",
    "shared_olink"    = "Shared proteins\n(N=48)",
    "transfer_xgb"    = "Transfer XGBoost\n(DIA-MS, N=72)",
    "transfer_logreg" = "Transfer LogReg\n(DIA-MS, N=72)"
  )
  platform_map <- c(
    "full_olink"      = "Within-platform",
    "shared_olink"    = "Within-platform",
    "transfer_xgb"    = "Cross-platform",
    "transfer_logreg" = "Cross-platform"
  )
  cp_plot[, label    := label_map[experiment]]
  cp_plot[is.na(label),    label    := experiment]
  cp_plot[, platform := platform_map[experiment]]
  cp_plot[is.na(platform), platform := "Cross-platform"]

  # Ordered factor for x axis (best to worst)
  cp_plot[, label := factor(label, levels = rev(cp_plot[order(-auroc), label]))]

  platform_colors <- c(
    "Within-platform" = masld_colors$down,
    "Cross-platform"  = "#E65100"   # deep orange
  )

  p_E <- ggplot(cp_plot, aes(x = label, y = auroc, fill = platform)) +
    geom_col(width = 0.55) +
    geom_errorbar(
      data = cp_plot[!is.na(sd_auroc)],
      aes(ymin = auroc - sd_auroc, ymax = auroc + sd_auroc),
      width = 0.2, linewidth = 0.4, color = "gray30"
    ) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               linewidth = 0.4, color = "gray50") +
    scale_fill_manual(values = platform_colors, name = NULL) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = number_format(accuracy = 0.1)
    ) +
    coord_flip() +
    labs(
      title = "Cross-platform transfer: MASLD vs other etiology (T5)",
      x     = NULL,
      y     = "AUROC"
    ) +
    theme_masld(base_size = 9) +
    theme(legend.position = "bottom")
}

save_panel(p_E, "E", width = fig_half_width + 0.5, height = 2.8)

# ---------------------------------------------------------------------------
# Panel F: Effect size concordance scatter (Olink vs DIA-MS)
# ---------------------------------------------------------------------------
cat("Building Panel F: Effect size concordance scatter...\n")

if (nrow(effect_conc) == 0 ||
    !all(c("olink_d", "diams_d") %in% names(effect_conc))) {
  cat("  WARNING: effect_conc missing required columns — skipping Panel F\n")
  p_F <- placeholder("Panel F: No data")
} else {
  # Significance classification
  has_pvals <- all(c("olink_pvalue", "diams_pvalue") %in% names(effect_conc))
  if (has_pvals) {
    effect_conc[, sig_class := fcase(
      olink_pvalue < 0.05 & diams_pvalue < 0.05, "Both significant",
      olink_pvalue < 0.05 | diams_pvalue < 0.05, "One significant",
      default                                     = "Neither"
    )]
  } else {
    effect_conc[, sig_class := "Neither"]
  }

  sig_colors <- c(
    "Both significant" = masld_colors$up,
    "One significant"  = "#E65100",
    "Neither"          = masld_colors$ns
  )

  # Compute Spearman rho from the data (supplement with cross_plat annotation)
  rho_val <- round(cor(effect_conc$olink_d, effect_conc$diams_d,
    method = "spearman", use = "complete.obs"), 3)
  pval <- tryCatch(
    cor.test(effect_conc$olink_d, effect_conc$diams_d,
             method = "spearman", exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  pval_label <- if (!is.na(pval) && pval < 0.001) {
    sprintf("rho = %.3f, p < 0.001", rho_val)
  } else if (!is.na(pval)) {
    sprintf("rho = %.3f, p = %.3f", rho_val, pval)
  } else {
    sprintf("rho = %.3f", rho_val)
  }

  lim_range <- range(c(effect_conc$olink_d, effect_conc$diams_d), na.rm = TRUE)
  lim_exp   <- extendrange(lim_range, f = 0.05)

  p_F <- ggplot(effect_conc, aes(x = olink_d, y = diams_d, color = sig_class)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3,
                linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "gray80") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "gray80") +
    geom_point(size = 1.8, alpha = 0.85) +
    annotate("text",
             x     = lim_exp[1] + diff(lim_exp) * 0.02,
             y     = lim_exp[2] - diff(lim_exp) * 0.04,
             label = pval_label, hjust = 0,
             size  = 2.5, color = "gray30") +
    scale_color_manual(values = sig_colors, name = "Significance") +
    coord_equal(xlim = lim_exp, ylim = lim_exp) +
    labs(
      title = "Effect size concordance: Olink vs DIA-MS (shared proteins)",
      x     = "Olink Cohen's d (F≥2 vs F<2)",
      y     = "DIA-MS Cohen's d"
    ) +
    theme_masld(base_size = 9) +
    theme(legend.position = "bottom")
}

save_panel(p_F, "F", width = fig_half_width + 0.5, height = 3.5)

# ---------------------------------------------------------------------------
# Panel G: Pathway NES heatmap — Hallmark, top 15 by mean |NES|
# ---------------------------------------------------------------------------
cat("Building Panel G: Pathway NES heatmap...\n")

if (nrow(pathways) == 0 || !"collection" %in% names(pathways)) {
  cat("  WARNING: pathways data missing or no 'collection' column — skipping Panel G\n")
  p_G <- placeholder("Panel G: No data")
} else {
  hall <- pathways[collection == "Hallmark"]

  if (nrow(hall) == 0) {
    cat("  WARNING: no Hallmark pathways found — using all collections\n")
    hall <- pathways
  }

  # Pivot to wide: pathway × etiology_class → NES
  if (!"etiology_class" %in% names(hall)) {
    cat("  WARNING: 'etiology_class' column missing in pathways — skipping Panel G\n")
    p_G <- placeholder("Panel G: Missing etiology_class")
  } else {
    hall_wide <- dcast(hall, pathway ~ etiology_class, value.var = "NES",
                       fun.aggregate = mean, na.rm = TRUE)

    nes_cols <- setdiff(names(hall_wide), "pathway")
    if (length(nes_cols) == 0) {
      p_G <- placeholder("Panel G: No NES columns after pivot")
    } else {
      # Compute mean |NES| across classes for ranking
      nes_mat <- as.matrix(hall_wide[, nes_cols, with = FALSE])
      hall_wide[, mean_abs_nes := rowMeans(abs(nes_mat), na.rm = TRUE)]
      hall_top <- hall_wide[order(-mean_abs_nes)][seq_len(min(15, .N))]

      # Clean pathway name (remove HALLMARK_ prefix, replace _ with space)
      hall_top[, pathway_clean := gsub("^HALLMARK_", "", pathway)]
      hall_top[, pathway_clean := gsub("_", " ", pathway_clean)]
      hall_top[, pathway_clean := stringr::str_to_title(pathway_clean)]

      # Melt back to long
      g_long <- melt(hall_top,
        id.vars     = c("pathway", "pathway_clean", "mean_abs_nes"),
        measure.vars = nes_cols,
        variable.name = "etiology_class",
        value.name    = "NES"
      )

      # Order pathways by mean |NES|
      g_long[, pathway_f := factor(pathway_clean,
        levels = hall_top[order(mean_abs_nes), pathway_clean])]

      nes_lim <- max(abs(g_long$NES), na.rm = TRUE)

      p_G <- ggplot(g_long, aes(x = etiology_class, y = pathway_f, fill = NES)) +
        geom_tile(color = "white", linewidth = 0.3) +
        scale_fill_gradient2(
          low      = masld_colors$down,
          mid      = "white",
          high     = masld_colors$up,
          midpoint = 0,
          limits   = c(-nes_lim, nes_lim),
          name     = "NES"
        ) +
        labs(
          title = "Pathway enrichment by etiology class (Hallmark, top 15 by mean |NES|)",
          x     = "Etiology class",
          y     = NULL
        ) +
        theme_masld(base_size = 9) +
        theme(
          axis.text.y     = element_text(size = 7),
          legend.position = "right"
        )
    }
  }
}

save_panel(p_G, "G", width = fig_half_width + 1, height = 4)

# ---------------------------------------------------------------------------
# Composite figure: patchwork layout
# ---------------------------------------------------------------------------
cat("\nAssembling composite figure...\n")

# Two rows:
#   Row 1: A | B | C
#   Row 2: D | E | F | G  (D spans 2 cols conceptually, but patchwork handles widths)

composite <- (p_A | p_B | p_C) /
             (p_D | p_E | p_F | p_G) +
  plot_annotation(
    tag_levels = "A",
    title      = "Plasma proteomics interpretation",
    theme      = theme_masld(base_size = 9) +
                 theme(plot.title = element_text(size = 10, face = "bold"))
  )

composite_path <- file.path(OUTDIR, "226e_plasma_interpretation.pdf")
pdf_dev <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
ggsave(
  composite_path,
  composite,
  width  = fig_full_width * 1.8,
  height = 9,
  device = pdf_dev
)
cat("  Saved composite:", composite_path, "\n")

cat("\n=== 226e COMPLETE ===\n")
cat("Finished:", as.character(Sys.time()), "\n")
