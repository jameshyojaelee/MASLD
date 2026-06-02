#!/usr/bin/env Rscript
# figS16_plasma_proteomics.R
# Supplementary Figure 16: Plasma proteomics biological characterization
# Panels:
#   A. CHI3L1 dominance lollipop (T1_binary top-20 SHAP)
#   B. Panel size curves — fibrosis (focal) vs etiology (diffuse)
#   C. Etiology-discriminating proteins: metabolic vs immune biology
#   D. Yang et al. forest plot — nested CV reveals optimism bias
# Output: figures/supplementary/figS10_prediction/figS16_plasma_proteomics.pdf (14 x 10 inches)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
})

set.seed(42)

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUTDIR  <- FIGS10_DIR
dir.create(file.path(OUTDIR, "panels"), showWarnings = FALSE, recursive = TRUE)
DATADIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram")

message("figS16 — Plasma proteomics: loading data...")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
load_csv <- function(fname) {
  path <- file.path(DATADIR, fname)
  if (!file.exists(path)) {
    warning("File not found, panel will be skipped: ", path)
    return(NULL)
  }
  fread(path)
}

target_facet_labels <- c(
  T1_binary        = "Fibrosis staging\n(F0-2 vs F3-4)",
  T5_etiology_binary = "Etiology discrimination\n(MASLD vs CVH)"
)

metric_axis_labels <- c(
  auroc    = "AUROC",
  macro_f1 = "Macro F1"
)

# ============================================================
# Panel A — CHI3L1 dominance lollipop (T1_binary, top-20)
# ============================================================
message("  Panel A: CHI3L1 dominance lollipop")

importance_raw <- load_csv("226a_protein_importance.csv")

p_a <- if (is.null(importance_raw)) {
  placeholder("Panel A: 226a_protein_importance.csv not found")
} else {
  t1_top20 <- importance_raw[
    target == "T1_binary"
  ][order(rank)][seq_len(min(20L, .N))]

  # Compute CHI3L1 fraction
  chi3l1_shap <- t1_top20[protein == "CHI3L1", mean_abs_shap]
  total_shap  <- sum(t1_top20$mean_abs_shap)
  chi3l1_pct  <- if (length(chi3l1_shap) > 0 && total_shap > 0) {
    round(100 * chi3l1_shap / total_shap)
  } else {
    NA_integer_
  }

  t1_top20[, protein := factor(protein, levels = rev(protein))]
  t1_top20[, is_chi3l1 := protein == "CHI3L1"]

  # Lollipop: horizontal for readability (proteins on y, SHAP on x)
  annot_label <- if (!is.na(chi3l1_pct)) {
    paste0("CHI3L1 = ", chi3l1_pct, "% of top-20 signal")
  } else {
    "CHI3L1 (dominant)"
  }

  chi3l1_x <- if (length(chi3l1_shap) > 0) chi3l1_shap else max(t1_top20$mean_abs_shap)

  ggplot(t1_top20, aes(x = mean_abs_shap, y = protein, color = is_chi3l1)) +
    geom_segment(aes(x = 0, xend = mean_abs_shap,
                     y = protein, yend = protein),
                 linewidth = 0.5) +
    geom_point(size = 2.2) +
    annotate("text",
             x     = chi3l1_x + max(t1_top20$mean_abs_shap) * 0.02,
             y     = "CHI3L1",
             label = annot_label,
             hjust = 0, size = 2.2, color = "#C2185B", fontface = "italic") +
    scale_color_manual(values = c("FALSE" = "#BDBDBD", "TRUE" = "#C2185B"),
                       guide  = "none") +
    scale_x_continuous(labels = scientific_format(digits = 2),
                       expand = expansion(mult = c(0, 0.35))) +
    labs(
      title = "A  Fibrosis signal dominated by a single protein",
      x     = "Mean |SHAP|",
      y     = "Plasma protein"
    ) +
    theme_masld(base_size = 7)
}

save_fig(p_a,
         file.path(OUTDIR, "panels", "figS16_panelA_chi3l1_dominance.pdf"),
         width = fig_half_width, height = 3.5)
message("    Panel A done.")

# ============================================================
# Panel B — Panel size curves (focal vs diffuse)
# ============================================================
message("  Panel B: Panel size curves")

curves_raw <- load_csv("226a_panel_curves.csv")

p_b <- if (is.null(curves_raw)) {
  placeholder("Panel B: 226a_panel_curves.csv not found")
} else {
  # Filter to T1_binary and T5_etiology_binary; compute metric per target
  # T1 uses auroc; T5 uses macro_f1 if available, else auroc
  t1_metric <- "auroc"
  t5_metric <- if ("macro_f1" %in% unique(curves_raw$metric_name)) "macro_f1" else "auroc"

  curves_filt <- curves_raw[
    target %in% c("T1_binary", "T5_etiology_binary") &
    ((target == "T1_binary"           & metric_name == t1_metric) |
     (target == "T5_etiology_binary"  & metric_name == t5_metric))
  ]

  # Summarise: mean + 2.5/97.5 percentile per target × panel_size × is_random
  curves_summ <- curves_filt[, .(
    mean_val   = mean(metric_value, na.rm = TRUE),
    lo         = quantile(metric_value, 0.025, na.rm = TRUE),
    hi         = quantile(metric_value, 0.975, na.rm = TRUE)
  ), by = .(target, panel_size, is_random)]

  # Labels
  curves_summ[, facet_label := target_facet_labels[target]]
  curves_summ[, line_type   := ifelse(is_random, "Random", "SHAP-ranked")]
  curves_summ[, y_label     := ifelse(target == "T1_binary",
                                       metric_axis_labels["auroc"],
                                       metric_axis_labels[t5_metric])]

  facet_order <- c("Fibrosis staging\n(F0-2 vs F3-4)",
                   "Etiology discrimination\n(MASLD vs CVH)")
  curves_summ[, facet_label := factor(facet_label, levels = facet_order)]
  curves_summ[, line_type   := factor(line_type,   levels = c("SHAP-ranked", "Random"))]

  panel_sizes_all <- sort(unique(curves_summ$panel_size))
  break_vals <- panel_sizes_all[panel_sizes_all %in% c(2, 5, 10, 20, 50, 100, 200, 1461)]

  ggplot(curves_summ,
         aes(x     = panel_size,
             y     = mean_val,
             color = line_type,
             fill  = line_type,
             group = line_type)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.2) +
    # Annotation: Random outperforms curated at n<200 in T5 facet
    annotate("text",
             x = 5, y = 0.45,
             label = "Random \u2265 SHAP\nat n < 200",
             size = 2.0, color = "#7B1FA2", hjust = 0, fontface = "italic") +
    facet_wrap(~ facet_label, scales = "free_y", nrow = 1) +
    scale_color_manual(values = c("SHAP-ranked" = "#1565C0", "Random" = "#9E9E9E"),
                       name = NULL) +
    scale_fill_manual(values  = c("SHAP-ranked" = "#1565C0", "Random" = "#9E9E9E"),
                      name = NULL) +
    scale_x_log10(breaks = break_vals,
                  labels = ifelse(break_vals == 1461, "All\n(1461)", as.character(break_vals))) +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    labs(
      title = "B  Fibrosis = focal; etiology = diffuse proteome encoding",
      x     = "Panel size (proteins)",
      y     = "Performance"
    ) +
    theme_masld(base_size = 7) +
    theme(legend.position = "bottom",
          axis.text.x     = element_text(size = 6))
}

save_fig(p_b,
         file.path(OUTDIR, "panels", "figS16_panelB_panel_curves.pdf"),
         width = fig_full_width, height = 3.2)
message("    Panel B done.")

# ============================================================
# Panel C — Etiology-discriminating proteins
# ============================================================
message("  Panel C: Etiology-discriminating proteins")

p_c <- if (is.null(importance_raw)) {
  placeholder("Panel C: 226a_protein_importance.csv not found")
} else {
  t5_top10 <- importance_raw[
    target == "T5_etiology_binary"
  ][order(rank)][seq_len(min(10L, .N))]

  # Classify proteins by known biology
  masld_metabolic_proteins <- c("NPY", "REN", "GCG", "PRKCQ", "LEP",
                                  "CPA2", "CES1", "GLP1R", "INSR",
                                  "ADIPOQ", "RETN", "NAMPT")
  immune_proteins           <- c("TPSAB1", "TNFRSF9", "MME", "VIM",
                                  "SFTPA2", "FCRL3", "HSP90B1", "TFF1",
                                  "PON2", "C1QA", "C1QB", "C3", "CFB",
                                  "IL6", "TNF", "CRP")

  t5_top10[, biology_class := fcase(
    protein %in% masld_metabolic_proteins, "MASLD-metabolic",
    protein %in% immune_proteins,          "CVH-immune",
    default = "Other"
  )]
  t5_top10[, protein := factor(protein,
                                levels = t5_top10[order(-rank), protein])]

  bio_colors <- c(
    "MASLD-metabolic" = "#F57F17",  # Warm amber
    "CVH-immune"      = "#1565C0",  # Deep blue
    "Other"           = "#9E9E9E"   # Gray
  )

  ggplot(t5_top10, aes(x = mean_abs_shap, y = protein, fill = biology_class)) +
    geom_col(width = 0.65) +
    scale_fill_manual(values = bio_colors, name = "Biology") +
    scale_x_continuous(labels = scientific_format(digits = 2),
                       expand = expansion(mult = c(0, 0.08))) +
    labs(
      title = "C  Etiology proteins reflect metabolic vs immune biology",
      x     = "Mean |SHAP|",
      y     = "Plasma protein"
    ) +
    theme_masld(base_size = 7) +
    theme(legend.position = "right")
}

save_fig(p_c,
         file.path(OUTDIR, "panels", "figS16_panelC_etiology_proteins.pdf"),
         width = fig_half_width, height = 3.2)
message("    Panel C done.")

# ============================================================
# Panel D — Yang et al. forest plot
# ============================================================
message("  Panel D: Yang et al. comparison forest plot")

forest_df <- data.table(
  model    = c(
    "Full Olink (1,461 proteins)",
    "NFASC+GDF15",
    "NFASC+GDF15",
    "FIB-4 (reference)"
  ),
  auroc    = c(0.790, 0.734, 0.870, 0.760),
  ci_low   = c(0.710, 0.664,    NA,    NA),
  ci_high  = c(0.870, 0.804,    NA,    NA),
  study    = c(
    "This study (nested CV)",
    "This study (nested CV)",
    "Yang et al. 2025 (single split)",
    "Literature"
  ),
  cv_type  = c(
    "10\u00d75-fold nested CV",
    "10\u00d75-fold nested CV",
    "Single 70/30 split",
    "Multiple cohorts"
  )
)

# Order: reference at bottom, Yang next, then ours
forest_df[, model_label := paste0(model, "\n[", cv_type, "]")]
forest_df[, row_y := .N:1]

study_colors <- c(
  "This study (nested CV)"         = "#1565C0",
  "Yang et al. 2025 (single split)" = "#C2185B",
  "Literature"                      = "#757575"
)

# Gap annotation between Yang 0.870 and our 0.734
gap_row_yang <- forest_df[study == "Yang et al. 2025 (single split)", row_y]
gap_row_ours <- forest_df[study == "This study (nested CV)" & grepl("NFASC", model), row_y]
mid_y <- mean(c(gap_row_yang, gap_row_ours))

p_d <- ggplot(forest_df,
              aes(x = auroc, y = row_y, color = study, shape = study)) +
  # Error bars (only where CI available)
  geom_errorbarh(
    data = forest_df[!is.na(ci_low)],
    aes(xmin = ci_low, xmax = ci_high),
    height = 0.25, linewidth = 0.5
  ) +
  geom_point(size = 2.8) +
  # Reference line at 0.5
  geom_vline(xintercept = 0.5, linetype = "dashed",
             color = "#BDBDBD", linewidth = 0.4) +
  # Optimism bias annotation: bracket-style segment between Yang and ours
  annotate("segment",
           x = 0.870, xend = 0.734,
           y = mid_y + 0.05, yend = mid_y + 0.05,
           arrow = arrow(ends = "both", length = unit(0.08, "cm"), type = "open"),
           color = "#F57F17", linewidth = 0.6) +
  annotate("text",
           x = (0.870 + 0.734) / 2,
           y = mid_y + 0.18,
           label = "~0.14 optimism bias",
           size = 2.1, color = "#F57F17", hjust = 0.5, fontface = "italic") +
  scale_color_manual(values = study_colors, name = NULL) +
  scale_shape_manual(
    values = c(
      "This study (nested CV)"          = 16,
      "Yang et al. 2025 (single split)" = 17,
      "Literature"                      = 15
    ), name = NULL) +
  scale_x_continuous(limits = c(0.55, 1.0),
                     labels = number_format(accuracy = 0.01)) +
  scale_y_continuous(
    breaks = forest_df$row_y,
    labels = forest_df$model_label
  ) +
  labs(
    title = "D  Nested CV reveals optimism in single-split evaluation",
    x     = "AUROC (fibrosis F0-2 vs F3-4)",
    y     = NULL
  ) +
  theme_masld(base_size = 7) +
  theme(
    legend.position = "bottom",
    legend.box      = "vertical",
    axis.text.y     = element_text(size = 6, hjust = 1)
  )

save_fig(p_d,
         file.path(OUTDIR, "panels", "figS16_panelD_yang_comparison.pdf"),
         width = fig_half_width, height = 3.2)
message("    Panel D done.")

# ============================================================
# Assembly — 2 × 2 with patchwork
# ============================================================
message("  Assembling final figure...")

final_plot <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(
    title = paste0(
      "Supplementary Figure 16. ",
      "Plasma proteomics biological characterization."
    ),
    theme = theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0)
    )
  )

save_fig(
  final_plot,
  file.path(OUTDIR, "figS16_plasma_proteomics.pdf"),
  width  = 14,
  height = 10
)

message("figS16 complete. Output: ", file.path(OUTDIR, "figS16_plasma_proteomics.pdf"))
