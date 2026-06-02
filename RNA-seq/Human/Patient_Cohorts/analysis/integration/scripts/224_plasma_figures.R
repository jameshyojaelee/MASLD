#!/usr/bin/env Rscript
# QUARANTINED 2026-06-01: plots supervised Olink plasma-classifier results that are WITHDRAWN (unrecoverable subject-label join). Do not regenerate for the manuscript. See docs/reviews/2026-06-01-proteomics-extreme-review.md.
# 224_plasma_figures.R
# Supplementary Figure S10 panels (c)-(f): Plasma Proteomics
#
# Generates four panels and assembles them into a composite PDF:
#   (c) Tissue-to-plasma funnel: 5,484 DEGs -> 812 plasma-detectable -> N dual-sig
#   (d) Plasma panel AUROC with random baseline violin overlay
#   (e) F2 switch scores by fibrosis stage (box/violin + KW + pairwise p-values)
#   (f) Cross-platform protein overlap (Olink vs DIA-MS Venn-style bar chart)
#
# Inputs (all in results/multiprogram/ unless noted):
#   tissue_plasma_bridge.csv      (staging_classifier results)
#   plasma_panel_results.csv
#   plasma_panel_summary.csv
#   plasma_random_baselines.csv
#   plasma_switch_scores.csv
#   plasma_switch_statistics.csv
#   cross_platform_validation_summary.csv  (from Script 223)
#
# Outputs -> figures/supplementary/prediction/panels/
#   panel_c_tissue_plasma_funnel.pdf
#   panel_d_plasma_auroc_random.pdf
#   panel_e_switch_scores_stage.pdf
#   panel_f_platform_overlap.pdf
#   fig_plasma_main.pdf   (composite 4-panel)
#
# SLURM: --partition=cpu --cpus-per-task=4 --mem=16G --time=48:00:00
# Env:   micromamba activate rnaseq

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
MP_DIR <- file.path(INT,  "results/multiprogram")
SC_DIR <- file.path(INT,  "results/staging_classifier")
FIGDIR <- file.path(BASE, "figures/supplementary/figS10_prediction")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(FIGDIR, "panels"), showWarnings = FALSE, recursive = TRUE)

source(file.path(BASE, "scripts/figures/publication_theme.R"))

cat("=== 224: Plasma Figure Panels (c)-(f) ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ── Helper: format p-value as label ───────────────────────────────────────────
fmt_pval <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return(sprintf("p=%.2g", p))
  if (p < 0.05)  return(sprintf("p=%.3f", p))
  return(sprintf("p=%.2f", p))
}

# ── PANEL (c): Tissue-to-plasma funnel ────────────────────────────────────────
cat("--- Panel (c): Tissue-to-plasma funnel ---\n")

bridge_file <- file.path(SC_DIR, "tissue_plasma_bridge.csv")
if (!file.exists(bridge_file)) {
  bridge_file <- file.path(MP_DIR, "plasma_selected_proteins.csv")
  cat("  Falling back to plasma_selected_proteins.csv for bridge data\n")
}

if (file.exists(bridge_file)) {
  bridge <- fread(bridge_file, na.strings = c("", "NA"))

  # Resolve column names for each source
  if ("human_symbol" %in% names(bridge)) {
    # tissue_plasma_bridge.csv from staging_classifier
    n_total_degs     <- sum(bridge$is_deg == TRUE, na.rm = TRUE)
    n_plasma_detect  <- sum(bridge$in_plasma == TRUE, na.rm = TRUE)
    # dual-significant: tissue DEG AND detected in plasma
    n_dual_sig       <- sum(bridge$is_deg == TRUE & bridge$in_plasma == TRUE, na.rm = TRUE)
    n_plasma_only    <- n_plasma_detect - n_dual_sig
    n_tissue_only    <- n_total_degs - n_dual_sig
  } else if ("is_tissue_deg" %in% names(bridge)) {
    # plasma_selected_proteins.csv
    n_total_degs    <- 5484L   # canonical dream DEG count
    n_plasma_detect <- sum(bridge$in_panel_a == "True" | bridge$in_panel_a == TRUE, na.rm = TRUE)
    n_dual_sig      <- sum((bridge$is_tissue_deg == TRUE | bridge$is_tissue_deg == "True") &
                             (bridge$in_panel_a == TRUE | bridge$in_panel_a == "True"),
                           na.rm = TRUE)
    n_plasma_only   <- n_plasma_detect - n_dual_sig
    n_tissue_only   <- n_total_degs - n_dual_sig
  } else {
    cat("  WARNING: bridge file has unexpected format; using defaults\n")
    n_total_degs    <- 5484L
    n_plasma_detect <- 812L
    n_dual_sig      <- 62L
    n_plasma_only   <- n_plasma_detect - n_dual_sig
    n_tissue_only   <- n_total_degs - n_dual_sig
  }
} else {
  cat("  WARNING: tissue_plasma_bridge.csv not found; using hardcoded counts\n")
  n_total_degs    <- 5484L
  n_plasma_detect <- 812L
  n_dual_sig      <- 62L
  n_plasma_only   <- n_plasma_detect - n_dual_sig
  n_tissue_only   <- n_total_degs - n_dual_sig
}

cat(sprintf("  Total DEGs: %d | Plasma-detectable: %d | Dual-sig: %d\n",
            n_total_degs, n_plasma_detect, n_dual_sig))

funnel_dt <- data.table(
  stage = factor(
    c("All DEGs\n(tissue)", "Plasma-\ndetectable", "Tissue+Plasma\nDEG"),
    levels = c("All DEGs\n(tissue)", "Plasma-\ndetectable", "Tissue+Plasma\nDEG")
  ),
  n      = c(n_total_degs, n_plasma_detect, n_dual_sig),
  fill   = c(masld_colors$down, masld_colors$masl, masld_colors$fibrosis)
)

panel_c <- ggplot(funnel_dt, aes(x = stage, y = n, fill = stage)) +
  geom_col(width = 0.6, color = NA) +
  geom_text(aes(label = scales::comma(n)), vjust = -0.4, size = 2.2,
            fontface = "bold") +
  scale_fill_manual(values = setNames(funnel_dt$fill, funnel_dt$stage),
                    guide = "none") +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "c  Tissue-to-plasma funnel",
    x     = NULL,
    y     = "Number of proteins"
  ) +
  theme_masld()

save_fig(panel_c,
         file.path(FIGDIR, "panels", "panel_c_tissue_plasma_funnel.pdf"),
         width  = 120 / 25.4,
         height = 100 / 25.4)
cat("  Saved: panel_c_tissue_plasma_funnel.pdf\n")

# ── PANEL (d): Plasma panel AUROC with random baseline ────────────────────────
cat("--- Panel (d): Plasma panel AUROC vs random baseline ---\n")

panel_results_file <- file.path(MP_DIR, "plasma_panel_results.csv")
random_file        <- file.path(MP_DIR, "plasma_random_baselines.csv")
summary_file       <- file.path(MP_DIR, "plasma_panel_summary.csv")

if (file.exists(panel_results_file) && file.exists(random_file)) {
  panel_res  <- fread(panel_results_file,  na.strings = c("", "NA"))
  random_res <- fread(random_file,         na.strings = c("", "NA"))

  # Use LOEO pooled AUROC as primary metric; fall back to CV mean
  auroc_col <- if ("auroc_loeo_pooled" %in% names(panel_res)) "auroc_loeo_pooled" else "auroc_cv_mean"
  rand_col  <- if ("auroc_loeo_pooled" %in% names(random_res)) "auroc_loeo_pooled" else "auroc_cv_mean"

  # Named panels with display labels
  panel_labels <- c(
    A_full_olink           = "A: All Olink\n(n=1,461)",
    B_tissue_informed_all  = "B: Tissue-\ninformed",
    C_top20                = "C: Top-20\npanel"
  )

  panels_dt <- panel_res[panel %in% names(panel_labels),
                          .(panel, auroc = get(auroc_col))]
  panels_dt[, label := panel_labels[panel]]
  panels_dt[, label := factor(label, levels = panel_labels)]

  # Random baseline distribution (all draws pooled; n=20 proteins)
  rand_vec <- random_res[[rand_col]]
  rand_dt  <- data.table(
    auroc = rand_vec,
    label = factor("Random\n(n=20)", levels = c(panel_labels, "Random\n(n=20)"))
  )

  # Empirical p-value for tissue-informed vs random
  pval_summary_file <- file.path(MP_DIR, "plasma_panel_summary.csv")
  pval_label <- ""
  if (file.exists(pval_summary_file)) {
    psumm <- fread(pval_summary_file, na.strings = c("", "NA"))
    # Use LOEO row if available
    prow <- psumm[grepl("LOEO|loeo", metric, ignore.case = TRUE)]
    if (nrow(prow) == 0) prow <- psumm[1]
    if ("empirical_pval" %in% names(psumm) && nrow(prow) > 0) {
      pv <- prow$empirical_pval[1]
      pval_label <- fmt_pval(as.numeric(pv))
    }
  }

  # Tissue-informed AUROC point for annotation
  ti_auroc <- panels_dt[panel == "B_tissue_informed_all", auroc]

  panel_d <- ggplot() +
    # Random baseline violin
    geom_violin(data = rand_dt,
                aes(x = label, y = auroc),
                fill  = masld_colors$ns, color = NA, alpha = 0.6,
                width = 0.7) +
    geom_boxplot(data = rand_dt,
                 aes(x = label, y = auroc),
                 width = 0.2, outlier.shape = NA,
                 fill = "white", color = "gray40", linewidth = 0.3) +
    # Named panels as points
    geom_point(data = panels_dt,
               aes(x = label, y = auroc),
               color = masld_colors$fibrosis, size = 2.5, shape = 18) +
    # P-value annotation for tissue-informed
    {
      if (nchar(pval_label) > 0 && length(ti_auroc) > 0) {
        annotate("text",
                 x     = "B: Tissue-\ninformed",
                 y     = ti_auroc + 0.04,
                 label = pval_label,
                 size  = 2, hjust = 0.5,
                 color = "gray30")
      } else NULL
    } +
    scale_y_continuous(limits = c(0.4, 1.0),
                       breaks = seq(0.4, 1.0, 0.1),
                       labels = scales::number_format(accuracy = 0.1)) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = "gray60", linewidth = 0.3) +
    labs(
      title = "d  Plasma panel AUROC vs random",
      x     = NULL,
      y     = "AUROC"
    ) +
    theme_masld()
} else {
  cat("  WARNING: panel results or random baselines not found; using placeholder\n")
  panel_d <- placeholder("d  AUROC panel\n(data pending)")
}

save_fig(panel_d,
         file.path(FIGDIR, "panels", "panel_d_plasma_auroc_random.pdf"),
         width  = 120 / 25.4,
         height = 100 / 25.4)
cat("  Saved: panel_d_plasma_auroc_random.pdf\n")

# ── PANEL (e): F2 switch scores by fibrosis stage ─────────────────────────────
cat("--- Panel (e): F2 switch scores by fibrosis stage ---\n")

switch_scores_file <- file.path(MP_DIR, "plasma_switch_scores.csv")
switch_stats_file  <- file.path(MP_DIR, "plasma_switch_statistics.csv")

if (file.exists(switch_scores_file)) {
  sw_scores <- fread(switch_scores_file, na.strings = c("", "NA"))
  cat(sprintf("  Loaded switch scores: %d subjects\n", nrow(sw_scores)))
  cat(sprintf("  Columns: %s\n", paste(names(sw_scores), collapse = ", ")))

  # Classify fibrosis_stage into 3 display groups
  sw_scores[, stage_group := fcase(
    fibrosis_stage %in% c("F0", "F1", "F2", "F0-2"), "F0-F2",
    fibrosis_stage %in% c("F3"),                      "F3",
    fibrosis_stage %in% c("F4"),                      "F4",
    default = NA_character_
  )]
  sw_scores <- sw_scores[!is.na(stage_group)]
  sw_scores[, stage_group := factor(stage_group, levels = c("F0-F2", "F3", "F4"))]

  cat(sprintf("  Stage groups: %s\n",
              paste(sw_scores[, .N, by = stage_group][, paste0(stage_group, "=", N)],
                    collapse = ", ")))

  # Pull KW p-value for switch_score from statistics file
  kw_pval  <- NA_real_
  kw_label <- ""
  if (file.exists(switch_stats_file)) {
    sw_stats <- fread(switch_stats_file, na.strings = c("", "NA"))
    kw_row   <- sw_stats[score == "switch_score" & test == "Kruskal-Wallis"]
    if (nrow(kw_row) > 0) {
      kw_pval  <- as.numeric(kw_row$padj[1])
      kw_label <- paste0("KW ", fmt_pval(kw_pval))
    }
    # Pairwise F0-2 vs F4
    pw_row   <- sw_stats[score == "switch_score" &
                           test == "Wilcoxon rank-sum" &
                           comparison == "F0-2 vs F4"]
    pw_label <- if (nrow(pw_row) > 0) {
      paste0("F0-2 vs F4 ", fmt_pval(as.numeric(pw_row$padj[1])))
    } else ""
  }

  stage_colors_3 <- c(
    "F0-F2" = masld_colors$control,
    "F3"    = masld_colors$masl,
    "F4"    = masld_colors$fibrosis
  )

  panel_e <- ggplot(sw_scores, aes(x = stage_group, y = switch_score,
                                    fill = stage_group)) +
    geom_violin(alpha = 0.4, color = NA, trim = FALSE, width = 0.8) +
    geom_boxplot(width = 0.25, outlier.shape = 21, outlier.size = 0.8,
                 outlier.alpha = 0.5, linewidth = 0.4, fill = "white") +
    scale_fill_manual(values = stage_colors_3, guide = "none") +
    labs(
      title = "e  Switch score by fibrosis stage",
      x     = "Fibrosis stage",
      y     = "Plasma switch score (NPX)"
    ) +
    {
      if (nchar(kw_label) > 0) {
        ymax <- max(sw_scores$switch_score, na.rm = TRUE)
        annotate("text", x = 2, y = ymax * 1.05,
                 label = kw_label, size = 2, hjust = 0.5, color = "gray30")
      } else NULL
    } +
    theme_masld()
} else {
  cat("  WARNING: plasma_switch_scores.csv not found; using placeholder\n")
  panel_e <- placeholder("e  Switch scores\n(data pending)")
}

save_fig(panel_e,
         file.path(FIGDIR, "panels", "panel_e_switch_scores_stage.pdf"),
         width  = 120 / 25.4,
         height = 100 / 25.4)
cat("  Saved: panel_e_switch_scores_stage.pdf\n")

# ── PANEL (e2): H5-inv random-null overlay (T0.12, added 2026-04-22) ──────────
# Histogram of 1,000 random 180-protein draws (KW H statistic), with the
# observed switch-panel H as a vertical line + empirical p annotation.
cat("--- Panel (e2): H5-inv random-null overlay ---\n")

null_rds_file <- file.path(MP_DIR, "plasma_switch_random_null_draws.rds")
null_csv_file <- file.path(MP_DIR, "plasma_switch_random_null.csv")

if (file.exists(null_rds_file) && file.exists(null_csv_file)) {
  null_obj   <- readRDS(null_rds_file)
  null_summ  <- fread(null_csv_file, na.strings = c("", "NA"))

  null_dt <- rbindlist(list(
    data.table(universe = "Full Olink\n(n=1,460)",    H = null_obj$null_full),
    data.table(universe = "Tissue-DEG\nsubset",        H = null_obj$null_deg),
    data.table(universe = "Non-DEG\nsubset",           H = null_obj$null_nondeg)
  ))
  null_dt[, universe := factor(universe,
           levels = c("Full Olink\n(n=1,460)", "Tissue-DEG\nsubset", "Non-DEG\nsubset"))]

  obs_H       <- null_obj$observed_kw
  emp_p_full  <- null_summ$emp_p_full_universe[1]
  emp_p_deg   <- null_summ$emp_p_deg_subset[1]
  emp_p_nondeg<- null_summ$emp_p_nondeg_subset[1]
  panel_size  <- null_obj$panel_size
  n_draws     <- null_obj$n_draws

  emp_label <- sprintf("obs H=%.1f\nemp_p (full)=%.3f\nemp_p (DEG)=%.3f\nemp_p (non-DEG)=%.3f\n(n=%d draws, k=%d)",
                       obs_H, emp_p_full, emp_p_deg, emp_p_nondeg,
                       n_draws, panel_size)

  panel_e_null <- ggplot(null_dt, aes(x = H, fill = universe)) +
    geom_histogram(bins = 40, alpha = 0.6, color = NA,
                   position = "identity") +
    geom_vline(xintercept = obs_H, linetype = "solid",
               color = masld_colors$fibrosis, linewidth = 0.6) +
    annotate("text",
             x = obs_H,
             y = Inf,
             label = emp_label,
             vjust = 1.2, hjust = 1.05, size = 1.8, color = "gray20") +
    scale_fill_manual(values = c(
      "Full Olink\n(n=1,460)" = masld_colors$ns,
      "Tissue-DEG\nsubset"    = masld_colors$masl,
      "Non-DEG\nsubset"       = masld_colors$control
    ), name = "Universe") +
    labs(
      title = "e2  H5-inv: size-matched random null",
      x     = "Kruskal-Wallis H (stage)",
      y     = "Random draws (count)"
    ) +
    theme_masld() +
    theme(legend.position = "right",
          legend.key.size = unit(3, "mm"))
} else {
  cat("  WARNING: null draws not found (expected after Script 222 rerun); placeholder\n")
  panel_e_null <- placeholder("e2  Random null\n(pending 222 rerun)")
}

save_fig(panel_e_null,
         file.path(FIGDIR, "panels", "panel_e2_switch_random_null.pdf"),
         width  = 140 / 25.4,
         height = 90 / 25.4)
cat("  Saved: panel_e2_switch_random_null.pdf\n")

# ── PANEL (f): Cross-platform protein overlap ─────────────────────────────────
cat("--- Panel (f): Cross-platform protein overlap ---\n")

xp_summary_file <- file.path(MP_DIR, "cross_platform_validation_summary.csv")

if (file.exists(xp_summary_file)) {
  xp_summ <- fread(xp_summary_file, na.strings = c("", "NA"))
  cat(sprintf("  Loaded cross-platform summary: %d rows\n", nrow(xp_summ)))

  get_val <- function(m) {
    v <- xp_summ[metric == m, value]
    if (length(v) == 0 || is.na(v[1])) NA_real_ else as.numeric(v[1])
  }

  n_olink  <- get_val("n_olink_proteins")
  n_ms     <- get_val("n_diams_proteins")
  n_ov     <- get_val("n_overlap_gene_symbol")
  rho_val  <- get_val("spearman_rho_mean_abundance")

  n_olink  <- if (is.na(n_olink))  1461 else n_olink
  n_ms     <- if (is.na(n_ms))      356 else n_ms
  n_ov     <- if (is.na(n_ov))       NA else n_ov
} else {
  cat("  WARNING: cross_platform_validation_summary.csv not found; using defaults\n")
  n_olink <- 1461; n_ms <- 356; n_ov <- NA_real_; rho_val <- NA_real_
}

# Build stacked bar representation of overlap (Venn in bar form)
n_ov_val    <- if (is.na(n_ov)) 0 else as.integer(n_ov)
olink_only  <- as.integer(n_olink) - n_ov_val
ms_only     <- as.integer(n_ms)    - n_ov_val

overlap_dt <- data.table(
  platform = c("Olink", "Olink", "DIA-MS", "DIA-MS"),
  category = c("Overlap", "Platform-only", "Overlap", "Platform-only"),
  n        = c(n_ov_val, olink_only, n_ov_val, ms_only)
)
overlap_dt[, category := factor(category, levels = c("Platform-only", "Overlap"))]

overlap_colors <- c(
  "Platform-only" = masld_colors$ns,
  "Overlap"       = masld_colors$fibrosis
)

rho_annot <- if (!is.na(rho_val)) {
  sprintf("Mean abundance\nSpearman rho=%.2f", rho_val)
} else ""

panel_f <- ggplot(overlap_dt, aes(x = platform, y = n, fill = category)) +
  geom_col(width = 0.5, color = "white", linewidth = 0.3) +
  geom_text(aes(label = scales::comma(n)),
            position = position_stack(vjust = 0.5),
            size = 2.2, color = "white", fontface = "bold") +
  scale_fill_manual(values = overlap_colors, name = "Proteins") +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.1))) +
  {
    if (nchar(rho_annot) > 0) {
      annotate("text",
               x = 1.5, y = max(n_olink, n_ms) * 1.02,
               label = rho_annot, size = 2, hjust = 0.5, color = "gray30")
    } else NULL
  } +
  labs(
    title = "f  Platform protein overlap",
    x     = "Platform",
    y     = "Number of proteins"
  ) +
  theme_masld() +
  theme(legend.position = "right")

save_fig(panel_f,
         file.path(FIGDIR, "panels", "panel_f_platform_overlap.pdf"),
         width  = 120 / 25.4,
         height = 100 / 25.4)
cat("  Saved: panel_f_platform_overlap.pdf\n")

# ── Composite: assemble all 4 panels ─────────────────────────────────────────
cat("--- Assembling composite fig_plasma_main.pdf ---\n")

composite <- (panel_c | panel_d) / (panel_e | panel_f) / (panel_e_null) +
  patchwork::plot_annotation(
    title = "Supplementary Figure S10 (c-f, e2): Plasma proteomics validation",
    theme = theme(
      plot.title = element_text(size = 7, face = "bold", hjust = 0,
                                family = "Helvetica")
    )
  )

save_fig(composite,
         file.path(FIGDIR, "fig_plasma_main.pdf"),
         width  = 180 / 25.4,
         height = 280 / 25.4)
cat("  Saved: fig_plasma_main.pdf\n")

cat("\nDone:", as.character(Sys.time()), "\n")
cat(sprintf("All outputs written to: %s\n", FIGDIR))
