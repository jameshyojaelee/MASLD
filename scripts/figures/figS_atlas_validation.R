#!/usr/bin/env Rscript
# ==========================================================================
# Supplementary Figure: Atlas validation & cross-study reproducibility
# 3-panel layout:
#   a: Evidence source architecture — 7-source coverage bar chart
#   b: LOO-CV stability — per-cohort recovery + gene robustness histogram
#   c: Cross-study AUROC heatmap (single-cohort + integrated)
# ==========================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT <- file.path(FIGS08_DIR, "figS_atlas_validation.pdf")
dir.create(file.path(FIGS08_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

# LOO-CV directory
LOO_DIR <- file.path(INT_RESULTS, "loo_cv")

# Pre-initialize placeholders
p_a <- placeholder("Panel a: multi_evidence_atlas.csv not found")
p_b <- placeholder("Panel b: LOO-CV data not found")
p_c <- placeholder("Panel c: AUROC data not found")

# ==========================================================================
# Panel a: Evidence source coverage bar chart
# Compute active coverage per source from the atlas columns
# ==========================================================================
atlas <- load_multi_evidence()
if (!is.null(atlas)) {

  n_genes <- nrow(atlas)

  # Recompute source-active flags using same logic as 45a_integrate_all_sources.R
  atlas[, src1_human_bulk := is_dream_deg(atlas)]
  atlas[, src2_mouse_bulk := !is.na(mouse_meta_padj) & mouse_meta_padj < 0.1]

  # S3: Genetic causal (mr_pval removed 2026-04-22 — MR ditched from paper)
  atlas[, src3_genetic := FALSE]
  atlas[!is.na(twas_pval) & twas_pval < 0.05,
        src3_genetic := TRUE]
  coloc_cols_s3 <- intersect(
    c("coloc_pp4", "sceqtl_coloc_best_pp4", "broadaway_coloc_pp4",
      "ukbb_alt_coloc_pp4", "ast_coloc_pp4", "ggt_coloc_pp4", "pdff_coloc_pp4"),
    names(atlas))
  for (cc in coloc_cols_s3) {
    atlas[!is.na(get(cc)) & get(cc) > 0.5, src3_genetic := TRUE]
  }
  if ("zenodo_nafld_coloc" %in% names(atlas))
    atlas[zenodo_nafld_coloc == TRUE, src3_genetic := TRUE]
  if ("ieqtl_disease_interaction" %in% names(atlas))
    atlas[ieqtl_disease_interaction == TRUE, src3_genetic := TRUE]
  if ("sceqtl_twas_best_fdr" %in% names(atlas))
    atlas[!is.na(sceqtl_twas_best_fdr) & sceqtl_twas_best_fdr < 0.05, src3_genetic := TRUE]

  # S4: Essentiality
  atlas[, src4_essentiality := !is.na(essentiality_chronos)]

  # S5: Epigenomic
  atlas[, src5_epigenomic := FALSE]
  if ("mouse_da_padj" %in% names(atlas))
    atlas[!is.na(mouse_da_padj) & mouse_da_padj < 0.1, src5_epigenomic := TRUE]
  if ("hepatocyte_da_padj" %in% names(atlas))
    atlas[!is.na(hepatocyte_da_padj) & hepatocyte_da_padj < 0.1, src5_epigenomic := TRUE]
  if ("scenic_grn_target" %in% names(atlas))
    atlas[scenic_grn_target == TRUE, src5_epigenomic := TRUE]
  if ("cross_species_promoter_conserved" %in% names(atlas))
    atlas[cross_species_promoter_conserved == TRUE, src5_epigenomic := TRUE]

  # S6: Spatial
  atlas[, src6_spatial := FALSE]
  if ("spatial_is_svg" %in% names(atlas))
    atlas[!is.na(spatial_is_svg) & spatial_is_svg == TRUE, src6_spatial := TRUE]

  # S7: Single-cell
  atlas[, src7_singlecell := FALSE]
  if ("sc_best_padj" %in% names(atlas))
    atlas[!is.na(sc_best_padj) & sc_best_padj < 0.05, src7_singlecell := TRUE]
  if ("liana_n_diff_interactions" %in% names(atlas))
    atlas[!is.na(liana_n_diff_interactions) & liana_n_diff_interactions > 0,
          src7_singlecell := TRUE]

  # Build summary table
  src_dt <- data.table(
    source = c("S1: Human bulk", "S2: Mouse bulk", "S3: Genetic causal",
               "S4: Essentiality", "S5: Epigenomic", "S6: Spatial",
               "S7: Single-cell"),
    n_active = c(sum(atlas$src1_human_bulk, na.rm = TRUE),
                 sum(atlas$src2_mouse_bulk, na.rm = TRUE),
                 sum(atlas$src3_genetic, na.rm = TRUE),
                 sum(atlas$src4_essentiality, na.rm = TRUE),
                 sum(atlas$src5_epigenomic, na.rm = TRUE),
                 sum(atlas$src6_spatial, na.rm = TRUE),
                 sum(atlas$src7_singlecell, na.rm = TRUE))
  )
  src_dt[, pct := n_active / n_genes * 100]
  src_dt[, source := factor(source, levels = rev(source))]

  # Color gradient: deeper magenta for higher coverage
  source_fills <- c(
    "S1: Human bulk"     = masld_colors$deg,
    "S2: Mouse bulk"     = masld_colors$down,
    "S3: Genetic causal" = masld_colors$mr,
    "S4: Essentiality"   = masld_colors$ns,
    "S5: Epigenomic"     = masld_colors$twas,
    "S6: Spatial"        = "#F48FB1",
    "S7: Single-cell"    = masld_colors$gwas
  )

  p_a <- ggplot(src_dt, aes(x = pct, y = source, fill = source)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = paste0(sprintf("%.1f", pct), "%  (n=", comma(n_active), ")")),
              hjust = -0.05, size = GEOM_TEXT_6PT, color = "gray20") +
    scale_fill_manual(values = source_fills, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.35)),
                       labels = function(x) paste0(x, "%")) +
    labs(x = paste0("Gene coverage (N = ", comma(n_genes), " genes)"),
         y = NULL) +
    theme_masld()
}

# ==========================================================================
# Panel b: LOO-CV stability
# Left: bar chart of per-cohort recovery rates
# Right: gene robustness histogram (inset-style, combined via patchwork)
# ==========================================================================
loo_summary_path <- file.path(LOO_DIR, "loo_cv_summary.csv")
loo_gene_path    <- file.path(LOO_DIR, "loo_cv_per_gene.csv")

if (file.exists(loo_summary_path) && file.exists(loo_gene_path)) {
  loo_summary <- fread(loo_summary_path)
  loo_genes   <- fread(loo_gene_path)

  # --- Left: per-cohort recovery bar chart ---
  mean_recovery <- mean(loo_summary$pct_full_recovered)

  loo_summary[, cohort := factor(held_out, levels = held_out[order(pct_full_recovered)])]

  p_b_bars <- ggplot(loo_summary, aes(x = pct_full_recovered, y = cohort)) +
    geom_col(fill = masld_colors$deg, width = 0.6) +
    geom_vline(xintercept = mean_recovery, linetype = "dashed",
               linewidth = 0.3, color = "gray40") +
    geom_text(aes(label = paste0(sprintf("%.1f", pct_full_recovered), "%")),
              hjust = -0.1, size = GEOM_TEXT_6PT, color = "gray20") +
    annotate("text", x = mean_recovery, y = 0.5,
             label = paste0("Mean: ", sprintf("%.1f", mean_recovery), "%"),
             hjust = -0.05, vjust = -0.5, size = GEOM_TEXT_6PT, color = "gray40") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2)),
                       labels = function(x) paste0(x, "%")) +
    labs(x = "DEG recovery (%)",
         y = NULL) +
    theme_masld()

  # --- Right: gene robustness histogram ---
  # Categorize: Robust (8/8), Stable (6-7/8), Fragile (0-5/8)
  # Only count genes that are significant in the full model
  loo_genes_sig <- loo_genes[full_sig == TRUE | full_sig == "TRUE"]
  if (nrow(loo_genes_sig) == 0) {
    # Fallback: use robustness column directly
    loo_genes_sig <- loo_genes
  }

  rob_summary <- loo_genes_sig[, .N, by = robustness]

  # Order: Robust > Stable > Fragile
  rob_levels <- c("Robust", "Stable", "Fragile")
  rob_summary <- rob_summary[robustness %in% rob_levels]
  rob_summary[, robustness := factor(robustness, levels = rob_levels)]

  rob_colors <- c(
    Robust  = masld_colors$conserved,
    Stable  = masld_colors$deg,
    Fragile = masld_colors$ns
  )

  p_b_hist <- ggplot(rob_summary, aes(x = robustness, y = N, fill = robustness)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = comma(N)), vjust = -0.3, size = GEOM_TEXT_6PT) +
    scale_fill_manual(values = rob_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                       labels = comma) +
    labs(x = NULL, y = "DEGs") +
    theme_masld()

  # Combine bar + histogram side by side
  p_b <- p_b_bars + p_b_hist +
    plot_layout(widths = c(3, 2))

} else if (file.exists(loo_summary_path)) {
  # Only summary available (no per-gene)
  loo_summary <- fread(loo_summary_path)
  mean_recovery <- mean(loo_summary$pct_full_recovered)
  loo_summary[, cohort := factor(held_out, levels = held_out[order(pct_full_recovered)])]

  p_b <- ggplot(loo_summary, aes(x = pct_full_recovered, y = cohort)) +
    geom_col(fill = masld_colors$deg, width = 0.6) +
    geom_vline(xintercept = mean_recovery, linetype = "dashed",
               linewidth = 0.3, color = "gray40") +
    geom_text(aes(label = paste0(sprintf("%.1f", pct_full_recovered), "%")),
              hjust = -0.1, size = GEOM_TEXT_6PT, color = "gray20") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2)),
                       labels = function(x) paste0(x, "%")) +
    labs(x = "DEG recovery (%)", y = NULL) +
    theme_masld()
}

# ==========================================================================
# Panel c: Cross-study AUROC heatmap
# Rows = training cohort (+ Integrated), Columns = test cohort
# ==========================================================================
auroc_mat_path <- file.path(LOO_DIR, "auroc_matrix.csv")
auroc_int_path <- file.path(LOO_DIR, "auroc_integrated.csv")

if (file.exists(auroc_mat_path)) {
  auroc_mat <- fread(auroc_mat_path)
  auroc_int <- if (file.exists(auroc_int_path)) fread(auroc_int_path) else NULL

  # Exclude self-prediction (train == test) — set to NA for cleaner display
  auroc_mat[train == test, auroc := NA_real_]

  # Add integrated row if available
  if (!is.null(auroc_int)) {
    int_rows <- data.table(
      train = "Integrated",
      test = auroc_int$held_out,
      auroc = auroc_int$auroc,
      n_test = auroc_int$n_test,
      n_disease = auroc_int$n_disease,
      n_control = auroc_int$n_control
    )
    # Fill missing columns
    if ("train_is_fibrosis_contrast" %in% names(auroc_mat))
      int_rows[, train_is_fibrosis_contrast := FALSE]
    auroc_combined <- rbind(auroc_mat, int_rows, fill = TRUE)
  } else {
    auroc_combined <- auroc_mat
  }

  # Order: cohorts alphabetically, Integrated at bottom
  train_levels <- sort(unique(auroc_combined$train[auroc_combined$train != "Integrated"]))
  if ("Integrated" %in% auroc_combined$train) {
    train_levels <- c(train_levels, "Integrated")
  }
  test_levels <- sort(unique(auroc_combined$test))

  auroc_combined[, train := factor(train, levels = rev(train_levels))]
  auroc_combined[, test := factor(test, levels = test_levels)]

  p_c <- ggplot(auroc_combined, aes(x = test, y = train, fill = auroc)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(is.na(auroc), "",
                                 sprintf("%.2f", auroc))),
              size = GEOM_TEXT_6PT, color = ifelse(
                is.na(auroc_combined$auroc), "white",
                ifelse(auroc_combined$auroc > 0.75, "white", "black")
              )) +
    scale_fill_gradient2(low = masld_colors$down,
                         mid = "white",
                         high = masld_colors$up,
                         midpoint = 0.7,
                         na.value = "gray95",
                         limits = c(0.4, 1),
                         oob = squish,
                         name = "AUROC") +
    labs(x = "Test cohort",
         y = "Training cohort") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          panel.grid = element_blank(),
          legend.key.height = unit(0.5, "cm"),
          legend.key.width = unit(0.25, "cm"))

  # Add separator line above Integrated row
  if ("Integrated" %in% train_levels) {
    p_c <- p_c +
      geom_hline(yintercept = 1.5, linewidth = 0.6, color = "gray30")
  }
}

# ==========================================================================
# Assemble: 2 rows — top: panel a (full width), bottom: panels b + c
# ==========================================================================
fig_atlas <- p_a / (p_b | p_c) +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(
    tag_levels = "a"
  ) &
  theme(plot.tag = element_text(size = 8, face = "plain"))

save_fig_tall(fig_atlas, OUT, height = 8)
message("[caption] Atlas architecture and cross-study validation")
message("Atlas validation supplementary figure saved to ", OUT)
