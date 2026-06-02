# DEPRECATED — output superseded by renumbered figures
##############################################################################
# Supplementary Figure: Evidence Convergence — Source Orthogonality
# 5 panels: (a) unique source coverage, (b) MI heatmap, (c) cumulative
#   enrichment curves, (d) source addition bootstrap, (e) Bayesian waterfall
# Moved from main Figure 5 to supplementary (2026-03-18).
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT <- file.path(FIGS08_DIR, "figS_convergence.pdf")
dir.create(file.path(FIGS08_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Source label mapping (used across all panels)
# ---------------------------------------------------------------------------
SOURCE_LABELS <- c(
  S1 = "Human bulk\nRNA-seq",
  S2 = "Genetic\ncausal",
  S3 = "Essentiality\n(DepMap)",
  S4 = "Epigenomic\n(ATAC/SCENIC+)",
  S5 = "Spatial\n(Visium)",
  S6 = "Single-cell"
)

SOURCE_LABELS_SHORT <- c(
  S1 = "S1: Human",
  S2 = "S2: Genetic",
  S3 = "S3: Essential.",
  S4 = "S4: Epigenomic",
  S5 = "S5: Spatial",
  S6 = "S6: scRNA"
)

# Source colors (6 distinct colors from the masld palette)
source_colors <- c(
  S1 = "#1565C0",   # Deep blue (human bulk)
  S2 = "#C2185B",   # Magenta (genetic causal)
  S3 = "#42A5F5",   # Light blue (essentiality)
  S4 = "#7B1FA2",   # Violet (epigenomic)
  S5 = "#00695C",   # Teal (spatial)
  S6 = "#E91E63"    # Bright magenta (single-cell)
)

# Pre-initialize all 5 panels as placeholders
p_a <- placeholder("(a) Source-unique coverage")
p_b <- placeholder("(b) MI heatmap")
p_c <- placeholder("(c) Cumulative enrichment")
p_d <- placeholder("(d) Source addition")
p_e <- placeholder("(e) Bayesian waterfall")

# ==========================================================================
# (a) Source-by-source unique coverage
# ==========================================================================
atlas <- load_multi_evidence()
if (!is.null(atlas)) {
  # Define which columns indicate source activity for each source (6 sources)
  # S1: Human bulk RNA-seq
  s1_active <- is_dream_deg(atlas)
  # S2: Genetic causal (TWAS, COLOC significant)
  # (mr_sig and sceqtl_twas_best_fdr removed 2026-04-22 — MR ditched from paper)
  s2_active <- (!is.na(atlas$twas_pval) & atlas$twas_pval < 0.05) |
               (!is.na(atlas$coloc_pp4) & atlas$coloc_pp4 > 0.5) |
               (!is.na(atlas$broadaway_coloc_pp4) & atlas$broadaway_coloc_pp4 > 0.5)
  # S3: Essentiality
  s3_active <- !is.na(atlas$essentiality_chronos) & atlas$essentiality_chronos < -0.5
  # S4: Epigenomic
  s4_active <- (!is.na(atlas$mouse_da_padj) & atlas$mouse_da_padj < 0.05) |
               (!is.na(atlas$hepatocyte_da_padj) & atlas$hepatocyte_da_padj < 0.05) |
               (!is.na(atlas$scenic_grn_target) & atlas$scenic_grn_target != "" & !is.na(atlas$scenic_grn_target)) |
               (!is.na(atlas$cross_species_promoter_conserved) & atlas$cross_species_promoter_conserved == TRUE)
  # S5: Spatial
  s5_active <- !is.na(atlas$spatial_is_svg) & atlas$spatial_is_svg == TRUE
  # S6: Single-cell
  s6_active <- (!is.na(atlas$sc_hepatocyte_padj) & atlas$sc_hepatocyte_padj < 0.05) |
               (!is.na(atlas$liana_n_diff_interactions) & atlas$liana_n_diff_interactions > 0)

  # Build source activity matrix (6 sources, no mouse)
  src_mat <- data.table(
    gene = atlas$human_symbol,
    S1 = s1_active, S2 = s2_active, S3 = s3_active,
    S4 = s4_active, S5 = s5_active, S6 = s6_active
  )

  # Replace NAs with FALSE
  for (col in paste0("S", 1:6)) {
    src_mat[is.na(get(col)), (col) := FALSE]
  }

  # Count how many sources each gene has
  src_mat[, n_sources := rowSums(.SD), .SDcols = paste0("S", 1:6)]

  # For each source, count genes that are ONLY captured by that source
  unique_counts <- data.table(
    source = paste0("S", 1:6),
    unique_genes = sapply(paste0("S", 1:6), function(s) {
      other_cols <- setdiff(paste0("S", 1:6), s)
      sum(src_mat[[s]] == TRUE & rowSums(src_mat[, ..other_cols]) == 0)
    }),
    total_active = sapply(paste0("S", 1:6), function(s) {
      sum(src_mat[[s]] == TRUE)
    })
  )
  unique_counts[, pct_unique := unique_genes / total_active * 100]
  unique_counts[is.nan(pct_unique), pct_unique := 0]
  unique_counts[, label := SOURCE_LABELS_SHORT[source]]
  unique_counts[, source := factor(source, levels = paste0("S", 1:6))]

  p_a <- ggplot(unique_counts, aes(x = source, y = unique_genes, fill = source)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = paste0(unique_genes, "\n(", round(pct_unique, 0), "%)")),
              vjust = -0.3, size = 1.8, lineheight = 0.8) +
    scale_x_discrete(labels = SOURCE_LABELS_SHORT) +
    scale_fill_manual(values = source_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(x = NULL, y = "Uniquely captured genes",
         title = "Source-exclusive gene coverage") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 5.5))
}

# ==========================================================================
# (b) Mutual information heatmap (lower triangle)
# ==========================================================================
mi_file <- file.path(ME, "mutual_information.csv")
if (file.exists(mi_file)) {
  mi_dt <- fread(mi_file)

  # Build 6x6 matrix (off-diagonal only for display)
  sources <- paste0("S", 1:6)
  mi_mat <- matrix(NA, 6, 6, dimnames = list(sources, sources))
  for (i in seq_len(nrow(mi_dt))) {
    sa <- mi_dt$Source_A[i]
    sb <- mi_dt$Source_B[i]
    if (sa %in% sources && sb %in% sources) {
      mi_mat[sa, sb] <- mi_dt$MI[i]
    }
  }

  # Lower triangle only (exclude diagonal)
  mi_plot_dt <- data.table()
  for (i in 2:6) {
    for (j in 1:(i - 1)) {
      mi_plot_dt <- rbind(mi_plot_dt, data.table(
        Source_A = sources[i], Source_B = sources[j],
        MI = mi_mat[sources[i], sources[j]]
      ))
    }
  }

  mi_plot_dt[, Source_A := factor(Source_A, levels = sources)]
  mi_plot_dt[, Source_B := factor(Source_B, levels = sources)]

  # Format MI values for display
  mi_plot_dt[, mi_label := ifelse(MI < 0.01, sprintf("%.3f", MI), sprintf("%.2f", MI))]

  p_b <- ggplot(mi_plot_dt, aes(x = Source_B, y = Source_A, fill = MI)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = mi_label), size = 2, color = "black") +
    scale_fill_gradient2(
      low = "white", mid = "#F48FB1", high = "#880E4F",
      midpoint = 0.03, name = "MI (bits)",
      limits = c(0, max(mi_plot_dt$MI, na.rm = TRUE))
    ) +
    scale_x_discrete(labels = SOURCE_LABELS_SHORT[sources[1:5]]) +
    scale_y_discrete(labels = SOURCE_LABELS_SHORT[sources[2:6]]) +
    labs(title = "Mutual information between sources",
         x = NULL, y = NULL) +
    theme_masld() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
      axis.text.y = element_text(size = 5.5),
      legend.key.width = unit(0.3, "cm"),
      legend.key.height = unit(0.6, "cm"),
      panel.border = element_blank(),
      axis.line = element_blank()
    ) +
    coord_fixed()
}

# ==========================================================================
# (c) Cumulative enrichment curves (DGIdb druggable genes)
# ==========================================================================
enr_file <- file.path(ME, "cumulative_enrichment.csv")
if (file.exists(enr_file)) {
  enr_dt <- fread(enr_file)

  # Focus on V1_DGIdb (druggable genome, clearest signal)
  enr_sub <- enr_dt[validation == "V1_DGIdb"]

  # Select key rankings to display
  rankings_show <- c("R_multi", "R_S1", "R_S3", "R_count")
  enr_sub <- enr_sub[ranking %in% rankings_show]

  ranking_labels <- c(
    R_multi  = "Convergence score",
    R_S1     = "S1 only (human bulk)",
    R_S3     = "S3 only (genetic)",
    R_count  = "Source count"
  )
  ranking_colors <- c(
    R_multi  = "#880E4F",   # Dark magenta (best)
    R_S1     = "#1565C0",   # Blue
    R_S3     = "#C2185B",   # Magenta
    R_count  = "#7B1FA2"    # Violet
  )

  enr_sub[, ranking_label := ranking_labels[ranking]]
  enr_sub[, ranking_label := factor(ranking_label,
    levels = c("Convergence score", "Source count", "S1 only (human bulk)", "S3 only (genetic)"))]

  p_c <- ggplot(enr_sub, aes(x = k, y = fold_enrichment,
                              color = ranking_label, linetype = ranking_label)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.2, shape = 16) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_x_log10(labels = comma) +
    scale_color_manual(values = c(
      "Convergence score" = "#880E4F",
      "Source count"       = "#7B1FA2",
      "S1 only (human bulk)" = "#1565C0",
      "S3 only (genetic)"  = "#C2185B"
    )) +
    scale_linetype_manual(values = c(
      "Convergence score" = "solid",
      "Source count"       = "dashed",
      "S1 only (human bulk)" = "dotted",
      "S3 only (genetic)"  = "dotdash"
    )) +
    labs(
      x = "Top-k genes ranked",
      y = "Fold enrichment\n(DGIdb druggable)",
      title = "Cumulative enrichment for druggable genes",
      color = "Ranking", linetype = "Ranking"
    ) +
    theme_masld() +
    theme(
      legend.position = c(0.7, 0.75),
      legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
      legend.key.width = unit(0.7, "cm"),
      legend.text = element_text(size = 5.5)
    )
}

# ==========================================================================
# (d) Source addition bootstrap (greedy incremental, DGIdb at top-200)
# ==========================================================================
comp_file <- file.path(ME, "complementarity_pvalues.csv")
if (file.exists(comp_file)) {
  comp_dt <- fread(comp_file)

  # Use V1_DGIdb validation set
  comp_dgidb <- comp_dt[validation == "V1_DGIdb"]

  # enrichment_all = with all 6 sources; enrichment_without = with 5 sources (dropping this one)
  # improvement = enrichment_all - enrichment_without (positive = source improves enrichment)
  # perm_pvalue = permutation test p-value for contribution

  comp_dgidb[, source_label := SOURCE_LABELS_SHORT[source]]
  comp_dgidb[, source := factor(source, levels = paste0("S", 1:6))]

  # Significance labels
  comp_dgidb[, sig_label := ifelse(perm_pvalue < 0.01, "**",
                             ifelse(perm_pvalue < 0.05, "*", "ns"))]

  p_d <- ggplot(comp_dgidb, aes(x = source, y = improvement, fill = source)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "gray30") +
    geom_text(aes(label = sig_label,
                  y = ifelse(improvement >= 0, improvement + 0.15, improvement - 0.15)),
              size = 2.5, fontface = "bold") +
    geom_text(aes(label = sprintf("p=%.3f", perm_pvalue)),
              y = -4.5, size = 1.8, color = "gray40") +
    scale_x_discrete(labels = SOURCE_LABELS_SHORT) +
    scale_fill_manual(values = source_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.1))) +
    labs(
      x = NULL,
      y = "Enrichment change\n(top-200 DGIdb)",
      title = "Source contribution to enrichment"
    ) +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 5.5))
}

# ==========================================================================
# (e) Convergence score waterfall for key drug targets
# ==========================================================================
bayes_file <- file.path(ME, "bayesian_posterior.csv")
if (file.exists(bayes_file)) {
  bayes_dt <- fread(bayes_file)

  # Extract key drug targets
  targets <- c("THRB", "NR1H4", "PPARA")
  target_dt <- bayes_dt[human_symbol %in% targets]

  if (nrow(target_dt) > 0) {
    # Melt delta columns to long format (6 sources, no mouse)
    delta_cols <- paste0("delta_S", 1:6, "_",
                         c("human_bulk", "genetic",
                           "essentiality", "epigenomic", "spatial", "singlecell"))
    delta_long <- melt(target_dt,
                       id.vars = c("human_symbol", "posterior_odds", "posterior_rank"),
                       measure.vars = delta_cols,
                       variable.name = "source_var", value.name = "delta_log_odds")

    # Extract source ID from variable name
    delta_long[, source := sub("delta_(S[0-9])_.*", "\\1", source_var)]
    delta_long[, source_label := SOURCE_LABELS_SHORT[source]]
    delta_long[, source := factor(source, levels = paste0("S", 1:6))]

    # Add rank annotation
    delta_long[, gene_label := paste0(human_symbol, "\n(rank ", posterior_rank, ")")]
    delta_long[, gene_label := factor(gene_label,
      levels = unique(gene_label[order(match(human_symbol, targets))]))]

    # Compute cumulative positions for waterfall stacking
    delta_long <- delta_long[order(gene_label, source)]
    delta_long[, cum_end := cumsum(delta_log_odds), by = gene_label]
    delta_long[, cum_start := cum_end - delta_log_odds]

    # Identify positive vs negative contributions
    delta_long[, direction := ifelse(delta_log_odds >= 0, "positive", "negative")]

    p_e <- ggplot(delta_long, aes(x = source, fill = source)) +
      geom_rect(aes(xmin = as.numeric(source) - 0.35,
                     xmax = as.numeric(source) + 0.35,
                     ymin = cum_start, ymax = cum_end),
                show.legend = FALSE) +
      geom_hline(yintercept = 0, linewidth = 0.3, color = "gray30") +
      facet_wrap(~gene_label, nrow = 1, scales = "free_y") +
      scale_x_discrete(labels = paste0("S", 1:6)) +
      scale_fill_manual(values = source_colors) +
      labs(
        x = "Evidence source",
        y = "Delta log-odds",
        title = "Per-source Bayesian contribution for key targets"
      ) +
      theme_masld() +
      theme(
        strip.text = element_text(size = 6, face = "bold"),
        axis.text.x = element_text(size = 5.5)
      )
  }
}

# ==========================================================================
# Compose figure with patchwork
# ==========================================================================
figS <- (p_a | p_b) /
         (p_c | p_d) /
         p_e +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

save_fig_tall(figS, OUT, height = 9)
message("Saved: ", OUT)
