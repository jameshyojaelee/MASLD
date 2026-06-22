#!/usr/bin/env Rscript
# =============================================================================
# Fig 1 LOO-CV panels: Leave-one-study-out cross-validation visualizations.
# Generates standalone PDFs:
#   A. LOO-CV stability (DEG count + Spearman rho)
#   C. Gene-level LOO robustness heatmap
# (The former Panel B "AUROC cross-study prediction heatmap"
#  [fig1_auroc_heatmap.pdf] was retired 2026-06-12 — stale panel.)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ─── Shared: cohort accession mapping ────────────────────────────────────────
STUDY_NAMES <- c(
  GSE213621  = "GSE213621",
  GSE135251  = "GSE135251",
  GSE130970  = "GSE130970",
  GSE162694  = "GSE162694",
  GSE174478  = "GSE174478",
  GSE193066  = "GSE193066",
  GSE240729  = "GSE240729",
  GSE126848  = "GSE126848"
)

LOO_DIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/loo_cv")

PANEL_DIR <- file.path(FIG1_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PANEL A: LOO-CV Stability
# =============================================================================
cat("\n── Panel A: LOO-CV stability ──\n")

loo_summary_file <- file.path(LOO_DIR, "loo_cv_summary.csv")
if (!file.exists(loo_summary_file)) {
  cat("WARNING: loo_cv_summary.csv not found. Skipping Panel A.\n")
} else {
  loo_summary <- fread(loo_summary_file)
  loo_summary[, author := STUDY_NAMES[held_out]]

  # Full model DEG count
  dream <- load_dream_results()
  full_n_degs <- sum(dream$bulk_padj < 0.1)

  # Load QC-passing sample counts (matches actual dream model input)
  qc_file <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv")
  if (file.exists(qc_file)) {
    qc <- fread(qc_file)
    cohort_n <- qc[pass_technical == TRUE, .(n_samples = .N), by = dataset]
  } else {
    # Fallback to unified metadata if QC report unavailable
    meta <- load_metadata()
    cohort_n <- meta[, .(n_samples = .N), by = dataset]
    cat("  NOTE: Using pre-QC sample counts (QC report not found)\n")
  }

  loo_summary <- merge(loo_summary, cohort_n, by.x = "held_out", by.y = "dataset", all.x = TRUE)

  # Sort by sample count
  loo_summary[, author := factor(author, levels = author[order(n_samples)])]

  # Left: DEG count bar chart
  p_degs <- ggplot(loo_summary, aes(x = author, y = n_degs)) +
    geom_col(fill = masld_colors$down, width = 0.7) +
    geom_hline(yintercept = full_n_degs, linetype = "dashed", color = masld_colors$up, linewidth = 0.5) +
    annotate("text", x = Inf, y = full_n_degs, label = paste0("Full model: ", comma(full_n_degs)),
             hjust = 1.1, vjust = -0.5, size = 2, color = masld_colors$up) +
    labs(x = "Held-out cohort", y = "DEGs (padj < 0.1)") +
    scale_y_continuous(labels = comma) +
    coord_flip() +
    theme_masld() +
    ggtitle(sprintf("Mean recovery: %.0f%%", mean(loo_summary$pct_full_recovered)))

  # Right: Spearman rho dot plot
  p_rho <- ggplot(loo_summary, aes(x = spearman_rho, y = author)) +
    geom_point(aes(size = n_samples), color = masld_colors$up, shape = 16) +
    geom_vline(xintercept = 1.0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_size_continuous(range = c(2, 5), name = "N samples") +
    labs(x = "Spearman rho vs full model", y = NULL) +
    theme_masld() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    ggtitle(sprintf("Mean rho = %.3f", mean(loo_summary$spearman_rho)))

  p_stability <- p_degs + p_rho + plot_layout(widths = c(2, 1.2))

  save_fig(p_stability, file.path(PANEL_DIR, "fig1_loo_stability.pdf"),
           width = fig_full_width, height = 3)
  cat("Saved:", file.path(PANEL_DIR, "fig1_loo_stability.pdf"), "\n")
}

# =============================================================================
# PANEL B (AUROC cross-study heatmap) RETIRED 2026-06-12.
# The standalone fig1_auroc_heatmap.pdf was a stale panel and is no longer
# generated here. (fig1_compact.R still renders its own inline AUROC sub-panel
# from auroc_matrix.csv -- that one is unaffected.)
# =============================================================================

# =============================================================================
# PANEL C: Gene-Level LOO Robustness
# =============================================================================
cat("\n── Panel C: Gene-level LOO robustness ──\n")

gene_stab_file <- file.path(LOO_DIR, "loo_cv_per_gene.csv")
if (!file.exists(gene_stab_file)) {
  cat("WARNING: loo_cv_per_gene.csv not found. Skipping Panel C.\n")
} else {
  gene_stab <- fread(gene_stab_file)
  dream <- load_dream_results()

  # Top 30 genes by |dream t-stat| that are also DEGs
  top_genes <- dream[bulk_padj < 0.1][order(-abs(t))][1:30]

  # Build significance grid: gene × LOO iteration
  COHORTS <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694",
               "GSE174478", "GSE193066", "GSE213621", "GSE240729")

  grid_list <- list()
  for (study in COHORTS) {
    loo_file <- file.path(LOO_DIR, paste0("dream_loo_", study, ".csv"))
    if (!file.exists(loo_file)) next
    loo <- fread(loo_file)
    loo_sig <- loo[gene %in% top_genes$gene, .(gene, sig = as.character(padj < 0.1))]
    loo_sig[, cohort := STUDY_NAMES[study]]
    grid_list[[study]] <- loo_sig
  }

  # Add full model column
  full_sig <- dream[gene %in% top_genes$gene, .(gene, sig = as.character(bulk_padj < 0.1))]
  full_sig[, cohort := "Full model"]
  grid_list[["full"]] <- full_sig

  grid <- rbindlist(grid_list)

  # Add symbols
  gene_map <- dream[gene %in% top_genes$gene, .(gene, symbol)]
  grid <- merge(grid, gene_map, by = "gene", all.x = TRUE)
  grid[is.na(symbol), symbol := gene]

  # Get robustness for sorting
  gene_rob <- gene_stab[gene %in% top_genes$gene, .(gene, n_loo_sig)]
  grid <- merge(grid, gene_rob, by = "gene", all.x = TRUE)

  # Order genes by robustness then |t-stat|
  gene_order <- merge(top_genes[, .(gene, t)], gene_rob, by = "gene", all.x = TRUE)
  gene_order[is.na(n_loo_sig), n_loo_sig := 0]
  gene_order <- gene_order[order(-n_loo_sig, -abs(t))]
  symbol_order <- gene_map$symbol[match(gene_order$gene, gene_map$gene)]

  grid[, symbol := factor(symbol, levels = rev(symbol_order))]

  # Order cohorts: alphabetical author names + "Full model" at right
  cohort_levels <- c(sort(unique(STUDY_NAMES[COHORTS])), "Full model")
  grid[, cohort := factor(cohort, levels = cohort_levels)]

  # Count robust genes
  n_robust <- sum(gene_rob$n_loo_sig == length(COHORTS), na.rm = TRUE)

  p_robust <- ggplot(grid, aes(x = cohort, y = symbol, fill = sig)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_manual(values = c("TRUE" = masld_colors$up, "FALSE" = masld_colors$ns),
                      labels = c("TRUE" = "Significant", "FALSE" = "Not sig."),
                      na.value = "gray95",
                      name = "padj < 0.1") +
    labs(x = NULL, y = NULL,
         title = sprintf("%d/%d top genes significant in all %d LOO iterations",
                         n_robust, nrow(top_genes), length(COHORTS))) +
    theme_masld() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  # fig1_gene_robustness.pdf removed — not used in publication
}

cat("\nAll LOO-CV panels complete!\n")
