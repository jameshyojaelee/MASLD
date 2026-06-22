#!/usr/bin/env Rscript
# 06_meta_analysis.R
# ---------------------------------------------------------------------------
# Random-effects meta-analysis across per-study DE results (Disease vs Control)
# using metafor::rma. Serves as a REFERENCE COMPARISON to the canonical dream
# mega-analysis (Script 05). Dream remains the primary/canonical method for
# all downstream analyses.
#
# Outputs:
#   results/integration/meta_analysis_results.csv
#   results/integration/dream_vs_meta_concordance.csv
#   results/integration/meta_influence_diagnostics.csv
#   results/integration/meta_heterogeneity_summary.pdf
#   results/integration/dream_vs_meta_scatter.pdf
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(metafor)
  library(yaml)
  library(ggplot2)
})

PROJECT_ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
BASE <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts")
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/integration")
dir.create(RDIR, recursive = TRUE, showWarnings = FALSE)

MIN_DATASETS <- 3L  # Require gene present in >= 3 datasets

# ==========================================================================
#  1. Load per-study DE results (same cohort selection as dream Script 05)
# ==========================================================================
cfg_path <- file.path(PROJECT_ROOT, "config/human_datasets.yaml")
datasets_cfg <- yaml.load_file(cfg_path)$datasets

# Select datasets with has_controls=true AND include_in_mega=true
datasets <- names(Filter(function(ds) {
  isTRUE(ds$has_controls) && isTRUE(ds$de$include_in_mega)
}, datasets_cfg))

# Keep only datasets for which per-study results exist
datasets <- Filter(function(ds) {
  file.exists(file.path(INT, "results/per_study", paste0(ds, "_de_results.csv")))
}, datasets)

cat("Datasets included in meta-analysis:", paste(datasets, collapse = ", "), "\n")
cat("Number of datasets:", length(datasets), "\n\n")

per_study <- rbindlist(lapply(datasets, function(ds) {
  f <- file.path(INT, "results/per_study", paste0(ds, "_de_results.csv"))
  dt <- fread(f)
  # Strip ENSEMBL version suffixes (.N) so gene IDs match dream_results.csv
  dt[, gene := sub("\\.[0-9]+$", "", gene)]
  # Use UNMODERATED SE (stdev.unscaled * sigma) saved by Script 02.
  # The old formula SE = |logFC / t| gives moderated SE (eBayes-shrunk),
  # which biases tau^2 downward in metafor's REML estimator.
  if ("SE_unmoderated" %in% names(dt)) {
    dt[, .(gene, logFC, t, SE = SE_unmoderated, dataset = ds)]
  } else {
    # Fallback for legacy per-study files that lack SE_unmoderated.
    warning(ds, ": SE_unmoderated column missing — falling back to moderated |logFC/t|. ",
            "Re-run Script 02 to produce unmoderated SE.", immediate. = TRUE)
    dt[, .(gene, logFC, t, SE = abs(logFC / t), dataset = ds)]
  }
}))

cat("Total gene-dataset pairs:", nrow(per_study), "\n")

# ==========================================================================
#  2. Random-effects meta-analysis per gene (REML)
# ==========================================================================
gene_counts <- per_study[!is.na(SE) & is.finite(SE) & SE > 0, .N, by = gene]
genes_eligible <- gene_counts[N >= MIN_DATASETS, gene]
cat("Genes in >=", MIN_DATASETS, "datasets:", length(genes_eligible), "\n")
cat("Running meta-analysis...\n")

meta_results <- rbindlist(lapply(genes_eligible, function(g) {
  sub <- per_study[gene == g & !is.na(SE) & is.finite(SE) & SE > 0]
  if (nrow(sub) < MIN_DATASETS) return(NULL)

  tryCatch({
    fit <- rma(yi = sub$logFC, sei = sub$SE, method = "REML")
    data.table(
      gene         = g,
      meta_logFC   = as.numeric(fit$beta),
      meta_SE      = fit$se,
      meta_pval    = fit$pval,
      meta_I2      = fit$I2,
      meta_tau2    = fit$tau2,
      Q_stat       = fit$QE,
      Q_pval       = fit$QEp,
      n_datasets   = nrow(sub),
      datasets_included = paste(sub$dataset, collapse = ";")
    )
  }, error = function(e) NULL)
}))

meta_results[, meta_padj := p.adjust(meta_pval, method = "BH")]
setorder(meta_results, meta_padj)

cat("\n===== META-ANALYSIS RESULTS =====\n")
cat("Total genes analyzed:", nrow(meta_results), "\n")
sig_01  <- meta_results[meta_padj < 0.1]
sig_lfc <- meta_results[meta_padj < 0.1 & abs(meta_logFC) >= 0.5]
cat("DEGs (padj < 0.1):                ", nrow(sig_01), "\n")
cat("DEGs (padj < 0.1, |logFC| >= 0.5): ", nrow(sig_lfc), "\n")
cat("  Up:  ", sum(sig_lfc$meta_logFC > 0), "\n")
cat("  Down:", sum(sig_lfc$meta_logFC < 0), "\n")
cat("\nHeterogeneity (I^2):\n")
cat("  Median:", round(median(meta_results$meta_I2, na.rm = TRUE), 1), "%\n")
cat("  IQR:   ", paste(round(quantile(meta_results$meta_I2, c(0.25, 0.75), na.rm = TRUE), 1), collapse = " - "), "%\n")
cat("  Genes with I^2 > 50%:", sum(meta_results$meta_I2 > 50, na.rm = TRUE),
    sprintf("(%.1f%%)", 100 * mean(meta_results$meta_I2 > 50, na.rm = TRUE)), "\n")
cat("  Genes with I^2 > 75%:", sum(meta_results$meta_I2 > 75, na.rm = TRUE),
    sprintf("(%.1f%%)", 100 * mean(meta_results$meta_I2 > 75, na.rm = TRUE)), "\n")

fwrite(meta_results, file.path(RDIR, "meta_analysis_results.csv"))
cat("\nSaved: meta_analysis_results.csv\n")

# ==========================================================================
#  3. Concordance with dream results
# ==========================================================================
cat("\n===== CONCORDANCE WITH DREAM =====\n")

dream_file <- file.path(RDIR, "dream_results.csv")
if (!file.exists(dream_file)) {
  cat("WARNING: dream_results.csv not found — skipping concordance analysis.\n")
  cat("Run Script 05 first, then re-run Script 06.\n")
} else {
  dream <- fread(dream_file)
  # dream columns: gene, logFC, padj (plus symbol, etc.)
  dream[, gene := sub("\\.[0-9]+$", "", gene)]
  dream_sub <- dream[, .(gene, dream_logFC = logFC, dream_padj = padj)]  # C2-OK-sensitivity
  meta_sub  <- meta_results[, .(gene, meta_logFC, meta_padj, meta_I2)]

  concordance <- merge(dream_sub, meta_sub, by = "gene", all = TRUE)
  concordance[, `:=`(
    dream_sig = !is.na(dream_padj) & dream_padj < 0.1,  # C2-OK-sensitivity
    meta_sig  = !is.na(meta_padj)  & meta_padj  < 0.1,
    direction_concordant = sign(dream_logFC) == sign(meta_logFC)  # C2-OK-sensitivity
  )]

  # Genes tested in both
  both <- concordance[!is.na(dream_logFC) & !is.na(meta_logFC)]  # C2-OK-sensitivity
  dream_only <- concordance[!is.na(dream_logFC) & is.na(meta_logFC)]  # C2-OK-sensitivity
  meta_only  <- concordance[is.na(dream_logFC) & !is.na(meta_logFC)]  # C2-OK-sensitivity

  # Spearman correlation of logFCs
  rho <- cor(both$dream_logFC, both$meta_logFC, method = "spearman", use = "complete.obs")  # C2-OK-sensitivity
  pearson_r <- cor(both$dream_logFC, both$meta_logFC, method = "pearson", use = "complete.obs")  # C2-OK-sensitivity

  # DEG overlap (padj < 0.1)
  dream_degs <- both[dream_sig == TRUE, gene]  # C2-OK-sensitivity
  meta_degs  <- both[meta_sig == TRUE, gene]
  overlap    <- intersect(dream_degs, meta_degs)
  jaccard    <- length(overlap) / length(union(dream_degs, meta_degs))

  # Direction concordance among genes significant in BOTH
  both_sig <- both[dream_sig == TRUE & meta_sig == TRUE]  # C2-OK-sensitivity
  dir_concord <- mean(both_sig$direction_concordant, na.rm = TRUE)

  cat("Genes in both analyses:          ", nrow(both), "\n")
  cat("Genes dream-only (not in meta):  ", nrow(dream_only), "\n")
  cat("Genes meta-only (not in dream):  ", nrow(meta_only), "\n")
  cat("\nlogFC correlation (common genes):\n")
  cat("  Spearman rho:", round(rho, 4), "\n")
  cat("  Pearson r:   ", round(pearson_r, 4), "\n")
  cat("\nDEG overlap (padj < 0.1, common genes):\n")
  cat("  Dream DEGs:", length(dream_degs), "\n")
  cat("  Meta DEGs: ", length(meta_degs), "\n")
  cat("  Overlap:   ", length(overlap), "\n")
  cat("  Jaccard:   ", round(jaccard, 4), "\n")
  cat("\nDirection concordance (genes sig in BOTH):\n")
  cat("  N:", nrow(both_sig), "\n")
  cat("  Concordant:", round(100 * dir_concord, 1), "%\n")

  fwrite(concordance, file.path(RDIR, "dream_vs_meta_concordance.csv"))
  cat("\nSaved: dream_vs_meta_concordance.csv\n")

  # ========================================================================
  #  4. Scatter plot: dream vs meta logFC
  # ========================================================================
  cat("\nGenerating concordance plots...\n")

  pdf(file.path(RDIR, "dream_vs_meta_scatter.pdf"), width = 12, height = 5)

  # Panel 1: logFC scatter
  p1 <- ggplot(both, aes(x = dream_logFC, y = meta_logFC)) +  # C2-OK-sensitivity
    geom_point(aes(color = meta_I2), size = 0.4, alpha = 0.5) +
    scale_color_viridis_c(name = expression(I^2), limits = c(0, 100)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "steelblue", linewidth = 0.5) +
    labs(
      title = "Disease vs Control: Dream vs Meta logFC",
      subtitle = sprintf("N = %s genes | Spearman rho = %.3f | Pearson r = %.3f",
                         format(nrow(both), big.mark = ","), rho, pearson_r),
      x = "Dream logFC", y = "Meta-analysis logFC"
    ) +
    theme_minimal(base_size = 11) +
    coord_fixed()

  # Panel 2: Venn-style DEG comparison (bar chart)
  venn_dt <- data.table(
    category = c("Dream only", "Both", "Meta only"),
    n = c(length(setdiff(dream_degs, meta_degs)),
          length(overlap),
          length(setdiff(meta_degs, dream_degs)))
  )
  venn_dt[, category := factor(category, levels = c("Dream only", "Both", "Meta only"))]

  p2 <- ggplot(venn_dt, aes(x = category, y = n, fill = category)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = format(n, big.mark = ",")), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c("Dream only" = "#4393C3", "Both" = "#7B3294", "Meta only" = "#D6604D"),
                      guide = "none") +
    labs(
      title = "DEG Overlap (padj < 0.1)",
      subtitle = sprintf("Jaccard = %.3f | Direction concordance = %.1f%%", jaccard, 100 * dir_concord),
      x = "", y = "Number of DEGs"
    ) +
    theme_minimal(base_size = 11) +
    ylim(0, max(venn_dt$n) * 1.15)

  library(gridExtra)
  grid.arrange(p1, p2, ncol = 2, widths = c(1.3, 1))
  dev.off()
  cat("Saved: dream_vs_meta_scatter.pdf\n")
}

# ==========================================================================
#  5. Leave-one-out influence diagnostics
# ==========================================================================
cat("\n===== LEAVE-ONE-OUT INFLUENCE ANALYSIS =====\n")

# Run LOO on top 500 DEGs (by meta_padj) to keep runtime manageable
top_genes <- head(meta_results[meta_padj < 0.1], 500)$gene
if (length(top_genes) == 0) {
  cat("No significant DEGs — skipping LOO analysis.\n")
} else {
  cat("Running LOO for top", length(top_genes), "DEGs...\n")

  loo_results <- rbindlist(lapply(top_genes, function(g) {
    sub <- per_study[gene == g & !is.na(SE) & is.finite(SE) & SE > 0]
    if (nrow(sub) < MIN_DATASETS) return(NULL)

    # Full model
    full_fit <- tryCatch(rma(yi = sub$logFC, sei = sub$SE, method = "REML"),
                         error = function(e) NULL)
    if (is.null(full_fit)) return(NULL)

    # Drop each study
    rbindlist(lapply(seq_len(nrow(sub)), function(i) {
      sub_loo <- sub[-i]
      if (nrow(sub_loo) < 2) return(NULL)
      loo_fit <- tryCatch(rma(yi = sub_loo$logFC, sei = sub_loo$SE, method = "REML"),
                          error = function(e) NULL)
      if (is.null(loo_fit)) return(NULL)
      data.table(
        gene           = g,
        dropped_study  = sub$dataset[i],
        full_logFC     = as.numeric(full_fit$beta),
        loo_logFC      = as.numeric(loo_fit$beta),
        full_padj      = p.adjust(full_fit$pval, method = "BH"),
        loo_pval       = loo_fit$pval,
        logFC_change   = as.numeric(loo_fit$beta) - as.numeric(full_fit$beta),
        pct_logFC_change = abs(as.numeric(loo_fit$beta) - as.numeric(full_fit$beta)) /
                           (abs(as.numeric(full_fit$beta)) + 1e-10) * 100
      )
    }))
  }))

  if (nrow(loo_results) > 0) {
    # Summarize: which study, when dropped, causes the largest average LFC change?
    study_influence <- loo_results[, .(
      mean_abs_lfc_change = mean(abs(logFC_change)),
      median_pct_change   = median(pct_logFC_change),
      n_genes_tested      = .N
    ), by = dropped_study]
    setorder(study_influence, -mean_abs_lfc_change)

    cat("\nStudy influence (mean |delta_logFC| when dropped):\n")
    print(study_influence)

    # Flag genes where dropping ANY study flips significance (pval crosses 0.05)
    loo_results[, sig_flip := (full_padj < 0.1) != (loo_pval < 0.1 / length(top_genes))]

    fwrite(loo_results, file.path(RDIR, "meta_influence_diagnostics.csv"))
    cat("\nSaved: meta_influence_diagnostics.csv\n")
  } else {
    cat("LOO analysis returned no results.\n")
  }
}

# ==========================================================================
#  6. Heterogeneity summary plot
# ==========================================================================
cat("\nGenerating heterogeneity summary plot...\n")

pdf(file.path(RDIR, "meta_heterogeneity_summary.pdf"), width = 12, height = 5)

# Panel 1: I^2 distribution
p_i2 <- ggplot(meta_results, aes(x = meta_I2)) +
  geom_histogram(bins = 50, fill = "#5E4FA2", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = c(25, 50, 75), linetype = "dashed", color = "grey40") +
  annotate("text", x = 12.5, y = Inf, label = "Low", vjust = 1.5, size = 3, color = "grey40") +
  annotate("text", x = 37.5, y = Inf, label = "Moderate", vjust = 1.5, size = 3, color = "grey40") +
  annotate("text", x = 62.5, y = Inf, label = "Substantial", vjust = 1.5, size = 3, color = "grey40") +
  annotate("text", x = 87.5, y = Inf, label = "Considerable", vjust = 1.5, size = 3, color = "grey40") +
  labs(
    title = expression(paste("Cross-Study Heterogeneity (", I^2, " Distribution)")),
    subtitle = sprintf("N = %s genes | Median I^2 = %.1f%% | %.1f%% genes with I^2 > 50%%",
                       format(nrow(meta_results), big.mark = ","),
                       median(meta_results$meta_I2, na.rm = TRUE),
                       100 * mean(meta_results$meta_I2 > 50, na.rm = TRUE)),
    x = expression(I^2 ~ "(%)"), y = "Number of genes"
  ) +
  theme_minimal(base_size = 11)

# Panel 2: I^2 vs meta_logFC — are high-heterogeneity genes systematically different?
p_i2_lfc <- ggplot(meta_results, aes(x = abs(meta_logFC), y = meta_I2)) +
  geom_point(aes(color = meta_padj < 0.1), size = 0.4, alpha = 0.4) +
  scale_color_manual(values = c("TRUE" = "#D6604D", "FALSE" = "grey70"),
                     labels = c("TRUE" = "padj < 0.1", "FALSE" = "NS"),
                     name = "") +
  labs(
    title = expression(paste(I^2, " vs Effect Size")),
    subtitle = "High heterogeneity does not imply large or small effects",
    x = "|Meta logFC|", y = expression(I^2 ~ "(%)")
  ) +
  theme_minimal(base_size = 11)

grid.arrange(p_i2, p_i2_lfc, ncol = 2)
dev.off()
cat("Saved: meta_heterogeneity_summary.pdf\n")

cat("\n=== Script 06 complete ===\n")
