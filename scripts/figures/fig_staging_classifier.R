#!/usr/bin/env Rscript
# fig_staging_classifier.R
# Publication figures for the MASLD staging classifier pipeline
#
# Generates multi-panel figures from staging_classifier results
# Output: figures/supplementary/figS10_prediction/figS_staging_classifier.pdf
#
# Usage: Rscript fig_staging_classifier.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
SDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier")
FIGDIR <- FIGS10_DIR
dir.create(file.path(FIGDIR, "panels"), showWarnings = FALSE, recursive = TRUE)

cat("=== Staging Classifier Figures ===\n")
cat("Input:", SDIR, "\n")
cat("Output:", FIGDIR, "\n\n")

# --- Publication theme ---
theme_pub <- theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "plain", size = 11),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    strip.text = element_text(size = 9, face = "plain"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 10, 5, 5)
  )

safe_fread <- function(path, desc = "") {
  if (file.exists(path)) {
    dt <- fread(path)
    cat("  Loaded", desc, ":", nrow(dt), "rows\n")
    return(dt)
  }
  cat("  MISSING:", desc, "-", path, "\n")
  return(NULL)
}

# ============================================================
# Load data
# ============================================================
cat("Loading results...\n")
final_report <- safe_fread(file.path(SDIR, "stage_classifier_final_report.csv"), "final report")
model_ranking <- safe_fread(file.path(SDIR, "model_ranking.csv"), "model ranking")
panel_curve <- safe_fread(file.path(SDIR, "panel_performance_curve.csv"), "panel curve")
published_bench <- safe_fread(file.path(SDIR, "published_panel_benchmark.csv"), "published bench")
specificity <- safe_fread(file.path(SDIR, "stage_specificity_scores.csv"), "specificity")
unique_sigs <- safe_fread(file.path(SDIR, "stage_unique_signatures.csv"), "unique sigs")
conformal <- safe_fread(file.path(SDIR, "conformal_prediction_sets.csv"), "conformal")
imputation <- safe_fread(file.path(SDIR, "imputation_confidence.csv"), "imputation confidence")
embedding_comp <- safe_fread(file.path(SDIR, "embedding_vs_raw_comparison.csv"), "embedding comparison")
attention_wt <- safe_fread(file.path(SDIR, "attention_weights_all.csv"), "attention weights")
plan1_sweep_fib <- safe_fread(file.path(SDIR, "plan1_sweep_fibrosis.csv"), "Plan1 sweep fibrosis")
plan1_sweep_binary <- safe_fread(file.path(SDIR, "plan1_sweep_binary.csv"), "Plan1 sweep binary")
cascade <- safe_fread(file.path(SDIR, "cascade_error_analysis.csv"), "cascade analysis")
ordinal_fib <- safe_fread(file.path(SDIR, "ordinal_loco_fibrosis.csv"), "ordinal fibrosis")

# ============================================================
# Panel A: Model comparison heatmap (Plan x Model x Target)
# ============================================================
cat("\n=== Panel A: Model comparison ===\n")
p_a <- NULL
if (!is.null(final_report) && nrow(final_report) > 0) {
  # Focus on key targets
  hm <- copy(final_report)
  hm[, display := fifelse(!is.na(auroc) & auroc > 0, auroc,
                   fifelse(!is.na(qwk) & qwk > 0, qwk, accuracy))]
  hm[, metric_type := fifelse(!is.na(auroc) & auroc > 0, "AUROC",
                       fifelse(!is.na(qwk) & qwk > 0, "QWK", "Accuracy"))]
  hm <- hm[!is.na(display) & display > 0]
  hm[, model_label := paste0(model, " (", plan, ")")]

  if (nrow(hm) > 0) {
    # Top 10 model-target combos
    top_hm <- hm[order(-display)][1:min(20, nrow(hm))]

    p_a <- ggplot(top_hm, aes(x = target, y = reorder(model_label, display), fill = display)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.3f", display)), size = 2.2, color = "black") +
      scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                           midpoint = 0.5, name = "Performance", limits = c(0, 1)) +
      labs(title = "A) LOCO Cross-Validation Performance",
           x = "", y = "") +
      theme_pub +
      theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
            axis.text.y = element_text(size = 7))

    ggsave(file.path(FIGDIR, "panels", "panel_a_model_heatmap.pdf"), p_a, width = 8, height = 6)
    cat("  Saved panel A\n")
  }
}

# ============================================================
# Panel B: Gene panel size vs performance
# ============================================================
cat("=== Panel B: Panel performance curve ===\n")
p_b <- NULL
if (!is.null(panel_curve) && nrow(panel_curve) > 0) {
  # Expect columns: n_genes (or panel_size), mean_auroc, sd_auroc
  # panel_curve has: panel_size, n_genes, mean_auroc, sd_auroc
  pc <- copy(panel_curve)
  # Use n_genes directly if present, otherwise panel_size
  if (!"n_genes" %in% names(pc) && "panel_size" %in% names(pc)) {
    setnames(pc, "panel_size", "n_genes")
  }
  if (!"auroc" %in% names(pc) && "mean_auroc" %in% names(pc)) {
    pc[, auroc := mean_auroc]
  } else if (!"auroc" %in% names(pc)) {
    pc[, auroc := NA_real_]
  }
  pc <- pc[!is.na(auroc) & !is.na(n_genes)]

  if (TRUE) {

    if (nrow(pc) > 0) {
      p_b <- ggplot(pc, aes(x = n_genes, y = auroc)) +
        geom_line(color = "#1a9850", linewidth = 1) +
        geom_point(size = 2.5, color = "#1a9850") +
        scale_x_log10(breaks = c(10, 15, 20, 25, 50, 100, 200, 500)) +
        coord_cartesian(ylim = c(0.5, 1.0)) +
        labs(title = "B) Gene Panel Size vs AUROC (Fibrosis ordinal)",
             x = "Number of genes", y = "Mean LOCO AUROC") +
        theme_pub

      # Add published panel benchmarks if available
      if (!is.null(published_bench) && nrow(published_bench) > 0) {
        pb <- copy(published_bench)
        pb_auroc_col <- intersect(names(pb), c("mean_auroc", "auroc"))[1]
        pb_ngene_col <- intersect(names(pb), c("n_genes", "panel_size"))[1]
        if (!is.na(pb_auroc_col) && !is.na(pb_ngene_col)) {
          setnames(pb, c(pb_ngene_col, pb_auroc_col), c("n_genes", "auroc"), skip_absent = TRUE)
          pb <- pb[!is.na(auroc)]
          if (nrow(pb) > 0) {
            p_b <- p_b +
              geom_point(data = pb, aes(x = n_genes, y = auroc),
                         shape = 17, size = 3, color = "#d73027") +
              geom_text(data = pb, aes(x = n_genes, y = auroc, label = panel_name),
                        nudge_y = 0.02, size = 2.5, color = "#d73027")
          }
        }
      }

      ggsave(file.path(FIGDIR, "panels", "panel_b_panel_curve.pdf"), p_b, width = 6, height = 4)
      cat("  Saved panel B\n")
    }
  }
}

# ============================================================
# Panel C: Stage-unique gene counts
# ============================================================
cat("=== Panel C: Stage-unique signatures ===\n")
p_c <- NULL
if (!is.null(unique_sigs) && nrow(unique_sigs) > 0) {
  counts <- unique_sigs[!is.na(unique_stage), .N, by = .(staging_axis, unique_stage)]
  counts[, stage_short := gsub("_vs_rest", "", unique_stage)]

  if (nrow(counts) > 0) {
    p_c <- ggplot(counts, aes(x = reorder(stage_short, -N), y = N, fill = staging_axis)) +
      geom_col(position = "dodge", width = 0.7) +
      scale_fill_manual(values = c("NAS" = "#4575b4", "Fibrosis" = "#d73027"), name = "Axis") +
      labs(title = "C) Stage-Unique Gene Signatures",
           x = "Stage", y = "Number of unique genes") +
      theme_pub +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(FIGDIR, "panels", "panel_c_unique_signatures.pdf"), p_c, width = 6, height = 4)
    cat("  Saved panel C\n")
  }
}

# ============================================================
# Panel D: Tau specificity distribution
# ============================================================
cat("=== Panel D: Tau specificity ===\n")
p_d <- NULL
if (!is.null(specificity) && "tau" %in% names(specificity)) {
  spec <- specificity[!is.na(tau)]
  if (nrow(spec) > 0) {
    p_d <- ggplot(spec, aes(x = tau, fill = staging_axis)) +
      geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
      scale_fill_manual(values = c("NAS" = "#4575b4", "Fibrosis" = "#d73027"), name = "Axis") +
      labs(title = "D) Gene Specificity Index (\u03C4)",
           x = "Tau specificity index", y = "Number of genes") +
      theme_pub

    ggsave(file.path(FIGDIR, "panels", "panel_d_tau_distribution.pdf"), p_d, width = 6, height = 4)
    cat("  Saved panel D\n")
  }
}

# ============================================================
# Panel E: Embedding vs Raw comparison
# ============================================================
cat("=== Panel E: Embedding vs Raw ===\n")
p_e <- NULL
if (!is.null(embedding_comp) && nrow(embedding_comp) > 0) {
  ec <- copy(embedding_comp)
  auroc_emb_col <- intersect(names(ec), c("emb_auroc", "embedding_auroc", "auroc_embedding"))[1]
  auroc_raw_col <- intersect(names(ec), c("raw_auroc", "gene_auroc", "auroc_raw"))[1]

  if (!is.na(auroc_emb_col) && !is.na(auroc_raw_col)) {
    setnames(ec, c(auroc_emb_col, auroc_raw_col), c("emb", "raw"), skip_absent = TRUE)
    ec <- ec[!is.na(emb) & !is.na(raw)]

    if (nrow(ec) > 0) {
      p_e <- ggplot(ec, aes(x = raw, y = emb)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
        geom_point(aes(color = target), size = 2.5, alpha = 0.7) +
        labs(title = "E) VAE Embeddings vs Raw Features",
             x = "Raw gene AUROC", y = "VAE embedding AUROC") +
        theme_pub +
        theme(legend.position = "bottom")

      ggsave(file.path(FIGDIR, "panels", "panel_e_embedding_vs_raw.pdf"), p_e, width = 5, height = 5)
      cat("  Saved panel E\n")
    }
  }
}

# ============================================================
# Panel F: Conformal prediction set sizes
# ============================================================
cat("=== Panel F: Conformal prediction ===\n")
p_f <- NULL
if (!is.null(conformal) && nrow(conformal) > 0) {
  set_size_col <- intersect(names(conformal), c("set_size", "prediction_set_size", "mean_set_size"))[1]
  if (!is.na(set_size_col)) {
    cf <- copy(conformal)
    setnames(cf, set_size_col, "set_size", skip_absent = TRUE)
    cf <- cf[!is.na(set_size)]

    if (nrow(cf) > 0) {
      p_f <- ggplot(cf, aes(x = factor(set_size))) +
        geom_bar(fill = "#4575b4", alpha = 0.8) +
        labs(title = "F) Conformal Prediction Set Sizes",
             x = "Prediction set size (# possible stages)",
             y = "Number of samples") +
        theme_pub

      ggsave(file.path(FIGDIR, "panels", "panel_f_conformal_sets.pdf"), p_f, width = 5, height = 4)
      cat("  Saved panel F\n")
    }
  }
}

# ============================================================
# Panel G: Imputation confidence by dataset
# ============================================================
cat("=== Panel G: Imputation confidence ===\n")
p_g <- NULL
if (!is.null(imputation) && nrow(imputation) > 0) {
  conf_col <- intersect(names(imputation), c("confidence_level", "confidence", "label"))[1]
  ds_col <- intersect(names(imputation), c("dataset", "cohort"))[1]

  if (!is.na(conf_col)) {
    imp <- copy(imputation)
    setnames(imp, conf_col, "confidence", skip_absent = TRUE)
    if (!is.na(ds_col)) setnames(imp, ds_col, "dataset", skip_absent = TRUE)

    if ("dataset" %in% names(imp)) {
      imp_counts <- imp[, .N, by = .(dataset, confidence)]
      p_g <- ggplot(imp_counts, aes(x = dataset, y = N, fill = confidence)) +
        geom_col(position = "stack") +
        scale_fill_manual(values = c("High" = "#1a9850", "Medium" = "#fee08b", "Low" = "#d73027"),
                          name = "Confidence") +
        labs(title = "G) Label Imputation Confidence",
             x = "Dataset", y = "Number of samples") +
        theme_pub +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } else {
      imp_counts <- imp[, .N, by = confidence]
      p_g <- ggplot(imp_counts, aes(x = confidence, y = N, fill = confidence)) +
        geom_col() +
        scale_fill_manual(values = c("High" = "#1a9850", "Medium" = "#fee08b", "Low" = "#d73027"),
                          name = "Confidence") +
        labs(title = "G) Label Imputation Confidence",
             x = "Confidence level", y = "Number of samples") +
        theme_pub
    }

    ggsave(file.path(FIGDIR, "panels", "panel_g_imputation_confidence.pdf"), p_g, width = 6, height = 4)
    cat("  Saved panel G\n")
  }
}

# ============================================================
# Composite figure
# ============================================================
cat("\n=== Composing multi-panel figure ===\n")

panels <- list(p_a, p_b, p_c, p_d, p_e, p_f, p_g)
panels <- panels[!sapply(panels, is.null)]
cat("  Panels available:", length(panels), "of 7\n")

if (length(panels) >= 2) {
  # Arrange in a grid
  n_panels <- length(panels)
  ncol_layout <- min(3, n_panels)
  nrow_layout <- ceiling(n_panels / ncol_layout)

  composite <- wrap_plots(panels, ncol = ncol_layout) +
    plot_annotation(
      title = "MASLD Staging Classifier: Three-Plan Comparison",
      subtitle = "Ordinal (Plan 1) | Hierarchical (Plan 2) | DL Embedding (Plan 3)",
      theme = theme(
        plot.title = element_text(face = "plain", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray30")
      )
    )

  composite_path <- file.path(FIGDIR, "figS_staging_classifier.pdf")
  ggsave(composite_path, composite,
         width = ncol_layout * 5, height = nrow_layout * 4.5)
  cat("  Saved composite:", composite_path, "\n")
}

# Copy UMAP from VAE to figures dir
umap_src <- file.path(SDIR, "latent_space_umap.png")
if (file.exists(umap_src)) {
  file.copy(umap_src, file.path(FIGDIR, "panel_vae_umap.png"), overwrite = TRUE)
  cat("  Copied VAE UMAP to figures dir\n")
}

# Copy comparison_figures.pdf too
comp_src <- file.path(SDIR, "comparison_figures.pdf")
if (file.exists(comp_src)) {
  file.copy(comp_src, file.path(FIGDIR, "comparison_figures_from_script75.pdf"), overwrite = TRUE)
  cat("  Copied Script 75 comparison figures\n")
}

cat("\n=== Figure generation complete ===\n")
cat("Individual panels:", FIGDIR, "\n")
cat("Composite figure:", file.path(FIGDIR, "figS_staging_classifier.pdf"), "\n")
