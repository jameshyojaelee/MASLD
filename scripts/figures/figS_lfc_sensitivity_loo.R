#!/usr/bin/env Rscript
##############################################################################
# figS_lfc_sensitivity_loo.R  (new 2026-04-27; relocated 2026-05-07)
# Supp panel F: LOOCV-based motivation of the |log2FC| = 0.5 cutoff.
#
# Originally lived as Fig 1G; relocated to figS_lfc_sensitivity (the LFC-cutoff
# sensitivity supplement) on 2026-05-07 after the Fig 1G slot was taken over by
# the patient-level LFC concordance heatmap (fig1g_patient_lfc_combined.R).
#
# For each leave-one-cohort-out (LOO) refit of the dream mega-analysis, we
# count DEGs (padj < 0.05) at a sweep of |log2FC| cutoffs and measure
# recovery against the full mega-analysis. The cutoff at which across-fold
# recovery stabilizes (elbow) is the empirically defensible cutoff for the
# integrated DEG set.
#
# y-axis: held-out cohort, ordered by % of patients held out
# x-axis: |log2FC| cutoff
# fill: % of full-model DEGs recovered in the LOO fold
# text: # DEGs in the LOO fold at that cutoff
#
# A side panel shows the across-fold CV of DEG counts vs LFC, with the
# auto-detected elbow marked.
#
# Source (2026-06-27): canonical C2 limma-voom-qw LOO refits
#   results/integration/loo_cv_C2/lvqw_loo_C2_<COHORT>.csv  (carry shrunk_logFC+lfsr)
#   results/integration/canonical_deg_results.csv           (full-analysis n_full)
#
# Emitted per scale {raw, shrunk, shrunk_padjgate} (suffix appended):
#   panels/loo_cv_stability<suffix>.pdf
#   loo_cv_stability_data<suffix>.csv
#   loo_cv_stability_pG<suffix>.rds
# raw = padj<0.05 & |logFC|>cut; shrunk = lfsr<0.05 & |shrunk_logFC|>cut
# (PRIMARY); shrunk_padjgate = padj<0.05 & |shrunk_logFC|>cut.
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(yaml)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

INT_DIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
# Canonical C2 (limma-voom-qw) LOO refits; per-fold files already carry
# shrunk_logFC + lfsr (no ashr re-run needed). Repointed off the retired
# pre-C2 loo_cv/dream_loo_* source 2026-06-27.
LOO_DIR   <- file.path(INT_DIR, "loo_cv_C2")
OUT_DIR   <- file.path(FIG_SUPP, "figS_methods_validation/lfc_sensitivity")
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# Scale definitions: each emits a suffixed PDF/CSV/RDS triple.
#   raw            -> filter padj < 0.05 & |logFC|       > cut
#   shrunk         -> filter lfsr < 0.05 & |shrunk_logFC|> cut   (PRIMARY)
#   shrunk_padjgate-> filter padj < 0.05 & |shrunk_logFC|> cut   (effect-only swap)
SCALES <- list(
  raw = list(
    suffix    = "_raw",
    eff_col   = "logFC",
    sig_col   = "padj",
    sig_cut   = 0.05,
    sig_label = "padj < 0.05",
    lfc_word  = "raw |log2FC|"),
  shrunk = list(
    suffix    = "_shrunk",
    eff_col   = "shrunk_logFC",
    sig_col   = "lfsr",
    sig_cut   = 0.05,
    sig_label = "lfsr < 0.05",
    lfc_word  = "ashr-shrunk |log2FC|"),
  shrunk_padjgate = list(
    suffix    = "_shrunk_padjgate",
    eff_col   = "shrunk_logFC",
    sig_col   = "padj",
    sig_cut   = 0.05,
    sig_label = "padj < 0.05",
    lfc_word  = "ashr-shrunk |log2FC|")
)

# Cohort first-author labels (full atlas)
ALL_STUDY_NAMES <- c(
  GSE126848 = "GSE126848", GSE130970 = "GSE130970", GSE135251 = "GSE135251",
  GSE162694 = "GSE162694",   GSE174478 = "GSE174478", GSE193066 = "GSE193066",
  GSE213621 = "GSE213621",   GSE240729 = "GSE240729"
)
# Mega-analysis cohort set is the canonical source of truth in
# config/human_datasets.yaml (`include_in_mega: true`).
ycfg_path <- file.path(BASE, "config/human_datasets.yaml")
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega),
                             yaml::read_yaml(ycfg_path)$datasets))
STUDY_NAMES <- ALL_STUDY_NAMES[intersect(names(ALL_STUDY_NAMES), mega_cohorts)]
stopifnot("All mega cohorts must have first-author labels" =
            length(STUDY_NAMES) == length(mega_cohorts))

LFC_GRID <- c(0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.75, 1.0)

# -----------------------------------------------------------------------------
# Cohort sizes -> % patients held out per fold
# -----------------------------------------------------------------------------
qc <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv"))
mega_cohorts <- names(STUDY_NAMES)
n_per <- qc[pass_technical == TRUE & dataset %in% mega_cohorts, .N, by = dataset]
n_total <- sum(n_per$N)
n_per[, pct_held := round(100 * N / n_total, 1)]
setkey(n_per, dataset)
cat(sprintf("Total mega-eligible patients: %d across %d cohorts\n",
            n_total, length(mega_cohorts)))

# -----------------------------------------------------------------------------
# Load full dream model + 8 LOO refits
# -----------------------------------------------------------------------------
full_dream <- fread(file.path(INT_DIR, "canonical_deg_results.csv"))
setnames(full_dream, "adj.P.Val", "padj", skip_absent = TRUE)

loo <- rbindlist(lapply(mega_cohorts, function(ds) {
  dt <- fread(file.path(LOO_DIR, sprintf("lvqw_loo_C2_%s.csv", ds)))
  setnames(dt, "adj.P.Val", "padj", skip_absent = TRUE)
  dt[, held_out := ds]
  dt
}))

# -----------------------------------------------------------------------------
# Per-scale driver: sweep -> stability -> two-panel figure (suffixed outputs)
# -----------------------------------------------------------------------------
run_scale <- function(sc) {
  eff <- sc$eff_col; sig <- sc$sig_col; sig_cut <- sc$sig_cut
  OUT_PDF <- file.path(PANEL_DIR, sprintf("loo_cv_stability%s.pdf", sc$suffix))
  OUT_RDS <- file.path(OUT_DIR,   sprintf("loo_cv_stability_pG%s.rds", sc$suffix))
  OUT_CSV <- file.path(OUT_DIR,   sprintf("loo_cv_stability_data%s.csv", sc$suffix))

  cat(sprintf("\n================ scale = %s (%s, %s) ================\n",
              sc$suffix, sc$lfc_word, sc$sig_label))

  # Sweep: per fold x LFC cutoff -> #DEGs and recovery vs full
  sweep <- rbindlist(lapply(LFC_GRID, function(cut) {
    full_set <- full_dream[get(sig) < sig_cut & abs(get(eff)) > cut, gene]
    rbindlist(lapply(mega_cohorts, function(ds) {
      fold <- loo[held_out == ds]
      fold_set <- fold[get(sig) < sig_cut & abs(get(eff)) > cut, gene]
      inter <- length(intersect(fold_set, full_set))
      union <- length(union(fold_set, full_set))
      data.table(
        held_out      = ds,
        author        = STUDY_NAMES[ds],
        pct_held      = n_per[ds, pct_held],
        lfc_cutoff    = cut,
        n_full        = length(full_set),
        n_fold        = length(fold_set),
        n_overlap     = inter,
        pct_recovered = if (length(full_set)) 100 * inter / length(full_set) else NA_real_,
        jaccard       = if (union) inter / union else NA_real_
      )
    }))
  }))

  fwrite(sweep, OUT_CSV)
  cat("Saved: ", OUT_CSV, "\n", sep = "")

  # Across-fold stability metric -> empirical CV minimum
  stab <- sweep[, .(
    mean_recovery = mean(pct_recovered, na.rm = TRUE),
    sd_recovery   = sd(pct_recovered,   na.rm = TRUE),
    mean_jaccard  = mean(jaccard,        na.rm = TRUE),
    sd_jaccard    = sd(jaccard,          na.rm = TRUE),
    cv_n_fold     = sd(n_fold) / mean(n_fold)
  ), by = lfc_cutoff][order(lfc_cutoff)]

  cat("\n--- Across-fold stability vs LFC cutoff ---\n")
  print(stab)

  # Empirical CV minimum: across-fold CV of DEG counts is U-shaped (decreasing
  # as LFC tightens, then climbing past the plateau due to small-N fold
  # variance). We mark the lowest LFC within 5% of the minimum CV. This is an
  # empirical readout, NOT a pre-selected canonical cutoff -- the full sweep is
  # shown so the plateau can be read directly.
  cv_min    <- min(stab$cv_n_fold, na.rm = TRUE)
  cv_thresh <- cv_min * 1.05
  elbow_lfc <- stab[cv_n_fold <= cv_thresh, min(lfc_cutoff)]
  message(sprintf(
    "[caption %s] Empirical CV minimum of N DEGs = %.1f%% at |%s| = %.2g; ",
    sc$suffix, 100 * cv_min, eff, stab[which.min(cv_n_fold), lfc_cutoff]),
    sprintf("lowest cutoff within 5%% of min CV = %.2g (%s).",
            elbow_lfc, sc$sig_label))
  cat(sprintf("\nMin across-fold CV(N DEGs) = %.3f at |LFC| = %.2g\n",
              cv_min, stab[which.min(cv_n_fold), lfc_cutoff]))
  cat(sprintf("Empirical CV min (lowest LFC within 5%% of min CV): |LFC| = %.2g\n",
              elbow_lfc))

  # Panel A: heatmap (folds x LFC), text = N DEGs
  sweep[, fold_label := sprintf("%s (%.1f%%)", author, pct_held)]
  fold_order <- sweep[lfc_cutoff == 0][order(pct_held), fold_label]
  sweep[, fold_label := factor(fold_label, levels = fold_order)]
  sweep[, lfc_label := factor(sprintf("%.2g", lfc_cutoff),
                              levels = sprintf("%.2g", LFC_GRID))]

  elbow_xnum <- which(levels(sweep$lfc_label) == sprintf("%.2g", elbow_lfc))
  n_folds <- length(levels(sweep$fold_label))

  p_hm <- ggplot(sweep, aes(x = lfc_label, y = fold_label, fill = n_fold)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = format(n_fold, big.mark = ",")),
              size = GEOM_TEXT_6PT, color = "black") +
    annotate("rect",
             xmin = elbow_xnum - 0.5, xmax = elbow_xnum + 0.5,
             ymin = 0.5,             ymax = n_folds + 0.5,
             color = "#FFB300", fill = NA, linewidth = 0.9) +
    scale_fill_gradient(low = "#F5F9FB", high = "#9CC2D6", trans = "log10",
                        breaks = c(300, 1000, 3000, 10000),
                        labels = scales::label_comma(),
                        name = "N DEGs\nin fold") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL,
         y = NULL) +
    theme_masld() + theme_pub() +
    theme(panel.grid    = element_blank(),
          axis.ticks    = element_blank(),
          axis.text     = element_text(size = PUB_AXIS_TEXT, color = "black"),
          legend.title  = element_text(size = PUB_LEGEND_TIT, face = "plain",
                                       color = "black"),
          legend.text   = element_text(size = PUB_LEGEND, color = "black"))

  # Panel B: across-fold CV-of-DEG-count curve (sharp empirical minimum)
  stab[, lfc_label := factor(sprintf("%.2g", lfc_cutoff),
                             levels = sprintf("%.2g", LFC_GRID))]
  elbow_xpos <- which(levels(stab$lfc_label) == sprintf("%.2g", elbow_lfc))

  p_curve <- ggplot(stab, aes(x = lfc_label, y = 100 * cv_n_fold, group = 1)) +
    geom_vline(xintercept = elbow_xpos,
               linetype = "dashed", color = "#FFB300", linewidth = 0.6) +
    geom_line(color = "#37474F", linewidth = 0.7) +
    geom_point(color = "#37474F", size = 1.8) +
    scale_y_continuous(labels = function(v) paste0(v, "%")) +
    labs(x = expression("|log"[2]*"FC| cutoff"),
         y = "Coefficient of variation\nof N DEGs") +
    theme_masld() + theme_pub() +
    theme(axis.text  = element_text(size = PUB_AXIS_TEXT, color = "black"),
          axis.title = element_text(size = PUB_AXIS_TITLE, color = "black"))

  # Compose (heatmap + stability curve, stacked)
  p_g <- p_hm / p_curve + plot_layout(heights = c(2.4, 1))

  ggsave(OUT_PDF, p_g, width = fig_col_width, height = 3.0, device = cairo_pdf)
  saveRDS(p_g, OUT_RDS)
  cat("Saved: ", OUT_PDF, "\n", sep = "")
  cat("Saved: ", OUT_RDS, "\n", sep = "")

  # Headline cells at empirical CV minimum
  cat(sprintf("\n--- LOO recovery at empirical CV min |LFC| = %.2g ---\n", elbow_lfc))
  print(sweep[lfc_cutoff == elbow_lfc,
              .(author, pct_held, n_fold, n_full,
                pct_recovered = round(pct_recovered, 1),
                jaccard = round(jaccard, 3))][order(pct_held)])
  invisible(NULL)
}

for (nm in names(SCALES)) run_scale(SCALES[[nm]])

cat("\nDone.\n")
