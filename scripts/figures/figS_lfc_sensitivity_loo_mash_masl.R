#!/usr/bin/env Rscript
##############################################################################
# figS_lfc_sensitivity_loo_mash_masl.R  (2026-05-13)
#
# Panel-F-style LOO-CV stability plot for the 5 new contrasts:
#   mash_vs_masl primary | mash_vs_masl strict
#   mash_vs_healthy primary | mash_vs_healthy strict
#   masl_vs_healthy
#
# Mirrors the methodology of figS_lfc_sensitivity_loo.R (which justified
# |LFC|>0.5 for the legacy Disease-vs-Control Tier 1):
#   - heatmap (held-out cohort x |LFC| threshold), fill+text = N DEGs in fold
#   - across-fold CV(N DEGs) curve, elbow detected as smallest |LFC| with
#     CV within 5% of the global minimum CV
#
# Run via:
#   Rscript figS_lfc_sensitivity_loo_mash_masl.R
# Optionally restrict to one contrast:
#   CONTRAST_TAG=mash_vs_masl Rscript figS_lfc_sensitivity_loo_mash_masl.R
#
# Outputs (per contrast tag):
#   figures/supplementary/figS_methods_validation/lfc_sensitivity/panels/F_loo_cv_stability_<tag>.pdf
#   figures/supplementary/figS_methods_validation/lfc_sensitivity/loo_cv_stability_<tag>_data.csv
##############################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(viridis)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

INT_DIR  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
DSIG_DIR <- file.path(INT_DIR, "results/disease_signatures")
LOO_ROOT <- file.path(INT_DIR, "results/integration/loo_cv")
OUT_DIR  <- file.path(FIG_SUPP, "figS_methods_validation/lfc_sensitivity")
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

STUDY_NAMES <- c(
  GSE126848   = "GSE126848", GSE130970 = "GSE130970", GSE135251 = "GSE135251",
  GSE162694   = "GSE162694", GSE174478 = "GSE174478", GSE193066 = "GSE193066",
  GSE213621   = "GSE213621", GSE240729 = "GSE240729",
  GSE167523   = "GSE167523",PRJNA512027 = "PRJNA512027")

PADJ_CUTOFF <- 0.05
LFC_GRID    <- c(0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.75, 1.0)
TIER1_LFC   <- 0.5  # uniform Tier 1 cutoff across all contrasts (decision 2026-05-13)

# Contrast registry --------------------------------------------------------
# Canonical MASH definition (2026-05-13): strict NAS >= 5 only.
# Primary (Borderline-grouped) variants kept as data CSVs but figure regen skips.
CONTRASTS <- list(
  mash_vs_masl_strict = list(
    full_csv  = file.path(DSIG_DIR, "mash_vs_masl_dream_strict.csv"),
    loo_dir   = file.path(LOO_ROOT, "mash_vs_masl_strict"),
    cohorts   = c("GSE126848","GSE130970","GSE135251","GSE162694",
                  "GSE167523","GSE174478","GSE193066"),
    title     = "LOOCV stability — MASH vs MASL"),
  mash_vs_healthy_strict = list(
    full_csv  = file.path(DSIG_DIR, "mash_vs_healthy_dream_strict.csv"),
    loo_dir   = file.path(LOO_ROOT, "mash_vs_healthy_strict"),
    cohorts   = c("GSE126848","GSE130970","GSE135251","GSE162694"),
    title     = "LOOCV stability — MASH vs Healthy"),
  masl_vs_healthy = list(
    full_csv  = file.path(DSIG_DIR, "masl_vs_healthy_dream.csv"),
    loo_dir   = file.path(LOO_ROOT, "masl_vs_healthy"),
    cohorts   = c("GSE126848","GSE130970","GSE135251","GSE162694"),
    title     = "LOOCV stability — MASL vs Healthy")
)

# Cohort sample counts (post pass_technical) -------------------------------
qc <- fread(file.path(INT_DIR, "qc/sample_qc_report.csv"))
meta <- as.data.table(readRDS(file.path(INT_DIR, "results/integration/meta_matched.rds")))
meta <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

# Per-contrast samples = cohort rows with either arm of the contrast
contrast_n_samples <- function(cohort, tag) {
  # diagnosis filter mirrors what the LOO/contrast script uses
  if (grepl("mash_vs_masl", tag)) {
    keep <- if (grepl("strict", tag)) c("NAFL","NASH") else c("NAFL","NASH","Borderline")
  } else if (grepl("mash_vs_healthy", tag)) {
    keep <- if (grepl("strict", tag)) c("Control","NASH") else c("Control","NASH","Borderline")
  } else { # masl_vs_healthy
    keep <- c("Control","NAFL")
  }
  meta[dataset == cohort & diagnosis_harmonized %in% keep, .N]
}

# Plot one contrast --------------------------------------------------------
plot_contrast <- function(tag) {
  spec <- CONTRASTS[[tag]]
  cat(sprintf("\n=== %s ===\n", tag))
  if (!file.exists(spec$full_csv))
    stop("Missing full dream csv: ", spec$full_csv)

  full <- fread(spec$full_csv)
  setnames(full, "adj.P.Val", "padj", skip_absent = TRUE)

  loo_list <- lapply(spec$cohorts, function(ds) {
    p <- file.path(spec$loo_dir, sprintf("dream_loo_%s.csv", ds))
    if (!file.exists(p)) { warning("Missing LOO fold: ", p); return(NULL) }
    dt <- fread(p); setnames(dt, "adj.P.Val", "padj", skip_absent = TRUE)
    dt[, held_out := ds]
    dt
  })
  loo <- rbindlist(loo_list[!sapply(loo_list, is.null)])

  cohorts_present <- unique(loo$held_out)
  n_per <- data.table(dataset = cohorts_present,
                      N = sapply(cohorts_present, contrast_n_samples, tag = tag))
  n_per[, pct_held := round(100 * N / sum(N), 1)]
  setkey(n_per, dataset)
  cat(sprintf("  Cohorts: %d, total samples: %d\n", nrow(n_per), sum(n_per$N)))

  # Sweep n_fold and recovery
  sweep <- rbindlist(lapply(LFC_GRID, function(cut) {
    full_set <- full[padj < PADJ_CUTOFF & abs(logFC) > cut, gene]
    rbindlist(lapply(cohorts_present, function(ds) {
      fold <- loo[held_out == ds]
      fold_set <- fold[padj < PADJ_CUTOFF & abs(logFC) > cut, gene]
      inter <- length(intersect(fold_set, full_set))
      union <- length(union(fold_set, full_set))
      data.table(
        held_out   = ds,
        author     = STUDY_NAMES[ds],
        pct_held   = n_per[ds, pct_held],
        lfc_cutoff = cut,
        n_full     = length(full_set),
        n_fold     = length(fold_set),
        pct_recovered = if (length(full_set)) 100 * inter / length(full_set) else NA_real_,
        jaccard    = if (union) inter / union else NA_real_)
    }))
  }))

  out_csv <- file.path(OUT_DIR, sprintf("loo_cv_stability_%s_data.csv", tag))
  fwrite(sweep, out_csv); cat(sprintf("  Saved: %s\n", out_csv))

  # Across-fold CV stability (computed for sweep curve + caption stability evidence)
  stab <- sweep[, .(cv_n_fold = sd(n_fold) / mean(n_fold)), by = lfc_cutoff][order(lfc_cutoff)]
  cv_min       <- min(stab$cv_n_fold, na.rm = TRUE)
  cv_min_lfc   <- stab[which.min(cv_n_fold), lfc_cutoff]
  cv_at_tier1  <- stab[lfc_cutoff == TIER1_LFC, cv_n_fold]
  n_at_tier1   <- sweep[lfc_cutoff == TIER1_LFC, unique(n_full)]
  cat(sprintf("  min CV = %.4f at LFC=%.2g | CV at uniform Tier 1 |LFC|>%.2g = %.4f (n_full=%d)\n",
              cv_min, cv_min_lfc, TIER1_LFC, cv_at_tier1, n_at_tier1))

  # Plot panels
  sweep[, fold_label := sprintf("%s (%.1f%%)", author, pct_held)]
  fold_order <- sweep[lfc_cutoff == 0][order(pct_held), fold_label]
  sweep[, fold_label := factor(fold_label, levels = fold_order)]
  sweep[, lfc_label  := factor(sprintf("%.2g", lfc_cutoff),
                               levels = sprintf("%.2g", LFC_GRID))]
  # Highlight the canonical uniform Tier 1 cutoff |LFC|>0.5 (not the elbow).
  tier1_x  <- which(levels(sweep$lfc_label) == sprintf("%.2g", TIER1_LFC))
  n_folds  <- length(levels(sweep$fold_label))

  p_hm <- ggplot(sweep, aes(x = lfc_label, y = fold_label, fill = n_fold)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = format(n_fold, big.mark = ",")),
              size = 2.0, color = "gray15") +
    annotate("rect",
             xmin = tier1_x - 0.5, xmax = tier1_x + 0.5,
             ymin = 0.5, ymax = n_folds + 0.5,
             color = "#FFB300", fill = NA, linewidth = 0.9) +
    scale_fill_gradient(low = "#F5F9FB", high = "#9CC2D6", trans = "log10",
                        labels = scales::label_comma(),
                        name = "N DEGs\nin fold") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL,
         title = spec$title,
         subtitle = sprintf("N DEGs in LOO fold (padj < %.2f); Tier 1 = |LFC| > %.2g (uniform across contrasts)",
                            PADJ_CUTOFF, TIER1_LFC)) +
    theme_masld() + theme_pub() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = PUB_TITLE + 1, face = "bold"),
          plot.subtitle = element_text(size = PUB_SUBTITLE + 1, color = "gray30"),
          axis.text = element_text(size = PUB_AXIS_TEXT, color = "black"),
          legend.title = element_text(size = PUB_LEGEND_TIT, face = "bold"),
          legend.text  = element_text(size = PUB_LEGEND))

  stab[, lfc_label := factor(sprintf("%.2g", lfc_cutoff),
                             levels = sprintf("%.2g", LFC_GRID))]
  tier1_x_curve <- which(levels(stab$lfc_label) == sprintf("%.2g", TIER1_LFC))
  p_curve <- ggplot(stab, aes(x = lfc_label, y = 100 * cv_n_fold, group = 1)) +
    geom_vline(xintercept = tier1_x_curve,
               linetype = "dashed", color = "#FFB300", linewidth = 0.6) +
    geom_hline(yintercept = 100 * cv_min,
               linetype = "dotted", color = "gray50", linewidth = 0.4) +
    geom_line(color = "#37474F", linewidth = 0.7) +
    geom_point(color = "#37474F", size = 1.8) +
    scale_y_continuous(labels = function(v) paste0(v, "%")) +
    labs(x = expression("|log"[2]*"FC| cutoff"),
         y = "CV of N DEGs across folds",
         caption = sprintf("Min CV = %.3f at |LFC|>%.2g (dotted); uniform Tier 1 = |LFC|>%.2g (dashed gold).",
                           cv_min, cv_min_lfc, TIER1_LFC)) +
    theme_masld() + theme_pub() +
    theme(axis.text  = element_text(size = PUB_AXIS_TEXT + 1, color = "black"),
          axis.title = element_text(size = PUB_AXIS_TITLE + 1),
          plot.caption = element_text(size = PUB_LEGEND, color = "gray40"))

  p <- p_hm / p_curve + plot_layout(heights = c(2.4, 1))
  out_pdf <- file.path(PANEL_DIR, sprintf("F_loo_cv_stability_%s.pdf", tag))
  ggsave(out_pdf, p, width = fig_col_width, height = 3.4, device = cairo_pdf)
  cat(sprintf("  Saved: %s\n", out_pdf))
}

# Driver -------------------------------------------------------------------
tag_filter <- Sys.getenv("CONTRAST_TAG", "")
tags <- if (nchar(tag_filter)) tag_filter else names(CONTRASTS)
for (t in tags) plot_contrast(t)
cat("\nDone.\n")
