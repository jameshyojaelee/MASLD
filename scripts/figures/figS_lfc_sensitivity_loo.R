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
# Output:
#   figures/supplementary/figS_methods_validation/lfc_sensitivity/panels/F_loo_cv_stability.pdf
#   figures/supplementary/figS_methods_validation/lfc_sensitivity/loo_cv_stability_data.csv
#   figures/supplementary/figS_methods_validation/lfc_sensitivity/F_loo_cv_stability_pG.rds
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
LOO_DIR   <- file.path(INT_DIR, "loo_cv")
OUT_DIR   <- file.path(FIG_SUPP, "figS_methods_validation/lfc_sensitivity")
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PDF <- file.path(PANEL_DIR, "F_loo_cv_stability.pdf")
OUT_RDS <- file.path(OUT_DIR,   "F_loo_cv_stability_pG.rds")
OUT_CSV <- file.path(OUT_DIR,   "loo_cv_stability_data.csv")

# Cohort first-author labels (full atlas)
ALL_STUDY_NAMES <- c(
  GSE126848 = "Suppli", GSE130970 = "Hoang", GSE135251 = "Govaere",
  GSE162694 = "Bril",   GSE174478 = "Kawamura", GSE193066 = "Hoshida",
  GSE213621 = "Chen",   GSE240729 = "Verschuren"
)
# Mega-analysis cohort set is the canonical source of truth in
# config/human_datasets.yaml (`include_in_mega: true`).
ycfg_path <- file.path(BASE, "config/human_datasets.yaml")
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega),
                             yaml::read_yaml(ycfg_path)$datasets))
STUDY_NAMES <- ALL_STUDY_NAMES[intersect(names(ALL_STUDY_NAMES), mega_cohorts)]
stopifnot("All mega cohorts must have first-author labels" =
            length(STUDY_NAMES) == length(mega_cohorts))

PADJ_CUTOFF <- 0.05  # primary integrated DEG significance threshold
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
  dt <- fread(file.path(LOO_DIR, sprintf("dream_loo_%s.csv", ds)))
  setnames(dt, "adj.P.Val", "padj", skip_absent = TRUE)
  dt[, held_out := ds]
  dt
}))

# -----------------------------------------------------------------------------
# Sweep: per fold x LFC cutoff -> #DEGs and recovery vs full
# -----------------------------------------------------------------------------
sweep <- rbindlist(lapply(LFC_GRID, function(cut) {
  full_set <- full_dream[padj < PADJ_CUTOFF & abs(logFC) > cut, gene]
  rbindlist(lapply(mega_cohorts, function(ds) {
    fold <- loo[held_out == ds]
    fold_set <- fold[padj < PADJ_CUTOFF & abs(logFC) > cut, gene]
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

# -----------------------------------------------------------------------------
# Across-fold stability metric -> elbow detection
# -----------------------------------------------------------------------------
stab <- sweep[, .(
  mean_recovery = mean(pct_recovered, na.rm = TRUE),
  sd_recovery   = sd(pct_recovered,   na.rm = TRUE),
  mean_jaccard  = mean(jaccard,        na.rm = TRUE),
  sd_jaccard    = sd(jaccard,          na.rm = TRUE),
  cv_n_fold     = sd(n_fold) / mean(n_fold)
), by = lfc_cutoff][order(lfc_cutoff)]

cat("\n--- Across-fold stability vs LFC cutoff ---\n")
print(stab)

# Stability elbow: across-fold CV of DEG counts has a U-shape (decreasing as
# LFC tightens, then climbing again past ~0.5 due to small-N fold variance).
# The "smallest LFC achieving near-minimum CV" is the empirically defensible
# cutoff -- the lowest stringency at which LOO refits all return comparable
# DEG sets.
cv_min <- min(stab$cv_n_fold, na.rm = TRUE)
cv_thresh <- cv_min * 1.05
elbow_lfc <- stab[cv_n_fold <= cv_thresh, min(lfc_cutoff)]
cat(sprintf("\nMin across-fold CV(N DEGs) = %.3f at |LFC| = %.2g\n",
            cv_min, stab[which.min(cv_n_fold), lfc_cutoff]))
cat(sprintf("Stability elbow (lowest LFC within 5%% of min CV): |LFC| = %.2g\n",
            elbow_lfc))

# -----------------------------------------------------------------------------
# Panel A: heatmap (folds x LFC), fill = pct_recovered, text = N DEGs
# -----------------------------------------------------------------------------
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
            size = 2.2, color = "gray15") +
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
       y = NULL,
       title = "LOOCV stability vs |log2FC| cutoff",
       subtitle = sprintf("N DEGs in fold (padj < %.2f)", PADJ_CUTOFF)) +
  theme_masld() + theme_pub() +
  theme(panel.grid    = element_blank(),
        axis.ticks    = element_blank(),
        plot.title    = element_text(size = PUB_TITLE + 1, face = "bold"),
        plot.subtitle = element_text(size = PUB_SUBTITLE + 1, color = "gray30"),
        axis.text     = element_text(size = PUB_AXIS_TEXT + 1, color = "black"),
        legend.title  = element_text(size = PUB_LEGEND_TIT + 1, face = "bold"),
        legend.text   = element_text(size = PUB_LEGEND + 1))

# -----------------------------------------------------------------------------
# Panel B: across-fold stability curve -- CV of DEG counts (this is the
# metric that actually has a sharp elbow; mean recovery saturates flat).
# -----------------------------------------------------------------------------
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
  theme(axis.text  = element_text(size = PUB_AXIS_TEXT + 1, color = "black"),
        axis.title = element_text(size = PUB_AXIS_TITLE + 1))

# -----------------------------------------------------------------------------
# Compose panel G (heatmap + stability curve, stacked)
# -----------------------------------------------------------------------------
p_g <- p_hm / p_curve + plot_layout(heights = c(2.4, 1))

ggsave(OUT_PDF, p_g, width = fig_col_width, height = 3.0, device = cairo_pdf)
saveRDS(p_g, OUT_RDS)
cat("Saved: ", OUT_PDF, "\n", sep = "")
cat("Saved: ", OUT_RDS, "\n", sep = "")

# Headline cells at elbow
cat(sprintf("\n--- LOO recovery at elbow |LFC| = %.2g ---\n", elbow_lfc))
print(sweep[lfc_cutoff == elbow_lfc,
            .(author, pct_held, n_fold, n_full, pct_recovered = round(pct_recovered, 1),
              jaccard = round(jaccard, 3))][order(pct_held)])

cat("\nDone.\n")
