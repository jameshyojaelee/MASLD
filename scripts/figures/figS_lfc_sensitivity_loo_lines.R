#!/usr/bin/env Rscript
##############################################################################
# figS_lfc_sensitivity_loo_lines.R  (new 2026-05-07)
# Supp panel G: line/dot version of the LOOCV stability plot. Same sweep as
# panel F but the per-cohort N DEGs are drawn as colored lines+points across
# |log2FC| cutoffs (one color per held-out cohort) instead of a heatmap.
#
# Bottom panel (across-fold CV of DEG counts) is unchanged from panel F.
#
# Output:
#   figures/supplementary/figS_lfc_sensitivity/panels/G_loo_cv_lines.pdf
#   figures/supplementary/figS_lfc_sensitivity/G_loo_cv_lines_pG.rds
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(yaml)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

INT_DIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
LOO_DIR   <- file.path(INT_DIR, "loo_cv")
OUT_DIR   <- file.path(FIG_SUPP, "figS_lfc_sensitivity")
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PDF <- file.path(PANEL_DIR, "G_loo_cv_lines.pdf")
OUT_RDS <- file.path(OUT_DIR,   "G_loo_cv_lines_pG.rds")

# Cohort first-author labels (full atlas)
ALL_STUDY_NAMES <- c(
  GSE126848 = "Suppli", GSE130970 = "Hoang", GSE135251 = "Govaere",
  GSE162694 = "Bril",   GSE174478 = "Kawamura", GSE193066 = "Hoshida",
  GSE213621 = "Chen",   GSE240729 = "Verschuren"
)
ycfg_path <- file.path(BASE, "config/human_datasets.yaml")
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega),
                             yaml::read_yaml(ycfg_path)$datasets))
STUDY_NAMES <- ALL_STUDY_NAMES[intersect(names(ALL_STUDY_NAMES), mega_cohorts)]
stopifnot("All mega cohorts must have first-author labels" =
            length(STUDY_NAMES) == length(mega_cohorts))

PADJ_CUTOFF <- 0.05
LFC_GRID <- c(0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.75, 1.0)

# -----------------------------------------------------------------------------
# Cohort sizes
# -----------------------------------------------------------------------------
qc <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv"))
mega_cohorts <- names(STUDY_NAMES)
n_per <- qc[pass_technical == TRUE & dataset %in% mega_cohorts, .N, by = dataset]
n_total <- sum(n_per$N)
n_per[, pct_held := round(100 * N / n_total, 1)]
setkey(n_per, dataset)

# -----------------------------------------------------------------------------
# Load full dream + LOO refits
# -----------------------------------------------------------------------------
full_dream <- fread(file.path(INT_DIR, "dream_results.csv"))
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

# -----------------------------------------------------------------------------
# Across-fold stability + elbow
# -----------------------------------------------------------------------------
stab <- sweep[, .(
  cv_n_fold = sd(n_fold) / mean(n_fold)
), by = lfc_cutoff][order(lfc_cutoff)]

cv_min   <- min(stab$cv_n_fold, na.rm = TRUE)
cv_thresh <- cv_min * 1.05
elbow_lfc <- stab[cv_n_fold <= cv_thresh, min(lfc_cutoff)]

# -----------------------------------------------------------------------------
# Top panel: per-cohort N DEGs vs |LFC| cutoff (line + points)
# -----------------------------------------------------------------------------
# Order cohorts by % held out (smallest -> largest) so the legend reads
# Suppli, Hoang, Bril, Govaere, Chen
fold_order <- unique(sweep[order(pct_held), .(author, pct_held)])
sweep[, fold_label := factor(sprintf("%s (%.1f%%)", author, pct_held),
                             levels = sprintf("%s (%.1f%%)",
                                              fold_order$author, fold_order$pct_held))]

# Colorblind-safe palette (Okabe-Ito), 5 entries
COHORT_COLS <- c("#0072B2", "#009E73", "#D55E00", "#CC79A7", "#E69F00")
names(COHORT_COLS) <- levels(sweep$fold_label)

p_top <- ggplot(sweep, aes(x = lfc_cutoff, y = n_fold,
                           color = fold_label, group = fold_label)) +
  geom_vline(xintercept = elbow_lfc,
             linetype = "dashed", color = "#FFB300", linewidth = 0.6) +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  geom_point(size = 2.0) +
  annotate("text", x = elbow_lfc,
           y = max(sweep$n_fold, na.rm = TRUE),
           label = sprintf("elbow |LFC| = %.2g", elbow_lfc),
           hjust = -0.05, vjust = 1.1, size = 2.6, color = "#B26A00") +
  scale_color_manual(values = COHORT_COLS, name = "Held-out cohort") +
  scale_x_continuous(breaks = LFC_GRID,
                     labels = sprintf("%.2g", LFC_GRID)) +
  scale_y_continuous(trans = "log10",
                     labels = label_comma(),
                     breaks = c(100, 300, 1000, 3000, 10000, 30000)) +
  labs(x = NULL,
       y = "N DEGs in LOO fold\n(padj < 0.05, log10)",
       title = "LOOCV DEG counts vs |log2FC| cutoff",
       subtitle = sprintf("One line per held-out cohort · auto-elbow at |LFC| = %.2g",
                          elbow_lfc)) +
  theme_masld() + theme_pub() +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(size = PUB_TITLE + 1, face = "bold"),
        plot.subtitle = element_text(size = PUB_SUBTITLE + 1, color = "gray30"),
        axis.text     = element_text(size = PUB_AXIS_TEXT + 1, color = "black"),
        axis.title    = element_text(size = PUB_AXIS_TITLE + 1),
        legend.title  = element_text(size = PUB_LEGEND_TIT + 1, face = "bold"),
        legend.text   = element_text(size = PUB_LEGEND + 1),
        legend.position = "right")

# -----------------------------------------------------------------------------
# Bottom panel: across-fold CV (unchanged from panel F)
# -----------------------------------------------------------------------------
p_curve <- ggplot(stab, aes(x = lfc_cutoff, y = 100 * cv_n_fold)) +
  geom_vline(xintercept = elbow_lfc,
             linetype = "dashed", color = "#FFB300", linewidth = 0.6) +
  geom_line(color = "#37474F", linewidth = 0.7) +
  geom_point(color = "#37474F", size = 1.8) +
  scale_x_continuous(breaks = LFC_GRID,
                     labels = sprintf("%.2g", LFC_GRID)) +
  scale_y_continuous(labels = function(v) paste0(v, "%")) +
  labs(x = expression("|log"[2]*"FC| cutoff"),
       y = "Coefficient of variation\nof N DEGs") +
  theme_masld() + theme_pub() +
  theme(axis.text  = element_text(size = PUB_AXIS_TEXT + 1, color = "black"),
        axis.title = element_text(size = PUB_AXIS_TITLE + 1))

# -----------------------------------------------------------------------------
# Compose
# -----------------------------------------------------------------------------
p_g <- p_top / p_curve + plot_layout(heights = c(2.4, 1))

ggsave(OUT_PDF, p_g, width = 7.0, height = 5.5, device = cairo_pdf)
saveRDS(p_g, OUT_RDS)
cat("Saved: ", OUT_PDF, "\n", sep = "")
cat("Saved: ", OUT_RDS, "\n", sep = "")
cat(sprintf("Auto-elbow at |LFC| = %.2g (min across-fold CV = %.2f%%)\n",
            elbow_lfc, 100 * cv_min))

cat("\nDone.\n")
