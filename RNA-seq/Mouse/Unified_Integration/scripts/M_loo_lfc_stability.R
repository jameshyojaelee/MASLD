#!/usr/bin/env Rscript
##############################################################################
# M_loo_lfc_stability.R  (new 2026-04-27)
#
# Mouse analogue of human Fig 1G LOO-CV |log2FC| stability analysis.
#
# For each leave-one-diet-model-out (LOO) refit of the pooled mouse dream
# mega-analysis (5 folds: MCD, HFD, CDAHFD, FPC, LIDPAD) we count DEGs at
# padj<0.05 across a grid of |log2FC| cutoffs and measure recovery /
# Jaccard / across-fold CV vs the full pooled model. The cutoff at which
# across-fold CV(N DEGs) reaches a stable plateau is the empirically
# defensible mouse threshold.
#
# Mouse fold-changes are typically much larger than human (single-diet
# induction with strict controls), so the LFC grid spans a wider range
# than the human (0-1) sweep.
#
# Outputs:
#   results/loo_lfc_stability_data.csv             (full sweep table)
#   figures/supplementary/figS06_cross_species/mouse_loo_lfc_stability.pdf
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(viridis)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration")
RDIR <- file.path(INT, "results")
LOO_DIR <- file.path(RDIR, "integration", "loo_cv")

# Reuse human figure styling if available (graceful fallback)
theme_path <- file.path(BASE, "scripts/figures/publication_theme.R")
load_path  <- file.path(BASE, "scripts/figures/load_figure_data.R")
if (file.exists(theme_path)) source(theme_path)
if (file.exists(load_path))  source(load_path)
if (!exists("theme_masld")) theme_masld <- function() theme_minimal(base_size = 9)
if (!exists("FIGS06_DIR")) {
  FIGS06_DIR <- file.path(BASE, "figures/supplementary/figS06_cross_species")
}

dir.create(FIGS06_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PDF <- file.path(FIGS06_DIR, "mouse_loo_lfc_stability.pdf")
OUT_CSV <- file.path(RDIR, "loo_lfc_stability_data.csv")

# Diet labels (held-out unit; analogous to human first-author cohort labels)
DIET_LABELS <- c(MCD = "MCD", HFD = "HFD", CDAHFD = "CDAHFD",
                 FPC = "FPC", LIDPAD = "LIDPAD")
diet_models <- names(DIET_LABELS)

PADJ_CUTOFF <- 0.05
# Mouse LFC grid: spans a wider range than human (0-1) since single-diet
# induction with strict controls produces larger fold changes.
LFC_GRID <- c(0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5)

# -----------------------------------------------------------------------------
# Cohort sizes -> % samples held out per fold
# (analogous to human "% patients held out")
# -----------------------------------------------------------------------------
qc <- fread(file.path(INT, "qc/sample_qc_report.csv"))
meta <- readRDS(file.path(RDIR, "meta_matched.rds"))
pass <- qc[pass_qc == TRUE, sample_id]
meta_pass <- meta[sample_id %in% pass]

# All samples that would be removed when holding out a diet:
# the diet's Disease samples + Controls drawn from the SAME datasets.
# Because mouse pooled dream uses all controls + all disease at once,
# "held out" for fold X = samples with diet_model == X.
n_per <- meta_pass[, .N, by = diet_model][diet_model %in% diet_models]
n_total <- nrow(meta_pass)
n_per[, pct_held := round(100 * N / n_total, 1)]
setkey(n_per, diet_model)

cat(sprintf("Total mega-eligible samples: %d across %d diet models\n",
            n_total, length(diet_models)))
print(n_per)

# -----------------------------------------------------------------------------
# Load full dream model + LOO refits
# -----------------------------------------------------------------------------
full_dream <- fread(file.path(RDIR, "meta_analysis/dream_pooled_results.csv"))
setnames(full_dream, "adj.P.Val", "padj", skip_absent = TRUE)

# Identify which folds are available on disk
loo_files <- file.path(LOO_DIR, sprintf("dream_loo_%s.csv", diet_models))
have <- file.exists(loo_files)
if (!any(have)) {
  stop("No LOO refit files found in ", LOO_DIR,
       "\n  Run sbatch RNA-seq/Mouse/Unified_Integration/scripts/run_mouse_loo_dream.sbatch first.")
}
if (!all(have)) {
  cat("\nWARNING: Missing LOO folds: ",
      paste(diet_models[!have], collapse = ", "), "\n",
      "Stability analysis will use only the ", sum(have), " folds present.\n",
      sep = "")
  diet_models <- diet_models[have]
  loo_files   <- loo_files[have]
}

loo <- rbindlist(lapply(seq_along(diet_models), function(i) {
  ds <- diet_models[i]
  dt <- fread(loo_files[i])
  setnames(dt, "adj.P.Val", "padj", skip_absent = TRUE)
  dt[, held_out := ds]
  dt
}))

# -----------------------------------------------------------------------------
# Sweep: per fold x LFC cutoff -> #DEGs and recovery vs full pooled
# -----------------------------------------------------------------------------
sweep <- rbindlist(lapply(LFC_GRID, function(cut) {
  full_set <- full_dream[padj < PADJ_CUTOFF & abs(logFC) > cut, gene]
  rbindlist(lapply(diet_models, function(ds) {
    fold <- loo[held_out == ds]
    fold_set <- fold[padj < PADJ_CUTOFF & abs(logFC) > cut, gene]
    inter <- length(intersect(fold_set, full_set))
    union <- length(union(fold_set, full_set))
    data.table(
      held_out      = ds,
      label         = DIET_LABELS[ds],
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
cat("\nSaved: ", OUT_CSV, "\n", sep = "")

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

cat("\n--- Across-fold stability vs LFC cutoff (mouse) ---\n")
print(stab)

# Stability elbow: smallest LFC achieving CV within 5% of the minimum.
cv_min   <- min(stab$cv_n_fold, na.rm = TRUE)
cv_thr   <- cv_min * 1.05
elbow_lfc <- stab[cv_n_fold <= cv_thr, min(lfc_cutoff)]
plateau  <- stab[cv_n_fold <= cv_thr, range(lfc_cutoff)]

cat(sprintf("\nMin across-fold CV(N DEGs) = %.3f at |LFC| = %.2g\n",
            cv_min, stab[which.min(cv_n_fold), lfc_cutoff]))
cat(sprintf("Stability elbow (lowest LFC within 5%% of min CV): |LFC| = %.2g\n",
            elbow_lfc))
cat(sprintf("Stability plateau (all LFC within 5%% of min CV): |LFC| = [%.2g, %.2g]\n",
            plateau[1], plateau[2]))

# -----------------------------------------------------------------------------
# Panel A: heatmap (fold x LFC), fill = pct_recovered, text = N DEGs
# -----------------------------------------------------------------------------
sweep[, fold_label := sprintf("%s (%.1f%%)", label, pct_held)]
fold_order <- sweep[lfc_cutoff == LFC_GRID[1]][order(pct_held), fold_label]
sweep[, fold_label := factor(fold_label, levels = fold_order)]
sweep[, lfc_label := factor(sprintf("%.2g", lfc_cutoff),
                            levels = sprintf("%.2g", LFC_GRID))]

p_hm <- ggplot(sweep, aes(x = lfc_label, y = fold_label, fill = pct_recovered)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = format(n_fold, big.mark = ",")),
            size = 2.5, color = "white") +
  geom_tile(data = sweep[lfc_cutoff == elbow_lfc],
            color = "#FFB300", fill = NA, linewidth = 0.9) +
  scale_fill_viridis(option = "mako", limits = c(50, 100),
                     oob = scales::squish,
                     name = "% of full-model\nDEGs recovered") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "Mouse LOO-CV stability of integrated DEGs vs |log2FC| cutoff",
       subtitle = sprintf(
         "Cell label = N DEGs in fold (padj < %.2f) · row = held-out diet model (%% samples held out) · auto-elbow at |LFC| = %.2g",
         PADJ_CUTOFF, elbow_lfc)) +
  theme_masld() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.subtitle = element_text(size = 7.5, color = "gray40"))

# -----------------------------------------------------------------------------
# Panel B: across-fold stability curve (CV of DEG counts)
# -----------------------------------------------------------------------------
stab[, lfc_label := factor(sprintf("%.2g", lfc_cutoff),
                           levels = sprintf("%.2g", LFC_GRID))]
elbow_xpos <- which(levels(stab$lfc_label) == sprintf("%.2g", elbow_lfc))

p_curve <- ggplot(stab, aes(x = lfc_label, y = 100 * cv_n_fold, group = 1)) +
  geom_line(color = "#37474F", linewidth = 0.7) +
  geom_point(color = "#37474F", size = 1.8) +
  geom_vline(xintercept = elbow_xpos,
             linetype = "dashed", color = "#FFB300", linewidth = 0.6) +
  annotate("text",
           x = elbow_xpos,
           y = max(100 * stab$cv_n_fold, na.rm = TRUE) * 0.95,
           label = sprintf("stability elbow |LFC| = %.2g", elbow_lfc),
           hjust = -0.05, vjust = 1, size = 2.6, color = "#B26A00") +
  scale_y_continuous(labels = function(v) paste0(v, "%")) +
  labs(x = expression("|log"[2]*"FC| cutoff (mouse)"),
       y = "Across-fold CV\nof N DEGs",
       caption = sprintf(
         "CV = sd(N DEGs) / mean(N DEGs) across %d LOO folds (held-out diet model). Mean recovery range %.1f%%-%.1f%%.",
         length(diet_models),
         min(stab$mean_recovery), max(stab$mean_recovery))) +
  theme_masld() +
  theme(plot.caption = element_text(size = 6.5, color = "gray40", hjust = 0))

# Compose stacked panels
p_g <- p_hm / p_curve + plot_layout(heights = c(2.4, 1))
ggsave(OUT_PDF, p_g, width = 7, height = 5.5, device = cairo_pdf)
cat("Saved: ", OUT_PDF, "\n", sep = "")

# Headline cells at elbow
cat(sprintf("\n--- LOO recovery at elbow |LFC| = %.2g ---\n", elbow_lfc))
print(sweep[lfc_cutoff == elbow_lfc,
            .(label, pct_held, n_fold, n_full,
              pct_recovered = round(pct_recovered, 1),
              jaccard       = round(jaccard, 3))][order(pct_held)])

cat("\nDone.\n")
