#!/usr/bin/env Rscript
# KEY MESSAGE: The integrated DEG set is stable to leave-one-cohort-out — refitting
# the pooled model on 4 cohorts recovers a mean 81.6% of the full DEGs (Jaccard 0.74),
# effect sizes correlate at Spearman rho 0.95, and direction is ~99.9% concordant.
# The only weak fold is the largest cohort (Chen, 42.3% of samples) held out.
# ============================================================================
# loo_cv_stability.R  — Figure 3 LOO-CV stability panel
#   full 5-cohort integrated fit vs each leave-one-cohort-out 4-cohort refit.
#   Source: .../results/integration/loo_cv_C2/loo_cv_C2_summary.csv (canonical C2).
# Output: FIG2_DIR/panels/figs3c_loo_cv_stability.pdf  (FIG2_DIR = figures/main/fig3_RNAseq)
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

s <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/loo_cv_C2/loo_cv_C2_summary.csv"))

# Held-out cohort -> first-author label + fraction of mega-analysis samples held out
auth    <- c(GSE126848 = "GSE126848",  GSE130970 = "GSE130970", GSE135251 = "GSE135251",
             GSE162694 = "GSE162694",    GSE213621 = "GSE213621")
pctheld <- c(GSE126848 = 6.5, GSE130970 = 9.0, GSE135251 = 25.4,
             GSE162694 = 16.8, GSE213621 = 42.3)
s[, author   := auth[held_out_cohort]]
s[, pct_held := pctheld[held_out_cohort]]
s[, label    := sprintf("%s (%.0f%%)", author, pct_held)]

# Long form: the four full-vs-leave-one-out stability metrics the manuscript cites.
RHO <- "Spearman ρ"
m <- rbindlist(list(
  s[, .(label, pct_held, metric = "Recovery (%)",  value = pct_full_recovered_01)],
  s[, .(label, pct_held, metric = "Jaccard",       value = full_vs_train_jaccard_01)],
  s[, .(label, pct_held, metric = RHO,             value = full_vs_train_spearman)],
  s[, .(label, pct_held, metric = "Direction (%)", value = full_vs_train_direction)]
))
metric_levels <- c("Recovery (%)", "Jaccard", RHO, "Direction (%)")
m[, metric := factor(metric, levels = metric_levels)]

# Order cohorts by recovery so the gradient (and the weak Chen fold) reads top->down
ord <- s[order(pct_full_recovered_01), label]
m[, label := factor(label, levels = ord)]

# Per-metric mean (the headline numbers) + a formatted label
means <- m[, .(mu = mean(value)), by = metric]
means[, lab := ifelse(metric %in% c("Recovery (%)", "Direction (%)"),
                      sprintf("mean %.1f", mu), sprintf("mean %.2f", mu))]
means[, ytop := length(ord) + 0.55]

p <- ggplot(m, aes(value, label)) +
  geom_vline(data = means, aes(xintercept = mu),
             linetype = "dashed", color = "gray55", linewidth = 0.3) +
  geom_point(aes(fill = pct_held), shape = 21, size = 2.6, stroke = 0.25, color = "black") +
  geom_text(data = means, aes(x = mu, y = ytop, label = lab),
            inherit.aes = FALSE, size = 2.3, fontface = "plain", color = "gray25") +
  facet_wrap(~ metric, nrow = 1, scales = "free_x") +
  scale_fill_gradient(low = masld_colors$masl, high = masld_colors$mash,
                      name = "% samples\nheld out") +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.0))) +
  labs(x = NULL, y = NULL,
       title = "Integrated DEGs are stable to leave-one-cohort-out refits") +
  theme_masld(base_size = 9) +
  theme(plot.title    = element_text(size = 10, face = "plain", margin = margin(b = 4)),
        strip.text    = element_text(size = 9, face = "plain"),
        axis.text.y   = element_text(size = 8.5),
        axis.text.x   = element_text(size = 7.5),
        panel.spacing = unit(0.5, "lines"),
        legend.position = "right",
        legend.key.height = unit(0.5, "cm"),
        legend.title  = element_text(size = 7.5),
        legend.text   = element_text(size = 7))

save_fig(p, file.path(PANEL_DIR, "figs3c_loo_cv_stability.pdf"),
         width = fig_full_width * 1.15, height = 2.4)

# Side table for caption transparency
fwrite(s[, .(held_out_cohort, author, pct_held,
             recovery_pct = pct_full_recovered_01,
             jaccard = full_vs_train_jaccard_01,
             spearman_rho = full_vs_train_spearman,
             direction_pct = full_vs_train_direction)],
       file.path(DATA_DIR, "loo_cv_stability.csv"))

cat("Wrote figs3c_loo_cv_stability.pdf\n")
cat(sprintf("means: recovery=%.1f%%  jaccard=%.3f  rho=%.3f  direction=%.1f%%\n",
            means[metric == "Recovery (%)", mu], means[metric == "Jaccard", mu],
            means[metric == RHO, mu], means[metric == "Direction (%)", mu]))
