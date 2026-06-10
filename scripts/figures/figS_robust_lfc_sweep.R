#!/usr/bin/env Rscript
##############################################################################
# figS_robust_lfc_sweep.R  (rewritten 2026-05-08; bootstrap symmetric panels
# added; K=50 dropped per request)
#
# KEY MESSAGE: As |LFC| tightens, the surviving DEGs become more stable on
#              every resampling metric (CPSS pi-hat, bootstrap freq, K-fold
#              recurrence, per-iter rho), at the cost of fewer DEGs. The
#              canonical |LFC|>0.5 cutoff lands at 1,885 DEGs with 98.1%
#              passing CPSS pi-hat>=0.7 and 47.6% passing bootstrap freq>=0.9.
#
# Panels:
#   A — DEG count by |logFC| cutoff
#   D — K-fold-10 recurrence (dropped K=50 per PI request)
#   B — mean CPSS pi-hat
#   B'— mean bootstrap selection frequency
#   C — % canonical DEGs with CPSS pi-hat >= 0.7
#   C'— % canonical DEGs with bootstrap freq >= 0.9
#   E — per-iter Spearman rho vs full dream (CPSS, Bootstrap, K-fold-10)
##############################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
AUDIT   <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
OUT_DIR <- FIGS_ROBUST_DIR  # consolidated under figS_methods_validation/ (2026-06-04)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF <- file.path(OUT_DIR, "figS_robust_lfc_sweep.pdf")

theme_robust <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
    axis.line   = element_line(colour = "black", linewidth = 0.3),
    axis.ticks  = element_line(colour = "black", linewidth = 0.3),
    legend.position = "top",
    legend.title  = element_blank(),
    legend.margin = margin(b = -3),
    legend.text   = element_text(size = 9),
    plot.title    = element_text(face = "bold", size = 11,
                                 margin = margin(b = 2)),
    plot.subtitle = element_text(size = 9, colour = "grey30",
                                 margin = margin(b = 6)),
    plot.tag      = element_text(face = "bold", size = 12),
    plot.margin   = margin(8, 10, 8, 10),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 10)
  )

THR_LEVS    <- c("0", "0.1", "0.3", "0.5", "0.75", "1", "2")
BAR_FILL    <- "#4F8FBA"   # steel-blue for all non-canonical cutoffs
CUTOFF_FILL <- "#C96450"   # muted terracotta; matches Liang et al figure palette
THR_COLS    <- setNames(rep(BAR_FILL, length(THR_LEVS)), THR_LEVS)
THR_COLS["0.5"] <- CUTOFF_FILL

# -----------------------------------------------------------------------------
# Load
# -----------------------------------------------------------------------------
sweep <- fread(file.path(AUDIT, "pillar_A_lfc_threshold_sweep.csv"))
setorder(sweep, threshold)
sweep[, x_lab := factor(sprintf("%g", threshold), levels = THR_LEVS)]

iter <- fread(file.path(AUDIT, "pillar_A_lfc_sweep_iter_long.csv"))
iter <- iter[source != "kfold_K50"]   # drop K=50 per PI request
iter[, x_lab := factor(sprintf("%g", threshold), levels = THR_LEVS)]
iter[, source_lbl := factor(source,
       levels = c("cpss", "bootstrap", "kfold_K10"),
       labels = c("CPSS", "Bootstrap", "K-fold (K=10)"))]

# -----------------------------------------------------------------------------
# A — DEG count
# -----------------------------------------------------------------------------
p_a <- ggplot(sweep, aes(x = x_lab, y = n_DEG, fill = x_lab)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = comma(n_DEG)), vjust = -0.4, size = 3) +
  scale_fill_manual(values = THR_COLS, guide = "none") +
  scale_y_continuous(labels = comma,
                     expand = expansion(mult = c(0, 0.18))) +
  labs(tag = "A",
       title = "DEG count",
       subtitle = "padj < 0.05 + |log2FC| > x",
       x = "|log2FC| cutoff", y = "DEGs") +
  theme_robust

# -----------------------------------------------------------------------------
# D — K-fold-10 recurrence (K=50 dropped)
# -----------------------------------------------------------------------------
p_d <- ggplot(sweep, aes(x = x_lab, y = 100 * mean_kfold10 / 10,
                          fill = x_lab)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * mean_kfold10 / 10)),
            vjust = -0.4, size = 3) +
  scale_fill_manual(values = THR_COLS, guide = "none") +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25)) +
  labs(tag = "D",
       title = "K-fold (K=10) recurrence",
       subtitle = "Mean % of folds where each canonical DEG is re-selected",
       x = "|log2FC| cutoff",
       y = "% folds re-selected") +
  theme_robust

# -----------------------------------------------------------------------------
# B — mean CPSS pi-hat
# -----------------------------------------------------------------------------
p_b <- ggplot(sweep, aes(x = x_lab, y = mean_pi_hat, fill = x_lab)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", mean_pi_hat)),
            vjust = -0.4, size = 3) +
  scale_fill_manual(values = THR_COLS, guide = "none") +
  scale_y_continuous(limits = c(0, 1.10), breaks = seq(0, 1, 0.25)) +
  labs(tag = "B",
       title = expression("CPSS — mean "*hat(pi)),
       subtitle = "100 paired half-samples",
       x = "|log2FC| cutoff",
       y = expression("Mean "*hat(pi))) +
  theme_robust

# -----------------------------------------------------------------------------
# B' — mean bootstrap selection frequency
# -----------------------------------------------------------------------------
p_bp <- ggplot(sweep, aes(x = x_lab, y = mean_boot, fill = x_lab)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", mean_boot)),
            vjust = -0.4, size = 3) +
  scale_fill_manual(values = THR_COLS, guide = "none") +
  scale_y_continuous(limits = c(0, 1.10), breaks = seq(0, 1, 0.25)) +
  labs(tag = "B'",
       title = "Bootstrap — mean selection frequency",
       subtitle = "B = 1000 iterations; each draws 80% per cohort (cohort-stratified subsampling)",
       x = "|log2FC| cutoff",
       y = "Mean P(selected)") +
  theme_robust

# -----------------------------------------------------------------------------
# C — % canonical DEGs with CPSS pi-hat >= 0.7
# -----------------------------------------------------------------------------
p_c <- ggplot(sweep, aes(x = x_lab, y = pct_pi_ge_07, fill = x_lab)) +
  geom_hline(yintercept = 70, linetype = "dashed",
             colour = "grey25", linewidth = 0.4) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct_pi_ge_07)),
            vjust = -0.4, size = 3) +
  scale_fill_manual(values = THR_COLS, guide = "none") +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25)) +
  labs(tag = "C",
       title = expression("CPSS — % DEGs with "*hat(pi) >= 0.7),
       subtitle = "Dashed line: Shah-Samworth 2013 stability cutoff (70%)",
       x = "|log2FC| cutoff",
       y = "% canonical DEGs") +
  theme_robust

# -----------------------------------------------------------------------------
# C' — % canonical DEGs with bootstrap freq >= 0.9
# -----------------------------------------------------------------------------
p_cp <- ggplot(sweep, aes(x = x_lab, y = pct_boot_ge_09, fill = x_lab)) +
  geom_hline(yintercept = 90, linetype = "dashed",
             colour = "grey25", linewidth = 0.4) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct_boot_ge_09)),
            vjust = -0.4, size = 3) +
  scale_fill_manual(values = THR_COLS, guide = "none") +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25)) +
  labs(tag = "C'",
       title = "Bootstrap — % DEGs with freq >= 0.9",
       subtitle = "Dashed line: stringent reproducibility cutoff (90%)",
       x = "|log2FC| cutoff",
       y = "% canonical DEGs") +
  theme_robust

# -----------------------------------------------------------------------------
# E — per-iter Spearman rho (CPSS, Bootstrap, K-fold-10)
# -----------------------------------------------------------------------------
p_e <- ggplot(iter[!is.na(rho)],
              aes(x = x_lab, y = rho, fill = x_lab)) +
  geom_violin(scale = "width", trim = TRUE,
              alpha = 0.6, colour = "grey20", linewidth = 0.25) +
  geom_boxplot(width = 0.18, outlier.size = 0.4,
               fill = "white", colour = "grey20", linewidth = 0.25) +
  geom_hline(yintercept = 0.85, linetype = "dashed",
             colour = "grey25", linewidth = 0.3) +
  facet_wrap(~ source_lbl, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = THR_COLS, guide = "none") +
  labs(tag = "E",
       title = expression("Per-iter Spearman "*rho*" vs full dream"),
       subtitle = "Each violin = resampling iterations at that cutoff and source",
       x = "|log2FC| cutoff",
       y = expression(rho)) +
  theme_robust

# -----------------------------------------------------------------------------
# Compose
#   Row 1: A | D                 (count + K-fold-10)
#   Row 2: B | B'                (mean stability metric, CPSS vs bootstrap)
#   Row 3: C | C'                (% above stability cutoff, CPSS vs bootstrap)
#   Row 4: E                     (per-iter rho violins, 3 sources)
# -----------------------------------------------------------------------------
composite <- (p_a | p_d) / (p_b | p_bp) / (p_c | p_cp) / p_e +
  plot_layout(heights = c(1, 1, 1, 1.05)) +
  plot_annotation(
    caption = "Red bar = canonical |log₂FC| > 0.5 cutoff used throughout the paper."
  )

ggsave(OUT_PDF, composite, width = 12, height = 12, device = cairo_pdf)
cat("Saved:", OUT_PDF, "\n")
