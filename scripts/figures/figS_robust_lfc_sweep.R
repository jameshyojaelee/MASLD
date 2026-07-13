#!/usr/bin/env Rscript
##############################################################################
# figS_robust_lfc_sweep.R
#
# Pillar-A resampling-robustness vs |log2FC| cutoff, rendered at TWO scales:
#   raw    : DEGs gated padj<0.05 & |log2FC|>x          -> *_raw.pdf
#   shrunk : DEGs gated lfsr<0.05 & |shrunk log2FC|>x   -> *_shrunk.pdf  (PRIMARY)
#
# Each scale reads its own sweep summaries (produced by
# RNA-seq/.../scripts/aggregate_pillar_A_shrunk_sweep.R):
#   pillar_A_lfc_threshold_sweep_<scale>.csv
#   pillar_A_lfc_sweep_iter_long_<scale>.csv
#
# Panels (descriptive, prefix-free filenames; panel identity in title + tag):
#   DEG count | K-fold(K=10) recurrence | CPSS mean pi-hat |
#   Bootstrap mean selection freq | CPSS % DEGs pi>=0.7 |
#   Bootstrap % DEGs freq>=0.9 | per-iter Spearman rho
# House style: black text only, no subtitles, short titles, caption via message().
##############################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
AUDIT   <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
OUT_DIR <- FIGS_ROBUST_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

theme_robust <- theme_minimal(base_size = 6) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
    axis.line   = element_line(colour = "black", linewidth = 0.3),
    axis.ticks  = element_line(colour = "black", linewidth = 0.3),
    axis.text   = element_text(colour = "black"),
    axis.title  = element_text(colour = "black"),
    legend.position = "top",
    legend.title  = element_blank(),
    legend.text   = element_text(size = 6, colour = "black"),
    plot.title    = element_text(face = "plain", size = 6, colour = "black",
                                 margin = margin(b = 4)),
    plot.subtitle = element_blank(),
    plot.tag      = element_text(face = "plain", size = 6, colour = "black"),
    plot.margin   = margin(8, 10, 8, 10),
    strip.background = element_blank(),
    strip.text       = element_text(face = "plain", size = 6, colour = "black")
  )

THR_LEVS    <- c("0", "0.1", "0.3", "0.5", "0.75", "1", "2")
BAR_FILL    <- "#4F8FBA"
CUTOFF_FILL <- "#C96450"   # canonical |log2FC|>0.5 cutoff highlight
THR_COLS    <- setNames(rep(BAR_FILL, length(THR_LEVS)), THR_LEVS)
THR_COLS["0.5"] <- CUTOFF_FILL

render_scale <- function(scale) {
  eff_lab  <- if (scale == "raw") "raw |log2FC|" else "ashr-shrunk |log2FC|"
  x_lab    <- if (scale == "raw") "|log2FC| cutoff" else "|shrunk log2FC| cutoff"

  sweep_f <- file.path(AUDIT, sprintf("pillar_A_lfc_threshold_sweep_%s.csv", scale))
  iter_f  <- file.path(AUDIT, sprintf("pillar_A_lfc_sweep_iter_long_%s.csv", scale))
  if (!file.exists(sweep_f) || !file.exists(iter_f)) {
    message("[skip] missing inputs for scale=", scale); return(invisible())
  }
  sweep <- fread(sweep_f); setorder(sweep, threshold)
  sweep[, x_lab := factor(sprintf("%g", threshold), levels = THR_LEVS)]

  iter <- fread(iter_f)
  iter <- iter[source != "kfold_K50"]          # drop K=50 per PI request
  iter[, x_lab := factor(sprintf("%g", threshold), levels = THR_LEVS)]
  iter[, source_lbl := factor(source,
         levels = c("cpss", "bootstrap", "kfold_K10"),
         labels = c("CPSS", "Bootstrap", "K-fold (K=10)"))]

  p_count <- ggplot(sweep, aes(x = x_lab, y = n_DEG, fill = x_lab)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = comma(n_DEG)), vjust = -0.4, size = GEOM_TEXT_6PT, colour = "black") +
    scale_fill_manual(values = THR_COLS, guide = "none") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.18))) +
    labs(tag = "A",
         x = x_lab, y = "DEGs") + theme_robust

  p_kfold <- ggplot(sweep, aes(x = x_lab, y = 100 * mean_kfold10 / 10, fill = x_lab)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.0f%%", 100 * mean_kfold10 / 10)),
              vjust = -0.4, size = GEOM_TEXT_6PT, colour = "black") +
    scale_fill_manual(values = THR_COLS, guide = "none") +
    scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25)) +
    labs(tag = "D",
         x = x_lab, y = "% folds re-selected") + theme_robust

  p_pi <- ggplot(sweep, aes(x = x_lab, y = mean_pi_hat, fill = x_lab)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.3f", mean_pi_hat)),
              vjust = -0.4, size = GEOM_TEXT_6PT, colour = "black") +
    scale_fill_manual(values = THR_COLS, guide = "none") +
    scale_y_continuous(limits = c(0, 1.10), breaks = seq(0, 1, 0.25)) +
    labs(tag = "B",
         x = x_lab, y = expression("Mean "*hat(pi))) + theme_robust

  p_boot <- ggplot(sweep, aes(x = x_lab, y = mean_boot, fill = x_lab)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.3f", mean_boot)),
              vjust = -0.4, size = GEOM_TEXT_6PT, colour = "black") +
    scale_fill_manual(values = THR_COLS, guide = "none") +
    scale_y_continuous(limits = c(0, 1.10), breaks = seq(0, 1, 0.25)) +
    labs(tag = "B'",
         x = x_lab, y = "Mean P(selected)") + theme_robust

  p_pi_pct <- ggplot(sweep, aes(x = x_lab, y = pct_pi_ge_07, fill = x_lab)) +
    geom_hline(yintercept = 70, linetype = "dashed", colour = "grey25", linewidth = 0.4) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", pct_pi_ge_07)),
              vjust = -0.4, size = GEOM_TEXT_6PT, colour = "black") +
    scale_fill_manual(values = THR_COLS, guide = "none") +
    scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25)) +
    labs(tag = "C",
         x = x_lab, y = "% DEGs") + theme_robust

  p_boot_pct <- ggplot(sweep, aes(x = x_lab, y = pct_boot_ge_09, fill = x_lab)) +
    geom_hline(yintercept = 90, linetype = "dashed", colour = "grey25", linewidth = 0.4) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", pct_boot_ge_09)),
              vjust = -0.4, size = GEOM_TEXT_6PT, colour = "black") +
    scale_fill_manual(values = THR_COLS, guide = "none") +
    scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25)) +
    labs(tag = "C'",
         x = x_lab, y = "% DEGs") + theme_robust

  p_rho <- ggplot(iter[!is.na(rho)], aes(x = x_lab, y = rho, fill = x_lab)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.6,
                colour = "grey20", linewidth = 0.25) +
    geom_boxplot(width = 0.18, outlier.size = 0.4, fill = "white",
                 colour = "grey20", linewidth = 0.25) +
    geom_hline(yintercept = 0.85, linetype = "dashed", colour = "grey25", linewidth = 0.3) +
    facet_wrap(~ source_lbl, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = THR_COLS, guide = "none") +
    labs(tag = "E",
         x = x_lab, y = expression(rho)) + theme_robust

  # ---- individual panels (prefix-free descriptive names) ----
  panels <- list(
    deg_count                 = p_count,
    kfold10_recurrence        = p_kfold,
    cpss_mean_pi              = p_pi,
    bootstrap_mean_freq       = p_boot,
    cpss_pct_pi_ge07          = p_pi_pct,
    bootstrap_pct_freq_ge09   = p_boot_pct,
    periter_rho               = p_rho)
  for (nm in names(panels)) {
    w <- if (nm == "periter_rho") fig_full_width else 4.2
    ggsave(file.path(OUT_DIR, sprintf("robust_lfc_sweep_%s_%s.pdf", nm, scale)),
           panels[[nm]], width = w, height = 3.4, device = cairo_pdf, useDingbats = FALSE)
  }

  # ---- composite ----
  composite <- (p_count | p_kfold) / (p_pi | p_boot) / (p_pi_pct | p_boot_pct) / p_rho +
    plot_layout(heights = c(1, 1, 1, 1.05))
  out_pdf <- file.path(OUT_DIR, sprintf("figS_robust_lfc_sweep_%s.pdf", scale))
  ggsave(out_pdf, composite, width = fig_full_width, height = fig_full_width, device = cairo_pdf, useDingbats = FALSE)

  message(sprintf("Caption [%s]: Pillar-A resampling robustness of disease-vs-control DEGs as the %s cutoff tightens (x). Bars = the |log2FC|>0.5 canonical cutoff highlighted (%s). Dashed lines: CPSS Shah-Samworth 0.7 pi-hat stability cutoff; bootstrap 0.9 reproducibility cutoff; per-iter rho 0.85. K=50 fold panel dropped per PI request. Panels: A = DEG count, B = CPSS mean pi-hat, B' = bootstrap mean selection frequency, C = CPSS %% DEGs pi-hat>=0.7, C' = bootstrap %% DEGs freq>=0.9, D = K-fold (K=10) recurrence, E = per-iter Spearman rho vs full fit.",
                  scale, eff_lab, CUTOFF_FILL))
  cat("Saved:", out_pdf, "\n")
}

for (s in c("raw", "shrunk")) render_scale(s)
