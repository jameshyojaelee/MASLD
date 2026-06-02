#!/usr/bin/env Rscript
# plot_lfc_sweep.R
# ---------------------------------------------------------------------------
# B3 — Plot LFC × padj grid sweep results.
#
# Inputs:
#   RNA-seq/results/audit_sensitivity/lfc_sweep/lfc_sweep_results.csv
#
# Outputs:
#   figures/supplementary/figS_lfc_sweep.pdf  (multi-panel)
#     Panel A: heatmap of mean_LOO_recovery vs LFC × padj
#     Panel B: line plot of mean_LOO_recovery vs LFC (one line per padj)
#                 with error bars (sd across 5 folds)
#     Panel C: line plot of n_DEG vs LFC (log scale) per padj
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)

IN_F   <- file.path(BASE,
  "RNA-seq/results/audit_sensitivity/lfc_sweep/lfc_sweep_results.csv")
OUT_F  <- file.path(BASE, "figures/supplementary/figS_lfc_sweep.pdf")

stopifnot(file.exists(IN_F))
res <- fread(IN_F)
res[, padj_lbl := factor(sprintf("padj < %.2f", padj),
                         levels = sprintf("padj < %.2f",
                                          sort(unique(padj))))]

# Try project theme; fall back to theme_bw
theme_pub <- tryCatch({
  src_f <- file.path(BASE, "scripts/figures/publication_theme.R")
  if (file.exists(src_f)) {
    source(src_f, local = TRUE)
    if (exists("theme_publication"))
      theme_publication() else theme_bw(base_size = 10)
  } else theme_bw(base_size = 10)
}, error = function(e) theme_bw(base_size = 10))

# Panel A: heatmap of mean recovery
pA <- ggplot(res, aes(x = factor(lfc), y = factor(padj),
                      fill = mean_LOO_recovery)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", mean_LOO_recovery)),
            size = 2.4, color = "black") +
  scale_fill_gradientn(
    colours = c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"),
    name    = "mean LOO\nrecovery",
    limits  = c(min(res$mean_LOO_recovery, na.rm = TRUE), 1)
  ) +
  scale_y_discrete(limits = rev) +
  labs(
    title = "A. 5-fold LOO recovery across LFC × padj grid",
    x = expression("|log"[2]*" fold-change| threshold"),
    y = "padj threshold"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

# Panel B: line plot of mean recovery vs LFC, per padj
pB <- ggplot(res, aes(x = lfc, y = mean_LOO_recovery,
                      color = padj_lbl, group = padj_lbl)) +
  geom_ribbon(aes(ymin = mean_LOO_recovery - sd_LOO_recovery,
                  ymax = mean_LOO_recovery + sd_LOO_recovery,
                  fill = padj_lbl),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  annotate("text", x = 0.5, y = min(res$mean_LOO_recovery, na.rm = TRUE),
           label = "current\n|LFC| = 0.5", hjust = -0.05, vjust = 0,
           size = 2.8, color = "grey30") +
  labs(
    title = "B. Mean LOO recovery vs LFC threshold (plateau)",
    x = expression("|log"[2]*" fold-change| threshold"),
    y = "Mean LOO recovery (5 folds)",
    color = "", fill = ""
  ) +
  theme_pub +
  theme(legend.position = "top")

# Panel C: n_DEG (log scale)
pC <- ggplot(res, aes(x = lfc, y = n_DEG,
                      color = padj_lbl, group = padj_lbl)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",")) +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  labs(
    title = "C. DEG count vs LFC threshold",
    x = expression("|log"[2]*" fold-change| threshold"),
    y = "# DEGs (log scale)",
    color = ""
  ) +
  theme_pub +
  theme(legend.position = "none")

combined <- pA / (pB | pC) + plot_layout(heights = c(1.1, 1))

dir.create(dirname(OUT_F), recursive = TRUE, showWarnings = FALSE)
ggsave(OUT_F, combined, width = 10, height = 8.5, device = "pdf")
cat("Saved:", OUT_F, "\n")
