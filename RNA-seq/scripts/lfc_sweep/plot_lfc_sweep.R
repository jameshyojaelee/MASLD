#!/usr/bin/env Rscript
# plot_lfc_sweep.R
# ---------------------------------------------------------------------------
# B3 — Plot LFC × significance grid sweep results, per scale.
#
# Inputs (either / both):
#   RNA-seq/results/audit_sensitivity/lfc_sweep/lfc_sweep_results_raw.csv
#   RNA-seq/results/audit_sensitivity/lfc_sweep/lfc_sweep_results_shrunk.csv
#
# Outputs (one multi-panel PDF per scale present):
#   figures/supplementary/figS_lfc_sweep_raw.pdf
#   figures/supplementary/figS_lfc_sweep_shrunk.pdf
#     Panel 1: heatmap of mean_LOO_recovery vs LFC × significance threshold
#     Panel 2: mean_LOO_recovery vs LFC (one line per sig threshold) ± sd
#     Panel 3: n_DEG vs LFC (log scale) per sig threshold
#
# The empirical CV / plateau is shown across the full sweep; no single cutoff
# is baked in as canonical.
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

IN_DIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity/lfc_sweep")

# Try project theme; fall back to theme_bw
theme_pub <- tryCatch({
  src_f <- file.path(BASE, "scripts/figures/publication_theme.R")
  if (file.exists(src_f)) {
    source(src_f, local = TRUE)
    if (exists("theme_publication")) theme_publication() else theme_bw(base_size = 10)
  } else theme_bw(base_size = 10)
}, error = function(e) theme_bw(base_size = 10))

SCALES <- list(
  raw    = list(suffix = "raw",    title = "raw |log2FC|"),
  shrunk = list(suffix = "shrunk", title = "ashr-shrunk |log2FC|")
)

for (scale in names(SCALES)) {
  sc <- SCALES[[scale]]
  IN_F  <- file.path(IN_DIR, sprintf("lfc_sweep_results_%s.csv", sc$suffix))
  OUT_F <- file.path(BASE, sprintf("figures/supplementary/figS_lfc_sweep_%s.pdf", sc$suffix))
  if (!file.exists(IN_F)) {
    cat("SKIP (missing input):", IN_F, "\n"); next
  }

  res <- fread(IN_F)
  sig_type <- if ("sig_type" %in% names(res)) res$sig_type[1] else "padj"
  res[, sig_lbl := factor(sprintf("%s < %.2f", sig_type, sig_thr),
                          levels = sprintf("%s < %.2f", sig_type, sort(unique(sig_thr))))]

  # Panel 1: heatmap of mean recovery
  pA <- ggplot(res, aes(x = factor(lfc), y = factor(sig_thr),
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
      title = sprintf("5-fold LOO recovery across grid (%s)", sc$title),
      x = sprintf("%s threshold", sc$title),
      y = sprintf("%s threshold", sig_type)
    ) +
    theme_pub +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )

  # Panel 2: line plot of mean recovery vs LFC, per sig threshold
  pB <- ggplot(res, aes(x = lfc, y = mean_LOO_recovery,
                        color = sig_lbl, group = sig_lbl)) +
    geom_ribbon(aes(ymin = mean_LOO_recovery - sd_LOO_recovery,
                    ymax = mean_LOO_recovery + sd_LOO_recovery,
                    fill = sig_lbl),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.6) +
    scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
    scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
    labs(
      title = sprintf("Mean LOO recovery vs threshold (%s)", sc$title),
      x = sprintf("%s threshold", sc$title),
      y = "Mean LOO recovery (5 folds)",
      color = "", fill = ""
    ) +
    theme_pub +
    theme(legend.position = "top")

  # Panel 3: n_DEG (log scale)
  pC <- ggplot(res, aes(x = lfc, y = n_DEG,
                        color = sig_lbl, group = sig_lbl)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.6) +
    scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
    scale_y_log10(labels = function(x) format(x, big.mark = ",")) +
    labs(
      title = sprintf("DEG count vs threshold (%s)", sc$title),
      x = sprintf("%s threshold", sc$title),
      y = "# DEGs (log scale)",
      color = ""
    ) +
    theme_pub +
    theme(legend.position = "none")

  combined <- pA / (pB | pC) + plot_layout(heights = c(1.1, 1))

  dir.create(dirname(OUT_F), recursive = TRUE, showWarnings = FALSE)
  ggsave(OUT_F, combined, width = 10, height = 8.5, device = "pdf", useDingbats = FALSE)
  cat("Saved:", OUT_F, "\n")
  message(sprintf(
    "Caption [%s]: 5-fold leave-one-cohort-out recovery of integrated DEGs across the %s x %s grid. Heatmap = mean LOO recovery; line panels = mean recovery +/- sd and DEG count vs threshold. Full sweep shown; no cutoff is canonical.",
    scale, sc$title, sig_type))
}
