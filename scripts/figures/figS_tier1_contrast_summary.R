#!/usr/bin/env Rscript
# figS_tier1_contrast_summary.R  (2026-05-13)
#
# Compact diverging-bar summary of Tier 1 DEG counts (padj<0.05, |logFC|>0.5)
# across the 4 canonical contrasts: Disease-vs-Control (legacy), MASH-vs-MASL,
# MASH-vs-Healthy, MASL-vs-Healthy. Strict MASH definition, PRJNA512027 excluded.
#
# Output: figures/supplementary/figS_methods_validation/lfc_sensitivity/panels/contrast_tier1_summary.pdf

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

DSIG   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures")
INTRES <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(FIGS_LFCSENS_DIR, "panels")  # consolidated under figS_methods_validation/ (2026-06-04)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

PADJ_CUT <- 0.05
LFC_CUT  <- 0.5

specs <- list(
  list(tag="DvC",      csv=file.path(INTRES, "canonical_deg_results.csv"),
       label="Disease vs Control",                        n_samples=847, n_cohorts=5),
  list(tag="MM",       csv=file.path(DSIG,   "mash_vs_masl_dream_strict.csv"),
       label="MASH vs MASL",                              n_samples=493, n_cohorts=7),
  list(tag="MH",       csv=file.path(DSIG,   "mash_vs_healthy_dream_strict.csv"),
       label="MASH vs Healthy",                           n_samples=261, n_cohorts=4),
  list(tag="MlH",      csv=file.path(DSIG,   "masl_vs_healthy_dream.csv"),
       label="MASL vs Healthy",                           n_samples=150, n_cohorts=4))

count_tier1 <- function(spec) {
  d <- fread(spec$csv)
  setnames(d, "adj.P.Val", "padj", skip_absent = TRUE)
  setnames(d, "logFC",     "lfc",  skip_absent = TRUE)
  sig <- d[!is.na(padj) & !is.na(lfc) & padj < PADJ_CUT & abs(lfc) > LFC_CUT]
  list(up = sum(sig$lfc > 0), down = sum(sig$lfc < 0))
}

counts <- rbindlist(lapply(specs, function(s) {
  k <- count_tier1(s)
  data.table(tag = s$tag, label = s$label,
             n_samples = s$n_samples, n_cohorts = s$n_cohorts,
             up = k$up, down = k$down, total = k$up + k$down)
}))
counts[, label_with_n := sprintf("%s\nn=%d (%d cohorts)", label, n_samples, n_cohorts)]
# Order from largest total to smallest, top to bottom (so DvC reference at bottom)
setorder(counts, total)
counts[, label_with_n := factor(label_with_n, levels = counts$label_with_n)]

# Long form for stacking — up positive, down negative (diverging)
plot_dt <- rbind(
  counts[, .(label_with_n, direction = "Up",   value = up,   count = up)],
  counts[, .(label_with_n, direction = "Down", value = -down, count = down)])
plot_dt[, direction := factor(direction, levels = c("Down", "Up"))]

# Asymmetric x-limits — pad just enough for the label on each side, so the bar
# panel hugs the actual data range and no whitespace on the left.
max_up   <- max(counts$up)
max_down <- max(counts$down)
x_lim_pos <- max_up   * 1.22
x_lim_neg <- max_down * 1.55     # extra pad on left for the "↓ N" label

p <- ggplot(plot_dt, aes(x = value, y = label_with_n, fill = direction)) +
  geom_col(width = 0.64) +
  geom_vline(xintercept = 0, color = "gray40", linewidth = 0.3) +
  # Up/down totals with arrows, outside bar end
  geom_text(data = counts,
            aes(x = up, y = label_with_n,
                label = sprintf("%s ↑", format(up, big.mark = ","))),
            inherit.aes = FALSE,
            hjust = -0.12, size = 3.4, fontface = "bold", color = masld_colors$up) +
  geom_text(data = counts,
            aes(x = -down, y = label_with_n,
                label = sprintf("↓ %s", format(down, big.mark = ","))),
            inherit.aes = FALSE,
            hjust = 1.12, size = 3.4, fontface = "bold", color = masld_colors$down) +
  scale_fill_manual(values = c("Up" = masld_colors$up, "Down" = masld_colors$down),
                    name = NULL,
                    breaks = c("Up", "Down"),
                    labels = c("Upregulated", "Downregulated")) +
  scale_x_continuous(limits = c(-x_lim_neg, x_lim_pos),
                     breaks = pretty_breaks(n = 5),
                     labels = function(v) format(abs(v), big.mark = ",")) +
  labs(x = "Tier 1 DEGs (padj < 0.05 & |log2FC| > 0.5)",
       y = NULL,
       title = "DEG counts across contrasts") +
  theme_masld() + theme_pub() +
  theme(legend.position    = "top",
        legend.justification = c(0, 1),
        legend.margin      = margin(0, 0, 0, 0),
        legend.text        = element_text(size = PUB_LEGEND + 2),
        plot.title         = element_text(size = PUB_TITLE + 3, face = "bold"),
        axis.text.y        = element_text(size = PUB_AXIS_TEXT + 3, lineheight = 0.95),
        axis.text.x        = element_text(size = PUB_AXIS_TEXT + 2),
        axis.title.x       = element_text(size = PUB_AXIS_TITLE + 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.margin        = margin(4, 8, 4, 4))

out_pdf <- file.path(OUT_DIR, "DEG_4_contrasts.pdf")
save_fig(p, out_pdf, width = fig_col_width * 0.95, height = 2.4)
message("Saved: ", out_pdf)

# Companion CSV
fwrite(counts, file.path(OUT_DIR, "DEG_4_contrasts_data.csv"))
message("Saved: ", file.path(OUT_DIR, "DEG_4_contrasts_data.csv"))
