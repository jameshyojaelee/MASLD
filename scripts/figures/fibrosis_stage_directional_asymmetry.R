#!/usr/bin/env Rscript
# ============================================================================
# fibrosis_stage_directional_asymmetry.R
# Supp Fig S02 — up:down DEG ratio collapses as fibrosis advances.
#
# Per fibrosis-stage contrast (F1..F4 vs F0) count significant up- and
# down-regulated genes (padj < 0.05 & |logFC| > 0.5; canonical Tier-1) and
# plot the up:down ratio as a lollipop. Early fibrosis is induction-dominated
# (~4.8:1); advanced fibrosis approaches parity as downregulation
# (parenchymal collapse) overtakes induction (~1.6:1).
#
# Source: fibrosis_stage_dream.csv (limma_voom_qw per-stage-vs-F0 contrasts)
#
# Output: figures/supplementary/figS02_progression/fibrosis_stage_directional_asymmetry.pdf
#   sized 70 x 65 mm
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- FIGS02_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
DATA_DIR <- file.path(OUT_DIR, "data")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF  <- file.path(OUT_DIR, "fibrosis_stage_directional_asymmetry.pdf")

# ---------------------------------------------------------------------------
# Load DE results (one row per gene per stage contrast)
# ---------------------------------------------------------------------------
de <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/fibrosis_stage_dream.csv"))

# Canonical Tier-1 DEG threshold: padj < 0.05 & |logFC| > 0.5 (binary stage
# contrasts dilute fold changes; the LFC floor isolates the directional
# asymmetry signal). Direction from logFC sign.
de <- de[!is.na(padj) & !is.na(logFC)]
sig <- de[padj < 0.05 & abs(logFC) > 0.5]
sig[, dir := ifelse(logFC > 0, "up", "down")]

# ---------------------------------------------------------------------------
# Per-stage up/down counts + ratio
# ---------------------------------------------------------------------------
tab <- sig[, .(
  n_up   = sum(dir == "up"),
  n_down = sum(dir == "down")
), by = contrast]
tab[, ratio := n_up / n_down]
tab[, stage := factor(contrast, levels = c("F1_vs_F0", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0"),
                      labels = c("F1", "F2", "F3", "F4"))]
setorder(tab, stage)

fwrite(tab, file.path(DATA_DIR, "fibrosis_stage_directional_asymmetry.csv"))

cat("[hero] up:down DEG ratio per stage (padj<0.05 & |logFC|>0.5):\n")
for (i in seq_len(nrow(tab))) {
  cat(sprintf("    %s: %d up / %d down = %.2f:1\n",
              tab$stage[i], tab$n_up[i], tab$n_down[i], tab$ratio[i]))
}

# ---------------------------------------------------------------------------
# Lollipop
# ---------------------------------------------------------------------------
p <- ggplot(tab, aes(x = stage, y = ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#9E9E9E", linewidth = 0.3) +
  geom_segment(aes(xend = stage, y = 1, yend = ratio),
               color = "#6D6D6D", linewidth = 0.5) +
  geom_point(aes(color = ratio), size = 3.4) +
  geom_text(aes(label = sprintf("%.1f:1", ratio)),
            vjust = -1.0, size = 2.1, fontface = "bold") +
  scale_color_gradient(low = "#C9265E", high = "#1565C0", guide = "none") +
  scale_y_continuous(name = "Up : down DEG ratio",
                     expand = expansion(mult = c(0.02, 0.15))) +
  labs(x = "Fibrosis stage (vs F0)",
       title = "Down-regulation overtakes induction\nas fibrosis advances") +
  theme_masld(base_size = 7) +
  theme(
    plot.title  = element_text(size = 7.0, face = "bold", margin = margin(b = 6)),
    axis.text.x = element_text(size = 6.5, face = "bold")
  )

ggsave(OUT_PDF, p,
       width  = 70 / 25.4,
       height = 65 / 25.4,
       units  = "in",
       device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
