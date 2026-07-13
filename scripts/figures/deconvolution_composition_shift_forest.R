#!/usr/bin/env Rscript
# ============================================================================
# deconvolution_composition_shift_forest.R
# Supp Fig S03 panel — bulk deconvolution composition shifts (forest)
#
# Bulk deconvolution shows significant hepatocyte DEPLETION and coordinated
# significant EXPANSION of fibroblasts / macrophages / T / endothelial cells
# along BOTH disease (Control_vs_Disease) and fibrosis (F_low_vs_F_high) axes.
#
# logit_diff +/- 95% CI per cell type, colored by BH-padj significance,
# faceted by contrast. CI derived from logit_diff and t: SE = logit_diff / t,
# 95% CI = logit_diff +/- 1.96 * SE.
#
# Output: figures/supplementary/figS03_deconvolution/
#           deconvolution_composition_shift_forest.pdf  (130 x 95 mm)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- FIGS03_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF <- file.path(OUT_DIR, "deconvolution_composition_shift_forest.pdf")

# ---------------------------------------------------------------------------
# Load + derive 95% CI from logit_diff and t-stat (no explicit CI columns)
# ---------------------------------------------------------------------------
d <- fread(file.path(BASE, "RNA-seq/results/celltype_attribution/composition_shifts.csv"))
d[, se := logit_diff / t]
d[, ci_lo := logit_diff - 1.96 * se]
d[, ci_hi := logit_diff + 1.96 * se]
d[, sig := padj < 0.05]

contr_levels <- c("Control_vs_Disease", "F_low_vs_F_high")
contr_labels <- c(Control_vs_Disease = "Disease vs Control",
                  F_low_vs_F_high     = "Fibrosis low vs high")
d <- d[contrast %in% contr_levels]
d[, contrast := factor(contrast, levels = contr_levels,
                       labels = contr_labels[contr_levels])]

# Order cell types by disease-axis logit_diff (consistent across facets)
ord <- d[contrast == contr_labels["Control_vs_Disease"]][order(logit_diff)]
d[, celltype := factor(celltype, levels = ord$celltype)]

# ---------------------------------------------------------------------------
# Hero numbers
# ---------------------------------------------------------------------------
hep <- d[celltype == "Hepatocytes"]
for (i in seq_len(nrow(hep))) {
  cat(sprintf("[hero] Hepatocytes %s: logit_diff = %.4f, padj = %.2e (%s)\n",
              as.character(hep$contrast[i]), hep$logit_diff[i], hep$padj[i],
              ifelse(hep$sig[i], "SIG", "n.s.")))
}
sig_both <- d[sig == TRUE, .N, by = celltype][N == 2]$celltype
sig_both <- setdiff(sig_both, "Hepatocytes")
cat(sprintf("[hero] non-hep cell types SIG on BOTH axes: %s\n",
            paste(sort(as.character(sig_both)), collapse = ", ")))

fwrite(d, file.path(OUT_DIR, "deconvolution_composition_shift_forest_data.csv"))

# ---------------------------------------------------------------------------
# Forest plot
# ---------------------------------------------------------------------------
sig_cols <- c(`TRUE` = "#C9265E", `FALSE` = "#9E9E9E")

p <- ggplot(d, aes(x = logit_diff, y = celltype, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#9E9E9E", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, linewidth = 0.4) +
  geom_point(size = 1.6) +
  facet_wrap(~ contrast, nrow = 1) +
  scale_color_manual(values = sig_cols, name = "BH padj < 0.05",
                     labels = c(`TRUE` = "Significant", `FALSE` = "n.s.")) +
  labs(x = "Composition shift (logit difference)", y = NULL) +
  theme_masld(base_size = 6) +
  theme(
    axis.text.y     = element_text(size = 6),
    strip.text      = element_text(size = 6, face = "plain"),
    legend.position = "bottom",
    legend.text     = element_text(size = 6),
    legend.title    = element_text(size = 6),
    legend.key.size = unit(0.2, "cm")
  )

message("[caption] Bulk deconvolution composition shifts")

ggsave(OUT_PDF, p,
       width  = 130 / 25.4,
       height = 95 / 25.4,
       units  = "in",
       device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
