#!/usr/bin/env Rscript
# ============================================================================
# scrna_clr_abundance_forest.R
# Supp Fig S03 panel — donor-level CLR compositional mixed model (forest)
#
# A donor-level compositional (CLR) linear mixed model finds only
# CHOLANGIOCYTES (padj ~8.5e-4) and ENDOTHELIAL cells (padj ~0.018)
# significantly EXPAND in MASLD; hepatocytes are NOT significant.
#
# beta +/- 95% CI (ci_lo / ci_hi) per cell type, colored by BH-padj
# significance. STANDALONE — must NOT be cross-plotted against bulk
# deconvolution (the scRNA arm is FACS-confounded).
#
# Output: figures/supplementary/figS03_deconvolution/
#           scrna_clr_abundance_forest.pdf  (95 x 90 mm)
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
OUT_PDF <- file.path(OUT_DIR, "scrna_clr_abundance_forest.pdf")

# ---------------------------------------------------------------------------
# Load (explicit CI columns supplied by the model)
# ---------------------------------------------------------------------------
d <- fread(file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/compositional/lmm_binary_results.csv"))
d[, sig := significant == TRUE | padj < 0.05]
setorder(d, beta)
d[, cell_type := factor(cell_type, levels = d$cell_type)]

# ---------------------------------------------------------------------------
# Hero numbers
# ---------------------------------------------------------------------------
sig_cells <- d[sig == TRUE][order(padj)]
for (i in seq_len(nrow(sig_cells))) {
  cat(sprintf("[hero] SIG: %s beta = %.3f, padj = %.2e (%s)\n",
              sig_cells$cell_type[i], sig_cells$beta[i], sig_cells$padj[i],
              sig_cells$direction[i]))
}
hep <- d[cell_type == "Hepatocytes"]
cat(sprintf("[hero] Hepatocytes beta = %.3f, padj = %.3f (NS)\n",
            hep$beta, hep$padj))

fwrite(d, file.path(OUT_DIR, "scrna_clr_abundance_forest_data.csv"))

# ---------------------------------------------------------------------------
# Forest plot (standalone; emphasize significant cells)
# ---------------------------------------------------------------------------
sig_cols <- c(`TRUE` = "#C9265E", `FALSE` = "#9E9E9E")

p <- ggplot(d, aes(x = beta, y = cell_type, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#9E9E9E", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, linewidth = 0.4) +
  geom_point(aes(size = sig)) +
  scale_color_manual(values = sig_cols, name = "BH padj < 0.05",
                     labels = c(`TRUE` = "Significant", `FALSE` = "n.s.")) +
  scale_size_manual(values = c(`TRUE` = 2.1, `FALSE` = 1.4), guide = "none") +
  labs(x = "CLR abundance shift (MASLD vs healthy, beta)", y = NULL,
       title = "Donor-level compositional (CLR) mixed model") +
  theme_masld(base_size = 7) +
  theme(
    plot.title      = element_text(size = 7.3, face = "bold", margin = margin(b = 6)),
    axis.text.y     = element_text(size = 6),
    legend.position = "bottom",
    legend.text     = element_text(size = 5.5),
    legend.title    = element_text(size = 6),
    legend.key.size = unit(0.2, "cm")
  )

ggsave(OUT_PDF, p,
       width  = 95 / 25.4,
       height = 90 / 25.4,
       units  = "in",
       device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
