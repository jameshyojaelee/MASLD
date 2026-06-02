#!/usr/bin/env Rscript
# fig3_panel_ieqtl_concordance.R
#
# ieQTL x bulk DEG directional concordance scatter — filtered to the
# high-confidence intersection: ieQTL FDR < 0.01 + bulk padj < 0.05 + |LFC| > 0.5
# Binomial test vs 50% null annotated on plot.
#
# Output: figures/main/fig3_regulatory_architecture/panels/fig3_panel_ieqtl_concordance.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG3_DIR, "panels")
OUT_PDF   <- file.path(PANEL_DIR, "fig3_panel_ieqtl_concordance.pdf")

IEQTL_FDR_CUT <- 0.01
DEG_PADJ_CUT  <- 0.05
DEG_LFC_CUT   <- 0.5

# ── data ──────────────────────────────────────────────────────────────────────
ieqtl <- fread(file.path(BASE,
  "RNA-seq/results/causal_inference/sceqtl/ieqtl_disease_genes.csv"))

d <- ieqtl[interaction_fdr < IEQTL_FDR_CUT &
            dream_padj      < DEG_PADJ_CUT  &
            abs(dream_logFC) > DEG_LFC_CUT]

cat(sprintf("[filter] %d genes pass ieQTL FDR<%.2f + padj<%.2f + |LFC|>%.1f\n",
            nrow(d), IEQTL_FDR_CUT, DEG_PADJ_CUT, DEG_LFC_CUT))

# per-gene: keep most significant ieQTL cell type
d <- d[, .SD[which.min(interaction_pval)], by = gene]
cat(sprintf("[dedup]  %d unique genes\n", nrow(d)))

# concordance
d[, concordant := (interaction_beta > 0 & dream_logFC > 0) |
                  (interaction_beta < 0 & dream_logFC < 0)]
n_total  <- nrow(d)
n_conc   <- sum(d$concordant)
n_disc   <- n_total - n_conc
pct_conc <- 100 * n_conc / n_total

# binomial test
binom_p <- binom.test(n_conc, n_total, p = 0.5, alternative = "greater")$p.value
p_label <- if (binom_p < 0.001) "p < 0.001" else sprintf("p = %.3f", binom_p)

cat(sprintf("[concordance] %d/%d (%.1f%%) concordant, binomial %s\n",
            n_conc, n_total, pct_conc, p_label))

# cell type clean labels
ct_map <- c(hepatocyte        = "Hepatocyte",
            endothelial_cell  = "Endothelial",
            cholangiocyte     = "Cholangiocyte",
            stellate_cell     = "Stellate")
d[, ct_clean := ct_map[cell_type]]
d[is.na(ct_clean), ct_clean := cell_type]

ct_pal <- c(
  Hepatocyte     = masld_colors$fibrosis,   # deep magenta
  Endothelial    = "#00695C",               # teal
  Cholangiocyte  = "#1565C0",              # deep blue
  Stellate       = "#F57F17"               # amber
)

# label known genes + top significant
known <- c("FASN", "HKDC1", "THRB", "HSD17B13", "PNPLA3",
           "FABP1", "RORA", "COL1A1", "SPP1", "TM6SF2")
setorder(d, interaction_pval)
top_n  <- head(d, 12)
label_dt <- d[gene %in% known | gene %in% top_n$gene]
label_dt <- label_dt[, .SD[which.min(interaction_pval)], by = gene]

# axis limits with padding
x_lim <- max(abs(d$dream_logFC), na.rm = TRUE) * 1.15
y_lim <- max(abs(d$interaction_beta), na.rm = TRUE) * 1.15

# quadrant annotation data
quad_dt <- data.table(
  x     = c( x_lim, -x_lim, -x_lim,  x_lim) * 0.92,
  y     = c( y_lim,  y_lim, -y_lim, -y_lim) * 0.90,
  hjust = c(1, 0, 0, 1),
  vjust = c(1, 1, 0, 0),
  label = c(
    sprintf("Concordant\nn = %d", d[interaction_beta > 0 & dream_logFC > 0, .N]),
    sprintf("Discordant\nn = %d", d[interaction_beta > 0 & dream_logFC < 0, .N]),
    sprintf("Discordant\nn = %d", d[interaction_beta < 0 & dream_logFC > 0, .N]),
    sprintf("Concordant\nn = %d", d[interaction_beta < 0 & dream_logFC < 0, .N])
  ),
  color = c(masld_colors$up, "gray50", "gray50", masld_colors$up)
)

# ── plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(d, aes(x = dream_logFC, y = interaction_beta, color = ct_clean)) +
  # quadrant shading
  annotate("rect", xmin = 0, xmax =  x_lim, ymin = 0, ymax =  y_lim,
           fill = masld_colors$up, alpha = 0.04) +
  annotate("rect", xmin = -x_lim, xmax = 0, ymin = -y_lim, ymax = 0,
           fill = masld_colors$up, alpha = 0.04) +
  # reference lines
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_vline(xintercept =  DEG_LFC_CUT, linewidth = 0.3,
             color = "gray50", linetype = "dashed") +
  geom_vline(xintercept = -DEG_LFC_CUT, linewidth = 0.3,
             color = "gray50", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(data = label_dt,
                  aes(label = gene),
                  size = 2.8, color = "black",
                  min.segment.length = 0.2,
                  segment.color = "gray60",
                  segment.size  = 0.3,
                  box.padding   = 0.4,
                  max.overlaps  = 20) +
  # quadrant labels
  geom_text(data = quad_dt,
            aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust,
                color = NULL),
            color = quad_dt$color,
            size = 2.8, fontface = "bold", inherit.aes = FALSE) +
  # concordance annotation
  annotate("text",
           x = -x_lim * 0.98, y = y_lim * 0.98,
           hjust = 0, vjust = 1,
           label = sprintf("%d/%d (%.0f%%) concordant\nBinomial %s",
                           n_conc, n_total, pct_conc, p_label),
           size = 3, color = "gray20") +
  scale_color_manual(values = ct_pal, name = "Cell type") +
  scale_x_continuous(limits = c(-x_lim, x_lim)) +
  scale_y_continuous(limits = c(-y_lim, y_lim)) +
  labs(x = "Bulk RNA-seq log₂FC  (disease vs control)",
       y = "ieQTL interaction β  (MASLD disease effect on eQTL)",
       title    = "ieQTL × DEG directional concordance",
       subtitle = sprintf("ieQTL FDR < %.2f  |  padj < %.2f  |  |log₂FC| > %.1f",
                          IEQTL_FDR_CUT, DEG_PADJ_CUT, DEG_LFC_CUT)) +
  theme_masld(base_size = 11) +
  theme(legend.position = "right",
        legend.key.size = unit(0.35, "cm"),
        legend.text     = element_text(size = 9),
        plot.title      = element_text(size = 12, face = "bold"),
        plot.subtitle   = element_text(size = 8, color = "gray40"))

ggsave(OUT_PDF, p, width = 6, height = 5.5, device = cairo_pdf)
cat(sprintf("Saved: %s\n", OUT_PDF))
