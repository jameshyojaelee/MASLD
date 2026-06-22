#!/usr/bin/env Rscript
# Regenerate fig3c with fixed PIP highlighting
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrastr)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

fig_half_width <- 89/25.4
save_panel <- function(p, fname, width = fig_half_width, height = 4.0) {
  out <- file.path(FIG3_DIR, "panels", fname)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  ggsave(out, p, width = width, height = height, device = cairo_pdf)
  cat("Saved:", out, "\n")
}

# --- Load data ---
FM_BASE <- file.path(BASE, "GWAS/finemapping")
loci_tbl <- fread(file.path(FM_BASE, "results/susiex/shared_loci.csv"))
LZ_LOCUS <- "locus_ALT_chr19_41353107"
lz <- loci_tbl[locus_id == LZ_LOCUS][1]
LZ_CHR <- lz$chr; LZ_WSTART <- lz$window_start; LZ_WEND <- lz$window_end
EUR_LEAD <- lz$eur_lead_pos; EAS_LEAD <- lz$eas_lead_pos
TRAIT <- sub("locus_([^_]+)_.*", "\\1", LZ_LOCUS)  # "ALT"

lz_eur_ss <- fread(file.path(FM_BASE, sprintf("data/sumstats/UKBB_%s_reformatted_hg19.tsv", TRAIT)),
                   select = c("chromosome","position","pval"))
lz_eur_ss <- lz_eur_ss[chromosome == LZ_CHR & position >= LZ_WSTART & position <= LZ_WEND]
lz_eur_ss[, mlog10p := -log10(pmax(as.numeric(pval), 1e-300))]

lz_eas_ss <- fread(file.path(FM_BASE, sprintf("data/sumstats/BBJ_%s_reformatted_hg19.tsv", TRAIT)),
                   select = c("chromosome","position","pval"))
lz_eas_ss <- lz_eas_ss[chromosome == LZ_CHR & position >= LZ_WSTART & position <= LZ_WEND]
lz_eas_ss[, mlog10p := -log10(pmax(as.numeric(pval), 1e-300))]

lz_pips <- fread(file.path(FM_BASE, "results/mesusie/mesusie_variant_summary.csv"))
lz_pips <- lz_pips[locus_id == LZ_LOCUS]

GWS_LINE <- -log10(5e-8)
COL_EUR <- "#1565C0"; COL_EAS <- "#F4511E"

lz_common <- list(
  scale_x_continuous(limits = c(LZ_WSTART, LZ_WEND), expand = c(0.01, 0)),
  theme_masld(base_size = 6),
  theme(panel.grid.minor = element_blank(), plot.title = element_blank(),
        plot.margin = margin(t = 1, r = 6, b = 1, l = 2))
)
lz_no_x <- theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                 axis.ticks.x = element_blank())

# GWAS panels
p3c_1 <- ggplot(lz_eur_ss, aes(position, mlog10p)) +
  rasterize(geom_point(color = "grey80", size = 0.5, alpha = 0.5), dpi = 500) +
  geom_hline(yintercept = GWS_LINE, linetype = "longdash", linewidth = 0.3) +
  geom_point(data = lz_eur_ss[position == EUR_LEAD], fill = COL_EUR,
             color = "black", size = 2, shape = 21, stroke = 0.35) +
  labs(y = expression(UKBB~ALT~-log[10](italic(p)))) +
  lz_common + lz_no_x

p3c_2 <- ggplot(lz_eas_ss, aes(position, mlog10p)) +
  rasterize(geom_point(color = "grey80", size = 0.5, alpha = 0.5), dpi = 500) +
  geom_hline(yintercept = GWS_LINE, linetype = "longdash", linewidth = 0.3) +
  geom_point(data = lz_eas_ss[position == EAS_LEAD], fill = COL_EAS,
             color = "black", size = 2, shape = 21, stroke = 0.35) +
  labs(y = expression(BBJ~ALT~-log[10](italic(p)))) +
  lz_common + lz_no_x

# MESuSiE ancestry-SPECIFIC PIPs (without pip_shared).
# This locus is EAS-driven: EUR pip is near-zero, EAS pip peaks at ~0.56.
# Showing the asymmetry IS the story — do not add pip_shared.
top_eur <- lz_pips[which.max(pip_eur)]
top_eas <- lz_pips[which.max(pip_eas)]
cat(sprintf("Top PIP EUR (specific): %.3f at %d | Top PIP EAS (specific): %.3f at %d\n",
            top_eur$pip_eur, top_eur$pos, top_eas$pip_eas, top_eas$pos))
cat(sprintf("  (pip_shared max: %.3f — dominates combined pip but hidden here to show ancestry contrast)\n",
            max(lz_pips$pip_shared)))

# Clean PIP panels: grey background, lollipop + color for ancestry-specific PIP > 0.02
PIP_HI <- 0.02

p3c_3 <- ggplot(lz_pips, aes(pos, pip_eur)) +
  geom_point(color = "grey75", size = 0.8, alpha = 0.5) +
  geom_segment(data = lz_pips[pip_eur > PIP_HI],
               aes(x = pos, xend = pos, y = 0, yend = pip_eur),
               color = COL_EUR, linewidth = 0.5, alpha = 0.5) +
  geom_point(data = lz_pips[pip_eur > PIP_HI], fill = COL_EUR,
             color = "black", size = 1.8, shape = 21, stroke = 0.3) +
  geom_point(data = top_eur, fill = COL_EUR,
             color = "black", size = 3, shape = 21, stroke = 0.5) +
  scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
  labs(y = "PIP EUR") +
  lz_common + lz_no_x

p3c_4 <- ggplot(lz_pips, aes(pos, pip_eas)) +
  geom_point(color = "grey75", size = 0.8, alpha = 0.5) +
  geom_segment(data = lz_pips[pip_eas > PIP_HI],
               aes(x = pos, xend = pos, y = 0, yend = pip_eas),
               color = COL_EAS, linewidth = 0.5, alpha = 0.5) +
  geom_point(data = lz_pips[pip_eas > PIP_HI], fill = COL_EAS,
             color = "black", size = 1.8, shape = 21, stroke = 0.3) +
  geom_point(data = top_eas, fill = COL_EAS,
             color = "black", size = 3, shape = 21, stroke = 0.5) +
  geom_text(data = top_eas, aes(label = sprintf("PIP = %.2f", pip_eas)),
            hjust = -0.15, vjust = 0.4, size = 2.2, color = COL_EAS, fontface = "bold") +
  scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
  labs(x = sprintf("chr%s position (hg19, bp)", LZ_CHR), y = "PIP EAS") +
  lz_common

p3c <- (p3c_1 / p3c_2 / p3c_3 / p3c_4) +
  plot_layout(heights = c(3, 3, 2, 2)) +
  plot_annotation(
    title    = sprintf("ALT locus (chr%s:%s\u2013%s Mb)", LZ_CHR,
                       format(round(LZ_WSTART/1e6, 2), nsmall=2),
                       format(round(LZ_WEND/1e6, 2), nsmall=2)),
    subtitle = sprintf(
      "Cross-ancestry fine-mapping \u00B7 EUR lead %s bp, EAS lead %s bp (\u0394 = %d kb)\nAncestry-specific PIP: EUR max %.2f | EAS max %.2f  (EAS-driven locus)",
      format(EUR_LEAD, big.mark = ","), format(EAS_LEAD, big.mark = ","),
      abs(round((EUR_LEAD - EAS_LEAD) / 1e3)),
      top_eur$pip_eur, top_eas$pip_eas),
    theme = theme(plot.title = element_text(size = 7, face = "bold"),
                  plot.subtitle = element_text(size = 5.5, color = "grey30"))
  )

save_panel(p3c, "rora_locus_zoom.pdf", width = fig_half_width, height = 4.0)
cat("Done.\n")
