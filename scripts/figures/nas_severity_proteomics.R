#!/usr/bin/env Rscript
# KEY MESSAGE: the protein-level disease signal is GRADED, not just present-or-
# absent. Of the proteins that are differentially abundant in disease (the
# proteomics DE that validates the transcriptome), a large fraction further
# change with histologic severity, every one in the same direction it took in
# disease. The metric is computed entirely WITHIN the proteomics DE (liver and
# plasma) -- no external gene panel -- so the denominator is the validated
# disease proteome itself.
#
# NUMBERS (computed live from the canonical DE table, self-consistent):
#   Liver  (PXD051911):  of 218 disease-vs-control sig proteins, 141 (65%) are
#          also sig across NAS grade (NAS>=4 vs <4); 141/141 same direction.
#   Plasma (PXD052937):  of 94 disease-vs-control sig proteins, 35 (37%) also
#          separate MASH from MASL; 35/35 same direction.
#   Severity-driven exemplars (NAS log2FC): GPNMB, THY1, SCD, FASN, AKR1B10, ANXA2.
# Source: Analysis/Proteomics/results/protein_differential_results_v3.csv
#   (limma, BH-adjusted padj; disease-vs-control + NAS-high-vs-low + MASH-vs-MASL).
#
# CAVEAT (emitted to stdout for the legend): liver disease DE is n=58 (solid);
# plasma disease DE rests on 7 controls, so the plasma 94/35 inherit that
# under-powering. Severity axes differ by compartment (liver = NAS grade; plasma
# = MASH vs MASL, its only severity contrast).
#
# Output: figures/main/fig4_validation/nas_severity_proteomics.pdf
# Env:    rnaseq

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

diff <- fread(file.path(BASE, "Analysis/Proteomics/results/protein_differential_results_v3.csv"))
diff <- diff[!is.na(gene) & gene != ""]          # drop blank-symbol rows (off-by-one guard)
sig_genes <- function(ds) diff[dataset == ds & !is.na(padj) & padj < 0.05, unique(gene)]
lfc_of    <- function(ds) { x <- diff[dataset == ds & !is.na(logFC)]; x <- x[!duplicated(gene)]; setNames(x$logFC, x$gene) }

# ── Per-compartment: disease-altered proteins that also track severity ────────
liv_dis <- sig_genes("PXD051911")
liv_nas <- sig_genes("PXD051911_nas_high_vs_low")
liv_ov  <- intersect(liv_dis, liv_nas)
pl_dis  <- sig_genes("PXD052937")
pl_mvm  <- sig_genes("PXD052937_mash_vs_masl")
pl_ov   <- intersect(pl_dis, pl_mvm)

# direction concordance (disease vs severity) within each overlap
lf_ld <- lfc_of("PXD051911");  lf_ln <- lfc_of("PXD051911_nas_high_vs_low")
lf_pd <- lfc_of("PXD052937");  lf_pm <- lfc_of("PXD052937_mash_vs_masl")
liv_conc <- sum(sign(lf_ld[liv_ov]) == sign(lf_ln[liv_ov]), na.rm = TRUE)
pl_conc  <- sum(sign(lf_pd[pl_ov])  == sign(lf_pm[pl_ov]),  na.rm = TRUE)

message(sprintf("[nas_severity] LIVER: %d/%d (%.0f%%) disease-sig proteins also NAS-sig; %d/%d same direction",
                length(liv_ov), length(liv_dis), 100*length(liv_ov)/length(liv_dis), liv_conc, length(liv_ov)))
message(sprintf("[nas_severity] PLASMA: %d/%d (%.0f%%) disease-sig proteins also MASH-vs-MASL-sig; %d/%d same direction",
                length(pl_ov), length(pl_dis), 100*length(pl_ov)/length(pl_dis), pl_conc, length(pl_ov)))
message("[nas_severity] NAS = NAFLD Activity Score (0-8: steatosis + lobular inflammation + hepatocyte ballooning); contrast = High (NAS>=4, n=20) vs Low (NAS<4, n=38) liver biopsies, BMI+age adjusted (PXD051911).")
message("[nas_severity] CAVEAT: liver disease DE n=58 (solid); plasma disease DE 7 controls (underpowered). Severity axis: liver=NAS>=4-vs-<4, plasma=MASH-vs-MASL (its only severity contrast).")

# ════════════════════════════════════════════════════════════════════════════
# Panel A — fraction of disease-altered proteins that also track severity
# ════════════════════════════════════════════════════════════════════════════
rate_df <- data.table(
  lab  = c("Plasma  (MASH vs MASL)", "Liver  (NAS >=4 vs <4)"),
  num  = c(length(pl_ov), length(liv_ov)),
  den  = c(length(pl_dis), length(liv_dis))
)
rate_df[, rate := num / den]
rate_df[, lab := factor(lab, levels = lab)]      # plasma bottom, liver top
rate_cols <- c(masld_colors$down, masld_colors$up)

pA <- ggplot(rate_df, aes(x = rate, y = lab, color = lab)) +
  geom_segment(aes(x = 0, xend = rate, yend = lab), linewidth = 1.1) +
  geom_point(size = 2.8) +
  geom_text(aes(label = sprintf("%.0f%%  (%d/%d)", rate*100, num, den)),
            hjust = 0, nudge_x = 0.02, size = PUB_GEOM_TEXT + 0.4,
            fontface = "bold", color = "black") +
  scale_color_manual(values = rate_cols, guide = "none") +
  scale_x_continuous(limits = c(0, 0.85), breaks = c(0, 0.25, 0.5, 0.75),
                     labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Disease-altered proteins that\nalso track severity (padj < 0.05)", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_masld() + theme_pub() +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(face = "bold", size = PUB_AXIS_TEXT),
        plot.margin = margin(5.5, 34, 5.5, 5.5))

# ════════════════════════════════════════════════════════════════════════════
# Panel B — top severity-driven proteins (liver disease-sig AND NAS-sig)
# ════════════════════════════════════════════════════════════════════════════
# Highlight a few recognizable severity-driven disease proteins that are in this
# set (disease-DE AND NAS-DE): GPNMB (MASH-severity marker), PLIN2 (lipid
# droplet), ANXA2, and the matrix proteins COL14A1 / LAMA5. These match the genes
# named in the fig4 text.
anchors <- c("GPNMB", "PLIN2", "ANXA2", "COL14A1", "LAMA5")
nas_tab <- diff[dataset == "PXD051911_nas_high_vs_low" & gene %in% liv_ov &
                !is.na(logFC) & !is.na(padj)]
ord  <- nas_tab[order(-logFC)]
show <- unique(rbind(ord[gene %in% anchors], ord))[seq_len(min(14L, nrow(ord)))]
show <- show[order(logFC)]
show[, gene := factor(gene, levels = gene)]

lab_max <- max(show$logFC)
pB <- ggplot(show, aes(x = logFC, y = gene)) +
  geom_segment(aes(x = 0, xend = logFC, yend = gene), color = masld_colors$up, linewidth = 0.9) +
  geom_point(aes(size = -log10(padj)), color = masld_colors$up) +
  scale_size_continuous(range = c(1.2, 3.2), name = "-log10 padj", breaks = c(2, 4, 6)) +
  scale_x_continuous(limits = c(0, lab_max * 1.12), breaks = c(0, 0.5, 1.0, 1.5),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "log2FC, NAS >=4 vs <4", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_masld() + theme_pub() +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = PUB_AXIS_TEXT, color = "black"),
        legend.position = c(0.82, 0.28),
        legend.background = element_blank(),
        plot.margin = margin(5.5, 8, 5.5, 5.5))

p <- pA + pB + plot_layout(widths = c(1, 1))

out <- file.path(FIG4_DIR, "_supp", "nas_severity_proteomics.pdf")
cairo_pdf(out, width = fig_full_width * 0.82, height = fig_half_width * 0.62, family = "Helvetica")
print(p)
dev.off()
message("Saved: ", out)
