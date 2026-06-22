#!/usr/bin/env Rscript
##############################################################################
# fig3a_hybrid.R
#
# Hybrid replacement for fig3a (the twin Manhattan PIP/PP.H4 panel).
#
#   Top strip   : compact Manhattan of locus-level max PIP across the 20 target
#                 GWAS (14 EUR + 3 BBJ EAS), with alternating chromosome bands.
#                 No GWAS color encoding here — this strip exists to convey
#                 "we surveyed the genome, here are the high-PIP loci".
#
#   Bottom plot : per-gene max(PIP) vs max(PP.H4) scatter.
#                 - one point per gene
#                 - color  = best-evidence GWAS category
#                 - size   = number of GWAS where gene reaches PP.H4 > 0.5
#                 - canonical MASLD genes + top corner genes labelled
#                 - subtle "high-confidence corner" shading
#                 - marginal rugs along top and right axes
#
# Output: figures/main/fig2_genetics/panels/coloc_deg_hybrid.pdf
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
  library(ggrastr)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG3_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# GWAS portfolio (20 validated)
# ---------------------------------------------------------------------------
EUR_17 <- c(
  "2019_31311600_NAFLD_EUR", "2020_32298765_NAFLD_EUR",
  "2021_34128465_PDFF_EUR",  "2021_34841290_NAFLD_EUR",
  "2021_34957434_PDFF_EUR",  "2022_36402844_PDFF_EUR",
  "2023_36280732_NAFLD_deCode_EUR", "2023_36280732_NAFLD_Intermountain_EUR",
  "2023_36280732_NAFLD_UKBB_EUR", "FinnGen_NAFLD", "FinnGen_NASH",
  "UKBB_ALT", "UKBB_AST", "UKBB_GGT"
)
FINNGEN_3       <- c("FinnGen_NAFLD", "FinnGen_NASH")
BBJ_ENZYME_3    <- c("BBJ_ALT", "BBJ_AST", "BBJ_GGT")
BBJ_5           <- BBJ_ENZYME_3
TARGET_22       <- c(EUR_17, BBJ_5)
# 23-GWAS COLOC refactor (2026-06-06) dropped cirrhosis/HCC GWAS; the displayed
# fine-mapping panel is the 17 MASLD/enzyme/PDFF studies: 14 EUR (incl. 2 FinnGen) + 3 BBJ EAS.
stopifnot(length(TARGET_22) == 17L)

gwas_category <- function(study) {
  fcase(
    study %in% BBJ_ENZYME_3,  "BBJ liver enzyme (EAS)",
    study %in% FINNGEN_3,     "FinnGen R12",
    grepl("NAFLD|NASH|Cirrhosis|HCC", study), "EUR disease",
    grepl("PDFF", study),                     "EUR PDFF",
    grepl("UKBB_(ALT|AST|GGT)$", study),      "EUR liver enzyme",
    default = "Other"
  )
}

# Convergence palette: where does this gene's causal signal show up?
# "Cross-ancestry" (EUR + EAS both converge on the same gene) is the
# headline category and gets the most saturated color.
anc_colors <- c(
  "Cross-ancestry" = "#C9265E",   # Liang deep magenta — EUR + EAS converge
  "EUR-only"       = "#5E2C5F",   # deep wine — EUR convergence only
  "EAS-only"       = "#EF6C00"    # warm orange — EAS convergence only
)
anc_levels <- names(anc_colors)

# ===========================================================================
# Data loading
# ===========================================================================
cat("[fig3a-hybrid] Loading susie_coloc_all_gwas.csv ...\n")
sc <- fread(file.path(BASE,
            "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"),
            select = c("gene", "gwas_name", "PP.H4.abf", "PP.H4.susie",
                       "top_snp_PP", "top_snp", "chr"))
sc <- sc[gwas_name %in% TARGET_22]
sc[, pp4      := pmax(PP.H4.susie, PP.H4.abf, na.rm = TRUE)]
sc[, ancestry := fifelse(gwas_name %in% BBJ_5, "EAS", "EUR")]

cat("[fig3a-hybrid] Loading combined_finemapping.csv (canonical) + ",
    "combined_finemapping_coloc_targeted.csv (COLOC-targeted re-run) ...\n")
fmap_canonical <- fread(file.path(BASE,
              "GWAS/finemapping/results/combined_finemapping.csv"),
              select = c("study", "chromosome", "position",
                         "susie_pip", "recommended_pip", "susie_cs"))
fmap_targeted <- fread(file.path(BASE,
              "GWAS/finemapping/results/combined_finemapping_coloc_targeted.csv"),
              select = c("study", "chromosome", "position",
                         "susie_pip", "recommended_pip", "susie_cs"))
fmap <- rbindlist(list(fmap_canonical, fmap_targeted), use.names = TRUE)
fmap <- fmap[study %in% TARGET_22 & chromosome %in% 1:22]
# de-duplicate on (study, chromosome, position) — when a variant is in both
# tables, prefer the row with the higher recommended_pip (canonical and
# targeted use the same SuSiE+CARMA pipeline, so duplicates should be rare;
# higher PIP wins in case of any drift).
fmap <- fmap[order(study, chromosome, position, -recommended_pip)]
fmap_unique <- unique(fmap, by = c("study", "chromosome", "position"))
cat(sprintf("  fmap_canonical=%d  fmap_targeted=%d  combined_unique=%d\n",
            nrow(fmap_canonical), nrow(fmap_targeted), nrow(fmap_unique)))

# ===========================================================================
# Per-gene scatter : x = SuSiE fine-mapping PIP of the COLOC top variant
#
# For each (gene, gwas) COLOC row, parse top_snp into (chr, pos) and
# inner-join against combined_finemapping.csv on (study, chr, pos). Only
# variants that *both* COLOC picked AND fine-mapping output cover are kept,
# so x is the actual SuSiE PIP of the colocalizing variant.
#
# Coverage: canonical fine-mapping (Step 01) clumped at p<5e-8, so sub-
# significant COLOC peaks were originally uncovered. The 2026-05-05 targeted
# re-run added 456 loci across 19 GWAS at COLOC-positive sites (PP.H4>=0.5)
# with a p<=1e-5 sanity gate, recovering THRB/MLIP and ~669 other genes.
# ===========================================================================
sc_pos <- sc[!is.na(top_snp_PP) & !is.na(pp4) & !is.na(top_snp)]
sc_pos[, snp_chr := as.integer(tstrsplit(top_snp, ":", fixed = TRUE)[[1]])]
sc_pos[, snp_pos := as.integer(tstrsplit(top_snp, ":", fixed = TRUE)[[2]])]

# Inner-join COLOC top_snps to fine-mapping: keep only matched variants
sc_match <- merge(sc_pos, fmap_unique,
                  by.x = c("gwas_name", "snp_chr", "snp_pos"),
                  by.y = c("study",     "chromosome", "position"))
cat(sprintf("  [scatter] COLOC rows with fine-mapping coverage: %d / %d\n",
            nrow(sc_match), nrow(sc_pos)))

# Best-pp4 matched row per gene; x = recommended_pip from fine-mapping
gene_best <- sc_match[sc_match[, .I[which.max(pp4)], by = gene]$V1,
                      .(best_gene = gene, best_gwas = gwas_name,
                        top_var_pp = recommended_pip, max_pp4 = pp4,
                        in_cs = !is.na(susie_cs) & susie_cs > 0)]

# Per-gene cross-ancestry status (best PP.H4 ≥ 0.5 in EUR / EAS / both)
gene_anc <- sc_pos[pp4 >= 0.5,
                   .(has_eur = any(gwas_name %in% EUR_17),
                     has_eas = any(gwas_name %in% BBJ_5)),
                   by = gene]
gene_anc[, ancestry_class := fcase(
  has_eur & has_eas, "Cross-ancestry",
  has_eur,           "EUR",
  has_eas,           "EAS",
  default            = NA_character_
)]

# Number of GWAS supporting the gene at PP.H4 ≥ 0.5
gene_n_gwas <- sc_pos[pp4 >= 0.5,
                      .(n_gwas = uniqueN(gwas_name)), by = gene]

gene_scatter <- merge(gene_best,
                      gene_anc[, .(best_gene = gene, ancestry_class)],
                      by = "best_gene", all.x = TRUE)
gene_scatter <- merge(gene_scatter, gene_n_gwas,
                      by.x = "best_gene", by.y = "gene", all.x = TRUE)
gene_scatter[is.na(n_gwas), n_gwas := 0L]

# Ancestry color palette: 3 categories
anc_colors <- c(
  "Cross-ancestry" = "#7B3294",   # purple — EUR + EAS converge
  "EUR"            = "#3B7DA5",   # cool blue
  "EAS"            = "#E07B39"    # warm orange
)

# Restrict to the high-confidence corner (top_var_pp >= 0.5 AND PP.H4 >= 0.5)
gene_corner <- gene_scatter[top_var_pp >= 0.5 & max_pp4 >= 0.5 &
                              !is.na(ancestry_class)]
gene_corner[, ancestry_class := factor(ancestry_class,
                                        levels = c("Cross-ancestry", "EUR", "EAS"))]
cat(sprintf("  [scatter] corner genes (top_var_pp>=0.5 & pp4>=0.5): %d\n",
            nrow(gene_corner)))

# Mild jitter at the 1.0 wall so coincident genes don't overplot perfectly
set.seed(42)
gene_corner[top_var_pp > 0.995, top_var_pp := 1 - runif(.N, 0, 0.018)]
gene_corner[max_pp4    > 0.995, max_pp4    := 1 - runif(.N, 0, 0.018)]

# ---- Build the scatter (zoomed to corner: x,y in [0.45, 1.02])
AX_LO <- 0.45
AX_HI <- 1.02

# Counts per ancestry class for legend labels
n_per_anc <- gene_corner[, .N, by = ancestry_class]
setkey(n_per_anc, ancestry_class)
anc_legend_labels <- sprintf("%s (n=%d)",
                              levels(gene_corner$ancestry_class),
                              n_per_anc[levels(gene_corner$ancestry_class), N])

# Size legend breaks: keep only values that occur in the data
size_breaks <- intersect(c(1, 2, 4, 6, 10),
                         seq(min(gene_corner$n_gwas), max(gene_corner$n_gwas)))
if (length(size_breaks) < 3) size_breaks <- pretty(gene_corner$n_gwas, n = 4)

# Minimal label set: 4 anchors illustrating quadrants (must be in 20-gene
# strict corner — i.e. fine-mapping coverage AND PIP, PP.H4 both >= 0.5).
#   - GGT1      : top-right  (canonical liver enzyme; PIP=1.0, PP.H4=0.997)
#   - HNF1A-AS1 : top-right  (HNF1A regulator; PIP=1.0, PP.H4=0.986)
#   - NT5C      : top-left   (high PP.H4 0.965, moderate PIP 0.634)
#   - TM6SF2    : bottom-right (canonical MASLD risk gene; PIP=0.957, PP.H4=0.520)
label_genes <- c("GGT1", "HNF1A-AS1", "NT5C", "TM6SF2")
gene_labels <- gene_corner[best_gene %in% label_genes]

p_sc <- ggplot() +
  geom_hline(yintercept = 0.9, linetype = "dashed",
             linewidth = 0.2, color = "gray55") +
  geom_vline(xintercept = 0.9, linetype = "dashed",
             linewidth = 0.2, color = "gray55") +
  geom_point(data = gene_corner,
             aes(x = top_var_pp, y = max_pp4,
                 color = ancestry_class, size = n_gwas),
             shape = 16, alpha = 0.88, stroke = 0) +
  geom_point(data = gene_corner,
             aes(x = top_var_pp, y = max_pp4, size = n_gwas),
             color = "white", stroke = 0.35, shape = 1, alpha = 0.45) +
  geom_text_repel(
    data = gene_labels,
    aes(x = top_var_pp, y = max_pp4, label = best_gene),
    color = "gray10",
    size = 2.6, fontface = "italic",
    box.padding = 0.5, point.padding = 0.3,
    segment.size = 0.2, segment.color = "gray45",
    min.segment.length = 0, force = 3, force_pull = 0.35,
    max.overlaps = Inf, seed = 42
  ) +
  scale_color_manual(values = anc_colors, name = "Gene replication",
                     labels = anc_legend_labels, drop = FALSE) +
  scale_size_continuous(name = "GWAS supporting\n(PP.H4 >= 0.5)",
                        range = c(1.6, 4.4),
                        breaks = size_breaks) +
  scale_x_continuous(limits = c(AX_LO, AX_HI),
                     breaks = c(0.5, 0.7, 0.9, 1.0),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(AX_LO, AX_HI),
                     breaks = c(0.5, 0.7, 0.9, 1.0),
                     expand = c(0, 0)) +
  labs(x = "SuSiE fine-mapping PIP  (COLOC lead variant)",
       y = "COLOC PP.H4  (best of SuSiE / ABF)",
       title = "High-confidence variant-to-gene COLOC genes") +
  theme_masld() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#EEEAE2", linewidth = 0.16),
    panel.border     = element_rect(color = "gray50", fill = NA,
                                    linewidth = 0.35),
    axis.text        = element_text(size = 6.5),
    axis.title       = element_text(size = 7),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.key.size  = unit(0.28, "cm"),
    legend.text      = element_text(size = 6),
    legend.title     = element_text(size = 6.5, face = "bold"),
    plot.title       = element_text(size = 8.5, face = "bold")
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.6, alpha = 0.9,
                                                   shape = 16)),
         size  = guide_legend(override.aes = list(color = "gray30",
                                                   alpha = 0.85)))

out <- file.path(PANEL_DIR, "coloc_deg_hybrid.pdf")
# RETIRED 2026-06-12 (not a Fig 2 panel): coloc_deg_hybrid.pdf
# save_fig(p_sc, out, width = fig_full_width * 0.6, height = 3.6)
# cat("[fig3a-hybrid] Saved:", out, "\n")

# Sidecar CSV: one row per corner gene
corner_export <- gene_corner[, .(best_gene, top_var_pp, max_pp4,
                                  best_gwas, ancestry_class, n_gwas)
                            ][order(-max_pp4, -top_var_pp)]
# RETIRED 2026-06-12 (not a Fig 2 panel): coloc_deg_hybrid_corner_genes.csv
# fwrite(corner_export,
#        file.path(PANEL_DIR, "coloc_deg_hybrid_corner_genes.csv"))
# cat("[fig3a-hybrid] Wrote sidecar CSV with", nrow(corner_export), "corner genes\n")

# ===========================================================================
# Version B : ALL genes plotted (axes from 0 to 1), labels limited to the
# top 5 per quadrant by combined evidence + 6 must-have anchors.
#
# A gene's color = "Cross-ancestry" if PP.H4 >= 0.5 in BOTH EUR and EAS GWAS;
# otherwise the ancestry of the best-pp4 GWAS (so low-evidence genes still
# get colored by which ancestry's GWAS gave their strongest, if weak, signal).
# ===========================================================================
gene_full <- copy(gene_scatter)
gene_best_anc <- sc_pos[sc_pos[, .I[which.max(pp4)], by = gene]$V1,
                        .(best_gene = gene,
                          best_anc  = fifelse(gwas_name %in% BBJ_5, "EAS", "EUR"))]
gene_full <- merge(gene_full, gene_best_anc, by = "best_gene", all.x = TRUE)
cross_genes <- gene_anc[ancestry_class == "Cross-ancestry", gene]
gene_full[, ancestry_class := fifelse(best_gene %in% cross_genes,
                                       "Cross-ancestry", best_anc)]
gene_full[, ancestry_class := factor(ancestry_class,
                                      levels = c("Cross-ancestry", "EUR", "EAS"))]
gene_full <- gene_full[!is.na(ancestry_class) & !is.na(top_var_pp) & !is.na(max_pp4)]

# Mild jitter only at the upper wall
set.seed(42)
gene_full[top_var_pp > 0.995, top_var_pp := 1 - runif(.N, 0, 0.018)]
gene_full[max_pp4    > 0.995, max_pp4    := 1 - runif(.N, 0, 0.018)]

# Quadrants split at 0.5 on each axis
gene_full[, quadrant := fcase(
  top_var_pp >= 0.5 & max_pp4 >= 0.5, "TR",
  top_var_pp <  0.5 & max_pp4 >= 0.5, "TL",
  top_var_pp <  0.5 & max_pp4 <  0.5, "BL",
  top_var_pp >= 0.5 & max_pp4 <  0.5, "BR"
)]

# Top-5 per quadrant by combined evidence (max_pp4 + top_var_pp)
quadrant_top5 <- gene_full[, .SD[order(-(max_pp4 + top_var_pp))][1:min(5, .N)],
                            by = quadrant]

# Must-have anchors (always labelled even if outside top-5).
# THRB, MLIP, RORC dropped: their COLOC top_snp has no fine-mapping coverage,
# so they cannot be plotted on the SuSiE-PIP axis.
must_have    <- c("PNPLA3", "TM4SF4", "RORA")
labels_b     <- gene_full[best_gene %in% unique(c(must_have,
                                                  quadrant_top5$best_gene))]

cat(sprintf("  [scatter B] %d genes plotted, %d labelled\n",
            nrow(gene_full), nrow(labels_b)))

# Size legend breaks for full data
size_breaks_b <- intersect(c(1, 2, 4, 6, 10),
                           seq(min(gene_full$n_gwas), max(gene_full$n_gwas)))
if (length(size_breaks_b) < 3) size_breaks_b <- pretty(gene_full$n_gwas, n = 4)

p_sc_b <- ggplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             linewidth = 0.2, color = "gray55") +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             linewidth = 0.2, color = "gray55") +
  # Background haze: genes with either axis < 0.5 → tiny gray dots
  ggrastr::rasterise(
    geom_point(data = gene_full[!(top_var_pp >= 0.5 & max_pp4 >= 0.5)],
               aes(x = top_var_pp, y = max_pp4),
               color = "gray70", shape = 16, size = 0.45, alpha = 0.45,
               stroke = 0),
    dpi = 350) +
  # Corner genes (both axes ≥ 0.5): small colored dots, ancestry-encoded
  ggrastr::rasterise(
    geom_point(data = gene_full[top_var_pp >= 0.5 & max_pp4 >= 0.5],
               aes(x = top_var_pp, y = max_pp4, color = ancestry_class),
               shape = 16, size = 0.95, alpha = 0.9, stroke = 0),
    dpi = 350) +
  geom_text_repel(
    data = labels_b,
    aes(x = top_var_pp, y = max_pp4, label = best_gene),
    color = "gray10",
    size = 2.4, fontface = "italic",
    box.padding   = 0.45,
    point.padding = 0.25,
    segment.size  = 0.18,
    segment.color = "gray45",
    min.segment.length = 0,
    force        = 5,
    force_pull   = 0.3,
    max.overlaps = Inf,
    max.iter     = 50000,
    max.time     = 5,
    seed         = 42
  ) +
  scale_color_manual(values = anc_colors, name = "Gene replication",
                     drop = FALSE) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  coord_cartesian(xlim = c(-0.02, 1.02), ylim = c(-0.02, 1.02),
                  expand = FALSE, clip = "off") +
  labs(x = "SuSiE fine-mapping PIP  (COLOC lead variant)",
       y = "COLOC PP.H4  (best of SuSiE / ABF)",
       title = "All COLOC genes with fine-mapping coverage") +
  theme_masld() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#EEEAE2", linewidth = 0.16),
    panel.border     = element_rect(color = "gray50", fill = NA,
                                    linewidth = 0.35),
    axis.text        = element_text(size = 6.5),
    axis.title       = element_text(size = 7),
    legend.position  = "right",
    legend.key.size  = unit(0.28, "cm"),
    legend.text      = element_text(size = 6),
    legend.title     = element_text(size = 6.5, face = "bold"),
    plot.title       = element_text(size = 8.5, face = "bold"),
    plot.margin      = margin(t = 10, r = 14, b = 8, l = 14)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.4, alpha = 0.95,
                                                   shape = 16)))

out_b <- file.path(PANEL_DIR, "coloc_deg_hybrid_b.pdf")
# RETIRED 2026-06-12 (not a Fig 2 panel): coloc_deg_hybrid_b.pdf
# save_fig(p_sc_b, out_b, width = fig_full_width * 0.85, height = 5.0)
# cat("[fig3a-hybrid] Saved version B (full-range, quadrant labels):", out_b, "\n")
