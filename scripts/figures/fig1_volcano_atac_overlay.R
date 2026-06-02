#!/usr/bin/env Rscript
# fig1_volcano_atac_overlay.R — Augmented Fig 1 volcano with ATAC validation tiers.
#
# Mirrors the data loading + thresholds of fig1_volcano.R but adds a secondary
# colour/shape encoding to highlight ATAC + COLOC validation status for Tier 1
# DEGs. Liang aesthetic: hollow circles for the cloud, filled gold disks for
# the 4-way validated TFs, threshold lines as thin gray dashes.
#
# Tier table (atac_tier):
#   FourWay     gold      #FFB300  THRB / HNF4A / RORA / MLXIPL / MAX (filled)
#   COLOC_ATAC  magenta   #C9265E  DEGs that are COLOC (PP.H4.susie>0.5) and
#                                  have a GWAS-ATAC variant within +/-50 kb
#   COLOC_only  light mag #FCE4EC  DEGs that are COLOC but no ATAC evidence
#   Conserved   teal      #00695C  Mouse+human concordant (is_conserved)
#   Other_DEG   gray      #9E9E9E  Remaining Tier 1 DEGs
#   n.s.        v light   #EEEEEE  background (padj>=0.05 or |LFC|<=0.5)
#
# Output: figures/main/fig1_atlas_overview/panels/fig1_volcano_atac_overlay.pdf
# Sibling of fig1_volcano.pdf — does NOT replace the canonical volcano.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG1_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

PADJ_CUT  <- 0.05
LFC_CUT   <- 0.5
COLOC_CUT <- 0.5         # SuSiE PP.H4
ATAC_WIN  <- 50000       # +/- 50 kb of any GWAS-ATAC variant in a peak

# 4-way validated TFs (hardcoded headline panel from GWAS-ATAC pipeline).
FOUR_WAY_TFS <- c("THRB", "HNF4A", "RORA", "MLXIPL", "MAX")

# ----------------------------------------------------------------------------
# Load: dream DEG table
# ----------------------------------------------------------------------------
message("Loading dream mega-analysis results...")
dream <- load_dream_results()
volc  <- dream[!is.na(dream_padj) & !is.na(dream_logFC),
               .(symbol, dream_logFC, dream_padj)]
volc[, neglog10p := -log10(dream_padj)]
volc[, is_sig := dream_padj < PADJ_CUT & abs(dream_logFC) > LFC_CUT]

# ----------------------------------------------------------------------------
# Load: gene-level SuSiE-COLOC (PP.H4.susie > 0.5)
# ----------------------------------------------------------------------------
message("Loading gene-level SuSiE-COLOC...")
coloc_f <- file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
coloc   <- fread(coloc_f)
# coerce numeric (the file has occasional empty strings)
coloc[, coloc_best_susie_pp4 := suppressWarnings(as.numeric(coloc_best_susie_pp4))]
coloc[, coloc_best_pp4       := suppressWarnings(as.numeric(coloc_best_pp4))]
# A gene "is COLOC" if EITHER SuSiE or ABF posterior > cutoff (matches atlas
# semantics — both posteriors preserved per Phase L docs).
coloc_genes <- unique(coloc[
  (!is.na(coloc_best_susie_pp4) & coloc_best_susie_pp4 > COLOC_CUT) |
  (!is.na(coloc_best_pp4)       & coloc_best_pp4       > COLOC_CUT),
  gene])
message(sprintf("  %d genes with COLOC PP4>%.2f (SuSiE or ABF)",
                length(coloc_genes), COLOC_CUT))

# ----------------------------------------------------------------------------
# Load: GWAS-ATAC variant annotation → genes near peak variants (+/- 50 kb)
# ----------------------------------------------------------------------------
message("Loading GWAS-ATAC variant annotation...")
atac_f <- file.path(BASE,
                    "GWAS/finemapping/results/gwas_atac/gwas_atac_variant_annotation.csv")
atac   <- fread(atac_f)
atac[, distance_to_tss := suppressWarnings(as.numeric(distance_to_tss))]
# Linked gene (SCENIC+/Cicero) is the primary mechanism call; nearest_gene
# within 50 kb of a peak-resident variant is a positional fallback.
atac_near_genes <- unique(c(
  atac[linked_gene  != "" & !is.na(linked_gene),  linked_gene],
  atac[nearest_gene != "" & !is.na(nearest_gene) &
       !is.na(distance_to_tss) & abs(distance_to_tss) <= ATAC_WIN,
       nearest_gene]
))
atac_near_genes <- atac_near_genes[atac_near_genes != ""]
message(sprintf("  %d unique genes within +/-%d kb of a GWAS-ATAC peak variant",
                length(atac_near_genes), ATAC_WIN %/% 1000))

# ----------------------------------------------------------------------------
# Load: conserved (mouse+human) from multi-evidence atlas
# ----------------------------------------------------------------------------
message("Loading multi-evidence atlas (is_conserved)...")
me <- load_multi_evidence()
conserved_genes <- if (!is.null(me) && "is_conserved" %in% names(me)) {
  unique(me[is_conserved %in% c(TRUE, "TRUE"), symbol])
} else character(0)
message(sprintf("  %d conserved-concordant genes (mouse+human)",
                length(conserved_genes)))

# ----------------------------------------------------------------------------
# Assign atac_tier per gene (Tier 1 DEGs only get the non-grey colours)
# ----------------------------------------------------------------------------
volc[, atac_tier := "n.s."]
volc[is_sig == TRUE, atac_tier := "Other_DEG"]
volc[is_sig == TRUE & symbol %in% conserved_genes, atac_tier := "Conserved"]
volc[is_sig == TRUE & symbol %in% coloc_genes,     atac_tier := "COLOC_only"]
volc[is_sig == TRUE & symbol %in% coloc_genes & symbol %in% atac_near_genes,
     atac_tier := "COLOC_ATAC"]
# 4-way TFs are master-regulator hits: relax LFC requirement (small effects in
# bulk are expected for TFs); require dream_padj < PADJ_CUT only.
volc[symbol %in% FOUR_WAY_TFS &
     !is.na(dream_padj) & dream_padj < PADJ_CUT, atac_tier := "FourWay"]

# Order so background plots first, headline tiers plot on top
tier_levels <- c("n.s.", "Other_DEG", "Conserved", "COLOC_only", "COLOC_ATAC", "FourWay")
volc[, atac_tier := factor(atac_tier, levels = tier_levels)]
setorder(volc, atac_tier)

tier_counts <- volc[, .N, by = atac_tier]
message("\nPer-tier counts:")
for (i in seq_len(nrow(tier_counts))) {
  message(sprintf("  %-12s %s", tier_counts$atac_tier[i],
                  comma(tier_counts$N[i])))
}

# 4-way TF audit (some may be n.s. in the dream table — flag them)
four_way_status <- volc[symbol %in% FOUR_WAY_TFS,
                        .(symbol, dream_logFC, dream_padj, atac_tier)]
message("\n4-way TF Tier-1 audit:")
for (tf in FOUR_WAY_TFS) {
  row <- four_way_status[symbol == tf]
  if (nrow(row) == 0) {
    message(sprintf("  %-8s NOT in dream table", tf))
  } else {
    message(sprintf("  %-8s logFC=%+.2f  padj=%.2e  tier=%s",
                    tf, row$dream_logFC, row$dream_padj, as.character(row$atac_tier)))
  }
}

# ----------------------------------------------------------------------------
# Labels: 4-way TFs (all that are present in dream) + top 5 conserved by |LFC|
# ----------------------------------------------------------------------------
label_four <- volc[symbol %in% FOUR_WAY_TFS & symbol != "" & !grepl("^ENSG", symbol)]
cons_pool  <- volc[atac_tier == "Conserved" & symbol != "" & !grepl("^ENSG", symbol)]
setorder(cons_pool, -dream_logFC)   # largest +LFC
top_cons_up   <- head(cons_pool, 3)
setorder(cons_pool, dream_logFC)    # largest -LFC
top_cons_down <- head(cons_pool, 2)
label_cons  <- rbind(top_cons_up, top_cons_down)
label_df    <- unique(rbind(label_four, label_cons), by = "symbol")
message(sprintf("\nLabelling %d genes (%d 4-way TFs + %d conserved top-LFC)",
                nrow(label_df), nrow(label_four), nrow(label_cons)))

# ----------------------------------------------------------------------------
# Aesthetics — Liang palette
# ----------------------------------------------------------------------------
tier_colors <- c(
  "FourWay"    = "#FFB300",   # gold
  "COLOC_ATAC" = "#C9265E",   # warm magenta
  "COLOC_only" = "#FCE4EC",   # light magenta (light enough; bordered black)
  "Conserved"  = "#00695C",   # deep teal
  "Other_DEG"  = "#9E9E9E",   # neutral gray
  "n.s."       = "#EEEEEE"    # very light gray
)
tier_labels <- c(
  "FourWay"    = "4-way validated TF (n=5)",
  "COLOC_ATAC" = "COLOC + ATAC peak",
  "COLOC_only" = "COLOC only",
  "Conserved"  = "Conserved (mouse+human)",
  "Other_DEG"  = "Other Tier 1 DEG",
  "n.s."       = "n.s."
)

x_lim <- max(abs(volc$dream_logFC), na.rm = TRUE) * 1.04
y_lim <- max(volc$neglog10p,        na.rm = TRUE) * 1.05

# Per-tier styling: shape, point size, stroke, fill-alpha control via aes mapping.
# Map both colour (ring) and fill (interior) to atac_tier. shape is selected
# per tier via scale_shape_manual. n.s. and Other_DEG use shape 16 (solid),
# Conserved uses shape 1 (hollow), COLOC_only/COLOC_ATAC/FourWay use shape 21
# (filled with ring).
tier_shapes <- c(
  "FourWay"    = 21,
  "COLOC_ATAC" = 21,
  "COLOC_only" = 21,
  "Conserved"  = 1,
  "Other_DEG"  = 1,
  "n.s."       = 16
)
# Ring colour (used by shape 21 stroke). For shape 1 / 16 ggplot uses `colour`.
ring_colors <- c(
  "FourWay"    = "black",
  "COLOC_ATAC" = "black",
  "COLOC_only" = "black",
  "Conserved"  = "#00695C",
  "Other_DEG"  = "#9E9E9E",
  "n.s."       = "#EEEEEE"
)

p <- ggplot() +
  # threshold lines
  geom_hline(yintercept = -log10(PADJ_CUT),
             linetype = "dashed", colour = "gray70", linewidth = 0.25) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT),
             linetype = "dashed", colour = "gray70", linewidth = 0.25) +
  # n.s. cloud — rasterised
  rasterize_layer(
    geom_point(data = volc[atac_tier == "n.s."],
               aes(x = dream_logFC, y = neglog10p, colour = atac_tier),
               shape = 16, size = 0.35, alpha = 0.5, inherit.aes = FALSE)
  ) +
  # Other DEGs — hollow gray
  rasterize_layer(
    geom_point(data = volc[atac_tier == "Other_DEG"],
               aes(x = dream_logFC, y = neglog10p, colour = atac_tier),
               shape = 1, size = 0.75, stroke = 0.25, alpha = 0.85,
               inherit.aes = FALSE)
  ) +
  # Conserved — hollow teal
  geom_point(data = volc[atac_tier == "Conserved"],
             aes(x = dream_logFC, y = neglog10p, colour = atac_tier),
             shape = 1, size = 0.95, stroke = 0.30, alpha = 0.9,
             inherit.aes = FALSE) +
  # COLOC only — filled light magenta with black ring
  geom_point(data = volc[atac_tier == "COLOC_only"],
             aes(x = dream_logFC, y = neglog10p, fill = atac_tier),
             shape = 21, size = 1.30, stroke = 0.30, colour = "black",
             alpha = 0.95, inherit.aes = FALSE) +
  # COLOC + ATAC — filled magenta with black ring (headline tier)
  geom_point(data = volc[atac_tier == "COLOC_ATAC"],
             aes(x = dream_logFC, y = neglog10p, fill = atac_tier),
             shape = 21, size = 1.65, stroke = 0.35, colour = "black",
             inherit.aes = FALSE) +
  # 4-way TFs — large filled gold
  geom_point(data = volc[atac_tier == "FourWay"],
             aes(x = dream_logFC, y = neglog10p, fill = atac_tier),
             shape = 21, size = 2.50, stroke = 0.45, colour = "black",
             inherit.aes = FALSE) +
  # Dummy invisible layer to ensure all 5 legend tiers register on BOTH scales.
  geom_point(data = data.table(
               dream_logFC = rep(NA_real_, 5),
               neglog10p   = rep(NA_real_, 5),
               atac_tier   = factor(c("FourWay","COLOC_ATAC","COLOC_only",
                                      "Conserved","Other_DEG"),
                                    levels = tier_levels)),
             aes(x = dream_logFC, y = neglog10p,
                 colour = atac_tier, fill = atac_tier),
             shape = 21, size = 1.4, stroke = 0.35, na.rm = TRUE,
             inherit.aes = FALSE, show.legend = TRUE) +
  # Labels
  geom_text_repel(data = label_df,
                  aes(x = dream_logFC, y = neglog10p, label = symbol),
                  inherit.aes = FALSE,
                  size = 2.4, colour = "black", fontface = "italic",
                  segment.size = 0.2, segment.colour = "gray45",
                  box.padding = 0.45, point.padding = 0.25,
                  min.segment.length = 0,
                  max.overlaps = Inf, force = 4, seed = 42,
                  show.legend = FALSE) +
  scale_colour_manual(values = ring_colors, guide = "none") +
  scale_fill_manual(values = tier_colors,
                    breaks = c("FourWay", "COLOC_ATAC", "COLOC_only",
                               "Conserved", "Other_DEG"),
                    labels = c("4-way validated TF (n=5)",
                               "COLOC + ATAC peak",
                               "COLOC only",
                               "Conserved (mouse+human)",
                               "Other Tier 1 DEG"),
                    name = NULL,
                    guide = guide_legend(override.aes = list(
                      shape  = c(21, 21, 21, 1, 1),
                      fill   = c("#FFB300", "#C9265E", "#FCE4EC", NA, NA),
                      colour = c("black", "black", "black", "#00695C", "#9E9E9E"),
                      size   = c(2.2, 1.6, 1.3, 1.1, 0.9),
                      stroke = c(0.45, 0.35, 0.30, 0.40, 0.30)
                    ))) +
  scale_x_continuous(limits = c(-x_lim, x_lim),
                     expand = expansion(mult = 0),
                     breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(limits = c(0, y_lim),
                     expand = expansion(mult = c(0, 0)),
                     breaks = pretty_breaks(n = 5)) +
  labs(
    x = expression(bold("log"[2]*" fold change")),
    y = expression(bold(-log[10]~"(adjusted P)"))
  ) +
  theme_masld(base_size = 7) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(linewidth = 0.18, colour = "gray92"),
    axis.title         = element_text(size = 8, face = "bold", colour = "black"),
    legend.position    = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background  = element_blank(),
    legend.key         = element_blank(),
    legend.key.size    = unit(0.30, "cm"),
    legend.text        = element_text(size = 6),
    legend.spacing.y   = unit(0.05, "cm"),
    plot.margin        = margin(4, 6, 2, 4)
  )

out_pdf <- file.path(PANEL_DIR, "fig1_volcano_atac_overlay.pdf")
save_fig(p, out_pdf, width = fig_half_width * 1.15, height = 3.2)

# Companion CSV — per-tier gene table for downstream provenance.
fwrite(volc[atac_tier != "n.s.",
            .(symbol, dream_logFC, dream_padj, neglog10p, atac_tier)],
       file.path(PANEL_DIR, "fig1_volcano_atac_overlay_tiers.csv"))
fwrite(tier_counts,
       file.path(PANEL_DIR, "fig1_volcano_atac_overlay_counts.csv"))

if (file.exists(out_pdf)) {
  message(sprintf("\nOutput: %s (%s)", out_pdf,
                  utils:::format.object_size(file.size(out_pdf), "auto")))
}
