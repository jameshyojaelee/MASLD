#!/usr/bin/env Rscript
##############################################################################
# fig3g_cross_ancestry_pp4.R
#
# Fig 3g — Cross-ancestry SuSiE-COLOC concordance.
#
# KEY MESSAGE: A core set of MASLD eQTL targets colocalize in BOTH European
#              and East Asian GWAS at PP.H4 ≥ 0.5 — the most portable
#              regulatory signals in the atlas. Single-ancestry replicators
#              (EUR-only / EAS-only) flank the corner.
#
# Per gene, computes:
#   x = max PP.H4 across 14 EUR GWAS
#   y = max PP.H4 across 3 EAS GWAS
# Each axis is gene-specific (per-gene COLOC output) — no locus-sharing.
#
# Output (PDF only, flat in FIG3_DIR):
#   cross_ancestry_pp4.pdf
#   cross_ancestry_pp4_table.csv  (sidecar: cross-ancestry replicators)
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(ggrepel)
  library(ggrastr)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ---------------------------------------------------------------------------
# GWAS portfolio — 14 EUR + 3 EAS
# ---------------------------------------------------------------------------
EUR_14 <- c(
  "2019_31311600_NAFLD_EUR", "2020_32298765_NAFLD_EUR",
  "2021_34128465_PDFF_EUR",  "2021_34841290_NAFLD_EUR",
  "2021_34957434_PDFF_EUR",  "2022_36402844_PDFF_EUR",
  "2023_36280732_NAFLD_deCode_EUR", "2023_36280732_NAFLD_Intermountain_EUR",
  "2023_36280732_NAFLD_UKBB_EUR", "FinnGen_NAFLD", "FinnGen_NASH",
  "UKBB_ALT", "UKBB_AST", "UKBB_GGT"
)
BBJ_3 <- c("BBJ_ALT", "BBJ_AST", "BBJ_GGT")
TARGET_17 <- c(EUR_14, BBJ_3)

# Colorblind-safe categorical colors (Set2-derived)
COL_EUR  <- "#3B7DA5"   # cool blue
COL_EAS  <- "#E07B39"   # warm orange
COL_BOTH <- "#7B3294"   # purple

# ===========================================================================
# Load data
# ===========================================================================
sc <- fread(file.path(BASE,
              "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"),
            select = c("gene", "gwas_name", "PP.H4.abf", "PP.H4.susie"))
sc <- sc[gwas_name %in% TARGET_17]
sc[, ancestry := fifelse(gwas_name %in% BBJ_3, "EAS", "EUR")]
sc[, pp4_any  := pmax(PP.H4.abf, PP.H4.susie, na.rm = TRUE)]

gl <- fread(file.path(BASE,
              "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"),
            select = c("gene", "is_mhc"))

# ===========================================================================
# Build per-gene EUR vs EAS table
# ===========================================================================
gene_anc <- sc[!is.na(pp4_any),
               .(max_pp4 = max(pp4_any, na.rm = TRUE)),
               by = .(gene, ancestry)]
gene_xa <- dcast(gene_anc, gene ~ ancestry, value.var = "max_pp4", fill = 0)
gene_xa <- gene_xa[!is.na(EUR) & !is.na(EAS)]

# Replication breadth: number of GWAS phenotypes hitting PP.H4 >= 0.5
n_phen <- sc[!is.na(pp4_any) & pp4_any >= 0.5,
             .(n_phen = uniqueN(gwas_name)), by = gene]
gene_xa <- merge(gene_xa, n_phen, by = "gene", all.x = TRUE)
gene_xa[is.na(n_phen), n_phen := 0L]

# Drop MHC (LD complexity makes COLOC unreliable there)
gene_xa <- merge(gene_xa, gl, by = "gene", all.x = TRUE)
gene_xa <- gene_xa[is_mhc == FALSE | is.na(is_mhc)]
gene_xa[, is_mhc := NULL]

gene_xa[, tier := fcase(
  EUR >= 0.5 & EAS >= 0.5, "both",
  EUR >= 0.5,              "EUR-only",
  EAS >= 0.5,              "EAS-only",
  default                  = "background"
)]

cat(sprintf("Genes tested in both ancestries: %d\n", nrow(gene_xa)))
cat(sprintf("  Cross-ancestry (both >= 0.5):  %d\n",
            sum(gene_xa$tier == "both")))
cat(sprintf("  EUR-only (PP.H4 >= 0.5):       %d\n",
            sum(gene_xa$tier == "EUR-only")))
cat(sprintf("  EAS-only (PP.H4 >= 0.5):       %d\n",
            sum(gene_xa$tier == "EAS-only")))

# Mild jitter at the 1.0 wall (visual clarity only)
gene_xa[EUR > 0.995, EUR := 1 - runif(.N, 0, 0.018)]
gene_xa[EAS > 0.995, EAS := 1 - runif(.N, 0, 0.018)]

# Top cross-ancestry replicators for labeling
gene_xa[, score := EUR + EAS]
both_lab <- gene_xa[tier == "both"][order(-score)][
  1:min(8, sum(gene_xa$tier == "both"))]

tier_colors <- c(
  "both"       = COL_BOTH,
  "EUR-only"   = COL_EUR,
  "EAS-only"   = COL_EAS,
  "background" = "#BFB9AF"
)

# ===========================================================================
# Plot
# ===========================================================================
p <- ggplot() +
  # Quadrant shading for the cross-ancestry corner
  geom_rect(aes(xmin = 0.5, xmax = 1.02, ymin = 0.5, ymax = 1.02),
            fill = "#FFF1CF", alpha = 0.55, inherit.aes = FALSE) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             linewidth = 0.2, color = "gray60") +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             linewidth = 0.2, color = "gray60") +
  # Background — no replication anywhere
  rasterize(
    geom_point(data = gene_xa[tier == "background"],
               aes(x = EUR, y = EAS),
               color = "#B5B0A8", size = 0.4, alpha = 0.25, stroke = 0),
    dpi = 600
  ) +
  # Single-ancestry replicators
  geom_point(data = gene_xa[tier %in% c("EUR-only", "EAS-only")],
             aes(x = EUR, y = EAS, color = tier),
             size = 1.4, alpha = 0.85, stroke = 0) +
  # Cross-ancestry replicators — sized by replication breadth
  geom_point(data = gene_xa[tier == "both"],
             aes(x = EUR, y = EAS, color = tier, size = n_phen),
             alpha = 0.92, stroke = 0) +
  # Subtle white outline on corner points for crispness
  geom_point(data = gene_xa[tier == "both"],
             aes(x = EUR, y = EAS, size = n_phen),
             color = "white", stroke = 0.4, shape = 1, alpha = 0.5) +
  # Labels — top 8 cross-ancestry replicators
  geom_label_repel(
    data = both_lab,
    aes(x = EUR, y = EAS, label = gene),
    color = "gray10", fill = alpha("white", 0.85),
    size = 2.5, fontface = "italic",
    label.size = 0.18, label.padding = unit(0.12, "lines"),
    label.r = unit(0.08, "lines"),
    box.padding = 0.6, point.padding = 0.3,
    segment.size = 0.22, segment.color = "gray35",
    min.segment.length = 0, force = 4, max.overlaps = Inf, seed = 42
  ) +
  scale_color_manual(values = tier_colors, name = "Replication",
                     breaks = c("both", "EUR-only", "EAS-only"),
                     labels = c("both",
                                sprintf("EUR-only (n=%d)",
                                        sum(gene_xa$tier == "EUR-only")),
                                sprintf("EAS-only (n=%d)",
                                        sum(gene_xa$tier == "EAS-only")))) +
  scale_size_continuous(range = c(1.2, 3.2),
                        name  = "# phenotypes\nPP.H4 >= 0.5",
                        breaks = c(1, 3, 5, 8)) +
  scale_x_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25),
                     expand = c(0, 0)) +
  labs(x = "max PP.H4  (14 EUR GWAS)",
       y = "max PP.H4  (3 EAS GWAS)",
       title = "Cross-ancestry colocalization concordance",
       subtitle = sprintf("%d genes replicate in both EUR and EAS at PP.H4 >= 0.5",
                          sum(gene_xa$tier == "both"))) +
  theme_masld() +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "gray50", fill = NA, linewidth = 0.35),
    axis.text        = element_text(size = 6.5),
    axis.title       = element_text(size = 7),
    legend.position  = "right",
    legend.title     = element_text(size = 6.5, face = "bold"),
    legend.text      = element_text(size = 6),
    legend.key.size  = unit(0.28, "cm"),
    plot.title       = element_text(size = 8.5, face = "bold"),
    plot.subtitle    = element_text(size = 7, color = "gray35")
  ) +
  guides(color = guide_legend(order = 1, override.aes = list(size = 2.2)),
         size  = guide_legend(order = 2,
                              override.aes = list(color = COL_BOTH)))

# RETIRED 2026-06-12 (cross_ancestry_pp4.pdf no longer a Fig 2 panel; sidecar CSV below is kept):
# out_pdf <- file.path(FIG3_DIR, "cross_ancestry_pp4.pdf")
# save_fig(p, out_pdf, width = fig_full_width * 0.6, height = 3.6)
# cat("Saved:", out_pdf, "\n")

# Sidecar CSV: cross-ancestry replicators for caption / supplement
out_tbl <- gene_xa[tier == "both",
                   .(gene, EUR_max_pp4 = round(EUR, 4),
                     EAS_max_pp4 = round(EAS, 4),
                     n_phen, score)
                  ][order(-score)]
fwrite(out_tbl,
       file.path(FIG3_DIR, "cross_ancestry_pp4_table.csv"))
cat("Wrote sidecar table with", nrow(out_tbl), "cross-ancestry replicators\n")
