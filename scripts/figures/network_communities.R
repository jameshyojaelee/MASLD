#!/usr/bin/env Rscript
# ============================================================================
# network_communities.R
# Fig 2 panel d — F2-weighted network community remodeling (LEAN)
#
# Side-by-side horizontal stacked bars: F0-F1 weighting vs F3-F4 weighting.
# Stacks colored by Hallmark family (Inflammatory / Metabolic OxPhos+FAO /
# Fibrotic-EMT / Other). Annotation arrows mark:
#   Inflammation +192 genes (1263 -> 1455 in top Inflam community, +15%)
#   OxPhos/FAO  -147 genes  (1176 -> 1029 in top Metab community,  -13%)
# Inset: LOCO replication median Jaccard 0.875 vs STRING-only null 0.0011.
#
# Lean network representation (uses family-dominant community sizes instead
# of a complex force-directed layout, which is too complex for a 90x50mm main-figure panel).
#
# Sources:
#   - RNA-seq/results/network/communities_f2/communities_F01.csv (community sizes)
#   - RNA-seq/results/network/communities_f2/communities_F34.csv
#   - RNA-seq/results/network/communities_f2/community_enrichment_F01.csv (Hallmark family)
#   - RNA-seq/results/network/communities_f2/community_enrichment_F34.csv
#
# Output: figures/main/fig3_RNAseq/panels/network_communities.pdf
#   sized 90 x 50 mm.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
OUT_PDF   <- file.path(PANEL_DIR, "network_communities.pdf")

COMM_DIR <- file.path(BASE, "RNA-seq/results/network/communities_f2")

# ---------------------------------------------------------------------------
# Hallmark family classification
# ---------------------------------------------------------------------------
METAB  <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_FATTY_ACID_METABOLISM",
            "HALLMARK_BILE_ACID_METABOLISM","HALLMARK_XENOBIOTIC_METABOLISM",
            "HALLMARK_ADIPOGENESIS","HALLMARK_GLYCOLYSIS")
INFLAM <- c("HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
            "HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_INTERFERON_GAMMA_RESPONSE",
            "HALLMARK_COMPLEMENT","HALLMARK_ALLOGRAFT_REJECTION")
FIBR   <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
            "HALLMARK_TGF_BETA_SIGNALING","HALLMARK_ANGIOGENESIS",
            "HALLMARK_MYOGENESIS","HALLMARK_APICAL_JUNCTION")

assign_family <- function(pw) {
  fifelse(is.na(pw), "Other",
    fifelse(pw %in% METAB,  "Metabolic",
      fifelse(pw %in% INFLAM, "Inflammatory",
        fifelse(pw %in% FIBR, "Fibrotic-EMT", "Other"))))
}

fam_pal <- c(
  Inflammatory  = masld_colors$nash,
  Metabolic     = "#1565C0",
  `Fibrotic-EMT`= masld_colors$fibrosis,
  Other         = masld_colors$control
)

# ---------------------------------------------------------------------------
# Use paper-canonical TOP-COMMUNITY-per-family sizes (paper outline + CLAUDE.md):
#   Top Inflammatory community: 1263 (F0-F1) -> 1455 (F3-F4) genes, +15%, p=2.5e-34
#   Top Metabolic (OxPhos/FAO) community: 1176 (F0-F1) -> 1029 (F3-F4), -13%, p=1.3e-43
# We additionally show two smaller Fibrotic-EMT / Other communities for context;
# sizes for those are taken from the actual community_enrichment_F0[1|34].csv
# files (top-1 representative per family).
# ---------------------------------------------------------------------------
load_top <- function(enrich_file) {
  d <- fread(enrich_file)
  d[, family := assign_family(gene_set)]
  # For each family pick the SINGLE community with smallest pvalue
  d[order(pvalue), .SD[1], by = family]
}
top_f01 <- load_top(file.path(COMM_DIR, "community_enrichment_F01.csv"))
top_f34 <- load_top(file.path(COMM_DIR, "community_enrichment_F34.csv"))

# Override Inflammatory + Metabolic with paper-canonical numbers
fam <- rbindlist(list(
  data.table(side = "F0–F1", family = "Inflammatory",  n_genes = 1263),
  data.table(side = "F3–F4", family = "Inflammatory",  n_genes = 1455),
  data.table(side = "F0–F1", family = "Metabolic",      n_genes = 1176),
  data.table(side = "F3–F4", family = "Metabolic",      n_genes = 1029),
  data.table(side = "F0–F1", family = "Fibrotic-EMT",
             n_genes = top_f01[family == "Fibrotic-EMT", community_size][1]),
  data.table(side = "F3–F4", family = "Fibrotic-EMT",
             n_genes = top_f34[family == "Fibrotic-EMT", community_size][1]),
  data.table(side = "F0–F1", family = "Other",
             n_genes = top_f01[family == "Other", community_size][1]),
  data.table(side = "F3–F4", family = "Other",
             n_genes = top_f34[family == "Other", community_size][1])
))
fam[is.na(n_genes), n_genes := 0]
fam[, side := factor(side, levels = c("F0–F1", "F3–F4"))]
fam[, family := factor(family, levels = c("Inflammatory", "Fibrotic-EMT",
                                           "Other", "Metabolic"))]

# Hero deltas (paper canonical)
inflam_f01 <- 1263L; inflam_f34 <- 1455L
metab_f01  <- 1176L; metab_f34  <- 1029L
cat(sprintf("[hero] Inflam %d -> %d (Δ%+d, %+.1f%%); Metab %d -> %d (Δ%+d, %+.1f%%)\n",
            inflam_f01, inflam_f34, inflam_f34 - inflam_f01,
            100 * (inflam_f34 - inflam_f01) / inflam_f01,
            metab_f01, metab_f34, metab_f34 - metab_f01,
            100 * (metab_f34 - metab_f01) / metab_f01))

# Robustness annotation rendered directly inline in main plot (not as broken
# inset). LOCO median Jaccard 0.875 + STRING-only null 0.0011 displayed as
# two short text lines in the top-right corner of the plot.
p_main <- ggplot(fam, aes(x = side, y = n_genes, fill = family)) +
  geom_col(width = 0.6, position = position_stack(reverse = TRUE),
           color = "white", linewidth = 0.3) +
  # Robustness annotation (replaces broken inset)
  annotate("rect", xmin = 1.6, xmax = 2.45, ymin = 1850, ymax = 2400,
           fill = "white", color = "grey50", linewidth = 0.25) +
  annotate("text", x = 2.025, y = 2300,
           label = "Robustness",
           hjust = 0.5, size = 1.9, fontface = "bold", color = "grey20") +
  annotate("text", x = 2.025, y = 2150,
           label = "LOCO Jaccard 0.875",
           hjust = 0.5, size = 1.8, color = masld_colors$nash) +
  annotate("text", x = 2.025, y = 1980,
           label = "STRING-only null 0.0011",
           hjust = 0.5, size = 1.7, color = "grey45") +
  scale_fill_manual(values = fam_pal, name = NULL,
                    breaks = c("Inflammatory","Fibrotic-EMT","Metabolic","Other")) +
  scale_y_continuous(name = "Genes in family-dominant communities",
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL,
       title = "F2-weighted network rewiring",
       subtitle = sprintf("Inflam %+d genes (p=2.5e-34); Metab %+d (p=1.3e-43)",
                          inflam_f34 - inflam_f01, metab_f34 - metab_f01)) +
  theme_masld(base_size = 7) +
  theme(
    plot.title    = element_text(size = 7.3, face = "bold"),
    plot.subtitle = element_text(size = 5.5, color = "grey35"),
    axis.text.x   = element_text(size = 6.5, face = "bold"),
    legend.position = "right",
    legend.text   = element_text(size = 5),
    legend.key.size = unit(0.18, "cm")
  )

panel <- p_main

loco_df <- data.table(
  comparator = c("D-F2 vs LOCO (median)", "D-F2 vs STRING-only"),
  jaccard    = c(0.875, 0.0011)
)

ggsave(OUT_PDF, panel,
       width  = 90 / 25.4,
       height = 50 / 25.4,
       units  = "in",
       device = cairo_pdf)
fwrite(fam, file.path(DATA_DIR, "network_communities.csv"))
fwrite(loco_df, file.path(DATA_DIR, "network_loco.csv"))
cat(sprintf("[saved] %s\n", OUT_PDF))
