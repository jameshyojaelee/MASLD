#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Fig 4e — "disease reorganizes liver space", ONE merged panel:
#
#  LEFT:  SVG dynamics bar chart — how many genes gained / maintained / lost
#         spatial structure in MASLD (GSE192741, donor-blocked Moran's I).
#
#  RIGHT: Moran's I scatter — Healthy (x) vs Steatotic (y) for all SVGs,
#         colored by category; hepatocyte-specific examples labeled.
#         Identity diagonal (y = x) marks no change; above = gained, below = lost.
#
# Data: Analysis/Spatial/results/svg/differential_svgs.csv
# Output: figures/main/fig4_validation/panels/fig4e_spatial_reorganization.pdf
# (relettered to 4e 2026-07-08: Fig4 rebuilt around "physical corroboration", then
# plasma-translation + spatial-CCC demoted to supp, shifting this panel 4f->4e —
# protein/chromatin layers now precede this panel in the lineup)
# Env:    rnaseq
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
SP <- file.path(BASE, "Analysis/Spatial/results")

# ── [a] SVG dynamics bar chart ────────────────────────────────────────────────
svg <- fread(file.path(SP, "svg/differential_svgs.csv"))
setnames(svg, 1, "gene")
tl <- function(x) tolower(as.character(x)) %in% c("true", "1", "t")
svg[, `:=`(h = tl(svg_healthy), m = tl(svg_masld))]
n_gained <- svg[m & !h, .N]; n_maint <- svg[h & m, .N]; n_lost <- svg[h & !m, .N]

svg_df <- data.table(
  cls  = factor(c("Gained", "Maintained", "Lost"),
                levels = c("Gained", "Maintained", "Lost")),
  n    = c(n_gained, n_maint, n_lost),
  fill = c(masld_colors$up, masld_colors$ns, masld_colors$down))

p_a <- ggplot(svg_df, aes(cls, n, fill = fill)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = formatC(n, big.mark = ",", format = "d")),
            vjust = -0.35, size = GEOM_TEXT_6PT) +
  scale_fill_identity() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = NULL, y = "spatially variable genes") +
  theme_masld_compact() +
  theme(axis.text   = element_text(color = house_ink),
        axis.text.x = element_text(color = house_ink),
        plot.margin = margin(2, 8, 0, 3))

# ── [b] Moran's I scatter: Healthy vs Steatotic ───────────────────────────────
# Clean category labels and colors
svg[, cat_lab := fcase(
  category == "disease_emergent_SVG", "Gained",
  category == "disease_lost_SVG",     "Lost",
  default = "Stable")]
svg[, cat_lab := factor(cat_lab, levels = c("Gained", "Stable", "Lost"))]

# Semantic hues kept (Liang magenta = disease-up, deep blue = down) for cross-panel
# consistency; prettiness comes from uniform small dots + tuned alpha + a lighter
# neutral stable gray, not from re-hueing.
cat_colors <- c(
  Gained = masld_colors$up,
  Stable = "#D9D9D9",
  Lost   = masld_colors$down)
cat_alpha  <- c(Gained = 0.80, Stable = 0.28, Lost = 0.80)
PT_SIZE    <- 0.5   # all dots identical size (no size mapping)

# Genes to label — chosen so every label is EITHER a true reorganizer (crossed the
# SVG threshold, i.e. a colored Gained/Lost dot) OR one of two maintained-zonation
# landmarks with strong cross-panel narrative. Rationale: emergent/lost SVGs are by
# construction low Moran's I (they pile into the bottom-left origin), so two high-I
# stable anchors give the labels spatial spread AND mark where zonation persists.
#   Maintained-zonation landmarks (Stable category, gray):
#     CYP3A4 (0.567→0.546, on diagonal)              — pericentral drug metabolism; own Fig4 zonation panel
#     FABP1  (0.337→0.212, below diagonal, weakened) — fatty-acid binding; Fig2E locus-zoom
#   Gained spatial structure (disease_emergent_SVG, magenta):
#     LDLR   (0.016→0.163) — LDL uptake
#     HMGCS1 (0.031→0.115) — cholesterol synthesis
#     COL1A1 (0.031→0.063) — fibrillar collagen; fibroblast ECM program (CCC in Fig 3)
#   Lost spatial structure (disease_lost_SVG, blue):
#     COL1A2 (0.121→0.037) — collagen partner of COL1A1 (the collagen pair gains/loses together)
hep_label <- c("CYP3A4","FABP1","LDLR","HMGCS1","COL1A1","COL1A2")
svg[, label := fifelse(gene %in% hep_label, gene, NA_character_)]

# Per-gene label nudge targets: fan labels into empty space — origin reorganizer
# cluster to the left gutter, the two stable landmarks up-right / lower-right for spread.
nudge_tbl <- data.table(
  gene = c("CYP3A4","FABP1","LDLR","HMGCS1","COL1A1","COL1A2"),
  nx   = c(-0.06,   0.06, -0.11, -0.13, -0.11,  0.05),
  ny   = c( 0.11,  -0.12,  0.10, -0.02, -0.10, -0.08))
svg[nudge_tbl, on = "gene", `:=`(nx = i.nx, ny = i.ny)]
lab_dt <- svg[!is.na(label)]

p_b <- ggplot(svg, aes(morans_I_healthy, morans_I_masld,
                        color = cat_lab, alpha = cat_lab)) +
  # identity line (no change)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60", linewidth = 0.3) +
  geom_point(shape = 16, size = PT_SIZE) +
  geom_text_repel(
    data = lab_dt,
    aes(label = label),
    nudge_x       = lab_dt$nx,
    nudge_y       = lab_dt$ny,
    size          = GEOM_TEXT_6PT,
    fontface      = "italic",
    color         = house_ink, # all labels identical: near-black ink, 6pt Helvetica italic
    alpha         = 1,         # override inherited alpha=cat_lab (fixes faded Stable labels)
    box.padding   = 0.35,
    point.padding = 0.20,
    segment.size  = 0.22,
    segment.color = "grey50",
    max.overlaps  = Inf,
    seed          = 1,
    force         = 2,
    force_pull    = 0.1,
    show.legend   = FALSE,
    min.segment.length = 0) +
  scale_color_manual(values = cat_colors, name = NULL) +
  scale_alpha_manual(values = cat_alpha,  name = NULL, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.16, 0.05))) +
  labs(x = "Moran's I (Healthy)", y = "Moran's I (Steatotic)") +
  coord_equal() +
  # legend swatches enlarged (plotted points are tiny) and at full opacity
  guides(color = guide_legend(override.aes = list(size = 1.4, alpha = 1))) +
  theme_masld_compact() +
  theme(legend.position  = c(0.16, 0.82),
        legend.key.size  = unit(0.22, "cm"),
        legend.background = element_rect(fill = NA, color = NA),
        axis.text = element_text(color = house_ink),
        plot.margin = margin(2, 4, 0, 4))

# ── compose ───────────────────────────────────────────────────────────────────
# NOTE: parenthesize (p_a | p_b) — `+` binds tighter than `|`, so without the
# parens plot_layout() would attach to p_b alone and the width ratio would never
# apply. p_b uses coord_equal() (square), so it is height-limited: give it the
# larger width share AND a taller canvas so the Moran's I scatter renders large.
fig <- (p_a | p_b) + plot_layout(widths = c(0.85, 1))
out <- file.path(FIG4_DIR, "panels", "fig4e_spatial_reorganization.pdf")
# Sized to the Fig4 layout slot (panel E, 3.46x1.88in; 2026-07-08).
save_fig(fig, out, width = 3.46, height = 1.88)
message("Saved: ", out)
message(sprintf("[spatial] SVG gained=%d maintained=%d lost=%d",
                n_gained, n_maint, n_lost))
