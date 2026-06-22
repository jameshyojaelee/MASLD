#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Fig 4b — "disease reorganizes liver space", ONE merged panel (replaces the
# separate spatial_svg_dynamics + cosmx_il32_colocalization panels).
#
#  [a] the TRANSCRIPTOME reorganizes (Visium GSE192741, donor-blocked): of the
#      disease spatially-variable genes, 94 GAIN spatial structure in disease,
#      271 MAINTAIN it, 99 LOSE it; the reorganized signal tilts to the
#      pericentral detoxification zone (pericentral : periportal among zoned SVGs).
#  [b] the CELLS reposition (Govaere CosMx single-cell spatial): in all 3 MASH
#      slides IL32+ hepatocytes sit CLOSER to macrophages and their nearest
#      macrophages carry MORE CD74 — shown as a compact direction tile (value in
#      cell; magenta = concordant), Normal slide faded. NO lollipop
#      (see memory/feedback-no-lollipop).
#
# Data: Analysis/Spatial/results/svg/differential_svgs.csv,
#       Analysis/Spatial/results/zonation/deg_zonation_classification.csv,
#       Analysis/Spatial/results/govaere2026/il32_colocalization/il32_coloc_per_slide.csv
# Output: figures/main/fig4_validation/fig4b_spatial_reorganization.pdf
# Env:    rnaseq
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
SP <- file.path(BASE, "Analysis/Spatial/results")

# ── [a] SVG dynamics + pericentral tilt ──────────────────────────────────────
svg <- fread(file.path(SP, "svg/differential_svgs.csv"))
setnames(svg, 1, "gene")
tl <- function(x) tolower(as.character(x)) %in% c("true", "1", "t")
svg[, `:=`(h = tl(svg_healthy), m = tl(svg_masld))]
n_gained <- svg[m & !h, .N]; n_maint <- svg[h & m, .N]; n_lost <- svg[h & !m, .N]

zon <- fread(file.path(SP, "zonation/deg_zonation_classification.csv"))
# Zonation tilt among the zoned DEGs (the canonical 8:4 count; full classification
# file — Pan-lobular excluded). Matches the manuscript.
zcol <- intersect(c("zonation_class"), names(zon))[1]
n_peri <- zon[get(zcol) == "Pericentral-enriched", .N]
n_port <- zon[get(zcol) == "Periportal-enriched", .N]

svg_df <- data.table(
  cls = factor(c("Gained", "Maintained", "Lost"), levels = c("Gained", "Maintained", "Lost")),
  n   = c(n_gained, n_maint, n_lost),
  fill = c(masld_colors$up, masld_colors$ns, masld_colors$down))

p_a <- ggplot(svg_df, aes(cls, n, fill = fill)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = formatC(n, big.mark = ",", format = "d")),
            vjust = -0.35, size = PUB_GEOM_TEXT + 0.2, fontface = "bold") +
  scale_fill_identity() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = NULL, y = "spatially variable genes") +
  theme_masld() + theme_pub() +
  theme(plot.title = element_text(size = PUB_TITLE, face = "bold", color = "black"),
        plot.subtitle = element_text(size = PUB_AXIS_TEXT - 0.5, color = "grey35"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(face = "bold", color = "black"),
        plot.margin = margin(2, 4, 0, 3))

# ── [b] CosMx IL32→macrophage repositioning (direction tile) ─────────────────
cm <- fread(file.path(SP, "govaere2026/il32_colocalization/il32_coloc_per_slide.csv"))
cm[, is_mash := tolower(as.character(is_mash)) == "true"]
cm <- cm[order(-is_mash, slide)]
cm[, slide_lab := fifelse(is_mash, slide, paste0(slide, " (Normal)"))]
cm[, slide_lab := factor(slide_lab, levels = rev(slide_lab))]

tile <- rbind(
  cm[, .(slide_lab, is_mash, metric = "IL32-hi hepatocyte\ncloser to macrophage",
         value = delta_high_minus_low_um, concord = delta_high_minus_low_um < 0,
         txt = sprintf("%.2f µm", delta_high_minus_low_um))],
  cm[, .(slide_lab, is_mash, metric = "neighbour macrophage\nCD74 (ρ)",
         value = rho_il32_vs_knn_cd74, concord = rho_il32_vs_knn_cd74 > 0,
         txt = sprintf("%+.2f", rho_il32_vs_knn_cd74))])
tile[, metric := factor(metric, levels = c("IL32-hi hepatocyte\ncloser to macrophage",
                                           "neighbour macrophage\nCD74 (ρ)"))]
# fill: concordant direction in MASH = magenta; not = grey; Normal = faded
tile[, fillcol := fifelse(!is_mash, "#E8E8E8",
                  fifelse(concord, masld_colors$up, masld_colors$down))]
n_mash <- cm[is_mash == TRUE, .N]
n_dc   <- cm[is_mash & delta_high_minus_low_um < 0, .N]
n_cc   <- cm[is_mash & rho_il32_vs_knn_cd74 > 0, .N]

p_b <- ggplot(tile, aes(metric, slide_lab, fill = fillcol)) +
  geom_tile(color = "white", linewidth = 1.0) +
  geom_text(aes(label = txt, color = ifelse(is_mash, "white", "grey55")),
            size = PUB_GEOM_TEXT, fontface = "bold") +
  scale_fill_identity() + scale_color_identity() +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_masld() + theme_pub() +
  theme(plot.title = element_text(size = PUB_TITLE, face = "bold", color = "black"),
        plot.subtitle = element_text(size = PUB_AXIS_TEXT - 0.5, color = "grey35"),
        axis.text.x.top = element_text(size = PUB_AXIS_TEXT, lineheight = 0.85, color = "black"),
        axis.text.y = element_text(size = PUB_AXIS_TEXT, face = "bold", color = "black"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        plot.margin = margin(1, 4, 2, 3))

# ── compose ──────────────────────────────────────────────────────────────────
fig <- p_a / p_b + plot_layout(heights = c(1, 1.05))
out <- file.path(FIG4_DIR, "fig4b_spatial_reorganization.pdf")
save_fig(fig, out, width = fig_full_width * 0.58, height = 2.6)
message("Saved: ", out)
message(sprintf("[spatial] SVG gained=%d maintained=%d lost=%d ; pericentral:periportal=%d:%d ; CosMx %d/%d closer, %d/%d CD74+",
                n_gained, n_maint, n_lost, n_peri, n_port, n_dc, n_mash, n_cc, n_mash))
