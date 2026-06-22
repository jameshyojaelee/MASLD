#!/usr/bin/env Rscript
# KEY MESSAGE: Spatial-transcriptomics summary for Figure 4b (GSE192741,
# Guilliams et al.; 5 slides, 6,546 spots, 2 healthy + 3 steatotic donors,
# donor-blocked). Two compact panels, CSV-only (no h5ad — light render):
#
#   (a) SVG dynamics across condition: of the disease-set SVGs,
#       94 GAINED spatial structure in steatotic liver (disease-emergent),
#       271 MAINTAINED it, 99 LOST it.
#   (b) DEG zonation tilt: pericentral-enriched DEGs (8) outnumber
#       periportal-enriched DEGs (4), concentrating the spatial disease
#       signal in the pericentral detoxification zone.
#
# Stats live in the figure legend (PI directive), emitted to stdout — not on
# the panel beyond the count annotations the bars themselves carry.
#
# Data sources (numbers verified from disk, not hardcoded blindly):
#   SVG gain/maintain/loss:
#     Analysis/Spatial/results/svg/differential_svgs.csv
#       gained     = svg_masld & !svg_healthy  (== category disease_emergent_SVG) = 94
#       maintained = svg_masld &  svg_healthy                                      = 271
#       lost       = !svg_masld &  svg_healthy (== category disease_lost_SVG)      = 99
#   Zonation counts/genes:
#     Analysis/Spatial/results/zonation/deg_zonation_classification.csv
#       zonation_class == "Pericentral-enriched" = 8 genes
#       zonation_class == "Periportal-enriched"  = 4 genes
#
# Output: figures/main/fig4_validation/spatial_svg_summary.pdf
# Env:    rnaseq

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# Zone palette (matches fig4f_cyp3a4_zonation.R): pericentral warm orange,
# periportal blue.
PERICENTRAL_COL <- "#E65100"
PERIPORTAL_COL  <- "#1565C0"

# ── Data 1: SVG dynamics ─────────────────────────────────────────────────────
svg <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/svg/differential_svgs.csv"),
  row.names = 1, stringsAsFactors = FALSE)

sh <- svg$svg_healthy %in% c("True", TRUE)
sm <- svg$svg_masld   %in% c("True", TRUE)

n_gained     <- sum( sm & !sh)   # disease-emergent: SVG in steatotic only
n_maintained <- sum( sm &  sh)   # SVG in both
n_lost       <- sum(!sm &  sh)   # SVG in healthy only

# Sanity check against the curated category column (covers gained + lost).
stopifnot(n_gained == sum(svg$category == "disease_emergent_SVG"))
stopifnot(n_lost   == sum(svg$category == "disease_lost_SVG"))

svg_df <- data.frame(
  state = factor(c("Gained", "Maintained", "Lost"),
                 levels = c("Gained", "Maintained", "Lost")),
  n     = c(n_gained, n_maintained, n_lost)
)
svg_fills <- c(Gained = PERICENTRAL_COL, Maintained = "#9E9E9E", Lost = PERIPORTAL_COL)

p_svg <- ggplot(svg_df, aes(x = state, y = n, fill = state)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = n), vjust = -0.35, size = PUB_GEOM_TEXT + 0.4,
            fontface = "bold", color = "black") +
  scale_fill_manual(values = svg_fills, guide = "none") +
  scale_y_continuous(limits = c(0, max(svg_df$n) * 1.18),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Disease-set SVGs (n)",
       title = "Spatial SVG dynamics") +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(size = PUB_AXIS_TEXT + 0.5))

# ── Data 2: DEG zonation tilt ────────────────────────────────────────────────
zon <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/zonation/deg_zonation_classification.csv"),
  stringsAsFactors = FALSE)

peri_central <- sort(zon$gene[zon$zonation_class == "Pericentral-enriched"])
peri_portal  <- sort(zon$gene[zon$zonation_class == "Periportal-enriched"])

# Horizontal paired bars: pericentral (8, warm) over periportal (4, blue).
# Gene names ride just past the bar end, wrapped over two short lines so the
# long pericentral list never collides with the bar or the panel edge.
wrap_genes <- function(g, per_line = 4) {
  idx <- ceiling(seq_along(g) / per_line)
  paste(vapply(split(g, idx),
               function(x) paste(x, collapse = ", "),
               character(1)),
        collapse = "\n")
}
peri_central_lab <- wrap_genes(peri_central, per_line = 4)   # 2 lines (8 genes)
peri_portal_lab  <- wrap_genes(peri_portal,  per_line = 4)   # 1 line  (4 genes)

zon_df <- data.frame(
  zone = factor(c("Pericentral", "Periportal"),
                levels = c("Periportal", "Pericentral")),
  n    = c(length(peri_central), length(peri_portal)),
  lab  = c(peri_central_lab, peri_portal_lab)
)

# Generous headroom on the count axis so the wrapped gene captions sit clear
# of the bars and inside the panel.
ymax <- max(zon_df$n) + 14
p_zon <- ggplot(zon_df, aes(x = zone, y = n, fill = zone)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = n), hjust = -0.35, size = PUB_GEOM_TEXT + 0.4,
            fontface = "bold", color = "black") +
  geom_text(aes(label = lab, color = zone), y = max(zon_df$n) + 1.4,
            hjust = 0, vjust = 0.5, size = PUB_GEOM_TEXT - 0.1,
            fontface = "italic", lineheight = 0.9) +
  scale_fill_manual(values = c(Pericentral = PERICENTRAL_COL,
                               Periportal  = PERIPORTAL_COL), guide = "none") +
  scale_color_manual(values = c(Pericentral = PERICENTRAL_COL,
                                Periportal  = PERIPORTAL_COL), guide = "none") +
  scale_y_continuous(limits = c(0, ymax),
                     breaks = c(0, 4, 8),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_flip() +
  labs(x = NULL, y = "Zonation-classified DEGs (n)",
       title = "DEG zonation tilt") +
  theme_masld() + theme_pub() +
  theme(axis.text.y = element_text(size = PUB_AXIS_TEXT + 0.5))

# ── Assemble (compact, half-to-two-thirds width, ~2.2in tall) ────────────────
p_out <- (p_svg | p_zon) + plot_layout(widths = c(1, 1.15))

# Stats to stdout (legend material; PI directive — keep off the panel).
message(sprintf(
  "[spatial_svg_summary legend] GSE192741 disease-set SVGs: %d gained / %d maintained / %d lost (SVG totals: %d steatotic, %d healthy). Zonation tilt: %d pericentral-enriched (%s) vs %d periportal-enriched (%s) DEGs.",
  n_gained, n_maintained, n_lost, sum(sm), sum(sh),
  length(peri_central), peri_central_lab,
  length(peri_portal),  peri_portal_lab))

out <- file.path(FIG4_DIR, "spatial_svg_summary.pdf")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(out, width = fig_col_width, height = 2.2)
print(p_out)
dev.off()
message("Saved: ", out)
