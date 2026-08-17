#!/usr/bin/env Rscript
# KEY MESSAGE: Disease reorganizes the spatial transcriptome. Of the disease-set
# spatially variable genes (SVGs) in GSE192741 (Guilliams et al.; 5 slides,
# 6,546 spots, 2 healthy + 3 steatotic donors, donor-blocked):
#     94 GAIN spatial structure in steatotic liver (disease-emergent),
#    271 MAINTAIN it, and  99 LOSE it.
# The reorganized disease signal concentrates in the pericentral detoxification
# zone: 8 pericentral-enriched DEGs outnumber 4 periportal-enriched DEGs (8:4).
#
# CAVEAT (ledger rows 34/35/49): GSE192741 is donor-blocked at n=5 (2 healthy /
# 3 steatotic). This panel reports spatial ORGANIZATION (whether a gene is
# spatially structured), NOT disease DIRECTION (up/down). Spot-level permutation
# is anticonservative at this donor count. SERPINE1-type organization gains are
# reported descriptively (cf. Kim 2024 J Hepatol).
#
# Per PI directive, statistics live in the figure legend (emitted to stdout via
# message()); the panel carries only the count/gene annotations the bars need.
#
# Data sources (read from disk; never hardcoded from prose):
#   SVG gain/maintain/loss  -> Analysis/Spatial/results/svg/differential_svgs.csv
#       gained     = svg_masld & !svg_healthy  (== category disease_emergent_SVG) = 94
#       maintained = svg_masld &  svg_healthy                                      = 271
#       lost       = !svg_masld &  svg_healthy (== category disease_lost_SVG)      = 99
#   Zonation tilt           -> Analysis/Spatial/results/zonation/deg_zonation_classification.csv
#       zonation_class == "Pericentral-enriched" = 8 genes
#       zonation_class == "Periportal-enriched"  = 4 genes
#
# Output: figures/main/fig5_molecular_context/spatial_svg_dynamics.pdf
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

# Colors. SVG dynamics is an ORGANIZATION axis, not a disease-direction axis, but
# we borrow the semantic palette as a visual cue: Gained = disease magenta,
# Maintained = neutral gray, Lost = blue. Zonation is a spatial-ZONE axis:
# pericentral warm orange (matches the CYP3A4 zonation panel) vs periportal blue.
GAINED_COL      <- masld_colors$up        # "#C9265E" disease-emergent / gained
MAINTAINED_COL  <- masld_colors$ns        # "#9E9E9E" neutral
LOST_COL        <- masld_colors$down      # "#1565C0" lost
PERICENTRAL_COL <- "#E65100"              # warm orange (CYP3A4 zonation panel)
PERIPORTAL_COL  <- masld_colors$down      # "#1565C0"

# ── Data 1: SVG dynamics ─────────────────────────────────────────────────────
svg <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/svg/differential_svgs.csv"),
  row.names = 1, stringsAsFactors = FALSE)

sh <- svg$svg_healthy %in% c("True", TRUE)
sm <- svg$svg_masld   %in% c("True", TRUE)

n_gained     <- sum( sm & !sh)   # disease-emergent: SVG in steatotic only
n_maintained <- sum( sm &  sh)   # SVG in both
n_lost       <- sum(!sm &  sh)   # SVG in healthy only

# Sanity-check the boolean crosstab against the curated category column.
stopifnot(n_gained == sum(svg$category == "disease_emergent_SVG"))
stopifnot(n_lost   == sum(svg$category == "disease_lost_SVG"))

svg_df <- data.frame(
  state = factor(c("Gained", "Maintained", "Lost"),
                 levels = c("Gained", "Maintained", "Lost")),
  n     = c(n_gained, n_maintained, n_lost)
)
svg_fills <- c(Gained = GAINED_COL, Maintained = MAINTAINED_COL, Lost = LOST_COL)

p_svg <- ggplot(svg_df, aes(x = state, y = n, fill = state)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = n), vjust = -0.45, size = PUB_GEOM_TEXT + 0.6,
            fontface = "plain", color = "black") +
  scale_fill_manual(values = svg_fills, guide = "none") +
  scale_y_continuous(limits = c(0, max(svg_df$n) * 1.18),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Disease-set SVGs (n)",
       title = "Spatial structure gained / maintained / lost") +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(size = PUB_AXIS_TEXT + 0.5))

# ── Data 2: DEG zonation tilt ────────────────────────────────────────────────
zon <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/zonation/deg_zonation_classification.csv"),
  stringsAsFactors = FALSE)

peri_central <- sort(zon$gene[zon$zonation_class == "Pericentral-enriched"])
peri_portal  <- sort(zon$gene[zon$zonation_class == "Periportal-enriched"])

# Direct gene labels ride just past each bar, wrapped so the 8-gene pericentral
# list never collides with the bar or the panel edge.
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

# Generous headroom on the count axis so the wrapped gene captions sit clear of
# the bars and inside the panel.
ymax <- max(zon_df$n) + 14
p_zon <- ggplot(zon_df, aes(x = zone, y = n, fill = zone)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = n), hjust = -0.4, size = PUB_GEOM_TEXT + 0.6,
            fontface = "plain", color = "black") +
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
       title = "Disease signal tilts pericentral (8:4)") +
  theme_masld() + theme_pub() +
  theme(axis.text.y = element_text(size = PUB_AXIS_TEXT + 0.5))

# ── Assemble (compact, ~column width, ~2.2in tall; NO composite title) ───────
p_out <- (p_svg | p_zon) + plot_layout(widths = c(1, 1.15))

# Legend material to stdout (PI directive — keep off the panel). Includes the
# donor-block + organization-not-direction caveat (ledger rows 34/35/49).
message(sprintf(
  "[spatial_svg_dynamics legend] Disease reorganizes the spatial transcriptome (GSE192741; 5 slides / 6,546 spots / 2 healthy + 3 steatotic donors, donor-blocked). Of the disease-set SVGs: %d GAINED spatial structure (disease-emergent), %d MAINTAINED, %d LOST (SVG totals: %d steatotic, %d healthy). The reorganized signal concentrates pericentrally: %d pericentral-enriched DEGs (%s) vs %d periportal-enriched DEGs (%s) = 8:4 tilt toward the detoxification zone. CAVEAT: this reports spatial ORGANIZATION (squidpy Moran's I, per-condition), not disease DIRECTION; n=5 donors and spot-level permutation is anticonservative (cf. Kim 2024 J Hepatol for SERPINE1-type organization gains).",
  n_gained, n_maintained, n_lost, sum(sm), sum(sh),
  length(peri_central), paste(peri_central, collapse = ", "),
  length(peri_portal),  paste(peri_portal,  collapse = ", ")))

out <- file.path(FIG4_DIR, "_supp", "spatial_svg_dynamics.pdf")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(out, width = fig_col_width, height = 2.2)
print(p_out)
dev.off()
message("Saved: ", out)
