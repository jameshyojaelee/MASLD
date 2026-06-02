#!/usr/bin/env Rscript
# ============================================================================
# fig2_progression.R
# Single compositor for Fig 2 (Progression). Imports 7 panel PDFs and tiles
# them at 180 x 215 mm via patchwork + magick. Each panel script is run
# independently; this compositor ONLY tiles the rendered PDFs.
#
# Panels (a-f):
#   a — Four-transition DEG cascade + LOCO inset (fig2_panel_cascade_degs.pdf)
#   b — F1-F3 inflection composite (fig2_panel_inflection_composite.pdf)
#   c — Cell-type proportion trajectory, full width (fig2_panel_celltype_proportion.pdf)
#   d — Hotspot autocorrelation cascade, full width (fig2_hotspot_cascade.pdf)
#   e — scRNA-seq UMAP with MASH vs Control density (fig2e.pdf)
#   f — Stage-progressive CCC chord (fig2_ccc_chord.pdf)
#
# Output: figures/main/fig2_progression_sex/fig2_progression.pdf
# ============================================================================

suppressPackageStartupMessages({
  library(magick)
  library(grid)
  library(ggplot2)
  library(patchwork)
  library(ggplotify)
})

# magick::image_read_pdf needs pdftools; magick::image_read alone uses the
# underlying ImageMagick PDF delegate (ghostscript). Use the latter to avoid
# the pdftools dependency in the rnaseq env.

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
OUT_PDF   <- file.path(FIG2_DIR, "fig2_progression.pdf")

panels <- list(
  a = "fig2_panel_cascade_degs.pdf",
  b = "fig2_panel_inflection_composite.pdf",
  c = "fig2_panel_celltype_proportion.pdf",
  d = "fig2_hotspot_cascade.pdf",
  e = "fig2e.pdf",
  f = "fig2_ccc_chord.pdf"
)

# Verify all exist
for (label in names(panels)) {
  p <- file.path(PANEL_DIR, panels[[label]])
  if (!file.exists(p)) stop("Missing panel: ", p)
}
cat("[ok] all 6 panels found on disk\n")

# Load each panel as a ggplot via magick::image_ggplot (handles bitmap class)
read_panel <- function(path, label) {
  img <- magick::image_read(path, density = 300)
  g <- magick::image_ggplot(img, interpolate = TRUE)
  g + labs(tag = label) +
    theme(plot.tag = element_text(size = 10, face = "bold",
                                   color = "black"),
          plot.tag.position = c(0.02, 0.97),
          plot.margin = margin(2, 2, 2, 2))
}

g_a <- read_panel(file.path(PANEL_DIR, panels$a), "a")
g_b <- read_panel(file.path(PANEL_DIR, panels$b), "b")
g_c <- read_panel(file.path(PANEL_DIR, panels$c), "c")
g_d <- read_panel(file.path(PANEL_DIR, panels$d), "d")
g_e <- read_panel(file.path(PANEL_DIR, panels$e), "e")
g_f <- read_panel(file.path(PANEL_DIR, panels$f), "f")

# Layout: 4 rows
#   Row 1 (height 65 mm): a | b
#   Row 2 (height 55 mm): c (full — hep loss + non-hep stack)
#   Row 3 (height 45 mm): d (full — hotspot cascade)
#   Row 4 (height 55 mm): e | f
fig <- (g_a | g_b) /
       g_c /
       g_d /
       (g_e | g_f) +
       plot_layout(heights = c(65, 55, 45, 55))

cat(sprintf("[save] %s\n", OUT_PDF))
ggsave(OUT_PDF, fig,
       width  = 180 / 25.4,
       height = 220 / 25.4,
       units  = "in",
       device = cairo_pdf)
cat("[done]\n")
