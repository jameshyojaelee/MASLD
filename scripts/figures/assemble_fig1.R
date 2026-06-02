##############################################################################
# assemble_fig1.R
# Composite assembly for Figure 1 (Atlas Overview).
# Loads canonical individual panel PDFs via magick, lays them out on a single
# cowplot canvas with bold A/B/C labels (top-left, Helvetica/Arial ~14pt).
#
# Inputs : figures/main/fig1_atlas_overview/panels/*.pdf
# Output : figures/main/fig1_atlas_overview/fig1_composite.pdf  (flat at fig root)
#
# Liang aesthetic: 7.2" wide (Cell/Nature double column), dense 3x3 grid,
# white background, ~0.1" inter-panel padding, no overall title.
#
# NEW (2026-05-18): includes fig1_volcano_atac_overlay.pdf as panel (i).
##############################################################################

suppressPackageStartupMessages({
  library(cowplot)
  library(ggplot2)
  library(grid)
  library(png)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG1_DIR, "panels")
OUT_PDF   <- file.path(FIG1_DIR, "fig1_composite.pdf")

# ---------------------------------------------------------------------------
# Canonical panel selection (see docs/manuscript/05_figure_legends.md Fig 1).
# Panel (a) study design is an Illustrator schematic with no PDF in panels/;
# rendered as a labelled white tile (placeholder) so the grid keeps labels A-I
# in sync.  Legacy / alternative files explicitly excluded:
#   - fig1_volcano_vector.pdf  (1.1 MB rasterised duplicate of fig1_volcano.pdf)
#   - fig1h_optionb.pdf        (alternative framing; discovery scatter chosen)
# ---------------------------------------------------------------------------
panels <- list(
  list(label = "a", file = NA_character_,
       title = "Study design\n(schematic)"),
  list(label = "b", file = "fig1b.pdf",                       title = NULL),
  list(label = "c", file = "fig1c.pdf",                       title = NULL),
  list(label = "d", file = "fig1d.pdf",                       title = NULL),
  list(label = "e", file = "fig1_volcano.pdf",                title = NULL),
  list(label = "f", file = "fig1_volcano_atac_overlay.pdf",   title = NULL),
  list(label = "g", file = "fig1f.pdf",                       title = NULL),
  list(label = "h", file = "fig1g.pdf",                       title = NULL),
  list(label = "i", file = "fig1h_integration_discovery.pdf", title = NULL)
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
RENDER_DENSITY <- 400  # dpi for PDF -> raster import (high for print clarity)

# Rasterise a single-page PDF to PNG using pdftoppm (preferred; no ghostscript
# required) and read back as a raster array.  Returns NULL if the conversion
# fails, so the caller can fall back to a placeholder.
read_panel <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  tmp_prefix <- tempfile()
  cmd_args <- c("-r", RENDER_DENSITY, "-png", "-f", "1", "-l", "1",
                shQuote(path), shQuote(tmp_prefix))
  status <- suppressWarnings(
    system2("pdftoppm", args = cmd_args, stdout = TRUE, stderr = TRUE)
  )
  png_path <- paste0(tmp_prefix, "-1.png")
  if (!file.exists(png_path)) png_path <- paste0(tmp_prefix, "-01.png")
  if (!file.exists(png_path)) {
    # Fallback to ImageMagick convert
    png_path <- paste0(tmp_prefix, ".png")
    system2("convert", args = c("-density", RENDER_DENSITY, shQuote(path),
                                "-flatten", shQuote(png_path)),
            stdout = FALSE, stderr = FALSE)
  }
  if (!file.exists(png_path)) return(NULL)
  img <- png::readPNG(png_path)
  unlink(png_path)
  img
}

placeholder_grob <- function(text) {
  # ggdraw with a thin black border + centered text
  ggdraw() +
    draw_label(text, fontface = "plain", size = 9,
               colour = "grey30", lineheight = 1.05, hjust = 0.5, vjust = 0.5) +
    theme(panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.3))
}

panel_to_grob <- function(p) {
  if (is.na(p$file)) {
    return(placeholder_grob(p$title %||% sprintf("Panel %s", toupper(p$label))))
  }
  fpath <- file.path(PANEL_DIR, p$file)
  img <- read_panel(fpath)
  if (is.null(img)) {
    return(placeholder_grob(sprintf("Missing: %s", p$file)))
  }
  # img is a numeric raster array (H x W x channels); wrap as a raster grob.
  raster_grob <- grid::rasterGrob(img, interpolate = TRUE,
                                  width = unit(1, "npc"),
                                  height = unit(1, "npc"))
  ggdraw() + cowplot::draw_grob(raster_grob)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# ---------------------------------------------------------------------------
# Build panel grobs + check missing files
# ---------------------------------------------------------------------------
missing_files <- character()
grob_list <- list()
for (i in seq_along(panels)) {
  p <- panels[[i]]
  if (!is.na(p$file) && !file.exists(file.path(PANEL_DIR, p$file))) {
    missing_files <- c(missing_files, p$file)
  }
  grob_list[[i]] <- panel_to_grob(p)
}

if (length(missing_files) > 0) {
  message("WARNING: missing panel PDFs (rendering placeholders):\n  - ",
          paste(missing_files, collapse = "\n  - "))
}

# ---------------------------------------------------------------------------
# Layout: 3x3 grid, 7.2" wide (Cell/Nature double-column standard).
# Inter-panel padding ~0.05" via cowplot rel_widths/heights default + theme.
# Per-row aspect: keep panels square-ish; total height ~ 7.2" so each panel
# occupies ~2.4" x 2.4".  Adjust height for label headroom.
# ---------------------------------------------------------------------------
W_IN <- 7.2
H_IN <- 7.2   # 3 rows of ~2.4" tiles

labels <- vapply(panels, function(p) toupper(p$label), character(1))

composite <- cowplot::plot_grid(
  plotlist     = grob_list,
  ncol         = 3,
  nrow         = 3,
  labels       = labels,
  label_size   = 14,
  label_fontface = "bold",
  label_fontfamily = "sans",  # Helvetica/Arial-equivalent on Linux
  label_x      = 0.005,
  label_y      = 0.995,
  hjust        = 0,
  vjust        = 1,
  rel_widths   = c(1, 1, 1),
  rel_heights  = c(1, 1, 1)
)

# Wrap with a white background canvas
composite_bg <- cowplot::ggdraw() +
  draw_grob(grid::rectGrob(gp = grid::gpar(fill = "white", col = NA))) +
  draw_plot(composite, x = 0, y = 0, width = 1, height = 1)

# ---------------------------------------------------------------------------
# Save (cairo_pdf for proper font embedding)
# ---------------------------------------------------------------------------
ggsave(filename = OUT_PDF,
       plot     = composite_bg,
       width    = W_IN,
       height   = H_IN,
       units    = "in",
       device   = cairo_pdf,
       bg       = "white")

cat("Saved composite:", OUT_PDF, "\n")
cat(sprintf("Dimensions: %.1f x %.1f inches (%d panels, 3x3 grid)\n",
            W_IN, H_IN, length(panels)))
if (length(missing_files) > 0) {
  cat("Panels rendered as placeholders:\n  - ",
      paste(missing_files, collapse = "\n  - "), "\n", sep = "")
}
