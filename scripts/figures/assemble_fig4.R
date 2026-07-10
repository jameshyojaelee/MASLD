# =====================================================================
# assemble_fig4.R
# ---------------------------------------------------------------------
# Compose Fig 4 (Proteomics + Spatial + Plasma + Pharma validation) from
# the 6 canonical sub-panel PDFs in figures/main/fig4_validation/panels/.
#
# Canonical letter order (per docs/manuscript/05_figure_legends.md
# Figure 4 + figures/README.md Fig 4 section):
#   A  proteomics DE volcano                figS4b.pdf (was fig4a.pdf)
#   (fig4b.pdf / fig4c.pdf panels RETIRED 2026-07-07 — removed from this composite)
#   D  plasma fibrosis classifier (F4)      fig4d.pdf
#   E  drug regulon disruption (Venn+forest) fig4_panel_E_drug_regulon_disruption.pdf  [NEW]
#   F  drug-target pharma panels            fig4_pharma_panels.pdf
#
# Layout: 3 rows x 2 cols (Liang aesthetic).
# Output: figures/main/fig4_validation/fig4_composite.pdf
# =====================================================================

suppressPackageStartupMessages({
  library(cowplot)
  library(ggplot2)
  library(magick)
})

source("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/load_figure_data.R")

PANEL_DIR <- file.path(FIG4_DIR, "panels")
OUT_PDF   <- file.path(FIG4_DIR, "fig4_composite.pdf")

panel_files <- c(
  A = "figS4b.pdf",
  D = "fig4d.pdf",
  E = "fig4_panel_E_drug_regulon_disruption.pdf",
  F = "fig4_pharma_panels.pdf"
)

# Verify all panels exist before reading.
for (lbl in names(panel_files)) {
  fp <- file.path(PANEL_DIR, panel_files[[lbl]])
  if (!file.exists(fp)) stop("Missing panel ", lbl, ": ", fp)
}

read_panel <- function(fp, density = 300) {
  # Use ghostscript directly (rnaseq env lacks `pdftools` which magick's
  # image_read_pdf() requires). Render page 1 to a hi-res PNG, then read.
  tmp_png <- tempfile(fileext = ".png")
  cmd <- sprintf(
    "gs -dQUIET -dBATCH -dNOPAUSE -dFirstPage=1 -dLastPage=1 -sDEVICE=pngalpha -r%d -sOutputFile=%s %s",
    density, shQuote(tmp_png), shQuote(fp)
  )
  rc <- system(cmd)
  if (rc != 0 || !file.exists(tmp_png)) stop("gs render failed: ", fp)
  magick::image_read(tmp_png)
}

panel_to_gg <- function(img, label) {
  # Render the magick image inside a ggdraw canvas and stamp the bold
  # A/B/C... label in the top-left of the panel area.
  info <- magick::image_info(img)
  aspect <- info$height / info$width
  p <- cowplot::ggdraw() +
       cowplot::draw_image(img, x = 0, y = 0, width = 1, height = 1) +
       cowplot::draw_label(
         label,
         x = 0.005, y = 0.995, hjust = 0, vjust = 1,
         fontface = "plain", size = 6
       )
  attr(p, "aspect") <- aspect
  p
}

cat("Reading panels...\n")
imgs <- lapply(panel_files, function(f) read_panel(file.path(PANEL_DIR, f)))
names(imgs) <- names(panel_files)

panels_gg <- mapply(panel_to_gg, imgs, names(imgs), SIMPLIFY = FALSE)

# ---- layout sizing -----------------------------------------------------
# 3 rows x 2 cols. Composite 7.2" wide (Cell/Nature single-figure width).
# Each cell is therefore ~3.6" wide. Pick row heights from per-panel
# aspect ratios so nothing is squashed.
COMPOSITE_W <- 7.2
N_COL       <- 2
N_ROW       <- 2
CELL_W      <- COMPOSITE_W / N_COL

panel_order <- c("A","D","E","F")
aspects <- sapply(panel_order, function(k) attr(panels_gg[[k]], "aspect"))
# Per-row height = max of the two panels in that row * cell width
row_heights <- sapply(seq_len(N_ROW), function(r) {
  in_row <- panel_order[((r - 1) * N_COL + 1):(r * N_COL)]
  max(aspects[in_row]) * CELL_W
})
COMPOSITE_H <- sum(row_heights) + 0.15  # tiny pad

cat(sprintf("Composite: %.2f x %.2f inches (rows: %s)\n",
            COMPOSITE_W, COMPOSITE_H,
            paste(sprintf("%.2f", row_heights), collapse = " / ")))

composite <- cowplot::plot_grid(
  plotlist  = panels_gg[panel_order],
  ncol      = N_COL,
  nrow      = N_ROW,
  rel_heights = row_heights,
  align     = "none"
)

ggplot2::ggsave(
  filename = OUT_PDF,
  plot     = composite,
  width    = COMPOSITE_W,
  height   = COMPOSITE_H,
  units    = "in",
  device   = cairo_pdf
)

cat("Wrote: ", OUT_PDF, "\n", sep = "")
