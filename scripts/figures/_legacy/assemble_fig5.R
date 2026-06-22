# Fig 5 composite assembler — multi-evidence convergence + TF drug-target prioritization.
#
# Combines the canonical Fig 5 panel set from individual panel PDFs into a single
# multi-panel composite (Liang aesthetic). The four canonical panels (May 2026):
#   (a) Sources-active histogram + per-source recovery               (fig5a.pdf)
#   (b) Multi-modal convergence matrix, top 50 × 7 sources           (fig5b.pdf)
#   (c) Refined TF convergence scatter (4-way validated TFs)         (fig5_tf_convergence_scatter.pdf)
#   (d) 4-way TF survival lollipop                                   (fig5_tf_4way_survival_lollipop.pdf)
#
# Layout: 2 rows × 2 columns at 7.2 inch composite width.
# Output: figures/main/fig5_convergence/fig5_composite.pdf (FLAT in FIG5_DIR root).

suppressPackageStartupMessages({
  library(magick)
  library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# FIG5_DIR is the canonical path constant from load_figure_data.R
PANELS_DIR <- file.path(FIG5_DIR, "panels")
OUT_PDF    <- file.path(FIG5_DIR, "fig5_composite.pdf")

panel_files <- c(
  a = file.path(PANELS_DIR, "fig5a.pdf"),
  b = file.path(PANELS_DIR, "fig5b.pdf"),
  c = file.path(PANELS_DIR, "fig5_tf_convergence_scatter.pdf"),
  d = file.path(PANELS_DIR, "fig5_tf_4way_survival_lollipop.pdf")
)
missing <- panel_files[!file.exists(panel_files)]
if (length(missing) > 0)
  stop("Missing panel PDFs:\n  ", paste(missing, collapse = "\n  "))

# Render each panel PDF to a high-DPI raster via ghostscript (consistent with
# assemble_fig3.R; magick::image_read_pdf relies on a working PDF backend in
# the rnaseq env that is not always installed).
TMP <- tempdir()
imgs <- lapply(names(panel_files), function(lab) {
  src <- panel_files[[lab]]
  png_out <- file.path(TMP, paste0("fig5_", lab, ".png"))
  cmd <- sprintf(
    "gs -dQUIET -dBATCH -dNOPAUSE -dFirstPage=1 -dLastPage=1 -sDEVICE=pngalpha -r300 -sOutputFile=%s %s >/dev/null 2>&1",
    png_out, src
  )
  system(cmd)
  if (!file.exists(png_out))
    stop("ghostscript rasterization failed for: ", src)
  image_read(png_out)
})
names(imgs) <- names(panel_files)

# Normalize each panel to a common cell width; rows share width via image_resize.
CELL_W <- 1800  # pixels per panel column at 300 DPI -> ~6 inch per cell; comfortably > 3.6"

imgs_lab <- lapply(imgs, function(im) image_resize(im, paste0(CELL_W, "x")))

# 2x2 layout (labels added via grid after raster placement to avoid magick font issues)
row1 <- image_append(c(imgs_lab$a, imgs_lab$b))
row2 <- image_append(c(imgs_lab$c, imgs_lab$d))

# Match row widths
RW <- max(image_info(row1)$width, image_info(row2)$width)
row1 <- image_resize(row1, paste0(RW, "x"))
row2 <- image_resize(row2, paste0(RW, "x"))

composite <- image_append(c(row1, row2), stack = TRUE)

# Final downscale to 7.2 inches wide × proportional height.
# We render the final PDF page at 7.2" × (height/width * 7.2)".
ci <- image_info(composite)
final_width_in  <- 7.2
final_height_in <- final_width_in * ci$height / ci$width

# Use cairo_pdf via grid for a vector-wrapped page (raster panels embedded).
cairo_pdf(OUT_PDF, width = final_width_in, height = final_height_in,
          family = "Helvetica")
grid.newpage()
# Convert magick image to raster grob spanning the full page.
ras <- as.raster(composite)
grid.raster(ras, width = unit(1, "npc"), height = unit(1, "npc"))

# Bold panel labels (A-D) in the top-left of each cell. Cells are arranged
# in a 2x2 grid; npc coordinates: column edges 0 / 0.5 / 1, row edges 1 / 0.5 / 0.
label_positions <- list(
  A = c(0.005, 0.985),
  B = c(0.505, 0.985),
  C = c(0.005, 0.485),
  D = c(0.505, 0.485)
)
for (lab in names(label_positions)) {
  xy <- label_positions[[lab]]
  grid.text(lab, x = unit(xy[1], "npc"), y = unit(xy[2], "npc"),
            just = c("left", "top"),
            gp = gpar(fontsize = 11, fontface = "bold", fontfamily = "Helvetica"))
}
dev.off()

cat(sprintf(
  "Wrote: %s (%.2f x %.2f in; composite raster %d x %d px)\n",
  OUT_PDF, final_width_in, final_height_in, ci$width, ci$height
))
