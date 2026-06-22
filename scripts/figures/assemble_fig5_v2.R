# Fig 5 composite assembler v2 (2026-06-10) — convergence capstone.
#
# Assembles the C2-current Fig 5 from individual panel PDFs into a single
# multi-panel composite (Liang aesthetic). Writes a NEW versioned output
# (fig5_composite_v2.pdf) so the live fig5_composite.pdf is NOT clobbered.
#
# Panels (all verified C2-current — read canonical_deg_results.csv / bulk_*):
#   5a  Multi-Modal Convergence Matrix      (panels/fig5a.pdf)              — centerpiece
#   5b  Per-modality DEG coverage + n-active histogram (panels/fig5b.pdf)
#   5c  Drug-target orthogonality landscape  (panels/fig5c_drug_target_orthogonality.pdf)
#         (swapped in 2026-06-18, relocated from Fig 3; expression vs genetic
#          anchoring + clinical calibration. The old evidence-layer correlation
#          heatmap is demoted to figS08_subtyping_convergence.)
#   5d  Held-out convergence-score calibration (panels/fig5d_calibration.pdf) — NEW
#
# CUT 2026-06-19: the former "panel 5e" (TF convergence scatter + 4-way survival
# lollipop) was REMOVED from the main figure. It used a banned lollipop chart,
# oversold a "4-way TF convergence" claim the Currin caQTL data REFUTES (only 2
# direction-concordant 4-way variants exist, both in ZNF701 — not a disease gene;
# see figS05_scatac_caqtl_concordance), and its scatter x-axis (SCENIC+ regulon
# activity) was built on the null n=18 multiome SCENIC+ output. Source scripts +
# panel PDFs moved to panels/_legacy/.
#
# Layout (composite width 7.2 in):
#   Row 1 : 5a                         (full width, tall centerpiece)
#   Row 2 : 5b | 5c | 5d               (three even cells)
#
# Panel letters are baked into each panel's own title (a/b/c/d).
#
# Output: figures/main/fig5_convergence/fig5_composite_v2.pdf

suppressPackageStartupMessages({
  library(magick)
  library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANELS_DIR <- file.path(FIG5_DIR, "panels")
OUT_PDF    <- file.path(FIG5_DIR, "fig5_composite_v2.pdf")

panel_files <- c(
  a  = file.path(PANELS_DIR, "fig5a.pdf"),
  b  = file.path(PANELS_DIR, "fig5b.pdf"),
  c  = file.path(PANELS_DIR, "fig5c_drug_target_orthogonality.pdf"),
  d  = file.path(PANELS_DIR, "fig5d_calibration.pdf")
)
missing <- panel_files[!file.exists(panel_files)]
if (length(missing) > 0)
  stop("Missing panel PDFs:\n  ", paste(missing, collapse = "\n  "))

# --- Rasterize each panel PDF to high-DPI PNG via ghostscript (consistent with
#     assemble_fig5.R / assemble_fig3.R; avoids magick PDF-backend issues). ---
TMP <- tempdir()
rasterize <- function(lab, src, r = 300) {
  png_out <- file.path(TMP, paste0("fig5v2_", lab, ".png"))
  cmd <- sprintf(
    "gs -dQUIET -dBATCH -dNOPAUSE -dFirstPage=1 -dLastPage=1 -sDEVICE=pngalpha -r%d -sOutputFile=%s %s >/dev/null 2>&1",
    r, png_out, src)
  system(cmd)
  if (!file.exists(png_out))
    stop("ghostscript rasterization failed for: ", src)
  image_read(png_out)
}
imgs <- Map(rasterize, names(panel_files), panel_files)

# --- Assemble rows. Pad cells to a common height before horizontal append so
#     panels of different native aspect ratios align cleanly. ---
pad_to_height <- function(im, h) {
  info <- image_info(im)
  if (info$height == h) return(im)
  image_extent(im, geometry = paste0(info$width, "x", h),
               gravity = "north", color = "white")
}
resize_w <- function(im, w) image_resize(im, paste0(w, "x"))

# Row 1: 5a alone, scaled to full composite width.
FULL_W <- 2400  # px at 300 DPI -> 8 in working canvas; downscaled to 7.2 at end
row1 <- resize_w(imgs$a, FULL_W)

# Row 2: 5b | 5c | 5d  (three even thirds).
cell2_w <- round(FULL_W / 3)
r2 <- lapply(list(imgs$b, imgs$c, imgs$d), resize_w, w = cell2_w)
r2h <- max(sapply(r2, function(x) image_info(x)$height))
r2 <- lapply(r2, pad_to_height, h = r2h)
row2 <- image_append(do.call(c, r2))
row2 <- resize_w(row2, FULL_W)

# Stack rows (5a centerpiece + 5b|5c|5d). Panel 5e was cut 2026-06-19.
composite <- image_append(c(
  image_extent(row1, paste0(FULL_W, "x", image_info(row1)$height),
               gravity = "north", color = "white"),
  row2), stack = TRUE)
composite <- image_background(composite, "white", flatten = TRUE)

# --- Final page at 7.2 in wide. ---
ci <- image_info(composite)
final_w_in <- 7.2
final_h_in <- final_w_in * ci$height / ci$width

cairo_pdf(OUT_PDF, width = final_w_in, height = final_h_in, family = "Helvetica")
grid.newpage()
grid.raster(as.raster(composite), width = unit(1, "npc"), height = unit(1, "npc"))
dev.off()

cat(sprintf(
  "Wrote: %s (%.2f x %.2f in; composite raster %d x %d px)\n",
  OUT_PDF, final_w_in, final_h_in, ci$width, ci$height))
