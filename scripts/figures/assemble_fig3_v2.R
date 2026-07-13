# Fig 3 (Multi-Ancestry Regulatory Architecture) composite assembler — v2.
# 2026-05-18: 9-panel canonical layout, 3 x 3 grid, integrates the 3 NEW ATAC
# panels (K exemplar -> J regulon heatmap -> L funnel) per legend reorder.
#
# Canonical panel selection (matches docs/manuscript/05_figure_legends.md Fig 3,
# pruned to 9 panels for compactness — legend panels (g) high-PIP TF list +
# (i) D4 cell-type heritability remain in the figure narrative as inline
# references; legend (j)/(k) are promoted into this composite as h/g):
#   a: COLOC Manhattan          (pip_coloc_manhattan.pdf)
#   b: Top COLOC heatmap        (ancestry_coloc_counts.pdf)
#   c: sc-eQTL cell-type dot    (sceqtl_celltype_dotplot.pdf)
#   d: Causal method coverage   (causal_coverage.pdf)
#   e: ieQTL-DEG concordance    (ieqtl_concordance.pdf)
#   f: Multi-ancestry COLOC     (multiancestry_coloc.pdf)
#   g: GWAS-ATAC chain exemplar (fig3_panel_K_gwas_atac_chain.pdf)   [BOTTOM-UP; not produced in this dir — dangling]
#   h: Disease regulons x CT    (disease_master_regulators_celltype.pdf) [TOP-DOWN]
#   i: Evidence funnel 24 -> 5  (evidence_funnel.pdf)
#
# Legacy duplicates excluded:
#   fig3a_hybrid.pdf, fig3a_hybrid_b.pdf, fig3a_coloc_manhattan.pdf
#   fig3c.pdf (May-17 alt)
#   fig3d.pdf (May-17 alt)
#   fig3e.pdf, fig3e_tf_heatmap.pdf
#   fig3f.pdf
#   fig3g_coloc_class_split.pdf, fig3g_cross_ancestry_pp4.pdf, fig3g_twas_concordance.pdf
#   fig3h_deg_coloc_funnel.pdf, fig3i_deg_coloc_scatter.pdf, fig3i_sctwas_dotplot.pdf
#   fig3_intro_*.pdf (intro variants superseded)

suppressPackageStartupMessages({
  library(magick)
  library(grid)
  library(gridExtra)
})

# DISABLED 2026-06-12: the Fig 2 (genetics) composite is retired — only the individual
# panels are kept (per user). This assembler no longer generates genetics_composite.pdf.
message("[assemble_fig3_v2] disabled — genetics_composite.pdf (Fig 2 composite) retired (2026-06-12)")
quit(save = "no", status = 0)

source("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/load_figure_data.R")

PANEL_DIR <- file.path(FIG3_DIR, "panels")
OUT_PDF   <- file.path(FIG3_DIR, "genetics_composite.pdf")

# Ordered (letter, file) tuples
panels <- list(
  a = "pip_coloc_manhattan.pdf",
  b = "ancestry_coloc_counts.pdf",
  c = "sceqtl_celltype_dotplot.pdf",
  d = "causal_coverage.pdf",
  e = "ieqtl_concordance.pdf",
  f = "multiancestry_coloc.pdf",
  g = "fig3_panel_K_gwas_atac_chain.pdf",
  h = "disease_master_regulators_celltype.pdf",
  i = "evidence_funnel.pdf"
)

# Sanity-check all panels exist
for (lab in names(panels)) {
  fp <- file.path(PANEL_DIR, panels[[lab]])
  if (!file.exists(fp)) stop("Missing panel ", lab, ": ", fp)
}

# Load each panel via Ghostscript rasterization at 300 dpi (pdftools not available in rnaseq env)
TMP <- tempdir()
load_panel <- function(fp, tag) {
  png_path <- file.path(TMP, paste0(tag, ".png"))
  cmd <- sprintf(
    "gs -dQUIET -dBATCH -dNOPAUSE -dFirstPage=1 -dLastPage=1 -sDEVICE=pngalpha -r300 -sOutputFile=%s %s",
    shQuote(png_path), shQuote(fp))
  status <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (status != 0 || !file.exists(png_path)) stop("gs rasterization failed for ", fp)
  image_read(png_path)
}
imgs <- mapply(function(p, lab) load_panel(file.path(PANEL_DIR, p), tag = paste0("p_", lab)),
               panels, names(panels), SIMPLIFY = FALSE)
names(imgs) <- names(panels)

# Composite dimensions: 7.2 inches wide; 3-col grid -> ~2.4 inches per col
# Maintain individual panel aspect ratio; standardize to fixed cell width.
COMP_W_IN <- 7.2
N_COL     <- 3
N_ROW     <- 3
CELL_W_IN <- COMP_W_IN / N_COL          # 2.4"
PAD_IN    <- 0.05                        # minimal padding between cells

# Resize each image to a fixed cell width (in pixels at 300 dpi); preserve aspect
DPI <- 300
CELL_W_PX <- round(CELL_W_IN * DPI)

resize_to_width <- function(img, w_px) {
  image_resize(img, paste0(w_px, "x"))
}
imgs_rs <- lapply(imgs, resize_to_width, w_px = CELL_W_PX)

# Determine per-row max height (so rows are uniform but columns natural)
heights_px <- sapply(imgs_rs, function(im) image_info(im)$height)
row_max_h <- sapply(0:(N_ROW-1), function(r) {
  idx <- (r * N_COL + 1):((r + 1) * N_COL)
  max(heights_px[idx])
})

# Pad each image to its row's max height (white background, top-aligned)
pad_to_height <- function(img, target_h) {
  info <- image_info(img)
  if (info$height >= target_h) return(img)
  # Add bottom padding
  image_extent(img,
               geometry = paste0(info$width, "x", target_h),
               gravity  = "north",
               color    = "white")
}

idx_lin <- 1
imgs_padded <- imgs_rs
for (r in 0:(N_ROW-1)) {
  for (c in 0:(N_COL-1)) {
    imgs_padded[[idx_lin]] <- pad_to_height(imgs_rs[[idx_lin]], row_max_h[r + 1])
    idx_lin <- idx_lin + 1
  }
}

# Build rows: append columns horizontally with a thin white spacer
PAD_PX <- round(PAD_IN * DPI)
spacer_h <- function(w_px, h_px) image_blank(width = w_px, height = h_px, color = "white")

append_with_pad <- function(im_list, pad_w) {
  # interleave with horizontal padding
  out <- im_list[[1]]
  for (k in 2:length(im_list)) {
    pad <- spacer_h(pad_w, image_info(out)$height)
    out <- image_append(c(out, pad, im_list[[k]]), stack = FALSE)
  }
  out
}

rows <- vector("list", N_ROW)
for (r in 0:(N_ROW-1)) {
  row_imgs <- imgs_padded[(r*N_COL+1):((r+1)*N_COL)]
  rows[[r+1]] <- append_with_pad(row_imgs, PAD_PX)
}

# Standardize row widths to the max row width
max_rw <- max(sapply(rows, function(r) image_info(r)$width))
rows_unif <- lapply(rows, function(r) {
  info <- image_info(r)
  if (info$width == max_rw) return(r)
  image_extent(r,
               geometry = paste0(max_rw, "x", info$height),
               gravity  = "west",
               color    = "white")
})

# Stack rows with vertical padding
stack_with_pad <- function(im_list, pad_h) {
  out <- im_list[[1]]
  for (k in 2:length(im_list)) {
    pad <- image_blank(width = image_info(out)$width, height = pad_h, color = "white")
    out <- image_append(c(out, pad, im_list[[k]]), stack = TRUE)
  }
  out
}
composite_img <- stack_with_pad(rows_unif, PAD_PX)

# Overlay bold A-I labels at top-left of each cell.
# Determine cell origins on the composite image (post-resize/pad).
cell_w_px_final <- image_info(rows_unif[[1]])$width
# Per-column x origin: compute cumulative widths from imgs_padded row 1
row1_widths <- sapply(imgs_padded[1:N_COL], function(im) image_info(im)$width)
col_x_origin <- cumsum(c(0, row1_widths[-N_COL] + PAD_PX))
# Per-row y origin
row_heights_final <- sapply(rows_unif, function(r) image_info(r)$height)
row_y_origin <- cumsum(c(0, row_heights_final[-N_ROW] + PAD_PX))

LABELS <- c("A","B","C","D","E","F","G","H","I")

# Write the composite as a PDF using grid + cairo_pdf (vector wrapper),
# with A-I labels overlaid as native PDF text (not raster) for crispness.
final_info <- image_info(composite_img)
final_w_in <- final_info$width / DPI
final_h_in <- final_info$height / DPI

cat(sprintf("[fig3] composite raster: %d x %d px -> %.2f x %.2f in\n",
            final_info$width, final_info$height, final_w_in, final_h_in))

# Cell origins in "inches from top-left"
col_x_origin_in <- col_x_origin / DPI
row_y_origin_in <- row_y_origin / DPI

cairo_pdf(OUT_PDF, width = final_w_in, height = final_h_in)
grid.newpage()
grid.raster(as.raster(composite_img),
            x = unit(0.5, "npc"), y = unit(0.5, "npc"),
            width = unit(1, "npc"), height = unit(1, "npc"),
            interpolate = TRUE)

# Overlay native PDF text labels A-I (top-left of each cell, with small inset)
INSET_IN <- 0.04
idx_lin <- 1
for (r in 0:(N_ROW-1)) {
  for (c in 0:(N_COL-1)) {
    x_in <- col_x_origin_in[c + 1] + INSET_IN
    # grid y-axis runs bottom-up; convert
    y_in_from_top <- row_y_origin_in[r + 1] + INSET_IN
    y_in <- final_h_in - y_in_from_top
    grid.text(LABELS[idx_lin],
              x = unit(x_in, "in"),
              y = unit(y_in, "in"),
              just = c("left", "top"),
              gp = gpar(fontsize = 12, fontface = "plain",
                        fontfamily = "Helvetica"))
    idx_lin <- idx_lin + 1
  }
}
dev.off()

cat(sprintf("[fig3] wrote %s\n", OUT_PDF))
cat(sprintf("[fig3] panels: a=%s b=%s c=%s d=%s e=%s f=%s g=%s h=%s i=%s\n",
            panels$a, panels$b, panels$c, panels$d, panels$e,
            panels$f, panels$g, panels$h, panels$i))
