#!/usr/bin/env Rscript
# =====================================================================
# assemble_fig4_validation.R
# ---------------------------------------------------------------------
# Compose the CURRENT Fig 4 (validation) from the standalone panel PDFs in
# figures/main/fig4_validation/ (flat, descriptive names — no fig-number
# prefix; user organizes main/supp in Illustrator). This is a convenience
# montage, not the final layout.
#
#   a  proteomics mRNA-protein concordance   proteomics_concordance.pdf
#   b  spatial zonation overview (Visium)    spatial_zonation.pdf
#   c  NR triad benchmark (THRB/RORA/NR1H4)  nuclear_receptor_triad.pdf
#   d  HKDC1 NRF2 circuit                     hkdc1_nrf2_circuit.pdf
#   e  SERPINE1 fibroblast signal             serpine1_ligand.pdf
#   f  CYP3A4 zonation                        cyp3a4_zonation.pdf
#      FADS2 cross-ancestry (supp case)       fads2_cross_ancestry.pdf
#   g  cross-modal evidence matrix            target_evidence_matrix.pdf
#
# Layout (preserves each panel's native aspect ratio):
#   Row 1: a | b   Row 2: c   Row 3: d | e   Row 4: CYP3A4
#   Row 5: FADS2   Row 6: cross-modal evidence matrix (synthesis)
#
# Output: figures/main/fig4_validation/validation_composite.pdf
# Env:    rnaseq (uses ghostscript + magick, like assemble_fig4.R; rnaseq
#         lacks pdftools so we raster each page via gs then read with magick)
# =====================================================================

suppressPackageStartupMessages({
  library(cowplot)
  library(ggplot2)
  library(magick)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG4_DIR, "panels")  # panels moved under panels/ (matches fig2/fig3)
OUT_PDF   <- file.path(FIG4_DIR, "validation_composite.pdf")

# Descriptive, no figure-number prefix (user organizes main/supp in Illustrator).
panel_files <- c(
  a = "proteomics_concordance.pdf",
  b = "spatial_zonation.pdf",
  b2 = "spatial_svg_summary.pdf",
  c = "nuclear_receptor_triad.pdf",
  d = "hkdc1_nrf2_circuit.pdf",
  e = "serpine1_ligand.pdf",
  f1 = "cyp3a4_zonation.pdf",
  f2 = "fads2_cross_ancestry.pdf",
  g  = "target_evidence_matrix.pdf"
)

for (lbl in names(panel_files)) {
  fp <- file.path(PANEL_DIR, panel_files[[lbl]])
  if (!file.exists(fp)) stop("Missing panel ", lbl, ": ", fp)
}

read_panel <- function(fp, density = 300) {
  tmp_png <- tempfile(fileext = ".png")
  cmd <- sprintf(
    "gs -dQUIET -dBATCH -dNOPAUSE -dFirstPage=1 -dLastPage=1 -sDEVICE=pngalpha -r%d -sOutputFile=%s %s",
    density, shQuote(tmp_png), shQuote(fp))
  rc <- system(cmd)
  if (rc != 0 || !file.exists(tmp_png)) stop("gs render failed: ", fp)
  magick::image_read(tmp_png)
}

panel_to_gg <- function(img, label) {
  info   <- magick::image_info(img)
  aspect <- info$height / info$width
  p <- cowplot::ggdraw() +
    cowplot::draw_image(img, x = 0, y = 0, width = 1, height = 1) +
    cowplot::draw_label(label, x = 0.004, y = 0.996, hjust = 0, vjust = 1,
                        fontface = "plain", size = 6)
  attr(p, "aspect") <- aspect
  p
}

cat("Reading panels...\n")
imgs <- lapply(panel_files, function(f) read_panel(file.path(PANEL_DIR, f)))
names(imgs) <- names(panel_files)

# Panel-letter stamps (CYP3A4 + FADS2 share letter f as f / f cont.)
letters_map <- c(a = "a", b = "b", b2 = "", c = "c", d = "d", e = "e",
                 f1 = "f", f2 = "", g = "g")
panels_gg <- mapply(function(im, k) panel_to_gg(im, letters_map[[k]]),
                    imgs, names(imgs), SIMPLIFY = FALSE)
names(panels_gg) <- names(panel_files)

asp <- sapply(panels_gg, function(p) attr(p, "aspect"))

# ---- Layout -----------------------------------------------------------
# Composite full two-column width. Build each row at fixed width and a
# height that preserves the tallest panel's aspect in that row.
COMPOSITE_W <- 7.2  # inches

# Row 1: a | b  (give b ~62% of the width — it's the wider spatial panel)
r1_w  <- c(a = 0.38, b = 0.62) * COMPOSITE_W
r1_h  <- max(asp[["a"]] * r1_w[["a"]], asp[["b"]] * r1_w[["b"]])
row1  <- cowplot::plot_grid(panels_gg$a, panels_gg$b, ncol = 2,
                            rel_widths = r1_w, align = "none")

# Row 1b: spatial SVG summary (full width — pairs with the spatial overview)
r1b_h <- asp[["b2"]] * COMPOSITE_W
row1b <- panels_gg$b2

# Row 2: c  (full width, short)
r2_h  <- asp[["c"]] * COMPOSITE_W
row2  <- panels_gg$c

# Row 3: d | e  (HKDC1 wider/short, SERPINE1 narrower/tall)
r3_w  <- c(d = 0.60, e = 0.40) * COMPOSITE_W
r3_h  <- max(asp[["d"]] * r3_w[["d"]], asp[["e"]] * r3_w[["e"]])
row3  <- cowplot::plot_grid(panels_gg$d, panels_gg$e, ncol = 2,
                            rel_widths = r3_w, align = "none")

# Row 4: CYP3A4 (full)
r4_h  <- asp[["f1"]] * COMPOSITE_W
row4  <- panels_gg$f1

# Row 5: FADS2 (full)
r5_h  <- asp[["f2"]] * COMPOSITE_W
row5  <- panels_gg$f2

# Row 6: cross-modal evidence matrix (synthesis: 6 case-study genes x 6 modalities)
r6_h  <- asp[["g"]] * COMPOSITE_W
row6  <- panels_gg$g

row_heights <- c(r1_h, r1b_h, r2_h, r3_h, r4_h, r5_h, r6_h)
COMPOSITE_H <- sum(row_heights) + 0.2

cat(sprintf("Composite: %.2f x %.2f in (rows: %s)\n",
            COMPOSITE_W, COMPOSITE_H,
            paste(sprintf("%.2f", row_heights), collapse = " / ")))

composite <- cowplot::plot_grid(
  row1, row1b, row2, row3, row4, row5, row6,
  ncol = 1, rel_heights = row_heights, align = "none")

ggplot2::ggsave(OUT_PDF, composite,
                width = COMPOSITE_W, height = COMPOSITE_H,
                units = "in", device = cairo_pdf, limitsize = FALSE)
cat("Wrote: ", OUT_PDF, "\n", sep = "")
