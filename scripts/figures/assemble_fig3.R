# Alternate Fig 3 composite assembler — magick-based stitch from existing PDFs.
# Layout (matches fig3_regulatory_architecture_v2.R composite, 2026-04-29 cascade):
#   Row 1: fig3a (full width)
#   Row 2: fig3b | fig3c     (per-ancestry COLOC bars | RORA locus LD-zoom)
#   Row 3: fig3d | fig3e     (GWAS-ATAC scatter | high-PIP TF heatmap)
#   Row 4: fig3f             (drug-target logFC × COLOC PP4)
suppressPackageStartupMessages(library(magick))

DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/main/fig3_regulatory_architecture"
TMP <- tempdir()

panels <- c("fig3a","fig3b","fig3c","fig3d","fig3e","fig3f")
for (p in panels) {
  system(sprintf("gs -dQUIET -dBATCH -dNOPAUSE -dFirstPage=1 -dLastPage=1 -sDEVICE=pngalpha -r150 -sOutputFile=%s/%s.png %s/panels/%s.pdf >/dev/null 2>&1",
                 TMP, p, DIR, p))
}
imgs <- lapply(panels, function(p) image_read(file.path(TMP, paste0(p, ".png"))))
names(imgs) <- panels

W <- 1600
imgs_rs <- lapply(imgs, function(im) image_resize(im, paste0(W, "x")))

row1 <- imgs_rs$fig3a
row2 <- image_append(c(imgs_rs$fig3b, imgs_rs$fig3c))
row3 <- image_append(c(imgs_rs$fig3d, imgs_rs$fig3e))
row4 <- imgs_rs$fig3f

RW <- max(sapply(list(row1, row2, row3, row4), function(r) image_info(r)$width))
rows_rs <- lapply(list(row1, row2, row3, row4), function(r) image_resize(r, paste0(RW, "x")))
combined <- image_append(do.call(c, rows_rs), stack = TRUE)

out <- file.path(DIR, "fig3_regulatory_architecture.pdf")
image_write(combined, path = out, format = "pdf")
info <- image_info(combined)
cat("Wrote:", out, " (", info$width, "x", info$height, ")\n")
