#!/usr/bin/env Rscript
# fig2_locus_legend_strip.R
# Shared horizontal legend for the Fig2F / Fig2G locus-zoom pair (the legend elements
# were stripped from Fig2G and consolidated here). One row, left-to-right:
#   r² gradient bar | Lead SNP diamond | EUR-in-CS / not-in-CS / coloc-shared-variant
#   markers | log2FC gradient bar | n.s. swatch
# 6pt Helvetica, cairo_pdf (embeds Helvetica), fixed page (no bbox trim). ~6.5 x 0.45 in.
# Output: figures/main/fig2_genetics/panels/Fig2H_locus_legend.pdf
suppressPackageStartupMessages({ library(plotgardener); library(grid) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT  <- file.path(BASE, "figures/main/fig2_genetics/panels/Fig2H_locus_legend.pdf")
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

PAGE_W <- 6.61; PAGE_H <- 0.44   # = Fig2F (3.10) + Fig2G (3.50) combined width
YC <- PAGE_H / 2                       # single-row vertical centre
t2g <- function(y) PAGE_H - y          # top-origin -> grid bottom-origin

# Palettes — identical to the locus-zoom panel (figS09_locus_zoom_disease.R)
ld_palette  <- colorRampPalette(c("#E6E6E5", "#1B2F5B"))                                   # r² to lead: grey -> navy
deg_palette <- colorRampPalette(c("#1565C0", "#90CAF9", "#FFFFFF", "#F48FB1", "#C9265E"))  # log2FC: blue -> white -> magenta

cairo_pdf(OUT, width = PAGE_W, height = PAGE_H, family = "Helvetica")
pageCreate(width = PAGE_W, height = PAGE_H, default.units = "inches",
           showGuides = FALSE, xgrid = 0, ygrid = 0)

# --- drawing helpers (advance an x cursor, in inches from the left) --------------
txt_w <- function(s) convertWidth(grobWidth(textGrob(s, gp = gpar(fontsize = 6, fontfamily = "Helvetica"))),
                                  "inches", valueOnly = TRUE)
put_txt <- function(s, x) {
  grid.text(s, unit(x, "in"), unit(t2g(YC), "in"), just = c("left", "center"),
            gp = gpar(fontsize = 6, fontfamily = "Helvetica", col = "black"))
  x + txt_w(s)
}
put_marker <- function(x, pch, fill, col = "black", size = 0.10, lwd = 0.7) {
  grid.points(unit(x + size/2, "in"), unit(t2g(YC), "in"), pch = pch,
              size = unit(size, "in"), gp = gpar(col = col, fill = fill, lwd = lwd))
  x + size
}
put_gradbar <- function(x, palfun, w = 0.5, h = 0.11) {
  cols <- palfun(50)
  for (k in 0:49)
    grid.rect(unit(x + k * (w/50), "in"), unit(t2g(YC), "in"),
              width = unit(w/50 + 0.004, "in"), height = unit(h, "in"),
              just = c("left", "center"), gp = gpar(col = NA, fill = cols[k + 1]))
  x + w
}
GAP <- 0.16   # gap between legend groups
cur <- 0.10

# 1) r² to lead gradient
cur <- put_txt("r² to lead", cur) + 0.06
cur <- put_txt("0", cur) + 0.03
cur <- put_gradbar(cur, ld_palette) + 0.03
cur <- put_txt("1", cur) + GAP

# 2) Lead SNP diamond
cur <- put_marker(cur, 23, "#C2185B", size = 0.11, lwd = 0.8) + 0.05
cur <- put_txt("Lead SNP", cur) + GAP

# 3) SuSiE PIP markers
cur <- put_marker(cur, 19, "#1565C0", col = "#1565C0", size = 0.09) + 0.05
cur <- put_txt("EUR (in CS)", cur) + 0.11
cur <- put_marker(cur, 21, "white", col = "#9E9E9E", size = 0.07) + 0.05
cur <- put_txt("not in CS", cur) + 0.11
cur <- put_marker(cur, 23, "#FFD600", size = 0.10, lwd = 0.7) + 0.05
cur <- put_txt("coloc shared variant", cur) + GAP

# 4) log2FC gradient
cur <- put_txt("log2FC", cur) + 0.06
cur <- put_txt("−0.3", cur) + 0.03
cur <- put_gradbar(cur, deg_palette) + 0.03
cur <- put_txt("+0.3", cur) + GAP

# 5) n.s. swatch
grid.rect(unit(cur + 0.05, "in"), unit(t2g(YC), "in"), width = unit(0.10, "in"),
          height = unit(0.10, "in"), just = c("left", "center"), gp = gpar(col = NA, fill = "#9E9E9E"))
cur <- cur + 0.12
cur <- put_txt("n.s.", cur)

invisible(dev.off())
cat(sprintf("[saved] %s  (%.2f x %.2f in; content ends at %.2f in)\n", OUT, PAGE_W, PAGE_H, cur))
