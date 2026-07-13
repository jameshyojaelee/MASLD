# Alternate Fig 3 composite assembler — magick-based stitch from existing PDFs.
# Layout (matches fig3_regulatory_architecture_v2.R composite, 2026-04-29 cascade):
#   Row 1: pip_coloc_manhattan         (full width)
#   Row 2: ancestry_coloc_counts | rora_locus_zoom      (per-ancestry COLOC bars | RORA locus LD-zoom)
#   Row 3: pip_vs_disruption_scatter | drug_target_validation   (GWAS-ATAC scatter | drug-target validation)
#   Row 4: susiex_cross_ancestry       (finemapping logFC × COLOC PP4)
suppressPackageStartupMessages(library(magick))

# DISABLED 2026-06-12: regulatory_architecture.pdf composite retired from Fig 2 (per user).
# Two of its input panels (pip_coloc_manhattan, pip_vs_disruption_scatter) were also retired.
# Script preserved for reference; it no longer generates output.
message("[assemble_fig3] disabled — regulatory_architecture.pdf retired (2026-06-12)")
quit(save = "no", status = 0)

DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/main/fig2_genetics"
TMP <- tempdir()

panels <- c("pip_coloc_manhattan","ancestry_coloc_counts","rora_locus_zoom",
            "pip_vs_disruption_scatter","drug_target_validation","susiex_cross_ancestry")
for (p in panels) {
  system(sprintf("gs -dQUIET -dBATCH -dNOPAUSE -dFirstPage=1 -dLastPage=1 -sDEVICE=pngalpha -r150 -sOutputFile=%s/%s.png %s/panels/%s.pdf >/dev/null 2>&1",
                 TMP, p, DIR, p))
}
imgs <- lapply(panels, function(p) image_read(file.path(TMP, paste0(p, ".png"))))
names(imgs) <- panels

W <- 1600
imgs_rs <- lapply(imgs, function(im) image_resize(im, paste0(W, "x")))

row1 <- imgs_rs[["pip_coloc_manhattan"]]
row2 <- image_append(c(imgs_rs[["ancestry_coloc_counts"]], imgs_rs[["rora_locus_zoom"]]))
row3 <- image_append(c(imgs_rs[["pip_vs_disruption_scatter"]], imgs_rs[["drug_target_validation"]]))
row4 <- imgs_rs[["susiex_cross_ancestry"]]

RW <- max(sapply(list(row1, row2, row3, row4), function(r) image_info(r)$width))
rows_rs <- lapply(list(row1, row2, row3, row4), function(r) image_resize(r, paste0(RW, "x")))
combined <- image_append(do.call(c, rows_rs), stack = TRUE)

out <- file.path(DIR, "regulatory_architecture.pdf")
image_write(combined, path = out, format = "pdf")
info <- image_info(combined)
cat("Wrote:", out, " (", info$width, "x", info$height, ")\n")
