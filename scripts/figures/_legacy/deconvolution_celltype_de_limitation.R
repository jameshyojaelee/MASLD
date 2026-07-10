#!/usr/bin/env Rscript
# deconvolution_celltype_de_limitation.R
# Fig S3 (supplementary diagnostic) — why bulk deconvolution cannot supply a
# multi-lineage cell-type-DE y-axis for the Fig 3E disconnect panel.
# Per lineage: mean tissue proportion (MuSiC) vs whether TOAST could estimate a
# trustworthy disease effect. Only lineages above bulk's ~1% detection floor are
# testable, and of those only hepatocytes (92% of reads) give a plausible beta;
# the rare non-parenchymal lineages inflate (Myeloid beta ~3.5 = noise). Motivates
# using scRNA pseudobulk (direct) rather than deconvolution for the disconnect panel.
#
# Output: FIG2_DIR/panels/figs3_deconvolution_celltype_de_limitation.pdf

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG2_DIR, "panels")
lab_size  <- 6 / ggplot2::.pt

# per-lineage tissue proportion (aggregate MuSiC fine -> 7 lineages)
prop <- fread(file.path(BASE, "RNA-seq/results/celltype_attribution/persample_celltype_proportions.csv"))
LMAP <- list(Hepatocytes="Hepatocytes", Cholangiocytes="Cholangiocytes",
             Endothelial="Endothelial cells", Mesenchymal="Fibroblasts",
             Myeloid=c("Macrophages","Kupffer cells","Mono+mono derived cells",
                       "Monocytes & Monocyte-derived cells","cDC1s","cDC2s","pDCs",
                       "Mig. cDCs","Basophils","Neutrophils"),
             `T / NK`=c("T cells","Circulating NK/NKT","Resident NK","NK cells","ILC1s"),
             `B / Plasma`=c("B cells","Plasma cells"))
ab <- rbindlist(lapply(names(LMAP), function(ln) {
  cols <- intersect(LMAP[[ln]], names(prop))
  data.table(lineage = ln, prop = mean(rowSums(as.matrix(prop[, ..cols]), na.rm=TRUE), na.rm=TRUE))
}))

# TOAST betas (which lineages were estimable + their effect magnitude)
toast <- tryCatch(fread(file.path(BASE,
  "RNA-seq/results/celltype_attribution/coarse_lineage/toast_coarse_lineage_summary.csv")),
  error = function(e) data.table(cell_type=character(), mean_abs_beta=numeric()))
tmap <- c(Hepatocytes="Hepatocytes", Myeloid="Myeloid", T_NK="T / NK",
          Endothelial="Endothelial", Mesenchymal="Mesenchymal",
          Cholangiocytes="Cholangiocytes", B_Plasma="B / Plasma")
toast[, lineage := tmap[cell_type]]
ab <- merge(ab, toast[, .(lineage, toast_beta = mean_abs_beta)], by = "lineage", all.x = TRUE)
ab[, status := fifelse(is.na(toast_beta), "not estimable from bulk (rare / collinear)",
              fifelse(lineage == "Hepatocytes", "estimable & plausible",
                      "estimable but noise-inflated"))]
ab[, lineage := factor(lineage, levels = ab[order(prop)]$lineage)]
scol <- c("not estimable from bulk (rare / collinear)"="#BDBDBD",
          "estimable & plausible"="#1565C0", "estimable but noise-inflated"="#C9265E")

p <- ggplot(ab, aes(x = prop, y = lineage, fill = status)) +
  geom_vline(xintercept = 0.01, linetype = "dashed", linewidth = 0.3, colour = "gray50") +
  geom_col(width = 0.68) +
  geom_text(aes(label = ifelse(is.na(toast_beta), "",
                sprintf("β=%.2f", toast_beta))), hjust = -0.15, size = lab_size, colour = "gray20") +
  annotate("text", x = 0.0105, y = 0.6, label = "deconvolution\ndetection floor ~1%",
           hjust = 0, vjust = 0, size = lab_size, colour = "gray45", lineheight = 0.9) +
  scale_x_log10(labels = scales::percent_format(accuracy = 0.1),
                expand = expansion(mult = c(0.02, 0.22))) +
  scale_fill_manual(values = scol, name = NULL) +
  labs(x = "Mean tissue proportion (bulk deconvolution, log scale)", y = NULL) +
  theme_masld_compact() +
  theme(legend.position = "top", legend.key.size = unit(0.28,"cm"))

message("CAPTION Fig S3 (deconvolution limitation): Bulk deconvolution cannot supply a ",
        "multi-lineage cell-type-DE axis. Bars = mean tissue proportion per lineage (MuSiC); ",
        "dashed line = the ~1% detection floor below which cell-type-specific DE is not ",
        "estimable from bulk. Only Hepatocytes (92%) yield a plausible TOAST disease effect ",
        "(mean|β|=0.22); Myeloid (5%) inflates to |β|=3.5 (43% of genes spuriously 'significant'); ",
        "the remaining lineages (<1%) are not estimable at all. Hence Fig 3E uses scRNA pseudobulk ",
        "(direct per-lineage measurement) rather than deconvolution.")
out <- file.path(PANEL_DIR, "figs3_deconvolution_celltype_de_limitation.pdf")
save_fig(p, out, width = fig_half_width + 0.8, height = 2.4)
cat("[figs3-deconv-limit] Saved:", out, "\n"); print(ab[order(-prop)])
