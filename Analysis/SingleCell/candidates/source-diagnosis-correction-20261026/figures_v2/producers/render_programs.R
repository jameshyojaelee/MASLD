#!/usr/bin/env Rscript
# KEY MESSAGE: When each dataset keeps its own residual variance and the two uncorrected
# GSE189600 donors are removed (H59), no prespecified program passes the 117-program family;
# the panel shows the donor-level stage pattern and per-program GLS stage effects without a call.
# Candidate for Fig 4C (24 rows) and Fig S4C (117 rows). Copied and adapted from
# Repair 2026-10-15: forest-axis ticks limited to the plotted range (the S4C "-3" label sat
# outside the forest and overlapped the rotated "Steatohepatitis" column name); 4 mm
# annotation padding. Values and row choice are unchanged; outputs go to repair/.
# scripts/figures/singlecell_module_heatmap.R (draw_modules); repo script not modified.
suppressPackageStartupMessages({library(data.table); library(ComplexHeatmap); library(circlize); library(grid)})
grDevices::pdf.options(useDingbats = FALSE)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
FIG  <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/figures_v2"
OUT <- FIG
source(file.path(BASE, "scripts/figures/publication_theme.R"))
STAGES <- c("Healthy", "Steatosis", "Steatohepatitis")
CTS    <- c("hepatocytes", "fibroblasts", "macrophages", "cholangiocytes", "tcells")
CT_LABEL <- c(hepatocytes = "Hepatocyte", fibroblasts = "Fibroblast", macrophages = "Macrophage",
              cholangiocytes = "Cholangiocyte", tcells = "T cell")
CT_SHORT <- c(hepatocytes = "Hep", fibroblasts = "Fib", macrophages = "Mac", cholangiocytes = "Chol", tcells = "T")
CT_FULL  <- c(hepatocytes = "Hepatocytes", fibroblasts = "Fibroblasts", macrophages = "Macrophages",
              cholangiocytes = "Cholangiocytes", tcells = "T cells")
CT_COL <- setNames(unname(ct_palette[CT_FULL[CTS]]), CT_LABEL[CTS])
G6 <- function(...) gpar(fontsize = 6, fontfamily = "Helvetica", ...)

W <- fread(file.path(FIG, "source_tables/fig4c_h59_registry_117.tsv"))
crit <- fread(file.path(BASE, "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/programs/null_critical_abs_t.tsv"))
W[, program := paste(cell_type, module, sep = "_")]
W <- merge(W, crit, by = "program", all.x = TRUE)
stopifnot(nrow(W) == 117L, !anyNA(W$null_abs_t_crit), all(is.finite(W$q)))
W[, `:=`(fw_lo = beta - null_abs_t_crit * SE, fw_hi = beta + null_abs_t_crit * SE)]
W[, cell_type := factor(cell_type, levels = CTS)]
fwrite(W, file.path(OUT, "source_tables/fig4c_h59_registry_117_with_familywise_interval.tsv"), sep = "\t")
cat(sprintf("[check] nominal CI excludes 0: %d/117; 117-program interval excludes 0: %d/117\n",
            sum(W$lo > 0 | W$hi < 0), sum(W$fw_lo > 0 | W$fw_hi < 0)))

LIM  <- as.numeric(quantile(abs(as.matrix(W[, ..STAGES])), 0.98, na.rm = TRUE))
COLF <- colorRamp2(c(-LIM, 0, LIM), c("#1565C0", "white", "#C9265E"))
GAP_MM <- 0.6
ht_opt(ROW_ANNO_PADDING = unit(4, "mm"))

draw_modules <- function(Wsub, out_pdf, pw, ph, row_mm) {
  Wsub <- copy(Wsub); setorder(Wsub, cell_type, -beta)
  mat <- as.matrix(Wsub[, ..STAGES])
  nm <- sub(" \\([^)]*\\)$", "", Wsub$name)
  rownames(mat) <- sprintf("%d · %s", Wsub$module, nm)
  present <- CTS[CTS %in% as.character(Wsub$cell_type)]
  row_split <- factor(CT_LABEL[as.character(Wsub$cell_type)], levels = CT_LABEL[present])
  CTv <- as.character(Wsub$cell_type)
  bracket_fun <- function(index, k, n) {
    lab <- unname(CT_SHORT[CTv[index[1]]])
    grid.segments(unit(0.78, "npc"), unit(0.04, "npc"), unit(0.78, "npc"), unit(0.96, "npc"), gp = gpar(col = "grey25", lwd = 0.6))
    grid.segments(unit(0.78, "npc"), unit(c(0.04, 0.96), "npc"), unit(0.98, "npc"), unit(c(0.04, 0.96), "npc"), gp = gpar(col = "grey25", lwd = 0.6))
    grid.text(lab, x = unit(0.02, "npc"), y = unit(0.5, "npc"), just = "left", gp = G6())
  }
  left_anno <- rowAnnotation(
    lineage = AnnotationFunction(fun = bracket_fun, which = "row", width = unit(8.5, "mm"), var_import = list(CTv = CTv)),
    ct = CT_LABEL[CTv], col = list(ct = CT_COL), show_legend = FALSE, show_annotation_name = FALSE,
    simple_anno_size = unit(1.6, "mm"), gap = unit(0.4, "mm"))
  BV <- Wsub$beta; LOv <- Wsub$lo; HIv <- Wsub$hi; FLv <- Wsub$fw_lo; FHv <- Wsub$fw_hi
  xr <- range(c(FLv, FHv, 0), na.rm = TRUE); xr <- xr + c(-0.06, 0.06) * diff(xr)
  ci_fun <- function(index, k, n) {
    m <- length(index); yy <- m:1
    pushViewport(viewport(xscale = xr, yscale = c(0.5, m + 0.5)))
    grid.segments(unit(0, "native"), unit(0, "npc"), unit(0, "native"), unit(1, "npc"), gp = gpar(col = "#9E9E9E", lwd = 0.4))
    for (i in seq_len(m)) {
      j <- index[i]; y <- yy[i]
      grid.segments(unit(FLv[j], "native"), unit(y, "native"), unit(FHv[j], "native"), unit(y, "native"),
                    gp = gpar(col = "grey82", lwd = 2.2, lineend = "butt"))
      grid.segments(unit(LOv[j], "native"), unit(y, "native"), unit(HIv[j], "native"), unit(y, "native"),
                    gp = gpar(col = "grey15", lwd = 0.6))
      fill <- if (BV[j] > 0) "#C9265E" else "#1565C0"
      grid.points(unit(BV[j], "native"), unit(y, "native"), pch = 21, size = unit(1.1, "mm"),
                  gp = gpar(fill = fill, col = "grey15", lwd = 0.3))
    }
    if (k == n) { at <- pretty(xr, 4); at <- at[at >= xr[1] & at <= xr[2]]
      grid.lines(unit(xr, "native"), unit(c(0, 0), "npc"), gp = gpar(col = "black", lwd = 0.5))
      grid.xaxis(at = at, gp = gpar(fontsize = 6, fontfamily = "Helvetica", lwd = 0.5)) }
    popViewport()
  }
  right_anno <- rowAnnotation(
    "GLS stage β" = AnnotationFunction(fun = ci_fun, which = "row", width = unit(1.7, "cm"),
      var_import = list(BV = BV, LOv = LOv, HIv = HIv, FLv = FLv, FHv = FHv, xr = xr)),
    annotation_name_gp = G6(), annotation_name_rot = 0, annotation_name_side = "bottom",
    annotation_name_offset = unit(4, "mm"))
  ht <- Heatmap(mat, name = "Module score\n(centered)", col = COLF, na_col = "grey92",
    cluster_rows = FALSE, cluster_columns = FALSE, row_split = row_split, row_gap = unit(GAP_MM, "mm"),
    height = unit(row_mm * nrow(mat), "mm"), width = unit(1.35, "cm"), row_title = NULL,
    row_names_side = "left", row_names_gp = G6(), column_names_gp = G6(), column_names_rot = 45,
    column_names_side = "bottom", column_title = NULL, left_annotation = left_anno, right_annotation = right_anno,
    rect_gp = gpar(col = "white", lwd = 0.3),
    heatmap_legend_param = list(title_gp = G6(), labels_gp = G6(), grid_width = unit(2.5, "mm")))
  int_lgd <- Legend(labels = c("nominal 95% CI", "117-program interval"), title = "GLS stage β", type = "lines",
                    legend_gp = gpar(col = c("grey15", "grey82"), lwd = c(0.6, 2.2)),
                    title_gp = G6(), labels_gp = G6(), grid_width = unit(3, "mm"))
  cairo_pdf(out_pdf, width = pw, height = ph)
  draw(ht, merge_legends = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right",
       annotation_legend_list = list(int_lgd), padding = unit(c(5, 1, 5, 1), "mm"))
  dev.off()
  cat(sprintf("[saved] %s (%d rows, %.2f x %.2f in)\n", out_pdf, nrow(mat), pw, ph))
}

main <- W[rk <= 3 & cell_type != "tcells"]
stopifnot(nrow(main) == 24L)
draw_modules(main, file.path(OUT, "fig4c_hotspot_stage_heatmap_H59_source_qualified.pdf"), pw = 4.0, ph = 3.0, row_mm = 2.0)
nblk <- 5L; ph_full <- (2.4 * 117 + GAP_MM * (nblk - 1)) / 25.4 + 1.05
draw_modules(W, file.path(OUT, "figs4c_hotspot_stage_heatmap_full_H59_source_qualified.pdf"), pw = 4.85, ph = ph_full, row_mm = 2.4)
writeLines(capture.output(sessionInfo()), file.path(OUT, "producers/sessionInfo_02r_render.txt"))
