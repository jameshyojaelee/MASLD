#!/usr/bin/env Rscript
# ============================================================================
# singlecell_module_heatmap_donorlevel.R  — DONOR-LEVEL variant of Fig 3G.
#
# Run-vs-donor pseudoreplication fix for the Hotspot module remodeling heatmap.
# Identical rendering to singlecell_module_heatmap.R EXCEPT both run-level
# statistics are recomputed at the biological-donor level:
#   (a) MODULE SELECTION reads disease_stage_q from all_modules_donor.tsv
#       (donor-level recompute by 511_donorlevel_disease_stage.R), NOT the
#       run-level all_modules.tsv.
#   (b) The in-script disease-slope beta+95%CI collapses each per-cell-type
#       donor_scores.tsv from RUN level (atlas `sample` = a sequencing run /
#       sort-fraction library) to one score per TRUE biological donor (mean
#       across the donor's runs, via lib_donor_collapse::build_srr_to_donor_map)
#       BEFORE lmer(score ~ ord + (1|dataset)).
# Writes to _DONORLEVEL filenames; the canonical fig3g PDF is left untouched.
# The fig3h LIANA/CCC panel is NOT regenerated here (separate donor-collapse track).
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ComplexHeatmap); library(circlize); library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
source(file.path(BASE, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))

HS        <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
ALL_MOD   <- file.path(HS, "donor_collapse/all_modules_donor.tsv")   # donor-level selection stats
META_F    <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv")
PANEL_DIR <- file.path(FIG4_SC_DIR, "panels", "_legacy_direct_render")
DATA_DIR  <- file.path(PANEL_DIR, "data")
OUT_MAIN  <- file.path(PANEL_DIR, "fig3g_singlecell_module_heatmap_DONORLEVEL.pdf")
OUT_SUPP  <- file.path(PANEL_DIR, "figs3_singlecell_module_heatmap_full_DONORLEVEL.pdf")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

STAGES   <- c("Healthy", "Steatosis", "Steatohepatitis")
CTS      <- c("hepatocytes", "fibroblasts", "macrophages", "cholangiocytes", "tcells")
CT_LABEL <- c(hepatocytes = "Hepatocyte", fibroblasts = "Fibroblast",
              macrophages = "Macrophage", cholangiocytes = "Cholangiocyte", tcells = "T cell")
CT_SHORT <- c(hepatocytes = "Hep", fibroblasts = "Fib", macrophages = "Mac",
              cholangiocytes = "Chol", tcells = "T")
CT_FULL  <- c(hepatocytes = "Hepatocytes", fibroblasts = "Fibroblasts",
              macrophages = "Macrophages", cholangiocytes = "Cholangiocytes", tcells = "T cells")
CT_COL   <- setNames(unname(ct_palette[CT_FULL[CTS]]), CT_LABEL[CTS])
G6 <- function(...) gpar(fontsize = 6, fontfamily = "Helvetica", ...)
MAIN_MODULE_WIDTH_IN  <- 3.54
MAIN_MODULE_HEIGHT_IN <- 2.54

# donor-collapse machinery ---------------------------------------------------
srr_to_donor <- build_srr_to_donor_map(BASE)
collapse_scores <- function(dt) {           # sample,module,score -> donor,module,score (mean over runs)
  dt <- copy(dt)
  dt[, donor := ifelse(sample %in% names(srr_to_donor), srr_to_donor[sample], sample)]
  dt[, .(score = mean(score, na.rm = TRUE)), by = .(module, donor)]
}

# ---------------------------------------------------------------------------
# (1) disease-significant modules (donor-level stats)
# ---------------------------------------------------------------------------
am  <- fread(ALL_MOD)
sig <- am[disease_stage_q < 0.05 & cell_type %in% CTS,
          .(cell_type, module, beta = disease_stage_beta, q = disease_stage_q,
            name = fifelse(is.na(module_name) | module_name == "", best_match_program, module_name))]
NSIG <- table(factor(sig$cell_type, levels = CTS))

# ---------------------------------------------------------------------------
# (2) module x stage matrix — dataset-centered median on DONOR-collapsed scores
# ---------------------------------------------------------------------------
meta <- fread(META_F, select = c("sample", "disease_stage_coarse", "dataset", "exclude_stage_analysis"))
meta <- meta[exclude_stage_analysis != TRUE & disease_stage_coarse %in% STAGES]
meta[, donor := ifelse(sample %in% names(srr_to_donor), srr_to_donor[sample], sample)]
meta_donor <- unique(meta[, .(donor, disease_stage_coarse, dataset)], by = "donor")  # invariant within donor

trend <- rbindlist(lapply(CTS, function(ct) {
  sm <- sig[cell_type == ct, module]; if (!length(sm)) return(NULL)
  ds <- merge(collapse_scores(fread(file.path(HS, ct, "donor_scores.tsv"))[module %in% sm]),
              meta_donor, by = "donor")
  ds[, sc := score - mean(score, na.rm = TRUE), by = .(module, dataset)]
  ds[, .(med = median(sc, na.rm = TRUE)), by = .(module, disease_stage_coarse)][, cell_type := ct][]
}))
trend[, stage := factor(disease_stage_coarse, levels = STAGES)]

use_lmer <- requireNamespace("lmerTest", quietly = TRUE)
if (use_lmer) suppressPackageStartupMessages(library(lmerTest))
slope <- rbindlist(lapply(CTS, function(ct) {
  sm <- sig[cell_type == ct, module]; if (!length(sm)) return(NULL)
  ds <- merge(collapse_scores(fread(file.path(HS, ct, "donor_scores.tsv"))[module %in% sm]),
              meta_donor, by = "donor")
  ds[, ord := as.integer(factor(disease_stage_coarse, levels = STAGES)) - 1L]
  rbindlist(lapply(sm, function(mm) {
    d  <- ds[module == mm]
    co <- tryCatch({
      if (use_lmer && length(unique(d$dataset)) > 1L) {
        summary(lmer(score ~ ord + (1 | dataset), data = d))$coefficients["ord", c("Estimate", "Std. Error")]
      } else summary(lm(score ~ ord, data = d))$coefficients["ord", c("Estimate", "Std. Error")]
    }, error = function(e) c(NA_real_, NA_real_))
    data.table(cell_type = ct, module = mm, beta = co[[1]], SE = co[[2]])
  }))
}))

W <- merge(dcast(trend, cell_type + module ~ stage, value.var = "med"),
           sig[, .(cell_type, module, name, q)], by = c("cell_type", "module"))
W <- merge(W, slope, by = c("cell_type", "module"))
W[, cell_type := factor(cell_type, levels = CTS)]
W[, dir := ifelse(beta > 0, "up", "down")]
W[, rk  := frank(q, ties.method = "first"), by = .(cell_type, dir)]
W[, `:=`(lo = beta - 1.96 * SE, hi = beta + 1.96 * SE)]
setorder(W, cell_type, -beta)
fwrite(W[, .(cell_type, module, name, beta, SE, lo, hi, q, dir, rk, Healthy, Steatosis, Steatohepatitis)],
       file.path(DATA_DIR, "fig3g_module_matrix_DONORLEVEL.csv"))

LIM  <- as.numeric(quantile(abs(as.matrix(W[, ..STAGES])), 0.98, na.rm = TRUE))
COLF <- colorRamp2(c(-LIM, 0, LIM), c("#1565C0", "white", "#C9265E"))

# ---------------------------------------------------------------------------
# module-heatmap renderer (verbatim from singlecell_module_heatmap.R)
# ---------------------------------------------------------------------------
ROW_MM  <- 2.4; GAP_MM  <- 0.6
draw_modules <- function(Wsub, out_pdf, block_titles = TRUE,
                         page_width_in = NULL, page_height_in = NULL,
                         row_mm = ROW_MM, overhead_in = NULL) {
  Wsub <- copy(Wsub); setorder(Wsub, cell_type, -beta)
  mat  <- as.matrix(Wsub[, ..STAGES])
  nm <- sub(" \\([^)]*\\)$", "", Wsub$name)
  rownames(mat) <- sprintf("%s-%d · %s", CT_SHORT[as.character(Wsub$cell_type)], Wsub$module, nm)
  nr      <- nrow(mat)
  present <- CTS[CTS %in% as.character(Wsub$cell_type)]
  nblk    <- length(present)
  row_split <- factor(CT_LABEL[as.character(Wsub$cell_type)], levels = CT_LABEL[present])
  left_anno  <- rowAnnotation(ct = CT_LABEL[as.character(Wsub$cell_type)],
                              col = list(ct = CT_COL), show_legend = FALSE,
                              show_annotation_name = FALSE, simple_anno_size = unit(1.6, "mm"))
  BV <- Wsub$beta; LOv <- Wsub$lo; HIv <- Wsub$hi
  xr <- range(c(LOv, HIv, 0), na.rm = TRUE); xr <- xr + c(-0.04, 0.04) * diff(xr)
  ci_fun <- function(index, k, n) {
    m <- length(index); yy <- m:1
    pushViewport(viewport(xscale = xr, yscale = c(0.5, m + 0.5)))
    grid.segments(unit(0, "native"), unit(0, "npc"), unit(0, "native"), unit(1, "npc"),
                  gp = gpar(col = "grey75", lwd = 0.4))
    for (i in seq_len(m)) {
      j <- index[i]; b <- BV[j]; l <- LOv[j]; h <- HIv[j]; y <- yy[i]
      fill <- if (!is.na(b) && b > 0) "#C9265E" else "#1565C0"
      grid.rect(x = unit(min(0, b), "native"), y = unit(y, "native"),
                width = unit(abs(b), "native"), height = unit(0.55, "native"),
                just = c("left", "centre"), gp = gpar(fill = fill, col = NA))
      if (!is.na(l) && !is.na(h)) {
        grid.segments(unit(l, "native"), unit(y, "native"), unit(h, "native"), unit(y, "native"),
                      gp = gpar(col = "grey20", lwd = 0.5))
        grid.segments(unit(l, "native"), unit(y - 0.18, "native"), unit(l, "native"), unit(y + 0.18, "native"), gp = gpar(col = "grey20", lwd = 0.5))
        grid.segments(unit(h, "native"), unit(y - 0.18, "native"), unit(h, "native"), unit(y + 0.18, "native"), gp = gpar(col = "grey20", lwd = 0.5))
      }
    }
    if (k == n) grid.xaxis(at = pretty(xr, 3), gp = gpar(fontsize = 6, fontfamily = "Helvetica"))
    popViewport()
  }
  right_anno <- rowAnnotation(
    "disease β (95% CI)" = AnnotationFunction(fun = ci_fun, which = "row",
        width = unit(1.5, "cm"), var_import = list(BV = BV, LOv = LOv, HIv = HIv, xr = xr)),
    annotation_name_gp = G6(), annotation_name_rot = 0, annotation_name_side = "bottom",
    annotation_name_offset = unit(4, "mm"))
  ht <- Heatmap(mat, name = "Module score\n(centered)", col = COLF, na_col = "grey92",
    cluster_rows = FALSE, cluster_columns = FALSE, row_split = row_split, row_gap = unit(GAP_MM, "mm"),
    height = unit(row_mm * nr, "mm"), width = unit(1.35, "cm"),
    row_title = if (block_titles) sprintf("%s\n(%d sig)", unname(CT_LABEL[present]), NSIG[present]) else NULL,
    row_title_gp = G6(), row_title_rot = 0, row_names_side = "left", row_names_gp = G6(),
    column_names_gp = G6(), column_names_rot = 45, column_names_side = "bottom",
    column_title = NULL, left_annotation = left_anno, right_annotation = right_anno,
    rect_gp = gpar(col = "white", lwd = 0.3),
    heatmap_legend_param = list(title_gp = G6(), labels_gp = G6(), grid_width = unit(2.5, "mm")))
  if (is.null(overhead_in)) overhead_in <- if (block_titles) 1.05 else 0.56
  ph <- (row_mm * nr + GAP_MM * (nblk - 1)) / 25.4 + overhead_in
  pw <- if (block_titles) 4.4 else 3.7
  if (!is.null(page_width_in))  pw <- page_width_in
  if (!is.null(page_height_in)) ph <- page_height_in
  cairo_pdf(out_pdf, width = pw, height = ph)
  draw(ht, merge_legends = TRUE, heatmap_legend_side = "right",
       annotation_legend_side = "right", padding = unit(c(5, 1, 1, 1), "mm"))
  dev.off()
  cat(sprintf("[saved] %s  (%d modules, %.1f x %.1f in)\n", out_pdf, nr, pw, ph))
}

draw_modules(W,          OUT_SUPP, block_titles = TRUE)
draw_modules(W[rk <= 3], OUT_MAIN, block_titles = FALSE,
             page_width_in = MAIN_MODULE_WIDTH_IN, page_height_in = MAIN_MODULE_HEIGHT_IN,
             row_mm = 2.0, overhead_in = 0.56)

n_sig_total <- nrow(am[disease_stage_q < 0.05 & cell_type %in% CTS])
n_tested    <- nrow(am[!is.na(disease_stage_q) & cell_type %in% CTS])
cat(sprintf("\n[DONOR-LEVEL] %d of %d modules disease-significant (q<0.05) across the 5 Fig3G cell types\n",
            n_sig_total, n_tested))
cat("per cell type (n sig): "); print(NSIG)
