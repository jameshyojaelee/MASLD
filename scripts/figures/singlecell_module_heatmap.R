#!/usr/bin/env Rscript
# KEY MESSAGE: The whole single-cell atlas remodels across MASLD — disease-
# significant Hotspot co-expression modules across 5 cell types shift coordinately
# Healthy → Steatosis → Steatohepatitis. Three outputs from one matrix (shared
# colour scale so they never drift):
#   MAIN(3I) fig3g_singlecell_module_heatmap.pdf — COMPACT: top-3 up + top-3 down
#         per cell type (17 modules) — the multi-cellular hand-off at a glance.
#   MAIN(3J) fig3h_ccc_lr_heatmap.pdf           — companion communication heatmap
#         (13 headline L-R pairs, row z-scored) as its own panel.
#   SUPP  figs3_singlecell_module_heatmap_full.pdf — the FULL 22-module breadth.
# ============================================================================
# singlecell_module_heatmap.R  — Fig 3G module remodeling heatmap(s)
#
# RIGOR: column axis = the VETTED coarse `disease_stage_coarse` (Healthy /
# Steatosis / Steatohepatitis; Cirrhosis excluded — protocol contamination). NOT
# crn_transition_scores.tsv (F0-F4) — that is the cross-cohort inferred F-stage the
# stat-audit flagged as circular. Values via the dataset-centered-median method
# used by hotspot_cascade.R / the NPC track. Cell-type colours = shared ct_palette
# (identical to the fig3g CCC chord).
#
# DONOR-LEVEL (canonical 2026-07-12): pseudoreplication fix. Module SELECTION reads
# the donor-collapsed disease_stage_q from donor_collapse/all_modules_donor.tsv, and
# the in-script disease-slope beta collapses each per-cell-type donor_scores.tsv from
# sequencing-run to TRUE biological donor (mean over runs, lib_donor_collapse) BEFORE
# lmer. 22 disease-significant modules across the 5 cell types (run-level was 54);
# compact main shows 17 (top-3 up + top-3 down per cell type). The fig3h L-R panel
# reads the donor-collapsed ccc_trajectories_data.csv (repointed in ccc_v3_panels.R).
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ComplexHeatmap); library(circlize); library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))   # ct_palette (shared with the CCC chord)
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
source(file.path(BASE, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))  # build_srr_to_donor_map

HS        <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
ALL_MOD   <- file.path(HS, "donor_collapse/all_modules_donor.tsv")   # DONOR-LEVEL selection stats (pseudoreplication fix)
META_F    <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv")
PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
OUT_MAIN  <- file.path(PANEL_DIR, "fig3g_singlecell_module_heatmap.pdf")            # compact (main; 3I in the A-J layout 2026-07-02)
OUT_SUPP  <- file.path(PANEL_DIR, "figs3_singlecell_module_heatmap_full.pdf")       # full (supp)
OUT_LIANA <- file.path(PANEL_DIR, "fig3h_ccc_lr_heatmap.pdf")                       # LR heatmap (main; 3J in the A-J layout 2026-07-02)
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

STAGES   <- c("Healthy", "Steatosis", "Steatohepatitis")
CTS      <- c("hepatocytes", "fibroblasts", "macrophages", "cholangiocytes", "tcells")
CT_LABEL <- c(hepatocytes = "Hepatocyte", fibroblasts = "Fibroblast",
              macrophages = "Macrophage", cholangiocytes = "Cholangiocyte", tcells = "T cell")
CT_SHORT <- c(hepatocytes = "Hep", fibroblasts = "Fib", macrophages = "Mac",
              cholangiocytes = "Chol", tcells = "T")
# Cell-type strip colours from ct_palette (shared with the supp CCC chord figs3_ccc_chord_all_lr).
CT_FULL  <- c(hepatocytes = "Hepatocytes", fibroblasts = "Fibroblasts",
              macrophages = "Macrophages", cholangiocytes = "Cholangiocytes", tcells = "T cells")
CT_COL   <- setNames(unname(ct_palette[CT_FULL[CTS]]), CT_LABEL[CTS])
G6 <- function(...) gpar(fontsize = 6, fontfamily = "Helvetica", ...)
MAIN_MODULE_WIDTH_IN  <- 3.54
MAIN_MODULE_HEIGHT_IN <- 2.54
MAIN_LIANA_WIDTH_IN   <- 3.54
MAIN_LIANA_HEIGHT_IN  <- 1.95

# donor-collapse machinery (pseudoreplication fix) --------------------------
# atlas `sample` = a sequencing run / sort-fraction library; collapse to true
# biological donor (mean over runs) so per-donor stats are not pseudoreplicated.
srr_to_donor <- build_srr_to_donor_map(BASE)
collapse_scores <- function(dt) {           # sample,module,score -> donor,module,score (mean over runs)
  dt <- copy(dt)
  dt[, donor := ifelse(sample %in% names(srr_to_donor), srr_to_donor[sample], sample)]
  dt[, .(score = mean(score, na.rm = TRUE)), by = .(module, donor)]
}

# ---------------------------------------------------------------------------
# (1) disease-significant modules (+ beta + biological name)  [DONOR-LEVEL q]
# ---------------------------------------------------------------------------
am  <- fread(ALL_MOD)
sig <- am[disease_stage_q < 0.05 & cell_type %in% CTS,
          .(cell_type, module, beta = disease_stage_beta, q = disease_stage_q,
            name = fifelse(is.na(module_name) | module_name == "", best_match_program, module_name))]
# name is the canonical module_name (merged from module_names.tsv into all_modules.tsv;
# verified 0 mismatches) with best_match_program only as a fallback — no hardcoded dict.
# DEFECT 3: count of disease-significant (q<0.05) modules per cell type — annotated onto
# each block title so light blocks (e.g. T cells, 1 sig) read as complete, not truncated.
NSIG <- table(factor(sig$cell_type, levels = CTS))

# ---------------------------------------------------------------------------
# (2) module x stage matrix — dataset-centered median (rigorous coarse axis)
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
# DEFECT 4: disease-slope beta + SE recomputed IN-SCRIPT on the SAME 3 displayed stages
# (Healthy/Steatosis/Steatohepatitis; Cirrhosis EXCLUDED) via
# module_score ~ stage_ordinal + (1|dataset)  [lmerTest; lm fallback if 1 dataset].
# Previously the bar beta/SE were read from phenotype_correlations.tsv, whose disease_stage
# axis is a 4-level ordinal INCLUDING cirrhosis — so the bars and the (cirrhosis-free) colour
# columns disagreed. `meta` is already filtered to the 3 stages above, so it is reused here.
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
W <- merge(W, slope, by = c("cell_type", "module"))          # 3-stage beta + SE drive the bars
W[, cell_type := factor(cell_type, levels = CTS)]
W[, dir := ifelse(beta > 0, "up", "down")]
# DEFECT 2: rank within cell type × direction by disease_stage_q (most significant first) —
# surfaces the most robustly stage-associated modules (was ranked by |beta|). Selection stays
# top-3 up + top-3 down among disease_stage_q<0.05; bars keep the gain->loss signed-beta sort.
W[, rk  := frank(q, ties.method = "first"), by = .(cell_type, dir)]
W[, `:=`(lo = beta - 1.96 * SE, hi = beta + 1.96 * SE)]
setorder(W, cell_type, -beta)
fwrite(W[, .(cell_type, module, name, beta, SE, lo, hi, q, dir, rk, Healthy, Steatosis, Steatohepatitis)],
       file.path(DATA_DIR, "fig3g_module_matrix.csv"))

# SHARED colour limit (full matrix) so compact + full panels share the exact scale
LIM  <- as.numeric(quantile(abs(as.matrix(W[, ..STAGES])), 0.98, na.rm = TRUE))
COLF <- colorRamp2(c(-LIM, 0, LIM), c("#1565C0", "white", "#C9265E"))

# ---------------------------------------------------------------------------
# module-heatmap renderer (reused for full + compact)
# ---------------------------------------------------------------------------
ROW_MM  <- 2.4                                    # tight row height (snug for 6 pt)
GAP_MM  <- 0.6
# block_titles: full/supp panel shows cell-type block titles (self-contained);
# the compact MAIN drops them (relies on the colored strip + row-label prefixes +
# the SHARED cell-type legend in the fig3g CCC chord, placed side-by-side).
draw_modules <- function(Wsub, out_pdf, block_titles = TRUE,
                         page_width_in = NULL, page_height_in = NULL,
                         row_mm = ROW_MM, overhead_in = NULL) {
  Wsub <- copy(Wsub); setorder(Wsub, cell_type, -beta)
  mat  <- as.matrix(Wsub[, ..STAGES])
  # Row labels: hyphenate the ID prefix (Hep-17) and drop ONLY a trailing gene
  # parenthetical from the name (the full name + gene stay in the cascade/geneset panels).
  nm <- sub(" \\([^)]*\\)$", "", Wsub$name)
  rownames(mat) <- sprintf("%s-%d · %s", CT_SHORT[as.character(Wsub$cell_type)], Wsub$module, nm)
  nr      <- nrow(mat)
  present <- CTS[CTS %in% as.character(Wsub$cell_type)]
  nblk    <- length(present)
  row_split <- factor(CT_LABEL[as.character(Wsub$cell_type)], levels = CT_LABEL[present])
  left_anno  <- rowAnnotation(ct = CT_LABEL[as.character(Wsub$cell_type)],
                              col = list(ct = CT_COL), show_legend = FALSE,
                              show_annotation_name = FALSE, simple_anno_size = unit(1.6, "mm"))
  # β bars WITH 95% CI whiskers (β ± 1.96·SE from the dataset-adjusted lmer). Custom
  # annotation: anno_barplot can't draw error bars, so draw bar + whisker per row in grid.
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
    # title UNDER the numeric β axis (both x-axis titles on the bottom row, aligned with the
    # heatmap's bottom "Disease stage (coarse)"). The grid.xaxis ticks are drawn at the bottom
    # of the last slice (k == n); annotation_name_offset pushes the title clear of those tick
    # labels so they don't overlap.
    annotation_name_gp = G6(), annotation_name_rot = 0, annotation_name_side = "bottom",
    annotation_name_offset = unit(4, "mm"))
  ht <- Heatmap(mat, name = "Module score\n(centered)", col = COLF, na_col = "grey92",
    cluster_rows = FALSE, cluster_columns = FALSE, row_split = row_split, row_gap = unit(GAP_MM, "mm"),
    height = unit(row_mm * nr, "mm"), width = unit(1.35, "cm"),        # tight body — no row stretch
    row_title = if (block_titles) sprintf("%s\n(%d sig)", unname(CT_LABEL[present]), NSIG[present]) else NULL,  # DEFECT 3: append sig-module count (NULL -> no block titles)
    row_title_gp = G6(), row_title_rot = 0, row_names_side = "left", row_names_gp = G6(),
    column_names_gp = G6(), column_names_rot = 45, column_names_side = "bottom",  # stage labels moved to bottom
    column_title = NULL,                                                          # "Disease stage (coarse)" dropped — the 3 stage names are self-explanatory
    left_annotation = left_anno, right_annotation = right_anno,
    rect_gp = gpar(col = "white", lwd = 0.3),
    heatmap_legend_param = list(title_gp = G6(), labels_gp = G6(),
                                grid_width = unit(2.5, "mm")))
  if (is.null(overhead_in)) overhead_in <- if (block_titles) 1.05 else 0.56
  ph <- (row_mm * nr + GAP_MM * (nblk - 1)) / 25.4 + overhead_in      # page height = body + overhead (col labels + barplot axis + bottom beta-axis title)
  pw <- if (block_titles) 4.4 else 3.7                               # drop the block-title column when compact
  if (!is.null(page_width_in))  pw <- page_width_in
  if (!is.null(page_height_in)) ph <- page_height_in
  cairo_pdf(out_pdf, width = pw, height = ph)
  draw(ht, merge_legends = TRUE, heatmap_legend_side = "right",
       annotation_legend_side = "right", padding = unit(c(5, 1, 1, 1), "mm"))  # extra bottom pad for the moved β-axis title
  dev.off()
  cat(sprintf("[saved] %s  (%d modules, %.1f x %.1f in)\n", out_pdf, nr, pw, ph))
}

draw_modules(W,          OUT_SUPP, block_titles = TRUE)    # FULL 54 -> supp (block titles + (n sig) help navigate)
draw_modules(W[rk <= 3], OUT_MAIN, block_titles = FALSE,
             page_width_in = MAIN_MODULE_WIDTH_IN,
             page_height_in = MAIN_MODULE_HEIGHT_IN,
             row_mm = 2.0,
             overhead_in = 0.56)  # compact 24 -> main: clean left margin; fixed to Fig 3 size contract

# ---------------------------------------------------------------------------
# (3) LIANA communication heatmap — standalone MAIN panel
# ---------------------------------------------------------------------------
lia <- fread(file.path(DATA_DIR, "ccc_trajectories_data.csv"))[disease_stage_coarse %in% STAGES]
LW  <- dcast(lia, headline_label ~ factor(disease_stage_coarse, levels = STAGES), value.var = "mean_score")
# Tidy L-R labels: abbreviate integrin heterodimer receptors (ITGA1_ITGB1 -> ITGA1/B1),
# collapse any other X_Y complex to X/Y, and squeeze the double space before the
# (Sender->Receiver) tag to a single space. The (Sender->Receiver) tag is kept — it
# disambiguates repeated L-R pairs (e.g. C4BPA->BMPR2 for both Hep->Fib and Hep->Endo).
lr_lab <- LW$headline_label
lr_lab <- gsub("ITGA(\\w+)_ITGB(\\w+)", "ITGA\\1/B\\2", lr_lab)
lr_lab <- gsub("_", "/", lr_lab)
lr_lab <- gsub(" +\\(", " (", lr_lab)
lmat <- t(scale(t(as.matrix(LW[, ..STAGES])))); rownames(lmat) <- lr_lab   # row z-score
fwrite(data.table(headline_label = lr_lab, lmat), file.path(DATA_DIR, "fig3h_liana_matrix.csv"))
# order by net change (Steatohepatitis − Healthy) so strengthen→weaken reads top→bottom
lmat <- lmat[order(-(lmat[, "Steatohepatitis"] - lmat[, "Healthy"])), , drop = FALSE]
lmat <- t(lmat)                                     # transpose -> 3 stages (rows) x 13 L-R pairs (cols) = WIDE
ZCOL <- colorRamp2(c(-1.5, 0, 1.5), c("#1565C0", "white", "#C9265E"))
ht_l <- Heatmap(lmat, name = "L-R z", col = ZCOL, na_col = "grey92",
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_title = NULL,
  row_names_side = "left", row_names_gp = G6(),
  column_names_gp = G6(), column_names_rot = 45, column_names_side = "bottom",
  column_title = NULL,
  height = unit(3 * 0.50, "cm"), width = unit(13 * 0.42, "cm"),
  rect_gp = gpar(col = "white", lwd = 0.5),
  heatmap_legend_param = list(title_gp = G6(), labels_gp = G6(), grid_width = unit(2.5, "mm")))
cairo_pdf(OUT_LIANA, width = MAIN_LIANA_WIDTH_IN, height = MAIN_LIANA_HEIGHT_IN)
draw(ht_l, heatmap_legend_side = "right", padding = unit(c(2, 1, 2, 1), "mm"))
dev.off()
cat(sprintf("[saved] %s  (%d L-R pairs)\n", OUT_LIANA, ncol(lmat)))

n_sig_total <- nrow(sig)          # donor-level disease-significant modules across the 5 cell types
n_compact   <- nrow(W[rk <= 3])   # compact main selection (top-3 up + top-3 down per cell type)
message(sprintf(paste0(
  "[Fig 3G/3H caption] Single-cell remodeling across MASLD (DONOR-LEVEL). fig3g (MAIN): per cell type the ",
  "3 disease-up + 3 disease-down Hotspot modules with the SMALLEST donor-level disease-stage FDR q ",
  "(ranked by disease_stage_q, most significant first; ", n_compact, " compact of ", n_sig_total, " donor-significant; full set in ",
  "figs3_..._full) split by cell type × coarse stage Healthy→Steatohepatitis; colour = ",
  "dataset-centered median module score over the 3 stages Healthy/Steatosis/Steatohepatitis ",
  "(Cirrhosis excluded; scale shared with the full supp panel). Right bars = disease slope ",
  "beta +/- 95%% CI computed IN-SCRIPT on the SAME 3 stages (beta +/- 1.96.SE from score ~ ",
  "stage_ordinal + (1|dataset), lmerTest; all shown modules FDR q<0.05); block titles annotate ",
  "the count of disease-significant modules per cell type. fig3h (MAIN): 13 headline ligand-receptor ",
  "pairs (row z-scored) on the same stage axis. Cell-type colours from ct_palette (shared with the supp CCC chord)."
)))
