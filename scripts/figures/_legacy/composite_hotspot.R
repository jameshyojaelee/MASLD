#!/usr/bin/env Rscript
# composite_hotspot.R — CANDIDATE composite (② Hotspot modules)
# LEFT : module-score co-variation matrix across donors (pie glyphs), grouped by cell type.
# RIGHT: single-cell disease slope vs bulk logFC per module (z-scored; sc<->bulk concordance).
# Output: FIG2_DIR/panels/composite_hotspot.pdf
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
source(file.path(BASE, "scripts/figures/composite_helpers.R"))
PANEL_DIR <- file.path(FIG2_DIR, "panels"); DATA_DIR <- file.path(PANEL_DIR, "data")
HS <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")

am <- fread(file.path(HS, "all_modules.tsv"))
am[, key := paste0(cell_type, "__", module)]
am <- am[!is.na(disease_stage_beta) & !is.na(mean_bulk_logFC)]
# top disease-progressing modules, up to 6 per cell type
sel <- am[disease_stage_q < 0.05][order(-abs(disease_stage_beta))][, head(.SD, 6), by = cell_type]
sel[, lab := ifelse(grepl("\\(", module_name), sub(".*\\(([^)]*)\\).*", "\\1", module_name), module_name)]
sel[, lab := make.unique(lab, sep = ".")]
disp <- c(hepatocytes = "Hepatocytes", macrophages = "Macrophages", fibroblasts = "Fibroblasts",
          cholangiocytes = "Cholangiocytes", tcells = "T cells")
sel[, ct := disp[cell_type]]
gl <- unname(disp[names(disp) %in% sel$cell_type])
sel[, ct := factor(ct, levels = gl)]; setorder(sel, ct)

# module-score correlation across donors
ds <- fread(file.path(HS, "donor_scores_all.tsv"))[module %in% sel$key]
w <- dcast(ds, sample ~ module, value.var = "score")
wm <- as.matrix(w[, -1]); rownames(wm) <- w$sample
cormat <- cor(wm, use = "pairwise.complete.obs", method = "spearman")
cormat <- cormat[sel$key, sel$key]; dimnames(cormat) <- list(sel$lab, sel$lab)
cormat[is.na(cormat)] <- 0
group_of <- setNames(as.character(sel$ct), sel$lab)

sel[, sc_z := as.numeric(scale(disease_stage_beta))][, bulk_z := as.numeric(scale(mean_bulk_logFC))]
ord <- order_by_group_then_clust(cormat, group_of, gl)
cat(sprintf("[hotspot] %d modules, %d cell types\n", nrow(sel), length(gl)))

mat <- pie_glyph_matrix(cormat, ord, group_of, group_levels = gl, italic_items = TRUE)
# MATRIX-ONLY: the sc-slope-vs-bulk-logFC lollipop was dropped — z-scoring two metrics with a ~6x SD
# mismatch (sc disease_stage_beta vs the tiny mean_bulk_logFC) made the "gap" a scaling artifact, not a
# real single-cell-vs-bulk discordance signal.
p_hs <- mat$plot
save_fig(p_hs, file.path(PANEL_DIR, "composite_hotspot.pdf"), width = 6.0, height = 5.4)
message(sprintf(paste0("CAPTION (Hotspot module co-variation): donor-level Spearman co-variation of %d ",
  "disease Hotspot modules (pie fill proportional to |correlation|, blue +/red -), grouped by cell type; ",
  "modules labeled by marker gene. CAVEAT: donor-level module-score correlations partly reflect shared ",
  "cell-type composition across donors (a donor richer in a cell type scores high on all of its modules), ",
  "not only functional co-regulation."), nrow(sel)))
fwrite(sel[, .(cell_type, module, module_name, lab, disease_stage_beta = round(disease_stage_beta,3),
        mean_bulk_logFC = round(mean_bulk_logFC,3))], file.path(DATA_DIR, "composite_hotspot.csv"))
