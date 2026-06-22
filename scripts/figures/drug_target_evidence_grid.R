#!/usr/bin/env Rscript
# ============================================================================
# drug_target_evidence_grid.R
# ALTERNATIVE to the Fig 5c scatter (drug_target_orthogonality_landscape.R):
# a dense multi-evidence GRID for the MASLD drug pipeline. Same data-driven
# target set (mash_clinical_pipeline.tsv + discovery anchors), but instead of
# 2 axes it packs ~7 evidence dimensions + clinical status/phase/mechanism.
#
# Rows  = drug targets, sorted by colocalization PP.H4 (anchoring gradient:
#         THRB/RORA/HKDC1 top -> un-anchored pipeline bottom).
# Left  annotations: clinical status (approved/active/discontinued/preclinical),
#         mechanism class, max clinical phase.
# Cols  : COLOC PP.H4 | Disease/scRNA-hep/Mouse/Protein log2FC | essentiality |
#         cell-type breadth.
#
# INTEGRITY: PP.H4 read ONLY from gene_level_coloc.csv (<=1). NOT a hand list.
# Output: figures/main/fig5_convergence/panels/fig5c_drug_target_grid.pdf
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})
set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
source(file.path(BASE, "scripts/publication_color_themes.R"))

PANEL_DIR <- file.path(FIG5_DIR, "panels")
OUT_PDF   <- file.path(PANEL_DIR, "fig5c_drug_target_grid.pdf")

# ---------------------------------------------------------------------------
# Target set (identical sourcing to the scatter panel)
# ---------------------------------------------------------------------------
pipe <- fread(file.path(BASE, "data/external/druggability/mash_clinical_pipeline.tsv"),
              fill = TRUE)
pipe[, mech_class := fcase(
  grepl("THR-beta", mechanism),                                            "Thyroid (THR-β)",
  grepl("FGF21|GLP-1|GCGR|GIP|klotho", mechanism, ignore.case = TRUE),     "Incretin / FGF21",
  grepl("PPAR|pyruvate carrier",       mechanism, ignore.case = TRUE),     "PPAR",
  grepl("DGAT2|SCD1|ACC inhibitor|FASN|ketohexokinase", mechanism, ignore.case = TRUE), "DNL / lipogenesis",
  grepl("FXR",                         mechanism),                         "FXR",
  grepl("ROR",                         mechanism),                         "RORα",
  grepl("CCR2|CCR5|ASK1|LOXL2",        mechanism, ignore.case = TRUE),     "Fibro-inflammatory",
  default =                                                                "Other metabolic / genetic"
)]
status_rank <- c(approved = 3L, active = 2L, discontinued = 1L)
pipe[, status_r := status_rank[status]]
setorder(pipe, -status_r, -max_phase)
targets <- pipe[, .(status = status[1], max_phase = max(max_phase, na.rm = TRUE),
                    mech_class = mech_class[1]), by = .(human_symbol = gene_symbol)]
disc <- data.table(human_symbol = c("HKDC1", "AKR1B10"),
                   status = "preclinical", max_phase = NA_real_,
                   mech_class = "Preclinical (validated)")
targets <- rbind(targets, disc[!human_symbol %in% targets$human_symbol])

# ---------------------------------------------------------------------------
# Evidence matrix
# ---------------------------------------------------------------------------
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
ev <- atlas[!duplicated(human_symbol),
            .(human_symbol, bulk_shrunk_logFC, sc_hepatocyte_logFC, mouse_meta_logFC,
              best_protein_logFC, sc_n_celltypes_sig, essentiality_chronos)]
coloc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
stopifnot(max(coloc$coloc_best_pp4, na.rm = TRUE) <= 1)
cu <- coloc[gene != "" & !is.na(gene) & !duplicated(gene), .(human_symbol = gene, coloc_best_pp4)]

d <- merge(targets, ev, by = "human_symbol", all.x = TRUE)
d <- merge(d, cu, by = "human_symbol", all.x = TRUE)
d[is.na(coloc_best_pp4), coloc_best_pp4 := 0]          # not a tested eGene / no signal
miss_x <- d[is.na(bulk_shrunk_logFC), human_symbol]    # GLP1R: not hepatically expressed
if (length(miss_x))
  cat(sprintf("[warn] dropped (no hepatic expression): %s\n", paste(miss_x, collapse = ", ")))
d <- d[!is.na(bulk_shrunk_logFC)]
setorder(d, -coloc_best_pp4)
genes <- d$human_symbol
cat(sprintf("[grid] %d targets x 7 evidence dims\n", nrow(d)))

clamp <- function(x, lo, hi) { y <- pmin(pmax(x, lo), hi); y[is.na(x)] <- NA; y }

mat_lfc <- clamp(as.matrix(d[, .(Disease = bulk_shrunk_logFC, `scRNA hep` = sc_hepatocyte_logFC,
                                 Mouse = mouse_meta_logFC, Protein = best_protein_logFC)]), -2, 2)
rownames(mat_lfc) <- genes
mat_pp4 <- matrix(d$coloc_best_pp4,       ncol = 1, dimnames = list(genes, "COLOC\nPP.H4"))
mat_ess <- matrix(d$essentiality_chronos, ncol = 1, dimnames = list(genes, "Chronos\n(DepMap)"))

# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------
col_lfc     <- colorRamp2(c(-2, 0, 2), c(sea_gradient[9], "white", sunset_gradient[7]))
col_pp4     <- colorRamp2(c(0, 0.5, 1), c("white", sunset_gradient[4], sunset_gradient[11]))
col_ess     <- colorRamp2(c(-1, -0.5, 0), c(purple_gradient[10], purple_gradient[5], "white"))  # more essential = darker

status_cols <- c(approved = "#2E7D32", active = "#90A4AE",
                 discontinued = "#C62828", preclinical = "#5E35B1")
mech_cols <- c(
  `Thyroid (THR-β)`           = "#00695C",
  `Incretin / FGF21`          = "#F39C12",
  `PPAR`                      = "#1565C0",
  `DNL / lipogenesis`         = "#C9265E",
  `FXR`                       = "#7B1FA2",
  `RORα`                      = "#AD1457",
  `Fibro-inflammatory`        = "#795548",
  `Other metabolic / genetic` = "#9E9E9E",
  `Preclinical (validated)`   = "#212121"
)
col_phase <- colorRamp2(c(1, 4), c("#E1BEE7", "#4A148C"))

ht_opt(
  heatmap_row_names_gp    = gpar(fontsize = 6.5, fontface = "italic"),
  heatmap_column_names_gp = gpar(fontsize = 6),
  heatmap_column_title_gp = gpar(fontsize = 7, fontface = "bold"),
  legend_title_gp         = gpar(fontsize = 6.5, fontface = "bold"),
  legend_labels_gp        = gpar(fontsize = 6))

# ---------------------------------------------------------------------------
# Left annotation: clinical status + mechanism + phase
# ---------------------------------------------------------------------------
ha <- rowAnnotation(
  Status    = d$status,
  Mechanism = as.character(d$mech_class),
  Phase     = d$max_phase,
  col = list(Status = status_cols, Mechanism = mech_cols, Phase = col_phase),
  na_col = "grey92",
  annotation_name_gp   = gpar(fontsize = 6),
  annotation_name_rot  = 45,
  simple_anno_size     = unit(3.2, "mm"),
  annotation_legend_param = list(
    Status    = list(title = "Clinical status", nrow = 4,
                     at = c("approved","active","discontinued","preclinical")),
    Mechanism = list(title = "Mechanism", ncol = 1),
    Phase     = list(title = "Max phase", at = c(1,2,3,4))))

h_pp4 <- Heatmap(mat_pp4, name = "PP.H4", col = col_pp4, column_title = "Genetic",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  left_annotation = ha, na_col = "grey92", border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5), width = unit(5, "mm"),
  cell_fun = function(j, i, x, y, w, h, fill)
    grid.text(sprintf("%.2f", mat_pp4[i, j]), x, y, gp = gpar(fontsize = 4.4,
              col = ifelse(mat_pp4[i, j] > 0.6, "white", "grey25"))),
  heatmap_legend_param = list(at = c(0, 0.5, 1)))

h_lfc <- Heatmap(mat_lfc, name = "log2FC", col = col_lfc, column_title = "Fold change",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  na_col = "grey92", border = TRUE, rect_gp = gpar(col = "white", lwd = 0.5),
  width = unit(5 * ncol(mat_lfc), "mm"),
  heatmap_legend_param = list(at = c(-2, 0, 2)))

h_ess <- Heatmap(mat_ess, name = "Chronos", col = col_ess, column_title = "Essential",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE,
  row_names_side = "right", na_col = "grey92", border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5), width = unit(5, "mm"),
  heatmap_legend_param = list(at = c(-1, -0.5, 0)))

ht <- h_pp4 + h_lfc + h_ess

cat(sprintf("Saving %s ...\n", OUT_PDF))
cairo_pdf(OUT_PDF, width = 5.4, height = 6.0, family = "Helvetica")
draw(ht,
  column_title = "c   Drug-target evidence grid (sorted by genetic anchoring)",
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  ht_gap = unit(1.6, "mm"),
  merge_legends = TRUE, heatmap_legend_side = "right",
  annotation_legend_side = "right", padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()
ht_opt(RESET = TRUE)
cat("Done.\n")
