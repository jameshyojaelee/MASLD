#!/usr/bin/env Rscript
# KEY MESSAGE: Cross-modality disease master-regulator TFs vs F0 baseline across F1-F4 stages; HNF4A, RORA, THRB anchor the largest motif-disruption load.
# NOTE: "motif-disruption load" is a PUTATIVE computational prediction (FIMO/
# motifbreakR on fine-mapped variants), NOT functionally validated -- the Currin
# caQTL test does not confirm it (motif x caQTL 46.3%, median 50%; only 2 four-way
# variants, in ZNF701; see figS05_scatac_caqtl_concordance). Report as PREDICTED
# regulatory disruption, never "validated".
# Mirrors figS05_gwas_atac_tf_transition_heatmap.R but uses F0 as the common reference instead of the prior stage.
# Source: fibrosis_stage_dream.csv (dream mega F{1-4}_vs_F0) + DoRothEA A/B/C + decoupleR run_wmean.
# TF panel: cross-modality disease master regulators = hepatocyte SCENIC+ regulon TFs that are
#   bulk MASLD DEGs (n=846) OR COLOC hits. The prior FDR-gated SCENIC+ disease_regulons.csv is
#   empty (the donor-level n=18 regulon-activity DE honestly finds 0; the old non-zero set was
#   cell-level pseudoreplication), so single-cohort regulon-DE is underpowered and not used to
#   define the panel here.

suppressPackageStartupMessages({
  library(data.table)
  library(decoupleR)
  library(dorothea)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# Cross-modality disease master-regulator set (default); env-overridable.
# Old FDR-gated disease_regulons.csv is empty (pseudoreplication-prone donor-level n=18 DE).
REGULON_F <- Sys.getenv("REGULON_FILE",
                        unset = file.path(ATAC_DIR, "scenic_plus", "disease_master_regulators.csv"))
DREAM_F   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/fibrosis_stage_dream.csv")
SYMBOL_F  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
MOTIF_F   <- file.path(BASE, "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv")
ACT_CSV   <- file.path(BASE, "RNA-seq/results/stratified_causal/fstage_vs_F0_tf_activity.csv")

stages <- c("F1_vs_F0", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0")

if (file.exists(ACT_CSV)) {
  message("Reusing cached activity table: ", ACT_CSV)
  tf_all <- fread(ACT_CSV)
} else {
  message("Computing TF activity vs F0 baseline...")
  dream_map <- fread(SYMBOL_F, select = c("gene", "symbol"))
  dream_map <- unique(dream_map[!is.na(symbol) & symbol != ""])

  stage_dt <- fread(DREAM_F)
  stage_dt <- merge(stage_dt, dream_map, by = "gene", all.x = TRUE)
  stage_dt <- stage_dt[!is.na(symbol) & symbol != ""]
  stage_dt <- stage_dt[order(-abs(t))][!duplicated(paste0(symbol, "_", contrast))]

  regulons <- dorothea_hs
  regulons <- regulons[regulons$confidence %in% c("A", "B", "C"), ]
  regulons <- as.data.frame(regulons[, c("tf", "target", "mor")])

  tf_results_list <- list()
  for (st in stages) {
    sub <- stage_dt[contrast == st]
    tstat_vec <- setNames(sub$t, sub$symbol)
    mat <- matrix(tstat_vec, ncol = 1, dimnames = list(names(tstat_vec), st))
    res <- run_wmean(mat = mat, net = regulons, .source = "tf", .target = "target",
                     .mor = "mor", times = 1000, minsize = 5)
    res_wmean <- res[res$statistic == "norm_wmean", ]
    res_wmean$contrast <- st
    tf_results_list[[st]] <- as.data.table(res_wmean)
  }
  tf_all <- rbindlist(tf_results_list)
  setnames(tf_all, "source", "tf")
  fwrite(tf_all[, .(tf, contrast, score, p_value, statistic)], ACT_CSV)
  message("Saved: ", ACT_CSV)
}

regulons_panel <- fread(REGULON_F)
# Master-regulator file keys TFs by `tf_name`; legacy regulon files used `regulon_id`.
if ("tf_name" %in% names(regulons_panel)) {
  tfs <- unique(regulons_panel$tf_name)
} else {
  tfs <- unique(sub("_regulon$", "", regulons_panel$regulon_id))
}
motif <- fread(MOTIF_F)

act <- tf_all[tf %in% tfs]
act[, fdr := p.adjust(p_value, method = "BH")]

score_mat <- matrix(NA_real_, nrow = length(tfs), ncol = length(stages),
                    dimnames = list(tfs, stages))
fdr_mat <- matrix(NA_real_, nrow = length(tfs), ncol = length(stages),
                  dimnames = list(tfs, stages))
for (i in seq_along(tfs)) {
  for (j in seq_along(stages)) {
    row <- act[tf == tfs[i] & contrast == stages[j]]
    if (nrow(row) == 1) {
      score_mat[i, j] <- row$score
      fdr_mat[i, j]   <- row$fdr
    }
  }
}

motif_total   <- motif[tf_name %in% tfs, .N, by = tf_name]
motif_regulon <- motif[tf_name %in% tfs & motif_in_disease_regulon == TRUE, .N, by = tf_name]
disrupt <- merge(data.table(tf_name = tfs), motif_total, by = "tf_name", all.x = TRUE)
disrupt <- merge(disrupt, motif_regulon, by = "tf_name", all.x = TRUE, suffixes = c("_total", "_regulon"))
setnames(disrupt, c("N_total", "N_regulon"), c("n_total", "n_regulon"), skip_absent = TRUE)
disrupt[is.na(n_total),   n_total   := 0L]
disrupt[is.na(n_regulon), n_regulon := 0L]
setkey(disrupt, tf_name)

ord <- order(-disrupt[tfs, n_regulon], -disrupt[tfs, n_total])
tfs_ord     <- tfs[ord]
score_ord   <- score_mat[ord, , drop = FALSE]
fdr_ord     <- fdr_mat[ord, , drop = FALSE]
disrupt_ord <- disrupt[tfs_ord]

sig_fun <- function(j, i, x, y, w, h, fill) {
  f <- fdr_ord[i, j]
  if (is.na(f)) return(invisible(NULL))
  lab <- if (f < 0.01) "**" else if (f < 0.10) "*" else ""
  if (lab != "")
    grid.text(lab, x, y, gp = gpar(fontsize = 6, fontface = "plain", col = "black"))
}

vmax <- max(abs(score_ord), na.rm = TRUE)
col_fun <- colorRamp2(c(-vmax, 0, vmax),
                      c("#2166AC", "#F7F7F7", "#B2182B"))

right_anno <- rowAnnotation(
  `Motif\ndisruptions` = anno_barplot(
    disrupt_ord[, n_total],
    gp = gpar(fill = "#9E9E9E", col = NA),
    bar_width = 0.75,
    width = unit(2.0, "cm")
  ),
  annotation_name_gp = gpar(fontsize = 6),
  annotation_name_rot = 0
)

col_labels <- c("F1", "F2", "F3", "F4")

ht <- Heatmap(
  score_ord,
  name = "TF activity\n(decoupleR z)",
  col = col_fun,
  na_col = "#F0F0F0",
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left", row_names_gp = gpar(fontsize = 6),
  column_labels = col_labels, column_names_rot = 0,
  column_names_gp = gpar(fontsize = 6),
  column_names_centered = TRUE,
  cell_fun = sig_fun,
  right_annotation = right_anno,
  rect_gp = gpar(col = "white", lwd = 0.6),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 6),
    legend_height = unit(2.5, "cm")
  ),
  width = unit(5.5, "cm"),
  height = unit(11, "cm")
)

out_pdf <- file.path(FIGS05_DIR, "figS05_gwas_atac_tf_fstage_heatmap.pdf")
sig_legend <- Legend(
  title      = "Significance",
  labels     = c("* FDR < 0.10", "** FDR < 0.01"),
  type       = "points",
  pch        = NA,
  legend_gp  = gpar(col = "white"),
  labels_gp  = gpar(fontsize = 6),
  title_gp   = gpar(fontsize = 6, fontface = "plain")
)

message("[caption] Cross-modality disease master regulators across fibrosis stages (vs F0)")
pdf(out_pdf, width = 6.0, height = 8.0)
draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_list = list(sig_legend),
     padding = unit(c(4, 4, 4, 4), "mm"))
dev.off()

out_csv <- file.path(FIGS05_DIR, "figS05_gwas_atac_tf_fstage_data.csv")
out_dt <- as.data.table(score_ord, keep.rownames = "tf")
out_dt_long <- melt(out_dt, id.vars = "tf", variable.name = "contrast", value.name = "score")
fdr_dt <- as.data.table(fdr_ord, keep.rownames = "tf")
fdr_long <- melt(fdr_dt, id.vars = "tf", variable.name = "contrast", value.name = "fdr")
out_long <- merge(out_dt_long, fdr_long, by = c("tf", "contrast"))
out_long <- merge(out_long, disrupt_ord, by.x = "tf", by.y = "tf_name", all.x = TRUE)
fwrite(out_long, out_csv)

message("Wrote: ", out_pdf)
message("Wrote: ", out_csv)
