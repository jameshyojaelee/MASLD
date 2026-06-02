# KEY MESSAGE: 24 GWAS-ATAC disease regulon TFs show transition-specific activity shifts across the 4 fibrosis transitions; HNF4A, RORA, THRB anchor the largest motif-disruption load.
# Output: figS05_epigenomic_spatial/ (sits alongside existing figS_gwas_atac_regulons*).
# Data gaps: 5/24 TFs (THRB, RORA, KLF15, MLXIPL, NR1H4, RORC) have no rows in transition_tf_activity.csv -- rendered as NA cells.

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

REGULON_F  <- file.path(ATAC_DIR, "scenic_plus", "disease_regulons.csv")
TRANS_F    <- file.path(BASE, "RNA-seq/results/stratified_causal/transition_tf_activity.csv")
MOTIF_F    <- file.path(BASE, "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv")

regulons <- fread(REGULON_F)
tfs      <- sub("_regulon$", "", regulons$regulon_id)
trans    <- fread(TRANS_F)
motif    <- fread(MOTIF_F)

trans <- trans[tf %in% tfs]
trans[, fdr := p.adjust(p_value, method = "BH")]

transition_levels <- c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4")
score_mat <- matrix(NA_real_, nrow = length(tfs), ncol = length(transition_levels),
                    dimnames = list(tfs, transition_levels))
fdr_mat <- matrix(NA_real_, nrow = length(tfs), ncol = length(transition_levels),
                  dimnames = list(tfs, transition_levels))
for (i in seq_along(tfs)) {
  for (j in seq_along(transition_levels)) {
    row <- trans[tf == tfs[i] & transition == transition_levels[j]]
    if (nrow(row) == 1) {
      score_mat[i, j] <- row$score
      fdr_mat[i, j]   <- row$fdr
    }
  }
}

motif_total <- motif[tf_name %in% tfs, .N, by = tf_name]
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


vmax <- max(abs(score_ord), na.rm = TRUE)
col_fun <- colorRamp2(c(-vmax, 0, vmax),
                      c("#2166AC", "#F7F7F7", "#B2182B"))

right_anno <- rowAnnotation(
  `Motif\ndisruptions` = anno_barplot(
    disrupt_ord[, .(n_regulon, n_total - n_regulon)],
    gp = gpar(fill = c("#9E9E9E", "#D9D9D9"), col = NA),
    bar_width = 0.75,
    width = unit(2.0, "cm")
  ),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_name_rot = 0
)

col_labels <- c("F0->F1", "F1->F2", "F2->F3", "F3->F4")

ht <- Heatmap(
  score_ord,
  name = "TF activity\n(decoupleR z)",
  col = col_fun,
  na_col = "#F0F0F0",
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left", row_names_gp = gpar(fontsize = 9),
  column_labels = col_labels, column_names_rot = 0,
  column_names_gp = gpar(fontsize = 10),
  column_names_centered = TRUE,
  right_annotation = right_anno,
  rect_gp = gpar(col = "white", lwd = 0.6),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8),
    legend_height = unit(2.5, "cm")
  ),
  width = unit(5.5, "cm"),
  height = unit(11, "cm")
)

out_pdf <- file.path(FIGS05_DIR, "figS05_gwas_atac_tf_transition_heatmap.pdf")
pdf(out_pdf, width = 6.0, height = 8.0)
draw(ht,
     column_title = "GWAS-ATAC disease regulons across fibrosis transitions",
     column_title_gp = gpar(fontsize = 11, fontface = "bold"),
     heatmap_legend_side = "right",
     padding = unit(c(4, 4, 4, 4), "mm"))

dev.off()

out_csv <- file.path(FIGS05_DIR, "figS05_gwas_atac_tf_transition_data.csv")
out_dt <- as.data.table(score_ord, keep.rownames = "tf")
out_dt_long <- melt(out_dt, id.vars = "tf", variable.name = "transition", value.name = "score")
fdr_dt <- as.data.table(fdr_ord, keep.rownames = "tf")
fdr_long <- melt(fdr_dt, id.vars = "tf", variable.name = "transition", value.name = "fdr")
out_long <- merge(out_dt_long, fdr_long, by = c("tf", "transition"))
out_long <- merge(out_long, disrupt_ord, by.x = "tf", by.y = "tf_name", all.x = TRUE)
fwrite(out_long, out_csv)

message("Wrote: ", out_pdf)
message("Wrote: ", out_csv)
