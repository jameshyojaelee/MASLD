# figS05 B: motif quality tier distribution + per-TF survival across tiers
# Left:  stacked bar of 1,412 motifs by grade A/B/C
# Right: heatmap rows = disease regulon TFs, cols = tier_all / tier_AB / tier_Aonly
#        fill = disruption count if tf retained, NA if dropped at that tier
#        side-bar marks drug-target TFs in gold (THRB / FXR=NR1H4 / PPARG / ESR1)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(ComplexHeatmap)
  library(circlize)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

QUAL_F <- file.path(BASE, "GWAS/finemapping/results/gwas_atac/motif_quality_tier.csv")
TIER_F <- file.path(BASE, "GWAS/finemapping/results/gwas_atac/motif_tier_comparison.csv")

qual <- fread(QUAL_F)
tier <- fread(TIER_F)

# ---- Left panel: stacked bar of 1,412 motifs by grade ----
grade_counts <- qual[, .N, by = grade]
grade_counts[, grade := factor(grade, levels = c("A", "B", "C"))]
grade_counts <- grade_counts[order(grade)]
total_n <- sum(grade_counts$N)

grade_colors <- c(A = "#C9265E", B = "#1565C0", C = "#9E9E9E")
grade_counts[, lab := paste0(grade, " (n=", N, ")")]

# stacked single bar (one column)
grade_counts[, x := "All motifs"]

p_left <- ggplot(grade_counts, aes(x = x, y = N, fill = grade)) +
  geom_col(width = 0.45, color = "black", linewidth = 0.3) +
  geom_text(aes(label = lab),
            position = position_stack(vjust = 0.5),
            size = 2.4, fontface = "bold", color = "white") +
  scale_fill_manual(values = grade_colors, name = "Grade",
                    breaks = c("A", "B", "C"),
                    labels = c("A (high IC, SELEX)",
                               "B (curated motif)",
                               "C (low confidence)")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  labs(x = NULL,
       y = paste0("Motifs (n=", total_n, ")"),
       title = "Motif quality tiers") +
  theme_masld(base_size = 8) +
  theme_pub() +
  theme(
    axis.title.y    = element_text(face = "bold", color = "black"),
    axis.text.x     = element_text(face = "bold", color = "black"),
    plot.title      = element_text(face = "bold", color = "black"),
    legend.position = "right",
    legend.title    = element_text(face = "bold")
  )

# ---- Right panel: TF survival heatmap across tier_all / tier_AB / tier_Aonly ----
dis_tfs <- tier[in_disease_regulon == TRUE]
# matrix: rows = TFs, cols = tier policies; fill = disruption count if kept else NA
mat <- as.matrix(dis_tfs[, .(
  tier_all   = ifelse(keep_all,   n_disruptions, NA_real_),
  tier_AB    = ifelse(keep_AB,    n_disruptions, NA_real_),
  tier_Aonly = ifelse(keep_Aonly, n_disruptions, NA_real_)
)])
rownames(mat) <- dis_tfs$tf_name
colnames(mat) <- c("tier_all", "tier_AB", "tier_Aonly")
# order by total disruptions
mat <- mat[order(-dis_tfs$n_disruptions), , drop = FALSE]
tf_order <- rownames(mat)

# drug-target gold side bar: known MASLD/MASH drug targets
drug_targets <- c("THRB", "NR1H4", "PPARG", "ESR1", "PPARA", "FXR")
drug_anno <- ifelse(tf_order %in% drug_targets, "Drug target", "Other")

vmax <- max(mat, na.rm = TRUE)
col_fun <- colorRamp2(c(0, vmax / 2, vmax),
                      c("#F7F7F7", "#F4A674", "#C9265E"))

right_anno <- rowAnnotation(
  Drug = anno_simple(
    drug_anno,
    col = c("Drug target" = "#FFB300", "Other" = "white"),
    border = FALSE,
    width = unit(0.25, "cm")
  ),
  show_annotation_name = FALSE
)

ht <- Heatmap(
  mat,
  name = "Disruptions\n(retained)",
  col = col_fun,
  na_col = "#E0E0E0",
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 7,
                     fontface = ifelse(tf_order %in% drug_targets,
                                       "bold.italic", "plain")),
  column_labels = c("all", "A+B", "A only"),
  column_names_rot = 0, column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 7),
  cell_fun = function(j, i, x, y, w, h, fill) {
    v <- mat[i, j]
    if (is.na(v)) {
      grid.text("-", x, y,
                gp = gpar(fontsize = 6, col = "gray40"))
    } else {
      grid.text(v, x, y,
                gp = gpar(fontsize = 6,
                          col = if (v > vmax * 0.6) "white" else "black",
                          fontface = "bold"))
    }
  },
  right_annotation = right_anno,
  rect_gp = gpar(col = "white", lwd = 0.5),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 6),
    legend_height = unit(1.5, "cm")
  ),
  width  = unit(2.8, "cm"),
  height = unit(0.30 * nrow(mat), "cm"),
  column_title = "Disease regulon TF survival",
  column_title_gp = gpar(fontsize = 8, fontface = "bold")
)

# capture heatmap to a grob via grid graphics
ht_grob <- grid.grabExpr({
  draw(ht,
       heatmap_legend_side = "right",
       padding = unit(c(2, 2, 2, 2), "mm"))
}, width = 4, height = 4)

# combine left ggplot + right heatmap with patchwork
p_combined <- p_left + wrap_elements(ht_grob) + plot_layout(widths = c(1, 1.4))

out_pdf <- file.path(FIGS05_DIR, "figS05_b_motif_quality_tier.pdf")
dir.create(FIGS05_DIR, showWarnings = FALSE, recursive = TRUE)
ggsave(out_pdf, p_combined, width = 8, height = 4, useDingbats = FALSE)
message("Wrote: ", out_pdf)

# data export
out_csv <- sub("\\.pdf$", "_tf_survival.csv", out_pdf)
fwrite(dis_tfs, out_csv)
fwrite(grade_counts, sub("\\.pdf$", "_grade_counts.csv", out_pdf))
message("Wrote: ", out_csv)
