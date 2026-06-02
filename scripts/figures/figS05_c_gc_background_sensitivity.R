# figS05 C: GC background sensitivity for motif disruption calls
# Rows = disease regulon TFs, cols = uniform / genome / peak GC backgrounds
# Side-bar = in-disease-regulon flag (warm magenta vs gray)
# 4-way validated TFs (THRB, HNF4A, RORA, MLXIPL, MAX) labelled bold italic.

suppressPackageStartupMessages({
  library(data.table)
  library(grid)
  library(ComplexHeatmap)
  library(circlize)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

BG_F   <- file.path(BASE, "GWAS/finemapping/results/gwas_atac/motif_background_comparison.csv")
TIER_F <- file.path(BASE, "GWAS/finemapping/results/gwas_atac/motif_tier_comparison.csv")
REG_F  <- file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv")
REG_LENIENT_F <- file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons_lenient.csv")

bg   <- fread(BG_F)
tier <- fread(TIER_F)

# disease regulon TF set: prefer lenient (broader); else strict
if (file.exists(REG_LENIENT_F) && file.info(REG_LENIENT_F)$size > 100) {
  reg <- fread(REG_LENIENT_F)
} else {
  reg <- fread(REG_F)
}
dis_tfs_full <- unique(reg$tf_name)
# Backup: also pull from tier file flag
dis_tfs_full <- unique(c(dis_tfs_full, tier[in_disease_regulon == TRUE, tf_name]))

# 4-way validated TFs (call out)
four_way <- c("THRB", "HNF4A", "RORA", "MLXIPL", "MAX")

# Per the design brief, target ~18 TFs: disease-regulon TFs that survive in bg
keep_tfs <- intersect(bg$tf_name, dis_tfs_full)
keep_tfs <- union(keep_tfs, intersect(bg$tf_name, four_way))

mat_dt <- bg[tf_name %in% keep_tfs, .(tf_name, uniform, genome, peak, total)]
setorder(mat_dt, -total)
mat_dt <- mat_dt[1:min(.N, 18)]
mat <- as.matrix(mat_dt[, .(uniform, genome, peak)])
rownames(mat) <- mat_dt$tf_name
tf_order <- rownames(mat)

# disease regulon side-bar
in_dr <- ifelse(tf_order %in% tier[in_disease_regulon == TRUE, tf_name],
                "Disease regulon", "Other")

vmax <- max(mat, na.rm = TRUE)
col_fun <- colorRamp2(c(0, vmax / 2, vmax),
                      c("#F7F7F7", "#F4A674", "#C9265E"))

right_anno <- rowAnnotation(
  Regulon = anno_simple(
    in_dr,
    col = c("Disease regulon" = "#C9265E", "Other" = "#9E9E9E"),
    border = FALSE,
    width = unit(0.25, "cm")
  ),
  annotation_name_gp = gpar(fontsize = 6, fontface = "bold"),
  annotation_name_rot = 0
)

# row name face: 4-way validated TFs in bold italic
row_face <- ifelse(tf_order %in% four_way, "bold.italic", "plain")

ht <- Heatmap(
  mat,
  name = "Disruptions",
  col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 7, fontface = row_face),
  column_labels = c("Uniform", "Genome", "Peak"),
  column_names_rot = 0, column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 7, fontface = "bold"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    v <- mat[i, j]
    grid.text(v, x, y,
              gp = gpar(fontsize = 6,
                        col = if (v > vmax * 0.6) "white" else "black",
                        fontface = "bold"))
  },
  right_annotation = right_anno,
  rect_gp = gpar(col = "white", lwd = 0.5),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 6),
    legend_height = unit(2, "cm")
  ),
  width  = unit(2.6, "cm"),
  height = unit(0.30 * nrow(mat), "cm")
)

out_pdf <- file.path(FIGS05_DIR, "figS05_c_gc_background_sensitivity.pdf")
dir.create(FIGS05_DIR, showWarnings = FALSE, recursive = TRUE)
pdf(out_pdf, width = 5, height = 6, useDingbats = FALSE)
draw(ht,
     column_title = "GC background sensitivity (disease regulon TFs)",
     column_title_gp = gpar(fontsize = 9, fontface = "bold"),
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     padding = unit(c(4, 4, 4, 4), "mm"))
grid.text("Bold italic = 4-way validated TF (THRB / HNF4A / RORA / MLXIPL / MAX)",
          x = 0.5, y = 0.03,
          gp = gpar(fontsize = 6, col = "gray30"))
dev.off()
message("Wrote: ", out_pdf)

out_csv <- sub("\\.pdf$", ".csv", out_pdf)
fwrite(mat_dt, out_csv)
message("Wrote: ", out_csv)
