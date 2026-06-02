
library(ComplexHeatmap)
library(circlize)

source("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/visualization/functions/theme_publication.R")

plot_heatmap <- function(mat, metadata, top_annotation_col="condition", title="Heatmap") {
  # Scale rows (Z-score)
  mat_scaled <- t(scale(t(mat)))
  
  # Create annotation
  if (!is.null(metadata) && top_annotation_col %in% colnames(metadata)) {
    # Use global condition_colors if available, ensuring match
    # Filter to levels present in metadata to ensure exact matching if needed, but ComplexHeatmap handles it
    
    ha = HeatmapAnnotation(df = metadata[, top_annotation_col, drop=FALSE],
                           col = list(condition = condition_colors),
                           annotation_name_gp = gpar(fontsize = 7, fontfamily="Helvetica"),
                           simple_anno_size = unit(3, "mm"))
  } else {
    ha = NULL
  }
  
  # Create Heatmap
  hm <- Heatmap(mat_scaled,
    name = "Z-score",
    top_annotation = ha,
    column_title = title,
    column_title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily="Helvetica"),
    row_names_gp = gpar(fontsize = 6, fontfamily="Helvetica"),
    column_names_gp = gpar(fontsize = 6, fontfamily="Helvetica"),
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    # Color Ramp: Blue -> White -> Red (Standard) or Purple/Magenta based?
    # Guidelines don't specify heatmap ramp, but standard is usually Blue-Red or Viridis.
    # Sanjana palette has Blue/Magenta/Orange. 
    # Let's keep Blue-White-Red as it's standard for Z-score (Cold to Hot)
    col = colorRamp2(c(-2, 0, 2), c("#4baeef", "white", "#e14b9d")), # Blue to Magenta! Matches palette better.
    
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily="Helvetica"),
      labels_gp = gpar(fontsize = 6, fontfamily="Helvetica")
    )
  )
  
  return(hm)
}
