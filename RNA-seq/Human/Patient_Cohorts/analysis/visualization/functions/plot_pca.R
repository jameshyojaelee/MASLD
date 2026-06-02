
library(ggplot2)
library(ggrepel)

source("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/visualization/functions/theme_publication.R")

plot_pca <- function(vst_data, metadata, color_by="condition", shape_by=NULL, title="PCA Plot") {
  # Calculate PCA
  pca <- prcomp(t(vst_data), center=TRUE, scale.=TRUE)
  percentVar <- round(100 * summary(pca)$importance[2, 1:2])
  
  # Prepare data for plotting
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], metadata)
  
  # Base plot
  p <- ggplot(d, aes_string(x="PC1", y="PC2", color=color_by, shape=shape_by)) +
    geom_point(size=2, alpha=0.9) + # Slightly smaller points
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle(title) +
    theme_publication() +
    scale_color_manual(values = condition_colors) + # Use strict palette
    theme(aspect.ratio=1)
    
  return(p)
}
