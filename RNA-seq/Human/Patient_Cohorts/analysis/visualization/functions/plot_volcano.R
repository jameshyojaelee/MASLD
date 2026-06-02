
library(EnhancedVolcano)
library(ggplot2)

source("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/visualization/functions/theme_publication.R")

plot_volcano <- function(res, title, subtitle=NULL, fc_cutoff=1, p_cutoff=0.05, top_n=10) {
  # Validate input
  if (!all(c("log2FoldChange", "padj", "symbol") %in% colnames(res))) {
    stop("Input dataframe must contain log2FoldChange, padj, and symbol columns")
  }

  # Clean NA values
  res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
  
  # Identify top genes for labeling
  top_genes <- res[order(res$padj), "symbol"][1:top_n]

  # Colors based on guidelines
  # NS: Grey, LFC: Blue, P: Blue, Both: Magenta (Primary)
  colors <- c(sanjana_colors[["Grey"]], sanjana_colors[["Blue"]], sanjana_colors[["Blue"]], sanjana_colors[["Magenta"]])
  names(colors) <- c("NS", "Log2 FC", "P", "FC & P")

  # Set default subtitle if NULL
  if (is.null(subtitle)) {
      subtitle <- paste0("FDR < ", p_cutoff, ", |LFC| > ", fc_cutoff)
  }

  # Create plot
  p <- EnhancedVolcano(res,
    lab = res$symbol,
    x = 'log2FoldChange',
    y = 'padj',
    title = title,
    subtitle = subtitle,
    pCutoff = p_cutoff,
    FCcutoff = fc_cutoff,
    pointSize = 1.5, # Smaller points for cleaner aesthetic
    labSize = 2.5, # ~7pt
    col = c(sanjana_colors[["Grey"]], sanjana_colors[["Blue"]], sanjana_colors[["Blue"]], sanjana_colors[["Magenta"]]),
    colAlpha = 0.8,
    legendPosition = 'right',
    legendLabSize = 7,
    legendIconSize = 2.0,
    drawConnectors = TRUE,
    widthConnectors = 0.3,
    selectLab = top_genes,
    gridlines.major = FALSE,
    gridlines.minor = FALSE
  ) + 
  theme_publication() +
  theme(legend.position = "right") # Ensure legend stays
  
  return(p)
}
