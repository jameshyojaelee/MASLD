
# 02_integrated_plots.R
# Generates visualizations for Human Integrated Analysis (Meta-analysis, Consensus, GSEA)
# Updated: Feb 2026 (Sanjana Lab Theme)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(EnhancedVolcano)
  library(ComplexHeatmap)
  library(circlize)
  library(forestplot)
  library(grid)
})

# Paths
BASE_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT_DIR <- file.path(BASE_DIR, "analysis/integration")
RESULTS_DIR <- file.path(BASE_DIR, "analysis/visualization/plots/integrated")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Source plotting functions & theme
FUNC_DIR <- file.path(BASE_DIR, "analysis/visualization/functions")
source(file.path(FUNC_DIR, "theme_publication.R"))
source(file.path(FUNC_DIR, "plot_volcano.R"))

# Load Data
message("Loading integrated results...")
meta_res_file <- file.path(INT_DIR, "results/integration/meta_analysis_results.csv") 
gsea_file <- file.path(INT_DIR, "results/integration/gsea_results.csv") 

# Check files
if (!file.exists(meta_res_file)) stop("Meta-analysis results not found.")

meta_res <- fread(meta_res_file)
setDF(meta_res)

# 1. Consensus Volcano Plot
message("Generating Consensus Volcano Plot...")
vol_data <- meta_res
# Map columns
# Map columns
if ("meta_logFC" %in% names(vol_data)) {
  vol_data$log2FoldChange <- vol_data$meta_logFC
  vol_data$padj <- vol_data$meta_padj
}

# Load Mapping
map_file <- file.path(INT_DIR, "results/integration/human_mouse_ortholog_comparison.csv")
if (file.exists(map_file)) {
    map_df <- fread(map_file)
    id_map <- map_df[, c("gene_base", "human_symbol")]
    id_map <- id_map[!duplicated(id_map$gene_base), ]
    
    # Map
    if ("gene" %in% names(vol_data)) {
        vol_data$gene_base <- gsub("\\..*", "", vol_data$gene)
        vol_data <- merge(vol_data, id_map, by="gene_base", all.x=TRUE)
        vol_data$symbol <- ifelse(is.na(vol_data$human_symbol) | vol_data$human_symbol == "", vol_data$gene, vol_data$human_symbol)
    }
} else {
    if (!"symbol" %in% names(vol_data) && "gene" %in% names(vol_data)) vol_data$symbol <- vol_data$gene
}

p_vol <- plot_volcano(vol_data, 
                      title="Meta-Analysis Consensus: Disease vs Control", 
                      subtitle = paste0("Combined P-values across cohorts"),
                      fc_cutoff=1, p_cutoff=0.05)
save_pdf(p_vol, file.path(RESULTS_DIR, "consensus_volcano.pdf"), width=6, height=6) 

# 2. Forest Plots for Top Genes
message("Generating Forest Plots for Top 10 Genes...")
top_genes <- vol_data[order(vol_data$padj), "symbol"][1:10]

# Load per-study results
per_study_files <- list.files(file.path(INT_DIR, "results/per_study"), pattern="_de_results.csv", full.names=TRUE)
study_data <- list()
for (f in per_study_files) {
  ds_name <- gsub("_de_results.csv", "", basename(f))
  dt <- fread(f)
  setDF(dt)
  if (!"symbol" %in% names(dt) && "gene" %in% names(dt)) dt$symbol <- dt$gene
  study_data[[ds_name]] <- dt
}

# Setup PDF for Forest Plots (multi-page or single grid?)
# Guidelines: Editable text. Forestplot package usually uses grid graphics.
pdf(file.path(RESULTS_DIR, "forest_plots_top10.pdf"), width=8, height=4, useDingbats = FALSE)
for (g in top_genes) {
  if (is.na(g)) next
  
  # ... (Data collection logic same as before) ...
  studies <- names(study_data)
  betas <- c()
  se <- c()
  
  for (s in studies) {
    row <- study_data[[s]][study_data[[s]]$symbol == g, ]
    if (nrow(row) > 0) {
      if ("logFC" %in% names(row)) b <- row$logFC[1] else b <- NA
      # Helper for SE
      if ("lfcSE" %in% names(row)) err <- row$lfcSE[1] 
      else if ("t" %in% names(row) && !is.na(b)) err <- abs(b / row$t[1]) 
      else err <- NA
      
      betas <- c(betas, b)
      se <- c(se, err)
    } else {
      betas <- c(betas, NA)
      se <- c(se, NA)
    }
  }
  
  # Add meta-analysis result
  meta_row <- vol_data[vol_data$symbol == g, ]
  betas <- c(betas, meta_row$meta_logFC)
  se <- c(se, meta_row$meta_SE)
  
  labels <- c(studies, "Meta-Analysis")
  
  # Filter NAs
  valid <- !is.na(betas) & !is.na(se)
  if (sum(valid) < 2) next
  
  # Forest Plot with Theme Fonts?
  # Forestplot allows fpTxtGp for gpar settings
  txt_gp <- fpTxtGp(label = gpar(fontfamily = "Helvetica", cex = 0.8),
                    ticks = gpar(fontfamily = "Helvetica", cex = 0.7),
                    xlab  = gpar(fontfamily = "Helvetica", cex = 0.8),
                    title = gpar(fontfamily = "Helvetica", fontface = "bold", cex = 1.0))
  
  forestplot(labeltext = labels[valid],
             mean = betas[valid],
             lower = betas[valid] - 1.96*se[valid],
             upper = betas[valid] + 1.96*se[valid],
             title = paste0("Forest Plot: ", g),
             xlab = "log2 Fold Change",
             zero = 0,
             boxsize = 0.2,
             lineheight = unit(8, "mm"),
             colgap = unit(4, "mm"),
             lwd.ci = 1.5, 
             ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             txt_gp = txt_gp, # Apply font theme
             col = fpColors(box = sanjana_colors[["Magenta"]], line = "black", summary = sanjana_colors[["Purple"]]))
}
dev.off()

# 3. GSEA Dotplot
if (file.exists(gsea_file)) {
  message("Generating GSEA Dotplot...")
  gsea_res <- fread(gsea_file)
  setDF(gsea_res)
  
  if ("NES" %in% names(gsea_res) && "padj" %in% names(gsea_res)) {
    gsea_res <- gsea_res[gsea_res$padj < 0.05, ]
    up <- gsea_res[gsea_res$NES > 0, ]
    down <- gsea_res[gsea_res$NES < 0, ]
    
    # Take top 10
    top_up <- up[order(up$padj), ][1:min(10, nrow(up)), ]
    top_down <- down[order(down$padj), ][1:min(10, nrow(down)), ]
    
    plot_data <- rbind(top_up, top_down)
    plot_data$pathway <- gsub("HALLMARK_", "", plot_data$pathway) # Clean names
    plot_data$pathway <- gsub("_", " ", plot_data$pathway)
    plot_data$pathway <- factor(plot_data$pathway, levels=plot_data$pathway[order(plot_data$NES)])
    
    # Custom Gradient: Blue (Down) -> Grey -> Red (Up)? Or Sanjana Palette?
    # Guidelines don't specify gradient for dotplots. 
    # Use Blue -> Red standard for NES, compatible with palette.
    # Alternatively: Green/Orange? 
    # Let's use sanjana_colors[["Blue"]] -> sanjana_colors[["Magenta"]] (or Red)
    
    p_gsea <- ggplot(plot_data, aes(x=NES, y=pathway, color=padj, size=-log10(padj))) +
      geom_point() +
      scale_color_gradient(low=sanjana_colors[["Magenta"]], high=sanjana_colors[["Blue"]]) + # Low p (high sig) = Magenta? No, usually color by NES or p-val.
      # Let's color by NES if we want directionality, or p-val for sig.
      # Standard is usually color by p-val (red=sig, blue=less sig) or direction.
      # Let's stick to standard Blue-Red gradient for p-val (low p = Red/Magenta).
      scale_color_gradient(low=sanjana_colors[["Magenta"]], high=sanjana_colors[["Blue"]]) +
      labs(title="Top Enriched Pathways (GSEA)", x="Normalized Enrichment Score (NES)", y="", color="Adj P-val") +
      geom_vline(xintercept=0, linetype="dashed", color="grey50") +
      theme_publication()
      
    save_pdf(p_gsea, file.path(RESULTS_DIR, "gsea_dotplot.pdf"), width=7, height=6)
  }
}

message("Done.")
