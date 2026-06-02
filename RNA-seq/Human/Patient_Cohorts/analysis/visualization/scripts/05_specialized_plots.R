
# 05_specialized_plots.R
# Generates Specialized visualizations (WGCNA Module Preservation)
# Updated: Feb 2026 (Sanjana Lab Theme)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

# Paths
BASE_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
WGCNA_DIR <- file.path(BASE_DIR, "Analysis/Cross_Species_Concordance/results")
RESULTS_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/visualization/plots/specialized"
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Source plotting functions & theme
FUNC_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/visualization/functions"
source(file.path(FUNC_DIR, "theme_publication.R"))

# --- 1. WGCNA Module Preservation ---
message("Generating WGCNA Module Preservation Plot...")

pres_file <- file.path(WGCNA_DIR, "wgcna_preservation_stats.csv")

if (file.exists(pres_file)) {
    pres_df <- fread(pres_file)
    setDF(pres_df)
    
    # Columns: diet, module, module_size, Zsummary, preservation
    # diet likely 'MCD'
    
    # Filter for valid Zsummary
    pres_df <- pres_df[!is.na(pres_df$Zsummary), ]
    
    # Define Preservation Categories explicitly for plotting if needed
    # Zsummary > 10: Strong
    # 2 < Zsummary < 10: Moderate
    # Zsummary < 2: None
    
    pres_df$PreservationCategory <- cut(pres_df$Zsummary, 
                                        breaks=c(-Inf, 2, 10, Inf), 
                                        labels=c("Not Preserved", "Moderately Preserved", "Strongly Preserved"))
    
    # Colors
    pres_colors <- c(
        "Strongly Preserved" = sanjana_colors[["Green"]],
        "Moderately Preserved" = sanjana_colors[["Blue"]],
        "Not Preserved" = sanjana_colors[["Grey"]]
    )
    
    # Plot Barplot
    p_pres <- ggplot(pres_df, aes(x=reorder(as.factor(module), Zsummary), y=Zsummary, fill=PreservationCategory)) +
        geom_bar(stat="identity", width=0.7) +
        geom_hline(yintercept=2, linetype="dashed", color=sanjana_colors[["Orange"]]) + # Thresholds
        geom_hline(yintercept=10, linetype="dashed", color=sanjana_colors[["Green"]]) +
        coord_flip() +
        theme_publication() +
        scale_fill_manual(values=pres_colors) +
        labs(title="WGCNA Module Preservation (Human -> Mouse MCD)",
             x="Module ID",
             y="Zsummary Statistic",
             fill="Preservation Status") +
        theme(legend.position="bottom")
        
    save_pdf(p_pres, file.path(RESULTS_DIR, "wgcna_module_preservation.pdf"), width=6, height=8)
}

message("Done.")
