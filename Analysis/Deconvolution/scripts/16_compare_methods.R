#!/usr/bin/env Rscript
# Compare Deconvolution Results: MuSiC vs InstaPrism
# Input: Results directories for MuSiC and InstaPrism
# Output: Correlation plots, stacked bar plots comparison

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(pheatmap)
  library(gridExtra)
})

# Configure PDF device
pdf.options(useDingbats = FALSE)

# Parameters — resolve from env var or default to correct project-relative path
project_root <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RESULTS_BASE <- file.path(project_root, "Analysis/Deconvolution/results")
OUTPUT_DIR <- file.path(project_root, "Analysis/Deconvolution/comparison_plots")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# All 8 datasets with both MuSiC and BayesPrism results.
# GSE126848, PRJNA512027, GSE167523 require BayesPrism backfill
# (run_bayesprism_missing_human.sbatch) before this script will produce
# meaningful output for those three. Missing BayesPrism files are silently
# skipped via warning() in load_props().
datasets <- c(
  "inhouse_MCD", "GSE156918", "GSE205974",   # mouse
  "GSE130970", "GSE135251",                   # human (original)
  "GSE126848", "PRJNA512027", "GSE167523"     # human (backfill)
)

# Function to load proportions
load_props <- function(dataset, method) {
  if (method == "music") {
    file_path <- file.path(RESULTS_BASE, dataset, paste0(dataset, "_music_prop_weighted.tsv"))
  } else if (method == "bayesprism") {
    file_path <- file.path(RESULTS_BASE, dataset, paste0(dataset, "_bayesprism_proportions.tsv"))
  } else if (method == "cibersortx") {
    file_path <- file.path(RESULTS_BASE, dataset, "cibersortx", "CIBERSORTx_Results.txt")
  }

  if (!file.exists(file_path)) {
    warning(paste("Missing file:", file_path))
    return(NULL)
  }
  
  props <- read.delim(file_path, row.names = 1, check.names = FALSE)
  props$Sample <- rownames(props)
  props$Method <- method
  props$Dataset <- dataset
  
  # Normalize to sum to 1 (just in case)
  # props_mat <- as.matrix(props[, !colnames(props) %in% c("Sample", "Method", "Dataset")])
  # props_mat <- sweep(props_mat, 1, rowSums(props_mat), "/")
  # props[, !colnames(props) %in% c("Sample", "Method", "Dataset")] <- props_mat
  
  pivot_longer(props, cols = -c(Sample, Method, Dataset), 
               names_to = "CellType", values_to = "Proportion")
}

# Load all data (CIBERSORTx silently skipped if not yet run)
all_data <- list()
for (ds in datasets) {
  music <- load_props(ds, "music")
  bp <- load_props(ds, "bayesprism")
  cbx <- load_props(ds, "cibersortx")
  all_data[[ds]] <- bind_rows(music, bp, cbx)
}
combined_df <- bind_rows(all_data)

# standardize cell type names if needed (e.g. MuSiC might enforce make.names)

# Plot 1: Correlation between methods for each dataset
cat("Generating correlation plots...\n")
pdf(file.path(OUTPUT_DIR, "method_correlation.pdf"), width = 10, height = 8)
for (ds in datasets) {
  df_ds <- combined_df %>% filter(Dataset == ds)
  
  # Spread to wide for correlation
  df_wide <- df_ds %>% 
    pivot_wider(names_from = Method, values_from = Proportion) %>%
    na.omit() # Remove cell types present in one but not other
    
  if(nrow(df_wide) == 0) next
  
  overall_cor <- cor(df_wide$music, df_wide$bayesprism, method = "pearson")
  
  p <- ggplot(df_wide, aes(x = music, y = bayesprism, color = CellType)) +
    geom_point(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    theme_minimal(base_family = "Helvetica", base_size = 7) +
    labs(title = paste0("MuSiC vs InstaPrism: ", ds),
         subtitle = paste("Pearson R =", round(overall_cor, 3)),
         x = "MuSiC Proportion", y = "InstaPrism Proportion") +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold", size = 8),
          axis.text = element_text(size = 6),
          axis.title = element_text(size = 8))
  print(p)
}
dev.off()

# Plot 2: Detailed Bar Plots comparison
cat("Generating stacked bar plots comparison...\n")
pdf(file.path(OUTPUT_DIR, "method_comparison_barplots.pdf"), width = 12, height = 8)
for (ds in datasets) {
  df_ds <- combined_df %>% filter(Dataset == ds)
  
  # Order samples for consistency
  samples <- unique(df_ds$Sample)
  
  # Create plot
  p <- ggplot(df_ds, aes(x = Sample, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~Method, ncol = 1, scales = "free_x") +
    theme_minimal(base_family = "Helvetica", base_size = 7) +
    labs(title = paste("Deconvolution Comparison:", ds), x = NULL) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 8),
          axis.text = element_text(size = 6),
          axis.title = element_text(size = 8))
  print(p)
}
dev.off()

cat("Analysis complete. Results in:", OUTPUT_DIR, "\n")
