
# generate_advanced_plots.R
# Generate comprehensive plots for ssGSEA and DoRothEA analysis
# strict dataset order enforced: Hoang, Govaere, InHouse_MCD, Ext_MCD_1, Ext_MCD_2

# Environment setup
conda_lib <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/pathway_analysis/.mamba/pathway_analysis/lib/R/library"
if (dir.exists(conda_lib)) {
    .libPaths(c(conda_lib))
}

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(pheatmap)
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
})

# Configuration
BASE_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/pathway_analysis"
OUTPUT_DIR <- file.path(BASE_DIR, "results/refactored")
PLOT_DIR <- file.path(BASE_DIR, "plots/refactored")
THEME_PATH <- file.path(BASE_DIR, "scripts/publication_theme.R")

if (file.exists(THEME_PATH)) {
    source(THEME_PATH)
} else {
    MAIN_COLOR <- "red"
}

# Ensure Plot Dir
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Dataset Order
DATASET_ORDER <- c("Hoang", "Govaere", "InHouse_MCD", "Ext_MCD_1", "Ext_MCD_2")

# 1. Load data
cat("Loading results...\n")
ssgsea_df <- readRDS(file.path(OUTPUT_DIR, "all_ssgsea_scores.rds"))
dorothea_list <- readRDS(file.path(OUTPUT_DIR, "all_dorothea_scores.rds"))

# Enforce factor levels for ssgsea_df
ssgsea_df <- ssgsea_df %>%
    mutate(Dataset = factor(Dataset, levels = DATASET_ORDER)) %>%
    filter(!is.na(Dataset)) # Safety filter

# Process DoRothEA list to DF
tf_df <- bind_rows(lapply(names(dorothea_list), function(n) {
    if (is.null(dorothea_list[[n]])) return(NULL)
    # decoupleR output is usually: condition, source, score, statistic...
    dorothea_list[[n]] %>% 
        mutate(Dataset = n, TF = toupper(source)) %>% # Normalize TF names
        dplyr::select(Sample = condition, TF, Score = score, Dataset)
}))

if (!is.null(tf_df)) {
    tf_df <- tf_df %>%
        mutate(Dataset = factor(Dataset, levels = DATASET_ORDER)) %>%
        filter(!is.na(Dataset))
}

# ----------------------------------------------------------------------------
# Plot 0: Fixed ssGSEA Validation Boxplot
# ----------------------------------------------------------------------------
cat("Generating ssGSEA Validation Boxplot...\n")

target_pathways <- c(
    "HALLMARK_FATTY_ACID_METABOLISM",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_APOPTOSIS"
)

# Clean names for plotting
clean_names <- gsub("HALLMARK_", "", target_pathways)
names(clean_names) <- target_pathways

p_box <- ssgsea_df %>%
    filter(Pathway %in% target_pathways) %>%
    mutate(PathwayName = factor(Pathway, levels = target_pathways, labels = clean_names)) %>%
    ggplot(aes(x = Dataset, y = Score, fill = Dataset)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, lwd=0.3) +
    geom_jitter(width = 0.2, size = 0.3, alpha = 0.3) +
    facet_wrap(~PathwayName, scales = "fixed", ncol = 3) +
    theme_minimal(base_size = 10) +
    theme(
        axis.text.x = element_blank(), # Remove X text for compactness (Legend suffices)
        axis.ticks.x = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm")
    ) +
    labs(title = "Comprehensive Pathway Validation", x = NULL, y = "ssGSEA Score")

ggsave(file.path(PLOT_DIR, "ssgsea_validation_boxplot.pdf"), p_box, width = 10, height = 7)




# ----------------------------------------------------------------------------
# Plot 1: Heatmap of Top Variable ssGSEA Pathways (Mean per Dataset)
# ----------------------------------------------------------------------------
cat("Generating ssGSEA Heatmap...\n")

# Calculate mean score per Dataset
pathway_means <- ssgsea_df %>%
    group_by(Pathway, Dataset) %>%
    summarise(MeanScore = mean(Score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Dataset, values_from = MeanScore) %>%
    column_to_rownames("Analysis/downstream_analysis/pathway_analysis")

    # pathway_means <- na.omit(pathway_means) # Too strict
    pathway_means <- pathway_means[rowSums(is.na(pathway_means)) < ncol(pathway_means), , drop = FALSE]
    
    
    if (nrow(pathway_means) < 5) {
        warning("Too few shared pathways across datasets for heatmap.")
    } else {
        # Select top 30 most variable pathways across datasets (ignoring NAs for var calc)
        vars <- apply(pathway_means, 1, var, na.rm = TRUE)
        top_pathways <- names(sort(vars, decreasing = TRUE))[1:min(30, length(vars), nrow(pathway_means))]
        plot_mat <- pathway_means[top_pathways, , drop = FALSE]
        
        # Ensure column order matches desired order (if present in data)
        existing_cols <- intersect(DATASET_ORDER, colnames(plot_mat))
        plot_mat <- plot_mat[, existing_cols, drop=FALSE]

        p1_file <- file.path(PLOT_DIR, "ssgsea_pathway_dataset_heatmap.pdf")
        pdf(p1_file, width = 8, height = 10)
        pheatmap::pheatmap(as.matrix(plot_mat), 
                 scale = "row", 
                 cluster_cols = FALSE, # Preserve requested order
                 main = "Top Variable ssGSEA Pathways (Unified Rank Space)",
                 color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                 border_color = "grey80",
                 fontsize_row = 8)
        dev.off()
    }

# ----------------------------------------------------------------------------
# Plot 2: Heatmap of Top Variable TFs (DoRothEA)
# ----------------------------------------------------------------------------
cat("Generating DoRothEA TF Heatmap...\n")

if (!is.null(tf_df)) {
    tf_means <- tf_df %>%
        group_by(TF, Dataset) %>%
        summarise(MeanActivity = mean(Score, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = Dataset, values_from = MeanActivity) %>%
        column_to_rownames("TF")
    
    # Handle NAs
    # tf_means <- na.omit(tf_means) # Too strict
    tf_means <- tf_means[rowSums(is.na(tf_means)) < ncol(tf_means), , drop=FALSE]
    
    if (nrow(tf_means) < 5) {
         warning("Too few shared TFs for heatmap.")
    } else {
        # Calculate Variance
        tf_vars <- apply(tf_means, 1, var)
        top_tfs <- names(sort(tf_vars, decreasing = TRUE))[1:min(30, length(tf_vars), nrow(tf_means))]
        tf_plot_mat <- tf_means[top_tfs, ]
        
        # Ensure column order
        existing_cols <- intersect(DATASET_ORDER, colnames(tf_plot_mat))
        tf_plot_mat <- tf_plot_mat[, existing_cols, drop=FALSE]
        
        p2_file <- file.path(PLOT_DIR, "dorothea_tf_dataset_heatmap.pdf")
        pdf(p2_file, width = 8, height = 10)
        pheatmap::pheatmap(as.matrix(tf_plot_mat), 
                 scale = "row",
                 cluster_cols = FALSE, # Preserve requested order
                 main = "Top Variable TF Activities (Mean per Dataset)",
                 color = colorRampPalette(c("purple", "white", "orange"))(100),
                 border_color = "grey80",
                 fontsize_row = 8)
        dev.off()
    }
}

# ----------------------------------------------------------------------------
# Plot 3: Combined Dot Plot of Core Inflammatory Pathways
# ----------------------------------------------------------------------------
cat("Generating Pathway Dot Plot...\n")

core_pathways <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE", 
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)

summ_stats <- ssgsea_df %>%
    filter(Pathway %in% core_pathways) %>%
    group_by(Dataset, Pathway) %>%
    summarise(
        Mean = mean(Score, na.rm=TRUE),
        SE = sd(Score, na.rm=TRUE) / sqrt(n()),
        .groups = "drop"
    )

p3 <- ggplot(summ_stats, aes(x = Dataset, y = Pathway)) +
    geom_point(aes(size = abs(Mean), color = Mean)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = "grey90")
    ) +
    labs(title = "Core Pathway Activity Summary", 
         size = "Abs Mean Score", 
         color = "Mean Score",
         y = NULL, x = NULL)

ggsave(file.path(PLOT_DIR, "ssgsea_core_pathways_dotplot.pdf"), p3, width = 8, height = 6)

# ----------------------------------------------------------------------------
# Plot 4: Top Active TFs Consensus
# ----------------------------------------------------------------------------
cat("Generating TF Consensus Plot...\n")

if (!is.null(tf_df)) {
    # Find TFs that are top active in most datasets
    top_tfs_per_dataset <- tf_df %>%
        group_by(Dataset, TF) %>%
        summarise(Mean = mean(Score, na.rm=TRUE), .groups="drop") %>%
        group_by(Dataset) %>%
        slice_max(order_by = abs(Mean), n = 10) %>%
        pull(TF) %>%
        unique()
    
    # Filter for these TFs
    tf_subset <- tf_df %>%
        filter(TF %in% top_tfs_per_dataset) %>%
        group_by(Dataset, TF) %>%
        summarise(Mean = mean(Score, na.rm=TRUE), .groups="drop")
    
    p4 <- ggplot(tf_subset, aes(x = Dataset, y = TF)) +
        geom_tile(aes(fill = Mean), color = "white") +
        scale_fill_gradient2(low = "purple", mid = "white", high = "orange") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Top Active TFs (Consensus)", x = NULL, y = NULL, fill = "Activity")
        
    ggsave(file.path(PLOT_DIR, "dorothea_top_tfs_consensus.pdf"), p4, width = 10, height = 8)
}

cat("Done. Plots saved to", PLOT_DIR, "\n")
