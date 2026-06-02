#!/usr/bin/env Rscript
# generate_ultimate_plot.R - Create the "Ultimate" Enrichment Landscape Visualization
#
# This script integrates GSEA and ORA results to produce a single, comprehensive
# summary figure ("The Landscape Plot") that captures:
# 1. Biological Theme (Facet by Collection)
# 2. Directionality & Magnitude (x-axis = GSEA NES)
# 3. Statistical Confidence (Color = p.adjust)
# 4. Multi-Method Validation (Point Shape/Mark = ORA Confirmation)
#
# Usage: Rscript scripts/generate_ultimate_plot.R

# Critical: Exclude user library to avoid conflicts
.libPaths(.libPaths()[!grepl("jameslee", .libPaths())])

suppressPackageStartupMessages({
    library(tidyverse)
    library(cowplot) # Replacement for patchwork
    library(scales)
    # library(ggnewscale) # Missing
})

# ============================================================================
# Configuration
# ============================================================================
PATHWAY_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/pathway_analysis"
ORA_DIR <- file.path(PATHWAY_DIR, "results/ora")
GSEA_DIR <- file.path(PATHWAY_DIR, "results/gsea")
PLOT_DIR <- file.path(PATHWAY_DIR, "plots/summary")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load publication theme
source(file.path(PATHWAY_DIR, "scripts/publication_theme.R"))

cat("=== Generating The ULTIMATE Enrichment Landscape ===\n\n")

# ============================================================================
# 1. Load & Integrate Data
# ============================================================================
cat("[1/4] Loading and integrating results...\n")

# --- Load GSEA Results (The Backbone) ---
load(file.path(GSEA_DIR, "all_gsea_results.RData"))
# Combine list into partial dataframe
gsea_df <- bind_rows(lapply(names(gsea_results), function(n) {
    gsea_results[[n]] %>% 
        as_tibble() %>% 
        select(pathway, NES, padj, size) %>% 
        mutate(collection = n)
}))

# --- Load ORA Results (The Validator) ---
load(file.path(ORA_DIR, "all_ora_results.RData"))
# Helper to extract significant ORA terms
extract_sig_ora <- function(res, name) {
    if (is.null(res)) return(tibble())
    # Input is now a tibble
    res %>% 
        filter(p.adjust < 0.1) %>% 
        select(ID) %>% 
        mutate(collection = name, is_ora_sig = TRUE)
}

ora_sig_df <- bind_rows(
    extract_sig_ora(ora_hallmark, "hallmark"),
    extract_sig_ora(ora_gobp, "go_bp"),
    extract_sig_ora(ora_kegg, "kegg"),
    extract_sig_ora(ora_reactome, "reactome")
)

# --- Merge ---
# We focus on the collections we want to plot: Hallmark, GO:BP, KEGG, Reactome
target_collections <- c("hallmark", "go_bp", "kegg", "reactome")

landscape_df <- gsea_df %>%
    filter(collection %in% target_collections) %>%
    # Join with ORA significance (Left Join: GSEA is the base)
    left_join(ora_sig_df, by = c("pathway" = "ID", "collection")) %>%
    mutate(
        is_ora_sig = replace_na(is_ora_sig, FALSE),
        # Clean collection names for plotting
        collection_label = case_when(
            collection == "hallmark" ~ "Hallmark",
            collection == "go_bp" ~ "GO: Biological Process",
            collection == "kegg" ~ "KEGG",
            collection == "reactome" ~ "Reactome"
        )
    )

cat(sprintf("  Integrated %d GSEA pathways (across 4 collections)\n", nrow(landscape_df)))

# ============================================================================
# 2. Select Top Pathways
# ============================================================================
cat("[2/4] Selecting top representative pathways...\n")

# Strategy: Top N significant pathways by padj, balanced by direction if possible
# We want ~10-15 per facet for readability

get_top_pathways <- function(df, n_top = 10) {
    df %>%
        filter(padj < 0.05) %>%
        group_by(collection) %>%
        arrange(padj) %>%
        slice_head(n = n_top) %>%
        ungroup()
}

plot_data <- get_top_pathways(landscape_df, n_top = 12)

# Order pathways by NES within each collection for the plot
plot_data <- plot_data %>%
    group_by(collection_label) %>%
    mutate(pathway_ordered = reorder(pathway, NES)) %>%
    ungroup()

# Clean pathway names (remove prefixes, replace underscores)
clean_names <- function(x) {
    x %>% 
        str_remove("^HALLMARK_") %>%
        str_remove("^GOBP_") %>%
        str_remove("^KEGG_") %>%
        str_remove("^REACTOME_") %>%
        str_replace_all("_", " ") %>%
        str_to_title() %>%
        str_trunc(45)
}

plot_data$pathway_clean <- clean_names(plot_data$pathway)
# Reorder factor based on NES
plot_data <- plot_data %>%
    arrange(collection_label, NES) %>%
    mutate(pathway_clean = factor(pathway_clean, levels = unique(pathway_clean)))

cat(sprintf("  Selected %d pathways for visualization\n", nrow(plot_data)))

# ============================================================================
# 3. Generate The Landscape Plot
# ============================================================================
cat("[3/4] Designing the visualization...\n")

p_landscape <- ggplot(plot_data, aes(x = NES, y = pathway_clean)) +
    # 0. Reference line
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
    
    # 1. Lollipop Stem
    geom_segment(aes(x = 0, xend = NES, y = pathway_clean, yend = pathway_clean), 
                 color = "grey40", linewidth = 0.5) +
    
    # 2. Lollipop Head (The Point)
    # Color by significance (gradient), Size by set size? Or constant?
    # Let's use Color = NES direction/magnitude (redundant but impactful) OR Significance
    # Let's stick to the Sanjana standard: Msgenta/Pink gradient for significance
    geom_point(aes(color = -log10(padj), size = size), alpha = 0.9) +
    
    # 3. ORA Validation Marker (The 'Halo' or distinct shape)
    # We add a second layer for ORA-confirmed hits
    geom_point(data = filter(plot_data, is_ora_sig), 
               aes(x = NES, y = pathway_clean), 
               shape = 21, color = "black", fill = NA, size = 6, stroke = 1.2) +
    
    # 4. Faceting
    facet_wrap(~collection_label, scales = "free_y", ncol = 2) +
    
    # 5. Scales & Theme
    scale_color_gradient(low = MAIN_COLOR_LIGHT, high = MAIN_COLOR_DARK, 
                        name = "-log10(FDR)") +
    scale_size_continuous(range = c(3, 8), name = "Gene Set Size") +
    
    # 6. Labels and Limits
    labs(
        title = "Integrated Pathway Enrichment Landscape",
        subtitle = "GSEA Normalized Enrichment Scores (NES) with ORA Confirmation",
        x = "Normalized Enrichment Score (NES)",
        y = "",
        caption = "Points encircled in black are independently validated by ORA (p < 0.1)\nSize = Number of genes in pathway | Color = GSEA Statistical Confidence"
    ) +
    
    PUB_THEME +
    theme(
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "grey95", color = NA),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 10, face = "bold", margin = margin(t=10)),
        panel.spacing.x = unit(1.5, "lines"), # More space between cols
        legend.position = "right",
        legend.box = "vertical"
    )

# ============================================================================
# 4. Save Output
# ============================================================================
cat("[4/4] Saving 'Ultimate' plot...\n")

output_file <- file.path(PLOT_DIR, "ultimate_enrichment_landscape.pdf")
# ggsave(output_file, p_landscape, width = 16, height = 12, dpi = 300)
# Use cowplot/grid save if needed, but since it's just ggplot facet, ggsave works fine.
ggsave(output_file, p_landscape, width = 16, height = 12, dpi = 300)

cat(sprintf("  SUCCESS: Saved to %s\n", basename(output_file)))

# Optional: Also save a CSV of the data used for the plot for transparency
write_csv(plot_data, file.path(PLOT_DIR, "ultimate_plot_source_data.csv"))
cat("  Saved source data CSV.\n")

cat("\n=== Done. ===\n")
