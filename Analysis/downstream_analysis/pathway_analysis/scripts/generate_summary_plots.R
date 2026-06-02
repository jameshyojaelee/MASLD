#!/usr/bin/env Rscript
# generate_summary_plots.R - Create comprehensive summary visualizations
# Combines results from ORA and GSEA across all collections

# Critical: Exclude user library to avoid Matrix/libRlapack conflicts
.libPaths(.libPaths()[!grepl("jameslee", .libPaths())])

suppressPackageStartupMessages({
    library(tidyverse)
    library(pheatmap)
    library(RColorBrewer)
    # library(ggupset) # Missing
    library(cowplot)   # Replacement for patchwork
    library(scales)
})

# ============================================================================
# Configuration
# ============================================================================
PATHWAY_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/pathway_analysis"
ORA_DIR <- file.path(PATHWAY_DIR, "results/ora")
GSEA_DIR <- file.path(PATHWAY_DIR, "results/gsea")
PLOT_DIR <- file.path(PATHWAY_DIR, "plots/summary")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load publication color theme
source("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/pathway_analysis/scripts/publication_theme.R")

cat("=== Generating Summary Visualizations ===\n\n")

# ============================================================================
# Load Results
# ============================================================================
cat("[1/6] Loading analysis results...\n")

# Load ORA results
load(file.path(ORA_DIR, "all_ora_results.RData"))
cat("  Loaded ORA results\n")

# Load GSEA results  
load(file.path(GSEA_DIR, "all_gsea_results.RData"))
cat("  Loaded GSEA results\n")

# ============================================================================
# Plot 1: Multi-Panel Dot Plot Comparison
# ============================================================================
cat("[2/6] Creating multi-panel summary...\n")

# Load actual background universe size for correct GeneRatio calculation
UNIVERSE_FILE <- file.path(PATHWAY_DIR, "data/processed/background_universe.txt")
if (file.exists(UNIVERSE_FILE)) {
    UNIVERSE_SIZE <- length(readLines(UNIVERSE_FILE))
} else {
    UNIVERSE_SIZE <- 17000  # Fallback estimate
    warning("Could not find background_universe.txt, using fallback size")
}
cat(sprintf("  Background universe size: %d genes\n", UNIVERSE_SIZE))

# Helper function to extract ORA data for plotting
extract_ora_top <- function(ora_result, name, top_n = 10) {
    if (is.null(ora_result) || nrow(ora_result) == 0) {
        return(NULL)
    }
    
    # Input is now a tibble from fgsea::fora
    ora_result %>%
        filter(p.adjust < 0.1) %>%
        arrange(p.adjust) %>%
        head(top_n) %>%
        mutate(
            collection = name,
            Description = str_wrap(ID, width = 40), # ID is the name in fora outputs
            # Use actual universe size for GeneRatio
            GeneRatio_numeric = overlap / UNIVERSE_SIZE,
            Count = overlap
        )
}


# Extract top results from each collection
ora_combined <- bind_rows(
    extract_ora_top(ora_hallmark, "Hallmark"),
    extract_ora_top(ora_gobp, "GO:BP"),
    extract_ora_top(ora_gomf, "GO:MF"),
    extract_ora_top(ora_gocc, "GO:CC"),
    extract_ora_top(ora_kegg, "KEGG"),
    extract_ora_top(ora_reactome, "Reactome")
)

if (nrow(ora_combined) > 0) {
    # Create faceted dot plot with magenta theme
    p_multipanel <- ggplot(ora_combined, 
                           aes(x = GeneRatio_numeric, y = reorder(Description, GeneRatio_numeric))) +
        geom_point(aes(size = Count, color = -log10(p.adjust))) +
        facet_wrap(~collection, scales = "free_y", ncol = 2) +
        scale_color_gradient(low = MAIN_COLOR_LIGHT, high = MAIN_COLOR_DARK, name = "-log10(padj)") +
        scale_size_continuous(range = c(2, 8), name = "Gene Count") +
        labs(
            title = "Pathway Enrichment Summary: ORA Results",
            x = "Gene Ratio",
            y = ""
        ) +
        PUB_THEME +
        theme(
            strip.text = element_text(face = "bold", size = 12),
            axis.text.y = element_text(size = 7)
        )
    
    ggsave(file.path(PLOT_DIR, "pathway_summary_multipanel_ora.pdf"), 
           p_multipanel, width = 16, height = 14, dpi = 300)
    cat("  Saved: pathway_summary_multipanel_ora.pdf\n")
}

# ============================================================================
# Plot 2: GSEA NES Comparison Heatmap
# ============================================================================
cat("[3/6] Creating GSEA heatmap...\n")

# Extract top GSEA pathways per collection
gsea_top <- lapply(names(gsea_results), function(name) {
    gsea_results[[name]] %>%
        filter(padj < 0.05) %>%
        arrange(padj) %>%
        head(8) %>%
        mutate(collection = name) %>%
        select(pathway, NES, padj, collection)
}) %>% bind_rows()

if (nrow(gsea_top) > 0) {
    # Create heatmap-style summary with magenta/blue theme
    p_gsea_summary <- ggplot(gsea_top, 
                             aes(x = collection, y = reorder(pathway, NES), fill = NES)) +
        geom_tile(color = "white", size = 0.5) +
        geom_text(aes(label = sprintf("%.2f", NES)), size = 2.5) +
        scale_fill_gradient2(low = palette1[6], mid = "white", high = MAIN_COLOR, 
                             midpoint = 0, name = "NES") +
        labs(
            title = "Top GSEA Pathways by Collection",
            x = "Gene Set Collection",
            y = ""
        ) +
        PUB_THEME +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
            axis.text.y = element_text(size = 6),
            panel.grid = element_blank()
        )
    
    ggsave(file.path(PLOT_DIR, "gsea_nes_heatmap.pdf"), 
           p_gsea_summary, width = 12, height = 14, dpi = 300)
    cat("  Saved: gsea_nes_heatmap.pdf\n")
}

# ============================================================================
# Plot 3: Pathway Count Summary
# ============================================================================
cat("[4/6] Creating pathway count summary...\n")

# Count significant pathways per collection (ORA)
# Input is now a tibble, so we check nrow(filter(...)) directly
ora_counts <- tibble(
    collection = c("Hallmark", "GO:BP", "GO:MF", "GO:CC", "KEGG", "Reactome"),
    ora_count = c(
        if(!is.null(ora_hallmark)) sum(ora_hallmark$p.adjust < 0.1) else 0,
        if(!is.null(ora_gobp)) sum(ora_gobp$p.adjust < 0.1) else 0,
        if(!is.null(ora_gomf)) sum(ora_gomf$p.adjust < 0.1) else 0,
        if(!is.null(ora_gocc)) sum(ora_gocc$p.adjust < 0.1) else 0,
        if(!is.null(ora_kegg)) sum(ora_kegg$p.adjust < 0.1) else 0,
        if(!is.null(ora_reactome)) sum(ora_reactome$p.adjust < 0.1) else 0
    )
)

# Count significant pathways per collection (GSEA)
gsea_counts <- tibble(
    collection = names(gsea_results),
    gsea_count = sapply(gsea_results, function(x) sum(x$padj < 0.05, na.rm = TRUE))
) %>%
    mutate(collection = case_when(
        collection == "hallmark" ~ "Hallmark",
        collection == "go_bp" ~ "GO:BP",
        collection == "go_mf" ~ "GO:MF",
        collection == "go_cc" ~ "GO:CC",
        collection == "kegg" ~ "KEGG",
        collection == "reactome" ~ "Reactome",
        TRUE ~ collection
    ))

# Merge and reshape
counts_df <- left_join(ora_counts, gsea_counts, by = "collection") %>%
    pivot_longer(cols = c(ora_count, gsea_count), 
                 names_to = "method", values_to = "count") %>%
    mutate(method = ifelse(method == "ora_count", "ORA", "GSEA"))

p_counts <- ggplot(counts_df, aes(x = collection, y = count, fill = method)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("ORA" = MAIN_COLOR, "GSEA" = palette1[6])) +
    labs(
        title = "Significant Pathways by Method and Collection",
        x = "Gene Set Collection",
        y = "Number of Significant Pathways",
        fill = "Method"
    ) +
    PUB_THEME +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

ggsave(file.path(PLOT_DIR, "pathway_counts_comparison.pdf"), 
       p_counts, width = 10, height = 6, dpi = 300)
cat("  Saved: pathway_counts_comparison.pdf\n")

# ============================================================================
# Plot 4: Top 5 from Each Collection - Combined
# ============================================================================
cat("[5/6] Creating top pathways summary...\n")

# Get top 5 from each collection for ORA
ora_top5 <- ora_combined %>%
    group_by(collection) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    mutate(
        label = paste0(collection, ": ", str_trunc(Description, 50)),
        neg_log_p = -log10(p.adjust)
    )

if (nrow(ora_top5) > 0) {
    p_top5 <- ggplot(ora_top5, 
                     aes(x = neg_log_p, y = reorder(label, neg_log_p), fill = collection)) +
        geom_col() +
        scale_fill_manual(values = COLLECTION_COLORS) +
        labs(
            title = "Top 5 Pathways from Each Collection (ORA)",
            x = "-log10(Adjusted P-value)",
            y = "",
            fill = "Collection"
        ) +
        PUB_THEME +
        theme(
            axis.text.y = element_text(size = 7)
        )
    
    ggsave(file.path(PLOT_DIR, "top5_per_collection_ora.pdf"), 
           p_top5, width = 14, height = 10, dpi = 300)
    cat("  Saved: top5_per_collection_ora.pdf\n")
}

# ============================================================================
# Plot 5: Save Summary Statistics
# ============================================================================
cat("[6/6] Saving summary statistics...\n")

# Create summary table
summary_stats <- counts_df %>%
    pivot_wider(names_from = method, values_from = count) %>%
    mutate(total = ORA + GSEA) %>%
    arrange(desc(total))

write_csv(summary_stats, file.path(PLOT_DIR, "pathway_analysis_summary.csv"))
cat("  Saved: pathway_analysis_summary.csv\n")

# Print summary
cat("\n=== Analysis Summary ===\n")
print(summary_stats)

cat("\n=== Summary Visualization Complete ===\n")
cat(sprintf("All summary plots saved to: %s\n", PLOT_DIR))
