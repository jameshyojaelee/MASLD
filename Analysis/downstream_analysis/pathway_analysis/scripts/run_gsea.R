#!/usr/bin/env Rscript
# run_gsea.R - Gene Set Enrichment Analysis for all gene set collections (Integrative + Per-Dataset)
# Uses fgsea for fast preranked GSEA

# Critical: Exclude user library to avoid Matrix/libRlapack conflicts
.libPaths(.libPaths()[!grepl("jameslee", .libPaths())])

suppressPackageStartupMessages({
    library(tidyverse)
    library(fgsea)
    # library(clusterProfiler) # Missing in env
    # library(enrichplot)      # Missing in env
    library(msigdbr)
    library(org.Hs.eg.db)
    library(ggplot2)
    library(pheatmap)
})

# ============================================================================
# Configuration
# ============================================================================
PATHWAY_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/pathway_analysis"
DATA_DIR <- file.path(PATHWAY_DIR, "data/processed")
GENESET_DIR <- file.path(PATHWAY_DIR, "data/genesets")
OUTPUT_DIR <- file.path(PATHWAY_DIR, "results/gsea")
PLOT_DIR <- file.path(PATHWAY_DIR, "plots/individual")

# Ensure output directories exist
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Parameters
P_CUTOFF <- 0.05
MIN_SIZE <- 15
MAX_SIZE <- 500
N_PERM <- 10000

cat("=== Gene Set Enrichment Analysis (GSEA) ===\n\n")

# ============================================================================
# Load Gene Sets
# ============================================================================
cat("[1/6] Loading gene sets...\n")

# Helper to convert term2gene to list format for fgsea
t2g_to_list <- function(t2g_df) {
    split(t2g_df$gene_symbol, t2g_df$gs_name)
}

# Load all collections
genesets <- list(
    hallmark = t2g_to_list(read_csv(file.path(GENESET_DIR, "hallmark_term2gene.csv"), show_col_types = FALSE)),
    go_bp = t2g_to_list(read_csv(file.path(GENESET_DIR, "go_bp_term2gene.csv"), show_col_types = FALSE)),
    go_mf = t2g_to_list(read_csv(file.path(GENESET_DIR, "go_mf_term2gene.csv"), show_col_types = FALSE)),
    go_cc = t2g_to_list(read_csv(file.path(GENESET_DIR, "go_cc_term2gene.csv"), show_col_types = FALSE)),
    kegg = t2g_to_list(read_csv(file.path(GENESET_DIR, "kegg_term2gene.csv"), show_col_types = FALSE)),
    reactome = t2g_to_list(read_csv(file.path(GENESET_DIR, "reactome_term2gene.csv"), show_col_types = FALSE))
)

for (name in names(genesets)) {
    cat(sprintf("  %s: %d gene sets\n", name, length(genesets[[name]])))
}

# ============================================================================
# Run GSEA for All Dataset Files
# ============================================================================
cat("[2/6] Running GSEA for all ranked lists...\n")

# Find all .rnk files
rnk_files <- list.files(DATA_DIR, pattern = "\\.rnk$", full.names = TRUE)
cat(sprintf("  Found %d ranked list files: %s\n", length(rnk_files), paste(basename(rnk_files), collapse=", ")))

# Storage for results
all_gsea_runs <- list() # dataset -> collection -> result
integrative_results <- list() # collection -> result (Legacy support)
integrative_stats <- NULL

for (f in rnk_files) {
    # Determine dataset name
    fname <- basename(f)
    if (fname == "gsea_ranked_full.rnk") {
        dataset_name <- "Integrative"
    } else {
        dataset_name <- str_remove(fname, "gsea_ranked_") %>% str_remove("\\.rnk")
    }
    
    cat(sprintf("\n  --- Processing Dataset: %s ---\n", dataset_name))
    
    # Load Ranks
    ranked_df <- read_tsv(f, col_names = c("gene_symbol", "rank"), show_col_types = FALSE)
    stats <- ranked_df$rank
    names(stats) <- ranked_df$gene_symbol
    stats <- sort(stats, decreasing = TRUE)
    
    if (dataset_name == "Integrative") {
        integrative_stats <- stats
    }
    
    # Run for each collection
    dataset_results <- list()
    for (name in names(genesets)) {
        cat(sprintf("    Running %s...\n", name))
        
        result <- fgsea(
            pathways = genesets[[name]],
            stats = stats,
            minSize = MIN_SIZE,
            maxSize = MAX_SIZE,
            nPermSimple = N_PERM
        )
        
        result$collection <- name
        result$dataset <- dataset_name
        
        # Save individual CSV
        out_csv <- file.path(OUTPUT_DIR, paste0(dataset_name, "_", name, "_gsea.csv"))
        write_csv(result %>% as_tibble() %>% dplyr::select(-leadingEdge) %>% arrange(padj), out_csv)
        
        dataset_results[[name]] <- result
    }
    
    all_gsea_runs[[dataset_name]] <- dataset_results
    
    # Store "Integrative" results in the legacy format format
    if (dataset_name == "Integrative") {
        integrative_results <- dataset_results
    }
}

# Use Integrative results as the "main" gsea_results for downstream script compatibility
gsea_results <- integrative_results
stats <- integrative_stats

# ============================================================================
# Generate Dataset Comparison Heatmap (Hallmark)
# ============================================================================
cat("[3/6] Generating Cross-Dataset Comparison...\n")

# Combine Hallmark NES
hallmark_nes_df <- bind_rows(lapply(names(all_gsea_runs), function(dname) {
    res <- all_gsea_runs[[dname]]$hallmark
    as_tibble(res) %>% dplyr::select(pathway, NES, padj) %>% mutate(Dataset = dname)
}))

if (nrow(hallmark_nes_df) > 0) {
    # Pivot for heatmap
    nes_mat <- hallmark_nes_df %>%
        dplyr::select(pathway, Dataset, NES) %>%
        pivot_wider(names_from = Dataset, values_from = NES) %>%
        column_to_rownames("pathway") %>%
        as.matrix()
    
    # Filter for variability if too large
    if(nrow(nes_mat) > 2) {
        # Save heatmap
        pdf(file.path(PLOT_DIR, "gsea_hallmark_cross_dataset_heatmap.pdf"), width = 10, height = 12)
        pheatmap(nes_mat, 
                 main = "GSEA Hallmark NES Across Datasets",
                 color = colorRampPalette(c("blue", "white", "red"))(100),
                 display_numbers = FALSE,
                 fontsize_row = 8)
        dev.off()
        cat("    Saved: gsea_hallmark_cross_dataset_heatmap.pdf\n")
    }
}

# ============================================================================
# Generate Standard Plots (Integrative Only)
# ============================================================================
cat("[4/6] Generating Standard GSEA plots (Integrative)...\n")

# Load publication color theme
source("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/pathway_analysis/scripts/publication_theme.R")

# Function to create GSEA table plot with publication theme
create_gsea_table <- function(result, name, top_n = 20) {
    sig_results <- result %>%
        as_tibble() %>%
        filter(padj < P_CUTOFF) %>%
        arrange(padj) %>%
        head(top_n)
    
    if (nrow(sig_results) == 0) {
        return(NULL)
    }
    
    p <- ggplot(sig_results, aes(x = reorder(pathway, NES), y = NES, fill = NES)) +
        geom_col() +
        coord_flip() +
        scale_fill_gradient2(low = palette1[6], mid = "white", high = MAIN_COLOR, 
                            midpoint = 0, name = "NES") +
        labs(
            title = paste0("GSEA (Integrative): ", toupper(name)),
            x = "",
            y = "Normalized Enrichment Score (NES)"
        ) +
        PUB_THEME +
        theme(axis.text.y = element_text(size = 8))
    
    return(p)
}

# Generate plots for each collection (Integrative)
for (name in names(gsea_results)) {
    result <- gsea_results[[name]]
    
    # NES bar plot
    p <- create_gsea_table(result, name)
    if (!is.null(p)) {
        ggsave(file.path(PLOT_DIR, paste0(name, "_gsea_nes_barplot.pdf")), 
               p, width = 10, height = 8, dpi = 300)
    }
    
    # Classic enrichment curve (Top Up/Down)
    top_up <- result %>% filter(padj < P_CUTOFF & NES > 0) %>% arrange(padj) %>% head(2)
    top_down <- result %>% filter(padj < P_CUTOFF & NES < 0) %>% arrange(padj) %>% head(2)
    
    if (nrow(top_up) > 0 || nrow(top_down) > 0) {
        pdf(file.path(PLOT_DIR, paste0(name, "_gsea_enrichment_curves.pdf")), width = 10, height = 6)
        for (pw in c(top_up$pathway, top_down$pathway)) {
            print(plotEnrichment(genesets[[name]][[pw]], stats) +
                  labs(title = str_wrap(pw, width = 60)) + theme_minimal())
        }
        dev.off()
    }
}

# ============================================================================
# Create Summary Comparison (Integrative)
# ============================================================================
cat("[5/6] Creating GSEA summary comparison (Integrative)...\n")

all_gsea <- bind_rows(lapply(gsea_results, as_tibble)) %>%
    filter(padj < P_CUTOFF) %>%
    dplyr::select(pathway, NES, padj, collection) %>%
    arrange(padj)

write_csv(all_gsea, file.path(OUTPUT_DIR, "all_gsea_significant_summary.csv"))

top_per_collection <- all_gsea %>%
    group_by(collection) %>%
    slice_head(n = 10) %>%
    ungroup()

write_csv(top_per_collection, file.path(OUTPUT_DIR, "top10_per_collection_gsea.csv"))

# ============================================================================
# Save All Results
# ============================================================================
cat("[6/6] Saving results...\n")

# Save Integrative results (legacy name)
save(gsea_results, stats, genesets, 
     file = file.path(OUTPUT_DIR, "all_gsea_results.RData"))

# Save ALL results (new)
save(all_gsea_runs, file = file.path(OUTPUT_DIR, "multi_dataset_gsea_results.RData"))

cat("  Saved all_gsea_results.RData (Integrative)\n")
cat("  Saved multi_dataset_gsea_results.RData (All Datasets)\n")
cat("\n=== GSEA Analysis Complete ===\n")
