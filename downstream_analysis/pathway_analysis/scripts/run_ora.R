#!/usr/bin/env Rscript
# run_ora.R - Over-Representation Analysis using fgsea::fora
# Replacement for clusterProfiler/ReactomePA version due to missing dependencies

# Critical: Exclude user library to avoid Matrix/libRlapack conflicts
.libPaths(.libPaths()[!grepl("jameslee", .libPaths())])

suppressPackageStartupMessages({
    library(tidyverse)
    library(fgsea) # Used for 'fora' function
    library(msigdbr)
    library(ggplot2)
    library(scales)
})

# ============================================================================
# Configuration
# ============================================================================
PATHWAY_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/downstream_analysis/pathway_analysis"
DATA_DIR <- file.path(PATHWAY_DIR, "data/processed")
GENESET_DIR <- file.path(PATHWAY_DIR, "data/genesets")
OUTPUT_DIR <- file.path(PATHWAY_DIR, "results/ora")
PLOT_DIR <- file.path(PATHWAY_DIR, "plots/individual")

# Ensure output directories exist
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Parameters
P_CUTOFF <- 0.05
Q_CUTOFF <- 0.1

cat("=== Over-Representation Analysis (ORA) via fgsea::fora ===\n\n")

# ============================================================================
# Load Gene List & Universe
# ============================================================================
cat("[1/6] Loading lists...\n")

gene_symbols <- readLines(file.path(DATA_DIR, "core_deg_symbols.txt"))
cat(sprintf("  Loaded %d DEG symbols (Query)\n", length(gene_symbols)))

background_universe <- readLines(file.path(DATA_DIR, "background_universe.txt"))
cat(sprintf("  Loaded %d background genes (Universe)\n", length(background_universe)))

# ============================================================================
# Load Gene Sets
# ============================================================================
cat("[2/6] Loading gene sets...\n")

# Helper to convert term2gene to list format for fgsea
t2g_to_list <- function(t2g_df) {
    split(t2g_df$gene_symbol, t2g_df$gs_name)
}

genesets <- list(
    hallmark = t2g_to_list(read_csv(file.path(GENESET_DIR, "hallmark_term2gene.csv"), show_col_types = FALSE)),
    go_bp = t2g_to_list(read_csv(file.path(GENESET_DIR, "go_bp_term2gene.csv"), show_col_types = FALSE)),
    go_mf = t2g_to_list(read_csv(file.path(GENESET_DIR, "go_mf_term2gene.csv"), show_col_types = FALSE)),
    go_cc = t2g_to_list(read_csv(file.path(GENESET_DIR, "go_cc_term2gene.csv"), show_col_types = FALSE)),
    kegg = t2g_to_list(read_csv(file.path(GENESET_DIR, "kegg_term2gene.csv"), show_col_types = FALSE)),
    reactome = t2g_to_list(read_csv(file.path(GENESET_DIR, "reactome_term2gene.csv"), show_col_types = FALSE))
)

# ============================================================================
# Run ORA (fora)
# ============================================================================
cat("[3/6] Running ORA...\n")

ora_results <- list()

for (name in names(genesets)) {
    cat(sprintf("  Running %s...\n", name))
    
    # fgsea::fora(pathways, genes, universe, minSize=1, maxSize=Inf)
    res <- fora(genesets[[name]], gene_symbols, background_universe, minSize = 10, maxSize = 500)
    
    # Add FDR/Padj if not present (fora returns pval and padj usually)
    # Rename columns to match expected output structure vaguely or standard style
    res <- res %>%
        as_tibble() %>%
        rename(ID = pathway, p.adjust = padj) %>%
        arrange(p.adjust) %>%
        mutate(
            Description = ID, # Placeholder if no description map
            GeneRatio = paste0(overlap, "/", length(gene_symbols)), # Approx ratio
            Count = overlap
        )
    
    ora_results[[name]] <- res
    
    # Save CSV
    write_csv(res, file.path(OUTPUT_DIR, paste0(name, "_ora_results.csv")))
}

# ============================================================================
# Generate Plots
# ============================================================================
cat("[4/6] Generating Plots...\n")
source("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/downstream_analysis/pathway_analysis/scripts/publication_theme.R")

plot_ora <- function(res, title, filename) {
    top_res <- res %>%
        filter(p.adjust < Q_CUTOFF) %>%
        head(20)
    
    if (nrow(top_res) == 0) return(NULL)
    
    p <- ggplot(top_res, aes(x = overlap, y = reorder(ID, overlap))) +
        geom_point(aes(color = -log10(p.adjust), size = overlap)) +
        scale_color_gradient(low = MAIN_COLOR_LIGHT, high = MAIN_COLOR_DARK) +
        labs(title = title, x = "Gene Count", y = "") +
        PUB_THEME +
        theme(axis.text.y = element_text(size = 8))
    
    ggsave(file.path(PLOT_DIR, filename), p, width = 10, height = 8)
}

plot_ora(ora_results$hallmark, "Hallmark ORA", "hallmark_ora_dotplot.pdf")
plot_ora(ora_results$go_bp, "GO:BP ORA", "go_bp_ora_dotplot.pdf")
plot_ora(ora_results$go_mf, "GO:MF ORA", "go_mf_ora_dotplot.pdf")
plot_ora(ora_results$go_cc, "GO:CC ORA", "go_cc_ora_dotplot.pdf")
plot_ora(ora_results$kegg, "KEGG ORA", "kegg_ora_dotplot.pdf")
plot_ora(ora_results$reactome, "Reactome ORA", "reactome_ora_dotplot.pdf")

# ============================================================================
# Save RData for Summary Script
# ============================================================================
cat("[5/6] Saving RData...\n")

# Map list items to variables expected by summary script
ora_hallmark <- ora_results$hallmark
ora_gobp <- ora_results$go_bp
ora_gomf <- ora_results$go_mf
ora_gocc <- ora_results$go_cc
ora_kegg <- ora_results$kegg
ora_reactome <- ora_results$reactome

# Note: These are now TIBBLES, not enrichResult objects. 
# Summary script must be updated to handle this!
save(ora_hallmark, ora_gobp, ora_gomf, ora_gocc, ora_kegg, ora_reactome,
     file = file.path(OUTPUT_DIR, "all_ora_results.RData"))

cat("Done.\n")
