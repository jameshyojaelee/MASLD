# WGCNA Pipeline - Step 06: Core DEG Enrichment Analysis
# Tests if "Core DEGs" are significantly enriched in specific WGCNA modules.

library(WGCNA)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# ==============================================================================
# CONFIGURATION
# ==============================================================================
BASE_DIR <- "../results"
OUTPUT_DIR <- "../results/06_enrichment"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

CORE_DEGS_FILE <- "../../../final_core_degs.csv"
DATASETS <- c("Patient", "Other_MCD", "InHouse_MCD")

# ==============================================================================
# LOAD CORE DEGS
# ==============================================================================
cat("Loading Core DEGs from:", CORE_DEGS_FILE, "\n")
core_degs_df <- read.csv(CORE_DEGS_FILE)
# Assume 'Human_Genes' and 'Mouse_Genes' columns or similar.
# Let's inspect columns or assume standard 'Human_ID', 'Mouse_ID', 'Gene'
# We will create a list of core genes for Human and Mouse.

if ("Human_Name" %in% colnames(core_degs_df)) {
    core_genes_human <- unique(core_degs_df$Human_Name)
} else {
    warning("Could not find 'Human_Name' in core DEGs. Using 'Gene' if available.")
    core_genes_human <- unique(core_degs_df$Gene)
}

if ("Mouse_Symbol" %in% colnames(core_degs_df)) {
    core_genes_mouse <- unique(core_degs_df$Mouse_Symbol)
} else {
    core_genes_mouse <- unique(core_degs_df$Gene) # Fallback
}

cat("Found", length(core_genes_human), "Human core genes and", 
    length(core_genes_mouse), "Mouse core genes.\n")

# ==============================================================================
# ENRICHMENT ANALYSIS
# ==============================================================================

enrichment_results <- list()

for (ds in DATASETS) {
    cat("Analyzing", ds, "...\n")
    
    # Load Network Data
    net_file <- file.path(BASE_DIR, "03_network_construction", paste0("WGCNA_network_", ds, ".RData"))
    if (!file.exists(net_file)) {
        cat("  Network file missing for", ds, "- skipping.\n")
        next
    }
    load(net_file) # Loads moduleColors, moduleLabels, net
    
    # Get Genes in Network
    # We need expression data to match genes to colors
    load(file.path(BASE_DIR, "01_preprocessing", paste0(ds, ".RData")))
    if (ds == "Patient") {
       network_genes <- colnames(datExpr_patient)
       target_core <- core_genes_human
    } else {
       if (ds == "Other_MCD") network_genes <- colnames(datExpr_other)
       if (ds == "InHouse_MCD") network_genes <- colnames(datExpr_inhouse)
       target_core <- core_genes_mouse
    }
    
    # Standardize gene names (Upper case)
    network_genes_map <- toupper(network_genes)
    target_core_map <- toupper(target_core)
    
    # Identify Core Genes in Network
    is_core <- network_genes_map %in% target_core_map
    cat("  ", sum(is_core), "Core DEGs found in", ds, "network.\n")
    
    # Hypergeometric Test per Module
    modules <- unique(moduleColors)
    pvals <- numeric(length(modules))
    names(pvals) <- modules
    
    universe <- length(network_genes)
    m <- sum(is_core)         # total core genes in universe
    n <- universe - m         # total non-core
    
    for (mod in modules) {
        # k: size of module
        mod_genes_idx <- which(moduleColors == mod)
        k <- length(mod_genes_idx)
        
        # q: intersection (core genes in module)
        q <- sum(is_core[mod_genes_idx])
        
        # phyper(q-1, m, n, k, lower.tail=FALSE)
        pvals[mod] <- phyper(q - 1, m, n, k, lower.tail = FALSE)
    }
    
    enrichment_results[[ds]] <- pvals
    
    # Plot Enrichment Barplot
    df_plot <- data.frame(Module = names(pvals), Pval = -log10(pvals))
    df_plot <- df_plot %>% arrange(desc(Pval))
    
    p <- ggplot(df_plot, aes(x = reorder(Module, Pval), y = Pval, fill = Module)) +
        geom_bar(stat = "identity") +
        scale_fill_identity() +
        coord_flip() +
        geom_hline(yintercept = -log10(0.05), linetype="dashed", color="red") +
        labs(title = paste("Core DEG Enrichment in", ds), y = "-log10(P-value)", x = "Module") +
        theme_minimal()
    
    ggsave(file.path(OUTPUT_DIR, paste0("Enrichment_", ds, ".pdf")), p, width = 6, height = 8)
}

save(enrichment_results, file = file.path(OUTPUT_DIR, "enrichment_results.RData"))
cat("Enrichment analysis complete.\n")
