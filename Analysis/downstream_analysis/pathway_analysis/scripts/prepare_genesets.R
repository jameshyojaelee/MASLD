#!/usr/bin/env Rscript
# prepare_genesets.R - Extract subset gene sets from MSigDB GMT file
# Extracts Hallmark, GO, KEGG, and Reactome collections

suppressPackageStartupMessages({
    library(tidyverse)
    library(msigdbr)
})

# ============================================================================
# Configuration
# ============================================================================
PATHWAY_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/pathway_analysis"
GMT_FILE <- file.path(PATHWAY_DIR, "msigdb.v2025.1.Hs.symbols.gmt")
OUTPUT_DIR <- file.path(PATHWAY_DIR, "data/genesets")

cat("=== Extracting Gene Set Collections from MSigDB ===\n\n")

# ============================================================================
# Method 1: Use msigdbr package (primary method)
# ============================================================================
cat("[1/3] Loading gene sets via msigdbr...\n")

# Define collections to extract
collections <- list(
    hallmark = list(collection = "H", subcollection = NULL, name = "Hallmark"),
    go_bp = list(collection = "C5", subcollection = "GO:BP", name = "GO Biological Process"),
    go_mf = list(collection = "C5", subcollection = "GO:MF", name = "GO Molecular Function"),
    go_cc = list(collection = "C5", subcollection = "GO:CC", name = "GO Cellular Component"),
    kegg = list(collection = "C2", subcollection = "CP:KEGG_MEDICUS", name = "KEGG Medicus"),
    reactome = list(collection = "C2", subcollection = "CP:REACTOME", name = "Reactome")
)

# Function to convert msigdbr output to GMT format
write_gmt <- function(df, output_file) {
    # Group by gene set name
    gmt_lines <- df %>%
        group_by(gs_name) %>%
        summarize(
            genes = paste(gene_symbol, collapse = "\t"),
            .groups = "drop"
        ) %>%
        mutate(
            line = paste(gs_name, "NA", genes, sep = "\t")
        )
    
    writeLines(gmt_lines$line, output_file)
    return(nrow(gmt_lines))
}

# Extract each collection
for (coll_id in names(collections)) {
    coll <- collections[[coll_id]]
    cat(sprintf("  Extracting %s (Collection: %s", coll$name, coll$collection))
    if (!is.null(coll$subcollection)) {
        cat(sprintf(", Subcollection: %s", coll$subcollection))
    }
    cat(")...\n")

    # Get gene sets from msigdbr (v10+ API: collection/subcollection)
    if (is.null(coll$subcollection)) {
        gs_df <- msigdbr(species = "Homo sapiens", collection = coll$collection)
    } else {
        gs_df <- msigdbr(species = "Homo sapiens", collection = coll$collection, subcollection = coll$subcollection)
    }
    
    # Write to GMT
    output_file <- file.path(OUTPUT_DIR, paste0(coll_id, ".gmt"))
    n_sets <- write_gmt(gs_df, output_file)
    cat(sprintf("    -> Saved %d gene sets to %s\n", n_sets, basename(output_file)))
}

# ============================================================================
# Method 2: Parse GMT file directly (backup/validation)
# ============================================================================
cat("\n[2/3] Validating against original GMT file...\n")

# Read GMT file
gmt_lines <- readLines(GMT_FILE)
cat(sprintf("  Total gene sets in MSigDB: %d\n", length(gmt_lines)))

# Count by prefix
prefixes <- c(
    "HALLMARK_" = "Hallmark",
    "GOBP_" = "GO:BP",
    "GOMF_" = "GO:MF", 
    "GOCC_" = "GO:CC",
    "KEGG_" = "KEGG",
    "REACTOME_" = "Reactome"
)

for (prefix in names(prefixes)) {
    count <- sum(grepl(paste0("^", prefix), gmt_lines))
    cat(sprintf("  %s: %d gene sets\n", prefixes[prefix], count))
}

# ============================================================================
# Create TERM2GENE and TERM2NAME mappings for clusterProfiler
# ============================================================================
cat("\n[3/3] Creating TERM2GENE mappings for clusterProfiler...\n")

for (coll_id in names(collections)) {
    coll <- collections[[coll_id]]
    
    # Get gene sets from msigdbr (v10+ API)
    if (is.null(coll$subcollection)) {
        gs_df <- msigdbr(species = "Homo sapiens", collection = coll$collection)
    } else {
        gs_df <- msigdbr(species = "Homo sapiens", collection = coll$collection, subcollection = coll$subcollection)
    }
    
    # TERM2GENE: two columns (term, gene)
    term2gene <- gs_df %>%
        select(gs_name, gene_symbol) %>%
        distinct()
    
    # TERM2NAME: two columns (term, term name) for readable labels
    term2name <- gs_df %>%
        select(gs_name, gs_description) %>%
        distinct()
    
    # Save
    t2g_file <- file.path(OUTPUT_DIR, paste0(coll_id, "_term2gene.csv"))
    t2n_file <- file.path(OUTPUT_DIR, paste0(coll_id, "_term2name.csv"))
    
    write_csv(term2gene, t2g_file)
    write_csv(term2name, t2n_file)
    
    cat(sprintf("  %s: TERM2GENE (%d rows), TERM2NAME (%d terms)\n", 
                coll_id, nrow(term2gene), nrow(term2name)))
}

cat("\n=== Gene Set Preparation Complete ===\n")
cat(sprintf("\nOutput directory: %s\n", OUTPUT_DIR))
cat("Files created:\n")
cat("  - {collection}.gmt        : GMT format gene sets\n")
cat("  - {collection}_term2gene.csv : For clusterProfiler enricher()\n")
cat("  - {collection}_term2name.csv : Human-readable term names\n")
