#!/usr/bin/env Rscript
# run_refactored_ssgsea.R - ssGSEA on Normalized Counts
#
# This script implements the refactored workflow:
# 1. Load DESeq2 Normalized Counts (Fundamental Data)
# 2. Map IDs (Ensembl -> Symbol)
# 3. Run ssGSEA (Rank-based single-sample enrichment)
# 4. DoRothEA (TF Activity)

# Critical: Exclude user library (prevents R version mismatch)
# We handle library paths via micromamba environment, but explicitly setting it ensures safety.
conda_lib <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/downstream_analysis/pathway_analysis/.mamba/pathway_analysis/lib/R/library"
if (dir.exists(conda_lib)) {
    .libPaths(c(conda_lib))
}

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(GSVA)
    library(decoupleR)
    library(OmnipathR)
    library(ComplexHeatmap)
    library(circlize)
    library(msigdbr)
    # Load AnnotationDbi last or explicitly use ::select to avoid conflicts
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
    library(AnnotationDbi)
})

# ============================================================================
# Configuration
# ============================================================================
BASE_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PATHWAY_DIR <- file.path(BASE_DIR, "downstream_analysis/pathway_analysis")
RNASEQ_DIR <- file.path(BASE_DIR, "RNA-seq")
OUTPUT_DIR <- file.path(PATHWAY_DIR, "results/refactored")
PLOT_DIR <- file.path(PATHWAY_DIR, "plots/refactored")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load theme if available
theme_path <- file.path(PATHWAY_DIR, "scripts/publication_theme.R")
if (file.exists(theme_path)) {
    source(theme_path)
} else {
    MAIN_COLOR <- "red" # Fallback
}

cat("=== Refactored ssGSEA & DoRothEA Workflow ===\n\n")

# ============================================================================
# Helpers
# ============================================================================

# Wrapper for ssGSEA
run_ssgsea_wrapper <- function(mat, geneset_list) {
    # ssGSEA on counts (rank based)
    param <- ssgseaParam(exprData = mat, geneSets = geneset_list)
    gsva(param, verbose = FALSE)
}

# Helper to map Ensembl to Symbol and aggregate max expression
process_counts_smart <- function(file_path, species = "human") {
    cat(sprintf("Processing %s...\n", basename(file_path)))
    if (!file.exists(file_path)) return(NULL)
    
    counts <- fread(file_path)
    
    # Handle headers
    if (colnames(counts)[1] == "V1" || colnames(counts)[1] == "") colnames(counts)[1] <- "gene_id"
    if (!"gene_id" %in% colnames(counts)) colnames(counts)[1] <- "gene_id" 
    
    ids <- counts$gene_id
    
    # Check if already symbols (no ENS prefix) or if a symbol column exists
    is_ensembl <- any(grepl("^ENS", ids[1:10]))
    
    # Check for common symbol column names
    symbol_col <- NULL
    possible_names <- c("gene_name", "gene_symbol", "symbol", "Symbol", "SYMBOL")
    for (p in possible_names) {
        if (p %in% colnames(counts)) {
            symbol_col <- p
            break
        }
    }
    
    if (!is.null(symbol_col)) {
        cat(sprintf("  Using existing symbol column: %s\n", symbol_col))
        counts$Symbol <- counts[[symbol_col]]
    } else if (is_ensembl) {
        cat("  Mapping Ensembl IDs to Symbols...\n")
        # Remove version
        counts$gene_id <- sub("\\..*", "", counts$gene_id)
        
        db <- if (species == "human") org.Hs.eg.db else org.Mm.eg.db
        # Map using vector
        mapping <- mapIds(db, keys = counts$gene_id, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
        counts$Symbol <- mapping
    } else {
        cat("  Using gene_id as Symbol (assuming already symbols)...\n")
        counts$Symbol <- counts$gene_id
    }
    
    # Aggregation
    mat_clean <- counts %>%
        as_tibble() %>%
        filter(!is.na(Symbol) & Symbol != "") %>%
        dplyr::select(Symbol, everything(), -gene_id) %>%
        # Fix: Force numeric conversion on sample columns to prevent pivot errors
        mutate(across(-Symbol, as.numeric)) %>%
        # Fix: Remove columns that became all NA (e.g. gene_name, gene_type)
        dplyr::select(Symbol, where(~ any(!is.na(.)))) %>%
        pivot_longer(-Symbol, names_to = "Sample", values_to = "Count") %>%
        group_by(Symbol, Sample) %>%
        summarise(Count = max(Count), .groups = "drop") %>%
        pivot_wider(names_from = Sample, values_from = Count)
    
    mat <- as.matrix(mat_clean %>% dplyr::select(-Symbol))
    rownames(mat) <- mat_clean$Symbol
    
    # Log-transform to stabilize variance (critical for DoRothEA, neutral for ssGSEA ranks)
    mat <- log1p(mat)
    
    return(mat)
}

# ============================================================================
# 1. Load Data
# ============================================================================
cat("[1/4] Loading Datasets...\n")

files <- list(
    "Hoang" = file.path(RNASEQ_DIR, "patient_RNAseq/results/GSE130970/deseq2_results/20260121_094851/deseq2_normalized_counts.csv"),
    "Govaere" = file.path(RNASEQ_DIR, "patient_RNAseq/results/GSE135251/deseq2_results/20260121_094851/deseq2_normalized_counts.csv"),
    "InHouse_MCD" = file.path(RNASEQ_DIR, "in-house_MCD_RNAseq/normalized_counts_all_samples.csv"),
    "Ext_MCD_1" = file.path(RNASEQ_DIR, "other_MCD_RNAseq/GSE156918/analysis_mcd_vs_control/normalized_counts.csv"),
    "Ext_MCD_2" = file.path(RNASEQ_DIR, "other_MCD_RNAseq/GSE205974/analysis_mcd_vs_control/normalized_counts.csv")
)

mats <- list()
mats[["Hoang"]] <- process_counts_smart(files[["Hoang"]], "human")
mats[["Govaere"]] <- process_counts_smart(files[["Govaere"]], "human")
mats[["InHouse_MCD"]] <- process_counts_smart(files[["InHouse_MCD"]], "mouse")
mats[["Ext_MCD_1"]] <- process_counts_smart(files[["Ext_MCD_1"]], "mouse")
mats[["Ext_MCD_2"]] <- process_counts_smart(files[["Ext_MCD_2"]], "mouse")

# ============================================================================
# 1.5. Gene Universe Info (Diagnostic Only - No Subsetting)
# ============================================================================
cat("[1.5/4] Analyzing Gene Universes (Diagnostic Only)...\n")

# Report gene counts and overlap for information only
# NOTE: We do NOT subset to common genes. ssGSEA is a rank-based method that
# normalizes within each sample. Subsetting would:
#   1. Remove informative genes that may be pathway members
#   2. Artificially constrain analysis to a reduced gene space
#   3. Not improve comparability (scores are already sample-normalized)

for (n in names(mats)) {
    if (!is.null(mats[[n]])) {
        cat(sprintf("  %s: %d genes\n", n, nrow(mats[[n]])))
    }
}

# Calculate intersection for diagnostic purposes only
gene_lists <- lapply(names(mats), function(n) {
    if (is.null(mats[[n]])) return(NULL)
    r <- rownames(mats[[n]])
    # Convert to upper for intersection calculation only
    if (n %in% c("InHouse_MCD", "Ext_MCD_1", "Ext_MCD_2")) return(toupper(r))
    return(r)
})

common_genes <- Reduce(intersect, gene_lists)
cat(sprintf("  Common intersection (for reference): %d genes\n", length(common_genes)))
cat("  NOTE: Using FULL transcriptome per dataset (no subsetting)\n")

# ============================================================================
# 1.6. Load Hallmark Gene Sets (Species-Specific)
# ============================================================================
hallmark_human_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_mouse_df <- msigdbr(species = "Mus musculus", collection = "H")

# Convert to lists
hallmark_human_list <- split(hallmark_human_df$gene_symbol, hallmark_human_df$gs_name)
hallmark_mouse_list <- split(hallmark_mouse_df$gene_symbol, hallmark_mouse_df$gs_name)

# ============================================================================
# 2. Run ssGSEA
# ============================================================================
cat("[2/4] Running ssGSEA...\n")

ssgsea_out_file <- file.path(OUTPUT_DIR, "all_ssgsea_scores.rds")

ssgsea_list <- list()
for (n in names(mats)) {
    if (!is.null(mats[[n]])) {
        cat(sprintf("  Running %s (%d samples)...\n", n, ncol(mats[[n]])))
        
        mat_for_gsea <- mats[[n]]
        
        # Determine species and select appropriate gene set + fix rownames
        if (n %in% c("InHouse_MCD", "Ext_MCD_1", "Ext_MCD_2")) {
             # Mouse: Use Mouse Hallmark List
             # Ensure rownames are Title Case (standard Mouse Symbol)
             # process_counts_smart returns standard symbols from OrgDb, which are usually Title Case for Mouse.
             # No toupper() needed.
             target_geneset <- hallmark_mouse_list
             cat("    Species: Mouse (using Mus musculus Hallmark)\n")
        } else {
             # Human: Use Human Hallmark List
             target_geneset <- hallmark_human_list
             cat("    Species: Human (using Homo sapiens Hallmark)\n")
        }
        
        res <- run_ssgsea_wrapper(mat_for_gsea, target_geneset)
        ssgsea_list[[n]] <- res
    }
}

# Combine into long format
ssgsea_df <- bind_rows(lapply(names(ssgsea_list), function(n) {
    res <- ssgsea_list[[n]]
    as_tibble(res, rownames = "Pathway") %>%
        pivot_longer(-Pathway, names_to = "Sample", values_to = "Score") %>%
        mutate(Dataset = n)
}))

# Infer Condition (Heuristic)

# Assign Condition based on dataset type and metadata
# We separate the data frame by dataset to apply specific logic
ssgsea_list_df <- list()

for (n in unique(ssgsea_df$Dataset)) {
    df_sub <- ssgsea_df %>% filter(Dataset == n)
    
    if (n == "Hoang") {
        # Load metadata for GSE130970
        meta_file <- file.path(BASE_DIR, "RNA-seq/patient_RNAseq/data/samplesheets/GSE130970_samplesheet.csv")
        if (file.exists(meta_file)) {
            meta <- read_csv(meta_file, show_col_types = FALSE)
            # Map based on fibrosis stage (0 = Control, >=1 = Disease)
            meta <- meta %>%
                mutate(Cond = if_else(fibrosis_stage == 0, "Control", "Disease")) %>%
                dplyr::select(sample, Cond)
            
            df_sub <- df_sub %>%
                left_join(meta, by = c("Sample" = "sample")) %>%
                mutate(Condition = coalesce(Cond, "Unknown")) %>%
                dplyr::select(-Cond)
        } else {
            df_sub <- df_sub %>% mutate(Condition = "Unknown")
        }
        
    } else if (n == "Govaere") {
        # Load metadata for GSE135251
        meta_file <- file.path(BASE_DIR, "RNA-seq/patient_RNAseq/data/samplesheets/GSE135251_samplesheet.csv")
        if (file.exists(meta_file)) {
            meta <- read_csv(meta_file, show_col_types = FALSE)
            # Map based on group_in_paper
            meta <- meta %>%
                mutate(Cond = if_else(grepl("control", group_in_paper, ignore.case=TRUE), "Control", "Disease")) %>%
                dplyr::select(sample, Cond)
            
            df_sub <- df_sub %>%
                left_join(meta, by = c("Sample" = "sample")) %>%
                mutate(Condition = coalesce(Cond, "Unknown")) %>%
                dplyr::select(-Cond)
        } else {
            df_sub <- df_sub %>% mutate(Condition = "Unknown")
        }
        
    } else {
        # Fallback / Mouse Regex Heuristic
        df_sub <- df_sub %>%
            mutate(
                Condition = case_when(
                    grepl("MCD", Sample, ignore.case = TRUE) ~ "Disease",
                    grepl("Chow|Control|Ctrl", Sample, ignore.case = TRUE) ~ "Control",
                    grepl("Fibrosis|F3|F4|NASH", Sample, ignore.case = TRUE) ~ "Disease",
                    grepl("(CF|CM)$", Sample) ~ "Control", 
                    grepl("(MF|MM)$", Sample) ~ "Disease",
                    TRUE ~ "Unknown" 
                )
            )
    }
    ssgsea_list_df[[n]] <- df_sub
}

ssgsea_df <- bind_rows(ssgsea_list_df)


saveRDS(ssgsea_df, ssgsea_out_file)

# ============================================================================
# 3. Viz
# ============================================================================
cat("[3/4] Visualization...\n")

target_pathways <- c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_FATTY_ACID_METABOLISM")

p_box <- ssgsea_df %>%
    filter(Pathway %in% target_pathways) %>%
    # filter(Condition != "Unknown") %>% # Show all
    ggplot(aes(x = Dataset, y = Score, fill = Dataset)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
    facet_wrap(~Pathway, scales = "fixed") +
    # scale_fill_manual(values = c("Control" = "grey70", "Disease" = MAIN_COLOR, "Unknown" = "grey90")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ssGSEA Scores (Unstratified)")

ggsave(file.path(PLOT_DIR, "ssgsea_validation_boxplot.pdf"), p_box, width = 12, height = 6)

# ============================================================================
# 4. DoRothEA
# ============================================================================
cat("[4/4] DoRothEA Integration...\n")

run_dorothea <- function(mat, species) {
    organism <- if (species == "human") "human" else "mouse"
    tryCatch({
        net <- get_dorothea(organism = organism, levels = c("A", "B", "C"))
        acts <- run_wmean(mat, net, .source = "source", .target = "target", .mor = "mor", times = 100, minsize = 5)
        return(acts)
    }, error = function(e) {
        warning(sprintf("DoRothEA failed for %s: %s", species, e$message))
        return(NULL)
    })
}

dorothea_list <- list()
for (n in names(mats)) {
    if (!is.null(mats[[n]])) {
        cat(sprintf("  Running DoRothEA for %s...\n", n))
        spec <- if (n %in% c("Hoang", "Govaere")) "human" else "mouse" 
        res <- run_dorothea(mats[[n]], spec)
        dorothea_list[[n]] <- res
    }
}

saveRDS(dorothea_list, file.path(OUTPUT_DIR, "all_dorothea_scores.rds"))

cat("Done.\n")
