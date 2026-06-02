#!/usr/bin/env Rscript
# meta_04_human_comparison.R
# Compare combined mouse MCD DEGs with human patient DEGs

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
})

root <- normalizePath(".")
out_dir <- file.path(root, "meta_analysis")
plot_dir <- file.path(root, "plots")

message("=== Cross-Species Comparison: Mouse MCD vs Human Patient ===\n")

# Load mouse gene sets
gene_sets <- readRDS(file.path(out_dir, "upregulated_gene_sets.rds"))
all_degs <- readRDS(file.path(out_dir, "all_degs_list.rds"))

# Load final_core_degs.csv for ortholog mapping
library_path <- file.path(root, "..", "..", "final_core_degs.csv")
if (!file.exists(library_path)) {
  library_path <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/results/library/final_core_degs.csv"
}

library_df <- read_csv(library_path, col_types = cols())
message("Loaded library with ", nrow(library_df), " genes")

# Create ortholog mapping (mouse -> human)
ortholog_map <- library_df %>%
  select(mouse_gene_id, human_ortholog_symbols) %>%
  filter(!is.na(human_ortholog_symbols), human_ortholog_symbols != "") %>%
  distinct()

message("Ortholog mappings available: ", nrow(ortholog_map))

# Load human patient DEG files
human_files <- tribble(
  ~dataset_id, ~path,
  "GSE130970_NAS", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/archive/old_results/target_selection/padj_0.1_upregulated/GSE130970_NAS_high_padj0.1_lfc0.75_upregulated.csv",
  "GSE135251_NAS", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/archive/old_results/target_selection/padj_0.1_upregulated/GSE135251_NAS_high_padj0.1_lfc1.20_upregulated.csv",
  "GSE130970_Fibrosis", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/archive/old_results/target_selection/padj_0.1_upregulated/GSE130970_Fibrosis_padj0.1_lfc0.25_upregulated.csv"
)

# Load human DEGs
human_degs <- list()
for (i in seq_len(nrow(human_files))) {
  path <- human_files$path[i]
  if (file.exists(path)) {
    df <- read_csv(path, col_types = cols())
    # Get gene symbols column
    if ("gene_symbol" %in% colnames(df)) {
      human_degs[[human_files$dataset_id[i]]] <- df$gene_symbol
    } else if ("gene_name" %in% colnames(df)) {
      human_degs[[human_files$dataset_id[i]]] <- df$gene_name
    } else {
      # Try first column
      human_degs[[human_files$dataset_id[i]]] <- df[[1]]
    }
    message("Loaded ", human_files$dataset_id[i], ": ", length(human_degs[[human_files$dataset_id[i]]]), " genes")
  } else {
    message("File not found: ", path)
  }
}

# Map mouse combined DEGs to human orthologs
combined_mouse <- gene_sets$Combined_lenient

# Get human orthologs for combined mouse DEGs
combined_with_orthologs <- tibble(mouse_gene_id = combined_mouse) %>%
  left_join(ortholog_map, by = "mouse_gene_id") %>%
  filter(!is.na(human_ortholog_symbols))

combined_human_orthologs <- unique(unlist(str_split(combined_with_orthologs$human_ortholog_symbols, ";")))
combined_human_orthologs <- combined_human_orthologs[!is.na(combined_human_orthologs) & combined_human_orthologs != ""]

message("\nCombined mouse upregulated DEGs: ", length(combined_mouse))
message("With human orthologs: ", nrow(combined_with_orthologs))
message("Unique human ortholog symbols: ", length(combined_human_orthologs))

# Also do for stringent
combined_mouse_stringent <- gene_sets$Combined_stringent
combined_stringent_with_orthologs <- tibble(mouse_gene_id = combined_mouse_stringent) %>%
  left_join(ortholog_map, by = "mouse_gene_id") %>%
  filter(!is.na(human_ortholog_symbols))
combined_stringent_human <- unique(unlist(str_split(combined_stringent_with_orthologs$human_ortholog_symbols, ";")))
combined_stringent_human <- combined_stringent_human[!is.na(combined_stringent_human) & combined_stringent_human != ""]

# Compute overlaps with human datasets
cross_species_summary <- list()
for (human_ds in names(human_degs)) {
  human_genes <- human_degs[[human_ds]]
  
  # Lenient
  overlap_lenient <- intersect(combined_human_orthologs, human_genes)
  
  # Stringent
  overlap_stringent <- intersect(combined_stringent_human, human_genes)
  
  cross_species_summary[[length(cross_species_summary) + 1]] <- tibble(
    comparison = paste0("Combined vs ", human_ds),
    mouse_degs = length(combined_human_orthologs),
    human_degs = length(human_genes),
    overlap_lenient = length(overlap_lenient),
    pct_mouse_lenient = round(100 * length(overlap_lenient) / length(combined_human_orthologs), 1),
    overlap_stringent = length(overlap_stringent),
    pct_human_stringent = round(100 * length(overlap_stringent) / length(human_genes), 1)
  )
  
  message(sprintf("\n%s: %d mouse orthologs overlap with %d human DEGs → %d concordant (%.1f%%)",
                  human_ds, length(combined_human_orthologs), length(human_genes), 
                  length(overlap_lenient), 100 * length(overlap_lenient) / length(combined_human_orthologs)))
}

cross_species_df <- bind_rows(cross_species_summary)
write_csv(cross_species_df, file.path(out_dir, "cross_species_overlap.csv"))

# Create intersection summary for cross-species comparison
if (length(human_degs) > 0) {
  # Use lenient mouse orthologs
  upset_list <- c(
    list(Combined_Mouse = combined_human_orthologs),
    human_degs
  )
  
  # Calculate key intersections
  all_intersect <- Reduce(intersect, upset_list)
  mouse_only <- setdiff(combined_human_orthologs, unique(unlist(human_degs)))
  
  cross_intersect <- tibble(
    intersection = c("All datasets", "Mouse only", "Human NAS overlap"),
    count = c(
      length(all_intersect),
      length(mouse_only),
      length(intersect(combined_human_orthologs, unique(c(human_degs$GSE130970_NAS, human_degs$GSE135251_NAS))))
    )
  )
  
  write_csv(cross_intersect, file.path(out_dir, "cross_species_intersections.csv"))
  message("\nCross-species intersections:")
  print(cross_intersect)
}

# Save concordant genes (in both mouse combined and human NAS)
all_human_nas <- unique(c(human_degs$GSE130970_NAS, human_degs$GSE135251_NAS))
concordant <- intersect(combined_human_orthologs, all_human_nas)
write_csv(tibble(gene_symbol = concordant), file.path(out_dir, "concordant_genes_human_mouse.csv"))
message("\nConcordant genes (mouse + human NAS): ", length(concordant))

message("\n=== Cross-Species Summary ===")
print(cross_species_df)
