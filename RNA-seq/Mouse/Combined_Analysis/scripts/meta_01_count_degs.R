#!/usr/bin/env Rscript
# meta_01_count_degs.R
# Count upregulated DEGs across all MCD datasets at different thresholds

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
})

root <- normalizePath(".")
out_dir <- file.path(root, "meta_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("=== Counting Upregulated DEGs Across Datasets ===\n")

# Define datasets and their DEG result files
datasets <- tribble(
  ~dataset_id, ~dataset_label, ~path,
  "Combined", "Combined (batch-corrected)", 
    file.path(root, "analysis_mcd_vs_control", "deseq2_mcd_vs_control.tsv"),
  "InHouse", "In-house pooled", 
    file.path(root, "..", "InHouse_MCD", "results", "combined_pooled", "deseq2_results.tsv"),
  "GSE156918", "GSE156918 (external)", 
    file.path(root, "..", "Public_MCD", "GSE156918", "analysis_mcd_vs_control", "deseq2_mcd_vs_control.tsv"),
  "GSE205974", "GSE205974 (external)", 
    file.path(root, "..", "Public_MCD", "GSE205974", "analysis_mcd_vs_control", "deseq2_mcd_vs_control.tsv")
)

# Load all DEG results
load_degs <- function(path, dataset_id) {
  if (!file.exists(path)) {
    warning("File not found: ", path)
    return(NULL)
  }
  df <- read_tsv(path, col_types = cols())
  
  # Standardize column names
  if ("gene" %in% colnames(df)) {
    df <- df %>% rename(gene_id = gene)
  }
  if (!"gene_id" %in% colnames(df)) {
    # Try to find gene column
    gene_col <- colnames(df)[1]
    df <- df %>% rename(gene_id = !!sym(gene_col))
  }
  
  df %>%
    mutate(
      dataset = dataset_id,
      # Strip version suffix from Ensembl ID for consistent matching
      gene_id = str_replace(gene_id, "\\.\\d+$", "")
    ) %>%
    select(dataset, gene_id, any_of(c("gene_name", "baseMean", "log2FoldChange", "padj")))
}

all_degs <- list()
for (i in seq_len(nrow(datasets))) {
  result <- load_degs(datasets$path[i], datasets$dataset_id[i])
  if (!is.null(result)) {
    all_degs[[datasets$dataset_id[i]]] <- result
    message("Loaded ", datasets$dataset_id[i], ": ", nrow(result), " genes")
  }
}

# Function to count upregulated DEGs at different thresholds
count_upregulated <- function(df, padj_cutoff = 0.1, lfc_cutoff = 0, tpm_proxy_cutoff = 0) {
  # Note: Using baseMean as proxy for expression level (not TPM, but correlated)
  df %>%
    filter(
      !is.na(padj), padj < padj_cutoff,
      !is.na(log2FoldChange), log2FoldChange > lfc_cutoff,
      !is.na(baseMean), baseMean > tpm_proxy_cutoff
    ) %>%
    nrow()
}

# Get upregulated gene lists
get_upregulated <- function(df, padj_cutoff = 0.1, lfc_cutoff = 0, basemean_cutoff = 0) {
  df %>%
    filter(
      !is.na(padj), padj < padj_cutoff,
      !is.na(log2FoldChange), log2FoldChange > lfc_cutoff,
      !is.na(baseMean), baseMean > basemean_cutoff
    ) %>%
    pull(gene_id)
}

# Count DEGs at different thresholds
thresholds <- tribble(
  ~threshold_id, ~padj, ~lfc, ~basemean, ~description,
  "lenient", 0.1, 0, 0, "padj<0.1, LFC>0",
  "stringent", 0.1, 0.58, 10, "padj<0.1, LFC>=0.58, baseMean>10"
)

# Create summary table
summary_list <- list()
for (ds_id in names(all_degs)) {
  df <- all_degs[[ds_id]]
  for (j in seq_len(nrow(thresholds))) {
    th <- thresholds[j, ]
    n <- count_upregulated(df, th$padj, th$lfc, th$basemean)
    summary_list[[length(summary_list) + 1]] <- tibble(
      dataset = ds_id,
      threshold = th$threshold_id,
      n_upregulated = n
    )
  }
}

summary_df <- bind_rows(summary_list) %>%
  pivot_wider(names_from = threshold, values_from = n_upregulated)

# Add labels
summary_df <- summary_df %>%
  left_join(datasets %>% select(dataset_id, dataset_label), by = c("dataset" = "dataset_id"))

message("\n=== Upregulated DEG Counts ===")
print(summary_df)

# Save summary
write_csv(summary_df, file.path(out_dir, "deg_counts_summary.csv"))
message("\nSaved: ", file.path(out_dir, "deg_counts_summary.csv"))

# Extract upregulated gene sets for overlap analysis
gene_sets <- list()
for (ds_id in names(all_degs)) {
  df <- all_degs[[ds_id]]
  
  # Lenient threshold
  gene_sets[[paste0(ds_id, "_lenient")]] <- get_upregulated(df, 0.1, 0, 0)
  
  # Stringent threshold
  gene_sets[[paste0(ds_id, "_stringent")]] <- get_upregulated(df, 0.1, 0.8, 10)
}

# Save gene sets as RDS for downstream scripts
saveRDS(gene_sets, file.path(out_dir, "upregulated_gene_sets.rds"))
saveRDS(all_degs, file.path(out_dir, "all_degs_list.rds"))
message("Saved gene sets and DEG lists for downstream analysis")

message("\n=== Summary ===")
for (ds_id in names(all_degs)) {
  lenient_n <- length(gene_sets[[paste0(ds_id, "_lenient")]])
  stringent_n <- length(gene_sets[[paste0(ds_id, "_stringent")]])
  message(sprintf("%s: %d lenient, %d stringent upregulated genes", 
                  ds_id, lenient_n, stringent_n))
}
