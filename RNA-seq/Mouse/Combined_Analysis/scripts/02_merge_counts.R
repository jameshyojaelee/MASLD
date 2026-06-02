#!/usr/bin/env Rscript
# 02_merge_counts.R
# Merge featureCounts matrices from all MCD datasets

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

root <- normalizePath(".")
metadata_path <- file.path(root, "metadata", "samples.tsv")
out_dir <- file.path(root, "counts")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("=== Merging Count Matrices ===")

# Read combined metadata
samples <- read_tsv(metadata_path, col_types = cols())
message("Samples in metadata: ", nrow(samples))

# Helper function to read featureCounts and extract sample columns
read_featurecounts <- function(path, sample_ids) {
  message("Reading: ", path)
  
  raw <- read_tsv(path, comment = "#", col_types = cols())
  
  # Clean column names - extract sample ID from BAM path
  colnames(raw) <- colnames(raw) %>%
    basename() %>%
    str_remove(".Aligned.sortedByCoord.out.bam")
  
  # Find available samples
  available <- intersect(sample_ids, colnames(raw))
  message("  Found ", length(available), " of ", length(sample_ids), " samples")
  
  if (length(available) == 0) {
    stop("No matching samples found in ", path)
  }
  
  # Keep Geneid + sample columns
  counts <- raw %>%
    select(Geneid, all_of(available))
  
  return(counts)
}

# Get unique source files
source_files <- samples %>%
  select(source_counts, batch) %>%
  distinct()

# Read each count matrix
count_list <- list()
for (i in seq_len(nrow(source_files))) {
  src <- source_files$source_counts[i]
  batch <- source_files$batch[i]
  
  batch_samples <- samples %>%
    filter(source_counts == src, batch == !!batch) %>%
    pull(sample_id)
  
  count_list[[batch]] <- read_featurecounts(src, batch_samples)
}

# Merge all count matrices by Geneid
message("\nMerging count matrices...")
merged <- count_list[[1]]
for (i in 2:length(count_list)) {
  merged <- inner_join(merged, count_list[[i]], by = "Geneid")
}

message("Final matrix: ", nrow(merged), " genes × ", ncol(merged) - 1, " samples")

# Verify all samples present
sample_cols <- setdiff(colnames(merged), "Geneid")
missing <- setdiff(samples$sample_id, sample_cols)
if (length(missing) > 0) {
  warning("Missing samples: ", paste(missing, collapse = ", "))
}

# Reorder columns to match metadata order
sample_order <- intersect(samples$sample_id, sample_cols)
merged <- merged %>%
  select(Geneid, all_of(sample_order))

# Write merged counts
out_file <- file.path(out_dir, "combined_gene_counts.tsv")
write_tsv(merged, out_file)
message("\nWrote merged counts to: ", out_file)

# Summary stats
message("\n=== Summary ===")
message("Total genes: ", nrow(merged))
message("Total samples: ", length(sample_order))
message("Samples by batch:")
print(table(samples$batch[samples$sample_id %in% sample_order]))
