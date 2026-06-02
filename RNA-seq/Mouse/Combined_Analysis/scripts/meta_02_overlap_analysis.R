#!/usr/bin/env Rscript
# meta_02_overlap_analysis.R
# Compute overlap metrics between combined and individual MCD datasets

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

root <- normalizePath(".")
out_dir <- file.path(root, "meta_analysis")
plot_dir <- file.path(root, "plots")

message("=== Overlap Analysis: Combined vs Individual ===\n")

# Load gene sets from previous step
gene_sets <- readRDS(file.path(out_dir, "upregulated_gene_sets.rds"))

# Define dataset pairs for comparison
datasets <- c("Combined", "InHouse", "GSE156918", "GSE205974")

# Function to compute Jaccard index
jaccard <- function(a, b) {
  intersect_n <- length(intersect(a, b))
  union_n <- length(union(a, b))
  if (union_n == 0) return(NA)
  intersect_n / union_n
}

# Compute Jaccard matrix for each threshold
for (threshold in c("lenient", "stringent")) {
  message("\n--- ", toupper(threshold), " threshold ---")
  
  # Build Jaccard matrix
  jaccard_mat <- matrix(NA, nrow = length(datasets), ncol = length(datasets),
                        dimnames = list(datasets, datasets))
  
  for (i in seq_along(datasets)) {
    for (j in seq_along(datasets)) {
      set_i <- gene_sets[[paste0(datasets[i], "_", threshold)]]
      set_j <- gene_sets[[paste0(datasets[j], "_", threshold)]]
      jaccard_mat[i, j] <- jaccard(set_i, set_j)
    }
  }
  
  message("Jaccard similarity matrix:")
  print(round(jaccard_mat, 3))
  
  # Save
  write.csv(jaccard_mat, file.path(out_dir, paste0("jaccard_matrix_", threshold, ".csv")))
  
  # Identify gained and lost genes (combined vs union of individuals)
  combined_genes <- gene_sets[[paste0("Combined_", threshold)]]
  individual_union <- unique(c(
    gene_sets[[paste0("InHouse_", threshold)]],
    gene_sets[[paste0("GSE156918_", threshold)]],
    gene_sets[[paste0("GSE205974_", threshold)]]
  ))
  
  gained <- setdiff(combined_genes, individual_union)
  lost <- setdiff(individual_union, combined_genes)
  
  message("\nGained in combined (not in any individual): ", length(gained))
  message("Lost in combined (in individual but not combined): ", length(lost))
  
  # Save gained/lost genes
  if (length(gained) > 0) {
    write_csv(tibble(gene_id = gained), 
              file.path(out_dir, paste0("gained_genes_", threshold, ".csv")))
  }
  if (length(lost) > 0) {
    write_csv(tibble(gene_id = lost), 
              file.path(out_dir, paste0("lost_genes_", threshold, ".csv")))
  }
  
  # Generate intersection counts instead of UpSet plot
  message("Computing intersection counts...")
  
  # Prepare list
  upset_list <- list(
    Combined = gene_sets[[paste0("Combined_", threshold)]],
    InHouse = gene_sets[[paste0("InHouse_", threshold)]],
    GSE156918 = gene_sets[[paste0("GSE156918_", threshold)]],
    GSE205974 = gene_sets[[paste0("GSE205974_", threshold)]]
  )
  
  # Calculate key intersections
  all_4 <- length(Reduce(intersect, upset_list))
  combined_only <- length(setdiff(upset_list$Combined, 
                                   union(union(upset_list$InHouse, upset_list$GSE156918), 
                                         upset_list$GSE205974)))
  
  intersect_counts <- tibble(
    intersection = c("All 4 datasets", "Combined only", "InHouse only", "GSE156918 only", "GSE205974 only"),
    count = c(
      all_4,
      combined_only,
      length(setdiff(upset_list$InHouse, union(union(upset_list$Combined, upset_list$GSE156918), upset_list$GSE205974))),
      length(setdiff(upset_list$GSE156918, union(union(upset_list$Combined, upset_list$InHouse), upset_list$GSE205974))),
      length(setdiff(upset_list$GSE205974, union(union(upset_list$Combined, upset_list$InHouse), upset_list$GSE156918)))
    )
  )
  
  write_csv(intersect_counts, file.path(out_dir, paste0("intersection_counts_", threshold, ".csv")))
  message("Key intersections:")
  print(intersect_counts)
  
  # Create bar chart of set sizes
  set_sizes <- tibble(
    dataset = names(upset_list),
    size = sapply(upset_list, length)
  )
  
  p <- ggplot(set_sizes, aes(x = reorder(dataset, -size), y = size)) +
    geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.8) +
    geom_text(aes(label = size), vjust = -0.5, size = 4) +
    labs(
      title = paste0("Upregulated DEGs by Dataset (", threshold, " threshold)"),
      x = "Dataset", y = "Number of Genes"
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(plot_dir, paste0("deg_bar_chart_", threshold, ".pdf")), p, width = 7, height = 5)
  message("Saved: ", file.path(plot_dir, paste0("deg_bar_chart_", threshold, ".pdf")))
}

# Create overlap summary table
overlap_summary <- tibble(
  comparison = c("Combined ∩ InHouse", "Combined ∩ GSE156918", "Combined ∩ GSE205974",
                 "Gained in Combined", "Lost in Combined"),
  lenient = c(
    length(intersect(gene_sets$Combined_lenient, gene_sets$InHouse_lenient)),
    length(intersect(gene_sets$Combined_lenient, gene_sets$GSE156918_lenient)),
    length(intersect(gene_sets$Combined_lenient, gene_sets$GSE205974_lenient)),
    length(setdiff(gene_sets$Combined_lenient, 
                   union(union(gene_sets$InHouse_lenient, gene_sets$GSE156918_lenient), 
                         gene_sets$GSE205974_lenient))),
    length(setdiff(union(union(gene_sets$InHouse_lenient, gene_sets$GSE156918_lenient), 
                         gene_sets$GSE205974_lenient), 
                   gene_sets$Combined_lenient))
  ),
  stringent = c(
    length(intersect(gene_sets$Combined_stringent, gene_sets$InHouse_stringent)),
    length(intersect(gene_sets$Combined_stringent, gene_sets$GSE156918_stringent)),
    length(intersect(gene_sets$Combined_stringent, gene_sets$GSE205974_stringent)),
    length(setdiff(gene_sets$Combined_stringent, 
                   union(union(gene_sets$InHouse_stringent, gene_sets$GSE156918_stringent), 
                         gene_sets$GSE205974_stringent))),
    length(setdiff(union(union(gene_sets$InHouse_stringent, gene_sets$GSE156918_stringent), 
                         gene_sets$GSE205974_stringent), 
                   gene_sets$Combined_stringent))
  )
)

write_csv(overlap_summary, file.path(out_dir, "overlap_summary.csv"))
message("\n=== Overlap Summary ===")
print(overlap_summary)
