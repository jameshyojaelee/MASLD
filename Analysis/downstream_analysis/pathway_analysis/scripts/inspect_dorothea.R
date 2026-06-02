
# inspect_dorothea.R
library(tidyverse)
conda_lib <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/pathway_analysis/.mamba/pathway_analysis/lib/R/library"
.libPaths(c(conda_lib))

dorothea_list <- readRDS("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/pathway_analysis/results/refactored/all_dorothea_scores.rds")

cat("Datasets in DoRothEA list:\n")
print(names(dorothea_list))

for (n in names(dorothea_list)) {
  res <- dorothea_list[[n]]
  if (is.null(res)) {
      cat(sprintf("%s: NULL\n", n))
  } else {
      cat(sprintf("%s: %d rows, %d unique TFs\n", n, nrow(res), length(unique(res$source))))
      cat("First 5 TFs:\n")
      print(head(unique(res$source), 5))
  }
}
