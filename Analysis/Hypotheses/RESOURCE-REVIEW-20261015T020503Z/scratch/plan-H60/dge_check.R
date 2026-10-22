suppressPackageStartupMessages(library(edgeR))
REL <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/"
dge <- readRDS(paste0(REL, "merged_dge.rds"))
cat("class", class(dge), "norm.factors range", range(dge$samples$norm.factors), "n != 1", sum(dge$samples$norm.factors != 1), "\n")
cat("sample cols", paste(names(dge$samples), collapse=","), "\n")
