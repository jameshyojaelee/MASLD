# plan-H51: structure of merged_dge.rds (normalization factors, gene annotation). No model fit.
suppressMessages({library(edgeR)})
d <- readRDS("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/merged_dge.rds")
print(class(d)); print(dim(d)); print(names(d)); print(summary(d$samples$norm.factors)); print(head(colnames(d$genes)))
print(head(rownames(d), 3))
