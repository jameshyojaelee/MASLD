suppressPackageStartupMessages({library(data.table);library(edgeR)})
arm <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
dge <- readRDS(file.path(arm,"merged_dge.rds")); meta <- as.data.table(readRDS(file.path(arm,"meta_matched.rds")))
print(dim(dge)); print(head(dge$samples)); print(summary(dge$samples$norm.factors))
print(names(meta)); print(meta[,.N,by=.(dataset,group_binary)])
print(meta[dataset=="GSE126848",.N,by=.(group_binary,condition)])
print(all(colnames(dge)==meta$sample_id))
