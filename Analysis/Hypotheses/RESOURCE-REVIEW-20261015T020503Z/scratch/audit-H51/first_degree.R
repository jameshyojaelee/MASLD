suppressMessages(library(data.table))
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
m <- as.data.table(readRDS(file.path(root, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/meta_matched.rds")))
p <- fread(file.path(root, "figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/full_filtered_related_pairs.tsv")); setnames(p, sub("^#", "", names(p)))
q <- p[IID1 %in% m$sample_id & IID2 %in% m$sample_id]
print(q[, .(IID1, IID2, KINSHIP)])
print(m[sample_id %in% c(q$IID1, q$IID2), .(sample_id, dataset, group_binary, condition, fibrosis_stage)])
