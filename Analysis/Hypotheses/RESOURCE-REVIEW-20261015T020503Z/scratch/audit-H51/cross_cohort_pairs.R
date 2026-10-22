suppressMessages(library(data.table))
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
m <- as.data.table(readRDS(file.path(root, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/meta_matched.rds")))
idr <- file.path(root, "figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified")
p <- fread(file.path(idr, "full_filtered_related_pairs.tsv")); setnames(p, sub("^#", "", names(p)))
r <- fread(file.path(idr, "resource_identity_roster.tsv"))
cat("roster cols:", paste(names(r), collapse=","), "\n")
ds <- setNames(r$dataset, r$sample_id)
p[, `:=`(d1 = ds[IID1], d2 = ds[IID2])]
p[, in1 := IID1 %in% m$sample_id][, in2 := IID2 %in% m$sample_id]
p[, cls := fifelse(KINSHIP > 0.354, "dup/MZ", fifelse(KINSHIP > 0.177, "1st", fifelse(KINSHIP > 0.0884, "2nd", "lower")))]
cat("all pairs:", nrow(p), "\n")
print(p[d1 != d2, .N, by = .(pair = paste(pmin(d1,d2), pmax(d1,d2), sep="|"), cls, both_in_774 = in1 & in2)][order(pair, cls)])
