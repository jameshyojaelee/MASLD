# H81 n-by-group counts. Read-only, no fit.
suppressMessages({library(edgeR); library(data.table)})
R <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
cnt <- readRDS("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_20261007/assembled_774_primary/counts_raw.rds")
m <- as.data.table(readRDS(file.path(R, "meta_matched.rds")))
d <- readRDS(file.path(R, "merged_dge.rds"))
cat("counts", dim(cnt), storage.mode(cnt), "non-integer", sum(cnt != round(cnt)), "negative", sum(cnt < 0), "\n")
cat("meta rows", nrow(m), "unique sample_id", uniqueN(m$sample_id), "unique analysis_unit_id", uniqueN(m$analysis_unit_id), "\n")
cat("setequal(meta, counts cols)", setequal(m$sample_id, colnames(cnt)), " dge dims", dim(d), " identical(dge cols, meta order)", identical(colnames(d), m$sample_id), "\n")
print(dcast(m[, .N, by = .(dataset, group_binary)], dataset ~ group_binary, value.var = "N"))
print(dcast(m[, .N, by = .(dataset, inferred_sex)], dataset ~ inferred_sex, value.var = "N", fill = 0))
print(m[, .(control = sum(group_binary == "Control"), disease = sum(group_binary == "Disease"), n = .N)])
