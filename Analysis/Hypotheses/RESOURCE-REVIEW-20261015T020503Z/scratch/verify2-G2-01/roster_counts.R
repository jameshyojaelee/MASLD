suppressPackageStartupMessages(library(data.table))
REL <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
m <- as.data.table(readRDS(file.path(REL, "meta_matched.rds")))
cat("meta rows", nrow(m), "\n"); print(names(m))
st <- m[dataset %in% c("GSE130970","GSE135251","GSE162694") & fibrosis_stage %in% 0:4 & !is.na(inferred_sex)]
st[, fib_group := paste0("F", fibrosis_stage)]
print(dcast(st[, .N, by=.(dataset, fib_group)], dataset ~ fib_group, value.var="N"))
print(st[fib_group=="F0", .N, by=.(dataset, group_binary)])
print(st[fib_group %in% c("F0","F1"), .N, by=.(dataset, fib_group, inferred_sex)][order(dataset, fib_group, inferred_sex)])
man <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/stage_sample_manifest.tsv")
cat("released stage manifest rows", nrow(man), " my roster rows", nrow(st), " setequal:", setequal(man$sample_id, st$sample_id), "\n")
