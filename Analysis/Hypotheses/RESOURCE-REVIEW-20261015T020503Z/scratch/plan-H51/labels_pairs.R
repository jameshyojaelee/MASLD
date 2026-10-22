# plan-H51: label agreement with SraRunTable and same-genotype pair label discordance. Metadata only.
suppressMessages({library(data.table)})
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
rel <- file.path(root, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
m <- as.data.table(readRDS(file.path(rel, "meta_matched.rds")))
sra <- fread(file.path(root, "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE213621/metadata/SraRunTable.csv"))
cat("SraRunTable rows:", nrow(sra), "\n"); print(sra[, .N, by=fibrotic_stage])
cat("Instrument x stage:\n"); print(sra[, .N, by=.(Instrument, LibraryLayout, AvgSpotLen>0)])
cat("cell_type/source_name values:\n"); print(sra[, .N, by=.(cell_type, source_name)])
g <- m[dataset == "GSE213621"]
x <- merge(g[, .(sample_id, condition)], sra[, .(Run, fibrotic_stage)], by.x="sample_id", by.y="Run")
cat("roster GSE213621 matched to SRA:", nrow(x), " of ", nrow(g), "\n")
print(table(x$condition, x$fibrotic_stage))
# Same-genotype pairs
idr <- file.path(root, "figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified")
p <- fread(file.path(idr, "full_filtered_related_pairs.tsv")); setnames(p, sub("^#", "", names(p)))
r <- fread(file.path(idr, "resource_identity_roster.tsv"))
cat("pairs cols:", paste(names(p), collapse=","), "\n")
ds <- setNames(r$dataset, r$sample_id)
lab <- setNames(sra$fibrotic_stage, sra$Run)
p[, `:=`(d1 = ds[IID1], d2 = ds[IID2])]
q <- p[d1 == "GSE213621" & d2 == "GSE213621" & KINSHIP > 0.45]
q[, `:=`(l1 = lab[IID1], l2 = lab[IID2])]
q[, in774_1 := IID1 %in% m$sample_id][, in774_2 := IID2 %in% m$sample_id]
cat("GSE213621 within-cohort pairs with kinship > 0.45:", nrow(q), "\n")
q2 <- q[!is.na(l1) & !is.na(l2)]
cat("  with both labels:", nrow(q2), "; both disease:", nrow(q2[l1 != "Control" & l2 != "Control"]), "\n")
q2[, combo := apply(cbind(l1, l2), 1, function(z) paste(sort(z), collapse="|"))]
print(q2[, .N, by=.(same = l1 == l2, combo)][order(same, combo)])
cat("  pairs with both members in the 774 roster:", nrow(q2[in774_1 & in774_2]), "\n")
cat("  discordant pairs with at least one member in the 774 roster:", nrow(q2[l1 != l2 & (in774_1 | in774_2)]), "\n")
cat("  discordant pairs by roster membership:\n"); print(q2[l1 != l2, .N, by=.(in774_1, in774_2)])
