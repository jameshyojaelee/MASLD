# H66 plan: count donors per lineage after donor collapse, with the minimum-cells rule. No gene detection is computed.
suppressPackageStartupMessages(library(data.table))
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
source(file.path(BASE, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))
m <- build_srr_to_donor_map(BASE)
cc <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/lineage_cell_counts.tsv"))
meta <- fread(file.path(BASE, "Analysis/SingleCell/metadata/snrna_donor_sex.csv"))[, .(sample, dataset, condition)]
cc <- merge(cc, meta, by = "sample", all.x = TRUE)
cat("runs", nrow(cc), "runs without dataset", sum(is.na(cc$dataset)), "\n")
cat("GSE289173 present:", any(grepl("GSE289173", cc$dataset)), "\n")
cc[, donor := ifelse(sample %in% names(m), m[sample], sample)]
lin <- grep("^n_", names(cc), value = TRUE)
d <- cc[, lapply(.SD, sum), by = .(donor, dataset), .SDcols = lin]
cat("donors", nrow(d), "\n"); print(d[, .N, by = dataset][order(-N)])
res <- rbindlist(lapply(lin, function(l) data.table(lineage = sub("^n_", "", l),
  runs_with_cells = sum(cc[[l]] > 0), donors_ge1 = sum(d[[l]] >= 1), donors_ge20 = sum(d[[l]] >= 20), donors_ge50 = sum(d[[l]] >= 50))))
print(res)
# pseudobulk column match and integer check (Hepatocytes, Cholangiocytes, Fibroblasts headers + first 200 rows)
for (L in c("Hepatocytes", "Cholangiocytes", "Fibroblasts")) {
  pb <- fread(file.path(BASE, sprintf("Analysis/SingleCell/results_gpu_v2/pseudobulk/%s_pseudobulk.csv", L)), nrows = 200)
  cols <- names(pb)[-1]
  v <- as.matrix(pb[, -1])
  cat(L, "cols", length(cols), "in cell-count table", sum(cols %in% cc$sample), "integer-valued", all(abs(v - round(v)) < 1e-6), "\n")
}
