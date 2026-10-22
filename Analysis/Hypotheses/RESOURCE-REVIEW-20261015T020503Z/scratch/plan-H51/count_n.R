# plan-H51: count units by cohort, bin, sex and accession-quartile stratum. Metadata only; no expression fit.
suppressMessages({library(data.table)})
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
rel <- file.path(root, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
m <- as.data.table(readRDS(file.path(rel, "meta_matched.rds")))
cat("meta rows:", nrow(m), " cols:", paste(names(m), collapse=","), "\n")
g <- m[dataset == "GSE213621"]
cat("\nGSE213621 condition x fibrosis_stage x inferred_sex:\n")
print(dcast(g[, .N, by=.(condition, fibrosis_stage, inferred_sex)], condition + fibrosis_stage ~ inferred_sex, value.var="N", fill=0))
gd <- g[group_binary == "Disease"]
cat("GSE213621 disease units:", nrow(gd), " NA sex:", sum(is.na(gd$inferred_sex)), "\n")
gd[, acc := as.numeric(sub("SRR", "", sample_id))]
gd[, acc_rank := frank(acc)]
gd[, acc_q := cut(acc_rank, quantile(acc_rank, 0:4/4), include.lowest=TRUE, labels=paste0("Q",1:4))]
cat("\nGSE213621 disease: accession quartile x sex x condition (strata for the permutation):\n")
print(dcast(gd[, .N, by=.(acc_q, inferred_sex, condition)], acc_q + inferred_sex ~ condition, value.var="N", fill=0))
cat("stratum sizes (acc_q x sex):\n"); print(gd[, .N, by=.(acc_q, inferred_sex)][order(acc_q, inferred_sex)])
bin_of <- function(fs) fifelse(fs <= 1L, "F0F1", fifelse(fs == 2L, "F2", "F3F4"))
d <- m[dataset %in% c("GSE130970","GSE135251","GSE162694") & group_binary == "Disease"]
cat("\nDiscovery disease units:", nrow(d), " with stage 0-4:", sum(d$fibrosis_stage %in% 0:4), " NA sex:", sum(is.na(d$inferred_sex)), "\n")
d <- d[fibrosis_stage %in% 0:4 & !is.na(inferred_sex)]
d[, bin := bin_of(fibrosis_stage)]
print(dcast(d[, .N, by=.(dataset, bin)], dataset ~ bin, value.var="N"))
cat("by stage:\n"); print(dcast(d[, .N, by=.(dataset, fibrosis_stage)], dataset ~ fibrosis_stage, value.var="N", fill=0))
cat("totals by bin:\n"); print(d[, .N, by=bin][order(bin)])
cat("\nLOCO: units in the two-cohort hit-selection fits and held-out cohort:\n")
for (h in c("GSE130970","GSE135251","GSE162694")) {
  s <- d[dataset != h]; o <- d[dataset == h]
  cat(h, " selection units", nrow(s), "(", paste(s[, .N, by=bin][order(bin)]$N, collapse="/"), ")",
      " held-out units", nrow(o), "(", paste(o[, .N, by=bin][order(bin)]$N, collapse="/"), ")\n")
}
cat("\nControls in 3E F0/F1/F2 reference (source controls with a stage):\n")
print(m[dataset %in% c("GSE130970","GSE135251","GSE162694") & group_binary == "Control", .N, by=.(dataset, fibrosis_stage)][order(dataset, fibrosis_stage)])
