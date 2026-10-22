#!/usr/bin/env Rscript
# Dynamic-range eligibility of the four covariate lineages (p05 >= 1e-4 in every required
# cohort) on the rebuilt composition, against the accepted table, for the three discovery
# cohorts that were rebuilt (GSE130970, GSE135251, GSE162694). GSE174478 and GSE193066 were not rebuilt.
# Environment: COMP_SUMMARY_DIR (existing summary dir; one new file is added).
suppressPackageStartupMessages(library(data.table))
set.seed(20260923)
out <- Sys.getenv("COMP_SUMMARY_DIR", "")
stopifnot(nzchar(out), dir.exists(out), !file.exists(file.path(out, "covariate_dynamic_range.tsv")))
V6 <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10"
old <- fread(file.path(V6, "workstreams/BULK-PROGRAM-MAP-v9/precontract/composition_acceptance/accepted_composition.tsv.gz"))
new <- fread(file.path(out, "composition_new_774.tsv.gz"))
cohorts <- c("GSE130970", "GSE135251", "GSE162694")
lin <- c("Cholangiocytes", "Fibroblasts", "Macrophages", "T cells")
one <- function(x, version) rbindlist(lapply(lin, function(l) {
  p <- x[dataset %in% cohorts, .(p05 = quantile(get(l), 0.05, names = FALSE)), by = dataset]
  data.table(version = version, lineage = l, p05_GSE130970 = p[dataset == "GSE130970", p05],
             p05_GSE135251 = p[dataset == "GSE135251", p05], p05_GSE162694 = p[dataset == "GSE162694", p05],
             min_p05 = min(p$p05), eligible_p05_ge_1e4 = min(p$p05) >= 1e-4)
}))
res <- rbind(one(old, "old_844"), one(new, "new_774"))
fwrite(res, file.path(out, "covariate_dynamic_range.tsv"), sep = "\t")
options(width = 200, digits = 4); print(res)
