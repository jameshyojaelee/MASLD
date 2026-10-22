# plan-H51: do the rebuilt fraction files cover every analysis unit? Coverage only.
suppressMessages(library(data.table))
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
m <- as.data.table(readRDS(file.path(root, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/meta_matched.rds")))
for (ds in c("GSE130970","GSE135251","GSE162694","GSE213621")) {
  f <- read.delim(file.path(root, "Analysis/Deconvolution/corrected_774_20261009T141331Z/results", ds, paste0(ds, "_bayesprism_proportions.tsv")), check.names = FALSE)
  u <- m[dataset == ds & group_binary == "Disease", sample_id]
  cat(ds, " fraction rows", nrow(f), " cols", ncol(f), " disease units", length(u), " covered", sum(u %in% rownames(f)),
      " max|rowsum-1|", signif(max(abs(rowSums(f) - 1)), 3), "\n")
}
