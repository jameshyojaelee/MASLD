# CGA2: per-cohort TMM logCPM for the 774-unit corrected substrate (primary cohorts only).
suppressMessages(library(edgeR))
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- readRDS("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/merged_dge.rds")
for (ds in c("GSE162694", "GSE213621")) {
  x <- d[, d$samples$dataset == ds]
  keep <- rowMeans(cpm(x) > 1) >= 0.2
  x <- calcNormFactors(x[keep, , keep.lib.sizes = FALSE], method = "TMM")
  m <- cpm(x, log = TRUE, prior.count = 1)
  rownames(m) <- sub("\\.\\d+$", "", rownames(m))
  write.table(data.frame(g = rownames(m), m, check.names = FALSE), gzfile(file.path(out, paste0(ds, "_logcpm.tsv.gz"))),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(ds, dim(m), "\n")
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
