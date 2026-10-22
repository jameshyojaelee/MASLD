# H81 audit: does build_dge (review_refit_v2.R lines 22-31) on the unthinned, unpermuted counts reproduce merged_dge.rds?
# Read-only. No fit, no randomness.
suppressPackageStartupMessages({library(edgeR); library(data.table)})
R <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
cnt <- readRDS("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_20261007/assembled_774_primary/counts_raw.rds")
meta <- as.data.table(readRDS(file.path(R, "meta_matched.rds")))
md <- readRDS(file.path(R, "merged_dge.rds"))
m <- cnt[, meta$sample_id, drop = FALSE]
dge <- DGEList(counts = m)
dge$samples$group_binary <- factor(meta$group_binary, levels = c("Control", "Disease"))
keep <- filterByExpr(dge, design = model.matrix(~ 0 + group_binary, data = dge$samples))
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE], method = "RLE")
o <- match(colnames(dge), colnames(md))
cat("genes", nrow(dge), "identical gene set/order", identical(rownames(dge), rownames(md)), "\n")
cat("max |count diff|", max(abs(dge$counts - md$counts[, o])), "\n")
cat("max |lib.size diff|", max(abs(dge$samples$lib.size - md$samples$lib.size[o])), "\n")
cat("max |norm.factors diff|", max(abs(dge$samples$norm.factors - md$samples$norm.factors[o])), "\n")
