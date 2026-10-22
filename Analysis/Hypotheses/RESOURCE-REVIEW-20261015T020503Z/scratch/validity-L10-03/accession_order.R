# Blinded check: does GSE213621 accession order (a proxy for processing order) track bin labels,
# and is there an expression axis along accession order? The observed bin-expression contrast is NOT fitted.
suppressMessages({library(data.table); library(edgeR)})
set.seed(20261015)
rel <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
L <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/validity-L10-03"
d <- readRDS(file.path(rel, "merged_dge.rds"))
m <- as.data.table(readRDS(file.path(rel, "meta_matched.rds")))
g <- m[dataset == "GSE213621" & group_binary == "Disease"]
g[, acc := as.integer(sub("SRR", "", sample_id))]
g[, acc_rank := frank(acc)]
g[, bin := c("F0F1","F2","F3F4")[fibrosis_stage]]
cat("Disease units:", nrow(g), "\n")
cat("Labels only: Spearman(accession rank, bin 1/2/3) =", signif(cor(g$acc_rank, g$fibrosis_stage, method = "spearman"), 3),
    " p =", signif(cor.test(g$acc_rank, g$fibrosis_stage, method = "spearman", exact = FALSE)$p.value, 3), "\n")
g[, acc_q := cut(acc_rank, quantile(acc_rank, 0:4/4), include.lowest = TRUE, labels = paste0("Q", 1:4))]
print(dcast(g[, .N, by = .(acc_q, bin)], acc_q ~ bin, value.var = "N"))
cat("Chi-square accession quartile x bin p =", signif(chisq.test(table(g$acc_q, g$bin))$p.value, 3), "\n")
cat("Sex by accession quartile:\n"); print(table(g$acc_q, g$inferred_sex))
# Expression axes without labels
dg <- d[, g$sample_id]
lc <- cpm(dg, log = TRUE, prior.count = 2)
keep <- rowMeans(lc) > 2
lc <- lc[keep, ]
v <- apply(lc, 1, var)
top <- lc[order(-v)[1:2000], ]
pc <- prcomp(t(top), center = TRUE, scale. = FALSE)
ve <- pc$sdev^2 / sum(pc$sdev^2)
lib <- colSums(dg$counts)
p <- sweep(dg$counts, 2, lib, "/")
top50 <- apply(p, 2, function(x) sum(sort(x, decreasing = TRUE)[1:50]))
res <- rbindlist(lapply(1:6, function(k) data.table(PC = k, var_expl = round(ve[k], 3),
  rho_accession = round(cor(pc$x[, k], g$acc_rank, method = "spearman"), 3),
  p_accession = signif(cor.test(pc$x[, k], g$acc_rank, method = "spearman", exact = FALSE)$p.value, 3),
  rho_top50 = round(cor(pc$x[, k], top50, method = "spearman"), 3),
  rho_sexM = round(cor(pc$x[, k], as.integer(g$inferred_sex == "M"), method = "spearman"), 3))))
print(res)
cat("Spearman(accession rank, top50 share) =", signif(cor(g$acc_rank, top50, method = "spearman"), 3), "\n")
cat("Spearman(accession rank, library size) =", signif(cor(g$acc_rank, lib, method = "spearman"), 3), "\n")
fwrite(res, file.path(L, "accession_order_pcs.tsv"), sep = "\t")
fwrite(g[, .(sample_id, acc_rank, acc_q, bin, inferred_sex)], file.path(L, "accession_order_bins.tsv"), sep = "\t")
