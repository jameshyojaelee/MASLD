suppressPackageStartupMessages({library(edgeR); library(data.table)})
set.seed(20261015)
P <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/"
NEW <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz"
d <- readRDS(paste0(P, "merged_dge.rds"))
cat("class", class(d), "dim", dim(d), "\n")
s <- d$samples
cnt <- d$counts
lib <- colSums(cnt)
top100 <- apply(cnt, 2, function(x) sum(sort(x, decreasing = TRUE)[1:100])) / lib
dz <- fread(paste0(P, "model_design.tsv"))
w <- fread(paste0(P, "lvqw_sample_weights.tsv"))
fr <- fread(NEW)
ids <- dz$sample_id
stopifnot(all(ids %in% colnames(cnt)), all(ids %in% fr$sample_id))
fr <- fr[match(ids, fr$sample_id)]
ds <- fr$dataset
grp <- dz$group_binaryDisease
loglib <- log10(lib[ids]); t100 <- top100[ids]
ct <- setdiff(names(fr), c("sample_id", "dataset"))
X <- as.matrix(fr[, ..ct]); X <- pmax(X, 1e-6); X <- X / rowSums(X)
lf <- log(X)
cat("\nWithin-cohort Spearman of log fraction with log10 library size / top-100 share (pooled after ranking within cohort)\n")
rk <- function(v) ave(v, ds, FUN = function(z) rank(z) / length(z))
for (c in ct) {
  r1 <- cor(rk(lf[, c]), rk(loglib)); r2 <- cor(rk(lf[, c]), rk(t100))
  cat(sprintf("  %-25s lib %6.3f  top100 %6.3f\n", c, r1, r2))
}
cat("\nWithin-cohort disease minus control, median log10 library size and top-100 share\n")
for (g in unique(ds)) {
  i <- ds == g
  cat(sprintf("  %s n=%d  lib D-C %6.3f  top100 D-C %6.3f\n", g, sum(i),
              median(loglib[i & grp == 1]) - median(loglib[i & grp == 0]),
              median(t100[i & grp == 1]) - median(t100[i & grp == 0])))
}
# R2 of 15 ILR-equivalent (clr of 16, rank 15) on top-100 share within dataset
clr <- lf - rowMeans(lf)
fit0 <- lm(t100 ~ factor(ds)); fit1 <- lm(t100 ~ factor(ds) + clr[, -16])
cat(sprintf("\npartial R2 of top-100 share on 15 composition coordinates | dataset: %.3f\n",
            1 - sum(resid(fit1)^2) / sum(resid(fit0)^2)))
fit0 <- lm(loglib ~ factor(ds)); fit1 <- lm(loglib ~ factor(ds) + clr[, -16])
cat(sprintf("partial R2 of log10 library size on 15 composition coordinates | dataset: %.3f\n",
            1 - sum(resid(fit1)^2) / sum(resid(fit0)^2)))
# near-threshold calls
deg <- fread(paste0(P, "deg_results.csv"))
calls <- deg[padj < 0.05 & abs(logFC) > 0.5]
cat("\ncalls", nrow(calls), "; |logFC| in (0.5, 0.6]:", sum(abs(calls$logFC) <= 0.6),
    "; (0.5, 0.55]:", sum(abs(calls$logFC) <= 0.55), "\n")
# expected share failing |logFC|>0.5 after adjustment under NO composition effect (SE diff = SE*sqrt(1.157^2-1))
sed <- calls$SE * sqrt(1.157^2 - 1)
pfail <- pnorm((0.5 - abs(calls$logFC)) / sed)
cat(sprintf("expected calls failing |logFC|>0.5 after adjustment by noise alone (rebuilt inflation 1.157): %.0f of %d\n",
            sum(pfail), nrow(calls)))
