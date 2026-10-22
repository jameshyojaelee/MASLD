suppressMessages({library(data.table); library(edgeR)})
rel <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
d <- readRDS(file.path(rel, "merged_dge.rds"))
m <- as.data.table(readRDS(file.path(rel, "meta_matched.rds")))
stopifnot(identical(colnames(d), m$sample_id) || all(m$sample_id %in% colnames(d)))
cnt <- d$counts[, m$sample_id]
lib <- colSums(cnt)
p <- sweep(cnt, 2, lib, "/")
ent <- apply(p, 2, function(x){x <- x[x>0]; -sum(x*log2(x))})
top50 <- apply(p, 2, function(x) sum(sort(x, decreasing=TRUE)[1:50]))
detected <- colSums(cnt >= 10)
m[, `:=`(libsize = lib, entropy = ent, top50 = top50, detected = detected)]
w <- fread(file.path(rel, "lvqw_sample_weights.tsv"))
cat("weights cols:", names(w), "\n")
bin_of <- function(fs) fifelse(fs <= 1L, "F0F1", fifelse(fs == 2L, "F2", "F3F4"))
g <- m[dataset == "GSE213621" & group_binary == "Disease"]
g[, bin := factor(c("F0F1","F2","F3F4")[fibrosis_stage])]
disc <- m[dataset %in% c("GSE130970","GSE135251","GSE162694") & group_binary == "Disease" & fibrosis_stage %in% 0:4]
disc[, bin := factor(bin_of(fibrosis_stage))]
both <- rbind(g[, .(dataset, bin, libsize, entropy, top50, detected, inferred_sex)],
              disc[, .(dataset, bin, libsize, entropy, top50, detected, inferred_sex)])
print(both[, .(n=.N, med_lib_M = round(median(libsize)/1e6,2), med_entropy = round(median(entropy),3),
              med_top50 = round(median(top50),4), med_detected = median(detected)), by = .(dataset, bin)][order(dataset, bin)])
for (ds in unique(both$dataset)) {
  x <- both[dataset == ds]
  cat(ds, " Kruskal p: entropy", signif(kruskal.test(entropy ~ bin, x)$p.value,3),
      " top50", signif(kruskal.test(top50 ~ bin, x)$p.value,3),
      " detected", signif(kruskal.test(detected ~ bin, x)$p.value,3),
      " libsize", signif(kruskal.test(libsize ~ bin, x)$p.value,3), "\n")
  cat("   Spearman(entropy, bin) ", signif(cor(x$entropy, as.integer(x$bin), method="spearman"),3),
      " Spearman(top50, bin) ", signif(cor(x$top50, as.integer(x$bin), method="spearman"),3), "\n")
}
fwrite(both, file.path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/validity-L10-03", "complexity_by_bin.tsv"), sep = "\t")
