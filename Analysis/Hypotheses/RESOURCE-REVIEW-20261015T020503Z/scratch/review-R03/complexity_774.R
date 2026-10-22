suppressMessages({library(edgeR); library(data.table)})
rel <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
out <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R03"
d <- readRDS(file.path(rel, "merged_dge.rds"))
m <- as.data.table(readRDS(file.path(rel, "meta_matched.rds")))
w <- fread(file.path(rel, "lvqw_sample_weights.tsv"))
stopifnot(setequal(colnames(d), m$sample_id), nrow(m) == ncol(d))
cnt <- d$counts; storage.mode(cnt) <- "double"
p <- sweep(cnt, 2, colSums(cnt), "/")
ent <- apply(p, 2, function(x) { x <- x[x > 0]; -sum(x * log(x)) })
gini <- apply(cnt, 2, function(x) { x <- sort(as.numeric(x)); n <- length(x); sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x)) })
topn <- function(k) apply(p, 2, function(x) sum(sort(x, decreasing = TRUE)[seq_len(k)]))
cx <- data.table(sample_id = colnames(d), lib_size = colSums(cnt), detected = colSums(cnt > 0),
  entropy = ent, gini = gini, top10 = topn(10), top100 = topn(100), norm_factor = d$samples$norm.factors)
cx <- merge(cx, m[, .(sample_id, dataset, group_binary, fibrosis_stage, nas_score, inferred_sex)], by = "sample_id")
cx <- merge(cx, w, by = "sample_id")
fwrite(cx, file.path(out, "complexity_774.tsv"), sep = "\t")
res <- list()
for (ds in sort(unique(cx$dataset))) for (v in c("lib_size","detected","entropy","gini","top10","top100","norm_factor","sample_weight")) {
  s <- cx[dataset == ds]; a <- s[group_binary == "Control"][[v]]; b <- s[group_binary == "Disease"][[v]]
  wt <- suppressWarnings(wilcox.test(b, a))
  res[[length(res)+1]] <- data.table(dataset = ds, metric = v, n_control = length(a), n_disease = length(b),
    median_control = median(a), median_disease = median(b), auc_disease_gt_control = as.numeric(wt$statistic) / (length(a)*length(b)), wilcox_p = wt$p.value)
}
res <- rbindlist(res); res[, bh_q_within_metric_family := p.adjust(wilcox_p, "BH")]
fwrite(res, file.path(out, "complexity_by_group.tsv"), sep = "\t")
print(res[metric %in% c("lib_size","detected","entropy","top100","norm_factor","sample_weight")], digits = 3)
# stage association within disease, exact-stage cohorts
st <- list()
for (ds in c("GSE130970","GSE135251","GSE162694")) for (v in c("lib_size","detected","entropy","gini","top10","top100","norm_factor","sample_weight")) {
  s <- cx[dataset == ds & !is.na(fibrosis_stage)]
  ct <- suppressWarnings(cor.test(as.numeric(s$fibrosis_stage), s[[v]], method = "spearman"))
  st[[length(st)+1]] <- data.table(dataset = ds, metric = v, n = nrow(s), rho_vs_fibrosis = unname(ct$estimate), p = ct$p.value)
}
st <- rbindlist(st); st[, bh_q := p.adjust(p, "BH")]
fwrite(st, file.path(out, "complexity_vs_fibrosis.tsv"), sep = "\t")
print(st, digits = 3)
