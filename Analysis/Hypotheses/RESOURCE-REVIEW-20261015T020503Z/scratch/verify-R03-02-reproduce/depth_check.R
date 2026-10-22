suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)})
set.seed(20261015)
base <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
arm <- file.path(base,"RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
ctx <- file.path(base,"figures/candidates/bulk-context-774-20261014T140904Z_22236311")
dge_all <- readRDS(file.path(arm,"merged_dge.rds"))
meta <- as.data.table(readRDS(file.path(arm,"meta_matched.rds")))
cat("dims", dim(dge_all), "\n")
cat("lib.size == colSums(counts)? ", isTRUE(all.equal(as.numeric(dge_all$samples$lib.size), as.numeric(colSums(dge_all$counts)))), "\n")
cat("norm.factors range", range(dge_all$samples$norm.factors), "\n")
meta[, lib := dge_all$samples$lib.size[match(sample_id, colnames(dge_all))]]
stopifnot(!anyNA(meta$lib))
# depth by cohort and group
dep <- meta[, .(n=.N, median_lib_M=median(lib)/1e6, mean_log2=mean(log2(lib))), by=.(dataset,group_binary)]
print(dep[order(dataset,group_binary)])
g <- meta[dataset=="GSE126848"]; setorder(g, sample_id)
dl <- g[group_binary=="Disease",lib]; cl <- g[group_binary=="Control",lib]
cat("GSE126848 median lib (M): disease", median(dl)/1e6, "control", median(cl)/1e6, "\n")
pairs <- outer(dl, cl, ">"); cat("disease>control pairs", sum(pairs), "of", length(pairs), " AUC", mean(pairs)+0.5*mean(outer(dl,cl,"==")), "\n")
print(wilcox.test(dl, cl, exact=TRUE))
cat("disease libs (M):", sort(round(dl/1e6,2)), "\n")
cat("control libs (M):", sort(round(cl/1e6,2)), "\n")
# lean vs obese controls
cc <- g[group_binary=="Control"]
print(table(cc$condition))
print(wilcox.test(cc[condition=="Control",lib], cc[condition=="Control_Obese",lib], exact=TRUE))
# quality weights per-cohort and pooled
qw <- fread(file.path(ctx,"GSE126848_quality_weights.tsv"))
qw <- merge(qw, g[,.(sample_id,group_binary)], by="sample_id")
print(qw[, .(median_w=median(quality_weight), mean_w=mean(quality_weight)), by=group_binary])
pw <- fread(file.path(arm,"lvqw_sample_weights.tsv"))
pw <- merge(pw, g[,.(sample_id,group_binary)], by="sample_id")
print(pw[, .(n=.N, median_pooled_w=median(sample_weight), mean_pooled_w=mean(sample_weight)), by=group_binary])
cat("cor(per-cohort weight, log2 lib):", cor(qw$quality_weight, log2(g$lib[match(qw$sample_id,g$sample_id)]), method="spearman"), "\n")
# voom as in producing script
info <- data.frame(group=factor(g$group_binary,levels=c("Control","Disease")),row.names=g$sample_id)
design <- model.matrix(~group,info)
dge <- dge_all[,match(g$sample_id,colnames(dge_all))]
v <- voomWithQualityWeights(dge,design,plot=FALSE)
cat("weights reproduce file:", isTRUE(all.equal(v$targets$sample.weights, qw$quality_weight[match(g$sample_id,qw$sample_id)])), "\n")
E <- v$E
ctrl <- g$group_binary=="Control"
x <- log2(g$lib[ctrl]); xc <- x-mean(x)
Ec <- E[,ctrl]
slope <- as.numeric((Ec - rowMeans(Ec)) %*% xc / sum(xc^2))
names(slope) <- rownames(E)
gap <- mean(log2(dl)) - mean(log2(cl))
cat("GSE126848 depth gap (mean log2 disease - control):", gap, "\n")
res <- fread(file.path(ctx,"cohort_disease_all_gene_results.tsv"))
out <- list()
for (ds in unique(res$dataset)) {
  r <- res[dataset==ds]; r <- r[match(names(slope), gene_id_versioned)]
  stopifnot(!anyNA(r$logFC))
  calls <- r$FDR<0.05 & abs(r$logFC)>0.5
  pred <- slope*gap
  m <- meta[dataset==ds]
  dgap <- m[group_binary=="Disease", mean(log2(lib))] - m[group_binary=="Control", mean(log2(lib))]
  out[[ds]] <- data.table(dataset=ds, n_genes=nrow(r), n_calls=sum(calls),
    rho_slope_logFC_all=cor(slope, r$logFC, method="spearman"),
    rho_slope_logFC_calls=cor(slope[calls], r$logFC[calls], method="spearman"),
    rho_pred126848gap_logFC_all=cor(pred, r$logFC, method="spearman"),
    sign_agree_pred_calls=mean(sign(pred[calls])==sign(r$logFC[calls])),
    depth_gap_log2=dgap)
}
print(rbindlist(out), digits=4)
# slope from a nested model without groups confounding: slope significance count
fwrite(data.table(gene=names(slope), slope=slope), "slope_controls.tsv.gz", sep="\t")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
