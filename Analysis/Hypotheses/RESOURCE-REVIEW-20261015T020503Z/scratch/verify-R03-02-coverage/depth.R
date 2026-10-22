suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)})
set.seed(20261015)
arm <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
cand <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/bulk-context-774-20261014T140904Z_22236311"
dge_all <- readRDS(file.path(arm,"merged_dge.rds")); meta <- as.data.table(readRDS(file.path(arm,"meta_matched.rds")))
cat("lib.size == colSums(counts)? ", isTRUE(all.equal(as.numeric(dge_all$samples$lib.size), as.numeric(colSums(dge_all$counts)))), "\n")
res <- fread(file.path(cand,"cohort_disease_all_gene_results.tsv"))
qw <- fread(file.path(cand,"GSE126848_quality_weights.tsv"))
sel <- meta[dataset=="GSE126848"]; setorder(sel,sample_id)
d <- dge_all[,match(sel$sample_id,colnames(dge_all))]
lib <- d$samples$lib.size; grp <- sel$group_binary
cat("median lib.size (M) control/disease:", median(lib[grp=="Control"])/1e6, median(lib[grp=="Disease"])/1e6, "\n")
cat("range control:", range(lib[grp=="Control"])/1e6, " disease:", range(lib[grp=="Disease"])/1e6, "\n")
w <- wilcox.test(lib[grp=="Disease"], lib[grp=="Control"]); print(w)
pairs <- outer(lib[grp=="Disease"], lib[grp=="Control"], ">"); cat("disease>control pairs:", sum(pairs), "of", length(pairs), " AUC:", mean(pairs), "\n")
q <- qw[match(sel$sample_id, sample_id), quality_weight]
cat("median quality weight control/disease:", median(q[grp=="Control"]), median(q[grp=="Disease"]), "\n")
cat("cor(log lib, qweight) all 30 spearman:", cor(log2(lib), q, method="spearman"), "\n")
ctl <- sel$condition[grp=="Control"]; lc <- lib[grp=="Control"]
print(wilcox.test(lc[ctl=="Control"], lc[ctl=="Control_Obese"]))
cat("median lean/obese control lib (M):", median(lc[ctl=="Control"])/1e6, median(lc[ctl=="Control_Obese"])/1e6, "\n")
# effective lib with norm factors
eff <- lib*d$samples$norm.factors
cat("median effective lib (M) control/disease:", median(eff[grp=="Control"])/1e6, median(eff[grp=="Disease"])/1e6, "\n")
cat("median norm.factor control/disease:", median(d$samples$norm.factors[grp=="Control"]), median(d$samples$norm.factors[grp=="Disease"]), "\n")
# Depth slope in controls; per cohort
slope_fun <- function(ds){
  s <- meta[dataset==ds]; dd <- dge_all[,match(s$sample_id,colnames(dge_all))]
  E <- cpm(dd, log=TRUE, prior.count=0.5)  # uses lib.size*norm.factors
  L <- log2(dd$samples$lib.size); g <- s$group_binary
  Ec <- E[, g=="Control", drop=FALSE]; x <- L[g=="Control"]; xc <- x-mean(x)
  b <- as.numeric((Ec - rowMeans(Ec)) %*% xc / sum(xc^2))
  gap <- mean(L[g=="Disease"]) - mean(L[g=="Control"])
  r <- res[dataset==ds][match(rownames(dd), gene_id_versioned)]
  call <- r$FDR<0.05 & abs(r$logFC)>0.5
  data.table(dataset=ds, n_ctl=sum(g=="Control"), n_dis=sum(g=="Disease"), gap_log2=gap,
    med_ctl_M=median(dd$samples$lib.size[g=="Control"])/1e6, med_dis_M=median(dd$samples$lib.size[g=="Disease"])/1e6,
    rho_slope_logFC_all=cor(b, r$logFC, method="spearman"),
    rho_slope_logFC_calls=cor(b[call], r$logFC[call], method="spearman"),
    n_calls=sum(call),
    sign_agree_calls=mean(sign(b[call]*gap)==sign(r$logFC[call])))
}
out <- rbindlist(lapply(c("GSE126848","GSE130970","GSE135251","GSE162694","GSE213621"), slope_fun))
print(out)
# predicted logFC = slope*gap; Spearman with observed equals sign(gap)*rho
# Sensitivity: add log2 depth to GSE126848 model
info <- data.frame(group=factor(grp,levels=c("Control","Disease")), ld=log2(lib)-mean(log2(lib)))
cat("cor(group, log2 depth):", cor(as.numeric(info$group=="Disease"), info$ld), "\n")
des2 <- model.matrix(~group+ld, info)
v2 <- voomWithQualityWeights(d, des2, plot=FALSE); f2 <- eBayes(lmFit(v2, des2)); t2 <- topTable(f2, coef="groupDisease", number=Inf, sort.by="none")
t2$FDR <- p.adjust(t2$P.Value,"BH")
cat("depth-adjusted GSE126848 calls (FDR<0.05 & |logFC|>0.5):", sum(t2$FDR<0.05 & abs(t2$logFC)>0.5), "\n")
r1 <- res[dataset=="GSE126848"][match(rownames(d), gene_id_versioned)]
c1 <- r1$FDR<0.05 & abs(r1$logFC)>0.5; c2 <- t2$FDR<0.05 & abs(t2$logFC)>0.5
cat("original calls:", sum(c1), " retained after adj:", sum(c1&c2), " new:", sum(c2&!c1), "\n")
cat("spearman logFC orig vs adj:", cor(r1$logFC, t2$logFC, method="spearman"), "\n")
# reproduce original 6147 from same pipeline as check
des1 <- model.matrix(~group, info); v1 <- voomWithQualityWeights(d, des1, plot=FALSE); f1 <- eBayes(lmFit(v1,des1)); tt <- topTable(f1, coef="groupDisease", number=Inf, sort.by="none")
cat("reproduced calls:", sum(p.adjust(tt$P.Value,"BH")<0.05 & abs(tt$logFC)>0.5), "\n")
# unique-to-GSE126848 genes in the union
rr <- res[, .(call=FDR<0.05 & abs(logFC)>0.5), by=.(dataset,gene_id_versioned)]
m <- rr[call==TRUE, .(n=.N, ds=paste(sort(dataset),collapse=",")), by=gene_id_versioned]
cat("union:", nrow(m), " single-cohort:", sum(m$n==1), " unique to GSE126848:", sum(m$n==1 & m$ds=="GSE126848"), "\n")
print(m[n==1, .N, by=ds])
