suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)})
set.seed(20261015)
arm <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
meta <- as.data.table(readRDS(file.path(arm,"meta_matched.rds"))); dge_all <- readRDS(file.path(arm,"merged_dge.rds"))
sel <- meta[dataset=="GSE126848"]; setorder(sel,sample_id); d <- dge_all[,match(sel$sample_id,colnames(dge_all))]
info <- data.frame(group=factor(sel$group_binary,levels=c("Control","Disease")))
des <- model.matrix(~group,info)
fit <- function(dd){v<-voomWithQualityWeights(dd,des,plot=FALSE);t<-topTable(eBayes(lmFit(v,des)),coef="groupDisease",number=Inf,sort.by="none");t$FDR<-p.adjust(t$P.Value,"BH");t}
t0 <- fit(d); c0 <- t0$FDR<0.05 & abs(t0$logFC)>0.5
keep <- filterByExpr(d, design=des); cat("genes kept by per-cohort filterByExpr:", sum(keep), "of", nrow(d), "\n")
t1 <- fit(d[keep,,keep.lib.sizes=TRUE]); c1 <- t1$FDR<0.05 & abs(t1$logFC)>0.5
cat("calls, per-cohort filtered family:", sum(c1), " up/down:", sum(c1 & t1$logFC>0), sum(c1 & t1$logFC<0), "\n")
cat("original calls inside filtered gene set:", sum(c0[keep]), " original calls outside it:", sum(c0[!keep]), " up among outside:", sum(c0[!keep] & t0$logFC[!keep]>0), "\n")
