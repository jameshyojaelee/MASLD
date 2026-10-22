suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)})
arm <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
cand <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/bulk-context-774-20261014T140904Z_22236311"
meta <- as.data.table(readRDS(file.path(arm,"meta_matched.rds"))); dge_all <- readRDS(file.path(arm,"merged_dge.rds"))
pw <- fread(file.path(arm,"lvqw_sample_weights.tsv")); pw <- merge(pw, meta[,.(sample_id,dataset,group_binary)], by="sample_id")
cat("Pooled-fit weights, median by cohort x group:\n"); print(dcast(pw[,.(med=round(median(sample_weight),3)),by=.(dataset,group_binary)], dataset~group_binary))
# expression stratification of depth-slope prediction in GSE126848 (voom fixed-0.5 log-CPM)
res <- fread(file.path(cand,"cohort_disease_all_gene_results.tsv"))
sel <- meta[dataset=="GSE126848"]; setorder(sel,sample_id); d <- dge_all[,match(sel$sample_id,colnames(dge_all))]
g <- sel$group_binary; L <- log2(d$samples$lib.size)
Eraw <- log2(t((t(d$counts)+0.5)/(d$samples$lib.size*d$samples$norm.factors+1))*1e6)
Escaled <- cpm(d,log=TRUE,prior.count=0.5)
slope <- function(E){Ec<-E[,g=="Control"];x<-L[g=="Control"];xc<-x-mean(x);as.numeric((Ec-rowMeans(Ec))%*%xc/sum(xc^2))}
b1 <- slope(Eraw); b2 <- slope(Escaled); gap <- mean(L[g=="Disease"])-mean(L[g=="Control"])
r <- res[dataset=="GSE126848"][match(rownames(d),gene_id_versioned)]; call <- r$FDR<0.05 & abs(r$logFC)>0.5
ave <- rowMeans(Eraw); medcount <- apply(d$counts,1,median)
bin <- cut(ave, quantile(ave,c(0,.25,.5,.75,1)), include.lowest=TRUE, labels=c("Q1_low","Q2","Q3","Q4_high"))
out <- data.table(bin=bin, call=call, up=r$logFC>0, b1=b1, b2=b2, lfc=r$logFC)[, .(n_genes=.N, n_calls=sum(call), frac_calls_up=round(mean(up[call]),3),
   rho_pred_voomprior=round(cor(b1*gap,lfc,method="spearman"),3), rho_pred_scaledprior=round(cor(b2*gap,lfc,method="spearman"),3),
   signagree_calls_voomprior=round(mean(sign(b1[call]*gap)==sign(lfc[call])),3)), by=bin][order(bin)]
print(out)
cat("calls up/down overall:", sum(call & r$logFC>0), sum(call & r$logFC<0), "\n")
cat("median count across 30 libs among calls, quartiles:\n"); print(quantile(medcount[call],c(.1,.25,.5,.75)))
cat("calls with median count < 10:", sum(call & medcount<10), "\n")
