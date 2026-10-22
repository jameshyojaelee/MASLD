suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)})
set.seed(20260923)
base <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out <- commandArgs(TRUE)[1]
arm <- file.path(base,"RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
dge <- readRDS(file.path(arm,"merged_dge.rds")); meta <- as.data.table(readRDS(file.path(arm,"meta_matched.rds")))
cat("dge", dim(dge), "\n"); print(head(dge$samples)); print(table(meta$dataset, meta$group_binary))
meta <- meta[match(colnames(dge), meta$sample_id)]; stopifnot(all(meta$sample_id==colnames(dge)))
res <- list()
## ---------- R03-02: GSE126848 ----------
cr <- fread(file.path(base,"figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv"))
rel <- cr[dataset=="GSE126848"]
ix <- which(meta$dataset=="GSE126848"); d <- dge[,ix]; g <- factor(meta$group_binary[ix],levels=c("Control","Disease"))
L <- d$samples$lib.size; Leff <- L*d$samples$norm.factors
X <- model.matrix(~g)
v <- voomWithQualityWeights(d,X); f <- eBayes(lmFit(v,X)); tt <- topTable(f,coef=2,number=Inf,sort.by="none")
call0 <- tt$adj.P.Val<.05 & abs(tt$logFC)>.5
m <- match(rownames(tt), rel$gene_id_versioned)
relcall <- rel$FDR[m]<.05 & abs(rel$logFC[m])>.5
z <- rowSums(d$counts)==0
X1 <- model.matrix(~g+log2(L)); v1 <- voomWithQualityWeights(d,X1); t1 <- topTable(eBayes(lmFit(v1,X1)),coef=2,number=Inf,sort.by="none")
call1 <- t1$adj.P.Val<.05 & abs(t1$logFC)>.5
k <- filterByExpr(d,group=g); d2 <- calcNormFactors(d[k,,keep.lib.sizes=FALSE],method="RLE")
v2 <- voomWithQualityWeights(d2,X); t2 <- topTable(eBayes(lmFit(v2,X)),coef=2,number=Inf,sort.by="none"); call2 <- t2$adj.P.Val<.05 & abs(t2$logFC)>.5
wt <- wilcox.test(L[g=="Disease"],L[g=="Control"],exact=TRUE)
sw <- v$targets$sample.weights
# weighted log2 effective lib ratio
wl <- function(gg) sum(sw[g==gg]*log2(Leff[g==gg]+1))/sum(sw[g==gg])
pw <- fread(file.path(arm,"lvqw_sample_weights.tsv"))[sample_id %in% colnames(d)]
pw$grp <- g[match(pw$sample_id,colnames(d))]
res$r0302 <- data.table(metric=c("n_ctrl","n_dis","n_genes","med_lib_ctrl","med_lib_dis","pairs_dis_deeper","n_pairs","wilcox_p","cor_grp_log2lib",
 "base_calls","base_up","base_down","match_release_calls","release_calls","logFC_maxabsdiff_vs_release","zero_read_calls","zero_call_logFC_min","zero_call_logFC_max","zero_all_up",
 "weighted_log2_efflib_ratio_ctrl_minus_dis","depthadj_calls","depthadj_kept","depthadj_spearman","fbe_genes","fbe_calls","fbe_up","fbe_down",
 "qw_med_ctrl_percohort","qw_med_dis_percohort","pooledqw_med_ctrl","pooledqw_med_dis"),
 value=c(sum(g=="Control"),sum(g=="Disease"),nrow(d),median(L[g=="Control"]),median(L[g=="Disease"]),
 sum(outer(L[g=="Disease"],L[g=="Control"],">")),sum(g=="Disease")*sum(g=="Control"),wt$p.value,cor(as.numeric(g),log2(L)),
 sum(call0),sum(call0&tt$logFC>0),sum(call0&tt$logFC<0),sum(call0==relcall),sum(relcall),max(abs(tt$logFC-rel$logFC[m])),
 sum(call0&z),min(tt$logFC[call0&z]),max(tt$logFC[call0&z]),all(tt$logFC[call0&z]>0),
 wl("Control")-wl("Disease"),sum(call1),sum(call0&call1),cor(tt$logFC,t1$logFC,method="spearman"),sum(k),sum(call2),sum(call2&t2$logFC>0),sum(call2&t2$logFC<0),
 median(sw[g=="Control"]),median(sw[g=="Disease"]),median(pw$sample_weight[pw$grp=="Control"]),median(pw$sample_weight[pw$grp=="Disease"])))
print(res$r0302)
fwrite(res$r0302,file.path(out,"r0302_check.tsv"),sep="\t")
fwrite(data.table(gene=rownames(tt),call0,zero=z,call1),file.path(out,"r0302_gene_calls.tsv.gz"),sep="\t")
## ---------- R03-09: pooled fixed-weight complexity ----------
deg <- fread(file.path(arm,"deg_results.csv")); base_call <- deg$padj<.05 & abs(deg$logFC)>.5
cat("release pooled calls", sum(base_call), sum(base_call&deg$logFC>0), "\n")
swr <- fread(file.path(arm,"lvqw_sample_weights.tsv")); swr <- swr$sample_weight[match(colnames(dge),swr$sample_id)]
des <- fread(file.path(arm,"model_design.tsv")); stopifnot(all(des$sample_id==colnames(dge)))
X0 <- as.matrix(des[,-1]); 
cnt <- dge$counts; lib <- colSums(cnt)
p <- t(t(cnt)/lib)
top100 <- apply(cnt,2,function(x) sum(sort(x,decreasing=TRUE)[1:100]))/lib
ent <- apply(p,2,function(x){x<-x[x>0]; -sum(x*log(x))})
det <- colSums(cnt>0)
zc <- function(x) ave(x, meta$dataset, FUN=function(y) (y-mean(y))/sd(y))
cov <- data.table(loglib=zc(log(dge$samples$lib.size)), top100=zc(top100), entropy=zc(ent), detected=zc(det))
runfit <- function(Xm, mode){
  if(mode=="A"){ v <- voom(dge,Xm,weights=swr) } else { v <- voom(dge,Xm); v$weights <- t(swr*t(v$weights)) }
  f <- eBayes(lmFit(v,Xm)); tt <- topTable(f,coef="group_binaryDisease",number=Inf,sort.by="none"); tt }
mods <- list(baseline=NULL, loglib="loglib", top100="top100", entropy="entropy", joint=c("loglib","top100","detected"))
rows <- list()
for(mode in c("A","B")) for(nm in names(mods)){
  Xm <- X0; if(!is.null(mods[[nm]])) Xm <- cbind(Xm, as.matrix(cov[,mods[[nm]],with=FALSE]))
  tt <- runfit(Xm,mode); stopifnot(all(rownames(tt)==deg$gene))
  cl <- tt$adj.P.Val<.05 & abs(tt$logFC)>.5
  rows[[length(rows)+1]] <- data.table(mode=mode,model=nm,calls=sum(cl),up=sum(cl&tt$logFC>0),down=sum(cl&tt$logFC<0),kept=sum(cl&base_call),
     spearman=cor(tt$logFC,deg$logFC,method="spearman"))
  print(rows[[length(rows)]])
}
r9 <- rbindlist(rows); fwrite(r9,file.path(out,"r0309_check.tsv"),sep="\t")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
