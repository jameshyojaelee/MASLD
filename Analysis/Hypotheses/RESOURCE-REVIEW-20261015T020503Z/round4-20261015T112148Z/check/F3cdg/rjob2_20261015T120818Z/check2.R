suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)})
set.seed(20260923)
base <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out <- commandArgs(TRUE)[1]
arm <- file.path(base,"RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
dge <- readRDS(file.path(arm,"merged_dge.rds")); meta <- as.data.table(readRDS(file.path(arm,"meta_matched.rds")))
meta <- meta[match(colnames(dge), meta$sample_id)]; stopifnot(all(meta$sample_id==colnames(dge)))
cat("limma", as.character(packageVersion("limma")), "\n")
print(body(voomWithQualityWeights))
cs <- colSums(dge$counts)
cat("lib.size==colSums:", all(abs(dge$samples$lib.size-cs)<1), " max rel diff", max(abs(dge$samples$lib.size/cs-1)), "\n")
## R03-02: what is the zero-read log2FC?
ix <- which(meta$dataset=="GSE126848"); d <- dge[,ix]; g <- factor(meta$group_binary[ix],levels=c("Control","Disease"))
L <- d$samples$lib.size; Leff <- L*d$samples$norm.factors
X <- model.matrix(~g); v <- voomWithQualityWeights(d,X); f <- lmFit(v,X)
z <- which(rowSums(d$counts)==0)
lc <- log2(Leff+1)
cat("zero gene logFC range", range(f$coefficients[z,2]), "\n")
cat("unweighted mean log2(Leff+1) ctrl-dis", mean(lc[g=="Control"])-mean(lc[g=="Disease"]), "\n")
cat("log2 ratio of median Leff", log2(median(Leff[g=="Control"])/median(Leff[g=="Disease"])), "\n")
cat("log2 ratio of median L", log2(median(L[g=="Control"])/median(L[g=="Disease"])), "\n")
w <- v$weights[z[1],]
cat("voom-weighted mean diff (gene 1)", sum(w[g=="Control"]*lc[g=="Control"])/sum(w[g=="Control"]) - sum(w[g=="Disease"]*lc[g=="Disease"])/sum(w[g=="Disease"]), "\n")
## R03-09
deg <- fread(file.path(arm,"deg_results.csv")); base_call <- deg$padj<.05 & abs(deg$logFC)>.5
swr <- fread(file.path(arm,"lvqw_sample_weights.tsv")); swr <- swr$sample_weight[match(colnames(dge),swr$sample_id)]
des <- fread(file.path(arm,"model_design.tsv")); stopifnot(all(des$sample_id==colnames(dge)))
X0 <- as.matrix(des[,-1]); rownames(X0) <- des$sample_id
cnt <- dge$counts
mkcov <- function(lib){
  top100 <- apply(cnt,2,function(x) sum(sort(x,decreasing=TRUE)[1:100]))/lib
  p <- t(t(cnt)/colSums(cnt)); ent <- apply(p,2,function(x){x<-x[x>0]; -sum(x*log(x))})
  det <- colSums(cnt>0)
  zc <- function(x) ave(x, meta$dataset, FUN=function(y) (y-mean(y))/sd(y))
  data.table(loglib=zc(log(lib)), top100=zc(top100), entropy=zc(ent), detected=zc(det)) }
covs <- list(libsize=mkcov(dge$samples$lib.size), efflib=mkcov(dge$samples$lib.size*dge$samples$norm.factors))
runfit <- function(Xm, mode){
  if(mode=="C"){ v <- voom(dge,Xm,weights=swr); v$weights <- t(swr*t(v$weights)) }
  if(mode=="D"){ v <- voomWithQualityWeights(dge,Xm) }
  f <- eBayes(lmFit(v,Xm)); topTable(f,coef="group_binaryDisease",number=Inf,sort.by="none") }
mods <- list(baseline=NULL, loglib="loglib", top100="top100", entropy="entropy", joint=c("loglib","top100","detected"))
rows <- list()
for(mode in c("C","D")) for(cv in names(covs)) for(nm in names(mods)){
  if(cv=="efflib" && nm %in% c("baseline")) next
  Xm <- X0; if(!is.null(mods[[nm]])) Xm <- cbind(Xm, as.matrix(covs[[cv]][,mods[[nm]],with=FALSE]))
  tt <- runfit(Xm,mode); stopifnot(all(rownames(tt)==deg$gene))
  cl <- tt$adj.P.Val<.05 & abs(tt$logFC)>.5
  r <- data.table(mode=mode,cov=cv,model=nm,calls=sum(cl),up=sum(cl&tt$logFC>0),down=sum(cl&tt$logFC<0),kept=sum(cl&base_call),
     spearman=cor(tt$logFC,deg$logFC,method="spearman"), maxabs_vs_release=max(abs(tt$logFC-deg$logFC)))
  print(r); rows[[length(rows)+1]] <- r
}
fwrite(rbindlist(rows),file.path(out,"r0309_check2.tsv"),sep="\t")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
