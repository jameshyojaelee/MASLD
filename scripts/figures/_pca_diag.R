suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma);library(matrixStats)})
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(BASE,"scripts/figures/load_figure_data.R"))
MEGA=c("GSE126848","GSE130970","GSE135251","GSE162694","GSE213621")
dge=load_merged_dge(); s=as.data.table(dge$samples,keep.rownames="sample_id")
k=s$dataset%in%MEGA; dge=dge[,k]; s=s[k]
grp=factor(s$group_binary,levels=c("Control","Disease")); ds=factor(s$dataset)
lc=edgeR::cpm(dge,log=TRUE,prior.count=1)
xi=rownames(dge)[grep("^ENSG00000229807",rownames(dge))];dy=rownames(dge)[grep("^ENSG00000067048",rownames(dge))]
sx=s$sex; es=t(lc[c(xi[1],dy[1]),,drop=FALSE]);km=kmeans(es,2,nstart=20,iter.max=50)
fem=as.integer(names(which.max(tapply(es[,1],km$cluster,mean))));inf=ifelse(km$cluster==fem,"F","M")
mi=is.na(s$sex)|s$sex=="";sx[mi]=inf[mi];sx[is.na(sx)|sx==""]="Unknown"
# HVG within-cohort top-2000
wc=lc;for(d in unique(s$dataset)){i=which(s$dataset==d);wc[,i]=lc[,i]-rowMeans(lc[,i,drop=FALSE])}
idx_hvg=order(matrixStats::rowVars(wc),decreasing=TRUE)[1:2000]
# DEG top-200 by |t|
deg=fread(file.path(BASE,"RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"))
deg[,base:=sub("\\..*","",gene)];rn=sub("\\..*","",rownames(dge))
deg=deg[base%in%rn&is.finite(t)];deg[,a:=abs(t)];setorder(deg,-a)
idx_deg=match(head(deg$base,200),rn);idx_deg=idx_deg[!is.na(idx_deg)]
report=function(idx,name){
  X=limma::removeBatchEffect(lc[idx,],batch=ds,covariates=as.numeric(sx=="F"),design=model.matrix(~grp))
  pr=prcomp(t(X),center=TRUE,scale.=FALSE); pve=100*pr$sdev^2/sum(pr$sdev^2)
  C=cor(t(X)); mr=mean(abs(C[upper.tri(C)]))
  totvar=sum(matrixStats::rowVars(X))
  # how aligned is PC1 with disease? point-biserial cor(PC1, disease)
  pcb=abs(cor(pr$x[,1], as.numeric(grp=="Disease")))
  cat(sprintf("\n== %s (%d genes) ==\n",name,length(idx)))
  cat(sprintf("  PVE PC1-10: %s\n",paste(sprintf("%.1f",pve[1:10]),collapse=" ")))
  cat(sprintf("  PC1 captures %.1f%%  |  mean |gene-gene r| = %.3f  |  total var = %.0f  |  |cor(PC1,disease)| = %.2f\n",
              pve[1],mr,totvar,pcb))
}
report(idx_hvg,"HVG-2000 (variance-selected)")
report(idx_deg,"DEG-200 (|t|-selected)")
cat("\nDIAG_DONE\n")
