#!/usr/bin/env Rscript
# T2 held-out / within-cohort reproducibility of the F1 rare substate (review fix:
# the substate was discovered on pooled F1, so test that it reproduces when each
# cohort is analysed INDEPENDENTLY — a single-cohort artifact would fail this).
# For each F1 cohort with ≥20 donors: own DV-feature PCA → bimodality (BIC+dip+sep).
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
INT <- "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration"
dge <- readRDS(file.path(INT, "merged_dge.rds"))
mm  <- as.data.table(readRDS(file.path(INT, "meta_matched.rds")))
meta <- mm[match(colnames(dge), sample_id)]; meta[, sid := colnames(dge)]
meta[, fstage := suppressWarnings(as.integer(gsub("[^0-9]","",as.character(fibrosis_stage))))]
dv <- fread("RNA-seq/results/heterogeneity_program/phase0/dv_atlas_columns.tsv")
dvg <- dv[dv_sig == TRUE | abs(dv_t) > 3, gene_ensembl]
rowVars <- function(m){mu<-rowMeans(m);rowSums((m-mu)^2)/(ncol(m)-1)}
fit_gmm_2 <- function(x,mi=200,tol=1e-6){n<-length(x);o<-sort(x);mu1<-mean(o[1:floor(n/2)]);mu2<-mean(o[(floor(n/2)+1):n])
 s1<-max(sd(o[1:floor(n/2)]),1e-3);s2<-max(sd(o[(floor(n/2)+1):n]),1e-3);p1<-.5;lo<--Inf
 for(i in 1:mi){d1<-p1*dnorm(x,mu1,s1);d2<-(1-p1)*dnorm(x,mu2,s2);g<-d1/(d1+d2+1e-300);p1<-mean(g)
  mu1<-sum(g*x)/sum(g);mu2<-sum((1-g)*x)/sum(1-g);s1<-max(sqrt(sum(g*(x-mu1)^2)/sum(g)),1e-3);s2<-max(sqrt(sum((1-g)*(x-mu2)^2)/sum(1-g)),1e-3)
  ll<-sum(log(p1*dnorm(x,mu1,s1)+(1-p1)*dnorm(x,mu2,s2)+1e-300));if(abs(ll-lo)<tol)break;lo<-ll}
 list(sep=abs(mu1-mu2)/sqrt((s1^2+s2^2)/2),pi=min(p1,1-p1),ll=ll)}
dipp <- function(x,nmc=499){d<-function(z){g<-seq(min(z),max(z),length.out=200);max(abs(ecdf(z)(g)-pnorm(g,mean(z),sd(z))))}
 o<-d(x);(sum(replicate(nmc,d(rnorm(length(x),mean(x),sd(x))))>=o)+1)/(nmc+1)}

f1 <- which(meta$fstage == 1)
res <- list()
for (co in names(which(table(meta$dataset[f1]) >= 20))) {
  idx <- f1[meta$dataset[f1] == co]
  d <- dge[, idx, keep.lib.sizes=FALSE]; d <- d[rownames(d) %in% dvg,,keep.lib.sizes=FALSE]; d <- calcNormFactors(d)
  E <- voom(d, design=NULL)$E; E <- E[order(rowVars(E),decreasing=TRUE)[1:min(1000,nrow(E))],]
  sc <- prcomp(t(E), scale.=TRUE)$x[,1]
  g <- fit_gmm_2(sc); n <- length(sc); ll1 <- sum(dnorm(sc,mean(sc),sd(sc),log=TRUE))
  bic_diff <- (-2*ll1+2*log(n)) - (-2*g$ll+5*log(n))
  dp <- dipp(sc)
  res[[co]] <- data.table(cohort=co, n=n, bic_diff=round(bic_diff,1), dip_p=signif(dp,2),
                          sep=round(g$sep,2), pi_rare=round(g$pi,3),
                          bimodal=bic_diff>0 & dp<0.05 & g$sep>1.5 & g$pi>=0.05)
}
r <- rbindlist(res); print(r)
cat(sprintf("\nF1 substate reproduces independently in %d / %d cohorts (≥2 = REPRODUCIBLE, not a pooling artifact)\n",
            r[bimodal==TRUE,.N], nrow(r)))
