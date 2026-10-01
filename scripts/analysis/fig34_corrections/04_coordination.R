#!/usr/bin/env Rscript
# Retrospective uncertainty for fixed communities. Participants are resampled,
# never edges. Gene calibration and membership remain frozen within each cohort.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(Matrix)})
out <- commandArgs(TRUE)[1]; stopifnot(!dir.exists(out)); dir.create(out, recursive=TRUE)
root <- getwd(); sd <- file.path(root,"scripts/manuscript/program_observability_map")
for (f in c("config.R","analysis_lib.R","calibration_lib.R","systems_contract.R","t4_interaction_lib.R")) source(file.path(sd,f))
set.seed(20260929); B <- 1000L
td <- file.path(workstream_root,"systems/T4_interaction")
comm <- fread(file.path(td,"program_communities.tsv")); features <- comm$feature_id
stopifnot(nrow(comm)==113L,uniqueN(comm$community)==4L)
assert_frozen_programs(program_registry,program_membership)
w <- t4_membership_weights(fread(program_membership),features)
pairs <- t4_pair_table(features,w)
lab <- fread(file.path(root,"Analysis/SingleCell/results_gpu_v2/program_label_vs_content/plc-20260907T105813Z/program_label_vs_content.tsv"))
ct <- setNames(lab$content_lineage,lab$program_uid); cm <- setNames(comm$community,comm$feature_id)
pairs[, `:=`(content_a=ct[program_a],content_b=ct[program_b],community_a=cm[program_a],community_b=cm[program_b])]
pairs[, cross_content := content_a != content_b]
fwrite(pairs,file.path(out,"fixed_pair_classes.tsv"),sep="\t")
coord <- fread(file.path(td,"coordinated_edges.tsv"))
key <- function(a,b) paste(pmin(a,b),pmax(a,b),sep="||")
coord[, ekey:=key(program_a,program_b)]; pairs[, ekey:=key(program_a,program_b)]
counts <- pairs[, .(n_possible=.N,n_discovery_supported=sum(ekey %in% coord$ekey)),by=.(community_a,community_b,cross_content)]
fwrite(counts,file.path(out,"edge_census.tsv"),sep="\t")
pick <- which(pairs$community_a==3 & pairs$community_b==3)
stopifnot(length(pick)==105L,sum(pairs$cross_content[pick])==56L)
anno <- fread(gene_annotation,select=c("gene_id","gene_name")); blocks <- list()
for (hold in c(FALSE,TRUE)) {
  dge <- readRDS(if(hold) holdout_dge else nonholdout_dge)
  md <- as.data.table(readRDS(if(hold) holdout_meta else nonholdout_meta))
  if(hold) {
    cw <- fread(gse193066_crosswalk)
    md <- merge(md,cw[,.(sample_id=run_id,participant_token,is_first_biopsy)],by="sample_id")
    md <- md[is_first_biopsy==TRUE]; stopifnot(!anyDuplicated(md$participant_token))
  } else md <- md[dataset %in% discovery_cohorts]
  md <- md[!is.na(fibrosis_stage)&!is.na(nas_score)]
  sm <- collapse_symbols(edgeR::cpm(dge,log=TRUE,prior.count=1),rownames(dge),anno)
  comp <- fread(if(hold) holdout_composition else nonholdout_composition)
  for(co in unique(md$dataset)) {
    b <- t4_build_block(sm,md[dataset==co],w,features,.8)
    op <- t4_disjoint_operators(b,pairs,features,.8)
    P <- t(b$Z)%*%b$W; A <- P[,pairs$i[pick],drop=FALSE]; D <- P[,pairs$j[pick],drop=FALSE]
    for(j in seq_along(pick)) {
      k <- match(pick[j],op$pair_index)
      if(!is.na(k)) {
        A[,j] <- (A[,j]-as.vector(t(b$Z)%*%op$Va[,k]))/(1-op$shared_a[k])
        D[,j] <- (D[,j]-as.vector(t(b$Z)%*%op$Vb[,k]))/(1-op$shared_b[k])
        if(!op$eligible[k]) { A[,j]<-NA;D[,j]<-NA }
      }
    }
    cf <- as.matrix(comp[match(b$samples,sample_id),setdiff(names(comp),c("sample_id","dataset")),with=FALSE])
    stopifnot(all(is.finite(cf)))
    clr <- log(pmax(cf,1e-6));clr <- clr-rowMeans(clr)
    cp <- svd(scale(clr,scale=FALSE),nu=3,nv=0)$u[,1:3,drop=FALSE]
    blocks[[co]] <- list(P=P,A=A,D=D,X=b$X,XC=cbind(b$X,cp),n=b$n,df=b$df)
    # Independently constructed score-space residualization must reproduce
    # the original gene-space procedure on the observed participants.
    ob <- t4_one_pass(list(b),pairs,list(op))
    attr(blocks[[co]],"expected") <- ob$z_partial_k1[pick]
  }
  rm(dge,sm);gc()
}
one <- function(b,ii=seq_len(b$n)) {
  np <- ncol(b$A); answer <- matrix(NA_real_,np,5,dimnames=list(NULL,c("unadjusted","histology","histology_PC1","histology_PC3","histology_composition3")))
  P <- b$P[ii,,drop=FALSE]; A <- b$A[ii,,drop=FALSE];D <- b$D[ii,,drop=FALSE]
  answer[,1] <- t4_atanh(t4_colwise_cor(A,D))
  for(comp in c(FALSE,TRUE)) {
    X <- if(comp)b$XC[ii,,drop=FALSE] else b$X[ii,,drop=FALSE]
    qrX<-qr(X); if(qrX$rank<ncol(X))next
    h<-rowSums(qr.Q(qrX)^2);if(any(h>=1-1e-10))next
    res <- function(y) {
      ok<-colSums(is.finite(y))==nrow(y);z<-matrix(NA_real_,nrow(y),ncol(y))
      if(any(ok))z[,ok]<-qr.resid(qrX,y[,ok,drop=FALSE])/sqrt(1-h)
      z
    }
    a<-res(A);d<-res(D)
    if(comp) {answer[,5]<-t4_atanh(t4_colwise_cor(a,d));next}
    answer[,2]<-t4_atanh(t4_colwise_cor(a,d))
    ps<-scale(res(P));ps[!is.finite(ps)]<-0
    U<-svd(ps,nu=3,nv=0)$u
    for(k in c(1,3)) {
      u<-U[,seq_len(k),drop=FALSE]
      answer[,if(k==1)3 else 4]<-t4_atanh(t4_colwise_cor(a-u%*%crossprod(u,a),d-u%*%crossprod(u,d)))
    }
  }
  answer
}
checks<-rbindlist(lapply(names(blocks),function(co) data.table(dataset=co,n=blocks[[co]]$n,max_abs_z_difference=max(abs(one(blocks[[co]])[,3]-attr(blocks[[co]],"expected")),na.rm=TRUE))))
stopifnot(all(checks$max_abs_z_difference<1e-8));fwrite(checks,file.path(out,"original_statistic_reproduction.tsv"),sep="\t")
observed<-lapply(blocks,one); cohorts<-names(blocks)
summary_stats<-function(ans) {
  rows<-list()
  groups<-c(setNames(lapply(cohorts,identity),cohorts),list(discovery=intersect(discovery_cohorts,cohorts)))
  for(g in names(groups)) {
    cs<-groups[[g]];wt<-sapply(blocks[cs],function(b)b$df-3)
    av<-Reduce(`+`,Map(function(a,w)a*w,ans[cs],wt))/sum(wt)
    for(cross in c(FALSE,TRUE)) for(j in seq_len(ncol(av))) {
      v<-av[pairs$cross_content[pick]==cross,j]
      rows[[length(rows)+1L]]<-data.table(dataset=g,cross_content=cross,model=colnames(av)[j],mean_fisher_z=if(any(is.finite(v)))mean(v,na.rm=TRUE) else NA_real_,positive_fraction=if(any(is.finite(v)))mean(v>0,na.rm=TRUE) else NA_real_,n_eligible_pairs=sum(is.finite(v)),n=sum(sapply(blocks[cs],`[[`,"n")))
    }
  };rbindlist(rows)
}
obs<-summary_stats(observed); draws<-vector("list",B)
for(i in seq_len(B)) {
  aa<-lapply(blocks,function(b)one(b,sample.int(b$n,b$n,replace=TRUE)))
  draws[[i]]<-summary_stats(aa)[,replicate:=i]
  if(i%%100==0)message("Participant bootstrap ",i,"/",B)
}
dd<-rbindlist(draws);lim<-dd[,.(z_low=quantile(mean_fisher_z,.025,na.rm=TRUE),z_high=quantile(mean_fisher_z,.975,na.rm=TRUE),positive_low=quantile(positive_fraction,.025,na.rm=TRUE),positive_high=quantile(positive_fraction,.975,na.rm=TRUE),n_bootstrap=sum(is.finite(mean_fisher_z))),by=.(dataset,cross_content,model)]
rr<-merge(obs,lim,by=c("dataset","cross_content","model"));rr[,`:=`(mean_r=tanh(mean_fisher_z),r_low=tanh(z_low),r_high=tanh(z_high),seed=20260929,interval="pointwise_conditional_on_frozen_scores",retrospective=TRUE)]
fwrite(rr,file.path(out,"participant_bootstrap_summary.tsv"),sep="\t")
edge<-rbindlist(lapply(cohorts,function(co){ z<-as.data.table(observed[[co]]);cbind(data.table(dataset=co,program_a=pairs$program_a[pick],program_b=pairs$program_b[pick],cross_content=pairs$cross_content[pick]),z)}))
fwrite(edge,file.path(out,"C3_cohort_pair_fisher_z.tsv"),sep="\t")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
writeLines("Retrospective fixed C3 summaries; additive categorical fibrosis/NAS plus sex. Composition sensitivity adds three cohort-fitted CLR components; PC1/PC3 are cohort-fitted residual-score components. Intervals resample participants within cohort; no edge-level inference. Bootstrap conditional on observed gene calibration and composition basis. No cell-to-cell regulatory interpretation.",file.path(out,"analysis_definition.txt"))
writeLines("Complete",file.path(out,"COMPLETE"))
