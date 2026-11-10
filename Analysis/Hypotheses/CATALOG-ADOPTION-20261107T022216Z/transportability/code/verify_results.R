#!/usr/bin/env Rscript
# Independent base-R restricted-likelihood and complete-family checks.
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
OUT <- '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/transportability'
T<-fread(file.path(OUT,'results/gene_transportability.tsv.gz'))
L<-fread(file.path(OUT,'cohort_effects_long.tsv.gz'))
stopifnot(nrow(T)==24115L,nrow(L)==120575L,sum(T$primary_call)==1854L,
 sum(L$native_testable)==112052L,sum(!L$native_testable)==8523L,
 all(L[native_testable==FALSE,dataset]=='GSE126848'),all(is.na(L[native_testable==FALSE,moderated_SE])),
 max(abs(p.adjust(T$p,'BH',n=24115L)-T$BH_q))<1e-12)
idx<-unique(round(seq(1,nrow(T),length.out=25)))
check<-rbindlist(lapply(idx,function(i){
 g<-T$gene[i];z<-L[gene_id_versioned==g & native_testable==TRUE]; y<-z$logFC;v<-z$moderated_SE^2;k<-length(y)
 objective<-function(tau){w<-1/(v+tau);b<-sum(w*y)/sum(w);.5*(sum(log(v+tau))+log(sum(w))+sum(w*(y-b)^2))}
 upper<-max(1,var(y)*20);tau<-optimize(objective,c(0,upper),tol=1e-10)$minimum
 if(objective(0)<=objective(tau))tau<-0
 w<-1/(v+tau);b<-sum(w*y)/sum(w);scale<-max(1,sum(w*(y-b)^2)/(k-1));se<-sqrt(scale/sum(w))
 data.table(gene=g,k=k,independent_tau2=tau,reported_tau2=T$tau2[i],tau_difference=abs(tau-T$tau2[i]),
 independent_beta=b,reported_beta=T$beta[i],beta_difference=abs(b-T$beta[i]),
 independent_SE=se,reported_SE=T$SE[i],SE_difference=abs(se-T$SE[i]))
}))
stopifnot(all(check$tau_difference<1e-4*(1+check$reported_tau2)),all(check$beta_difference<1e-4),all(check$SE_difference<1e-4))
fwrite(check,file.path(OUT,'independent_REML_gene_checks.tsv'),sep='\t')
write_json(list(pass=TRUE,n_checked_genes=nrow(check),BH_full_family=24115L,fixed_calls=1854L,
 max_tau2_difference=max(check$tau_difference),max_beta_difference=max(check$beta_difference),
 max_SE_difference=max(check$SE_difference)),file.path(OUT,'INDEPENDENT_RESULTS_PASS.json'),auto_unbox=TRUE,pretty=TRUE,digits=NA)
cat('PASS: independent restricted likelihood on 25 fixed-position genes; coverage and full-family BH\n')
