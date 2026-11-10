#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(metafor);library(jsonlite);library(parallel)})
ROOT <- '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
OUT <- file.path(ROOT,'Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/transportability')
setDTthreads(1);options(digits=17);set.seed(20261107)
stopifnot(file.exists(file.path(OUT,'SE_COVERAGE_PASS.json')))
RESULT <- file.path(OUT,'results');stopifnot(!dir.exists(RESULT));dir.create(RESULT)
M<-fread(file.path(OUT,'input_manifest.tsv'));P<-setNames(M$path,M$input)
sha<-function(p)sub(' .*','',system2('sha256sum',shQuote(p),stdout=TRUE))
stopifnot(all(M$sha256==vapply(M$path,sha,'')))
spec_sha<-strsplit(readLines(file.path(OUT,'ANALYSIS_SPEC.sha256')),'  ',fixed=TRUE)[[1]][1]
stopifnot(identical(sha(file.path(OUT,'ANALYSIS_SPEC.md')),spec_sha))
L<-fread(file.path(OUT,'cohort_effects_long.tsv.gz'))
fit_one<-function(y,se){
 k<-length(y);warn<-character()
 fit<-tryCatch(withCallingHandlers(rma.uni(yi=y,sei=se,method='REML',test='adhoc',control=list(maxiter=1000,stepadj=.5)),
  warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 if(inherits(fit,'error'))return(list(k=k,status='fit_failed',failure=conditionMessage(fit),p=1))
 w<-1/(se^2+fit$tau2);b<-sum(w*y)/sum(w);scale<-max(1,sum(w*(y-b)^2)/(k-1));s<-sqrt(scale/sum(w))
 p<-2*pt(-abs(b/s),df=k-1)
 stopifnot(abs(b-as.numeric(fit$b))<1e-8,abs(s-as.numeric(fit$se))<1e-8,abs(p-as.numeric(fit$pval))<1e-8)
 piok<-k==5L && is.finite(fit$tau2) && fit$tau2>0
 piwidth<-if(piok)qt(.975,k-2)*sqrt(fit$tau2+s^2)else NA_real_
 list(k=k,status='estimated',failure='',beta=b,SE=s,ci_low=b-qt(.975,k-1)*s,ci_high=b+qt(.975,k-1)*s,p=p,
  tau2=fit$tau2,I2=fit$I2,Q=fit$QE,Q_p=fit$QEp,kh_scale=scale,max_weight=max(w)/sum(w),
  pi_low=if(piok)b-piwidth else NA_real_,pi_high=if(piok)b+piwidth else NA_real_,
  pi_state=if(piok)'exploratory_k5_positive_tau2_exchangeability_assumed'else'not_reported_k_less_5_or_tau2_boundary',
  n_positive=sum(y>0),n_negative=sum(y<0),warnings=paste(unique(warn),collapse=';'))
}
# Planted summary-data controls exercise fixed scientific invariants.
f0<-fit_one(rep(.6,5),rep(.1,5));stopifnot(abs(f0$beta-.6)<1e-10,f0$kh_scale>=1,f0$SE>=.1/sqrt(5),f0$tau2<1e-8)
fp<-fit_one(c(.4,.5,.6,.7,.8),rep(.1,5));fm<-fit_one(-c(.4,.5,.6,.7,.8),rep(.1,5))
f4<-fit_one(c(.4,.5,.6,.7),rep(.1,4));stopifnot(abs(fp$beta+fm$beta)<1e-10,abs(fp$p-fm$p)<1e-10,is.na(f4$pi_low))
write_json(list(equal_effect=TRUE,KH_floor=TRUE,sign_reversal=TRUE,missing_cohort_PI_suppressed=TRUE),file.path(RESULT,'planted_controls.json'),auto_unbox=TRUE,pretty=TRUE)
genes<-fread(P[['pooled']])$gene
X<-split(L[native_testable==TRUE,.(gene_id_versioned,logFC,moderated_SE)],by='gene_id_versioned',keep.by=FALSE)
stopifnot(length(X)==24115L,setequal(names(X),genes))
cores<-min(2L,as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','1')))
ans<-vector('list',length(genes));chunk<-1000L
for(start in seq(1L,length(genes),by=chunk)){
 idx<-start:min(length(genes),start+chunk-1L)
 ans[idx]<-mclapply(genes[idx],function(g){z<-X[[g]];c(list(gene=g),fit_one(z$logFC,z$moderated_SE))},mc.cores=cores,mc.preschedule=TRUE)
 cat('REML completed',max(idx),'/',length(genes),'\n');flush.console()
}
T<-rbindlist(ans,fill=TRUE);stopifnot(nrow(T)==24115L,!anyDuplicated(T$gene))
saveRDS(T,file.path(RESULT,'random_effects_before_annotations.rds'))
T[is.na(p)|!is.finite(p),`:=`(p=1,status='fit_failed')]
T[,BH_q:=p.adjust(p,'BH',n=24115L)]
A<-fread(P[['pooled']])[,.(gene,symbol,pooled_logFC=logFC,pooled_SE=SE,pooled_q=padj,primary_call=padj<.05 & abs(logFC)>.5)]
T<-merge(A,T,by='gene',all.x=TRUE,sort=FALSE)
C<-fread(P[['classes']]);setnames(C,'gene_id_versioned','gene')
setnames(C,setdiff(names(C),'gene'),paste0('existing_cohort_calls_',setdiff(names(C),'gene')))
T<-merge(T,C,by='gene',all.x=TRUE,sort=FALSE)
H<-fread(P[['ancestry']])[,.(gene,ancestry_pc_status=M1_status,ancestry_pc_sensitive_descriptive=D1_rule_M1,
 ancestry_pc_adjusted_logFC=M1_logFC,ancestry_pc_adjusted_q=M1_padj,ancestry_EUR_status=ME_status)]
T<-merge(T,H,by='gene',all.x=TRUE,sort=FALSE)
SS<-fread(P[['ss_sidecar']])[,.(gene=gene_id_versioned,steatosite_stage_state,
 ss_discovery_F3F4_logFC=D0_F3F4_log2FC,ss_F3F4_logFC=steatosite_F3F4_log2FC,ss_F3F4_q=steatosite_F3F4_BH_q,
 ss_F3F4_MDE=steatosite_F3F4_MDE_g)]
ss_blank_native_ids<-sum(is.na(SS$gene)|SS$gene=='')
SS<-SS[gene %in% genes]
stopifnot(!anyDuplicated(SS$gene))
T<-merge(T,SS,by='gene',all.x=TRUE,sort=FALSE)
T[is.na(steatosite_stage_state),steatosite_stage_state:='not_joined_to_existing_Catalog_stage_sidecar']
SG<-fread(P[['ss_gene']])[coef=='binF3F4',.(gene,ss_F3F4_SE=SE)]
T<-merge(T,SG,by='gene',all.x=TRUE,sort=FALSE)
T[,`:=`(same_sign_as_pooled=sign(beta)==sign(pooled_logFC),RE_q_supported=BH_q<.05)]
T[,sensitivity_state:=fifelse(status!='estimated','unavailable',
 fifelse(BH_q<.05 & same_sign_as_pooled,'same_direction_RE_support',
 fifelse(BH_q<.05 & !same_sign_as_pooled,'opposite_direction_RE_support','RE_uncertain')))]
T[,`:=`(analysis_role='retrospective_random_effects_sensitivity_not_primary_call_replacement',
 SE_type='empirical_Bayes_moderated_input_SE',BH_family_size=24115L,
 new_blind_confirmation='not_started_no_qualified_unexposed_cohort_established')]
stopifnot(nrow(T)==24115L,sum(T$primary_call)==1854L,!anyNA(T$status),all(T$k %in% 4:5))
fwrite(T,file.path(RESULT,'gene_transportability.tsv.gz'),sep='\t')
fwrite(T[primary_call==TRUE],file.path(RESULT,'fixed_1854_transportability.tsv.gz'),sep='\t')
side<-T[,.(gene,symbol,primary_call,transport_k=k,transport_RE_beta=beta,transport_RE_CI_low=ci_low,
 transport_RE_CI_high=ci_high,transport_RE_BH_q=BH_q,transport_tau2=tau2,transport_I2=I2,
 transport_state=sensitivity_state,ancestry_pc_status,ancestry_pc_sensitive_descriptive,steatosite_stage_state,
 ss_F3F4_logFC,ss_F3F4_q,analysis_role)]
fwrite(side,file.path(RESULT,'catalog_transportability_proposed_sidecar.tsv.gz'),sep='\t')
cohort<-merge(L,A[,.(gene_id_versioned=gene,primary_call,pooled_logFC)],by='gene_id_versioned')
cs<-cohort[primary_call==TRUE,.(fixed_calls=.N,testable=sum(native_testable),
 same_direction=sum(native_testable & sign(logFC)==sign(pooled_logFC),na.rm=TRUE),
 called_same_direction=sum(canonical_cohort_call & sign(logFC)==sign(pooled_logFC),na.rm=TRUE),
 called_opposite_direction=sum(canonical_cohort_call & sign(logFC)!=sign(pooled_logFC),na.rm=TRUE)),by=dataset]
fwrite(cs,file.path(RESULT,'fixed_signature_cohort_summary.tsv'),sep='\t')
fwrite(fread(P[['ss_family']]),file.path(RESULT,'existing_steatosite_grouped_stage_family.tsv'),sep='\t')
summary<-list(n_genes=nrow(T),n_fixed_primary_calls=sum(T$primary_call),fit_states=as.list(table(T$status)),
 ss_catalog_blank_native_ids_not_joined=ss_blank_native_ids,
 fixed_sensitivity_states=as.list(table(T[primary_call==TRUE,sensitivity_state])),
 all_RE_BH_calls=sum(T$BH_q<.05),fixed_RE_BH_same_direction=sum(T$primary_call & T$BH_q<.05 & T$same_sign_as_pooled),
 fixed_Q_nominal_heterogeneity=sum(T$primary_call & T$Q_p<.05),fixed_median_I2=median(T[primary_call==TRUE,I2]),
 k_counts=as.list(table(T$k)),PI_role='exploratory conditional interval only, k5 and positive tau2; no calibrated external prediction',
 new_blind_confirmation='not_started_no_qualified_unexposed_cohort_established')
write_json(summary,file.path(RESULT,'summary.json'),auto_unbox=TRUE,pretty=TRUE,digits=NA)
writeLines(capture.output(sessionInfo()),file.path(RESULT,'sessionInfo.txt'))
write_json(list(complete=TRUE,n_genes=nrow(T),input_manifest_sha256=sha(file.path(OUT,'input_manifest.tsv'))),file.path(RESULT,'COMPLETE.json'),auto_unbox=TRUE,pretty=TRUE)
print(summary)
