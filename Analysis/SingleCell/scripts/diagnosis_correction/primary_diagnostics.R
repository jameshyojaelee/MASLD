# Recompute H59 uncertainty and permutation checks for the corrected donor roster.
suppressPackageStartupMessages({library(data.table);library(nlme);library(parallel)})
setDTthreads(1)
root<-Sys.getenv("MASLD_PROJECT_ROOT")
task<-file.path(root,"Analysis/SingleCell/candidates/source-diagnosis-correction-20261026")
source(file.path(root,"Analysis/Hypotheses/H59-20261015T112432Z/scripts/h59_lib.R"))
out<-file.path(task,"programs_diagnostics");stopifnot(!dir.exists(out));dir.create(out)
r<-fread(file.path(task,"programs/results_primary_family.tsv"))
p<-fread(file.path(task,"programs/donor_program_scores_primary.tsv"))[dataset%in%FIT_DATASETS]
null<-readRDS(file.path(task,"programs/permutation_nulls.rds"))
stopifnot(nrow(null$t)==19999L,ncol(null$t)==117L,all(is.finite(null$t)),identical(colnames(null$t),r$program))
set.seed(20260923L)
K<-matrix(runif(19999L*length(null$donor_order)),nrow=19999L)
checks<-rbindlist(lapply(seq_len(nrow(r)),function(i) {
  d<-prog_data(p,r$cell_type[i],r$module[i],FIT_DATASETS)
  id<-perm_stage(d,matrix(seq_along(null$donor_order),nrow=1L),null$donor_order)
  perm<-perm_stage(d,K[1:100,,drop=FALSE],null$donor_order)
  stopifnot(identical(as.numeric(id[1,]),as.numeric(d$stage)),!anyDuplicated(d$donor))
  for (ds in FIT_DATASETS) for (stage in 0:2) {
    idx<-which(d$dataset==ds)
    stopifnot(all(rowSums(perm[,idx,drop=FALSE]==stage)==sum(d$stage[idx]==stage)))
  }
  X<-model.matrix(~stage+dataset,data.frame(stage=d$stage,dataset=factor(d$dataset)))
  w<-ifelse(d$dataset==FIT_DATASETS[1],r$sigma_GSE202379[i],r$sigma_GSE244832[i])^-2
  independent<-solve(crossprod(X*w,X),crossprod(X*w,d$score))[2]
  stopifnot(abs(independent-r$gls_beta[i])<1e-8)
  sg<-sigma_cells(d);sx<-sxx_by(d);se<-1/sqrt(sum(sx/sg^2))
  data.table(program=r$program[i],n=nrow(d),independent_weighted_beta=independent,
    sd_ratio=sg[2]/sg[1],mde_lone_call=(qnorm(1-.05/117/2)+qnorm(.8))*se,
    GLS_nominal_null_rate=mean(2*pt(abs(null$t[,i]),nrow(d)-3,lower.tail=FALSE)<.05))
}))
fwrite(checks,file.path(out,"model_checks.tsv"),sep="\t")
reconstructed<-vapply(seq_len(nrow(r)),function(i) {
  if (!r$stage_count_eligible[i]) return(1)
  (1+sum(abs(null$t[,i])>=abs(r$gls_t[i])*(1-1e-12)))/20000
},0)
stopifnot(max(abs(reconstructed-r$perm_p))<1e-12,max(abs(bh117(reconstructed)-r$perm_q))<1e-12)
# Each null draw is compared with the other draws, including ties correctly.
PN<-sapply(seq_len(ncol(null$t)),function(j) rank(-abs(null$t[,j]),ties.method="max")/nrow(null$t))
PN[,!r$stage_count_eligible]<-1
calls<-apply(PN[1:2000,,drop=FALSE],1,function(x) sum(bh117(x)<.05))
nc1<-data.table(n_draws=2000L,fraction_any_call=mean(calls>0),max_calls=max(calls))
fwrite(nc1,file.path(out,"NC1_null_family.tsv"),sep="\t")
stopifnot(nc1$fraction_any_call<=.10)
# Same known-truth simulation as H59, with actual corrected hepatocyte stages.
hep<-prog_data(p,"hepatocytes",20L,FIT_DATASETS)
sdv<-ifelse(hep$dataset==FIT_DATASETS[1],.33,.95)
set.seed(20260923L)
sim<-rbindlist(lapply(c(0,.30),function(b) {
  ys<-lapply(seq_len(2000L),function(k) b*hep$stage+rnorm(nrow(hep))*sdv)
  z<-do.call(rbind,mclapply(ys,function(y) fit_gls(y,hep$stage,hep$dataset),mc.cores=16L))
  data.table(true_beta=b,n_draws=2000L,n_failed=sum(!is.finite(z[,"beta"])),
    coverage95=mean(abs(z[,"beta"]-b)<=qt(.975,nrow(hep)-3)*z[,"se"]),
    bias=mean(z[,"beta"]-b))
}))
fwrite(sim,file.path(out,"PC3_known_truth.tsv"),sep="\t")
# Native fixed alternatives, not effects chosen from the corrected results.
detection<-rbindlist(lapply(seq_len(nrow(NAMED)),function(i) {
  d<-prog_data(p,NAMED$cell_type[i],NAMED$module[i],FIT_DATASETS)
  perm<-perm_stage(d,K[1:1000,,drop=FALSE],null$donor_order)
  crit<-r[program==paste(NAMED$cell_type[i],NAMED$module[i],sep="_"),null_critical_abs_t]
  rbindlist(lapply(c(.20,.30,NAMED$printed_beta[i]),function(delta) {
    tt<-unlist(mclapply(seq_len(1000L),function(k) fit_gls(d$score+delta*perm[k,],perm[k,],d$dataset)[["t"]],mc.cores=16L))
    stopifnot(all(is.finite(tt)))
    data.table(program=paste(NAMED$cell_type[i],NAMED$module[i],sep="_"),delta=delta,n_draws=1000L,detection=mean(abs(tt)>crit))
  }))
}))
fwrite(detection,file.path(out,"PC4_detection.tsv"),sep="\t")
deletion<-rbindlist(lapply(c(8L,20L),function(m) {
  d<-prog_data(p,"hepatocytes",m,FIT_DATASETS)
  rbindlist(lapply(seq_len(nrow(d)),function(i) {
    x<-d[-i];f<-fit_gls(x$score,x$stage,x$dataset)
    data.table(module=m,deleted_donor=d$donor[i],beta=f[["beta"]],se=f[["se"]])
  }))
}))
fwrite(deletion,file.path(out,"single_donor_deletions.tsv"),sep="\t")
f<-fread(file.path(task,"metadata/documented_F_source.tsv"))[!is.na(F_stage_documented)]
fwrite(f,file.path(out,"documented_fstage_donor_map.tsv"),sep="\t")
fwrite(data.table(n_documented=nrow(f),n_F0=sum(f$F_stage_documented==0),n_F1=sum(f$F_stage_documented==1),
  n_F0_F1=sum(f$F_stage_documented%in%0:1),minimum=15,decision=ifelse(sum(f$F_stage_documented%in%0:1)<15,"STOP_fewer_than_15_F0_F1_donors","eligible")),file.path(out,"H7_stop_rule.tsv"),sep="\t")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
stopifnot(all(sim$n_failed==0),min(sim$coverage95)>=.92)
writeLines("Independent weighted coefficients, permutation p-values, BH117, identity/category preservation, null-family calibration and model-CI simulation passed.",file.path(out,"COMPLETED"))
