# Independently collect and check the saved complete permutation families.
suppressPackageStartupMessages({library(data.table);library(nlme)})
root<-Sys.getenv("MASLD_PROJECT_ROOT")
task<-file.path(root,"Analysis/SingleCell/candidates/source-diagnosis-correction-20261026")
source(file.path(root,"Analysis/Hypotheses/H59-20261015T112432Z/scripts/h59_lib.R"))
p<-fread(file.path(task,"programs/donor_program_scores_primary.tsv"))[dataset%in%FIT_DATASETS]
meta<-unique(p[,.(donor,dataset,disease_stage_coarse,stage)])
an7<-file.path(root,"Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN7/20261021T124200Z/out")
checks<-list()
for (mode in c("primary","qc","ambient","tcell20")) {
  path<-file.path(task,if(mode=="primary") "programs" else if(mode=="ambient") "programs_ambient_v2" else paste0("programs_",mode))
  r<-fread(file.path(path,"results_primary_family.tsv"))
  null<-readRDS(file.path(path,"permutation_nulls.rds"))
  stopifnot(identical(colnames(null$t),r$program),nrow(null$t)==19999L,all(is.finite(null$t)))
  fit<-copy(p)
  if(mode=="qc") {
    x<-fread(file.path(an7,"donor_scores_v2qc_cells.tsv"))
    fit<-merge(x[,.(donor,cell_type,module=as.integer(module),score)],meta,by="donor")
  }
  if(mode=="ambient") {
    x<-fread(inp("ambient_scores"))
    fit<-merge(x[is.finite(corrected),.(donor,cell_type,module=as.integer(module),score=corrected)],meta,by="donor")
    fit<-fit[cell_type!="tcells"]
  }
  if(mode=="tcell20") {
    x<-fread(file.path(an7,"donor_scores_all_cells.tsv"))
    fit<-merge(fit,x[,.(donor,cell_type,module=as.integer(module),n_cells)],by=c("donor","cell_type","module"))
    fit<-fit[cell_type=="tcells"&n_cells>=20L]
  }
  check<-rbindlist(lapply(seq_len(nrow(r)),function(i) {
    d<-prog_data(fit,r$cell_type[i],r$module[i],FIT_DATASETS)
    f<-fit_gls(d$score,d$stage,d$dataset)
    stopifnot(!anyDuplicated(d$donor),nrow(d)==r$n_fit[i],abs(f[["beta"]]-r$gls_beta[i])<1e-8,
      abs(f[["se"]]-r$gls_se[i])<1e-8)
    pval<-(1+sum(abs(null$t[,i])>=abs(f[["t"]])*(1-1e-12)))/20000
    p3<-if(all(is.finite(null$hc3_t[,i]))) (1+sum(abs(null$hc3_t[,i])>=abs(f[["t3"]])*(1-1e-12)))/20000 else NA_real_
    if(!r$stage_count_eligible[i]) {pval<-1;p3<-1}
    stopifnot(abs(pval-r$perm_p[i])<1e-12,isTRUE(all.equal(p3,r$perm_p_t3[i],tolerance=1e-12)))
    data.table(program=r$program[i],n=nrow(d),permutation_p=pval,hc3_permutation_p=p3)
  }))
  fwrite(fit,file.path(path,"fitted_program_scores.tsv"),sep="\t")
  fwrite(check,file.path(path,"independent_null_checks.tsv"),sep="\t")
  checks[[mode]]<-data.table(mode=mode,programs=nrow(r),draws=19999L,all_null_finite=TRUE)
}
fwrite(rbindlist(checks),file.path(task,"program_null_verification.tsv"),sep="\t")
writeLines(capture.output(sessionInfo()),file.path(task,"sessionInfo_program_verification.txt"))
