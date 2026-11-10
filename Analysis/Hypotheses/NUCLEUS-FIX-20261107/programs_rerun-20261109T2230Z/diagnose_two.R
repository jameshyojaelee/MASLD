# Diagnose the two programs that failed in job 22846014. Read-only on prior outputs.
suppressPackageStartupMessages({library(data.table);library(nlme)})
root <- '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
task <- file.path(root,'Analysis/Hypotheses/NUCLEUS-FIX-20261107')
here <- file.path(task,'programs_rerun-20261109T2230Z')
source(file.path(root,'Analysis/Hypotheses/H59-20261015T112432Z/scripts/h59_lib.R'))
current <- fread(file.path(task,'fixed_programs/donor_scores.tsv'))
donors <- sort(unique(current$donor)); draws <- 19999L
set.seed(20260923L); keys <- matrix(runif(draws*length(donors)),nrow=draws)
for (pr in list(c('cholangiocytes',18L), c('fibroblasts',25L))) {
  ct <- pr[1]; mod <- as.integer(pr[2])
  d <- prog_data(current,ct,mod,FIT_DATASETS)
  cat('==',ct,mod,'n=',nrow(d),'\n'); print(d[,.N,by=.(dataset,stage)])
  obs <- fit_gls(d$score,d$stage,d$dataset); print(obs)
  permuted <- perm_stage(d,keys,donors)
  tt <- vapply(seq_len(draws),function(b) fit_gls(d$score,permuted[b,],d$dataset)[['t']],0.)
  bad <- which(!is.finite(tt))
  cat('non-finite permutation t:',length(bad),'of',draws,'; first:',head(bad,10),'\n')
  if (length(bad)) {
    b <- bad[1]; df <- data.frame(y=d$score,stage=permuted[b,],dataset=factor(d$dataset))
    e <- tryCatch(gls(y~stage+dataset,data=df,weights=varIdent(form=~1|dataset),method='REML',control=GLS_CTRL),error=function(e) e)
    cat('gls error on draw',b,':',conditionMessage(e),'\n')
  }
  saveRDS(list(d=d,obs=obs,t=tt,bad=bad),file.path(here,paste0('diag_',ct,'_',mod,'.rds')))
}
