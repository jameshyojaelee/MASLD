# Rerun of scripts/refit_fixed_programs.R for the two programs whose job 22846014
# workers failed (cholangiocytes_18, fibroblasts_25), then assemble all 117.
# Root cause: on a few permuted stage vectors nlme::gls (nlminb) stops with
# "false convergence (8)"; fit_gls returns NA and stopifnot(all(is.finite(null)))
# errors. Fix: same REML GLS model (stage + dataset, per-dataset variance), refit
# a failed fit once with the optim optimizer. A fit that still fails counts as at
# least as extreme as observed (conservative). Memberships, weights, donor
# scores, eligibility, seed, keys and draws are unchanged. The 115 completed
# program caches are reused as-is (all their fits converged with nlminb).
suppressPackageStartupMessages({library(data.table);library(nlme)})
root <- '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
task <- file.path(root,'Analysis/Hypotheses/NUCLEUS-FIX-20261107')
old_out <- file.path(task,'fixed_programs')
out <- file.path(task,'programs_rerun-20261109T2230Z')
source(file.path(root,'Analysis/Hypotheses/H59-20261015T112432Z/scripts/h59_lib.R'))
stopifnot(!file.exists(file.path(out,'results.tsv')))
reg <- fread(inp('registry'))[,.(cell_type,module=as.integer(module),stability_median)]
reg <- unique(reg);setorder(reg,cell_type,module);stopifnot(nrow(reg)==117)
# Donor scores written by job 22846014 from the unique-count projection (unchanged inputs).
current <- fread(file.path(old_out,'donor_scores.tsv'))
previous <- fread(file.path(old_out,'historical_projection_donor_scores.tsv'))
donors <- sort(unique(current$donor)); draws <- 19999L
set.seed(20260923L);keys <- matrix(runif(draws*length(donors)),nrow=draws)

GLS_OPTIM <- glsControl(maxIter=200,msMaxIter=200,opt='optim')
fit_gls_fallback <- function(y,stage,dataset) {
  r <- fit_gls(y,stage,dataset)
  if (is.finite(r[['t']])) return(c(r,fallback=0))
  df <- data.frame(y=y,stage=stage,dataset=factor(dataset))
  f <- tryCatch(gls(y~stage+dataset,data=df,weights=varIdent(form=~1|dataset),
                    method='REML',control=GLS_OPTIM),error=function(e) NULL)
  if (is.null(f)) return(c(beta=NA,se=NA,t=NA,t3=NA,s1=NA,s2=NA,fallback=2))
  beta <- unname(coef(f)['stage']); se <- sqrt(vcov(f)['stage','stage'])
  c(beta=beta,se=se,t=beta/se,t3=NA,s1=NA,s2=NA,fallback=1)
}

one <- function(i) {
  ct <- reg$cell_type[i];mod <- reg$module[i]
  cache <- file.path(old_out,paste0(ct,'_',mod,'.rds'))
  if(file.exists(cache)) { r <- readRDS(cache); r$row[,`:=`(n_optim_fallback=0L,n_failed_draws=0L)]; return(r) }
  d <- prog_data(current,ct,mod,FIT_DATASETS)
  old <- prog_data(previous,ct,mod,FIT_DATASETS)
  eligible <- all(vapply(FIT_DATASETS,function(ds) {
    tab<-table(d[dataset==ds,stage]);sum(tab)>=8 && sum(tab>=2)>=2
  },TRUE))
  observed <- if(eligible) fit_gls(d$score,d$stage,d$dataset) else rep(NA_real_,6)
  oldfit <- if(nrow(old)>5 && uniqueN(old$dataset)==2) fit_gls(old$score,old$stage,old$dataset) else rep(NA_real_,6)
  null <- rep(NA_real_,draws); p <- 1; nfb <- 0L; nfail <- 0L
  if(eligible && is.finite(observed[3])) {
    permuted <- perm_stage(d,keys,donors)
    fits <- vapply(seq_len(draws),function(b) fit_gls_fallback(d$score,permuted[b,],d$dataset)[c('t','fallback')],c(t=0,fallback=0))
    null <- fits['t',]; nfb <- sum(fits['fallback',]==1); nfail <- sum(fits['fallback',]==2)
    null[!is.finite(null)] <- Inf   # conservative: unresolved draw counts as extreme
    p <- (1+sum(abs(null)>=abs(observed[3])*(1-1e-12)))/(draws+1)
  } else eligible <- FALSE
  row <- data.table(cell_type=ct,module=mod,n=nrow(d),beta=observed[1],se=observed[2],t=observed[3],
    sigma_GSE202379=observed[5],sigma_GSE244832=observed[6],
    previous_projection_beta=oldfit[1],previous_projection_n=nrow(old),
    testable=eligible,p=p,stability_median=reg$stability_median[i],
    effect_unit='fixed weighted mean log1p(CP10k) per source stage step',
    n_optim_fallback=nfb,n_failed_draws=nfail)
  result<-list(row=row,null=null);saveRDS(result,file.path(out,paste0(ct,'_',mod,'.rds')));result
}
results <- lapply(seq_len(117),one)
r<-rbindlist(lapply(results,`[[`,'row'));r[,q:=p.adjust(p,'BH',n=117)]
r[,selected:=testable & q<.05 & stability_median>=.5]
fwrite(r,file.path(out,'results.tsv'),sep='\t')
saveRDS(list(seed=20260923L,donor_order=donors,draws=draws,null=lapply(results,`[[`,'null')),file.path(out,'nulls.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines(sprintf('%d/117 selected; %d testable; min q %.6g. Fixed-weight projection; not native Hotspot graph/PCA. Optim fallback draws: %d; unresolved draws: %d.',
  sum(r$selected),sum(r$testable),min(r$q),sum(r$n_optim_fallback),sum(r$n_failed_draws)),file.path(out,'COMPLETED'))
