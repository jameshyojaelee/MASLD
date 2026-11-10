#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(susieR);library(coloc)})
run <- 'Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z'
out <- file.path(run,paste0('gram_refit_',Sys.getenv('SLURM_JOB_ID')))
stopifnot(!dir.exists(out));dir.create(out)
source('GWAS/finemapping/src/coloc_signal_eligibility.R')
# Use the same checked fitting/PSD/convergence routine as the preceding six arms.
expressions <- parse(file.path(run,'code/refit_mprip.R'))
routine <- Filter(function(x) is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name('fit_one')),as.list(expressions))
stopifnot(length(routine)==1L);eval(routine[[1]])
diagnostics <- list();results <- list();fits <- new.env()
ids <- fread(file.path(run,'gram_inputs/genotype_LD_order.tsv'))$SNP
con <- file(file.path(run,'gram_inputs/genotype_LD_float64.bin'),'rb')
R <- matrix(readBin(con,'double',n=length(ids)^2,size=8,endian='little'),length(ids));close(con)
g <- fread(file.path(run,'mprip_sources/MPRIP_raw_LD_order.tsv'))
e <- fread(file.path(run,'mprip_refit_inputs/eqtl_aligned.tsv'))
g <- g[SNP %in% ids];e <- e[SNP %in% ids]
stopifnot(all(c('17:17164840:C:CTAAA','17:17164840:CTAAA:C') %in% ids))
for(arm in c('literal_deletion','frequency_inferred_insertion','repeat_excluded')) {
  gevent <- if(arm=='literal_deletion') '17:17164840:CTAAA:C' else '17:17164840:C:CTAAA'
  gd <- g[position!=17164840 | (arm!='repeat_excluded' & SNP==gevent)]
  ed <- e[position!=17164840 | (arm!='repeat_excluded' & SNP=='17:17164840:C:CTAAA')]
  ix <- match(gd$SNP,ids);sg <- fit_one(gd,R[ix,ix],406803L,paste0(arm,'__gwas'))
  eqtag <- if(arm=='repeat_excluded') 'repeat_excluded__eqtl' else 'source_insertion__eqtl'
  if(!exists(eqtag,fits,inherits=FALSE)) {
    ix <- match(ed$SNP,ids);assign(eqtag,fit_one(ed,R[ix,ix],1183L,eqtag),fits)
  }
  se <- get(eqtag,fits)
  ans <- data.table(interpretation=arm,reference_participants=388L,n_gwas=nrow(gd),n_eqtl=nrow(ed),
    gwas_converged=isTRUE(sg$converged),eqtl_converged=isTRUE(se$converged),
    gwas_cs=length(sg$sets$cs),eqtl_cs=length(se$sets$cs),pairs=0L,eligible_pairs=0L,max_H4=NA_real_,
    state='untestable_nonconverged_or_no_CS')
  if(isTRUE(sg$converged)&&isTRUE(se$converged)&&length(sg$sets$cs)&&length(se$sets$cs)) {
    set.seed(42);cr <- coloc.susie(sg,se,p1=1e-4,p2=1e-4,p12=5e-6)
    pt <- coloc_susie_pair_table(sg,se,cr);stopifnot(pt$matches_coloc)
    fwrite(pt$pairs,file.path(out,paste0(arm,'__signal_pairs.tsv')),sep='\t')
    saveRDS(cr,file.path(out,paste0(arm,'__coloc.rds')))
    ans[,`:=`(pairs=nrow(pt$pairs),eligible_pairs=sum(pt$pairs$eligible),state=coloc_susie_method(pt$kept))]
    if(nrow(pt$kept)) ans[,max_H4:=max(pt$kept$PP.H4.abf)]
  }
  results[[length(results)+1L]] <- ans
  fwrite(rbindlist(results),file.path(out,'MPRIP_genotype_LD_sensitivity.tsv'),sep='\t');print(ans)
}
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
