#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(susieR);library(coloc)})
root <- normalizePath('.')
run <- file.path(root,'Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z')
out <- file.path(run,paste0('mprip_refit_',Sys.getenv('SLURM_JOB_ID')))
stopifnot(!dir.exists(out));dir.create(out)
fm <- file.path(root,'GWAS/finemapping')
source(file.path(fm,'src/finemapping_functions.R'))
source(file.path(fm,'src/perf/get_ld_per_locus_cached.R'))
source(file.path(fm,'src/coloc_signal_eligibility.R'))
old <- readRDS(file.path(root,'Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/genetics/pilot_22811698/exports/MVP_ALT_EUR__ENSG00000133030__inputs_and_fits.rds'))
g <- as.data.table(old$aligned_summary)
e <- fread(file.path(run,'mprip_refit_inputs/eqtl_aligned.tsv'))
src <- fread(file.path(run,'mprip_sources/MPRIP_original_eqtl.tsv'))
original_ss <- data.frame(chromosome=src$CHR,position=src$POS,allele1=src$NEA,
                         allele2=src$EA,beta=src$Beta,standard_error=src$SE)
ld <- get_ld_per_locus(original_ss,'ENSG00000133030',17L,
                      min(src$POS)-500000L,max(src$POS)+500000L,'EUR')
stopifnot(!is.null(ld))
es <- as.data.table(ld[[1]])
stopifnot(setequal(es$SNP,e$SNP),!anyDuplicated(es$SNP))
e <- e[match(es$SNP,e$SNP)]
stopifnot(identical(es$SNP,e$SNP),
          isTRUE(all.equal(es$beta,-e$beta,tolerance=1e-12)),
          identical(paste0('17:',e$position),colnames(old$eqtl_fit$lbf_variable)),
          setequal(g[duplicated(position)|duplicated(position,fromLast=TRUE),position],17164840L),
          setequal(e[duplicated(position)|duplicated(position,fromLast=TRUE),position],17164840L))
read_matrix <- function(path,n) {
  con <- file(path,'rb');on.exit(close(con))
  x <- readBin(con,'double',n=n*n,size=8,endian='little')
  stopifnot(length(x)==n*n);matrix(x,n,n)
}
graw <- read_matrix(file.path(run,'mprip_sources/MPRIP_raw_single_window_LD_float64.bin'),nrow(g))
eraw_order <- fread(file.path(run,'mprip_refit_inputs/eqtl_aligned.tsv'))$SNP
eraw <- read_matrix(file.path(run,'mprip_refit_inputs/eqtl_raw_LD_float64.bin'),nrow(e))
ix <- match(e$SNP,eraw_order);eraw <- eraw[ix,ix]
gold <- old$R_regularized-diag(.001,nrow(g))
eold <- as.matrix(ld[[2]])
rm(ld);.LD_CACHE$entry <- NULL;gc(FALSE)
diagnostics <- list();results <- list();fit_cache <- new.env()
fit_one <- function(d,R,n,tag,allow_duplicate=FALSE) {
  stopifnot(nrow(R)==nrow(d),isTRUE(all.equal(unname(R),unname(t(R)),tolerance=1e-12)),
            all(is.finite(R)),all(is.finite(d$beta/d$standard_error)),all(d$standard_error>0))
  ids <- if(allow_duplicate) paste0('17:',d$position) else d$SNP
  if(!allow_duplicate) stopifnot(!anyDuplicated(ids),!anyDuplicated(d$position))
  ev <- eigen(R,symmetric=TRUE,only.values=TRUE)$values
  Rreg <- R+diag(.001,nrow(R));dimnames(Rreg) <- list(ids,ids)
  z <- setNames(d$beta/d$standard_error,ids)
  diagnostic <- data.table(arm=tag,n_variants=length(ids),n=n,
    minimum_eigenvalue=min(ev),minimum_regularized_eigenvalue=min(ev)+.001,
    eigenvalues_below_minus_1e_6=sum(ev < -1e-6),
    lambda_s=as.numeric(estimate_s_rss(z=z,R=Rreg,n=n)),
    status=if(min(ev)+.001 < -1e-6) 'non_PSD_regularized_LD' else 'valid_matrix')
  diagnostics[[length(diagnostics)+1L]] <<- diagnostic
  fwrite(rbindlist(diagnostics),file.path(out,'LD_diagnostics.tsv'),sep='\t')
  if(diagnostic$status!='valid_matrix' && !allow_duplicate) return(NULL)
  set.seed(42)
  fit <- susie_rss(z=z,R=Rreg,n=n,L=10,estimate_residual_variance=FALSE,
    residual_variance=1,check_R=FALSE,min_abs_corr=.5,max_iter=500)
  stopifnot(identical(colnames(fit$lbf_variable),ids))
  saveRDS(fit,file.path(out,paste0(tag,'_fit.rds')))
  fwrite(d,file.path(out,paste0(tag,'_aligned.tsv')),sep='\t')
  cat('Finished',tag,'converged',fit$converged,'CS',length(fit$sets$cs),'\n')
  fit
}
# Reproduce the historical eQTL input orientation and complete variant set.
echeck <- copy(e);echeck[,beta:=-beta]
replica <- fit_one(echeck,eold,1183L,'historical_eqtl_reconstruction',TRUE)
checks <- data.table(n_eqtl=nrow(e),n_gwas=nrow(g),
  source_eQTL_orientation_is_global_flip=TRUE,
  historical_eqtl_converged=replica$converged,
  historical_eqtl_max_alpha_difference=max(abs(replica$alpha-old$eqtl_fit$alpha)),
  raw_vs_block_gwas_max_difference=max(abs(graw-gold)),
  raw_vs_block_eqtl_max_difference=max(abs(eraw-eold)))
fwrite(checks,file.path(out,'reconstruction_checks.tsv'),sep='\t')
for(panel in c('rounded_blocks','raw_single_window')) {
  gr <- if(panel=='rounded_blocks') gold else graw
  er <- if(panel=='rounded_blocks') eold else eraw
  for(arm in c('literal_deletion','frequency_inferred_insertion','repeat_excluded')) {
    g_event <- if(arm=='literal_deletion') '17:17164840:CTAAA:C' else '17:17164840:C:CTAAA'
    gkeep <- g$position!=17164840L | (arm!='repeat_excluded' & g$SNP==g_event)
    ekeep <- e$position!=17164840L | (arm!='repeat_excluded' & e$SNP=='17:17164840:C:CTAAA')
    tag <- paste(panel,arm,sep='__')
    sg <- fit_one(g[gkeep],gr[gkeep,gkeep],old$n_eff,paste0(tag,'__gwas'))
    eq_tag <- paste(panel,if(arm=='repeat_excluded') 'repeat_excluded' else 'source_insertion','eqtl',sep='__')
    if(!exists(eq_tag,fit_cache,inherits=FALSE)) {
      assign(eq_tag,fit_one(e[ekeep],er[ekeep,ekeep],1183L,eq_tag),fit_cache)
    }
    se <- get(eq_tag,fit_cache)
    ans <- data.table(panel=panel,interpretation=arm,n_gwas=sum(gkeep),n_eqtl=sum(ekeep),
      gwas_converged=isTRUE(sg$converged),eqtl_converged=isTRUE(se$converged),
      gwas_cs=length(sg$sets$cs),eqtl_cs=length(se$sets$cs),
      pairs=0L,eligible_pairs=0L,max_H4=NA_real_,
      state=if(is.null(sg)||is.null(se)) 'untestable_non_PSD_LD' else 'untestable_nonconverged_or_no_CS')
    if(isTRUE(sg$converged)&&isTRUE(se$converged)&&length(sg$sets$cs)&&length(se$sets$cs)) {
      set.seed(42)
      cr <- coloc.susie(sg,se,p1=1e-4,p2=1e-4,p12=5e-6)
      pt <- coloc_susie_pair_table(sg,se,cr)
      stopifnot(pt$matches_coloc)
      fwrite(pt$pairs,file.path(out,paste0(tag,'__signal_pairs.tsv')),sep='\t')
      saveRDS(cr,file.path(out,paste0(tag,'__coloc.rds')))
      ans[,`:=`(pairs=nrow(pt$pairs),eligible_pairs=sum(pt$pairs$eligible),state=coloc_susie_method(pt$kept))]
      if(nrow(pt$kept)) ans[,max_H4:=max(pt$kept$PP.H4.abf)]
    }
    results[[length(results)+1L]] <- ans
    fwrite(rbindlist(results),file.path(out,'MPRIP_sensitivity_summary.tsv'),sep='\t')
    print(ans)
  }
}
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
stopifnot(nrow(rbindlist(results))==6L)
