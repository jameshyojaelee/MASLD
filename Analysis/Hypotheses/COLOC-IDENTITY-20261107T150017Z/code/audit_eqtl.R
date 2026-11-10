#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
run <- normalizePath('Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z')
root <- normalizePath('.')
out <- file.path(run,'audit')
coll <- fread(file.path(out,'LD_ambiguous_positions.tsv'))
lookup <- split(coll$pos,paste(coll$panel,coll$chr,sep=':'))
paths <- unlist(lapply(1:22,function(ch) list.files(file.path(root,'GWAS/finemapping/results/eqtl_susie_polyfun',paste0('chr',ch)),pattern='_susie.rds$',full.names=TRUE)))
inspect <- function(path) {
  fit <- readRDS(path); ids <- colnames(fit$lbf_variable)
  ch <- as.integer(sub('chr','',basename(dirname(path))))
  gene <- sub('_susie.rds$','',basename(path)); pos <- as.integer(sub('^[0-9]+:','',ids))
  ans <- data.table(ensembl=gene,chr=ch,n_variants=length(ids),n_duplicate_ids=sum(duplicated(ids)),
    converged=isTRUE(fit$converged),min_position=min(pos),max_position=max(pos))
  for(panel in unique(coll$panel)) ans[[paste0('reference_collision_positions_',panel)]] <- sum(unique(pos) %in% lookup[[paste(panel,ch,sep=':')]])
  ans
}
rows <- parallel::mclapply(paths,inspect,mc.cores=4L)
stopifnot(!any(vapply(rows,inherits,logical(1),'try-error')))
d <- rbindlist(rows);stopifnot(!anyDuplicated(d$ensembl))
fwrite(d,file.path(out,'all_saved_eqtl_fit_identity_checks.tsv'),sep='\t')
master <- fread(file.path(out,'all_50_study_records.tsv.gz'))
joined <- merge(master,d,by=c('ensembl','chr'),all.x=TRUE,sort=FALSE)
joined[,actual_panel:=fifelse(ancestry=='EUR','polyfun_eur',paste0('1kg_',tolower(ancestry),'_gram'))]
joined[,reference_collision_positions:=NA_integer_]
for(panel in unique(joined$actual_panel)) {
 col <- paste0('reference_collision_positions_',panel)
 if(col %in% names(joined)) joined[actual_panel==panel,reference_collision_positions:=get(col)]
}
fwrite(joined,file.path(out,'all_50_study_identity_risk.tsv.gz'),sep='\t')
fwrite(joined[PP.H4.susie>.5],file.path(out,'positive_gene_study_identity_risk.tsv'),sep='\t')
fwrite(joined[,.(gene_study_records=.N,positive=sum(PP.H4.susie>.5,na.rm=TRUE),
  missing_eqtl_fit=sum(is.na(n_variants)),duplicate_eqtl_fit=sum(n_duplicate_ids>0,na.rm=TRUE),
  positive_reference_exposed=sum(PP.H4.susie>.5 & reference_collision_positions>0,na.rm=TRUE)),by=.(gwas_name,ancestry,actual_panel)],file.path(out,'study_risk_summary.tsv'),sep='\t')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo_eqtl_audit.txt'))
print(d[,.(fits=.N,duplicate_fits=sum(n_duplicate_ids>0),variants=sum(n_variants))])
