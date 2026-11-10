# Fixed memberships and L1 weights on unique integer counts; no discovery.
suppressPackageStartupMessages({library(data.table);library(nlme);library(parallel)})
root <- '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
task <- file.path(root,'Analysis/Hypotheses/NUCLEUS-FIX-20261107')
source(file.path(root,'Analysis/Hypotheses/H59-20261015T112432Z/scripts/h59_lib.R'))
out <- file.path(task,'fixed_programs'); dir.create(out,showWarnings=FALSE)
stopifnot(!file.exists(file.path(out,'results.tsv')))
meta <- fread(file.path(task,'unique_counts/metadata/donor_metadata.tsv'))
reg <- fread(inp('registry'))[,.(cell_type,module=as.integer(module),stability_median)]
reg <- unique(reg);setorder(reg,cell_type,module);stopifnot(nrow(reg)==117)
lineage <- c(hepatocytes='Hepatocytes',fibroblasts='Fibroblasts',macrophages='Macrophages',cholangiocytes='Cholangiocytes',tcells='T cells')
prepare <- function(path) {
  x <- fread(path); x <- x[lineage == unname(lineage_map[cell_type])]
  x <- x[(cell_type=='tcells' | n_cells>=30)]
  x <- x[cell_type=='hepatocytes' | !dataset%in%c('GSE136103','Liver_Atlas')]
  x <- x[,.(score=sum(score_sum)/sum(n_cells),n_cells=sum(n_cells)),by=.(donor,dataset,cell_type,module)]
  x <- merge(x,meta[,.(donor,disease_stage_coarse,exclude_stage_analysis)],by='donor')
  x[,stage:=match(disease_stage_coarse,STAGES)-1L]
  x[!exclude_stage_analysis & !is.na(stage) & dataset%in%FIT_DATASETS]
}
lineage_map <- lineage
current <- prepare(file.path(task,'unique_counts/fixed_weight_projection_by_capture.tsv.gz'))
previous <- prepare(file.path(task,'unique_counts/historical_fixed_weight_projection_by_run.tsv.gz'))
fwrite(current,file.path(out,'donor_scores.tsv'),sep='\t')
fwrite(previous,file.path(out,'historical_projection_donor_scores.tsv'),sep='\t')
donors <- sort(unique(current$donor)); draws <- 19999L
set.seed(20260923L);keys <- matrix(runif(draws*length(donors)),nrow=draws)
one <- function(i) {
  ct <- reg$cell_type[i];mod <- reg$module[i]
  cache <- file.path(out,paste0(ct,'_',mod,'.rds'))
  if(file.exists(cache)) return(readRDS(cache))
  d <- prog_data(current,ct,mod,FIT_DATASETS)
  old <- prog_data(previous,ct,mod,FIT_DATASETS)
  eligible <- all(vapply(FIT_DATASETS,function(ds) {
    tab<-table(d[dataset==ds,stage]);sum(tab)>=8 && sum(tab>=2)>=2
  },TRUE))
  observed <- if(eligible) fit_gls(d$score,d$stage,d$dataset) else rep(NA_real_,6)
  oldfit <- if(nrow(old)>5 && uniqueN(old$dataset)==2) fit_gls(old$score,old$stage,old$dataset) else rep(NA_real_,6)
  null <- rep(NA_real_,draws); p <- 1
  if(eligible && is.finite(observed[3])) {
    permuted <- perm_stage(d,keys,donors)
    null <- vapply(seq_len(draws),function(b) fit_gls(d$score,permuted[b,],d$dataset)[['t']],0.)
    stopifnot(all(is.finite(null)))
    p <- (1+sum(abs(null)>=abs(observed[3])*(1-1e-12)))/(draws+1)
  } else eligible <- FALSE
  row <- data.table(cell_type=ct,module=mod,n=nrow(d),beta=observed[1],se=observed[2],t=observed[3],
    sigma_GSE202379=observed[5],sigma_GSE244832=observed[6],
    previous_projection_beta=oldfit[1],previous_projection_n=nrow(old),
    testable=eligible,p=p,stability_median=reg$stability_median[i],
    effect_unit='fixed weighted mean log1p(CP10k) per source stage step')
  result<-list(row=row,null=null);saveRDS(result,cache);result
}
results <- mclapply(seq_len(117),one,mc.cores=as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','16')),mc.preschedule=FALSE)
stopifnot(!any(vapply(results,inherits,TRUE,'try-error')))
r<-rbindlist(lapply(results,`[[`,'row'));r[,q:=p.adjust(p,'BH',n=117)]
r[,selected:=testable & q<.05 & stability_median>=.5]
fwrite(r,file.path(out,'results.tsv'),sep='\t')
saveRDS(list(seed=20260923L,donor_order=donors,draws=draws,null=lapply(results,`[[`,'null')),file.path(out,'nulls.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines(sprintf('%d/117 selected; %d testable. Fixed-weight projection; not native Hotspot graph/PCA.',sum(r$selected),sum(r$testable)),file.path(out,'COMPLETED'))
