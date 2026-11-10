#!/usr/bin/env Rscript
# Separate-library RNA corroboration, using donor pseudobulks only.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3L,nzchar(Sys.getenv('SLURM_JOB_ID')))
source(file.path(args[1L],'scripts/compare_da.R'))
data.table::setDTthreads(1L)
input <- args[2L]; output <- args[3L]
stopifnot(!dir.exists(output))
dir.create(output,recursive=TRUE)
receipt <- fromJSON(file.path(input,'completion.json'))
stopifnot(identical(receipt$status,'completed'),!isTRUE(receipt$technical_pilot))
con <- gzfile(file.path(input,'pseudobulk.mtx.gz'),'rb')
counts <- as(readMM(con),'CsparseMatrix'); close(con)
donors <- fread(file.path(input,'donors.tsv'))
genes <- fread(file.path(input,'genes.tsv'))
stopifnot(nrow(counts)==nrow(genes),ncol(counts)==nrow(donors),!anyDuplicated(genes$gene_id))
rownames(counts) <- genes$gene_id
results <- list(); designs <- list(); normalized <- list()
for (lineage_id in c('hepatocyte','stellate','macrophage')) {
  idx <- which(donors$lineage==lineage_id & donors$eligible & donors$condition %chin% c('NORMAL','MASH'))
  d <- copy(donors[idx]); x <- t(counts[,idx,drop=FALSE])
  eligible <- all(table(factor(d$condition,levels=c('NORMAL','MASH'))) >= 4L)
  if (!eligible) {
    designs[[length(designs)+1L]] <- data.table(lineage=lineage_id,adjustment='condition_only',status='untestable',
      reason='fewer_than_four_donors_per_condition_with_100_qualified_nuclei',n_donors=nrow(d))
    next
  }
  # RNA-seq native expression filtering; one fixed gene family for the matched models.
  design <- make_da_design(d)$design
  y <- DGEList(counts=t(x)); keep <- filterByExpr(y,design=design)
  x <- x[,keep,drop=FALSE]
  for (covariate in c('condition_only','batch_jb')) {
    cv <- if (covariate=='condition_only') character() else 'batch_jb'
    fit <- tryCatch(fit_da(x,d,'voom_qw',cv),error=function(e)e)
    if (inherits(fit,'error')) {
      designs[[length(designs)+1L]] <- data.table(lineage=lineage_id,adjustment=covariate,status='untestable',reason=conditionMessage(fit),n_donors=nrow(d))
      next
    }
    value <- copy(fit$results); setnames(value,'peak_coordinate','gene_id')
    value[, `:=`(lineage=lineage_id,cohort='GSE244832',adjustment=covariate,n_donors=nrow(d),
      effect_unit='log2_normalized_RNA_UMI_count',bh_family='all_RNA_filter_pass_genes_within_lineage_model',
      interpretation='separate_library_donor_level_RNA_corroboration_not_paired_nucleus_prediction')]
    results[[length(results)+1L]] <- merge(value,genes,by='gene_id',all.x=TRUE)
    designs[[length(designs)+1L]] <- cbind(data.table(lineage=lineage_id,adjustment=covariate,status='testable',reason='',n_donors=nrow(d)),
      fit$qc[1L,.(design_rank,residual_df,n_tested_peaks)])
  }
  y <- calcNormFactors(DGEList(counts=t(x)),method='TMM')
  logcpm <- cpm(y,log=TRUE,prior.count=0.5)
  dt <- as.data.table(logcpm,keep.rownames='gene_id')
  setnames(dt,names(dt)[-1L],d$donor_id)
  dt <- melt(dt,id.vars='gene_id',variable.name='donor_id',value.name='logCPM')
  dt[,lineage:=lineage_id]
  normalized[[length(normalized)+1L]] <- dt
}
fwrite(rbindlist(results,fill=TRUE),file.path(output,'rna_effects.tsv.gz'),sep='\t')
fwrite(rbindlist(designs,fill=TRUE),file.path(output,'model_decisions.tsv'),sep='\t')
fwrite(rbindlist(normalized,fill=TRUE),file.path(output,'donor_logCPM.tsv.gz'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(output,'sessionInfo.txt'))
write_json(list(status='completed',candidate_only=TRUE,same_nucleus_pairing=FALSE),
  file.path(output,'completion.json'),auto_unbox=TRUE,pretty=TRUE)
