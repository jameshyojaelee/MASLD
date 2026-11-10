#!/usr/bin/env Rscript
# Descriptive cross-assay donor associations, no edge-level significance claims.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3L,nzchar(Sys.getenv('SLURM_JOB_ID')))
a <- normalizePath(args[1L]); candidate <- args[2L]; output <- args[3L]
source(file.path(a,'scripts/compare_da.R')); data.table::setDTthreads(1L)
stopifnot(!dir.exists(output)); dir.create(output,recursive=TRUE)
cfg <- fromJSON(file.path(a,'config/analysis.json'),simplifyVector=FALSE)
pairing <- fread(cfg$cohort_metadata$GSE244832)
annotation_sources <- file.path(candidate,'annotation_sources.json')
link_root <- if(file.exists(annotation_sources)) fromJSON(annotation_sources)$disease_links else 'disease_links'
rna <- fread(file.path(candidate,'rna_models/donor_logCPM.tsv.gz'))
roster <- fread(file.path(candidate,'rna_pseudobulk/donors.tsv'))
decisions <- list(); tables <- list()
for (lineage_id in c('hepatocyte','stellate')) {
  x <- load_da_counts(cfg$counts_root,'GSE244832',lineage_id)
  d <- x$donors
  d[,batch_mm:=pairing$batch_mm[match(donor_id,pairing$donor_id)]]
  rr <- roster[roster$lineage==lineage_id & roster$eligible & roster$cross_assay_match_eligible & roster$condition %chin% c('NORMAL','MASH')]
  retained <- d$donor_id %in% rr$donor_id & pairing$confidence[match(d$donor_id,pairing$donor_id)]=='HIGH'
  d <- d[retained]
  condition_counts <- table(factor(d$condition,levels=c('NORMAL','MASH')))
  design <- if (all(condition_counts >= 4L)) tryCatch(make_da_design(d,'batch_mm')$design,error=function(e)e) else simpleError('fewer_than_four_verified_matches_per_condition')
  if (inherits(design,'error') || nrow(design)-qr(design)$rank < 4L) {
    reason <- if(inherits(design,'error')) conditionMessage(design) else 'fewer_than_four_residual_degrees_of_freedom'
    decisions[[length(decisions)+1L]] <- data.table(lineage=lineage_id,status='untestable',reason,n_donors=nrow(d))
    next
  }
  peaks <- fread(file.path(candidate,'disease_synthesis/results.tsv.gz'))
  coordinates <- peaks[cohort=='GSE244832' & lineage==lineage_id & adjustment=='condition_only',peak_coordinate]
  counts <- x$counts[,coordinates,drop=FALSE]
  y <- calcNormFactors(DGEList(counts=t(counts)),method='TMM')
  peak_values <- t(cpm(y,log=TRUE,prior.count=0.5))
  peak_values <- peak_values[match(d$donor_id,x$donors$donor_id),,drop=FALSE]
  rg <- rna[rna$lineage==lineage_id & rna$donor_id %in% d$donor_id]
  expression <- dcast(rg,gene_id~donor_id,value.var='logCPM')
  gene_ids <- expression$gene_id
  gene_values <- t(as.matrix(expression[,d$donor_id,with=FALSE]))
  stopifnot(!anyNA(gene_values))
  q <- qr(design)
  pr <- qr.resid(q,peak_values); gr <- qr.resid(q,gene_values)
  pn <- sqrt(colSums(pr^2)); gn <- sqrt(colSums(gr^2))
  links <- fread(file.path(candidate,link_root,paste0('GSE244832_',lineage_id),'links.tsv.gz'))
  links <- links[method=='cicero']
  links[, `:=`(gene_index=match(gene_id,gene_ids),peak_index=match(peak_coordinate,coordinates))]
  links[,status:=fifelse(is.na(gene_index),'RNA_gene_not_evaluable',fifelse(is.na(peak_index),'peak_not_in_joint_primary_universe','evaluable'))]
  links[,partial_correlation:=NA_real_]
  idx <- which(links$status=='evaluable')
  for (chunk in split(idx,ceiling(seq_along(idx)/5000))) {
    gi <- links$gene_index[chunk]; pi <- links$peak_index[chunk]
    denominator <- gn[gi]*pn[pi]
    value <- colSums(gr[,gi,drop=FALSE]*pr[,pi,drop=FALSE])/denominator
    value[denominator <= 1e-12] <- NA_real_
    links$partial_correlation[chunk] <- value
    links$status[chunk[!is.finite(value)]] <- 'no_residual_variation'
  }
  links[, `:=`(gene_index=NULL,peak_index=NULL,n_verified_donors=nrow(d),residual_df=nrow(d)-q$rank,
    adjustment='condition_plus_recorded_batch_mm',inference='descriptive_only_no_P_values_or_edge_calls',
    pairing='same_source_declared_donor_separate_RNA_ATAC_libraries_not_same_nucleus')]
  tables[[length(tables)+1L]] <- links
  decisions[[length(decisions)+1L]] <- data.table(lineage=lineage_id,status='descriptive_estimable',reason='',
    n_donors=nrow(d),residual_df=nrow(d)-q$rank,n_evaluable_links=sum(links$status=='evaluable'))
}
if(length(tables)) fwrite(rbindlist(tables,fill=TRUE),file.path(output,'donor_partial_correlations.tsv.gz'),sep='\t')
fwrite(rbindlist(decisions,fill=TRUE),file.path(output,'decisions.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(output,'sessionInfo.txt'))
write_json(list(status='completed',hypothesis_testing=FALSE,scientific_adoption=FALSE),file.path(output,'completion.json'),auto_unbox=TRUE,pretty=TRUE)
