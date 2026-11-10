#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
run <- 'Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z'
fm <- 'GWAS/finemapping/results'
p <- fread(file.path(run,'audit/confirmed_positive_identity_checks.tsv'))
p <- p[confirmed_duplicate_eqtl_fit==TRUE]
primary <- fread(file.path(run,'audit/primary_436_identity_checks.tsv'))$ensembl
pairs <- fread(file.path(fm,'susie_coloc_r2_20260924T053244Z/susie_coloc_signal_pairs.tsv.gz'))
pairs <- pairs[PP.H4>.5 & eligible==TRUE & ensembl %in% p$ensembl]
genes <- unique(p[,.(ensembl,chr,gene)])
rows <- lapply(seq_len(nrow(genes)),function(i) {
  gene_info <- genes[i]
  fit <- readRDS(file.path(fm,'eqtl_susie_polyfun',paste0('chr',gene_info$chr),paste0(gene_info$ensembl,'_susie.rds')))
  ids <- colnames(fit$lbf_variable)
  dup <- which(duplicated(ids)|duplicated(ids,fromLast=TRUE))
  pp <- pairs[ensembl==gene_info$ensembl]
  rbindlist(lapply(seq_len(nrow(pp)),function(j) {
    pair <- pp[j];k <- pair$idx2
    cs <- fit$sets$cs[[which(fit$sets$cs_index==k)]]
    data.table(gene=gene_info$gene,ensembl=gene_info$ensembl,chr=gene_info$chr,
      primary=gene_info$ensembl %in% primary,gwas=pair$gwas,idx1=pair$idx1,idx2=k,
      historical_H4=pair$PP.H4,duplicated_fit_positions=paste(unique(ids[dup]),collapse=';'),
      eqtl_CS_contains_duplicate_position=any(cs %in% dup),
      eqtl_alpha_mass_on_duplicate_positions=sum(fit$alpha[k,dup]))
  }))
})
d <- rbindlist(rows)
fwrite(d,file.path(run,'audit/affected_positive_signals.tsv'),sep='\t')
fwrite(d[,.(signal_pairs=.N,genes=uniqueN(ensembl),
  pairs_with_duplicate_CS_members=sum(eqtl_CS_contains_duplicate_position),
  genes_with_duplicate_CS_members=uniqueN(ensembl[eqtl_CS_contains_duplicate_position])),by=primary],
  file.path(run,'audit/affected_signal_summary.tsv'),sep='\t')
# No threshold on this mass grants validity: even a low-posterior duplicate can
# affect conditional regression. This describes where the defect is visible.
print(fread(file.path(run,'audit/affected_signal_summary.tsv')))
