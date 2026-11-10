#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma);library(jsonlite)})
ROOT <- '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
OUT <- file.path(ROOT,'Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/transportability')
setDTthreads(1); options(digits=17); set.seed(20261107)
M <- fread(file.path(OUT,'input_manifest.tsv'))
sha <- function(p) sub(' .*','',system2('sha256sum',shQuote(p),stdout=TRUE))
stopifnot(all(M$sha256==vapply(M$path,sha,'')))
P <- setNames(M$path,M$input)
# SHA entry is relative to OUT; compare explicitly without changing process cwd.
spec_sha <- strsplit(readLines(file.path(OUT,'ANALYSIS_SPEC.sha256')),'  ',fixed=TRUE)[[1]][1]
stopifnot(identical(sha(file.path(OUT,'ANALYSIS_SPEC.md')),spec_sha))
if(file.exists(file.path(OUT,'SE_COVERAGE_PASS.json'))) stop('Already prepared; refusing overwrite')
pooled <- fread(P[['pooled']]); stopifnot(nrow(pooled)==24115L,!anyDuplicated(pooled$gene))
genes <- pooled$gene
C <- fread(P[['cohort']]); N <- fread(P[['native126']]); meta<-as.data.table(readRDS(P[['metadata']]))
roster<-fread(P[['roster']]); dge<-readRDS(P[['dge']])
stopifnot(nrow(C)==5L*24115L,nrow(N)==15592L,nrow(meta)==774L,nrow(roster)==774L,
 !anyDuplicated(meta$analysis_unit_id),!anyDuplicated(roster$analysis_unit_id),
 setequal(meta$sample_id,roster$sample_id),setequal(meta$sample_id,colnames(dge)),
 identical(rownames(dge),genes),!anyDuplicated(C[,.(gene_id_versioned,dataset)]),
 setequal(N$gene_id_versioned,genes[genes %in% N$gene_id_versioned]))
R <- merge(meta,roster,by='sample_id',suffixes=c('','_roster'))
stopifnot(all(R$analysis_unit_id==R$analysis_unit_id_roster))
fwrite(R,file.path(OUT,'immutable_unit_metadata.tsv'),sep='\t')
datasets<-sort(unique(C$dataset));stopifnot(length(datasets)==5L)
L <- copy(C); L[,native_testable:=TRUE]
idx<-which(L$dataset=='GSE126848'); j<-match(L$gene_id_versioned[idx],N$gene_id_versioned)
for(v in c('logFC','t','P.Value','FDR')) set(L,i=idx,j=v,value=N[[v]][j])
L[idx,native_testable:=!is.na(j)]
L[,moderated_SE:=logFC/t]
stopifnot(all(is.finite(L[native_testable==TRUE,moderated_SE])),all(L[native_testable==TRUE,moderated_SE]>0))
L[,native_family_size:=ifelse(dataset=='GSE126848',15592L,24115L)]
L[,canonical_cohort_call:=native_testable & FDR<.05 & abs(logFC)>.5]
L[is.na(canonical_cohort_call),canonical_cohort_call:=FALSE]
checks<-list()
for(ds in datasets){
 s<-meta[dataset==ds];setorder(s,sample_id)
 varying<-uniqueN(s$inferred_sex)>1L
 inf<-data.frame(inferred_sex=factor(s$inferred_sex),group=factor(s$group_binary,levels=c('Control','Disease')))
 X<-model.matrix(if(varying)~inferred_sex+group else ~group,inf)
 d<-dge[,match(s$sample_id,colnames(dge))]
 if(ds=='GSE126848'){
  keep<-filterByExpr(d,group=inf$group);d<-d[keep,,keep.lib.sizes=FALSE];d<-calcNormFactors(d,method='RLE')
  stopifnot(setequal(rownames(d),N$gene_id_versioned))
 }
 stopifnot(qr(X)$rank==ncol(X))
 v<-voomWithQualityWeights(d,X,plot=FALSE);f<-eBayes(lmFit(v,X));ix<-match('groupDisease',colnames(X))
 b<-f$coefficients[,ix];se<-f$stdev.unscaled[,ix]*sqrt(f$s2.post)
 ref<-L[dataset==ds & native_testable==TRUE][match(rownames(d),gene_id_versioned)]
 z<-data.table(dataset=ds,n_genes=nrow(d),n_units=nrow(s),n_control=sum(inf$group=='Control'),n_disease=sum(inf$group=='Disease'),
  max_abs_effect_difference=max(abs(b-ref$logFC)),max_abs_SE_difference=max(abs(se-ref$moderated_SE)),
  max_abs_t_difference=max(abs(f$t[,ix]-ref$t)))
 checks[[ds]]<-z; print(z)
 stopifnot(z$max_abs_effect_difference<1e-8,z$max_abs_SE_difference<1e-8,z$max_abs_t_difference<1e-7,
   all(ref$n_control==z$n_control),all(ref$n_disease==z$n_disease))
}
V<-rbindlist(checks);fwrite(V,file.path(OUT,'exact_refit_moderated_SE_checks.tsv'),sep='\t')
L[,unit_definition:='conservative_genetic_identity_cluster_not_verified_person']
L[,SE_definition:='eBayes_moderated_SE_logFC_divided_by_t_verified_against_stdev_unscaled_sqrt_s2_post']
fwrite(L,file.path(OUT,'cohort_effects_long.tsv.gz'),sep='\t')
S<-meta[,.(n_units=.N,n_control=sum(group_binary=='Control'),n_disease=sum(group_binary=='Disease')),by=dataset]
fwrite(S,file.path(OUT,'cohort_unit_counts.tsv'),sep='\t')
write_json(list(pass=TRUE,n_genes=length(genes),n_gene_cohort_slots=nrow(L),n_eligible_slots=sum(L$native_testable),
 n_untestable_GSE126848=sum(!L$native_testable),n_units=nrow(meta),n_fixed_calls=sum(pooled$padj<.05 & abs(pooled$logFC)>.5),
 SE_type='moderated',input_manifest_sha256=sha(file.path(OUT,'input_manifest.tsv'))),
 file.path(OUT,'SE_COVERAGE_PASS.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(OUT,'sessionInfo_preflight.txt'))
