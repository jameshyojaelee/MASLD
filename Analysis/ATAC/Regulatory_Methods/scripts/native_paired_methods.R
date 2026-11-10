#!/usr/bin/env Rscript
# Native TRIPOD fit on one training-only stage; Python prepares donor splits.
# Native p-values are not donor-calibrated. Never use them as disease evidence.
args <- commandArgs(TRUE)
stopifnot(nzchar(Sys.getenv("SLURM_JOB_ID")))
suppressPackageStartupMessages(library(TRIPOD))
run_stage <- function(a) {
stopifnot(length(a)==2L,!dir.exists(a[2]))
stopifnot(any(grepl("TRIPOD_full_axis_CPM1e6_no_log",readLines(file.path(a[1],"stage.json")))))
read <- function(f) read.delim(file.path(a[1],f),check.names=FALSE)
meta <- read("metacells.tsv"); rna <- as.matrix(read("rna.tsv"))
atac <- as.matrix(read("atac.tsv")); motifs <- read("motifs.tsv")
target <- read("target.tsv")$gene_id
stopifnot(length(target)==1L, length(unique(meta$donor_id))>=3L,
          all(meta$split=="training"), nrow(meta)<=1000L,
          nrow(rna)==nrow(meta), nrow(atac)==nrow(meta), ncol(atac)<=25L,
          all(is.finite(rna)), all(is.finite(atac)), all(rna>=0),all(atac>=0))
tfs <- sort(intersect(setdiff(unique(motifs$tf_gene_id),target),colnames(rna)))
tfs <- tfs[colSums(rna[,tfs,drop=FALSE])>0]
stopifnot(length(tfs)>0L,length(tfs)<=10L)
pm <- matrix(0,ncol(atac),length(tfs),dimnames=list(colnames(atac),tfs))
for(i in seq_len(nrow(motifs))) {
  if(motifs$peak_id[i] %in% rownames(pm) && motifs$tf_gene_id[i] %in% tfs)
    pm[motifs$peak_id[i],motifs$tf_gene_id[i]] <- 1
}
stopifnot(sum(pm)>0L,sum(pm)<=250L)
x <- list(gene.name=target,Yg=rna[,target],Xt=atac,Yj=rna[,tfs,drop=FALSE],
          peakxmotif.g=pm,TF.g=tfs,nonzero.peakxmotif.g=which(pm!=0,arr.ind=TRUE))
set.seed(42)
fit <- TRIPOD::fitModel(x,model.name="TRIPOD",match.by="Xt",log=NULL,verbose=FALSE)
ij <- which(pm!=0,arr.ind=TRUE)
trios <- data.frame(gene_id=target,peak_id=rownames(pm)[ij[,1]],
  tf_gene_id=colnames(pm)[ij[,2]],coefficient_level1=fit$coef.level1[ij],
  native_nominal_p_level1=fit$pval.level1[ij],coefficient_level2=fit$coef.level2[ij],
  native_nominal_p_level2=fit$pval.level2[ij],
  evidence_state="training_only_motif_compatible_association")
dir.create(a[2],recursive=TRUE)
write.table(trios,file.path(a[2],"trios.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
valid <- trios[is.finite(trios$coefficient_level1),]
links <- data.frame(gene_id=character(),peak_id=character(),score=double())
if(nrow(valid)) {
  links <- aggregate(abs(valid$coefficient_level1),valid[c("gene_id","peak_id")],max)
  names(links)[3] <- "score"
  links <- links[order(-links$score,links$peak_id),]
}
links$rank <- seq_len(nrow(links))
links$method <- rep("TRIPOD_Xt_L1",nrow(links))
links$score_unit <- rep("max_absolute_native_level1_coefficient_across_TFs",nrow(links))
links$evidence_state <- rep("training_only_motif_compatible_association",nrow(links))
write.table(links,file.path(a[2],"links.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file=file.path(a[2],"sessionInfo.txt"))
writeLines('{"status":"native_training_fit_completed_not_evaluated"}',file.path(a[2],"native_fit.json"))
}
if(length(args)==2L) {
  run_stage(args)
} else if(length(args)==1L) {
  jobs <- read.delim(args[1],stringsAsFactors=FALSE)
  stopifnot(all(c("stage","output") %in% names(jobs)),nrow(jobs)<=200L)
  for(i in seq_len(nrow(jobs))) {
    tryCatch(run_stage(c(jobs$stage[i],jobs$output[i])),error=function(e) {
      message("STAGE FAILED: ",jobs$stage[i],": ",conditionMessage(e))
    })
  }
} else stop("Require STAGE NEWOUTPUT or batch.tsv")
