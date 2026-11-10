# Same donor-level composition model, applied to unique physical cells.
suppressPackageStartupMessages({library(data.table);library(limma)})
setDTthreads(1);set.seed(20261026L)
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
task <- file.path(root,"Analysis/Hypotheses/NUCLEUS-FIX-20261107")
out <- file.path(task,"composition");stopifnot(!dir.exists(out));dir.create(out)
lib <- fread(file.path(task,"unique_counts/metadata/library_metadata_source_F.tsv"))
lib <- lib[tolower(as.character(binary_eligible))=="true" & preparation_method%in%c("nuclei","unsorted")]
ct <- names(lib)[match("n_B_cells",names(lib)):match("n_pDCs",names(lib))]
stopifnot(length(ct)==16L)
counts <- lib[,lapply(.SD,sum),by=.(donor,dataset,condition_harmonized),.SDcols=ct]
counts[,total:=rowSums(.SD),.SDcols=ct];counts<-counts[total>=50]
fit_proportions <- function(counts, columns, label) {
  counts<-counts[rowSums(as.matrix(counts[,..columns]))>=50]
  mat<-as.matrix(counts[,..columns]);total<-rowSums(mat);p<-t(mat/total)
  design<-model.matrix(~dataset+condition,data=data.frame(dataset=factor(counts$dataset),
    condition=factor(counts$condition_harmonized,levels=c("Healthy","MASLD"))))
  stopifnot(qr(design)$rank==ncol(design),sum(counts$condition_harmonized=="Healthy")>=3,
    sum(counts$condition_harmonized=="MASLD")>=3)
  # The arcsine-square-root donor proportion model is native to the scRNA
  # propeller approach. Dataset is retained as an adjustment term.
  fit<-eBayes(lmFit(asin(sqrt(p)),design),robust=TRUE,trend=FALSE)
  k<-which(colnames(design)=="conditionMASLD")
  tt<-as.data.table(topTable(fit,coef=k,number=Inf,sort.by="none"),keep.rownames="feature")
  tt[,`:=`(feature=columns,donor_n=nrow(counts),n_Healthy=sum(counts$condition_harmonized=="Healthy"),
    n_MASLD=sum(counts$condition_harmonized=="MASLD"),
    mean_Healthy=rowMeans(p[,counts$condition_harmonized=="Healthy",drop=FALSE]),
    mean_MASLD=rowMeans(p[,counts$condition_harmonized=="MASLD",drop=FALSE]),
    effect_unit="arcsine-square-root proportion",multiple_testing=paste0("BH over ",length(columns)," ",label))]
  tt[,difference_percentage_points:=100*(mean_MASLD-mean_Healthy)]
  fwrite(tt,file.path(out,paste0(label,"_results.tsv")),sep="\t")
  fwrite(counts,file.path(out,paste0(label,"_donor_counts.tsv")),sep="\t")
}
fit_proportions(counts,ct,"lineages")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
writeLines("Captured-lineage proportions after unique-count and unanimous-label eligibility; not whole-tissue composition. Historical subtype clusters were not transferred.",file.path(out,"interpretation.txt"))
