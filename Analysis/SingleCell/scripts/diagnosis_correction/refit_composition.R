# Donor-level transformed-proportion models for scRNA composition; fixed labels.
suppressPackageStartupMessages({library(data.table);library(limma)})
setDTthreads(1);set.seed(20261026L)
root <- Sys.getenv("MASLD_PROJECT_ROOT")
task <- file.path(root,"Analysis/SingleCell/candidates/source-diagnosis-correction-20261026")
out <- file.path(task,"composition")
stopifnot(!dir.exists(out));dir.create(out)
source(file.path(root,"Analysis/SingleCell/scripts/lib_source_diagnosis.R"))
lib <- load_source_diagnosis(root,"library")
# Whole-liver composition cannot be estimated from CD45-enriched/depleted fractions.
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
hep<-fread(file.path(root,"Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_subtype_metadata.csv"),select=c("sample","hepatocyte_subtype"))
j<-match(hep$sample,lib$sample)
hep[,`:=`(donor=lib$donor[j],dataset=lib$dataset[j],condition_harmonized=lib$condition_harmonized[j])]
hep<-hep[!is.na(donor)]
mapping<-fread(file.path(root,"Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/meta_subtype_mapping.csv"))
hep[,meta_subtype:=mapping$meta_subtype[match(as.character(hepatocyte_subtype),as.character(mapping$subtype))]]
hep[is.na(meta_subtype),meta_subtype:="Other"]
for (field in c("hepatocyte_subtype","meta_subtype")) {
  tab<-hep[,.N,by=c("donor","dataset","condition_harmonized",field)]
  wide<-dcast(tab,as.formula(paste("donor+dataset+condition_harmonized~",field)),value.var="N",fill=0)
  cols<-setdiff(names(wide),c("donor","dataset","condition_harmonized"))
  fit_proportions(wide,cols,field)
}
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
writeLines("Donor proportions; dataset-adjusted arcsine-square-root linear models with empirical Bayes moderation. BH separately over 16 lineages, frozen raw hepatocyte clusters, and frozen meta-clusters. Association is not lineage specificity or longitudinal progression.",file.path(out,"METHOD.txt"))
