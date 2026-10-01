#!/usr/bin/env Rscript
# KEY MESSAGE: protein abundance covariation is distinct from metabolic function.
suppressPackageStartupMessages({library(data.table); library(limma); library(ggplot2)})
args <- commandArgs(TRUE); stopifnot(length(args)==1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out,recursive=TRUE)
set.seed(20260929); setDTthreads(2); pdf.options(useDingbats=FALSE)
write <- function(x,name) fwrite(x,file.path(out,name),sep="\t",na="NA")
triad <- c("GNMT","MAT1A","CYP2C19")
background_sets <- c("HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_BILE_ACID_METABOLISM","HALLMARK_FATTY_ACID_METABOLISM")
writeLines(c("seed=20260929", "bootstrap=2000", "biological_unit=participant",
 "primary_family=3_pairs", "adjustment_family=12_pair_by_model_tests",
 "background_family=3_pairs_sensitivity_only", "common_complete_cases_across_models=true",
 "histology_covariates=recorded_numeric_grades", "selection_conditioned_example=true",
 paste0("background_sets=",paste(background_sets,collapse=";"))),file.path(out,"analysis_definition.txt"))
raw <- fread("data/PXD051911/liver_protein_quant.txt")
md <- fread("data/PXD051911/meta_data.txt")
samples <- setdiff(names(raw),c("ProteinAccessions","Genes","ProteinDescriptions"))
md <- unique(md[liver_proteomics_filename %in% samples],by="liver_proteomics_filename")
stopifnot(nrow(md)==58L)
md[, `:=`(sample_id=liver_proteomics_filename,
 group=ifelse(saf_diagnosis=="No_MASLD","Control","MASLD"),
 batch=factor(ifelse(grepl("^2019",liver_proteomics_filename),"2019","2020")),
 age=as.numeric(alder), bmi_num=as.numeric(bmi), sex=factor(gender),
 Fibrosis=as.numeric(sub("^F","",kleiner_fibrosis_grade)),
 Steatosis=as.numeric(steatosis_score), Ballooning=as.numeric(hepatocellular_ballooning_score),
 Inflammation=as.numeric(lobular_inflammation_score), NAS=as.numeric(nafld_activity_score))]
setorder(md, sample_id)
raw <- raw[!is.na(Genes)&Genes!=""&!grepl(";",Genes,fixed=TRUE)]
x <- as.matrix(raw[,md$sample_id,with=FALSE]); storage.mode(x)<-"double"
x[!is.finite(x)|x<=0]<-NA_real_;x<-log2(x);rownames(x)<-raw$Genes
x<-x[rowMeans(is.na(x))<=.5,,drop=FALSE]
x<-normalizeBetweenArrays(x,method="quantile")
indices<-split(seq_len(nrow(x)),rownames(x))
gene<-t(vapply(indices,function(i) if(length(i)==1) x[i,] else apply(x[i,,drop=FALSE],2,median,na.rm=TRUE),numeric(ncol(x))))
colnames(gene)<-md$sample_id;gene[!is.finite(gene)]<-NA_real_
stopifnot(all(triad %in% rownames(gene)))
for(g in triad)md[,(g):=gene[g,sample_id]]
# An outcome-blind, fixed biological background, not a selected list of DE proteins.
gmt<-strsplit(readLines("Analysis/downstream_analysis/pathway_analysis/data/genesets/hallmark.gmt"),"\t")
bg<-setdiff(unique(unlist(lapply(gmt,function(v) if(v[1]%in%background_sets)v[-c(1,2)] else character()))),triad)
bg<-intersect(bg,rownames(gene));stopifnot(length(bg)>=10)
write(data.table(gene=bg,background="three_fixed_Hallmark_metabolic_sets"),"metabolic_background.tsv")
zz<-t(scale(t(gene[bg,,drop=FALSE])))
md[, metabolic_background:=apply(zz,2,median,na.rm=TRUE)]
cols<-c(triad,"batch","age","bmi_num","sex","Fibrosis","NAS","Steatosis","Ballooning","Inflammation","metabolic_background")
d<-droplevels(md[group=="MASLD" & complete.cases(md[,..cols])])
stopifnot(!anyDuplicated(d$sample_id),nrow(d)>20)
write(d[,c("sample_id","group",cols),with=FALSE],"common_participants.tsv")
models<-list(baseline=c("batch","age","bmi_num","sex"),
 fibrosis=c("batch","age","bmi_num","sex","Fibrosis"),
 nas=c("batch","age","bmi_num","sex","NAS"),
 components=c("batch","age","bmi_num","sex","Fibrosis","Steatosis","Ballooning","Inflammation"),
 metabolic_background=c("batch","age","bmi_num","sex","metabolic_background"))
pairs<-combn(triad,2,simplify=FALSE)
calc<-function(z,pair,vars){
 xx<-model.matrix(reformulate(vars),data=z)
 if(qr(xx)$rank<ncol(xx)||nrow(xx)<=ncol(xx)+2)return(c(rho=NA,p=NA,df=NA))
 yy<-sapply(pair,function(g)rank(z[[g]],ties.method="average"))
 rr<-qr.resid(qr(xx),yy);r<-cor(rr[,1],rr[,2]);df<-nrow(xx)-ncol(xx)-1
 c(rho=r,p=2*pt(-abs(r*sqrt(df/(1-r*r))),df),df=df)
}
observed<-rbindlist(lapply(names(models),function(model)rbindlist(lapply(pairs,function(pair){
 v<-calc(d,pair,models[[model]])
 data.table(pair=paste(pair,collapse="–"),model=model,n=nrow(d),rho=v["rho"],p=v["p"],df=v["df"])
}))))
draws<-array(NA_real_,dim=c(2000,length(models),3))
for(b in seq_len(2000)){
 ii<-sample.int(nrow(d),nrow(d),replace=TRUE);z<-droplevels(d[ii])
 for(j in seq_along(models))for(k in seq_along(pairs))draws[b,j,k]<-tryCatch(calc(z,pairs[[k]],models[[j]])["rho"],error=function(e)NA_real_)
}
for(j in seq_along(models))for(k in seq_along(pairs)){
 v<-draws[,j,k];delta<-v-draws[,1,k];idx<-which(observed$model==names(models)[j]&observed$pair==paste(pairs[[k]],collapse="–"))
 observed[idx,`:=`(low=quantile(v,.025,na.rm=TRUE),high=quantile(v,.975,na.rm=TRUE),
 n_bootstrap=sum(is.finite(v)),delta_vs_baseline=rho-observed[model=="baseline"&pair==paste(pairs[[k]],collapse="–"),rho],
 delta_low=quantile(delta,.025,na.rm=TRUE),delta_high=quantile(delta,.975,na.rm=TRUE))]
}
observed[,q_three_pairs:=p.adjust(p,"BH",n=3),by=model]
observed[model!="metabolic_background",q_twelve:=p.adjust(p,"BH",n=12)]
observed[,`:=`(seed=20260929L,interval="pointwise_participant_percentile_bootstrap",selection_conditioned=TRUE)]
write(observed,"protein_pair_covariance.tsv")
write(data.table(metric=c("all_liver_participants","common_MASLD_participants","background_genes"),value=c(nrow(md),nrow(d),length(bg))),"denominators.tsv")
# Reconstruct the historical baseline as a check; mismatches are reported, not hidden.
expected<-data.table(pair=c("GNMT–MAT1A","GNMT–CYP2C19","MAT1A–CYP2C19"),reported_rho=c(.472,.406,.472))
check<-merge(observed[model=="baseline",.(pair,n,rho)],expected,by="pair")
check[,difference_from_rounded_report:=rho-reported_rho];write(check,"historical_reproduction.tsv")
p<-ggplot(observed[model!="metabolic_background"],aes(rho,pair,color=model,shape=model))+
 geom_vline(xintercept=0,color="#9E9E9E",linewidth=.3)+
 geom_errorbar(aes(xmin=low,xmax=high),orientation="y",width=.15,position=position_dodge(width=.6),linewidth=.3)+
 geom_point(position=position_dodge(width=.6),size=1)+
 labs(x="Adjusted protein rank correlation (95% bootstrap interval)",y=NULL,color=NULL,shape=NULL)+
 theme_classic(base_size=6)+theme(text=element_text(size=6),legend.position="bottom")
ggsave(file.path(out,"S3_metabolic_protein_covariance.pdf"),p,width=4.6,height=2.1,device=pdf,useDingbats=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
writeLines("Candidate: common-participant covariance, no functional or mediation inference.",file.path(out,"COMPLETE"))
