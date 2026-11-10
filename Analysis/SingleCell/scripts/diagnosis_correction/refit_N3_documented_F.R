# Numeric fibrosis sensitivity only; retain completed primary/permutation results.
suppressPackageStartupMessages({library(data.table);library(limma);library(edgeR);library(jsonlite)})
setDTthreads(1);set.seed(20260923L)
root<-Sys.getenv("MASLD_PROJECT_ROOT")
task<-file.path(root,"Analysis/SingleCell/candidates/source-diagnosis-correction-20261026")
# Use the native RNA-count fit functions without evaluating the full producer.
expr<-parse(file.path(task,"refit_scripts/n3_model.R"))
for (e in expr) if (is.call(e)&&identical(e[[1]],as.name("<-"))&&is.symbol(e[[2]])&&as.character(e[[2]])%in%c("fit_voom","coef_table","check_rank")) eval(e)
for (tag in c("raw","corrected")) {
  base<-file.path(task,"N3",tag)
  native<-file.path(base,if(tag=="raw") "results" else "results_corrected")
  out<-file.path(base,"results_source_F");stopifnot(!dir.exists(out));dir.create(out)
  excluded<-c("S2_GSE202379_documented_F_targets.tsv","target_classes.tsv","summary.json","donors_in_model.tsv")
  for (p in list.files(native,full.names=TRUE)) if (!basename(p)%in%excluded) stopifnot(file.symlink(p,file.path(out,basename(p))))
  dm<-fread(file.path(base,"donor_meta_source_F.tsv"))
  inc<-dm[included==TRUE & dataset=="GSE202379" & !is.na(F_documented)]
  inc[,sex:=factor(sex_call,levels=c("Female","Male"))]
  design<-model.matrix(~sex+F_documented,inc);check_rank(design,"numeric SAF sensitivity")
  cnt<-fread(file.path(base,"hep_donor_counts.tsv.gz"));tl<-fread(file.path(base,"target_list.tsv"))
  mat<-as.matrix(cnt[,inc$donor,with=FALSE]);rownames(mat)<-cnt$symbol;storage.mode(mat)<-"double"
  fit<-fit_voom(mat,design)
  tab<-coef_table(fit$fit,"F_documented")[symbol%in%tl$symbol]
  tab[,`:=`(q_bh=p.adjust(P.Value,"BH"),n_donors=nrow(inc))]
  fwrite(tab,file.path(out,"S2_GSE202379_documented_F_targets.tsv"),sep="\t")
  cls<-fread(file.path(native,"target_classes.tsv"))
  drop<-intersect(c("S2_logFC_per_F","S2_q_bh"),names(cls));cls[,(drop):=NULL]
  cls<-merge(cls,tab[,.(symbol,S2_logFC_per_F=logFC,S2_q_bh=q_bh)],by="symbol",all.x=TRUE)
  fwrite(cls,file.path(out,"target_classes.tsv"),sep="\t")
  summary<-read_json(file.path(native,"summary.json"),simplifyVector=TRUE)
  summary$S2_n_donors<-nrow(inc)
  summary$S2_source<-"Numeric deposited SAF only; free-text healthy control and end stage do not establish measured F0/F4"
  write_json(summary,file.path(out,"summary.json"),auto_unbox=TRUE,pretty=TRUE,digits=10,na="null")
  roster<-fread(file.path(native,"donors_in_model.tsv"))
  roster[,F_documented:=dm$F_documented[match(donor,dm$donor)]]
  fwrite(roster,file.path(out,"donors_in_model.tsv"),sep="\t")
  fwrite(inc,file.path(out,"numeric_SAF_donors.tsv"),sep="\t")
  writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo_documented_F.txt"))
}
writeLines("Numeric SAF sensitivity refit; primary fits and permutation results unchanged.",file.path(task,"N3_DOCUMENTED_F_COMPLETED"))
