# KEY MESSAGE: Corrected donor labels change the estimated stage effects and their uncertainty.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
setDTthreads(1);grDevices::pdf.options(useDingbats=FALSE)
root<-Sys.getenv("MASLD_PROJECT_ROOT")
task<-file.path(root,"Analysis/SingleCell/candidates/source-diagnosis-correction-20261026")
out<-file.path(task,"figures_v3")
stopifnot(!dir.exists(out));dir.create(out)
dir.create(file.path(out,"source_tables"));dir.create(file.path(out,"producers"))
r<-fread(file.path(task,"programs/results_primary_family.tsv"))
old<-file.path(root,"figures/candidates/review-r3-fixes-20261015T112148Z/F4abcd")
template<-fread(file.path(old,"source_tables/fig4c_h59_registry_117.tsv"))[,.(cell_type,module,name,rk)]
w<-merge(template,r,by=c("cell_type","module"))
stopifnot(nrow(w)==117L,!anyNA(w$gls_beta))
w[,`:=`(beta=gls_beta,SE=gls_se,lo=gls_ci_low,hi=gls_ci_high,q=perm_q)]
scores<-fread(file.path(task,"programs/donor_program_scores_primary.tsv"))
scores<-scores[dataset%in%c("GSE202379","GSE244832") & is.finite(score)]
scores[,centered:=score-mean(score),by=.(cell_type,module)]
means<-dcast(scores[,.(centered=mean(centered)),by=.(cell_type,module,disease_stage_coarse)],
  cell_type+module~disease_stage_coarse,value.var="centered")
w<-merge(w,means,by=c("cell_type","module"))
fwrite(w,file.path(out,"source_tables/fig4c_h59_registry_117.tsv"),sep="\t")
critical<-r[,.(program,null_abs_t_crit=null_critical_abs_t)]
critical_path<-file.path(task,"programs/null_critical_abs_t.tsv")
if(file.exists(critical_path)) stopifnot(isTRUE(all.equal(fread(critical_path),critical))) else fwrite(critical,critical_path,sep="\t")
# Keep the displayed program identities fixed; never rerank the main-panel
# membership after seeing corrected effects.
producer<-file.path(old,"repair/producers/02r_render_fig4c_h59.R")
s<-paste(readLines(producer),collapse="\n")
replace_once<-function(a,b) {
  stopifnot(length(gregexpr(a,s,fixed=TRUE)[[1]])==1L,grepl(a,s,fixed=TRUE))
  s<<-sub(a,b,s,fixed=TRUE)
}
replace_once('FIG  <- file.path(BASE, "figures/candidates/review-r3-fixes-20261015T112148Z/F4abcd")',paste0('FIG  <- "',out,'"'))
replace_once('OUT  <- file.path(FIG, "repair")','OUT <- FIG')
replace_once('Analysis/Hypotheses/H59-20261015T112432Z/out/null_critical_abs_t.tsv',
             'Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/programs/null_critical_abs_t.tsv')
replace_once('sum(W$selected_H59) == 0L','all(is.finite(W$q))')
s<-gsub('_candidate_repair.pdf','_source_qualified.pdf',s,fixed=TRUE)
s<-gsub('"117-program"','"117-test interval"',s,fixed=TRUE)
script<-file.path(out,"producers/render_programs.R");writeLines(s,script)
source(script)
ccc<-fread(file.path(task,"CCC/fig4d_program_ccc_associations.tsv"))
ccc[,label:=headline_label]
ccc[,label:=factor(label,levels=rev(unique(label)))]
ccc[,`:=`(low=beta-qt(.975,residual_df)*se,high=beta+qt(.975,residual_df)*se)]
fwrite(ccc,file.path(out,"source_tables/ccc_26.tsv"),sep="\t")
p<-ggplot(ccc[testable==TRUE],aes(beta,label))+geom_vline(xintercept=0,color="#9E9E9E",linewidth=.25)+
  geom_errorbar(aes(xmin=low,xmax=high),orientation="y",width=.2,linewidth=.3)+
  geom_point(aes(fill=evidence_state),shape=21,size=1.2,stroke=.25)+
  scale_fill_manual(values=c(supported="#C9265E",indeterminate="#9E9E9E"))+
  facet_wrap(~program_label,ncol=1)+labs(x="Adjusted program association (95% model CI)",y=NULL,fill=NULL,
    caption="20 estimable associations of 26 planned; BH correction across 26 tests")+
  theme_classic(base_size=6)+theme(text=element_text(size=6,face="plain"),legend.position="bottom")
ggsave(file.path(out,"ccc_26_source_qualified.pdf"),p,width=5.5,height=6,device=cairo_pdf)
context<-fread(file.path(task,"OLS_ambient/hotspot/context_effects.tsv"))[cell_type=="hepatocytes" & module%in%c(8L,20L)]
pooled<-fread(file.path(task,"OLS_ambient/hotspot/primary_effects.tsv"))[cell_type=="hepatocytes" & module%in%c(8L,20L)]
pooled[,`:=`(analysis_type="pooled",dataset_context="",held_out_dataset="")]
context<-rbindlist(list(pooled,context),fill=TRUE)
context[,label:=fcase(analysis_type=="pooled","Pooled sensitivity",analysis_type=="cohort",dataset_context,
  default=paste("Omit",held_out_dataset))]
context[,label:=paste0(label," (n=",n_donors,")")]
context[,program_label:=ifelse(module==8L,"Stromal ECM","Ductular injury")]
context[,`:=`(low=beta-qt(.975,residual_df)*hc3_se,high=beta+qt(.975,residual_df)*hc3_se)]
fwrite(context,file.path(out,"source_tables/OLS_source_records_sensitivity.tsv"),sep="\t")
p<-ggplot(context[estimable==TRUE & is.finite(low) & is.finite(high)],aes(beta,label))+
  geom_vline(xintercept=0,color="#9E9E9E",linewidth=.25)+
  geom_errorbar(aes(xmin=low,xmax=high),orientation="y",width=.2,linewidth=.3)+
  geom_point(size=1,color="#9E9E9E")+facet_wrap(~program_label,ncol=1)+
  labs(x="Equal-variance stage beta (nominal HC3 95% CI)",y=NULL,
    caption="Sensitivity includes GSE189600 records of unresolved independence; not primary evidence")+
  theme_classic(base_size=6)+theme(text=element_text(size=6,face="plain"))
ggsave(file.path(out,"OLS_source_records_sensitivity.pdf"),p,width=5.5,height=5,device=cairo_pdf)
writeLines(capture.output(sessionInfo()),file.path(out,"producers/sessionInfo.txt"))
