# Corrected targeted MASH-minus-steatosis contrasts, identical available cells.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
setDTthreads(1)
root<-Sys.getenv("MASLD_PROJECT_ROOT")
task<-file.path(root,"Analysis/SingleCell/candidates/source-diagnosis-correction-20261026")
out<-file.path(task,"figures_v3/S4B_content_lineage_source_qualified.pdf")
stopifnot(!file.exists(out))
d<-fread(file.path(task,"content_lineage/lineage_cohort_effects.tsv"))[
  minimum_cells==50 & annotation_filter=="all" & contrast=="primary_steatohepatitis_vs_steatosis" &
  arm%in%c("raw","ambient_corrected_original_exposure")]
stopifnot(nrow(d)==8L,all(d$estimable))
d[,arm_label:=ifelse(arm=="raw","Raw, same cells","Ambient corrected")]
d[,label:=paste0(dataset," (",n_reference," steatosis, ",n_disease," MASH)")]
fwrite(d,file.path(task,"figures_v3/source_tables/S4B_content_lineage_source_qualified.tsv"),sep="\t")
p<-ggplot(d,aes(beta,label,color=arm_label))+
  geom_vline(xintercept=0,color="#9E9E9E",linewidth=.25)+
  geom_errorbar(aes(xmin=hc3_low,xmax=hc3_high),orientation="y",width=.15,linewidth=.3,
    position=position_dodge(width=.45))+
  geom_point(size=1.3,position=position_dodge(width=.45))+
  facet_wrap(~program,ncol=1)+scale_color_manual(values=c("Ambient corrected"="#C9265E","Raw, same cells"="#9E9E9E"))+
  labs(x="MASH − steatosis projected score (nominal HC3 95% CI)",y=NULL,color=NULL)+
  theme_classic(base_size=6)+theme(text=element_text(size=6,face="plain"),legend.position="bottom")
ggsave(out,p,width=5.5,height=3,device=cairo_pdf)
