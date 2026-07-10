#!/usr/bin/env Rscript
# figS_supervised_disease_pca_top200DEG_clinical.R
# Supervised disease-separation PCA (top-200 DEGs by |t|, batch+sex corrected),
# AUGMENTED with NAS-score and fibrosis-stage colourings.
# 5 subplots: Cohort | Disease | Sex | NAS | Fibrosis.  2D PDF + 3D HTML.
# Outputs -> panels/pca_top200/  (new *_clinical* files; does NOT overwrite the 3-panel version)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(edgeR); library(limma); library(plotly); library(htmlwidgets)
})
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
OUT <- file.path(BASE,
  "figures/supplementary/figS_methods_validation/multimethod_validation/panels/pca_top200")
CTRL <- "#9E9E9E"; NDEG <- 200L; NAg <- "grey85"
MEGA <- c("GSE126848","GSE130970","GSE135251","GSE162694","GSE213621")
cs   <- c(GSE126848="GSE126848",GSE130970="GSE130970",GSE135251="GSE135251",GSE162694="GSE162694",GSE213621="GSE213621")
cohort_pal  <- c(GSE126848="#1F77B4",GSE130970="#FF7F0E",GSE135251="#2CA02C",GSE162694="#D62728",GSE213621="#9467BD")
disease_pal <- c(Control=CTRL,Disease="#C0392B")
sex_pal     <- c(F="#C0392B",M="#1F77B4",Unknown="grey80")
fib_pal     <- c(F0="#FEE0B6",F1="#FDB863",F2="#E08214",F3="#B35806",F4="#7F3B08",Unknown=NAg)

# --- data + sex inference (same recipe as unsupervised panels) ---------------
dge  <- load_merged_dge(); stopifnot(!is.null(dge))
samp <- as.data.table(dge$samples, keep.rownames="sample_id")
keep <- samp$dataset %in% MEGA; dge <- dge[,keep]; samp <- samp[keep]
grp  <- factor(samp$group_binary, levels=c("Control","Disease")); ds <- factor(samp$dataset)
logcpm_all <- edgeR::cpm(dge, log=TRUE, prior.count=1)
xi <- rownames(dge)[grep("^ENSG00000229807",rownames(dge))]; dy <- rownames(dge)[grep("^ENSG00000067048",rownames(dge))]
sx <- samp$sex
if (length(xi)>0 && length(dy)>0){
  es<-t(logcpm_all[c(xi[1],dy[1]),,drop=FALSE]); km<-kmeans(es,2,nstart=20,iter.max=50)
  fem<-as.integer(names(which.max(tapply(es[,1],km$cluster,mean)))); inf<-ifelse(km$cluster==fem,"F","M")
  mi<-is.na(samp$sex)|samp$sex==""; sx[mi]<-inf[mi]
}
sx[is.na(sx)|sx==""]<-"Unknown"

# --- NAS + fibrosis from unified metadata ------------------------------------
meta <- fread(file.path(BASE,"RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"),
              select=c("sample_id","fibrosis_stage","nas_score"))
samp <- merge(samp, meta, by="sample_id", all.x=TRUE, sort=FALSE)
# GSE213621 (Chen) fibrosis is COARSE (F0F1/F3F4 grouped, not true Kleiner) -> show as
# Unknown rather than mislabel as precise F1/F3 (see load_figure_data.R cohort guards).
samp$fibrosis_stage[samp$dataset == "GSE213621"] <- NA
cat(sprintf("NAS non-NA=%d  fibrosis non-NA=%d  of %d\n",
            sum(!is.na(samp$nas_score)), sum(!is.na(samp$fibrosis_stage)), nrow(samp)))

# --- supervised gene selection: top-200 by canonical DEG |t| -----------------
deg <- fread(file.path(BASE,"RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"))
deg[, base:=sub("\\..*","",gene)]; rn<-sub("\\..*","",rownames(dge))
deg <- deg[base %in% rn & is.finite(t)]; deg[, abst:=abs(t)]; setorder(deg,-abst)
idx <- match(head(deg$base,NDEG), rn); idx <- idx[!is.na(idx)]
cat(sprintf("Selected %d DEGs by |t|\n", length(idx)))

X <- limma::removeBatchEffect(logcpm_all[idx,], batch=ds,
       covariates=as.numeric(sx=="F"), design=model.matrix(~grp))
pr <- prcomp(t(X), center=TRUE, scale.=FALSE); pve<-round(100*pr$sdev^2/sum(pr$sdev^2),1)
if (mean(pr$x[grp=="Disease",1]) < mean(pr$x[grp=="Control",1])) pr$x[,1]<- -pr$x[,1]
dt <- data.table(sample_id=samp$sample_id,
  cohort=factor(cs[samp$dataset],levels=unname(cs)),
  disease=factor(samp$group_binary,levels=c("Control","Disease")), sex=sx,
  nas=as.numeric(samp$nas_score),
  fibrosis=factor(ifelse(is.na(samp$fibrosis_stage),"Unknown",paste0("F",samp$fibrosis_stage)),
                  levels=c("F0","F1","F2","F3","F4","Unknown")),
  PC1=pr$x[,1], PC2=pr$x[,2], PC3=pr$x[,3])
auc <- {r<-rank(dt$PC1); nD<-sum(dt$disease=="Disease"); nC<-sum(dt$disease=="Control")
        (sum(r[dt$disease=="Disease"])-nD*(nD+1)/2)/(nD*nC)}
cat(sprintf("PVE PC1=%.1f%% PC2=%.1f%%  PC1 disease AUC(in-sample)=%.3f\n",pve[1],pve[2],auc))
fwrite(dt, file.path(OUT,"supervised_disease_top200DEG_clinical_data.csv"))

# --- 2D: 5 panels -----------------------------------------------------------
bt <- function() theme_masld(base_size=7)+theme(axis.text=element_blank(),axis.ticks=element_blank(),
  panel.grid=element_blank(),legend.position="right",legend.title=element_text(size=6,face="plain"),
  legend.text=element_text(size=5.5),legend.key.size=unit(0.22,"cm"),plot.title=element_text(size=7.5,face="plain"))
xy <- function() labs(x=sprintf("PC1 (%.1f%%)",pve[1]),y=sprintf("PC2 (%.1f%%)",pve[2]))
pC <- ggplot(dt,aes(PC1,PC2,colour=cohort))+geom_point(size=.55,alpha=.7)+scale_colour_manual(values=cohort_pal,name="Cohort")+xy()+labs(title="by Cohort")+bt()
pD <- ggplot(dt,aes(PC1,PC2,colour=disease))+geom_point(size=.55,alpha=.7)+scale_colour_manual(values=disease_pal,name="Disease")+xy()+labs(title="by Disease")+bt()
pS <- ggplot(dt,aes(PC1,PC2,colour=sex))+geom_point(size=.55,alpha=.7)+scale_colour_manual(values=sex_pal,name="Sex")+xy()+labs(title="by Sex")+bt()
pN <- ggplot(dt[order(!is.na(nas))],aes(PC1,PC2,colour=nas))+geom_point(size=.55,alpha=.75)+
      scale_colour_viridis_c(name="NAS",na.value=NAg,limits=c(0,8))+xy()+labs(title="by NAS score")+
      guides(colour=guide_colourbar(barwidth=.4,barheight=2.6))+bt()
pF <- ggplot(dt[order(fibrosis=="Unknown")],aes(PC1,PC2,colour=fibrosis))+geom_point(size=.55,alpha=.75)+
      scale_colour_manual(values=fib_pal,name="Fibrosis")+xy()+labs(title="by Fibrosis stage")+bt()
fig2d <- (pC|pD|pS|pN|pF)+plot_annotation(
  title=sprintf("Supervised PCA (top-%d DEGs by |t|, batch+sex corrected, n=%d) coloured by clinical severity",length(idx),ncol(dge)),
  subtitle=sprintf("PC1 disease AUC=%.2f (in-sample, double-dipped — DEGs selected on these same samples; descriptive). NAS n=%d, fibrosis n=%d (grey = not available, e.g. GSE126848).",
                   auc, sum(!is.na(dt$nas)), sum(dt$fibrosis!="Unknown")),
  theme=theme(plot.title=element_text(size=9,face="plain"),plot.subtitle=element_text(size=6,colour="grey35")))
ggsave(file.path(OUT,"supervised_disease_top200DEG_clinical.pdf"),fig2d,width=16.5,height=3.5,device=cairo_pdf)
cat("Wrote supervised_disease_top200DEG_clinical.pdf\n")

# --- 3D: 5 scenes -----------------------------------------------------------
scn <- function() list(
  xaxis=list(title=sprintf("PC1 %.1f%%",pve[1]),titlefont=list(size=8),tickfont=list(size=6)),
  yaxis=list(title=sprintf("PC2 %.1f%%",pve[2]),titlefont=list(size=8),tickfont=list(size=6)),
  zaxis=list(title=sprintf("PC3 %.1f%%",pve[3]),titlefont=list(size=8),tickfont=list(size=6)),
  camera=list(eye=list(x=1.5,y=1.5,z=0.8)))
disc <- function(cby,pal) plot_ly(dt,x=~PC1,y=~PC2,z=~PC3,color=dt[[cby]],colors=pal,
  type="scatter3d",mode="markers",marker=list(size=2.4,opacity=.72,line=list(width=0))) |> layout(scene=scn())
cont <- function(){
  ok<-is.finite(dt$nas)
  plot_ly() |>
   add_trace(x=dt$PC1[!ok],y=dt$PC2[!ok],z=dt$PC3[!ok],type="scatter3d",mode="markers",
     name="NA",showlegend=FALSE,marker=list(size=2.2,color=NAg,opacity=.45)) |>
   add_trace(x=dt$PC1[ok],y=dt$PC2[ok],z=dt$PC3[ok],type="scatter3d",mode="markers",showlegend=FALSE,
     marker=list(size=2.6,color=dt$nas[ok],colorscale="Viridis",cmin=0,cmax=8,opacity=.8,
                 colorbar=list(title="NAS",len=.45,thickness=8))) |> layout(scene=scn())
}
fig3d <- subplot(disc("cohort",cohort_pal),disc("disease",disease_pal),disc("sex",sex_pal),
                 cont(),disc("fibrosis",fib_pal),
                 nrows=1,shareX=FALSE,shareY=FALSE,titleX=TRUE,titleY=TRUE) |>
  layout(title=list(text=sprintf("Supervised 3D PCA (top-%d DEGs) — Cohort | Disease | Sex | NAS | Fibrosis<br><sub>in-sample PC1 disease AUC=%.2f (double-dipped, descriptive)</sub>",length(idx),auc),
                    font=list(size=12)), margin=list(l=0,r=0,t=60,b=0),
         annotations=list(
           list(text="Cohort", x=0.09,y=1.04,xref="paper",yref="paper",showarrow=FALSE,font=list(size=10)),
           list(text="Disease",x=0.30,y=1.04,xref="paper",yref="paper",showarrow=FALSE,font=list(size=10)),
           list(text="Sex",    x=0.50,y=1.04,xref="paper",yref="paper",showarrow=FALSE,font=list(size=10)),
           list(text="NAS",    x=0.70,y=1.04,xref="paper",yref="paper",showarrow=FALSE,font=list(size=10)),
           list(text="Fibrosis",x=0.91,y=1.04,xref="paper",yref="paper",showarrow=FALSE,font=list(size=10))))
htmlwidgets::saveWidget(fig3d, file.path(OUT,"supervised_disease_top200DEG_clinical_3d.html"),
                        selfcontained=TRUE, title="supervised_disease_top200DEG_clinical_3d")
cat("Wrote supervised_disease_top200DEG_clinical_3d.html\nCLINICAL_DONE\n")
