#!/usr/bin/env Rscript
# figS_pca_clinical.R  — PCA (846 mega samples, batch+sex corrected) coloured by
# Cohort | Disease | Sex | NAS | Fibrosis.  2D PDF + 3D HTML.
# Env-parameterized gene selection:
#   PCA_SEL_MODE = "deg" (top-N canonical DEGs by |t|)  | "hvg" (top-N within-cohort HVGs)
#   PCA_NSEL     = N genes
#   PCA_OUTSUB   = output subdir under panels/      (e.g. "pca_top200" or "pca")
#   PCA_BASEFN   = output filename base             (e.g. "supervised_disease_top200DEG_clinical")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(edgeR); library(limma); library(matrixStats); library(plotly); library(htmlwidgets)
})
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
MODE   <- Sys.getenv("PCA_SEL_MODE", "deg")
NSEL   <- as.integer(Sys.getenv("PCA_NSEL", "200"))
OUTSUB <- Sys.getenv("PCA_OUTSUB", "pca_top200")
BASEFN <- Sys.getenv("PCA_BASEFN", sprintf("supervised_disease_top%dDEG_clinical", NSEL))
OUT <- file.path(BASE, "figures/supplementary/figS_methods_validation/multimethod_validation/panels", OUTSUB)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
CTRL<-"#9E9E9E"; NAg<-"grey85"
MEGA <- c("GSE126848","GSE130970","GSE135251","GSE162694","GSE213621")
cs   <- c(GSE126848="GSE126848",GSE130970="GSE130970",GSE135251="GSE135251",GSE162694="GSE162694",GSE213621="GSE213621")
cohort_pal  <- c(GSE126848="#1F77B4",GSE130970="#FF7F0E",GSE135251="#2CA02C",GSE162694="#D62728",GSE213621="#9467BD")
disease_pal <- c(Control=CTRL,Disease="#C0392B")
sex_pal     <- c(F="#C0392B",M="#1F77B4",Unknown="grey80")
fib_pal     <- c(F0="#FEE0B6",F1="#FDB863",F2="#E08214",F3="#B35806",F4="#7F3B08",Unknown=NAg)

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

meta <- fread(file.path(BASE,"RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"),
              select=c("sample_id","fibrosis_stage","nas_score"))
samp <- merge(samp, meta, by="sample_id", all.x=TRUE, sort=FALSE)
# GSE213621 (Chen) fibrosis is COARSE (study reports F0F1/F2/F3F4 groups; harmonize
# maps F0F1->1, F3F4->3) — NOT true Kleiner. Show as "Unknown" in the fibrosis panel
# rather than mislabel it as precise F1/F3 (Chen stays in the PCA for cohort/disease/
# NAS coloring; its NAS is already absent -> grey). See load_figure_data.R cohort guards.
samp$fibrosis_stage[samp$dataset == "GSE213621"] <- NA

# --- gene selection ---------------------------------------------------------
if (MODE=="hvg") {
  wc <- logcpm_all
  for (d in unique(samp$dataset)){ i<-which(samp$dataset==d); wc[,i]<-logcpm_all[,i]-rowMeans(logcpm_all[,i,drop=FALSE]) }
  idx <- order(matrixStats::rowVars(wc), decreasing=TRUE)[seq_len(NSEL)]
  sel_lab <- sprintf("Unsupervised PCA (top-%d HVGs by within-cohort variance, batch+sex corrected)", NSEL)
} else {
  deg <- fread(file.path(BASE,"RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"))
  deg[, base:=sub("\\..*","",gene)]; rn<-sub("\\..*","",rownames(dge))
  deg <- deg[base %in% rn & is.finite(t)]; deg[, abst:=abs(t)]; setorder(deg,-abst)
  idx <- match(head(deg$base,NSEL), rn); idx <- idx[!is.na(idx)]
  sel_lab <- sprintf("Supervised PCA (top-%d DEGs by |t|, batch+sex corrected)", length(idx))
}
cat(sprintf("MODE=%s  selected %d genes  OUT=%s\n", MODE, length(idx), OUT))

X  <- limma::removeBatchEffect(logcpm_all[idx,], batch=ds, covariates=as.numeric(sx=="F"), design=model.matrix(~grp))
pr <- prcomp(t(X), center=TRUE, scale.=FALSE); pve<-round(100*pr$sdev^2/sum(pr$sdev^2),1)
if (mean(pr$x[grp=="Disease",1]) < mean(pr$x[grp=="Control",1])) pr$x[,1]<- -pr$x[,1]
dt <- data.table(sample_id=samp$sample_id,
  cohort=factor(cs[samp$dataset],levels=unname(cs)),
  disease=factor(samp$group_binary,levels=c("Control","Disease")), sex=sx,
  nas=as.numeric(samp$nas_score),
  fibrosis=factor(ifelse(is.na(samp$fibrosis_stage),"Unknown",paste0("F",samp$fibrosis_stage)),
                  levels=c("F0","F1","F2","F3","F4","Unknown")),
  PC1=pr$x[,1],PC2=pr$x[,2],PC3=pr$x[,3])
auc <- {r<-rank(dt$PC1); nD<-sum(dt$disease=="Disease"); nC<-sum(dt$disease=="Control")
        (sum(r[dt$disease=="Disease"])-nD*(nD+1)/2)/(nD*nC)}
cat(sprintf("PVE PC1=%.1f%% PC2=%.1f%%  PC1 disease AUC=%.3f\n",pve[1],pve[2],auc))
fwrite(dt, file.path(OUT, paste0(BASEFN,"_data.csv")))

bt <- function() theme_masld(base_size=6)+theme(axis.text=element_blank(),axis.ticks=element_blank(),
  panel.grid=element_blank(),legend.position="right",legend.title=element_text(size=6,face="plain"),
  legend.text=element_text(size=6),legend.key.size=unit(0.22,"cm"))
xy <- function() labs(x=sprintf("PC1 (%.1f%%)",pve[1]),y=sprintf("PC2 (%.1f%%)",pve[2]))
pC<-ggplot(dt,aes(PC1,PC2,colour=cohort))+geom_point(size=.55,alpha=.7)+scale_colour_manual(values=cohort_pal,name="Cohort")+xy()+bt()
pD<-ggplot(dt,aes(PC1,PC2,colour=disease))+geom_point(size=.55,alpha=.7)+scale_colour_manual(values=disease_pal,name="Disease")+xy()+bt()
pS<-ggplot(dt,aes(PC1,PC2,colour=sex))+geom_point(size=.55,alpha=.7)+scale_colour_manual(values=sex_pal,name="Sex")+xy()+bt()
# NAS / fibrosis: SHAPE encodes group (open=control, filled=disease), COLOUR the
# score (gray when unscored). Scored control = open coloured ring; scored disease
# = filled coloured dot; unscored = gray (open if control, filled if disease).
pN<-ggplot(dt[order(!is.na(nas))],aes(PC1,PC2,colour=nas,shape=disease))+
    geom_point(size=.7,stroke=.3,alpha=.8)+
    scale_colour_gradient(name="NAS",low="#fad0ce",high="#741816",limits=c(0,8),na.value=NAg)+
    scale_shape_manual(name=NULL,values=c(Control=1,Disease=16))+
    guides(colour=guide_colourbar(order=1,barwidth=.4,barheight=2.6),
           shape=guide_legend(order=2,override.aes=list(size=1.8,colour="grey30")))+
    xy()+bt()
pF<-ggplot(dt[order(fibrosis=="Unknown",decreasing=TRUE)],aes(PC1,PC2,colour=fibrosis,shape=disease))+
    geom_point(size=.7,stroke=.3,alpha=.85)+
    scale_colour_manual(values=fib_pal,name="Fibrosis")+
    scale_shape_manual(name=NULL,values=c(Control=1,Disease=16))+
    guides(colour=guide_legend(order=1,override.aes=list(size=1.8,shape=16)),
           shape=guide_legend(order=2,override.aes=list(size=1.8,colour="grey30")))+
    xy()+bt()
sub <- if (MODE=="hvg")
  sprintf("Variance-selected genes (unsupervised) — disease/NAS/fibrosis are NOT high-variance axes, so they do not separate. PC1 disease AUC=%.2f. NAS n=%d, fibrosis n=%d (grey=N/A).",auc,sum(!is.na(dt$nas)),sum(dt$fibrosis!="Unknown")) else
  sprintf("PC1 disease AUC=%.2f (in-sample, double-dipped — DEGs selected on these same samples; descriptive). NAS n=%d, fibrosis n=%d (grey=N/A).",auc,sum(!is.na(dt$nas)),sum(dt$fibrosis!="Unknown"))
message(sprintf("[caption] %s, n=%d. %s", sel_lab, ncol(dge), sub))
fig2d <- (pC|pS|pD|pN|pF)
ggsave(file.path(OUT, paste0(BASEFN,".pdf")), fig2d, width=fig_full_width, height=1.5, device=cairo_pdf)
cat("Wrote", paste0(BASEFN,".pdf"), "\n")

scn<-function() list(xaxis=list(title=sprintf("PC1 %.1f%%",pve[1]),titlefont=list(size=8),tickfont=list(size=6)),
  yaxis=list(title=sprintf("PC2 %.1f%%",pve[2]),titlefont=list(size=8),tickfont=list(size=6)),
  zaxis=list(title=sprintf("PC3 %.1f%%",pve[3]),titlefont=list(size=8),tickfont=list(size=6)),
  camera=list(eye=list(x=1.5,y=1.5,z=0.8)))
disc<-function(cby,pal) plot_ly(dt,x=~PC1,y=~PC2,z=~PC3,color=dt[[cby]],colors=pal,type="scatter3d",mode="markers",
  marker=list(size=2.4,opacity=.72,line=list(width=0))) |> layout(scene=scn())
cont<-function(){ ctl<-dt$disease=="Control"; mis<-dt$disease=="Disease"&!is.finite(dt$nas); ok<-is.finite(dt$nas)
  plot_ly() |>
   add_trace(x=dt$PC1[mis],y=dt$PC2[mis],z=dt$PC3[mis],type="scatter3d",mode="markers",name="N/A",showlegend=FALSE,
     marker=list(size=2.2,color=NAg,opacity=.4)) |>
   add_trace(x=dt$PC1[ok],y=dt$PC2[ok],z=dt$PC3[ok],type="scatter3d",mode="markers",showlegend=FALSE,
     marker=list(size=2.6,color=dt$nas[ok],colorscale="Reds",cmin=0,cmax=8,opacity=.85,
                 colorbar=list(title="NAS",len=.45,thickness=8))) |>
   add_trace(x=dt$PC1[ctl],y=dt$PC2[ctl],z=dt$PC3[ctl],type="scatter3d",mode="markers",name="Control",showlegend=FALSE,
     marker=list(size=2.6,color="rgba(0,0,0,0)",opacity=.85,line=list(color="#9E9E9E",width=1.2))) |>
   layout(scene=scn()) }
fig3d <- subplot(disc("cohort",cohort_pal),disc("sex",sex_pal),disc("disease",disease_pal),cont(),disc("fibrosis",fib_pal),
  nrows=1,shareX=FALSE,shareY=FALSE,titleX=TRUE,titleY=TRUE) |>
  layout(title=list(text=sprintf("%s — Cohort | Sex | Disease | NAS | Fibrosis<br><sub>PC1 disease AUC=%.2f</sub>",sel_lab,auc),font=list(size=11)),
         margin=list(l=0,r=0,t=60,b=0),
         annotations=list(
           list(text="Cohort",x=0.09,y=1.04,xref="paper",yref="paper",showarrow=FALSE,font=list(size=10)),
           list(text="Sex",x=0.30,y=1.04,xref="paper",yref="paper",showarrow=FALSE,font=list(size=10)),
           list(text="Disease",x=0.50,y=1.04,xref="paper",yref="paper",showarrow=FALSE,font=list(size=10)),
           list(text="NAS",x=0.70,y=1.04,xref="paper",yref="paper",showarrow=FALSE,font=list(size=10)),
           list(text="Fibrosis",x=0.91,y=1.04,xref="paper",yref="paper",showarrow=FALSE,font=list(size=10))))
htmlwidgets::saveWidget(fig3d, file.path(OUT, paste0(BASEFN,"_3d.html")), selfcontained=TRUE, title=BASEFN)
cat("Wrote", paste0(BASEFN,"_3d.html"), "\nCLINICAL_DONE\n")
