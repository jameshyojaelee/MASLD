#!/usr/bin/env Rscript
# 22_replication_summary.R -- PAPER-READY replication framing (not pooling/meta).
# ============================================================================
# Per the reviewer-risk decision: present GSE281367 as INDEPENDENT REPLICATION of
# the GSE244832 accessibility signal, NOT as a fused pooled/meta count. This makes
# the strongest, least-attackable claim ("two independent studies agree") and adds
# no cross-study modeling assumption.
#
# For each cell type, fit condition per cohort SEPARATELY (limma voom-QW, each
# cohort's own depth/dispersion) and report:
#   - n significant FDR<0.05 in EACH cohort alone (GSE244832 is underpowered: 5 ctrl)
#   - directional concordance: of peaks the DEEP cohort (GSE281367) calls significant,
#     what fraction move the SAME direction in GSE244832 (the correct way to frame the
#     underpowered old cohort -- "concordant-but-underpowered", not "non-replicating")
#   - genome-wide logFC Pearson r (all jointly-testable peaks)
#   - overlap of the two significant sets + sign-concordance on that overlap
# Reportable = hep, stellate (abundant). mac/chol shown but FLAGGED underpowered
# (few hundred-few thousand testable peaks; do NOT headline).
# Env: rnaseq.  Fast reader (readLines) -- runs in ~1-2 min total.
# ============================================================================
suppressPackageStartupMessages({ library(edgeR); library(limma) })
ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
G244 <- file.path(ROOT, "Analysis/ATAC/Human_Multiome/results/snapatac2")
G281 <- file.path(ROOT, "Analysis/ATAC/Human_External/pseudobulk")
OUT  <- file.path(ROOT, "Analysis/ATAC/Human_External/results")
msg  <- function(...) cat(sprintf(...), "\n", sep="")

MASLD_SET  <- c("MASLD","MASH","MASL","NASH","NAFLD","NAFL","1")
NORMAL_SET <- c("NORMAL","HEALTHY","CONTROL","0")
to01 <- function(x){ xu<-toupper(trimws(as.character(x))); o<-rep(NA_real_,length(xu))
  o[xu %in% MASLD_SET]<-1; o[xu %in% NORMAL_SET]<-0; o }

fast_read_counts <- function(path){
  con <- gzfile(path,"r"); L <- readLines(con); close(con)
  peaks <- strsplit(L[1],"\t",fixed=TRUE)[[1]][-1]; body <- L[-1]; nd <- length(body)
  donors <- character(nd); mat <- matrix(0,nrow=nd,ncol=length(peaks))
  for (i in seq_len(nd)){ f<-strsplit(body[i],"\t",fixed=TRUE)[[1]]
    donors[i]<-f[1]; mat[i,]<-as.numeric(f[-1]) }
  rownames(mat)<-donors; colnames(mat)<-peaks; as.data.frame(mat,check.names=FALSE)
}

# per-cohort voom-QW fit -> per-peak logFC + FDR (that cohort's own depth/dispersion)
fit_cohort <- function(counts, cond01){
  mat <- t(as.matrix(counts)); storage.mode(mat)<-"double"
  grp <- factor(ifelse(cond01==1,"MASLD","NORMAL"), levels=c("NORMAL","MASLD"))
  design <- model.matrix(~ grp)
  dge <- DGEList(counts=mat); keep <- filterByExpr(dge, design)
  dge <- dge[keep,,keep.lib.sizes=FALSE]; dge <- calcNormFactors(dge,"TMM")
  v <- voomWithQualityWeights(dge, design); fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef=2, number=Inf, sort.by="none")
  data.frame(feature=rownames(tt), logFC=tt$logFC, FDR=tt$adj.P.Val, stringsAsFactors=FALSE)
}

CTS <- c("hep","stellate","macrophage","cholangiocyte")
rows <- list()
for (ct in CTS){
  cf244<-file.path(G244,sprintf("%s_pseudobulk_counts.tsv.gz",ct)); cd244<-file.path(G244,sprintf("%s_pseudobulk_coldata.tsv",ct))
  cf281<-file.path(G281,sprintf("%s_pseudobulk_counts_GSE281367.tsv.gz",ct)); cd281<-file.path(G281,sprintf("%s_pseudobulk_coldata_GSE281367.tsv",ct))
  if (!all(file.exists(cf244,cd244,cf281,cd281))){ msg("[SKIP] %s",ct); next }

  a<-fast_read_counts(cf244); a<-a[,!duplicated(colnames(a)),drop=FALSE]
  ca<-read.delim(cd244); rownames(ca)<-ca$donor_id; ca<-ca[rownames(a),]; conda<-to01(ca$condition)
  b<-fast_read_counts(cf281); b<-b[,!duplicated(colnames(b)),drop=FALSE]
  cb<-read.delim(cd281); rownames(cb)<-cb$donor_id; cb<-cb[rownames(b),]; condb<-to01(cb$condition)

  fa<-fit_cohort(a,conda); fb<-fit_cohort(b,condb)
  m<-merge(fa,fb,by="feature",suffixes=c("_244","_281"))          # jointly-testable peaks
  # per-peak substrate for the replication figure (logFC + FDR in EACH cohort)
  write.csv(m[,c("feature","logFC_244","FDR_244","logFC_281","FDR_281")],
            file.path(OUT, sprintf("replication_peaks_%s.csv", ct)), row.names=FALSE)
  n244<-sum(fa$FDR<0.05,na.rm=TRUE); n281<-sum(fb$FDR<0.05,na.rm=TRUE)
  # peaks significant in the DEEP cohort (281): directional concordance in the old cohort
  sig281 <- m[m$FDR_281<0.05,]
  concord_dir <- if(nrow(sig281)) mean(sign(sig281$logFC_244)==sign(sig281$logFC_281),na.rm=TRUE) else NA
  # overlap: significant in BOTH, same sign
  both <- m[m$FDR_244<0.05 & m$FDR_281<0.05,]
  n_both <- nrow(both); concord_both <- if(n_both) mean(sign(both$logFC_244)==sign(both$logFC_281),na.rm=TRUE) else NA
  r_gw <- suppressWarnings(cor(m$logFC_244,m$logFC_281,use="complete.obs"))

  msg("\n=== %s ===  jointly-testable peaks: %d", ct, nrow(m))
  msg("  GSE244832 alone (n=18, 5 ctrl)  significant FDR<0.05 : %d", n244)
  msg("  GSE281367 alone (n=12, 6/6)     significant FDR<0.05 : %d", n281)
  msg("  of GSE281367-sig peaks, same direction in GSE244832  : %.1f%% (n=%d)", 100*concord_dir, nrow(sig281))
  msg("  significant in BOTH (independent replication)         : %d (sign-concord %.3f)", n_both, concord_both)
  msg("  genome-wide logFC Pearson r                          : %.3f", r_gw)
  rows[[ct]] <- data.frame(cell_type=ct, jointly_testable=nrow(m),
    n_sig_GSE244832_alone=n244, n_sig_GSE281367_alone=n281,
    pct_281sig_dir_concordant_in_244=round(100*concord_dir,1),
    n_sig_both_replicated=n_both, sign_concord_both=round(concord_both,3),
    logFC_r_genomewide=round(r_gw,3), stringsAsFactors=FALSE)
}
if (length(rows)){
  S<-do.call(rbind,rows)
  S$reportable <- ifelse(S$jointly_testable>=10000 & S$logFC_r_genomewide>=0.3, "YES","underpowered-flag")
  write.csv(S, file.path(OUT,"replication_summary.csv"), row.names=FALSE)
  msg("\n=== REPLICATION SUMMARY (independent per-cohort + concordance) ==="); print(S)
  msg("\nReport hep+stellate as independent replication. mac/chol = underpowered (do not headline).")
  msg("Wrote %s", file.path(OUT,"replication_summary.csv"))
}
