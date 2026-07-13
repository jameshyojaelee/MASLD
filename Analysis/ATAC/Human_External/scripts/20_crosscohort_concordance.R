#!/usr/bin/env Rscript
# 20_crosscohort_concordance.R -- is the n=30 pooled DA gain REAL or INFLATED?
# ============================================================================
# The pooled ~cohort+condition model gives a large jump in sig peaks (e.g. hep
# 47->11,648). Two explanations:
#   (A) genuine power: both cohorts independently carry the SAME-direction MASH
#       effect, and pooling + balanced controls (5->11) surfaces it.
#   (B) inflation: edgeR pools dispersion across cohorts; if GSE281367 (deeper,
#       more cells) has lower within-group variance, the shared dispersion can be
#       dragged down and over-call significance not backed by both cohorts.
# Decisive test: for the POOLED-SIGNIFICANT peaks, do GSE244832-alone and
# GSE281367-alone agree in DIRECTION and correlate in logFC? High concordance +
# positive logFC correlation => (A). Near-0.5 sign, ~0 correlation => (B).
#
# Uses the per-peak tables 18 already saved (pooled + 244832-only) and refits
# GSE281367-alone on the common peaks. Env: rnaseq (edgeR).
# ============================================================================
suppressPackageStartupMessages({ library(edgeR) })
ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
G281 <- file.path(ROOT, "Analysis/ATAC/Human_External/pseudobulk")
RES  <- file.path(ROOT, "Analysis/ATAC/Human_External/results")
msg  <- function(...) cat(sprintf(...), "\n", sep="")

MASLD_SET  <- c("MASLD","MASH","MASL","NASH","NAFLD","NAFL","1")
NORMAL_SET <- c("NORMAL","HEALTHY","CONTROL","0")
to01 <- function(x){ xu<-toupper(trimws(as.character(x))); o<-rep(NA_real_,length(xu))
  o[xu %in% MASLD_SET]<-1; o[xu %in% NORMAL_SET]<-0; o }

edger_lfc <- function(counts_donor_x_peak, cond01){
  mat <- t(as.matrix(counts_donor_x_peak)); storage.mode(mat) <- "double"
  dge <- DGEList(counts=mat); dge <- calcNormFactors(dge, method="TMM")
  design <- model.matrix(~ cond01)
  dge <- estimateDisp(dge, design); fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef="cond01")
  tt <- topTags(qlf, n=Inf, sort.by="none")$table
  data.frame(feature=rownames(mat), logFC=tt$logFC, PValue=tt$PValue, FDR=tt$FDR)
}

CTS <- c("hep","stellate","macrophage","cholangiocyte")
out <- list()
for (ct in CTS){
  pooled_f <- file.path(RES, sprintf("pooled_da_%s_edger.csv", ct))
  base_f   <- file.path(RES, sprintf("pooled_da_%s_GSE244832only.csv", ct))
  cnt281_f <- file.path(G281, sprintf("%s_pseudobulk_counts_GSE281367.tsv.gz", ct))
  cold281_f<- file.path(G281, sprintf("%s_pseudobulk_coldata_GSE281367.tsv", ct))
  if (!all(file.exists(pooled_f, base_f, cnt281_f, cold281_f))){ msg("[SKIP] %s (missing inputs)", ct); next }

  pooled <- read.csv(pooled_f); base244 <- read.csv(base_f)
  cnt281 <- read.delim(cnt281_f, row.names=1, check.names=FALSE)
  col281 <- read.delim(cold281_f, stringsAsFactors=FALSE)
  rownames(col281) <- as.character(col281$donor_id); col281 <- col281[rownames(cnt281),,drop=FALSE]
  cond281 <- to01(col281$condition)
  fit281 <- edger_lfc(cnt281, cond281)   # GSE281367-alone logFC on its own peaks

  # merge the three on common features
  m <- Reduce(function(x,y) merge(x,y,by="feature"),
              list(setNames(pooled[,c("feature","logFC","FDR")], c("feature","lfc_pool","fdr_pool")),
                   setNames(base244[,c("feature","logFC","PValue")], c("feature","lfc_244","p_244")),
                   setNames(fit281[,c("feature","logFC","PValue")], c("feature","lfc_281","p_281"))))
  sig <- m[!is.na(m$fdr_pool) & m$fdr_pool < 0.05, ]

  # cross-cohort agreement among POOLED-SIG peaks
  concord_sig <- mean(sign(sig$lfc_244) == sign(sig$lfc_281), na.rm=TRUE)
  r_sig  <- suppressWarnings(cor(sig$lfc_244, sig$lfc_281, method="pearson", use="complete.obs"))
  rho_sig<- suppressWarnings(cor(sig$lfc_244, sig$lfc_281, method="spearman", use="complete.obs"))
  # of pooled-sig, how many also nominally replicate (p<0.05 same sign) in EACH cohort alone
  rep244 <- mean(sig$p_244 < 0.05 & sign(sig$lfc_244)==sign(sig$lfc_pool), na.rm=TRUE)
  rep281 <- mean(sig$p_281 < 0.05 & sign(sig$lfc_281)==sign(sig$lfc_pool), na.rm=TRUE)
  rep_both <- mean(sig$p_244<0.05 & sig$p_281<0.05 &
                   sign(sig$lfc_244)==sign(sig$lfc_281), na.rm=TRUE)
  # genome-wide (all common peaks) logFC correlation, as a reference
  r_all <- suppressWarnings(cor(m$lfc_244, m$lfc_281, method="pearson", use="complete.obs"))

  msg("\n=== %s ===  pooled-sig peaks: %d", ct, nrow(sig))
  msg("  cross-cohort sign-concordance (pooled-sig)   : %.3f", concord_sig)
  msg("  cross-cohort logFC Pearson r  (pooled-sig)   : %.3f", r_sig)
  msg("  cross-cohort logFC Spearman   (pooled-sig)   : %.3f", rho_sig)
  msg("  logFC Pearson r (ALL common peaks, reference): %.3f", r_all)
  msg("  pooled-sig also nominal p<.05 same-sign 244  : %.1f%%", 100*rep244)
  msg("  pooled-sig also nominal p<.05 same-sign 281  : %.1f%%", 100*rep281)
  msg("  pooled-sig replicating in BOTH (p<.05 & sign): %.1f%%", 100*rep_both)
  out[[ct]] <- data.frame(cell_type=ct, n_pooled_sig=nrow(sig),
    sign_concord_sig=round(concord_sig,3), r_sig=round(r_sig,3), rho_sig=round(rho_sig,3),
    r_all=round(r_all,3), pct_rep_244=round(100*rep244,1), pct_rep_281=round(100*rep281,1),
    pct_rep_both=round(100*rep_both,1))
}
if (length(out)){
  S <- do.call(rbind, out)
  write.csv(S, file.path(RES, "crosscohort_concordance_summary.csv"), row.names=FALSE)
  msg("\n=== CONCORDANCE SUMMARY ==="); print(S)
  msg("interpretation: sign_concord_sig >~0.8 AND r_sig >~0.4 => real shared signal;")
  msg("                near 0.5 / ~0 => pooled significance is dispersion-inflation.")
}
