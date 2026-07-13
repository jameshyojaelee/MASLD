#!/usr/bin/env Rscript
# 21_meta_da.R -- DEPTH-ROBUST cross-cohort DA via per-cohort voom-QW + IVW meta.
# ============================================================================
# Naive count-pooling (18_pooled_da: ~cohort+condition edgeR) inflates: it forces
# ONE dispersion across two cohorts of ~6x-different pseudobulk depth (GSE281367
# median lib 56M vs GSE244832 9.4M), dragging dispersion down and over-calling
# (hep p<1e-4 tail 0.09%->2.03%; 47->11,648 sig). Correct approach = combine
# EFFECT SIZES, not counts:
#   1. Fit condition per cohort with limma-voomWithQualityWeights (~condition).
#      voom precision weights + array quality weights respect each cohort's depth
#      & per-sample quality (this is the project-canonical limma_voom_qw logic).
#   2. Inverse-variance fixed-effect meta per peak: beta_meta = Σw_i b_i / Σw_i,
#      w_i = 1/SE_i^2; SE_meta = 1/sqrt(Σw_i); z, p, BH-FDR. Cochran's Q + I^2 for
#      heterogeneity. Report the HONEST count: meta-sig AND same-sign in both.
# Peaks = shared GSE244832 coordinate set; effective n = donors. Env: rnaseq.
# ============================================================================
suppressPackageStartupMessages({ library(edgeR); library(limma) })
ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
G244 <- file.path(ROOT, "Analysis/ATAC/Human_Multiome/results/snapatac2")
G281 <- file.path(ROOT, "Analysis/ATAC/Human_External/pseudobulk")
OUT  <- file.path(ROOT, "Analysis/ATAC/Human_External/results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(sprintf(...), "\n", sep = "")

MASLD_SET  <- c("MASLD","MASH","MASL","NASH","NAFLD","NAFL","1")
NORMAL_SET <- c("NORMAL","HEALTHY","CONTROL","0")
to01 <- function(x){ xu<-toupper(trimws(as.character(x))); o<-rep(NA_real_,length(xu))
  o[xu %in% MASLD_SET]<-1; o[xu %in% NORMAL_SET]<-0; o }

# per-cohort voom-QW fit -> per-peak logFC + moderated SE (for that cohort's depth)
cohort_fit <- function(counts_donor_x_peak, cond01){
  mat <- t(as.matrix(counts_donor_x_peak)); storage.mode(mat) <- "double"  # peaks x donors
  grp <- factor(ifelse(cond01==1,"MASLD","NORMAL"), levels=c("NORMAL","MASLD"))
  design <- model.matrix(~ grp)                       # coef 2 = MASLD vs NORMAL
  dge <- DGEList(counts = mat)
  keep <- filterByExpr(dge, design)                   # cohort-specific expressed peaks
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  v   <- voomWithQualityWeights(dge, design)          # depth + per-donor quality weights
  fit <- eBayes(lmFit(v, design))
  co  <- colnames(design)[2]
  se  <- sqrt(fit$s2.post) * fit$stdev.unscaled[, co] # moderated SE of the coef
  data.frame(feature = rownames(fit$coefficients),
             logFC = fit$coefficients[, co], SE = se, stringsAsFactors = FALSE)
}

CTS <- c("hep","stellate","macrophage","cholangiocyte")
summ <- list()
for (ct in CTS){
  cf244 <- file.path(G244, sprintf("%s_pseudobulk_counts.tsv.gz", ct))
  cd244 <- file.path(G244, sprintf("%s_pseudobulk_coldata.tsv", ct))
  cf281 <- file.path(G281, sprintf("%s_pseudobulk_counts_GSE281367.tsv.gz", ct))
  cd281 <- file.path(G281, sprintf("%s_pseudobulk_coldata_GSE281367.tsv", ct))
  if (!all(file.exists(cf244,cd244,cf281,cd281))){ msg("[SKIP] %s", ct); next }

  a <- read.delim(cf244, row.names=1, check.names=FALSE); ca <- read.delim(cd244);
  rownames(ca)<-ca$donor_id; ca<-ca[rownames(a),]; conda<-to01(ca$condition)
  b <- read.delim(cf281, row.names=1, check.names=FALSE); cb <- read.delim(cd281)
  rownames(cb)<-cb$donor_id; cb<-cb[rownames(b),]; condb<-to01(cb$condition)

  # DEDUPE peaks: the GSE244832 cell_type_peak_sets_v2 BEDs carry ~24% coordinate-
  # identical intervals (hep 209,042 -> 158,481 unique). Identical coordinates ->
  # identical tile->peak sums -> identical count columns, so keep-first is lossless.
  # Must dedupe BEFORE fitting so the BH-FDR denominator is the #unique peaks, and
  # so the effect-size merge is 1:1 (not a many-to-many blow-up).
  a <- a[, !duplicated(colnames(a)), drop=FALSE]
  b <- b[, !duplicated(colnames(b)), drop=FALSE]
  msg("  unique peaks: 244=%d  281=%d", ncol(a), ncol(b))

  fa <- cohort_fit(a, conda); fb <- cohort_fit(b, condb)
  m  <- merge(fa, fb, by="feature", suffixes=c("_244","_281"))   # peaks testable in BOTH

  # inverse-variance fixed-effect meta
  w244 <- 1/m$SE_244^2; w281 <- 1/m$SE_281^2; W <- w244 + w281
  beta <- (w244*m$logFC_244 + w281*m$logFC_281)/W
  se   <- sqrt(1/W); z <- beta/se; p <- 2*pnorm(-abs(z)); fdr <- p.adjust(p, "BH")
  # Cochran's Q heterogeneity
  Q  <- w244*(m$logFC_244-beta)^2 + w281*(m$logFC_281-beta)^2
  I2 <- pmax(0, (Q-1)/Q)*100
  same_sign <- sign(m$logFC_244)==sign(m$logFC_281)

  res <- data.frame(feature=m$feature, logFC_meta=beta, SE_meta=se, z=z, PValue=p, FDR=fdr,
                    logFC_244=m$logFC_244, logFC_281=m$logFC_281, same_sign=same_sign,
                    Q=Q, I2=I2)
  write.csv(res, file.path(OUT, sprintf("meta_da_%s.csv", ct)), row.names=FALSE)

  n_meta   <- sum(res$FDR < 0.05, na.rm=TRUE)
  n_honest <- sum(res$FDR < 0.05 & res$same_sign, na.rm=TRUE)   # meta-sig AND concordant
  n_hetero <- sum(res$FDR < 0.05 & res$I2 > 50, na.rm=TRUE)     # sig but heterogeneous
  msg("\n=== %s ===  peaks testable in both: %d", ct, nrow(res))
  msg("  meta FDR<0.05                       : %d", n_meta)
  msg("  meta FDR<0.05 AND same-sign (HONEST): %d", n_honest)
  msg("  of those, high-heterogeneity (I2>50): %d", n_hetero)
  msg("  cross-cohort logFC Pearson r        : %.3f",
      suppressWarnings(cor(m$logFC_244, m$logFC_281, use="complete.obs")))
  summ[[ct]] <- data.frame(cell_type=ct, n_peaks=nrow(res),
    n_meta_sig=n_meta, n_honest_concordant=n_honest, n_hetero=n_hetero,
    r_logfc=round(suppressWarnings(cor(m$logFC_244,m$logFC_281,use="complete.obs")),3))
}
if (length(summ)){
  S <- do.call(rbind, summ)
  # attach n=18 baseline for the honest before/after
  S$n_sig_n18 <- c(hep=47, stellate=0, macrophage=0, cholangiocyte=0)[S$cell_type]
  write.csv(S, file.path(OUT, "meta_da_power_summary.csv"), row.names=FALSE)
  msg("\n=== HONEST POWER SUMMARY (meta-analysis, depth-robust) ==="); print(S)
  msg("n_honest_concordant = the defensible DA count (meta-sig, same direction in both cohorts).")
}
