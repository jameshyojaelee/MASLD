#!/usr/bin/env Rscript
# 18_pooled_da_gse281367.R  --  Layer A A4 (pooled donor-level DA, power boost)
# ============================================================================
# Pool GSE244832 (18 donors) + GSE281367 (12 donors) at the DONOR PSEUDOBULK level
# and re-run the Squair-2021 edgeR-QLF differential-accessibility test with a
# COHORT covariate, at n=30. Reports n sig at n=18 (GSE244832-only baseline) vs
# n=30 (pooled) per cell type, so the power gain (esp. control arm 5 -> 11) is
# explicit. Also a within-GSE281367 replication of GSE244832's top DA peaks.
#
# Pseudobulk pooling (NOT cell-level integration) is the defensible cross-lab
# approach: the peak coordinate set is shared (GSE244832 cell_type_peak_sets_v2
# BEDs; GSE281367 counted over the SAME peaks), and cohort enters as a fixed
# covariate to absorb the batch/protocol difference. edgeR eBayes shares
# information across peaks; effective n = DONORS (never cells).
#
# Design:  ~ cohort + condition   (coef = condition = MASLD vs NORMAL)
#   MASLD = {MASL, MASH}; NORMAL control. GSE244832: 13 MASLD / 5 NORMAL.
#   GSE281367: 6 MASH / 6 NORMAL. Pooled: 19 MASLD / 11 NORMAL.
#
# Inputs (GSE244832, exist): results/snapatac2/{hep,stellate,macrophage,
#   cholangiocyte}_pseudobulk_counts.tsv.gz + _coldata.tsv
# Inputs (GSE281367, from 17-equiv pseudobulk on the SAME peaks):
#   Analysis/ATAC/Human_External/pseudobulk/{ct}_pseudobulk_counts_GSE281367.tsv.gz + _coldata
# Outputs: Analysis/ATAC/Human_External/results/pooled_da_{ct}_edger.csv
#          Analysis/ATAC/Human_External/results/pooled_da_power_summary.csv
# Env: rnaseq (edgeR).  Run baseline-only (no GSE281367 yet) to smoke-test.
# ============================================================================

suppressPackageStartupMessages({ library(edgeR) })

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

# cell types with GSE244832 pseudobulk on disk
CTS <- c(hep="hep", stellate="stellate", macrophage="macrophage", cholangiocyte="cholangiocyte")

load_counts <- function(path_counts, path_cold, cohort){
  if (!file.exists(path_counts) || !file.exists(path_cold)) return(NULL)
  cnt <- read.delim(path_counts, row.names = 1, check.names = FALSE)  # donors x peaks
  col <- read.delim(path_cold, stringsAsFactors = FALSE)
  rownames(col) <- as.character(col$donor_id)
  col <- col[rownames(cnt), , drop = FALSE]
  col$cond01 <- to01(col$condition)
  col$cohort <- cohort
  # prefix donor ids with cohort to guarantee uniqueness when pooling
  rownames(cnt) <- paste0(cohort, ":", rownames(cnt))
  rownames(col) <- rownames(cnt)
  list(counts = cnt, col = col)
}

edger_da <- function(counts_donor_x_peak, cond01, cohort=NULL){
  mat <- t(as.matrix(counts_donor_x_peak)); storage.mode(mat) <- "double"  # peaks x donors
  dge <- DGEList(counts = mat); dge <- calcNormFactors(dge, method = "TMM")
  if (is.null(cohort) || length(unique(cohort)) < 2) {
    design <- model.matrix(~ cond01); coefn <- "cond01"
  } else {
    cohort <- factor(cohort); design <- model.matrix(~ cohort + cond01); coefn <- "cond01"
  }
  dge <- estimateDisp(dge, design); fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = coefn)
  tt <- topTags(qlf, n = Inf, sort.by = "none")$table
  data.frame(feature = rownames(mat), logFC = tt$logFC, PValue = tt$PValue, FDR = tt$FDR)
}

summ <- list()
for (ct in names(CTS)) {
  tag <- CTS[[ct]]
  msg("\n=== cell type: %s ===", ct)
  a <- load_counts(file.path(G244, sprintf("%s_pseudobulk_counts.tsv.gz", tag)),
                   file.path(G244, sprintf("%s_pseudobulk_coldata.tsv", tag)), "GSE244832")
  if (is.null(a)) { msg("  [SKIP] no GSE244832 counts for %s", ct); next }
  b <- load_counts(file.path(G281, sprintf("%s_pseudobulk_counts_GSE281367.tsv.gz", tag)),
                   file.path(G281, sprintf("%s_pseudobulk_coldata_GSE281367.tsv", tag)), "GSE281367")

  # baseline: GSE244832 only (n=18)
  base <- edger_da(a$counts, a$col$cond01)
  n18  <- sum(base$FDR < 0.05, na.rm = TRUE)
  write.csv(base, file.path(OUT, sprintf("pooled_da_%s_GSE244832only.csv", ct)), row.names = FALSE)
  msg("  GSE244832-only: %d donors (%d MASLD/%d NORMAL) -> %d peaks tested, %d sig FDR<0.05",
      nrow(a$col), sum(a$col$cond01==1), sum(a$col$cond01==0), nrow(base), n18)

  row <- data.frame(cell_type=ct, n_donors_244832=nrow(a$col),
                    n_normal_244832=sum(a$col$cond01==0), n_sig_n18=n18,
                    n_donors_pooled=NA, n_normal_pooled=NA, n_sig_n30=NA,
                    repl_sign_concord=NA, stringsAsFactors=FALSE)

  if (!is.null(b)) {
    # align on common peaks (shared coordinate columns)
    common <- intersect(colnames(a$counts), colnames(b$counts))
    msg("  common peaks GSE244832 ∩ GSE281367: %d", length(common))
    A <- a$counts[, common, drop=FALSE]; B <- b$counts[, common, drop=FALSE]
    pooled_cnt <- rbind(A, B)
    pooled_col <- rbind(a$col[, c("cond01","cohort")], b$col[, c("cond01","cohort")])
    pooled <- edger_da(pooled_cnt, pooled_col$cond01, cohort = pooled_col$cohort)
    n30 <- sum(pooled$FDR < 0.05, na.rm = TRUE)
    write.csv(pooled, file.path(OUT, sprintf("pooled_da_%s_edger.csv", ct)), row.names = FALSE)
    msg("  POOLED (~cohort+condition): %d donors (%d MASLD/%d NORMAL) -> %d sig FDR<0.05  [n18=%d]",
        nrow(pooled_col), sum(pooled_col$cond01==1), sum(pooled_col$cond01==0), n30, n18)

    # within-GSE281367 replication of GSE244832's top DA peaks (sign concordance)
    top244 <- base[order(base$FDR)[1:min(500, nrow(base))], ]
    b_only <- tryCatch(edger_da(b$counts[, intersect(colnames(b$counts), top244$feature), drop=FALSE],
                                b$col$cond01), error=function(e) NULL)
    if (!is.null(b_only)) {
      m <- merge(top244[,c("feature","logFC")], b_only[,c("feature","logFC")], by="feature", suffixes=c("_244","_281"))
      concord <- mean(sign(m$logFC_244) == sign(m$logFC_281), na.rm=TRUE)
      row$repl_sign_concord <- round(concord, 3)
      msg("  within-GSE281367 replication of top-500 GSE244832 DA peaks: sign-concordance %.3f (n=%d)",
          concord, nrow(m))
    }
    row$n_donors_pooled <- nrow(pooled_col); row$n_normal_pooled <- sum(pooled_col$cond01==0); row$n_sig_n30 <- n30
  } else {
    msg("  [GSE281367 pseudobulk not present yet -> baseline-only; pooled test will run once it exists]")
  }
  summ[[ct]] <- row
}

if (length(summ)) {
  S <- do.call(rbind, summ)
  write.csv(S, file.path(OUT, "pooled_da_power_summary.csv"), row.names = FALSE)
  msg("\n=== POWER SUMMARY (n18 baseline vs n30 pooled) ==="); print(S)
  msg("Wrote %s", file.path(OUT, "pooled_da_power_summary.csv"))
}
