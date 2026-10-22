# R11-04 repair: does the first-round 845 come from fitting peaks in sorted order instead of the release `common` order?
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(Matrix)})
set.seed(20260923)
CANDIDATE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"
ld <- function(cohort, L) {
  root <- file.path(CANDIDATE, "counts", cohort); con <- gzfile(file.path(root, paste0(L, ".mtx.gz")), "rb")
  M <- as(readMM(con), "dgCMatrix"); close(con)
  d <- fread(file.path(root, paste0(L, ".donors.tsv"))); p <- fread(file.path(root, paste0(L, ".peaks.tsv")))
  rownames(M) <- d$donor_id; colnames(M) <- p$peak_coordinate
  k <- d$contrast_eligible == TRUE & d$condition %chin% c("NORMAL", "MASH"); list(M = M[k, ], d = d[k])
}
fb <- function(M) filterByExpr(DGEList(counts = t(M)), design = matrix(1, nrow(M), 1))
fit <- function(M, cond) {
  design <- model.matrix(~ factor(cond, levels = c("NORMAL", "MASH")))
  y <- calcNormFactors(DGEList(counts = t(M)), method = "TMM")
  f <- eBayes(lmFit(voomWithQualityWeights(y, design, plot = FALSE), design), robust = TRUE)
  setNames(f$coefficients[, 2], rownames(f$coefficients))
}
L <- "hepatocyte"
a <- ld("GSE244832", L); b <- ld("GSE281367", L)
common <- intersect(colnames(a$M)[fb(a$M)], colnames(b$M)[fb(b$M)])
saved <- fread(cmd = paste("zcat", file.path(CANDIDATE, "da/da_peak_results.tsv.gz")))[lineage == L]
disc <- saved[qvalue_gse281367 < 0.05]
for (ord in c("release_common_order", "sorted_da_peak_results_order")) {
  pk <- if (ord == "release_common_order") common else saved$peak_coordinate
  lfc <- fit(a$M[, pk], a$d$condition)
  cat(ord, ": same sign", sum(sign(lfc[disc$peak_coordinate]) == sign(disc$logFC_gse281367)),
      "; max |logFC - saved|", max(abs(lfc[saved$peak_coordinate] - saved$logFC_gse244832)),
      "; identical order to common:", identical(pk, common), "\n")
}
writeLines(capture.output(sessionInfo()), "sessionInfo_order_check.txt")
