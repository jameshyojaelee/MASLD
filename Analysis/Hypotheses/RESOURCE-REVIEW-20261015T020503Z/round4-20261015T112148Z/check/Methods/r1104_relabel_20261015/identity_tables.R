#!/usr/bin/env Rscript
# Write rerun primary fits (GSE244832, GSE281367) for hepatocyte and stellate so
# they can be compared with the saved release da_peak_results.tsv.gz.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(Matrix)})
set.seed(20260923)
OUT <- commandArgs(trailingOnly = TRUE)[1]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CAND <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"
rd <- function(p) { con <- gzfile(p, "rb"); on.exit(close(con)); as(readMM(con), "CsparseMatrix") }
ld <- function(cohort, lin) {
  r <- file.path(CAND, "counts", cohort)
  m <- rd(file.path(r, paste0(lin, ".mtx.gz"))); d <- fread(file.path(r, paste0(lin, ".donors.tsv"))); p <- fread(file.path(r, paste0(lin, ".peaks.tsv")))
  rownames(m) <- d$donor_id; colnames(m) <- p$peak_coordinate
  k <- d$contrast_eligible == TRUE & d$condition %chin% c("NORMAL", "MASH")
  list(counts = m[k, , drop = FALSE], donors = d[k])
}
flt <- function(cn) { y <- DGEList(counts = t(cn)); filterByExpr(y, design = matrix(1, ncol(y), 1L)) }
fitv <- function(cn, cond) {
  design <- model.matrix(~ factor(cond, levels = c("NORMAL", "MASH")))
  y <- calcNormFactors(DGEList(counts = t(cn)), method = "TMM")
  v <- voomWithQualityWeights(y, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design), robust = TRUE)
  tt <- topTable(fit, coef = 2L, number = Inf, sort.by = "none")
  data.table(peak_coordinate = rownames(tt), logFC = tt$logFC, pvalue = tt$P.Value, qvalue = p.adjust(tt$P.Value, "BH"),
             sample_weights = paste(signif(v$targets$sample.weights, 8), collapse = ","))
}
res <- list()
for (lin in c("stellate", "hepatocyte")) {
  a <- ld("GSE244832", lin); b <- ld("GSE281367", lin)
  common <- intersect(colnames(a$counts)[flt(a$counts)], colnames(b$counts)[flt(b$counts)])
  fa <- fitv(as.matrix(a$counts[, common, drop = FALSE]), a$donors$condition)
  fb <- fitv(as.matrix(b$counts[, common, drop = FALSE]), b$donors$condition)
  j <- merge(fa, fb, by = "peak_coordinate", suffixes = c("_gse244832", "_gse281367"))
  j[, lineage := lin]
  res[[lin]] <- j
}
fwrite(rbindlist(res), file.path(OUT, paste0("rerun_primary_", Sys.getenv("TAG"), ".tsv.gz")), sep = "\t")
writeLines(c(Sys.getlocale(), capture.output(sessionInfo())), file.path(OUT, paste0("sessionInfo_", Sys.getenv("TAG"), ".txt")))
