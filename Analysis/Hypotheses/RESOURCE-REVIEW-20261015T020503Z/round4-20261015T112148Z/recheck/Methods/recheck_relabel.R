# Blind recheck of R11-04 / R11-08: independent full rerun of the GSE244832 donor-relabeling null.
# Model follows the release script Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/04_run_da.R
# (eligible NORMAL/MASH donors, condition-blind filterByExpr in each cohort, common peaks in release
# intersect order, TMM, voomWithQualityWeights, robust eBayes, BH within cohort). GSE281367 fit is fixed
# at the saved release values. Usage: Rscript recheck_relabel.R <lineage> <outdir>
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(Matrix); library(parallel)})
set.seed(20260923)
args <- commandArgs(TRUE); L <- args[1]; OUTDIR <- args[2]
CAND <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"

load_cohort <- function(cohort) {
  root <- file.path(CAND, "counts", cohort)
  con <- gzfile(file.path(root, paste0(L, ".mtx.gz")), "rb")
  m <- as(readMM(con), "dgCMatrix"); close(con)
  d <- fread(file.path(root, paste0(L, ".donors.tsv")))
  p <- fread(file.path(root, paste0(L, ".peaks.tsv")))
  dimnames(m) <- list(d$donor_id, p$peak_coordinate)
  k <- d$contrast_eligible == TRUE & d$condition %in% c("NORMAL", "MASH")
  list(m = m[k, , drop = FALSE], d = d[k])
}
blind_keep <- function(m) {
  y <- DGEList(counts = t(m)); filterByExpr(y, design = matrix(1, ncol(y), 1))
}
a <- load_cohort("GSE244832"); b <- load_cohort("GSE281367")
common <- intersect(colnames(a$m)[blind_keep(a$m)], colnames(b$m)[blind_keep(b$m)])
X <- a$m[, common, drop = FALSE]
lab0 <- a$d$condition

rel <- fread(cmd = paste("zcat", file.path(CAND, "da/da_peak_results.tsv.gz")))[lineage == L]
rel <- rel[match(common, peak_coordinate)]
stopifnot(!anyNA(rel$peak_coordinate))
disc <- rel$qvalue_gse281367 < 0.05
s281 <- sign(rel$logFC_gse281367)

fit1 <- function(lab) {
  g <- factor(lab, levels = c("NORMAL", "MASH"))
  des <- model.matrix(~ g)
  y <- calcNormFactors(DGEList(counts = t(X)), method = "TMM")
  v <- voomWithQualityWeights(y, des, plot = FALSE)
  f <- eBayes(lmFit(v, des), robust = TRUE)
  lfc <- f$coefficients[, 2]
  q <- p.adjust(f$p.value[, 2], method = "BH")
  list(lfc = lfc, q = q)
}
summ <- function(f) {
  same <- sign(f$lfc) == s281
  c(same = sum(same & disc), sig244 = sum(f$q < 0.05),
    supp = sum(f$q < 0.05 & disc & same), disco = sum(f$q < 0.05 & disc & !same))
}
f0 <- fit1(lab0)
id <- data.table(lineage = L, n = length(lab0), n_mash = sum(lab0 == "MASH"), n_normal = sum(lab0 == "NORMAL"),
  n_joint = length(common), n_disc = sum(disc),
  maxdiff_lfc = max(abs(f0$lfc - rel$logFC_gse244832)), maxdiff_q = max(abs(f0$q - rel$qvalue_gse244832)),
  t(summ(f0)))
print(id)
fwrite(id, file.path(OUTDIR, paste0("identity_", L, ".tsv")), sep = "\t")

cmb <- combn(length(lab0), sum(lab0 == "MASH"))
res <- mclapply(seq_len(ncol(cmb)), function(j) {
  lab <- rep("NORMAL", length(lab0)); lab[cmb[, j]] <- "MASH"
  c(j = j, obs = as.integer(all(lab == lab0)), summ(fit1(lab)))
}, mc.cores = 16)
out <- as.data.table(do.call(rbind, res)); out[, lineage := L]
fwrite(out, file.path(OUTDIR, paste0("null_", L, ".tsv")), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, paste0("sessionInfo_", L, ".txt")))
