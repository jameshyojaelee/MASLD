# Chunked version of ../7130/perm_limma_full.R: one lineage, relabelings j0..j1 of choose(14,9)=2002 (combn order), same fit.
# Exhaustive GSE244832 donor-label permutation null for ATAC v3 directional replication.
# Copied from RESOURCE-REVIEW-20261015T020503Z/jobs/1111/perm_limma.R; change: enumerate all
# choose(14,9)=2002 relabelings (identity included) instead of 64 random draws; exact p = #(stat>=obs)/2002.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(Matrix); library(parallel)})
set.seed(20260923)
C <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"
res <- fread(cmd = paste("zcat", file.path(C, "da/da_peak_results.tsv.gz")))
out <- list()
args <- commandArgs(TRUE); L <- args[1]; j0 <- as.integer(args[2]); j1 <- as.integer(args[3])
{
  r <- res[lineage == L]
  con <- gzfile(file.path(C, "counts/GSE244832", paste0(L, ".mtx.gz")), "rb")
  M <- as(as(readMM(con), "CsparseMatrix"), "dgCMatrix"); close(con)
  donors <- fread(file.path(C, "counts/GSE244832", paste0(L, ".donors.tsv")))
  peaks <- fread(file.path(C, "counts/GSE244832", paste0(L, ".peaks.tsv")))
  rownames(M) <- donors$donor_id; colnames(M) <- peaks$peak_coordinate
  keep <- donors$contrast_eligible == TRUE & donors$condition %chin% c("NORMAL", "MASH")
  X <- M[keep, r$peak_coordinate, drop = FALSE]; d <- donors[keep]
  disc <- r[qvalue_gse281367 < 0.05]
  s281 <- sign(disc$logFC_gse281367)
  y <- calcNormFactors(DGEList(counts = t(as.matrix(X))), method = "TMM")
  fitstat <- function(cond) {
    group <- factor(cond, levels = c("NORMAL", "MASH"))
    design <- model.matrix(~ group)
    v <- voomWithQualityWeights(y, design, plot = FALSE)
    fit <- eBayes(lmFit(v, design), robust = TRUE)
    lfc <- fit$coefficients[disc$peak_coordinate, 2]
    sum(sign(lfc) == s281)
  }
  obs <- fitstat(d$condition)
  nm <- sum(d$condition == "MASH")
  combs <- combn(nrow(d), nm)
  obs_idx <- which(apply(combs, 2, function(cc) setequal(cc, which(d$condition == "MASH"))))
  jj <- j0:min(j1, ncol(combs))
  null <- unlist(mclapply(jj, function(j) {
    cond <- rep("NORMAL", nrow(d)); cond[combs[, j]] <- "MASH"; fitstat(cond)
  }, mc.cores = 8, mc.preschedule = TRUE))
  fwrite(data.table(lineage = L, relabeling = jj, is_observed = jj == obs_idx, same_sign = null, observed = obs),
         paste0("null_chunk_", L, "_", j0, "_", j1, ".tsv"), sep = "\t")
}
writeLines(capture.output(sessionInfo()), paste0("sessionInfo_", j0, ".txt"))
