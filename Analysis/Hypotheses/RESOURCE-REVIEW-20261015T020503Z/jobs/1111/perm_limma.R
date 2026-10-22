# Exact-pipeline donor-label permutation null for ATAC v3 directional replication.
# Reproduces 04_run_da.R fit_voom (TMM + voomWithQualityWeights + eBayes robust) on GSE244832,
# restricted to the jointly tested peaks, then counts sign agreement with GSE281367 BH discoveries.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(Matrix); library(parallel)})
set.seed(20261015)
C <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"
res <- fread(cmd = paste("zcat", file.path(C, "da/da_peak_results.tsv.gz")))
t0 <- Sys.time(); budget <- 5.5 * 60
out <- list()
for (L in c("stellate", "hepatocyte")) {
  r <- res[lineage == L]
  con <- gzfile(file.path(C, "counts/GSE244832", paste0(L, ".mtx.gz")), "rb")
  M <- as(readMM(con), "dgCMatrix"); close(con)
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
  cat(L, ": donors", nrow(d), "MASH", sum(d$condition == "MASH"), "; discoveries", nrow(disc),
      "; reproduced same-sign", obs, "(saved", sum(sign(disc$logFC_gse244832) == s281), ")\n")
  combs <- combn(nrow(d), sum(d$condition == "MASH"))
  ord <- sample(ncol(combs))
  per_lineage_budget <- if (L == "stellate") 1.2 * 60 else budget - as.numeric(difftime(Sys.time(), t0, units = "secs"))
  t1 <- Sys.time(); null <- integer(0); i <- 0L; chunk <- 16L
  while (i < length(ord) && as.numeric(difftime(Sys.time(), t1, units = "secs")) < per_lineage_budget) {
    idx <- ord[(i + 1L):min(i + chunk, length(ord))]
    vals <- unlist(mclapply(idx, function(j) {
      cond <- rep("NORMAL", nrow(d)); cond[combs[, j]] <- "MASH"; fitstat(cond)
    }, mc.cores = 4))
    null <- c(null, vals); i <- i + length(idx)
  }
  cat(L, ": permutations completed", length(null), "of", ncol(combs),
      "; null median", median(null), "sd", round(sd(null), 1), "p95", quantile(null, 0.95), "max", max(null),
      "; draws >= observed", sum(null >= obs), "; p =", round((sum(null >= obs) + 1) / (length(null) + 1), 4),
      "; binomial sd under independence", round(sqrt(nrow(disc) * 0.25), 1), "\n")
  out[[L]] <- data.table(lineage = L, n_disc = nrow(disc), observed = obs, n_perm = length(null),
                         null_median = median(null), null_sd = sd(null), null_p95 = quantile(null, 0.95),
                         null_max = max(null), n_ge_obs = sum(null >= obs),
                         p_perm = (sum(null >= obs) + 1) / (length(null) + 1))
  fwrite(data.table(lineage = L, null = null), paste0("null_", L, ".tsv"), sep = "\t")
}
fwrite(rbindlist(out), "perm_summary.tsv", sep = "\t")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
