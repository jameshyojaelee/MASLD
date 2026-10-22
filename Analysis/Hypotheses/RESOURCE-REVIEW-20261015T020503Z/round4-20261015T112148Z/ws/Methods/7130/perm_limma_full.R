# Exhaustive GSE244832 donor-label permutation null for ATAC v3 directional replication.
# Copied from RESOURCE-REVIEW-20261015T020503Z/jobs/1111/perm_limma.R; change: enumerate all
# choose(14,9)=2002 relabelings (identity included) instead of 64 random draws; exact p = #(stat>=obs)/2002.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(Matrix); library(parallel)})
set.seed(20260923)
C <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"
res <- fread(cmd = paste("zcat", file.path(C, "da/da_peak_results.tsv.gz")))
out <- list()
for (L in c("stellate", "hepatocyte")) {
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
  null <- unlist(mclapply(seq_len(ncol(combs)), function(j) {
    cond <- rep("NORMAL", nrow(d)); cond[combs[, j]] <- "MASH"; fitstat(cond)
  }, mc.cores = 8, mc.preschedule = TRUE))
  cat(L, ": donors", nrow(d), "MASH", nm, "; discoveries", nrow(disc), "; observed", obs,
      "; relabelings", length(null), "; identity index", obs_idx, "value", null[obs_idx],
      "; n>=obs", sum(null >= obs), "; exact p", sum(null >= obs) / length(null), "\n")
  out[[L]] <- data.table(lineage = L, n_donors = nrow(d), n_mash = nm, n_disc = nrow(disc), observed = obs,
                         saved_same_sign = sum(sign(disc$logFC_gse244832) == s281),
                         n_relabelings = length(null), null_median = median(null), null_sd = sd(null),
                         null_q95 = quantile(null, 0.95), null_max = max(null), n_ge_obs = sum(null >= obs),
                         p_exact = sum(null >= obs) / length(null),
                         binomial_p_independent_peaks = binom.test(obs, nrow(disc), 0.5, alternative = "greater")$p.value)
  fwrite(data.table(lineage = L, relabeling = seq_along(null), is_observed = seq_along(null) == obs_idx, same_sign = null),
         paste0("null_full_", L, ".tsv"), sep = "\t")
}
o <- rbindlist(out)
o[, bh_q_donor_perm_2lineages := p.adjust(p_exact, "BH")]
fwrite(o, "perm_full_summary.tsv", sep = "\t")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
