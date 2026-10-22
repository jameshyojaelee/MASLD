#!/usr/bin/env Rscript
# Blind check of R11-04: exact GSE244832 donor-relabeling null for the
# directional-replication count (GSE281367 q < 0.05 discovery peaks, same logFC
# sign in GSE244832). Model functions copied from
# Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/04_run_da.R
# (copy in this directory). Reads the release counts directly.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(Matrix)
  library(parallel)
})
set.seed(20260923)

args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]
NCORES <- as.integer(args[2])
DO_NULL <- length(args) < 3 || args[3] != "identity_only"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

CAND <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"

read_sparse_gz <- function(path) {
  con <- gzfile(path, "rb")
  on.exit(close(con), add = TRUE)
  as(readMM(con), "dgCMatrix")
}
load_counts <- function(cohort, lineage) {
  root <- file.path(CAND, "counts", cohort)
  counts <- read_sparse_gz(file.path(root, paste0(lineage, ".mtx.gz")))
  donors <- fread(file.path(root, paste0(lineage, ".donors.tsv")))
  peaks <- fread(file.path(root, paste0(lineage, ".peaks.tsv")))
  stopifnot(nrow(counts) == nrow(donors), ncol(counts) == nrow(peaks))
  rownames(counts) <- donors$donor_id
  colnames(counts) <- peaks$peak_coordinate
  list(counts = counts, donors = donors, peaks = peaks)
}
eligible_primary <- function(x) {
  keep <- x$donors$contrast_eligible == TRUE & x$donors$condition %chin% c("NORMAL", "MASH")
  list(counts = x$counts[keep, , drop = FALSE], donors = x$donors[keep])
}
condition_blind_filter <- function(counts) {
  y <- DGEList(counts = t(counts))
  design <- matrix(1, nrow = ncol(y), ncol = 1L)
  filterByExpr(y, design = design)
}
# Same model as fit_voom(): TMM, voomWithQualityWeights, lmFit, robust eBayes.
# Returns full table on request, otherwise only logFC (eBayes does not alter
# coefficients, so logFC is identical with or without it).
fit_voom <- function(counts, condition, full = FALSE) {
  group <- factor(condition, levels = c("NORMAL", "MASH"))
  design <- model.matrix(~ group)
  y <- DGEList(counts = t(counts))
  y <- calcNormFactors(y, method = "TMM")
  v <- voomWithQualityWeights(y, design, plot = FALSE)
  lf <- lmFit(v, design)
  if (!full) return(setNames(lf$coefficients[, 2L], rownames(lf$coefficients)))
  fit <- eBayes(lf, robust = TRUE)
  tt <- topTable(fit, coef = 2L, number = Inf, sort.by = "none")
  data.table(peak_coordinate = rownames(tt), logFC = tt$logFC,
             pvalue = tt$P.Value, qvalue = p.adjust(tt$P.Value, method = "BH"))
}

saved <- fread(cmd = paste("zcat", file.path(CAND, "da/da_peak_results.tsv.gz")))
summ <- list(); nulls <- list(); idrows <- list()

for (lineage in c("stellate", "hepatocyte")) {
  message("[check] ", lineage, " ", Sys.time())
  x244 <- eligible_primary(load_counts("GSE244832", lineage))
  x281 <- eligible_primary(load_counts("GSE281367", lineage))
  keep244 <- condition_blind_filter(x244$counts)
  keep281 <- condition_blind_filter(x281$counts)
  common <- intersect(colnames(x244$counts)[keep244], colnames(x281$counts)[keep281])
  c244 <- as.matrix(x244$counts[, common, drop = FALSE])
  c281 <- as.matrix(x281$counts[, common, drop = FALSE])
  f281 <- fit_voom(c281, x281$donors$condition, full = TRUE)
  f244 <- fit_voom(c244, x244$donors$condition, full = TRUE)
  sv <- saved[saved$lineage == lineage]
  disc <- f281[qvalue < 0.05, peak_coordinate]
  disc_saved <- sv[qvalue_gse281367 < 0.05, peak_coordinate]
  sign281 <- sign(f281$logFC[match(disc, f281$peak_coordinate)])
  obs <- sum(sign(f244$logFC[match(disc, f244$peak_coordinate)]) == sign281)
  m <- match(sv$peak_coordinate, f244$peak_coordinate)
  same_saved <- sv[qvalue_gse281367 < 0.05 & sign(logFC_gse244832) == sign(logFC_gse281367), peak_coordinate]
  same_rerun <- disc[sign(f244$logFC[match(disc, f244$peak_coordinate)]) == sign281]
  idrows[[lineage]] <- data.table(
    lineage = lineage, locale = Sys.getlocale("LC_COLLATE"),
    n_joint = length(common), n_joint_saved = nrow(sv),
    n_donor_244 = nrow(c244), n_normal_244 = sum(x244$donors$condition == "NORMAL"),
    n_mash_244 = sum(x244$donors$condition == "MASH"),
    n_disc_rerun = length(disc), n_disc_saved = length(disc_saved),
    disc_identical = setequal(disc, disc_saved),
    n_sig244_rerun = sum(f244$qvalue < 0.05),
    n_sig244_saved = sum(sv$qvalue_gse244832 < 0.05),
    max_abs_dlogFC_244 = max(abs(f244$logFC[m] - sv$logFC_gse244832)),
    max_abs_dlogFC_281 = max(abs(f281$logFC[match(sv$peak_coordinate, f281$peak_coordinate)] - sv$logFC_gse281367)),
    same_sign_rerun = obs, same_sign_saved = length(same_saved),
    rerun_subset_of_saved = all(same_rerun %in% same_saved),
    n_saved_not_rerun = sum(!same_saved %in% same_rerun)
  )
  print(idrows[[lineage]])
  if (!DO_NULL) next

  n <- nrow(c244)
  n_norm <- sum(x244$donors$condition == "NORMAL")
  combos <- combn(n, n_norm)
  obs_norm <- which(x244$donors$condition == "NORMAL")
  id_index <- which(apply(combos, 2, function(cc) setequal(cc, obs_norm)))
  stopifnot(length(id_index) == 1L, ncol(combos) == choose(n, n_norm))
  res <- mclapply(seq_len(ncol(combos)), function(j) {
    cond <- rep("MASH", n); cond[combos[, j]] <- "NORMAL"
    lfc <- fit_voom(c244, cond)
    sum(sign(lfc[disc]) == sign281)
  }, mc.cores = NCORES, mc.preschedule = TRUE)
  ss <- unlist(res)
  stopifnot(length(ss) == ncol(combos), ss[id_index] == obs)
  nulls[[lineage]] <- data.table(lineage = lineage, relabeling = seq_along(ss),
                                 is_observed = seq_along(ss) == id_index, same_sign = ss)
  summ[[lineage]] <- data.table(
    lineage = lineage, n_disc = length(disc), observed = obs,
    n_relabelings = length(ss), null_mean = mean(ss), null_median = median(ss),
    null_sd = sd(ss), binom_sd = sqrt(length(disc) / 4),
    null_q95 = unname(quantile(ss, 0.95)), null_max = max(ss),
    n_ge_observed = sum(ss >= obs), p_exact = mean(ss >= obs),
    n_ge_saved = sum(ss >= length(same_saved)), p_exact_vs_saved = mean(ss >= length(same_saved))
  )
  print(summ[[lineage]])
}

fwrite(rbindlist(idrows), file.path(OUT, "identity_check.tsv"), sep = "\t")
if (DO_NULL) {
  s <- rbindlist(summ)
  s[, bh_q := p.adjust(p_exact, "BH")]
  s[, bh_q_vs_saved := p.adjust(p_exact_vs_saved, "BH")]
  fwrite(s, file.path(OUT, "relabel_summary.tsv"), sep = "\t")
  fwrite(rbindlist(nulls), file.path(OUT, "relabel_null.tsv.gz"), sep = "\t")
  print(s)
}
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
message("[check] done ", Sys.time())
