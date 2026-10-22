# R11-04 repair: exact GSE244832 donor-relabeling null for the ATAC v3 directional statistic,
# using the release model in Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/04_run_da.R.
# Functions read_sparse_gz, load_counts, eligible_primary, condition_blind_filter and the fit inside
# fit_voom are copied from that script (lines 43-126). The first-round rerun (ws/Methods/7131) subset
# counts by the sorted da_peak_results order instead of the release `common` order and did not
# reproduce the saved logFC; this script keeps the release order and checks identity first.
# Usage: Rscript relabel_exact.R <lineage> <j0> <j1>
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(Matrix); library(parallel)})
set.seed(20260923)
CANDIDATE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"
args <- commandArgs(TRUE); L <- args[1]; j0 <- as.integer(args[2]); j1 <- as.integer(args[3])

read_sparse_gz <- function(path) {
  con <- gzfile(path, "rb"); on.exit(close(con), add = TRUE)
  as(readMM(con), "dgCMatrix")
}
load_counts <- function(cohort, lineage) {
  root <- file.path(CANDIDATE, "counts", cohort)
  counts <- read_sparse_gz(file.path(root, paste0(lineage, ".mtx.gz")))
  donors <- fread(file.path(root, paste0(lineage, ".donors.tsv")))
  peaks <- fread(file.path(root, paste0(lineage, ".peaks.tsv")))
  stopifnot(nrow(counts) == nrow(donors), ncol(counts) == nrow(peaks), !anyDuplicated(peaks$peak_coordinate))
  rownames(counts) <- donors$donor_id; colnames(counts) <- peaks$peak_coordinate
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
fit_primary <- function(counts, condition) {
  group <- factor(condition, levels = c("NORMAL", "MASH"))
  design <- model.matrix(~ group)
  if (qr(design)$rank != ncol(design)) stop("DA design is not full rank")
  y <- calcNormFactors(DGEList(counts = t(counts)), method = "TMM")
  voom <- voomWithQualityWeights(y, design, plot = FALSE)
  fit <- eBayes(lmFit(voom, design), robust = TRUE)
  table <- topTable(fit, coef = 2L, number = Inf, sort.by = "none")
  data.table(peak_coordinate = rownames(table), logFC = table$logFC,
             qvalue = p.adjust(table$P.Value, method = "BH"))
}

raw244 <- load_counts("GSE244832", L); raw281 <- load_counts("GSE281367", L)
x244 <- eligible_primary(raw244); x281 <- eligible_primary(raw281)
keep244 <- condition_blind_filter(x244$counts); keep281 <- condition_blind_filter(x281$counts)
common <- intersect(colnames(x244$counts)[keep244], colnames(x281$counts)[keep281])
counts244 <- x244$counts[, common, drop = FALSE]
cond0 <- x244$donors$condition

saved <- fread(cmd = paste("zcat", file.path(CANDIDATE, "da/da_peak_results.tsv.gz")))[lineage == L]
stopifnot(nrow(saved) == length(common), setequal(saved$peak_coordinate, common))
saved <- saved[match(common, peak_coordinate)]
disc <- saved$qvalue_gse281367 < 0.05
s281 <- sign(saved$logFC_gse281367)

stat <- function(f) {
  same <- sign(f$logFC) == s281
  c(same_sign_disc = sum(same[disc]),
    n_sig_gse244832 = sum(f$qvalue < 0.05),
    n_supported = sum(f$qvalue < 0.05 & disc & same),
    n_discordant = sum(f$qvalue < 0.05 & disc & !same))
}

f0 <- fit_primary(counts244, cond0)
ident <- data.table(lineage = L, n_donors = length(cond0), n_mash = sum(cond0 == "MASH"),
                    n_joint_peaks = length(common), n_disc_gse281367 = sum(disc),
                    max_abs_logFC_diff_vs_release = max(abs(f0$logFC - saved$logFC_gse244832)),
                    max_abs_q_diff_vs_release = max(abs(f0$qvalue - saved$qvalue_gse244832)),
                    saved_same_sign_disc = sum(sign(saved$logFC_gse244832[disc]) == s281[disc]),
                    saved_n_supported = sum(saved$evidence_state == "supported"),
                    t(stat(f0)))
print(ident)
if (j0 == 1L) fwrite(ident, paste0("identity_", L, ".tsv"), sep = "\t")

nm <- sum(cond0 == "MASH")
combs <- combn(length(cond0), nm)
obs_idx <- which(apply(combs, 2, function(cc) setequal(cc, which(cond0 == "MASH"))))
jj <- j0:min(j1, ncol(combs))
res <- mclapply(jj, function(j) {
  cond <- rep("NORMAL", length(cond0)); cond[combs[, j]] <- "MASH"
  stat(fit_primary(counts244, cond))
}, mc.cores = 8, mc.preschedule = TRUE)
out <- data.table(lineage = L, relabeling = jj, is_observed = jj == obs_idx, do.call(rbind, res))
fwrite(out, paste0("null_", L, "_", j0, "_", j1, ".tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), paste0("sessionInfo_", L, "_", j0, ".txt"))
