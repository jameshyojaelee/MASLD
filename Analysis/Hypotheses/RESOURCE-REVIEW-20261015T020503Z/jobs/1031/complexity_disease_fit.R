#!/usr/bin/env Rscript
# Review R03, job 1031. Assay: bulk liver RNA-seq fragment counts (limma-voom, native to counts).
# Question: are pooled and per-cohort disease-vs-control calls sensitive to library depth/complexity?
# Unit: one library per genetic-identity cluster (774). Testing family: BH over 24,115 genes per fit.
# Reads the adopted release; writes only to this job directory.
suppressPackageStartupMessages({library(edgeR); library(limma); library(data.table)})
set.seed(20261015)
REL <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
CX <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R03/complexity_774.tsv"
OUT <- getwd()
t0 <- Sys.time()
d <- readRDS(file.path(REL, "merged_dge.rds"))
m <- as.data.table(readRDS(file.path(REL, "meta_matched.rds")))[match(colnames(d), sample_id)]
cx <- fread(CX)[match(colnames(d), sample_id)]
sw <- fread(file.path(REL, "lvqw_sample_weights.tsv"))[match(colnames(d), sample_id)]
stored <- fread(file.path(REL, "deg_results.csv"))
stopifnot(identical(m$sample_id, colnames(d)), identical(cx$sample_id, colnames(d)), identical(sw$sample_id, colnames(d)))
zc <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v)) / sd(v))
info <- data.frame(dataset = factor(m$dataset), inferred_sex = factor(m$inferred_sex),
  group_binary = factor(m$group_binary, levels = c("Control", "Disease")),
  z_loglib = zc(log2(cx$lib_size), m$dataset), z_top100 = zc(cx$top100, m$dataset),
  z_entropy = zc(cx$entropy, m$dataset), z_detected = zc(cx$detected, m$dataset), row.names = colnames(d))
call_of <- function(lfc, q) q < 0.05 & abs(lfc) > 0.5
tab <- function(fit, label) {
  tt <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
  data.table(model = label, gene = rownames(tt), logFC = tt$logFC, padj = tt$adj.P.Val)
}
fit_qw <- function(dge, design) eBayes(lmFit(voomWithQualityWeights(dge, design), design))
fit_fixed <- function(dge, design, w) {
  v <- voom(dge, design, weights = w)
  v$weights <- t(w * t(v$weights))
  eBayes(lmFit(v, design))
}
res <- list()
# B. fixed stored sample weights, baseline and covariate models
forms <- list(
  pooled_fixedw_baseline = ~ dataset + inferred_sex + group_binary,
  pooled_fixedw_loglib = ~ dataset + inferred_sex + z_loglib + group_binary,
  pooled_fixedw_top100 = ~ dataset + inferred_sex + z_top100 + group_binary,
  pooled_fixedw_entropy = ~ dataset + inferred_sex + z_entropy + group_binary,
  pooled_fixedw_loglib_top100_detected = ~ dataset + inferred_sex + z_loglib + z_top100 + z_detected + group_binary)
for (nm in names(forms)) {
  ds <- model.matrix(forms[[nm]], info); stopifnot(qr(ds)$rank == ncol(ds))
  res[[nm]] <- tab(fit_fixed(d, ds, sw$sample_weight), nm)
  cat(nm, format(Sys.time() - t0), "\n")
  fwrite(rbindlist(res), file.path(OUT, "partial_results.tsv.gz"), sep = "\t")
}
# C. per-cohort GSE126848 and GSE135251 with full recipe
for (coh in c("GSE126848", "GSE135251")) {
  idx <- which(info$dataset == coh); dd <- d[, idx]; ii <- droplevels(info[idx, ])
  zl <- (log2(cx$lib_size[idx]) - mean(log2(cx$lib_size[idx]))) / sd(log2(cx$lib_size[idx]))
  zt <- (cx$top100[idx] - mean(cx$top100[idx])) / sd(cx$top100[idx])
  ii$zl <- zl; ii$zt <- zt
  sexterm <- if (nlevels(ii$inferred_sex) > 1) "inferred_sex + " else ""
  for (cov in c("", "zl + ", "zt + ", "zl + zt + ")) {
    f <- as.formula(paste0("~ ", sexterm, cov, "group_binary"))
    ds <- model.matrix(f, ii); stopifnot(qr(ds)$rank == ncol(ds))
    nm <- paste0(coh, "_qw_", ifelse(cov == "", "baseline", gsub(" \\+ ", "_", sub(" \\+ $", "", cov))))
    res[[nm]] <- tab(fit_qw(dd, ds), nm)
    cat(nm, format(Sys.time() - t0), "\n")
  }
  # group-vs-depth separation in this cohort
  auc <- wilcox.test(cx$lib_size[idx][ii$group_binary == "Disease"], cx$lib_size[idx][ii$group_binary == "Control"])$statistic /
    (sum(ii$group_binary == "Disease") * sum(ii$group_binary == "Control"))
  cat(coh, "AUC disease libsize > control:", auc, " cor(zl, disease):", cor(zl, as.numeric(ii$group_binary == "Disease")), "\n")
}
# A. guard: full recipe on 774 (run last; skipped if elapsed time is too long)
dsg0 <- model.matrix(~ dataset + inferred_sex + group_binary, info)
if (as.numeric(difftime(Sys.time(), t0, units = "secs")) < 260) {
  res$qw_base <- tab(fit_qw(d, dsg0), "pooled_qw_baseline")
  cat("guard qw time", format(Sys.time() - t0), "\n")
} else cat("qw guard skipped for time\n")
all <- rbindlist(res)
fwrite(all, file.path(OUT, "all_gene_results.tsv.gz"), sep = "\t")
ref_pooled <- res$pooled_fixedw_baseline[, .(gene, ref_lfc = logFC, ref_call = call_of(logFC, padj))]
summ <- all[, .(n_genes = .N, calls = sum(call_of(logFC, padj)), up = sum(call_of(logFC, padj) & logFC > 0),
               down = sum(call_of(logFC, padj) & logFC < 0)), by = model]
cmp <- function(mod, ref) {
  x <- merge(all[model == mod, .(gene, logFC, call = call_of(logFC, padj))], ref, by = "gene")
  data.table(model = mod, ref_calls = sum(x$ref_call), kept_same_sign = sum(x$ref_call & x$call & sign(x$logFC) == sign(x$ref_lfc)),
             new_calls = sum(!x$ref_call & x$call), spearman_lfc_all = cor(x$logFC, x$ref_lfc, method = "spearman"))
}
c1 <- rbindlist(lapply(grep("^pooled", unique(all$model), value = TRUE), cmp, ref = ref_pooled))
c2 <- rbindlist(lapply(c("GSE126848", "GSE135251"), function(coh) {
  r <- res[[paste0(coh, "_qw_baseline")]][, .(gene, ref_lfc = logFC, ref_call = call_of(logFC, padj))]
  rbindlist(lapply(grep(paste0("^", coh), unique(all$model), value = TRUE), cmp, ref = r))
}))
st <- stored[, .(gene, ref_lfc = logFC, ref_call = call_of(logFC, padj))]
guard <- rbindlist(lapply(intersect(c("pooled_fixedw_baseline", "pooled_qw_baseline"), unique(all$model)), cmp, ref = st))
fwrite(summ, file.path(OUT, "call_counts.tsv"), sep = "\t")
fwrite(rbind(c1, c2), file.path(OUT, "retention_vs_baseline.tsv"), sep = "\t")
fwrite(guard, file.path(OUT, "guard_vs_stored_release.tsv"), sep = "\t")
print(summ); print(rbind(c1, c2)); print(guard)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("total time", format(Sys.time() - t0), "\n")
