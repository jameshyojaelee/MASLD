#!/usr/bin/env Rscript
# verify2-G2-01, job 4191. Bulk liver RNA-seq fragment counts; limma voomWithQualityWeights + eBayes (count-native),
# matching scripts/figures/build_pi_stage_extensions.R fit_contrast(): design ~ dataset + inferred_sex + group,
# full-774 DGEList norm factors carried into the subset, BH over all 24,115 genes per fit.
# Unit: genetic-identity cluster. Contrast: F1 vs F0 on the three-cohort stage roster (205 units).
suppressPackageStartupMessages({library(edgeR); library(limma); library(data.table)})
set.seed(20261015)
REL <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
OUT <- getwd(); t0 <- Sys.time(); el <- function() round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
d <- readRDS(file.path(REL, "merged_dge.rds"))
m <- as.data.table(readRDS(file.path(REL, "meta_matched.rds")))
st <- m[dataset %in% c("GSE130970","GSE135251","GSE162694") & fibrosis_stage %in% 0:1 & !is.na(inferred_sex) & !is.na(dataset)]
st[, grp := factor(paste0("F", fibrosis_stage), levels = c("F0","F1"))]
setorder(st, sample_id)
stopifnot(nrow(st) == 205)
res <- list(); genes <- list()
fitit <- function(label, sub) {
  info <- data.frame(dataset = droplevels(factor(sub$dataset)), sex = droplevels(factor(sub$inferred_sex)), grp = sub$grp, row.names = sub$sample_id)
  X <- if (nlevels(info$dataset) > 1) model.matrix(~ dataset + sex + grp, info) else model.matrix(~ sex + grp, info)
  stopifnot(qr(X)$rank == ncol(X))
  fit <- eBayes(lmFit(voomWithQualityWeights(d[, match(sub$sample_id, colnames(d))], X, plot = FALSE), X))
  tt <- topTable(fit, coef = "grpF1", number = Inf, sort.by = "none")
  q <- p.adjust(tt$P.Value, "BH")
  res[[length(res) + 1]] <<- data.table(label = label, cohorts = paste(sort(unique(sub$dataset)), collapse = ","),
    n_units = nrow(sub), n_F0 = sum(sub$grp == "F0"), n_F1 = sum(sub$grp == "F1"), design = paste(colnames(X), collapse = "+"),
    n_genes = nrow(tt), bh_q05 = sum(q < 0.05), up = sum(q < 0.05 & tt$logFC > 0), down = sum(q < 0.05 & tt$logFC < 0), t_s = el())
  genes[[length(genes) + 1]] <<- data.table(label = label, gene = rownames(tt), logFC = tt$logFC, P = tt$P.Value, q = q)
  fwrite(rbindlist(res), file.path(OUT, "f1_vs_f0_counts.tsv"), sep = "\t")
  cat(label, nrow(sub), "units; BH q<0.05:", sum(q < 0.05), "; t", el(), "s\n")
}
fitit("loco_drop_GSE130970", st[dataset != "GSE130970"])
fitit("pooled_all_three", st)
fitit("loco_drop_GSE135251", st[dataset != "GSE135251"])
fitit("loco_drop_GSE162694", st[dataset != "GSE162694"])
for (g in c("GSE130970","GSE162694","GSE135251")) { if (el() > 420) { cat("time budget; skip alone", g, "\n"); next }; fitit(paste0("alone_", g), st[dataset == g]) }
fwrite(rbindlist(genes), file.path(OUT, "f1_vs_f0_genes.tsv.gz"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("done", el(), "s\n")
