#!/usr/bin/env Rscript
# Review R03, job 1032. Assay: bulk liver RNA-seq fragment counts; limma voomWithQualityWeights (native to counts).
# Figure 3E stage contrasts (F1-F4 vs F0) on the 387-unit, three-cohort stage roster.
# Checks: (1) guard counts, (2) cohort x stage interaction F-test, (3) leave-one-cohort-out,
# (4) library-complexity covariates, (5) GSE213621 binned-stage contrast as held-out comparison.
# Unit: genetic-identity cluster. Family: BH over all 24,115 genes per fit. Writes only to this job dir.
suppressPackageStartupMessages({library(edgeR); library(limma); library(data.table)})
set.seed(20261015)
REL <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
STG <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/stage_sample_manifest.tsv"
CX <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R03/complexity_774.tsv"
OUT <- getwd(); t0 <- Sys.time()
el <- function() as.numeric(difftime(Sys.time(), t0, units = "secs"))
d <- readRDS(file.path(REL, "merged_dge.rds"))
m <- as.data.table(readRDS(file.path(REL, "meta_matched.rds")))
s <- fread(STG)
cx <- fread(CX)
s <- merge(s, cx[, .(sample_id, lib_size, top100, detected)], by = "sample_id")
zc <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v)) / sd(v))
s[, `:=`(z_loglib = zc(log2(lib_size), dataset), z_top100 = zc(top100, dataset), z_detected = zc(detected, dataset))]
counts_tab <- list(); genes_tab <- list()
record <- function(label, contrast, fit, coef, n) {
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  q <- if ("adj.P.Val" %in% names(tt)) tt$adj.P.Val else p.adjust(tt$P.Value, "BH")
  lfc <- if ("logFC" %in% names(tt)) tt$logFC else NA_real_
  counts_tab[[length(counts_tab) + 1]] <<- data.table(label = label, contrast = contrast, n_units = n,
    n_genes = nrow(tt), bh_sig = sum(q < 0.05), up = sum(q < 0.05 & lfc > 0, na.rm = TRUE), down = sum(q < 0.05 & lfc < 0, na.rm = TRUE),
    elapsed_s = round(el()))
  genes_tab[[length(genes_tab) + 1]] <<- data.table(label = label, contrast = contrast, gene = rownames(tt), logFC = lfc, q = q)
  fwrite(rbindlist(counts_tab), file.path(OUT, "stage_check_counts.tsv"), sep = "\t")
  fwrite(rbindlist(genes_tab), file.path(OUT, "stage_check_genes.tsv.gz"), sep = "\t")
  cat(label, contrast, "n", n, "bh_sig", sum(q < 0.05), "t", round(el()), "\n")
}
fit_one <- function(sub, form) {
  sub <- sub[order(sample_id)]
  info <- data.frame(dataset = droplevels(factor(sub$dataset)), inferred_sex = droplevels(factor(sub$inferred_sex)),
    group = factor(sub$grp, levels = c("REF", "CMP")), z_loglib = sub$z_loglib, z_top100 = sub$z_top100,
    z_detected = sub$z_detected, row.names = sub$sample_id)
  X <- model.matrix(form, info); stopifnot(qr(X)$rank == ncol(X))
  dge <- d[, match(sub$sample_id, colnames(d))]
  list(fit = eBayes(lmFit(voomWithQualityWeights(dge, X, plot = FALSE), X)), X = X, n = nrow(sub))
}
plan <- list(c("F4", "base"), c("F4", "interaction"), c("F1", "base"), c("F1", "interaction"),
             c("GSE213621", "binned"),
             c("F4", "loco"), c("F1", "loco"), c("F4", "complexity"), c("F1", "complexity"),
             c("F3", "base"), c("F3", "interaction"), c("F2", "base"), c("F2", "interaction"))
for (p in plan) {
  if (el() > 400) { cat("time budget reached; skipping", p, "\n"); next }
  k <- p[1]; what <- p[2]
  if (k == "GSE213621") {
    mm <- m[dataset == "GSE213621" & group_binary == "Disease" & fibrosis_stage %in% c(1, 3)]
    sub <- data.table(sample_id = mm$sample_id, dataset = "GSE213621", inferred_sex = mm$inferred_sex,
      grp = ifelse(mm$fibrosis_stage == 3, "CMP", "REF"), z_loglib = 0, z_top100 = 0, z_detected = 0)
    r <- fit_one(sub, ~ inferred_sex + group)
    record("GSE213621_F3F4_vs_F0F1_source_bins", "GSE213621_binned", r$fit, "groupCMP", r$n); next
  }
  sub <- copy(s[fib_group %in% c("F0", k)]); sub[, grp := ifelse(fib_group == "F0", "REF", "CMP")]
  ctr <- paste0(k, "_vs_F0")
  if (what == "base") { r <- fit_one(sub, ~ dataset + inferred_sex + group); record("baseline", ctr, r$fit, "groupCMP", r$n) }
  if (what == "complexity") { r <- fit_one(sub, ~ dataset + inferred_sex + z_loglib + z_top100 + z_detected + group); record("complexity_covariates", ctr, r$fit, "groupCMP", r$n) }
  if (what == "interaction") {
    r <- fit_one(sub, ~ inferred_sex + dataset + dataset:group)
    cn <- grep(":groupCMP$", colnames(r$X), value = TRUE)
    for (cc in cn) record(paste0("per_cohort_", sub("dataset", "", sub(":groupCMP", "", cc))), ctr, r$fit, cc, r$n)
    # heterogeneity: contrasts of per-cohort stage effects against the first cohort
    C <- matrix(0, ncol(r$X), length(cn) - 1, dimnames = list(colnames(r$X), paste0("diff_", seq_len(length(cn) - 1))))
    for (j in 2:length(cn)) { C[cn[j], j - 1] <- 1; C[cn[1], j - 1] <- -1 }
    fc <- eBayes(contrasts.fit(r$fit, C))
    record("cohort_by_stage_interaction_Ftest", ctr, fc, NULL, r$n)
  }
  if (what == "loco") for (drop in sort(unique(sub$dataset))) {
    if (el() > 420) break
    ss <- sub[dataset != drop]
    r <- fit_one(ss, ~ dataset + inferred_sex + group); record(paste0("loco_minus_", drop), ctr, r$fit, "groupCMP", r$n)
  }
}
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("done", round(el()), "s\n")
