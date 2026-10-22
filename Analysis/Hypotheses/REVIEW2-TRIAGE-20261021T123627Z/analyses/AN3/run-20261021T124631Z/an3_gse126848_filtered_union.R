#!/usr/bin/env Rscript
# AN3 task 1 (PLAN A3, review B-B2): producer run of the GSE126848 within-cohort filterByExpr refit and the
# Figure 3C cross-cohort union with that refit substituted for the adopted GSE126848 calls.
# The refit code is copied from
#   figures/candidates/review-r3-fixes-20261015T112148Z/F3cdg/R03-02_gse126848_depth/gse126848_depth_sensitivity.R
# (same inputs, same seed 20260923, same call rule BH FDR < 0.05 and |log2FC| > 0.5). That script printed the
# filterByExpr call count but did not save the per-gene refit, so the union could not be rebuilt from files.
# The union rule is the one in panel-deep-fixes-20261017T035633Z/fig3/code/00_build_source_tables.py (3C).
# Also writes the single-cohort dependence table for the 1,854 pooled calls (review B-B3).
# Unit of replication: genetic-identity clusters (30 in GSE126848: 23 control, 7 disease), not verified donors.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma)})
set.seed(20260923)
base <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out <- commandArgs(TRUE)[1]
stopifnot(dir.exists(out))
arm <- file.path(base, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
cohort_fits <- file.path(base, "figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv")
prior_union <- file.path(base, "figures/candidates/review-r3-fixes-20261015T112148Z/F3cdg/R03-02_gse126848_depth/fig3c_union_under_gse126848_sensitivities.tsv")

## ---- 1. GSE126848 refit (verbatim model code from gse126848_depth_sensitivity.R) ----
dge_all <- readRDS(file.path(arm, "merged_dge.rds"))
meta <- as.data.table(readRDS(file.path(arm, "meta_matched.rds")))
sel <- meta[dataset == "GSE126848"]; setorder(sel, sample_id)
stopifnot(nrow(sel) == 30)
grp <- factor(sel$group_binary, levels = c("Control", "Disease"))
dge <- dge_all[, match(sel$sample_id, colnames(dge_all))]
call <- function(t) t$adj.P.Val < .05 & abs(t$logFC) > .5
fit <- function(d, X) {v <- voomWithQualityWeights(d, X, plot = FALSE); f <- eBayes(lmFit(v, X))
  list(t = topTable(f, coef = "groupDisease", number = Inf, sort.by = "none"), w = v$targets$sample.weights)}
X0 <- model.matrix(~grp); colnames(X0) <- c("(Intercept)", "groupDisease")
b <- fit(dge, X0); c0 <- call(b$t)
keep <- filterByExpr(dge, group = grp)
d2 <- dge[keep, , keep.lib.sizes = FALSE]; d2 <- calcNormFactors(d2, method = "RLE")
f2 <- fit(d2, X0); c2 <- call(f2$t)
refit <- data.table(gene_id_versioned = rownames(d2), logFC = f2$t$logFC, t = f2$t$t, P.Value = f2$t$P.Value,
                    FDR = f2$t$adj.P.Val, call = c2)
fwrite(refit, file.path(out, "gse126848_filterByExpr_refit_genes.tsv.gz"), sep = "\t")
base_genes <- data.table(gene_id_versioned = rownames(dge), base_logFC = b$t$logFC, base_FDR = b$t$adj.P.Val,
                         base_call = c0, total_reads = rowSums(dge$counts), filterByExpr_keep = keep)

## ---- 2. Check the base refit against the adopted per-cohort fit used by Figure 3C ----
r <- fread(cohort_fits)
adopted126 <- r[dataset == "GSE126848"]
m <- merge(base_genes, adopted126[, .(gene_id_versioned, logFC, FDR)], by = "gene_id_versioned")
check <- data.table(
  metric = c("genes_in_dge", "genes_in_adopted_GSE126848_fit", "genes_matched", "max_abs_dlogFC_base_vs_adopted",
             "base_calls", "adopted_calls", "base_calls_equal_adopted_calls"),
  value = c(nrow(base_genes), nrow(adopted126), nrow(m), max(abs(m$base_logFC - m$logFC)), sum(c0),
            sum(adopted126$FDR < 0.05 & abs(adopted126$logFC) > 0.5),
            as.integer(identical(sort(base_genes[base_call == TRUE, gene_id_versioned]),
                                 sort(adopted126[FDR < 0.05 & abs(logFC) > 0.5, gene_id_versioned])))))
print(check)

## ---- 3. GSE126848 call accounting ----
bg <- merge(base_genes, refit[, .(gene_id_versioned, refit_logFC = logFC, refit_call = call)], by = "gene_id_versioned", all.x = TRUE)
acct <- data.table(
  metric = c("filterByExpr_genes", "filterByExpr_calls", "filterByExpr_up", "filterByExpr_down",
             "adopted_calls", "adopted_up", "adopted_calls_on_genes_removed_by_filter", "adopted_calls_removed_by_filter_up",
             "adopted_calls_zero_reads_all30", "adopted_calls_lt10_reads", "refit_calls_shared_with_adopted",
             "refit_calls_shared_same_sign", "refit_calls_not_in_adopted"),
  value = c(sum(keep), sum(c2), sum(c2 & f2$t$logFC > 0), sum(c2 & f2$t$logFC < 0),
            sum(bg$base_call), sum(bg$base_call & bg$base_logFC > 0),
            sum(bg$base_call & !bg$filterByExpr_keep), sum(bg$base_call & !bg$filterByExpr_keep & bg$base_logFC > 0),
            sum(bg$base_call & bg$total_reads == 0), sum(bg$base_call & bg$total_reads < 10),
            sum(bg$base_call & bg$refit_call %in% TRUE),
            sum(bg$base_call & bg$refit_call %in% TRUE & sign(bg$base_logFC) == sign(bg$refit_logFC)),
            sum(bg$refit_call %in% TRUE & !bg$base_call)))
print(acct)
fwrite(rbind(check, acct), file.path(out, "gse126848_filterByExpr_accounting.tsv"), sep = "\t")

## ---- 4. Figure 3C union under each GSE126848 version (00_build_source_tables.py rule) ----
calls_all <- r[FDR < 0.05 & abs(logFC) > 0.5, .(gene_id_versioned, dataset, sgn = sign(logFC))]
stopifnot(nrow(calls_all) == 20054L)
other <- calls_all[dataset != "GSE126848"]
filt126 <- refit[call == TRUE, .(gene_id_versioned, dataset = "GSE126848", sgn = sign(logFC))]
classify <- function(cl, label) {
  g <- cl[, .(k = .N, nup = sum(sgn > 0), ndn = sum(sgn < 0), has126 = any(dataset == "GSE126848")), by = gene_id_versioned]
  g[, cls := fifelse(k == 1L, "unique", fifelse(nup == 0L | ndn == 0L, paste0("same", k), "mixed"))]
  multi <- g[k >= 2L]
  data.table(label = label, memberships = nrow(cl), union = nrow(g), unique = sum(g$cls == "unique"), multi = nrow(multi),
             mixed = sum(multi$cls == "mixed"), mixed_pct = round(100 * sum(multi$cls == "mixed") / nrow(multi), 2),
             mixed_with_GSE126848 = sum(multi$cls == "mixed" & multi$has126),
             same2 = sum(g$cls == "same2"), same3 = sum(g$cls == "same3"), same4 = sum(g$cls == "same4"),
             same5 = sum(g$cls == "same5"))
}
u <- rbind(classify(calls_all, "adopted calls (five cohorts)"),
           classify(rbind(other, filt126), "GSE126848 within-cohort filterByExpr refit"),
           classify(other, "four cohorts without GSE126848"))
print(u)
fwrite(u, file.path(out, "fig3c_union_filterByExpr_refit.tsv"), sep = "\t")
# side-by-side with the saved masked / depth-adjusted / four-cohort rows (prior table read, not recomputed)
pu <- fread(prior_union)
pu[, mixed_pct := round(100 * mixed / multi, 2)]
side <- rbind(pu[, .(label, source = "R03-02 saved table", memberships, union, unique, multi, mixed, mixed_pct)],
              u[label == "GSE126848 within-cohort filterByExpr refit",
                .(label, source = "AN3 this run", memberships, union, unique, multi, mixed, mixed_pct)])
print(side)
fwrite(side, file.path(out, "fig3c_union_all_gse126848_versions.tsv"), sep = "\t")
# per-gene classes for the redraw (four-cohort and filtered versions)
g_filt <- rbind(other, filt126)[, .(k = .N, nup = sum(sgn > 0), ndn = sum(sgn < 0),
                                     cohorts = paste(sort(dataset), collapse = ";")), by = gene_id_versioned]
fwrite(g_filt, file.path(out, "fig3c_gene_classes_filterByExpr_refit.tsv.gz"), sep = "\t")

## ---- 5. Single-cohort dependence of the 1,854 pooled calls (B-B3) ----
pooled <- fread(file.path(arm, "deg_results.csv"))
pc <- pooled[padj < 0.05 & abs(logFC) > 0.5, .(gene_id_versioned = gene, pooled_sgn = sign(logFC))]
stopifnot(nrow(pc) == 1854L)
dep <- rbindlist(lapply(sort(unique(r$dataset)), function(ds) {
  cc <- calls_all[dataset == ds]
  mm <- merge(pc, cc, by = "gene_id_versioned")
  data.table(cohort = ds, cohort_calls = nrow(cc), pooled_calls = nrow(pc), overlap = nrow(mm),
             overlap_same_sign = sum(mm$pooled_sgn == mm$sgn), overlap_opposite_sign = sum(mm$pooled_sgn != mm$sgn))
}))
dep[, same_sign_frac_of_pooled := round(overlap_same_sign / pooled_calls, 4)]
print(dep)
fwrite(dep, file.path(out, "pooled_1854_called_by_each_cohort_alone.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
