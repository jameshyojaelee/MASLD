#!/usr/bin/env Rscript
# T2 supporting validation (GATE-K: supporting only, never headline; n~15/arm) —
# Is the T1 heterogeneity axis clinically anchored to the fibrosis TRAJECTORY?
# GSE193066 = 58 paired biopsies (progressor/stable/regressor). We rank its paired
# DE by t-stat and test (fgsea) whether the T1 gene partitions are enriched. The
# decisive contrast: the 481 VARIANCE-ONLY genes (DV-significant but NOT mean-shift
# DEGs — entirely invisible to standard DEG analysis) vs the mean-shift sets. If
# variance-only genes are enriched in the progression trajectory, within-condition
# heterogeneity carries clinical-fate information the mean-shift analysis misses.
suppressPackageStartupMessages({ library(data.table); library(fgsea) })
set.seed(42)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
P  <- file.path(BASE, "RNA-seq/results/heterogeneity_program/phase0")
REV <- file.path(BASE, "RNA-seq/results/reversal")

dv <- fread(file.path(P, "dv_results.csv"))
dv <- dv[!is.na(symbol) & symbol != ""]
# expression-matched negative control (review M1): variance_only genes are
# high-expression → easier to detect as DE in ANY contrast. An AveExpr-decile-
# matched control from `neither` is the fair null: if variance_only is enriched
# but the matched control is NOT, the enrichment is biology, not detection bias.
dv[, expr_decile := cut(AveExpr, quantile(AveExpr, 0:10/10, na.rm=TRUE), include.lowest=TRUE, labels=FALSE)]
vo <- dv[partition=="variance_only" & is.finite(expr_decile)]
matched <- unlist(lapply(names(table(vo$expr_decile)), function(d) {
  cand <- dv[partition=="neither" & is.finite(expr_decile) & expr_decile==as.integer(d), symbol]
  sample(cand, min(table(vo$expr_decile)[[d]], length(cand))) }))
sets <- list(
  variance_only     = unique(vo$symbol),                                   # 481 — orthogonal heterogeneity
  exprmatched_ctrl  = unique(matched),                                     # AveExpr-matched `neither` null
  mean_and_variance = unique(dv[partition=="mean_and_variance", symbol]),  # 87
  mean_only         = unique(dv[partition=="mean_only", symbol]),          # 1261 — mean-shift comparator
  dv_sig_all        = unique(dv[dv_sig==TRUE, symbol]))                    # 568
cat("gene-set sizes:\n"); for (s in names(sets)) cat(sprintf("  %-18s %d\n", s, length(sets[[s]])))

rankings <- c(progression_specific = "de_progression_specific",
              reversal_vs_progression = "de_reversal_vs_progression",
              progress_change = "de_progress_change",
              regress_change  = "de_regress_change")
out <- list()
for (rn in names(rankings)) {
  de <- fread(file.path(REV, paste0(rankings[[rn]], ".csv")))
  de <- de[!is.na(symbol) & symbol != "" & is.finite(t)]
  de <- de[order(-abs(t))][!duplicated(symbol)]            # one stat per symbol
  rk <- setNames(de$t, de$symbol)
  fg <- fgsea(pathways = sets, stats = rk, minSize = 10, maxSize = 2000, nPermSimple = 10000)
  fg[, ranking := rn]; out[[rn]] <- fg[, .(ranking, pathway, size, NES, pval, padj,
                                           leadingEdge = sapply(leadingEdge, length))]
}
R <- rbindlist(out); setorder(R, pathway, ranking)
fwrite(R[, .(ranking, pathway, size, NES=round(NES,2), pval=signif(pval,2),
             padj=signif(padj,2), n_leadingEdge=leadingEdge)],
       file.path(P, "t2_gse193066_dv_enrichment.tsv"), sep="\t")
cat("\n── GSE193066 paired-trajectory enrichment of T1 partitions ──\n")
print(R[, .(ranking, pathway, size, NES=round(NES,2), padj=signif(padj,2))])
# Honest interpretation (review C1/M1):
#  - the FAIR claim = variance_only is enriched in the progression-SPECIFIC axis,
#    AND it survives the expression-matched null (if exprmatched_ctrl is ns there).
#  - we do NOT claim "bidirectional rewind": progression_specific & reversal_vs_progression
#    share the progress_change term with opposite signs, so their sign flip is mechanical
#    (~74% of all genes flip); the independent test is r(progress_change, regress_change),
#    which is ~+0.13 for variance_only (no independent reversal). Report progression-anchor only.
vo_prog  <- R[pathway=="variance_only"    & ranking=="progression_specific"]
ctl_prog <- R[pathway=="exprmatched_ctrl" & ranking=="progression_specific"]
cat(sprintf("\nProgression-specific axis (the fair, GATE-K supporting claim):\n  variance_only     NES=%.2f padj=%.1e\n  exprmatched_ctrl  NES=%.2f padj=%.1e\n",
    vo_prog$NES, vo_prog$padj, ctl_prog$NES, ctl_prog$padj))
cat(sprintf("=> variance_only enriched in progression%s the matched control => %s.\n",
    if (vo_prog$padj<0.05 && ctl_prog$padj>=0.05) " but NOT" else " AND (check)",
    if (vo_prog$padj<0.05 && ctl_prog$padj>=0.05) "real, not a detection artifact" else "INSPECT — confound possible"))
