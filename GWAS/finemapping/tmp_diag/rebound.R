suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
# which gene crossed, and by how much?
k <- fread(file.path(FM,"results/susie_coloc_abtest_keep/MVP_ALT_EUR/susie_coloc_chr22.csv"))
d <- fread(file.path(FM,"results/susie_coloc_abtest_drop/MVP_ALT_EUR/susie_coloc_chr22.csv"))
g <- grep("^gene", names(k), value=TRUE)[1]
m <- merge(k[, .(gg=get(g), su_k=PP.H4.susie, abf_k=PP.H4.abf)],
           d[, .(gg=get(g), su_d=PP.H4.susie, abf_d=PP.H4.abf)], by="gg")
cat("=== the SuSiE threshold crossing ===\n")
print(m[is.finite(su_k) & is.finite(su_d) & xor(su_k>0.5, su_d>0.5),
        .(gene=gg, susie_KEEP=round(su_k,4), susie_DROP=round(su_d,4),
          delta=round(su_k-su_d,4), abf_KEEP=round(abf_k,4), abf_DROP=round(abf_d,4))])
u <- m[is.finite(su_k) & is.finite(su_d)]
cat(sprintf("\n=== SuSiE delta distribution (n=%d, the best-covered study) ===\n", nrow(u)))
u[, dd := abs(su_k - su_d)]
cat(sprintf("  median %.6f | q90 %.5f | q99 %.4f | max %.4f | n with |d|>0.01: %d (%.1f%%)\n",
    median(u$dd), quantile(u$dd,.90), quantile(u$dd,.99), max(u$dd),
    sum(u$dd>0.01), 100*mean(u$dd>0.01)))
# genome-wide at-risk bound using the CORRECTED max|d|
a <- fread(file.path(FM,"results/susie_coloc/susie_coloc_all_gwas.csv"), showProgress=FALSE)
gg <- grep("gene|symbol", names(a), value=TRUE, ignore.case=TRUE)[1]
cat("\n=== corrected genome-wide bound ===\n")
for (nm in c("PP.H4.susie","PP.H4.abf")) {
  tol <- if (nm=="PP.H4.susie") 0.2139 else 0.1293
  x <- a[is.finite(get(nm))]
  near <- x[abs(get(nm)-0.5) <= tol]
  cat(sprintf("  %-13s genes >0.5 %4d | within +/-%.4f: %5d rows -> AT-RISK GENES %4d (%.1f%% of the >0.5 set)\n",
      nm, uniqueN(x[get(nm)>0.5][[gg]]), tol, nrow(near), uniqueN(near[[gg]]),
      100*uniqueN(near[[gg]])/uniqueN(x[get(nm)>0.5][[gg]])))
}
# realistic rate, not worst case: observed crossing rate among tested pairs
cat(sprintf("\n  OBSERVED crossing rate: 1 of 109 pooled SuSiE comparisons = %.2f%%\n", 100/109))
cat(sprintf("  applied to %s SuSiE-tested rows -> ~%.0f rows might cross\n",
    format(nrow(a[is.finite(PP.H4.susie)]), big.mark=","), nrow(a[is.finite(PP.H4.susie)])*1/109))
