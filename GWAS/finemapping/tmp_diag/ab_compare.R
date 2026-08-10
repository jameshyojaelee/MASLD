# KEEP arm (palindromes retained, new code) vs CANONICAL (palindromes dropped,
# old code) on the SAME study x chr, same eQTL fits, same LD panels.
suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
for (s in c("MVP_NAFLD_EUR", "MVP_ALT_EUR", "MVP_ALT_AFR")) {
  fk <- file.path(FM, "results/susie_coloc_abtest_keep", s, "susie_coloc_chr22.csv")
  fc <- file.path(FM, "results/susie_coloc", s, "susie_coloc_chr22.csv")
  if (!file.exists(fk) || !file.exists(fc)) { cat(sprintf("\n### %s: keep=%s canon=%s -- skip\n", s, file.exists(fk), file.exists(fc))); next }
  k <- fread(fk); c0 <- fread(fc)
  gcol <- grep("^gene", names(k), value = TRUE)[1]
  cat(sprintf("\n### %s chr22  KEEP %d genes | CANON %d genes\n", s, nrow(k), nrow(c0)))
  if ("n_snps" %in% names(k))
    cat(sprintf("  median n_snps  KEEP %.0f vs CANON %.0f  (+%.1f%% variants)\n",
        median(k$n_snps, na.rm=TRUE), median(c0$n_snps, na.rm=TRUE),
        100*(median(k$n_snps,na.rm=TRUE)/median(c0$n_snps,na.rm=TRUE)-1)))
  if ("n_palindromic_kept" %in% names(k))
    cat(sprintf("  palindromes retained: median %.0f/gene | AF-dropped total %d\n",
        median(k$n_palindromic_kept, na.rm=TRUE), sum(k$n_palindromic_af_dropped, na.rm=TRUE)))
  m <- merge(k[, .(g=get(gcol), abf_k=PP.H4.abf, su_k=PP.H4.susie)],
             c0[, .(g=get(gcol), abf_c=PP.H4.abf, su_c=PP.H4.susie)], by="g")
  cat(sprintf("  genes compared: %d\n", nrow(m)))
  ab <- m[is.finite(abf_k) & is.finite(abf_c)]
  cat(sprintf("  ABF   : cor %.5f | median |delta| %.5f | max |delta| %.4f\n",
      cor(ab$abf_k, ab$abf_c), median(abs(ab$abf_k-ab$abf_c)), max(abs(ab$abf_k-ab$abf_c))))
  cat(sprintf("  ABF >0.5 : KEEP %d vs CANON %d | crossings: gained %d, lost %d\n",
      sum(ab$abf_k>0.5), sum(ab$abf_c>0.5), sum(ab$abf_k>0.5 & ab$abf_c<=0.5), sum(ab$abf_k<=0.5 & ab$abf_c>0.5)))
  su <- m[is.finite(su_k) & is.finite(su_c)]
  if (nrow(su)) {
    cat(sprintf("  SuSiE : n=%d | cor %.5f | median |delta| %.5f | max |delta| %.4f\n",
        nrow(su), cor(su$su_k, su$su_c), median(abs(su$su_k-su$su_c)), max(abs(su$su_k-su$su_c))))
    cat(sprintf("  SuSiE >0.5 : KEEP %d vs CANON %d | gained %d, lost %d\n",
        sum(su$su_k>0.5), sum(su$su_c>0.5), sum(su$su_k>0.5 & su$su_c<=0.5), sum(su$su_k<=0.5 & su$su_c>0.5)))
  } else cat(sprintf("  SuSiE : no gene has a SuSiE result in BOTH arms (KEEP %d, CANON %d)\n",
        sum(is.finite(m$su_k)), sum(is.finite(m$su_c))))
}
