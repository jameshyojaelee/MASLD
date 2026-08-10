# DEFINITIVE A/B: palindromes RETAINED vs DROPPED, same eQTL fits, same LD,
# same code revision, same chromosome. The only difference is the palindrome
# change, so any delta here is attributable to it alone.
suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
all <- list()
for (s in c("MVP_NAFLD_EUR", "MVP_ALT_EUR", "MVP_ALT_AFR")) {
  fk <- file.path(FM, "results/susie_coloc_abtest_keep", s, "susie_coloc_chr22.csv")
  fd <- file.path(FM, "results/susie_coloc_abtest_drop", s, "susie_coloc_chr22.csv")
  if (!file.exists(fk) || !file.exists(fd)) { cat(sprintf("### %s -- missing arm\n", s)); next }
  k <- fread(fk); d <- fread(fd)
  g <- grep("^gene", names(k), value = TRUE)[1]
  m <- merge(k[, .(gg = get(g), abf_k = PP.H4.abf, su_k = PP.H4.susie, n_k = n_snps)],
             d[, .(gg = get(g), abf_d = PP.H4.abf, su_d = PP.H4.susie, n_d = n_snps)], by = "gg")
  cat(sprintf("\n### %s chr22 -- %d genes in both arms\n", s, nrow(m)))
  cat(sprintf("  variants/gene: KEEP %.0f vs DROP %.0f  (+%.1f%%)\n",
      median(m$n_k, na.rm=TRUE), median(m$n_d, na.rm=TRUE),
      100*(median(m$n_k,na.rm=TRUE)/median(m$n_d,na.rm=TRUE)-1)))
  a <- m[is.finite(abf_k) & is.finite(abf_d)]
  cat(sprintf("  ABF   n=%3d cor %.6f | med|d| %.6f | max|d| %.4f | >0.5: %d vs %d | gained %d lost %d\n",
      nrow(a), cor(a$abf_k,a$abf_d), median(abs(a$abf_k-a$abf_d)), max(abs(a$abf_k-a$abf_d)),
      sum(a$abf_k>0.5), sum(a$abf_d>0.5),
      sum(a$abf_k>0.5 & a$abf_d<=0.5), sum(a$abf_k<=0.5 & a$abf_d>0.5)))
  u <- m[is.finite(su_k) & is.finite(su_d)]
  if (nrow(u)) cat(sprintf("  SuSiE n=%3d cor %.6f | med|d| %.6f | max|d| %.4f | >0.5: %d vs %d | gained %d lost %d\n",
      nrow(u), cor(u$su_k,u$su_d), median(abs(u$su_k-u$su_d)), max(abs(u$su_k-u$su_d)),
      sum(u$su_k>0.5), sum(u$su_d>0.5),
      sum(u$su_k>0.5 & u$su_d<=0.5), sum(u$su_k<=0.5 & u$su_d>0.5)))
  cat(sprintf("  SuSiE coverage: KEEP %d genes vs DROP %d genes with a SuSiE result\n",
      sum(is.finite(m$su_k)), sum(is.finite(m$su_d))))
  all[[s]] <- m
}
M <- rbindlist(all, idcol = "study")
a <- M[is.finite(abf_k) & is.finite(abf_d)]; u <- M[is.finite(su_k) & is.finite(su_d)]
cat(sprintf("\n=== POOLED (%d study x gene rows) ===\n", nrow(M)))
cat(sprintf("  ABF   : n=%d cor %.6f max|d| %.4f | threshold crossings: %d\n", nrow(a),
    cor(a$abf_k,a$abf_d), max(abs(a$abf_k-a$abf_d)),
    sum(xor(a$abf_k>0.5, a$abf_d>0.5))))
cat(sprintf("  SuSiE : n=%d cor %.6f max|d| %.4f | threshold crossings: %d\n", nrow(u),
    cor(u$su_k,u$su_d), max(abs(u$su_k-u$su_d)),
    sum(xor(u$su_k>0.5, u$su_d>0.5))))
