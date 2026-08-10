suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
fs <- c(list.files(file.path(FM,"results/susie_coloc_rerun"), "susie_coloc_chr(1|2|6|19)\\.csv$",
                   recursive=TRUE, full.names=TRUE))
for (f in fs) {
  study <- basename(dirname(f)); chrf <- basename(f)
  cf <- file.path(FM,"results/susie_coloc", study, chrf)
  r <- fread(f)
  cat(sprintf("\n%s / %s\n  rerun rows: %d\n", study, chrf, nrow(r)))
  if (!file.exists(cf)) { cat("  no canonical counterpart\n"); next }
  c0 <- fread(cf)
  m <- merge(r[, .(gg=gene, a_r=PP.H4.abf, s_r=PP.H4.susie, n_r=n_snps)],
             c0[, .(gg=gene, a_c=PP.H4.abf, s_c=PP.H4.susie, n_c=n_snps)], by="gg")
  a <- m[is.finite(a_r) & is.finite(a_c)]
  cat(sprintf("  canonical rows: %d | compared: %d\n", nrow(c0), nrow(a)))
  cat(sprintf("  ABF  cor %.5f | med|d| %.6f | max|d| %.4f | crossings %d\n",
      cor(a$a_r,a$a_c), median(abs(a$a_r-a$a_c)), max(abs(a$a_r-a$a_c)),
      sum(xor(a$a_r>0.5, a$a_c>0.5))))
  u <- m[is.finite(s_r) & is.finite(s_c)]
  if (nrow(u)) cat(sprintf("  SuSiE n=%d cor %.5f | max|d| %.4f | crossings %d\n",
      nrow(u), cor(u$s_r,u$s_c), max(abs(u$s_r-u$s_c)), sum(xor(u$s_r>0.5, u$s_c>0.5))))
  cat(sprintf("  variants/gene: %+.1f%% (palindromes retained)\n",
      100*(median(m$n_r,na.rm=TRUE)/median(m$n_c,na.rm=TRUE)-1)))
  cat(sprintf("  SuSiE coverage: rerun %d vs canonical %d genes\n",
      sum(is.finite(m$s_r)), sum(is.finite(m$s_c))))
}
