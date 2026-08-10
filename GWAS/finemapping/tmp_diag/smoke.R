# The smoke test that should have preceded 1,100 tasks: does a COMPLETED rerun
# task agree with canonical where it must, and differ only where expected?
suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
done <- list.files(file.path(FM,"results/susie_coloc_rerun"), pattern="susie_coloc_chr.*csv$",
                   recursive=TRUE, full.names=TRUE)
cat(sprintf("completed rerun files: %d\n\n", length(done)))
ok <- 0; bad <- 0
for (f in head(done, 8)) {
  study <- basename(dirname(f)); chrf <- basename(f)
  cf <- file.path(FM,"results/susie_coloc", study, chrf)
  if (!file.exists(cf)) { cat(sprintf("  %-32s %-20s no canonical counterpart\n", study, chrf)); next }
  r <- fread(f); c0 <- fread(cf)
  g <- grep("^gene$", names(r), value=TRUE)[1]
  m <- merge(r[, .(gg=get(g), abf_r=PP.H4.abf, su_r=PP.H4.susie, n_r=n_snps)],
             c0[, .(gg=get(g), abf_c=PP.H4.abf, su_c=PP.H4.susie, n_c=n_snps)], by="gg")
  a <- m[is.finite(abf_r) & is.finite(abf_c)]
  cr <- if (nrow(a) > 2) cor(a$abf_r, a$abf_c) else NA_real_
  md <- if (nrow(a)) median(abs(a$abf_r-a$abf_c)) else NA_real_
  cross <- if (nrow(a)) sum(xor(a$abf_r>0.5, a$abf_c>0.5)) else NA_integer_
  vinc <- if (nrow(m)) 100*(median(m$n_r,na.rm=TRUE)/median(m$n_c,na.rm=TRUE)-1) else NA_real_
  # P1: cor > 0.99 and median|d| < 0.001 ; variants must INCREASE (palindromes kept)
  pass <- is.finite(cr) && cr > 0.99 && md < 0.001 && is.finite(vinc) && vinc > 5
  if (pass) ok <- ok + 1 else bad <- bad + 1
  cat(sprintf("  %-30s %-16s n=%4d cor %.5f med|d| %.6f cross %2d variants %+.1f%%  %s\n",
      study, chrf, nrow(a), cr, md, cross, vinc, if (pass) "PASS" else "*** FAIL"))
}
cat(sprintf("\nPASS %d | FAIL %d\n", ok, bad))
cat("Expected: cor>0.99, median|d|<0.001, variants +~15%% (palindromes retained).\n")
