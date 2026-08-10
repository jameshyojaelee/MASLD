# median|d| is identical everywhere but correlation drops on chr18/21/22, so the
# difference lives in the TAIL. Find the genes driving it.
suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
S <- "2019_31311600_NAFLD_EUR"
for (ch in c(22, 18, 13)) {
  f  <- file.path(FM,"results/susie_coloc_rerun",S,sprintf("susie_coloc_chr%d.csv",ch))
  cf <- file.path(FM,"results/susie_coloc",S,sprintf("susie_coloc_chr%d.csv",ch))
  r <- fread(f); c0 <- fread(cf)
  m <- merge(r[, .(gg=gene, a_r=PP.H4.abf, n_r=n_snps, top_r=top_snp)],
             c0[, .(gg=gene, a_c=PP.H4.abf, n_c=n_snps, top_c=top_snp)], by="gg")
  m <- m[is.finite(a_r) & is.finite(a_c)][, d := abs(a_r-a_c)]
  cat(sprintf("\n=== chr%d  n=%d  cor %.5f  med|d| %.6f  max|d| %.4f ===\n",
      ch, nrow(m), cor(m$a_r,m$a_c), median(m$d), max(m$d)))
  cat(sprintf("  genes with |d| > 0.01 : %d (%.2f%%)   > 0.05 : %d\n",
      sum(m$d>0.01), 100*mean(m$d>0.01), sum(m$d>0.05)))
  cat("  worst 4:\n")
  print(m[order(-d)][1:4, .(gene=gg, abf_rerun=round(a_r,4), abf_canon=round(a_c,4),
        d=round(d,4), n_rerun=n_r, n_canon=n_c, top_same=(top_r==top_c))])
}
