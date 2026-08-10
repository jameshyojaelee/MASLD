# Is the low correlation on chr18/21/22 a DEFECT or an artifact of my metric?
# Correlation is scale-free: if PP.H4 has almost no variance on a chromosome
# (no signal anywhere), a fixed tiny absolute shift destroys r while |d| stays
# identical. Testable prediction: the "FAIL" chromosomes have much LOWER SD.
suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
S <- "2019_31311600_NAFLD_EUR"
cat(sprintf("%-8s %6s %9s %10s %10s %10s %9s\n","chr","n","cor","med|d|","SD canon","max canon","verdict"))
res <- list()
for (ch in c(13,16,20,18,21,22)) {
  f <- file.path(FM,"results/susie_coloc_rerun",S,sprintf("susie_coloc_chr%d.csv",ch))
  cf <- file.path(FM,"results/susie_coloc",S,sprintf("susie_coloc_chr%d.csv",ch))
  if (!file.exists(f) || !file.exists(cf)) next
  r <- fread(f); c0 <- fread(cf)
  m <- merge(r[, .(gg=gene, a_r=PP.H4.abf)], c0[, .(gg=gene, a_c=PP.H4.abf)], by="gg")
  m <- m[is.finite(a_r) & is.finite(a_c)]
  cr <- cor(m$a_r,m$a_c); md <- median(abs(m$a_r-m$a_c)); sdc <- sd(m$a_c); mx <- max(m$a_c)
  cat(sprintf("chr%-5d %6d %9.5f %10.6f %10.6f %10.4f %9s\n", ch, nrow(m), cr, md, sdc, mx,
      if (cr>0.99) "passed" else "FAILED"))
  res[[length(res)+1]] <- data.table(chr=ch, cr=cr, sdc=sdc, mx=mx)
}
d <- rbindlist(res)
cat(sprintf("\ncor(r, SD of canonical PP4) across these 6 chromosomes = %.4f\n", cor(d$cr, d$sdc)))
cat(sprintf("mean SD where r>0.99 : %.6f\nmean SD where r<=0.99: %.6f  (%.1fx lower)\n",
    mean(d[cr>0.99]$sdc), mean(d[cr<=0.99]$sdc), mean(d[cr>0.99]$sdc)/mean(d[cr<=0.99]$sdc)))
cat("\nIf the low-r chromosomes have far lower SD, correlation is the wrong metric\n")
cat("here and median|d| is the valid one. If SD is comparable, it is a real defect.\n")
