# CORRECTED, NON-CIRCULAR. For a palindrome comp(a1)==a2, so:
#   {direct, complement_swapped}  => gwas allele1 sits on ld_allele1
#   {swapped, complement_direct}  => gwas allele1 sits on ld_allele2
# The relabel (direct<->complement_swapped, swapped<->complement_direct) maps
# WITHIN a letter-position class, so this mapping is invariant to it.
FM_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
src <- readLines(file.path(FM,"src/06_susie_coloc.R"))
i0 <- grep("^PANEL_AF_ROOT <- ", src)
i1 <- min(grep("^\\}$", src)[grep("^\\}$", src) > grep("^strand_verdict_vs_panel <- ", src)])
eval(parse(text=paste(src[i0:i1], collapse="\n")))
v <- fread(file.path(FM,"runs/uniform35_v4_2026-08-04/aggregate/susie_variant_results.tsv.gz"),
  select=c("study_name","chromosome","position","match_class","palindromic",
           "strand_resolution","strand_source","af","panel_af_a1","pip","primary_eligible"),
  showProgress=FALSE)
p <- unique(v[palindromic==TRUE], by=c("study_name","chromosome","position"))
p[, ref := fifelse(match_class %in% c("direct","complement_swapped"), panel_af_a1, 1-panel_af_a1)]
p[, recomputed := strand_verdict_vs_panel(af, ref)]
five <- c("FinnGen_NAFLD","FinnGen_NASH","2023_36280732_NAFLD_deCode_EUR",
          "2023_36280732_NAFLD_UKBB_EUR","2023_36280732_NAFLD_Intermountain_EUR")
s <- p[, .(n=.N,
      rec_opp   = sum(recomputed %in% "opposite"),
      rec_same  = sum(recomputed %in% "same"),
      rec_abst  = sum(is.na(recomputed)),
      run_opp   = sum(strand_resolution %in% "opposite"),
      pct_rec_opp = round(100*mean(recomputed %in% "opposite"),4),
      pct_run_opp = round(100*mean(strand_resolution %in% "opposite"),4)),
   by=.(study_name)][]
s[, affected := study_name %in% five]
setorder(s, -affected, -pct_rec_opp)
cat("=== recomputed (non-circular) vs run-recorded, palindromic variants ===\n")
print(s, nrows=40)
cat("\n=== agreement recomputed vs recorded ===\n")
q <- p[!is.na(recomputed) & !is.na(strand_resolution)]
print(q[, .(n=.N, agree=sum(recomputed==strand_resolution),
            pct=round(100*mean(recomputed==strand_resolution),3)), by=study_name][order(pct)][1:8])
cat("\n=== EUR: affected vs other ===\n")
e <- s[study_name %in% unique(p$study_name)]
a <- s[affected==TRUE]; b <- s[affected==FALSE & study_name %like% "EUR|FinnGen|PDFF|UKBB"]
cat(sprintf("  affected  : %d opp / %d pal = %.4f%%\n", sum(a$rec_opp), sum(a$n), 100*sum(a$rec_opp)/sum(a$n)))
cat(sprintf("  other EUR : %d opp / %d pal = %.4f%%\n", sum(b$rec_opp), sum(b$n), 100*sum(b$rec_opp)/sum(b$n)))
if (sum(a$n)>0 && sum(b$n)>0) {
  ft <- fisher.test(matrix(c(sum(a$rec_opp), sum(a$n)-sum(a$rec_opp),
                             sum(b$rec_opp), sum(b$n)-sum(b$rec_opp)), nrow=2))
  cat(sprintf("  Fisher OR=%.2f p=%.3g\n", ft$estimate, ft$p.value)) }
