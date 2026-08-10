#!/usr/bin/env Rscript
# ALL re-oriented palindromes in the frozen v4 release (not just PIP>=0.10):
# margin distribution, survival at AF_CERT_OVERRIDE_MARGIN=0.50, anchor overlap.
suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
OUT <- file.path(FM,"tmp_diag/all190"); RR <- file.path(FM,"runs/uniform35_v4_2026-08-04")
v <- fread(file.path(RR,"aggregate/susie_variant_results.tsv.gz"),
  select=c("locus_index","locus_id","study_name","ancestry","chromosome","position","rsid",
           "effect_allele","other_allele","beta","se","pval","match_class","palindromic",
           "strand_resolution","strand_source","af","panel_af_a1","pip","cs_id",
           "primary_eligible","secondary_eligible"), showProgress=FALSE)
t <- v[palindromic==TRUE & match_class %like% "^complement_"]
cat("re-oriented palindromic variants (ALL PIP):", nrow(t), "\n")
cat("strand_source tally:\n"); print(t[, .N, by=strand_source])
# relabel-invariant letter-position grouping
t[, af_on_effect := fifelse(match_class %in% c("direct","complement_swapped"), af, 1-af)]
t[, d_same := abs(af_on_effect - panel_af_a1)]
t[, d_opp  := abs(af_on_effect - (1-panel_af_a1))]
t[, margin := d_same - d_opp]
t[, survives_050 := margin >= 0.50]
fwrite(t, file.path(OUT,"all_reoriented.tsv"), sep="\t")

cat("\n=== (1) MARGIN DISTRIBUTION, all", nrow(t), "===\n")
print(summary(t$margin))
cat("\n  deciles:\n"); print(round(quantile(t$margin, probs=seq(0,1,0.1), na.rm=TRUE),4))
cat("\n  survival at thresholds:\n")
for (th in c(0.10,0.20,0.30,0.40,0.50,0.60,0.70)) {
  cat(sprintf("    margin >= %.2f : %4d of %d kept (%.1f%%)  -> %4d dropped\n",
      th, sum(t$margin>=th,na.rm=TRUE), nrow(t), 100*mean(t$margin>=th,na.rm=TRUE),
      sum(t$margin<th,na.rm=TRUE)))
}
cat("\n  by PIP stratum:\n")
t[, pip_stratum := fifelse(pip>=0.10, "PIP>=0.10 (the 28)", "PIP<0.10 (the rest)")]
print(t[, .(n=.N, median_margin=round(median(margin),4), surv_050=sum(survives_050),
            pct_surv=round(100*mean(survives_050),1)), by=pip_stratum])

cat("\n=== (4) ANCHOR OVERLAP ===\n")
ap <- fread(file.path(FM,"results/seqfunc/haplotype_saturation/v3/anchors/anchor_pip_support.tsv"))
setnames(ap, c("chromosome_hg19","position_hg19"), c("chromosome","position"))
# (a) is the variant itself an anchor-supporting row?
hit <- merge(t[, .(chromosome, position, study_name, rsid, pip, margin, survives_050, primary_eligible, locus_id)],
             ap[, .(chromosome, position, study_name, anchor_rank, variant_id_hg38, anchor_pip=pip)],
             by=c("chromosome","position","study_name"))
cat("  (a) re-oriented variants that ARE anchor-supporting rows:", nrow(hit), "\n")
if (nrow(hit)) print(hit[order(anchor_rank)])
# (b) re-oriented variants sitting in a locus that contributed ANY anchor
locs <- unique(ap[, .(locus_id, anchor_rank)])
inloc <- merge(t[, .(locus_id, study_name, chromosome, position, rsid, pip, margin, survives_050, primary_eligible)],
               locs, by="locus_id", allow.cartesian=TRUE)
cat("\n  (b) re-oriented variants in an anchor-contributing LOCUS:", uniqueN(inloc[, .(chromosome,position,study_name)]),
    "across", uniqueN(inloc$locus_id), "loci\n")
if (nrow(inloc)) {
  s <- unique(inloc, by=c("chromosome","position","study_name","anchor_rank"))
  print(s[order(anchor_rank, -pip)][1:min(40,.N)])
  fwrite(s, file.path(OUT,"anchor_locus_overlap.tsv"), sep="\t")
}
cat("\nwrote:", OUT, "\n")
