# Re-oriented palindromes carrying high PIP in the FROZEN v4 release.
# These are variants whose EFFECT DIRECTION was inverted by the resolver.
# Reads RECORDED fields only (match_class / strand_source / pip) -- no verdict
# is recomputed, so the circularity trap in match_class grouping does not apply.
suppressPackageStartupMessages(library(data.table))
V4 <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/runs/uniform35_v4_2026-08-04"
v <- fread(file.path(V4, "aggregate", "susie_variant_results.tsv.gz"),
           select = c("study_name","trait","tier","ancestry","chromosome","position","rsid",
                      "effect_allele","other_allele","beta","pip","cs_id","palindromic",
                      "match_class","strand_resolution","strand_source","af","panel_af_a1",
                      "primary_eligible","locus_id"))
cat(sprintf("v4 variant rows: %s\n\n", format(nrow(v), big.mark=",")))
cat("=== strand_source distribution among PALINDROMES ===\n")
print(v[palindromic == TRUE, .N, by = strand_source][order(-N)])
cat("\n=== match_class among palindromes (complement_* == RE-ORIENTED) ===\n")
print(v[palindromic == TRUE, .N, by = match_class][order(-N)])
ro <- v[palindromic == TRUE & grepl("^complement_", match_class)]
cat(sprintf("\nRE-ORIENTED palindromes total: %s\n", format(nrow(ro), big.mark=",")))
cat(sprintf("  of which strand_source == 'af_over_certificate': %s\n",
            format(sum(ro$strand_source == "af_over_certificate"), big.mark=",")))
for (thr in c(0.10, 0.50, 0.90, 0.999)) {
  h <- ro[pip >= thr]
  cat(sprintf("  PIP >= %-5.3f : %4d variants | %3d loci | %2d studies | primary_eligible %d\n",
      thr, nrow(h), uniqueN(h$locus_id), uniqueN(h$study_name), sum(h$primary_eligible %in% TRUE)))
}
cat("\n=== the high-PIP re-oriented set (PIP >= 0.10), most confident first ===\n")
h <- ro[pip >= 0.10][order(-pip)]
h[, af_delta_same := abs(af - panel_af_a1)]
h[, af_delta_opp  := abs(af - (1 - panel_af_a1))]
print(h[, .(study_name, rsid, chromosome, position, ea = effect_allele, oa = other_allele,
            beta = round(beta,4), pip = round(pip,3), src = strand_source,
            af = round(af,4), panel = round(panel_af_a1,4),
            d_same = round(af_delta_same,3), d_opp = round(af_delta_opp,3),
            elig = primary_eligible)], nrows = 40)
cat("\n=== decision margin |d_same - d_opp| for these calls (bigger = more confident) ===\n")
h[, margin := abs(af_delta_same - af_delta_opp)]
print(h[, .(n=.N, min_margin=round(min(margin),3), median_margin=round(median(margin),3),
            max_margin=round(max(margin),3), n_below_0.2=sum(margin<0.2))])
fwrite(h, "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/tmp_diag/reoriented_highpip.tsv", sep="\t")
