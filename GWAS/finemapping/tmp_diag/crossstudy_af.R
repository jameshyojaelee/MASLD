# Are the 28 re-oriented AF calls right? Independent on-disk check.
# For each position, pull `af` from EVERY study in the release that carries it.
# The AF call said "the study frequency matches 1-panel, not panel". If that is
# correct, OTHER independent studies at the same position should agree with the
# study, not the panel. If instead every study agrees with the panel and only
# these disagree, the study is flipped. If ALL studies disagree with the panel,
# the PANEL record is the suspect (e.g. a different variant at that position).
suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
V4 <- file.path(FM, "runs/uniform35_v4_2026-08-04")
h  <- fread(file.path(FM, "tmp_diag/reoriented_highpip.tsv"))
v  <- fread(file.path(V4, "aggregate/susie_variant_results.tsv.gz"),
            select = c("study_name","ancestry","chromosome","position","rsid","effect_allele",
                       "other_allele","match_class","strand_resolution","strand_source",
                       "af","panel_af_a1","pip","palindromic","primary_eligible"),
            showProgress = FALSE)
setkey(v, chromosome, position)
# Orient every study's af onto the SAME reference axis as panel_af_a1, using the
# relabel-invariant class grouping: {direct, complement_swapped} -> gwas a1 sits
# on ld_allele1; {swapped, complement_direct} -> it sits on ld_allele2.
v[, af_on_a1 := fifelse(match_class %in% c("direct","complement_swapped"), af, 1 - af)]
cat(sprintf("=== cross-study AF at the %d re-oriented high-PIP positions ===\n\n", nrow(h)))
res <- rbindlist(lapply(seq_len(nrow(h)), function(i) {
  r <- h[i]
  o <- v[.(r$chromosome, r$position)][!is.na(af_on_a1)]
  if (!nrow(o)) return(NULL)
  pan <- r$panel_af_a1
  # distance of each study's oriented af to the panel, and to its complement
  o[, d_pan := abs(af_on_a1 - pan)]
  o[, d_flip := abs(af_on_a1 - (1 - pan))]
  data.table(rsid = r$rsid, chr = r$chromosome, pos = r$position,
             pip = r$pip, elig = r$primary_eligible, panel = round(pan, 4),
             n_studies = nrow(o),
             n_agree_panel = sum(o$d_pan < o$d_flip),
             n_agree_flip  = sum(o$d_flip < o$d_pan),
             med_af_on_a1 = round(median(o$af_on_a1), 4),
             spread = round(diff(range(o$af_on_a1)), 4),
             n_reoriented = sum(grepl("^complement_", o$match_class)))
}), fill = TRUE)
setorder(res, -pip, rsid)
print(res, nrows = 40)
cat("\n=== interpretation ===\n")
cat(sprintf("positions where ALL studies agree with the FLIP (panel record suspect): %d / %d\n",
            sum(res$n_agree_panel == 0), nrow(res)))
cat(sprintf("positions where studies SPLIT (some agree panel, some flip)          : %d / %d\n",
            sum(res$n_agree_panel > 0 & res$n_agree_flip > 0), nrow(res)))
cat(sprintf("positions where ALL studies agree with the PANEL (call would be wrong): %d / %d\n",
            sum(res$n_agree_flip == 0), nrow(res)))
cat(sprintf("\nmedian across-study spread in oriented af: %.4f (small => studies concordant)\n",
            median(res$spread, na.rm = TRUE)))
fwrite(res, file.path(FM, "tmp_diag/crossstudy_af.tsv"), sep = "\t")
