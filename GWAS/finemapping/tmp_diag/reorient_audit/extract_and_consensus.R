#!/usr/bin/env Rscript
# EXTERNAL VALIDATION of the palindromic variants the resolver RE-ORIENTED
# (match_class complement_*) at high PIP.
#
# Logic: the re-orientation fired because study S's reported AF disagreed with
# the 1kg EUR panel under the letters reading. Independent EUR cohorts break the
# tie WITHOUT using the panel at all:
#   * other EUR studies AGREE with S  -> the PANEL is the outlier -> re-orientation likely WRONG
#   * other EUR studies AGREE with the panel (disagree with S) -> S's AF is genuinely
#     inverted -> re-orientation RIGHT
suppressPackageStartupMessages(library(data.table))
FM  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
OUT <- file.path(FM, "tmp_diag/reorient_audit")
src <- readLines(file.path(FM,"src/06_susie_coloc.R"))
FM_DIR <- FM
i0 <- grep("^PANEL_AF_ROOT <- ", src); cl <- grep("^\\}$", src)
i1 <- min(cl[cl > grep("^strand_verdict_vs_panel <- ", src)])
eval(parse(text=paste(src[i0:i1], collapse="\n")))
cat("helpers sourced; panel_freq_of has", length(formals(panel_freq_of)), "args\n")

v <- fread(file.path(FM,"runs/uniform35_v4_2026-08-04/aggregate/susie_variant_results.tsv.gz"),
  select=c("locus_id","study_name","ancestry","chromosome","position","rsid","effect_allele",
           "other_allele","beta","se","pval","match_class","palindromic","strand_resolution",
           "strand_source","af","panel_af_a1","pip","cs_id","primary_eligible","secondary_eligible"),
  showProgress=FALSE)
tgt <- v[palindromic==TRUE & match_class %in% c("complement_direct","complement_swapped") & pip >= 0.10]
setorder(tgt, -pip)
cat("re-oriented palindromic variants with PIP>=0.10:", nrow(tgt), "\n")
# study AF oriented onto the recorded effect_allele (= ld_allele1)
tgt[, af_on_effect := fifelse(match_class %in% c("direct","complement_swapped"), af, 1-af)]
tgt[, d_same := abs(af_on_effect - panel_af_a1)]
tgt[, d_opp  := abs(af_on_effect - (1-panel_af_a1))]
tgt[, margin := d_same - d_opp]
fwrite(tgt, file.path(OUT,"reoriented_highpip.tsv"), sep="\t")

# ---- independent EUR cohorts ----
afsrc <- fread(file.path(FM,"config/gwas_af_sources.tsv"))
reg   <- fread(file.path(FM,"config/gwas_registry.tsv"))
eur <- reg[ancestry=="EUR", study_name]
cand <- afsrc[study_name %in% eur & af_coverage > 0.5 & file.exists(af_sumstats_path)]
cat("independent EUR cohorts with usable AF:", nrow(cand), "\n")
pos <- unique(tgt[, .(chromosome, position)])
acc <- list()
for (i in seq_len(nrow(cand))) {
  S <- cand$study_name[i]
  g <- tryCatch(fread(cand$af_sumstats_path[i],
        select=c("chromosome","position","allele1","allele2","af"), showProgress=FALSE),
        error=function(e) NULL)
  if (is.null(g)) { cat("  skip", S, "\n"); next }
  h <- merge(pos, g, by=c("chromosome","position"))
  if (nrow(h)) { h[, src := S]; acc[[length(acc)+1L]] <- h }
  cat("  ", S, ": ", nrow(h), " of ", nrow(pos), " target positions\n", sep=""); flush.console()
}
o <- rbindlist(acc)
o[, `:=`(a1=toupper(allele1), a2=toupper(allele2))]
m <- merge(o, tgt[, .(chromosome, position, effect_allele, other_allele, reoriented_study=study_name,
                      reoriented_af_on_effect=af_on_effect, panel_af_a1, pip, margin)],
           by=c("chromosome","position"), allow.cartesian=TRUE)
# orient each cohort's AF onto the recorded effect_allele, by letters
m[, cohort_af_on_effect := fifelse(a1==effect_allele, af, fifelse(a2==effect_allele, 1-af, NA_real_))]
m <- m[!is.na(cohort_af_on_effect) & src != reoriented_study]
fwrite(m, file.path(OUT,"cohort_af_at_reoriented.tsv"), sep="\t")

cons <- m[, .(n_cohorts=.N,
              median_cohort_af=median(cohort_af_on_effect),
              min_cohort_af=min(cohort_af_on_effect), max_cohort_af=max(cohort_af_on_effect)),
          by=.(chromosome, position, effect_allele, other_allele, reoriented_study,
               reoriented_af_on_effect, panel_af_a1, pip, margin)]
cons[, dist_cohorts_to_study := abs(median_cohort_af - reoriented_af_on_effect)]
cons[, dist_cohorts_to_panel := abs(median_cohort_af - panel_af_a1)]
cons[, verdict := fifelse(n_cohorts < 2, "insufficient_cohorts",
                   fifelse(dist_cohorts_to_study < dist_cohorts_to_panel - 0.10, "REORIENT_LIKELY_WRONG (cohorts side with the study)",
                    fifelse(dist_cohorts_to_panel < dist_cohorts_to_study - 0.10, "REORIENT_SUPPORTED (cohorts side with the panel)",
                            "ambiguous")))]
setorder(cons, -pip, verdict)
fwrite(cons, file.path(OUT,"consensus_verdicts.tsv"), sep="\t")
cat("\n=== CONSENSUS VERDICTS ===\n"); print(cons[, .N, by=verdict])
cat("\n=== per variant ===\n")
print(cons[, .(chromosome, position, effect_allele, reoriented_study, study_af=round(reoriented_af_on_effect,4),
               panel_af=round(panel_af_a1,4), cohort_med=round(median_cohort_af,4), n_cohorts,
               margin=round(margin,4), pip=round(pip,3), verdict)], nrows=60)
cat("\nwrote:", OUT, "\n")
