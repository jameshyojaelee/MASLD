# Corrected: compare ONLY within the same ancestry (same LD panel => same
# panel_af_a1 reference), and record whether the agreeing studies are
# INDEPENDENT of the study that was re-oriented. Nearly all re-orientations are
# in MVP strata, which share a cohort and a processing pipeline, so agreement
# among MVP studies is NOT independent corroboration.
suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
V4 <- file.path(FM, "runs/uniform35_v4_2026-08-04")
h  <- unique(fread(file.path(FM, "tmp_diag/reoriented_highpip.tsv")),
             by = c("study_name","chromosome","position"))
v  <- fread(file.path(V4, "aggregate/susie_variant_results.tsv.gz"),
            select = c("study_name","ancestry","chromosome","position","rsid","match_class",
                       "strand_resolution","strand_source","af","panel_af_a1","pip",
                       "palindromic","primary_eligible"), showProgress = FALSE)
v[, af_on_a1 := fifelse(match_class %in% c("direct","complement_swapped"), af, 1 - af)]
v[, cohort := fifelse(grepl("^MVP_", study_name), "MVP",
               fifelse(grepl("^FinnGen", study_name), "FinnGen",
                fifelse(grepl("^PanUKBB|^UKBB", study_name), "UKBB",
                 fifelse(grepl("^BBJ", study_name), "BBJ", "other"))))]
setkey(v, chromosome, position, ancestry)
out <- rbindlist(lapply(seq_len(nrow(h)), function(i) {
  r <- h[i]
  me <- v[.(r$chromosome, r$position)][study_name == r$study_name][1]
  if (is.na(me$ancestry)) return(NULL)
  o <- v[.(r$chromosome, r$position)][ancestry == me$ancestry & !is.na(af_on_a1)]
  pan <- me$panel_af_a1
  o[, agrees_flip := abs(af_on_a1 - (1-pan)) < abs(af_on_a1 - pan)]
  ind <- o[cohort != me$cohort]          # independent = different cohort family
  data.table(study = r$study_name, rsid = r$rsid, anc = me$ancestry, pip = round(r$pip,3),
             elig = r$primary_eligible, panel = round(pan,4), my_af = round(me$af_on_a1,4),
             n_same_anc = nrow(o), n_indep = nrow(ind),
             indep_agree_flip = sum(ind$agrees_flip), indep_agree_panel = sum(!ind$agrees_flip),
             indep_cohorts = paste(sort(unique(ind$cohort)), collapse=","))
}), fill = TRUE)
setorder(out, -elig, -pip)
print(out, nrows = 40)
cat("\n=== VERDICT PER VARIANT ===\n")
out[, verdict := fifelse(n_indep == 0, "NO INDEPENDENT DATA",
                  fifelse(indep_agree_panel == 0, "independent studies BACK the flip",
                   fifelse(indep_agree_flip == 0, "*** independent studies CONTRADICT the flip ***",
                    "independent studies SPLIT")))]
print(out[, .N, by = verdict][order(-N)])
cat("\n--- restricted to the 7 primary_eligible (anchor-pool feeders) ---\n")
print(out[elig == TRUE, .(study, rsid, pip, n_indep, indep_agree_flip, indep_agree_panel, verdict)])
fwrite(out, file.path(FM, "tmp_diag/crossstudy_af2.tsv"), sep = "\t")
