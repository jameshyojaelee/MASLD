# EXTERNAL adjudication of the re-oriented variants.
# The resolver saw study_af vs panel_af and concluded "opposite", flipping beta.
# gnomAD is independent of BOTH. For the study's own effect allele:
#   if gnomAD agrees with the STUDY  -> the study was fine; the FLIP WAS WRONG.
#   if gnomAD agrees with the PANEL  -> the study was flipped; the flip was RIGHT.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
h <- unique(fread(file.path(FM,"tmp_diag/reoriented_highpip.tsv")),
            by=c("study_name","chromosome","position"))
POP <- c(EUR="gnomADg:nfe", AFR="gnomADg:afr", AMR="gnomADg:amr",
         EAS="gnomADg:eas", SAS="gnomADg:sas")
anc_of <- function(s) sub(".*_", "", s)
get_freq <- function(rs, allele, pop) {
  f <- file.path(FM,"tmp_diag/gnomad", paste0(rs,".json"))
  if (!file.exists(f)) return(NA_real_)
  j <- tryCatch(fromJSON(f), error=function(e) NULL); if (is.null(j$populations)) return(NA_real_)
  p <- as.data.table(j$populations)
  al <- allele; x <- p[population == pop & allele == al]
  if (!nrow(x)) return(NA_real_)
  as.numeric(x$frequency[1])
}
out <- rbindlist(lapply(seq_len(nrow(h)), function(i) {
  r <- h[i]; anc <- r$ancestry; pop <- POP[[anc]]
  if (is.null(pop) || is.na(anc)) return(NULL)
  g_ea <- get_freq(r$rsid, r$effect_allele, pop)
  data.table(study=r$study_name, rsid=r$rsid, anc=anc, pip=round(r$pip,3),
             elig=r$primary_eligible, ea=r$effect_allele, oa=r$other_allele,
             study_af=round(r$af,4), panel_af=round(r$panel_af_a1,4),
             gnomad_ea=round(g_ea,4))
}), fill=TRUE)
out[, d_study := abs(study_af - gnomad_ea)]
out[, d_panel := abs(panel_af - gnomad_ea)]
out[, verdict := fifelse(is.na(gnomad_ea), "no gnomAD data",
                  fifelse(d_study < d_panel, "*** FLIP WAS WRONG (gnomAD backs the STUDY) ***",
                                             "flip was RIGHT (gnomAD backs the panel)"))]
setorder(out, -elig, -pip)
print(out[, .(study, rsid, anc, pip, elig, ea, study_af, panel_af, gnomad_ea,
              d_study=round(d_study,3), d_panel=round(d_panel,3), verdict)], nrows=40)
cat("\n=== SUMMARY ===\n"); print(out[, .N, by=verdict][order(-N)])
cat("\n--- the 7 primary_eligible (anchor-pool feeders) ---\n")
print(out[elig==TRUE, .(study, rsid, pip, study_af, panel_af, gnomad_ea, verdict)])
fwrite(out, file.path(FM,"tmp_diag/gnomad_verdict.tsv"), sep="\t")
