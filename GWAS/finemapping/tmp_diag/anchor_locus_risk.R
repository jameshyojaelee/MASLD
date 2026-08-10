# Can a WRONG-signed re-oriented variant distort the credible set of a locus
# that CONTRIBUTED an anchor, even where the variant is not itself the anchor?
# Join the anchor-contributing loci to the external verdicts.
suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
o <- fread(file.path(FM, "tmp_diag/all190/anchor_locus_overlap.tsv"))
v <- fread(file.path(FM, "tmp_diag/all190/ensembl_verdicts_190.tsv"))
v[, is_wrong := grepl("WRONG", verdict, ignore.case = TRUE)]
# key on study + position: rsid alone duplicates across studies
m <- merge(o, v[, .(study, pos, verdict, is_wrong)],
           by.x = c("study_name","position"), by.y = c("study","pos"), all.x = TRUE)
cat(sprintf("re-oriented variants inside anchor-contributing loci: %d rows, %d distinct loci, %d anchors\n\n",
            nrow(m), uniqueN(m$locus_id), uniqueN(m$anchor_rank)))
cat("=== verdict breakdown ===\n"); print(m[, .N, by = verdict][order(-N)])
w <- m[is_wrong %in% TRUE]
cat(sprintf("\nWRONG re-orientations inside anchor loci: %d\n", nrow(w)))
if (nrow(w)) {
  cat(sprintf("  their PIPs: max %.4f | median %.4f | n with PIP>=0.10: %d | n with PIP>=0.50: %d\n",
      max(w$pip), median(w$pip), sum(w$pip >= 0.10), sum(w$pip >= 0.50)))
  cat(sprintf("  how many SURVIVE the 0.50 fix: %d of %d\n",
      sum(w$survives_050 %in% c(TRUE,"TRUE")), nrow(w)))
  cat("\n  the wrong ones with PIP >= 0.01 (potential credible-set influence):\n")
  print(w[pip >= 0.01][order(-pip), .(anchor_rank, locus_id, study_name, rsid,
        pip = round(pip,4), margin = round(margin,3), survives_050)])
}
cat("\n=== AFTER applying the 0.50 fix, what wrong variants remain in anchor loci? ===\n")
rem <- w[survives_050 %in% c(TRUE,"TRUE")]
if (!nrow(rem)) cat("  NONE — the fix removes every wrong re-orientation from every anchor-contributing locus\n") else
  print(rem[order(-pip), .(anchor_rank, study_name, rsid, pip = round(pip,4), margin = round(margin,3))])
