# Corrected: ensembl_verdicts_190.tsv is self-contained (its own margin, pip,
# primary, verdict), so no merge is needed. My previous version keyed on rsid
# alone because the two files name their columns differently, which duplicated
# any variant appearing in several studies.
suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
v <- fread(file.path(FM, "tmp_diag/all190/ensembl_verdicts_190.tsv"))
cat(sprintf("rows %d | distinct study x variant %d\n", nrow(v),
            uniqueN(v[, .(study, chr, pos)])))
print(v[, .N, by = verdict][order(-N)])
v[, is_wrong := grepl("WRONG", verdict, ignore.case = TRUE)]
adj <- v[!is.na(margin) & (is_wrong | grepl("CORRECT", verdict, ignore.case = TRUE))]
cat(sprintf("\nadjudicated with a margin: %d (%d wrong / %d correct)\n\n",
            nrow(adj), sum(adj$is_wrong), sum(!adj$is_wrong)))
cat("=== ANCHOR RISK: max PIP among WRONG calls surviving each threshold ===\n")
cat("   (anchors are the top 500 by PIP, so a wrong call at PIP~0 is unselectable)\n")
for (t in c(0.30,0.40,0.50,0.60,0.70,0.80)) {
  s <- adj[margin >= t]; w <- s[is_wrong %in% TRUE]
  cat(sprintf("  >= %.2f : kept %3d | wrong %2d | purity %5.1f%% | max PIP(wrong) %.4f | wrong PIP>=0.10: %d | wrong primary: %d\n",
      t, nrow(s), nrow(w), 100*mean(!s$is_wrong), if (nrow(w)) max(w$pip, na.rm=TRUE) else 0,
      sum(w$pip >= 0.10, na.rm=TRUE), sum(w$primary %in% c(TRUE,"TRUE"), na.rm=TRUE)))
}
cat("\n=== COST: correct re-orientations lost at each threshold ===\n")
for (t in c(0.30,0.40,0.50,0.60,0.70,0.80)) {
  l <- adj[margin < t & is_wrong %in% FALSE]
  cat(sprintf("  >= %.2f : loses %3d correct | PIP>=0.10: %2d | PIP>=0.90: %2d\n",
      t, nrow(l), sum(l$pip >= 0.10, na.rm=TRUE), sum(l$pip >= 0.90, na.rm=TRUE)))
}
