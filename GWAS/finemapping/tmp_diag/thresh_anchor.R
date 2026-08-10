# The threshold decision should be made on ANCHOR RISK, not on overall purity:
# anchors drive expensive saturation experiments, and they are selected top-500
# by max PIP. A WRONG call at PIP ~ 0 is eligible-but-unselectable and cannot
# become an anchor; a WRONG call at high PIP can.
suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
D  <- file.path(FM, "tmp_diag/all190")
a <- fread(file.path(D, "all_reoriented.tsv"))
v <- fread(file.path(D, "ensembl_verdicts_190.tsv"))
cat("all_reoriented cols :", paste(names(a), collapse=", "), "\n")
cat("verdict cols        :", paste(names(v), collapse=", "), "\n\n")
kc <- intersect(names(a), names(v)); kc <- kc[kc %in% c("rsid","chromosome","position","study_name")]
m <- merge(a, v, by = kc, all.x = TRUE, suffixes = c("", ".v"))
wc <- grep("verdict|wrong|correct", names(m), value = TRUE, ignore.case = TRUE)
cat("verdict-ish columns:", paste(wc, collapse=", "), "\n")
mc <- grep("^margin$|af_margin", names(m), value = TRUE)[1]
pc <- grep("^pip$", names(m), value = TRUE)[1]
cat(sprintf("using margin col '%s', pip col '%s'\n\n", mc, pc))
wcol <- wc[1]
m[, is_wrong := grepl("WRONG", get(wcol), ignore.case = TRUE)]
cat("=== THE DECIDING NUMBER: max PIP among WRONG calls that SURVIVE each threshold ===\n")
for (t in c(0.30,0.40,0.50,0.60,0.70,0.80)) {
  s <- m[is.finite(get(mc)) & get(mc) >= t]
  w <- s[is_wrong %in% TRUE]
  cat(sprintf("  margin >= %.2f : kept %3d | wrong %3d | max PIP among wrong %.4f | wrong with PIP>=0.10: %d\n",
      t, nrow(s), nrow(w), if (nrow(w)) max(w[[pc]], na.rm=TRUE) else 0, sum(w[[pc]] >= 0.10, na.rm=TRUE)))
}
cat("\n=== correct re-orientations LOST at each threshold (the cost side) ===\n")
for (t in c(0.30,0.40,0.50,0.60,0.70,0.80)) {
  lost <- m[is.finite(get(mc)) & get(mc) < t & is_wrong %in% FALSE]
  cat(sprintf("  margin >= %.2f : loses %3d correct | of which PIP>=0.10: %d | PIP>=0.90: %d\n",
      t, nrow(lost), sum(lost[[pc]] >= 0.10, na.rm=TRUE), sum(lost[[pc]] >= 0.90, na.rm=TRUE)))
}
