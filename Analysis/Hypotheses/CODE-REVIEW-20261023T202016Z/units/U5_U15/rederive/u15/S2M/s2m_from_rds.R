# S2M from the raw meSuSiE fits: count credible sets and their ancestry category over the
# ALT/AST/NAFLD loci; flag pairs of credible sets within a locus that share variants. No randomness.
suppressPackageStartupMessages(library(data.table))
R <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/mesusie_mvp"
OUT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U5_U15/rederive/u15/S2M"
fs <- unlist(lapply(c("ALT", "AST", "NAFLD"), function(t) list.files(file.path(R, t), "_mesusie\\.rds$", full.names = TRUE)))
cat("rds files:", length(fs), "\n")
rows <- rbindlist(lapply(fs, function(f) {
  x <- readRDS(f)
  cs <- x$credible_sets
  if (!length(cs)) return(data.table(locus_id = x$locus_id, converged = x$converged, cs = NA_integer_, category = NA_character_, size = NA_integer_, overlaps_other = NA))
  ovl <- sapply(seq_along(cs), function(i) any(sapply(seq_along(cs)[-i], function(j) length(intersect(cs[[i]], cs[[j]])) > 0)))
  data.table(locus_id = x$locus_id, converged = x$converged, cs = seq_along(cs), category = as.character(x$cs_category),
             size = lengths(cs), overlaps_other = ovl)
}))
cat("loci:", uniqueN(rows$locus_id), "; converged:", paste(names(table(rows[!duplicated(locus_id), converged])), table(rows[!duplicated(locus_id), converged])), "\n")
cs <- rows[!is.na(cs)]
cat("credible sets:", nrow(cs), "in", uniqueN(cs$locus_id), "loci\n")
cs[, sharing := fifelse(lengths(strsplit(category, "_")) >= 2, "Model-shared (>=2 ancestries)", paste0(category, "-restricted"))]
tab <- cs[, .N, by = sharing][, pct := round(100 * N / sum(N), 1)][order(-N)]
print(tab)
cat("credible sets sharing variants with another set in the same locus:\n"); print(cs[overlaps_other == TRUE])
cat("EAS-only sets:", cs[category == "EAS", .N], "\n")
fwrite(cs, file.path(OUT, "S2M_from_rds_credible_sets.csv")); fwrite(tab, file.path(OUT, "S2M_from_rds.csv"))
