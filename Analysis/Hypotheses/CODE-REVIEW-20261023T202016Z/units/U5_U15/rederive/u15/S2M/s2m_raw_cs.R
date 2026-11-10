# Inspect the raw meSuSiE credible sets at the two loci where the per-variant table holds one CS
# fewer than the locus summary. No randomness.
R <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/mesusie_mvp"
for (f in c("ALT/locus_MVP_ALT_chr14_94838142_mesusie.rds", "AST/locus_MVP_AST_chr10_101157438_mesusie.rds")) {
  x <- readRDS(file.path(R, f))
  cat("==", f, " n_cs:", x$n_cs, "\n")
  cs <- x$credible_sets
  for (i in seq_along(cs)) {
    ids <- x$snp_ids[cs[[i]]]
    cat(sprintf("CS %d (%s) category %s purity %s size %d: %s\n", i, names(cs)[i], x$cs_category[i],
                if (is.data.frame(x$cs_purity)) paste(round(unlist(x$cs_purity[i, , drop = TRUE]), 3), collapse = "/") else "NA", length(ids), paste(head(ids, 6), collapse = " ")))
  }
  ov <- outer(seq_along(cs), seq_along(cs), Vectorize(function(a, b) length(intersect(cs[[a]], cs[[b]]))))
  cat("pairs of CS sharing variants:\n"); print(which(ov > 0 & upper.tri(ov), arr.ind = TRUE))
}
