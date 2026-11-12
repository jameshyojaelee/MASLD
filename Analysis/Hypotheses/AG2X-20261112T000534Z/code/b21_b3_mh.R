# Stage B, B3 (PRESPEC_AG2X_stageB.md section 5): exact conditional Mantel-Haenszel odds ratio, strata = locus,
# sequence-supported vs other members of the same credible sets. Holm over the 3 arms. Nothing is adopted.
suppressPackageStartupMessages(library(data.table))
set.seed(20261111)
TAB <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/stageB/tables"
d <- fread(file.path(TAB, "stageB_B3_input.tsv"))
res <- rbindlist(lapply(split(d, d$arm), function(x) {
  # strata with variation in both supported and outcome carry information; keep strata with >= 2 rows
  x <- x[, if (.N >= 2 && uniqueN(supported) == 2) .SD, by = locus_id]
  out <- data.table(arm = x$arm[1], n_members = nrow(x), n_loci = uniqueN(x$locus_id),
                    supported_with_outcome = sum(x$supported == 1 & x$outcome == 1), supported_n = sum(x$supported == 1),
                    other_with_outcome = sum(x$supported == 0 & x$outcome == 1), other_n = sum(x$supported == 0),
                    or_exact = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, p = NA_real_)
  if (uniqueN(x$locus_id) >= 1 && uniqueN(x$outcome) == 2) {
    tab <- xtabs(~ factor(supported, 0:1) + factor(outcome, 0:1) + factor(locus_id), data = x)
    if (dim(tab)[3] == 1) {
      ft <- fisher.test(tab[, , 1])
      out[, `:=`(or_exact = unname(ft$estimate), ci_lo = ft$conf.int[1], ci_hi = ft$conf.int[2], p = ft$p.value)]
    } else {
      mh <- mantelhaen.test(tab, exact = TRUE)
      out[, `:=`(or_exact = unname(mh$estimate), ci_lo = mh$conf.int[1], ci_hi = mh$conf.int[2], p = mh$p.value)]
    }
  }
  out
}))
res[, p_holm := p.adjust(p, "holm")]
res[, decision := fifelse(is.na(p_holm), "not testable", fifelse(p_holm < 0.05, "pass", "fail"))]
res[, candidate_not_adopted := TRUE]
fwrite(res, file.path(TAB, "stageB_B3.tsv"), sep = "\t")
print(res)
writeLines(capture.output(sessionInfo()), file.path(TAB, "sessionInfo_b3.txt"))
