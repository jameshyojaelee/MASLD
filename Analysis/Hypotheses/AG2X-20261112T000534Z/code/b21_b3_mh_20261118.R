# Stage B 2026-11-18 B3 exact Mantel-Haenszel + pooled Fisher summary on out_20261118/stageB
suppressPackageStartupMessages(library(data.table))
set.seed(20261111)
args <- commandArgs(trailingOnly = TRUE)
TAB <- if (length(args) >= 1) args[1] else
  "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/out_20261118/stageB"
d <- fread(file.path(TAB, "stageB_B3_input.tsv"))
res <- rbindlist(lapply(split(d, d$arm), function(x) {
  x_all <- copy(x)
  ft_all <- if (uniqueN(x_all$supported) == 2 && uniqueN(x_all$outcome) == 2) {
    fisher.test(table(factor(x_all$supported, 0:1), factor(x_all$outcome, 0:1)))
  } else NULL
  x <- x[, if (.N >= 2 && uniqueN(supported) == 2) .SD, by = locus_id]
  out <- data.table(
    arm = x_all$arm[1],
    n_measured_total = nrow(x_all),
    supported_with_outcome_total = sum(x_all$supported == 1 & x_all$outcome == 1),
    supported_n_total = sum(x_all$supported == 1),
    other_with_outcome_total = sum(x_all$supported == 0 & x_all$outcome == 1),
    other_n_total = sum(x_all$supported == 0),
    or_fisher_pooled = if (!is.null(ft_all)) unname(ft_all$estimate) else NA_real_,
    p_fisher_pooled = if (!is.null(ft_all)) ft_all$p.value else NA_real_,
    n_members = nrow(x),
    n_loci = uniqueN(x$locus_id),
    supported_with_outcome = sum(x$supported == 1 & x$outcome == 1),
    supported_n = sum(x$supported == 1),
    other_with_outcome = sum(x$supported == 0 & x$outcome == 1),
    other_n = sum(x$supported == 0),
    or_exact = NA_real_,
    ci_lo = NA_real_,
    ci_hi = NA_real_,
    p = NA_real_
  )
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
