# Stage B2 2026-11-18 B3 exact conditional MH + pooled Fisher by locus, Holm over 3.
suppressPackageStartupMessages(library(data.table))
set.seed(20261111)
args <- commandArgs(trailingOnly = TRUE)
TAB <- if (length(args) >= 1) args[1] else
  "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/out_20261118/stageB2"
ARMS <- c("B3a_caqtl_matching_sign", "B3b_hepatocyte_allelic_imbalance", "B3c_mpra_allelic_activity")
d <- fread(file.path(TAB, "stageB2_B3_input.tsv"))
res <- rbindlist(lapply(ARMS, function(arm_name) {
  x_all <- d[arm == arm_name]
  ft_all <- if (nrow(x_all) > 0 && uniqueN(x_all$supported) == 2 && uniqueN(x_all$outcome) == 2) {
    fisher.test(table(factor(x_all$supported, 0:1), factor(x_all$outcome, 0:1)))
  } else NULL
  x <- x_all[, if (.N >= 2 && uniqueN(supported) == 2) .SD, by = locus_id]
  out <- data.table(
    arm = arm_name,
    n_measured_total = nrow(x_all),
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
res[, p_holm := p.adjust(p, "holm", n = 3)]
res[, decision := fifelse(is.na(p_holm), "not testable", fifelse(p_holm < 0.05, "pass", "fail"))]
res[, candidate_not_adopted := TRUE]
fwrite(res, file.path(TAB, "stageB2_B3.tsv"), sep = "\t")
print(res)
writeLines(capture.output(sessionInfo()), file.path(TAB, "sessionInfo_b2_b3.txt"))
