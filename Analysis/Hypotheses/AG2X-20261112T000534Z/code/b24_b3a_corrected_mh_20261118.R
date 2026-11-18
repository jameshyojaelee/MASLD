# Stage B 2026-11-18 B3a orientation-corrected & measured-only MH + pooled Fisher test
suppressPackageStartupMessages(library(data.table))
set.seed(20261111)
args <- commandArgs(trailingOnly = TRUE)
TAB <- if (length(args) >= 1) args[1] else
  "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/out_20261118/stageB"
to01 <- function(v) as.integer(tolower(as.character(v)) %in% c("1", "1.0", "true", "t"))
x_all <- fread(file.path(TAB, "stageB_B3a_corrected_input.tsv"))
x_all[, supported := to01(supported)]
x_all[, outcome := to01(outcome)]
ft_all <- if (nrow(x_all) > 0 && uniqueN(x_all$supported) == 2 && uniqueN(x_all$outcome) == 2) {
  fisher.test(table(factor(x_all$supported, 0:1), factor(x_all$outcome, 0:1)))
} else NULL
x <- x_all[, if (.N >= 2 && uniqueN(supported) == 2) .SD, by = locus_id]
or_exact <- NA_real_; ci_lo <- NA_real_; ci_hi <- NA_real_; p_exact <- NA_real_
if (uniqueN(x$locus_id) >= 1 && uniqueN(x$outcome) == 2) {
  tab <- xtabs(~ factor(supported, 0:1) + factor(outcome, 0:1) + factor(locus_id), data = x)
  if (dim(tab)[3] == 1) {
    mh <- fisher.test(tab[, , 1])
  } else {
    mh <- mantelhaen.test(tab, exact = TRUE)
  }
  or_exact <- unname(mh$estimate)
  ci_lo <- mh$conf.int[1]
  ci_hi <- mh$conf.int[2]
  p_exact <- mh$p.value
}
a <- data.table(
  arm = x_all$arm[1],
  n_measured_total = nrow(x_all),
  supported_with_outcome_total = sum(x_all$supported == 1 & x_all$outcome == 1),
  supported_n_total = sum(x_all$supported == 1),
  other_with_outcome_total = sum(x_all$supported == 0 & x_all$outcome == 1),
  other_n_total = sum(x_all$supported == 0),
  or_fisher_pooled = if (!is.null(ft_all)) unname(ft_all$estimate) else NA_real_,
  ci_lo_fisher_pooled = if (!is.null(ft_all)) ft_all$conf.int[1] else NA_real_,
  ci_hi_fisher_pooled = if (!is.null(ft_all)) ft_all$conf.int[2] else NA_real_,
  p_fisher_pooled = if (!is.null(ft_all)) ft_all$p.value else NA_real_,
  n_members = nrow(x),
  n_loci = uniqueN(x$locus_id),
  supported_with_outcome = sum(x$supported == 1 & x$outcome == 1),
  supported_n = sum(x$supported == 1),
  other_with_outcome = sum(x$supported == 0 & x$outcome == 1),
  other_n = sum(x$supported == 0),
  or_exact = or_exact,
  ci_lo = ci_lo,
  ci_hi = ci_hi,
  p = p_exact
)
orig <- fread(file.path(TAB, "stageB_B3.tsv"))[
  arm != "B3a_caqtl_matching_sign",
  .(arm, n_measured_total, supported_with_outcome_total, supported_n_total,
    other_with_outcome_total, other_n_total, or_fisher_pooled,
    ci_lo_fisher_pooled = NA_real_, ci_hi_fisher_pooled = NA_real_, p_fisher_pooled,
    n_members, n_loci, supported_with_outcome, supported_n,
    other_with_outcome, other_n, or_exact, ci_lo, ci_hi, p)
]
res <- rbind(a, orig)
res[, p_holm := p.adjust(p, "holm")]
res[, decision := fifelse(is.na(p_holm), "not testable", fifelse(p_holm < 0.05, "pass", "fail"))]
res[, `:=`(posthoc_orientation_correction = arm == x_all$arm[1], candidate_not_adopted = TRUE)]
fwrite(res, file.path(TAB, "stageB_B3_corrected_posthoc.tsv"), sep = "\t")
print(res)
