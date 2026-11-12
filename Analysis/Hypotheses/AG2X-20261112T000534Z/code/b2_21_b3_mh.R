# Stage B2, B3 (PRESPEC_AG2X_stageB2.md section 4): exact conditional MH by locus, Holm over 3. Nothing adopted.
# sequence-supported vs other members of the same credible sets. Holm over the 3 arms. Nothing is adopted.
suppressPackageStartupMessages(library(data.table))
set.seed(20261111)
args <- commandArgs(trailingOnly = TRUE)
TAB <- if (length(args) >= 1) args[1] else
  "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/stageB2/tables"
ARMS <- c("B3a_caqtl_matching_sign", "B3b_hepatocyte_allelic_imbalance", "B3c_mpra_allelic_activity")
d <- fread(file.path(TAB, "stageB2_B3_input.tsv"))
res <- rbindlist(lapply(ARMS, function(arm_name) {
  x <- d[arm == arm_name]
  # strata with variation in both supported and outcome carry information; keep strata with >= 2 rows
  x <- x[, if (.N >= 2 && uniqueN(supported) == 2) .SD, by = locus_id]
  out <- data.table(arm = arm_name, n_members = nrow(x), n_loci = uniqueN(x$locus_id),
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
res[, p_holm := p.adjust(p, "holm", n = 3)]  # Holm over the 3 prespecified arms, also if one is not testable
res[, decision := fifelse(is.na(p_holm), "not testable", fifelse(p_holm < 0.05, "pass", "fail"))]
res[, candidate_not_adopted := TRUE]
fwrite(res, file.path(TAB, "stageB2_B3.tsv"), sep = "\t")
print(res)
writeLines(capture.output(sessionInfo()), file.path(TAB, "sessionInfo_b2_b3.txt"))
