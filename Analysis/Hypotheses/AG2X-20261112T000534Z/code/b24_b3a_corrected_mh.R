# POST HOC (see b23_b3a_corrected.py): exact MH for the orientation-corrected B3a; Holm with the unchanged B3b, B3c.
suppressPackageStartupMessages(library(data.table))
set.seed(20261111)
TAB <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/stageB/tables"
x <- fread(file.path(TAB, "stageB_B3a_corrected_input.tsv"))
x <- x[, if (.N >= 2 && uniqueN(supported) == 2) .SD, by = locus_id]
tab <- xtabs(~ factor(supported, 0:1) + factor(outcome, 0:1) + factor(locus_id), data = x)
mh <- mantelhaen.test(tab, exact = TRUE)
a <- data.table(arm = x$arm[1], n_members = nrow(x), n_loci = uniqueN(x$locus_id),
                supported_with_outcome = sum(x$supported == 1 & x$outcome == 1), supported_n = sum(x$supported == 1),
                other_with_outcome = sum(x$supported == 0 & x$outcome == 1), other_n = sum(x$supported == 0),
                or_exact = unname(mh$estimate), ci_lo = mh$conf.int[1], ci_hi = mh$conf.int[2], p = mh$p.value)
orig <- fread(file.path(TAB, "stageB_B3.tsv"))[arm != "B3a_caqtl_matching_sign",
                                                .(arm, n_members, n_loci, supported_with_outcome, supported_n,
                                                  other_with_outcome, other_n, or_exact, ci_lo, ci_hi, p)]
res <- rbind(a, orig)
res[, p_holm := p.adjust(p, "holm")]
res[, decision := fifelse(p_holm < 0.05, "pass", "fail")]
res[, `:=`(posthoc_orientation_correction = arm == x$arm[1], candidate_not_adopted = TRUE)]
fwrite(res, file.path(TAB, "stageB_B3_corrected_posthoc.tsv"), sep = "\t")
print(res)
