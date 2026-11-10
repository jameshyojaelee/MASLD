# U8a: recompute permutation p, BH q (n = 117) and critical |t| from the stored H59 null
# with the independently refitted observed t (out/observed/observed_gls_117.tsv).
here <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U8a/rederive"
source(file.path(here, "u8a_lib.R"))
out <- file.path(here, "out/stored_null"); dir.create(out, showWarnings = FALSE)
obs <- fread(file.path(here, "out/observed/observed_gls_117.tsv"))
TG <- as.matrix(fread(file.path(H59, "out/null_gls_t.tsv.gz")))
say("stored null: %d draws x %d programs; non-finite %d", nrow(TG), ncol(TG), sum(!is.finite(TG)))
stopifnot(setequal(colnames(TG), obs$program))
obs[, B_ok := sapply(program, function(p) sum(is.finite(TG[, p])))]
obs[, perm_p := sapply(seq_len(.N), function(i) perm_p(t[i], TG[, program[i]]))]
obs[, perm_q := p.adjust(perm_p, "BH", n = 117)]
obs[, crit := sapply(program, function(p) quantile(abs(TG[, p]), 1 - 0.05 / 117, type = 7, na.rm = TRUE))]
obs[, selected := perm_q < 0.05 & stability_median >= 0.5]
mx <- apply(abs(TG), 1, max)
obs[, maxT_p := sapply(t, function(x) (1 + sum(mx >= abs(x))) / (1 + length(mx)))]

pr <- fread(file.path(H59, "out/outcome/results_primary_family.tsv"))
cr <- fread(file.path(H59, "out/null_critical_abs_t.tsv"))
cmp <- merge(obs[, .(program, module_name, t, perm_p, perm_q, crit, maxT_p, selected)],
             pr[, .(program, p_perm_p = perm_p, p_perm_q = perm_q, p_maxT = maxT_fwer_p, p_sel = selected_H59)], by = "program")
cmp <- merge(cmp, cr[, .(program, p_crit = null_abs_t_crit)], by = "program")
say("vs producer: max|dp| %.3g, max|dq| %.3g, max|dcrit| %.3g, max|dmaxT| %.3g, selection mismatches %d",
    max(abs(cmp$perm_p - cmp$p_perm_p)), max(abs(cmp$perm_q - cmp$p_perm_q)), max(abs(cmp$crit - cmp$p_crit)),
    max(abs(cmp$maxT_p - cmp$p_maxT)), sum(cmp$selected != cmp$p_sel))
setorder(cmp, perm_p)
print(cmp[1:8])
print(cmp[program %in% c("hepatocytes_20", "hepatocytes_8")])
say("selected %d; smallest q %.6f; smallest maxT p %.5f", sum(cmp$selected), min(cmp$perm_q), min(cmp$maxT_p))
# ties near observed |t|: does the (1 - 1e-12) tolerance matter?
obs[, p_strict := sapply(seq_len(.N), function(i) { x <- TG[, program[i]]; (1 + sum(abs(x) >= abs(t[i]) * (1 - 1e-6))) / (1 + length(x)) })]
say("max |p(tol 1e-6) - p| = %.3g", max(abs(obs$p_strict - obs$perm_p)))
fwrite(cmp, file.path(out, "stored_null_recomputed.tsv"), sep = "\t")
