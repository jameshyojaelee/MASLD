# U8a: compare regenerated nulls (02_null.R) with the stored H59 null and across seeds.
# Run after jobs 22383479 (seed 20260923), 22383480 (seed 1), 22383481 (seed 20261023) finish.
here <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U8a/rederive"
source(file.path(here, "u8a_lib.R"))
seeds <- c(20260923L, 1L, 20261023L)
pr <- fread(file.path(H59, "out/outcome/results_primary_family.tsv"))
cr <- fread(file.path(H59, "out/null_critical_abs_t.tsv"))

# same-seed null matrix vs stored
A <- as.matrix(fread(file.path(here, "out/null_seed_20260923/null_gls_t.tsv.gz")))
S <- as.matrix(fread(file.path(H59, "out/null_gls_t.tsv.gz")))[, colnames(A)]
say("same seed GLS null: %d x %d; max |t_own - t_stored| = %.3g; draws differing > 1e-4: %d",
    nrow(A), ncol(A), max(abs(A - S)), sum(abs(A - S) > 1e-4))
A62 <- as.matrix(fread(file.path(here, "out/null_seed_20260923/null_ols62_t.tsv.gz")))
S62 <- as.matrix(fread(file.path(H59, "out/null_ols62_t.tsv.gz")))[, colnames(A62)]
say("same seed OLS62 null: max |diff| = %.3g", max(abs(A62 - S62)))

res <- rbindlist(lapply(seeds, function(s) {
  x <- fread(file.path(here, sprintf("out/null_seed_%d/perm_results_117.tsv", s)))
  x[, seed := s][, rank := frank(perm_p, ties.method = "min")][]
}))
fwrite(res, file.path(here, "out/seed_comparison_117.tsv"), sep = "\t")
w <- dcast(res, program + module_name ~ seed, value.var = c("perm_p", "perm_q", "selected", "rank"))
say("selected per seed: %s", paste(res[, sum(selected), by = seed][, paste(seed, V1)], collapse = "; "))
say("q < 0.05 per seed: %s", paste(res[, sum(perm_q < 0.05), by = seed][, paste(seed, V1)], collapse = "; "))
say("smallest q per seed: %s", paste(res[, sprintf("%.6f", min(perm_q)), by = seed][, paste(seed, V1)], collapse = "; "))
top <- res[rank <= 5][order(seed, perm_p), .(seed, program, module_name, perm_p, perm_q, rank)]
print(top)
say("top-5 sets identical across seeds: %s",
    identical(sort(res[seed == seeds[1] & rank <= 5, program]), sort(res[seed == seeds[2] & rank <= 5, program])) &&
    identical(sort(res[seed == seeds[1] & rank <= 5, program]), sort(res[seed == seeds[3] & rank <= 5, program])))
named <- res[program %in% c("hepatocytes_20", "hepatocytes_8"), .(seed, program, perm_p, perm_q, perm_q62, perm_q60, crit)]
print(named)
same <- merge(res[seed == 20260923L, .(program, perm_p, perm_q, crit, perm_q62, maxT_p)],
              pr[, .(program, p_p = perm_p, p_q = perm_q, p_q62 = ols62_perm_q, p_maxT = maxT_fwer_p)], by = "program")
same <- merge(same, cr[, .(program, p_crit = null_abs_t_crit)], by = "program")
say("same seed vs producer: max|dp| %.3g max|dq| %.3g max|dcrit| %.3g max|dq62| %.3g max|dmaxT| %.3g",
    max(abs(same$perm_p - same$p_p)), max(abs(same$perm_q - same$p_q)), max(abs(same$crit - same$p_crit)),
    max(abs(same$perm_q62 - same$p_q62)), max(abs(same$maxT_p - same$p_maxT)))
for (s in seeds) { f <- file.path(here, sprintf("out/null_seed_%d/pc4.tsv", s)); say("PC4 seed %d:", s); print(fread(f)) }
