#!/usr/bin/env Rscript
# 03b: is the axis instability real, or an artefact of how I summarised branches?
#
# The concern. When Slingshot returns a different NUMBER of lineages under
# different settings, the curve-weighted average pseudotime is a different object
# each time -- an average over 2 branches is not comparable to an average over 5.
# So part of the observed disagreement between settings could be my summarisation
# choice rather than instability of the method.
#
# The test. Stratify the setting-pairs by whether the two settings returned the
# SAME number of lineages. If agreement is high within a stratum and low across
# strata, the instability is really branch-count selection, and the honest claim
# becomes "the method does not identify how many trajectories exist" rather than
# "the ordering is unstable". If agreement is low even within a stratum, the
# ordering is unstable in its own right and my summary is not to blame.
#
# Either way the underlying finding survives, but they are different claims and
# the paper must make the right one.

suppressPackageStartupMessages({ library(data.table) })
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")

rows <- list()
for (arm in c("A1", "A2")) {
  sp <- file.path(OUT, "arms", sprintf("%s_axis_stability.tsv", arm))
  gp <- file.path(OUT, "arms", sprintf("%s_grid_settings.tsv", arm))
  if (!file.exists(sp)) { log_step("skipping ", arm, ": no stability matrix"); next }

  s <- fread(sp); M <- as.matrix(s[, -1]); rownames(M) <- s$setting
  key <- data.table(setting = rownames(M))
  if (file.exists(gp)) {
    g <- fread(gp)
    g[, setting := sprintf("k%d_npc%d_seed%d", k, npc, seed)]
    key <- merge(key, g[, .(setting, n_lineages, k, npc, seed)], by = "setting", all.x = TRUE)
  } else {
    key[, `:=`(n_lineages = NA_integer_, k = NA_integer_, npc = NA_integer_, seed = NA_integer_)]
  }
  setkey(key, setting)

  idx <- which(upper.tri(M), arr.ind = TRUE)
  pairs <- data.table(
    arm = arm,
    a = rownames(M)[idx[, 1]], b = rownames(M)[idx[, 2]],
    rho = M[upper.tri(M)])
  pairs[, lin_a := key[.(a), n_lineages]]
  pairs[, lin_b := key[.(b), n_lineages]]
  pairs[, k_a := key[.(a), k]][, k_b := key[.(b), k]]
  pairs[, same_lineage_count := lin_a == lin_b]
  pairs[, same_k := k_a == k_b]
  rows[[arm]] <- pairs
}
assert_true(length(rows) > 0L, "No stability matrices found for A1 or A2")
allp <- rbindlist(rows)
write_tsv_once(allp, file.path(OUT, "arms", "lineage_stratified_pairs.tsv"))

summ <- allp[, .(n_pairs = .N,
                 median_abs_rho = median(abs(rho), na.rm = TRUE),
                 frac_below_0.3 = mean(abs(rho) < 0.3, na.rm = TRUE),
                 frac_above_0.8 = mean(abs(rho) >= 0.8, na.rm = TRUE)),
             by = .(arm, same_lineage_count)][order(arm, -same_lineage_count)]
byk <- allp[, .(n_pairs = .N,
                median_abs_rho = median(abs(rho), na.rm = TRUE)),
            by = .(arm, same_k)][order(arm, -same_k)]
write_tsv_once(summ, file.path(OUT, "arms", "lineage_stratified_summary.tsv"))
write_tsv_once(byk, file.path(OUT, "arms", "k_stratified_summary.tsv"))

cat("\nAgreement stratified by whether both settings returned the same lineage count:\n")
print(summ)
cat("\nAgreement stratified by whether both settings used the same k:\n")
print(byk)

verdict <- summ[same_lineage_count == TRUE]
if (nrow(verdict)) {
  cat("\nInterpretation guide:\n")
  for (i in seq_len(nrow(verdict))) {
    m <- verdict$median_abs_rho[i]
    cat(sprintf("  %s: within-stratum median |rho| = %.3f -> %s\n", verdict$arm[i], m,
                if (m >= 0.8) "instability is branch-count selection, not the ordering itself"
                else "the ordering is unstable even holding branch count fixed"))
  }
}
log_step("LINEAGE_STRATIFIED_COMPLETE")
