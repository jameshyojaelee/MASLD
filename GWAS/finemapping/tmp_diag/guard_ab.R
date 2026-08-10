#!/usr/bin/env Rscript
# A/B the palindrome AF guard: OLD (uncalibrated) vs NEW (common.R calibration).
# CONTROL = non-palindromic variants.  Their orientation is fixed by letters
# alone, so EVERY flag on them is by definition a false positive.  A guard that
# is detecting strand should flag ~0% of controls; one that is detecting
# ancestry frequency drift flags controls and palindromes at the same rate.
suppressPackageStartupMessages(library(data.table))
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
EQ   <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/broadaway_eqtl"
afsrc <- fread(file.path(BASE, "config/gwas_af_sources.tsv"))
args <- commandArgs(trailingOnly = TRUE)
STUDIES <- strsplit(args[1], ",")[[1]]; CHRS <- as.integer(strsplit(args[2], ",")[[1]])
MAFMAX <- 0.40; MAXDIST <- 0.15; MARGIN <- 0.10

res <- list()
for (CHR in CHRS) {
  eq <- fread(file.path(EQ, sprintf("chr%d_marginal_summary_results.tsv", CHR)))
  eq[, merge_key := paste(CHR, POS, sep = ":")]
  ev <- unique(eq[, .(merge_key, eqtl_ea = EA, eqtl_nea = NEA, eqtl_eaf = EAF)], by = "merge_key")
  for (S in STUDIES) {
    h <- afsrc[study_name == S]
    if (!nrow(h) || !file.exists(h$af_sumstats_path[[1]])) next
    g <- fread(h$af_sumstats_path[[1]])[chromosome == CHR]
    if (!"af" %in% names(g) || !nrow(g)) next
    g[, merge_key := paste(chromosome, position, sep = ":")]
    gv <- unique(g[, .(merge_key, gwas_a1 = allele1, gwas_a2 = allele2, gwas_af = af)], by = "merge_key")
    m <- merge(ev, gv, by = "merge_key")
    m[, allele_match := (eqtl_ea == gwas_a1 & eqtl_nea == gwas_a2)]
    m[, allele_flip  := (eqtl_ea == gwas_a2 & eqtl_nea == gwas_a1)]
    m <- m[allele_match | allele_flip]
    m[, is_pal := (gwas_a1 %in% c("A","T") & gwas_a2 %in% c("A","T")) |
                  (gwas_a1 %in% c("C","G") & gwas_a2 %in% c("C","G"))]
    m[, eqtl_or := fifelse(allele_flip, 1 - eqtl_eaf, eqtl_eaf)]
    m <- m[is.finite(gwas_af) & is.finite(eqtl_or)]
    if (!nrow(m)) next
    m[, d_same := abs(gwas_af - eqtl_or)]; m[, d_opp := abs(gwas_af - (1 - eqtl_or))]
    m[, old_flag := d_same > MAXDIST & d_opp < d_same]                       # superseded rule
    m[, testable := pmin(gwas_af, 1-gwas_af) <= MAFMAX &                     # new rule
                    pmin(eqtl_or, 1-eqtl_or) <= MAFMAX &
                    pmin(d_same, d_opp) <= MAXDIST]
    m[, new_flag := testable & (d_same - d_opp) >= MARGIN]
    for (grp in c(FALSE, TRUE)) { x <- m[is_pal == grp]; if (!nrow(x)) next
      res[[length(res)+1]] <- data.table(study = S, ancestry = h$ancestry[[1]], chr = CHR,
        grp = if (grp) "palindromic" else "CONTROL(nonpal)", n = nrow(x),
        old_pct = 100*mean(x$old_flag), new_pct = 100*mean(x$new_flag)) }
  }
}
r <- rbindlist(res)
agg <- r[, .(n = sum(n), old_pct = 100*sum(n*old_pct/100)/sum(n),
             new_pct = 100*sum(n*new_pct/100)/sum(n)), by = .(study, ancestry, grp)]
setorder(agg, ancestry, study, grp)
cat(sprintf("\nchromosomes: %s   MAF<=%.2f  maxdist<=%.2f  margin>=%.2f\n\n",
            paste(CHRS, collapse=","), MAFMAX, MAXDIST, MARGIN))
cat(sprintf("%-22s %-4s %-18s %9s %9s %9s\n","study","anc","group","n","OLD %","NEW %"))
for (i in seq_len(nrow(agg))) cat(sprintf("%-22s %-4s %-18s %9d %9.3f %9.3f\n",
  agg$study[i], agg$ancestry[i], agg$grp[i], agg$n[i], agg$old_pct[i], agg$new_pct[i]))
ctl <- agg[grp == "CONTROL(nonpal)"]
cat(sprintf("\n=== FALSE-POSITIVE RATE ON CONTROLS (every flag here is wrong) ===\n"))
cat(sprintf("  OLD rule: mean %.3f%%   worst %.3f%% (%s)\n", mean(ctl$old_pct), max(ctl$old_pct), ctl$study[which.max(ctl$old_pct)]))
cat(sprintf("  NEW rule: mean %.3f%%   worst %.3f%% (%s)\n", mean(ctl$new_pct), max(ctl$new_pct), ctl$study[which.max(ctl$new_pct)]))
cat(sprintf("  reduction: %.1fx\n", mean(ctl$old_pct)/max(mean(ctl$new_pct), 1e-9)))
