# SCALE of the palindromic sidecar defect, measured with NO external network.
# Non-circular design: compare the local panel against a CONSENSUS of >=3
# INDEPENDENT EUR GWAS. The studies are independent of the panel and (across
# cohort families) of each other, so if several agree with each other and
# disagree with the panel, the panel is the outlier.
# CONTROL = non-palindromic variants, where letters fix orientation, so the
# panel and the studies must agree there whatever the strand handling.
suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
V4 <- file.path(FM, "runs/uniform35_v4_2026-08-04/aggregate/susie_variant_results.tsv.gz")
v <- fread(V4, select = c("study_name","ancestry","chromosome","position","match_class",
                          "palindromic","af","panel_af_a1"), showProgress = FALSE)
v <- v[ancestry == "EUR" & is.finite(af) & is.finite(panel_af_a1) & !is.na(match_class)]
v[, cohort := fifelse(grepl("^MVP_", study_name), "MVP",
              fifelse(grepl("^FinnGen", study_name), "FinnGen",
               fifelse(grepl("^UKBB|^PanUKBB", study_name), "UKBB",
                fifelse(grepl("36280732", study_name), "GCST36280732", "other"))))]
# panel_af_a1 is already the panel frequency of the study's allele1, so `af` and
# `panel_af_a1` are on the same axis and directly comparable.
v[, d := abs(af - panel_af_a1)]
key <- c("chromosome","position")
agg <- v[, .(n_studies = .N, n_cohorts = uniqueN(cohort),
             med_d = median(d), max_spread = diff(range(af)),
             palindromic = palindromic[1]), by = key]
# require >=3 studies from >=2 independent cohort families, and tight agreement
# AMONG the studies (so a real per-study error cannot masquerade as consensus)
ok <- agg[n_studies >= 3 & n_cohorts >= 2 & max_spread <= 0.05]
cat(sprintf("positions with >=3 studies from >=2 cohorts and study-spread <=0.05: %s\n",
            format(nrow(ok), big.mark = ",")))
ok[, panel_outlier := med_d > 0.15]
tab <- ok[, .(n = .N, panel_outlier = sum(panel_outlier),
              pct = round(100*mean(panel_outlier), 3)), by = palindromic]
cat("\n=== does the panel disagree with the study consensus? ===\n")
print(tab[order(palindromic)])
if (nrow(tab) == 2) {
  p <- tab[palindromic == TRUE]; c_ <- tab[palindromic == FALSE]
  ft <- fisher.test(matrix(c(p$panel_outlier, p$n - p$panel_outlier,
                             c_$panel_outlier, c_$n - c_$panel_outlier), nrow = 2))
  cat(sprintf("\npalindromic %.3f%% vs NON-palindromic control %.3f%%  ratio %.1fx  Fisher p = %.3g\n",
              p$pct, c_$pct, p$pct/max(c_$pct,1e-9), ft$p.value))
  cat(sprintf("estimated affected palindromic positions (EUR, this measurable subset): %d\n",
              p$panel_outlier))
}
fwrite(ok, file.path(FM, "tmp_diag/sidecar_scale.tsv"), sep = "\t")
