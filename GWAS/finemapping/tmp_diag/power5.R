# Is "0 opposites in finemap" actually inconsistent with the COLOC-measured
# 0.25% flip rate? Only if finemap had power to see it. Calibrate the two
# pipelines using CONTROL studies neither of us thinks is defective.
suppressPackageStartupMessages(library(data.table))
# finemap recorded rates (from reconcile3) vs my COLOC sweep rates, control studies
ctrl <- data.table(
  study   = c("MVP_NAFLD_EUR","MVP_ALT_EUR","MVP_AST_EUR","2021_34957434_PDFF_EUR"),
  fm_pct  = c(0.0151, 0.0157, 0.0180, 0.0130),
  coloc_pct = c(0.0543, 0.0496, 0.0458, 0.0430))
ctrl[, ratio := coloc_pct / fm_pct]
cat("=== calibration on control studies (COLOC rate / finemap rate) ===\n"); print(ctrl)
k <- median(ctrl$ratio)
cat(sprintf("\nmedian offset k = %.2f  (COLOC sees ~%.1fx the rate finemap does)\n\n", k, k))
# The five, with their AF-RESOLVED palindrome counts (abstentions cannot flag)
five <- data.table(
  study    = c("FinnGen_NAFLD","FinnGen_NASH","2023_36280732_NAFLD_UKBB_EUR",
               "2023_36280732_NAFLD_deCode_EUR","2023_36280732_NAFLD_Intermountain_EUR"),
  resolved = c(3468, 2741, 3517, 891, 497),
  observed = c(0,0,0,0,0),
  coloc_pct= c(0.2542, 0.2453, 0.2491, 0.2565, 0.2477))
five[, expected_fm := resolved * (coloc_pct/k) / 100]
five[, p_zero := ppois(0, expected_fm)]
cat("=== if the COLOC-measured flip rate were real, what would finemap have seen? ===\n")
print(five[, .(study, resolved, coloc_pct, expected = round(expected_fm,2),
               observed, p_obs_zero = round(p_zero,4))])
tot_e <- sum(five$expected_fm)
cat(sprintf("\nPOOLED across the five: expected %.2f opposites, observed 0, Poisson p = %.4f\n",
            tot_e, ppois(0, tot_e)))
cat(sprintf("VERDICT: %s\n", if (ppois(0, tot_e) > 0.05)
  "0 is NOT significantly below expectation -- the pipelines do NOT demonstrably disagree" else
  "0 IS significantly below expectation -- the pipelines genuinely disagree"))
