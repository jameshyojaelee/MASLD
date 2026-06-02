#!/usr/bin/env Rscript
# Canonical loci per ancestry under the policy:
#   UKBB EUR  → PolyFun LD
#   All other → 1KG LD (TOP-LD shown as alternative)
#
# UKBB EUR studies: UKBB_ALT, UKBB_AST, UKBB_GGT, 2023_36280732_NAFLD_UKBB_EUR,
#   2021_34128465_PDFF_EUR, 2021_34957434_PDFF_EUR, 2022_36402844_PDFF_EUR
#   (all three PDFF studies are UKBB MRI per CLAUDE.md).
suppressPackageStartupMessages({ library(data.table) })

ukbb_studies <- c(
  "UKBB_ALT", "UKBB_AST", "UKBB_GGT",
  "2023_36280732_NAFLD_UKBB_EUR",
  "2021_34128465_PDFF_EUR",
  "2021_34957434_PDFF_EUR",
  "2022_36402844_PDFF_EUR")

reg <- fread("GWAS/finemapping/config/gwas_registry.tsv",
             select = c("study_name","ancestry"))
setnames(reg, "study_name", "study")
reg[, panel_choice := fifelse(study %in% ukbb_studies, "PolyFun", "1KG")]

# 1KG converged loci per study
k <- fread("GWAS/finemapping/results/combined_finemapping_1kg.csv",
           select = c("study","locus","susie_converged"))
k_per <- k[susie_converged == TRUE, .(n_loci = uniqueN(locus)), by = study]

# TOP-LD converged loci per study (alternative)
t <- fread("GWAS/finemapping/results/combined_finemapping_topld.csv",
           select = c("study","locus","susie_converged"))
t_per <- t[susie_converged == TRUE, .(n_loci = uniqueN(locus)), by = study]

# PolyFun re-aggregated to lead-SNP windows (verified script v2)
p_per <- fread("scripts/figures/sketches_fig3_intro/data/polyfun_loci_per_study_v2.csv")
setnames(p_per, "n_loci_polyfun", "n_loci")

# Canonical: PolyFun for UKBB EUR studies, 1KG for the rest
canonical <- copy(reg)
canonical <- merge(canonical, p_per, by = "study", all.x = TRUE, suffixes = c("","_p"))
setnames(canonical, "n_loci", "n_loci_polyfun")
canonical <- merge(canonical, k_per, by = "study", all.x = TRUE)
setnames(canonical, "n_loci", "n_loci_1kg")
canonical <- merge(canonical, t_per, by = "study", all.x = TRUE)
setnames(canonical, "n_loci", "n_loci_topld")

canonical[, n_loci_canonical := fifelse(panel_choice == "PolyFun",
                                          n_loci_polyfun,
                                          n_loci_1kg)]

# Optional: TOP-LD-as-alternative for non-UKBB
canonical[, n_loci_canonical_topld_alt := fifelse(panel_choice == "PolyFun",
                                                   n_loci_polyfun,
                                                   n_loci_topld)]

fwrite(canonical, "scripts/figures/sketches_fig3_intro/data/canonical_loci_per_study.csv")

cat("=== Canonical per-study counts ===\n")
print(canonical[, .(study, ancestry, panel_choice, n_loci_polyfun,
                    n_loci_1kg, n_loci_topld, n_loci_canonical)])

cat("\n=== Per-ancestry totals (canonical: PolyFun for UKBB EUR, 1KG otherwise) ===\n")
agg <- canonical[, .(n_loci = sum(n_loci_canonical, na.rm = TRUE),
                     n_studies = uniqueN(study[!is.na(n_loci_canonical) &
                                                n_loci_canonical > 0])),
                 by = ancestry]
print(agg)
cat("\nGRAND TOTAL canonical:", sum(agg$n_loci), "loci across",
    sum(agg$n_studies), "studies\n")

cat("\n=== Per-ancestry totals (alternative: PolyFun for UKBB EUR, TOP-LD otherwise) ===\n")
agg2 <- canonical[, .(n_loci = sum(n_loci_canonical_topld_alt, na.rm = TRUE),
                      n_studies = uniqueN(study[!is.na(n_loci_canonical_topld_alt) &
                                                 n_loci_canonical_topld_alt > 0])),
                  by = ancestry]
print(agg2)
cat("\nGRAND TOTAL alternative:", sum(agg2$n_loci), "\n")

cat("\n=== Side-by-side comparison: PolyFun vs 1KG on UKBB EUR studies ===\n")
print(canonical[panel_choice == "PolyFun",
                .(study, n_loci_polyfun, n_loci_1kg, n_loci_topld)])
