#!/usr/bin/env Rscript
# 207b_sex_coloc_refresh.R — Refresh of Script 207 with current v2 sex_class +
# SuSiE-canonical PP4 + threshold sweep [0.1, 0.3, 0.5, 0.7, 0.8, 0.9] +
# global BH correction across panel.
#
# Reasons for refresh (audit findings):
#  - Existing sex_coloc_enrichment.csv uses STALE v1 sex_class (n=1,678 F_biased
#    vs current v2 n=1,978). Script 207 was run before the 2026-04-03 v2 rewrite.
#  - Script 207 uses ABF (`PP.H4.abf`) for per-GWAS-group breakdowns
#    (line 93: `max(PP.H4.abf)`) but SuSiE has been canonical since 2026-04-22.
#  - Script 207 tested PP4>0.9 only. Audit reveals OR=1.52, p=0.014 at PP4>0.5
#    (M0 quick analysis).
#
# Output: RNA-seq/results/stratified_causal/sex_coloc_refresh_v2.csv

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- Load canonical sex_class (v3 mashr Bayesian preferred; v2 fallback) ----
# v3 CSV provides v2-compatible `sex_class` alias column; downstream filters unchanged.
sex_v3_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_v3/sex_deg_classification_v3.csv")
sex_v2_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_deg_classification.csv")
sex_class_file <- if (file.exists(sex_v3_path)) sex_v3_path else sex_v2_path
cat(sprintf("sex_class source: %s\n", basename(dirname(sex_class_file))))
sex_class <- fread(sex_class_file)
sex_class[, ensembl_base := sub("[.].*", "", gene)]

gencode <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))
sex_class <- merge(sex_class, gencode[, .(ensembl_base, gene_name)],
                   by = "ensembl_base", all.x = TRUE)

cat(sprintf("sex_class loaded: %d genes\n", nrow(sex_class)))
cat("sex_class breakdown (v2 canonical):\n")
print(sex_class[, .N, by = sex_class][order(-N)])

# ---- Load SuSiE-canonical COLOC + per-GWAS pairs ----------------------------
coloc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc[, pp4 := fifelse(is.na(coloc_best_susie_pp4), coloc_best_pp4, coloc_best_susie_pp4)]

coloc_all <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))

# GWAS family classification
GWAS_DISEASE <- c(
  "2019_31311600_NAFLD_EUR", "2020_32298765_NAFLD_EUR",
  "2020_32514122_Cirrhosis_EAS", "2020_32514122_HCC_EAS",
  "2021_34841290_NAFLD_EUR", "2023_36280732_NAFLD_deCode_EUR",
  "2023_36280732_NAFLD_Intermountain_EUR", "2023_36280732_NAFLD_UKBB_EUR",
  "FinnGen_HCC", "FinnGen_NAFLD", "FinnGen_NASH",
  "Ghouse_Cirrhosis", "Ghouse_HCC"
)
GWAS_ENZYME <- c(
  "BBJ_ALT", "BBJ_AST", "BBJ_GGT",
  "UKBB_ALT", "UKBB_AST", "UKBB_GGT",
  "PanUKBB_AFR_ALT", "PanUKBB_AFR_AST", "PanUKBB_AFR_GGT",
  "PanUKBB_CSA_ALT", "PanUKBB_CSA_AST", "PanUKBB_CSA_GGT",
  "2021_34128465_PDFF_EUR", "2021_34957434_PDFF_EUR",
  "2022_36402844_PDFF_EUR"
)

# Per-gene max PP.H4.susie within each stratum (using SuSiE canonical, NOT ABF)
coloc_all[, pp4 := fifelse(is.na(PP.H4.susie), PP.H4.abf, PP.H4.susie)]
get_pp4 <- function(strat_gwas) {
  d <- coloc_all[gwas_name %in% strat_gwas]
  d[, .(pp4 = max(pp4, na.rm = TRUE)), by = .(gene, ensembl)]
}
pp4_all     <- coloc[, .(gene, pp4)]
pp4_disease <- get_pp4(GWAS_DISEASE)
pp4_enzyme  <- get_pp4(GWAS_ENZYME)

cat(sprintf("\nPP4 sources: all=%d, disease=%d, enzyme=%d genes\n",
            nrow(pp4_all), nrow(pp4_disease), nrow(pp4_enzyme)))

# ---- Fisher panel across thresholds × sex_class × GWAS family ---------------
strata <- list(all = pp4_all, disease = pp4_disease, enzyme = pp4_enzyme)
sex_levels <- c("Female_biased", "Male_biased", "Divergent")
thresholds <- c(0.1, 0.3, 0.5, 0.7, 0.8, 0.9)

results <- list()
for (gwas_strat in names(strata)) {
  d <- merge(sex_class[, .(gene_name, sex_class)],
             strata[[gwas_strat]], by.x = "gene_name", by.y = "gene",
             all.x = TRUE)
  d[is.na(pp4), pp4 := 0]

  for (sl in sex_levels) {
    for (cut in thresholds) {
      a <- sum(d$sex_class == sl & d$pp4 > cut, na.rm = TRUE)
      b <- sum(d$sex_class == sl & d$pp4 <= cut, na.rm = TRUE)
      cc <- sum(d$sex_class != sl & d$pp4 > cut, na.rm = TRUE)
      dd <- sum(d$sex_class != sl & d$pp4 <= cut, na.rm = TRUE)
      tab <- matrix(c(a, b, cc, dd), nrow = 2)
      if (a + cc < 1 || b + dd < 1) next
      ft <- fisher.test(tab)
      results[[length(results) + 1]] <- data.table(
        gwas_stratum = gwas_strat,
        sex_class = sl,
        pp4_cut = cut,
        n_class = a + b,
        n_class_coloc = a,
        n_bg = cc + dd,
        n_bg_coloc = cc,
        odds_ratio = ft$estimate,
        ci_low = ft$conf.int[1],
        ci_high = ft$conf.int[2],
        pvalue = ft$p.value
      )
    }
  }
}
res <- rbindlist(results)
res[, padj_bh := p.adjust(pvalue, "BH")]
setorder(res, gwas_stratum, sex_class, pp4_cut)

cat("\n=== Fisher panel results ===\n")
print(res)

cat(sprintf("\nTotal tests: %d\n", nrow(res)))
cat(sprintf("Significant at padj_bh<0.1: %d\n", sum(res$padj_bh < 0.1)))
cat(sprintf("Significant at padj_bh<0.05: %d\n", sum(res$padj_bh < 0.05)))
cat(sprintf("Nominal p<0.05: %d\n", sum(res$pvalue < 0.05)))

cat("\n=== Top 10 by p-value ===\n")
print(res[order(pvalue)][1:10])

fwrite(res, file.path(OUT, "sex_coloc_refresh_v2.csv"))
cat(sprintf("\nWrote %s\n", file.path(OUT, "sex_coloc_refresh_v2.csv")))

# Quick sanity: also do "any sex_dimorphic" vs "Concordant" overall
cat("\n=== Sanity: ANY sex_dimorphic vs Concordant ===\n")
sex_class[, sex_dim := sex_class != "Concordant"]
for (gwas_strat in names(strata)) {
  d <- merge(sex_class[, .(gene_name, sex_dim)],
             strata[[gwas_strat]], by.x = "gene_name", by.y = "gene", all.x = TRUE)
  d[is.na(pp4), pp4 := 0]
  for (cut in thresholds) {
    tab <- table(d$sex_dim, d$pp4 > cut)
    if (all(dim(tab) == c(2, 2))) {
      ft <- fisher.test(tab)
      cat(sprintf("  %s PP4>%.1f: sex_dim_coloc=%d/%d (%.2f%%); bg_coloc=%d/%d (%.2f%%); OR=%.2f p=%.3g\n",
                  gwas_strat, cut,
                  tab["TRUE","TRUE"], sum(tab["TRUE",]), 100*tab["TRUE","TRUE"]/sum(tab["TRUE",]),
                  tab["FALSE","TRUE"], sum(tab["FALSE",]), 100*tab["FALSE","TRUE"]/sum(tab["FALSE",]),
                  ft$estimate, ft$p.value))
    }
  }
}
