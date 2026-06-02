#!/usr/bin/env Rscript
# 207_sex_stratified_coloc.R — Sex-Stratified COLOC Enrichment
#
# Projects existing COLOC/TWAS results onto sex-stratified DEG classification
# to determine whether genetically causal genes (PP.H4 > 0.9) are enriched
# among female-specific, male-specific, or sex-divergent DEGs.
#
# Analysis:
#   1. Fisher's exact tests: COLOC genes enriched in each sex class?
#   2. GSEA-style: rank genes by PP.H4, test enrichment of sex gene sets (fgsea)
#   3. Per-GWAS phenotype: enzyme GWAS (ALT/AST/GGT) vs disease GWAS
#      (NAFLD/cirrhosis/HCC) show different sex enrichment?
#   4. Continuous sex-causal score: PP4 × |female_logFC - male_logFC|
#
# Inputs:
#   - GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv (599 genes PP.H4>0.9)
#   - GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv (per-GWAS PP.H4)
#   - RNA-seq/Human/.../results/integration/sex_deg_classification.csv (sex DEGs)
#   - RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv (TWAS z-scores)
#   - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv (gene ID mapping)
#
# Outputs:
#   - RNA-seq/results/stratified_causal/sex_coloc_enrichment.csv
#   - RNA-seq/results/stratified_causal/sex_coloc_gsea.csv
#   - RNA-seq/results/stratified_causal/sex_causal_scores.csv
#
# SLURM: cpu partition, 4 CPUs, 32GB, 48h

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(fgsea)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

outdir <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# 1. Load gene ID mapping from atlas
# ===========================================================================
cat("=== Loading gene ID mapping from multi-evidence atlas ===\n")
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("ensembl_id", "human_symbol"))
cat("  Atlas:", nrow(atlas), "genes\n")

# ===========================================================================
# 2. Load COLOC gene-level results
# ===========================================================================
cat("\n=== Loading COLOC gene-level results ===\n")
coloc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc <- coloc[gene != "" & !is.na(gene)]
# T0.4 (2026-04-22): prefer SuSiE PP4 over legacy ABF. CLAUDE.md canonical: SuSiE is
# canonical across 28 GWAS; ABF retained only as fallback for SuSiE-unconverged genes.
# Keep `coloc_best_pp4` as the in-script handle (referenced 14x below) but sourced from
# SuSiE where available. `coloc_best_abf_pp4` preserved for provenance.
if ("coloc_best_susie_pp4" %in% names(coloc)) {
  coloc[, coloc_best_abf_pp4 := coloc_best_pp4]
  coloc[, coloc_best_pp4 := fifelse(!is.na(coloc_best_susie_pp4),
                                    coloc_best_susie_pp4, coloc_best_pp4)]
  cat("  T0.4: coloc_best_pp4 now sourced from SuSiE with ABF fallback\n")
}
cat("  COLOC:", nrow(coloc), "genes total\n")
cat("  PP.H4 > 0.9:", sum(coloc$coloc_best_pp4 > 0.9), "\n")
cat("  PP.H4 > 0.8:", sum(coloc$coloc_best_pp4 > 0.8), "\n")

# ===========================================================================
# 3. Load per-GWAS COLOC results (for phenotype stratification)
# ===========================================================================
cat("\n=== Loading per-GWAS COLOC results ===\n")
coloc_all <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
cat("  Per-GWAS COLOC:", nrow(coloc_all), "gene-GWAS entries\n")
cat("  Unique GWAS:", length(unique(coloc_all$gwas_name)), "\n")

# Classify GWAS into phenotype groups
coloc_all[, gwas_group := fcase(
  grepl("ALT|AST|GGT", gwas_name, ignore.case = TRUE), "enzyme",
  grepl("PDFF", gwas_name, ignore.case = TRUE), "imaging",
  grepl("NAFLD|NASH|Anstee", gwas_name, ignore.case = TRUE), "disease",
  grepl("Cirrhosis|cirrhosis", gwas_name), "cirrhosis",
  grepl("HCC|hcc", gwas_name, ignore.case = TRUE), "hcc",
  grepl("Obesity|obesity", gwas_name, ignore.case = TRUE), "obesity",
  default = "other"
)]
cat("  GWAS group distribution:\n")
print(coloc_all[, .N, by = gwas_group][order(-N)])

# Best PP.H4 per gene per GWAS group
coloc_by_group <- coloc_all[, .(
  pp4_best = max(PP.H4.abf, na.rm = TRUE),
  n_gwas_h4_05 = sum(PP.H4.abf > 0.5, na.rm = TRUE),
  best_gwas = gwas_name[which.max(PP.H4.abf)]
), by = .(gene, ensembl, gwas_group)]

# ===========================================================================
# 4. Load sex-stratified DEG classification
# ===========================================================================
cat("\n=== Loading sex DEG classification ===\n")
sex_file <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_deg_classification.csv")
sex_degs <- fread(sex_file)
cat("  Sex DEGs:", nrow(sex_degs), "genes\n")

# Strip version from Ensembl IDs for joining
sex_degs[, ensembl_id := sub("\\..*", "", gene)]
setnames(sex_degs, "gene", "ensembl_versioned")

cat("  Sex class distribution:\n")
print(sex_degs[, .N, by = sex_class][order(-N)])

# Auto-detect classification scheme (v1: significance-based; v2: interaction-based)
is_v2 <- any(sex_degs$sex_class %in% c("Female_biased", "Male_biased", "Concordant"))
if (is_v2) {
  cat("  Detected v2 (interaction-based) classification\n")
  SEX_CLASSES_ALL <- c("Female_biased", "Male_biased", "Divergent", "Concordant")
  SEX_CLASSES_SIG <- c("Female_biased", "Male_biased", "Divergent")
  SEX_NS          <- "Concordant"
} else {
  cat("  Detected v1 (significance-based) classification\n")
  SEX_CLASSES_ALL <- c("Female_specific", "Male_specific", "Divergent", "Shared", "Not_significant")
  SEX_CLASSES_SIG <- c("Female_specific", "Male_specific", "Divergent")
  SEX_NS          <- "Not_significant"
}

# ===========================================================================
# 5. Load TWAS results
# ===========================================================================
cat("\n=== Loading TWAS results ===\n")
twas <- fread(file.path(BASE, "RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv"))
twas_best <- twas[, .(twas_z = zscore[which.min(pvalue)],
                       twas_p = min(pvalue)),
                  by = .(gene_name)]
cat("  TWAS:", nrow(twas_best), "genes\n")

# ===========================================================================
# 6. Build merged table: COLOC + sex classification
# ===========================================================================
cat("\n=== Building merged gene table ===\n")

# Map COLOC gene symbols to Ensembl IDs via atlas
coloc_mapped <- merge(coloc, atlas, by.x = "gene", by.y = "human_symbol", all.x = TRUE)
# Also try direct Ensembl match for genes that didn't map via symbol
coloc_mapped[is.na(ensembl_id), ensembl_id := ensembl]

# Join with sex DEG classification
merged <- merge(coloc_mapped, sex_degs, by = "ensembl_id", all = TRUE)

# Fill in PP.H4 as 0 for genes not in COLOC
merged[is.na(coloc_best_pp4), coloc_best_pp4 := 0]

# Define universe: genes present in both COLOC universe and sex classification
# (i.e., genes that were tested in both analyses)
merged <- merged[!is.na(sex_class) & sex_class != ""]
cat("  Merged table:", nrow(merged), "genes in universe\n")
cat("  COLOC genes (PP.H4>0.9) in universe:", sum(merged$coloc_best_pp4 > 0.9), "\n")

# ===========================================================================
# 7. Fisher's exact tests: COLOC enrichment per sex class
# ===========================================================================
cat("\n=== Fisher's Exact Tests: COLOC × Sex DEG Class ===\n")

fisher_results <- list()
coloc_threshold <- 0.9
coloc_genes <- merged$coloc_best_pp4 > coloc_threshold

sex_classes <- SEX_CLASSES_ALL

for (sc in sex_classes) {
  in_class <- merged$sex_class == sc

  # 2×2 contingency table
  a <- sum(coloc_genes & in_class)   # COLOC+ & class+
  b <- sum(coloc_genes & !in_class)  # COLOC+ & class-
  c <- sum(!coloc_genes & in_class)  # COLOC- & class+
  d <- sum(!coloc_genes & !in_class) # COLOC- & class-

  mat <- matrix(c(a, b, c, d), nrow = 2)
  ft <- fisher.test(mat)

  fisher_results[[sc]] <- data.table(
    sex_class = sc,
    n_class = sum(in_class),
    n_coloc_in_class = a,
    n_coloc_total = sum(coloc_genes),
    pct_coloc_in_class = round(100 * a / max(sum(in_class), 1), 2),
    pct_coloc_overall = round(100 * sum(coloc_genes) / nrow(merged), 2),
    odds_ratio = round(ft$estimate, 3),
    ci_lower = round(ft$conf.int[1], 3),
    ci_upper = round(ft$conf.int[2], 3),
    pvalue = ft$p.value,
    direction = ifelse(ft$estimate > 1, "enriched", "depleted")
  )

  cat(sprintf("  %-20s n=%5d  coloc_in_class=%3d  OR=%.2f [%.2f-%.2f]  p=%.2e  %s\n",
              sc, sum(in_class), a, ft$estimate,
              ft$conf.int[1], ft$conf.int[2],
              ft$p.value,
              ifelse(ft$estimate > 1, "ENRICHED", "depleted")))
}

fisher_dt <- rbindlist(fisher_results)
fisher_dt[, padj := p.adjust(pvalue, method = "BH")]

# Also test at stricter threshold (PP.H4 > 0.8)
cat("\n  --- Sensitivity: PP.H4 > 0.8 ---\n")
coloc_strict <- merged$coloc_best_pp4 > 0.8
for (sc in SEX_CLASSES_SIG) {
  in_class <- merged$sex_class == sc
  a <- sum(coloc_strict & in_class)
  b <- sum(coloc_strict & !in_class)
  c <- sum(!coloc_strict & in_class)
  d <- sum(!coloc_strict & !in_class)
  mat <- matrix(c(a, b, c, d), nrow = 2)
  ft <- fisher.test(mat)
  cat(sprintf("  %-20s PP4>0.8: n=%3d  OR=%.2f  p=%.2e\n",
              sc, a, ft$estimate, ft$p.value))
}

# ===========================================================================
# 8. GSEA-style: Rank all genes by PP.H4, test enrichment of sex gene sets
# ===========================================================================
cat("\n=== GSEA-style: sex gene sets ranked by PP.H4 ===\n")

# Build ranking vector: gene symbol → PP.H4 (use atlas to map back)
# Need gene symbols for fgsea pathway labels
merged_with_symbol <- merge(merged, atlas, by = "ensembl_id", all.x = TRUE)

# Use gene symbol for ranking; fall back to ensembl_id
merged_with_symbol[, rank_id := fifelse(
  !is.na(human_symbol) & human_symbol != "",
  human_symbol,
  ensembl_id
)]

# Deduplicate: keep gene with highest PP.H4 if duplicate symbols
merged_with_symbol <- merged_with_symbol[order(-coloc_best_pp4)]
merged_with_symbol <- merged_with_symbol[!duplicated(rank_id)]

# Create named ranking vector (PP.H4 as statistic)
ranks <- setNames(merged_with_symbol$coloc_best_pp4, merged_with_symbol$rank_id)

# Add small jitter to break ties at 0 (fgsea needs variability)
set.seed(42)
ranks <- ranks + runif(length(ranks), 0, 1e-10)

# Sort descending
ranks <- sort(ranks, decreasing = TRUE)

# Define gene sets: sex classes → gene symbols
pathways <- list()
for (sc in SEX_CLASSES_SIG) {
  ens_ids <- merged_with_symbol[sex_class == sc, rank_id]
  if (length(ens_ids) >= 5) {
    pathways[[sc]] <- ens_ids
  }
}

cat("  Gene set sizes:", paste(sapply(pathways, length), collapse = ", "), "\n")
cat("  Ranking vector length:", length(ranks), "\n")

# Run fgsea
gsea_res <- fgsea(pathways = pathways,
                  stats = ranks,
                  minSize = 5,
                  maxSize = 10000,
                  nPermSimple = 10000)

gsea_dt <- as.data.table(gsea_res)
gsea_dt[, leadingEdge := sapply(leadingEdge, function(x) paste(head(x, 10), collapse = ";"))]

cat("\n  GSEA results:\n")
for (i in seq_len(nrow(gsea_dt))) {
  cat(sprintf("  %-20s NES=%.3f  padj=%.2e  size=%d\n",
              gsea_dt$pathway[i], gsea_dt$NES[i],
              gsea_dt$padj[i], gsea_dt$size[i]))
}

# ===========================================================================
# 9. Per-GWAS phenotype: sex enrichment by GWAS category
# ===========================================================================
cat("\n=== Per-GWAS phenotype: sex enrichment ===\n")

# For each GWAS group, get genes with PP.H4 > 0.9 and test sex enrichment
gwas_groups <- c("enzyme", "disease", "cirrhosis", "hcc", "imaging")

# Get per-group gene lists
per_group_genes <- coloc_by_group[pp4_best > 0.9, .(gene = unique(gene)), by = gwas_group]

# Map gene symbols to ensembl for joining
per_group_genes <- merge(per_group_genes, atlas,
                         by.x = "gene", by.y = "human_symbol", all.x = TRUE)

gwas_sex_results <- list()

for (grp in gwas_groups) {
  grp_ens <- per_group_genes[gwas_group == grp & !is.na(ensembl_id), ensembl_id]
  if (length(grp_ens) < 5) {
    cat(sprintf("  Skipping %s — only %d genes with PP.H4>0.9\n", grp, length(grp_ens)))
    next
  }

  for (sc in SEX_CLASSES_SIG) {
    in_class <- merged$sex_class == sc
    in_coloc_grp <- merged$ensembl_id %in% grp_ens

    a <- sum(in_coloc_grp & in_class)
    b <- sum(in_coloc_grp & !in_class)
    c <- sum(!in_coloc_grp & in_class)
    d <- sum(!in_coloc_grp & !in_class)

    mat <- matrix(c(a, b, c, d), nrow = 2)
    ft <- fisher.test(mat)

    gwas_sex_results[[paste(grp, sc, sep = "_")]] <- data.table(
      gwas_group = grp,
      sex_class = sc,
      n_coloc_group = sum(in_coloc_grp),
      n_coloc_in_class = a,
      odds_ratio = round(ft$estimate, 3),
      ci_lower = round(ft$conf.int[1], 3),
      ci_upper = round(ft$conf.int[2], 3),
      pvalue = ft$p.value
    )
  }
}

gwas_sex_dt <- rbindlist(gwas_sex_results)
gwas_sex_dt[, padj := p.adjust(pvalue, method = "BH")]

cat("\n  Per-GWAS phenotype × sex class:\n")
for (i in seq_len(nrow(gwas_sex_dt))) {
  r <- gwas_sex_dt[i]
  cat(sprintf("  %-12s × %-20s n_overlap=%3d  OR=%.2f  p=%.2e  padj=%.2e\n",
              r$gwas_group, r$sex_class, r$n_coloc_in_class,
              r$odds_ratio, r$pvalue, r$padj))
}

# ===========================================================================
# 10. Per-GWAS phenotype GSEA: rank by per-group PP.H4
# ===========================================================================
cat("\n=== Per-GWAS phenotype GSEA ===\n")

gwas_gsea_results <- list()

for (grp in gwas_groups) {
  # Get per-gene best PP.H4 for this GWAS group
  grp_pp4 <- coloc_by_group[gwas_group == grp, .(gene, pp4 = pp4_best)]
  grp_pp4 <- merge(grp_pp4, atlas, by.x = "gene", by.y = "human_symbol", all.x = TRUE)

  # Map to merged table IDs
  grp_merged <- merge(
    merged_with_symbol[, .(rank_id, ensembl_id)],
    grp_pp4[, .(ensembl_id, pp4)],
    by = "ensembl_id", all.x = TRUE
  )
  grp_merged[is.na(pp4), pp4 := 0]
  grp_merged <- grp_merged[!duplicated(rank_id)]

  grp_ranks <- setNames(grp_merged$pp4 + runif(nrow(grp_merged), 0, 1e-10),
                         grp_merged$rank_id)
  grp_ranks <- sort(grp_ranks, decreasing = TRUE)

  if (length(grp_ranks) < 100) next

  grp_gsea <- fgsea(pathways = pathways,
                    stats = grp_ranks,
                    minSize = 5,
                    maxSize = 10000,
                    nPermSimple = 10000)

  grp_gsea_dt <- as.data.table(grp_gsea)
  grp_gsea_dt[, gwas_group := grp]
  grp_gsea_dt[, leadingEdge := sapply(leadingEdge, function(x) paste(head(x, 10), collapse = ";"))]
  gwas_gsea_results[[grp]] <- grp_gsea_dt
}

gwas_gsea_dt <- rbindlist(gwas_gsea_results, fill = TRUE)

cat("  Per-GWAS GSEA results:\n")
for (i in seq_len(nrow(gwas_gsea_dt))) {
  r <- gwas_gsea_dt[i]
  cat(sprintf("  %-12s × %-20s NES=%.3f  padj=%.2e\n",
              r$gwas_group, r$pathway, r$NES, r$padj))
}

# ===========================================================================
# 11. Continuous sex-causal score: PP4 × |female_logFC - male_logFC|
# ===========================================================================
cat("\n=== Computing sex-causal scores ===\n")

scores <- merged[, .(
  ensembl_id,
  sex_class,
  coloc_best_pp4,
  coloc_best_gwas,
  logFC_F,
  logFC_M,
  padj_F,
  padj_M,
  lfc_diff
)]

# Compute continuous score
scores[, sex_causal_score := coloc_best_pp4 * abs(logFC_F - logFC_M)]

# Also compute signed version (positive = female-biased causal)
scores[, sex_causal_signed := coloc_best_pp4 * (logFC_F - logFC_M)]

# Add gene symbols via atlas
scores <- merge(scores, atlas, by = "ensembl_id", all.x = TRUE)

# Add TWAS z-score
scores <- merge(scores, twas_best[, .(gene_name, twas_z, twas_p)],
                by.x = "human_symbol", by.y = "gene_name", all.x = TRUE)

# Rank
scores <- scores[order(-sex_causal_score)]
scores[, rank := .I]

cat("  Genes with sex_causal_score > 0:", sum(scores$sex_causal_score > 0), "\n")
cat("  Top 20 sex-causal genes:\n")
top20 <- scores[rank <= 20, .(
  rank, human_symbol, sex_class, coloc_best_pp4,
  logFC_F = round(logFC_F, 3),
  logFC_M = round(logFC_M, 3),
  sex_causal_score = round(sex_causal_score, 3)
)]
print(top20)

# ===========================================================================
# 12. Summary statistics
# ===========================================================================
cat("\n=== Summary Statistics ===\n")

# Proportion of COLOC genes in each sex class
n_coloc <- sum(merged$coloc_best_pp4 > 0.9)
for (sc in sex_classes) {
  n <- sum(merged$coloc_best_pp4 > 0.9 & merged$sex_class == sc)
  pct <- round(100 * n / max(n_coloc, 1), 1)
  cat(sprintf("  COLOC PP4>0.9 in %-20s: %3d / %d (%.1f%%)\n", sc, n, n_coloc, pct))
}

# Sex-causal score distribution among COLOC genes
coloc_scored <- scores[coloc_best_pp4 > 0.9]
cat("\n  Sex-causal score among COLOC genes (PP4>0.9):\n")
print(summary(coloc_scored$sex_causal_score))

# Correlation: PP4 vs |LFC diff|
cor_test <- cor.test(scores$coloc_best_pp4, abs(scores$logFC_F - scores$logFC_M),
                     method = "spearman", exact = FALSE)
cat(sprintf("\n  Spearman cor(PP4, |LFC_F - LFC_M|): rho=%.3f  p=%.2e\n",
            cor_test$estimate, cor_test$p.value))

# ===========================================================================
# 13. Combine enrichment results (Fisher + GSEA + per-GWAS)
# ===========================================================================
cat("\n=== Saving results ===\n")

# Enrichment CSV: Fisher's exact + per-GWAS stratification
enrichment_out <- list(
  overall_fisher = fisher_dt,
  per_gwas_fisher = gwas_sex_dt
)
enrichment_combined <- rbindlist(list(
  fisher_dt[, .(test_type = "overall_fisher", sex_class, n_class,
                n_coloc_in_class, odds_ratio, ci_lower, ci_upper,
                pvalue, padj, direction)],
  gwas_sex_dt[, .(test_type = paste0("gwas_", gwas_group), sex_class,
                  n_class = n_coloc_group,
                  n_coloc_in_class, odds_ratio, ci_lower, ci_upper,
                  pvalue, padj, direction = ifelse(odds_ratio > 1, "enriched", "depleted"))]
), fill = TRUE)

fwrite(enrichment_combined,
       file.path(outdir, "sex_coloc_enrichment.csv"))
cat("  sex_coloc_enrichment.csv:", nrow(enrichment_combined), "rows\n")

# GSEA CSV: overall + per-GWAS
gsea_combined <- rbindlist(list(
  gsea_dt[, .(gwas_group = "overall", pathway, pval, padj, NES, size, leadingEdge)],
  gwas_gsea_dt[, .(gwas_group, pathway, pval, padj, NES, size, leadingEdge)]
), fill = TRUE)

fwrite(gsea_combined,
       file.path(outdir, "sex_coloc_gsea.csv"))
cat("  sex_coloc_gsea.csv:", nrow(gsea_combined), "rows\n")

# Sex-causal scores: all genes
fwrite(scores,
       file.path(outdir, "sex_causal_scores.csv"))
cat("  sex_causal_scores.csv:", nrow(scores), "rows\n")

cat("\nDone.\n")
