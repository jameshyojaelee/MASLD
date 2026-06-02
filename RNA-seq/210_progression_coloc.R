#!/usr/bin/env Rscript
# 210_progression_coloc.R — Progression-stratified COLOC enrichment
#
# Part of the Stratified Causal Architecture Pipeline (Module B).
# Tests whether genetically causal genes (COLOC PP.H4 > 0.9) are enriched
# at specific fibrosis stage transitions, correlates stage-specificity (tau)
# with COLOC strength, classifies genes as onset/progression/pan-stage,
# and compares cirrhosis-specific vs enzyme COLOC.
#
# Inputs:
#   - COLOC gene-level:       GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv
#   - COLOC per-GWAS:         GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv
#   - Stage-unique signatures: .../results/staging_classifier/stage_unique_signatures.csv
#   - Fibrosis OVR DEGs:      .../results/staging_classifier/one_vs_rest_fibrosis_dream.csv
#   - Multi-evidence atlas:   RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
#
# Outputs (in RNA-seq/results/stratified_causal/):
#   - progression_coloc_enrichment.csv    Per-transition Fisher enrichment
#   - progression_coloc_by_transition.csv Per-gene COLOC x transition annotation
#   - onset_vs_progression_genes.csv      Onset / progression / pan-stage classification
#
# DEG THRESHOLD SYSTEM (Two-Tier):
#   Tier 1 (Primary): padj < 0.05, |logFC| > 0.5 — main dream DEGs (Script 05b)
#   Tier 2 (Progression): padj < 0.05, no LFC filter — binary and adjacent contrasts
#   Rationale: Binary contrasts pool multiple stages, diluting per-gene fold changes.
#              Adjacent transitions have lower N (200-400 vs 1,444), making LFC estimates
#              noisier. Cross-contrast comparisons use rank-based enrichment (fgsea)
#              to avoid confounding power with biology.
#   See: figures/supplementary/sensitivity/figS_deg_threshold_landscape.pdf
#
# SLURM: cpu partition, 4 CPUs, 32GB, 48h

suppressPackageStartupMessages({
  library(data.table)
})

cat("=== 210_progression_coloc.R ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

outdir <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# 1. Load COLOC gene-level results
# ===========================================================================
cat("=== Loading COLOC gene-level results ===\n")
coloc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc <- coloc[gene != "" & !is.na(gene)]
# T0.4 (2026-04-22): prefer SuSiE PP4 over legacy ABF.
if ("coloc_best_susie_pp4" %in% names(coloc)) {
  coloc[, coloc_best_abf_pp4 := coloc_best_pp4]
  coloc[, coloc_best_pp4 := fifelse(!is.na(coloc_best_susie_pp4),
                                    coloc_best_susie_pp4, coloc_best_pp4)]
  cat("  T0.4: coloc_best_pp4 now sourced from SuSiE with ABF fallback\n")
}
cat("  Total genes in COLOC:", nrow(coloc), "\n")
cat("  PP.H4 > 0.9:", sum(coloc$coloc_best_pp4 > 0.9), "\n")
cat("  PP.H4 > 0.8:", sum(coloc$coloc_best_pp4 > 0.8), "\n")

# ===========================================================================
# 2. Load per-GWAS COLOC for cirrhosis vs enzyme comparison
# ===========================================================================
cat("\n=== Loading per-GWAS COLOC ===\n")
coloc_all <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
cat("  Total gene-GWAS entries:", nrow(coloc_all), "\n")
cat("  Unique GWAS:", length(unique(coloc_all$gwas_name)), "\n")

# ===========================================================================
# 3. Load multi-evidence atlas for Ensembl <-> symbol mapping
# ===========================================================================
cat("\n=== Loading atlas for gene ID mapping ===\n")
atlas <- fread(file.path(BASE,
  "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
  select = c("ensembl_id", "human_symbol"))
# Create mapping: strip version from ensembl_id for matching
atlas[, ensembl_base := sub("\\.[0-9]+$", "", ensembl_id)]
symbol_map <- atlas[human_symbol != "" & !is.na(human_symbol),
                    .(ensembl_base, symbol = human_symbol)]
# Deduplicate (keep first)
symbol_map <- symbol_map[!duplicated(ensembl_base)]
cat("  Gene ID mapping:", nrow(symbol_map), "Ensembl -> symbol pairs\n")

# ===========================================================================
# 4. Load stage-unique signatures (fibrosis axis only)
# ===========================================================================
cat("\n=== Loading stage-unique signatures ===\n")
stage_uniq <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier",
  "stage_unique_signatures.csv"))
cat("  Total stage-unique genes:", nrow(stage_uniq), "\n")

# Filter to fibrosis axis
stage_fib <- stage_uniq[staging_axis == "Fibrosis"]
stage_fib[, ensembl_base := sub("\\.[0-9]+$", "", gene)]
# Map to symbols
stage_fib <- merge(stage_fib, symbol_map, by = "ensembl_base", all.x = TRUE)
cat("  Fibrosis-axis stage-unique genes:", nrow(stage_fib), "\n")
cat("  Mapped to symbols:", sum(!is.na(stage_fib$symbol)), "\n")
cat("  Per-stage counts:\n")
print(table(stage_fib$unique_stage))

# ===========================================================================
# 5. Load fibrosis OVR dream DEGs
# ===========================================================================
cat("\n=== Loading fibrosis OVR dream DEGs ===\n")
fib_ovr <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier",
  "one_vs_rest_fibrosis_dream.csv"))
fib_ovr[, ensembl_base := sub("\\.[0-9]+$", "", gene)]
# Map to symbols
fib_ovr <- merge(fib_ovr, symbol_map, by = "ensembl_base", all.x = TRUE)
cat("  Total OVR rows:", nrow(fib_ovr), "\n")
cat("  Stages:", paste(sort(unique(fib_ovr$stage_label)), collapse = ", "), "\n")

# Define significant OVR DEGs per stage (Tier 2: padj < 0.05, no LFC filter)
fib_ovr[, is_sig := padj < 0.05]
cat("  Significant DEGs per stage:\n")
print(fib_ovr[is_sig == TRUE, .N, by = stage_label][order(stage_label)])

# ===========================================================================
# 6. Define stage transitions and their DEG sets
# ===========================================================================
cat("\n=== Defining stage-transition DEG sets ===\n")

# Define transitions as contrasts between adjacent stages.
# For onset (F0->F1): genes significant in F1_vs_rest but NOT in F0_vs_rest,
#   OR significant in F0_vs_rest with opposite direction in F1.
# More robustly: genes that show differential expression at each stage boundary.
# We define transition DEGs as genes significantly upregulated at stage N+1
# (i.e., significant in F(N+1)_vs_rest with positive logFC) OR significantly
# downregulated at stage N (i.e., significant in F(N)_vs_rest with negative logFC).
# For simplicity and interpretability, we use the per-stage significant DEGs.

stages <- c("F0_vs_rest", "F1_vs_rest", "F2_vs_rest", "F3_vs_rest", "F4_vs_rest")

# Get significant gene sets per stage (symbols)
sig_by_stage <- list()
for (stg in stages) {
  sig_genes <- fib_ovr[stage_label == stg & is_sig == TRUE & !is.na(symbol),
                       unique(symbol)]
  sig_by_stage[[stg]] <- sig_genes
  cat("  ", stg, ":", length(sig_genes), "significant genes with symbols\n")
}

# Define transition gene sets:
# F0->F1 transition = genes sig at F1 but NOT sig at F0 (newly activated)
# F1->F2 transition = genes sig at F2 but NOT sig at F1
# F2->F3 transition = genes sig at F3 but NOT sig at F2
# F3->F4 transition = genes sig at F4 but NOT sig at F3
transitions <- list(
  "F0_to_F1" = setdiff(sig_by_stage[["F1_vs_rest"]], sig_by_stage[["F0_vs_rest"]]),
  "F1_to_F2" = setdiff(sig_by_stage[["F2_vs_rest"]], sig_by_stage[["F1_vs_rest"]]),
  "F2_to_F3" = setdiff(sig_by_stage[["F3_vs_rest"]], sig_by_stage[["F2_vs_rest"]]),
  "F3_to_F4" = setdiff(sig_by_stage[["F4_vs_rest"]], sig_by_stage[["F3_vs_rest"]])
)

for (tr in names(transitions)) {
  cat("  Transition", tr, ":", length(transitions[[tr]]), "newly significant genes\n")
}

# Also define cumulative sets (all genes significant at each stage)
cat("\n  Also computing cumulative stage DEG sets:\n")
cumulative_by_stage <- list()
for (stg in stages) {
  cumulative_by_stage[[stg]] <- sig_by_stage[[stg]]
  cat("    Cumulative", stg, ":", length(cumulative_by_stage[[stg]]), "\n")
}

# ===========================================================================
# 7. Per-transition Fisher's exact enrichment
# ===========================================================================
cat("\n=== Per-transition Fisher's exact enrichment ===\n")

# Universe: all genes present in both COLOC and OVR data
universe <- intersect(coloc$gene, unique(fib_ovr$symbol[!is.na(fib_ovr$symbol)]))
cat("  Universe size (genes in both COLOC and OVR):", length(universe), "\n")

coloc_sig <- coloc[coloc_best_pp4 > 0.9, gene]
coloc_sig_in_univ <- intersect(coloc_sig, universe)
cat("  COLOC PP.H4>0.9 in universe:", length(coloc_sig_in_univ), "\n")

enrichment_results <- list()

# Test enrichment for each transition
for (tr in names(transitions)) {
  tr_genes <- intersect(transitions[[tr]], universe)

  # 2x2 contingency:
  #           COLOC+   COLOC-
  # Trans+    a        b
  # Trans-    c        d
  a <- length(intersect(tr_genes, coloc_sig_in_univ))
  b <- length(setdiff(tr_genes, coloc_sig_in_univ))
  c <- length(setdiff(coloc_sig_in_univ, tr_genes))
  d <- length(universe) - a - b - c

  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft <- fisher.test(mat, alternative = "two.sided")

  enrichment_results[[tr]] <- data.table(
    transition = tr,
    n_transition_genes = length(tr_genes),
    n_coloc_genes = length(coloc_sig_in_univ),
    n_overlap = a,
    odds_ratio = ft$estimate,
    ci_lower = ft$conf.int[1],
    ci_upper = ft$conf.int[2],
    pvalue = ft$p.value,
    universe_size = length(universe)
  )

  cat(sprintf("  %s: overlap=%d, OR=%.2f [%.2f-%.2f], p=%.4g\n",
              tr, a, ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
}

# Also test cumulative per-stage sets
for (stg in stages) {
  stg_genes <- intersect(sig_by_stage[[stg]], universe)

  a <- length(intersect(stg_genes, coloc_sig_in_univ))
  b <- length(setdiff(stg_genes, coloc_sig_in_univ))
  c <- length(setdiff(coloc_sig_in_univ, stg_genes))
  d <- length(universe) - a - b - c

  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft <- fisher.test(mat, alternative = "two.sided")

  enrichment_results[[stg]] <- data.table(
    transition = stg,
    n_transition_genes = length(stg_genes),
    n_coloc_genes = length(coloc_sig_in_univ),
    n_overlap = a,
    odds_ratio = ft$estimate,
    ci_lower = ft$conf.int[1],
    ci_upper = ft$conf.int[2],
    pvalue = ft$p.value,
    universe_size = length(universe)
  )

  cat(sprintf("  %s: overlap=%d, OR=%.2f [%.2f-%.2f], p=%.4g\n",
              stg, a, ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
}

enrichment_dt <- rbindlist(enrichment_results)
enrichment_dt[, padj := p.adjust(pvalue, method = "BH")]

cat("\n  BH-adjusted results:\n")
print(enrichment_dt[, .(transition, n_overlap, odds_ratio = round(odds_ratio, 2),
                         pvalue = signif(pvalue, 3), padj = signif(padj, 3))])

# ===========================================================================
# 8. Tau x COLOC correlation (Spearman)
# ===========================================================================
cat("\n=== Tau x COLOC correlation ===\n")

# Merge fibrosis stage-unique genes with COLOC PP.H4
stage_fib_coloc <- merge(
  stage_fib[!is.na(symbol), .(symbol, tau, unique_stage, most_specific_stage, direction)],
  coloc[, .(symbol = gene, coloc_best_pp4, coloc_best_gwas,
            coloc_n_gwas_h4_05, coloc_n_gwas_h4_08)],
  by = "symbol"
)
cat("  Stage-unique fibrosis genes with COLOC data:", nrow(stage_fib_coloc), "\n")

if (nrow(stage_fib_coloc) > 10) {
  cor_test <- cor.test(stage_fib_coloc$tau, stage_fib_coloc$coloc_best_pp4,
                       method = "spearman", exact = FALSE)
  cat(sprintf("  Spearman rho = %.4f, p = %.4g (tau vs PP.H4)\n",
              cor_test$estimate, cor_test$p.value))

  # Also test separately for high-tau genes (tau > 0.7)
  high_tau <- stage_fib_coloc[tau > 0.7]
  if (nrow(high_tau) > 10) {
    cor_high <- cor.test(high_tau$tau, high_tau$coloc_best_pp4,
                         method = "spearman", exact = FALSE)
    cat(sprintf("  High tau (>0.7) subset: rho = %.4f, p = %.4g (n=%d)\n",
                cor_high$estimate, cor_high$p.value, nrow(high_tau)))
  }
} else {
  cat("  WARNING: Too few genes for correlation test\n")
}

# ===========================================================================
# 9. Classify genes as onset / progression / pan-stage
# ===========================================================================
cat("\n=== Classifying genes as onset / progression / pan-stage ===\n")

# Build per-gene annotation: which stages is each gene significant at?
gene_stage_profile <- fib_ovr[is_sig == TRUE & !is.na(symbol),
  .(stages_sig = paste(sort(unique(stage_label)), collapse = ";"),
    n_stages_sig = uniqueN(stage_label),
    earliest_stage = min(stage_value),
    latest_stage = max(stage_value),
    max_abs_logfc = max(abs(logFC)),
    best_stage = stage_label[which.max(abs(logFC))],
    best_logfc = logFC[which.max(abs(logFC))]),
  by = symbol]

cat("  Genes with at least 1 significant OVR stage:", nrow(gene_stage_profile), "\n")

# Merge with COLOC
gene_classified <- merge(
  gene_stage_profile,
  coloc[, .(symbol = gene, coloc_best_pp4, coloc_best_gwas,
            coloc_n_gwas_h4_05)],
  by = "symbol"
)
cat("  Genes with stage profile + COLOC data:", nrow(gene_classified), "\n")

# Classification logic:
# - onset: significant at F0 or F1 but NOT at F3 or F4 (early-stage specific)
# - progression: significant at F3 or F4 but NOT at F0 or F1 (late-stage specific)
# - pan-stage: significant at both early (F0/F1) and late (F3/F4) stages
# - mid-stage: significant only at F2 (neither early nor late)
# Only classify COLOC-positive genes (PP.H4 > 0.9)
gene_classified[, has_early := grepl("F0_vs_rest|F1_vs_rest", stages_sig)]
gene_classified[, has_late  := grepl("F3_vs_rest|F4_vs_rest", stages_sig)]
gene_classified[, has_mid   := grepl("F2_vs_rest", stages_sig)]

gene_classified[, progression_class := fifelse(
  has_early & has_late, "pan_stage",
  fifelse(has_early & !has_late, "onset",
  fifelse(!has_early & has_late, "progression",
  fifelse(has_mid, "mid_stage", "unclassified")))
)]

# Mark COLOC status
gene_classified[, is_coloc := coloc_best_pp4 > 0.9]

cat("\n  Progression class distribution (all genes with stage profile):\n")
print(gene_classified[, .N, by = progression_class][order(-N)])

cat("\n  Among COLOC PP.H4 > 0.9 genes:\n")
print(gene_classified[is_coloc == TRUE, .N, by = progression_class][order(-N)])

# ===========================================================================
# 10. Cirrhosis-specific vs enzyme COLOC comparison
# ===========================================================================
cat("\n=== Cirrhosis vs enzyme COLOC comparison ===\n")

# Define GWAS phenotype categories
# ARCHIVED 2026-04-09: Whitfield 2023 (36653562) — provenance unverified, Cirrhosis/HCC duplicate Ghouse cases
# cirrhosis_gwas <- c("2023_36653562_Cirrhosis_EUR", "Ghouse_Cirrhosis")
cirrhosis_gwas <- c("Ghouse_Cirrhosis")
enzyme_gwas    <- c("UKBB_ALT", "UKBB_AST", "UKBB_GGT")
nafld_gwas     <- grep("NAFLD|NASH", unique(coloc_all$gwas_name), value = TRUE)
hcc_gwas       <- grep("HCC", unique(coloc_all$gwas_name), value = TRUE)
pdff_gwas      <- grep("PDFF", unique(coloc_all$gwas_name), value = TRUE)

cat("  Cirrhosis GWAS:", paste(cirrhosis_gwas, collapse = ", "), "\n")
cat("  Enzyme GWAS:", paste(enzyme_gwas, collapse = ", "), "\n")
cat("  NAFLD/NASH GWAS:", length(nafld_gwas), "studies\n")
cat("  HCC GWAS:", length(hcc_gwas), "studies\n")
cat("  PDFF GWAS:", length(pdff_gwas), "studies\n")

# Get max PP.H4 per gene per phenotype category
get_max_pp4 <- function(gwas_names) {
  coloc_all[gwas_name %in% gwas_names,
            .(pp4 = max(PP.H4.abf, na.rm = TRUE)),
            by = .(gene)]
}

cirrhosis_pp4 <- get_max_pp4(cirrhosis_gwas)
setnames(cirrhosis_pp4, "pp4", "cirrhosis_pp4")
enzyme_pp4 <- get_max_pp4(enzyme_gwas)
setnames(enzyme_pp4, "pp4", "enzyme_pp4")
nafld_pp4 <- get_max_pp4(nafld_gwas)
setnames(nafld_pp4, "pp4", "nafld_pp4")
hcc_pp4 <- get_max_pp4(hcc_gwas)
setnames(hcc_pp4, "pp4", "hcc_pp4")
pdff_pp4 <- get_max_pp4(pdff_gwas)
setnames(pdff_pp4, "pp4", "pdff_pp4")

# Merge phenotype-specific PP.H4
pheno_coloc <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE),
                      list(cirrhosis_pp4, enzyme_pp4, nafld_pp4, hcc_pp4, pdff_pp4))

# Replace NA with 0 (gene not tested in that GWAS category)
for (col in c("cirrhosis_pp4", "enzyme_pp4", "nafld_pp4", "hcc_pp4", "pdff_pp4")) {
  pheno_coloc[is.na(get(col)), (col) := 0]
  # Handle -Inf from max of empty set
  pheno_coloc[is.infinite(get(col)), (col) := 0]
}

cat("\n  Phenotype-specific COLOC (PP.H4 > 0.9):\n")
cat("    Cirrhosis:", sum(pheno_coloc$cirrhosis_pp4 > 0.9), "\n")
cat("    Enzyme:", sum(pheno_coloc$enzyme_pp4 > 0.9), "\n")
cat("    NAFLD/NASH:", sum(pheno_coloc$nafld_pp4 > 0.9), "\n")
cat("    HCC:", sum(pheno_coloc$hcc_pp4 > 0.9), "\n")
cat("    PDFF:", sum(pheno_coloc$pdff_pp4 > 0.9), "\n")

# Classify: cirrhosis-specific = cirrhosis PP.H4>0.9 & enzyme PP.H4<0.3
pheno_coloc[, gwas_phenotype_class := fifelse(
  cirrhosis_pp4 > 0.9 & enzyme_pp4 < 0.3, "cirrhosis_specific",
  fifelse(enzyme_pp4 > 0.9 & cirrhosis_pp4 < 0.3, "enzyme_specific",
  fifelse(cirrhosis_pp4 > 0.9 & enzyme_pp4 > 0.9, "shared_cirrhosis_enzyme",
  fifelse(nafld_pp4 > 0.9 | hcc_pp4 > 0.9 | pdff_pp4 > 0.9, "other_phenotype",
  "not_significant")))
)]

cat("\n  GWAS phenotype classification:\n")
print(pheno_coloc[, .N, by = gwas_phenotype_class][order(-N)])

# ===========================================================================
# 11. Cross-reference: progression class x GWAS phenotype
# ===========================================================================
cat("\n=== Cross-referencing progression class x GWAS phenotype ===\n")

gene_full <- merge(
  gene_classified,
  pheno_coloc,
  by.x = "symbol", by.y = "gene",
  all.x = TRUE
)

# For genes with both classification and GWAS phenotype data
cat("  Genes with both progression and phenotype classification:",
    nrow(gene_full[!is.na(gwas_phenotype_class)]), "\n")

# Cross-table for COLOC-positive genes
if (nrow(gene_full[is_coloc == TRUE & !is.na(gwas_phenotype_class)]) > 0) {
  cat("\n  Progression class x GWAS phenotype (COLOC PP.H4>0.9):\n")
  cross_tab <- gene_full[is_coloc == TRUE & !is.na(gwas_phenotype_class),
                         .N, by = .(progression_class, gwas_phenotype_class)]
  cross_tab_wide <- dcast(cross_tab, progression_class ~ gwas_phenotype_class,
                          value.var = "N", fill = 0)
  print(cross_tab_wide)
}

# ===========================================================================
# 12. Sanity checks: known MASLD genes
# ===========================================================================
cat("\n=== Sanity check: known MASLD drug targets ===\n")
known_genes <- c("THRB", "NR1H4", "PPARA", "PNPLA3", "TM6SF2", "HSD17B13",
                 "GCKR", "MARC1", "GPAM")

for (g in known_genes) {
  gc <- gene_classified[symbol == g]
  pc <- pheno_coloc[gene == g]

  if (nrow(gc) > 0) {
    cat(sprintf("  %-12s class=%-13s PP4=%.3f stages=%s\n",
                g, gc$progression_class[1], gc$coloc_best_pp4[1],
                gc$stages_sig[1]))
  } else if (nrow(pc) > 0) {
    cat(sprintf("  %-12s (no stage profile) cirr=%.3f enz=%.3f nafld=%.3f\n",
                g, pc$cirrhosis_pp4[1], pc$enzyme_pp4[1], pc$nafld_pp4[1]))
  } else {
    cat(sprintf("  %-12s NOT FOUND in either dataset\n", g))
  }
}

# ===========================================================================
# 13. Compute progression-causal score
# ===========================================================================
cat("\n=== Computing progression-causal score ===\n")

# progression_causal = PP4 x tau x |best_logFC|
# For genes in stage_fib_coloc (stage-unique fibrosis genes with COLOC)
stage_fib_coloc_score <- merge(
  stage_fib_coloc,
  gene_stage_profile[, .(symbol, max_abs_logfc, best_stage, best_logfc)],
  by = "symbol",
  all.x = TRUE
)

stage_fib_coloc_score[, progression_causal_score :=
  coloc_best_pp4 * tau * fifelse(is.na(max_abs_logfc), 0, max_abs_logfc)]

cat("  Genes with progression causal score:", nrow(stage_fib_coloc_score), "\n")
cat("  Score > 0.1:", sum(stage_fib_coloc_score$progression_causal_score > 0.1, na.rm = TRUE), "\n")
cat("  Score > 0.01:", sum(stage_fib_coloc_score$progression_causal_score > 0.01, na.rm = TRUE), "\n")

# Top progression causal genes
cat("\n  Top 15 by progression causal score:\n")
top_prog <- stage_fib_coloc_score[order(-progression_causal_score)][1:min(15, .N)]
print(top_prog[, .(symbol, tau = round(tau, 3), pp4 = round(coloc_best_pp4, 3),
                    logfc = round(best_logfc, 3),
                    score = round(progression_causal_score, 4),
                    stage = unique_stage, gwas = coloc_best_gwas)])

# ===========================================================================
# 14. Build per-gene transition annotation table
# ===========================================================================
cat("\n=== Building per-gene transition annotation ===\n")

# For each gene in the universe, annotate with: which transitions it participates in,
# its COLOC PP.H4, progression class, and phenotype class
per_gene_transition <- gene_classified[, .(
  symbol,
  n_stages_sig,
  earliest_stage,
  latest_stage,
  max_abs_logfc,
  best_stage,
  best_logfc,
  coloc_best_pp4,
  coloc_best_gwas,
  progression_class,
  is_coloc,
  has_early,
  has_late,
  stages_sig
)]

# Add transition membership
per_gene_transition[, in_F0_F1 := symbol %in% transitions[["F0_to_F1"]]]
per_gene_transition[, in_F1_F2 := symbol %in% transitions[["F1_to_F2"]]]
per_gene_transition[, in_F2_F3 := symbol %in% transitions[["F2_to_F3"]]]
per_gene_transition[, in_F3_F4 := symbol %in% transitions[["F3_to_F4"]]]

# Merge phenotype-specific PP.H4
per_gene_transition <- merge(
  per_gene_transition,
  pheno_coloc,
  by.x = "symbol", by.y = "gene",
  all.x = TRUE
)

cat("  Total genes in transition table:", nrow(per_gene_transition), "\n")
cat("  COLOC+ genes:", sum(per_gene_transition$is_coloc, na.rm = TRUE), "\n")

# ===========================================================================
# 15. Build onset vs progression output table
# ===========================================================================
cat("\n=== Building onset vs progression table ===\n")

onset_prog <- gene_classified[progression_class %in%
  c("onset", "progression", "pan_stage", "mid_stage")]

# Add phenotype data
onset_prog <- merge(
  onset_prog,
  pheno_coloc,
  by.x = "symbol", by.y = "gene",
  all.x = TRUE
)

# Add tau from stage-unique signatures (if available)
tau_map <- stage_fib[!is.na(symbol), .(symbol, tau, unique_stage)]
tau_map <- tau_map[!duplicated(symbol)]
onset_prog <- merge(onset_prog, tau_map, by = "symbol", all.x = TRUE)

# Compute progression_causal_score for this table too
onset_prog[, progression_causal_score :=
  coloc_best_pp4 * fifelse(is.na(tau), 0.5, tau) * max_abs_logfc]

cat("  Total classified genes:", nrow(onset_prog), "\n")
cat("  By class:\n")
print(onset_prog[, .N, by = progression_class][order(-N)])
cat("  COLOC+ by class:\n")
print(onset_prog[is_coloc == TRUE, .N, by = progression_class][order(-N)])

# ===========================================================================
# 16. Save outputs
# ===========================================================================
cat("\n=== Saving outputs ===\n")

# 1. Progression COLOC enrichment
fwrite(enrichment_dt, file.path(outdir, "progression_coloc_enrichment.csv"))
cat("  progression_coloc_enrichment.csv:", nrow(enrichment_dt), "rows\n")

# 2. Per-gene transition annotation
fwrite(per_gene_transition[order(-coloc_best_pp4)],
       file.path(outdir, "progression_coloc_by_transition.csv"))
cat("  progression_coloc_by_transition.csv:", nrow(per_gene_transition), "rows\n")

# 3. Onset vs progression genes
fwrite(onset_prog[order(-progression_causal_score)],
       file.path(outdir, "onset_vs_progression_genes.csv"))
cat("  onset_vs_progression_genes.csv:", nrow(onset_prog), "rows\n")

# ===========================================================================
# 17. Summary
# ===========================================================================
cat("\n=== Summary ===\n")
cat("Universe:", length(universe), "genes\n")
cat("COLOC PP.H4>0.9 in universe:", length(coloc_sig_in_univ), "\n")

# Best enrichment
best_enr <- enrichment_dt[which.min(pvalue)]
cat(sprintf("Strongest enrichment: %s (OR=%.2f, p=%.3g)\n",
            best_enr$transition, best_enr$odds_ratio, best_enr$pvalue))

# Onset vs progression counts (COLOC-positive)
cat("COLOC+ onset genes:", nrow(onset_prog[is_coloc == TRUE & progression_class == "onset"]), "\n")
cat("COLOC+ progression genes:", nrow(onset_prog[is_coloc == TRUE & progression_class == "progression"]), "\n")
cat("COLOC+ pan-stage genes:", nrow(onset_prog[is_coloc == TRUE & progression_class == "pan_stage"]), "\n")

# Cirrhosis vs enzyme
cat("Cirrhosis-specific COLOC:", sum(pheno_coloc$gwas_phenotype_class == "cirrhosis_specific", na.rm = TRUE), "\n")
cat("Enzyme-specific COLOC:", sum(pheno_coloc$gwas_phenotype_class == "enzyme_specific", na.rm = TRUE), "\n")
cat("Shared cirrhosis+enzyme:", sum(pheno_coloc$gwas_phenotype_class == "shared_cirrhosis_enzyme", na.rm = TRUE), "\n")

cat("\nAll outputs in:", outdir, "\n")
cat("Completed:", format(Sys.time()), "\n")
