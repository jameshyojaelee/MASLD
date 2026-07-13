#!/usr/bin/env Rscript
# 27_multi_evidence_score.R
# Strategy 7: Multi-evidence gene scoring
# Integrates 7 evidence layers into a single 0-100 composite score per gene

library(data.table)
library(ggplot2)

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUTDIR <- file.path(BASE, "RNA-seq/results/multi_evidence")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Strategy 7: Multi-Evidence Gene Scoring ===\n\n")

# ----------------------------------------------------------------
# Helper: rank-percentile normalization (0-1) within non-zero genes
# Genes with zero/NA/missing evidence stay at 0.
# Non-zero genes are ranked among themselves to span (0, 1].
# This prevents sparse layers (97% zeros) from collapsing to 0.
# ----------------------------------------------------------------
rank_norm_nonzero <- function(x) {
  out <- rep(0, length(x))
  nonzero <- !is.na(x) & x != 0
  n_nz <- sum(nonzero)
  if (n_nz > 1) {
    out[nonzero] <- rank(x[nonzero], ties.method = "average") / n_nz
  } else if (n_nz == 1) {
    out[nonzero] <- 1.0
  }
  out
}

# Legacy full-rank normalization (for reference)
rank_norm <- function(x) {
  out <- rep(0, length(x))
  valid <- !is.na(x)
  if (sum(valid) > 1) {
    out[valid] <- rank(x[valid], ties.method = "average") / sum(valid)
  }
  out
}

# ================================================================
# Layer 1: Human DE evidence (meta-analysis)
# ================================================================
cat("Loading Layer 1: Human consensus DEGs...\n")
consensus <- fread(file.path(RDIR, "consensus_degs.csv"))
consensus[, ensembl_clean := sub("\\..*", "", gene)]

# Use the ortholog comparison to get HGNC symbols
ortho <- fread(file.path(RDIR, "human_mouse_ortholog_comparison.csv"))
ortho[, ensembl_clean := sub("\\..*", "", gene_base)]
symbol_map <- unique(ortho[!is.na(human_symbol) & human_symbol != "", .(ensembl_clean, human_symbol)])
symbol_map <- symbol_map[!duplicated(ensembl_clean)]

consensus <- merge(consensus, symbol_map, by = "ensembl_clean", all.x = TRUE)

# Human DE score: -log10(bulk_padj) * sign(bulk_logFC), capped at 50
consensus[, human_de_raw := ifelse(!is.na(bulk_padj) & bulk_padj > 0,
                                    -log10(bulk_padj) * sign(bulk_logFC), 0)]
consensus[human_de_raw > 50, human_de_raw := 50]
consensus[human_de_raw < -50, human_de_raw := -50]

cat("  Genes with HGNC symbol:", sum(!is.na(consensus$human_symbol)), "/", nrow(consensus), "\n")

# ================================================================
# Layer 2: Mouse DE evidence
# ================================================================
cat("Loading Layer 2: Mouse consensus DEGs...\n")
mouse <- fread(file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/mouse_consensus_degs.csv"))
mouse[, ensembl_clean := sub("\\..*", "", gene)]

# Map mouse Ensembl to human symbol via ortholog comparison
mouse_map <- unique(ortho[!is.na(human_symbol) & human_symbol != "",
                           .(mouse_gene_id, human_symbol)])
mouse_map <- mouse_map[!duplicated(mouse_gene_id)]

mouse <- merge(mouse, mouse_map, by.x = "ensembl_clean", by.y = "mouse_gene_id", all.x = TRUE)
mouse[, mouse_de_raw := ifelse(!is.na(dream_padj) & dream_padj > 0,  # C2-OK-sensitivity: MOUSE pipeline (mouse_consensus_degs.csv) still uses dream_* — not migrated
                                -log10(dream_padj) * sign(dream_logFC), 0)]  # C2-OK-sensitivity: mouse consensus columns
mouse[mouse_de_raw > 50, mouse_de_raw := 50]
mouse[mouse_de_raw < -50, mouse_de_raw := -50]

mouse_layer <- mouse[!is.na(human_symbol), .(human_symbol,
                                               mouse_de_raw,
                                               mouse_tier = tier)]
mouse_layer <- mouse_layer[!duplicated(human_symbol)]
cat("  Mouse genes with human symbol:", nrow(mouse_layer), "\n")

# ================================================================
# Layer 3: Cross-species concordance
# ================================================================
cat("Loading Layer 3: Cross-species concordance atlas...\n")
concordance <- fread(file.path(BASE, "Analysis/Cross_Species_Concordance/results/concordance_atlas_unified.csv"))
conc_layer <- concordance[, .(human_symbol, n_concordant, primary_category)]
conc_layer <- conc_layer[!duplicated(human_symbol)]
cat("  Genes in concordance atlas:", nrow(conc_layer), "\n")

# ================================================================
# Layer 4: Causal inference (multi-GWAS replication)
# ================================================================
cat("Loading Layer 4: Causal inference (multi-GWAS)...\n")

# TWAS replication across GWAS
twas_combined_file <- file.path(BASE, "RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv")
if (file.exists(twas_combined_file)) {
  twas_combined <- fread(twas_combined_file)
  # Count how many GWAS each gene is nominally significant in
  twas_rep <- twas_combined[!is.na(pvalue) & pvalue < 0.05, .(
    n_gwas_twas   = .N,
    best_twas_p   = min(pvalue, na.rm = TRUE),
    best_twas_z   = {idx <- which.min(pvalue); if (length(idx) > 0) zscore[idx] else NA_real_}
  ), by = gene]
  # Map Ensembl gene IDs to HGNC symbols (S-PrediXcan writes Ensembl in 'gene', symbol in 'gene_name')
  twas_rep[, gene_clean := sub("\\..*", "", gene)]
  if ("gene_name" %in% names(twas_combined)) {
    sym_lookup <- unique(twas_combined[, .(gene_clean = sub("\\..*", "", gene), gene_name)])
    twas_rep <- merge(twas_rep, sym_lookup, by = "gene_clean", all.x = TRUE)
    twas_rep[, human_symbol := gene_name]
  } else {
    # Fallback: use symbol_map via Ensembl IDs
    twas_rep <- merge(twas_rep, symbol_map[, .(ensembl_clean, human_symbol)],
                      by.x = "gene_clean", by.y = "ensembl_clean", all.x = TRUE)
  }
  cat("  TWAS genes nominally significant in >=1 GWAS:", nrow(twas_rep), "\n")
} else {
  cat("  Multi-GWAS TWAS file not found, using empty table\n")
  twas_rep <- data.table(gene = character(0), human_symbol = character(0),
                          n_gwas_twas = integer(0),
                          best_twas_p = numeric(0), best_twas_z = numeric(0))
}

# MR layer REMOVED 2026-04-22 — MR ditched from paper.
# TWAS-only causal score (MR contribution previously gave +2 for mr_sig).
# Keyed on human_symbol (HGNC).
causal <- twas_rep[, .(human_symbol, n_gwas_twas, best_twas_p, best_twas_z)]
causal[is.na(n_gwas_twas), n_gwas_twas := 0L]
causal[, causal_score := (n_gwas_twas) + ifelse(is.na(best_twas_z), 0, abs(best_twas_z)/10)]
causal[is.na(causal_score), causal_score := 0]

# Build causal_layer for merging (expects 'human_symbol' key)
causal_layer <- causal[, .(human_symbol, causal_score, n_gwas_twas, best_twas_p, best_twas_z)]
causal_layer <- causal_layer[!duplicated(human_symbol)]
cat("  Genes with any causal evidence:", sum(causal_layer$causal_score > 0, na.rm = TRUE), "\n")

# ================================================================
# Layer 5: Sex-stratified robustness
# ================================================================
cat("Loading Layer 5: Sex-stratified analysis...\n")
sex_class_file <- file.path(RDIR, "sex_deg_classification.csv")
if (file.exists(sex_class_file)) {
  sex_class_dt <- fread(sex_class_file)
  sex_class_dt[, ensembl_clean := sub("\\..*", "", gene)]
  sex_class_dt <- merge(sex_class_dt, symbol_map, by = "ensembl_clean", all.x = TRUE)

  # Robustness: Concordant/Shared=1.0, biased/specific=0.5, Divergent/Not_significant=0
  # Handles both old names (Shared, Female_specific, Male_specific) and new (Concordant, Female_biased, Male_biased)
  sex_class_dt[, sex_robust := fcase(
    sex_class %in% c("Shared", "Concordant"),                                      1.0,
    sex_class %in% c("Male_specific", "Female_specific", "Male_biased", "Female_biased"), 0.5,
    default = 0.0
  )]
  sex_layer <- sex_class_dt[!is.na(human_symbol),
                              .(human_symbol, sex_robust,
                                sex_differential = as.logical(sex_differential),
                                sex_class)]
  sex_layer <- sex_layer[!duplicated(human_symbol)]
  cat("  Genes with sex-stratified data:", nrow(sex_layer), "\n")
  cat("  Sex-differential genes:", sum(sex_layer$sex_differential, na.rm = TRUE), "\n")
} else {
  cat("  Sex classification file not found, using empty table\n")
  sex_layer <- data.table(human_symbol = character(0), sex_robust = numeric(0),
                           sex_differential = logical(0), sex_class = character(0))
}

# ================================================================
# Layer 6: Pathway membership (GSEA leading edges)
# ================================================================
cat("Loading Layer 6: GSEA pathway membership...\n")
gsea <- fread(file.path(RDIR, "gsea_results.csv"))
sig_pathways <- gsea[padj < 0.05]
cat("  Significant pathways:", nrow(sig_pathways), "\n")

# Count how many enriched pathways each gene is a leading edge member of
if (nrow(sig_pathways) > 0) {
  le_list <- strsplit(sig_pathways$leadingEdge, ";")
  le_dt <- data.table(
    ensembl_id = unlist(le_list),
    pathway = rep(sig_pathways$pathway, sapply(le_list, length))
  )
  le_dt[, ensembl_clean := sub("\\..*", "", ensembl_id)]
  pathway_counts <- le_dt[, .(n_pathways = uniqueN(pathway)), by = ensembl_clean]
  pathway_counts <- merge(pathway_counts, symbol_map, by = "ensembl_clean", all.x = TRUE)
  pathway_layer <- pathway_counts[!is.na(human_symbol), .(human_symbol, n_pathways)]
  pathway_layer <- pathway_layer[!duplicated(human_symbol)]
} else {
  pathway_layer <- data.table(human_symbol = character(0), n_pathways = integer(0))
}
cat("  Genes in leading edges:", nrow(pathway_layer), "\n")

# ================================================================
# Layer 7: Essentiality (inverse — non-essential preferred)
# ================================================================
cat("Loading Layer 7: Essentiality...\n")
ess <- fread(file.path(BASE, "Analysis/downstream_analysis/essentiality/results/Unified_Essentiality_Combined.csv"))
ess_layer <- ess[, .(human_symbol = gene_symbol,
                      essentiality_score = Mean_Essentiality_Score)]
ess_layer <- ess_layer[!duplicated(human_symbol)]
cat("  Genes with essentiality data:", nrow(ess_layer), "\n")

# ================================================================
# Merge all layers
# ================================================================
cat("\nMerging layers...\n")

# Start with all unique human symbols from consensus
scored <- unique(consensus[!is.na(human_symbol) & human_symbol != "",
                           .(human_symbol, ensembl_clean, human_de_raw,
                             bulk_logFC, bulk_padj)])
scored <- scored[!duplicated(human_symbol)]

scored <- merge(scored, mouse_layer, by = "human_symbol", all.x = TRUE)
scored <- merge(scored, conc_layer, by = "human_symbol", all.x = TRUE)
scored <- merge(scored, causal_layer, by = "human_symbol", all.x = TRUE)
scored <- merge(scored, sex_layer, by = "human_symbol", all.x = TRUE)
scored <- merge(scored, pathway_layer, by = "human_symbol", all.x = TRUE)
scored <- merge(scored, ess_layer, by = "human_symbol", all.x = TRUE)

cat("Total genes after merge:", nrow(scored), "\n")

# ================================================================
# Normalize each layer to 0-1 via rank percentile WITHIN NON-ZERO genes
# This ensures sparse layers (e.g., L4 with 97.9% zeros) still spread
# their non-zero genes across the full 0-1 range.
# ================================================================
cat("Normalizing layers...\n")

# Layer 1: Human DE (rank by absolute DE strength, preserving direction in raw columns)
scored[, L1_human_de := rank_norm_nonzero(abs(human_de_raw))]

# Layer 2: Mouse DE (rank by absolute DE strength)
scored[, L2_mouse_de := rank_norm_nonzero(abs(mouse_de_raw))]

# Layer 3: Concordance (more concordant diets = better)
scored[, L3_concordance := rank_norm_nonzero(n_concordant)]

# Layer 4: Causal score (downweighted — underpowered eQTLs)
scored[, L4_causal := rank_norm_nonzero(causal_score)]

# Layer 5: Sex robustness (already 0/0.5/1)
scored[is.na(sex_robust), sex_robust := 0]
scored[, L5_sex := sex_robust]

# Layer 6: Pathway count
scored[, L6_pathway := rank_norm_nonzero(n_pathways)]

# Layer 7: Essentiality (missing data imputed to baseline 0; more positive = safe)
scored[, has_essentiality := !is.na(essentiality_score)]
scored[is.na(essentiality_score), essentiality_score := 0]
# Dense ranking preserving ties to smoothly map everything to 0-1
scored[, L7_safety := (rank(essentiality_score, ties.method = "average") - 1) / (.N - 1)]

# ================================================================
# Count active layers per gene (how many layers contribute non-zero evidence)
# ================================================================
layer_cols_all <- c("L1_human_de", "L2_mouse_de", "L3_concordance",
                    "L4_causal", "L5_sex", "L6_pathway", "L7_safety")

# NA contributions = 0
for (col in layer_cols_all) {
  scored[is.na(get(col)), (col) := 0]
}

scored[, layers_active := (L1_human_de > 0) + (L2_mouse_de > 0) +
         (L3_concordance > 0) + (L4_causal > 0) + (L5_sex > 0) +
         (L6_pathway > 0) + as.numeric(has_essentiality)]

cat("Layers active distribution:\n")
print(scored[, .N, by = layers_active][order(layers_active)])

# ================================================================
# Weighted composite score (0-100)
# L4_causal downweighted to 5% (underpowered GTEx N=208)
# Weight redistributed: +5% each to L1, L2, L3, L5
# ================================================================
weights <- c(L1 = 0.30, L2 = 0.20, L3 = 0.20, L4 = 0.05,
             L5 = 0.10, L6 = 0.10, L7 = 0.05)

cat("Weights:", paste(names(weights), weights, sep = "=", collapse = ", "), "\n")
cat("NOTE: L4_causal downweighted from 20% to 5% pending Broadaway eQTL integration\n")

scored[, composite_score := 100 * (
  weights["L1"] * L1_human_de +
  weights["L2"] * L2_mouse_de +
  weights["L3"] * L3_concordance +
  weights["L4"] * L4_causal +
  weights["L5"] * L5_sex +
  weights["L6"] * L6_pathway +
  weights["L7"] * L7_safety
)]

# ================================================================
# Priority-based classification
# Priority 1: Top 1% AND multi-layer (layers_active >= 3)
# Priority 2: Top 5%
# Priority 3: Top 10%
# ================================================================
q99 <- quantile(scored$composite_score, 0.99)
q95 <- quantile(scored$composite_score, 0.95)
q90 <- quantile(scored$composite_score, 0.90)

scored[, priority := fcase(
  composite_score >= q99 & layers_active >= 3, "Priority_1",
  composite_score >= q95, "Priority_2",
  composite_score >= q90, "Priority_3",
  default = "Unranked"
)]

cat("\nPriority distribution:\n")
print(scored[, .N, by = priority][order(priority)])

# ================================================================
# Output ranked table
# ================================================================
scored <- scored[order(-composite_score)]
scored[, rank := .I]

cat("\n=== Top 30 genes by composite score ===\n")
print(head(scored[, .(rank, human_symbol, composite_score, priority, layers_active,
                       L1_human_de, L2_mouse_de, L3_concordance,
                       L4_causal, L5_sex, L6_pathway, L7_safety)], 30))

fwrite(scored, file.path(OUTDIR, "multi_evidence_scored_genes.csv"))
cat("\nSaved:", file.path(OUTDIR, "multi_evidence_scored_genes.csv"), "\n")

# ================================================================
# Layer correlation matrix
# ================================================================
cat("\nComputing layer correlations...\n")
layer_cols <- c("L1_human_de", "L2_mouse_de", "L3_concordance",
                "L4_causal", "L5_sex", "L6_pathway", "L7_safety")
cor_mat <- cor(scored[, ..layer_cols], use = "pairwise.complete.obs")
cat("Layer correlation matrix:\n")
print(round(cor_mat, 3))

# Save correlation heatmap
pdf(file.path(OUTDIR, "layer_correlation_heatmap.pdf"), width = 8, height = 7)
cor_dt <- as.data.table(as.data.frame(as.table(cor_mat)))
setnames(cor_dt, c("Layer1", "Layer2", "Correlation"))
ggplot(cor_dt, aes(x = Layer1, y = Layer2, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Multi-Evidence Layer Correlation Matrix")
dev.off()
cat("Saved: layer_correlation_heatmap.pdf\n")

# ================================================================
# Score distribution
# ================================================================
pdf(file.path(OUTDIR, "score_distribution.pdf"), width = 8, height = 5)
ggplot(scored, aes(x = composite_score)) +
  geom_histogram(bins = 50, fill = "#2166AC", alpha = 0.7) +
  geom_vline(xintercept = quantile(scored$composite_score, 0.95),
             linetype = "dashed", color = "red") +
  annotate("text", x = quantile(scored$composite_score, 0.95) + 1,
           y = Inf, vjust = 2, hjust = 0, label = "Top 5%", color = "red") +
  theme_bw(base_size = 12) +
  labs(title = "Multi-Evidence Composite Score Distribution",
       x = "Composite Score (0-100)", y = "Count")
dev.off()
cat("Saved: score_distribution.pdf\n")

# ================================================================
# Radar chart for top genes
# ================================================================
top_genes <- head(scored, 10)$human_symbol
layer_cols_radar <- c("L1_human_de", "L2_mouse_de", "L3_concordance",
                      "L4_causal", "L5_sex", "L6_pathway", "L7_safety")
radar_dt <- melt(scored[human_symbol %in% top_genes,
                         c("human_symbol", layer_cols_radar), with = FALSE],
                  id.vars = "human_symbol",
                  variable.name = "layer", value.name = "score")

pdf(file.path(OUTDIR, "top10_radar_profiles.pdf"), width = 10, height = 8)
ggplot(radar_dt, aes(x = layer, y = score, group = human_symbol, color = human_symbol)) +
  geom_polygon(fill = NA, linewidth = 0.8) +
  geom_point(size = 1.5) +
  coord_polar() +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(size = 8),
        legend.position = "right") +
  labs(title = "Evidence Profiles: Top 10 Genes", color = "Gene")
dev.off()
cat("Saved: top10_radar_profiles.pdf\n")

# ================================================================
# Summary stats
# ================================================================
cat("\n=== Summary ===\n")
cat("Total scored genes:", nrow(scored), "\n")
cat("Score range:", round(min(scored$composite_score), 2), "-",
    round(max(scored$composite_score), 2), "\n")
cat("Median score:", round(median(scored$composite_score), 2), "\n")
cat("Top 5% threshold:", round(quantile(scored$composite_score, 0.95), 2), "\n")
cat("Top 1% threshold:", round(quantile(scored$composite_score, 0.99), 2), "\n")
cat("Genes above top 5%:", sum(scored$composite_score >= quantile(scored$composite_score, 0.95)), "\n")

# Priority summary
cat("\nPriority summary:\n")
cat("  Priority 1 (top 1%, multi-layer):", sum(scored$priority == "Priority_1"), "\n")
cat("  Priority 2 (top 5%):", sum(scored$priority == "Priority_2"), "\n")
cat("  Priority 3 (top 10%):", sum(scored$priority == "Priority_3"), "\n")
cat("  Mean layers_active in Priority 1:", round(mean(scored[priority == "Priority_1"]$layers_active), 1), "\n")

# Layer coverage
cat("\nLayer coverage (non-zero genes):\n")
n_total <- nrow(scored)
cat(sprintf("  L1 Human DE:   %d (%0.1f%%)\n", sum(scored$L1_human_de > 0), 100*sum(scored$L1_human_de > 0)/n_total))
cat(sprintf("  L2 Mouse DE:   %d (%0.1f%%)\n", sum(scored$L2_mouse_de > 0), 100*sum(scored$L2_mouse_de > 0)/n_total))
cat(sprintf("  L3 Concordance:%d (%0.1f%%)\n", sum(scored$L3_concordance > 0), 100*sum(scored$L3_concordance > 0)/n_total))
cat(sprintf("  L4 Causal:     %d (%0.1f%%)\n", sum(scored$L4_causal > 0), 100*sum(scored$L4_causal > 0)/n_total))
cat(sprintf("  L5 Sex robust: %d (%0.1f%%)\n", sum(scored$L5_sex > 0), 100*sum(scored$L5_sex > 0)/n_total))
cat(sprintf("  L6 Pathway:    %d (%0.1f%%)\n", sum(scored$L6_pathway > 0), 100*sum(scored$L6_pathway > 0)/n_total))
cat(sprintf("  L7 Safety:     %d (%0.1f%%)\n", sum(scored$L7_safety > 0), 100*sum(scored$L7_safety > 0)/n_total))

# Validation: check known MASLD expression genes
known_genes <- c("COL1A1", "ACTA2", "TGFB1", "CIDEC", "FAP", "LGALS3", "FASN")
cat("\nValidation — known MASLD expression genes:\n")
for (g in known_genes) {
  row <- scored[human_symbol == g]
  if (nrow(row) > 0) {
    cat(sprintf("  %s: rank=%d, score=%.1f, priority=%s, layers=%d\n",
                g, row$rank[1], row$composite_score[1], row$priority[1], row$layers_active[1]))
  } else {
    cat(sprintf("  %s: not in scored genes\n", g))
  }
}

cat("\nDone.\n")
