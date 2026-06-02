#!/usr/bin/env Rscript
# 245_resilient_genetics.R
#
# Healthy-control characterization pipeline, step 6 of 6.
#
# Genetic anchor for the resilient transcriptome: do "resilient-up" DEGs (genes
# higher in resilient than obese-MASLD) overlap with protective MASLD genetic
# loci? Three layers:
#   (i)   Hand-curated 3-gene allowlist {PNPLA3, HSD17B13, MTARC1} — boxplot
#         of bulk expression across 3 groups (vis only).
#   (ii)  Hypergeometric test of resilient-up DEGs vs Bayesian-evidence
#         "Protective-LOF" set (407 genes) from
#         convergence_evidence_genetic_down_coherent.csv.
#   (iii) Signed-Z TWAS concordance: are resilient-up genes enriched for
#         negative TWAS Z (= lower-expression-protective) at MASLD GWAS
#         genes? Computed against twas_multi_gwas_combined.csv.
#
# Inputs:
#   resilient_de_pairwise_GSE126848.csv (within-cohort, primary)
#   resilient_de_pairwise_crosscohort.csv (cross-cohort, sensitivity)
#   convergence_evidence_genetic_down_coherent.csv
#   twas_multi_gwas_combined.csv
#   merged_dge.rds (for boxplot expression)
#
# Output: resilient_genetic_enrichment.csv
#
# Spec: docs/superpowers/specs/2026-04-27-healthy-control-audit-design.md
# Env: rnaseq

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HCDIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity/healthy_control_audit")
DEDIR <- file.path(HCDIR, "resilience_de")
GENDIR <- file.path(HCDIR, "gwas_overlap")
TWAS  <- file.path(BASE, "RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv")
PROT_LOF <- file.path(BASE, "RNA-seq/results/multi_evidence/convergence_evidence_genetic_down_coherent.csv")
ATLAS <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
DGE_RDS <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")

dir.create(GENDIR, recursive = TRUE, showWarnings = FALSE)

# Hand-curated protective allowlist (literature: Ghodsian 2022, Abul-Husn 2018,
# BasuRay 2017 for PNPLA3 I148M-protective discussion; Innes 2023 for HSD17B13;
# Vujkovic 2022 for MTARC1)
HAND_PROTECTIVE <- c("PNPLA3", "HSD17B13", "MTARC1")

# Atlas gene mapping
atlas <- fread(ATLAS, select = c("ensembl_id", "human_symbol"))
atlas[, ensembl_clean := sub("\\..*", "", ensembl_id)]
gene_map <- unique(atlas[, .(ensembl_clean, symbol = human_symbol)])

# Bayesian Protective-LOF gene set
prot_lof <- fread(PROT_LOF)
cat(sprintf("Protective-LOF gene set: %d entries\n", nrow(prot_lof)))
prot_lof_genes <- prot_lof[concordance_state == "Protective-LOF" |
                            (!is.na(concordance_state) & grepl("Protective", concordance_state)),
                            human_symbol]
cat(sprintf("  Distinct symbols: %d\n", length(unique(prot_lof_genes))))

# TWAS multi-GWAS combined
twas <- fread(TWAS)
cat(sprintf("\nTWAS combined: %d rows\n", nrow(twas)))
# Aggregate per gene_name across GWAS by best |zscore|
twas_agg <- twas[!is.na(zscore),
                  .(twas_best_z = zscore[which.max(abs(zscore))],
                    twas_best_gwas = gwas[which.max(abs(zscore))],
                    n_gwas_tested = .N,
                    n_gwas_sig = sum(fdr < 0.05, na.rm = TRUE)),
                  by = .(gene, gene_name)]
cat(sprintf("  Aggregated: %d genes\n", nrow(twas_agg)))

# Load DE results
de_files <- c(
  GSE126848_within = file.path(DEDIR, "resilient_de_pairwise_GSE126848.csv"),
  cross_cohort     = file.path(DEDIR, "resilient_de_pairwise_crosscohort.csv")
)

# Resilient-up gene set: (resilient > obese_MASLD) significantly + same direction in (resilient > lean_healthy)
results <- list()
for (label in names(de_files)) {
  fp <- de_files[[label]]
  if (!file.exists(fp)) {
    cat(sprintf("\n[skip] %s not yet present (Script 242 not finished?)\n", label))
    next
  }
  cat(sprintf("\n=== %s ===\n", label))
  d <- fread(fp)

  # Two contrasts in the file
  d_v_mas  <- d[contrast == "groupresilient - groupobese_MASLD"]
  d_v_lean <- d[contrast == "groupresilient - grouplean_healthy"]

  if (nrow(d_v_mas) == 0 || nrow(d_v_lean) == 0) {
    cat(sprintf("  [warn] missing contrasts in %s\n", label))
    next
  }

  # Map ensembl -> symbol
  d_v_mas[,  ensembl_clean := sub("\\..*", "", gene_id)]
  d_v_lean[, ensembl_clean := sub("\\..*", "", gene_id)]
  d_v_mas <- merge(d_v_mas, gene_map, by = "ensembl_clean", all.x = TRUE)
  d_v_lean <- merge(d_v_lean, gene_map, by = "ensembl_clean", all.x = TRUE)

  # Resilient-up: padj<0.1 in vs-MASLD, logFC>0, and direction-consistent in vs-lean
  setkey(d_v_lean, ensembl_clean)
  d_v_mas[, lean_logFC := d_v_lean[ensembl_clean, logFC]]
  resil_up <- d_v_mas[!is.na(symbol) & symbol != "" &
                       adj.P.Val < 0.1 & logFC > 0 &
                       !is.na(lean_logFC) & lean_logFC >= 0,
                       unique(symbol)]
  cat(sprintf("  Resilient-up genes (FDR<0.1, vs-MASLD positive, vs-lean non-negative): %d\n",
              length(resil_up)))

  if (length(resil_up) < 5) {
    cat("  [skip] Too few resilient-up genes\n")
    next
  }

  # Universe of tested genes
  universe <- intersect(unique(d_v_mas$symbol), unique(c(prot_lof_genes, twas_agg$gene_name)))
  universe <- universe[!is.na(universe) & universe != ""]
  cat(sprintf("  Universe: %d testable genes\n", length(universe)))

  # (i) Hand-curated allowlist hits
  hand_hits <- intersect(HAND_PROTECTIVE, resil_up)
  cat(sprintf("  Hand-curated protective allowlist hits: %d (%s)\n",
              length(hand_hits), paste(hand_hits, collapse = ",")))

  # (ii) Hypergeometric: resilient-up vs Protective-LOF set
  k <- sum(resil_up %in% prot_lof_genes)        # white balls drawn
  m <- sum(prot_lof_genes %in% universe)        # white balls in urn
  N <- length(universe)                          # total balls in urn
  n <- length(intersect(resil_up, universe))    # balls drawn
  if (m > 0 && n > 0 && N > m) {
    hyper_p <- phyper(k - 1, m, N - m, n, lower.tail = FALSE)
  } else {
    hyper_p <- NA
  }
  cat(sprintf("  Hypergeom Protective-LOF: k=%d (intersect), m=%d (set in universe), n=%d (drawn), N=%d\n",
              k, m, n, N))
  cat(sprintf("  Hypergeom p = %.3g\n", hyper_p))

  # (iii) Signed-Z TWAS concordance: resilient-up vs MASLD-GWAS
  # Hypothesis: resilient-up genes have negative TWAS Z (lower-expression = protective)
  # A negative correlation between resilient_logFC and twas_z = consistent protection
  d_v_mas_with_twas <- merge(d_v_mas, twas_agg[, .(symbol = gene_name, twas_best_z)],
                              by = "symbol", all.x = TRUE)
  d_resil_with_twas <- d_v_mas_with_twas[symbol %in% resil_up & !is.na(twas_best_z)]
  cat(sprintf("  Resilient-up with TWAS data: %d genes\n", nrow(d_resil_with_twas)))

  # Concordance-rate sign analysis
  if (nrow(d_resil_with_twas) >= 10) {
    n_neg_z <- sum(d_resil_with_twas$twas_best_z < 0)
    n_pos_z <- sum(d_resil_with_twas$twas_best_z > 0)
    binom_p <- binom.test(n_neg_z, n_neg_z + n_pos_z, p = 0.5)$p.value
    cat(sprintf("  Resilient-up TWAS sign: %d neg / %d pos (binom p=%.3g)\n",
                n_neg_z, n_pos_z, binom_p))
    # Background: all DE-tested genes
    d_bg <- d_v_mas_with_twas[!is.na(twas_best_z)]
    bg_neg <- sum(d_bg$twas_best_z < 0)
    bg_pos <- sum(d_bg$twas_best_z > 0)
    cat(sprintf("  Background TWAS sign: %d neg / %d pos\n", bg_neg, bg_pos))
    # Fisher: enrichment of negative-Z in resilient-up vs background
    ft <- fisher.test(matrix(c(n_neg_z, n_pos_z, bg_neg - n_neg_z, bg_pos - n_pos_z), nrow = 2))
    cat(sprintf("  Fisher OR (neg-Z enrichment) = %.2f, p = %.3g\n",
                ft$estimate, ft$p.value))
    fish_or <- ft$estimate
    fish_p <- ft$p.value
  } else {
    n_neg_z <- NA; n_pos_z <- NA; binom_p <- NA; fish_or <- NA; fish_p <- NA
  }

  results[[label]] <- data.table(
    analysis = label,
    n_resilient_up = length(resil_up),
    hand_curated_hits = paste(hand_hits, collapse = ","),
    hypergeom_protective_lof_p = hyper_p,
    hypergeom_k = k, hypergeom_m = m, hypergeom_n = n, hypergeom_N = N,
    twas_neg_z = n_neg_z, twas_pos_z = n_pos_z, twas_binom_p = binom_p,
    twas_fisher_or = fish_or, twas_fisher_p = fish_p
  )

  # Per-gene table for the figure: resilient-up genes with TWAS + protective annotations
  resil_table <- d_v_mas_with_twas[symbol %in% resil_up,
                                     .(symbol, ensembl_id = gene_id,
                                       resilient_logFC = logFC,
                                       resilient_padj = adj.P.Val,
                                       lean_logFC,
                                       twas_best_z,
                                       in_protective_lof = symbol %in% prot_lof_genes,
                                       in_hand_curated = symbol %in% HAND_PROTECTIVE)]
  resil_table[, analysis := label]
  fwrite(resil_table[order(-resilient_logFC)],
         file.path(GENDIR, sprintf("resilient_up_genes_%s.csv", label)))
}

if (length(results) > 0) {
  out <- rbindlist(results, use.names = TRUE, fill = TRUE)
  fwrite(out, file.path(GENDIR, "resilient_genetic_enrichment.csv"))
  cat("\nWrote: resilient_genetic_enrichment.csv\n")
  print(out)
}

# --- Hand-curated 3-gene expression boxplot data (for panel h.i) ---
cat("\nExtracting expression of PNPLA3 / HSD17B13 / MTARC1 across 3 groups...\n")
defs <- fread(file.path(HCDIR, "controls_definitions.csv"))
dge <- readRDS(DGE_RDS)
# Match symbols to ensembl in counts
ens_targets <- gene_map[symbol %in% HAND_PROTECTIVE]
cat("  Matching ensembl IDs:\n")
print(ens_targets)

# Counts -> log2(CPM+1)
log_cpm <- cpm(dge, log = TRUE, prior.count = 2)
# Find the rows in dge for our targets
gene_rows <- which(sub("\\..*", "", rownames(log_cpm)) %in% ens_targets$ensembl_clean)
if (length(gene_rows) > 0) {
  expr_subset <- log_cpm[gene_rows, , drop = FALSE]
  rownames(expr_subset) <- gene_map[match(sub("\\..*", "", rownames(log_cpm)[gene_rows]),
                                            ensembl_clean), symbol]
  cat(sprintf("  Found %d / %d genes in counts: %s\n",
              length(gene_rows), length(HAND_PROTECTIVE),
              paste(rownames(expr_subset), collapse = ",")))

  # Build long-form for boxplot
  expr_dt <- as.data.table(expr_subset, keep.rownames = "symbol")
  expr_long <- melt(expr_dt, id.vars = "symbol", variable.name = "sample_id",
                     value.name = "log2_cpm")
  defs_3way <- defs[is_resilient | is_lean_healthy | is_obese_MASLD,
                     .(sample_id, dataset, group3 = fcase(
                       is_resilient,    "resilient",
                       is_lean_healthy, "lean_healthy",
                       is_obese_MASLD,  "obese_MASLD"))]
  expr_long <- merge(expr_long, defs_3way, by = "sample_id")
  fwrite(expr_long, file.path(GENDIR, "hand_curated_expression_long.csv"))
  cat(sprintf("  Wrote: hand_curated_expression_long.csv (%d rows)\n", nrow(expr_long)))

  # Quick summary
  smry <- expr_long[, .(mean_log2cpm = mean(log2_cpm),
                          median_log2cpm = median(log2_cpm),
                          n = .N),
                     by = .(symbol, group3)]
  cat("\nMean log2 CPM by group:\n")
  print(dcast(smry, symbol ~ group3, value.var = "mean_log2cpm"))
} else {
  cat("  [warn] None of PNPLA3/HSD17B13/MTARC1 found in counts!\n")
}

cat("\nDone.\n")
