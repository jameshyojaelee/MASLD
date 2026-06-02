#!/usr/bin/env Rscript
# 143_progression_regulon_convergence.R
# ---------------------------------------------------------------------------
# Test whether progression-specific TF regulons overlap with progression
# drug targets and GWAS loci — the progression-specific convergence analysis.
#
# Follows patterns from Script 42 (regulon-drug-GWAS convergence).
# Key question: does TF-Drug-GWAS convergence hold for *progression*
# contrasts (like the 3.5x enrichment observed for disease onset)?
#
# Steps:
#   1. Load progression TF activity results (from Script 142)
#   2. Load drug target genes (DGIdb, LINCS, OpenTargets, CGP reversal)
#   3. Load GWAS loci (COLOC, MR from multi-evidence atlas)
#   4. For each progression-significant TF, test if its DoRothEA regulon
#      members are enriched among drug targets and GWAS genes (Fisher exact)
#   5. Permutation test for overall convergence enrichment
#
# Output: results/progression/
#   progression_regulon_convergence.csv   — per-TF enrichment results
#   progression_convergence_summary.csv   — overall enrichment fold/p-value
#
# Usage: Rscript 143_progression_regulon_convergence.R
# Compute: login node OK (~5 min)
# Requires: data.table, decoupleR (for DoRothEA)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(decoupleR)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
ODIR <- file.path(INT, "results/progression")
dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 143: Progression Regulon-Drug-GWAS Convergence ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# ── 1. Load progression TF activity results ──────────────────────────────
cat("--- Loading progression TF activity results ---\n")

contrasts <- c("c2", "c3", "c5")
tf_results <- list()

for (cc in contrasts) {
  f <- file.path(ODIR, sprintf("progression_tf_activity_%s.csv", cc))
  if (file.exists(f)) {
    dt <- fread(f)
    tf_results[[cc]] <- dt
    n_sig <- dt[wilcox_padj < 0.05, .N]
    cat(sprintf("  %s: %d TFs loaded, %d significant (padj<0.05)\n",
      toupper(cc), nrow(dt), n_sig))
  } else {
    cat(sprintf("  WARNING: %s not found — run Script 142 first\n", f))
  }
}

if (length(tf_results) == 0) {
  cat("ERROR: No TF activity results found. Run Script 142 first.\n")
  quit(save = "no", status = 1)
}

# Collect significant TFs across all contrasts
sig_tfs_all <- unique(unlist(lapply(tf_results, function(dt) {
  dt[wilcox_padj < 0.05, tf]
})))
cat(sprintf("  Unique significant TFs across contrasts: %d\n", length(sig_tfs_all)))

# ── 2. Load DoRothEA regulons ────────────────────────────────────────────
cat("\n--- Loading DoRothEA regulons (confidence A-C) ---\n")
tf_net <- as.data.table(get_dorothea(organism = "human", levels = c("A", "B", "C")))
cat("  DoRothEA:", nrow(tf_net), "interactions,",
    length(unique(tf_net$source)), "TFs\n")

# Build TF -> target map for significant TFs
tf_target_map <- tf_net[source %in% sig_tfs_all, .(tf = source, target)]
cat("  Regulon edges for significant TFs:", nrow(tf_target_map), "\n")

all_regulon_targets <- unique(tf_target_map$target)
all_regulon_genes <- unique(c(all_regulon_targets, sig_tfs_all))
cat("  All regulon-associated genes (targets + TFs):", length(all_regulon_genes), "\n")

# ── 3. Load multi-evidence atlas for GWAS + drug annotations ────────────
cat("\n--- Loading multi-evidence atlas ---\n")

atlas_path <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
if (!file.exists(atlas_path)) {
  stop("Multi-evidence atlas not found: ", atlas_path)
}
atlas <- fread(atlas_path)
cat("  Atlas:", nrow(atlas), "genes\n")

# Universe of gene symbols
universe <- atlas[!is.na(human_symbol) & human_symbol != "", human_symbol]
cat("  Universe size:", length(universe), "\n")

# ── 4. Build GWAS gene set ──────────────────────────────────────────────
cat("\n--- Building GWAS gene set ---\n")

# COLOC genes from atlas (any PP4 > 0.5)
coloc_cols <- intersect(c("broadaway_coloc_pp4", "best_liver_enzyme_pp4",
                           "ukbb_alt_coloc_pp4", "finngen_nafld_coloc_pp4",
                           "finngen_nash_coloc_pp4", "bbj_alt_coloc_pp4",
                           "bbj_ast_coloc_pp4", "bbj_ggt_coloc_pp4"),
                         names(atlas))
coloc_genes <- character(0)
for (cc in coloc_cols) {
  hits <- atlas[!is.na(get(cc)) & get(cc) > 0.5, human_symbol]
  coloc_genes <- union(coloc_genes, hits)
  cat(sprintf("  %s > 0.5: %d genes\n", cc, length(hits)))
}

# MR genes — REMOVED 2026-04-22 (MR ditched from paper)
mr_genes <- character(0)

# TWAS genes
twas_genes <- if ("twas_pval" %in% names(atlas)) {
  atlas[!is.na(twas_pval) & twas_pval < 0.05, human_symbol]
} else { character(0) }
cat(sprintf("  TWAS p<0.05: %d genes\n", length(twas_genes)))

gwas_genes <- unique(c(coloc_genes, mr_genes, twas_genes))
cat(sprintf("  Total GWAS causal genes: %d\n", length(gwas_genes)))

# Also load closest GWAS gene list
gwas_closest_path <- file.path(BASE, "GWAS/closest_gene_list.txt")
if (file.exists(gwas_closest_path)) {
  gwas_closest <- fread(gwas_closest_path, header = FALSE)$V1
  # Remove leading whitespace/counts if present
  gwas_closest <- trimws(gwas_closest)
  gwas_closest <- gwas_closest[gwas_closest != ""]
  gwas_genes <- unique(c(gwas_genes, gwas_closest))
  cat(sprintf("  After adding closest genes: %d total GWAS genes\n", length(gwas_genes)))
}

# ── 5. Build drug target gene set ───────────────────────────────────────
cat("\n--- Building drug target gene set ---\n")

drug_genes <- character(0)

# DGIdb
dgidb_path <- file.path(BASE, "RNA-seq/results/drug_repurposing/dgidb_drug_gene_interactions.csv")
if (file.exists(dgidb_path)) {
  dgidb <- fread(dgidb_path)
  dgidb_genes <- unique(dgidb$gene)
  drug_genes <- union(drug_genes, dgidb_genes)
  cat(sprintf("  DGIdb: %d drug-gene interactions, %d unique genes\n",
    nrow(dgidb), length(dgidb_genes)))
} else {
  cat("  WARNING: DGIdb file not found\n")
}

# LINCS reversal from atlas
lincs_genes <- atlas[!is.na(lincs_reversal) & lincs_reversal == TRUE, human_symbol]
drug_genes <- union(drug_genes, lincs_genes)
cat(sprintf("  LINCS reversal: %d genes\n", length(lincs_genes)))

# OpenTargets
ot_genes <- atlas[!is.na(opentargets_drug) & opentargets_drug == TRUE, human_symbol]
drug_genes <- union(drug_genes, ot_genes)
cat(sprintf("  OpenTargets: %d genes\n", length(ot_genes)))

# CGP reversal hits
cgp_path <- file.path(BASE, "RNA-seq/results/drug_repurposing/cgp_reversal_hits.csv")
if (file.exists(cgp_path)) {
  cgp <- fread(cgp_path)
  # Extract leading edge gene symbols from CGP
  if ("leadingEdge" %in% names(cgp)) {
    # Leading edges are ensembl IDs — map to symbols via atlas
    le_ensembl <- unique(unlist(strsplit(cgp$leadingEdge, ";")))
    le_map <- atlas[sub("\\..*", "", ensembl_id) %in% sub("\\..*", "", le_ensembl),
                     human_symbol]
    drug_genes <- union(drug_genes, le_map)
    cat(sprintf("  CGP leading edge: %d genes mapped\n", length(le_map)))
  }
} else {
  cat("  WARNING: CGP reversal hits not found\n")
}

# Convergent drug targets
conv_path <- file.path(BASE, "RNA-seq/results/drug_repurposing/convergent_drug_targets.csv")
if (file.exists(conv_path)) {
  conv <- fread(conv_path)
  conv_sym <- if ("symbol" %in% names(conv)) conv$symbol else character(0)
  drug_genes <- union(drug_genes, conv_sym)
  cat(sprintf("  Convergent drug targets: %d genes\n", length(conv_sym)))
}

cat(sprintf("  Total drug target genes: %d\n", length(drug_genes)))

# ── 6. Per-TF enrichment tests ──────────────────────────────────────────
cat("\n--- Per-TF enrichment tests (Fisher exact) ---\n")

run_per_tf_enrichment <- function(contrast_label, sig_tf_list) {
  if (length(sig_tf_list) == 0) {
    cat(sprintf("  %s: no significant TFs, skipping\n", contrast_label))
    return(NULL)
  }

  results <- rbindlist(lapply(sig_tf_list, function(tf_name) {
    targets <- tf_net[source == tf_name, target]
    targets_in_universe <- intersect(targets, universe)
    n_targets <- length(targets_in_universe)

    if (n_targets < 5) return(NULL)  # Skip TFs with too few mapped targets

    # GWAS enrichment
    n_gwas <- sum(targets_in_universe %in% gwas_genes)
    gwas_fisher <- fisher.test(matrix(c(
      n_gwas,
      n_targets - n_gwas,
      sum(universe %in% gwas_genes) - n_gwas,
      length(universe) - n_targets - sum(universe %in% gwas_genes) + n_gwas
    ), nrow = 2), alternative = "greater")

    # Drug enrichment
    n_drug <- sum(targets_in_universe %in% drug_genes)
    drug_fisher <- fisher.test(matrix(c(
      n_drug,
      n_targets - n_drug,
      sum(universe %in% drug_genes) - n_drug,
      length(universe) - n_targets - sum(universe %in% drug_genes) + n_drug
    ), nrow = 2), alternative = "greater")

    # Triple: targets in both GWAS and drug
    n_triple <- sum(targets_in_universe %in% gwas_genes &
                    targets_in_universe %in% drug_genes)

    data.table(
      tf = tf_name,
      contrast = contrast_label,
      n_targets = n_targets,
      n_gwas_targets = n_gwas,
      gwas_odds_ratio = gwas_fisher$estimate,
      gwas_pval = gwas_fisher$p.value,
      n_drug_targets = n_drug,
      drug_odds_ratio = drug_fisher$estimate,
      drug_pval = drug_fisher$p.value,
      n_triple_targets = n_triple,
      gwas_target_genes = paste(targets_in_universe[targets_in_universe %in% gwas_genes],
                                 collapse = ";"),
      drug_target_genes = paste(targets_in_universe[targets_in_universe %in% drug_genes],
                                 collapse = ";")
    )
  }))

  if (nrow(results) > 0) {
    results[, gwas_padj := p.adjust(gwas_pval, method = "BH")]
    results[, drug_padj := p.adjust(drug_pval, method = "BH")]
    setorder(results, gwas_pval)
  }
  results
}

# Run per contrast
convergence_results <- list()
for (cc in names(tf_results)) {
  sig_list <- tf_results[[cc]][wilcox_padj < 0.05, tf]
  label <- toupper(cc)
  res <- run_per_tf_enrichment(label, sig_list)
  if (!is.null(res) && nrow(res) > 0) {
    convergence_results[[cc]] <- res
    n_gwas_sig <- res[gwas_padj < 0.1, .N]
    n_drug_sig <- res[drug_padj < 0.1, .N]
    cat(sprintf("  %s: %d TFs tested, %d GWAS-enriched, %d drug-enriched (padj<0.1)\n",
      label, nrow(res), n_gwas_sig, n_drug_sig))
  }
}

# Combine all contrast results
all_convergence <- rbindlist(convergence_results)

if (nrow(all_convergence) > 0) {
  fwrite(all_convergence, file.path(ODIR, "progression_regulon_convergence.csv"))
  cat(sprintf("\n  Total entries: %d\n", nrow(all_convergence)))
}

# ── 7. Permutation test for overall convergence ────────────────────────
cat("\n--- Permutation test: progression regulon convergence ---\n")

n_perm <- 10000

# Observed: how many progression-significant TF regulon targets are GWAS genes
observed_gwas <- sum(all_regulon_genes %in% gwas_genes)
observed_drug <- sum(all_regulon_genes %in% drug_genes)
observed_triple <- sum(all_regulon_genes %in% gwas_genes &
                        all_regulon_genes %in% drug_genes)

perm_gwas <- numeric(n_perm)
perm_drug <- numeric(n_perm)
perm_triple <- numeric(n_perm)

n_regulon <- length(all_regulon_genes)

for (i in seq_len(n_perm)) {
  random_genes <- sample(universe, min(n_regulon, length(universe)))
  perm_gwas[i] <- sum(random_genes %in% gwas_genes)
  perm_drug[i] <- sum(random_genes %in% drug_genes)
  perm_triple[i] <- sum(random_genes %in% gwas_genes &
                          random_genes %in% drug_genes)
}

# Compute p-values and fold enrichment
p_gwas <- (sum(perm_gwas >= observed_gwas) + 1) / (n_perm + 1)
p_drug <- (sum(perm_drug >= observed_drug) + 1) / (n_perm + 1)
p_triple <- (sum(perm_triple >= observed_triple) + 1) / (n_perm + 1)

perm_summary <- data.table(
  test = c("GWAS_enrichment", "Drug_enrichment", "Triple_convergence"),
  context = "progression",
  observed = c(observed_gwas, observed_drug, observed_triple),
  expected_mean = c(mean(perm_gwas), mean(perm_drug), mean(perm_triple)),
  expected_sd = c(sd(perm_gwas), sd(perm_drug), sd(perm_triple)),
  fold_enrichment = c(
    observed_gwas / max(mean(perm_gwas), 0.001),
    observed_drug / max(mean(perm_drug), 0.001),
    observed_triple / max(mean(perm_triple), 0.001)
  ),
  perm_pvalue = c(p_gwas, p_drug, p_triple),
  n_permutations = n_perm,
  n_regulon_genes = n_regulon,
  n_gwas_universe = length(gwas_genes),
  n_drug_universe = length(drug_genes),
  n_universe = length(universe)
)

cat("\nPermutation results (progression regulon convergence):\n")
print(perm_summary[, .(test, observed, expected = round(expected_mean, 1),
  fold = round(fold_enrichment, 2), p = signif(perm_pvalue, 3))])

# Compare with onset (3.5x enrichment from Script 42)
cat("\n  Comparison with onset convergence (Script 42): 3.5x GWAS enrichment (p=0.014)\n")
cat(sprintf("  Progression GWAS enrichment: %.2fx (p=%s)\n",
  perm_summary[test == "GWAS_enrichment", fold_enrichment],
  signif(perm_summary[test == "GWAS_enrichment", perm_pvalue], 3)))

# ── 8. Per-contrast convergence summary ─────────────────────────────────
cat("\n--- Per-contrast convergence summary ---\n")

contrast_summaries <- rbindlist(lapply(names(tf_results), function(cc) {
  dt <- tf_results[[cc]]
  sig_tfs <- dt[wilcox_padj < 0.05, tf]

  # Get regulon targets for this contrast's sig TFs
  targets <- unique(tf_net[source %in% sig_tfs, target])
  targets_in_universe <- intersect(targets, universe)

  n_gwas <- sum(targets_in_universe %in% gwas_genes)
  n_drug <- sum(targets_in_universe %in% drug_genes)
  n_both <- sum(targets_in_universe %in% gwas_genes &
                 targets_in_universe %in% drug_genes)

  # Expected by chance
  exp_gwas <- length(targets_in_universe) * length(gwas_genes) / length(universe)
  exp_drug <- length(targets_in_universe) * length(drug_genes) / length(universe)

  data.table(
    contrast = toupper(cc),
    n_sig_tfs = length(sig_tfs),
    n_regulon_targets = length(targets_in_universe),
    n_gwas_overlap = n_gwas,
    expected_gwas = round(exp_gwas, 1),
    gwas_fold = round(n_gwas / max(exp_gwas, 0.001), 2),
    n_drug_overlap = n_drug,
    expected_drug = round(exp_drug, 1),
    drug_fold = round(n_drug / max(exp_drug, 0.001), 2),
    n_triple = n_both
  )
}))

cat("\nPer-contrast convergence:\n")
print(contrast_summaries)

# ── 9. Write summary output ─────────────────────────────────────────────
cat("\n--- Writing outputs ---\n")

summary_out <- rbind(perm_summary, fill = TRUE)
# Append per-contrast details as extra rows
if (nrow(contrast_summaries) > 0) {
  cs_rows <- rbindlist(lapply(seq_len(nrow(contrast_summaries)), function(i) {
    row <- contrast_summaries[i]
    data.table(
      test = paste0(row$contrast, "_GWAS_enrichment"),
      context = "progression",
      observed = row$n_gwas_overlap,
      expected_mean = row$expected_gwas,
      expected_sd = NA_real_,
      fold_enrichment = row$gwas_fold,
      perm_pvalue = NA_real_,
      n_permutations = NA_integer_,
      n_regulon_genes = row$n_regulon_targets,
      n_gwas_universe = NA_integer_,
      n_drug_universe = NA_integer_,
      n_universe = NA_integer_
    )
  }))
  summary_out <- rbind(summary_out, cs_rows, fill = TRUE)
}

fwrite(summary_out, file.path(ODIR, "progression_convergence_summary.csv"))
cat("  Written: progression_convergence_summary.csv\n")

if (nrow(all_convergence) > 0) {
  cat("  Written: progression_regulon_convergence.csv\n")
}

# ── Summary ──────────────────────────────────────────────────────────────
cat("\n=== RESULTS SUMMARY ===\n")
cat(sprintf("  Unique progression-significant TFs: %d\n", length(sig_tfs_all)))
cat(sprintf("  Regulon-associated genes: %d\n", length(all_regulon_genes)))
cat(sprintf("  GWAS enrichment: %.2fx (p=%s)\n",
  perm_summary[test == "GWAS_enrichment", fold_enrichment],
  signif(perm_summary[test == "GWAS_enrichment", perm_pvalue], 3)))
cat(sprintf("  Drug enrichment: %.2fx (p=%s)\n",
  perm_summary[test == "Drug_enrichment", fold_enrichment],
  signif(perm_summary[test == "Drug_enrichment", perm_pvalue], 3)))
cat(sprintf("  Triple convergence: %.2fx (p=%s)\n",
  perm_summary[test == "Triple_convergence", fold_enrichment],
  signif(perm_summary[test == "Triple_convergence", perm_pvalue], 3)))
cat(sprintf("  (Onset reference: 3.5x GWAS enrichment, p=0.014)\n"))
cat("\nOutputs written to:", ODIR, "\n")
cat("\nDone:", format(Sys.time()), "\n")
