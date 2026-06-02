#!/usr/bin/env Rscript
# 164_tas_permutation_null.R
# ---------------------------------------------------------------------------
# Permutation null testing for the 6 TAS validation tests from Script 147.
#
# For each of the 6 validation tests, run 1,000 permutations of the relevant
# labels/scores (preserving within-stage structure where appropriate) to build
# a null distribution. Compute empirical p-value and 95% null CI.
#
# Key question: are V2 (rho=0.679) and V6 (OR=5.64) above the permutation
# null, confirming they are not artifacts of the data structure?
#
# Input (all under results/progression/ unless noted):
#   - transition_activity_scores.csv   (1,444 samples x 7 TAS + assignments)
#   - transition_subtype_assignments.csv (dominant transition per patient)
#   - transition_subtype_validation.csv  (6 original test results)
#   - results/staging_classifier/modeling_metadata.csv (fibrosis, dataset)
#   - consensus_pseudotime.csv          (pseudotime)
#   - cibersortx_celltype_expression/bayesprism_proportions.csv (cell types)
#   - fate_probabilities.csv            (P(F4), subtype)
#   - transition_genesets.csv           (gene sets for V6)
#
# Output (results/prognosis_v2/):
#   - tas_permutation_null.csv
#
# SLURM: cpu, 8 CPUs, 16G, 48h
# Env:   micromamba activate rnaseq
#
# Usage:
#   sbatch --job-name=stg164_perm \
#          --partition=io --cpus-per-task=8 --mem=16G --time=48:00:00 \
#          --output=logs/164_tas_permutation_null_%j.out \
#          --error=logs/164_tas_permutation_null_%j.err \
#          --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && \
#                  micromamba activate rnaseq && \
#                  cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \
#                  Rscript 164_tas_permutation_null.R'"
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
  library(msigdbr)
})

set.seed(42)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR   <- file.path(INT, "results")
PROGDIR <- file.path(RDIR, "progression")
OUTDIR  <- file.path(RDIR, "prognosis_v2")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(INT, "scripts/logs"), showWarnings = FALSE, recursive = TRUE)

NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
N_PERM <- 1000L

cat("=== 164: TAS Permutation Null ===\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("Cores:", NCORES, "\n")
cat("Permutations:", N_PERM, "\n\n")

# ============================================================
# STEP 1: Load all data
# ============================================================
cat("--- Step 1: Load data ---\n")

tas <- fread(file.path(PROGDIR, "transition_activity_scores.csv"))
cat("  TAS:", nrow(tas), "samples x", ncol(tas), "columns\n")

assign_dt <- fread(file.path(PROGDIR, "transition_subtype_assignments.csv"))
cat("  Assignments:", nrow(assign_dt), "samples\n")

orig_val <- fread(file.path(PROGDIR, "transition_subtype_validation.csv"))
cat("  Original validation:", nrow(orig_val), "tests\n")
cat("  Tests:", paste(orig_val$test, collapse = ", "), "\n")

meta <- fread(file.path(RDIR, "staging_classifier/modeling_metadata.csv"))
cat("  Metadata:", nrow(meta), "samples\n")

ptime <- fread(file.path(PROGDIR, "consensus_pseudotime.csv"))
cat("  Pseudotime:", nrow(ptime), "samples\n")

cellprop <- fread(file.path(PROGDIR, "cibersortx_celltype_expression/bayesprism_proportions.csv"))
cat("  Cell-type proportions:", nrow(cellprop), "samples\n")

fate <- fread(file.path(PROGDIR, "fate_probabilities.csv"))
cat("  Fate probabilities:", nrow(fate), "samples\n")

genesets <- fread(file.path(PROGDIR, "transition_genesets.csv"))
cat("  Gene sets:", nrow(genesets), "genes across", uniqueN(genesets$transition), "transitions\n")

cat("\n")

# ============================================================
# STEP 2: Merge into unified working table
# ============================================================
cat("--- Step 2: Build unified working table ---\n")

# Start from TAS
dt <- copy(tas)

# Merge fibrosis stage and dataset from metadata
meta_sub <- meta[, .(sample_id, fibrosis_stage, dataset)]
dt <- merge(dt, meta_sub, by = "sample_id", all.x = TRUE)

# Merge pseudotime
pt_sub <- ptime[, .(sample_id, pseudotime_consensus)]
dt <- merge(dt, pt_sub, by = "sample_id", all.x = TRUE)

# Merge cell-type proportions (Stellate)
if ("Stellate" %in% names(cellprop)) {
  cp_sub <- cellprop[, .(sample_id, Stellate)]
} else {
  # Try case-insensitive match
  stellate_col <- grep("stellate", names(cellprop), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(stellate_col)) {
    cp_sub <- cellprop[, .(sample_id, Stellate = get(stellate_col))]
  } else {
    stop("Cannot find Stellate column in BayesPrism proportions")
  }
}
dt <- merge(dt, cp_sub, by = "sample_id", all.x = TRUE)

# Merge fate probabilities
fate_sub <- fate[, .(sample_id, fate_prob_F4)]
dt <- merge(dt, fate_sub, by = "sample_id", all.x = TRUE)

# Merge dominant assignments from assign_dt
asgn_sub <- assign_dt[, .(sample_id, fib_dominant, is_quiescent)]
dt <- merge(dt, asgn_sub, by = "sample_id", all.x = TRUE,
            suffixes = c("", ".assign"))

# Use the assignment file's fib_dominant if TAS file's column differs
if ("fib_dominant.assign" %in% names(dt)) {
  # Prefer the assignment file's version
  dt[, fib_dominant := fib_dominant.assign]
  dt[, fib_dominant.assign := NULL]
}

cat("  Unified table:", nrow(dt), "rows x", ncol(dt), "columns\n")
cat("  Fibrosis stages:", paste(sort(unique(dt$fibrosis_stage)), collapse = ", "), "\n")
cat("  Fib dominant:", paste(sort(unique(dt$fib_dominant)), collapse = ", "), "\n\n")

# ============================================================
# STEP 3: Extract observed statistics from original validation
# ============================================================
cat("--- Step 3: Extract observed statistics ---\n")

get_obs <- function(test_name) {
  row <- orig_val[test == test_name]
  if (nrow(row) == 0) {
    cat("  WARNING: test", test_name, "not found in original validation\n")
    return(NA_real_)
  }
  # Use the statistic column (which stores the test-specific metric)
  row$statistic
}

obs_v1 <- get_obs("V1_celltype_concordance")   # Wilcoxon W
obs_v2 <- get_obs("V2_pseudotime_ordering")    # Spearman rho
obs_v3 <- get_obs("V3_fate_probability")       # Wilcoxon W
obs_v4 <- get_obs("V4_within_stage_heterogeneity") # Mean H_norm
obs_v5 <- get_obs("V5_loco_cross_dataset")     # Fraction consistent
obs_v6 <- get_obs("V6_pathway_coherence")      # Fisher's OR

cat(sprintf("  V1 observed (W):    %.4f\n", obs_v1))
cat(sprintf("  V2 observed (rho):  %.4f\n", obs_v2))
cat(sprintf("  V3 observed (W):    %.4f\n", obs_v3))
cat(sprintf("  V4 observed (H):    %.4f\n", obs_v4))
cat(sprintf("  V5 observed (frac): %.4f\n", obs_v5))
cat(sprintf("  V6 observed (OR):   %.4f\n", obs_v6))
cat("\n")

# Also extract the effect_size column (more interpretable for some tests)
get_obs_effect <- function(test_name) {
  row <- orig_val[test == test_name]
  if (nrow(row) == 0) return(NA_real_)
  row$effect_size
}

obs_v1_es <- get_obs_effect("V1_celltype_concordance")
obs_v2_es <- get_obs_effect("V2_pseudotime_ordering")
obs_v3_es <- get_obs_effect("V3_fate_probability")
obs_v4_es <- get_obs_effect("V4_within_stage_heterogeneity")
obs_v5_es <- get_obs_effect("V5_loco_cross_dataset")
obs_v6_es <- get_obs_effect("V6_pathway_coherence")

# ============================================================
# Helper: safe boolean mask
# ============================================================
safe_mask <- function(...) {
  m <- Reduce(`&`, list(...))
  m[is.na(m)] <- FALSE
  m
}

# ============================================================
# STEP 4: Define permutation functions for each test
# ============================================================
cat("--- Step 4: Define permutation functions ---\n")

# ---- V1: Cell-type concordance ----
# Permute fib_dominant labels within each fibrosis stage (among F1/F2 samples).
# Recompute Wilcoxon of Stellate fraction between F1->F2 dominant vs others.
perm_v1 <- function(seed) {
  set.seed(seed)
  # Focus on stages 1 and 2 with non-unassigned dominant labels
  mask_f12 <- safe_mask(
    dt$fibrosis_stage %in% c(1, 2),
    dt$fib_dominant != "unassigned"
  )
  sub <- dt[mask_f12]
  if (nrow(sub) < 10) return(NA_real_)

  # Permute fib_dominant within each fibrosis stage
  sub[, fib_dominant_perm := {
    sample(.SD$fib_dominant)
  }, by = fibrosis_stage]

  is_f1f2 <- sub$fib_dominant_perm == "F1_to_F2"
  is_other <- !is_f1f2
  if (sum(is_f1f2) < 3 || sum(is_other) < 3) return(NA_real_)

  # Return the effect size (median difference) rather than W statistic,
  # as effect size is more comparable across permutations
  median(sub$Stellate[is_f1f2], na.rm = TRUE) -
    median(sub$Stellate[is_other], na.rm = TRUE)
}

# ---- V2: Pseudotime ordering ----
# Permute TAS_F2_to_F3 scores across samples within F2 stage.
# Recompute Spearman(TAS_F2_to_F3, pseudotime) within F2.
perm_v2 <- function(seed) {
  set.seed(seed)
  mask_f2 <- safe_mask(
    dt$fibrosis_stage == 2,
    !is.na(dt$pseudotime_consensus)
  )
  sub <- dt[mask_f2]
  if (nrow(sub) < 10) return(NA_real_)

  # Permute TAS scores within F2
  perm_tas <- sample(sub$TAS_F2_to_F3)
  rho <- cor(perm_tas, sub$pseudotime_consensus,
             method = "spearman", use = "complete.obs")
  return(rho)
}

# ---- V3: Fate probability ----
# Permute fib_dominant labels across F2/F3 samples.
# Recompute mean P(F4) difference between advanced dominant vs others.
perm_v3 <- function(seed) {
  set.seed(seed)
  mask_f23 <- safe_mask(
    dt$fibrosis_stage %in% c(2, 3),
    !is.na(dt$fate_prob_F4),
    dt$fib_dominant != "unassigned"
  )
  sub <- dt[mask_f23]
  if (nrow(sub) < 10) return(NA_real_)

  # Permute dominant labels (no within-stage constraint -- labels shuffled freely)
  sub$fib_dominant_perm <- sample(sub$fib_dominant)

  is_adv <- sub$fib_dominant_perm %in% c("F2_to_F3", "F3_to_F4")
  is_early <- !is_adv
  if (sum(is_adv) < 3 || sum(is_early) < 3) return(NA_real_)

  # Effect size: median P(F4) difference
  median(sub$fate_prob_F4[is_adv], na.rm = TRUE) -
    median(sub$fate_prob_F4[is_early], na.rm = TRUE)
}

# ---- V4: Within-stage heterogeneity ----
# Permute fib_dominant labels ACROSS all stages (destroying stage-adjacency
# constraint). Recompute per-stage entropy, return mean normalized entropy.
# Bug fix: permuting labels WITHIN a stage is a no-op for frequency-based
# entropy because table(sample(x)) == table(x). The correct null hypothesis
# is: "does the real stage-adjacency-constrained assignment produce different
# entropy than random cross-stage assignment?"
perm_v4 <- function(seed) {
  set.seed(seed)
  # Get all non-unassigned samples
  mask_all <- dt$fib_dominant != "unassigned"
  if (sum(mask_all) < 10) return(NA_real_)

  # Permute fib_dominant labels ACROSS all samples (not within stage)
  perm_labels_all <- sample(dt$fib_dominant[mask_all])

  entropy_vals <- c()
  idx_all <- which(mask_all)
  for (fstage in c("0", "1", "2", "3", "4")) {
    stage_mask <- as.character(dt$fibrosis_stage[idx_all]) == fstage
    if (sum(stage_mask) < 5) next
    perm_labels <- perm_labels_all[stage_mask]
    tab <- table(perm_labels)
    probs <- tab / sum(tab)
    probs <- probs[probs > 0]
    H <- -sum(probs * log2(probs))
    H_max <- log2(length(probs))
    H_norm <- ifelse(H_max > 0, H / H_max, 0)
    entropy_vals <- c(entropy_vals, H_norm)
  }
  if (length(entropy_vals) == 0) return(NA_real_)
  mean(entropy_vals, na.rm = TRUE)
}

# ---- V5: Cross-dataset consistency ----
# Permute DATASET labels (not TAS values) across samples, keeping TAS fixed.
# Recompute KS test consistency across the permuted "datasets".
# Bug fix: shuffling TAS across datasets homogenizes distributions, making the
# null MORE consistent than real data (inverted null). The correct null
# hypothesis is: "does the real dataset structure matter?" — so we permute
# which samples belong to which dataset while holding TAS constant.
perm_v5 <- function(seed) {
  set.seed(seed)
  fib_transitions <- c("F0_to_F1", "F1_to_F2")
  tas_cols <- paste0("TAS_", fib_transitions)

  # Permute dataset labels across all samples (TAS stays fixed)
  perm_dt <- copy(dt)
  perm_dt[, dataset_perm := sample(dataset)]

  datasets <- unique(perm_dt$dataset_perm)
  ks_pvals <- c()
  for (tr in tas_cols) {
    for (i in seq_along(datasets)) {
      for (j in seq(i + 1, length(datasets))) {
        if (j > length(datasets)) break
        d1_mask <- perm_dt$dataset_perm == datasets[i]
        d2_mask <- perm_dt$dataset_perm == datasets[j]
        if (sum(d1_mask) < 10 || sum(d2_mask) < 10) next
        ks <- suppressWarnings(ks.test(perm_dt[[tr]][d1_mask], perm_dt[[tr]][d2_mask]))
        ks_pvals <- c(ks_pvals, ks$p.value)
      }
    }
  }
  if (length(ks_pvals) == 0) return(NA_real_)
  mean(ks_pvals > 0.05, na.rm = TRUE)
}

# ---- V6: Pathway coherence ----
# Permute gene-set membership for the F0->F1 transition.
# Recompute Fisher's test for Hallmark enrichment, return best OR.
# Bug fix: random gene sets were drawn from the small transition-geneset
# universe (~1,343 genes) rather than the full expressed gene universe (~24K).
# This inflates overlap by chance and makes the null OR too high.
# Fix: load gene annotation to get all expressed gene symbols and use that as
# both the random-draw pool and the Fisher's test denominator.
cat("  Loading Hallmark gene sets for V6...\n")
hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

immune_hallmarks <- c("HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_INFLAMMATORY_RESPONSE",
                       "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                       "HALLMARK_COMPLEMENT", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
immune_hallmarks <- immune_hallmarks[immune_hallmarks %in% names(hallmark_sets)]

# Get the F0->F1 UP gene set
f0f1_up <- genesets[transition == "F0_to_F1" & direction == "UP", gene_symbol]

# Load full expressed gene universe from annotation + dream results
cat("  Loading full gene universe for V6...\n")
gene_annot <- fread(file.path(RDIR, "gene_annotation/human_ensg_to_symbol.tsv"))
dream_res  <- fread(file.path(RDIR, "integration/dream_results.csv"))
# Map dream ENSG IDs to symbols via annotation
dream_genes <- gene_annot[gene_id %in% dream_res$gene, unique(symbol)]
# Remove empty/NA symbols
dream_genes <- dream_genes[!is.na(dream_genes) & nchar(dream_genes) > 0]
# This is the full expressed gene universe (~24K symbols)
full_universe <- dream_genes
universe_size <- length(full_universe)

cat(sprintf("  F0->F1 UP genes: %d, full universe: %d genes (was %d)\n",
            length(f0f1_up), universe_size, uniqueN(genesets$gene_symbol)))

perm_v6 <- function(seed) {
  set.seed(seed)
  # Draw random gene set from the FULL expressed universe (not just transition genes)
  perm_genes <- sample(full_universe, size = length(f0f1_up), replace = FALSE)

  best_or <- NA_real_
  for (pw in immune_hallmarks) {
    pw_genes <- hallmark_sets[[pw]]
    a <- sum(perm_genes %in% pw_genes)                    # in both
    b <- length(perm_genes) - a                             # in perm set only
    c_val <- sum(pw_genes %in% full_universe) - a          # in pathway only (within universe)
    d <- universe_size - a - b - c_val                      # in neither
    if (a < 0 || b < 0 || c_val < 0 || d < 0) next
    ft <- fisher.test(matrix(c(a, b, c_val, d), 2, 2), alternative = "greater")
    if (is.na(best_or) || ft$estimate > best_or) {
      best_or <- ft$estimate
    }
  }
  return(ifelse(is.null(best_or), NA_real_, best_or))
}

cat("  All permutation functions defined.\n\n")

# ============================================================
# STEP 5: Run permutations
# ============================================================
cat("--- Step 5: Run 1000 permutations per test ---\n")

# Pre-generate seeds for reproducibility
perm_seeds <- sample.int(1e7, N_PERM)

run_perm <- function(perm_fn, test_name) {
  cat(sprintf("  Running %s (%d permutations, %d cores)...\n", test_name, N_PERM, NCORES))
  t0 <- Sys.time()
  null_dist <- unlist(mclapply(perm_seeds, perm_fn, mc.cores = NCORES))
  t1 <- Sys.time()
  elapsed <- round(difftime(t1, t0, units = "secs"), 1)
  n_valid <- sum(!is.na(null_dist))
  cat(sprintf("    Done in %.1f sec, %d valid permutations (%.1f%%)\n",
              elapsed, n_valid, 100 * n_valid / N_PERM))
  cat(sprintf("    Null range: [%.4f, %.4f], mean=%.4f, sd=%.4f\n",
              min(null_dist, na.rm = TRUE), max(null_dist, na.rm = TRUE),
              mean(null_dist, na.rm = TRUE), sd(null_dist, na.rm = TRUE)))
  return(null_dist)
}

null_v1 <- run_perm(perm_v1, "V1_celltype_concordance")
null_v2 <- run_perm(perm_v2, "V2_pseudotime_ordering")
null_v3 <- run_perm(perm_v3, "V3_fate_probability")
null_v4 <- run_perm(perm_v4, "V4_within_stage_heterogeneity")
null_v5 <- run_perm(perm_v5, "V5_loco_cross_dataset")
null_v6 <- run_perm(perm_v6, "V6_pathway_coherence")

cat("\n")

# ============================================================
# STEP 6: Compute empirical p-values and null CIs
# ============================================================
cat("--- Step 6: Compute empirical p-values ---\n")

compute_perm_summary <- function(null_dist, observed, test_name, direction = "greater") {
  # Remove NAs
  null_valid <- null_dist[!is.na(null_dist)]
  n <- length(null_valid)

  if (n == 0) {
    return(data.table(
      test_name = test_name,
      observed_stat = observed,
      perm_pvalue = NA_real_,
      null_ci_lo = NA_real_,
      null_ci_hi = NA_real_,
      null_mean = NA_real_,
      null_sd = NA_real_,
      n_perms = 0L
    ))
  }

  # Empirical p-value: (# permutations with stat >= observed) / n
  # Use (count + 1) / (n + 1) to avoid p=0 (conservative correction)
  if (direction == "greater") {
    count_extreme <- sum(null_valid >= observed)
  } else {
    count_extreme <- sum(null_valid <= observed)
  }
  perm_pval <- (count_extreme + 1) / (n + 1)

  # 95% null CI
  ci_lo <- quantile(null_valid, 0.025, na.rm = TRUE, names = FALSE)
  ci_hi <- quantile(null_valid, 0.975, na.rm = TRUE, names = FALSE)

  data.table(
    test_name = test_name,
    observed_stat = observed,
    perm_pvalue = perm_pval,
    null_ci_lo = ci_lo,
    null_ci_hi = ci_hi,
    null_mean = mean(null_valid),
    null_sd = sd(null_valid),
    n_perms = as.integer(n)
  )
}

# For V1: test is whether effect size (Stellate median diff) is greater than null
# Observed effect_size from original validation
result_v1 <- compute_perm_summary(null_v1, obs_v1_es, "V1_celltype_concordance", "greater")

# For V2: test is whether rho is greater than null
result_v2 <- compute_perm_summary(null_v2, obs_v2_es, "V2_pseudotime_ordering", "greater")

# For V3: test is whether effect size (P(F4) median diff) is greater than null
result_v3 <- compute_perm_summary(null_v3, obs_v3_es, "V3_fate_probability", "greater")

# For V4: test is whether entropy is greater than null
# Higher entropy = more heterogeneity = nontrivial assignment
result_v4 <- compute_perm_summary(null_v4, obs_v4_es, "V4_within_stage_heterogeneity", "greater")

# For V5: test is whether consistency fraction is greater than null
result_v5 <- compute_perm_summary(null_v5, obs_v5_es, "V5_loco_cross_dataset", "greater")

# For V6: test is whether OR is greater than null
result_v6 <- compute_perm_summary(null_v6, obs_v6_es, "V6_pathway_coherence", "greater")

results <- rbindlist(list(result_v1, result_v2, result_v3, result_v4, result_v5, result_v6))

cat("\n=== Permutation Null Summary ===\n")
cat(sprintf("%-30s %12s %12s %20s %12s\n",
            "Test", "Observed", "Perm.pval", "95% Null CI", "Above Null?"))
cat(paste(rep("-", 90), collapse = ""), "\n")
for (i in seq_len(nrow(results))) {
  r <- results[i]
  above <- ifelse(r$perm_pvalue < 0.05, "YES", "no")
  cat(sprintf("%-30s %12.4f %12.4f [%8.4f, %8.4f] %12s\n",
              r$test_name, r$observed_stat, r$perm_pvalue,
              r$null_ci_lo, r$null_ci_hi, above))
}
cat("\n")

# ============================================================
# STEP 7: Save output
# ============================================================
cat("--- Step 7: Save output ---\n")

fwrite(results, file.path(OUTDIR, "tas_permutation_null.csv"))
cat("  Saved:", file.path(OUTDIR, "tas_permutation_null.csv"), "\n")
cat("  Columns:", paste(names(results), collapse = ", "), "\n")
cat("  Rows:", nrow(results), "\n")

# Summary interpretation
cat("\n=== Interpretation ===\n")
for (i in seq_len(nrow(results))) {
  r <- results[i]
  if (is.na(r$perm_pvalue)) {
    cat(sprintf("  %s: INDETERMINATE (no valid permutations)\n", r$test_name))
  } else if (r$perm_pvalue < 0.001) {
    cat(sprintf("  %s: STRONG evidence above null (p<0.001, obs=%.4f >> null [%.4f,%.4f])\n",
                r$test_name, r$observed_stat, r$null_ci_lo, r$null_ci_hi))
  } else if (r$perm_pvalue < 0.05) {
    cat(sprintf("  %s: SIGNIFICANT above null (p=%.4f, obs=%.4f > null [%.4f,%.4f])\n",
                r$test_name, r$perm_pvalue, r$observed_stat, r$null_ci_lo, r$null_ci_hi))
  } else {
    cat(sprintf("  %s: NOT above null (p=%.4f, obs=%.4f within null [%.4f,%.4f])\n",
                r$test_name, r$perm_pvalue, r$observed_stat, r$null_ci_lo, r$null_ci_hi))
  }
}

cat("\n=== 164: COMPLETE ===\n")
cat("Finished:", as.character(Sys.time()), "\n")
