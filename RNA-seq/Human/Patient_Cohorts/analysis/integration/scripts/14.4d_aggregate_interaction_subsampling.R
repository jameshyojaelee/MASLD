#!/usr/bin/env Rscript
# 14.4d_aggregate_interaction_subsampling.R
# ---------------------------------------------------------------------------
# Aggregates results from 14.4c interaction-based subsampling iterations.
#
# Computes:
#   1. Per-fraction summary: Spearman rho (interaction logFC), Jaccard per class,
#      recovery rate of sex-dimorphic genes, class-switching matrix
#   2. Per-gene robustness: fraction of iterations retaining same class
#   3. Comparison to full-model classification from Script 26 v2
#
# Output: .../audit_sensitivity/sex_interaction_subsampling/
#   sex_interaction_subsampling_iter_metrics.csv
#   sex_interaction_subsampling_fraction_summary.csv
#   sex_interaction_subsampling_per_gene.csv
#   sex_interaction_subsampling_class_switching.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ITER_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_interaction_subsampling/iterations"
OUT_DIR  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_interaction_subsampling"
FULL_CLS <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_deg_classification.csv"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load full-model classification (Script 26 v2)
# ---------------------------------------------------------------------------
cat("Loading full-model classification...\n")
full <- fread(FULL_CLS)
cat("  Full model:", nrow(full), "genes\n")
cat("  Full model class distribution:\n")
print(full[, .N, by = sex_class][order(-N)])

SEX_CLASSES <- c("Female_biased", "Male_biased", "Divergent", "Concordant")

# ---------------------------------------------------------------------------
# Discover iteration files
# ---------------------------------------------------------------------------
cat("\nDiscovering iteration files...\n")
cls_files <- list.files(ITER_DIR, pattern = "^sub_class_v2_f\\d+_i\\d+\\.csv$", full.names = TRUE)
int_files <- list.files(ITER_DIR, pattern = "^sub_interaction_f\\d+_i\\d+\\.csv$", full.names = TRUE)

cat("  Classification files:", length(cls_files), "\n")
cat("  Interaction files:", length(int_files), "\n")

if (length(cls_files) == 0) stop("No iteration files found in ", ITER_DIR)

# Parse fraction and iteration from filenames
parse_tags <- function(fpath) {
  bn <- basename(fpath)
  frac <- as.integer(sub(".*_f(\\d+)_i\\d+\\.csv$", "\\1", bn))
  iter <- as.integer(sub(".*_f\\d+_i(\\d+)\\.csv$", "\\1", bn))
  data.table(file = fpath, frac = frac, iter = iter)
}

cls_meta <- rbindlist(lapply(cls_files, parse_tags))
int_meta <- rbindlist(lapply(int_files, parse_tags))

fractions <- sort(unique(cls_meta$frac))
cat("  Fractions found:", paste(fractions, collapse = ", "), "\n")
for (f in fractions) {
  cat(sprintf("    Frac %d%%: %d iterations\n", f, sum(cls_meta$frac == f)))
}

# ---------------------------------------------------------------------------
# Aggregate per-iteration metrics
# ---------------------------------------------------------------------------
cat("\nAggregating per-iteration metrics...\n")

iter_results <- list()

for (i in seq_len(nrow(cls_meta))) {
  cls_file <- cls_meta$file[i]
  frac_val <- cls_meta$frac[i]
  iter_val <- cls_meta$iter[i]

  sub_cls <- fread(cls_file)

  # Merge with full-model classification
  merged <- merge(
    full[, .(gene, full_class = sex_class, full_interaction_padj = interaction_padj,
             full_interaction_logFC = interaction_logFC)],
    sub_cls[, .(gene, sub_class = sex_class, sub_interaction_padj = interaction_padj,
                sub_interaction_logFC = interaction_logFC)],
    by = "gene"
  )

  # 1. Spearman rho of interaction logFC
  rho_int <- cor(merged$full_interaction_logFC, merged$sub_interaction_logFC,
                 method = "spearman", use = "complete.obs")

  # 2. Pearson rho of interaction logFC
  r_int <- cor(merged$full_interaction_logFC, merged$sub_interaction_logFC,
               method = "pearson", use = "complete.obs")

  # 3. Recovery rate: fraction of full-model sex-dimorphic genes still dimorphic
  full_dimorphic <- merged$full_interaction_padj < 0.05
  sub_dimorphic  <- merged$sub_interaction_padj < 0.05
  recovery <- sum(full_dimorphic & sub_dimorphic, na.rm = TRUE) /
              max(sum(full_dimorphic, na.rm = TRUE), 1)

  # 4. Class concordance
  class_match <- sum(merged$full_class == merged$sub_class, na.rm = TRUE) / nrow(merged)

  # 5. Per-class Jaccard
  jaccard_vals <- list()
  for (sc in SEX_CLASSES) {
    in_full <- merged$full_class == sc
    in_sub  <- merged$sub_class == sc
    isect <- sum(in_full & in_sub)
    union <- sum(in_full | in_sub)
    jaccard_vals[[paste0("jaccard_", sc)]] <- if (union > 0) isect / union else NA_real_
  }

  # 6. Direction concordance (among genes dimorphic in both)
  both_dimorphic <- full_dimorphic & sub_dimorphic
  if (sum(both_dimorphic) > 0) {
    dir_concord <- sum(sign(merged$full_interaction_logFC[both_dimorphic]) ==
                       sign(merged$sub_interaction_logFC[both_dimorphic])) /
                   sum(both_dimorphic)
  } else {
    dir_concord <- NA_real_
  }

  # 7. Counts per class in subsample
  n_per_class <- sub_cls[, .N, by = sex_class]

  iter_results[[i]] <- data.table(
    frac = frac_val,
    iter = iter_val,
    spearman_interaction_logFC = round(rho_int, 4),
    pearson_interaction_logFC  = round(r_int, 4),
    recovery_dimorphic         = round(recovery, 4),
    class_concordance          = round(class_match, 4),
    direction_concordance      = round(dir_concord, 4),
    n_dimorphic_sub            = sum(sub_dimorphic, na.rm = TRUE),
    n_Female_biased            = sub_cls[sex_class == "Female_biased", .N],
    n_Male_biased              = sub_cls[sex_class == "Male_biased", .N],
    n_Divergent                = sub_cls[sex_class == "Divergent", .N],
    n_Concordant               = sub_cls[sex_class == "Concordant", .N],
    jaccard_Female_biased      = jaccard_vals[["jaccard_Female_biased"]],
    jaccard_Male_biased        = jaccard_vals[["jaccard_Male_biased"]],
    jaccard_Divergent          = jaccard_vals[["jaccard_Divergent"]],
    jaccard_Concordant         = jaccard_vals[["jaccard_Concordant"]]
  )

  if (i %% 10 == 0 || i == nrow(cls_meta)) {
    cat(sprintf("  Processed %d / %d iterations\n", i, nrow(cls_meta)))
  }
}

iter_dt <- rbindlist(iter_results)
fwrite(iter_dt, file.path(OUT_DIR, "sex_interaction_subsampling_iter_metrics.csv"))
cat("Saved: sex_interaction_subsampling_iter_metrics.csv\n")

# ---------------------------------------------------------------------------
# Fraction summary
# ---------------------------------------------------------------------------
cat("\nComputing per-fraction summaries...\n")

frac_summary <- iter_dt[, .(
  n_iters = .N,
  spearman_mean = round(mean(spearman_interaction_logFC, na.rm = TRUE), 4),
  spearman_sd   = round(sd(spearman_interaction_logFC, na.rm = TRUE), 4),
  pearson_mean  = round(mean(pearson_interaction_logFC, na.rm = TRUE), 4),
  recovery_mean = round(mean(recovery_dimorphic, na.rm = TRUE), 4),
  recovery_sd   = round(sd(recovery_dimorphic, na.rm = TRUE), 4),
  class_concordance_mean = round(mean(class_concordance, na.rm = TRUE), 4),
  direction_concordance_mean = round(mean(direction_concordance, na.rm = TRUE), 4),
  n_dimorphic_mean = round(mean(n_dimorphic_sub)),
  jaccard_Female_biased_mean = round(mean(jaccard_Female_biased, na.rm = TRUE), 4),
  jaccard_Male_biased_mean   = round(mean(jaccard_Male_biased, na.rm = TRUE), 4),
  jaccard_Divergent_mean     = round(mean(jaccard_Divergent, na.rm = TRUE), 4)
), by = frac]

cat("\nFraction summaries:\n")
print(frac_summary)

fwrite(frac_summary, file.path(OUT_DIR, "sex_interaction_subsampling_fraction_summary.csv"))
cat("Saved: sex_interaction_subsampling_fraction_summary.csv\n")

# ---------------------------------------------------------------------------
# Per-gene robustness: selection probability + subsampling intervals
# ---------------------------------------------------------------------------
# Primary metric: Selection probability (Meinshausen & Bühlmann 2010,
# Stability Selection). For each gene, the fraction of 70% subsamples
# where the interaction term achieves padj < 0.05. This directly measures
# how consistently a gene is identified as sex-dimorphic under data
# perturbation and is the standard metric for assessing discovery stability.
#
# Secondary metric: Subsampling interval — the 2.5th-97.5th percentile
# of interaction logFC across subsamples. This is NOT a bootstrap CI
# (subsampling is without replacement at 70%, not with replacement at 100%).
# It characterizes effect size variability under data perturbation; if the
# interval excludes zero, the direction of the sex effect is consistent.
#
# Robustness tiers (based on selection probability):
#   Robust:   selection_prob >= 0.60 (Stability Selection standard)
#   Stable:   selection_prob >= 0.40
#   Moderate: selection_prob >= 0.20
#   Fragile:  selection_prob <  0.20
#
# For Concordant genes, robustness is inverted: how often does the gene
# remain non-significant (padj >= 0.05)?
# ---------------------------------------------------------------------------
cat("\nComputing per-gene robustness (selection probability + subsampling intervals)...\n")

primary_frac <- 70
primary_int_files <- int_meta[frac == primary_frac]
primary_cls_files <- cls_meta[frac == primary_frac]

if (nrow(primary_int_files) > 0) {
  # Load per-iteration interaction results
  int_list <- list()
  cls_list <- list()
  for (i in seq_len(nrow(primary_int_files))) {
    int_list[[i]] <- fread(primary_int_files$file[i])[, .(gene, sub_logFC = logFC, sub_padj = padj)]
    cls_list[[i]] <- fread(primary_cls_files$file[i])[, .(gene, sex_class)]
  }

  all_int <- rbindlist(int_list, idcol = "iteration")
  all_cls <- rbindlist(cls_list, idcol = "iteration")
  n_iters <- length(int_list)

  # Full-model reference
  full_ref <- full[, .(gene,
    full_class = sex_class,
    full_int_logFC = interaction_logFC,
    full_int_padj  = interaction_padj,
    full_dimorphic = (!is.na(interaction_padj) & interaction_padj < 0.05)
  )]

  # Per-gene aggregation: bootstrap CI + selection probability
  gene_stability <- all_int[, .(
    n_iters_present    = .N,
    median_sub_logFC   = median(sub_logFC, na.rm = TRUE),
    mean_sub_logFC     = mean(sub_logFC, na.rm = TRUE),
    sd_sub_logFC       = sd(sub_logFC, na.rm = TRUE),
    subsamp_lower      = quantile(sub_logFC, 0.025, na.rm = TRUE),
    subsamp_upper      = quantile(sub_logFC, 0.975, na.rm = TRUE),
    selection_prob     = sum(sub_padj < 0.05, na.rm = TRUE) / .N,
    pct_same_direction = sum(sign(sub_logFC) == sign(median(sub_logFC, na.rm = TRUE)),
                             na.rm = TRUE) / .N
  ), by = gene]

  # Subsampling interval excludes zero?
  gene_stability[, interval_excludes_zero := (subsamp_lower > 0 | subsamp_upper < 0)]

  # Merge with full-model reference
  gene_stability <- merge(gene_stability, full_ref, by = "gene")

  # Class-switching summary
  cls_agg <- all_cls[, .(
    pct_Female_biased = sum(sex_class == "Female_biased") / .N,
    pct_Male_biased   = sum(sex_class == "Male_biased") / .N,
    pct_Divergent     = sum(sex_class == "Divergent") / .N,
    pct_Concordant    = sum(sex_class == "Concordant") / .N,
    modal_class       = names(sort(table(sex_class), decreasing = TRUE))[1]
  ), by = gene]

  gene_stability <- merge(gene_stability, cls_agg, by = "gene", all.x = TRUE)

  # Robustness tiers based on selection probability
  # For dimorphic genes: how often does interaction stay significant?
  # For concordant genes: how often does it stay non-significant?
  gene_stability[, robustness := fifelse(
    full_dimorphic,
    fcase(
      selection_prob >= 0.60, "Robust",
      selection_prob >= 0.40, "Stable",
      selection_prob >= 0.20, "Moderate",
      default = "Fragile"
    ),
    # Concordant: inverted — robust means consistently non-significant
    fcase(
      (1 - selection_prob) >= 0.90, "Robust",
      (1 - selection_prob) >= 0.70, "Stable",
      (1 - selection_prob) >= 0.50, "Moderate",
      default = "Fragile"
    )
  )]

  # --- Report ---
  cat("\nPer-gene robustness at", primary_frac, "% subsampling (30 iterations):\n")
  cat("  Selection probability thresholds: Robust>=0.60, Stable>=0.40, Moderate>=0.20\n\n")

  for (sc in SEX_CLASSES) {
    sub <- gene_stability[full_class == sc]
    if (nrow(sub) > 0) {
      n_int <- sum(sub$interval_excludes_zero, na.rm = TRUE)
      cat(sprintf("  %-16s N=%5d  Robust=%4.0f%%  Stable=%4.0f%%  Moderate=%4.0f%%  Fragile=%4.0f%%\n",
                  sc, nrow(sub),
                  100 * mean(sub$robustness == "Robust"),
                  100 * mean(sub$robustness == "Stable"),
                  100 * mean(sub$robustness == "Moderate"),
                  100 * mean(sub$robustness == "Fragile")))
      cat(sprintf("    Subsampling interval excludes zero: %d / %d (%.0f%%)\n",
                  n_int, nrow(sub), 100 * n_int / nrow(sub)))
      cat(sprintf("    Mean selection probability: %.1f%%\n",
                  100 * mean(sub$selection_prob)))
      cat(sprintf("    Direction consistency: %.0f%%\n",
                  100 * mean(sub$pct_same_direction)))
    }
  }

  # Highlight top robust genes per dimorphic class
  cat("\n  Top robust sex-dimorphic genes (selection_prob >= 0.60):\n")
  robust_genes <- gene_stability[full_dimorphic == TRUE & selection_prob >= 0.60]
  robust_genes <- robust_genes[order(-selection_prob)]
  if (nrow(robust_genes) > 0) {
    cat(sprintf("    Total: %d genes (%d Female_biased, %d Male_biased, %d Divergent)\n",
                nrow(robust_genes),
                sum(robust_genes$full_class == "Female_biased"),
                sum(robust_genes$full_class == "Male_biased"),
                sum(robust_genes$full_class == "Divergent")))
    cat("    Top 20 by selection probability:\n")
    top20 <- head(robust_genes[, .(gene, full_class, full_int_logFC = round(full_int_logFC, 3),
                                    selection_prob = round(selection_prob, 2),
                                    subsamp_lower = round(subsamp_lower, 3),
                                    subsamp_upper = round(subsamp_upper, 3))], 20)
    print(top20)
  } else {
    cat("    None at >= 0.60 threshold\n")
  }

  fwrite(gene_stability, file.path(OUT_DIR, "sex_interaction_subsampling_per_gene.csv"))
  cat("\nSaved: sex_interaction_subsampling_per_gene.csv\n")
}

# ---------------------------------------------------------------------------
# Class-switching matrix (at primary fraction)
# ---------------------------------------------------------------------------
if (nrow(primary_files) > 0) {
  cat("\nClass-switching matrix (full → subsample mode):\n")

  switch_dt <- gene_stability[, .(gene, full_class, modal_class)]
  switch_mat <- table(switch_dt$full_class, switch_dt$modal_class)
  print(switch_mat)

  switch_long <- as.data.table(switch_mat)
  setnames(switch_long, c("from_full", "to_subsample_mode", "count"))
  fwrite(switch_long, file.path(OUT_DIR, "sex_interaction_subsampling_class_switching.csv"))
  cat("Saved: sex_interaction_subsampling_class_switching.csv\n")
}

cat("\nDone.\n")
