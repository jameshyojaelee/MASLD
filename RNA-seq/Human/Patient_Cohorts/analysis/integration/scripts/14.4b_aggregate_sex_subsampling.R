#!/usr/bin/env Rscript
# 14.4b_aggregate_sex_subsampling.R
# ---------------------------------------------------------------------------
# Aggregate sex subsampling stability + permutation results.
#
# Reads per-iteration CSVs produced by 14.4_sex_subsampling_iter.R and
# per-permutation CSVs from 14.4a_sex_permutation_iter.R, then computes:
#   A. Full-model reference counts
#   B. Per-iteration subsampling metrics (vs full model)
#   C. Fraction-level summary statistics
#   D. Per-gene stability across iterations
#   E. Class-switching transition matrix
#   F. Permutation null distribution + significance
#   G. Console interpretation summary
#
# Outputs (in AUDIT_DIR):
#   sex_subsampling_iter_metrics.csv    — one row per iteration (250 expected)
#   sex_subsampling_fraction_summary.csv — one row per fraction (5 rows)
#   sex_subsampling_per_gene.csv        — one row per gene (~34,453)
#   sex_subsampling_class_switching.csv — transition matrix per fraction
#   sex_permutation_summary.csv         — observed vs null per sex_class
# ---------------------------------------------------------------------------

t0 <- proc.time()
library(data.table)

# ===========================================================================
# Paths
# ===========================================================================
BASE        <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT_RESULTS <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
AUDIT_DIR   <- file.path(BASE, "RNA-seq/results/audit_sensitivity/sex_subsampling")
ITER_DIR    <- file.path(AUDIT_DIR, "iterations")
PERM_DIR    <- file.path(AUDIT_DIR, "permutations")

stopifnot(dir.exists(ITER_DIR))
stopifnot(dir.exists(PERM_DIR))

# ===========================================================================
# Part A: Load full-model references
# ===========================================================================
cat("===== Part A: Loading full-model references =====\n")

ref_male   <- fread(file.path(INT_RESULTS, "sex_stratified_results_male.csv"))
ref_female <- fread(file.path(INT_RESULTS, "sex_stratified_results_female.csv"))
ref_class  <- fread(file.path(INT_RESULTS, "sex_deg_classification.csv"))

cat("  Reference male genes:", nrow(ref_male), "\n")
cat("  Reference female genes:", nrow(ref_female), "\n")
cat("  Reference classified genes:", nrow(ref_class), "\n")

full_counts <- ref_class[, .N, by = sex_class]
cat("  Full-model sex_class distribution:\n")
print(full_counts)

# Pre-compute full-model DEG sets (padj < 0.05 & |LFC| > 0.25)
ref_sig_m <- ref_male[padj < 0.05 & abs(logFC) > 0.25, gene]
ref_sig_f <- ref_female[padj < 0.05 & abs(logFC) > 0.25, gene]
cat("  Full-model DEGs: male =", length(ref_sig_m), ", female =", length(ref_sig_f), "\n")

# ===========================================================================
# Part B: Compute per-iteration subsampling metrics
# ===========================================================================
cat("\n===== Part B: Per-iteration subsampling metrics =====\n")

class_files <- list.files(ITER_DIR, pattern = "^sub_class_f\\d+_i\\d+\\.csv$", full.names = TRUE)
cat("  Found", length(class_files), "classification files\n")

if (length(class_files) == 0) {
  stop("No iteration files found in ", ITER_DIR)
}

# Define sex_class levels for Jaccard computation
SEX_CLASSES <- c("Female_specific", "Male_specific", "Shared", "Divergent")

iter_metrics_list <- vector("list", length(class_files))
n_loaded <- 0
n_failed <- 0

for (fi in seq_along(class_files)) {
  cf <- class_files[fi]

  # Parse fraction and iteration from filename
  bn <- basename(cf)
  frac_match <- regmatches(bn, regexec("sub_class_f(\\d+)_i(\\d+)\\.csv", bn))[[1]]
  if (length(frac_match) < 3) {
    warning("Cannot parse filename: ", bn)
    n_failed <- n_failed + 1
    next
  }
  frac_val <- as.numeric(frac_match[2]) / 100
  iter_val <- as.integer(frac_match[3])
  frac_tag <- frac_match[2]
  iter_tag <- frac_match[3]

  # Load classification file and re-classify with current thresholds
  iter_class <- tryCatch(fread(cf), error = function(e) NULL)
  if (is.null(iter_class) || nrow(iter_class) == 0) {
    warning("Empty or unreadable: ", bn)
    n_failed <- n_failed + 1
    next
  }
  # Re-apply current thresholds (padj < 0.05 & |LFC| > 0.25)
  if (all(c("logFC_M", "padj_M", "logFC_F", "padj_F") %in% names(iter_class))) {
    iter_class[, sig_M := (!is.na(padj_M) & padj_M < 0.05 &
                           !is.na(logFC_M) & abs(logFC_M) > 0.25)]
    iter_class[, sig_F := (!is.na(padj_F) & padj_F < 0.05 &
                           !is.na(logFC_F) & abs(logFC_F) > 0.25)]
    iter_class[, same_sign := !is.na(logFC_M) & !is.na(logFC_F) &
                              sign(logFC_M) == sign(logFC_F)]
    iter_class[, sex_class := fcase(
      sig_M & !sig_F,                     "Male_specific",
      !sig_M & sig_F,                     "Female_specific",
      sig_M & sig_F & same_sign,          "Shared",
      sig_M & sig_F & !same_sign,         "Divergent",
      default = "Not_significant"
    )]
    iter_class[, c("sig_M", "sig_F", "same_sign") := NULL]
  }

  # Load corresponding male and female results
  male_file   <- file.path(ITER_DIR, sprintf("sub_male_f%s_i%s.csv", frac_tag, iter_tag))
  female_file <- file.path(ITER_DIR, sprintf("sub_female_f%s_i%s.csv", frac_tag, iter_tag))

  iter_male   <- tryCatch(fread(male_file),   error = function(e) NULL)
  iter_female <- tryCatch(fread(female_file), error = function(e) NULL)

  if (is.null(iter_male) || is.null(iter_female)) {
    warning("Missing male/female file for f", frac_tag, "_i", iter_tag)
    n_failed <- n_failed + 1
    next
  }

  # --- Per-stratum metrics: MALE ---
  common_m <- intersect(ref_male$gene, iter_male$gene)
  ref_lfc_m  <- ref_male$logFC[match(common_m, ref_male$gene)]
  iter_lfc_m <- iter_male$logFC[match(common_m, iter_male$gene)]

  spearman_rho_male <- cor(ref_lfc_m, iter_lfc_m, method = "spearman", use = "complete.obs")

  iter_sig_m  <- iter_male[padj < 0.05 & abs(logFC) > 0.25, gene]
  jaccard_male <- length(intersect(ref_sig_m, iter_sig_m)) / length(union(ref_sig_m, iter_sig_m))
  recovery_male <- if (length(ref_sig_m) > 0) sum(ref_sig_m %in% iter_sig_m) / length(ref_sig_m) else NA_real_
  n_degs_male <- length(iter_sig_m)

  both_sig_m <- intersect(ref_sig_m, iter_sig_m)
  if (length(both_sig_m) > 0) {
    ref_signs_m  <- sign(ref_male$logFC[match(both_sig_m, ref_male$gene)])
    iter_signs_m <- sign(iter_male$logFC[match(both_sig_m, iter_male$gene)])
    dir_conc_male <- mean(ref_signs_m == iter_signs_m)
  } else {
    dir_conc_male <- NA_real_
  }

  # --- Per-stratum metrics: FEMALE ---
  common_f <- intersect(ref_female$gene, iter_female$gene)
  ref_lfc_f  <- ref_female$logFC[match(common_f, ref_female$gene)]
  iter_lfc_f <- iter_female$logFC[match(common_f, iter_female$gene)]

  spearman_rho_female <- cor(ref_lfc_f, iter_lfc_f, method = "spearman", use = "complete.obs")

  iter_sig_f  <- iter_female[padj < 0.05 & abs(logFC) > 0.25, gene]
  jaccard_female <- length(intersect(ref_sig_f, iter_sig_f)) / length(union(ref_sig_f, iter_sig_f))
  recovery_female <- if (length(ref_sig_f) > 0) sum(ref_sig_f %in% iter_sig_f) / length(ref_sig_f) else NA_real_
  n_degs_female <- length(iter_sig_f)

  both_sig_f <- intersect(ref_sig_f, iter_sig_f)
  if (length(both_sig_f) > 0) {
    ref_signs_f  <- sign(ref_female$logFC[match(both_sig_f, ref_female$gene)])
    iter_signs_f <- sign(iter_female$logFC[match(both_sig_f, iter_female$gene)])
    dir_conc_female <- mean(ref_signs_f == iter_signs_f)
  } else {
    dir_conc_female <- NA_real_
  }

  # --- Per-class Jaccard ---
  jaccard_per_class <- setNames(rep(NA_real_, length(SEX_CLASSES)), paste0("jaccard_", SEX_CLASSES))
  for (cls in SEX_CLASSES) {
    ref_genes  <- ref_class[sex_class == cls, gene]
    iter_genes <- iter_class[sex_class == cls, gene]
    n_union <- length(union(ref_genes, iter_genes))
    if (n_union > 0) {
      jaccard_per_class[paste0("jaccard_", cls)] <- length(intersect(ref_genes, iter_genes)) / n_union
    } else {
      jaccard_per_class[paste0("jaccard_", cls)] <- NA_real_
    }
  }

  # --- Class switches ---
  merged <- merge(
    ref_class[, .(gene, ref_class = sex_class)],
    iter_class[, .(gene, iter_class = sex_class)],
    by = "gene"
  )
  n_class_switches <- sum(merged$ref_class != merged$iter_class)

  # Assemble row
  row <- data.table(
    frac                        = frac_val,
    iter                        = iter_val,
    spearman_rho_male           = spearman_rho_male,
    spearman_rho_female         = spearman_rho_female,
    jaccard_male                = jaccard_male,
    jaccard_female              = jaccard_female,
    recovery_male               = recovery_male,
    recovery_female             = recovery_female,
    direction_concordance_male  = dir_conc_male,
    direction_concordance_female = dir_conc_female,
    n_degs_male                 = n_degs_male,
    n_degs_female               = n_degs_female,
    jaccard_Female_specific     = jaccard_per_class[["jaccard_Female_specific"]],
    jaccard_Male_specific       = jaccard_per_class[["jaccard_Male_specific"]],
    jaccard_Shared              = jaccard_per_class[["jaccard_Shared"]],
    jaccard_Divergent           = jaccard_per_class[["jaccard_Divergent"]],
    n_class_switches            = n_class_switches
  )
  iter_metrics_list[[fi]] <- row
  n_loaded <- n_loaded + 1

  # Progress
  if (n_loaded %% 50 == 0) {
    cat("  Processed", n_loaded, "of", length(class_files), "files\n")
  }
}

iter_metrics <- rbindlist(iter_metrics_list, use.names = TRUE, fill = TRUE)
cat("  Successfully loaded:", n_loaded, "iterations\n")
cat("  Failed/skipped:", n_failed, "iterations\n")

# Expected: 5 fracs x 50 iters = 250
expected <- 250
if (n_loaded < expected) {
  cat("  WARNING: Expected", expected, "iterations, got", n_loaded, "\n")
}

# Save
iter_metrics_file <- file.path(AUDIT_DIR, "sex_subsampling_iter_metrics.csv")
fwrite(iter_metrics[order(frac, iter)], iter_metrics_file)
cat("  Saved:", iter_metrics_file, "\n")

# ===========================================================================
# Part C: Fraction summary
# ===========================================================================
cat("\n===== Part C: Fraction summary =====\n")

# Numeric metric columns (everything except frac and iter)
metric_cols <- setdiff(names(iter_metrics), c("frac", "iter"))

# Compute mean, sd, q2.5, q97.5 per fraction for each metric
frac_summary_list <- list()

for (f in sort(unique(iter_metrics$frac))) {
  sub <- iter_metrics[frac == f]
  row <- data.table(frac = f)
  for (mc in metric_cols) {
    vals <- sub[[mc]]
    vals <- vals[!is.na(vals)]
    row[, paste0(mc, "_mean")  := if (length(vals) > 0) mean(vals) else NA_real_]
    row[, paste0(mc, "_sd")    := if (length(vals) > 1) sd(vals) else NA_real_]
    row[, paste0(mc, "_q2.5")  := if (length(vals) > 0) quantile(vals, 0.025) else NA_real_]
    row[, paste0(mc, "_q97.5") := if (length(vals) > 0) quantile(vals, 0.975) else NA_real_]
  }
  frac_summary_list[[length(frac_summary_list) + 1]] <- row
}

fraction_summary <- rbindlist(frac_summary_list, use.names = TRUE, fill = TRUE)

frac_summary_file <- file.path(AUDIT_DIR, "sex_subsampling_fraction_summary.csv")
fwrite(fraction_summary, frac_summary_file)
cat("  Saved:", frac_summary_file, "\n")
cat("  Fractions:", paste(fraction_summary$frac, collapse = ", "), "\n")

# ===========================================================================
# Part D: Per-gene stability
# ===========================================================================
cat("\n===== Part D: Per-gene stability =====\n")

# For each gene, count how many iterations (per fraction) returned the same
# sex_class as the full model.

# Re-load all classification files into a single long table (gene, sex_class, frac, iter)
cat("  Loading all classification files for per-gene analysis...\n")
all_class_list <- vector("list", length(class_files))

for (fi in seq_along(class_files)) {
  cf <- class_files[fi]
  bn <- basename(cf)
  frac_match <- regmatches(bn, regexec("sub_class_f(\\d+)_i(\\d+)\\.csv", bn))[[1]]
  if (length(frac_match) < 3) next

  frac_val <- as.numeric(frac_match[2]) / 100
  iter_val <- as.integer(frac_match[3])

  dt <- tryCatch(fread(cf, select = c("gene", "sex_class")), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) next

  dt[, frac := frac_val]
  dt[, iter := iter_val]
  all_class_list[[fi]] <- dt

  if (fi %% 50 == 0) cat("  Re-loaded", fi, "of", length(class_files), "files\n")
}

all_class <- rbindlist(all_class_list, use.names = TRUE, fill = TRUE)
cat("  Total gene-iteration rows:", nrow(all_class), "\n")

# Build reference lookup
ref_lookup <- ref_class[, .(gene, full_model_class = sex_class)]

# Merge with reference
all_class <- merge(all_class, ref_lookup, by = "gene", all.x = TRUE)
all_class[, same_class := (sex_class == full_model_class)]

# Count same_class per gene per fraction
gene_frac_counts <- all_class[, .(n_same = sum(same_class, na.rm = TRUE),
                                   n_total = .N),
                               by = .(gene, frac)]

# Pivot to wide: one row per gene, columns for each fraction
fracs <- sort(unique(gene_frac_counts$frac))
frac_pct_tags <- as.integer(fracs * 100)  # e.g. 50, 60, 70, 80, 90

gene_stability <- ref_lookup[, .(gene, full_model_class)]

for (i in seq_along(fracs)) {
  f <- fracs[i]
  tag <- frac_pct_tags[i]
  sub <- gene_frac_counts[frac == f, .(gene, n_same)]
  setnames(sub, "n_same", paste0("n_same_class_", tag))
  gene_stability <- merge(gene_stability, sub, by = "gene", all.x = TRUE)
}

# Fill NA with 0 (genes not present in any iteration at that fraction)
for (tag in frac_pct_tags) {
  col <- paste0("n_same_class_", tag)
  gene_stability[is.na(get(col)), (col) := 0L]
}

# Primary stability metric at 70% fraction
if ("n_same_class_70" %in% names(gene_stability)) {
  # Determine n_iters at 70% from iter_metrics
  n_iters_70 <- nrow(iter_metrics[frac == 0.70])
  if (n_iters_70 == 0) n_iters_70 <- 50  # fallback
  gene_stability[, pct_same_class_70 := n_same_class_70 / n_iters_70]
} else {
  # If 70% fraction not present, use the middle fraction
  mid_tag <- frac_pct_tags[ceiling(length(frac_pct_tags) / 2)]
  n_iters_mid <- nrow(iter_metrics[frac == mid_tag / 100])
  if (n_iters_mid == 0) n_iters_mid <- 50
  gene_stability[, pct_same_class_70 := get(paste0("n_same_class_", mid_tag)) / n_iters_mid]
  cat("  NOTE: 70% fraction not found; using", mid_tag, "% as primary metric\n")
}

# Robustness classification
gene_stability[, robustness := fcase(
  pct_same_class_70 >= 0.90, "Robust",
  pct_same_class_70 >= 0.70, "Stable",
  pct_same_class_70 >= 0.50, "Moderate",
  default = "Fragile"
)]

gene_stability_file <- file.path(AUDIT_DIR, "sex_subsampling_per_gene.csv")
fwrite(gene_stability, gene_stability_file)
cat("  Saved:", gene_stability_file, "\n")
cat("  Total genes:", nrow(gene_stability), "\n")
cat("  Robustness distribution:\n")
print(gene_stability[, .N, by = robustness][order(-N)])

# ===========================================================================
# Part E: Class-switching transition matrix
# ===========================================================================
cat("\n===== Part E: Class-switching transition matrix =====\n")

# For each fraction, compute mean transition counts (from full-model class → subsample class)
ALL_CLASSES <- c("Female_specific", "Male_specific", "Shared", "Divergent", "Not_significant")

switching_list <- list()

for (f in fracs) {
  sub_iters <- all_class[frac == f]
  n_iters_f <- length(unique(sub_iters$iter))

  # Count transitions per iteration, then average
  trans <- sub_iters[, .N, by = .(iter, full_model_class, sex_class)]
  trans_mean <- trans[, .(mean_count = mean(N)), by = .(full_model_class, sex_class)]

  # Compute total per from_class for percentages
  totals <- trans_mean[, .(total = sum(mean_count)), by = full_model_class]
  trans_mean <- merge(trans_mean, totals, by = "full_model_class")
  trans_mean[, mean_pct := mean_count / total]
  trans_mean[, frac := f]

  setnames(trans_mean, c("full_model_class", "sex_class"), c("from_class", "to_class"))
  switching_list[[length(switching_list) + 1]] <- trans_mean[, .(frac, from_class, to_class, mean_count, mean_pct)]
}

switching <- rbindlist(switching_list, use.names = TRUE)
switching_file <- file.path(AUDIT_DIR, "sex_subsampling_class_switching.csv")
fwrite(switching[order(frac, from_class, -mean_count)], switching_file)
cat("  Saved:", switching_file, "\n")

# ===========================================================================
# Part F: Permutation null
# ===========================================================================
cat("\n===== Part F: Permutation null =====\n")

perm_files <- list.files(PERM_DIR, pattern = "^perm_summary_i\\d+\\.csv$", full.names = TRUE)
cat("  Found", length(perm_files), "permutation summary files\n")

if (length(perm_files) == 0) {
  cat("  WARNING: No permutation files found — skipping Part F\n")
  perm_summary <- NULL
} else {
  perm_summaries <- rbindlist(lapply(perm_files, function(f) {
    tryCatch(fread(f), error = function(e) NULL)
  }), use.names = TRUE, fill = TRUE)

  # Remove rows where everything is NA (failed iterations)
  perm_summaries <- perm_summaries[!is.na(n_Female_specific) | !is.na(n_Male_specific)]
  cat("  Valid permutation iterations:", nrow(perm_summaries), "\n")

  if (nrow(perm_summaries) > 0) {
    perm_results <- list()
    for (cls in SEX_CLASSES) {
      col_name <- paste0("n_", cls)
      observed <- full_counts[sex_class == cls, N]
      if (length(observed) == 0) observed <- 0L

      null_vals <- perm_summaries[[col_name]]
      null_vals <- null_vals[!is.na(null_vals)]

      if (length(null_vals) > 0) {
        null_mean <- mean(null_vals)
        null_sd   <- sd(null_vals)
        z_score   <- if (null_sd > 0) (observed - null_mean) / null_sd else NA_real_
        # Two-sided permutation p-value (proportion as extreme or more)
        p_value   <- (sum(null_vals >= observed) + 1) / (length(null_vals) + 1)
      } else {
        null_mean <- NA_real_
        null_sd   <- NA_real_
        z_score   <- NA_real_
        p_value   <- NA_real_
      }

      perm_results[[length(perm_results) + 1]] <- data.table(
        sex_class = cls,
        observed  = observed,
        null_mean = null_mean,
        null_sd   = null_sd,
        z_score   = z_score,
        p_value   = p_value
      )
    }
    perm_summary <- rbindlist(perm_results)
  } else {
    cat("  WARNING: All permutation iterations failed\n")
    perm_summary <- NULL
  }

  if (!is.null(perm_summary)) {
    perm_summary_file <- file.path(AUDIT_DIR, "sex_permutation_summary.csv")
    fwrite(perm_summary, perm_summary_file)
    cat("  Saved:", perm_summary_file, "\n")
  }
}

# ===========================================================================
# Part G: Console interpretation summary
# ===========================================================================
cat("\n\n")
cat("=====================================================\n")
cat("  SEX SUBSAMPLING AGGREGATION — SUMMARY\n")
cat("=====================================================\n\n")

cat("--- Iteration loading ---\n")
cat("  Loaded:", n_loaded, "/", expected, "expected iterations\n")
cat("  Failed:", n_failed, "\n\n")

cat("--- Fraction summary (key metrics at each subsample fraction) ---\n")
for (i in seq_len(nrow(fraction_summary))) {
  f <- fraction_summary$frac[i]
  cat(sprintf("  Frac %.0f%%:  rho_M=%.3f (%.3f-%.3f)  rho_F=%.3f (%.3f-%.3f)  Jacc_M=%.3f  Jacc_F=%.3f  Switches=%d\n",
              f * 100,
              fraction_summary[[paste0("spearman_rho_male_mean")]][i],
              fraction_summary[[paste0("spearman_rho_male_q2.5")]][i],
              fraction_summary[[paste0("spearman_rho_male_q97.5")]][i],
              fraction_summary[[paste0("spearman_rho_female_mean")]][i],
              fraction_summary[[paste0("spearman_rho_female_q2.5")]][i],
              fraction_summary[[paste0("spearman_rho_female_q97.5")]][i],
              fraction_summary[[paste0("jaccard_male_mean")]][i],
              fraction_summary[[paste0("jaccard_female_mean")]][i],
              round(fraction_summary[[paste0("n_class_switches_mean")]][i])))
}

cat("\n--- Gene robustness (at 70% subsample) ---\n")
print(gene_stability[, .N, by = robustness][order(-N)])

cat("\n===== INTERPRETATION =====\n")

# Primary interpretation at 70% fraction
if (0.70 %in% fraction_summary$frac) {
  rho_70_male   <- fraction_summary[frac == 0.70, spearman_rho_male_mean]
  rho_70_female <- fraction_summary[frac == 0.70, spearman_rho_female_mean]

  for (sex in c("Male", "Female")) {
    rho_val <- if (sex == "Male") rho_70_male else rho_70_female
    if (!is.na(rho_val)) {
      level <- if (rho_val > 0.95) {
        "NEGLIGIBLE instability"
      } else if (rho_val >= 0.85) {
        "MODERATE — male control scarcity contributes"
      } else {
        "SUBSTANTIAL — results require cautious interpretation"
      }
      cat(sprintf("  %s stratum (70%% subsample): rho=%.3f -> %s\n", sex, rho_val, level))
    }
  }
} else {
  cat("  70% fraction not available for interpretation\n")
}

cat("\n--- Permutation test ---\n")
if (!is.null(perm_summary)) {
  for (i in seq_len(nrow(perm_summary))) {
    cat(sprintf("  %s: observed=%d, null_mean=%.1f, z=%.2f, p=%.4f\n",
                perm_summary$sex_class[i],
                perm_summary$observed[i],
                perm_summary$null_mean[i],
                perm_summary$z_score[i],
                perm_summary$p_value[i]))
  }
} else {
  cat("  No permutation results available\n")
}

cat("\n--- Output files ---\n")
cat("  ", file.path(AUDIT_DIR, "sex_subsampling_iter_metrics.csv"), "\n")
cat("  ", file.path(AUDIT_DIR, "sex_subsampling_fraction_summary.csv"), "\n")
cat("  ", file.path(AUDIT_DIR, "sex_subsampling_per_gene.csv"), "\n")
cat("  ", file.path(AUDIT_DIR, "sex_subsampling_class_switching.csv"), "\n")
if (!is.null(perm_summary)) {
  cat("  ", file.path(AUDIT_DIR, "sex_permutation_summary.csv"), "\n")
}

elapsed <- (proc.time() - t0)["elapsed"]
cat("\nElapsed time:", round(elapsed / 60, 1), "minutes\n")
cat("Done:", as.character(Sys.time()), "\n")
