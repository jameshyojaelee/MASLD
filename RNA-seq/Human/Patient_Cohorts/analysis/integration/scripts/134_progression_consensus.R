#!/usr/bin/env Rscript
# 134_progression_consensus.R
# ---------------------------------------------------------------------------
# Progression Consensus Matrix: classify genes by their disease-stage
# specificity across ALL progression contrasts.
#
# Loads dream results from:
#   C1:  MASLD vs Control           (Script 05)
#   C2:  NAFL vs NASH               (Script 13)
#   C3:  F3-F4 vs F0-F2             (Script 130)
#   C4:  NAFL vs Control            (Script 130)
#   C5:  NAS >= 5 vs NAS < 5        (Script 130)
#   C6:  Extreme endpoints          (Script 130)
#   C8:  Cirrhosis binary           (Script 130)
#   C9:  F2 inflection point        (Script 130)
#   C7a: Steatosis ordinal          (Script 131)
#   C7b: Inflammation ordinal       (Script 131)
#   C7c: Ballooning ordinal         (Script 131)
#   C11: NASH vs Control            (Script 130)
#   C12: Within-NASH fibrosis prog  (Script 130)
#   C13: NASH vs NAFL fib-adjusted  (Script 130)
#   C17: Fibrosis dose-response     (Script 131)
#
# Builds a gene x contrast consensus matrix with logFC, padj, t-statistic,
# and significance flag per contrast, then classifies each gene:
#   onset_only            — sig in C1/C4/C11, not in any progression (C2/C3/C5/C12/C13)
#   progression_only      — sig in C2/C3/C5/C12/C13, not in onset (C1/C4/C11)
#   both_onset_and_progression — sig in >= 1 onset AND >= 1 progression
#   fibrosis_specific     — sig in C3/C8/C9/C17, not in C2/C5
#   inflammation_specific — sig in C7b, not in C7a/C3
#   extreme_only          — sig in C6 only
#   ubiquitous            — sig in >= 40% of contrasts
#
# Tau specificity index (same formula as Script 115):
#   tau = (N - sum(x_i / x_max)) / (N - 1)
#   where x_i = |t-statistic| per contrast, N = number of contrasts
#
# Output (to results/progression/):
#   progression_consensus_matrix.csv      — wide: gene x contrast stats + class
#   gene_progression_classification.csv   — gene, class, tau, peak, counts
#   progression_classification_summary.csv — count per class
#
# Usage: Rscript 134_progression_consensus.R
# SLURM: cpu, 4 CPU, 32GB RAM, ~15min
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results")
ODIR <- file.path(RDIR, "progression")
dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 134: Progression Consensus Matrix ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ============================================================
# 1. Load all contrast dream results
# ============================================================
cat("=== Loading contrast dream results ===\n")

# Registry: contrast_id -> file path, padj column name
contrast_registry <- list(
  # range = covariate span (1 for binary, N for ordinal 0-N)
  # For ordinal contrasts, total_logFC = logFC * range
  # Primary significance: padj < PADJ_THRESH AND |total_logFC| >= TOTAL_LFC_THRESH
  C1  = list(file = file.path(RDIR, "integration/dream_results.csv"),
    padj_col = "padj", label = "MASLD vs Control", range = 1),
  C2  = list(file = file.path(RDIR, "disease_signatures/nafl_vs_nash_dream.csv"),
    padj_col = "adj.P.Val", label = "NASH vs NAFL", range = 1),
  C3  = list(file = file.path(RDIR, "progression/c3_adv_vs_early_fib_dream.csv"),
    padj_col = "padj", label = "F3-F4 vs F0-F2", range = 1),
  C4  = list(file = file.path(RDIR, "progression/c4_nafl_vs_ctrl_dream.csv"),
    padj_col = "padj", label = "NAFL vs Control", range = 1),
  C5  = list(file = file.path(RDIR, "progression/c5_nas_ge5_vs_lt5_dream.csv"),
    padj_col = "padj", label = "NAS >= 5 vs < 5", range = 1),
  C6  = list(file = file.path(RDIR, "progression/c6_extreme_endpoints_dream.csv"),
    padj_col = "padj", label = "Extreme endpoints", range = 1),
  C8  = list(file = file.path(RDIR, "progression/c8_cirrhosis_dream.csv"),
    padj_col = "padj", label = "Cirrhosis binary", range = 1),
  C9  = list(file = file.path(RDIR, "progression/c9_f2_inflection_dream.csv"),
    padj_col = "padj", label = "F2 inflection", range = 1),
  C7a = list(file = file.path(RDIR, "progression/c7a_steatosis_ordinal_dream.csv"),
    padj_col = "padj", label = "Steatosis ordinal", range = 3),
  C7b = list(file = file.path(RDIR, "progression/c7b_inflammation_ordinal_dream.csv"),
    padj_col = "padj", label = "Inflammation ordinal", range = 3),
  C7c = list(file = file.path(RDIR, "progression/c7c_ballooning_ordinal_dream.csv"),
    padj_col = "padj", label = "Ballooning ordinal", range = 2),
  C11 = list(file = file.path(RDIR, "progression/c11_nash_vs_ctrl_dream.csv"),
    padj_col = "padj", label = "NASH vs Control", range = 1),
  C12 = list(file = file.path(RDIR, "progression/c12_early_vs_late_nash_dream.csv"),
    padj_col = "padj", label = "Within-NASH fibrosis progression", range = 1),
  C13 = list(file = file.path(RDIR, "progression/c13_nash_vs_nafl_fib_adj_dream.csv"),
    padj_col = "padj", label = "NASH vs NAFL fibrosis-adjusted", range = 1),
  C17 = list(file = file.path(RDIR, "progression/c17_fibrosis_ordinal_dream.csv"),
    padj_col = "padj", label = "Fibrosis dose-response ordinal", range = 4)
)

# Load each contrast, harmonising column names
loaded <- list()
for (cid in names(contrast_registry)) {
  entry <- contrast_registry[[cid]]
  if (!file.exists(entry$file)) {
    cat(sprintf("  WARNING: %s file not found — skipping: %s\n", cid, basename(entry$file)))
    next
  }

  dt <- fread(entry$file, select = c("gene", "logFC", "t", entry$padj_col))

  # Harmonise padj column name
  if (entry$padj_col != "padj") {
    setnames(dt, entry$padj_col, "padj")
  }

  # Deduplicate — keep first occurrence per gene (should be unique, but be safe)
  dt <- dt[!duplicated(gene)]

  loaded[[cid]] <- dt
  n_sig <- sum(dt$padj < 0.1, na.rm = TRUE)
  cat(sprintf("  %s (%s): %d genes, %d DEGs (padj<0.1)\n",
    cid, entry$label, nrow(dt), n_sig))
}

n_contrasts <- length(loaded)
contrast_ids <- names(loaded)
cat(sprintf("\nLoaded %d / %d contrasts\n\n", n_contrasts, length(contrast_registry)))

if (n_contrasts < 2) {
  stop("Need at least 2 contrasts to build consensus matrix. Aborting.")
}

# ============================================================
# 2. Gene annotations from C1 (dream_results.csv)
# ============================================================
cat("=== Loading gene annotations ===\n")
# Gene annotations (dream_results.csv lacks symbol/gene_type; use disease signatures)
annot_file <- file.path(RDIR, "disease_signatures/unified_disease_signatures.csv")
if (file.exists(annot_file)) {
  gene_annot <- fread(annot_file, select = c("gene", "symbol", "gene_type"))
} else {
  gene_annot <- fread(file.path(RDIR, "integration/dream_results.csv"),
    select = intersect(c("gene", "symbol", "gene_type"),
      names(fread(file.path(RDIR, "integration/dream_results.csv"), nrows = 0))))
}
gene_annot <- gene_annot[!duplicated(gene)]
cat(sprintf("  Gene annotations: %d genes\n\n", nrow(gene_annot)))

# ============================================================
# 3. Build consensus matrix (wide format: gene x contrast)
# ============================================================
cat("=== Building consensus matrix ===\n")

# Collect all unique genes across all contrasts
all_genes <- unique(unlist(lapply(loaded, function(dt) dt$gene)))
cat(sprintf("  Total unique genes across contrasts: %d\n", length(all_genes)))

# Initialise the wide matrix
consensus <- data.table(gene = all_genes)

for (cid in contrast_ids) {
  dt <- loaded[[cid]]
  suffix <- tolower(cid)

  # Rename columns for this contrast
  dt_sub <- dt[, .(gene, logFC, padj, t)]
  setnames(dt_sub, c("logFC", "padj", "t"),
    c(paste0("logFC_", suffix), paste0("padj_", suffix), paste0("t_", suffix)))

  # Merge into consensus (left join, NAs for missing genes)
  consensus <- merge(consensus, dt_sub, by = "gene", all.x = TRUE)
}

# ---- Total-effect normalization for significance ----
# For ordinal contrasts, logFC is a per-unit slope. To compare across contrast
# types, we define "total logFC" = logFC * range, where range = covariate span
# (1 for binary, 2 for ballooning 0-2, 3 for steatosis/inflammation 0-3, 4 for fibrosis 0-4).
# Primary significance: padj < 0.05 AND |total_logFC| >= 0.5
# Sensitivity: padj < 0.1 (no logFC filter) — used for GSEA ranking, tau computation

PADJ_THRESH       <- 0.05   # Primary padj threshold
TOTAL_LFC_THRESH  <- 0.5    # Minimum total log2FC across full covariate range
PADJ_SENSITIVITY  <- 0.1    # Sensitivity threshold (no logFC filter)

# Build range lookup from registry
contrast_ranges <- setNames(
  sapply(contrast_registry[contrast_ids], `[[`, "range"),
  contrast_ids
)
cat("  Covariate ranges:", paste(sprintf("%s=%d", names(contrast_ranges), contrast_ranges), collapse=", "), "\n")

# Add primary significance flags (padj < 0.05 AND |total_logFC| >= 0.5)
for (cid in contrast_ids) {
  suffix   <- tolower(cid)
  padj_col <- paste0("padj_", suffix)
  lfc_col  <- paste0("logFC_", suffix)
  sig_col  <- paste0("sig_", suffix)
  rng      <- contrast_ranges[cid]

  # Primary: padj + total-effect filter
  consensus[, (sig_col) := fifelse(
    !is.na(get(padj_col)) & get(padj_col) < PADJ_THRESH &
    !is.na(get(lfc_col))  & abs(get(lfc_col)) * rng >= TOTAL_LFC_THRESH,
    1L, 0L)]
}

# Also add sensitivity flags (padj < 0.1 only, for tau and GSEA)
for (cid in contrast_ids) {
  suffix   <- tolower(cid)
  padj_col <- paste0("padj_", suffix)
  sens_col <- paste0("sens_", suffix)
  consensus[, (sens_col) := fifelse(
    !is.na(get(padj_col)) & get(padj_col) < PADJ_SENSITIVITY, 1L, 0L)]
}

# Count significant contrasts per gene (PRIMARY threshold)
sig_cols <- paste0("sig_", tolower(contrast_ids))
consensus[, n_contrasts_sig := rowSums(.SD, na.rm = TRUE), .SDcols = sig_cols]

# Count sensitivity contrasts per gene
sens_cols <- paste0("sens_", tolower(contrast_ids))
consensus[, n_contrasts_sens := rowSums(.SD, na.rm = TRUE), .SDcols = sens_cols]

ubi_thresh <- ceiling(length(contrast_ids) * 0.4)

cat(sprintf("  Consensus matrix: %d genes x %d contrasts\n", nrow(consensus), n_contrasts))
cat(sprintf("  Thresholds: PRIMARY padj<%.2f + |total_logFC|>=%.1f; SENSITIVITY padj<%.1f\n",
  PADJ_THRESH, TOTAL_LFC_THRESH, PADJ_SENSITIVITY))
cat(sprintf("  Genes sig (primary) in 0 contrasts: %d\n", sum(consensus$n_contrasts_sig == 0)))
cat(sprintf("  Genes sig (primary) in >= 1 contrast: %d\n", sum(consensus$n_contrasts_sig >= 1)))
cat(sprintf("  Genes sig (primary) in >= %d contrasts (40%%): %d\n",
  ubi_thresh, sum(consensus$n_contrasts_sig >= ubi_thresh)))
cat(sprintf("  Genes sig (sensitivity) in >= 1 contrast: %d\n", sum(consensus$n_contrasts_sens >= 1)))

# ============================================================
# 4. Compute tau specificity index across contrasts
# ============================================================
cat("\n=== Computing tau specificity index ===\n")

# Build |t-statistic| matrix for tau computation
t_cols <- paste0("t_", tolower(contrast_ids))
t_mat <- as.matrix(consensus[, ..t_cols])
t_mat_abs <- abs(t_mat)
t_mat_abs[is.na(t_mat_abs)] <- 0  # treat missing as 0

N <- ncol(t_mat_abs)
x_max <- apply(t_mat_abs, 1, max)

# tau = (N - sum(x_i / x_max)) / (N - 1)
# For genes with x_max == 0 (never tested or all t=0), tau = 0
tau_vals <- ifelse(x_max == 0, 0,
  (N - rowSums(t_mat_abs / x_max)) / (N - 1))

# Peak contrast = which contrast has highest |t| for each gene
peak_idx <- apply(t_mat_abs, 1, which.max)
peak_contrast <- contrast_ids[peak_idx]

consensus[, tau := tau_vals]
consensus[, peak_contrast := peak_contrast]

cat(sprintf("  tau > 0.8 (highly specific):  %d genes\n", sum(tau_vals > 0.8)))
cat(sprintf("  tau > 0.6 (moderately specific): %d genes\n", sum(tau_vals > 0.6)))
cat(sprintf("  tau < 0.3 (ubiquitous):       %d genes\n", sum(tau_vals < 0.3)))

# ============================================================
# 5. Classify genes into progression categories
# ============================================================
cat("\n=== Classifying genes ===\n")

# Define contrast groups for classification
# Onset contrasts: C1 (MASLD vs Ctrl), C4 (NAFL vs Ctrl), C11 (NASH vs Ctrl)
onset_ids   <- intersect(c("C1", "C4", "C11"), contrast_ids)
# Progression contrasts: C2 (NAFL vs NASH), C3 (Adv vs Early Fib), C5 (NAS>=5),
#   C12 (within-NASH fib progression), C13 (NASH vs NAFL fib-adjusted)
prog_ids    <- intersect(c("C2", "C3", "C5", "C12", "C13"), contrast_ids)
# Fibrosis contrasts: C3, C8 (cirrhosis), C9 (F2 inflection), C17 (fibrosis ordinal)
fib_ids     <- intersect(c("C3", "C8", "C9", "C17"), contrast_ids)
# Inflammation contrast: C7b
inflam_ids  <- intersect(c("C7b"), contrast_ids)
# Non-inflammation fibrosis/steatosis: C7a, C3
non_inflam_ids <- intersect(c("C7a", "C3"), contrast_ids)

# Helper: check if gene is significant in ANY of a set of contrasts
is_sig_any <- function(contrast_set) {
  if (length(contrast_set) == 0) return(rep(FALSE, nrow(consensus)))
  cols <- paste0("sig_", tolower(contrast_set))
  cols <- intersect(cols, names(consensus))
  if (length(cols) == 0) return(rep(FALSE, nrow(consensus)))
  rowSums(consensus[, ..cols], na.rm = TRUE) > 0
}

# Compute flags
sig_onset     <- is_sig_any(onset_ids)
sig_prog      <- is_sig_any(prog_ids)
sig_fib       <- is_sig_any(fib_ids)
sig_inflam_c2 <- is_sig_any(c("C2"))       # NAFL vs NASH (inflammation component)
sig_inflam_c5 <- is_sig_any(c("C5"))       # NAS >= 5 (inflammation component)
sig_c7b       <- is_sig_any(c("C7b"))      # inflammation ordinal
sig_c7a       <- is_sig_any(c("C7a"))      # steatosis ordinal
sig_c3        <- is_sig_any(c("C3"))       # fibrosis
sig_c6_only   <- is_sig_any(c("C6"))       # extreme endpoints

# Classification priority (checked in order; first match wins)
# 1. ubiquitous (>= 40% of contrasts) — overrides all
# 2. both_onset_and_progression
# 3. onset_only
# 4. progression_only
# 5. fibrosis_specific
# 6. inflammation_specific
# 7. extreme_only
# 8. unclassified (sig in 1-4 contrasts but doesn't fit above)
# 9. not_significant (sig in 0 contrasts)

consensus[, gene_class := "not_significant"]

# Start from most specific, then overwrite with broader categories
# extreme_only: significant ONLY in C6
consensus[n_contrasts_sig == 1 & sig_c6_only,
  gene_class := "extreme_only"]

# inflammation_specific: sig in C7b but NOT in C7a or C3
consensus[sig_c7b & !sig_c7a & !sig_c3,
  gene_class := "inflammation_specific"]

# fibrosis_specific: sig in C3/C8/C9 but NOT in C2 or C5
consensus[sig_fib & !sig_inflam_c2 & !sig_inflam_c5,
  gene_class := "fibrosis_specific"]

# progression_only: sig in C2/C3/C5 but NOT in onset (C1/C4)
consensus[sig_prog & !sig_onset,
  gene_class := "progression_only"]

# onset_only: sig in C1/C4 but NOT in any progression (C2/C3/C5)
consensus[sig_onset & !sig_prog,
  gene_class := "onset_only"]

# both_onset_and_progression: sig in >= 1 onset AND >= 1 progression
consensus[sig_onset & sig_prog,
  gene_class := "both_onset_and_progression"]

# ubiquitous: proportional threshold — overrides everything above
consensus[n_contrasts_sig >= ceiling(length(contrast_ids) * 0.4),
  gene_class := "ubiquitous"]

# Print classification summary
cat("\n  Gene classification summary:\n")
class_summary <- consensus[, .N, by = gene_class][order(-N)]
for (i in seq_len(nrow(class_summary))) {
  cat(sprintf("    %-30s %6d genes\n", class_summary$gene_class[i], class_summary$N[i]))
}

# ============================================================
# 6. Merge gene annotations and save outputs
# ============================================================
cat("\n=== Saving outputs ===\n")

# --- Output 1: Full consensus matrix ---
consensus_out <- merge(consensus, gene_annot, by = "gene", all.x = TRUE)

# Reorder columns: gene, symbol, gene_type, then per-contrast stats, then classification
meta_cols <- c("gene", "symbol", "gene_type")
stat_cols <- sort(c(
  paste0("logFC_", tolower(contrast_ids)),
  paste0("padj_", tolower(contrast_ids)),
  paste0("sig_", tolower(contrast_ids)),
  paste0("t_", tolower(contrast_ids))
))
class_cols <- c("gene_class", "n_contrasts_sig", "tau", "peak_contrast")

# Ensure all expected columns exist; use only those present
avail_cols <- intersect(c(meta_cols, stat_cols, class_cols), names(consensus_out))
remaining <- setdiff(names(consensus_out), avail_cols)
setcolorder(consensus_out, c(avail_cols, remaining))

out_f1 <- file.path(ODIR, "progression_consensus_matrix.csv")
fwrite(consensus_out, out_f1)
cat(sprintf("  Saved: %s (%d genes x %d cols)\n",
  basename(out_f1), nrow(consensus_out), ncol(consensus_out)))

# --- Output 2: Classification table (compact) ---
classification <- consensus_out[, .(gene, symbol, gene_type, gene_class, tau,
  peak_contrast, n_contrasts_sig)]

# Sort: ubiquitous first (by n_contrasts_sig desc), then by tau desc
classification <- classification[order(-n_contrasts_sig, -tau)]

out_f2 <- file.path(ODIR, "gene_progression_classification.csv")
fwrite(classification, out_f2)
cat(sprintf("  Saved: %s (%d genes)\n", basename(out_f2), nrow(classification)))

# --- Output 3: Classification summary ---
summary_dt <- consensus_out[, .(
  n_genes = .N,
  n_protein_coding = sum(gene_type == "protein_coding", na.rm = TRUE),
  n_lncRNA = sum(gene_type == "lncRNA", na.rm = TRUE),
  mean_tau = round(mean(tau, na.rm = TRUE), 3),
  median_n_contrasts = median(n_contrasts_sig, na.rm = TRUE)
), by = gene_class][order(-n_genes)]

out_f3 <- file.path(ODIR, "progression_classification_summary.csv")
fwrite(summary_dt, out_f3)
cat(sprintf("  Saved: %s (%d classes)\n", basename(out_f3), nrow(summary_dt)))

# ============================================================
# 7. Diagnostic summaries
# ============================================================
cat("\n=== Diagnostic Summary ===\n")

# Top ubiquitous genes (highest n_contrasts)
cat("\nTop 20 ubiquitous genes (by n_contrasts_sig):\n")
ubi <- consensus_out[gene_class == "ubiquitous"][order(-n_contrasts_sig, -abs(logFC_c1))]
if (nrow(ubi) > 0) {
  for (i in seq_len(min(20, nrow(ubi)))) {
    cat(sprintf("  %-15s tau=%.3f  n_sig=%d  peak=%s\n",
      ubi$symbol[i], ubi$tau[i], ubi$n_contrasts_sig[i], ubi$peak_contrast[i]))
  }
}

# Top onset-only genes
cat("\nTop 20 onset-only genes (by tau):\n")
onset <- consensus_out[gene_class == "onset_only"][order(-tau)]
if (nrow(onset) > 0) {
  for (i in seq_len(min(20, nrow(onset)))) {
    cat(sprintf("  %-15s tau=%.3f  n_sig=%d  peak=%s\n",
      onset$symbol[i], onset$tau[i], onset$n_contrasts_sig[i], onset$peak_contrast[i]))
  }
}

# Top progression-only genes
cat("\nTop 20 progression-only genes (by tau):\n")
prog <- consensus_out[gene_class == "progression_only"][order(-tau)]
if (nrow(prog) > 0) {
  for (i in seq_len(min(20, nrow(prog)))) {
    cat(sprintf("  %-15s tau=%.3f  n_sig=%d  peak=%s\n",
      prog$symbol[i], prog$tau[i], prog$n_contrasts_sig[i], prog$peak_contrast[i]))
  }
}

# Top fibrosis-specific genes
cat("\nTop 20 fibrosis-specific genes (by tau):\n")
fib_spec <- consensus_out[gene_class == "fibrosis_specific"][order(-tau)]
if (nrow(fib_spec) > 0) {
  for (i in seq_len(min(20, nrow(fib_spec)))) {
    cat(sprintf("  %-15s tau=%.3f  n_sig=%d  peak=%s\n",
      fib_spec$symbol[i], fib_spec$tau[i], fib_spec$n_contrasts_sig[i],
      fib_spec$peak_contrast[i]))
  }
}

# Top inflammation-specific genes
cat("\nTop 20 inflammation-specific genes (by tau):\n")
inf_spec <- consensus_out[gene_class == "inflammation_specific"][order(-tau)]
if (nrow(inf_spec) > 0) {
  for (i in seq_len(min(20, nrow(inf_spec)))) {
    cat(sprintf("  %-15s tau=%.3f  n_sig=%d  peak=%s\n",
      inf_spec$symbol[i], inf_spec$tau[i], inf_spec$n_contrasts_sig[i],
      inf_spec$peak_contrast[i]))
  }
}

# Cross-tabulation: gene_class vs peak_contrast
cat("\nGene class x Peak contrast cross-tabulation:\n")
ct <- consensus_out[gene_class != "not_significant", .N, by = .(gene_class, peak_contrast)]
ct_wide <- dcast(ct, gene_class ~ peak_contrast, value.var = "N", fill = 0)
print(ct_wide)

cat("\n=== 134: COMPLETE ===\n")
cat("Finished:", as.character(Sys.time()), "\n")
