#!/usr/bin/env Rscript
# 142_progression_tf_activity.R
# ---------------------------------------------------------------------------
# Compute differential TF activity for progression contrasts using decoupleR.
#
# For each contrast (C2: NAFL-vs-NASH, C3: Adv-vs-Early-Fib,
# C5: NAS>=5-vs-NAS<5), computes per-sample TF activity scores using
# decoupleR (ULM method) with DoRothEA regulons (confidence A-C), then
# runs Wilcoxon rank-sum tests for differential activity between groups.
#
# Follows patterns from Script 51 (decoupleR functional activity).
#
# Output: results/progression/
#   progression_tf_activity_c2.csv  — NAFL-vs-NASH differential TF activity
#   progression_tf_activity_c3.csv  — Fibrosis differential TF activity
#   progression_tf_activity_c5.csv  — NAS binary differential TF activity
#   progression_tf_activity_summary.csv — top TFs per contrast
#
# Usage: Rscript 142_progression_tf_activity.R
# Compute: login node OK (~10 min)
# Requires: decoupleR 2.12+, edgeR, data.table
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(decoupleR)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
ODIR <- file.path(INT, "results/progression")
dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 142: Progression TF Activity (decoupleR) ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# ── 1. Load expression data ─────────────────────────────────────────────
cat("--- Loading expression data ---\n")

dge_path <- file.path(INT, "results/integration/merged_dge.rds")
if (!file.exists(dge_path)) {
  stop("merged_dge.rds not found: ", dge_path)
}
dge <- readRDS(dge_path)
cat("  DGEList loaded:", nrow(dge), "genes x", ncol(dge), "samples\n")

# Also load raw counts for TMM-normalised logCPM
counts_path <- file.path(INT, "results/integration/merged_counts_raw.rds")
if (file.exists(counts_path)) {
  counts_raw <- readRDS(counts_path)
  cat("  Raw counts loaded:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")
} else {
  counts_raw <- NULL
  cat("  WARNING: merged_counts_raw.rds not found, using DGE object only\n")
}

# ── 2. Load metadata ────────────────────────────────────────────────────
cat("\n--- Loading metadata ---\n")

# QC-passing samples
qc_path <- file.path(INT, "qc/sample_qc_report.csv")
if (!file.exists(qc_path)) {
  stop("sample_qc_report.csv not found: ", qc_path)
}
qc <- fread(qc_path)
pass_ids <- qc[pass_technical == TRUE, sample_id]

# Modeling metadata for fibrosis/NAS/diagnosis
mm_path <- file.path(INT, "results/staging_classifier/modeling_metadata.csv")
if (!file.exists(mm_path)) {
  stop("modeling_metadata.csv not found: ", mm_path)
}
modeling_meta <- fread(mm_path)
cat("  Modeling metadata:", nrow(modeling_meta), "samples\n")

# Load unified metadata for dataset info
meta_path <- file.path(INT, "results/integration/meta_matched.rds")
if (file.exists(meta_path)) {
  meta <- readRDS(meta_path)
  meta <- meta[sample_id %in% pass_ids]
  # Merge in modeling metadata
  meta <- merge(meta, modeling_meta[, .(sample_id, fibrosis_stage, nas_score,
    diagnosis_harmonized, fib_ge3, nas_ge5)],
    by = "sample_id", all.x = TRUE, suffixes = c("", ".mm"))
  # Resolve duplicates — prefer modeling_metadata values
  for (col in c("fibrosis_stage", "nas_score", "diagnosis_harmonized")) {
    mm_col <- paste0(col, ".mm")
    if (mm_col %in% names(meta)) {
      na_mask <- is.na(meta[[col]]) | meta[[col]] == ""
      if (any(na_mask)) meta[[col]][na_mask] <- meta[[mm_col]][na_mask]
      meta[, (mm_col) := NULL]
    }
  }
} else {
  # Fallback: use modeling_metadata directly
  meta <- modeling_meta[sample_id %in% pass_ids]
}

# Restrict DGE to QC-passing samples present in both DGE and metadata
shared_ids <- intersect(intersect(colnames(dge), pass_ids), meta$sample_id)
cat("  QC-passing samples in expression data:", length(shared_ids), "\n")

# ── 3. Build normalised expression matrix ────────────────────────────────
cat("\n--- Building normalised expression matrix ---\n")

# Use TMM-normalised logCPM for per-sample TF activity
if (!is.null(counts_raw)) {
  shared_genes <- intersect(rownames(dge), rownames(counts_raw))
  dge_sub <- DGEList(counts = counts_raw[shared_genes, shared_ids])
} else {
  dge_sub <- dge[, shared_ids]
}
dge_sub <- calcNormFactors(dge_sub, method = "TMM")
logcpm <- cpm(dge_sub, log = TRUE, prior.count = 1)
cat("  logCPM matrix:", nrow(logcpm), "genes x", ncol(logcpm), "samples\n")

# Map gene IDs to symbols for decoupleR
# Gene annotations from dream results
dream_path <- file.path(INT, "results/integration/dream_results.csv")
# Gene annotations — dream_results.csv may lack symbol; try disease signatures first
annot_file <- file.path(INT, "results/disease_signatures/unified_disease_signatures.csv")
if (file.exists(annot_file)) {
  gene_annot <- fread(annot_file, select = c("gene", "symbol"))
  gene_annot <- unique(gene_annot, by = "gene")
  gene_annot[, ensembl_clean := sub("\\..*", "", gene)]
} else if (file.exists(dream_path)) {
  avail_cols <- names(fread(dream_path, nrows = 0))
  sel_cols <- intersect(c("gene", "symbol"), avail_cols)
  gene_annot <- fread(dream_path, select = sel_cols)
  if (!"symbol" %in% names(gene_annot)) gene_annot[, symbol := NA_character_]
  gene_annot[, ensembl_clean := sub("\\..*", "", gene)]
  gene_annot <- gene_annot[!is.na(symbol) & symbol != ""]
  gene_annot <- gene_annot[!duplicated(ensembl_clean)]
} else {
  stop("dream_results.csv not found")
}

# Map rownames (ensembl IDs) to symbols
row_ids <- sub("\\..*", "", rownames(logcpm))
idx <- match(row_ids, gene_annot$ensembl_clean)
valid <- !is.na(idx)
logcpm_sym <- logcpm[valid, ]
rownames(logcpm_sym) <- gene_annot$symbol[idx[valid]]

# Remove duplicated symbols (keep first)
dup_mask <- !duplicated(rownames(logcpm_sym))
logcpm_sym <- logcpm_sym[dup_mask, ]
cat("  Symbol-mapped matrix:", nrow(logcpm_sym), "genes x", ncol(logcpm_sym), "samples\n")

# ── 4. Load DoRothEA regulons ────────────────────────────────────────────
cat("\n--- Loading DoRothEA regulons (confidence A-C) ---\n")
tf_net <- as.data.table(get_dorothea(organism = "human", levels = c("A", "B", "C")))
cat("  DoRothEA:", nrow(tf_net), "interactions,",
    length(unique(tf_net$source)), "TFs\n")

# ── 5. Compute per-sample TF activity via ULM ───────────────────────────
cat("\n--- Computing per-sample TF activity (ULM) ---\n")

tf_acts <- as.data.table(run_ulm(
  mat = logcpm_sym,
  net = tf_net,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5
))

# Pivot to TF x sample matrix
tf_wide <- dcast(tf_acts, source ~ condition, value.var = "score")
tf_names <- tf_wide$source
tf_mat <- as.matrix(tf_wide[, -1, with = FALSE])
rownames(tf_mat) <- tf_names
cat("  TF activity matrix:", nrow(tf_mat), "TFs x", ncol(tf_mat), "samples\n")

# ── 6. Define contrasts ─────────────────────────────────────────────────
cat("\n--- Defining progression contrasts ---\n")

build_contrast <- function(contrast_id, meta_dt, group_col, levels_vec, description) {
  cat(sprintf("\n  %s: %s\n", contrast_id, description))
  # Subset to samples with valid group assignment
  sub <- meta_dt[!is.na(get(group_col)) & get(group_col) != ""]
  sub <- sub[get(group_col) %in% levels_vec]
  # Restrict to samples in TF activity matrix
  sub <- sub[sample_id %in% colnames(tf_mat)]
  sub[, contrast_group := factor(get(group_col), levels = levels_vec)]

  cat(sprintf("    Samples: %d (%s=%d, %s=%d)\n",
    nrow(sub), levels_vec[1], sum(sub$contrast_group == levels_vec[1]),
    levels_vec[2], sum(sub$contrast_group == levels_vec[2])))

  if (nrow(sub) < 10) {
    cat("    WARNING: Too few samples, skipping.\n")
    return(NULL)
  }
  sub
}

# C2: NASH vs NAFL
c2_meta <- build_contrast("C2", meta, "diagnosis_harmonized",
  c("NAFL", "NASH"), "NASH vs NAFL")

# C3: Advanced vs Early Fibrosis (F3-F4 vs F0-F2)
meta[, fib_binary := fifelse(fibrosis_stage %in% c(3, 4), "Advanced",
                     fifelse(fibrosis_stage %in% c(0, 1, 2), "Early", NA_character_))]
c3_meta <- build_contrast("C3", meta, "fib_binary",
  c("Early", "Advanced"), "Adv (F3-F4) vs Early (F0-F2) Fibrosis")

# C5: NAS >= 5 vs NAS < 5
meta[, nas_binary := fifelse(!is.na(nas_score) & nas_score >= 5, "NAS_high",
                    fifelse(!is.na(nas_score) & nas_score < 5, "NAS_low", NA_character_))]
c5_meta <- build_contrast("C5", meta, "nas_binary",
  c("NAS_low", "NAS_high"), "NAS>=5 vs NAS<5")

# ── 7. Differential TF activity per contrast ────────────────────────────
cat("\n--- Running differential TF activity ---\n")

run_diff_tf <- function(contrast_id, meta_sub, levels_vec) {
  if (is.null(meta_sub)) {
    cat(sprintf("  %s: SKIPPED (no valid metadata)\n", contrast_id))
    return(NULL)
  }

  group1_ids <- meta_sub[contrast_group == levels_vec[1], sample_id]
  group2_ids <- meta_sub[contrast_group == levels_vec[2], sample_id]

  # Restrict to TFs present in matrix
  results <- rbindlist(lapply(rownames(tf_mat), function(tf) {
    vals1 <- tf_mat[tf, group1_ids, drop = TRUE]
    vals2 <- tf_mat[tf, group2_ids, drop = TRUE]

    # Wilcoxon rank-sum test
    wt <- tryCatch(
      wilcox.test(vals2, vals1, conf.int = FALSE),
      error = function(e) NULL
    )

    # Also compute t-test for effect size
    tt <- tryCatch(
      t.test(vals2, vals1),
      error = function(e) NULL
    )

    data.table(
      tf = tf,
      mean_group1 = mean(vals1, na.rm = TRUE),
      mean_group2 = mean(vals2, na.rm = TRUE),
      diff = mean(vals2, na.rm = TRUE) - mean(vals1, na.rm = TRUE),
      wilcox_pval = if (!is.null(wt)) wt$p.value else NA_real_,
      ttest_pval = if (!is.null(tt)) tt$p.value else NA_real_,
      ttest_stat = if (!is.null(tt)) tt$statistic else NA_real_
    )
  }))

  results[, wilcox_padj := p.adjust(wilcox_pval, method = "BH")]
  results[, ttest_padj := p.adjust(ttest_pval, method = "BH")]
  results[, direction := fifelse(diff > 0, "activated_in_group2", "repressed_in_group2")]
  results[, abs_diff := abs(diff)]
  setorder(results, wilcox_pval)

  # Add contrast metadata
  results[, contrast := contrast_id]
  results[, group1 := levels_vec[1]]
  results[, group2 := levels_vec[2]]
  results[, n_group1 := length(group1_ids)]
  results[, n_group2 := length(group2_ids)]

  results
}

c2_results <- run_diff_tf("C2_NASH_vs_NAFL", c2_meta, c("NAFL", "NASH"))
c3_results <- run_diff_tf("C3_Adv_vs_Early_Fib", c3_meta, c("Early", "Advanced"))
c5_results <- run_diff_tf("C5_NAS_high_vs_low", c5_meta, c("NAS_low", "NAS_high"))

# ── 8. Write per-contrast results ───────────────────────────────────────
cat("\n--- Writing per-contrast results ---\n")

write_contrast <- function(dt, filename, label) {
  if (is.null(dt)) {
    cat(sprintf("  %s: SKIPPED\n", label))
    return()
  }
  fwrite(dt, file.path(ODIR, filename))
  n_sig <- dt[wilcox_padj < 0.05, .N]
  n_act <- dt[wilcox_padj < 0.05 & direction == "activated_in_group2", .N]
  n_rep <- dt[wilcox_padj < 0.05 & direction == "repressed_in_group2", .N]
  cat(sprintf("  %s: %d TFs tested, %d significant (padj<0.05): %d activated, %d repressed\n",
    label, nrow(dt), n_sig, n_act, n_rep))
  if (n_sig > 0) {
    cat("    Top 10 TFs:\n")
    print(dt[1:min(10, n_sig), .(tf, diff = round(diff, 3),
      wilcox_padj = signif(wilcox_padj, 3), direction)])
  }
}

write_contrast(c2_results, "progression_tf_activity_c2.csv", "C2 NAFL-vs-NASH")
write_contrast(c3_results, "progression_tf_activity_c3.csv", "C3 Adv-vs-Early Fibrosis")
write_contrast(c5_results, "progression_tf_activity_c5.csv", "C5 NAS>=5-vs-NAS<5")

# ── 9. Build summary across contrasts ───────────────────────────────────
cat("\n--- Building cross-contrast summary ---\n")

all_results <- rbindlist(Filter(Negate(is.null), list(c2_results, c3_results, c5_results)))

if (nrow(all_results) > 0) {
  # Summarise per TF: how many contrasts significant, consistency of direction
  tf_summary <- all_results[, .(
    n_contrasts_tested = .N,
    n_contrasts_sig = sum(wilcox_padj < 0.05, na.rm = TRUE),
    contrasts_sig = paste(contrast[wilcox_padj < 0.05], collapse = ";"),
    mean_diff = mean(diff, na.rm = TRUE),
    max_abs_diff = max(abs_diff, na.rm = TRUE),
    min_padj = min(wilcox_padj, na.rm = TRUE),
    direction_consistent = length(unique(sign(diff[wilcox_padj < 0.05]))) <= 1
  ), by = tf]

  setorder(tf_summary, -n_contrasts_sig, min_padj)

  fwrite(tf_summary, file.path(ODIR, "progression_tf_activity_summary.csv"))

  n_multi <- tf_summary[n_contrasts_sig >= 2, .N]
  n_all <- tf_summary[n_contrasts_sig == tf_summary[, max(n_contrasts_tested)], .N]

  cat(sprintf("  Total TFs assessed: %d\n", nrow(tf_summary)))
  cat(sprintf("  Significant in >=2 contrasts: %d\n", n_multi))
  cat(sprintf("  Significant in all contrasts: %d\n", n_all))

  if (n_multi > 0) {
    cat("\n  TFs significant in >=2 progression contrasts:\n")
    print(tf_summary[n_contrasts_sig >= 2,
      .(tf, n_sig = n_contrasts_sig, mean_diff = round(mean_diff, 3),
        min_padj = signif(min_padj, 3), consistent = direction_consistent)][1:min(20, n_multi)])
  }

  # Key drug-target TFs
  cat("\n  Key MASLD drug-target TF status:\n")
  key_tfs <- c("PPARA", "PPARG", "PPARD", "NR1H4", "THRB", "HNF4A", "NFE2L2",
               "SREBF1", "SREBF2", "NR1H3", "RXRA", "CEBPA", "FOXA1", "FOXA2")
  for (tf_name in key_tfs) {
    row <- tf_summary[tf == tf_name]
    if (nrow(row) > 0) {
      cat(sprintf("    %s: %d/%d contrasts sig, mean_diff=%.3f, min_padj=%.2e\n",
        tf_name, row$n_contrasts_sig, row$n_contrasts_tested,
        row$mean_diff, row$min_padj))
    }
  }
} else {
  cat("  WARNING: No contrast results to summarise\n")
}

# ── Summary ──────────────────────────────────────────────────────────────
cat("\n=== RESULTS SUMMARY ===\n")
for (r in list(list(c2_results, "C2"), list(c3_results, "C3"), list(c5_results, "C5"))) {
  dt <- r[[1]]; lbl <- r[[2]]
  if (!is.null(dt)) {
    cat(sprintf("  %s: %d TFs, %d sig (padj<0.05)\n", lbl, nrow(dt),
      dt[wilcox_padj < 0.05, .N]))
  } else {
    cat(sprintf("  %s: SKIPPED\n", lbl))
  }
}
cat("\nOutputs written to:", ODIR, "\n")
cat("  progression_tf_activity_c2.csv\n")
cat("  progression_tf_activity_c3.csv\n")
cat("  progression_tf_activity_c5.csv\n")
cat("  progression_tf_activity_summary.csv\n")
cat("\nDone:", format(Sys.time()), "\n")
