#!/usr/bin/env Rscript
# 145_expand_atlas_progression.R
# Expand multi-evidence atlas with ~25 progression columns
#
# Loads the existing multi-evidence atlas and merges per-gene progression
# results from Scripts 130-144.  All joins use the atlas's unversioned
# Ensembl ID (ensembl_id) matched against the versioned IDs produced by the
# dream pipeline (strip version from progression "gene" columns).
#
# Input:
#   RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
#   RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/*
#   RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression/*
#
# Output:
#   RNA-seq/results/multi_evidence/multi_evidence_atlas.csv        (overwritten)
#   RNA-seq/results/multi_evidence/multi_evidence_atlas_pre_progression.csv (backup)
#
# Usage:
#   Rscript 145_expand_atlas_progression.R
# ============================================================

suppressPackageStartupMessages(library(data.table))

cat("=== Script 145: Expand Atlas with Progression Columns ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ============================================================
# Paths
# ============================================================
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATLAS_DIR  <- file.path(BASE, "RNA-seq/results/multi_evidence")
ATLAS_FILE <- file.path(ATLAS_DIR, "multi_evidence_atlas.csv")
BACKUP_FILE <- file.path(ATLAS_DIR, "multi_evidence_atlas_pre_progression.csv")

INT_RESULTS <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results")
SIG_DIR  <- file.path(INT_RESULTS, "disease_signatures")
PROG_DIR <- file.path(INT_RESULTS, "progression")

# ============================================================
# 1. Load atlas
# ============================================================
cat("--- Loading atlas ---\n")
stopifnot(file.exists(ATLAS_FILE))
atlas <- fread(ATLAS_FILE)
# If a pre-progression backup exists, load THAT instead (avoids
# duplicate columns from re-running this script on an already-expanded atlas)
if (file.exists(BACKUP_FILE)) {
  cat("  Found pre-progression backup — loading clean atlas to avoid duplicates\n")
  atlas <- fread(BACKUP_FILE)
} else {
  # First run — back up the clean atlas
  cat(sprintf("  Creating backup: %s\n", basename(BACKUP_FILE)))
  fwrite(atlas, BACKUP_FILE)
}

n_orig_cols <- ncol(atlas)
cat(sprintf("  Atlas: %d genes x %d columns\n", nrow(atlas), n_orig_cols))

# ============================================================
# Helper: safe loader  (returns NULL if file missing)
# ============================================================
safe_load <- function(path, label) {
  if (!file.exists(path)) {
    cat(sprintf("  [SKIP] %s — file not found\n", label))
    return(NULL)
  }
  dt <- fread(path)
  cat(sprintf("  [OK]   %s — %d rows x %d cols\n", label, nrow(dt), ncol(dt)))
  dt
}

# ============================================================
# Helper: strip Ensembl version to match atlas ensembl_id
# The atlas `ensembl_id` is unversioned (ENSG00000000003).
# Progression files use versioned IDs (ENSG00000000003.17).
# ============================================================
strip_version <- function(dt, gene_col = "gene") {
  dt[, ensembl_id := sub("\\..*", "", get(gene_col))]
  dt
}

# ============================================================
# Helper: merge a small data.table onto atlas by ensembl_id
# Returns the updated atlas (by reference where possible)
# ============================================================
merge_onto_atlas <- function(atlas, layer, new_cols, label) {
  if (is.null(layer) || nrow(layer) == 0) {
    cat(sprintf("  [SKIP MERGE] %s — no data\n", label))
    return(atlas)
  }

  # De-duplicate by ensembl_id (keep first)
  layer <- layer[!duplicated(ensembl_id)]

  # Check for column name conflicts with existing atlas
  existing <- names(atlas)
  conflict <- intersect(new_cols, existing)
  if (length(conflict) > 0) {
    cat(sprintf("  WARNING: Column conflict for %s: %s → prefixing with prog_\n",
      label, paste(conflict, collapse = ", ")))
    for (cc in conflict) {
      new_name <- paste0("prog_", cc)
      setnames(layer, cc, new_name)
      new_cols[new_cols == cc] <- new_name
    }
  }

  # Select join key + new columns only
  keep_cols <- intersect(c("ensembl_id", new_cols), names(layer))
  layer_slim <- layer[, ..keep_cols]

  before <- ncol(atlas)
  atlas <- merge(atlas, layer_slim, by = "ensembl_id", all.x = TRUE)
  after <- ncol(atlas)
  n_non_na <- sum(!is.na(atlas[[new_cols[1]]]))
  cat(sprintf("  [MERGED] %s: +%d cols, %d/%d genes with data (%.1f%%)\n",
    label, after - before, n_non_na, nrow(atlas),
    100 * n_non_na / nrow(atlas)))

  atlas
}

# ============================================================
# 3. Load and merge each progression data source
# ============================================================
cat("\n--- Loading progression results ---\n\n")

# ------ 3a. NAFL vs NASH (disease_signatures/nafl_vs_nash_dream.csv) ------
# This uses adj.P.Val column naming (limma/dream convention with symbol join)
nn <- safe_load(file.path(SIG_DIR, "nafl_vs_nash_dream.csv"), "NAFL vs NASH dream")
if (!is.null(nn)) {
  strip_version(nn)
  # Handle both adj.P.Val and padj column naming
  padj_col <- if ("adj.P.Val" %in% names(nn)) "adj.P.Val" else "padj"
  nn[, nafl_vs_nash_logFC := logFC]
  nn[, nafl_vs_nash_padj  := get(padj_col)]
  nn[, nafl_vs_nash_tstat := t]
  atlas <- merge_onto_atlas(atlas, nn,
    c("nafl_vs_nash_logFC", "nafl_vs_nash_padj", "nafl_vs_nash_tstat"),
    "NAFL vs NASH")
}

# ------ 3b. C3: Advanced vs Early Fibrosis ------
c3 <- safe_load(file.path(PROG_DIR, "c3_adv_vs_early_fib_dream.csv"),
  "C3: Adv vs Early Fib")
if (!is.null(c3)) {
  strip_version(c3)
  padj_col <- if ("adj.P.Val" %in% names(c3)) "adj.P.Val" else "padj"
  c3[, adv_fib_logFC := logFC]
  c3[, adv_fib_padj  := get(padj_col)]
  atlas <- merge_onto_atlas(atlas, c3,
    c("adv_fib_logFC", "adv_fib_padj"), "C3 Adv vs Early Fib")
}

# ------ 3c. C4: NAFL vs Control ------
c4 <- safe_load(file.path(PROG_DIR, "c4_nafl_vs_ctrl_dream.csv"),
  "C4: NAFL vs Ctrl")
if (!is.null(c4)) {
  strip_version(c4)
  padj_col <- if ("adj.P.Val" %in% names(c4)) "adj.P.Val" else "padj"
  c4[, nafl_vs_ctrl_logFC := logFC]
  c4[, nafl_vs_ctrl_padj  := get(padj_col)]
  atlas <- merge_onto_atlas(atlas, c4,
    c("nafl_vs_ctrl_logFC", "nafl_vs_ctrl_padj"), "C4 NAFL vs Ctrl")
}

# ------ 3d. C5: NAS >= 5 vs NAS < 5 ------
c5 <- safe_load(file.path(PROG_DIR, "c5_nas_ge5_vs_lt5_dream.csv"),
  "C5: NAS >= 5 vs < 5")
if (!is.null(c5)) {
  strip_version(c5)
  padj_col <- if ("adj.P.Val" %in% names(c5)) "adj.P.Val" else "padj"
  c5[, nas_ge5_logFC := logFC]
  c5[, nas_ge5_padj  := get(padj_col)]
  atlas <- merge_onto_atlas(atlas, c5,
    c("nas_ge5_logFC", "nas_ge5_padj"), "C5 NAS >= 5")
}

# ------ 3e. C6: Extreme Endpoints ------
c6 <- safe_load(file.path(PROG_DIR, "c6_extreme_endpoints_dream.csv"),
  "C6: Extreme Endpoints")
if (!is.null(c6)) {
  strip_version(c6)
  padj_col <- if ("adj.P.Val" %in% names(c6)) "adj.P.Val" else "padj"
  c6[, extreme_logFC := logFC]
  c6[, extreme_padj  := get(padj_col)]
  atlas <- merge_onto_atlas(atlas, c6,
    c("extreme_logFC", "extreme_padj"), "C6 Extreme Endpoints")
}

# ------ 3f. C7a: Steatosis ordinal ------
c7a <- safe_load(file.path(PROG_DIR, "c7a_steatosis_ordinal_dream.csv"),
  "C7a: Steatosis ordinal")
if (!is.null(c7a)) {
  strip_version(c7a)
  c7a[, steatosis_ordinal_coef := logFC]
  atlas <- merge_onto_atlas(atlas, c7a,
    c("steatosis_ordinal_coef"), "C7a Steatosis ordinal")
}

# ------ 3g. C7b: Inflammation ordinal ------
c7b <- safe_load(file.path(PROG_DIR, "c7b_inflammation_ordinal_dream.csv"),
  "C7b: Inflammation ordinal")
if (!is.null(c7b)) {
  strip_version(c7b)
  c7b[, inflammation_ordinal_coef := logFC]
  atlas <- merge_onto_atlas(atlas, c7b,
    c("inflammation_ordinal_coef"), "C7b Inflammation ordinal")
}

# ------ 3h. C7c: Ballooning ordinal ------
c7c <- safe_load(file.path(PROG_DIR, "c7c_ballooning_ordinal_dream.csv"),
  "C7c: Ballooning ordinal")
if (!is.null(c7c)) {
  strip_version(c7c)
  c7c[, ballooning_ordinal_coef := logFC]
  atlas <- merge_onto_atlas(atlas, c7c,
    c("ballooning_ordinal_coef"), "C7c Ballooning ordinal")
}

# ------ 3i. C8: Cirrhosis (F4 vs F0-F3) ------
c8 <- safe_load(file.path(PROG_DIR, "c8_cirrhosis_dream.csv"),
  "C8: Cirrhosis")
if (!is.null(c8)) {
  strip_version(c8)
  padj_col <- if ("adj.P.Val" %in% names(c8)) "adj.P.Val" else "padj"
  c8[, cirrhosis_logFC := logFC]
  c8[, cirrhosis_padj  := get(padj_col)]
  atlas <- merge_onto_atlas(atlas, c8,
    c("cirrhosis_logFC", "cirrhosis_padj"), "C8 Cirrhosis")
}

# ------ 3j. C9: F2 Inflection ------
c9 <- safe_load(file.path(PROG_DIR, "c9_f2_inflection_dream.csv"),
  "C9: F2 Inflection")
if (!is.null(c9)) {
  strip_version(c9)
  padj_col <- if ("adj.P.Val" %in% names(c9)) "adj.P.Val" else "padj"
  c9[, f2_inflection_logFC := logFC]
  c9[, f2_inflection_padj  := get(padj_col)]
  atlas <- merge_onto_atlas(atlas, c9,
    c("f2_inflection_logFC", "f2_inflection_padj"), "C9 F2 Inflection")
}

# ------ 3k. Deconvolution attribution (C2: NAFL vs NASH) ------
deconv_c2 <- safe_load(file.path(PROG_DIR, "progression_deconv_attribution_c2.csv"),
  "Deconv attribution (C2)")
if (!is.null(deconv_c2)) {
  strip_version(deconv_c2)
  deconv_c2[, progression_deconv_class := category]
  atlas <- merge_onto_atlas(atlas, deconv_c2,
    c("progression_deconv_class"), "Deconv attribution C2")
}

# ------ 3l. Progression TWAS unified ------
# Try the unified file first (best p across contrasts); fall back to summary
twas_unified <- safe_load(file.path(PROG_DIR, "progression_twas_unified.csv"),
  "Progression TWAS unified")
if (!is.null(twas_unified)) {
  # The unified file has gene_base (unversioned Ensembl) + symbol
  if ("gene_base" %in% names(twas_unified)) {
    setnames(twas_unified, "gene_base", "ensembl_id")
  } else {
    strip_version(twas_unified)
  }
  twas_unified[, progression_twas_n_contrasts := n_contrasts_sig]
  twas_unified[, progression_twas_direction := direction]
  atlas <- merge_onto_atlas(atlas, twas_unified,
    c("progression_twas_n_contrasts", "progression_twas_direction"),
    "Progression TWAS unified")
} else {
  # Fall back: try individual per-contrast TWAS files and take best
  twas_files <- list.files(PROG_DIR, pattern = "^progression_twas_c[0-9]",
                           full.names = TRUE)
  if (length(twas_files) > 0) {
    cat(sprintf("  Falling back to %d individual TWAS files\n", length(twas_files)))
    twas_all <- rbindlist(lapply(twas_files, function(f) {
      dt <- fread(f)
      if (!"gene" %in% names(dt) && "gene_base" %in% names(dt))
        setnames(dt, "gene_base", "gene")
      dt
    }), fill = TRUE)
    strip_version(twas_all)
    # Take best (min) p-value per gene across contrasts
    pval_col <- intersect(c("prog_padj", "twas_pval", "pvalue"), names(twas_all))
    if (length(pval_col) > 0) {
      twas_summary <- twas_all[, .(
        progression_twas_pval = min(get(pval_col[1]), na.rm = TRUE)
      ), by = ensembl_id]
      twas_summary <- twas_summary[is.finite(progression_twas_pval)]
      atlas <- merge_onto_atlas(atlas, twas_summary,
        c("progression_twas_pval"), "Progression TWAS (individual)")
    }
  } else {
    cat("  [SKIP] No TWAS files found\n")
  }
}

# ------ 3m. Progression COLOC combined ------
coloc <- safe_load(file.path(PROG_DIR, "progression_coloc_combined.csv"),
  "Progression COLOC combined")
if (!is.null(coloc) && nrow(coloc) > 0 && !"note" %in% names(coloc)) {
  strip_version(coloc)
  # Take best PP.H4 per gene across all contrasts x GWAS
  coloc_summary <- coloc[, .(
    progression_coloc_best_pp4 = max(PP.H4, na.rm = TRUE),
    progression_coloc_n_contrasts = uniqueN(contrast[PP.H4 > 0.5])
  ), by = ensembl_id]
  coloc_summary <- coloc_summary[is.finite(progression_coloc_best_pp4)]
  atlas <- merge_onto_atlas(atlas, coloc_summary,
    c("progression_coloc_best_pp4", "progression_coloc_n_contrasts"),
    "Progression COLOC")
}

# ------ 3n. Drug reversal hits (per-gene from leading edge) ------
drug_hits <- safe_load(file.path(PROG_DIR, "progression_drug_reversal_hits.csv"),
  "Drug reversal hits")
if (!is.null(drug_hits) && nrow(drug_hits) > 0 && !"note" %in% names(drug_hits)) {
  # The hits file is per-pathway. Extract genes from leadingEdge column.
  # leadingEdge is semicolon-separated Ensembl IDs (may be versioned)
  if ("leadingEdge" %in% names(drug_hits)) {
    le_genes <- drug_hits[, .(
      ensembl_versioned = unlist(strsplit(leadingEdge, ";"))
    ), by = .(contrast_id, pathway)]
    le_genes[, ensembl_id := sub("\\..*", "", ensembl_versioned)]
    # Count: per gene, how many unique pathway hits involve it
    drug_gene_summary <- le_genes[, .(
      progression_drug_reversal_n = uniqueN(pathway)
    ), by = ensembl_id]
    atlas <- merge_onto_atlas(atlas, drug_gene_summary,
      c("progression_drug_reversal_n"), "Drug reversal (per gene)")
  } else {
    cat("  [SKIP] Drug hits file lacks leadingEdge column\n")
  }
}

# ------ 3o. Sex-stratified progression classification ------
sex_cls <- safe_load(file.path(PROG_DIR, "c10_sex_progression_classification.csv"),
  "Sex progression classification")
if (!is.null(sex_cls)) {
  strip_version(sex_cls)
  sex_cls[, sex_progression_class_new := sex_progression_class]
  atlas <- merge_onto_atlas(atlas, sex_cls,
    c("sex_progression_class_new"), "Sex progression class")
  # Rename to desired name (avoids conflict check issue with intermediate name)
  setnames(atlas, "sex_progression_class_new", "sex_progression_class",
           skip_absent = TRUE)
}

# ------ 3p. Gene progression classification (from consensus Script 134) ------
gene_cls <- safe_load(file.path(PROG_DIR, "gene_progression_classification.csv"),
  "Gene progression classification")
if (!is.null(gene_cls)) {
  strip_version(gene_cls)
  # Derive binary flags from the classification
  gene_cls[, progression_gene_class := gene_class]
  gene_cls[, is_progression_specific := gene_class %in%
    c("progression_only", "fibrosis_specific", "inflammation_specific")]
  gene_cls[, is_onset_specific := gene_class == "onset_only"]
  gene_cls[, progression_n_contrasts_sig := n_contrasts_sig]
  gene_cls[, progression_tau := tau]
  gene_cls[, progression_peak_contrast := peak_contrast]
  atlas <- merge_onto_atlas(atlas, gene_cls,
    c("progression_gene_class", "is_progression_specific", "is_onset_specific",
      "progression_n_contrasts_sig", "progression_tau", "progression_peak_contrast"),
    "Gene progression classification")
}

# ------ 3r. Progression concordance (Script 140) ------
# Use the classified file which has per-gene summary columns
concordance <- safe_load(file.path(PROG_DIR, "progression_concordance_classified.csv"),
  "Progression concordance (classified)")
if (!is.null(concordance)) {
  # The classified file has mouse_gene_id, not Ensembl gene ID.
  # Aggregate across signatures: take max n_concordant and best category per gene.
  if ("mouse_gene_id" %in% names(concordance) && !"gene" %in% names(concordance)) {
    # Need ortholog mapping to get human gene IDs
    cat("  [NOTE] Concordance uses mouse_gene_id — skipping atlas merge (use concordance_summary.csv for cross-species)\n")
  } else if ("gene" %in% names(concordance)) {
    strip_version(concordance)
    avail <- intersect(c("n_concordant", "category", "n_diets_sig"),
      names(concordance))
    if (length(avail) > 0) {
      # Aggregate per gene across signatures
      conc_agg <- concordance[, .(
        progression_concordance_n = max(n_concordant, na.rm = TRUE),
        progression_concordance_class = category[which.max(n_concordant)],
        progression_n_diets_sig = max(n_diets_sig, na.rm = TRUE)
      ), by = gene]
      atlas <- merge(atlas, conc_agg, by = "gene", all.x = TRUE)
      n_filled <- sum(!is.na(atlas$progression_concordance_n))
      cat(sprintf("  Added concordance columns (%d/%d genes filled)\n",
        n_filled, nrow(atlas)))
    }
  }
}

# ------ 3s. Progression TF activity (Script 142) ------
tf_summary <- safe_load(file.path(PROG_DIR, "progression_tf_activity_summary.csv"),
  "Progression TF activity summary")
if (!is.null(tf_summary)) {
  # TF activity is per-TF, not per-gene, so we check if it has a gene-level mapping
  if ("gene" %in% names(tf_summary) || "ensembl_id" %in% names(tf_summary)) {
    if ("gene" %in% names(tf_summary)) strip_version(tf_summary)
    avail <- intersect(c("n_contrasts_sig", "mean_t"), names(tf_summary))
    if (length(avail) > 0) {
      rename_map <- c(n_contrasts_sig = "progression_tf_n_contrasts",
                      mean_t = "progression_tf_mean_t")
      for (old_nm in avail) {
        new_nm <- rename_map[old_nm]
        if (!is.na(new_nm)) setnames(tf_summary, old_nm, new_nm)
      }
      new_cols <- unname(rename_map[avail])
      new_cols <- new_cols[!is.na(new_cols)]
      atlas <- merge_onto_atlas(atlas, tf_summary, new_cols,
        "Progression TF activity")
    }
  } else {
    cat("  [NOTE] TF activity summary is per-TF, not per-gene — not merged\n")
  }
}

# ------ 3t. C11: NASH vs Control ------
c11 <- safe_load(file.path(PROG_DIR, "c11_nash_vs_ctrl_dream.csv"),
  "C11: NASH vs Ctrl")
if (!is.null(c11)) {
  strip_version(c11)
  padj_col <- if ("adj.P.Val" %in% names(c11)) "adj.P.Val" else "padj"
  c11[, nash_vs_ctrl_logFC := logFC]
  c11[, nash_vs_ctrl_padj  := get(padj_col)]
  atlas <- merge_onto_atlas(atlas, c11,
    c("nash_vs_ctrl_logFC", "nash_vs_ctrl_padj"), "C11 NASH vs Ctrl")
}

# ------ 3u. C12: Early vs Late NASH ------
c12 <- safe_load(file.path(PROG_DIR, "c12_early_vs_late_nash_dream.csv"),
  "C12: Early vs Late NASH")
if (!is.null(c12)) {
  strip_version(c12)
  padj_col <- if ("adj.P.Val" %in% names(c12)) "adj.P.Val" else "padj"
  c12[, early_late_nash_logFC := logFC]
  c12[, early_late_nash_padj  := get(padj_col)]
  atlas <- merge_onto_atlas(atlas, c12,
    c("early_late_nash_logFC", "early_late_nash_padj"), "C12 Early vs Late NASH")
}

# ------ 3v. C13: NASH vs NAFL (fibrosis-adjusted) ------
c13 <- safe_load(file.path(PROG_DIR, "c13_nash_vs_nafl_fib_adj_dream.csv"),
  "C13: NASH vs NAFL (fib-adj)")
if (!is.null(c13)) {
  strip_version(c13)
  padj_col <- if ("adj.P.Val" %in% names(c13)) "adj.P.Val" else "padj"
  c13[, nash_vs_nafl_fibadj_logFC := logFC]
  c13[, nash_vs_nafl_fibadj_padj  := get(padj_col)]
  atlas <- merge_onto_atlas(atlas, c13,
    c("nash_vs_nafl_fibadj_logFC", "nash_vs_nafl_fibadj_padj"),
    "C13 NASH vs NAFL (fib-adj)")
}

# ------ 3w. C17: Fibrosis ordinal (slope per unit increase) ------
c17 <- safe_load(file.path(PROG_DIR, "c17_fibrosis_ordinal_dream.csv"),
  "C17: Fibrosis ordinal")
if (!is.null(c17)) {
  strip_version(c17)
  c17[, fibrosis_ordinal_coef := logFC]
  atlas <- merge_onto_atlas(atlas, c17,
    c("fibrosis_ordinal_coef"), "C17 Fibrosis ordinal")
}

# ============================================================
# 4. Check for column name conflicts (final validation)
# ============================================================
cat("\n--- Validation ---\n")
dup_cols <- names(atlas)[duplicated(names(atlas))]
if (length(dup_cols) > 0) {
  cat(sprintf("  ERROR: Duplicate columns found: %s\n", paste(dup_cols, collapse = ", ")))
  # Remove duplicates (keep first)
  atlas <- atlas[, .SD, .SDcols = unique(names(atlas))]
  cat("  Removed duplicates, keeping first occurrence.\n")
}

# ============================================================
# 5. Summary report
# ============================================================
n_new_cols <- ncol(atlas) - n_orig_cols
new_col_names <- setdiff(names(atlas), names(fread(ATLAS_FILE, nrows = 0)))

cat(sprintf("\n--- Summary ---\n"))
cat(sprintf("  Original columns: %d\n", n_orig_cols))
cat(sprintf("  New columns:      %d\n", n_new_cols))
cat(sprintf("  Final columns:    %d\n", ncol(atlas)))
cat(sprintf("  Atlas rows:       %d (unchanged)\n", nrow(atlas)))

if (length(new_col_names) > 0) {
  cat("\n  New columns added:\n")
  for (i in seq_along(new_col_names)) {
    col <- new_col_names[i]
    n_na <- sum(is.na(atlas[[col]]))
    pct_na <- round(100 * n_na / nrow(atlas), 1)
    cat(sprintf("    %2d. %-40s  NA: %5d (%5.1f%%)\n", i, col, n_na, pct_na))
  }
}

# ============================================================
# 6. Save expanded atlas
# ============================================================
cat(sprintf("\n--- Saving expanded atlas ---\n"))
fwrite(atlas, ATLAS_FILE)
cat(sprintf("  Saved: %s (%d x %d)\n", basename(ATLAS_FILE), nrow(atlas), ncol(atlas)))

cat(sprintf("\n=== Script 145 complete: %s ===\n", as.character(Sys.time())))
