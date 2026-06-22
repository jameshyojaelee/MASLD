#!/usr/bin/env Rscript
# 75_integrate_causal_overhaul.R
# ---------------------------------------------------------------------------
# Phase 5A: Integrate all new causal inference results into the multi-evidence
# atlas. This is the final step after all Phases 1-4 complete.
#
# New columns added (expanding atlas by ~25 columns):
#   OTTERS TWAS:     otters_broadaway_z, otters_broadaway_pval, otters_n_gwas_sig
#   cTWAS:           ctwas_pip, ctwas_pval
#   SuSiE-COLOC:     coloc_susie_best_pp4, coloc_susie_n_signals
#   Multi-trait:     hyprcoloc_posterior, hyprcoloc_n_traits, hyprcoloc_traits
#   Derived:         causal_methods_sig, causal_robustness
#
# 2026-04-22: Enhanced MR (mr_n_instruments/mr_steiger_pass/mr_presso_pval/
# mr_raps_pval) and Bidirectional MR (mr_bidirectional_pval) blocks removed.
# MR is permanently ditched from the paper; TWAS + COLOC + INTACT is the
# causal framework. Archived MR scripts: archive/mr_ditched_2026-04-22/.
#
# Follows the idempotent merge-by-human_symbol pattern.
# ---------------------------------------------------------------------------

library(data.table)

cat("=== Script 75: Integrate Causal Overhaul Results ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATLAS_FILE <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
CAUSAL_DIR <- file.path(BASE, "RNA-seq/results/causal_inference")

stopifnot(file.exists(ATLAS_FILE))

atlas <- fread(ATLAS_FILE)
n_col_start <- ncol(atlas)
cat("Loaded atlas:", nrow(atlas), "genes x", n_col_start, "cols\n\n")

# Gene ID column — atlas uses "human_symbol" as primary key
gene_col <- "human_symbol"
if (!gene_col %in% names(atlas)) {
  # Fallback to gene
  gene_col <- "gene"
}

# ==============================================================================
# Helper: safe merge (left join, remove old cols if re-running)
# ==============================================================================
safe_merge <- function(atlas, new_dt, by_col = gene_col, new_cols) {
  # Remove existing columns (idempotent re-run)
  existing <- intersect(new_cols, names(atlas))
  if (length(existing) > 0) {
    atlas[, (existing) := NULL]
  }
  merge(atlas, new_dt, by = by_col, all.x = TRUE)
}

# ==============================================================================
# 1-2. Enhanced MR + Bidirectional MR — REMOVED 2026-04-22
# ==============================================================================
# MR permanently ditched from the paper. TWAS + COLOC + INTACT is the sole
# causal framework. The Enhanced-MR (TwoSampleMR) and Bidirectional-MR blocks
# that previously added mr_ivw_pval / mr_ivw_beta / mr_n_instruments /
# mr_n_gwas_sig / mr_presso_pval / mr_bidirectional_pval columns have been
# excised. Archived upstream scripts: archive/mr_ditched_2026-04-22/ (Scripts 19,
# 19b, 203, 36).
cat("--- 1-2. Enhanced/Bidirectional MR: REMOVED (MR ditched 2026-04-22) ---\n")

# ==============================================================================
# 3. OTTERS TWAS (Phase 2A)
# ==============================================================================
cat("\n--- 3. OTTERS TWAS ---\n")

otters_files <- list.files(file.path(CAUSAL_DIR, "otters_broadaway"),
                           pattern = "otters_twas_combined\\.csv",
                           recursive = TRUE, full.names = TRUE)

if (length(otters_files) > 0) {
  otters_all <- rbindlist(lapply(otters_files, function(f) {
    dt <- tryCatch(fread(f), error = function(e) NULL)
    if (!is.null(dt) && nrow(dt) > 0) {
      gwas <- basename(dirname(f))
      dt[, gwas_source := gwas]
      return(dt)
    }
    NULL
  }), fill = TRUE)

  if (nrow(otters_all) > 0) {
    cat("  Loaded OTTERS TWAS:", nrow(otters_all), "entries\n")

    # Standardize gene column
    gene_sym_col <- intersect(c("gene_symbol", "genename", "GeneName"), names(otters_all))[1]
    z_col <- intersect(c("otters_acat_z", "best_z", "Zscore", "TWAS_Z", "FUSION_Z", "z"), names(otters_all))[1]
    p_col <- intersect(c("otters_acat_pval", "acat_pval", "best_pval", "P", "pvalue", "FUSION_PVAL"), names(otters_all))[1]

    if (!is.na(gene_sym_col) && !is.na(p_col)) {
      # Best OTTERS result per gene (lowest p-value)
      setorderv(otters_all, p_col)
      otters_best <- otters_all[, .SD[1], by = gene_sym_col]

      # Count GWAS where significant
      if ("fdr" %in% names(otters_all)) {
        n_gwas_sig <- otters_all[fdr < 0.05, .(otters_n_gwas_sig = uniqueN(gwas_source)),
                                 by = gene_sym_col]
      } else {
        n_gwas_sig <- data.table()
        n_gwas_sig[[gene_sym_col]] <- character(0)
        n_gwas_sig[, otters_n_gwas_sig := integer(0)]
      }

      otters_dt <- otters_best[, .SD, .SDcols = c(gene_sym_col)]
      otters_dt[, otters_broadaway_pval := otters_best[[p_col]]]
      if (!is.na(z_col)) {
        otters_dt[, otters_broadaway_z := otters_best[[z_col]]]
      }

      if (nrow(n_gwas_sig) > 0) {
        otters_dt <- merge(otters_dt, n_gwas_sig, by = gene_sym_col, all.x = TRUE)
      } else {
        otters_dt[, otters_n_gwas_sig := NA_integer_]
      }

      if (gene_col != gene_sym_col) setnames(otters_dt, gene_sym_col, gene_col)
      new_cols <- intersect(c("otters_broadaway_pval", "otters_broadaway_z", "otters_n_gwas_sig"),
                            names(otters_dt))
      atlas <- safe_merge(atlas, otters_dt, new_cols = new_cols)
      cat("  Added", length(new_cols), "OTTERS columns\n")
    }
  }
} else {
  cat("  No OTTERS TWAS results found\n")
}

# ==============================================================================
# 4. cTWAS (Phase 2B)
# ==============================================================================
cat("\n--- 4. cTWAS [DROPPED 2026-05-30] ---\n")
# cTWAS dropped from the atlas (review A09#1/#6): 64_ctwas_gtex.R normalized PIPs by the
# genome-wide ABF sum (not per-region), so ctwas_pip was permanently < 0.002 — it never
# cleared the > 0.5 gate and contributed zero genes, while the methods text claimed cTWAS
# evidence. Script 64 is archived under archive/ctwas_dropped_2026-05-30/. Drop any stale
# cTWAS columns left over from a previous atlas build.
ctwas_stale <- intersect(c("ctwas_pip", "ctwas_pval"), names(atlas))
if (length(ctwas_stale) > 0) {
  atlas[, (ctwas_stale) := NULL]
  cat("  Dropped stale cTWAS columns:", paste(ctwas_stale, collapse = ", "), "\n")
} else {
  cat("  cTWAS not present in atlas (dropped) — nothing to remove\n")
}

# ==============================================================================
# 5. SuSiE-COLOC (Phase 3A)
# ==============================================================================
cat("\n--- 5. SuSiE-COLOC ---\n")

# Primary source: gene-level COLOC summary (ABF + real SuSiE when available)
# NOTE: coloc_susie_best_pp4 name preserved for downstream compat (13+ scripts).
# Sources from real coloc_best_susie_pp4 when present, ABF coloc_best_pp4 otherwise.
susie_gene_file <- file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
if (file.exists(susie_gene_file)) {
  susie_gene <- fread(susie_gene_file)
  has_real_susie <- "coloc_best_susie_pp4" %in% names(susie_gene)

  # Use real SuSiE PP.H4 when available, ABF fallback otherwise
  if (has_real_susie) {
    susie_best <- susie_gene[gene != "", .(
      coloc_susie_best_pp4  = fifelse(!is.na(coloc_best_susie_pp4),
                                       coloc_best_susie_pp4, coloc_best_pp4),
      coloc_susie_n_signals = coloc_n_gwas_h4_05,
      coloc_susie_best_gwas = fifelse(!is.na(coloc_best_susie_gwas),
                                       coloc_best_susie_gwas, coloc_best_gwas)
    ), by = .(gene)]
    cat("  Using real SuSiE PP.H4 with ABF fallback\n")
  } else {
    susie_best <- susie_gene[gene != "", .(
      coloc_susie_best_pp4  = coloc_best_pp4,
      coloc_susie_n_signals = coloc_n_gwas_h4_05,
      coloc_susie_best_gwas = coloc_best_gwas
    ), by = .(gene)]
    cat("  NOTE: No real SuSiE columns — using ABF as coloc_susie_best_pp4\n")
  }
  setnames(susie_best, "gene", gene_col)
  susie_best <- susie_best[!duplicated(get(gene_col))]
  atlas <- safe_merge(atlas, susie_best,
                      new_cols = c("coloc_susie_best_pp4", "coloc_susie_n_signals", "coloc_susie_best_gwas"))
  cat("  COLOC gene-level:", nrow(susie_best), "genes\n")
  cat("    PP4 > 0.5:", sum(susie_best$coloc_susie_best_pp4 > 0.5, na.rm = TRUE), "\n")
} else {
  # Fallback: per-GWAS SuSiE results (ABF-based)
  susie_files <- list.files(CAUSAL_DIR, pattern = "susie_coloc.*\\.csv",
                            recursive = TRUE, full.names = TRUE)
  if (length(susie_files) > 0) {
    susie_all <- rbindlist(lapply(susie_files, function(f) {
      tryCatch(fread(f), error = function(e) NULL)
    }), fill = TRUE)
    if (nrow(susie_all) > 0) {
      pp4_col <- intersect(c("coloc_best_pp4", "PP.H4.abf", "best_pp4", "coloc_pp4"), names(susie_all))[1]
      gene_sym <- intersect(c("gene_symbol", "gene", "GeneSymbol"), names(susie_all))[1]
      if (!is.na(pp4_col) && !is.na(gene_sym)) {
        susie_best <- susie_all[, .(
          coloc_susie_best_pp4 = max(get(pp4_col), na.rm = TRUE),
          coloc_susie_n_signals = .N
        ), by = gene_sym]
        setnames(susie_best, gene_sym, gene_col)
        atlas <- safe_merge(atlas, susie_best, new_cols = c("coloc_susie_best_pp4", "coloc_susie_n_signals"))
        cat("  Added SuSiE-COLOC from per-GWAS fallback\n")
      }
    }
  } else {
    cat("  No SuSiE-COLOC results found\n")
  }
}

# ==============================================================================
# 6. HyPrColoc / Multi-Trait COLOC (Phase 3B)
# ==============================================================================
cat("\n--- 6. Multi-Trait COLOC ---\n")

hypr_file <- file.path(CAUSAL_DIR, "hyprcoloc/hyprcoloc_results.csv")
# Also check checkpoint file (for in-progress runs)
hypr_ckpt <- file.path(CAUSAL_DIR, "hyprcoloc/hyprcoloc_results_checkpoint.csv")
if (file.exists(hypr_ckpt) && file.size(hypr_ckpt) > file.size(hypr_file)) {
  cat("  Using checkpoint file (job still running)\n")
  hypr_file <- hypr_ckpt
}
if (file.exists(hypr_file)) {
  hypr <- fread(hypr_file)
  cat("  Loaded HyPrColoc:", nrow(hypr), "genes\n")

  if (nrow(hypr) > 0 && "gene_symbol" %in% names(hypr)) {
    # Carry the method flag (review A08#2): distinguishes a real HyPrColoc cluster
    # posterior from the ad-hoc coloc_abf_gmean_fallback. Only the real method counts
    # toward L4 below. Default to fallback when the column is absent (older runs).
    hypr_method <- if ("method" %in% names(hypr)) hypr$method else "coloc_abf_gmean_fallback"
    hypr_dt <- hypr[, .(
      gene_symbol,
      hyprcoloc_posterior = hyprcoloc_posterior,
      hyprcoloc_n_traits = hyprcoloc_n_traits,
      hyprcoloc_traits = hyprcoloc_traits_list,
      hyprcoloc_method = hypr_method
    )]
    if (gene_col != "gene_symbol") setnames(hypr_dt, "gene_symbol", gene_col)
    atlas <- safe_merge(atlas, hypr_dt,
                        new_cols = c("hyprcoloc_posterior", "hyprcoloc_n_traits", "hyprcoloc_traits", "hyprcoloc_method"))
    cat("  Added HyPrColoc columns\n")
  } else {
    cat("  HyPrColoc results empty, skipping\n")
  }
} else {
  cat("  No HyPrColoc results found\n")
}

# ==============================================================================
# 7. Cell-type GWAS-signal enrichment (donor-correct MAGMA-equivalent arm)
# ==============================================================================
# REPOINTED 2026-06-20 (mega-review A6.4 / A8): the previously consumed arm,
# scdrs/<gwas>/group_results.csv, runs a per-CELL one-sample/welch t-test over
# 200K–657K cells. Cells from the same donor are not independent, so that test
# pseudoreplicates and returns top_celltype_pval = 0.0 / fdr = 0.0 (see the
# scdrs_summary diagnostic below) — an artefact, not a donor-level signal.
#
# The donor-correct, MAGMA.Celltyping-equivalent enrichment is produced by
# 201_celltype_heritability.R: it asks whether genes specifically expressed in
# each cell type are enriched for the canonical GWAS COLOC/TWAS evidence
# (gene-set Wilcoxon, with Bonferroni control across cell types). That arm
# REVERSES the per-cell "hepatocyte is the top enriched cell type" claim:
# under the donor-correct arm hepatocytes carry the WEAKEST COLOC
# fold-enrichment of all 11 native lineages (1.16x, rank 11/11); the leading
# lineages are immune/stromal (Circulating NK/NKT 1.40x, T cells 1.38x,
# Macrophages 1.32x). We therefore report the donor-correct values here.
cat("\n--- 7. Cell-type GWAS-signal enrichment ---\n")

# 7a. Per-cell scDRS summary — DIAGNOSTIC ONLY (pseudoreplicated, p inflated).
scdrs_summary <- file.path(CAUSAL_DIR, "scdrs/scdrs_summary.csv")
if (file.exists(scdrs_summary)) {
  scdrs <- fread(scdrs_summary)
  cat("  [diagnostic only — per-cell t-test pseudoreplicates over",
      "~200K–657K cells; p=0.0 is an artefact, NOT consumed]\n")
  for (i in seq_len(nrow(scdrs))) {
    row <- scdrs[i]
    cat(sprintf("    %s: %d sig cells (%.1f%%), per-cell top CT: %s (p=%.2e, INVALID)\n",
                row$gwas, row$n_sig_cells, row$frac_sig_cells * 100,
                row$top_celltype, row$top_celltype_pval))
  }
}

# 7b. Donor-correct MAGMA-equivalent enrichment (201_celltype_heritability.R).
herit_file <- file.path(BASE, "RNA-seq/results/gwas_rna_integration",
                        "celltype_heritability_results.csv")
if (file.exists(herit_file)) {
  herit <- fread(herit_file)
  # Restrict to native per-cell-type sets; drop progression-contrast (*_vs_*)
  # and deconvolution (*_bayesprism) rows so the reported enrichment is the
  # clean cell-type heritability arm (producer fix tracked under A8).
  herit_ct <- herit[!grepl("_vs_|bayesprism", cell_type)]
  cat("  Donor-correct MAGMA-equivalent COLOC enrichment:",
      nrow(herit_ct), "cell types\n")
  n_bonf_sig <- sum(herit_ct$wilcox_p_coloc_bonf < 0.05, na.rm = TRUE)
  cat("  Bonferroni-significant (wilcox_p_coloc_bonf<0.05):",
      n_bonf_sig, "cell types\n")
  setorder(herit_ct, -fold_enrichment_coloc)
  for (i in seq_len(nrow(herit_ct))) {
    row <- herit_ct[i]
    cat(sprintf("    %s: COLOC fold-enrich %.2fx (wilcox_p_coloc_bonf=%.2g)\n",
                row$cell_type, row$fold_enrichment_coloc, row$wilcox_p_coloc_bonf))
  }
  hep <- herit_ct[cell_type == "Hepatocytes"]
  if (nrow(hep) == 1) {
    hep_rank <- which(herit_ct$cell_type == "Hepatocytes")
    cat(sprintf(paste0("  NOTE: hepatocytes rank %d/%d by COLOC fold-enrichment ",
                       "(%.2fx, bonf p=%.2g) — the per-cell 'top hepatocyte' ",
                       "claim does NOT replicate donor-correct.\n"),
                hep_rank, nrow(herit_ct), hep$fold_enrichment_coloc,
                hep$wilcox_p_coloc_bonf))
  }
} else {
  cat("  WARNING: donor-correct heritability file not found:", herit_file, "\n")
}

# ==============================================================================
# 8. Derived columns
# ==============================================================================
cat("\n--- 8. Derived columns ---\n")

# causal_methods_sig: count of significant causal methods per gene (p < 0.05)
# MR columns (mr_ivw_pval, mr_presso_pval, mr_raps_pval, mr_bidirectional_pval)
# removed 2026-04-22 — MR ditched from paper.
causal_pval_cols <- intersect(
  c("otters_broadaway_pval", "ctwas_pval",
    "twas_pval_ghodsian", "twas_pval_chen"),
  names(atlas)
)

# Also count COLOC PP.H4 > 0.5 and SuSiE PP.H4 > 0.5
coloc_cols <- intersect(
  c("coloc_pp4_ukbb_alt", "coloc_pp4_ukbb_ast", "coloc_pp4_ukbb_ggt",
    "coloc_pp4_pdff", "coloc_susie_best_pp4", "hyprcoloc_posterior"),
  names(atlas)
)

if (length(causal_pval_cols) > 0 || length(coloc_cols) > 0) {
  atlas[, causal_methods_sig := 0L]

  for (pcol in causal_pval_cols) {
    atlas[!is.na(get(pcol)) & get(pcol) < 0.05, causal_methods_sig := causal_methods_sig + 1L]
  }
  for (ccol in coloc_cols) {
    atlas[!is.na(get(ccol)) & get(ccol) > 0.5, causal_methods_sig := causal_methods_sig + 1L]
  }

  atlas[, causal_robustness := fifelse(
    causal_methods_sig >= 3, "Robust",
    fifelse(causal_methods_sig >= 1, "Suggestive", "None")
  )]

  cat("  Causal robustness distribution:\n")
  print(table(atlas$causal_robustness))
  cat("\n")
}

# ==============================================================================
# 8b. Update layers_active to include new overhaul columns
# ==============================================================================
cat("\n--- 8b. Update L4 layer count ---\n")

# Recalculate layers_active if it exists (27a computes it but misses 75's new columns)
if ("layers_active" %in% names(atlas)) {
  old_l4 <- sum(atlas$layers_active >= 1, na.rm = TRUE)

  # Check if any new causal columns add L4 evidence for genes that had layers_active=0 or missing L4
  # (mr_ivw_pval, mr_presso_pval removed 2026-04-22 — MR ditched)
  new_l4_pval <- intersect(c("otters_broadaway_pval"), names(atlas))
  # cTWAS dropped (review A09#1); HyPrColoc gated on real-method only, below (review A08#2)
  new_l4_pp <- intersect(c("coloc_susie_best_pp4"), names(atlas))

  # Build a boolean mask: does this gene gain L4 from overhaul columns?
  new_l4_gain <- rep(FALSE, nrow(atlas))
  for (.col in new_l4_pval) {
    new_l4_gain <- new_l4_gain | (!is.na(atlas[[.col]]) & atlas[[.col]] < 0.05)
  }
  for (.col in new_l4_pp) {
    new_l4_gain <- new_l4_gain | (!is.na(atlas[[.col]]) & atlas[[.col]] > 0.5)
  }
  # HyPrColoc contributes to L4 only when it is a real HyPrColoc cluster posterior, not the
  # coloc_abf_gmean_fallback ad-hoc geometric mean (review A08#2). Absent method => not counted.
  if ("hyprcoloc_posterior" %in% names(atlas)) {
    hypr_real <- if ("hyprcoloc_method" %in% names(atlas)) atlas$hyprcoloc_method == "hyprcoloc" else FALSE
    new_l4_gain <- new_l4_gain | (!is.na(atlas$hyprcoloc_posterior) &
                                  atlas$hyprcoloc_posterior > 0.5 & hypr_real)
  }

  # For genes that didn't already have any L4 signal, add 1 to layers_active
  # (We check the original L4-related columns to avoid double-counting)
  # (mr_pval removed 2026-04-22 — MR ditched)
  orig_l4_cols_pval <- intersect(c("twas_pval"), names(atlas))
  # FinnGen columns removed 2026-04-08: GWAS archived (duplicate of Whitfield 2023)
  orig_l4_cols_pp4 <- intersect(c("coloc_pp4", "broadaway_coloc_pp4",
    "ukbb_alt_coloc_pp4", "ast_coloc_pp4", "ggt_coloc_pp4", "pdff_coloc_pp4",
    "sceqtl_coloc_best_pp4", "bbj_alt_coloc_pp4", "bbj_ast_coloc_pp4",
    "bbj_ggt_coloc_pp4"), names(atlas))

  had_orig_l4 <- rep(FALSE, nrow(atlas))
  for (.col in orig_l4_cols_pval) {
    had_orig_l4 <- had_orig_l4 | (!is.na(atlas[[.col]]) & atlas[[.col]] < 0.05)
  }
  for (.col in orig_l4_cols_pp4) {
    had_orig_l4 <- had_orig_l4 | (!is.na(atlas[[.col]]) & atlas[[.col]] > 0.5)
  }

  # Genes gaining L4 = new evidence AND didn't have original L4
  gains_l4 <- new_l4_gain & !had_orig_l4
  n_gains <- sum(gains_l4, na.rm = TRUE)
  atlas[gains_l4, layers_active := layers_active + 1L]

  cat(sprintf("  Genes gaining L4 from overhaul: %d\n", n_gains))
  cat(sprintf("  L4 coverage: was ~%d, now ~%d\n", old_l4,
              sum(had_orig_l4 | new_l4_gain, na.rm = TRUE)))
}

# ==============================================================================
# 9. Save updated atlas
# ==============================================================================
cat("--- 9. Save ---\n")
cat("Atlas dimensions:", nrow(atlas), "x", ncol(atlas), "\n")
cat("New columns added:", ncol(atlas) - n_col_start, "\n")

fwrite(atlas, ATLAS_FILE)
cat("Saved:", ATLAS_FILE, "\n")

# Summary of new S3 coverage (MR columns removed 2026-04-22 — MR ditched)
s3_cols <- intersect(
  c("otters_broadaway_pval", "ctwas_pip", "coloc_susie_best_pp4",
    "hyprcoloc_posterior"),
  names(atlas)
)

for (col in s3_cols) {
  n_non_na <- sum(!is.na(atlas[[col]]))
  pct <- round(100 * n_non_na / nrow(atlas), 1)
  cat("  ", col, ":", n_non_na, "(", pct, "%)\n")
}

cat("\n=== Integration complete ===\n")
cat("End:", format(Sys.time()), "\n")
