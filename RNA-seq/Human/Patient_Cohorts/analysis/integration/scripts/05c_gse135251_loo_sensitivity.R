#!/usr/bin/env Rscript
# 05c_gse135251_loo_sensitivity.R
# ---------------------------------------------------------------------------
# GSE135251 leave-one-out sensitivity analysis for proteomics circularity.
#
# PROBLEM: GSE135251 (Govaere et al.) contributes ~160 samples to the
#   discovery dream mega-analysis AND appears in the Govaere proteomics
#   validation set (SomaScan). Any proteomics validation claims could
#   therefore be partially circular.
#
# THIS SCRIPT:
#   1. Reruns dream mega-analysis excluding GSE135251
#   2. Applies ashr adaptive shrinkage to the LOO results
#   3. Compares DEG results (full vs LOO): correlation, direction, Jaccard
#   4. Checks whether proteomics-validated genes survive LOO exclusion
#   5. Saves results to RNA-seq/results/audit_sensitivity/
#
# Runs AFTER: 05_dream_mega_analysis.R, 05b_ashr_shrinkage.R
# Pattern: based on dream_loo_cv.R + 05b_ashr_shrinkage.R
#
# Usage:
#   sbatch --job-name=gse135251_loo \
#     --partition=cpu --cpus-per-task=16 --mem=64G --time=48:00:00 \
#     --wrap="micromamba run -n rnaseq Rscript 05c_gse135251_loo_sensitivity.R"
# ---------------------------------------------------------------------------

#SBATCH --job-name=gse135251_loo
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/gse135251_loo_%j.out
#SBATCH --error=logs/gse135251_loo_%j.err

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(yaml)
  library(edgeR)
})

# Force injection into lme4 namespace BEFORE loading variancePartition
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({
      unlockBinding(fn, ns_lme4)
      assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
      lockBinding(fn, ns_lme4)
    }, silent = TRUE)
  }
}

suppressPackageStartupMessages({
  library(variancePartition)
  library(BiocParallel)
  library(ashr)
})

cat("============================================================\n")
cat("05c: GSE135251 LOO Sensitivity (Proteomics Circularity Audit)\n")
cat("============================================================\n")
cat("Time:", as.character(Sys.time()), "\n\n")

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")

# Output directory
OUT_DIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

HELD_OUT <- "GSE135251"

# =========================================================================
# 1. Load existing full-model results (ashr-shrunk)
# =========================================================================
cat("--- Loading full-model results ---\n")

full_ashr_f <- file.path(RDIR, "dream_results_ashr.csv")
full_raw_f  <- file.path(RDIR, "dream_results.csv")

if (!file.exists(full_ashr_f)) {
  stop("Full ashr results not found: ", full_ashr_f,
       "\nRun 05_dream_mega_analysis.R + 05b_ashr_shrinkage.R first.")
}
full <- fread(full_ashr_f)
cat("Full model:", nrow(full), "genes\n")
cat("Full model columns:", paste(names(full), collapse = ", "), "\n")

# Define full-model DEGs (standard thresholds)
PADJ_THRESH <- 0.05
LFC_THRESH  <- 0.5  # Migrated 0.3 -> 0.5 (LOO-CV stability; ~1.41x fold change)

full[, full_sig := padj < PADJ_THRESH & abs(logFC) > LFC_THRESH]
n_full_deg <- sum(full$full_sig, na.rm = TRUE)
cat("Full-model DEGs (padj<", PADJ_THRESH, ", |logFC|>", LFC_THRESH, "):",
    n_full_deg, "\n\n")

# =========================================================================
# 2. Check for pre-existing LOO result; if not, run dream
# =========================================================================
loo_raw_f <- file.path(RDIR, "loo_cv", paste0("dream_loo_", HELD_OUT, ".csv"))

if (file.exists(loo_raw_f)) {
  cat("--- Found pre-existing LOO dream result ---\n")
  cat("File:", loo_raw_f, "\n")
  loo_raw <- fread(loo_raw_f)
  cat("LOO raw genes:", nrow(loo_raw), "\n\n")
} else {
  cat("--- Running dream LOO (excluding", HELD_OUT, ") ---\n")

  # Load merged DGE
  dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

  # Exclude permanently excluded datasets + GSE135251
  ycfg <- yaml::read_yaml(file.path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design", "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
  keep_samples <- dge$samples$dataset %in% setdiff(mega_cohorts, HELD_OUT)
  dge_loo <- dge[, keep_samples]

  cat("Mega cohorts (yaml):", paste(mega_cohorts, collapse = ", "), "
")
  cat("Held out:", HELD_OUT, "\n")
  cat("Samples remaining:", ncol(dge_loo), "\n")
  cat("Genes:", nrow(dge_loo), "\n")

  # Load metadata for inferred sex
  meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
  stopifnot("All DGE samples must be in meta_matched" =
              all(colnames(dge_loo) %in% meta_new$sample_id))

  matched_sex <- meta_new$inferred_sex[match(colnames(dge_loo), meta_new$sample_id)]
  n_na_sex <- sum(is.na(matched_sex))
  if (n_na_sex > 0) {
    warning(n_na_sex, " samples have NA inferred_sex; these will be dropped by dream().")
  }

  info <- data.frame(
    group_binary = factor(dge_loo$samples$group_binary, levels = c("Control", "Disease")),
    dataset = droplevels(factor(dge_loo$samples$dataset)),
    inferred_sex = factor(matched_sex),
    stringsAsFactors = FALSE
  )
  rownames(info) <- colnames(dge_loo)

  cat("\nGroup distribution (LOO, excluding", HELD_OUT, "):\n")
  print(table(info$group_binary, info$dataset))
  cat("Sex distribution:\n")
  print(table(info$inferred_sex, useNA = "always"))

  # Dream model

  form <- ~ group_binary + inferred_sex + (1|dataset)
  cat("\nFormula:", deparse(form), "\n")

  ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
  cat("Using", ncpus, "CPU cores\n")
  param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

  cat("Running voomWithDreamWeights...\n")
  v <- suppressWarnings(voomWithDreamWeights(dge_loo, form, info, BPPARAM = param))

  cat("Running dream()...\n")
  fit <- suppressWarnings(dream(v, form, info, BPPARAM = param))
  # NOTE: do NOT call eBayes() after dream()

  res <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  loo_raw <- as.data.table(res)
  setnames(loo_raw, "adj.P.Val", "padj")

  # Save raw LOO result
  loo_dir <- file.path(RDIR, "loo_cv")
  dir.create(loo_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(loo_raw, loo_raw_f)
  cat("Saved raw LOO result:", loo_raw_f, "\n\n")
}

# =========================================================================
# 3. Apply ashr shrinkage to LOO results
# =========================================================================
cat("--- Applying ashr shrinkage to LOO results ---\n")

# Derive standard errors: SE = |logFC / t|
loo_raw[, se := abs(logFC / t)]
valid <- !is.na(loo_raw$se) & is.finite(loo_raw$se) & loo_raw$se > 0
n_bad <- sum(!valid)
if (n_bad > 0) {
  cat("WARNING:", n_bad, "genes with invalid SE — excluded from ashr\n")
}

cat("Running ashr on", sum(valid), "genes...\n")
ash_fit <- ash(
  betahat     = loo_raw$logFC[valid],
  sebetahat   = loo_raw$se[valid],
  mixcompdist = "halfuniform",
  method      = "shrink"
)

loo_raw[, shrunk_logFC := NA_real_]
loo_raw[, lfsr         := NA_real_]
loo_raw[valid, shrunk_logFC := ash_fit$result$PosteriorMean]
loo_raw[valid, lfsr         := ash_fit$result$lfsr]

# Define LOO DEGs with same thresholds
loo_raw[, loo_sig := padj < PADJ_THRESH & abs(logFC) > LFC_THRESH]
n_loo_deg <- sum(loo_raw$loo_sig, na.rm = TRUE)

cat("LOO DEGs (padj<", PADJ_THRESH, ", |logFC|>", LFC_THRESH, "):",
    n_loo_deg, "\n")
cat("  Up:", sum(loo_raw$loo_sig & loo_raw$logFC > 0, na.rm = TRUE), "\n")
cat("  Down:", sum(loo_raw$loo_sig & loo_raw$logFC < 0, na.rm = TRUE), "\n\n")

# =========================================================================
# 4. Compare full vs LOO: correlation, direction, Jaccard
# =========================================================================
cat("--- Comparing full model vs LOO (excluding", HELD_OUT, ") ---\n")

# Merge on gene
comp <- merge(
  full[, .(gene, full_logFC = logFC, full_shrunk = shrunk_logFC,
           full_lfsr = lfsr, full_padj = padj, full_sig)],
  loo_raw[, .(gene, loo_logFC = logFC, loo_shrunk = shrunk_logFC,
              loo_lfsr = lfsr, loo_padj = padj, loo_sig)],
  by = "gene"
)
cat("Genes in common:", nrow(comp), "\n")

# 4a. Spearman correlation (all genes)
rho_raw    <- cor(comp$full_logFC, comp$loo_logFC, method = "spearman", use = "complete.obs")
rho_shrunk <- cor(comp$full_shrunk, comp$loo_shrunk, method = "spearman", use = "complete.obs")
cat("Spearman rho (raw logFC):", round(rho_raw, 4), "\n")
cat("Spearman rho (shrunk logFC):", round(rho_shrunk, 4), "\n")

# 4b. Spearman correlation (among full-model DEGs only)
deg_comp <- comp[full_sig == TRUE]
rho_deg <- cor(deg_comp$full_shrunk, deg_comp$loo_shrunk, method = "spearman", use = "complete.obs")
cat("Spearman rho (shrunk, among full DEGs):", round(rho_deg, 4), "\n")

# 4c. Direction concordance
valid_dir <- !is.na(comp$full_shrunk) & !is.na(comp$loo_shrunk) &
             comp$full_shrunk != 0 & comp$loo_shrunk != 0
dir_match <- sum(sign(comp$full_shrunk[valid_dir]) == sign(comp$loo_shrunk[valid_dir]))
dir_total <- sum(valid_dir)
dir_pct   <- round(100 * dir_match / dir_total, 2)
cat("Direction concordance (all genes):", dir_match, "/", dir_total,
    "(", dir_pct, "%)\n")

# Direction concordance among DEGs
deg_dir <- deg_comp[!is.na(full_shrunk) & !is.na(loo_shrunk) &
                    full_shrunk != 0 & loo_shrunk != 0]
deg_dir_match <- sum(sign(deg_dir$full_shrunk) == sign(deg_dir$loo_shrunk))
deg_dir_pct   <- round(100 * deg_dir_match / nrow(deg_dir), 2)
cat("Direction concordance (DEGs only):", deg_dir_match, "/", nrow(deg_dir),
    "(", deg_dir_pct, "%)\n")

# 4d. Jaccard similarity of DEG sets
full_degs <- comp[full_sig == TRUE, gene]
loo_degs  <- comp[loo_sig == TRUE, gene]
jaccard   <- length(intersect(full_degs, loo_degs)) /
             length(union(full_degs, loo_degs))
cat("Jaccard similarity (DEG sets):", round(jaccard, 4), "\n")

# 4e. Recovery rate: fraction of full-model DEGs that remain DEGs in LOO
recovered      <- sum(full_degs %in% loo_degs)
recovery_rate  <- round(100 * recovered / length(full_degs), 2)
cat("Recovery rate:", recovered, "/", length(full_degs),
    "(", recovery_rate, "%)\n")

# 4f. DEGs gained in LOO (not in full model)
gained <- sum(loo_degs %in% full_degs == FALSE)
cat("DEGs gained in LOO:", gained, "\n")

# 4g. DEGs lost in LOO (in full but not LOO)
lost <- length(full_degs) - recovered
cat("DEGs lost in LOO:", lost, "\n\n")

# =========================================================================
# 5. Proteomics-specific validation
# =========================================================================
cat("--- Proteomics validation audit ---\n")

PROTEO_DIR <- file.path(BASE, "Analysis/Proteomics/results")
conc_f     <- file.path(PROTEO_DIR, "protein_transcript_concordance_v3.csv")

if (!file.exists(conc_f)) {
  cat("WARNING: Proteomics concordance file not found:", conc_f, "\n")
  cat("Skipping proteomics validation check.\n")
  proteo_audit <- data.table()
} else {
  conc <- fread(conc_f)
  cat("Loaded proteomics concordance:", nrow(conc), "gene-dataset entries\n")
  cat("Datasets:", paste(unique(conc$dataset), collapse = ", "), "\n")
  cat("Dream comparators:", paste(unique(conc$dream_comparator), collapse = ", "), "\n\n")

  # --- Gene symbol mapping ---
  # Full model uses Ensembl IDs; proteomics uses gene symbols.
  # Load the atlas for Ensembl -> symbol mapping.
  atlas_f <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
  if (file.exists(atlas_f)) {
    atlas_map <- fread(atlas_f, select = c("ensembl_id", "human_symbol"))
    atlas_map[, ensembl_clean := sub("\\..*", "", ensembl_id)]
    atlas_map <- unique(atlas_map[!is.na(human_symbol) & human_symbol != "",
                                  .(ensembl_clean, symbol = human_symbol)])
    cat("Gene map loaded:", nrow(atlas_map), "Ensembl->symbol\n")
  } else {
    cat("WARNING: Atlas not found for gene mapping. Trying ashr symbol column.\n")
    atlas_map <- NULL
  }

  # Map LOO results to symbols
  loo_mapped <- copy(loo_raw)
  if ("symbol" %in% names(full)) {
    # Use symbol column from full results (already mapped in 05b)
    sym_map <- unique(full[!is.na(symbol) & symbol != "", .(gene, symbol)])
    loo_mapped <- merge(loo_mapped, sym_map, by = "gene", all.x = TRUE)
  } else if (!is.null(atlas_map)) {
    loo_mapped[, ensembl_clean := sub("\\..*", "", gene)]
    loo_mapped <- merge(loo_mapped, atlas_map, by = "ensembl_clean", all.x = TRUE)
    loo_mapped[, ensembl_clean := NULL]
  } else {
    loo_mapped[, symbol := NA_character_]
  }

  # Also map full results if symbol not present
  full_mapped <- copy(full)
  if (!"symbol" %in% names(full_mapped)) {
    if (!is.null(atlas_map)) {
      full_mapped[, ensembl_clean := sub("\\..*", "", gene)]
      full_mapped <- merge(full_mapped, atlas_map, by = "ensembl_clean", all.x = TRUE)
      full_mapped[, ensembl_clean := NULL]
    }
  }

  # --- 5a. Disease vs control concordance: the key circularity check ---
  # These contrasts use dream_results.csv (Disease vs Control), which includes
  # GSE135251 in the full model. Check if concordance holds with LOO.
  cat("\n--- 5a. Disease-vs-Control concordance (affected by GSE135251 exclusion) ---\n")

  disease_conc <- conc[dream_comparator == "disease_vs_control"]
  if (nrow(disease_conc) > 0) {
    for (ds in unique(disease_conc$dataset)) {
      ds_conc <- disease_conc[dataset == ds]
      cat("\n  Dataset:", ds, "(", nrow(ds_conc), "genes in concordance table)\n")

      # Replace dream_logFC/dream_padj with LOO values  # C2-OK-sensitivity
      if (!is.null(loo_mapped) && "symbol" %in% names(loo_mapped)) {
        loo_slim <- loo_mapped[!is.na(symbol) & symbol != "",
                               .(symbol, loo_shrunk = shrunk_logFC, loo_lfsr = lfsr,
                                 loo_logFC = logFC, loo_padj = padj, loo_sig)]
        loo_slim <- loo_slim[!duplicated(symbol)]

        merged <- merge(ds_conc, loo_slim, by.x = "gene", by.y = "symbol", all.x = TRUE)

        n_mapped <- sum(!is.na(merged$loo_logFC))
        cat("    Mapped to LOO results:", n_mapped, "/", nrow(merged), "\n")

        if (n_mapped > 0) {
          mapped <- merged[!is.na(loo_logFC)]

          # Concordance with LOO logFC
          rho_loo <- cor(mapped$protein_logFC, mapped$loo_logFC,
                         method = "spearman", use = "complete.obs")
          # Concordance with original full-model logFC
          rho_full <- cor(mapped$protein_logFC, mapped$dream_logFC,  # C2-OK-sensitivity
                          method = "spearman", use = "complete.obs")

          dir_loo  <- sum(sign(mapped$protein_logFC) == sign(mapped$loo_logFC), na.rm = TRUE)
          dir_full <- sum(sign(mapped$protein_logFC) == sign(mapped$dream_logFC), na.rm = TRUE)  # C2-OK-sensitivity

          cat("    Protein-transcript Spearman rho (full model):", round(rho_full, 4), "\n")
          cat("    Protein-transcript Spearman rho (LOO, excl GSE135251):", round(rho_loo, 4), "\n")
          cat("    Direction concordance (full):", dir_full, "/", nrow(mapped),
              "(", round(100 * dir_full / nrow(mapped), 1), "%)\n")
          cat("    Direction concordance (LOO):", dir_loo, "/", nrow(mapped),
              "(", round(100 * dir_loo / nrow(mapped), 1), "%)\n")

          # Filtered concordance: |LFC| > 0.5 in both
          filt_full <- abs(mapped$protein_logFC) > 0.5 & abs(mapped$dream_logFC) > 0.5  # C2-OK-sensitivity
          filt_loo  <- abs(mapped$protein_logFC) > 0.5 & abs(mapped$loo_logFC) > 0.5

          if (sum(filt_full) > 0) {
            filt_full_conc <- sum(sign(mapped$protein_logFC[filt_full]) ==  # C2-OK-sensitivity
                                 sign(mapped$dream_logFC[filt_full]))  # C2-OK-sensitivity
            cat("    Filtered concordance (full, |LFC|>0.5 both):",
                filt_full_conc, "/", sum(filt_full),
                "(", round(100 * filt_full_conc / sum(filt_full), 1), "%)\n")
          }
          if (sum(filt_loo) > 0) {
            filt_loo_conc <- sum(sign(mapped$protein_logFC[filt_loo]) ==
                                sign(mapped$loo_logFC[filt_loo]))
            cat("    Filtered concordance (LOO, |LFC|>0.5 both):",
                filt_loo_conc, "/", sum(filt_loo),
                "(", round(100 * filt_loo_conc / sum(filt_loo), 1), "%)\n")
          }

          # Genes that were concordant-filtered in full model: do they survive LOO?
          full_concordant <- mapped[direction_concordant_filtered == TRUE]
          if (nrow(full_concordant) > 0) {
            loo_still_concordant <- sum(
              sign(full_concordant$protein_logFC) == sign(full_concordant$loo_logFC) &
              abs(full_concordant$loo_logFC) > 0.5,
              na.rm = TRUE
            )
            cat("    Genes concordant-filtered in full model:", nrow(full_concordant), "\n")
            cat("    Still concordant with LOO:", loo_still_concordant,
                "(", round(100 * loo_still_concordant / nrow(full_concordant), 1), "%)\n")
          }
        }
      } else {
        cat("    WARNING: Could not map gene symbols. Check atlas availability.\n")
      }
    }
  } else {
    cat("  No disease_vs_control entries in concordance table.\n")
  }

  # --- 5b. Other contrasts (fibrosis, NAFL vs NASH) ---
  # These use different dream comparators that are NOT from the main
  # disease vs control analysis, so the circularity is less direct.
  # Still worth noting for completeness.
  cat("\n--- 5b. Other contrasts (not directly affected by GSE135251 exclusion) ---\n")
  other_conc <- conc[dream_comparator != "disease_vs_control"]
  if (nrow(other_conc) > 0) {
    for (comp_label in unique(other_conc$dream_comparator)) {
      n_genes <- nrow(other_conc[dream_comparator == comp_label])
      n_conc  <- sum(other_conc[dream_comparator == comp_label, direction_concordant], na.rm = TRUE)
      cat("  ", comp_label, ":", n_genes, "genes,", n_conc, "direction-concordant",
          "(", round(100 * n_conc / n_genes, 1), "%)\n")
    }
    cat("  NOTE: These contrasts use disease-signature dream comparators, not the\n")
    cat("  main Disease vs Control mega-analysis. Circularity risk is lower.\n")
  }

  # --- 5c. DEG survival for proteomics-validated genes ---
  cat("\n--- 5c. Proteomics-validated gene survival in LOO ---\n")

  # Proteomics-validated = direction_concordant_filtered == TRUE in the
  # disease_vs_control contrast (the one affected by circularity)
  prot_validated <- disease_conc[direction_concordant_filtered == TRUE]
  cat("Proteomics-validated genes (concordant + |LFC|>0.5 both, disease_vs_control):",
      nrow(prot_validated), "\n")

  if (nrow(prot_validated) > 0 && !is.null(loo_mapped) && "symbol" %in% names(loo_mapped)) {
    loo_degs_sym <- loo_mapped[loo_sig == TRUE & !is.na(symbol), symbol]
    full_degs_sym <- full_mapped[full_sig == TRUE & !is.na(symbol), symbol]

    prot_in_full <- sum(prot_validated$gene %in% full_degs_sym)
    prot_in_loo  <- sum(prot_validated$gene %in% loo_degs_sym)

    cat("  In full-model DEGs:", prot_in_full, "/", nrow(prot_validated),
        "(", round(100 * prot_in_full / nrow(prot_validated), 1), "%)\n")
    cat("  In LOO DEGs:", prot_in_loo, "/", nrow(prot_validated),
        "(", round(100 * prot_in_loo / nrow(prot_validated), 1), "%)\n")

    # Detail: which genes are lost?
    lost_genes <- prot_validated$gene[prot_validated$gene %in% full_degs_sym &
                                     !(prot_validated$gene %in% loo_degs_sym)]
    if (length(lost_genes) > 0) {
      cat("  Proteomics-validated genes lost in LOO:", length(lost_genes), "\n")
      cat("    Genes:", paste(head(lost_genes, 20), collapse = ", "), "\n")
      if (length(lost_genes) > 20) cat("    ... and", length(lost_genes) - 20, "more\n")
    } else {
      cat("  No proteomics-validated genes lost in LOO. Validation is robust.\n")
    }

    # Also check: genes that become DEGs only in LOO (gained)
    gained_prot <- prot_validated$gene[!(prot_validated$gene %in% full_degs_sym) &
                                       prot_validated$gene %in% loo_degs_sym]
    if (length(gained_prot) > 0) {
      cat("  Proteomics-validated genes gained in LOO:", length(gained_prot), "\n")
    }
  }

  # --- 5d. Ranked enrichment comparison ---
  cat("\n--- 5d. Ranked enrichment: proteomics t-stats vs LOO dream ---\n")

  # For each proteomics dataset with disease_vs_control, check if dream DEG
  # enrichment in proteomics ranked list is preserved with LOO
  if (requireNamespace("fgsea", quietly = TRUE)) {
    library(fgsea)

    # Build gene sets from LOO and full model
    if (!is.null(loo_mapped) && "symbol" %in% names(loo_mapped)) {
      loo_deg_up   <- loo_mapped[loo_sig == TRUE & shrunk_logFC > 0 & !is.na(symbol), symbol]
      loo_deg_down <- loo_mapped[loo_sig == TRUE & shrunk_logFC < 0 & !is.na(symbol), symbol]
      loo_deg_all  <- loo_mapped[loo_sig == TRUE & !is.na(symbol), symbol]

      full_deg_up   <- full_mapped[full_sig == TRUE & shrunk_logFC > 0 & !is.na(symbol), symbol]
      full_deg_down <- full_mapped[full_sig == TRUE & shrunk_logFC < 0 & !is.na(symbol), symbol]
      full_deg_all  <- full_mapped[full_sig == TRUE & !is.na(symbol), symbol]

      gene_sets_full <- list(
        full_DEG_up   = full_deg_up,
        full_DEG_down = full_deg_down,
        full_DEG_all  = full_deg_all
      )
      gene_sets_loo <- list(
        loo_DEG_up   = loo_deg_up,
        loo_DEG_down = loo_deg_down,
        loo_DEG_all  = loo_deg_all
      )

      # Load proteomics DE results for disease_vs_control datasets
      prot_de_f <- file.path(PROTEO_DIR, "protein_differential_results_v3.csv")
      if (file.exists(prot_de_f)) {
        prot_de <- fread(prot_de_f)
        dvc_datasets <- c("PXD052937", "PXD051911")

        for (ds in dvc_datasets) {
          ds_de <- prot_de[dataset == ds & !is.na(t)]
          if (nrow(ds_de) == 0) next

          # Build ranked list
          ranked <- setNames(ds_de$t, ds_de$gene)
          ranked <- ranked[!duplicated(names(ranked))]
          ranked <- sort(ranked)

          if (length(ranked) < 100) next

          cat("\n  Proteomics dataset:", ds, "(", length(ranked), "proteins ranked)\n")

          # Test both full and LOO gene sets
          for (label in c("full", "loo")) {
            gs <- if (label == "full") gene_sets_full else gene_sets_loo
            gs_filt <- lapply(gs, function(g) intersect(g, names(ranked)))
            gs_filt <- gs_filt[sapply(gs_filt, length) >= 5]

            if (length(gs_filt) == 0) {
              cat("    ", label, ": insufficient overlap for fgsea\n")
              next
            }

            res_fgsea <- tryCatch({
              fgsea(pathways = gs_filt, stats = ranked, minSize = 5,
                    maxSize = 15000, nPermSimple = 10000)
            }, error = function(e) NULL)

            if (!is.null(res_fgsea)) {
              for (i in seq_len(nrow(res_fgsea))) {
                cat(sprintf("    %s %s: NES=%.2f, padj=%.2e (size=%d)\n",
                    label, res_fgsea$pathway[i], res_fgsea$NES[i],
                    res_fgsea$padj[i], res_fgsea$size[i]))
              }
            }
          }
        }
      }
    }
  } else {
    cat("  fgsea not available. Skipping ranked enrichment.\n")
  }
}

# =========================================================================
# 6. Save consolidated results
# =========================================================================
cat("\n--- Saving results ---\n")

# Main comparison table
summary_dt <- data.table(
  metric = c(
    "held_out", "n_full_degs", "n_loo_degs",
    "spearman_rho_raw_logFC", "spearman_rho_shrunk_logFC",
    "spearman_rho_shrunk_among_degs",
    "direction_concordance_all_pct", "direction_concordance_degs_pct",
    "jaccard_similarity", "recovery_rate_pct",
    "degs_gained_in_loo", "degs_lost_in_loo",
    "lfsr_threshold", "shrunk_lfc_threshold"
  ),
  value = c(
    HELD_OUT, n_full_deg, n_loo_deg,
    round(rho_raw, 4), round(rho_shrunk, 4),
    round(rho_deg, 4),
    dir_pct, deg_dir_pct,
    round(jaccard, 4), recovery_rate,
    gained, lost,
    PADJ_THRESH, LFC_THRESH
  )
)

out_f <- file.path(OUT_DIR, "gse135251_loo_sensitivity.csv")
fwrite(summary_dt, out_f)
cat("Saved summary:", out_f, "\n")

# Per-gene comparison table
per_gene <- comp[, .(gene,
                     full_shrunk_logFC = full_shrunk,
                     loo_shrunk_logFC  = loo_shrunk,
                     full_lfsr = full_lfsr,
                     loo_lfsr  = loo_lfsr,
                     full_deg  = full_sig,
                     loo_deg   = loo_sig,
                     direction_same = sign(full_shrunk) == sign(loo_shrunk),
                     status = fifelse(full_sig & loo_sig, "retained",
                              fifelse(full_sig & !loo_sig, "lost",
                              fifelse(!full_sig & loo_sig, "gained", "not_deg"))))]

per_gene_f <- file.path(OUT_DIR, "gse135251_loo_per_gene.csv")
fwrite(per_gene, per_gene_f)
cat("Saved per-gene comparison:", per_gene_f, "\n")

# =========================================================================
# 7. Print final summary
# =========================================================================
cat("\n")
cat("============================================================\n")
cat("        GSE135251 LOO SENSITIVITY — FINAL SUMMARY           \n")
cat("============================================================\n")
cat("\n")
cat("QUESTION: Is the proteomics validation circular because\n")
cat("  GSE135251 appears in both discovery and validation?\n")
cat("\n")
cat("DREAM MEGA-ANALYSIS ROBUSTNESS:\n")
cat("  Full model DEGs:       ", n_full_deg, "\n")
cat("  LOO DEGs (excl GSE135251):", n_loo_deg, "\n")
cat("  Spearman rho (shrunk):  ", round(rho_shrunk, 4), "\n")
cat("  Direction concordance:   ", deg_dir_pct, "%\n")
cat("  DEG recovery rate:       ", recovery_rate, "%\n")
cat("  Jaccard similarity:      ", round(jaccard, 4), "\n")
cat("\n")
cat("STATUS BREAKDOWN:\n")
cat("  Retained (DEG in both):", sum(per_gene$status == "retained"), "\n")
cat("  Lost (full only):      ", sum(per_gene$status == "lost"), "\n")
cat("  Gained (LOO only):     ", sum(per_gene$status == "gained"), "\n")
cat("\n")

if (recovery_rate >= 85) {
  cat("CONCLUSION: GSE135251 exclusion has MINIMAL impact on the dream\n")
  cat("  mega-analysis (recovery >=85%). Proteomics validation claims\n")
  cat("  are robust to this potential circularity.\n")
} else if (recovery_rate >= 70) {
  cat("CONCLUSION: GSE135251 exclusion has MODERATE impact. Core DEGs\n")
  cat("  are largely preserved but a substantial fraction is affected.\n")
  cat("  Proteomics validation claims should note this sensitivity.\n")
} else {
  cat("CONCLUSION: GSE135251 exclusion has SUBSTANTIAL impact. DEG\n")
  cat("  recovery is below 70%, suggesting the mega-analysis may be\n")
  cat("  over-reliant on this cohort. Proteomics validation claims\n")
  cat("  require careful qualification.\n")
}

cat("\n")
cat("REFERENCE: Existing LOO-CV summary shows GSE135251 recovery =",
    "89.4% (from loo_cv_summary.csv, padj<0.1 threshold)\n")
cat("============================================================\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("Done!\n")
