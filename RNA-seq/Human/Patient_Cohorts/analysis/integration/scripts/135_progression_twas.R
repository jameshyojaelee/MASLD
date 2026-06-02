#!/usr/bin/env Rscript
# 135_progression_twas.R
# ---------------------------------------------------------------------------
# Progression-specific TWAS re-analysis
#
# The existing TWAS pipeline (Script 19) runs S-PrediXcan to identify genes
# whose genetically-predicted expression is associated with MASLD risk.
# That analysis uses Disease-vs-Control contrasts implicitly (GWAS compares
# cases to controls).
#
# This script cross-references S-PrediXcan TWAS results with progression-
# specific dream DE statistics to identify genes that are BOTH:
#   (a) causally linked to MASLD via GWAS (TWAS z-score)
#   (b) differentially expressed during disease progression (dream t-stat)
#
# For each GWAS, the script either loads existing S-PrediXcan results or
# runs S-PrediXcan to generate them. Then for each progression contrast,
# it merges TWAS gene associations with progression t-statistics and
# computes a combined progression-causal score.
#
# Progression contrasts:
#   C2: NAFL vs NASH (from disease_signatures/nafl_vs_nash_dream.csv)
#   C3: Advanced vs Early Fibrosis (from progression/c3_adv_vs_early_fib_dream.csv
#       OR disease_signatures/adv_vs_early_fibrosis_dream.csv)
#   C5: NAS >= 5 vs NAS < 5 (from progression/c5_nas_ge5_vs_lt5_dream.csv)
#
# GWAS:
#   ghodsian, chen, finngen_nafld, finngen_nash, bbj_alt, bbj_ast, bbj_ggt
#
# Combined scoring:
#   - Fisher's combined p-value from TWAS p and progression p
#   - Product z-score = sign-concordant product of |TWAS z| * |progression t|
#   - Direction concordance flag (TWAS effect vs DE direction)
#
# Output: results/progression/progression_twas_{contrast}_{gwas}.csv
#         results/progression/progression_twas_summary.csv
#
# Usage:
#   Rscript 135_progression_twas.R
#   # Or run for a single GWAS:
#   GWAS_FILTER=ghodsian Rscript 135_progression_twas.R
#
# SLURM: bigmem partition recommended if S-PrediXcan needs to be run
#        (GWAS loading can use ~10-50 GB). If all TWAS results already
#        exist, cpu partition with 16 GB is sufficient.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(RSQLite)
})

cat("=== Script 135: Progression-Specific TWAS Re-Analysis ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# Configuration
# ==============================================================================
BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
setwd(BASE_DIR)

INT_DIR     <- file.path(BASE_DIR, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RESULTS_DIR <- file.path(INT_DIR, "results/progression")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

GWAS_DIR     <- file.path(BASE_DIR, "GWAS/MR_Data")
CAUSAL_DIR   <- file.path(BASE_DIR, "RNA-seq/results/causal_inference")
METAXCAN_DIR <- file.path(BASE_DIR, "MetaXcan/software")
SPREDIXCAN   <- file.path(METAXCAN_DIR, "SPrediXcan.py")

# PredictDB elastic net model (same as Script 19)
PREDICT_DB  <- file.path(GWAS_DIR, "en_Liver.db")
PREDICT_COV <- file.path(GWAS_DIR, "en_Liver.txt.gz")

# Gene annotation cache
GENE_CACHE <- file.path(INT_DIR, "results/gene_annotation/human_ensg_to_symbol.tsv")

# Optional: filter to a single GWAS for debugging
GWAS_FILTER <- Sys.getenv("GWAS_FILTER", "ALL")

# ==============================================================================
# GWAS registry
# ==============================================================================
# Each entry has:
#   file: path to S-PrediXcan formatted GWAS file
#   N: total sample size
#   twas_exists: whether Script 19 already produced TWAS results
#   skip_twas: whether S-PrediXcan is impossible (e.g., missing SE)

gwas_registry <- list(
  ghodsian = list(
    file = file.path(GWAS_DIR, "Ghodsian_2021_for_spredixcan.csv.gz"),
    N    = 778614L,
    label = "Ghodsian NAFLD (EHR, N=778K)"
  ),
  chen = list(
    file = file.path(GWAS_DIR, "Chen_2023_for_spredixcan.csv.gz"),
    N    = 691479L,
    label = "Chen NAFLD (meta-analysis, N=691K)"
  ),
  finngen_nafld = list(
    file = file.path(GWAS_DIR, "FinnGen/FinnGen_NAFLD_for_spredixcan.csv.gz"),
    N    = 438857L,
    label = "FinnGen NAFLD (R12, N=439K)"
  ),
  finngen_nash = list(
    file = file.path(GWAS_DIR, "FinnGen/FinnGen_NASH_for_spredixcan.csv.gz"),
    N    = 436066L,
    label = "FinnGen NASH (R12, N=436K)"
  ),
  bbj_alt = list(
    file = file.path(GWAS_DIR, "BBJ/BBJ_ALT_for_spredixcan.csv.gz"),
    N    = 261406L,
    label = "BBJ ALT (East Asian, N=261K)"
  ),
  bbj_ast = list(
    file = file.path(GWAS_DIR, "BBJ/BBJ_AST_for_spredixcan.csv.gz"),
    N    = 261270L,
    label = "BBJ AST (East Asian, N=261K)"
  ),
  bbj_ggt = list(
    file = file.path(GWAS_DIR, "BBJ/BBJ_GGT_for_spredixcan.csv.gz"),
    N    = 164006L,
    label = "BBJ GGT (East Asian, N=164K)"
  )
)

# Apply GWAS filter
if (GWAS_FILTER != "ALL") {
  gwas_names <- strsplit(GWAS_FILTER, ",")[[1]]
  gwas_registry <- gwas_registry[gwas_names[gwas_names %in% names(gwas_registry)]]
  if (length(gwas_registry) == 0) {
    stop("No valid GWAS found in GWAS_FILTER='", GWAS_FILTER,
         "'. Available: ", paste(names(gwas_registry), collapse = ", "))
  }
}

cat("GWAS to process:", paste(names(gwas_registry), collapse = ", "), "\n")

# ==============================================================================
# Progression contrast registry
# ==============================================================================
# Each contrast has a primary path (Script 130 output) and fallback
# (disease_signatures from older scripts). The t column is the dream
# moderated t-statistic, which serves as the progression z-score.

contrast_registry <- list(
  c2_nafl_vs_nash = list(
    paths = c(
      file.path(INT_DIR, "results/disease_signatures/nafl_vs_nash_dream.csv")
    ),
    label       = "NASH vs NAFL",
    description = "Steatohepatitis transition (positive logFC = UP in NASH)"
  ),
  c3_adv_vs_early_fib = list(
    paths = c(
      file.path(INT_DIR, "results/progression/c3_adv_vs_early_fib_dream.csv"),
      file.path(INT_DIR, "results/disease_signatures/adv_vs_early_fibrosis_dream.csv")
    ),
    label       = "Advanced vs Early Fibrosis",
    description = "F3-F4 vs F0-F2 fibrosis progression"
  ),
  c5_nas_ge5_vs_lt5 = list(
    paths = c(
      file.path(INT_DIR, "results/progression/c5_nas_ge5_vs_lt5_dream.csv")
    ),
    label       = "NAS >= 5 vs NAS < 5",
    description = "Clinical NASH threshold"
  )
)

# ==============================================================================
# Helper: Load gene symbol map
# ==============================================================================
load_gene_symbols <- function() {
  if (file.exists(GENE_CACHE)) {
    gm <- fread(GENE_CACHE)
    if ("gene_base" %in% names(gm)) {
      return(gm[!duplicated(gene_base), .(gene = gene_base, symbol)])
    }
  }
  cat("  WARNING: Gene annotation cache not found at", GENE_CACHE, "\n")
  return(NULL)
}

# ==============================================================================
# Helper: Load or run S-PrediXcan TWAS for a GWAS
# ==============================================================================
load_or_run_twas <- function(gwas_name, gwas_cfg) {
  # Check for existing results from Script 19
  existing_file <- file.path(CAUSAL_DIR, gwas_name, "twas_spredixcan_liver.csv")
  if (file.exists(existing_file) && file.size(existing_file) > 500) {
    cat("  Loading existing TWAS:", existing_file, "\n")
    res <- fread(existing_file)
    cat("    ", nrow(res), "genes loaded\n")
    return(res)
  }

  # Check if the GWAS file exists for S-PrediXcan
  if (!file.exists(gwas_cfg$file)) {
    cat("  WARNING: GWAS file not found:", gwas_cfg$file, "\n")
    cat("  Skipping TWAS for", gwas_name, "\n")
    return(NULL)
  }

  # Check if S-PrediXcan is available
  if (!file.exists(SPREDIXCAN)) {
    cat("  WARNING: S-PrediXcan not found:", SPREDIXCAN, "\n")
    return(NULL)
  }
  if (!file.exists(PREDICT_DB)) {
    cat("  WARNING: PredictDB not found:", PREDICT_DB, "\n")
    return(NULL)
  }

  # Run S-PrediXcan
  twas_out_dir <- file.path(CAUSAL_DIR, gwas_name)
  dir.create(twas_out_dir, recursive = TRUE, showWarnings = FALSE)
  twas_out <- file.path(twas_out_dir, "twas_spredixcan_liver.csv")

  cat("  Running S-PrediXcan for", gwas_name, "...\n")
  cmd <- paste(
    "micromamba run -n spatial python", shQuote(SPREDIXCAN),
    "--model_db_path", shQuote(PREDICT_DB),
    "--covariance", shQuote(PREDICT_COV),
    "--gwas_file", shQuote(gwas_cfg$file),
    "--snp_column snp",
    "--effect_allele_column effect_allele",
    "--non_effect_allele_column non_effect_allele",
    "--beta_column beta",
    "--se_column standard_error",
    "--pvalue_column pvalue",
    "--separator ','",
    "--gwas_N", gwas_cfg$N,
    "--output_file", shQuote(twas_out),
    "--overwrite --remove_ens_version --additional_output"
  )

  exit_code <- system(cmd)
  if (exit_code != 0) {
    cat("  WARNING: S-PrediXcan exited with code", exit_code, "\n")
  }

  if (file.exists(twas_out) && file.size(twas_out) > 500) {
    res <- fread(twas_out)
    res[, fdr := p.adjust(pvalue, method = "fdr")]
    fwrite(res, twas_out)
    cat("    TWAS complete:", nrow(res), "genes\n")
    cat("    TWAS FDR < 0.05:", nrow(res[fdr < 0.05]), "\n")
    return(res)
  }
  cat("  TWAS output not found or empty after S-PrediXcan run.\n")
  return(NULL)
}

# ==============================================================================
# Helper: Load progression dream results for a contrast
# ==============================================================================
load_progression_contrast <- function(contrast_id, contrast_cfg) {
  for (p in contrast_cfg$paths) {
    if (file.exists(p)) {
      cat("  Loading:", basename(p), "\n")
      dt <- fread(p)

      # Standardize column names
      # Some files have "adj.P.Val" instead of "padj"
      if ("adj.P.Val" %in% names(dt) && !"padj" %in% names(dt)) {
        setnames(dt, "adj.P.Val", "padj")
      }
      # Some files have "P.Value" instead of "pvalue"
      if ("P.Value" %in% names(dt) && !"pvalue" %in% names(dt)) {
        setnames(dt, "P.Value", "pvalue")
      }

      # Ensure gene column exists and strip Ensembl version
      if (!"gene" %in% names(dt)) {
        stop("  ERROR: No 'gene' column in ", p)
      }
      dt[, gene_base := sub("\\..*", "", gene)]

      cat("    ", nrow(dt), "genes,",
          sum(dt$padj < 0.1, na.rm = TRUE), "DEGs (padj<0.1)\n")
      return(dt)
    }
  }
  cat("  WARNING: No dream results found for", contrast_id, "\n")
  cat("    Searched:", paste(contrast_cfg$paths, collapse = "\n             "), "\n")
  cat("    Run Script 130 (progression_contrasts_dream.R) first.\n")
  return(NULL)
}

# ==============================================================================
# Helper: Fisher's combined p-value (two p-values)
# ==============================================================================
fisher_combine_p <- function(p1, p2) {
  # Fisher's method: -2 * sum(log(p)) ~ chi2(2k)
  # For exactly 2 p-values: df = 4
  valid <- !is.na(p1) & !is.na(p2) & p1 > 0 & p2 > 0
  result <- rep(NA_real_, length(p1))
  if (any(valid)) {
    chi2 <- -2 * (log(p1[valid]) + log(p2[valid]))
    result[valid] <- pchisq(chi2, df = 4, lower.tail = FALSE)
  }
  return(result)
}

# ==============================================================================
# Helper: Compute signed combined z-score
# ==============================================================================
# This is the Stouffer-like combination: z_combined = (z1 + z2) / sqrt(2)
# where we sign-align z1 to match the expected direction
signed_combined_z <- function(twas_z, prog_t) {
  valid <- !is.na(twas_z) & !is.na(prog_t)
  result <- rep(NA_real_, length(twas_z))
  if (any(valid)) {
    # TWAS z > 0 means genetically-upregulated expression increases MASLD risk
    # prog t > 0 means gene is upregulated during progression
    # Concordance: both positive or both negative
    result[valid] <- (twas_z[valid] + prog_t[valid]) / sqrt(2)
  }
  return(result)
}

# ==============================================================================
# MAIN: Load TWAS results for all GWAS
# ==============================================================================
cat("\n--- Step 1: Loading/running TWAS for each GWAS ---\n\n")

twas_results <- list()
for (gwas_name in names(gwas_registry)) {
  cat(sprintf("[GWAS: %s] %s\n", gwas_name, gwas_registry[[gwas_name]]$label))
  twas_results[[gwas_name]] <- load_or_run_twas(gwas_name, gwas_registry[[gwas_name]])
  cat("\n")
}

# Check how many GWAS have results
n_twas <- sum(!sapply(twas_results, is.null))
cat(sprintf("TWAS results available for %d / %d GWAS\n\n", n_twas, length(gwas_registry)))

if (n_twas == 0) {
  cat("ERROR: No TWAS results available. Cannot proceed.\n")
  cat("Run Script 19 first with GWAS_NAME=ghodsian to generate TWAS results.\n")
  quit(status = 1)
}

# ==============================================================================
# MAIN: Load progression contrast results
# ==============================================================================
cat("--- Step 2: Loading progression dream results ---\n\n")

prog_results <- list()
for (cid in names(contrast_registry)) {
  cat(sprintf("[Contrast: %s] %s\n", cid, contrast_registry[[cid]]$label))
  prog_results[[cid]] <- load_progression_contrast(cid, contrast_registry[[cid]])
  cat("\n")
}

n_contrasts <- sum(!sapply(prog_results, is.null))
cat(sprintf("Progression contrasts available: %d / %d\n\n",
            n_contrasts, length(contrast_registry)))

if (n_contrasts == 0) {
  cat("ERROR: No progression contrast results available.\n")
  cat("Run Script 130 (progression_contrasts_dream.R) and/or Script 13 first.\n")
  quit(status = 1)
}

# ==============================================================================
# Load gene symbols
# ==============================================================================
gene_symbols <- load_gene_symbols()

# ==============================================================================
# Load PredictDB gene metadata for annotation
# ==============================================================================
cat("--- Loading PredictDB model metadata ---\n")
if (file.exists(PREDICT_DB)) {
  con <- dbConnect(SQLite(), PREDICT_DB)
  db_extra <- as.data.table(dbReadTable(con, "extra"))
  dbDisconnect(con)
  # Strip Ensembl version from gene IDs
  db_extra[, gene_base := sub("\\..*", "", gene)]
  db_genes <- db_extra[, .(gene_base, genename, pred.perf.R2, pred.perf.pval,
                           n.snps.in.model)]
  cat("  PredictDB models:", nrow(db_genes), "genes\n\n")
} else {
  db_genes <- NULL
  cat("  WARNING: PredictDB not found; skipping model annotation\n\n")
}

# ==============================================================================
# MAIN: Cross-reference TWAS with progression contrasts
# ==============================================================================
cat("--- Step 3: Progression-TWAS cross-reference ---\n\n")

all_combined <- list()
summary_rows <- list()

for (cid in names(prog_results)) {
  prog <- prog_results[[cid]]
  if (is.null(prog)) next

  for (gwas_name in names(twas_results)) {
    twas <- twas_results[[gwas_name]]
    if (is.null(twas)) next

    cat(sprintf("  %s x %s: ", cid, gwas_name))

    # Prepare TWAS data: strip Ensembl version
    twas_dt <- copy(twas)
    twas_dt[, gene_base := sub("\\..*", "", gene)]

    # Select relevant TWAS columns
    twas_cols <- c("gene_base", "gene_name", "zscore", "effect_size", "pvalue",
                   "var_g", "pred_perf_r2", "pred_perf_pval", "n_snps_used",
                   "n_snps_in_model")
    twas_cols_avail <- intersect(twas_cols, names(twas_dt))
    twas_sub <- twas_dt[, ..twas_cols_avail]

    # Deduplicate by gene_base (keep first = most significant by default order)
    twas_sub <- twas_sub[!duplicated(gene_base)]

    # Prepare progression data
    prog_sub <- prog[, .(
      gene_base,
      prog_logFC = logFC,
      prog_t     = t,
      prog_pval  = pvalue,
      prog_padj  = padj,
      prog_aveexpr = AveExpr
    )]
    prog_sub <- prog_sub[!duplicated(gene_base)]

    # Merge on stripped Ensembl ID
    merged <- merge(twas_sub, prog_sub, by = "gene_base", all = FALSE)

    # Rename TWAS columns for clarity
    setnames(merged, c("zscore", "pvalue", "effect_size"),
             c("twas_z", "twas_pval", "twas_effect"),
             skip_absent = TRUE)

    if (nrow(merged) == 0) {
      cat("0 genes matched\n")
      next
    }

    # --- Compute combined scores ---

    # 1. Fisher's combined p-value
    merged[, fisher_pval := fisher_combine_p(twas_pval, prog_pval)]
    merged[, fisher_fdr  := p.adjust(fisher_pval, method = "fdr")]

    # 2. Signed combined z-score (Stouffer)
    merged[, combined_z := signed_combined_z(twas_z, prog_t)]
    merged[, combined_pval := 2 * pnorm(-abs(combined_z))]
    merged[, combined_fdr  := p.adjust(combined_pval, method = "fdr")]

    # 3. Direction concordance
    merged[, direction_concordant := sign(twas_z) == sign(prog_logFC)]

    # 4. Product score (geometric mean of -log10 p, sign-aware)
    merged[, product_score := sign(twas_z) * sign(prog_t) *
             sqrt(abs(twas_z) * abs(prog_t))]

    # 5. TWAS FDR (recalculate for this subset)
    merged[, twas_fdr := p.adjust(twas_pval, method = "fdr")]

    # --- Add gene symbol annotation ---
    if (!is.null(gene_symbols)) {
      merged <- merge(merged, gene_symbols, by.x = "gene_base", by.y = "gene",
                      all.x = TRUE)
      # Fill missing symbols from TWAS gene_name
      if ("gene_name" %in% names(merged)) {
        merged[is.na(symbol) | symbol == "", symbol := gene_name]
      }
    } else if ("gene_name" %in% names(merged)) {
      merged[, symbol := gene_name]
    }

    # --- Add PredictDB model quality ---
    if (!is.null(db_genes)) {
      merged <- merge(merged, db_genes[, .(gene_base, pred.perf.R2, n.snps.in.model)],
                      by = "gene_base", all.x = TRUE, suffixes = c("", ".db"))
    }

    # --- Add metadata columns ---
    merged[, contrast_id := cid]
    merged[, gwas_name := gwas_name]
    merged[, contrast_label := contrast_registry[[cid]]$label]
    gwas_lbl <- if (gwas_name %in% names(gwas_registry)) gwas_registry[[gwas_name]]$label else gwas_name
    merged[, gwas_label := gwas_lbl]

    # --- Sort by combined significance ---
    merged <- merged[order(combined_pval)]

    # --- Summary statistics ---
    n_total       <- nrow(merged)
    n_twas_sig    <- sum(merged$twas_fdr < 0.05, na.rm = TRUE)
    n_prog_sig    <- sum(merged$prog_padj < 0.1, na.rm = TRUE)
    n_both_sig    <- sum(merged$twas_fdr < 0.05 & merged$prog_padj < 0.1, na.rm = TRUE)
    n_fisher_sig  <- sum(merged$fisher_fdr < 0.05, na.rm = TRUE)
    n_combined_sig <- sum(merged$combined_fdr < 0.05, na.rm = TRUE)
    n_concordant  <- sum(merged$direction_concordant, na.rm = TRUE)
    pct_concordant <- round(100 * n_concordant / n_total, 1)

    cat(sprintf("%d genes, TWAS-sig=%d, Prog-sig=%d, Both=%d, ",
                n_total, n_twas_sig, n_prog_sig, n_both_sig))
    cat(sprintf("Fisher=%d, Combined=%d, Concordant=%s%%\n",
                n_fisher_sig, n_combined_sig, pct_concordant))

    # --- Save per-contrast-per-GWAS results ---
    out_file <- file.path(RESULTS_DIR,
      sprintf("progression_twas_%s_%s.csv", cid, gwas_name))
    fwrite(merged, out_file)

    # Store for combined analysis
    all_combined[[paste(cid, gwas_name, sep = "__")]] <- merged

    # Store summary row
    summary_rows[[length(summary_rows) + 1]] <- data.table(
      contrast_id    = cid,
      contrast_label = contrast_registry[[cid]]$label,
      gwas_name      = gwas_name,
      gwas_label     = gwas_lbl,
      n_genes_tested = n_total,
      n_twas_fdr05   = n_twas_sig,
      n_prog_padj01  = n_prog_sig,
      n_both_sig     = n_both_sig,
      n_fisher_fdr05 = n_fisher_sig,
      n_combined_fdr05 = n_combined_sig,
      pct_concordant = pct_concordant
    )
  }
}

# ==============================================================================
# Summary table
# ==============================================================================
cat("\n--- Step 4: Summary ---\n\n")

summary_dt <- rbindlist(summary_rows)
print(summary_dt)

fwrite(summary_dt, file.path(RESULTS_DIR, "progression_twas_summary.csv"))
cat("\nSaved: progression_twas_summary.csv\n")

# ==============================================================================
# Step 5: Cross-GWAS consensus for each contrast
# ==============================================================================
cat("\n--- Step 5: Cross-GWAS consensus analysis ---\n\n")

consensus_all <- list()

for (cid in names(prog_results)) {
  if (is.null(prog_results[[cid]])) next

  # Collect all GWAS results for this contrast
  cid_results <- list()
  for (gwas_name in names(twas_results)) {
    key <- paste(cid, gwas_name, sep = "__")
    if (key %in% names(all_combined)) {
      dt <- all_combined[[key]]
      cid_results[[gwas_name]] <- dt[, .(
        gene_base, symbol,
        twas_z, twas_pval, twas_fdr,
        prog_logFC, prog_t, prog_pval, prog_padj,
        fisher_fdr, combined_z, combined_fdr,
        direction_concordant
      )]
    }
  }

  if (length(cid_results) == 0) next

  cat(sprintf("  Consensus for %s across %d GWAS:\n", cid, length(cid_results)))

  # Count how many GWAS each gene is TWAS-significant in
  gene_counts <- rbindlist(lapply(names(cid_results), function(gn) {
    dt <- cid_results[[gn]]
    dt[, .(
      gene_base, symbol,
      twas_sig = twas_fdr < 0.05,
      fisher_sig = fisher_fdr < 0.05,
      combined_sig = combined_fdr < 0.05,
      concordant = direction_concordant,
      gwas = gn
    )]
  }))

  consensus <- gene_counts[, .(
    n_gwas_tested       = .N,
    n_gwas_twas_sig     = sum(twas_sig, na.rm = TRUE),
    n_gwas_fisher_sig   = sum(fisher_sig, na.rm = TRUE),
    n_gwas_combined_sig = sum(combined_sig, na.rm = TRUE),
    n_gwas_concordant   = sum(concordant, na.rm = TRUE),
    gwas_list_twas_sig  = paste(gwas[twas_sig], collapse = ";"),
    gwas_list_fisher_sig = paste(gwas[fisher_sig], collapse = ";")
  ), by = .(gene_base, symbol)]

  # Add progression statistics from the first available source
  first_key <- paste(cid, names(cid_results)[1], sep = "__")
  first_dt <- all_combined[[first_key]]
  prog_stats <- first_dt[, .(gene_base, prog_logFC, prog_t, prog_pval, prog_padj)]
  consensus <- merge(consensus, prog_stats, by = "gene_base", all.x = TRUE)

  # Sort by number of GWAS with significant combined evidence
  consensus <- consensus[order(-n_gwas_combined_sig, -n_gwas_twas_sig, prog_pval)]

  # Summary
  n_any_twas   <- sum(consensus$n_gwas_twas_sig > 0)
  n_multi_twas <- sum(consensus$n_gwas_twas_sig >= 2)
  n_any_comb   <- sum(consensus$n_gwas_combined_sig > 0)
  n_multi_comb <- sum(consensus$n_gwas_combined_sig >= 2)
  cat(sprintf("    Total genes: %d\n", nrow(consensus)))
  cat(sprintf("    TWAS-sig in >= 1 GWAS: %d\n", n_any_twas))
  cat(sprintf("    TWAS-sig in >= 2 GWAS: %d\n", n_multi_twas))
  cat(sprintf("    Combined-sig in >= 1 GWAS: %d\n", n_any_comb))
  cat(sprintf("    Combined-sig in >= 2 GWAS: %d\n", n_multi_comb))

  # Top hits
  top_hits <- consensus[n_gwas_combined_sig > 0][order(-n_gwas_combined_sig, prog_pval)]
  if (nrow(top_hits) > 0) {
    cat("    Top progression-causal genes:\n")
    show_n <- min(20, nrow(top_hits))
    for (i in seq_len(show_n)) {
      row <- top_hits[i]
      cat(sprintf("      %s (%s): %d/%d GWAS combined-sig, prog_logFC=%.3f, prog_padj=%.1e\n",
        row$symbol, row$gene_base,
        row$n_gwas_combined_sig, row$n_gwas_tested,
        row$prog_logFC, row$prog_padj))
    }
  }

  out_file <- file.path(RESULTS_DIR,
    sprintf("progression_twas_%s_consensus.csv", cid))
  fwrite(consensus, out_file)
  cat(sprintf("    Saved: %s\n\n", basename(out_file)))

  consensus[, contrast_id := cid]
  consensus_all[[cid]] <- consensus
}

# ==============================================================================
# Step 6: Unified progression-causal gene set
# ==============================================================================
cat("--- Step 6: Unified progression-causal gene set ---\n\n")

if (length(consensus_all) > 0) {
  # For each gene, count across how many contrasts it appears as progression-causal
  unified_list <- rbindlist(lapply(names(consensus_all), function(cid) {
    dt <- consensus_all[[cid]]
    dt[n_gwas_combined_sig > 0, .(
      gene_base, symbol, contrast_id = cid,
      n_gwas_combined_sig, prog_logFC, prog_padj
    )]
  }))

  if (nrow(unified_list) > 0) {
    unified <- unified_list[, .(
      n_contrasts_sig     = .N,
      contrasts           = paste(contrast_id, collapse = ";"),
      max_gwas_combined   = max(n_gwas_combined_sig),
      mean_prog_logFC     = mean(prog_logFC, na.rm = TRUE),
      min_prog_padj       = min(prog_padj, na.rm = TRUE),
      direction           = as.character(ifelse(all(prog_logFC > 0, na.rm = TRUE), "up",
                            ifelse(all(prog_logFC < 0, na.rm = TRUE), "down", "mixed")))
    ), by = .(gene_base, symbol)]

    unified <- unified[order(-n_contrasts_sig, -max_gwas_combined, min_prog_padj)]

    cat(sprintf("  Total progression-causal genes: %d\n", nrow(unified)))
    cat(sprintf("  Significant in >= 2 contrasts: %d\n",
                sum(unified$n_contrasts_sig >= 2)))
    cat(sprintf("  Significant in all %d contrasts: %d\n",
                n_contrasts, sum(unified$n_contrasts_sig == n_contrasts)))

    # Print top genes
    cat("\n  Top 20 progression-causal genes (by contrast replication):\n")
    show_n <- min(20, nrow(unified))
    for (i in seq_len(show_n)) {
      row <- unified[i]
      cat(sprintf("    %2d. %s: %d contrasts, max %d GWAS, %s, logFC=%.3f\n",
        i, row$symbol, row$n_contrasts_sig, row$max_gwas_combined,
        row$direction, row$mean_prog_logFC))
    }

    fwrite(unified, file.path(RESULTS_DIR, "progression_twas_unified.csv"))
    cat("\n  Saved: progression_twas_unified.csv\n")
  } else {
    cat("  No genes with combined significance across any contrast.\n")
  }
} else {
  cat("  No consensus results to unify.\n")
}

# ==============================================================================
# Step 7: Enrichment test — are TWAS-significant genes enriched for
#         progression DEGs? (hypergeometric test)
# ==============================================================================
cat("\n--- Step 7: Enrichment of TWAS-sig genes among progression DEGs ---\n\n")

enrichment_rows <- list()

for (cid in names(prog_results)) {
  prog <- prog_results[[cid]]
  if (is.null(prog)) next

  for (gwas_name in names(twas_results)) {
    twas <- twas_results[[gwas_name]]
    if (is.null(twas)) next

    # Build the contingency table
    twas_dt <- copy(twas)
    twas_dt[, gene_base := sub("\\..*", "", gene)]
    if (!"fdr" %in% names(twas_dt)) {
      twas_dt[, fdr := p.adjust(pvalue, method = "fdr")]
    }

    prog_sub <- prog[, .(gene_base, prog_padj = padj)]
    prog_sub <- prog_sub[!duplicated(gene_base)]

    # Universe = genes tested in both
    universe <- intersect(twas_dt$gene_base, prog_sub$gene_base)
    n_universe <- length(universe)

    twas_sig_genes <- twas_dt[gene_base %in% universe & fdr < 0.05, gene_base]
    prog_sig_genes <- prog_sub[gene_base %in% universe & prog_padj < 0.1, gene_base]

    n_twas_sig <- length(twas_sig_genes)
    n_prog_sig <- length(prog_sig_genes)
    n_overlap  <- length(intersect(twas_sig_genes, prog_sig_genes))

    # Hypergeometric test
    if (n_twas_sig > 0 && n_prog_sig > 0 && n_universe > 0) {
      hyper_pval <- phyper(n_overlap - 1, n_prog_sig,
                           n_universe - n_prog_sig,
                           n_twas_sig, lower.tail = FALSE)
      odds_ratio <- if (n_overlap > 0) {
        (n_overlap * (n_universe - n_twas_sig - n_prog_sig + n_overlap)) /
        ((n_twas_sig - n_overlap) * (n_prog_sig - n_overlap))
      } else 0

      enrichment_rows[[length(enrichment_rows) + 1]] <- data.table(
        contrast_id = cid,
        gwas_name   = gwas_name,
        n_universe  = n_universe,
        n_twas_sig  = n_twas_sig,
        n_prog_sig  = n_prog_sig,
        n_overlap   = n_overlap,
        odds_ratio  = round(odds_ratio, 2),
        hyper_pval  = hyper_pval,
        enriched    = hyper_pval < 0.05
      )
    }
  }
}

if (length(enrichment_rows) > 0) {
  enrichment_dt <- rbindlist(enrichment_rows)
  enrichment_dt[, hyper_fdr := p.adjust(hyper_pval, method = "fdr")]
  print(enrichment_dt)
  fwrite(enrichment_dt, file.path(RESULTS_DIR, "progression_twas_enrichment.csv"))
  cat("\nSaved: progression_twas_enrichment.csv\n")

  n_enriched <- sum(enrichment_dt$enriched)
  cat(sprintf("\n  Enriched (p < 0.05): %d / %d contrast-GWAS pairs\n",
              n_enriched, nrow(enrichment_dt)))
} else {
  cat("  No enrichment tests could be performed.\n")
}

# ==============================================================================
# Final report
# ==============================================================================
cat("\n")
cat(strrep("=", 70), "\n")
cat("  PROGRESSION TWAS ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n\n")

cat("Output files:\n")
out_files <- list.files(RESULTS_DIR, pattern = "^progression_twas_",
                        full.names = FALSE)
for (f in sort(out_files)) {
  cat(sprintf("  - %s\n", f))
}

cat(sprintf("\nTotal contrast x GWAS combinations: %d\n", length(all_combined)))
cat(sprintf("Total contrasts with consensus: %d\n", length(consensus_all)))

cat("\nEnd time:", format(Sys.time()), "\n")
cat("=== Script 135 complete ===\n")
