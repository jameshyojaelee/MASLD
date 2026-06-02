#!/usr/bin/env Rscript
# 35s_combine_susie_results.R
# ---------------------------------------------------------------------------
# Combine per-chromosome SuSiE-COLOC results into final output per GWAS.
# Run after all array jobs complete.
#
# Reads:  susie_coloc_results_chr{1-22}.csv from each GWAS results dir
# Writes: susie_coloc_results.csv (combined) + susie_summary.csv (aggregate)
#
# Usage:
#   Rscript RNA-seq/35s_combine_susie_results.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

BASE_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

GWAS_LIST <- c("UKBB_ALT", "UKBB_AST", "UKBB_GGT",
               "FINNGEN_NAFLD", "FINNGEN_NASH", "FINNGEN_HCC")

cat("=== Combining SuSiE-COLOC per-chromosome results ===\n")
cat("Time:", format(Sys.time()), "\n\n")

for (gwas_name in GWAS_LIST) {
  gwas_short <- tolower(gwas_name)
  results_dir <- file.path(BASE_DIR, "RNA-seq/results/causal_inference",
                           paste0("susie_broadaway_", gwas_short))

  cat("--- ", gwas_name, " ---\n")

  # Collect per-chromosome result files
  chr_files <- list()
  missing_chr <- integer(0)

  for (chr_num in 1:22) {
    # Try results file first, fall back to checkpoint
    res_file <- file.path(results_dir, sprintf("susie_coloc_results_chr%d.csv", chr_num))
    ckpt_file <- file.path(results_dir, sprintf("susie_coloc_checkpoint_chr%d.csv", chr_num))

    if (file.exists(res_file) && file.size(res_file) > 0) {
      dt <- fread(res_file)
      if (nrow(dt) > 0) {
        chr_files[[length(chr_files) + 1]] <- dt
        cat("  chr", chr_num, ":", nrow(dt), "genes (results)\n")
      } else {
        cat("  chr", chr_num, ": 0 genes (empty results)\n")
      }
    } else if (file.exists(ckpt_file) && file.size(ckpt_file) > 0) {
      dt <- fread(ckpt_file)
      if (nrow(dt) > 0) {
        chr_files[[length(chr_files) + 1]] <- dt
        cat("  chr", chr_num, ":", nrow(dt), "genes (checkpoint fallback)\n")
      } else {
        cat("  chr", chr_num, ": 0 genes (empty checkpoint)\n")
      }
    } else {
      missing_chr <- c(missing_chr, chr_num)
      cat("  chr", chr_num, ": MISSING\n")
    }
  }

  if (length(missing_chr) > 0) {
    cat("  WARNING: Missing chromosomes:", paste(missing_chr, collapse = ", "), "\n")
  }

  if (length(chr_files) == 0) {
    cat("  No results found for", gwas_name, ", skipping\n\n")
    next
  }

  # Combine all chromosomes
  combined <- rbindlist(chr_files, fill = TRUE)

  # Deduplicate (in case a gene appeared in both results and checkpoint)
  combined <- combined[!duplicated(gene)]

  setorder(combined, -PP.H4.abf)

  # Write combined results
  out_file <- file.path(results_dir, "susie_coloc_results.csv")
  fwrite(combined, out_file)

  # Summary stats
  n_total <- nrow(combined)
  susie_genes <- combined[!is.na(PP.H4.susie)]

  cat("\n  Combined:", n_total, "genes\n")
  cat("  ABF PP.H4 > 0.8:", nrow(combined[PP.H4.abf > 0.8]), "\n")
  cat("  ABF PP.H4 > 0.5:", nrow(combined[PP.H4.abf > 0.5]), "\n")
  cat("  ABF PP.H4 > 0.3:", nrow(combined[PP.H4.abf > 0.3]), "\n")

  if (nrow(susie_genes) > 0) {
    cat("  SuSiE PP.H4 > 0.8:", nrow(susie_genes[PP.H4.susie > 0.8]), "\n")
    cat("  SuSiE PP.H4 > 0.5:", nrow(susie_genes[PP.H4.susie > 0.5]), "\n")
    cat("  SuSiE PP.H4 > 0.3:", nrow(susie_genes[PP.H4.susie > 0.3]), "\n")

    if (nrow(susie_genes) > 10) {
      rho <- cor(susie_genes$PP.H4.abf, susie_genes$PP.H4.susie,
                 use = "complete.obs", method = "spearman")
      cat("  ABF vs SuSiE Spearman rho:", round(rho, 3), "\n")
    }
  }

  # Write aggregate summary
  summary_dt <- data.table(
    metric = c("genes_tested", "susie_available", "chromosomes_with_results",
               "chromosomes_missing",
               "ABF_PP.H4_gt_0.8", "ABF_PP.H4_gt_0.5", "ABF_PP.H4_gt_0.3",
               "SuSiE_PP.H4_gt_0.8", "SuSiE_PP.H4_gt_0.5", "SuSiE_PP.H4_gt_0.3"),
    value = c(n_total, nrow(susie_genes), 22 - length(missing_chr),
              length(missing_chr),
              nrow(combined[PP.H4.abf > 0.8]),
              nrow(combined[PP.H4.abf > 0.5]),
              nrow(combined[PP.H4.abf > 0.3]),
              nrow(susie_genes[PP.H4.susie > 0.8]),
              nrow(susie_genes[PP.H4.susie > 0.5]),
              nrow(susie_genes[PP.H4.susie > 0.3]))
  )
  fwrite(summary_dt, file.path(results_dir, "susie_summary.csv"))

  cat("  Saved:", out_file, "\n\n")
}

cat("=== Done ===\n")
cat("Time:", format(Sys.time()), "\n")
