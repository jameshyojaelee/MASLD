#!/usr/bin/env Rscript
# 04e_aggregate_coloc_targeted.R
# Aggregator for the targeted COLOC fine-mapping re-run.
#
# Reads SuSiE + CARMA per-locus outputs from output/<study>_coloc_targeted/
# and produces results/combined_finemapping_coloc_targeted.csv with the same
# canonical schema as combined_finemapping.csv (modulo SuSiEX/MESuSiE columns,
# which are joint cross-ancestry and not computed for the targeted re-run).
#
# Schema written:
#   chromosome, position, allele1, allele2, locus, study, ancestry,
#   susie_pip, susie_cs, susie_converged,
#   carma_pip, carma_cs, carma_outlier,
#   susie_reliable, susie_pip_clean,
#   concordant_pip, max_pip, both_in_cs, either_in_cs, variant_id,
#   recommended_pip
#
# IMPORTANT: depends on CARMA .txt.gz files being non-empty. If they are 0-byte
# (data.table::fwrite zlib bug), run src/04d_reexport_carma_txt.R first.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

BASE     <- Sys.getenv("MASLD_PROJECT_ROOT",
                       unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR   <- file.path(BASE, "GWAS/finemapping")
LEAD_DIR <- file.path(FM_DIR, "data/lead_snps_coloc_targeted")
OUT_PATH <- file.path(FM_DIR, "results/combined_finemapping_coloc_targeted.csv")

if (!dir.exists(LEAD_DIR)) {
  stop("LEAD_DIR not found: ", LEAD_DIR)
}

REGISTRY <- fread(file.path(FM_DIR, "config/gwas_registry.tsv"))

lead_files <- list.files(LEAD_DIR, pattern = "_leadSNPs\\.tsv$", full.names = TRUE)
if (length(lead_files) == 0) stop("No lead-SNP files found in ", LEAD_DIR)

cat("============================================================\n")
cat("Aggregating COLOC-targeted fine-mapping results\n")
cat("============================================================\n")

all_susie <- list()
all_carma <- list()

for (lf in lead_files) {
  study <- sub("_leadSNPs\\.tsv$", "", basename(lf))
  reg_row <- REGISTRY[study_name == study][1]
  if (nrow(reg_row) == 0) {
    cat(sprintf("[skip] %s : not in gwas_registry.tsv\n", study))
    next
  }
  ld_pop     <- toupper(reg_row$ancestry)
  window_mb  <- reg_row$window_mb
  ancestry   <- reg_row$ancestry
  trait_type <- reg_row$trait_type
  out_base   <- file.path(FM_DIR, "output",
                          paste0(study, "_coloc_targeted"),
                          paste0(ld_pop, "_", window_mb, "Mb"))

  cat("\n--- Study:", study, "---\n")

  # ----- SuSiE -----
  susie_dir <- file.path(out_base, "susie")
  susie_files <- if (dir.exists(susie_dir)) {
    list.files(susie_dir, pattern = "\\.tsv$", full.names = TRUE)
  } else character(0)

  if (length(susie_files) > 0) {
    n_ok <- 0L
    for (f in susie_files) {
      tryCatch({
        dt <- fread(f)
        if (nrow(dt) == 0) return(invisible())
        dt[, study := study]
        dt[, ancestry := ancestry]
        dt[, trait_type := trait_type]
        dt[, method := "susie"]
        is_nonconv <- grepl("_notconverged\\.tsv$", f) ||
          ("converged" %in% colnames(dt) && any(dt$converged == FALSE, na.rm = TRUE))
        if (is_nonconv) {
          dt[, PIP := NA_real_]
          dt[, CS := NA_integer_]
        }
        all_susie[[length(all_susie) + 1]] <- dt
        n_ok <- n_ok + 1L
      }, error = function(e) {
        cat("  WARN SuSiE ", basename(f), ":", conditionMessage(e), "\n")
      })
    }
    cat("  SuSiE: read", n_ok, "of", length(susie_files), "locus files\n")
  } else {
    cat("  SuSiE: no results\n")
  }

  # ----- CARMA -----
  carma_dir <- file.path(out_base, "CARMA")
  carma_files <- if (dir.exists(carma_dir)) {
    list.files(carma_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
  } else character(0)

  if (length(carma_files) > 0) {
    n_ok <- 0L; n_empty <- 0L
    for (f in carma_files) {
      if (file.size(f) == 0) { n_empty <- n_empty + 1L; next }
      tryCatch({
        dt <- fread(f)
        if (nrow(dt) == 0) { n_empty <- n_empty + 1L; return(invisible()) }
        dt[, study := study]
        dt[, ancestry := ancestry]
        dt[, trait_type := trait_type]
        dt[, method := "carma"]
        all_carma[[length(all_carma) + 1]] <- dt
        n_ok <- n_ok + 1L
      }, error = function(e) {
        n_empty <<- n_empty + 1L
        cat("  WARN CARMA ", basename(f), ":", conditionMessage(e), "\n")
      })
    }
    cat(sprintf("  CARMA: read %d loci (%d empty/unreadable; run 04d if many)\n",
                n_ok, n_empty))
  } else {
    cat("  CARMA: no results\n")
  }
}

# ----- Combine -----
cat("\n=== Combining results ===\n")

susie_combined <- if (length(all_susie) > 0) rbindlist(all_susie, fill = TRUE) else data.table()
carma_combined <- if (length(all_carma) > 0) rbindlist(all_carma, fill = TRUE) else data.table()
cat("SuSiE rows:", nrow(susie_combined), "  CARMA rows:", nrow(carma_combined), "\n")

if (nrow(susie_combined) == 0 && nrow(carma_combined) == 0) {
  stop("No SuSiE or CARMA results to aggregate.")
}

if (nrow(susie_combined) > 0 && nrow(carma_combined) > 0) {
  susie_slim <- susie_combined %>%
    select(chromosome, position, allele1, allele2, locus, study, ancestry,
           susie_pip = PIP, susie_cs = CS, susie_converged = converged)
  carma_slim <- carma_combined %>%
    select(chromosome, position, allele1, allele2, locus, study, ancestry,
           carma_pip = PIP, carma_cs = CS, carma_outlier = outlier)

  combined <- full_join(susie_slim, carma_slim,
                        by = c("chromosome", "position", "allele1", "allele2",
                               "locus", "study", "ancestry"))

  combined <- combined %>%
    mutate(
      susie_reliable  = (!is.na(susie_converged) & susie_converged == TRUE),
      susie_pip_clean = ifelse(susie_reliable, susie_pip, NA_real_),
      concordant_pip  = ifelse(!is.na(susie_pip_clean) & !is.na(carma_pip),
                               pmin(susie_pip_clean, carma_pip), NA_real_),
      max_pip         = ifelse(ancestry != "EUR" & !is.na(carma_pip),
                               carma_pip,
                               pmax(susie_pip, carma_pip, na.rm = TRUE)),
      both_in_cs      = (!is.na(susie_cs) & !is.na(carma_cs) &
                         susie_cs > 0 & carma_cs > 0),
      either_in_cs    = ((!is.na(susie_cs) & susie_cs > 0) |
                         (!is.na(carma_cs) & carma_cs > 0)),
      variant_id      = paste(chromosome, position, allele1, allele2, sep = ":"),
      recommended_pip = case_when(
        ancestry != "EUR" ~ carma_pip,
        !susie_reliable   ~ carma_pip,
        TRUE              ~ pmax(susie_pip, carma_pip, na.rm = TRUE)
      )
    )
} else if (nrow(susie_combined) > 0) {
  combined <- susie_combined %>%
    select(chromosome, position, allele1, allele2, locus, study, ancestry,
           susie_pip = PIP, susie_cs = CS, susie_converged = converged) %>%
    mutate(
      carma_pip = NA_real_, carma_cs = NA_integer_, carma_outlier = NA_integer_,
      susie_reliable  = (!is.na(susie_converged) & susie_converged == TRUE),
      susie_pip_clean = ifelse(susie_reliable, susie_pip, NA_real_),
      concordant_pip  = NA_real_,
      max_pip         = susie_pip,
      both_in_cs      = FALSE,
      either_in_cs    = (!is.na(susie_cs) & susie_cs > 0),
      variant_id      = paste(chromosome, position, allele1, allele2, sep = ":"),
      recommended_pip = ifelse(susie_reliable, susie_pip, NA_real_)
    )
} else {
  stop("CARMA-only mode not implemented; need SuSiE outputs.")
}

combined <- as.data.table(combined)
setorder(combined, study, chromosome, position)
dir.create(dirname(OUT_PATH), showWarnings = FALSE, recursive = TRUE)
fwrite(combined, OUT_PATH)
cat("\nWrote:", OUT_PATH, "(", nrow(combined), "rows )\n")

# Summary
cat("\n=== Summary ===\n")
cat("Total variant rows:", nrow(combined), "\n")
cat("Loci:", uniqueN(combined[, .(study, locus)]), "\n")
cat("Studies:", uniqueN(combined$study), "\n")
cat("Variants with susie_pip > 0.5 (converged):",
    sum(combined$susie_pip > 0.5 & combined$susie_reliable, na.rm = TRUE), "\n")
cat("Variants with carma_pip > 0.5:",
    sum(combined$carma_pip > 0.5, na.rm = TRUE), "\n")
cat("Variants with recommended_pip > 0.5:",
    sum(combined$recommended_pip > 0.5, na.rm = TRUE), "\n")
cat("Variants in either CS:",
    sum(combined$either_in_cs, na.rm = TRUE), "\n")
cat("Variants in both CS:",
    sum(combined$both_in_cs, na.rm = TRUE), "\n")

per_study <- combined[, .(
  n_variants = .N,
  n_loci     = uniqueN(locus),
  n_susie_cs = sum(susie_cs > 0, na.rm = TRUE),
  n_carma_cs = sum(carma_cs > 0, na.rm = TRUE),
  n_high_pip = sum(recommended_pip > 0.5, na.rm = TRUE)
), by = study]
cat("\nPer-study summary:\n")
print(per_study)
