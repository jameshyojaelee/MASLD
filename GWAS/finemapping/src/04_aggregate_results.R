#!/usr/bin/env Rscript
# 04_aggregate_results.R
# Aggregates per-locus SuSiE + CARMA results into master tables
# Usage: Rscript 04_aggregate_results.R

library(data.table)
library(dplyr)
library(tidyr)
library(readr)

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR <- file.path(BASE_DIR, "GWAS/finemapping")
setwd(FM_DIR)

REGISTRY <- fread("config/gwas_registry.tsv")
RESULTS_DIR <- file.path(FM_DIR, "results")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("Aggregating fine-mapping results\n")
cat("============================================================\n")

all_susie <- list()
all_carma <- list()
skipped_loci <- list()

for (i in 1:nrow(REGISTRY)) {
  study <- REGISTRY$study_name[i]
  ld_pop <- REGISTRY$ancestry[i]
  window_mb <- REGISTRY$window_mb[i]
  out_base <- file.path(FM_DIR, "output", study, paste0(ld_pop, "_", window_mb, "Mb"))

  cat("\n--- Study:", study, "---\n")

  # Collect SuSiE results
  susie_files <- list.files(file.path(out_base, "susie"),
                            pattern = "\\.tsv$", full.names = TRUE)
  if (length(susie_files) > 0) {
    for (f in susie_files) {
      tryCatch({
        dt <- fread(f)
        dt$study <- study
        dt$ancestry <- REGISTRY$ancestry[i]
        dt$trait_type <- REGISTRY$trait_type[i]
        dt$method <- "susie"

        # Non-converged SuSiE results can report PIP=1.0 for multiple variants
        # simultaneously (physically impossible — PIPs within a credible set
        # should sum to ~1 for a single causal signal). Null these out to
        # prevent spurious high-confidence calls from entering downstream tables.
        is_nonconv <- grepl("_notconverged\\.tsv$", f) ||
          ("converged" %in% colnames(dt) && any(dt$converged == FALSE, na.rm = TRUE))
        if (is_nonconv) {
          cat("  NOTE: Nullifying PIP/CS for non-converged locus:", basename(f), "\n")
          dt[, PIP := NA_real_]
          dt[, CS := NA_integer_]
        }

        all_susie[[length(all_susie) + 1]] <- dt
      }, error = function(e) {
        cat("  WARNING: Could not read", basename(f), ":", conditionMessage(e), "\n")
      })
    }
    cat("  SuSiE: read", length(susie_files), "locus files\n")
  } else {
    cat("  SuSiE: no results found\n")
  }

  # Collect CARMA results
  carma_files <- list.files(file.path(out_base, "CARMA"),
                            pattern = "\\.txt\\.gz$", full.names = TRUE)
  if (length(carma_files) > 0) {
    for (f in carma_files) {
      tryCatch({
        dt <- fread(f)
        dt$study <- study
        dt$ancestry <- REGISTRY$ancestry[i]
        dt$trait_type <- REGISTRY$trait_type[i]
        dt$method <- "carma"
        all_carma[[length(all_carma) + 1]] <- dt
      }, error = function(e) {
        cat("  WARNING: Could not read", basename(f), ":", conditionMessage(e), "\n")
      })
    }
    cat("  CARMA: read", length(carma_files), "locus files\n")
  } else {
    cat("  CARMA: no results found\n")
  }

  # Count skipped loci
  skip_files <- list.files(out_base, pattern = "SKIPPED_", full.names = TRUE)
  if (length(skip_files) > 0) {
    cat("  Skipped:", length(skip_files), "loci\n")
    skipped_loci[[length(skipped_loci) + 1]] <- data.frame(
      study = study, n_skipped = length(skip_files)
    )
  }
}

# Combine all results
cat("\n=== Combining results ===\n")

if (length(all_susie) > 0) {
  susie_combined <- rbindlist(all_susie, fill = TRUE)
  cat("SuSiE: ", nrow(susie_combined), "total variant-locus entries across",
      length(unique(susie_combined$study)), "studies\n")
  fwrite(susie_combined, file.path(RESULTS_DIR, "susie_all_results.csv"))
} else {
  cat("WARNING: No SuSiE results to combine\n")
  susie_combined <- data.table()
}

if (length(all_carma) > 0) {
  carma_combined <- rbindlist(all_carma, fill = TRUE)
  cat("CARMA:", nrow(carma_combined), "total variant-locus entries across",
      length(unique(carma_combined$study)), "studies\n")
  fwrite(carma_combined, file.path(RESULTS_DIR, "carma_all_results.csv"))
} else {
  cat("WARNING: No CARMA results to combine\n")
  carma_combined <- data.table()
}

# Cross-method concordance
if (nrow(susie_combined) > 0 && nrow(carma_combined) > 0) {
  cat("\n=== Cross-method concordance ===\n")

  # Join on variant position + study + locus
  susie_slim <- susie_combined %>%
    select(chromosome, position, allele1, allele2, locus, study, ancestry,
           susie_pip = PIP, susie_cs = CS, susie_converged = converged)
  carma_slim <- carma_combined %>%
    select(chromosome, position, allele1, allele2, locus, study, ancestry,
           carma_pip = PIP, carma_cs = CS, carma_outlier = outlier)

  combined <- full_join(susie_slim, carma_slim,
                        by = c("chromosome", "position", "allele1", "allele2",
                               "locus", "study", "ancestry"))

  # For non-converged SuSiE loci, set susie_pip to NA in concordance metrics
  # (PIP values from non-converged runs are unreliable)
  combined <- combined %>%
    mutate(
      susie_reliable = (!is.na(susie_converged) & susie_converged == TRUE),
      susie_pip_clean = ifelse(susie_reliable, susie_pip, NA_real_)
    )

  # Concordance metrics
  combined <- combined %>%
    mutate(
      # concordant_pip: only defined when BOTH methods ran and SuSiE converged
      # (avoid inflating single-method results via pmin(..., na.rm=TRUE)
      # returning the lone value)
      concordant_pip = ifelse(!is.na(susie_pip_clean) & !is.na(carma_pip),
                              pmin(susie_pip_clean, carma_pip), NA_real_),
      # For non-EUR studies, prefer CARMA PIP (1KG LD reference may cause SuSiE ghost CS)
      max_pip = ifelse(ancestry != "EUR" & !is.na(carma_pip),
                       carma_pip,
                       pmax(susie_pip, carma_pip, na.rm = TRUE)),
      both_in_cs = (!is.na(susie_cs) & !is.na(carma_cs) & susie_cs > 0 & carma_cs > 0),
      either_in_cs = ((!is.na(susie_cs) & susie_cs > 0) | (!is.na(carma_cs) & carma_cs > 0)),
      variant_id = paste(chromosome, position, allele1, allele2, sep = ":")
    )

  cat("Combined:", nrow(combined), "variant entries\n")

  # Compute recommended_pip:
  #   EUR converged  -> max(susie_pip, carma_pip)
  #   EUR not-converged -> carma_pip only (PIP already NA from ingestion step)
  #   non-EUR        -> carma_pip only (1KG LD panels have ghost CS artifact
  #                     rate due to small sample sizes; CARMA is preferred)
  combined <- combined %>%
    mutate(
      recommended_pip = case_when(
        ancestry != "EUR" ~ carma_pip,
        !susie_reliable   ~ carma_pip,
        TRUE              ~ pmax(susie_pip, carma_pip, na.rm = TRUE)
      )
    )

  # ===========================================================================
  # Integrate cross-ancestry joint fine-mapping (SuSiEX + MESuSiE)
  # ===========================================================================
  susiex_file <- file.path(RESULTS_DIR, "susiex", "susiex_variant_summary.csv")
  mesusie_file <- file.path(RESULTS_DIR, "mesusie", "mesusie_variant_summary.csv")

  if (file.exists(susiex_file)) {
    cat("\n=== Integrating SuSiEX results ===\n")
    susiex <- fread(susiex_file, select = c("chr", "pos", "a1", "a2", "OVRL_PIP"))
    setnames(susiex, c("chr", "pos", "a1", "a2", "OVRL_PIP"),
             c("chromosome", "position", "allele1", "allele2", "susiex_pip"))
    susiex[, chromosome := as.integer(chromosome)]
    susiex[, position  := as.integer(position)]
    # Deduplicate: keep max PIP per variant (same variant may appear across traits)
    susiex <- susiex[, .(susiex_pip = max(susiex_pip, na.rm = TRUE)),
                     by = .(chromosome, position, allele1, allele2)]
    cat("  SuSiEX: ", nrow(susiex), " unique variants loaded\n")
  } else {
    cat("  SuSiEX summary not found at:", susiex_file, "\n")
    susiex <- NULL
  }

  if (file.exists(mesusie_file)) {
    cat("\n=== Integrating MESuSiE results ===\n")
    mesusie <- fread(mesusie_file,
                     select = c("chr", "pos", "a1", "a2", "pip", "in_cs", "cs_category"))
    setnames(mesusie, c("chr", "pos", "a1", "a2", "pip"),
             c("chromosome", "position", "allele1", "allele2", "mesusie_pip"))
    mesusie[, chromosome := as.integer(chromosome)]
    mesusie[, position  := as.integer(position)]
    mesusie[, mesusie_in_shared_cs := (in_cs == TRUE & grepl("EUR_EAS|shared", cs_category, ignore.case = TRUE))]
    # Deduplicate: keep max PIP and any shared CS hit per variant
    mesusie <- mesusie[, .(mesusie_pip = max(mesusie_pip, na.rm = TRUE),
                           mesusie_in_shared_cs = any(mesusie_in_shared_cs, na.rm = TRUE)),
                       by = .(chromosome, position, allele1, allele2)]
    cat("  MESuSiE:", nrow(mesusie), "unique variants loaded\n")
  } else {
    cat("  MESuSiE summary not found at:", mesusie_file, "\n")
    mesusie <- NULL
  }

  # Left-join joint PIPs onto combined (handle allele flips by trying both orientations)
  combined <- as.data.table(combined)

  if (!is.null(susiex)) {
    # Direct match
    combined <- merge(combined, susiex, by = c("chromosome", "position", "allele1", "allele2"),
                      all.x = TRUE)
    # Allele-flipped match for unmatched rows
    unmatched <- is.na(combined$susiex_pip)
    if (any(unmatched)) {
      susiex_flip <- copy(susiex)
      setnames(susiex_flip, c("allele1", "allele2"), c("allele2", "allele1"))
      setnames(susiex_flip, "susiex_pip", "susiex_pip_flip")
      combined <- merge(combined, susiex_flip,
                        by = c("chromosome", "position", "allele1", "allele2"),
                        all.x = TRUE)
      combined[is.na(susiex_pip) & !is.na(susiex_pip_flip), susiex_pip := susiex_pip_flip]
      combined[, susiex_pip_flip := NULL]
    }
    cat("  SuSiEX matched:", sum(!is.na(combined$susiex_pip)), "/", nrow(combined), "variants\n")
  } else {
    combined[, susiex_pip := NA_real_]
  }

  if (!is.null(mesusie)) {
    # Direct match
    combined <- merge(combined, mesusie, by = c("chromosome", "position", "allele1", "allele2"),
                      all.x = TRUE)
    # Allele-flipped match for unmatched rows
    unmatched <- is.na(combined$mesusie_pip)
    if (any(unmatched)) {
      mesusie_flip <- copy(mesusie)
      setnames(mesusie_flip, c("allele1", "allele2"), c("allele2", "allele1"))
      setnames(mesusie_flip, c("mesusie_pip", "mesusie_in_shared_cs"),
               c("mesusie_pip_flip", "mesusie_shared_flip"))
      combined <- merge(combined, mesusie_flip,
                        by = c("chromosome", "position", "allele1", "allele2"),
                        all.x = TRUE)
      combined[is.na(mesusie_pip) & !is.na(mesusie_pip_flip),
               `:=`(mesusie_pip = mesusie_pip_flip,
                    mesusie_in_shared_cs = mesusie_shared_flip)]
      combined[, c("mesusie_pip_flip", "mesusie_shared_flip") := NULL]
    }
    cat("  MESuSiE matched:", sum(!is.na(combined$mesusie_pip)), "/", nrow(combined), "variants\n")
  } else {
    combined[, `:=`(mesusie_pip = NA_real_, mesusie_in_shared_cs = NA)]
  }

  # Derived columns
  combined[, in_joint_cs := fifelse(
    (!is.na(susiex_pip) & susiex_pip > 0.5) | (mesusie_in_shared_cs == TRUE),
    TRUE, FALSE, na = FALSE)]

  # Update recommended_pip: max of existing and joint PIPs
  combined[, recommended_pip := pmax(recommended_pip, susiex_pip, mesusie_pip, na.rm = TRUE)]

  # Update either_in_cs to include joint CS membership
  combined[, either_in_cs := (either_in_cs | in_joint_cs)]

  n_joint <- sum(combined$in_joint_cs, na.rm = TRUE)
  cat("  Joint CS variants:", n_joint, "\n")
  cat("  Updated either_in_cs:", sum(combined$either_in_cs, na.rm = TRUE), "variants\n")

  fwrite(combined, file.path(RESULTS_DIR, "combined_finemapping.csv"))

  # High-PIP variants (use recommended_pip which accounts for EAS artifacts/non-convergence)
  high_pip <- combined %>%
    filter(recommended_pip > 0.5) %>%
    arrange(desc(recommended_pip))
  fwrite(high_pip, file.path(RESULTS_DIR, "high_pip_variants.csv"))
  cat("High PIP (recommended >0.5):", nrow(high_pip), "variants\n")

  # Credible sets summary
  cs_variants <- combined %>%
    filter(either_in_cs) %>%
    arrange(study, locus, desc(max_pip))
  fwrite(cs_variants, file.path(RESULTS_DIR, "credible_sets.csv"))
  cat("In credible sets:", nrow(cs_variants), "variants\n")

  # Per-study summary
  study_summary <- combined %>%
    group_by(study) %>%
    summarise(
      n_variants = n(),
      n_loci = n_distinct(locus),
      n_susie_cs = sum(susie_cs > 0, na.rm = TRUE),
      n_carma_cs = sum(carma_cs > 0, na.rm = TRUE),
      n_both_cs = sum(both_in_cs, na.rm = TRUE),
      n_high_pip = sum(max_pip > 0.5, na.rm = TRUE),
      mean_concordant_pip = mean(concordant_pip, na.rm = TRUE),
      .groups = "drop"
    )
  fwrite(study_summary, file.path(RESULTS_DIR, "study_summary.csv"))
  cat("\nPer-study summary:\n")
  print(as.data.frame(study_summary))
}

cat("\n============================================================\n")
cat("Results written to:", RESULTS_DIR, "\n")
cat("============================================================\n")
