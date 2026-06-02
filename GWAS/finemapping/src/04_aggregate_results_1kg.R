#!/usr/bin/env Rscript
# 04_aggregate_results_1kg.R
# 1KG EUR LD-panel analog of 04_aggregate_results.R. Aggregates per-locus
# SuSiE + CARMA outputs from the 1KG EUR re-run into master tables.
#
# Differences vs 04_aggregate_results.R:
#   - For EUR-ancestry GWAS, reads from output/<study>/1kg_eur_0.5Mb/
#     instead of output/<study>/EUR_0.5Mb/. (The 1KG re-run only re-fit
#     the 17 EUR GWAS; non-EUR studies are unchanged.)
#   - For non-EUR GWAS (EAS / AFR / SAS), keeps the v1 path pattern
#     <ancestry>_<window>Mb/ — those outputs are NOT re-run and remain
#     identical to the v1 aggregator's view.
#   - CARMA results are read directly from the per-locus .rds files
#     because the 1kg_eur run wrote 0-byte .txt.gz files (zlib compression
#     bug in fwrite on this system; see 04d_reexport_carma_txt.R for
#     background). This script reconstructs the same tabular CARMA output
#     in-memory by joining .rds[[1]]$PIPs / Credible set / Outliers with
#     the LD-matched SuSiE .tsv (same variant order). For non-EUR studies
#     where the v1 .txt.gz files are valid, those are read as before.
#   - Outputs are suffixed with _1kg and write side-by-side with v1:
#       results/susie_all_results_1kg.csv
#       results/carma_all_results_1kg.csv
#       results/combined_finemapping_1kg.csv
#       results/high_pip_variants_1kg.csv
#       results/credible_sets_1kg.csv
#       results/study_summary_1kg.csv
#   - SuSiEX integration: prefers results/susiex_1kg/susiex_variant_summary_1kg.csv
#     (sister 1KG run); falls back to v1 results/susiex/susiex_variant_summary.csv
#     if 1KG version not yet produced. MESuSiE integration uses v1 path
#     unless a _1kg variant exists.
#
# Usage: micromamba run -n rnaseq Rscript GWAS/finemapping/src/04_aggregate_results_1kg.R

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

OUT_SUFFIX <- "_1kg"

# ---------------------------------------------------------------------------
# Helper: reconstruct CARMA output table from .rds + matching SuSiE .tsv
# Mirrors the primary path in 04d_reexport_carma_txt.R (PIPs, Credible set,
# Outliers fields). Returns a data.table or NULL on failure.
# ---------------------------------------------------------------------------
load_carma_from_rds <- function(rds_path, susie_dir) {
  base <- sub("\\.rds$", "", basename(rds_path))

  # Parse {study}_{LDpanel}_{window}Mb_{locus}.
  # NOTE: Some study_names (e.g. 2023_36280732_NAFLD_UKBB_EUR) themselves
  # contain "_UKBB_" so strsplit on the first occurrence is wrong.
  # We anchor the LD-panel + window from the right end of the filename.
  m <- regmatches(base,
                  regexec("_(UKBB|1KG_EUR|1KG_EAS|1KG_AFR|1KG_SAS)_([0-9.]+)Mb_([0-9]+\\.[0-9]+)$",
                          base))[[1]]
  if (length(m) < 4) return(NULL)
  ldpanel <- m[2]
  locus   <- m[4]

  # Locate the SuSiE .tsv (variant metadata in LD-matched order)
  susie_pattern <- paste0("^", ldpanel, "_", locus, "_cov0\\.95_")
  susie_candidates <- list.files(susie_dir, pattern = susie_pattern,
                                 full.names = TRUE)
  susie_candidates <- susie_candidates[grepl("\\.tsv$", susie_candidates)]
  if (length(susie_candidates) == 0) return(NULL)
  susie_path <- susie_candidates[1]

  rds      <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  susie_dt <- tryCatch(fread(susie_path), error = function(e) NULL)
  if (is.null(rds) || is.null(susie_dt)) return(NULL)

  n_pips  <- length(rds[[1]]$PIPs)
  n_susie <- nrow(susie_dt)
  if (n_pips != n_susie) {
    # Dimension mismatch — fall back to NULL; v1 04d has a heavier fallback
    # using the SuSiE RDS pip names + ss/ files but it is not needed for the
    # 1KG EUR run because all loci have matching dimensions (verified at
    # runtime; mismatches are reported and skipped).
    return(NULL)
  }

  result <- data.table(
    SNP        = susie_dt$SNP,
    chromosome = susie_dt$chromosome,
    position   = susie_dt$position,
    allele1    = susie_dt$allele1,
    allele2    = susie_dt$allele2,
    beta       = susie_dt$beta,
    se         = susie_dt$se,
    pval       = susie_dt$pval,
    locus      = susie_dt$locus,
    Z          = susie_dt$beta / susie_dt$se,
    PIP        = rds[[1]]$PIPs,
    CS         = 0L,
    outlier    = 0L
  )

  # Assign credible sets
  cs_obj <- rds[[1]]$`Credible set`
  if (length(cs_obj) >= 2 && length(cs_obj[[2]]) > 0) {
    for (l in seq_along(cs_obj[[2]])) {
      idx <- cs_obj[[2]][[l]]
      idx <- idx[idx >= 1 & idx <= nrow(result)]
      result$CS[idx] <- l
    }
  }

  # Assign outliers
  outlier_idx <- rds[[1]]$Outliers$Index
  if (length(outlier_idx) > 0) {
    outlier_idx <- outlier_idx[outlier_idx >= 1 & outlier_idx <= nrow(result)]
    result$outlier[outlier_idx] <- 1L
  }

  return(result)
}

cat("============================================================\n")
cat("Aggregating fine-mapping results (1KG EUR LD panel)\n")
cat("============================================================\n")

all_susie <- list()
all_carma <- list()
skipped_loci <- list()
n_carma_from_rds <- 0L
n_carma_from_txt <- 0L

for (i in 1:nrow(REGISTRY)) {
  study <- REGISTRY$study_name[i]
  ld_pop <- REGISTRY$ancestry[i]
  window_mb <- REGISTRY$window_mb[i]

  # KEY DIFFERENCE: EUR studies use the new 1kg_eur_<win>Mb path; non-EUR
  # studies keep the legacy <ancestry>_<win>Mb path (unchanged from v1).
  if (ld_pop == "EUR") {
    out_base <- file.path(FM_DIR, "output", study,
                          paste0("1kg_eur_", window_mb, "Mb"))
  } else {
    out_base <- file.path(FM_DIR, "output", study,
                          paste0(ld_pop, "_", window_mb, "Mb"))
  }

  cat("\n--- Study:", study, "(ancestry=", ld_pop, ", path=",
      basename(out_base), ") ---\n")

  if (!dir.exists(out_base)) {
    cat("  WARNING: output dir does not exist:", out_base, "\n")
    next
  }

  # ---- SuSiE ----
  susie_dir <- file.path(out_base, "susie")
  susie_files <- list.files(susie_dir, pattern = "\\.tsv$", full.names = TRUE)
  if (length(susie_files) > 0) {
    for (f in susie_files) {
      tryCatch({
        dt <- fread(f)
        dt$study     <- study
        dt$ancestry  <- REGISTRY$ancestry[i]
        dt$trait_type <- REGISTRY$trait_type[i]
        dt$method    <- "susie"

        is_nonconv <- grepl("_notconverged\\.tsv$", f) ||
          ("converged" %in% colnames(dt) && any(dt$converged == FALSE, na.rm = TRUE))
        if (is_nonconv) {
          cat("  NOTE: Nullifying PIP/CS for non-converged locus:", basename(f), "\n")
          dt[, PIP := NA_real_]
          dt[, CS  := NA_integer_]
        }

        all_susie[[length(all_susie) + 1]] <- dt
      }, error = function(e) {
        cat("  WARNING: Could not read", basename(f), ":",
            conditionMessage(e), "\n")
      })
    }
    cat("  SuSiE: read", length(susie_files), "locus files\n")
  } else {
    cat("  SuSiE: no results found\n")
  }

  # ---- CARMA ----
  # Prefer .txt.gz when valid (non-EUR retains v1 outputs). For EUR (1kg path)
  # the .txt.gz files are 0 bytes from a zlib bug, so we read .rds directly
  # and reconstruct the table in-memory using the matching SuSiE tsv.
  carma_dir <- file.path(out_base, "CARMA")
  if (dir.exists(carma_dir)) {
    txt_files <- list.files(carma_dir, pattern = "\\.txt\\.gz$",
                            full.names = TRUE)
    rds_files <- list.files(carma_dir, pattern = "\\.rds$",
                            full.names = TRUE)

    txt_used <- character(0)
    if (length(txt_files) > 0) {
      for (f in txt_files) {
        if (file.size(f) == 0) next  # skip 0-byte (broken) files
        tryCatch({
          dt <- fread(f)
          if (ncol(dt) == 0 || nrow(dt) == 0) {
            cat("  WARNING: Empty CARMA table:", basename(f), "\n")
            return(invisible())
          }
          dt$study      <- study
          dt$ancestry   <- REGISTRY$ancestry[i]
          dt$trait_type <- REGISTRY$trait_type[i]
          dt$method     <- "carma"
          all_carma[[length(all_carma) + 1]] <- dt
          n_carma_from_txt <<- n_carma_from_txt + 1L
          txt_used <<- c(txt_used, sub("\\.txt\\.gz$", "", basename(f)))
        }, error = function(e) {
          cat("  WARNING: Could not read", basename(f), ":",
              conditionMessage(e), "\n")
        })
      }
    }

    # For .rds without a usable .txt.gz, reconstruct in-memory.
    if (length(rds_files) > 0) {
      for (f in rds_files) {
        base_no_ext <- sub("\\.rds$", "", basename(f))
        if (base_no_ext %in% txt_used) next  # already covered by .txt.gz
        dt <- load_carma_from_rds(f, susie_dir)
        if (is.null(dt)) {
          cat("  WARNING: Could not reconstruct from RDS:", basename(f), "\n")
          next
        }
        dt$study      <- study
        dt$ancestry   <- REGISTRY$ancestry[i]
        dt$trait_type <- REGISTRY$trait_type[i]
        dt$method     <- "carma"
        all_carma[[length(all_carma) + 1]] <- dt
        n_carma_from_rds <- n_carma_from_rds + 1L
      }
    }

    cat("  CARMA: read", length(txt_files), "txt.gz +",
        length(rds_files), "rds files (txt_used=",
        length(txt_used), ", rds_reconstructed=",
        length(rds_files) - length(txt_used), ")\n")
  } else {
    cat("  CARMA: directory not found\n")
  }

  # ---- Skipped loci ----
  skip_files <- list.files(out_base, pattern = "SKIPPED_", full.names = TRUE)
  if (length(skip_files) > 0) {
    cat("  Skipped:", length(skip_files), "loci\n")
    skipped_loci[[length(skipped_loci) + 1]] <- data.frame(
      study = study, n_skipped = length(skip_files)
    )
  }
}

cat("\nCARMA ingestion summary:\n")
cat("  From .txt.gz:", n_carma_from_txt, "loci\n")
cat("  From .rds (reconstructed):", n_carma_from_rds, "loci\n")

# ---------------------------------------------------------------------------
# Combine
# ---------------------------------------------------------------------------
cat("\n=== Combining results ===\n")

if (length(all_susie) > 0) {
  susie_combined <- rbindlist(all_susie, fill = TRUE)
  cat("SuSiE: ", nrow(susie_combined), "total variant-locus entries across",
      length(unique(susie_combined$study)), "studies\n")
  fwrite(susie_combined,
         file.path(RESULTS_DIR, paste0("susie_all_results", OUT_SUFFIX, ".csv")))
} else {
  cat("WARNING: No SuSiE results to combine\n")
  susie_combined <- data.table()
}

if (length(all_carma) > 0) {
  carma_combined <- rbindlist(all_carma, fill = TRUE)
  cat("CARMA:", nrow(carma_combined), "total variant-locus entries across",
      length(unique(carma_combined$study)), "studies\n")
  fwrite(carma_combined,
         file.path(RESULTS_DIR, paste0("carma_all_results", OUT_SUFFIX, ".csv")))
} else {
  cat("WARNING: No CARMA results to combine\n")
  carma_combined <- data.table()
}

# ---------------------------------------------------------------------------
# Cross-method concordance
# ---------------------------------------------------------------------------
if (nrow(susie_combined) > 0 && nrow(carma_combined) > 0) {
  cat("\n=== Cross-method concordance ===\n")

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
      susie_reliable = (!is.na(susie_converged) & susie_converged == TRUE),
      susie_pip_clean = ifelse(susie_reliable, susie_pip, NA_real_)
    )

  combined <- combined %>%
    mutate(
      concordant_pip = ifelse(!is.na(susie_pip_clean) & !is.na(carma_pip),
                              pmin(susie_pip_clean, carma_pip), NA_real_),
      max_pip = ifelse(ancestry != "EUR" & !is.na(carma_pip),
                       carma_pip,
                       pmax(susie_pip, carma_pip, na.rm = TRUE)),
      both_in_cs = (!is.na(susie_cs) & !is.na(carma_cs) & susie_cs > 0 & carma_cs > 0),
      either_in_cs = ((!is.na(susie_cs) & susie_cs > 0) | (!is.na(carma_cs) & carma_cs > 0)),
      variant_id = paste(chromosome, position, allele1, allele2, sep = ":")
    )

  cat("Combined:", nrow(combined), "variant entries\n")

  # recommended_pip — same EAS-aware policy as v1
  combined <- combined %>%
    mutate(
      recommended_pip = case_when(
        ancestry != "EUR" ~ carma_pip,
        !susie_reliable   ~ carma_pip,
        TRUE              ~ pmax(susie_pip, carma_pip, na.rm = TRUE)
      )
    )

  # ---------------------------------------------------------------------------
  # Cross-ancestry joint fine-mapping: prefer 1kg variant if available
  # ---------------------------------------------------------------------------
  susiex_file_1kg <- file.path(RESULTS_DIR, "susiex_1kg",
                               "susiex_variant_summary_1kg.csv")
  susiex_file_v1  <- file.path(RESULTS_DIR, "susiex",
                               "susiex_variant_summary.csv")
  if (file.exists(susiex_file_1kg)) {
    susiex_file <- susiex_file_1kg
    cat("\n=== Integrating SuSiEX (1KG) results ===\n")
  } else if (file.exists(susiex_file_v1)) {
    susiex_file <- susiex_file_v1
    cat("\n=== Integrating SuSiEX (v1 fallback — 1KG variant not yet produced) ===\n")
  } else {
    susiex_file <- NULL
  }

  mesusie_file_1kg <- file.path(RESULTS_DIR, "mesusie",
                                "mesusie_variant_summary_1kg.csv")
  mesusie_file_v1  <- file.path(RESULTS_DIR, "mesusie",
                                "mesusie_variant_summary.csv")
  if (file.exists(mesusie_file_1kg)) {
    mesusie_file <- mesusie_file_1kg
  } else if (file.exists(mesusie_file_v1)) {
    mesusie_file <- mesusie_file_v1
  } else {
    mesusie_file <- NULL
  }

  if (!is.null(susiex_file)) {
    susiex <- fread(susiex_file, select = c("chr", "pos", "a1", "a2", "OVRL_PIP"))
    setnames(susiex, c("chr", "pos", "a1", "a2", "OVRL_PIP"),
             c("chromosome", "position", "allele1", "allele2", "susiex_pip"))
    susiex[, chromosome := as.integer(chromosome)]
    susiex[, position  := as.integer(position)]
    susiex <- susiex[, .(susiex_pip = max(susiex_pip, na.rm = TRUE)),
                     by = .(chromosome, position, allele1, allele2)]
    cat("  SuSiEX: ", nrow(susiex), " unique variants loaded from",
        basename(susiex_file), "\n")
  } else {
    cat("  SuSiEX summary not found (1KG or v1)\n")
    susiex <- NULL
  }

  if (!is.null(mesusie_file)) {
    cat("\n=== Integrating MESuSiE results from", basename(mesusie_file), "===\n")
    mesusie <- fread(mesusie_file,
                     select = c("chr", "pos", "a1", "a2", "pip", "in_cs", "cs_category"))
    setnames(mesusie, c("chr", "pos", "a1", "a2", "pip"),
             c("chromosome", "position", "allele1", "allele2", "mesusie_pip"))
    mesusie[, chromosome := as.integer(chromosome)]
    mesusie[, position  := as.integer(position)]
    mesusie[, mesusie_in_shared_cs := (in_cs == TRUE &
              grepl("EUR_EAS|shared", cs_category, ignore.case = TRUE))]
    mesusie <- mesusie[, .(mesusie_pip = max(mesusie_pip, na.rm = TRUE),
                           mesusie_in_shared_cs = any(mesusie_in_shared_cs, na.rm = TRUE)),
                       by = .(chromosome, position, allele1, allele2)]
    cat("  MESuSiE:", nrow(mesusie), "unique variants loaded\n")
  } else {
    cat("  MESuSiE summary not found (1KG or v1)\n")
    mesusie <- NULL
  }

  combined <- as.data.table(combined)

  if (!is.null(susiex)) {
    combined <- merge(combined, susiex, by = c("chromosome", "position", "allele1", "allele2"),
                      all.x = TRUE)
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
    cat("  SuSiEX matched:", sum(!is.na(combined$susiex_pip)), "/",
        nrow(combined), "variants\n")
  } else {
    combined[, susiex_pip := NA_real_]
  }

  if (!is.null(mesusie)) {
    combined <- merge(combined, mesusie, by = c("chromosome", "position", "allele1", "allele2"),
                      all.x = TRUE)
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
    cat("  MESuSiE matched:", sum(!is.na(combined$mesusie_pip)), "/",
        nrow(combined), "variants\n")
  } else {
    combined[, `:=`(mesusie_pip = NA_real_, mesusie_in_shared_cs = NA)]
  }

  combined[, in_joint_cs := fifelse(
    (!is.na(susiex_pip) & susiex_pip > 0.5) | (mesusie_in_shared_cs == TRUE),
    TRUE, FALSE, na = FALSE)]
  combined[, recommended_pip := pmax(recommended_pip, susiex_pip, mesusie_pip, na.rm = TRUE)]
  combined[, either_in_cs := (either_in_cs | in_joint_cs)]

  n_joint <- sum(combined$in_joint_cs, na.rm = TRUE)
  cat("  Joint CS variants:", n_joint, "\n")
  cat("  Updated either_in_cs:", sum(combined$either_in_cs, na.rm = TRUE), "variants\n")

  out_combined <- file.path(RESULTS_DIR, paste0("combined_finemapping", OUT_SUFFIX, ".csv"))
  fwrite(combined, out_combined)

  high_pip <- combined %>%
    filter(recommended_pip > 0.5) %>%
    arrange(desc(recommended_pip))
  fwrite(high_pip, file.path(RESULTS_DIR, paste0("high_pip_variants", OUT_SUFFIX, ".csv")))
  cat("High PIP (recommended >0.5):", nrow(high_pip), "variants\n")

  cs_variants <- combined %>%
    filter(either_in_cs) %>%
    arrange(study, locus, desc(max_pip))
  fwrite(cs_variants, file.path(RESULTS_DIR, paste0("credible_sets", OUT_SUFFIX, ".csv")))
  cat("In credible sets:", nrow(cs_variants), "variants\n")

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
  fwrite(study_summary, file.path(RESULTS_DIR, paste0("study_summary", OUT_SUFFIX, ".csv")))
  cat("\nPer-study summary:\n")
  print(as.data.frame(study_summary))

  # ---------------------------------------------------------------------------
  # Headline stats + comparison vs v1 combined_finemapping.csv
  # ---------------------------------------------------------------------------
  cat("\n============================================================\n")
  cat("HEADLINE STATS (1KG EUR)\n")
  cat("============================================================\n")
  cat("Total variant-locus entries:        ", nrow(combined), "\n")
  cat("Unique variants:                    ", length(unique(combined$variant_id)), "\n")
  cat("Unique loci (study x locus):        ",
      nrow(unique(as.data.frame(combined)[, c("study", "locus")])), "\n")
  cat("Unique loci (locus only):           ", length(unique(combined$locus)), "\n")
  cat("Studies:                            ", length(unique(combined$study)), "\n")
  cat("Variants with recommended_pip > 0.5:", sum(combined$recommended_pip > 0.5,
                                                  na.rm = TRUE), "\n")
  cat("Variants with recommended_pip > 0.9:", sum(combined$recommended_pip > 0.9,
                                                  na.rm = TRUE), "\n")
  cat("Variants in any CS (susie or carma):", sum(combined$either_in_cs, na.rm = TRUE), "\n")
  cat("n loci with at least one SuSiE CS variant:",
      length(unique(paste(combined$study[combined$susie_cs > 0 & !is.na(combined$susie_cs)],
                          combined$locus[combined$susie_cs > 0 & !is.na(combined$susie_cs)],
                          sep = "::"))), "\n")

  v1_path <- file.path(RESULTS_DIR, "combined_finemapping.csv")
  if (file.exists(v1_path)) {
    cat("\n--- Comparison vs v1 (combined_finemapping.csv) ---\n")
    v1 <- fread(v1_path, select = c("study", "locus", "variant_id",
                                    "recommended_pip", "either_in_cs"))

    cat(sprintf("                          v1            1kg            delta\n"))
    cat(sprintf("Total entries:            %-13d %-13d %+d\n",
                nrow(v1), nrow(combined), nrow(combined) - nrow(v1)))
    cat(sprintf("Unique variants:          %-13d %-13d %+d\n",
                length(unique(v1$variant_id)),
                length(unique(combined$variant_id)),
                length(unique(combined$variant_id)) - length(unique(v1$variant_id))))
    cat(sprintf("recommended_pip > 0.5:    %-13d %-13d %+d\n",
                sum(v1$recommended_pip > 0.5, na.rm = TRUE),
                sum(combined$recommended_pip > 0.5, na.rm = TRUE),
                sum(combined$recommended_pip > 0.5, na.rm = TRUE) -
                  sum(v1$recommended_pip > 0.5, na.rm = TRUE)))
    cat(sprintf("recommended_pip > 0.9:    %-13d %-13d %+d\n",
                sum(v1$recommended_pip > 0.9, na.rm = TRUE),
                sum(combined$recommended_pip > 0.9, na.rm = TRUE),
                sum(combined$recommended_pip > 0.9, na.rm = TRUE) -
                  sum(v1$recommended_pip > 0.9, na.rm = TRUE)))
    cat(sprintf("Variants in any CS:       %-13d %-13d %+d\n",
                sum(v1$either_in_cs, na.rm = TRUE),
                sum(combined$either_in_cs, na.rm = TRUE),
                sum(combined$either_in_cs, na.rm = TRUE) -
                  sum(v1$either_in_cs, na.rm = TRUE)))
  } else {
    cat("\n(v1 combined_finemapping.csv not found; skipping comparison)\n")
  }
}

cat("\n============================================================\n")
cat("Results written to:", RESULTS_DIR, "\n")
cat("Primary output:    ",
    file.path(RESULTS_DIR, paste0("combined_finemapping", OUT_SUFFIX, ".csv")), "\n")
cat("============================================================\n")
