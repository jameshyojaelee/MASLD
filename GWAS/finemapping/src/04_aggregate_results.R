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
  anc <- REGISTRY$ancestry[i]
  # FM output dir = the ld_panel subdir for NON-EUR (the fresh 1kg_<anc> runs from
  # run_fm_nonEUR_expanded.sh -> 03); EUR canonical FM lives in <EUR>_0.5Mb. Using
  # ancestry for all (prior behaviour) pointed non-EUR at <ANC>_0.5Mb, which is empty
  # (MVP) or stale (BBJ/PanUKBB) -> the entire fresh non-EUR/MVP run was silently missed.
  ld_pop <- if (anc == "EUR") anc else REGISTRY$ld_panel[i]
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
           susie_pip = PIP, susie_cs = CS, susie_converged = converged, lambda_s)
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
  # Union the legacy enzyme 2-way run (UKBB<->BBJ ALT/AST/GGT) with the within-MVP
  # N-way cross-ancestry run (EUR/AFR/AMR/EAS; 7 traits, 411 loci; 2026-07-05). The
  # _mvp file is schema-identical; the dedup below keeps max PIP per (variant, trait),
  # so ALT/AST take the best of either run and GGT (legacy-only) / the MVP-only
  # disease+biomarker traits are preserved.
  susiex_files <- c(file.path(RESULTS_DIR, "susiex", "susiex_variant_summary.csv"),
                    file.path(RESULTS_DIR, "susiex_mvp", "susiex_variant_summary_mvp.csv"))
  susiex_files <- susiex_files[file.exists(susiex_files)]

  if (length(susiex_files) > 0) {
    cat("\n=== Integrating SuSiEX results (", length(susiex_files), "source(s):",
        paste(basename(susiex_files), collapse = ", "), ") ===\n")
    susiex <- rbindlist(lapply(susiex_files, function(f)
      fread(f, select = c("chr", "pos", "a1", "a2", "OVRL_PIP", "trait_pair"))),
      use.names = TRUE)
    setnames(susiex, c("chr", "pos", "a1", "a2", "OVRL_PIP", "trait_pair"),
             c("chromosome", "position", "allele1", "allele2", "susiex_pip", "trait"))
    susiex[, chromosome := as.integer(chromosome)]
    susiex[, position  := as.integer(position)]
    susiex[, trait     := toupper(trait)]
    # Deduplicate WITHIN trait: a variant may recur across loci for the same
    # trait, but we must NOT collapse across traits (that broadcast a
    # liver-enzyme OVRL_PIP=1.0 onto disease-GWAS rows). Keep max PIP per
    # (variant, trait); the join below is keyed on trait so each PIP stays with
    # its own study.
    susiex <- susiex[, .(susiex_pip = max(susiex_pip, na.rm = TRUE)),
                     by = .(chromosome, position, allele1, allele2, trait)]
    cat("  SuSiEX: ", nrow(susiex), " unique (variant, trait) entries loaded\n")
  } else {
    cat("  SuSiEX summary not found (legacy or _mvp)\n")
    susiex <- NULL
  }

  mesusie_files <- c(file.path(RESULTS_DIR, "mesusie", "mesusie_variant_summary.csv"),
                     file.path(RESULTS_DIR, "mesusie_mvp", "mesusie_variant_summary_mvp.csv"))
  mesusie_files <- mesusie_files[file.exists(mesusie_files)]

  if (length(mesusie_files) > 0) {
    cat("\n=== Integrating MESuSiE results (", length(mesusie_files), "source(s):",
        paste(basename(mesusie_files), collapse = ", "), ") ===\n")
    mesusie <- rbindlist(lapply(mesusie_files, function(f)
      fread(f, select = c("chr", "pos", "a1", "a2", "pip", "in_cs",
                          "cs_category", "trait_pair"))), use.names = TRUE)
    setnames(mesusie, c("chr", "pos", "a1", "a2", "pip", "trait_pair"),
             c("chromosome", "position", "allele1", "allele2", "mesusie_pip", "trait"))
    mesusie[, chromosome := as.integer(chromosome)]
    mesusie[, position  := as.integer(position)]
    mesusie[, trait     := toupper(trait)]
    # Shared CS = credible set spanning >=2 ancestries. Legacy encodes this as
    # "EUR_EAS"/"shared"; the MVP N-way run encodes the ancestry combo directly
    # (EUR_AFR, EAS_AFR_AMR, ...). Any underscore-joined multi-ancestry category
    # (or the literal "shared") counts; a single-ancestry category (EUR, AFR, AMR)
    # is ancestry-specific.
    mesusie[, mesusie_in_shared_cs := (in_cs == TRUE & grepl("_|shared", cs_category, ignore.case = TRUE))]
    # Deduplicate WITHIN trait (same broadcast hazard as SuSiEX): keep max PIP
    # and any shared CS hit per (variant, trait). The join below is keyed on
    # trait, so a PIP never lands on a different-trait GWAS row.
    mesusie <- mesusie[, .(mesusie_pip = max(mesusie_pip, na.rm = TRUE),
                           mesusie_in_shared_cs = any(mesusie_in_shared_cs, na.rm = TRUE)),
                       by = .(chromosome, position, allele1, allele2, trait)]
    cat("  MESuSiE:", nrow(mesusie), "unique (variant, trait) entries loaded\n")
  } else {
    cat("  MESuSiE summary not found (legacy or _mvp)\n")
    mesusie <- NULL
  }

  # Left-join joint PIPs onto combined (handle allele flips by trying both orientations)
  combined <- as.data.table(combined)

  # ---------------------------------------------------------------------------
  # Trait key for the joint-PIP join.
  #
  # The joint fine-mappers (SuSiEX / MESuSiE) were only run on the liver-enzyme
  # cross-ancestry meta-loci, keyed by `trait_pair` in {ALT, AST, GGT}. Disease /
  # PDFF GWAS (NAFLD, NASH, PDFF) have NO joint fine-map. If we join joint PIPs to
  # `combined` on coordinates alone, a liver-enzyme variant's OVRL_PIP (often 1.0)
  # is broadcast onto every study's row at that position — including the
  # disease-GWAS rows for PNPLA3 / TM6SF2, manufacturing a fake high PIP at the
  # exact loci a genetics reviewer checks first. Keying the join on (coords, trait)
  # confines each joint PIP to the matching liver-enzyme study and leaves disease
  # rows as NA (no joint fine-map exists for them).
  #
  # Trait is the study_name token (UKBB_ALT, BBJ_GGT, MVP_Albumin_EUR,
  # MVP_NAFLD_AMR, ...), matched case-insensitively and upper-cased to match the
  # toupper'd trait_pair on the joint tables. The MVP N-way run (2026-07-05) added
  # the Albumin/ChronLiver/Cirrhosis/NAFLD/Platelet traits, so they are matched too
  # (a joint PIP still only attaches to same-trait studies). regmatches() drops
  # non-matches, so derive per-row explicitly.
  combined[, trait := vapply(study, function(s) {
    m <- regmatches(s, regexpr("ALT|AST|GGT|Albumin|ChronLiver|Cirrhosis|NAFLD|NASH|Platelet|PDFF",
                               s, ignore.case = TRUE))
    if (length(m) == 0) NA_character_ else toupper(m)
  }, character(1))]

  if (!is.null(susiex)) {
    # Direct match — keyed on trait so liver-enzyme PIPs only attach to the
    # matching liver-enzyme study (disease/PDFF rows have trait=NA -> no match).
    combined <- merge(combined, susiex,
                      by = c("chromosome", "position", "allele1", "allele2", "trait"),
                      all.x = TRUE)
    # Allele-flipped match for unmatched rows
    unmatched <- is.na(combined$susiex_pip)
    if (any(unmatched)) {
      susiex_flip <- copy(susiex)
      setnames(susiex_flip, c("allele1", "allele2"), c("allele2", "allele1"))
      setnames(susiex_flip, "susiex_pip", "susiex_pip_flip")
      combined <- merge(combined, susiex_flip,
                        by = c("chromosome", "position", "allele1", "allele2", "trait"),
                        all.x = TRUE)
      combined[is.na(susiex_pip) & !is.na(susiex_pip_flip), susiex_pip := susiex_pip_flip]
      combined[, susiex_pip_flip := NULL]
    }
    cat("  SuSiEX matched:", sum(!is.na(combined$susiex_pip)), "/", nrow(combined), "variants\n")
  } else {
    combined[, susiex_pip := NA_real_]
  }

  if (!is.null(mesusie)) {
    # Direct match — keyed on trait (same broadcast guard as SuSiEX).
    combined <- merge(combined, mesusie,
                      by = c("chromosome", "position", "allele1", "allele2", "trait"),
                      all.x = TRUE)
    # Allele-flipped match for unmatched rows
    unmatched <- is.na(combined$mesusie_pip)
    if (any(unmatched)) {
      mesusie_flip <- copy(mesusie)
      setnames(mesusie_flip, c("allele1", "allele2"), c("allele2", "allele1"))
      setnames(mesusie_flip, c("mesusie_pip", "mesusie_in_shared_cs"),
               c("mesusie_pip_flip", "mesusie_shared_flip"))
      combined <- merge(combined, mesusie_flip,
                        by = c("chromosome", "position", "allele1", "allele2", "trait"),
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
