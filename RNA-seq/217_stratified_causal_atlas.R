#!/usr/bin/env Rscript
# 217_stratified_causal_atlas.R — Stratified Causal Atlas Integration
#
# Merges all Module A-D stratified causal results (Scripts 207-216) into the
# multi-evidence atlas, adding ~15 new columns that classify each gene by
# sex, subtype, progression, spatial, and pharmacogenomic stratification.
#
# Run AFTER all Module A-D scripts (207-216) have completed.
#
# Inputs (from RNA-seq/results/stratified_causal/):
#   - sex_causal_scores.csv (Script 207)
#   - subtype_causal_scores.csv (Script 208)
#   - s2_genetic_risk.csv (Script 208)
#   - progression_coloc_enrichment.csv (Script 210)
#   - onset_vs_progression_genes.csv (Script 210)
#   - progression_driver_genetics.csv (Script 211)
#   - spatial_coloc_enrichment.csv (Script 213)
#   - zone_causal_scores.csv (Script 213)
#   - pharmacogenomic_stratified_targets.csv (Script 215)
#   - Multi-evidence atlas: RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
#
# Outputs:
#   - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv (updated in-place)
#   - RNA-seq/results/stratified_causal/atlas_integration_summary.csv
#
# SLURM: cpu partition, 4 CPUs, 32GB, 48h

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

sc_dir   <- file.path(BASE, "RNA-seq/results/stratified_causal")
atlas_fp <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")

# ===========================================================================
# Helper: safely load a CSV, returning NULL + warning if missing
# ===========================================================================
safe_load <- function(path, label) {
  if (!file.exists(path)) {
    cat(sprintf("  WARNING: %s not found — skipping (%s)\n", label, basename(path)))
    return(NULL)
  }
  dt <- fread(path)
  cat(sprintf("  Loaded %s: %d rows x %d cols\n", label, nrow(dt), ncol(dt)))
  dt
}

# ===========================================================================
# Helper: standardize gene column to "gene"
# ===========================================================================
standardize_gene_col <- function(dt, label) {
  if (is.null(dt)) return(NULL)
  # Try common gene column names
  candidates <- c("gene", "human_symbol", "gene_symbol", "Gene", "symbol", "target_gene")
  found <- intersect(names(dt), candidates)
  if (length(found) == 0L) {
    cat(sprintf("  WARNING: No gene column found in %s (cols: %s) — skipping\n",
                label, paste(head(names(dt), 5), collapse = ", ")))
    return(NULL)
  }
  if (found[1] != "gene") setnames(dt, found[1], "gene")
  dt
}

# ===========================================================================
# 1. Load multi-evidence atlas
# ===========================================================================
cat("=" |> rep(72) |> paste(collapse = ""), "\n")
cat("Script 217: Stratified Causal Atlas Integration\n")
cat("=" |> rep(72) |> paste(collapse = ""), "\n\n")

cat("Loading multi-evidence atlas...\n")
atlas <- fread(atlas_fp)
cat(sprintf("  Atlas: %d genes x %d columns\n", nrow(atlas), ncol(atlas)))

# Standardize gene column
gene_col <- intersect(names(atlas), c("human_symbol", "gene", "gene_symbol"))[1]
if (gene_col != "gene") {
  setnames(atlas, gene_col, "gene")
}

# ===========================================================================
# 2. Remove existing stratified_causal columns if re-running (idempotency)
# ===========================================================================
strat_cols <- c(
  "sex_coloc_class", "sex_causal_score",
  "subtype_coloc_enriched", "s2_causal_score",    # deprecated aliases
  "s2_genetic_risk_prob",                          # deprecated alias
  "dominant_program", "dominant_program_code",    # k-program primary
  "dominant_program_for_gene", "dominant_program_logFC",
  "progression_coloc_class", "progression_causal_score",
  "zone_coloc_enriched", "zone_causal_score",
  "drug_target_stratification",
  "progression_driver_validated",
  "spatial_coloc_overlap",
  "lincs_reversal_stratum",
  # LD-panel parallel-atlas columns (added 2026-04-23 during LD diversification).
  # These mirror the canonical coloc_best_* columns but are populated from
  # Phase 4 (1KG EUR) and Phase 2 (TOP-LD EUR) production re-runs.
  "coloc_susie_best_pp4_1kg", "coloc_susie_best_gwas_1kg",
  "coloc_n_gwas_susie_h4_05_1kg", "coloc_n_gwas_susie_h4_08_1kg",
  "coloc_susie_best_pp4_topld", "coloc_susie_best_gwas_topld",
  "coloc_n_gwas_susie_h4_05_topld", "coloc_n_gwas_susie_h4_08_topld",
  # Monogenic convergence columns from 53_monogenic_convergence.R (Okur et al. 2026)
  "monogenic_panel_okur2026", "monogenic_pl_okur2026", "monogenic_mafld_okur2026",
  "monogenic_category_okur2026", "monogenic_indication_okur2026"
)
# Per-program score column pattern (added dynamically below)
prog_col_prefix <- "program_"
existing <- intersect(strat_cols, names(atlas))
if (length(existing) > 0L) {
  cat(sprintf("  Removing %d existing stratified_causal columns: %s\n",
              length(existing), paste(existing, collapse = ", ")))
  atlas[, (existing) := NULL]
}
# Also pattern-strip per-program columns from any prior k-run (prog_*_causal_score
# etc accumulate across runs with different k; must clean to prevent dup-col merge error)
# T0.2 (2026-04-22): broadened regex to capture .x/.y Cartesian-merge suffixes that
# leaked through 27a's bare-name left-join (fixed below).
prog_cols <- grep("^prog_.*_(causal_score|logFC)(\\.[xy])?$", names(atlas), value = TRUE)
dom_cols  <- grep("^dominant_program_(for_gene|logFC)(\\.[xy])?$", names(atlas), value = TRUE)
prog_cols <- union(prog_cols, dom_cols)
if (length(prog_cols) > 0L) {
  cat(sprintf("  Removing %d prior per-program columns (pattern): %s\n",
              length(prog_cols), paste(head(prog_cols, 8), collapse = ", ")))
  atlas[, (prog_cols) := NULL]
}

# ===========================================================================
# 3. Load Module A: Sex-stratified causal (Script 207)
# ===========================================================================
cat("\n--- Module A: Sex & Subtype ---\n")

sex_scores <- safe_load(file.path(sc_dir, "sex_causal_scores.csv"), "sex causal scores")
sex_scores <- standardize_gene_col(sex_scores, "sex_causal_scores")

if (!is.null(sex_scores)) {
  # Expected columns: gene, sex_coloc_class, sex_causal_score
  keep_cols <- intersect(names(sex_scores),
                         c("gene", "sex_coloc_class", "sex_causal_score"))
  if (length(keep_cols) > 1L) {
    sex_scores_u <- unique(sex_scores[, ..keep_cols], by = "gene")
    atlas <- merge(atlas, sex_scores_u, by = "gene", all.x = TRUE)
    atlas <- unique(atlas, by = "gene")
    cat(sprintf("  Merged sex scores: %d genes with values\n",
                sum(!is.na(atlas$sex_causal_score))))
  } else {
    cat("  WARNING: sex_causal_scores.csv missing expected columns\n")
    # Initialize as NA
    atlas[, sex_coloc_class := NA_character_]
    atlas[, sex_causal_score := NA_real_]
  }
} else {
  # Derive sex_coloc_class from existing atlas columns if possible
  cat("  Attempting fallback derivation from existing atlas columns...\n")
  has_sex_class <- "sex_class" %in% names(atlas)
  has_coloc     <- "coloc_susie_best_pp4" %in% names(atlas)

  if (has_sex_class && has_coloc) {
    atlas[, sex_coloc_class := ifelse(
      is.na(coloc_susie_best_pp4) | coloc_susie_best_pp4 < 0.5,
      "not_causal",
      fifelse(sex_class %in% c("Female_Specific", "female_specific"), "female_enriched",
              fifelse(sex_class %in% c("Male_Specific", "male_specific"), "male_enriched",
                      fifelse(sex_class %in% c("Sex_Divergent", "sex_divergent"), "sex_divergent",
                              "shared")))
    )]
    # Compute sex_causal_score: PP4 * |female_logFC - male_logFC|
    if (all(c("dream_logFC_F", "dream_logFC_M") %in% names(atlas))) {
      atlas[, sex_causal_score := fifelse(
        !is.na(coloc_susie_best_pp4) & !is.na(dream_logFC_F) & !is.na(dream_logFC_M),
        coloc_susie_best_pp4 * abs(dream_logFC_F - dream_logFC_M),
        NA_real_
      )]
    } else {
      atlas[, sex_causal_score := NA_real_]
    }
    n_causal <- sum(atlas$sex_coloc_class != "not_causal", na.rm = TRUE)
    cat(sprintf("  Fallback sex derivation: %d causal genes classified\n", n_causal))
  } else {
    atlas[, sex_coloc_class := NA_character_]
    atlas[, sex_causal_score := NA_real_]
    cat("  No fallback possible — columns initialized as NA\n")
  }
}

# ===========================================================================
# 4. Load Module A: k-Program causal (Script 208, refactored 2026-04-20)
# ===========================================================================
# Read per-program long-format causal scores and pivot to per-program columns.
prog_causal <- safe_load(file.path(sc_dir, "program_causal_scores.csv"),
                         "per-program causal scores")
program_labels_dt <- tryCatch(
  fread(file.path(BASE, "RNA-seq/results/subtypes/program_labels.csv")),
  error = function(e) NULL)

if (!is.null(prog_causal) && !is.null(program_labels_dt)) {
  # Script 208's CSV has 'gene' (Ensembl ID) AND 'human_symbol' — drop Ensembl
  # col first to prevent duplicate-name collision when standardize_gene_col
  # renames human_symbol -> gene.
  if (all(c("gene", "human_symbol") %in% names(prog_causal))) {
    prog_causal[, gene := NULL]
  }
  prog_causal <- standardize_gene_col(prog_causal, "program_causal_scores")
  # Bring in biological label
  label_map <- program_labels_dt[, .(program_code, biological_label)]
  prog_causal <- merge(prog_causal, label_map, by = "program_code", all.x = TRUE)
  prog_causal[, label_clean := gsub("[^A-Za-z0-9]+", "_", biological_label)]

  # Wide-pivot signed causal score per program label
  wide_causal <- dcast(prog_causal[, .(gene, label_clean, program_causal_score)],
                       gene ~ label_clean, value.var = "program_causal_score",
                       fun.aggregate = function(x) x[1])
  causal_cols <- setdiff(names(wide_causal), "gene")
  setnames(wide_causal, causal_cols, paste0("prog_", causal_cols, "_causal_score"))

  wide_lfc <- dcast(prog_causal[, .(gene, label_clean, logFC)],
                    gene ~ label_clean, value.var = "logFC",
                    fun.aggregate = function(x) x[1])
  lfc_cols <- setdiff(names(wide_lfc), "gene")
  setnames(wide_lfc, lfc_cols, paste0("prog_", lfc_cols, "_logFC"))

  # Defensive dedup: dcast+aggregate should give unique gene rows, but merge
  # fails cartesian if atlas$gene already has duplicates from upstream merges.
  wide_causal <- unique(wide_causal, by = "gene")
  wide_lfc    <- unique(wide_lfc,    by = "gene")
  atlas <- unique(atlas, by = "gene")
  atlas <- merge(atlas, wide_causal, by = "gene", all.x = TRUE)
  atlas <- merge(atlas, wide_lfc,    by = "gene", all.x = TRUE)
  cat(sprintf("  Merged per-program causal: %d programs x gene-level scores\n",
              length(causal_cols)))

  # Dominant program per gene = program with max |logFC| in the causal table
  dom_per_gene <- prog_causal[, .SD[which.max(abs(logFC))], by = gene,
                              .SDcols = c("program_code","biological_label","logFC")]
  setnames(dom_per_gene, c("biological_label","logFC"),
           c("dominant_program_for_gene","dominant_program_logFC"))
  atlas <- merge(atlas, dom_per_gene[, .(gene, dominant_program_for_gene,
                                          dominant_program_logFC)],
                 by = "gene", all.x = TRUE)

  # DEPRECATED ALIAS: s2_causal_score = Fibrogenic program column (if exists)
  # T0.1 rename (2026-04-21): Fibrotic -> Fibrogenic; regex updated accordingly.
  fib_col <- names(atlas)[grepl("^prog_Fibrogenic_", names(atlas)) &
                           grepl("_causal_score$", names(atlas))]
  if (length(fib_col) >= 1) {
    atlas[, s2_causal_score := atlas[[fib_col[1]]]]
    atlas[, subtype_coloc_enriched := fifelse(
      !is.na(s2_causal_score) & s2_causal_score > 0,
      "Fibrogenic_program", NA_character_
    )]
  } else {
    atlas[, s2_causal_score := NA_real_]
    atlas[, subtype_coloc_enriched := NA_character_]
  }
} else {
  cat("  WARNING: program_causal_scores.csv or program_labels.csv missing — initializing NA\n")
  atlas[, s2_causal_score := NA_real_]
  atlas[, subtype_coloc_enriched := NA_character_]
}

# Per-sample genetic risk probabilities (Script 208 multinomial): gene-atlas level
# doesn't have sample predictions, so just NA this column (was gene-level mean
# in legacy but poorly defined). Kept for schema compatibility.
if (!"s2_genetic_risk_prob" %in% names(atlas)) {
  atlas[, s2_genetic_risk_prob := NA_real_]
}

# Ensure atlas remains unique-gene-keyed before subsequent merges
atlas <- unique(atlas, by = "gene")

# ===========================================================================
# 5. Load Module B: Progression-stratified causal (Scripts 210-211)
# ===========================================================================
cat("\n--- Module B: Progression ---\n")

prog_enrich <- safe_load(file.path(sc_dir, "progression_coloc_enrichment.csv"),
                         "progression COLOC enrichment")
prog_enrich <- standardize_gene_col(prog_enrich, "progression_coloc_enrichment")

onset_prog <- safe_load(file.path(sc_dir, "onset_vs_progression_genes.csv"),
                        "onset vs progression genes")
onset_prog <- standardize_gene_col(onset_prog, "onset_vs_progression_genes")

prog_drivers <- safe_load(file.path(sc_dir, "progression_driver_genetics.csv"),
                          "progression driver genetics")
prog_drivers <- standardize_gene_col(prog_drivers, "progression_driver_genetics")

# Merge onset/progression classification
if (!is.null(onset_prog)) {
  # Actual column names: progression_class (not progression_coloc_class)
  keep_cols <- intersect(names(onset_prog),
                         c("gene", "progression_class", "progression_coloc_class",
                           "progression_causal_score"))
  if (length(keep_cols) > 1L) {
    atlas <- merge(atlas, onset_prog[, ..keep_cols], by = "gene", all.x = TRUE)
    # Rename if source uses progression_class instead of progression_coloc_class
    if ("progression_class" %in% names(atlas) && !"progression_coloc_class" %in% names(atlas)) {
      setnames(atlas, "progression_class", "progression_coloc_class")
    }
    cat(sprintf("  Merged progression class: %d genes classified\n",
                sum(!is.na(atlas$progression_coloc_class))))
  } else {
    atlas[, progression_coloc_class := NA_character_]
    atlas[, progression_causal_score := NA_real_]
  }
} else {
  # Fallback: derive from existing atlas columns
  cat("  Attempting fallback derivation from existing atlas columns...\n")
  has_pp4 <- "coloc_susie_best_pp4" %in% names(atlas)
  has_tau <- "progression_tau" %in% names(atlas)
  has_peak <- "progression_peak_contrast" %in% names(atlas)

  if (has_pp4 && has_tau && has_peak) {
    atlas[, progression_coloc_class := {
      pc <- fifelse(
        is.na(coloc_susie_best_pp4) | coloc_susie_best_pp4 < 0.5,
        "not_causal",
        fifelse(grepl("F0|NAFL|ctrl", progression_peak_contrast, ignore.case = TRUE),
                "onset",
                fifelse(grepl("F3|F4|cirr", progression_peak_contrast, ignore.case = TRUE),
                        "progression",
                        "pan")))
      pc
    }]
    atlas[, progression_causal_score := fifelse(
      !is.na(coloc_susie_best_pp4) & !is.na(progression_tau),
      coloc_susie_best_pp4 * progression_tau,
      NA_real_
    )]
    n_onset <- sum(atlas$progression_coloc_class == "onset", na.rm = TRUE)
    n_prog  <- sum(atlas$progression_coloc_class == "progression", na.rm = TRUE)
    cat(sprintf("  Fallback: %d onset, %d progression, %d pan causal genes\n",
                n_onset, n_prog,
                sum(atlas$progression_coloc_class == "pan", na.rm = TRUE)))
  } else {
    atlas[, progression_coloc_class := NA_character_]
    atlas[, progression_causal_score := NA_real_]
    cat("  No fallback possible — columns initialized as NA\n")
  }
}

# Merge progression driver validation flag
if (!is.null(prog_drivers)) {
  # Actual column name is genetically_validated (from Script 211)
  keep_cols <- intersect(names(prog_drivers),
                         c("gene", "genetically_validated", "progression_driver_validated"))
  if (length(keep_cols) > 1L) {
    # Dedup defensively
    prog_drivers_u <- unique(prog_drivers[, ..keep_cols], by = "gene")
    atlas <- unique(atlas, by = "gene")
    atlas <- merge(atlas, prog_drivers_u, by = "gene", all.x = TRUE)
    if ("genetically_validated" %in% names(atlas) && !"progression_driver_validated" %in% names(atlas)) {
      setnames(atlas, "genetically_validated", "progression_driver_validated")
    }
  }
}
if (!"progression_driver_validated" %in% names(atlas)) {
  atlas[, progression_driver_validated := NA]
}

# ===========================================================================
# 6. Load Module C: Spatial / zonation causal (Script 213)
# ===========================================================================
cat("\n--- Module C: Spatial / Zonation ---\n")

spatial_enrich <- safe_load(file.path(sc_dir, "spatial_coloc_enrichment.csv"),
                            "spatial COLOC enrichment")
spatial_enrich <- standardize_gene_col(spatial_enrich, "spatial_coloc_enrichment")

zone_scores <- safe_load(file.path(sc_dir, "zone_causal_scores.csv"), "zone causal scores")
zone_scores <- standardize_gene_col(zone_scores, "zone_causal_scores")

if (!is.null(zone_scores)) {
  keep_cols <- intersect(names(zone_scores),
                         c("gene", "zone_coloc_class", "zone_coloc_enriched", "zone_causal_score"))
  if (length(keep_cols) > 1L) {
    atlas <- merge(atlas, zone_scores[, ..keep_cols], by = "gene", all.x = TRUE)
    # Harmonize column name: zone_coloc_class -> zone_coloc_enriched
    if ("zone_coloc_class" %in% names(atlas) && !"zone_coloc_enriched" %in% names(atlas)) {
      setnames(atlas, "zone_coloc_class", "zone_coloc_enriched")
    }
    cat(sprintf("  Merged zone scores: %d genes with values\n",
                sum(!is.na(atlas$zone_causal_score))))
  } else {
    atlas[, zone_coloc_enriched := NA_character_]
    atlas[, zone_causal_score := NA_real_]
  }
} else {
  # Fallback: derive from existing zonation_class + COLOC
  cat("  Attempting fallback derivation from existing atlas columns...\n")
  has_zone <- "zonation_class" %in% names(atlas)
  has_pp4  <- "coloc_susie_best_pp4" %in% names(atlas)

  if (has_zone && has_pp4) {
    atlas[, zone_coloc_enriched := fifelse(
      is.na(coloc_susie_best_pp4) | coloc_susie_best_pp4 < 0.5,
      NA_character_,
      fifelse(grepl("[Pp]eriportal", zonation_class), "periportal",
              fifelse(grepl("[Pp]ericentral", zonation_class), "pericentral",
                      fifelse(!is.na(zonation_class) & zonation_class != "" &
                                !grepl("^NA$", zonation_class),
                              "pan", NA_character_)))
    )]
    atlas[, zone_causal_score := fifelse(
      !is.na(coloc_susie_best_pp4) & !is.na(zonation_class) &
        zonation_class != "" & !grepl("^NA$", zonation_class),
      coloc_susie_best_pp4,
      NA_real_
    )]
    n_zone <- sum(!is.na(atlas$zone_coloc_enriched), na.rm = TRUE)
    cat(sprintf("  Fallback zone derivation: %d causal genes with zonation\n", n_zone))
  } else {
    atlas[, zone_coloc_enriched := NA_character_]
    atlas[, zone_causal_score := NA_real_]
    cat("  No fallback possible — columns initialized as NA\n")
  }
}

# Spatial COLOC overlap flag from enrichment file
if (!is.null(spatial_enrich)) {
  keep_cols <- intersect(names(spatial_enrich), c("gene", "spatial_coloc_overlap"))
  if (length(keep_cols) > 1L) {
    atlas <- merge(atlas, spatial_enrich[, ..keep_cols], by = "gene", all.x = TRUE)
  }
}
if (!"spatial_coloc_overlap" %in% names(atlas)) {
  atlas[, spatial_coloc_overlap := NA]
}

# ===========================================================================
# 7. Load Module D: Pharmacogenomic stratification (Script 215)
# ===========================================================================
cat("\n--- Module D: Pharmacogenomics ---\n")

pharma <- safe_load(file.path(sc_dir, "pharmacogenomic_stratified_targets.csv"),
                    "pharmacogenomic targets")
pharma <- standardize_gene_col(pharma, "pharmacogenomic_stratified_targets")

if (!is.null(pharma)) {
  keep_cols <- intersect(names(pharma),
                         c("gene", "drug_target_stratification", "lincs_reversal_stratum"))
  if (length(keep_cols) > 1L) {
    atlas <- merge(atlas, pharma[, ..keep_cols], by = "gene", all.x = TRUE)
    cat(sprintf("  Merged pharma: %d genes with stratification\n",
                sum(!is.na(atlas$drug_target_stratification))))
  } else {
    atlas[, drug_target_stratification := NA_character_]
    atlas[, lincs_reversal_stratum := NA_character_]
  }
} else {
  # Build drug_target_stratification from available info
  cat("  Attempting fallback derivation from existing atlas columns...\n")
  has_drug <- "dgidb_druggable" %in% names(atlas) | "opentargets_drug" %in% names(atlas)

  if (has_drug && "sex_coloc_class" %in% names(atlas) &&
      "progression_coloc_class" %in% names(atlas)) {
    # Build summary string: "sex=female_enriched;stage=onset;subtype=S2"
    atlas[, drug_target_stratification := {
      # Only annotate druggable genes with causal evidence
      is_drug <- (!is.na(dgidb_druggable) & dgidb_druggable == TRUE) |
                 (!is.na(opentargets_drug) & opentargets_drug != "")
      is_causal <- !is.na(coloc_susie_best_pp4) & coloc_susie_best_pp4 >= 0.5

      strat <- fifelse(
        is_drug & is_causal,
        paste0(
          "sex=", fifelse(!is.na(sex_coloc_class) & sex_coloc_class != "not_causal",
                          sex_coloc_class, "unspec"),
          ";stage=", fifelse(!is.na(progression_coloc_class) &
                              progression_coloc_class != "not_causal",
                            progression_coloc_class, "unspec"),
          fifelse(!is.na(subtype_coloc_enriched),
                  paste0(";subtype=", subtype_coloc_enriched), "")
        ),
        NA_character_
      )
      strat
    }]
    n_strat <- sum(!is.na(atlas$drug_target_stratification))
    cat(sprintf("  Fallback pharma derivation: %d druggable causal genes stratified\n",
                n_strat))
  } else {
    atlas[, drug_target_stratification := NA_character_]
    cat("  No fallback possible — column initialized as NA\n")
  }
  atlas[, lincs_reversal_stratum := NA_character_]
}

# ===========================================================================
# 8. Rename gene column back to human_symbol
# ===========================================================================
setnames(atlas, "gene", "human_symbol")

# ===========================================================================
# 9. Write updated atlas
# ===========================================================================
cat("\n--- Writing updated atlas ---\n")

# Backup current atlas
backup_fp <- sub("\\.csv$", sprintf("_backup_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")),
                 atlas_fp)
file.copy(atlas_fp, backup_fp)
cat(sprintf("  Backup: %s\n", backup_fp))

fwrite(atlas, atlas_fp)
cat(sprintf("  Updated atlas: %d genes x %d columns\n", nrow(atlas), ncol(atlas)))

# ===========================================================================
# 10. Summary statistics
# ===========================================================================
cat("\n--- Summary Statistics ---\n")

# Count new columns
new_cols <- intersect(strat_cols, names(atlas))
cat(sprintf("  New stratified columns added: %d\n", length(new_cols)))
for (col in new_cols) {
  n_nona <- sum(!is.na(atlas[[col]]))
  cat(sprintf("    %-35s: %d non-NA (%.1f%%)\n", col, n_nona,
              100 * n_nona / nrow(atlas)))
}

# Sex classification breakdown
if ("sex_coloc_class" %in% names(atlas)) {
  cat("\n  Sex-COLOC classification:\n")
  tbl <- table(atlas$sex_coloc_class, useNA = "ifany")
  for (nm in names(tbl)) {
    cat(sprintf("    %-25s: %d\n", ifelse(is.na(nm), "<NA>", nm), tbl[nm]))
  }
}

# Progression classification breakdown
if ("progression_coloc_class" %in% names(atlas)) {
  cat("\n  Progression-COLOC classification:\n")
  tbl <- table(atlas$progression_coloc_class, useNA = "ifany")
  for (nm in names(tbl)) {
    cat(sprintf("    %-25s: %d\n", ifelse(is.na(nm), "<NA>", nm), tbl[nm]))
  }
}

# Zone classification breakdown
if ("zone_coloc_enriched" %in% names(atlas)) {
  cat("\n  Zone-COLOC classification:\n")
  tbl <- table(atlas$zone_coloc_enriched, useNA = "ifany")
  for (nm in names(tbl)) {
    cat(sprintf("    %-25s: %d\n", ifelse(is.na(nm), "<NA>", nm), tbl[nm]))
  }
}

# Drug target stratification count
if ("drug_target_stratification" %in% names(atlas)) {
  n_strat <- sum(!is.na(atlas$drug_target_stratification))
  cat(sprintf("\n  Druggable causal genes with stratification: %d\n", n_strat))
}

# ===========================================================================
# 11. Write integration summary
# ===========================================================================
summary_dt <- data.table(
  column = new_cols,
  n_nona = vapply(new_cols, function(c) sum(!is.na(atlas[[c]])), integer(1)),
  pct_nona = vapply(new_cols, function(c) round(100 * sum(!is.na(atlas[[c]])) / nrow(atlas), 2),
                    numeric(1)),
  source = vapply(new_cols, function(c) {
    if (grepl("sex", c))        return("Module A: Sex (207)")
    if (grepl("subtype|s2", c)) return("Module A: Subtype (208)")
    if (grepl("prog", c))       return("Module B: Progression (210-211)")
    if (grepl("zone|spatial", c)) return("Module C: Spatial (213)")
    if (grepl("drug|lincs", c)) return("Module D: Pharma (215)")
    "Unknown"
  }, character(1))
)

summary_fp <- file.path(sc_dir, "atlas_integration_summary.csv")
fwrite(summary_dt, summary_fp)
cat(sprintf("\n  Integration summary: %s\n", summary_fp))

cat("\n")
cat("=" |> rep(72) |> paste(collapse = ""), "\n")
cat("Script 217 COMPLETE\n")
cat(sprintf("  Atlas: %s (%d genes x %d cols)\n", atlas_fp, nrow(atlas), ncol(atlas)))
cat(sprintf("  Summary: %s\n", summary_fp))
cat("=" |> rep(72) |> paste(collapse = ""), "\n")
