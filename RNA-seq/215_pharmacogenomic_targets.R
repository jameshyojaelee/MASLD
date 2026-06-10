#!/usr/bin/env Rscript
# 215_pharmacogenomic_targets.R — Pharmacogenomic target x stratification mapping
#
# Part of the Stratified Causal Architecture Pipeline (Module D).
# Maps 10 clinical drug targets onto stratification axes (COLOC, sex, stage,
# subtype, GWAS-ATAC), then stratifies 133 CGP reversal gene sets by sex,
# stage, and subtype overlap to produce a precision medicine matrix.
#
# Inputs:
#   - Clinical drug targets:   RNA-seq/results/drug_repurposing/clinical_drug_validation_table.csv
#   - COLOC gene-level:        GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv
#   - Sex DEG classification:  .../results/integration/sex_deg_classification.csv
#   - Stage-unique signatures: .../results/staging_classifier/stage_unique_signatures.csv
#   - NMF subtype assignments: RNA-seq/results/subtypes/nmf_assignments.csv
#   - Subtype markers:         RNA-seq/results/subtypes/subtype_markers.csv
#   - CGP reversal hits:       RNA-seq/results/drug_repurposing/cgp_reversal_hits.csv
#   - GWAS-ATAC motif:         GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv
#   - Gene-level GWAS-ATAC:    GWAS/finemapping/results/gwas_atac/gene_level_gwas_atac.csv
#   - Canonical DEGs:          .../results/integration/canonical_deg_results.csv
#   - Multi-evidence atlas:    RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
#
# Outputs (in RNA-seq/results/stratified_causal/):
#   - pharmacogenomic_stratified_targets.csv  Per-drug-target stratified evidence
#   - lincs_stratified_reversal.csv           Per-CGP-hit sex/stage/subtype reversal
#   - precision_medicine_matrix.csv           Long-format drug x patient-stratum matrix
#
# SLURM: cpu partition, 4 CPUs, 32GB, 48h

suppressPackageStartupMessages({
  library(data.table)
})

cat("=== 215_pharmacogenomic_targets.R ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

outdir <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

integ_dir <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
stag_dir <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier")

# ===========================================================================
# 1. Load multi-evidence atlas for Ensembl <-> symbol mapping
# ===========================================================================
cat("=== Loading atlas for gene ID mapping ===\n")
atlas <- fread(file.path(BASE,
  "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
  select = c("ensembl_id", "human_symbol"))
atlas[, ensembl_base := sub("\\.[0-9]+$", "", ensembl_id)]
symbol_map <- atlas[human_symbol != "" & !is.na(human_symbol),
                    .(ensembl_base, symbol = human_symbol)]
symbol_map <- symbol_map[!duplicated(ensembl_base)]
cat("  Gene ID mapping:", nrow(symbol_map), "Ensembl -> symbol pairs\n")

# Reverse map: symbol -> ensembl_base
sym2ens <- setNames(symbol_map$ensembl_base, symbol_map$symbol)

# ===========================================================================
# 2. Load clinical drug targets
# ===========================================================================
cat("\n=== Loading clinical drug targets ===\n")
drugs <- fread(file.path(BASE,
  "RNA-seq/results/drug_repurposing/clinical_drug_validation_table.csv"))
cat("  Drug-target rows:", nrow(drugs), "\n")
cat("  Unique drugs:", length(unique(drugs$drug)), "\n")
cat("  Unique target genes:", length(unique(drugs$target_gene)), "\n")
target_genes <- unique(drugs$target_gene)
cat("  Targets:", paste(target_genes, collapse = ", "), "\n")

# ===========================================================================
# 3. Load COLOC gene-level results
# ===========================================================================
cat("\n=== Loading COLOC gene-level results ===\n")
coloc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc <- coloc[gene != "" & !is.na(gene)]
# T0.4 (2026-04-22): prefer SuSiE PP4 over legacy ABF.
if ("coloc_best_susie_pp4" %in% names(coloc)) {
  coloc[, coloc_best_abf_pp4 := coloc_best_pp4]
  coloc[, coloc_best_pp4 := fifelse(!is.na(coloc_best_susie_pp4),
                                    coloc_best_susie_pp4, coloc_best_pp4)]
  cat("  T0.4: coloc_best_pp4 now sourced from SuSiE with ABF fallback\n")
}
cat("  Total genes in COLOC:", nrow(coloc), "\n")
cat("  PP.H4 > 0.9:", sum(coloc$coloc_best_pp4 > 0.9), "\n")
cat("  PP.H4 > 0.8:", sum(coloc$coloc_best_pp4 > 0.8), "\n")

# ===========================================================================
# 4. Load sex DEG classification
# ===========================================================================
cat("\n=== Loading sex DEG classification ===\n")
# v3 mashr Bayesian preferred; v2 fallback. v3 CSV provides v2-compatible `sex_class` alias.
sex_v3_path <- file.path(integ_dir, "sex_v3", "sex_deg_classification_v3.csv")
sex_v2_path <- file.path(integ_dir, "sex_deg_classification.csv")
sex_deg_file <- if (file.exists(sex_v3_path)) sex_v3_path else sex_v2_path
cat("  Source:", basename(dirname(sex_deg_file)), "/", basename(sex_deg_file), "\n")
sex_deg <- fread(sex_deg_file)
# v3 column compatibility shim: v3 lacks `lfc_diff` (use abs(logFC_F - logFC_M)) and
# names interaction logFC as `beta_interaction`.
if (!"lfc_diff" %in% names(sex_deg) && all(c("logFC_F", "logFC_M") %in% names(sex_deg))) {
  sex_deg[, lfc_diff := logFC_F - logFC_M]
}
if (!"interaction_logFC" %in% names(sex_deg) && "beta_interaction" %in% names(sex_deg)) {
  sex_deg[, interaction_logFC := beta_interaction]
}
sex_deg[, ensembl_base := sub("\\.[0-9]+$", "", gene)]
sex_deg <- merge(sex_deg, symbol_map, by = "ensembl_base", all.x = TRUE)
cat("  Total genes:", nrow(sex_deg), "\n")
cat("  Sex class distribution:\n")
print(sex_deg[, .N, by = sex_class][order(-N)])

# ===========================================================================
# 5. Load stage-unique signatures
# ===========================================================================
cat("\n=== Loading stage-unique signatures ===\n")
stage_uniq <- fread(file.path(stag_dir, "stage_unique_signatures.csv"))
stage_uniq[, ensembl_base := sub("\\.[0-9]+$", "", gene)]
stage_uniq <- merge(stage_uniq, symbol_map, by = "ensembl_base", all.x = TRUE)
cat("  Total stage-unique genes:", nrow(stage_uniq), "\n")
cat("  By axis:\n")
print(stage_uniq[, .N, by = staging_axis])

# ===========================================================================
# 6. Load NMF subtype markers (S1 vs S2 marker genes)
# ===========================================================================
cat("\n=== Loading NMF subtype markers ===\n")
sub_markers <- fread(file.path(BASE,
  "RNA-seq/results/subtypes/subtype_markers.csv"))
sub_markers[, ensembl_base := sub("\\.[0-9]+$", "", gene)]
sub_markers <- merge(sub_markers, symbol_map, by = "ensembl_base", all.x = TRUE)
cat("  Total subtype marker genes:", nrow(sub_markers), "\n")
cat("  By subtype:\n")
print(sub_markers[, .N, by = subtype])

# Extract S1 and S2 marker gene sets (Ensembl base IDs)
s1_markers_ens <- sub_markers[subtype == "S1", unique(ensembl_base)]
s2_markers_ens <- sub_markers[subtype == "S2", unique(ensembl_base)]
s1_markers_sym <- sub_markers[subtype == "S1" & !is.na(symbol), unique(symbol)]
s2_markers_sym <- sub_markers[subtype == "S2" & !is.na(symbol), unique(symbol)]
cat("  S1 markers (symbols):", length(s1_markers_sym), "\n")
cat("  S2 markers (symbols):", length(s2_markers_sym), "\n")

# ===========================================================================
# 7. Load CGP reversal hits
# ===========================================================================
cat("\n=== Loading CGP reversal hits ===\n")
cgp <- fread(file.path(BASE,
  "RNA-seq/results/drug_repurposing/cgp_reversal_hits.csv"))
cat("  Total CGP hits:", nrow(cgp), "\n")
cat("  CGP class distribution:\n")
print(cgp[, .N, by = cgp_class][order(-N)])

# ===========================================================================
# 8. Load GWAS-ATAC motif disruption
# ===========================================================================
cat("\n=== Loading GWAS-ATAC motif disruption ===\n")
motif <- fread(file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv"))
cat("  Total motif disruption entries:", nrow(motif), "\n")
cat("  Unique TFs disrupted:", length(unique(motif$tf_name)), "\n")

# Also load gene-level GWAS-ATAC
gene_atac <- fread(file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/gene_level_gwas_atac.csv"))
gene_atac[, ensembl_base := sub("\\.[0-9]+$", "", assigned_gene)]
cat("  Gene-level GWAS-ATAC entries:", nrow(gene_atac), "\n")

# ===========================================================================
# 9. Load dream results for the full disease t-statistic signature
# ===========================================================================
cat("\n=== Loading dream results ===\n")
dream <- fread(file.path(integ_dir, "canonical_deg_results.csv"))
cat("  Total dream genes:", nrow(dream), "\n")
# Canonical DEG table carries its own `symbol`; drop it so the symbol_map merge
# (ENSEMBL -> symbol) below does not collide with a `symbol.x`/`symbol.y` split.
if ("symbol" %in% names(dream)) dream[, symbol := NULL]
dream[, ensembl_base := sub("\\.[0-9]+$", "", gene)]
# Merge symbols
dream <- merge(dream, symbol_map, by = "ensembl_base", all.x = TRUE)
cat("  Mapped to symbols:", sum(!is.na(dream$symbol)), "\n")

# Build t-stat signature (used for LINCS stratification)
tstat_sig <- dream[, .(ensembl_base, symbol, t, logFC, padj)]

# ===========================================================================
# PART A: Per-drug-target stratified annotation
# ===========================================================================
cat("\n===================================================\n")
cat("PART A: Pharmacogenomic stratified target annotation\n")
cat("===================================================\n\n")

results_list <- list()

for (i in seq_len(nrow(drugs))) {
  tg <- drugs$target_gene[i]
  drug_name <- drugs$drug[i]
  cat("  Processing:", drug_name, "->", tg, "\n")

  row <- data.table(
    target_gene = tg,
    drug = drug_name,
    stage = drugs$stage[i],
    moa = drugs$moa[i],
    atlas_support = drugs$atlas_support[i]
  )

  # --- COLOC evidence ---
  coloc_row <- coloc[gene == tg]
  if (nrow(coloc_row) > 0) {
    pp4 <- coloc_row$coloc_best_pp4[1]
    row[, coloc_pp4 := pp4]
    row[, coloc_tier := fifelse(pp4 > 0.9, "Strong",
                        fifelse(pp4 > 0.8, "Moderate", "Weak"))]
    row[, coloc_best_gwas := coloc_row$coloc_best_gwas[1]]
    row[, coloc_n_gwas_h4_05 := coloc_row$coloc_n_gwas_h4_05[1]]
  } else {
    row[, c("coloc_pp4", "coloc_tier", "coloc_best_gwas", "coloc_n_gwas_h4_05") :=
      list(NA_real_, "Absent", NA_character_, 0L)]
  }

  # --- Direction-of-effect consistency (review A03; Nelson 2015 / Minikel 2024) ---
  # A therapeutic should OPPOSE the disease dysregulation: inhibit an up-in-disease target or
  # agonize a down-in-disease one. We use the disease expression direction (dream logFC) as the
  # available proxy for the genetic direction of effect; the gold standard is the lead-CS eQTL beta
  # aligned to the GWAS risk allele (NOT yet plumbed here — follow-up). Pairs where the drug pushes
  # the SAME direction as the disease (e.g. agonist of an up-regulated target) get FALSE. THRB
  # (down in disease; resmetirom is an agonist) validates as TRUE.
  .moa <- tolower(as.character(row$moa))
  drug_dir <- fifelse(grepl("agonist|activator|positive allosteric|inducer|stimulat", .moa), 1,
              fifelse(grepl("antagonist|inhibitor|blocker|negative allosteric|suppress|degrader", .moa), -1,
                      NA_real_))
  dis_lfc <- tstat_sig[symbol == tg, logFC][1]
  row[, drug_direction := drug_dir]
  row[, disease_logFC := dis_lfc]
  row[, direction_consistent := fifelse(is.na(drug_dir) | is.na(dis_lfc), NA, (drug_dir * dis_lfc) < 0)]

  # --- Sex specificity ---
  ens_id <- sym2ens[tg]
  if (!is.na(ens_id)) {
    sex_row <- sex_deg[ensembl_base == ens_id]
    if (nrow(sex_row) > 0) {
      row[, sex_class := sex_row$sex_class[1]]
      row[, sex_logFC_F := sex_row$logFC_F[1]]
      row[, sex_logFC_M := sex_row$logFC_M[1]]
      row[, sex_lfc_diff := sex_row$lfc_diff[1]]
    } else {
      row[, c("sex_class", "sex_logFC_F", "sex_logFC_M", "sex_lfc_diff") :=
        list("Not_found", NA_real_, NA_real_, NA_real_)]
    }
  } else {
    row[, c("sex_class", "sex_logFC_F", "sex_logFC_M", "sex_lfc_diff") :=
      list("Not_mapped", NA_real_, NA_real_, NA_real_)]
  }

  # --- Stage specificity ---
  if (!is.na(ens_id)) {
    stg_row <- stage_uniq[ensembl_base == ens_id]
    if (nrow(stg_row) > 0) {
      row[, stage_unique := stg_row$is_stage_unique[1]]
      row[, stage_tau := stg_row$tau[1]]
      row[, stage_most_specific := stg_row$most_specific_stage[1]]
      row[, stage_axis := stg_row$staging_axis[1]]
      row[, stage_max_lfc := stg_row$max_abs_lfc[1]]
      # Classify: early (F0/F1/NAS0-2), late (F3/F4/NAS5-7), or pan
      spec <- stg_row$most_specific_stage[1]
      if (grepl("F0|F1|NAS0|NAS1|NAS2", spec)) {
        row[, stage_class := "early"]
      } else if (grepl("F3|F4|NAS[5-8]", spec)) {
        row[, stage_class := "late"]
      } else {
        row[, stage_class := "mid"]
      }
    } else {
      row[, c("stage_unique", "stage_tau", "stage_most_specific",
              "stage_axis", "stage_max_lfc", "stage_class") :=
        list(FALSE, NA_real_, "not_stage_unique", NA_character_, NA_real_, "pan")]
    }
  } else {
    row[, c("stage_unique", "stage_tau", "stage_most_specific",
            "stage_axis", "stage_max_lfc", "stage_class") :=
      list(NA, NA_real_, "not_mapped", NA_character_, NA_real_, NA_character_)]
  }

  # --- Subtype specificity ---
  if (!is.na(ens_id)) {
    sub_row <- sub_markers[ensembl_base == ens_id]
    if (nrow(sub_row) > 0) {
      row[, subtype_marker := sub_row$subtype[1]]
      row[, subtype_direction := sub_row$direction[1]]
      row[, subtype_logFC := sub_row$logFC[1]]
      row[, subtype_padj := sub_row$adj.P.Val[1]]
    } else {
      row[, c("subtype_marker", "subtype_direction",
              "subtype_logFC", "subtype_padj") :=
        list("Neither", NA_character_, NA_real_, NA_real_)]
    }
  } else {
    row[, c("subtype_marker", "subtype_direction",
            "subtype_logFC", "subtype_padj") :=
      list("Not_mapped", NA_character_, NA_real_, NA_real_)]
  }

  # --- GWAS-ATAC motif disruption ---
  # Check if this target gene's TF motif is disrupted by any GWAS variant
  motif_hits <- motif[toupper(tf_name) == toupper(tg)]
  row[, n_motif_disruptions := nrow(motif_hits)]
  if (nrow(motif_hits) > 0) {
    row[, motif_max_pip := max(motif_hits$max_pip, na.rm = TRUE)]
    row[, motif_in_disease_regulon := any(motif_hits$motif_in_disease_regulon == TRUE)]
    row[, motif_best_effect := motif_hits[which.max(abs(as.numeric(alleleDiff))),
                                          effect][1]]
    row[, atac_evidence := fifelse(any(motif_hits$motif_in_disease_regulon == TRUE),
                                   "Strong", "Present")]
  } else {
    # Also check gene-level GWAS-ATAC (variant overlapping this gene's peaks)
    if (!is.na(ens_id)) {
      atac_row <- gene_atac[ensembl_base == ens_id]
      if (nrow(atac_row) > 0 && atac_row$gwas_variant_in_peak[1] == TRUE) {
        row[, c("motif_max_pip", "motif_in_disease_regulon",
                "motif_best_effect", "atac_evidence") :=
          list(atac_row$gwas_max_pip_in_peak[1], FALSE, NA_character_, "Weak")]
      } else {
        row[, c("motif_max_pip", "motif_in_disease_regulon",
                "motif_best_effect", "atac_evidence") :=
          list(NA_real_, FALSE, NA_character_, "Absent")]
      }
    } else {
      row[, c("motif_max_pip", "motif_in_disease_regulon",
              "motif_best_effect", "atac_evidence") :=
        list(NA_real_, FALSE, NA_character_, "Absent")]
    }
  }

  results_list[[i]] <- row
}

pharma_targets <- rbindlist(results_list, fill = TRUE)
cat("\n  Pharmacogenomic target summary:\n")
cat("    Targets annotated:", nrow(pharma_targets), "\n")
cat("    COLOC tiers:\n")
print(pharma_targets[, .N, by = coloc_tier][order(-N)])
cat("    Sex classes:\n")
print(pharma_targets[, .N, by = sex_class][order(-N)])
cat("    Stage classes:\n")
print(pharma_targets[, .N, by = stage_class][order(-N)])
cat("    Subtype markers:\n")
print(pharma_targets[, .N, by = subtype_marker][order(-N)])
cat("    ATAC evidence:\n")
print(pharma_targets[, .N, by = atac_evidence][order(-N)])

fwrite(pharma_targets,
       file.path(outdir, "pharmacogenomic_stratified_targets.csv"))
cat("\n  Saved: pharmacogenomic_stratified_targets.csv\n")

# ===========================================================================
# PART B: LINCS/CGP stratified reversal
# ===========================================================================
cat("\n===================================================\n")
cat("PART B: CGP reversal stratified by sex/stage/subtype\n")
cat("===================================================\n\n")

# Build gene sets for stratified scoring:
# (a) Female-specific DEGs (Ensembl base)
# Auto-detect v2 (interaction-based) vs v1 (stratified) class names
sex_cls_vals <- unique(sex_deg$sex_class)
use_v2_215 <- "Female_biased" %in% sex_cls_vals
fem_label_215 <- if (use_v2_215) "Female_biased" else "Female_specific"
mal_label_215 <- if (use_v2_215) "Male_biased"   else "Male_specific"
# In v2, "Concordant" replaces both "Shared" and "Not_significant"
shared_labels <- if (use_v2_215) "Concordant" else "Shared"

female_ens <- sex_deg[sex_class == fem_label_215, unique(ensembl_base)]
male_ens   <- sex_deg[sex_class == mal_label_215, unique(ensembl_base)]
shared_ens <- sex_deg[sex_class %in% shared_labels, unique(ensembl_base)]
cat("  Female-biased DEG set:", length(female_ens), "genes\n")
cat("  Male-biased DEG set:", length(male_ens), "genes\n")
cat("  Concordant DEG set:", length(shared_ens), "genes\n")

# (b) Early vs late stage genes (fibrosis axis)
early_stages <- c("F0_vs_rest", "F1_vs_rest")
late_stages  <- c("F3_vs_rest", "F4_vs_rest")
early_ens <- stage_uniq[staging_axis == "Fibrosis" &
                        most_specific_stage %in% early_stages, unique(ensembl_base)]
late_ens  <- stage_uniq[staging_axis == "Fibrosis" &
                        most_specific_stage %in% late_stages, unique(ensembl_base)]
cat("  Early-stage gene set:", length(early_ens), "genes\n")
cat("  Late-stage gene set:", length(late_ens), "genes\n")

# (c) S1 vs S2 marker genes (already loaded)
cat("  S1 marker gene set:", length(s1_markers_ens), "genes\n")
cat("  S2 marker gene set:", length(s2_markers_ens), "genes\n")

# Dream t-stat signature for weighting
tstat_lookup <- setNames(dream$t, dream$ensembl_base)

# ---- Process each CGP hit ----
lincs_results <- list()

for (j in seq_len(nrow(cgp))) {
  pathway_name <- cgp$pathway[j]
  nes <- cgp$NES[j]
  padj_val <- cgp$padj[j]
  cgp_class <- cgp$cgp_class[j]
  le_str <- cgp$leadingEdge[j]

  # Parse leading edge Ensembl IDs
  le_genes <- unique(unlist(strsplit(le_str, ";")))
  le_genes_base <- sub("\\.[0-9]+$", "", le_genes)
  n_le <- length(le_genes_base)

  if (n_le == 0) {
    next
  }

  # Reversal z-score (absolute NES; negative NES = reversed signature)
  reversal_z <- abs(nes)

  # --- Sex-stratified overlap ---
  n_female <- sum(le_genes_base %in% female_ens)
  n_male   <- sum(le_genes_base %in% male_ens)
  n_shared <- sum(le_genes_base %in% shared_ens)
  frac_female <- n_female / n_le
  frac_male   <- n_male / n_le
  frac_shared <- n_shared / n_le

  # Sex-specific reversal: enrichment of female vs male DEGs in leading edge
  # Fisher-like score: log2(frac_female / frac_male) with pseudocount
  sex_log2ratio <- log2((frac_female + 0.01) / (frac_male + 0.01))
  # Positive = female-biased reversal, negative = male-biased

  # Sex match factor: compound reverses MORE genes of a specific sex
  sex_female_match <- frac_female + frac_shared  # relevance for female patients
  sex_male_match   <- frac_male + frac_shared    # relevance for male patients

  # --- Stage-stratified overlap ---
  n_early <- sum(le_genes_base %in% early_ens)
  n_late  <- sum(le_genes_base %in% late_ens)
  frac_early <- n_early / n_le
  frac_late  <- n_late / n_le
  stage_log2ratio <- log2((frac_late + 0.01) / (frac_early + 0.01))
  # Positive = late-stage biased, negative = early-stage biased

  stage_early_match <- frac_early
  stage_late_match  <- frac_late

  # --- Subtype-stratified overlap ---
  n_s1 <- sum(le_genes_base %in% s1_markers_ens)
  n_s2 <- sum(le_genes_base %in% s2_markers_ens)
  frac_s1 <- n_s1 / n_le
  frac_s2 <- n_s2 / n_le
  subtype_log2ratio <- log2((frac_s2 + 0.01) / (frac_s1 + 0.01))
  # Positive = S2-biased, negative = S1-biased

  # --- Weighted t-stat reversal scores per stratum ---
  # For each leading edge gene, weight by stratum membership and t-stat
  le_tstats <- tstat_lookup[le_genes_base]
  le_tstats <- le_tstats[!is.na(le_tstats)]

  # S1-early score: mean |t| of leading edge genes that are S1-marker OR early-stage
  s1_early_genes <- le_genes_base[le_genes_base %in% s1_markers_ens |
                                  le_genes_base %in% early_ens]
  s1_late_genes  <- le_genes_base[le_genes_base %in% s1_markers_ens |
                                  le_genes_base %in% late_ens]
  s2_early_genes <- le_genes_base[le_genes_base %in% s2_markers_ens |
                                  le_genes_base %in% early_ens]
  s2_late_genes  <- le_genes_base[le_genes_base %in% s2_markers_ens |
                                  le_genes_base %in% late_ens]

  mean_t_or_zero <- function(genes) {
    ts <- tstat_lookup[genes]
    ts <- ts[!is.na(ts)]
    if (length(ts) > 0) mean(abs(ts)) else 0
  }

  s1_early_score <- reversal_z * (sex_male_match + stage_early_match + frac_s1) / 3
  s1_late_score  <- reversal_z * (sex_male_match + stage_late_match + frac_s1) / 3
  s2_early_score <- reversal_z * (sex_female_match + stage_early_match + frac_s2) / 3
  s2_late_score  <- reversal_z * (sex_female_match + stage_late_match + frac_s2) / 3

  # Also compute per-sex strata
  female_early_score <- reversal_z * (sex_female_match + stage_early_match) / 2
  female_late_score  <- reversal_z * (sex_female_match + stage_late_match) / 2
  male_early_score   <- reversal_z * (sex_male_match + stage_early_match) / 2
  male_late_score    <- reversal_z * (sex_male_match + stage_late_match) / 2

  lincs_results[[j]] <- data.table(
    pathway = pathway_name,
    NES = nes,
    padj = padj_val,
    cgp_class = cgp_class,
    n_le = n_le,
    reversal_z = reversal_z,
    # Sex overlap
    n_female_le = n_female,
    n_male_le = n_male,
    n_shared_le = n_shared,
    frac_female = frac_female,
    frac_male = frac_male,
    sex_log2ratio = sex_log2ratio,
    sex_bias = fifelse(sex_log2ratio > 0.5, "Female_biased",
               fifelse(sex_log2ratio < -0.5, "Male_biased", "Balanced")),
    # Stage overlap
    n_early_le = n_early,
    n_late_le = n_late,
    frac_early = frac_early,
    frac_late = frac_late,
    stage_log2ratio = stage_log2ratio,
    stage_bias = fifelse(stage_log2ratio > 0.5, "Late_biased",
                 fifelse(stage_log2ratio < -0.5, "Early_biased", "Pan_stage")),
    # Subtype overlap
    n_s1_le = n_s1,
    n_s2_le = n_s2,
    frac_s1 = frac_s1,
    frac_s2 = frac_s2,
    subtype_log2ratio = subtype_log2ratio,
    subtype_bias = fifelse(subtype_log2ratio > 0.5, "S2_biased",
                   fifelse(subtype_log2ratio < -0.5, "S1_biased", "Balanced")),
    # Per-stratum relevance scores
    s1_early_score = s1_early_score,
    s1_late_score = s1_late_score,
    s2_early_score = s2_early_score,
    s2_late_score = s2_late_score,
    female_early_score = female_early_score,
    female_late_score = female_late_score,
    male_early_score = male_early_score,
    male_late_score = male_late_score,
    # Best stratum
    best_stratum = c("S1-Early", "S1-Late", "S2-Early", "S2-Late",
                     "Female-Early", "Female-Late", "Male-Early", "Male-Late")[
      which.max(c(s1_early_score, s1_late_score, s2_early_score, s2_late_score,
                  female_early_score, female_late_score,
                  male_early_score, male_late_score))],
    best_stratum_score = max(c(s1_early_score, s1_late_score,
                               s2_early_score, s2_late_score,
                               female_early_score, female_late_score,
                               male_early_score, male_late_score))
  )
}

lincs_strat <- rbindlist(lincs_results, fill = TRUE)
cat("  LINCS stratified results:", nrow(lincs_strat), "compounds\n\n")

cat("  Sex bias distribution:\n")
print(lincs_strat[, .N, by = sex_bias][order(-N)])
cat("\n  Stage bias distribution:\n")
print(lincs_strat[, .N, by = stage_bias][order(-N)])
cat("\n  Subtype bias distribution:\n")
print(lincs_strat[, .N, by = subtype_bias][order(-N)])
cat("\n  Best stratum distribution:\n")
print(lincs_strat[, .N, by = best_stratum][order(-N)])

fwrite(lincs_strat,
       file.path(outdir, "lincs_stratified_reversal.csv"))
cat("\n  Saved: lincs_stratified_reversal.csv\n")

# ===========================================================================
# PART C: Precision medicine matrix (long format)
# ===========================================================================
cat("\n===================================================\n")
cat("PART C: Precision medicine matrix\n")
cat("===================================================\n\n")

# Create long-format: rows = drugs/CGP-sets, columns = strata, values = scores
strata_names <- c("S1-Early", "S1-Late", "S2-Early", "S2-Late",
                  "Female-Early", "Female-Late", "Male-Early", "Male-Late")
strata_cols  <- c("s1_early_score", "s1_late_score",
                  "s2_early_score", "s2_late_score",
                  "female_early_score", "female_late_score",
                  "male_early_score", "male_late_score")

pm_list <- list()
for (k in seq_len(nrow(lincs_strat))) {
  for (s in seq_along(strata_names)) {
    pm_list[[length(pm_list) + 1]] <- data.table(
      pathway = lincs_strat$pathway[k],
      cgp_class = lincs_strat$cgp_class[k],
      reversal_z = lincs_strat$reversal_z[k],
      stratum = strata_names[s],
      score = lincs_strat[[strata_cols[s]]][k]
    )
  }
}

pm_matrix <- rbindlist(pm_list)

# Also add clinical drug targets to the precision medicine matrix
# For each drug, compute stratum scores based on their target gene properties
for (i in seq_len(nrow(pharma_targets))) {
  tg_row <- pharma_targets[i]
  drug_name <- tg_row$drug

  # Stratum scores for clinical drugs based on target properties
  # Sex match: if target is female-specific -> higher female stratum score
  # Support v1 (Female_specific/Shared) and v2 (Female_biased/Concordant)
  sex_f_weight <- fifelse(tg_row$sex_class %in% c("Female_specific", "Female_biased",
                                                    "Shared", "Concordant"),
                          1.0, 0.3)
  sex_m_weight <- fifelse(tg_row$sex_class %in% c("Male_specific", "Male_biased",
                                                    "Shared", "Concordant"),
                          1.0, 0.3)

  # Stage match
  stage_early_weight <- fifelse(tg_row$stage_class %in% c("early", "pan"), 1.0, 0.3)
  stage_late_weight  <- fifelse(tg_row$stage_class %in% c("late", "pan"), 1.0, 0.3)

  # Subtype match
  sub_s1_weight <- fifelse(tg_row$subtype_marker %in% c("S1", "Neither"), 0.5, 0.3)
  sub_s2_weight <- fifelse(tg_row$subtype_marker %in% c("S2", "Neither"), 0.5, 0.3)

  # COLOC amplifier
  # FIX (review B5#3): test >0.9 before >0.8 so the 3-tier weight is reachable. Previously
  # any pp4>0.9 matched the >0.8 branch first (weight 1.5), making the 1.2 tier dead code and
  # over-upweighting moderate (0.8-0.9) genes. Tiers: >0.9->1.5, 0.8-0.9->1.2, else 0.8.
  coloc_weight <- fifelse(is.na(tg_row$coloc_pp4), 0.5,
                  fifelse(tg_row$coloc_pp4 > 0.9, 1.5,
                  fifelse(tg_row$coloc_pp4 > 0.8, 1.2, 0.8)))

  base_score <- coloc_weight

  drug_strata <- data.table(
    pathway = paste0("DRUG:", drug_name),
    cgp_class = "Clinical_Drug",
    reversal_z = NA_real_,
    stratum = strata_names,
    score = c(
      base_score * sub_s1_weight * stage_early_weight * sex_m_weight,   # S1-Early
      base_score * sub_s1_weight * stage_late_weight * sex_m_weight,    # S1-Late
      base_score * sub_s2_weight * stage_early_weight * sex_f_weight,   # S2-Early
      base_score * sub_s2_weight * stage_late_weight * sex_f_weight,    # S2-Late
      base_score * stage_early_weight * sex_f_weight,                    # Female-Early
      base_score * stage_late_weight * sex_f_weight,                     # Female-Late
      base_score * stage_early_weight * sex_m_weight,                    # Male-Early
      base_score * stage_late_weight * sex_m_weight                      # Male-Late
    )
  )
  pm_matrix <- rbind(pm_matrix, drug_strata)
}

cat("  Precision medicine matrix:", nrow(pm_matrix), "rows\n")
cat("  Unique compounds/drugs:", length(unique(pm_matrix$pathway)), "\n")
cat("  Strata:", length(unique(pm_matrix$stratum)), "\n")

fwrite(pm_matrix,
       file.path(outdir, "precision_medicine_matrix.csv"))
cat("  Saved: precision_medicine_matrix.csv\n")

# ===========================================================================
# Summary statistics
# ===========================================================================
cat("\n===================================================\n")
cat("SUMMARY\n")
cat("===================================================\n\n")

cat("Pharmacogenomic targets:\n")
cat("  Total drug-target pairs:", nrow(pharma_targets), "\n")
cat("  With Strong COLOC (PP.H4 > 0.9):",
    sum(pharma_targets$coloc_tier == "Strong", na.rm = TRUE), "\n")
cat("  With Moderate COLOC (PP.H4 > 0.8):",
    sum(pharma_targets$coloc_tier == "Moderate", na.rm = TRUE), "\n")
cat("  With motif disruption:",
    sum(pharma_targets$n_motif_disruptions > 0), "\n")
cat("  Sex-specific targets:",
    sum(pharma_targets$sex_class %in% c("Female_specific", "Male_specific",
                                        "Female_biased", "Male_biased")), "\n")
cat("  Stage-unique targets:",
    sum(pharma_targets$stage_unique == TRUE, na.rm = TRUE), "\n")

cat("\nCGP reversal stratification:\n")
cat("  Total CGP hits stratified:", nrow(lincs_strat), "\n")
cat("  Female-biased:", sum(lincs_strat$sex_bias == "Female_biased"), "\n")
cat("  Male-biased:", sum(lincs_strat$sex_bias == "Male_biased"), "\n")
cat("  Late-stage biased:", sum(lincs_strat$stage_bias == "Late_biased"), "\n")
cat("  Early-stage biased:", sum(lincs_strat$stage_bias == "Early_biased"), "\n")
cat("  S2-biased:", sum(lincs_strat$subtype_bias == "S2_biased"), "\n")
cat("  S1-biased:", sum(lincs_strat$subtype_bias == "S1_biased"), "\n")

# Top compounds by best_stratum_score
cat("\nTop 10 CGP hits by best stratum score:\n")
top10 <- lincs_strat[order(-best_stratum_score)][1:min(10, nrow(lincs_strat)),
  .(pathway, reversal_z, sex_bias, stage_bias, subtype_bias,
    best_stratum, best_stratum_score)]
print(top10, nrows = 10)

# Top liver/MASLD hits
liver_hits <- lincs_strat[cgp_class == "Liver/MASLD"][order(-best_stratum_score)]
cat("\nTop 10 Liver/MASLD CGP hits:\n")
print(liver_hits[1:min(10, nrow(liver_hits)),
  .(pathway, reversal_z, sex_bias, stage_bias, best_stratum,
    best_stratum_score)], nrows = 10)

cat("\n\nCompleted:", format(Sys.time()), "\n")
cat("=== 215_pharmacogenomic_targets.R DONE ===\n")
