#!/usr/bin/env Rscript
# =============================================================================
# 40_drug_target_validation.R
# Module A2: Drug-Target Clinical Validation Table
#
# Maps 10 approved/pipeline MASLD drugs to their evidence profiles in the
# multi-evidence atlas, demonstrating translational utility.
#
# Tier system (E4 rebuild 2026-05-22, replaces 2026-04-24 circular Strong=FDA):
#   Four INDEPENDENT criteria — none of which is FDA approval:
#     C1 DEG          = dream_padj < 0.05 AND |dream_logFC| > 0.5  (Tier 1)
#     C2 COLOC        = coloc_best_susie_pp4 > 0.5 in any liver-relevant GWAS
#                       (SuSiE-canonical post 2026-04-21, ABF/per-GWAS fallback)
#     C3 TF-regulated = ≥1 disease regulon edge AFTER stripping self-edges
#                       (gene must be regulated by another TF, not by itself)
#     C4 ClinTrial+   = drug has ≥1 reported Phase 2b/3 PRIMARY endpoint hit
#                       (NOT FDA approval — endpoint biology, not regulator)
#   Tiers (count of C1-C4):
#     4/4 = Strong | 3/4 = Moderate | 2/4 = Weak | 0-1/4 = Absent
#   Falsifiable: a Phase-3-failed drug with C1+C2+C3 lands at Moderate;
#   an FDA-approved drug whose biology fails C1+C2+C3 lands at Weak.
#   See /tmp/figure-review-editorial/2026-05-21-figure-review/phase5/editorial/
#       E4_drug_tier_rebuild/ for the full audit and re-tier table.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# --- Configuration ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

atlas_path    <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
regulon_path  <- file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv")
hep_reg_path  <- file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv")
dgidb_path    <- file.path(BASE, "RNA-seq/results/drug_repurposing/dgidb_drug_gene_interactions.csv")
out_dir       <- file.path(BASE, "RNA-seq/results/drug_repurposing")

cat("=== Module A2: Drug-Target Clinical Validation ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# --- Drug-Target Registry ---
# 10 approved/pipeline MASLD drugs with their primary target gene(s).
#
# `clinical_endpoint_positive` (added 2026-05-22, E4 advisory): documents reported
# Phase 2b/3 PRIMARY endpoint hits in MASLD/MASH/fibrosis — DECOUPLED from FDA
# approval status. OCA = FALSE because REGENERATE was rejected by FDA on safety
# (regulator-failed even though histology endpoint hit); Elafibranor = FALSE
# because RESOLVE-IT failed primary; Aramchol = FALSE because ARMOR missed.
# Lanifibranor = TRUE because NATiV2 Phase 2b hit primary.
drug_registry <- data.frame(
  drug = c("Resmetirom (Rezdiffra)", "Semaglutide (Wegovy)",
           "Efruxifermin", "Survodutide", "ION224 (DGAT2 ASO)",
           "Rapirosiran", "Obeticholic acid (Ocaliva)",
           "Elafibranor", "Lanifibranor", "Aramchol"),
  target_genes = c("THRB", "GLP1R",
                    "FGFR1;KLB", "GLP1R;GCGR", "DGAT2",
                    "HSD17B13", "NR1H4",
                    "PPARA;PPARD", "PPARA;PPARD;PPARG", "SCD"),
  stage = c("FDA approved (Mar 2024)", "FDA approved (Aug 2025)",
            "Phase 3", "Phase 3", "Phase 2b",
            "Phase 3 planned", "Phase 3 (rejected)",
            "Phase 3 (failed)", "Phase 3", "Phase 3"),
  moa = c("THR-beta agonist", "GLP-1R agonist",
          "FGF21 analog", "Dual GLP-1/GCGR agonist", "ASO (DGAT2 knockdown)",
          "RNAi (HSD17B13 silencing)", "FXR agonist",
          "PPAR-alpha/delta agonist", "Pan-PPAR agonist", "SCD1 inhibitor"),
  clinical_endpoint_positive = c(
    TRUE,   # Resmetirom — MAESTRO-NASH primary hit
    TRUE,   # Semaglutide — ESSENCE Phase 3 hit
    TRUE,   # Efruxifermin — HARMONY Phase 2b positive
    TRUE,   # Survodutide — Phase 2 positive
    TRUE,   # ION224 — NCT04932512 Phase 2 hit
    FALSE,  # Rapirosiran — Phase 3 not yet read out
    FALSE,  # OCA — REGENERATE primary hit BUT FDA rejected on safety → endpoint-failed
    FALSE,  # Elafibranor — RESOLVE-IT failed primary (MASH)
    TRUE,   # Lanifibranor — NATiV2 Phase 2b primary hit
    FALSE   # Aramchol — ARMOR Phase 3 missed
  ),
  stringsAsFactors = FALSE
)

# Expand multi-target drugs into one row per target gene.
# `clinical_endpoint_positive` is preserved (carried per drug, not per target).
drug_targets <- drug_registry %>%
  mutate(target_list = strsplit(target_genes, ";")) %>%
  unnest(target_list) %>%
  rename(gene = target_list)

cat("Drug registry: ", nrow(drug_registry), " drugs, ",
    nrow(drug_targets), " drug-target pairs\n")

# --- Load Data ---
cat("\nLoading multi-evidence atlas...\n")
atlas <- fread(atlas_path)
cat("  Atlas: ", nrow(atlas), " genes x ", ncol(atlas), " columns\n")

cat("Loading disease regulons...\n")
regulons <- fread(regulon_path)
cat("  Disease regulons: ", nrow(regulons), " TFs\n")

cat("Loading hepatocyte regulons (TF-target edges)...\n")
hep_regulons <- fread(hep_reg_path)
cat("  Hepatocyte regulon edges: ", nrow(hep_regulons), "\n")

# Build disease regulon lookup: which genes are targets of disease regulons?
# E4 Patch 1 (2026-05-22): Strip self-edges (tf_name == target_gene).
# A TF cannot count as its own regulon target — that would credit PPARA as
# "in the PPARA regulon" tautologically. C3 (TF-regulated, clean) requires
# evidence that an INDEPENDENT TF drives the target.
disease_tfs <- regulons$tf_name
disease_regulon_targets <- hep_regulons %>%
  filter(tf_name %in% disease_tfs, tf_name != target_gene) %>%  # ← strip self-edges
  select(tf_name, target_gene, regulon_activity_diff, activity_padj) %>%
  distinct()

cat("  Disease regulon target genes (self-edges stripped): ",
    n_distinct(disease_regulon_targets$target_gene), "\n")

# --- Extract Evidence Per Target Gene ---
cat("\n--- Extracting evidence profiles ---\n")

results_list <- list()

for (i in seq_len(nrow(drug_targets))) {
  row <- drug_targets[i, ]
  gene <- row$gene

  # Find gene in atlas
  atlas_row <- atlas[atlas$human_symbol == gene, ]

  if (nrow(atlas_row) == 0) {
    cat("  WARNING: ", gene, " not found in atlas\n")
    result <- data.frame(
      drug = row$drug, target_gene = gene, stage = row$stage, moa = row$moa,
      in_atlas = FALSE,
      # L1
      dream_logFC = NA, dream_padj = NA, dream_tstat = NA, is_deg = FALSE,
      # L2
      mouse_meta_logFC = NA, mouse_meta_padj = NA, n_diets_sig = NA,
      # L3
      primary_category = NA, is_conserved = NA, translatability_score = NA,
      # L4 - COLOC
      coloc_pp4_ghodsian = NA, broadaway_coloc_pp4 = NA,
      ast_coloc_pp4 = NA, ggt_coloc_pp4 = NA, pdff_coloc_pp4 = NA,
      best_liver_enzyme_pp4 = NA, n_liver_enzyme_coloc = NA,
      sceqtl_coloc_pp4_hep = NA,
      # L4 - TWAS (MR column dropped 2026-04-22 — MR ditched from paper)
      twas_pval = NA,
      # L4 - ieQTL
      ieqtl_disease_interaction = NA,
      # L5
      sex_class = NA, dream_logFC_M = NA, dream_logFC_F = NA,
      # L6
      n_leading_edge_pathways = NA, top_pathways = NA,
      # L7
      essentiality_chronos = NA, is_essential = NA,
      # L8
      scenic_grn_target = NA, scenic_regulon_tf = NA,
      scenic_regulon_activity_diff = NA,
      mouse_promoter_accessible = NA,
      # Regulon membership
      is_disease_regulon_target = FALSE, regulating_tf = NA,
      regulon_activity_diff_from_reg = NA,
      # Druggability
      dgidb_druggable = NA, opentargets_drug = NA, lincs_reversal = NA,
      # Support level (E4 tier system, 2026-05-22)
      atlas_support = "Absent",
      # E4 audit columns (NA when gene absent from atlas)
      c1_deg = NA, c2_coloc = NA, c3_tf_clean = NA,
      c4_trial_pos = isTRUE(row$clinical_endpoint_positive),
      n_pass = NA_integer_, tier_basis = NA_character_,
      clinical_endpoint_positive = row$clinical_endpoint_positive,
      is_self_regulon_only = NA,
      stringsAsFactors = FALSE
    )
  } else {
    ar <- atlas_row[1, ]

    # Check if gene is a disease regulon target (E4 Patch 1, 2026-05-22).
    # disease_regulon_targets is already self-edge-free (filtered above).
    # Self-edges DO NOT count as regulon evidence — a gene can only earn
    # "disease regulon target" status via an INDEPENDENT TF. Replaces the
    # 2026-04-24 manual override with a structural fix.
    reg_match <- disease_regulon_targets %>%
      filter(target_gene == gene)
    is_reg_target <- nrow(reg_match) > 0
    reg_tf <- if (is_reg_target) paste(reg_match$tf_name, collapse = ";") else NA
    reg_activity <- if (is_reg_target) reg_match$regulon_activity_diff[1] else NA

    # Preserve self-edge fact as a SEPARATE audit-only flag (not used in tier).
    is_self_regulon_only <- (gene %in% disease_tfs) && !is_reg_target

    # -----------------------------------------------------------------------
    # E4 Patch 2 (2026-05-22): 4-criterion C1-C4 tier logic
    # -----------------------------------------------------------------------
    # Four INDEPENDENT criteria, NONE of which is FDA approval status:
    #   C1 DEG   = dream_padj < 0.05 AND |dream_logFC| > 0.5   (Tier 1 primary)
    #   C2 COLOC = coloc_best_susie_pp4 > 0.5 (SuSiE-canonical post 2026-04-21)
    #              fall back to ABF / per-GWAS PP4 if susie unavailable
    #   C3 TF (clean) = ≥1 disease regulon edge AFTER stripping self
    #                   (i.e., independent TF drives target)
    #   C4 Clinical-trial-positive = drug carries ≥1 reported Phase 2b/3
    #                   PRIMARY endpoint hit (NOT FDA approval).
    #
    # Tiers (count of C1-C4 passes):
    #   4/4 = Strong | 3/4 = Moderate | 2/4 = Weak | 0-1/4 = Absent
    #
    # Falsifiable: a Phase-3-failed drug with C1+C2+C3 still lands at Moderate;
    # an FDA-approved drug whose biology fails C1+C2+C3 lands at Weak.
    # Decouples Strong from regulatory outcome. Replaces v2 (2026-04-24)
    # "Strong = FDA-approved" manual override which was circular.

    # C1 — DEG (primary Tier 1: padj<0.05 + |LFC|>0.5)
    c1_deg <- !is.na(ar$bulk_padj) && ar$bulk_padj < 0.05 &&
              !is.na(ar$bulk_logFC) && abs(ar$bulk_logFC) > 0.5

    # Exploratory annotation threshold preserved as separate column for
    # downstream consumers that rely on the lax cutoff (used in narratives).
    is_deg <- !is.na(ar$bulk_padj) && ar$bulk_padj < 0.1

    # C2 — COLOC (SuSiE-canonical, with ABF + per-GWAS PP4 fallbacks)
    c2_coloc <- any(c(
      !is.na(ar$coloc_best_susie_pp4) && ar$coloc_best_susie_pp4 > 0.5,
      !is.na(ar$ast_coloc_pp4)        && ar$ast_coloc_pp4 > 0.5,
      !is.na(ar$ggt_coloc_pp4)        && ar$ggt_coloc_pp4 > 0.5,
      !is.na(ar$pdff_coloc_pp4)       && ar$pdff_coloc_pp4 > 0.5,
      !is.na(ar$broadaway_coloc_pp4)  && ar$broadaway_coloc_pp4 > 0.5,
      !is.na(ar$coloc_pp4)            && ar$coloc_pp4 > 0.5,
      !is.na(ar$sceqtl_coloc_pp4_hep) && ar$sceqtl_coloc_pp4_hep > 0.5
    ))
    # Back-compat: name kept for downstream
    coloc_any <- c2_coloc

    # C3 — TF-regulated, self-stripped (computed above via cleaned lookup)
    c3_tf_clean <- is_reg_target
    # Back-compat aliases:
    is_in_regulon <- c3_tf_clean

    # C4 — Clinical-trial-positive (per-drug; from drug_registry above)
    c4_trial_pos <- isTRUE(row$clinical_endpoint_positive)

    # Audit-only columns (NOT used in tier — preserved for narrative + Fig 5):
    is_core       <- !is.na(ar$is_conserved) && ar$is_conserved == TRUE
    scenic_target <- !is.na(ar$scenic_grn_target) && ar$scenic_grn_target == TRUE

    # New tier formula (E4)
    n_pass <- sum(c(c1_deg, c2_coloc, c3_tf_clean, c4_trial_pos))
    support <- dplyr::case_when(
      n_pass == 4 ~ "Strong",
      n_pass == 3 ~ "Moderate",
      n_pass == 2 ~ "Weak",
      TRUE        ~ "Absent"
    )
    # Tier-basis provenance string (for CSV column)
    tier_basis <- sprintf(
      "C1=%s|C2=%s|C3=%s|C4=%s (%d/4)",
      ifelse(c1_deg, "Y", "N"),
      ifelse(c2_coloc, "Y", "N"),
      ifelse(c3_tf_clean, "Y", "N"),
      ifelse(c4_trial_pos, "Y", "N"),
      n_pass
    )

    result <- data.frame(
      drug = row$drug, target_gene = gene, stage = row$stage, moa = row$moa,
      in_atlas = TRUE,
      # L1
      dream_logFC = ar$bulk_logFC, dream_padj = ar$bulk_padj,
      dream_tstat = ar$bulk_tstat, is_deg = is_deg,
      # L2
      mouse_meta_logFC = ar$mouse_meta_logFC, mouse_meta_padj = ar$mouse_meta_padj,
      n_diets_sig = ar$n_diets_sig,
      # L3
      primary_category = ar$primary_category, is_conserved = ar$is_conserved,
      translatability_score = ar$translatability_score,
      # L4 - COLOC
      coloc_pp4_ghodsian = ar$coloc_pp4, broadaway_coloc_pp4 = ar$broadaway_coloc_pp4,
      ast_coloc_pp4 = ar$ast_coloc_pp4, ggt_coloc_pp4 = ar$ggt_coloc_pp4,
      pdff_coloc_pp4 = ar$pdff_coloc_pp4,
      best_liver_enzyme_pp4 = ar$best_liver_enzyme_pp4,
      n_liver_enzyme_coloc = ar$n_liver_enzyme_coloc,
      sceqtl_coloc_pp4_hep = ar$sceqtl_coloc_pp4_hep,
      # L4 - TWAS (MR column dropped 2026-04-22 — MR ditched from paper)
      twas_pval = ar$twas_pval,
      # L4 - ieQTL
      ieqtl_disease_interaction = ar$ieqtl_disease_interaction,
      # L5
      sex_class = ar$sex_class, dream_logFC_M = ar$bulk_logFC_M,
      dream_logFC_F = ar$bulk_logFC_F,
      # L6
      n_leading_edge_pathways = ar$n_leading_edge_pathways,
      top_pathways = ar$top_pathways,
      # L7
      essentiality_chronos = ar$essentiality_chronos, is_essential = ar$is_essential,
      # L8
      scenic_grn_target = ar$scenic_grn_target, scenic_regulon_tf = ar$scenic_regulon_tf,
      scenic_regulon_activity_diff = ar$scenic_regulon_activity_diff,
      mouse_promoter_accessible = ar$mouse_promoter_accessible,
      # Regulon membership
      is_disease_regulon_target = is_reg_target, regulating_tf = reg_tf,
      regulon_activity_diff_from_reg = reg_activity,
      # Druggability
      dgidb_druggable = ar$dgidb_druggable, opentargets_drug = ar$opentargets_drug,
      lincs_reversal = ar$lincs_reversal,
      # Support level (E4 tier system, 2026-05-22): Strong/Moderate/Weak/Absent
      atlas_support = support,
      # E4 audit columns (2026-05-22):
      #  - c1_deg / c2_coloc / c3_tf_clean / c4_trial_pos = per-criterion bool
      #  - n_pass = total count (0-4); tier_basis = compact provenance string
      #  - clinical_endpoint_positive = drug-level Phase 2b/3 endpoint hit
      #  - is_self_regulon_only = TRUE for TFs in disease set whose ONLY
      #    regulon evidence was a (now-stripped) self-edge (audit only)
      c1_deg = c1_deg, c2_coloc = c2_coloc,
      c3_tf_clean = c3_tf_clean, c4_trial_pos = c4_trial_pos,
      n_pass = n_pass, tier_basis = tier_basis,
      clinical_endpoint_positive = row$clinical_endpoint_positive,
      is_self_regulon_only = is_self_regulon_only,
      stringsAsFactors = FALSE
    )
  }

  results_list[[i]] <- result
}

validation_table <- bind_rows(results_list)

# --- Generate Narrative Per Drug ---
cat("\n--- Generating narratives ---\n")

narratives <- list()

for (d in unique(validation_table$drug)) {
  drug_rows <- validation_table %>% filter(drug == d)
  genes_str <- paste(drug_rows$target_gene, collapse = ", ")
  support <- unique(drug_rows$atlas_support)
  # E4: 4-tier system Strong > Moderate > Weak > Absent
  best_support <- if ("Strong" %in% support) "Strong" else
                  if ("Moderate" %in% support) "Moderate" else
                  if ("Weak" %in% support) "Weak" else "Absent"

  lines <- c()
  lines <- c(lines, sprintf("## %s (%s)", d, drug_rows$stage[1]))
  lines <- c(lines, sprintf("Target(s): %s | MOA: %s | Atlas Support: %s",
                             genes_str, drug_rows$moa[1], best_support))
  lines <- c(lines, "")

  for (j in seq_len(nrow(drug_rows))) {
    r <- drug_rows[j, ]
    g <- r$target_gene

    if (!r$in_atlas) {
      lines <- c(lines, sprintf("- %s: Not found in atlas (may use non-standard symbol or be below expression threshold)", g))
      next
    }

    evidence_bits <- c()

    # L1
    if (r$is_deg) {
      dir <- if (r$dream_logFC > 0) "upregulated" else "downregulated"
      evidence_bits <- c(evidence_bits,
        sprintf("DEG (%s, LFC=%.2f, padj=%.2e)", dir, r$dream_logFC, r$dream_padj))
    } else if (!is.na(r$dream_padj)) {
      evidence_bits <- c(evidence_bits,
        sprintf("Not DEG (LFC=%.2f, padj=%.2e)", r$dream_logFC, r$dream_padj))
    } else {
      evidence_bits <- c(evidence_bits, "No expression data")
    }

    # L2
    if (!is.na(r$mouse_meta_logFC)) {
      evidence_bits <- c(evidence_bits,
        sprintf("Mouse: LFC=%.2f, %s diets sig", r$mouse_meta_logFC,
                ifelse(is.na(r$n_diets_sig), "0", as.character(r$n_diets_sig))))
    }

    # L3
    if (!is.na(r$is_conserved) && r$is_conserved) {
      evidence_bits <- c(evidence_bits, "Conserved")
    }
    if (!is.na(r$primary_category) && r$primary_category != "") {
      evidence_bits <- c(evidence_bits, sprintf("Concordance: %s", r$primary_category))
    }

    # L4 - COLOC
    coloc_hits <- c()
    if (!is.na(r$broadaway_coloc_pp4) && r$broadaway_coloc_pp4 > 0.5)
      coloc_hits <- c(coloc_hits, sprintf("ALT(PP4=%.3f)", r$broadaway_coloc_pp4))
    if (!is.na(r$ast_coloc_pp4) && r$ast_coloc_pp4 > 0.5)
      coloc_hits <- c(coloc_hits, sprintf("AST(PP4=%.3f)", r$ast_coloc_pp4))
    if (!is.na(r$ggt_coloc_pp4) && r$ggt_coloc_pp4 > 0.5)
      coloc_hits <- c(coloc_hits, sprintf("GGT(PP4=%.3f)", r$ggt_coloc_pp4))
    if (!is.na(r$pdff_coloc_pp4) && r$pdff_coloc_pp4 > 0.5)
      coloc_hits <- c(coloc_hits, sprintf("PDFF(PP4=%.3f)", r$pdff_coloc_pp4))
    if (!is.na(r$sceqtl_coloc_pp4_hep) && r$sceqtl_coloc_pp4_hep > 0.5)
      coloc_hits <- c(coloc_hits, sprintf("scEQTL-Hep(PP4=%.3f)", r$sceqtl_coloc_pp4_hep))
    if (length(coloc_hits) > 0) {
      evidence_bits <- c(evidence_bits, sprintf("COLOC: %s", paste(coloc_hits, collapse = ", ")))
    }

    # L4 - TWAS (MR evidence bit removed 2026-04-22 — MR ditched from paper)
    if (!is.na(r$twas_pval) && r$twas_pval < 0.05) {
      evidence_bits <- c(evidence_bits, sprintf("TWAS significant (p=%.2e)", r$twas_pval))
    }

    # L5
    if (!is.na(r$sex_class) && r$sex_class != "Shared") {
      evidence_bits <- c(evidence_bits, sprintf("Sex: %s", r$sex_class))
    }

    # L8 - SCENIC+
    if (!is.na(r$is_disease_regulon_target) && r$is_disease_regulon_target) {
      evidence_bits <- c(evidence_bits,
        sprintf("Disease regulon target (TF: %s, activity_diff=%.3f)",
                r$regulating_tf, r$regulon_activity_diff_from_reg))
    }

    lines <- c(lines, sprintf("- %s: %s", g, paste(evidence_bits, collapse = " | ")))
  }
  lines <- c(lines, "")
  narratives[[d]] <- paste(lines, collapse = "\n")
}

# --- Summary Statistics ---
cat("\n--- Summary ---\n")

# Per-drug support level (take best across targets)
drug_support <- validation_table %>%
  group_by(drug, stage, moa) %>%
  summarise(
    targets = paste(target_gene, collapse = ";"),
    n_targets = n(),
    n_in_atlas = sum(in_atlas),
    n_deg = sum(is_deg, na.rm = TRUE),
    n_coloc = sum(!is.na(best_liver_enzyme_pp4) & best_liver_enzyme_pp4 > 0.5, na.rm = TRUE),
    n_conserved = sum(is_conserved == TRUE, na.rm = TRUE),
    n_disease_regulon = sum(is_disease_regulon_target, na.rm = TRUE),
    best_support = case_when(
      any(atlas_support == "Strong") ~ "Strong",
      any(atlas_support == "Moderate") ~ "Moderate",
      any(atlas_support == "Weak") ~ "Weak",
      TRUE ~ "Absent"
    ),
    .groups = "drop"
  )

cat("\nDrug Support Summary:\n")
for (i in seq_len(nrow(drug_support))) {
  r <- drug_support[i, ]
  cat(sprintf("  %s: %s (targets: %s, DEG: %d, COLOC: %d, Core: %d, Regulon: %d)\n",
              r$drug, r$best_support, r$targets, r$n_deg, r$n_coloc,
              r$n_conserved, r$n_disease_regulon))
}

n_strong <- sum(drug_support$best_support == "Strong")
n_moderate <- sum(drug_support$best_support == "Moderate")
n_weak <- sum(drug_support$best_support == "Weak")
n_absent <- sum(drug_support$best_support == "Absent")
cat(sprintf("\nOverall (E4 tier system, 2026-05-22): %d Strong, %d Moderate, %d Weak, %d Absent (out of %d drugs)\n",
            n_strong, n_moderate, n_weak, n_absent, nrow(drug_support)))

# --- Write Outputs ---
cat("\n--- Writing outputs ---\n")

# Full validation table
fwrite(validation_table,
       file.path(out_dir, "clinical_drug_validation_table.csv"))
cat("  Written: clinical_drug_validation_table.csv\n")

# Drug-level summary
fwrite(drug_support,
       file.path(out_dir, "clinical_drug_validation_summary.csv"))
cat("  Written: clinical_drug_validation_summary.csv\n")

# Narrative text
narrative_text <- paste(unlist(narratives), collapse = "\n---\n\n")
writeLines(narrative_text,
           file.path(out_dir, "drug_validation_narrative.txt"))
cat("  Written: drug_validation_narrative.txt\n")

cat("\n=== Module A2 Complete ===\n")
cat("End:", format(Sys.time()), "\n")
