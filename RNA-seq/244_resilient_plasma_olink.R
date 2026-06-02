#!/usr/bin/env Rscript
# 244_resilient_plasma_olink.R
# LEGACY SCRIPT — "F2-switch" framing retired 2026-05-09.
# Current thesis: multi-step cellular cascade across 4 CRN transitions (see paper_outline.md).
# Script name and output filenames (resilient_olink_f2switch*) retained for provenance.
#
# Healthy-control characterization pipeline, step 5 of 6.
#
# Plasma F-stage-boundary protein score per Olink subject. Tests whether the bulk
# stage-transition protein signature is detectable in plasma proteomics across
# Yang et al. 2025 cohort (Cell Reports Medicine; 1,461 proteins × 218 subjects).
#
# Subject-level disease metadata (40 healthy / 72 MASLD / 58 CVH / 14 ARLD /
# 34 HCC / 112 validation) is NOT yet harmonized in the project (per README
# Phase 2B/2C in 81/82 are pending). This script computes the F2-switch Z
# score per Olink subject and writes both raw scores and an ordering by
# signature strength. Mapping to clinical groups requires Yang et al.
# supplementary metadata to be downloaded separately.
#
# Outputs:
#   resilient_olink_f2switch.csv  – per-subject F2-switch Z-score
#   resilient_olink_f2switch_genes.csv – genes used for scoring (for reproducibility)
#
# Spec: docs/superpowers/specs/2026-04-27-healthy-control-audit-design.md
# Env: rnaseq

suppressPackageStartupMessages({
  library(data.table)
})

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HCDIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity/healthy_control_audit")
SWITCH <- file.path(BASE, "RNA-seq/results/stratified_causal/switch_gene_classification.csv")
OLINK  <- file.path(BASE, "Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt")
GEO_META <- file.path(BASE, "Analysis/Proteomics/results/gse276114_disease_metadata.csv")

dir.create(HCDIR, recursive = TRUE, showWarnings = FALSE)

cat("Loading switch_gene_classification.csv...\n")
sw <- fread(SWITCH)
cat(sprintf("  Total switch genes: %d\n", nrow(sw)))

# Map ensembl -> symbol via atlas
ATLAS <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
atlas <- fread(ATLAS, select = c("ensembl_id", "human_symbol"))
atlas[, ensembl_clean := sub("\\..*", "", ensembl_id)]
sw[, ensembl_clean := sub("\\..*", "", gene)]
sw <- merge(sw, unique(atlas[, .(ensembl_clean, symbol = human_symbol)]),
             by = "ensembl_clean", all.x = TRUE)
cat(sprintf("  Switch genes mapped to symbol: %d\n", sum(!is.na(sw$symbol))))

# F2-switch protein set: classification %in% switch_up/switch_down + |logFC| > 0.3
f2 <- sw[classification %in% c("switch_up", "switch_down") & abs(logFC) > 0.3 &
          !is.na(symbol) & symbol != "",
          .(gene = symbol, classification, logFC, padj, t, direction)]
cat(sprintf("  F2-switch protein-set candidates (|logFC|>0.3): %d\n", nrow(f2)))
print(table(f2$classification))

cat("\nLoading Olink data...\n")
olink <- fread(OLINK)
olink_genes <- olink$Assay
n_subjects <- ncol(olink) - 1
cat(sprintf("  Olink: %d proteins × %d subjects\n", length(olink_genes), n_subjects))

# Match F2-switch gene set to Olink panel
overlap <- intersect(f2$gene, olink_genes)
cat(sprintf("\nF2-switch ∩ Olink: %d genes\n", length(overlap)))
if (length(overlap) < 5) {
  cat("  [warn] Too few overlap genes; F2-switch score will be unreliable\n")
}

f2_olink <- f2[gene %in% overlap]
cat("F2-switch direction breakdown in overlap:\n")
print(table(f2_olink$direction))

# Subset Olink matrix to F2-switch genes
ol_mat <- as.matrix(olink[Assay %in% overlap, -"Assay"])
rownames(ol_mat) <- olink$Assay[olink$Assay %in% overlap]

# Z-score per protein across subjects (since NPX is already log2, we standardize)
ol_z <- t(scale(t(ol_mat)))  # per-protein Z

# Sign by F2-switch direction
sign_vec <- f2_olink[match(rownames(ol_z), gene), fifelse(direction == "up", 1, -1)]
ol_z_signed <- ol_z * sign_vec  # +Z = direction-consistent (= "more F2-switch-like")

# Per-subject F2-switch score = mean signed Z across F2-switch proteins
subject_scores <- data.table(
  subject = colnames(ol_z),
  f2_switch_z = colMeans(ol_z_signed, na.rm = TRUE),
  n_proteins_used = colSums(!is.na(ol_z_signed))
)
# Disease label mapping:
# Yang et al. 2025 discovery cohort = 218 subjects.
# GEO GSE276114 has 177 liver-biopsy RNA-seq samples ("Liver sample 1"..."Liver sample 177").
# Hypothesis: Olink Subject 1-177 = GEO liver samples (have biopsies, are CLD/HCC).
#             Subject 178-218 (n=41) = healthy plasma-only controls (no liver biopsy → absent from GEO).
# This matches the known 40 healthy + 177 CLD/HCC composition.
cat("\nMerging disease labels from GEO metadata (sequential numbering hypothesis)...\n")
geo <- fread(GEO_META)
geo[, subject := paste0("Subject ", sample_number)]
geo[, disease_broad := fcase(
  disease_group %in% c("F0-2"),  "CLD_early",
  disease_group %in% c("F3"),    "CLD_F3",
  disease_group %in% c("F4"),    "CLD_F4",
  default = NA_character_
)]
subject_scores[, subject_num := as.integer(sub("Subject ", "", subject))]
subject_scores <- merge(subject_scores,
                        geo[, .(subject, disease_group, disease, disease_broad)],
                        by = "subject", all.x = TRUE)
# Assign healthy to subjects beyond GEO range (178-218)
subject_scores[is.na(disease_group) & subject_num >= 178, disease_group := "Healthy"]
subject_scores[is.na(disease_broad)  & subject_num >= 178, disease_broad  := "Healthy"]

cat(sprintf("  Disease label coverage: %d / %d subjects labelled\n",
            sum(!is.na(subject_scores$disease_group)), nrow(subject_scores)))
print(table(subject_scores$disease_broad, useNA = "ifany"))

subject_scores <- subject_scores[order(-f2_switch_z)]
fwrite(subject_scores, file.path(HCDIR, "resilient_olink_f2switch.csv"))
cat(sprintf("\nWrote: resilient_olink_f2switch.csv (%d subjects)\n", nrow(subject_scores)))

# Score distribution by disease group
cat("\nScore distribution by disease group:\n")
subject_scores[, .(
  n = .N,
  mean_z = mean(f2_switch_z, na.rm = TRUE),
  median_z = median(f2_switch_z, na.rm = TRUE)
), by = disease_broad][order(mean_z)] |> print()

# Write the scoring gene set for reproducibility
fwrite(f2_olink, file.path(HCDIR, "resilient_olink_f2switch_genes.csv"))
cat(sprintf("\nWrote: resilient_olink_f2switch_genes.csv (%d genes used)\n", nrow(f2_olink)))

cat("\nDone.\n")
