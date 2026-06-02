#!/usr/bin/env Rscript
# 241_subclinical_screen.R
#
# Healthy-control characterization pipeline, step 2 of 6.
#
# For control-pool subjects (NAS<3 + fibrosis<=F1), compute a molecular
# subclinical score and stratify into clean / intermediate / suspect tertiles.
#
# Primary axis: NMF P_Pro_inflammatory_score (Spearman 0.472 with NAS in disease,
# validated in pre-flight). Secondary axis: pseudo-NAS-component score from
# nas_components_all.csv DE results projected onto each subject by the bulk
# atlas dream_logFC anchors (only used if NMF coverage <50%).
#
# Outputs subclinical_screen.csv: per-subject screen labels for the 92 controls
# plus 25 control-pool samples missing from NMF (handled as 'unscored').
#
# Spec: docs/superpowers/specs/2026-04-27-healthy-control-audit-design.md
# Env: rnaseq

suppressPackageStartupMessages({
  library(data.table)
})

BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HCDIR  <- file.path(BASE, "RNA-seq/results/audit_sensitivity/healthy_control_audit")
NMF    <- file.path(BASE, "RNA-seq/results/subtypes/nmf_assignments.csv")
META   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")

dir.create(HCDIR, recursive = TRUE, showWarnings = FALSE)

cat("Loading controls_definitions.csv...\n")
defs <- fread(file.path(HCDIR, "controls_definitions.csv"))
ctrl <- defs[in_control_pool == TRUE]
cat(sprintf("  Control pool: %d subjects\n", nrow(ctrl)))

cat("Loading nmf_assignments.csv...\n")
nmf <- fread(NMF, select = c("sample_id", "dominant_program",
                              "P_Pro_inflammatory_score",
                              "P_Innate_immune_score",
                              "P_Parenchymal_score",
                              "P_lncRNA_score",
                              "P_Hepatic_metabolic_score",
                              "P_Stellate_myofibroblast_score"))
cat(sprintf("  NMF assignments: %d subjects\n", nrow(nmf)))

# Join NMF to control pool
ctrl_nmf <- merge(ctrl, nmf, by = "sample_id", all.x = TRUE)
cat(sprintf("  Control pool with NMF coverage: %d / %d (%.1f%%)\n",
            sum(!is.na(ctrl_nmf$P_Pro_inflammatory_score)),
            nrow(ctrl_nmf),
            100 * sum(!is.na(ctrl_nmf$P_Pro_inflammatory_score)) / nrow(ctrl_nmf)))

# --- Validate axis: Spearman(P1, NAS) in disease set ---
cat("\nValidating P_Pro_inflammatory_score vs NAS (disease subjects, sanity check)\n")
m_all <- fread(META)
disease <- m_all[diagnosis_harmonized %in% c("Borderline", "NAFL", "NASH") &
                  !is.na(nas_score)]
disease_nmf <- merge(disease, nmf, by = "sample_id")
sp_p1 <- cor(disease_nmf$P_Pro_inflammatory_score, disease_nmf$nas_score,
             method = "spearman", use = "pairwise.complete.obs")
sp_p3 <- cor(disease_nmf$P_Parenchymal_score, disease_nmf$nas_score,
             method = "spearman", use = "pairwise.complete.obs")
cat(sprintf("  Spearman(P_Pro_inflammatory, NAS) = %.3f\n", sp_p1))
cat(sprintf("  Spearman(P_Parenchymal,        NAS) = %.3f\n", sp_p3))

# --- Tertile assignment on control pool ---
ctrl_scored <- ctrl_nmf[!is.na(P_Pro_inflammatory_score)]
ctrl_unscored <- ctrl_nmf[is.na(P_Pro_inflammatory_score)]

q_p1 <- quantile(ctrl_scored$P_Pro_inflammatory_score,
                 probs = c(1/3, 2/3), na.rm = TRUE)
cat(sprintf("\nTertile cuts on P_Pro_inflammatory_score (control pool, n=%d):\n",
            nrow(ctrl_scored)))
cat(sprintf("  Lower 1/3 cut: %.4f\n", q_p1[1]))
cat(sprintf("  Upper 2/3 cut: %.4f\n", q_p1[2]))

ctrl_scored[, p1_tertile := fifelse(P_Pro_inflammatory_score <= q_p1[1], "clean",
                                    fifelse(P_Pro_inflammatory_score >= q_p1[2], "suspect",
                                            "intermediate"))]
ctrl_scored[, p1_tertile := factor(p1_tertile, levels = c("clean","intermediate","suspect"))]

cat("\nTertile assignment summary (control pool):\n")
print(table(ctrl_scored$p1_tertile, ctrl_scored$dataset))

# Diagnosis_harmonized × tertile (does suspect tertile enrich for NAFL-labeled?)
cat("\nDiagnosis × tertile cross-tab:\n")
print(table(ctrl_scored$diagnosis_harmonized, ctrl_scored$p1_tertile))

# Suspect tertile enrichment for low-NAS NAFL (the suspected contamination class)
n_susp_nafl <- sum(ctrl_scored$p1_tertile == "suspect" &
                    ctrl_scored$diagnosis_harmonized == "NAFL")
n_clean_nafl <- sum(ctrl_scored$p1_tertile == "clean" &
                     ctrl_scored$diagnosis_harmonized == "NAFL")
n_susp_ctrl <- sum(ctrl_scored$p1_tertile == "suspect" &
                    ctrl_scored$diagnosis_harmonized == "Control")
n_clean_ctrl <- sum(ctrl_scored$p1_tertile == "clean" &
                     ctrl_scored$diagnosis_harmonized == "Control")
cat(sprintf("\nNAFL-labeled enrichment in suspect vs clean tertile:\n"))
cat(sprintf("  suspect: %d NAFL / %d Control\n", n_susp_nafl, n_susp_ctrl))
cat(sprintf("  clean:   %d NAFL / %d Control\n", n_clean_nafl, n_clean_ctrl))

if (n_susp_ctrl > 0 && n_clean_nafl > 0) {
  ft <- fisher.test(matrix(c(n_susp_nafl, n_susp_ctrl, n_clean_nafl, n_clean_ctrl),
                            nrow = 2))
  cat(sprintf("  Fisher OR (suspect-vs-clean for NAFL-vs-Control): %.2f, p=%.3g\n",
              ft$estimate, ft$p.value))
}

# --- Add unscored subjects back ---
ctrl_unscored[, p1_tertile := factor("unscored",
                                       levels = c("clean","intermediate","suspect","unscored"))]
ctrl_scored[, p1_tertile := factor(as.character(p1_tertile),
                                    levels = c("clean","intermediate","suspect","unscored"))]
out <- rbindlist(list(ctrl_scored, ctrl_unscored), use.names = TRUE, fill = TRUE)

# Write subclinical_screen.csv
out_cols <- c("sample_id", "dataset", "diagnosis_harmonized",
              "fibrosis_stage", "nas_score", "sex", "age",
              "obesity_proxy", "recruitment_note",
              "P_Pro_inflammatory_score", "P_Parenchymal_score",
              "P_Stellate_myofibroblast_score", "p1_tertile")
out_cols <- intersect(out_cols, names(out))
fwrite(out[, ..out_cols], file.path(HCDIR, "subclinical_screen.csv"))
cat(sprintf("\nWrote: subclinical_screen.csv (%d rows)\n", nrow(out)))

# --- Validation summary ---
val <- data.table(
  metric = c("control_pool_total",
             "nmf_covered",
             "spearman_p1_nas_disease",
             "spearman_p3_nas_disease",
             "tertile_cut_low",
             "tertile_cut_high",
             "n_clean", "n_intermediate", "n_suspect", "n_unscored"),
  value = c(nrow(ctrl_nmf),
            nrow(ctrl_scored),
            round(sp_p1, 4),
            round(sp_p3, 4),
            round(q_p1[1], 4),
            round(q_p1[2], 4),
            sum(out$p1_tertile == "clean"),
            sum(out$p1_tertile == "intermediate"),
            sum(out$p1_tertile == "suspect"),
            sum(out$p1_tertile == "unscored"))
)
fwrite(val, file.path(HCDIR, "subclinical_screen_validation.csv"))
cat("\nValidation summary:\n")
print(val)

cat("\nDone.\n")
