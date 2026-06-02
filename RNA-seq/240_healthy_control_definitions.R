#!/usr/bin/env Rscript
# 240_healthy_control_definitions.R
#
# Healthy-control characterization pipeline, step 1 of 6.
#
# Re-derive the control pool (NAS<3 + fibrosis<=F1) and the obese-MASLD pool
# (NAS>=3 OR fibrosis>=F2) from unified_metadata.csv. Emit:
#   controls_definitions.csv   – per-subject membership in {control_pool, obese_MASLD_pool, other}
#   bmi_curation_template.csv  – control-pool subjects with placeholder BMI for manual fill
#                                from primary-publication supplementary tables.
#
# Cohort-context BMI proxy is also encoded as a fallback when curated BMI is
# absent: GSE162694 = bariatric (obese by recruitment), GSE126848 = lean
# controls, others = mixed.
#
# Spec: docs/superpowers/specs/2026-04-27-healthy-control-audit-design.md
# (plan: ~/.claude/plans/for-masld-i-composed-horizon.md)
# Env: rnaseq

suppressPackageStartupMessages({
  library(data.table)
})

BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
META   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
OUTDIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity/healthy_control_audit")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Cohort-context BMI proxy table
# Source: original cohort recruitment criteria (curated from publications).
# obesity_proxy = "obese" => recruitment biases to BMI>=30
# obesity_proxy = "lean"  => recruitment biases to BMI<25
# obesity_proxy = "mixed" => mixed BMI distribution
cohort_proxy <- data.table(
  dataset           = c("GSE126848", "GSE130970", "GSE135251", "GSE162694",
                        "GSE167523", "GSE174478", "GSE193066", "GSE213621",
                        "GSE240729", "PRJNA512027"),
  obesity_proxy    = c("lean",      "mixed",     "obese",     "obese",
                        "mixed",     "mixed",     "mixed",     "mixed",
                        "mixed",     "mixed"),
  recruitment_note = c("lean controls + obese NAFLD donors (Suppli 2019)",
                        "consecutive bariatric biopsies",
                        "Newcastle EPoS / European NAFLD biopsy registry",
                        "bariatric surgery cohort (high BMI by design)",
                        "general hepatology referral",
                        "histologically staged biopsies",
                        "Pavlides Hepatology multi-site",
                        "Ravi STELLAR pre-treatment fibrosis cohort",
                        "Gawrieh / Liver Research Program biopsies",
                        "Govaere / Hardy NAFLD progression")
)

# Reference DOIs / supp-table notes (for BMI manual curation guidance)
# Update as the team identifies the exact supplementary table per cohort.
publication_refs <- data.table(
  dataset = cohort_proxy$dataset,
  doi     = c("10.1016/j.jhep.2018.11.030",  # Suppli — GSE126848
              "10.1097/HEP.0000000000000037", # Hoang/Lake — GSE130970 (verify)
              "10.1126/scitranslmed.aba4448", # Govaere 2020 Sci Transl Med — GSE135251
              "10.1136/gutjnl-2020-322602",   # GSE162694 candidate
              "",                              # GSE167523 – TBD
              "",                              # GSE174478 – TBD
              "10.1002/hep.32220",            # Pavlides — GSE193066 candidate
              "",                              # GSE213621 – TBD (STELLAR)
              "",                              # GSE240729 – TBD
              "10.1016/j.jhep.2019.02.012"),  # Hardy – PRJNA512027
  supplementary_table_hint = c(
    "Suppl Table 1 – patient demographics",
    "manuscript Table 1 / Suppl",
    "Suppl Table S1 – clinical features",
    "manuscript Suppl Table 1",
    "GEO Series Matrix or Suppl S2",
    "GEO Series Matrix",
    "manuscript Suppl Table",
    "STELLAR trial publication Suppl",
    "GEO Series Matrix or original paper Suppl",
    "Hardy 2019 Suppl Table"
  )
)

cat("Loading unified metadata...\n")
m <- fread(META)
cat(sprintf("  Total samples: %d\n", nrow(m)))

# Numeric coercion of fibrosis_stage (handles "0", "F0", factor, etc.)
fib_num <- suppressWarnings(as.numeric(as.character(m$fibrosis_stage)))

# --- Pool A: control pool (clinical-healthy) ---
pool_ctrl <- !is.na(fib_num) & fib_num <= 1 & !is.na(m$nas_score) & m$nas_score < 3
cat(sprintf("\n[Pool A] Control pool (fibrosis<=F1 AND NAS<3): %d subjects\n", sum(pool_ctrl)))
cat("  Per-cohort breakdown:\n")
print(table(m$dataset[pool_ctrl]))
cat("  Per-diagnosis_harmonized breakdown:\n")
print(table(m$diagnosis_harmonized[pool_ctrl], useNA = "ifany"))

# --- Pool B: obese-MASLD pool (advanced disease) ---
pool_mas <- (!is.na(m$nas_score) & m$nas_score >= 3) |
            (!is.na(fib_num) & fib_num >= 2)
cat(sprintf("\n[Pool B] Obese-MASLD pool (NAS>=3 OR fibrosis>=F2): %d subjects\n", sum(pool_mas)))
cat("  Per-cohort breakdown:\n")
print(table(m$dataset[pool_mas]))

# --- Build controls_definitions.csv ---
defs <- copy(m)
defs[, fib_num := fib_num]
defs[, in_control_pool := pool_ctrl]
defs[, in_obese_MASLD_pool := pool_mas]

# Pool label (priority: control < obese-MASLD; otherwise other/intermediate)
defs[, pool_label := fifelse(in_control_pool, "control_pool",
                             fifelse(in_obese_MASLD_pool, "obese_MASLD_pool",
                                     "intermediate"))]

# Join cohort proxy
defs <- merge(defs, cohort_proxy, by = "dataset", all.x = TRUE)

# Encode preliminary 3-way assignment using cohort proxy as fallback BMI.
# Original "histology-pool" proxy (kept for backwards-compatibility):
defs[, group3_proxy := fifelse(pool_label == "control_pool" & obesity_proxy == "obese", "resilient_proxy",
                       fifelse(pool_label == "control_pool" & obesity_proxy == "lean", "lean_healthy_proxy",
                       fifelse(pool_label == "control_pool" & obesity_proxy == "mixed", "mixed_healthy",
                       fifelse(pool_label == "obese_MASLD_pool", "obese_MASLD_proxy",
                               "out_of_scope"))))]

# --- Richer condition-based definitions (priority over histology pool) ---
# Resilient (gold-standard, n_explicit ~12+31=43):
#   GSE126848 Control_Obese (explicit obese-no-MASLD per Suppli 2019)
#   GSE162694 Control (bariatric, explicit normal histology per cohort)
# Resilient_broad (presumed obese-no-MASLD, +n ~137):
#   Add Control labels from bariatric/NAFLD-clinic cohorts:
#     GSE130970, GSE135251, GSE213621, PRJNA512027
# Lean_healthy_explicit (only confirmed lean cohort):
#   GSE126848 Control (n=14)
# Obese_MASLD: NAFL/NASH labels OR Fibrosis_F2+ from any cohort

# Bariatric / obesity-recruited cohorts (Controls from these are presumed obese)
bariatric_cohorts <- c("GSE126848",  # Control_Obese only — Suppli 2019
                        "GSE130970", "GSE135251", "GSE162694",
                        "GSE213621", "PRJNA512027")

defs[, is_bariatric_cohort := dataset %in% bariatric_cohorts]

defs[, group3_explicit := fifelse(
  # GSE126848-specific
  dataset == "GSE126848" & condition == "Control_Obese", "resilient_explicit",
  fifelse(dataset == "GSE126848" & condition == "Control", "lean_healthy_explicit",
  # Other bariatric cohorts: Controls = presumed obese (resilient_broad)
  fifelse(is_bariatric_cohort & condition == "Control" & dataset != "GSE126848",
          "resilient_broad",
  # Pool B (advanced disease) labelled obese MASLD
  fifelse(in_obese_MASLD_pool == TRUE, "obese_MASLD",
  # Other Controls (no clear obesity context) → ambiguous_healthy
  fifelse(condition == "Control", "ambiguous_healthy",
          "intermediate_or_other")))))]

# Combined "resilient_combined" flag (resilient_explicit + resilient_broad)
defs[, is_resilient := group3_explicit %in% c("resilient_explicit", "resilient_broad")]
defs[, is_lean_healthy := group3_explicit == "lean_healthy_explicit"]
defs[, is_obese_MASLD := group3_explicit == "obese_MASLD"]

# Reorder columns
out_cols <- c("sample_id", "dataset", "diagnosis_harmonized", "condition",
              "group_binary", "sex", "age", "fibrosis_stage", "fib_num",
              "nas_score", "obesity_proxy", "recruitment_note",
              "is_bariatric_cohort",
              "in_control_pool", "in_obese_MASLD_pool",
              "pool_label", "group3_proxy", "group3_explicit",
              "is_resilient", "is_lean_healthy", "is_obese_MASLD")
fwrite(defs[, ..out_cols],
       file.path(OUTDIR, "controls_definitions.csv"))
cat(sprintf("\nWrote: %s (%d rows, %d cols)\n",
            file.path(OUTDIR, "controls_definitions.csv"), nrow(defs), length(out_cols)))

# --- Group-3 proxy summary (histology-pool fallback) ---
cat("\n3-way proxy (histology-pool fallback) summary:\n")
print(table(defs$group3_proxy))

cat("\n3-way proxy by dataset:\n")
print(table(defs$dataset, defs$group3_proxy))

# --- Group-3 explicit summary (condition-based, primary) ---
cat("\n3-way explicit (condition-based, primary) summary:\n")
print(table(defs$group3_explicit))

cat("\n3-way explicit by dataset:\n")
print(table(defs$dataset, defs$group3_explicit))

cat("\nResilient pool composition (gold + broad):\n")
print(defs[is_resilient == TRUE, .N, by = .(dataset, condition, group3_explicit)][order(-N)])
cat(sprintf("\nResilient total: %d  (explicit=%d, broad=%d)\n",
            sum(defs$is_resilient),
            sum(defs$group3_explicit == "resilient_explicit"),
            sum(defs$group3_explicit == "resilient_broad")))
cat(sprintf("Lean_healthy_explicit (GSE126848 Control): %d\n",
            sum(defs$is_lean_healthy)))
cat(sprintf("Obese_MASLD: %d\n", sum(defs$is_obese_MASLD)))

# --- BMI curation template ---
ctrl <- defs[in_control_pool == TRUE, .(sample_id, dataset, sex, age,
                                          fibrosis_stage, nas_score,
                                          diagnosis_harmonized,
                                          obesity_proxy, recruitment_note)]
ctrl <- merge(ctrl, publication_refs, by = "dataset", all.x = TRUE)
ctrl[, `:=`(
  bmi             = NA_real_,
  bmi_source      = NA_character_,
  bmi_curated_by  = NA_character_,
  bmi_curated_on  = NA_character_,
  notes           = NA_character_
)]
setcolorder(ctrl, c("sample_id", "dataset", "diagnosis_harmonized",
                     "fibrosis_stage", "nas_score", "sex", "age",
                     "obesity_proxy", "recruitment_note",
                     "doi", "supplementary_table_hint",
                     "bmi", "bmi_source", "bmi_curated_by", "bmi_curated_on", "notes"))

fwrite(ctrl, file.path(OUTDIR, "bmi_curation_template.csv"))
cat(sprintf("\nWrote: %s (%d rows; manual fill needed for `bmi` column)\n",
            file.path(OUTDIR, "bmi_curation_template.csv"), nrow(ctrl)))

# --- Also emit a parallel BMI template for the obese-MASLD pool (smaller subset) ---
mas <- defs[in_obese_MASLD_pool == TRUE, .(sample_id, dataset, sex, age,
                                             fibrosis_stage, nas_score,
                                             diagnosis_harmonized,
                                             obesity_proxy, recruitment_note)]
mas <- merge(mas, publication_refs, by = "dataset", all.x = TRUE)
mas[, `:=`(bmi = NA_real_, bmi_source = NA_character_,
           bmi_curated_by = NA_character_, bmi_curated_on = NA_character_,
           notes = NA_character_)]
setcolorder(mas, c("sample_id", "dataset", "diagnosis_harmonized",
                    "fibrosis_stage", "nas_score", "sex", "age",
                    "obesity_proxy", "recruitment_note",
                    "doi", "supplementary_table_hint",
                    "bmi", "bmi_source", "bmi_curated_by", "bmi_curated_on", "notes"))
fwrite(mas, file.path(OUTDIR, "bmi_curation_template_MASLD.csv"))
cat(sprintf("\nWrote: %s (%d rows; obese-MASLD pool template for matched-BMI controls)\n",
            file.path(OUTDIR, "bmi_curation_template_MASLD.csv"), nrow(mas)))

# --- Write per-cohort manual-curation effort estimate ---
effort <- ctrl[, .(n_to_curate = .N), by = dataset]
effort <- merge(effort, publication_refs, by = "dataset", all.x = TRUE)
effort[, est_hours := pmax(2, ceiling(n_to_curate / 10))]
fwrite(effort, file.path(OUTDIR, "bmi_curation_effort_estimate.csv"))
cat("\nManual curation effort estimate (control pool only):\n")
print(effort[, .(dataset, n_to_curate, est_hours, doi)])
cat(sprintf("\nTotal estimated effort: %d hours across %d cohorts\n",
            sum(effort$est_hours), nrow(effort)))

# Save session info for provenance
sink(file.path(OUTDIR, "240_session_info.txt"))
cat("Date:", as.character(Sys.time()), "\n\n")
print(sessionInfo())
sink()

cat("\nDone.\n")
