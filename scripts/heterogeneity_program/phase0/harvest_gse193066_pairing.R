#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Phase 0d — GSE193066 paired-biopsy harvest → paired_193066.tsv
#
# The paired-trajectory call already exists (results/reversal/
# gse193066_paired_trajectory.csv: 58 patients × 2 timepoints, fib_delta +
# traj_class). This script reformats it into the canonical Phase-0 artifact
# `paired_193066.tsv`, keyed to the harmonized `sample_id` (= SRA run), normalises
# the trajectory labels to {progressor, stable, regressor}, and VALIDATES the
# join against unified_metadata.csv so downstream T2/T3 can score frozen
# signatures on the bulk samples and test progression/reversal.
#
# Rigor note (GATE-K): this arm is SUPPORTING ONLY (n≈15 progressor / 15
# regressor). No n=15 result is a headline; a held-out cohort carries the primary.
#
# Output: RNA-seq/results/heterogeneity_program/phase0/paired_193066.tsv
# Env:    rnaseq.  Lightweight (116-row reformat) — safe off-node.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages(library(data.table))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
traj_path <- file.path(BASE, "RNA-seq/results/reversal/gse193066_paired_trajectory.csv")
meta_path <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
out_path  <- file.path(BASE, "RNA-seq/results/heterogeneity_program/phase0/paired_193066.tsv")

stopifnot(file.exists(traj_path))
tr <- fread(traj_path)
cat(sprintf("Loaded trajectory: %d rows, %d patients\n", nrow(tr), uniqueN(tr$patient)))

# normalise trajectory labels → {progressor, stable, regressor}
lab_map <- c(progress = "progressor", regress = "regressor", stable = "stable",
             progressor = "progressor", regressor = "regressor")
tr[, traj_class := lab_map[as.character(traj_class)]]
stopifnot(!any(is.na(tr$traj_class)))

paired <- tr[, .(
  sample_id     = run,           # harmonized id = SRA run (matches unified_metadata)
  patient,
  timepoint,                     # 1 = 1st biopsy, 2 = 2nd biopsy
  fibrosis_stage = fib,
  nas_score      = nas,
  fib_1st, fib_2nd, fib_delta,
  nas_1st, nas_2nd,
  risk_1st,
  traj_class,                    # patient-level: progressor / stable / regressor
  paired
)]

# ── Validate the join to the harmonized metadata (GATE: keyed correctly) ──────
um <- fread(meta_path, select = c("sample_id", "dataset"))
um193 <- um[dataset == "GSE193066"]
n_match <- sum(paired$sample_id %in% um193$sample_id)
cat(sprintf("Join check: %d/%d paired samples found in unified_metadata GSE193066 (n=%d)\n",
            n_match, nrow(paired), nrow(um193)))
if (n_match < nrow(paired)) {
  miss <- paired$sample_id[!paired$sample_id %in% um193$sample_id]
  cat("  UNMATCHED sample_ids (first 10): ", paste(head(miss, 10), collapse=", "), "\n")
}

# patient-level trajectory summary (the supporting-validation strata)
pat <- unique(paired[, .(patient, traj_class)])
cat("\nPatient-level trajectory classes:\n"); print(pat[, .N, by = traj_class][order(-N)])

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
fwrite(paired, out_path, sep = "\t")
cat(sprintf("\nWrote %s  (%d rows × %d cols)\n", out_path, nrow(paired), ncol(paired)))
