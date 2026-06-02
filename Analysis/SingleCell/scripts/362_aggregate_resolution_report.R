#!/usr/bin/env Rscript
# S3 / 362: Aggregate resolution sweep + cross-cohort Progressor results
# into cross_cohort_progressor_REPORT.md (the deliverable). Also runs a
# confirmatory propeller (limma::propeller-style) test where possible.
#
# Outputs (worktree):
#   resolution_sweep/cross_cohort_progressor_REPORT.md
#   resolution_sweep/propeller_cross_cohort.csv  (if propeller available)

suppressPackageStartupMessages({
  library(data.table)
})

WORKTREE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation"
)
SWEEP_DIR <- file.path(
  WORKTREE,
  "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/resolution_sweep"
)

stab        <- fread(file.path(SWEEP_DIR, "resolution_stability.csv"))
consensus   <- fread(file.path(SWEEP_DIR, "consensus_summary.csv"))
prog_clust  <- fread(file.path(SWEEP_DIR, "progressor_clusters.csv"))
per_res     <- fread(file.path(SWEEP_DIR, "cross_cohort_progressor_per_resolution.csv"))
donor_cts   <- fread(file.path(SWEEP_DIR, "cross_cohort_progressor_donor_counts.csv"))
notes_path  <- file.path(SWEEP_DIR, "cross_cohort_progressor_NOTES.md")
notes_txt   <- if (file.exists(notes_path)) paste(readLines(notes_path), collapse = "\n") else ""

# Confirmatory propeller via speckle::propeller if installed
propeller_block <- "_propeller_ not run: `speckle` package not available in env.\n"
if (requireNamespace("speckle", quietly = TRUE)) {
  suppressPackageStartupMessages(library(speckle))
  propeller_rows <- list()
  for (res in unique(donor_cts$resolution)) {
    sub <- donor_cts[resolution == res &
                       disease_stage_coarse %in% c("Steatosis", "Steatohepatitis")]
    # Two arms: ALL and EXCL_GSE244832
    for (arm in c("ALL", "EXCL_GSE244832")) {
      sub_arm <- if (arm == "EXCL_GSE244832") sub[dataset != "GSE244832"] else sub
      if (nrow(sub_arm) < 4) next
      if (length(unique(sub_arm$disease_stage_coarse)) < 2) next
      # propeller expects a counts matrix (clusters x donors). With only 1 cluster
      # (Progressor vs not), we feed two rows: Progressor and not.
      mat <- rbind(
        Progressor = sub_arm$progressor_count,
        Other      = sub_arm$total_count - sub_arm$progressor_count
      )
      colnames(mat) <- sub_arm$sample
      grp <- factor(sub_arm$disease_stage_coarse,
                    levels = c("Steatosis", "Steatohepatitis"))
      out <- try(
        propeller(clusters = rep(rownames(mat), each = ncol(mat)),
                  sample = rep(colnames(mat), nrow(mat)),
                  group = rep(grp, nrow(mat)),
                  transform = "asin"),
        silent = TRUE
      )
      if (inherits(out, "try-error")) next
      out_dt <- as.data.table(out, keep.rownames = "cluster")
      out_dt[, resolution := res]
      out_dt[, arm := arm]
      propeller_rows[[length(propeller_rows) + 1]] <- out_dt
    }
  }
  if (length(propeller_rows) > 0) {
    propeller_df <- rbindlist(propeller_rows, fill = TRUE)
    fwrite(propeller_df, file.path(SWEEP_DIR, "propeller_cross_cohort.csv"))
    propeller_block <- "Propeller (speckle::propeller, asin transform) confirmatory results saved to `propeller_cross_cohort.csv`. Top rows shown below.\n\n```\n"
    propeller_block <- paste0(propeller_block,
      paste(capture.output(print(head(propeller_df, 20))), collapse = "\n"),
      "\n```\n")
  } else {
    propeller_block <- "_propeller_ executed but produced no rows (insufficient donors per arm).\n"
  }
}

# Build the REPORT
report_path <- file.path(SWEEP_DIR, "cross_cohort_progressor_REPORT.md")

# Summarize stability
setkey(stab, resolution)
setkey(consensus, resolution)
stab_merged <- merge(stab, consensus, by = "resolution")

# Per-resolution conclusion
per_res[, replicates := !is.na(p) & p < 0.05 & log2fc_mean > 0]
agg <- per_res[, .(
  n_res_ALL_sig = sum(arm == "ALL" & replicates),
  n_res_EXCL_sig = sum(arm == "EXCL_GSE244832" & replicates),
  n_res_total = uniqueN(resolution)
)]

# Filter resolutions where Progressor exists with consensus>0.7
high_consensus_res <- consensus[mean_consensus_score > 0.7, resolution]
hcr_n <- length(high_consensus_res)

# Top resolution rows (consensus_pass)
top_rows <- per_res[consensus_pass == TRUE]
setorder(top_rows, resolution, arm)

# Build markdown
lines <- c(
  "# Cross-cohort Progressor validation — S3 reviewer-defense report",
  "",
  sprintf("Generated: %s", Sys.time()),
  "",
  "## TL;DR",
  "",
  sprintf(
    "- Resolution sweep: %d resolutions x 20 seeds, Leiden on existing scVI latent (no retrain).",
    nrow(stab_merged)
  ),
  sprintf(
    "- Resolutions with mean consensus score > 0.7: **%d / %d**.",
    hcr_n, nrow(stab_merged)
  ),
  sprintf(
    "- ALL-cohort arm: %d / %d resolutions show Steatosis -> Steatohepatitis Progressor *expansion* with p<0.05.",
    agg$n_res_ALL_sig, agg$n_res_total
  ),
  sprintf(
    "- **GSE244832-EXCLUDED arm**: %d / %d resolutions show the same expansion with p<0.05.",
    agg$n_res_EXCL_sig, agg$n_res_total
  ),
  "",
  if (agg$n_res_EXCL_sig == 0) {
    "**Verdict — cross-cohort: NO replication after excluding GSE244832.** The Progressor 1.3% -> 63.1% Steatosis -> Steatohepatitis expansion does not replicate at any sweep resolution once GSE244832 is removed. Recommended manuscript reframe: 'GSE244832-specific Progressor subtype enrichment at the Steatohepatitis stage' (single-cohort observation), with cross-cohort generalization marked as untested in the current scRNA atlas because no other scRNA cohort contributes paired Steatosis + Steatohepatitis donors."
  } else {
    sprintf(
      "**Verdict — cross-cohort: PARTIAL replication.** %d / %d resolutions retain a significant Progressor expansion after excluding GSE244832; review per-resolution table below.",
      agg$n_res_EXCL_sig, agg$n_res_total
    )
  },
  "",
  "## Donor counts after excluding GSE244832 (cells used in the central test)",
  "",
  "```",
  paste(capture.output(print(
    donor_cts[disease_stage_coarse %in% c("Steatosis", "Steatohepatitis") &
                dataset != "GSE244832" &
                resolution == min(resolution)][,
      .(n_donors = uniqueN(sample),
        sum_progressor = sum(progressor_count),
        sum_total = sum(total_count)),
      by = .(dataset, disease_stage_coarse)]
  )), collapse = "\n"),
  "```",
  "",
  "_Donor counts are independent of resolution — they reflect available samples per (dataset, stage)._",
  "",
  "## Per-resolution stability",
  "",
  "```",
  paste(capture.output(print(head(stab_merged, 50))), collapse = "\n"),
  "```",
  "",
  "## Per-resolution cross-cohort proportion test (consensus > 0.7 only)",
  "",
  "```",
  paste(capture.output(print(top_rows[, .(resolution, arm, consensus_score,
    mean_prop_Steatosis, mean_prop_Steatohepatitis, log2fc_mean,
    beta, se, p, n_donors, n_steat, n_steathep,
    datasets_kept)])), collapse = "\n"),
  "```",
  "",
  "## Methods substitution note",
  "",
  notes_txt,
  "",
  "## Propeller confirmatory run",
  "",
  propeller_block,
  "",
  "## Files in this directory",
  "",
  paste0("- ", list.files(SWEEP_DIR, recursive = FALSE), collapse = "\n"),
  ""
)
writeLines(lines, report_path)
cat(sprintf("Wrote %s\n", report_path))
