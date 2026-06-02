#!/usr/bin/env Rscript
# tmm_directional_bias_audit.R (Agent B4 — 2026-05-20)
# ---------------------------------------------------------------------------
# Per-cohort directional-bias audit for TMM assumption violation.
#
# RATIONALE
#   TMM normalization assumes most genes are NOT differentially expressed,
#   and that DE genes are balanced between up- and down-regulated. A cohort
#   with severe global directional bias (e.g., >70% of expressed genes shifted
#   the same direction in cases vs controls) likely violates TMM's assumption.
#
# DESIGN
#   - For each cohort in the 5-cohort mega-analysis:
#     - Compute per-gene mean log-CPM in Disease vs Control.
#     - Report fraction up vs down at |delta| > 0.5 log2-CPM, > 1.0.
#     - Report median per-sample TMM norm.factor (deviation from 1.0).
#   - Flag cohorts with > 60% directional skew at |delta| > 0.5 as "TMM-risk".
#
# OUTPUT
#   RNA-seq/results/audit_sensitivity/tmm_vs_uq/per_cohort_directional_bias.csv
#   RNA-seq/results/audit_sensitivity/tmm_vs_uq/REPORT.md
# ---------------------------------------------------------------------------

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(yaml)
})

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(PROJ, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
OUT_DIR <- file.path(PROJ, "RNA-seq/results/audit_sensitivity/tmm_vs_uq")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Loading merged_dge.rds (TMM norm factors)...\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

ycfg <- yaml::read_yaml(file.path(PROJ, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
dge <- dge[, dge$samples$dataset %in% mega_cohorts]

cohorts <- sort(unique(dge$samples$dataset))
cat("Cohorts:", paste(cohorts, collapse = ", "), "\n")

rows <- list()
for (coh in cohorts) {
  idx <- which(dge$samples$dataset == coh)
  d_c <- dge[, idx]
  grp <- d_c$samples$group_binary
  n_ctrl <- sum(grp == "Control"); n_dis <- sum(grp == "Disease")
  if (n_ctrl < 2 || n_dis < 2) {
    cat("Skipping", coh, "(insufficient samples per group)\n")
    next
  }
  logcpm <- edgeR::cpm(d_c, log = TRUE, prior.count = 1)
  delta  <- rowMeans(logcpm[, grp == "Disease", drop = FALSE]) -
            rowMeans(logcpm[, grp == "Control", drop = FALSE])
  # Expressed-gene mask: mean log-CPM > 0 in either arm
  mu_any <- pmax(rowMeans(logcpm[, grp == "Disease", drop = FALSE]),
                 rowMeans(logcpm[, grp == "Control", drop = FALSE]))
  expr_mask <- mu_any > 0
  delta_e <- delta[expr_mask]

  frac_up_05   <- mean(delta_e >  0.5)
  frac_dn_05   <- mean(delta_e < -0.5)
  frac_up_10   <- mean(delta_e >  1.0)
  frac_dn_10   <- mean(delta_e < -1.0)
  skew_05 <- if ((frac_up_05 + frac_dn_05) > 0) {
    max(frac_up_05, frac_dn_05) / (frac_up_05 + frac_dn_05)
  } else NA_real_
  median_nf <- median(d_c$samples$norm.factors)
  mad_nf    <- mad(d_c$samples$norm.factors)

  rows[[coh]] <- data.table(
    cohort = coh, n_ctrl = n_ctrl, n_disease = n_dis,
    n_expressed_genes = sum(expr_mask),
    frac_up_lfc05  = frac_up_05,
    frac_dn_lfc05  = frac_dn_05,
    frac_up_lfc10  = frac_up_10,
    frac_dn_lfc10  = frac_dn_10,
    directional_skew_lfc05 = skew_05,
    median_norm_factor = median_nf,
    mad_norm_factor    = mad_nf,
    tmm_risk_flag = !is.na(skew_05) && skew_05 > 0.60
  )
}
audit_dt <- rbindlist(rows)
fwrite(audit_dt, file.path(OUT_DIR, "per_cohort_directional_bias.csv"))
cat("Wrote per_cohort_directional_bias.csv\n"); print(audit_dt)

# ---- REPORT.md ----
jacc_path <- file.path(OUT_DIR, "jaccard_summary.csv")
jacc_block <- if (file.exists(jacc_path)) {
  jdt <- fread(jacc_path)
  paste(c("| metric | value |", "|---|---|",
          apply(jdt, 1, function(r) sprintf("| %s | %s |", r[["metric"]],
                                            format(as.numeric(r[["value"]]), digits = 4)))),
        collapse = "\n")
} else "_(jaccard_summary.csv not yet produced; run 14_5_uq_normalization_sensitivity.R first)_"

audit_md <- paste(c(
  "| cohort | n_ctrl | n_dis | n_expr | up>0.5 | dn>0.5 | skew | med_NF | flag |",
  "|---|---|---|---|---|---|---|---|---|",
  apply(audit_dt, 1, function(r) sprintf(
    "| %s | %s | %s | %s | %.3f | %.3f | %.3f | %.3f | %s |",
    r[["cohort"]], r[["n_ctrl"]], r[["n_disease"]], r[["n_expressed_genes"]],
    as.numeric(r[["frac_up_lfc05"]]), as.numeric(r[["frac_dn_lfc05"]]),
    as.numeric(r[["directional_skew_lfc05"]]),
    as.numeric(r[["median_norm_factor"]]),
    ifelse(as.logical(r[["tmm_risk_flag"]]), "**RISK**", "ok")))),
  collapse = "\n")

report <- sprintf(
"# TMM vs UQ sensitivity (B4 — 2026-05-20)

## Background
TMM (Robinson & Oshlack 2010) assumes most genes are NOT differentially
expressed and that DE genes are balanced up/down. The MASLD mega-analysis
shows broad metabolic remodeling; this audit asks whether any cohort
violates that assumption strongly enough to alter the canonical DEG set.

## TMM vs UQ DEG concordance
%s

Interpretation: Jaccard > 0.85 on Tier-1 DEGs (padj<0.05, |LFC|>0.5) is the
de-facto bar for normalization-robustness in this project. Directional
concordance among shared DEGs > 0.99 is expected.

## Per-cohort directional-bias audit
%s

`skew` = max(frac_up, frac_dn) / (frac_up + frac_dn) at |delta log2-CPM| > 0.5
among expressed genes. Flag threshold = skew > 0.60.

`med_NF` = median TMM norm.factor per cohort. Strong systematic deviation
from 1.0 indicates global library composition shift that TMM is attempting
to correct.

## Verdict
- If `n_uq_deg_lfc05` is within ±10 percent of `n_tmm_deg_lfc05` AND
  `jaccard_lfc05` > 0.85 AND `pearson_lfc_all_genes` > 0.97: TMM is robust;
  declare in Methods that UQ gives equivalent results (supplementary table).
- If any cohort flags RISK: discuss in §1.B of the manuscript and provide
  per-cohort sensitivity (LOCO with that cohort dropped).
- Otherwise: TMM accepted as canonical with this audit as documentation.

(Generated by `RNA-seq/scripts/tmm_directional_bias_audit.R`.)
", jacc_block, audit_md)

writeLines(report, file.path(OUT_DIR, "REPORT.md"))
cat("Wrote REPORT.md\n")
