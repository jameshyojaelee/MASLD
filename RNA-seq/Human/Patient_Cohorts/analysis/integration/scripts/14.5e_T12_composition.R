#!/usr/bin/env Rscript
# 14.5e_T12_composition.R
# ---------------------------------------------------------------------------
# A7 T12: Composition decomposition for sex-interaction LFC
#
# Question: how much of the sex-interaction signal is driven by cell-type
# composition (hepatocyte / macrophage fractions) differing by sex?
#
# Approach:
#   1. Use BayesPrism (and MuSiC where available) per-sample cell fractions.
#   2. Compute, per gene, the disease-effect contrast in females and in males
#      from the existing stratified DE tables.
#   3. For each gene, regress:
#         interaction_LFC[gene] ~ (mean Δhep)_F * sex + (mean Δmac)_F * sex
#      Actually we do this at the *sample* / *cohort* level: build
#      delta_hep_by_sex (Disease_mean - Control_mean within each cohort x sex)
#      and quantify how much sex_class variance is explained by composition
#      shift.
#
# Concretely we:
#   a) Read unified_bayesprism_proportions.csv (sample-level fractions).
#   b) Merge with meta_matched + info to get sex, group, dataset.
#   c) For each (cohort, sex) cell, compute mean Δhep = mean(Disease) -
#      mean(Control) and Δmac similarly.
#   d) Build a sex composition shift contrast: (Δhep_F - Δhep_M),
#      (Δmac_F - Δmac_M) per cohort.
#   e) Cross-tabulate vs n_sex_dimorphic by cohort (loo_summary from T11
#      if available, else canonical only). Report Pearson r.
#   f) Gene-level analog: for each gene we don't have per-sample LFC, but we
#      can regress logFC_F - logFC_M (interaction proxy) against
#      gene-level "hepatocyte specificity" (whether gene is a hepatocyte
#      marker) and macrophage specificity. We use Bayesprism attribution
#      scores (per-gene Hepatocyte/Macrophage contribution) as proxies.
#
# Outputs:
#   audit_sensitivity/sex_composition_decomp/composition_regression.csv
#   audit_sensitivity/sex_composition_decomp/composition_shifts_by_cohort.csv
#   audit_sensitivity/sex_composition_decomp/summary.csv
# ---------------------------------------------------------------------------

set.seed(42)

suppressPackageStartupMessages({ library(data.table); library(yaml) })

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
DECONV <- file.path(INT, "results/deconvolution/bayesprism")
OUT  <- file.path(BASE, "RNA-seq/results/audit_sensitivity/sex_composition_decomp")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat("=== T12: Composition decomposition ===\n")

# --- Load sex classification + interaction dream ---
sex_cls <- fread(file.path(RDIR, "sex_deg_classification.csv"))
canon   <- fread(file.path(RDIR, "sex_interaction_dream.csv"))
cat("sex_cls rows:", nrow(sex_cls), "  canon rows:", nrow(canon), "\n")

# --- Load BayesPrism proportions ---
unif_path <- file.path(DECONV, "unified_bayesprism_proportions.csv")
attr_path <- file.path(DECONV, "bayesprism_attribution_scores.csv")
if (!file.exists(unif_path)) {
  stop("Missing: ", unif_path)
}
fract <- fread(unif_path)
cat("BayesPrism fractions:", nrow(fract), "samples,", ncol(fract), "cols\n")
cat("Columns:", paste(names(fract), collapse = ", "), "\n")

# Identify hep + mac columns flexibly
hep_col <- grep("(?i)hep", names(fract), value = TRUE)[1]
mac_col <- grep("(?i)mac|kupffer", names(fract), value = TRUE)[1]
cat("Hepatocyte col:", hep_col, "  Macrophage col:", mac_col, "\n")
if (is.na(hep_col) || is.na(mac_col))
  stop("Could not auto-detect hep/mac columns in unified_bayesprism_proportions.csv")

# Need sample_id and dataset
sid_col <- grep("(?i)sample|^id$", names(fract), value = TRUE)[1]
if (is.na(sid_col)) sid_col <- names(fract)[1]
cat("Sample id col:", sid_col, "\n")
setnames(fract, sid_col, "sample_id")

# --- Merge with metadata ---
ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
meta_match <- readRDS(file.path(RDIR, "meta_matched.rds"))
meta_dt <- as.data.table(meta_match)
fract_meta <- merge(fract[, .(sample_id, hep = get(hep_col), mac = get(mac_col))],
                    meta_dt[, .(sample_id, dataset, group_binary, inferred_sex)],
                    by = "sample_id")
fract_meta <- fract_meta[dataset %in% mega_cohorts]
fract_meta <- fract_meta[!is.na(inferred_sex)]
cat("Fractions x metadata merged:", nrow(fract_meta), "samples\n")
cat("Group x sex x dataset:\n"); print(table(fract_meta$group_binary, fract_meta$inferred_sex, fract_meta$dataset))

# --- Sex-stratified composition shifts per cohort ---
shifts <- fract_meta[, .(
  mean_hep = mean(hep, na.rm = TRUE),
  mean_mac = mean(mac, na.rm = TRUE),
  n        = .N
), by = .(dataset, group_binary, inferred_sex)]

# Pivot to wide: hep Disease vs Control per sex per cohort
shifts_w <- dcast(shifts, dataset + inferred_sex ~ group_binary,
                  value.var = c("mean_hep", "mean_mac"))
shifts_w[, delta_hep := mean_hep_Disease - mean_hep_Control]
shifts_w[, delta_mac := mean_mac_Disease - mean_mac_Control]
cat("\nDelta(hep) and Delta(mac) per cohort x sex:\n")
print(shifts_w)

# Sex-specific composition shift (F - M)
shifts_FM <- dcast(shifts_w[, .(dataset, inferred_sex, delta_hep, delta_mac)],
                   dataset ~ inferred_sex, value.var = c("delta_hep", "delta_mac"))
# Robust to F/Female and M/Male labels
fnm <- grep("delta_hep_F", names(shifts_FM), value = TRUE)[1]
mnm <- grep("delta_hep_M", names(shifts_FM), value = TRUE)[1]
fnm_m <- grep("delta_mac_F", names(shifts_FM), value = TRUE)[1]
mnm_m <- grep("delta_mac_M", names(shifts_FM), value = TRUE)[1]
shifts_FM[, sex_diff_delta_hep := get(fnm) - get(mnm)]
shifts_FM[, sex_diff_delta_mac := get(fnm_m) - get(mnm_m)]
fwrite(shifts_FM, file.path(OUT, "composition_shifts_by_cohort.csv"))
cat("Saved composition_shifts_by_cohort.csv\n"); print(shifts_FM)

# --- Gene-level composition regression ---
# For each gene we regress interaction_logFC ~ hep_attr + mac_attr (per-gene
# bayesprism attribution scores). If composition explains most variance,
# slope on attr will be large and R^2 will absorb most of the LFC variance.
if (file.exists(attr_path)) {
  attr_dt <- fread(attr_path)
  cat("BayesPrism attribution rows:", nrow(attr_dt),
      " cols:", paste(names(attr_dt), collapse = ", "), "\n")

  # Heuristic: detect gene + hep + mac columns
  gcol <- grep("(?i)gene|ensembl|symbol", names(attr_dt), value = TRUE)[1]
  if (is.na(gcol)) gcol <- names(attr_dt)[1]
  setnames(attr_dt, gcol, "gene")
  attr_hep <- grep("(?i)hep", names(attr_dt), value = TRUE)[1]
  attr_mac <- grep("(?i)mac|kupffer", names(attr_dt), value = TRUE)[1]
  cat("Attribution hep col:", attr_hep, "  mac col:", attr_mac, "\n")

  if (!is.na(attr_hep) && !is.na(attr_mac)) {
    setnames(attr_dt, attr_hep, "attr_hep")
    setnames(attr_dt, attr_mac, "attr_mac")
    keep_cols <- c("gene", "attr_hep", "attr_mac")
    attr_use <- attr_dt[, ..keep_cols]

    canon_lite <- canon[, .(gene, interaction_logFC = logFC, interaction_padj = padj)]
    mrg <- merge(canon_lite, attr_use, by = "gene")
    cat("Genes with both interaction LFC and bayesprism attr:", nrow(mrg), "\n")

    mrg <- mrg[is.finite(interaction_logFC) & is.finite(attr_hep) & is.finite(attr_mac)]
    fit_full <- lm(interaction_logFC ~ attr_hep + attr_mac, data = mrg)
    r2 <- summary(fit_full)$r.squared
    coef_tbl <- as.data.table(summary(fit_full)$coefficients, keep.rownames = "term")
    coef_tbl[, R2 := r2]
    fwrite(coef_tbl, file.path(OUT, "composition_regression.csv"))
    cat("Saved composition_regression.csv  R^2 =", round(r2, 4), "\n")

    # Per-sex_class composition signature
    cls <- sex_cls[, .(gene, sex_class)]
    mrg_cls <- merge(mrg, cls, by = "gene")
    summ <- mrg_cls[, .(
      n = .N,
      mean_attr_hep = mean(attr_hep, na.rm = TRUE),
      mean_attr_mac = mean(attr_mac, na.rm = TRUE),
      median_attr_hep = median(attr_hep, na.rm = TRUE),
      median_attr_mac = median(attr_mac, na.rm = TRUE)
    ), by = sex_class][order(-n)]
    fwrite(summ, file.path(OUT, "summary.csv"))
    cat("\nPer-class composition attribution:\n"); print(summ)
  } else {
    cat("Could not auto-detect attribution hep/mac cols.\n")
    fwrite(data.table(), file.path(OUT, "composition_regression.csv"))
    fwrite(data.table(), file.path(OUT, "summary.csv"))
  }
} else {
  cat("WARN: attribution scores file missing, gene-level skipped.\n")
  fwrite(data.table(), file.path(OUT, "composition_regression.csv"))
  fwrite(data.table(), file.path(OUT, "summary.csv"))
}

cat("Done:", as.character(Sys.time()), "\n")
