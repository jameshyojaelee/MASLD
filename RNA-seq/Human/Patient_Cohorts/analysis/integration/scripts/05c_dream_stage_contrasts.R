#!/usr/bin/env Rscript
# 05c_dream_stage_contrasts.R
# ---------------------------------------------------------------------------
# Stage-stratified bulk dream contrasts vs Healthy (Control).
#
# Three independent binary contrasts (one dream() fit each) used as the
# stage-matched anchor for sc CCC LR pair concordance (Script 348):
#
#   stage = Steatosis       : diagnosis_harmonized == "NAFL"  (NAS<3)
#                            OR nas_score in [0,2] when label is missing
#                            -> dream_results_stage_steatosis.csv
#
#   stage = Steatohepatitis : diagnosis_harmonized in {"NASH","Borderline"}
#                            (NAS>=3, borderline-grouped) AND
#                            (fibrosis_stage < 4 OR fibrosis_stage NA)
#                            -> dream_results_stage_sh.csv
#
#   stage = Cirrhosis       : fibrosis_stage == 4
#                            OR diagnosis_harmonized matches cirrhosis-coded
#                            condition labels
#                            -> dream_results_stage_cirrhosis.csv
#
# All three reuse the canonical dream formula from
# 05f_mash_vs_healthy_dream.R / 05g_masl_vs_healthy_dream.R:
#
#     ~ stage_status + sex_for_model + (1 | dataset_subbatch)
#
# Cohort filter: only cohorts with >=3 Control AND >=3 stage-target samples.
# (PRJNA512027 — with its L0/S0 library-prep confound — was permanently removed
# from the pipeline 2026-05-15.)
#
# Outputs (results/integration/, mirroring dream_results_ashr.csv schema):
#   logFC, AveExpr, t, P.Value, padj (=adj.P.Val), z.std, gene, symbol
#
# Usage:
#   micromamba run -n rnaseq Rscript 05c_dream_stage_contrasts.R
# SLURM: cpu --qos=interactive --time=04:00:00 --mem=128G --cpus-per-task=16
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
})

# Inject reformulas into lme4 namespace BEFORE loading variancePartition
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({
      unlockBinding(fn, ns_lme4)
      assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
      lockBinding(fn, ns_lme4)
    }, silent = TRUE)
  }
}

suppressPackageStartupMessages({
  library(variancePartition)
  library(BiocParallel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
PCDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts")

stopifnot(dir.exists(RDIR))

cat("=== 05c: Stage-stratified Dream contrasts vs Healthy ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# --- Load shared inputs ---
counts <- readRDS(file.path(RDIR, "merged_counts_raw.rds"))
meta   <- readRDS(file.path(RDIR, "meta_matched.rds"))
qc     <- fread(file.path(INT, "qc/sample_qc_report.csv"))
meta   <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

# Master review M-P0-8: enforce PRJNA512027 exclusion. The cohort was removed
# from human_datasets.yaml on 2026-05-15 for an L0/S0 library-prep batch
# confound that makes its disease-vs-control contrast unreliable, but the
# merged_counts_raw.rds / meta_matched.rds files were not rebuilt — so the
# samples leak back in here via the cohort-eligibility filter. Drop them
# explicitly to keep this script in sync with the rest of the pipeline.
EXCLUDED_COHORTS <- c("PRJNA512027")
n_before <- nrow(meta)
meta <- meta[!(dataset %in% EXCLUDED_COHORTS)]
n_after <- nrow(meta)
cat(sprintf("[cohort filter] dropped %d samples from %s (%d -> %d)\n",
            n_before - n_after, paste(EXCLUDED_COHORTS, collapse = ","),
            n_before, n_after))

# Ensembl -> symbol map (matches dream_results_ashr.csv `symbol` col)
sym_map <- fread(file.path(INT, "results/gene_annotation/human_ensg_to_symbol.tsv"))
setnames(sym_map, c("gene_id", "symbol", "gene_type", "gene_base"),
                  c("gene_id", "symbol", "gene_type", "gene_base"),
         skip_absent = TRUE)

cat("Loaded:", ncol(counts), "samples in count matrix,",
    nrow(meta), "after QC filter\n")
cat("diagnosis_harmonized table:\n")
print(meta[, .N, by = diagnosis_harmonized])
cat("fibrosis_stage table:\n")
print(meta[, .N, by = fibrosis_stage])

# Detect cirrhosis-coded condition labels (regex on condition field) ---------
cirrhosis_labels <- unique(meta$condition[
  grepl("cirrho|Cirrho|CIRRHO|F4", meta$condition, ignore.case = FALSE)])
cat("Condition labels detected as cirrhosis-coded:",
    paste(cirrhosis_labels, collapse = ", "), "\n\n")

# --- Stage definitions ---
# Steatosis (MASL): NAFL label OR NAS 0-2 (no Borderline; NAS<3)
meta[, steatosis_target := (
  diagnosis_harmonized == "NAFL" |
  (is.na(diagnosis_harmonized) & !is.na(nas_score) & nas_score < 3 &
   group_binary == "Disease")
)]

# Steatohepatitis (MASH borderline-grouped): NASH or Borderline label,
# excluding F4 (which is the Cirrhosis arm). Where label is missing,
# accept nas_score >= 3 (matches Kleiner harmonization).
meta[, sh_target := (
  (diagnosis_harmonized %in% c("NASH", "Borderline") |
   (is.na(diagnosis_harmonized) & !is.na(nas_score) & nas_score >= 3 &
    group_binary == "Disease")) &
  (is.na(fibrosis_stage) | fibrosis_stage < 4)
)]

# Cirrhosis: fibrosis_stage == 4 (highest-fidelity), OR explicit cirrhosis label
meta[, cirrhosis_target := (
  (!is.na(fibrosis_stage) & fibrosis_stage == 4) |
  (condition %in% cirrhosis_labels & !is.na(condition))
)]

# Healthy: Control (diagnosis_harmonized == "Control")
meta[, is_healthy := (diagnosis_harmonized == "Control")]

cat("Target counts per stage (before cohort filter):\n")
cat("  Steatosis      :", sum(meta$steatosis_target,  na.rm = TRUE), "\n")
cat("  Steatohepatitis:", sum(meta$sh_target,         na.rm = TRUE), "\n")
cat("  Cirrhosis      :", sum(meta$cirrhosis_target,  na.rm = TRUE), "\n")
cat("  Healthy        :", sum(meta$is_healthy,        na.rm = TRUE), "\n\n")

# ---------------------------------------------------------------------------
# Single-contrast dream wrapper. Mirrors 05f/05g structure exactly.
# ---------------------------------------------------------------------------
run_stage_dream <- function(target_flag, stage_label, out_fname) {
  cat("\n==============================================================\n")
  cat(sprintf("Stage contrast: %s vs Healthy\n", stage_label))
  cat("==============================================================\n")

  meta_s <- meta[get(target_flag) | is_healthy]
  meta_s[, stage_status := factor(
    fifelse(is_healthy, "Control", stage_label),
    levels = c("Control", stage_label))]

  # Per-cohort sample counts; keep cohorts with >=3 in BOTH arms
  cohort_counts <- meta_s[, .(
    n_ctrl   = sum(stage_status == "Control"),
    n_target = sum(stage_status == stage_label)),
    by = dataset]
  cat("Per-cohort counts:\n"); print(cohort_counts)
  eligible <- cohort_counts[n_ctrl >= 3 & n_target >= 3, dataset]
  cat("Eligible cohorts (>=3 each arm):",
      paste(sort(eligible), collapse = ", "), "\n")
  if (length(eligible) < 1) {
    cat("[WARN] no eligible cohorts; skipping", stage_label, "\n")
    return(invisible(NULL))
  }
  meta_s <- meta_s[dataset %in% eligible][order(sample_id)]
  cat("Final samples:", nrow(meta_s), "\n")

  idx <- colnames(counts) %in% meta_s$sample_id
  dge <- DGEList(counts = counts[, idx])
  dge$samples <- cbind(dge$samples,
    meta_s[match(colnames(dge), meta_s$sample_id),
           .(stage_status, dataset, sex, inferred_sex, age)])

  # dataset_subbatch == dataset (PRJNA L0/S0 subbatch logic permanently removed
  # 2026-05-15 with the cohort itself).
  dge$samples$dataset_subbatch <- as.character(dge$samples$dataset)

  # sex_for_model: annotated sex, fallback to inferred_sex
  dge$samples$sex_for_model <- dge$samples$sex
  na_sex <- is.na(dge$samples$sex_for_model) |
            dge$samples$sex_for_model == ""
  dge$samples$sex_for_model[na_sex] <-
    as.character(dge$samples$inferred_sex[na_sex])
  dge$samples$sex_for_model <- factor(dge$samples$sex_for_model)

  dge <- calcNormFactors(dge, method = "TMM")
  keep <- filterByExpr(dge, group = dge$samples$stage_status)
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  cat("Dream genes after filterByExpr:", nrow(dge), "\n")
  cat("Group distribution:\n")
  print(table(dge$samples$stage_status, dge$samples$dataset))
  cat("Sex levels:",
      paste(levels(dge$samples$sex_for_model), collapse = ", "), "\n")

  # Formula: drop sex if only one level
  if (nlevels(dge$samples$sex_for_model) > 1) {
    form <- ~ stage_status + sex_for_model + (1 | dataset_subbatch)
  } else {
    cat("[note] sex_for_model has only one level; dropping from formula\n")
    form <- ~ stage_status + (1 | dataset_subbatch)
  }
  # Drop random effect if only one cohort
  if (length(unique(dge$samples$dataset_subbatch)) < 2) {
    cat("[note] only one dataset_subbatch level; dropping random intercept\n")
    if (nlevels(dge$samples$sex_for_model) > 1) {
      form <- ~ stage_status + sex_for_model
    } else {
      form <- ~ stage_status
    }
  }
  cat("Formula:", deparse(form), "\n")

  n_cores <- min(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")), 32)
  cat("Using", n_cores, "CPU cores\n")
  BPPARAM <- MulticoreParam(n_cores, progressbar = FALSE)

  vobj <- voomWithDreamWeights(dge, form, dge$samples, BPPARAM = BPPARAM)
  cat("Running dream...\n")
  fit  <- dream(vobj, form, dge$samples, BPPARAM = BPPARAM)
  # Master review M-P0-8 follow-up: when dream() falls back to limma::lmFit
  # (single random-effect level, e.g. Cirrhosis arm with only GSE135251 after
  # PRJNA512027 exclusion), the resulting MArrayLM is NOT eBayes-shrunk and
  # topTable errors with "Need to run eBayes or treat first". eBayes is a
  # no-op on already-shrunk MArrayLM2 objects, so it's safe to call always.
  if (inherits(fit, "MArrayLM") && !inherits(fit, "MArrayLM2")) {
    cat("[note] dream fell back to limma::lmFit; running eBayes\n")
    fit <- limma::eBayes(fit)
  }
  coef_name <- paste0("stage_status", stage_label)
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)

  # Match dream_results_ashr.csv schema: rename adj.P.Val -> padj, add symbol
  setDT(tt)
  setnames(tt, "adj.P.Val", "padj")
  tt[, gene_base := sub("\\..*", "", gene)]
  if ("gene_base" %in% names(sym_map)) {
    tt <- merge(tt, sym_map[, .(gene_base, symbol)],
                by = "gene_base", all.x = TRUE)
  } else {
    # fallback: join by full gene_id
    tt <- merge(tt, sym_map[, .(gene_id, symbol)],
                by.x = "gene", by.y = "gene_id", all.x = TRUE)
  }
  tt[, gene_base := NULL]
  setcolorder(tt, intersect(c("gene", "symbol", "logFC", "AveExpr", "t",
                              "P.Value", "padj", "z.std", "B"),
                            names(tt)))

  sig    <- sum(tt$padj < 0.05, na.rm = TRUE)
  sig_up <- sum(tt$padj < 0.05 & tt$logFC > 0, na.rm = TRUE)
  sig_dn <- sum(tt$padj < 0.05 & tt$logFC < 0, na.rm = TRUE)
  sig_lfc<- sum(tt$padj < 0.05 & abs(tt$logFC) > 0.5, na.rm = TRUE)
  cat(sprintf("===== %s vs Healthy DREAM =====\n", stage_label))
  cat("Total genes tested:", nrow(tt), "\n")
  cat(sprintf("DEGs (padj<0.05): %d (Up: %d, Down: %d)\n",
              sig, sig_up, sig_dn))
  cat(sprintf("DEGs (padj<0.05, |logFC|>0.5): %d\n", sig_lfc))

  out_path <- file.path(RDIR, out_fname)
  fwrite(tt, out_path)
  cat("Saved:", out_path, "\n")
  invisible(tt)
}

# --- Run all three contrasts ---
run_stage_dream("steatosis_target",  "Steatosis",       "dream_results_stage_steatosis.csv")
run_stage_dream("sh_target",         "Steatohepatitis", "dream_results_stage_sh.csv")
run_stage_dream("cirrhosis_target",  "Cirrhosis",       "dream_results_stage_cirrhosis.csv")

cat("\n=== 05c complete:", as.character(Sys.time()), "===\n")
