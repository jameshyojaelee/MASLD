#!/usr/bin/env Rscript
# 14.4f_composition_adjusted_T13.R
# ---------------------------------------------------------------------------
# A8 / T13: Composition-adjusted sex-interaction dream refit.
#
# C1 critique #1: Script 26 explicitly excluded deconvolution covariates,
# justifying that "deconvolution covariates collapse DEG counts to near-zero."
# That collapse is the diagnostic finding itself: a large fraction of the
# 3,080 sex_dimorphic disease signal is likely cell-composition-driven
# (macrophage / Kupffer / NK fraction shifts in female MASLD livers), not
# cell-intrinsic transcriptional dimorphism.
#
# Adversarial test: refit the canonical dream interaction model with
# progressively stricter composition covariates and report
#   M0 (canonical):     ~ group_binary * inferred_sex + (1|dataset)
#   M1 (hep-adjusted):  ~ group_binary * inferred_sex + hep_fraction + (1|dataset)
#   M2 (fully adjusted):~ group_binary * inferred_sex + hep_fraction
#                          + macrophage_fraction + endothelial_fraction
#                          + cholangiocyte_fraction + (1|dataset)
#
# For each model:
#   - Count n_sex_dimorphic (interaction padj < 0.05)
#   - Recompute sex_class per gene (Female_biased / Male_biased / Divergent /
#     Concordant) using same rules as script 26
#   - Compare to M0: how many stable / lost / new
#
# Headline: subset of M0-sex_dimorphic genes that REMAIN sex_dimorphic in M2.
# These are the composition-INDEPENDENT sex-DEGs and the only ones we can
# claim are cell-intrinsic transcriptional dimorphism.
#
# Outputs (under RNA-seq/results/audit_sensitivity/sex_composition_adjusted/):
#   composition_adjusted_sex_classification.csv  — per gene, all 3 sex_class
#   summary.csv                                  — per model count summary
#   composition_independent_genes.csv            — headline list
#   interaction_lfc_M0.csv, _M1.csv, _M2.csv     — raw interaction LFC tables
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
  library(yaml)
})

# Force injection into lme4 namespace BEFORE loading variancePartition
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

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE   <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR   <- file.path(INT, "results/integration")
PROP_F <- file.path(BASE, "RNA-seq/results/celltype_attribution",
                    "persample_celltype_proportions.csv")
OUT    <- file.path(BASE, "RNA-seq/results/audit_sensitivity",
                    "sex_composition_adjusted")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat("=== A8 / T13: Composition-adjusted sex-interaction refit ===\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Output dir:", OUT, "\n")

# ---------------------------------------------------------------------------
# Parallel setup
# ---------------------------------------------------------------------------
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# ---------------------------------------------------------------------------
# Load DGE + metadata (replicates Script 26 lines 84-141)
# ---------------------------------------------------------------------------
cat("\n===== Loading data =====\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("yaml-declared mega cohorts (k =", length(mega_cohorts), "):",
    paste(mega_cohorts, collapse = ", "), "\n")

keep_samples <- dge$samples$dataset %in% mega_cohorts
dge_mega <- dge[, keep_samples]
cat("Samples for mega-analysis:", ncol(dge_mega), "\n")

qc_report <- fread(file.path(INT, "qc/sample_qc_report.csv"))
sex_pass_ids <- qc_report[pass_sex == TRUE, sample_id]
sex_fail <- !colnames(dge_mega) %in% sex_pass_ids
if (any(sex_fail)) {
  cat("Removing", sum(sex_fail), "samples that failed sex check\n")
  dge_mega <- dge_mega[, !sex_fail]
}

meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]

# ---------------------------------------------------------------------------
# Load cell-type proportions
# ---------------------------------------------------------------------------
cat("\n===== Loading cell-type proportions =====\n")
prop <- fread(PROP_F)
cat("Proportions table:", nrow(prop), "samples,", ncol(prop), "cols\n")

# Keep only relevant columns and rename to safe identifiers
prop_use <- prop[, .(
  sample_id   = sample_id,
  hep_frac    = `Hepatocytes`,
  mac_frac    = `Macrophages`,
  endo_frac   = `Endothelial cells`,
  chol_frac   = `Cholangiocytes`
)]
cat("hep NA:", sum(is.na(prop_use$hep_frac)),
    " mac NA:", sum(is.na(prop_use$mac_frac)),
    " endo NA:", sum(is.na(prop_use$endo_frac)),
    " chol NA:", sum(is.na(prop_use$chol_frac)), "\n")

# Replace NA with 0 (BayesPrism / MuSiC unmeasured cell types absent in panel)
# Conservative: do NOT impute — drop samples missing any of the 4 fractions
n_before <- nrow(prop_use)
prop_use <- prop_use[!is.na(hep_frac) & !is.na(mac_frac) &
                     !is.na(endo_frac) & !is.na(chol_frac)]
cat("Samples with all 4 fractions:", nrow(prop_use), "/", n_before, "\n")

# ---------------------------------------------------------------------------
# Build info table — base
# ---------------------------------------------------------------------------
info0 <- data.frame(
  sample_id    = colnames(dge_mega),
  group_binary = factor(dge_mega$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = factor(dge_mega$samples$dataset),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE
)
rownames(info0) <- info0$sample_id

# Merge with fractions; drop samples missing any covariate (sex or fractions)
info_mrg <- merge(info0, as.data.frame(prop_use), by = "sample_id", all.x = TRUE)
rownames(info_mrg) <- info_mrg$sample_id

has_all <- !is.na(info_mrg$inferred_sex) &
           !is.na(info_mrg$hep_frac) & !is.na(info_mrg$mac_frac) &
           !is.na(info_mrg$endo_frac) & !is.na(info_mrg$chol_frac)
cat("\nCovariate completeness:\n")
cat("  with sex      :", sum(!is.na(info_mrg$inferred_sex)), "\n")
cat("  with hep_frac :", sum(!is.na(info_mrg$hep_frac)), "\n")
cat("  with all 4    :", sum(has_all), "\n")

info_mrg <- info_mrg[has_all, , drop = FALSE]
dge_mega <- dge_mega[, colnames(dge_mega) %in% info_mrg$sample_id]
# Ensure same order
info_mrg <- info_mrg[match(colnames(dge_mega), info_mrg$sample_id), , drop = FALSE]
rownames(info_mrg) <- info_mrg$sample_id
stopifnot(identical(rownames(info_mrg), colnames(dge_mega)))

cat("\nFinal sample count:", ncol(dge_mega), "\n")
cat("Group x sex:\n"); print(table(info_mrg$group_binary, info_mrg$inferred_sex))
cat("Group x dataset:\n"); print(table(info_mrg$group_binary, info_mrg$dataset))

# Standardize composition covariates to zero-mean unit-SD for numerical
# stability in dream
for (cv in c("hep_frac", "mac_frac", "endo_frac", "chol_frac")) {
  v <- info_mrg[[cv]]
  info_mrg[[cv]] <- as.numeric(scale(v))
}
cat("\nStandardized composition covariates (mean ~0, sd ~1).\n")

# ---------------------------------------------------------------------------
# Model fitting helper
# ---------------------------------------------------------------------------
fit_model <- function(model_label, formula, dge_in, info_in, param_in) {
  cat("\n===== Fitting", model_label, "=====\n")
  cat("Formula:", deparse(formula), "\n")
  cat("N samples:", ncol(dge_in), "  N genes:", nrow(dge_in), "\n")

  cat("  voomWithDreamWeights...\n")
  v <- suppressWarnings(voomWithDreamWeights(dge_in, formula, info_in, BPPARAM = param_in))

  cat("  dream...\n")
  fit <- suppressWarnings(dream(v, formula, info_in, BPPARAM = param_in))

  all_coefs <- colnames(fit$coefficients)
  cat("  Available coefficients:", paste(all_coefs, collapse = ", "), "\n")
  int_coef <- grep("group_binary.*inferred_sex|inferred_sex.*group_binary",
                   all_coefs, value = TRUE)
  if (length(int_coef) == 0) stop("No interaction coefficient found.")
  int_coef <- int_coef[1]
  cat("  Interaction coef:", int_coef, "\n")

  res <- topTable(fit, coef = int_coef, number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  res_dt  <- as.data.table(res)
  setnames(res_dt, "adj.P.Val", "padj")

  n_dim <- sum(res_dt$padj < 0.05, na.rm = TRUE)
  cat("  n sex_dimorphic (interaction padj<0.05):", n_dim, "\n")

  list(table = res_dt, model = model_label, formula = deparse(formula),
       int_coef = int_coef)
}

# ---------------------------------------------------------------------------
# Per-stratum dream (replicates Script 26 PART 2) for sex_class assignment
# ---------------------------------------------------------------------------
run_stratum_dream <- function(sex_label, dge_full, info_full, param_in,
                              extra_cov = character(0)) {
  cat("\n--- Stratum:", sex_label,
      "; extra covariates:", paste(extra_cov, collapse=","), "---\n")

  info_s <- info_full[info_full$inferred_sex == sex_label, , drop = FALSE]
  dge_s  <- dge_full[, colnames(dge_full) %in% rownames(info_s)]

  # Drop single-group datasets
  ds_counts <- table(info_s$dataset, info_s$group_binary)
  keep_ds <- rownames(ds_counts)[apply(ds_counts, 1, function(x) all(x > 0))]
  if (length(keep_ds) < nrow(ds_counts)) {
    info_s <- info_s[info_s$dataset %in% keep_ds, , drop = FALSE]
    info_s$dataset <- droplevels(info_s$dataset)
    dge_s  <- dge_s[, colnames(dge_s) %in% rownames(info_s)]
  }
  if (nlevels(info_s$dataset) < 2) return(NULL)

  rhs <- "group_binary"
  if (length(extra_cov) > 0) rhs <- paste(c(rhs, extra_cov), collapse = " + ")
  rhs <- paste0(rhs, " + (1|dataset)")
  form <- as.formula(paste0("~ ", rhs))
  cat("  Formula:", deparse(form), "\n")

  v_s   <- suppressWarnings(voomWithDreamWeights(dge_s, form, info_s, BPPARAM = param_in))
  fit_s <- suppressWarnings(dream(v_s, form, info_s, BPPARAM = param_in))
  res_s <- topTable(fit_s, coef = "group_binaryDisease", number = Inf, sort.by = "none")
  res_s$gene <- rownames(res_s)
  res_s_dt <- as.data.table(res_s)
  setnames(res_s_dt, "adj.P.Val", "padj")
  res_s_dt
}

# ---------------------------------------------------------------------------
# Matched-dataset restriction (same as Script 26)
# ---------------------------------------------------------------------------
sex_levels <- levels(info_mrg$inferred_sex)
male_label   <- sex_levels[grepl("^M", sex_levels)][1]
female_label <- sex_levels[grepl("^F", sex_levels)][1]
cat("Sex levels:", paste(sex_levels, collapse = ", "), "\n")
cat("male=", male_label, " female=", female_label, "\n")

ds_m_ctrl <- unique(as.character(info_mrg$dataset[info_mrg$group_binary == "Control" & info_mrg$inferred_sex == male_label]))
ds_f_ctrl <- unique(as.character(info_mrg$dataset[info_mrg$group_binary == "Control" & info_mrg$inferred_sex == female_label]))
matched_ds <- intersect(ds_m_ctrl, ds_f_ctrl)
keep_matched <- info_mrg$dataset %in% matched_ds
dge_matched  <- dge_mega[, keep_matched]
info_matched <- info_mrg[keep_matched, , drop = FALSE]
info_matched$dataset <- droplevels(info_matched$dataset)
rownames(info_matched) <- info_matched$sample_id
cat("Matched samples:", ncol(dge_matched), "\n")

# ---------------------------------------------------------------------------
# Classification helper — mirror Script 26 PART 3 exactly
# ---------------------------------------------------------------------------
PADJ_THRESH       <- 0.05
LFC_THRESH        <- 0.25
INT_PADJ_THRESH   <- 0.05
DIVERGENT_MIN_LFC <- 0.1

classify_sex <- function(res_int, res_male, res_female) {
  cls <- merge(
    res_male[,   .(gene, logFC_M = logFC, padj_M = padj)],
    res_female[, .(gene, logFC_F = logFC, padj_F = padj)],
    by = "gene", all = TRUE
  )
  cls <- merge(cls,
    res_int[, .(gene, interaction_logFC = logFC, interaction_padj = padj)],
    by = "gene", all.x = TRUE
  )
  cls[, same_sign := !is.na(logFC_M) & !is.na(logFC_F) & sign(logFC_M) == sign(logFC_F)]
  cls[, sex_dimorphic := !is.na(interaction_padj) & interaction_padj < INT_PADJ_THRESH]

  cls[, sex_class := fcase(
    !sex_dimorphic, "Concordant",
    sex_dimorphic & !is.na(logFC_F) & !is.na(logFC_M) &
      sign(logFC_F) != sign(logFC_M) &
      abs(logFC_F) > DIVERGENT_MIN_LFC & abs(logFC_M) > DIVERGENT_MIN_LFC,
      "Divergent",
    sex_dimorphic & !is.na(logFC_F) & !is.na(logFC_M) &
      abs(logFC_F) >= abs(logFC_M), "Female_biased",
    sex_dimorphic & !is.na(logFC_F) & !is.na(logFC_M) &
      abs(logFC_M) > abs(logFC_F), "Male_biased",
    default = "Unclassified"
  )]
  cls
}

# ===========================================================================
# RUN THREE MODELS
# ===========================================================================
form_M0 <- ~ group_binary * inferred_sex + (1|dataset)
form_M1 <- ~ group_binary * inferred_sex + hep_frac + (1|dataset)
form_M2 <- ~ group_binary * inferred_sex + hep_frac + mac_frac + endo_frac + chol_frac + (1|dataset)

M0 <- fit_model("M0_canonical",    form_M0, dge_mega, info_mrg, param)
M1 <- fit_model("M1_hep_adjusted", form_M1, dge_mega, info_mrg, param)
M2 <- fit_model("M2_fully_adjusted", form_M2, dge_mega, info_mrg, param)

fwrite(M0$table, file.path(OUT, "interaction_lfc_M0.csv"))
fwrite(M1$table, file.path(OUT, "interaction_lfc_M1.csv"))
fwrite(M2$table, file.path(OUT, "interaction_lfc_M2.csv"))

# ---------------------------------------------------------------------------
# Stratified runs per model (for sex_class assignment)
# ---------------------------------------------------------------------------
cat("\n===== Stratified runs (M0) =====\n")
res_M_M0 <- run_stratum_dream(male_label,   dge_matched, info_matched, param)
res_F_M0 <- run_stratum_dream(female_label, dge_matched, info_matched, param)

cat("\n===== Stratified runs (M1 with hep_frac) =====\n")
res_M_M1 <- run_stratum_dream(male_label,   dge_matched, info_matched, param,
                              extra_cov = "hep_frac")
res_F_M1 <- run_stratum_dream(female_label, dge_matched, info_matched, param,
                              extra_cov = "hep_frac")

cat("\n===== Stratified runs (M2 with all 4 fractions) =====\n")
extras <- c("hep_frac", "mac_frac", "endo_frac", "chol_frac")
res_M_M2 <- run_stratum_dream(male_label,   dge_matched, info_matched, param,
                              extra_cov = extras)
res_F_M2 <- run_stratum_dream(female_label, dge_matched, info_matched, param,
                              extra_cov = extras)

cls_M0 <- classify_sex(M0$table, res_M_M0, res_F_M0)
cls_M1 <- classify_sex(M1$table, res_M_M1, res_F_M1)
cls_M2 <- classify_sex(M2$table, res_M_M2, res_F_M2)

# ---------------------------------------------------------------------------
# Build composition_adjusted_sex_classification.csv
# ---------------------------------------------------------------------------
out_cls <- merge(
  cls_M0[, .(gene,
             sex_class_M0       = sex_class,
             interaction_padj_M0 = interaction_padj,
             interaction_lfc_M0  = interaction_logFC,
             logFC_F_M0 = logFC_F, logFC_M_M0 = logFC_M)],
  cls_M1[, .(gene,
             sex_class_M1       = sex_class,
             interaction_padj_M1 = interaction_padj,
             interaction_lfc_M1  = interaction_logFC)],
  by = "gene", all = TRUE
)
out_cls <- merge(out_cls,
  cls_M2[, .(gene,
             sex_class_M2       = sex_class,
             interaction_padj_M2 = interaction_padj,
             interaction_lfc_M2  = interaction_logFC)],
  by = "gene", all = TRUE
)
out_cls[, retained_in_M2 := !is.na(sex_class_M0) & sex_class_M0 != "Concordant" &
                            !is.na(sex_class_M2) & sex_class_M2 != "Concordant"]
out_cls[, lost_in_M2     := !is.na(sex_class_M0) & sex_class_M0 != "Concordant" &
                            (is.na(sex_class_M2) | sex_class_M2 == "Concordant")]
out_cls[, new_in_M2      := (is.na(sex_class_M0) | sex_class_M0 == "Concordant") &
                            !is.na(sex_class_M2) & sex_class_M2 != "Concordant"]
out_cls[, sex_class_change := !is.na(sex_class_M0) & !is.na(sex_class_M2) &
                              sex_class_M0 != "Concordant" & sex_class_M2 != "Concordant" &
                              sex_class_M0 != sex_class_M2]

fwrite(out_cls, file.path(OUT, "composition_adjusted_sex_classification.csv"))
cat("\nSaved composition_adjusted_sex_classification.csv (",
    nrow(out_cls), "genes)\n")

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
summarize_model <- function(cls, label) {
  data.table(
    model            = label,
    n_genes_tested   = nrow(cls),
    n_sex_dimorphic  = sum(cls$sex_dimorphic, na.rm = TRUE),
    n_female_biased  = sum(cls$sex_class == "Female_biased", na.rm = TRUE),
    n_male_biased    = sum(cls$sex_class == "Male_biased", na.rm = TRUE),
    n_divergent      = sum(cls$sex_class == "Divergent", na.rm = TRUE),
    n_unclassified   = sum(cls$sex_class == "Unclassified", na.rm = TRUE),
    n_concordant     = sum(cls$sex_class == "Concordant", na.rm = TRUE)
  )
}

summary_dt <- rbindlist(list(
  summarize_model(cls_M0, "M0_canonical"),
  summarize_model(cls_M1, "M1_hep_adjusted"),
  summarize_model(cls_M2, "M2_fully_adjusted")
))
fwrite(summary_dt, file.path(OUT, "summary.csv"))
cat("\n===== SUMMARY =====\n"); print(summary_dt)

# ---------------------------------------------------------------------------
# Composition-INDEPENDENT genes (HEADLINE)
# Subset of M0 sex_dimorphic that REMAIN sex_dimorphic in M2 with
# coherent sex_class.
# ---------------------------------------------------------------------------
indep <- out_cls[
  !is.na(sex_class_M0) & sex_class_M0 %in% c("Female_biased","Male_biased","Divergent") &
  !is.na(sex_class_M2) & sex_class_M2 %in% c("Female_biased","Male_biased","Divergent")
]

# Mark whether the class call agrees between M0 and M2
indep[, sex_class_stable := sex_class_M0 == sex_class_M2]
indep_out <- indep[, .(gene,
                       sex_class_M0, sex_class_M1, sex_class_M2,
                       interaction_padj_M0, interaction_padj_M1, interaction_padj_M2,
                       interaction_lfc_M0, interaction_lfc_M1, interaction_lfc_M2,
                       logFC_F_M0, logFC_M_M0,
                       sex_class_stable)][
                       order(interaction_padj_M2)]

fwrite(indep_out, file.path(OUT, "composition_independent_genes.csv"))
cat("\n===== HEADLINE: Composition-independent sex-DEGs =====\n")
cat("Total M0 sex_dimorphic genes:    ", nrow(out_cls[!is.na(sex_class_M0) & sex_class_M0 != "Concordant"]), "\n")
cat("Composition-independent (M0 & M2):", nrow(indep_out), "\n")
cat("  ...with same sex_class in both: ", sum(indep_out$sex_class_stable), "\n")
cat("\nBreakdown by M2 sex_class:\n")
print(table(indep_out$sex_class_M2))
cat("\nM0->M2 class transitions (composition-independent subset):\n")
print(table(indep_out$sex_class_M0, indep_out$sex_class_M2,
            dnn = c("M0", "M2")))

cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
