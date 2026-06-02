#!/usr/bin/env Rscript
# 05_diagnostic_5cohort_re_vs_fe.R
# ---------------------------------------------------------------------------
# Diagnostic: with k=5 datasets after enforcing config/human_datasets.yaml
# `include_in_mega: true`, the (1|dataset) random-effect variance is poorly
# identified (5 levels). Compare three models on the same merged DGE:
#
#   M0 = 8-cohort random-effect (current production; loaded from disk)
#   M1 = 5-cohort random-effect: ~ group_binary + inferred_sex + (1|dataset)
#   M2 = 5-cohort fixed-effect : ~ group_binary + inferred_sex + dataset
#
# Compares: DEG counts at standard thresholds; Spearman rho of logFC across
# the shared gene universe; Jaccard / direction concordance among DEGs;
# distribution of group_binaryDisease SE; dataset-RE SD (M1 only).
#
# Output:
#   results/integration/diagnostic_5cohort/diagnostic_summary.txt
#   results/integration/diagnostic_5cohort/per_gene_3way.csv
#   results/integration/diagnostic_5cohort/dream_5cohort_re.csv
#   results/integration/diagnostic_5cohort/dream_5cohort_fe.csv
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

PROJECT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
BASE    <- file.path(PROJECT, "RNA-seq/Human/Patient_Cohorts")
INT     <- file.path(BASE, "analysis/integration")
RDIR    <- file.path(INT, "results/integration")
OUTDIR  <- file.path(RDIR, "diagnostic_5cohort")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# --- Determine the include_in_mega cohort set from yaml ----------------------
yaml_path <- file.path(PROJECT, "config/human_datasets.yaml")
ycfg <- yaml::read_yaml(yaml_path)$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("yaml-declared mega cohorts (k =", length(mega_cohorts), "):",
    paste(mega_cohorts, collapse = ", "), "\n")

# --- Load merged DGE ---------------------------------------------------------
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
cat("Full merged DGE:", ncol(dge), "samples,", nrow(dge), "genes\n")

# --- 8-cohort production reference (M0) --------------------------------------
m0_path <- file.path(RDIR, "dream_results.csv")
m0 <- fread(m0_path)
setnames(m0, "adj.P.Val", "padj", skip_absent = TRUE)
cat("Loaded 8-cohort production DEGs (M0) from", m0_path, ":", nrow(m0), "genes\n")

# --- Subset to 5 yaml-included cohorts --------------------------------------
keep <- dge$samples$dataset %in% mega_cohorts
dge5 <- dge[, keep]
cat("5-cohort DGE:", ncol(dge5), "samples,", nrow(dge5), "genes\n")

meta_new    <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge5), meta_new$sample_id)]
n_na_sex    <- sum(is.na(matched_sex))
if (n_na_sex > 0) warning(n_na_sex, " samples have NA inferred_sex.")

info5 <- data.frame(
  group_binary = factor(dge5$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = droplevels(factor(dge5$samples$dataset)),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE
)
rownames(info5) <- colnames(dge5)

cat("\n5-cohort group x dataset:\n");  print(table(info5$group_binary, info5$dataset))
cat("5-cohort sex:\n");               print(table(info5$inferred_sex, useNA = "always"))

# --- Parallel param ---------------------------------------------------------
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# --- Fit M1: 5-cohort random-effect ----------------------------------------
form_re   <- ~ group_binary + inferred_sex + (1 | dataset)
m1_path   <- file.path(OUTDIR, "dream_5cohort_re.csv")
m1_v_path <- file.path(OUTDIR, "voom_5cohort_re.rds")
if (file.exists(m1_path) && file.exists(m1_v_path)) {
  cat("\n=== M1: 5-cohort RE — loading cached results from", m1_path, "===\n")
  res_re <- fread(m1_path)
  v_re   <- readRDS(m1_v_path)
  cat("Loaded:", nrow(res_re), "genes\n")
} else {
  cat("\n=== M1: 5-cohort RE — formula:", deparse(form_re), "===\n")
  v_re   <- suppressWarnings(voomWithDreamWeights(dge5, form_re, info5, BPPARAM = param))
  fit_re <- suppressWarnings(dream(v_re, form_re, info5, BPPARAM = param))
  # dream() with random effects already moderates via Satterthwaite — do NOT eBayes
  res_re <- topTable(fit_re, coef = "group_binaryDisease", number = Inf, sort.by = "none")
  res_re$gene <- rownames(res_re)
  res_re <- as.data.table(res_re)
  setnames(res_re, "adj.P.Val", "padj")
  fwrite(res_re, m1_path)
  saveRDS(v_re, m1_v_path)
  cat("Saved 5-cohort RE results:", nrow(res_re), "genes\n")
}

# --- Fit M2: 5-cohort fixed-effect dataset --------------------------------
form_fe <- ~ group_binary + inferred_sex + dataset
cat("\n=== M2: 5-cohort FE — formula:", deparse(form_fe), "===\n")
v_fe   <- suppressWarnings(voomWithDreamWeights(dge5, form_fe, info5, BPPARAM = param))
fit_fe <- suppressWarnings(dream(v_fe, form_fe, info5, BPPARAM = param))
# dream() with NO random terms returns a plain MArrayLM (lmFit equivalent);
# eBayes() is required for moderated t-stats and topTable()-ready columns.
fit_fe <- eBayes(fit_fe)
res_fe <- topTable(fit_fe, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res_fe$gene <- rownames(res_fe)
res_fe <- as.data.table(res_fe)
setnames(res_fe, "adj.P.Val", "padj")
fwrite(res_fe, file.path(OUTDIR, "dream_5cohort_fe.csv"))
cat("Saved 5-cohort FE results:", nrow(res_fe), "genes\n")

# --- Build 3-way per-gene comparison table --------------------------------
m0_sub <- m0[, .(gene, logFC_M0 = logFC, padj_M0 = padj)]
m1_sub <- res_re[, .(gene, logFC_M1 = logFC, padj_M1 = padj)]
m2_sub <- res_fe[, .(gene, logFC_M2 = logFC, padj_M2 = padj)]
three  <- merge(merge(m0_sub, m1_sub, by = "gene"), m2_sub, by = "gene")
fwrite(three, file.path(OUTDIR, "per_gene_3way.csv"))
cat("\n3-way merge:", nrow(three), "genes\n")

# --- Comparison metrics ---------------------------------------------------
deg_set <- function(dt, lfc_col, padj_col, padj_thr = 0.05, lfc_thr = 0.3) {
  dt[ get(padj_col) < padj_thr & abs(get(lfc_col)) > lfc_thr, gene ]
}

s_M0 <- deg_set(three, "logFC_M0", "padj_M0")
s_M1 <- deg_set(three, "logFC_M1", "padj_M1")
s_M2 <- deg_set(three, "logFC_M2", "padj_M2")

jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))
dir_conc <- function(a, b, dt, ca, cb) {
  shared <- intersect(a, b)
  if (!length(shared)) return(NA_real_)
  sub <- dt[gene %in% shared]
  mean(sign(sub[[ca]]) == sign(sub[[cb]]))
}

rho_01 <- cor(three$logFC_M0, three$logFC_M1, method = "spearman")
rho_02 <- cor(three$logFC_M0, three$logFC_M2, method = "spearman")
rho_12 <- cor(three$logFC_M1, three$logFC_M2, method = "spearman")

dc_01 <- dir_conc(s_M0, s_M1, three, "logFC_M0", "logFC_M1")
dc_02 <- dir_conc(s_M0, s_M2, three, "logFC_M0", "logFC_M2")
dc_12 <- dir_conc(s_M1, s_M2, three, "logFC_M1", "logFC_M2")

# SE distribution comparison (5-cohort RE vs FE)
res_re[, se := abs(logFC) / abs(t)]
res_fe[, se := abs(logFC) / abs(t)]
se_re_med <- median(res_re$se, na.rm = TRUE)
se_fe_med <- median(res_fe$se, na.rm = TRUE)
se_re_mean <- mean(res_re$se,  na.rm = TRUE)
se_fe_mean <- mean(res_fe$se,  na.rm = TRUE)

# Dataset RE variance summary (M1)
re_var_summary <- tryCatch({
  vp <- as.data.frame(VarPartition <- variancePartition::fitExtractVarPartModel(
    v_re, form_re, info5, BPPARAM = param))
  c(med_dataset = median(vp$dataset, na.rm = TRUE),
    mean_dataset = mean(vp$dataset, na.rm = TRUE),
    pct_genes_high_RE = mean(vp$dataset > 0.20, na.rm = TRUE))
}, error = function(e) {
  cat("VarPart fit failed:", conditionMessage(e), "\n"); rep(NA_real_, 3)
})

# --- Write summary --------------------------------------------------------
summary_path <- file.path(OUTDIR, "diagnostic_summary.txt")
sink(summary_path)
cat("============================================================\n")
cat("Dream mega-analysis diagnostic — 8-cohort RE vs 5-cohort RE/FE\n")
cat("Generated:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

cat("Cohort sets:\n")
cat("  M0 (8-cohort RE, production): GSE126848, GSE130970, GSE135251, GSE162694,\n")
cat("                                 GSE174478, GSE193066, GSE213621, GSE240729\n")
cat("  M1/M2 (5-cohort, yaml-include_in_mega):", paste(mega_cohorts, collapse = ", "), "\n\n")

cat("Sample counts:\n")
cat(sprintf("  M0: %d samples (loaded from disk)\n", NA_integer_))  # not loaded here
cat(sprintf("  M1/M2: %d samples\n", ncol(dge5)))

cat("\nGroup × dataset (5-cohort):\n")
print(table(info5$group_binary, info5$dataset))

cat("\n----- DEG counts (padj < 0.05, |logFC| > 0.3) -----\n")
cat(sprintf("  M0  (8-cohort RE)     : %d DEGs\n", length(s_M0)))
cat(sprintf("  M1  (5-cohort RE)     : %d DEGs\n", length(s_M1)))
cat(sprintf("  M2  (5-cohort FE)     : %d DEGs\n", length(s_M2)))

cat("\n----- Spearman rho of logFC (across all", nrow(three), "shared genes) -----\n")
cat(sprintf("  rho(M0, M1)  : %.4f\n", rho_01))
cat(sprintf("  rho(M0, M2)  : %.4f\n", rho_02))
cat(sprintf("  rho(M1, M2)  : %.4f\n", rho_12))

cat("\n----- Jaccard of DEG sets -----\n")
cat(sprintf("  J(M0, M1)    : %.4f  (intersect %d / union %d)\n",
            jacc(s_M0, s_M1), length(intersect(s_M0, s_M1)), length(union(s_M0, s_M1))))
cat(sprintf("  J(M0, M2)    : %.4f  (intersect %d / union %d)\n",
            jacc(s_M0, s_M2), length(intersect(s_M0, s_M2)), length(union(s_M0, s_M2))))
cat(sprintf("  J(M1, M2)    : %.4f  (intersect %d / union %d)\n",
            jacc(s_M1, s_M2), length(intersect(s_M1, s_M2)), length(union(s_M1, s_M2))))

cat("\n----- Direction concordance among shared DEGs -----\n")
cat(sprintf("  M0 ∩ M1 sign agreement: %.4f\n", dc_01))
cat(sprintf("  M0 ∩ M2 sign agreement: %.4f\n", dc_02))
cat(sprintf("  M1 ∩ M2 sign agreement: %.4f\n", dc_12))

cat("\n----- group_binaryDisease SE distribution (5-cohort) -----\n")
cat(sprintf("  M1 RE  median SE : %.4f   mean SE : %.4f\n", se_re_med, se_re_mean))
cat(sprintf("  M2 FE  median SE : %.4f   mean SE : %.4f\n", se_fe_med, se_fe_mean))
cat(sprintf("  ratio RE/FE median: %.4f  (>>1 indicates RE inflation due to k=5 levels)\n",
            se_re_med / se_fe_med))

cat("\n----- Dataset RE variance fraction (M1, fitExtractVarPartModel) -----\n")
print(re_var_summary)

cat("\n----- Interpretation guide -----\n")
cat("  - rho(M1,M2) > 0.95  -> RE and FE give equivalent ranking; either is fine.\n")
cat("  - rho(M1,M2) < 0.90  -> RE form materially differs from FE; investigate.\n")
cat("  - SE ratio RE/FE within 1.0-1.2 -> RE is well-identified at k=5.\n")
cat("  - SE ratio RE/FE >> 1.5         -> RE absorbs too much variance; prefer FE.\n")
cat("  - rho(M0,M1) >> rho(M0,M2)      -> 8-cohort production was tracking dataset RE,\n")
cat("                                     not within-cohort biology; FE is the safer call.\n")
sink()

cat("\nWrote:", summary_path, "\n")
cat("Done.\n")
