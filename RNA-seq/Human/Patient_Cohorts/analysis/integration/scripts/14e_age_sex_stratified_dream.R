#!/usr/bin/env Rscript
# 14e_age_sex_stratified_dream.R
# Age + sex stratified dream mega-analyses for Figure 2 panels f-g
# Produces: per_sample_dysregulation.csv, age_stratified_dream.csv,
#           age_sex_interaction_dream.csv, age_continuous_dream.csv
#
# Usage: Rscript 14e_age_sex_stratified_dream.R
# SLURM: bigmem, 16 CPU, 260GB RAM, ~4h

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(yaml)
  library(edgeR)
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

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results")
SIGS <- file.path(RDIR, "disease_signatures")
dir.create(SIGS, showWarnings = FALSE, recursive = TRUE)

cat("=== 14e: Age + Sex Stratified Dream ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# --- Setup parallel backend ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# --- Load merged counts and metadata ---
cat("Loading merged counts...\n")
dge <- readRDS(file.path(RDIR, "integration/merged_dge.rds"))

cat("Loading matched metadata...\n")
meta_matched <- readRDS(file.path(RDIR, "integration/meta_matched.rds"))

# Load unified metadata for age column
meta_unified <- fread(file.path(INT, "metadata/unified_metadata.csv"))
cat("Unified metadata:", nrow(meta_unified), "samples,", ncol(meta_unified), "columns\n")

# Load QC report for pass_technical filter
qc <- fread(file.path(INT, "qc/sample_qc_report.csv"))
pass_samples <- qc[pass_technical == TRUE, sample_id]
cat("QC-passing samples (pass_technical):", length(pass_samples), "\n\n")

# ============================================================
# PART 0: Per-Sample Dysregulation Score (ALL QC-passing samples)
# ============================================================
cat("=== PART 0: Per-Sample Dysregulation Score ===\n")

# Load dream results for top DEGs
dream_file <- file.path(RDIR, "integration/dream_results.csv")
if (!file.exists(dream_file)) stop("dream_results.csv not found — run 05 first")
dream <- fread(dream_file)

# Normalize column names
if ("adj.P.Val" %in% names(dream) && !"padj" %in% names(dream))
  setnames(dream, "adj.P.Val", "padj")

# Top 200 DEGs by absolute t-statistic
dream_sig <- dream[padj < 0.1][order(-abs(t))]
top200 <- head(dream_sig$gene, 200)
cat("Top 200 DEGs selected for dysregulation score\n")

# Compute log2(CPM+1) for all QC-passing samples
all_qc_samples <- intersect(pass_samples, colnames(dge))
dge_all <- dge[, all_qc_samples]

# Get CPM
cpm_mat <- cpm(dge_all, log = FALSE)

# Filter to top 200 genes present in the matrix
top200_present <- intersect(top200, rownames(cpm_mat))
cat("Top DEGs present in count matrix:", length(top200_present), "of 200\n")
cpm_sub <- log2(cpm_mat[top200_present, ] + 1)

# Identify healthy/control samples for baseline
meta_all <- meta_matched[match(all_qc_samples, meta_matched$sample_id)]
ctrl_idx <- which(meta_all$group_binary == "Control")
cat("Control samples for baseline:", length(ctrl_idx), "\n")

# Compute healthy mean and SD per gene
ctrl_mean <- rowMeans(cpm_sub[, ctrl_idx, drop = FALSE], na.rm = TRUE)
ctrl_sd   <- apply(cpm_sub[, ctrl_idx, drop = FALSE], 1, sd, na.rm = TRUE)
# Avoid division by zero
ctrl_sd[ctrl_sd < 0.01] <- 0.01

# Z-score each gene relative to healthy, then mean |z| per sample
z_mat <- sweep(cpm_sub, 1, ctrl_mean, "-")
z_mat <- sweep(z_mat, 1, ctrl_sd, "/")
dysreg_score <- colMeans(abs(z_mat), na.rm = TRUE)

# Build output table
dysreg_dt <- data.table(
  sample_id = all_qc_samples,
  dysreg_score = dysreg_score
)

# Merge metadata
meta_for_merge <- meta_unified[, .(sample_id, dataset, condition, group_binary,
                                    age, sex, nas_score, fibrosis_stage,
                                    diagnosis_harmonized)]
dysreg_dt <- merge(dysreg_dt, meta_for_merge, by = "sample_id", all.x = TRUE)

# Add inferred_sex from meta_matched
sex_map <- meta_matched[, .(sample_id, inferred_sex)]
dysreg_dt <- merge(dysreg_dt, sex_map, by = "sample_id", all.x = TRUE)

fwrite(dysreg_dt, file.path(SIGS, "per_sample_dysregulation.csv"))
cat("Per-sample dysregulation saved:", nrow(dysreg_dt), "samples\n")
cat("  Mean dysreg (Control):", round(mean(dysreg_dt[group_binary == "Control", dysreg_score]), 3), "\n")
cat("  Mean dysreg (Disease):", round(mean(dysreg_dt[group_binary == "Disease", dysreg_score]), 3), "\n\n")

# ============================================================
# PART 1: Age-Stratified Dream (6 datasets with age)
# ============================================================
cat("=== PART 1: Age-Stratified Dream ===\n")

# Strict consistency: restrict age-stratified cohorts to mega-analysis set
# (config/human_datasets.yaml include_in_mega). Drops GSE167523/174478/193066.
# (PRJNA512027 permanently removed from pipeline 2026-05-15.)
ycfg <- yaml::read_yaml(file.path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design", "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
AGE_DATASETS <- intersect(c("GSE130970", "GSE162694", "GSE174478", "GSE193066",
                              "GSE167523"), mega_cohorts)
cat("Age-eligible mega cohorts (k =", length(AGE_DATASETS), "):", paste(AGE_DATASETS, collapse=", "), "\n")

meta_age <- meta_unified[
  dataset %in% AGE_DATASETS &
  !is.na(age) &
  sample_id %in% pass_samples
]
cat("Age-annotated QC-passing samples:", nrow(meta_age), "\n")

# Age groups
meta_age[, age_group := fifelse(age < 50, "Young", "Old")]
age_dist <- meta_age[, .N, by = .(age_group, group_binary)]
cat("\nAge group × disease distribution:\n")
print(age_dist)

# Check we have both groups in both conditions
if (nrow(age_dist) < 4) {
  cat("WARNING: Not all age_group × group_binary combinations present\n")
}

# Subset DGE
age_samples <- intersect(meta_age$sample_id, colnames(dge))
cat("Samples in DGE:", length(age_samples), "\n")
dge_age <- dge[, age_samples]

# Build info dataframe
matched_meta_age <- meta_age[match(colnames(dge_age), meta_age$sample_id)]
matched_sex_age  <- meta_matched$inferred_sex[match(colnames(dge_age), meta_matched$sample_id)]

info_age <- data.frame(
  group_binary = factor(matched_meta_age$group_binary),
  age_group    = factor(matched_meta_age$age_group),
  inferred_sex = factor(matched_sex_age),
  dataset      = factor(matched_meta_age$dataset),
  row.names    = colnames(dge_age),
  stringsAsFactors = FALSE
)
info_age$group_binary <- relevel(info_age$group_binary, ref = "Control")
info_age$age_group    <- relevel(info_age$age_group, ref = "Young")

# Filter by expression
keep_age <- filterByExpr(dge_age, group = info_age$group_binary)
dge_age <- dge_age[keep_age, , keep.lib.sizes = FALSE]
dge_age <- calcNormFactors(dge_age, method = "TMM")
cat("After expression filter:", nrow(dge_age), "genes\n")

# Dream: disease × age_group interaction
form_age <- ~ group_binary * age_group + inferred_sex + (1 | dataset)
cat("\nDream formula:", deparse(form_age), "\n")

cat("Running voomWithDreamWeights...\n")
v_age <- voomWithDreamWeights(dge_age, form_age, info_age, BPPARAM = param)

cat("Running dream...\n")
fit_age <- dream(v_age, form_age, info_age, BPPARAM = param)

# Extract interaction coefficient
int_coef <- "group_binaryDisease:age_groupOld"
if (int_coef %in% colnames(coef(fit_age))) {
  tt_age_int <- topTable(fit_age, coef = int_coef, number = Inf, sort.by = "none")
  tt_age_int$gene <- rownames(tt_age_int)
  tt_age_int$coefficient <- int_coef
  age_int_dt <- as.data.table(tt_age_int)
  setnames(age_int_dt, "adj.P.Val", "padj", skip_absent = TRUE)
  fwrite(age_int_dt, file.path(SIGS, "age_stratified_dream.csv"))
  n_sig <- sum(age_int_dt$padj < 0.05, na.rm = TRUE)
  cat("Age-stratified interaction DEGs (padj<0.05):", n_sig, "\n")
} else {
  cat("WARNING: Interaction coefficient not found. Available:",
      paste(colnames(coef(fit_age)), collapse = ", "), "\n")
  # Save empty result
  fwrite(data.table(), file.path(SIGS, "age_stratified_dream.csv"))
}

# ============================================================
# PART 2: Age × Sex Interaction Dream
# ============================================================
cat("\n=== PART 2: Age x Sex Interaction Dream ===\n")

# Check sufficient representation in all 8 cells (2 disease × 2 age × 2 sex)
cross_tab <- info_age[, c("group_binary", "age_group", "inferred_sex")]
cell_counts <- table(cross_tab$group_binary, cross_tab$age_group, cross_tab$inferred_sex)
cat("Cell counts (disease × age × sex):\n")
print(cell_counts)

min_cell <- min(cell_counts)
cat("Minimum cell count:", min_cell, "\n")

if (min_cell >= 5) {
  form_age_sex <- ~ group_binary * age_group * inferred_sex + (1 | dataset)
  cat("Dream formula:", deparse(form_age_sex), "\n")

  cat("Running voomWithDreamWeights...\n")
  v_age_sex <- voomWithDreamWeights(dge_age, form_age_sex, info_age, BPPARAM = param)

  cat("Running dream...\n")
  fit_age_sex <- dream(v_age_sex, form_age_sex, info_age, BPPARAM = param)

  # Extract 3-way interaction
  int3_coef <- "group_binaryDisease:age_groupOld:inferred_sexMale"
  if (int3_coef %in% colnames(coef(fit_age_sex))) {
    tt_3way <- topTable(fit_age_sex, coef = int3_coef, number = Inf, sort.by = "none")
    tt_3way$gene <- rownames(tt_3way)
    tt_3way$coefficient <- int3_coef
    age_sex_dt <- as.data.table(tt_3way)
    setnames(age_sex_dt, "adj.P.Val", "padj", skip_absent = TRUE)
    fwrite(age_sex_dt, file.path(SIGS, "age_sex_interaction_dream.csv"))
    n_sig3 <- sum(age_sex_dt$padj < 0.05, na.rm = TRUE)
    cat("Age x sex interaction DEGs (padj<0.05):", n_sig3, "\n")
  } else {
    cat("WARNING: 3-way interaction not found. Available:",
        paste(colnames(coef(fit_age_sex)), collapse = ", "), "\n")
    fwrite(data.table(), file.path(SIGS, "age_sex_interaction_dream.csv"))
  }
} else {
  cat("SKIPPING: Insufficient samples in some cells (min =", min_cell, ")\n")
  cat("Writing empty output\n")
  fwrite(data.table(), file.path(SIGS, "age_sex_interaction_dream.csv"))
}

# ============================================================
# PART 3: Continuous Age Dream
# ============================================================
cat("\n=== PART 3: Continuous Age Dream ===\n")

# Add continuous age to info
info_age_cont <- info_age
info_age_cont$age <- meta_age$age[match(rownames(info_age), meta_age$sample_id)]

# Check for NAs
n_na_age <- sum(is.na(info_age_cont$age))
if (n_na_age > 0) {
  cat("WARNING: Dropping", n_na_age, "samples with NA age\n")
  keep_rows <- !is.na(info_age_cont$age)
  info_age_cont <- info_age_cont[keep_rows, , drop = FALSE]
  dge_age_cont <- dge_age[, keep_rows]
} else {
  dge_age_cont <- dge_age
}

# Scale age for numerical stability
info_age_cont$age_scaled <- scale(info_age_cont$age)[, 1]
cat("Age range:", range(info_age_cont$age), "\n")
cat("Scaled age range:", round(range(info_age_cont$age_scaled), 2), "\n")

form_cont <- ~ group_binary + age_scaled + inferred_sex + (1 | dataset)
cat("Dream formula:", deparse(form_cont), "\n")

cat("Running voomWithDreamWeights...\n")
v_cont <- voomWithDreamWeights(dge_age_cont, form_cont, info_age_cont, BPPARAM = param)

cat("Running dream...\n")
fit_cont <- dream(v_cont, form_cont, info_age_cont, BPPARAM = param)

# Extract age coefficient
if ("age_scaled" %in% colnames(coef(fit_cont))) {
  tt_age_cont <- topTable(fit_cont, coef = "age_scaled", number = Inf, sort.by = "none")
  tt_age_cont$gene <- rownames(tt_age_cont)
  tt_age_cont$coefficient <- "age_scaled"
  age_cont_dt <- as.data.table(tt_age_cont)
  setnames(age_cont_dt, "adj.P.Val", "padj", skip_absent = TRUE)
  fwrite(age_cont_dt, file.path(SIGS, "age_continuous_dream.csv"))
  n_sig_cont <- sum(age_cont_dt$padj < 0.05, na.rm = TRUE)
  cat("Continuous age DEGs (padj<0.05):", n_sig_cont, "\n")
} else {
  cat("WARNING: age_scaled coefficient not found\n")
  fwrite(data.table(), file.path(SIGS, "age_continuous_dream.csv"))
}

cat("\n=== 14e completed:", as.character(Sys.time()), "===\n")
