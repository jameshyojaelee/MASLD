#!/usr/bin/env Rscript
# 15_chen_downsampling.R
# ---------------------------------------------------------------------------
# ONE iteration of Chen cohort downsampling sensitivity analysis.
#
# Env vars:
#   DOWNSAMPLE_SIZE  — target N for Chen (default: 118 = median of other cohorts)
#   DOWNSAMPLE_ITER  — iteration number >= 1 (default: 1)
#
# For each iteration:
#   1. Stratified subsample Chen (preserving Disease:Control ratio)
#   2. Run dream mega-analysis with the same formula as script 05
#   3. Save results CSV
#
# Output: .../audit_sensitivity/chen_influence/dream_downsample_i{II}.csv
# ---------------------------------------------------------------------------

t0 <- proc.time()

# ---------------------------------------------------------------------------
# Library loading + reformulas injection (MUST match script 05 / 14.4)
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT     <- file.path(BASE, "analysis/integration")
RDIR    <- file.path(INT, "results/integration")
OUT_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/chen_influence"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Parallel setup
# ---------------------------------------------------------------------------
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# ---------------------------------------------------------------------------
# Read env vars
# ---------------------------------------------------------------------------
DS_SIZE <- as.integer(Sys.getenv("DOWNSAMPLE_SIZE", "118"))
ITER    <- as.integer(Sys.getenv("DOWNSAMPLE_ITER", "1"))
SEED    <- 42000L + ITER

cat("DOWNSAMPLE_SIZE:", DS_SIZE, "\n")
cat("DOWNSAMPLE_ITER:", ITER, "\n")
cat("SEED:", SEED, "\n")

iter_tag <- sprintf("%02d", ITER)
out_file <- file.path(OUT_DIR, paste0("dream_downsample_i", iter_tag, ".csv"))

# ============================================================================
# Load data (same as script 05)
# ============================================================================
cat("\n===== Loading data =====\n")

dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

# Exclude datasets with irrecoverable confounds (MUST match script 05)
ycfg <- yaml::read_yaml(file.path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design", "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
keep_samples <- dge$samples$dataset %in% mega_cohorts
dge_mega <- dge[, keep_samples]
cat("Mega cohorts (yaml):", paste(mega_cohorts, collapse = ", "), "
")
cat("Samples before downsampling:", ncol(dge_mega), "\n")

# Load metadata with inferred sex
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]

# Build info data.frame
info <- data.frame(
  sample_id    = colnames(dge_mega),
  group_binary = factor(dge_mega$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = factor(dge_mega$samples$dataset),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE
)
rownames(info) <- colnames(dge_mega)

# Drop samples with NA sex
na_sex <- is.na(info$inferred_sex)
if (any(na_sex)) {
  cat("Dropping", sum(na_sex), "samples with NA sex\n")
  dge_mega <- dge_mega[, !na_sex]
  info     <- info[!na_sex, , drop = FALSE]
}

cat("\nPre-downsample distribution:\n")
print(table(info$dataset, info$group_binary))

# ============================================================================
# Stratified downsampling of Chen (GSE213621)
# ============================================================================
cat("\n===== Stratified downsampling of GSE213621 =====\n")

set.seed(SEED)

chen_idx <- which(info$dataset == "GSE213621")
chen_info <- info[chen_idx, ]
n_chen <- nrow(chen_info)
cat(sprintf("Chen samples: %d (Control: %d, Disease: %d)\n",
            n_chen,
            sum(chen_info$group_binary == "Control"),
            sum(chen_info$group_binary == "Disease")))

# Stratified sample preserving group_binary ratio
chen_ctrl <- which(chen_info$group_binary == "Control")
chen_dis  <- which(chen_info$group_binary == "Disease")

frac <- DS_SIZE / n_chen
n_ctrl_draw <- max(round(length(chen_ctrl) * frac), 10)  # at least 10 controls
n_dis_draw  <- DS_SIZE - n_ctrl_draw

cat(sprintf("Drawing %d samples (Control: %d, Disease: %d)\n",
            DS_SIZE, n_ctrl_draw, n_dis_draw))

drawn_ctrl <- sample(chen_ctrl, n_ctrl_draw)
drawn_dis  <- sample(chen_dis, n_dis_draw)
drawn_chen_local <- c(drawn_ctrl, drawn_dis)
drawn_chen_global <- chen_idx[drawn_chen_local]

# Keep all non-Chen samples + downsampled Chen
non_chen_idx <- which(info$dataset != "GSE213621")
keep_idx <- sort(c(non_chen_idx, drawn_chen_global))

dge_ds <- dge_mega[, keep_idx]
info_ds <- info[keep_idx, , drop = FALSE]
rownames(info_ds) <- colnames(dge_ds)

cat(sprintf("Samples after downsampling: %d\n", ncol(dge_ds)))
cat("\nPost-downsample distribution:\n")
print(table(info_ds$dataset, info_ds$group_binary))

# ============================================================================
# Dream mega-analysis (same formula as script 05)
# ============================================================================
cat("\n===== Running dream =====\n")

form <- ~ group_binary + inferred_sex + (1|dataset)
cat("Formula:", deparse(form), "\n")

cat("Running voomWithDreamWeights...\n")
v <- suppressWarnings(voomWithDreamWeights(dge_ds, form, info_ds, BPPARAM = param))

cat("Running dream()...\n")
fit <- suppressWarnings(dream(v, form, info_ds, BPPARAM = param))

# Extract results for Disease vs Control
res <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res_dt <- as.data.table(res, keep.rownames = "gene")

cat(sprintf("DEGs (padj<0.1): %d / %d\n",
            sum(res_dt$adj.P.Val < 0.1, na.rm = TRUE), nrow(res_dt)))

# ============================================================================
# Save
# ============================================================================
fwrite(res_dt, out_file)
cat("\nSaved:", out_file, "\n")

elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("Done in %.1f minutes\n", elapsed / 60))
