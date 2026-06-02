#!/usr/bin/env Rscript
# 14.4e_sex_bootstrap_T4.R
# ---------------------------------------------------------------------------
# T4 (Team A audit): Sample-level bootstrap CI on interaction LFC.
#
# Strategy: resample samples WITH REPLACEMENT, stratified by
# dataset × group_binary × inferred_sex (each stratum bootstrapped to its
# original size). Refit the dream interaction model
#   ~ group_binary * inferred_sex + (1|dataset)
# on each bootstrap replicate. Collect per-gene interaction logFC + adj.P.Val.
#
# Aggregate across B = N_BOOT replicates:
#   - bootstrap mean interaction logFC
#   - 95% percentile CI (2.5%, 97.5%)
#   - fraction of replicates with padj < 0.05  (call stability)
#   - SD across replicates
#
# Implementation note: duplicated sample names are de-duplicated by suffixing
# (sample_id__b1, __b2, ...) before passing to voomWithDreamWeights so the
# random-effect grouping factor remains well-defined. Since cell-wise replicates
# are identical libraries, duplicates inherit identical group/sex/dataset labels
# but get unique column names so the dream weight machinery doesn't collapse.
#
# Env vars:
#   BOOT_ITER   — integer >=1 (default 1)   array task ID
#   N_BOOT      — integer >=1 (default 5)   bootstrap replicates per array task
#   R_PARALLEL_SEED — integer (default 42)  base seed
#
# To produce B=200 total replicates with --array=1-40, set N_BOOT=5 (5×40=200).
#
# Outputs (one per array task):
#   .../audit_sensitivity/sex_bootstrap/iterations/
#     boot_int_t{II}.csv     — gene + per-replicate logFC + per-replicate padj
#                              wide format: logFC_b1..bN, padj_b1..bN
#
# Aggregation (separate script 14.4e_aggregate_bootstrap_T4.R) reads all
# task csvs and produces interaction_lfc_ci.csv.
# ---------------------------------------------------------------------------

t0 <- proc.time()

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(yaml)
  library(edgeR)
})

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
OUT_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_bootstrap/iterations"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Parallel + seed setup
# ---------------------------------------------------------------------------
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
seed_R_par <- as.integer(Sys.getenv("R_PARALLEL_SEED", "42"))
param <- if (ncpus > 1) MulticoreParam(workers = ncpus, RNGseed = seed_R_par) else SerialParam()

# ---------------------------------------------------------------------------
# Env vars
# ---------------------------------------------------------------------------
BOOT_ITER <- as.integer(Sys.getenv("BOOT_ITER", "1"))
N_BOOT    <- as.integer(Sys.getenv("N_BOOT", "5"))

cat("BOOT_ITER:", BOOT_ITER, "\n")
cat("N_BOOT (replicates per task):", N_BOOT, "\n")

if (is.na(BOOT_ITER) || BOOT_ITER < 1) stop("BOOT_ITER must be >= 1")
if (is.na(N_BOOT)   || N_BOOT   < 1) stop("N_BOOT must be >= 1")

iter_tag <- sprintf("%03d", BOOT_ITER)

# ============================================================================
# Load data
# ============================================================================
cat("\n===== Loading data =====\n")

dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

ycfg <- yaml::read_yaml(file.path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design", "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
keep_samples <- dge$samples$dataset %in% mega_cohorts
dge_mega <- dge[, keep_samples]
cat("Mega cohorts:", paste(mega_cohorts, collapse = ", "), "\n")

qc_report <- fread(file.path(INT, "qc/sample_qc_report.csv"))
sex_pass_ids <- qc_report[pass_sex == TRUE, sample_id]
sex_fail <- !colnames(dge_mega) %in% sex_pass_ids
if (any(sex_fail)) dge_mega <- dge_mega[, !sex_fail]

meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]

info <- data.frame(
  group_binary = factor(dge_mega$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = factor(dge_mega$samples$dataset),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE
)
rownames(info) <- colnames(dge_mega)

na_sex <- is.na(info$inferred_sex)
if (any(na_sex)) {
  dge_mega <- dge_mega[, !na_sex]
  info     <- info[!na_sex, , drop = FALSE]
}

cat("Group x sex distribution:\n")
print(table(info$group_binary, info$inferred_sex))

# stratum membership for bootstrap
info$stratum <- interaction(info$dataset, info$group_binary, info$inferred_sex, drop = TRUE)
strata <- split(seq_len(nrow(info)), info$stratum)
cat("N strata for stratified bootstrap:", length(strata), "\n")

# ============================================================================
# Bootstrap function
# ============================================================================
draw_bootstrap <- function(strata_idx_list, seed) {
  set.seed(seed)
  idx <- unlist(lapply(strata_idx_list, function(v) sample(v, length(v), replace = TRUE)),
                use.names = FALSE)
  idx
}

# Helper to refit dream interaction on a bootstrap draw
fit_one_boot <- function(idx_vec, dge_full, info_full, param) {
  # rename duplicated columns: __b1, __b2, ... so dge column names stay unique
  occ <- ave(seq_along(idx_vec), idx_vec, FUN = seq_along)
  new_names <- paste0(colnames(dge_full)[idx_vec], "__b", occ)
  dge_b <- dge_full[, idx_vec]
  colnames(dge_b) <- new_names
  info_b <- info_full[idx_vec, c("group_binary", "dataset", "inferred_sex"), drop = FALSE]
  rownames(info_b) <- new_names

  form_int <- ~ group_binary * inferred_sex + (1|dataset)
  v_int <- suppressWarnings(voomWithDreamWeights(dge_b, form_int, info_b, BPPARAM = param))
  fit_int <- suppressWarnings(dream(v_int, form_int, info_b, BPPARAM = param))
  all_coefs <- colnames(fit_int$coefficients)
  int_coef <- grep("group_binary.*inferred_sex|inferred_sex.*group_binary",
                   all_coefs, value = TRUE)
  if (length(int_coef) == 0) stop("No interaction coefficient found")
  int_coef <- int_coef[1]
  res <- topTable(fit_int, coef = int_coef, number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  res_dt <- as.data.table(res)
  setnames(res_dt, "adj.P.Val", "padj")
  res_dt[, .(gene, logFC, padj)]
}

# ============================================================================
# Run N_BOOT replicates in this array task
# ============================================================================
cat("\n===== Running", N_BOOT, "bootstrap replicates =====\n")

base_seed <- 1000L + BOOT_ITER * 10000L + seed_R_par

logfc_list <- list()
padj_list  <- list()
gene_ref   <- NULL
status     <- character(N_BOOT)

for (b in seq_len(N_BOOT)) {
  seed_b <- base_seed + b
  cat("\n--- Replicate", b, " seed =", seed_b, "---\n")
  idx <- draw_bootstrap(strata, seed_b)
  tryCatch({
    res_b <- fit_one_boot(idx, dge_mega, info, param)
    if (is.null(gene_ref)) gene_ref <- res_b$gene
    res_b <- res_b[match(gene_ref, gene), ]
    logfc_list[[b]] <- res_b$logFC
    padj_list[[b]]  <- res_b$padj
    status[b] <- "ok"
    cat("  sig padj<0.05:", sum(res_b$padj < 0.05, na.rm = TRUE), "\n")
  }, error = function(e) {
    warning("Replicate ", b, " failed: ", conditionMessage(e))
    status[b] <<- paste0("error: ", conditionMessage(e))
    logfc_list[[b]] <<- rep(NA_real_, length(gene_ref))
    padj_list[[b]]  <<- rep(NA_real_, length(gene_ref))
  })
}

# ============================================================================
# Save per-task wide outputs
# ============================================================================
cat("\n===== Saving =====\n")

if (is.null(gene_ref)) stop("All bootstrap replicates failed in this task.")

logfc_mat <- do.call(cbind, lapply(logfc_list, function(v) if (length(v) == length(gene_ref)) v else rep(NA_real_, length(gene_ref))))
padj_mat  <- do.call(cbind, lapply(padj_list,  function(v) if (length(v) == length(gene_ref)) v else rep(NA_real_, length(gene_ref))))

# unique global replicate IDs for this task
rep_ids <- sprintf("b%04d", (BOOT_ITER - 1) * N_BOOT + seq_len(N_BOOT))
colnames(logfc_mat) <- paste0("logFC_",  rep_ids)
colnames(padj_mat)  <- paste0("padj_",   rep_ids)

out_dt <- data.table(gene = gene_ref)
out_dt <- cbind(out_dt, as.data.table(logfc_mat), as.data.table(padj_mat))
fwrite(out_dt, file.path(OUT_DIR, sprintf("boot_int_t%s.csv", iter_tag)))
cat("Saved: boot_int_t", iter_tag, ".csv\n", sep = "")

# status log
status_dt <- data.table(task = BOOT_ITER, replicate = seq_len(N_BOOT), rep_id = rep_ids, status = status)
fwrite(status_dt, file.path(OUT_DIR, sprintf("boot_status_t%s.csv", iter_tag)))

cat(sprintf("\nDone. Elapsed: %.1f min (task=%s, %d reps)\n",
            (proc.time() - t0)[3] / 60, iter_tag, N_BOOT))
