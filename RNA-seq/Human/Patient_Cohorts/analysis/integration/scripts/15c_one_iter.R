#!/usr/bin/env Rscript
# 15c_one_iter.R
# ---------------------------------------------------------------------------
# ONE iteration of 15c power saturation, parameterized by SLURM array task ID
# (or by env vars TARGET_N + ITER). Writes one row CSV to be aggregated later.
#
# Array layout: SLURM_ARRAY_TASK_ID 1..60
#   target_n = c(200, 400, 600, 800, 1000, 1200)[((id - 1) %/% 10) + 1]
#   iter     = ((id - 1) %% 10) + 1
#
# Output: <OUT_DIR>/iter/iter_N{n}_i{iter}.csv  (single row per task)
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

# --- Paths ---
PROJECT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
BASE    <- file.path(PROJECT, "RNA-seq/Human/Patient_Cohorts")
INT     <- file.path(BASE, "analysis/integration")
RDIR    <- file.path(INT, "results/integration")
OUT_DIR <- file.path(PROJECT, "RNA-seq/results/audit_sensitivity/power_saturation")
ITER_DIR <- file.path(OUT_DIR, "iter")
dir.create(ITER_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Resolve task ID -> (target_n, iter) ---
TARGET_NS <- c(200L, 400L, 600L, 800L, 1000L, 1200L)
N_ITER    <- 10L
PADJ_THRESH <- 0.1

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
if (task_id == 0L) {
  # Allow manual env override for debugging
  TGT  <- as.integer(Sys.getenv("TARGET_N", ""))
  ITER <- as.integer(Sys.getenv("ITER", ""))
  if (is.na(TGT) || is.na(ITER)) {
    stop("Set SLURM_ARRAY_TASK_ID 1..60, or both TARGET_N and ITER env vars.")
  }
  target_n <- TGT
  iter     <- ITER
} else {
  if (task_id < 1L || task_id > length(TARGET_NS) * N_ITER) {
    stop("SLURM_ARRAY_TASK_ID out of range: ", task_id)
  }
  target_n <- TARGET_NS[((task_id - 1L) %/% N_ITER) + 1L]
  iter     <- ((task_id - 1L) %% N_ITER) + 1L
}

seed <- 42000L + target_n + iter
out_csv <- file.path(ITER_DIR, sprintf("iter_N%d_i%d.csv", target_n, iter))
cat(sprintf("Task %d -> target_n=%d iter=%d seed=%d\n",
            task_id, target_n, iter, seed))
cat("Output:", out_csv, "\n")

# Skip if already exists (resumable)
if (file.exists(out_csv) && nzchar(Sys.getenv("FORCE", ""))) {
  cat("File exists and FORCE not set — skipping.\n")
  quit(save = "no", status = 0)
}

set.seed(seed)

# --- Parallel ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("CPUs:", ncpus, "\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# --- Load reference (full dream + positive controls) ---
full_dream <- fread(file.path(RDIR, "dream_results.csv"))
full_deg_genes <- full_dream[padj < PADJ_THRESH, gene]
full_lfc <- setNames(full_dream$logFC, full_dream$gene)

pc <- fread(file.path(PROJECT, "results/library/positive_control.csv"))
pc_genes <- intersect(pc$`Gene symbol`, full_dream$gene)

# --- Load DGE + apply yaml-aware mega filter ---
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
ycfg <- yaml::read_yaml(file.path(PROJECT, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
keep_samples <- dge$samples$dataset %in% mega_cohorts
dge_mega <- dge[, keep_samples]

meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]

info <- data.frame(
  sample_id    = colnames(dge_mega),
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

N_FULL   <- ncol(dge_mega)
datasets <- levels(info$dataset)
ds_indices <- lapply(datasets, function(d) which(info$dataset == d))
names(ds_indices) <- datasets
ds_sizes <- sapply(ds_indices, length)

if (target_n > N_FULL) {
  cat("Skipping: target N exceeds available samples\n")
  quit(save = "no", status = 0)
}

# --- Stratified sampling (same logic as monolithic 15c) ---
frac <- target_n / N_FULL
draw_per_ds <- setNames(pmax(round(ds_sizes * frac), 2L), datasets)

total_allocated <- sum(draw_per_ds)
if (total_allocated != target_n) {
  diff <- total_allocated - target_n
  ds_order <- names(sort(ds_sizes, decreasing = TRUE))
  i <- 1L
  while (diff != 0) {
    d <- ds_order[((i - 1L) %% length(ds_order)) + 1L]
    if (diff > 0 && draw_per_ds[d] > 2L) {
      draw_per_ds[d] <- draw_per_ds[d] - 1L
      diff <- diff - 1L
    } else if (diff < 0 && draw_per_ds[d] < ds_sizes[d]) {
      draw_per_ds[d] <- draw_per_ds[d] + 1L
      diff <- diff + 1L
    }
    i <- i + 1L
    if (i > length(ds_order) * abs(diff) + 100L) break
  }
}

sampled_idx <- integer(0)
for (d in datasets) {
  d_idx <- ds_indices[[d]]
  d_info <- info[d_idx, ]
  n_draw <- min(draw_per_ds[d], length(d_idx))

  ctrl_idx <- d_idx[d_info$group_binary == "Control"]
  dis_idx  <- d_idx[d_info$group_binary == "Disease"]

  if (length(ctrl_idx) == 0 || length(dis_idx) == 0) {
    sampled_idx <- c(sampled_idx, sample(d_idx, n_draw))
    next
  }
  ctrl_frac <- length(ctrl_idx) / length(d_idx)
  n_ctrl <- max(round(n_draw * ctrl_frac), 1L)
  n_dis  <- n_draw - n_ctrl
  n_ctrl <- min(n_ctrl, length(ctrl_idx))
  n_dis  <- min(n_dis, length(dis_idx))
  sampled_idx <- c(sampled_idx,
                   sample(ctrl_idx, n_ctrl),
                   sample(dis_idx, n_dis))
}

dge_sub <- dge_mega[, sampled_idx]
info_sub <- info[sampled_idx, , drop = FALSE]
rownames(info_sub) <- colnames(dge_sub)
info_sub$dataset <- droplevels(info_sub$dataset)
actual_n <- ncol(dge_sub)

if (length(unique(info_sub$group_binary)) < 2) {
  cat("SKIP (missing group level)\n")
  quit(save = "no", status = 0)
}
if (nlevels(info_sub$dataset) < 2) {
  cat("SKIP (single dataset)\n")
  quit(save = "no", status = 0)
}

# --- Run dream ---
form <- ~ group_binary + inferred_sex + (1 | dataset)
v_sub <- suppressWarnings(voomWithDreamWeights(dge_sub, form, info_sub, BPPARAM = param))
fit_sub <- suppressWarnings(dream(v_sub, form, info_sub, BPPARAM = param))
res_sub <- topTable(fit_sub, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res_sub$gene <- rownames(res_sub)
res_dt <- as.data.table(res_sub)

# --- Metrics ---
sub_deg_genes <- res_dt[adj.P.Val < PADJ_THRESH, gene]
n_degs <- length(sub_deg_genes)
shared_genes <- intersect(res_dt$gene, full_dream$gene)
sub_lfc <- setNames(res_dt$logFC, res_dt$gene)
rho <- cor(full_lfc[shared_genes], sub_lfc[shared_genes],
           method = "spearman", use = "complete.obs")

shared_degs <- intersect(sub_deg_genes, full_deg_genes)
dir_concordance <- if (length(shared_degs) > 0) {
  mean(sign(full_lfc[shared_degs]) == sign(sub_lfc[shared_degs]))
} else NA_real_

union_degs <- length(union(sub_deg_genes, full_deg_genes))
jaccard <- if (union_degs > 0) length(shared_degs) / union_degs else 0
pc_recovery <- sum(pc_genes %in% sub_deg_genes) / length(pc_genes)

row <- data.table(
  target_n      = target_n,
  actual_n      = actual_n,
  iteration     = iter,
  seed          = seed,
  n_degs        = n_degs,
  n_shared_degs = length(shared_degs),
  spearman_rho  = rho,
  dir_concordance = dir_concordance,
  jaccard       = jaccard,
  pc_recovery   = pc_recovery,
  n_genes_tested = nrow(res_dt)
)
fwrite(row, out_csv)

elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("Done: N=%d iter=%d DEGs=%d rho=%.3f jaccard=%.3f PC=%.1f%% (%.0f min)\n",
            target_n, iter, n_degs, rho, jaccard, pc_recovery * 100, elapsed/60))
