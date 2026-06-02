#!/usr/bin/env Rscript
# 14.5a_T3_cohort_perm.R
# ---------------------------------------------------------------------------
# A7 T3: Cohort-stratified permutation null for sex-interaction model
#
# Permute `inferred_sex` labels WITHIN each cohort x disease group cell
# (preserves marginals). For each permutation, refit the dream interaction
# model and tabulate n_sex_dimorphic (interaction padj < 0.05).
#
# Chunked array: SLURM_ARRAY_TASK_ID selects perm range (default B=500 over 10 chunks).
#
# Outputs:
#   .../audit_sensitivity/sex_cohort_permutation/chunk_<id>_null.csv (per array task)
#   .../audit_sensitivity/sex_cohort_permutation/null_distribution.csv (aggregated by post-merge)
#   .../audit_sensitivity/sex_cohort_permutation/empirical_pvalue.csv (post-merge)
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

# Inject reformulas funcs into lme4 namespace before loading variancePartition
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
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
OUT  <- file.path(BASE, "RNA-seq/results/audit_sensitivity/sex_cohort_permutation")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# --- Parallel setup ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("CPUs:", ncpus, "\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# --- Array task ID ---
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
n_perms_per_chunk <- as.integer(Sys.getenv("N_PERMS_PER_CHUNK", "50"))
B_total <- as.integer(Sys.getenv("B_TOTAL", "500"))
cat("Chunk task_id:", task_id, " perms/chunk:", n_perms_per_chunk, "\n")

# --- Load data (same setup as Script 26) ---
cat("\nLoading data...\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("Mega cohorts (k =", length(mega_cohorts), "):", paste(mega_cohorts, collapse = ", "), "\n")

keep_samples <- dge$samples$dataset %in% mega_cohorts
dge_mega <- dge[, keep_samples]

qc_report <- fread(file.path(INT, "qc/sample_qc_report.csv"))
sex_pass_ids <- qc_report[pass_sex == TRUE, sample_id]
sex_fail <- !colnames(dge_mega) %in% sex_pass_ids
if (any(sex_fail)) dge_mega <- dge_mega[, !sex_fail]
cat("Samples after sex check filter:", ncol(dge_mega), "\n")

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
  info <- info[!na_sex, , drop = FALSE]
}
cat("Final samples:", nrow(info), "\n")
cat("group x sex table:\n"); print(table(info$group_binary, info$inferred_sex))

# --- Pre-filter genes (filterByExpr) so each perm fit uses same gene set ---
keep_g <- filterByExpr(dge_mega, group = info$group_binary)
dge_filt <- dge_mega[keep_g, , keep.lib.sizes = FALSE]
dge_filt <- calcNormFactors(dge_filt, method = "TMM")
cat("Genes after filterByExpr:", nrow(dge_filt), "\n")

# --- Permutation function: shuffle sex WITHIN each cohort x disease cell ---
permute_sex_within_cells <- function(info_df, seed) {
  set.seed(seed)
  out <- info_df
  cells <- split(seq_len(nrow(out)), list(out$dataset, out$group_binary), drop = TRUE)
  for (cell_idx in cells) {
    out$inferred_sex[cell_idx] <- sample(out$inferred_sex[cell_idx])
  }
  out
}

# --- Run a single permutation, return n_sex_dimorphic ---
run_perm <- function(perm_id, info_obs, dge_obj, BPPARAM) {
  info_p <- permute_sex_within_cells(info_obs, seed = 1000L + perm_id)
  form <- ~ group_binary * inferred_sex + (1 | dataset)
  v <- tryCatch(
    suppressWarnings(voomWithDreamWeights(dge_obj, form, info_p, BPPARAM = BPPARAM)),
    error = function(e) { cat("voom err perm", perm_id, ":", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(v)) return(NA_integer_)
  fit <- tryCatch(
    suppressWarnings(dream(v, form, info_p, BPPARAM = BPPARAM)),
    error = function(e) { cat("dream err perm", perm_id, ":", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(fit)) return(NA_integer_)
  all_coefs <- colnames(fit$coefficients)
  int_coef <- grep("group_binary.*inferred_sex|inferred_sex.*group_binary",
                   all_coefs, value = TRUE)
  if (length(int_coef) == 0) return(NA_integer_)
  int_coef <- int_coef[1]
  res <- topTable(fit, coef = int_coef, number = Inf, sort.by = "none")
  n_sig_05 <- sum(res$adj.P.Val < 0.05, na.rm = TRUE)
  n_sig_10 <- sum(res$adj.P.Val < 0.10, na.rm = TRUE)
  c(n_sig_pAdj05 = n_sig_05, n_sig_pAdj10 = n_sig_10)
}

# --- Compute perm range for this chunk ---
start_perm <- (task_id - 1L) * n_perms_per_chunk + 1L
end_perm   <- min(task_id * n_perms_per_chunk, B_total)
cat("\nRunning perms", start_perm, "to", end_perm, "\n")

out_chunk <- data.table(perm_id = integer(),
                        n_sig_pAdj05 = integer(),
                        n_sig_pAdj10 = integer())

for (pid in start_perm:end_perm) {
  t0 <- Sys.time()
  res <- run_perm(pid, info, dge_filt, param)
  dt_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  Perm %3d: pAdj<0.05 n=%d  pAdj<0.10 n=%d  (%.1fs)\n",
              pid, res["n_sig_pAdj05"], res["n_sig_pAdj10"], dt_s))
  out_chunk <- rbind(out_chunk, data.table(
    perm_id = pid,
    n_sig_pAdj05 = as.integer(res["n_sig_pAdj05"]),
    n_sig_pAdj10 = as.integer(res["n_sig_pAdj10"])
  ))
  # Save after each perm for resumability
  fwrite(out_chunk, file.path(OUT, sprintf("chunk_%02d_null.csv", task_id)))
}

cat("\nChunk", task_id, "complete. Saved", nrow(out_chunk), "perms\n")
cat("Done:", as.character(Sys.time()), "\n")
