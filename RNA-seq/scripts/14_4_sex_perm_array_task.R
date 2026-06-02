#!/usr/bin/env Rscript
# Sex-interaction permutation null — single-task worker (NO BiocParallel).
# Each SLURM array task does N_PER_TASK permutations in a serial for-loop.
# Total = array_size × N_PER_TASK = 10,000.
#
# Env: TASK_ID (1..100), N_PER_TASK (default 100)
# Out: results/audit_sensitivity/sex_v6_permutation_null/perm_task_{TASK_ID}.rds

set.seed(as.integer(Sys.getenv("TASK_ID", "1")) + 54321L)

suppressPackageStartupMessages({
  library(reformulas); library(lme4); library(data.table); library(edgeR); library(yaml)
})
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) try({
    unlockBinding(fn, ns_lme4)
    assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
    lockBinding(fn, ns_lme4)
  }, silent = TRUE)
}
suppressPackageStartupMessages(library(limma))

TASK_ID    <- as.integer(Sys.getenv("TASK_ID", "1"))
N_PER_TASK <- as.integer(Sys.getenv("N_PER_TASK", "100"))

PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(PROJ, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
OUT_DIR <- file.path(PROJ, "RNA-seq/results/audit_sensitivity/sex_v6_permutation_null")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_FILE <- file.path(OUT_DIR, sprintf("perm_task_%03d.rds", TASK_ID))

if (file.exists(OUT_FILE)) {
  cat("[task", TASK_ID, "] Already done — skipping\n")
  quit(save = "no", status = 0)
}

cat("[task", TASK_ID, "] Loading data...\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
meta <- readRDS(file.path(RDIR, "meta_matched.rds"))

ycfg <- yaml::read_yaml(file.path(PROJ, "config/human_datasets.yaml"))$datasets
mega <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
qc <- fread(file.path(INT, "qc/sample_qc_report.csv"))
pass <- qc[pass_technical == TRUE & dataset %in% mega, sample_id]
keep <- intersect(pass, colnames(dge))
dge <- dge[, keep]

sex_v <- meta$inferred_sex[match(colnames(dge), meta$sample_id)]
info <- data.frame(
  group_binary = factor(dge$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = factor(dge$samples$dataset),
  inferred_sex = factor(sex_v),
  stringsAsFactors = FALSE)
rownames(info) <- colnames(dge)

design <- model.matrix(~ group_binary * inferred_sex, data = info)
int_col <- grep("group_binary.*:inferred_sex|inferred_sex.*:group_binary",
                colnames(design), value = TRUE)
stopifnot(length(int_col) == 1L)

cat("[task", TASK_ID, "] voom + duplicateCorrelation...\n")
v <- voom(dge, design, plot = FALSE)
dupCor <- duplicateCorrelation(v, design, block = info$dataset)
cat("  consensus cor:", round(dupCor$consensus, 4), "\n")

cat("[task", TASK_ID, "] Observed lmFit...\n")
fit_obs <- lmFit(v, design, block = info$dataset, correlation = dupCor$consensus)
fit_obs <- eBayes(fit_obs)
t_obs <- fit_obs$t[, int_col]
abs_t_obs <- abs(t_obs)
ngene <- length(t_obs)

strata <- interaction(info$dataset, info$group_binary, drop = TRUE)
strata_idx <- split(seq_along(strata), strata)

cat("[task", TASK_ID, "] Running", N_PER_TASK, "permutations serially...\n")
n_extreme <- integer(ngene)
t0 <- Sys.time()

for (k in seq_len(N_PER_TASK)) {
  info_p <- info
  new_sex <- info_p$inferred_sex
  for (idx in strata_idx) {
    new_sex[idx] <- sample(info_p$inferred_sex[idx])
  }
  info_p$inferred_sex <- new_sex
  design_p <- model.matrix(~ group_binary * inferred_sex, data = info_p)
  if (!int_col %in% colnames(design_p)) next

  fit_p <- tryCatch(
    lmFit(v, design_p, block = info_p$dataset, correlation = dupCor$consensus),
    error = function(e) NULL)
  if (is.null(fit_p)) next
  fit_p <- eBayes(fit_p)
  t_p <- fit_p$t[, int_col]
  n_extreme <- n_extreme + as.integer(abs(t_p) >= abs_t_obs)

  if (k %% 10 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("[task %d] %d/%d done (%.1f min elapsed, %.1f sec/perm)\n",
                TASK_ID, k, N_PER_TASK, elapsed, elapsed * 60 / k))
  }
}

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("[task %d] COMPLETE — %d perms in %.1f min (%.1f sec/perm)\n",
            TASK_ID, N_PER_TASK, elapsed, elapsed * 60 / N_PER_TASK))

saveRDS(list(task_id = TASK_ID, n_perm = N_PER_TASK,
             n_extreme = n_extreme, gene_ids = names(t_obs),
             t_obs = t_obs, elapsed_min = elapsed),
        OUT_FILE)
cat("Saved:", OUT_FILE, "\n")
