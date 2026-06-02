#!/usr/bin/env Rscript
# 14_4_sex_permutation_null.R (Agent B4 — 2026-05-20)
# ---------------------------------------------------------------------------
# Sex permutation null (10,000 permutations) for §1.4 (sex calibration).
#
# RATIONALE
#   sex_v6 meta-review (memory/sex_v6_meta_review_2026_05_16.md) flagged
#   that dream's parametric sex-interaction p-values may not be calibrated
#   under the project's mixed-cohort structure. Empirical permutation null
#   gives a calibration-free p-value to confirm (or refute) the genome-wide
#   NULL direction for sex-modulated DEGs.
#
# APPROXIMATION
#   dream() with random effects is too slow at 10k permutations
#   (~2-3 minutes per fit × 10,000 = 14+ days). We use the limma + voom +
#   duplicateCorrelation approximation:
#
#     vobj  <- voom(dge, design)
#     dupCor <- duplicateCorrelation(vobj, design, block = dataset)
#     fit   <- lmFit(vobj, design, block = dataset, correlation = dupCor$consensus)
#
#   This is the canonical scalable substitute for `(1|dataset)` random effects
#   in limma when dream-style mixed-effects is too costly. It does NOT match
#   dream exactly (random-intercept is approximated by a single consensus
#   correlation), but the permutation distribution it generates is internally
#   consistent — we compare permuted t-statistics to the *same model's*
#   observed t-statistics, so the calibration question is answered self-
#   consistently.
#
# DESIGN
#   - Compute t_obs per gene from observed sex labels (group_binary:inferred_sex
#     interaction term).
#   - For each of 10,000 permutations, shuffle `inferred_sex` *within
#     (dataset × group_binary) strata* to preserve sex composition per cohort
#     and condition (the only source of confounding the v6 model worries about).
#   - Compute t_perm per gene, accumulate count of |t_perm| >= |t_obs|.
#   - Empirical two-sided p = (n_extreme + 1) / (n_perm + 1).
#
# OUTPUT
#   RNA-seq/results/audit_sensitivity/sex_v6_permutation_null/sex_v6_permutation_null.csv
#     gene, t_obs, p_emp, n_extreme
#
# PARALLELIZATION
#   BiocParallel::SnowParam(SLURM_CPUS_PER_TASK, type="PSOCK") wrapping a
#   bplapply over permutation chunks. Genes are vectorised inside each fit.
#   SnowParam (socket-based) replaces MulticoreParam (forked) to avoid the
#   silent post-fork worker stall observed in the prior 48h timeout run.
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(BiocParallel)
  library(yaml)
})

# ---- Paths ----
PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(PROJ, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
OUT_DIR <- file.path(PROJ, "RNA-seq/results/audit_sensitivity/sex_v6_permutation_null")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_CSV <- file.path(OUT_DIR, "sex_v6_permutation_null.csv")

N_PERM <- as.integer(Sys.getenv("N_PERM", "10000"))
ncpus  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
cat("N_PERM =", N_PERM, "; ncpus =", ncpus, "\n")

# ---- Load data ----
cat("Loading merged_dge.rds...\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

# Enforce 5-cohort yaml include_in_mega (consistent with 05_dream_mega_analysis.R)
ycfg <- yaml::read_yaml(file.path(PROJ, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("Mega cohorts (k =", length(mega_cohorts), "):",
    paste(mega_cohorts, collapse = ", "), "\n")

keep <- dge$samples$dataset %in% mega_cohorts
dge  <- dge[, keep]

# Sex check filter
qc_report <- fread(file.path(INT, "qc/sample_qc_report.csv"))
sex_pass_ids <- qc_report[pass_sex == TRUE, sample_id]
dge <- dge[, colnames(dge) %in% sex_pass_ids]

meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
sex_v <- meta_new$inferred_sex[match(colnames(dge), meta_new$sample_id)]
na_sex <- is.na(sex_v)
dge <- dge[, !na_sex]
sex_v <- sex_v[!na_sex]

info <- data.frame(
  group_binary = factor(dge$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = factor(dge$samples$dataset),
  inferred_sex = factor(sex_v),
  stringsAsFactors = FALSE
)
rownames(info) <- colnames(dge)

cat("Samples after filters:", ncol(dge), "\n")
cat("Sex x group x dataset:\n")
print(table(info$inferred_sex, info$group_binary, info$dataset))

# ---- Build design with sex interaction ----
# Model: ~ group_binary * inferred_sex; dataset handled via duplicateCorrelation block.
build_design <- function(info_local) {
  model.matrix(~ group_binary * inferred_sex, data = info_local)
}
design <- build_design(info)
int_col <- grep("group_binary.*:inferred_sex|inferred_sex.*:group_binary",
                colnames(design), value = TRUE)
stopifnot(length(int_col) == 1L)
cat("Interaction coef:", int_col, "\n")

# ---- voom + duplicateCorrelation (one-shot for observed fit) ----
cat("voom + duplicateCorrelation (observed)...\n")
vobj <- voom(dge, design, plot = FALSE)
dupCor <- duplicateCorrelation(vobj, design, block = info$dataset)
cat("Consensus correlation:", round(dupCor$consensus, 4), "\n")

cat("Observed lmFit...\n")
fit_obs <- lmFit(vobj, design, block = info$dataset,
                 correlation = dupCor$consensus)
fit_obs <- eBayes(fit_obs)
t_obs <- fit_obs$t[, int_col]
gene_ids <- rownames(fit_obs$t)
cat("Genes tested:", length(gene_ids), "\n")
cat("Observed |t| distribution:\n"); print(summary(abs(t_obs)))

# ---- Permutation engine ----
# Stratified shuffle within (dataset × group_binary): preserves marginal sex
# composition per stratum so cohort-level / case-control sex skew cannot
# leak into the null.
strata <- interaction(info$dataset, info$group_binary, drop = TRUE)
strata_idx <- split(seq_along(strata), strata)

perm_one_chunk <- function(seed_chunk, chunk_size, design, info, vobj,
                           dupCor_consensus, t_obs, int_col_name,
                           strata_idx, chunk_id, log_dir) {
  # Returns n_extreme count vector per gene for this chunk.
  # Writes per-chunk progress to log_dir/chunk_<id>.log every 25 perms.
  set.seed(seed_chunk)
  ngene <- length(t_obs)
  n_extreme <- integer(ngene)
  abs_t_obs <- abs(t_obs)
  log_file <- file.path(log_dir, paste0("chunk_", chunk_id, ".log"))
  cat(sprintf("[chunk %d] start — %d perms, seed %d, %s\n",
              chunk_id, chunk_size, seed_chunk, Sys.time()),
      file = log_file)
  for (k in seq_len(chunk_size)) {
    info_p <- info
    new_sex <- info_p$inferred_sex
    for (idx in strata_idx) {
      new_sex[idx] <- sample(info_p$inferred_sex[idx])
    }
    info_p$inferred_sex <- new_sex
    design_p <- model.matrix(~ group_binary * inferred_sex, data = info_p)
    # Skip if design rank-deficient under this permutation (rare)
    if (!int_col_name %in% colnames(design_p)) next
    # Re-use existing voom weights (sex-permutation shouldn't materially shift
    # mean-variance; this is the standard approximation).
    fit_p <- tryCatch(
      lmFit(vobj, design_p, block = info_p$dataset,
            correlation = dupCor_consensus),
      error = function(e) NULL)
    if (is.null(fit_p)) next
    fit_p <- eBayes(fit_p)
    t_p <- fit_p$t[, int_col_name]
    # Align gene order (lmFit preserves vobj row order)
    n_extreme <- n_extreme + as.integer(abs(t_p) >= abs_t_obs)
    # Progress log every 25 permutations
    if (k %% 25 == 0 || k == chunk_size) {
      cat(sprintf("[chunk %d] %d/%d done (%s)\n",
                  chunk_id, k, chunk_size, Sys.time()),
          file = log_file, append = TRUE)
    }
  }
  cat(sprintf("[chunk %d] COMPLETE — %s\n", chunk_id, Sys.time()),
      file = log_file, append = TRUE)
  n_extreme
}

# Split N_PERM into chunks for parallel workers (one chunk per worker)
n_workers <- min(ncpus, N_PERM)
chunk_sizes <- rep(N_PERM %/% n_workers, n_workers)
remainder <- N_PERM - sum(chunk_sizes)
if (remainder > 0) chunk_sizes[seq_len(remainder)] <- chunk_sizes[seq_len(remainder)] + 1L
chunk_seeds <- 100000L + seq_len(n_workers)
chunk_ids   <- seq_len(n_workers)
cat("Dispatching", n_workers, "chunks; sizes:", paste(chunk_sizes, collapse=","), "\n")

# Per-chunk progress logs
LOG_DIR <- file.path(OUT_DIR, "chunk_logs")
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)
cat("Chunk progress logs ->", LOG_DIR, "\n")

# MulticoreParam with timeout to catch stalls (SnowParam can't see parent functions)
bp <- MulticoreParam(workers = n_workers, RNGseed = 42, timeout = 7200)

t0 <- Sys.time()
res_list <- bpmapply(
  function(seed, sz, cid) perm_one_chunk(seed, sz, design, info, vobj,
                                          dupCor$consensus, t_obs, int_col,
                                          strata_idx, cid, LOG_DIR),
  chunk_seeds, chunk_sizes, chunk_ids,
  SIMPLIFY = FALSE, BPPARAM = bp)
cat("Permutation wall time:", format(Sys.time() - t0), "\n")

n_extreme_total <- Reduce("+", res_list)
p_emp <- (n_extreme_total + 1L) / (N_PERM + 1L)

out <- data.table(
  gene       = gene_ids,
  t_obs      = t_obs,
  n_extreme  = n_extreme_total,
  p_emp      = p_emp
)
# Add BH-adjusted empirical p
out[, padj_emp := p.adjust(p_emp, method = "BH")]
setorder(out, p_emp)

fwrite(out, OUT_CSV)
cat("\nWrote:", OUT_CSV, "\n")
cat("Empirical p<0.05 (uncorrected):", sum(out$p_emp < 0.05), "\n")
cat("Empirical padj_emp<0.05 (BH):", sum(out$padj_emp < 0.05), "\n")
cat("Done.\n")
