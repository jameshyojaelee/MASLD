#!/usr/bin/env Rscript
# bootstrap_multimethod.R
# ===========================================================================
# Multi-method DE validation harness — Pillar A: cohort-stratified subsampling
# (50%, WITHOUT replacement) run on IDENTICAL splits across dream / DESeq2 /
# (optionally) metafor, so per-method selection-frequency is directly
# comparable on the same resampled data.
#
# Subsampling, not with-replacement bootstrap: ties from replacement distort
# rank-based selection (De Bin 2016 Biometrics). The split is produced by
# de_validation_helpers::stratified_half_indices(), seeded IDENTICALLY to the
# existing dream_bootstrap_iter.R (set.seed(42 + ITER*7919) BEFORE any RNG), so
# this iteration's subsample is byte-identical to that pipeline's iter_{ITER}.
#
# Env:  ITER (integer >= 1; one array task per iteration)
#       VALIDATION_METHODS (comma list; default "dream,deseq2" — metafor gated)
#       SLURM_CPUS_PER_TASK (thread count for BiocParallel)
# Out:  results/integration/multimethod_validation/bootstrap/iter_{ITER:04d}.csv
#       Resumable: if the iter file already exists, message + quit(status=0).
# ===========================================================================

t0 <- proc.time()

# --- 1. Source the shared (TESTED) helpers; do not edit them ----------------
PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HELPERS <- file.path(PROJECT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts",
  "multimethod_validation/de_validation_helpers.R")
source(HELPERS)

ITER <- as.integer(Sys.getenv("ITER", "1"))
if (is.na(ITER) || ITER < 1) stop("ITER must be an integer >= 1")
cat("=== multimethod bootstrap iter", ITER, "===\n")

out_dir  <- file.path(RDIR, "multimethod_validation/bootstrap")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# BOOTSTRAP_OUT_SUFFIX lets a separate single-method run (VALIDATION_METHODS=metafor,
# BOOTSTRAP_OUT_SUFFIX=_metafor) write iter_metafor_{N}.csv on the SAME seeds without
# colliding with / being skipped by the existing dream+deseq2 iter_{N}.csv files.
# Deterministic resample (set.seed(42+ITER*7919)) => identical subsample => merges.
OUT_SUFFIX <- Sys.getenv("BOOTSTRAP_OUT_SUFFIX", "")
out_file <- file.path(out_dir, sprintf("iter%s_%04d.csv", OUT_SUFFIX, ITER))
if (file.exists(out_file)) {
  cat("[multimethod iter", ITER, "] Already done — skipping\n")
  quit(save = "no", status = 0)
}

methods <- ACTIVE_METHODS()
bp      <- bp_param()
cat("ACTIVE_METHODS =", paste(methods, collapse = ", "), "\n")

# Columns of every iter file (long; one row per gene per method).
OUT_COLS <- c("iter", "method", "gene", "logFC", "stat", "padj",
              "n_sub", "n_datasets", "ds_ok")

# --- 2. Seed BEFORE any RNG, then load data ---------------------------------
# CRITICAL: the seed must be set before *any* random draw so the subsample is
# byte-identical to dream_bootstrap_iter.R and across all methods this iter.
set.seed(42 + ITER * 7919)

d <- load_mega_data()
counts <- d$counts
meta   <- d$meta

# Drop NA inferred_sex rows (none expected for mega cohorts, but guard so the
# subsample index space matches what each method runner will accept). Subset
# counts + meta TOGETHER so they stay aligned.
na_sex <- is.na(meta$inferred_sex)
if (any(na_sex)) {
  cat("Dropping", sum(na_sex), "NA inferred_sex sample(s)\n")
  counts <- counts[, !na_sex, drop = FALSE]
  meta   <- meta[!na_sex]
}
stopifnot(all(meta$sample_id == colnames(counts)))

# --- 3. Cohort-stratified 50% subsample (consumes the seed ONCE) ------------
# stratified_half_indices() is the FIRST and ONLY RNG consumer here; no other
# random number is drawn before it, preserving cross-method/cross-pipeline
# reproducibility.
sub_idx <- stratified_half_indices(meta)
counts_s <- counts[, sub_idx, drop = FALSE]
meta_s   <- meta[sub_idx]
cat("Subsample n =", ncol(counts_s), "/", ncol(counts), "\n")

# --- 4. Failure guard: need >= 2 cohorts retaining BOTH groups --------------
ds_ok <- datasets_with_both_groups(meta_s)
write_empty_and_exit <- function() {
  # Mirror dream_bootstrap_iter.R: empty-but-valid file, exit 0 (resumable).
  empty <- rbindlist(lapply(methods, function(m) {
    data.table(iter = ITER, method = m,
               gene = character(0), logFC = numeric(0),
               stat = numeric(0), padj = numeric(0),
               n_sub = ncol(counts_s), n_datasets = length(ds_ok),
               ds_ok = paste(ds_ok, collapse = ";"))
  }))
  # Guarantee a valid header even if `methods` were somehow empty.
  if (nrow(empty) == 0) {
    empty <- data.table(iter = integer(0), method = character(0),
                        gene = character(0), logFC = numeric(0),
                        stat = numeric(0), padj = numeric(0),
                        n_sub = integer(0), n_datasets = integer(0),
                        ds_ok = character(0))
  }
  setcolorder(empty, OUT_COLS)
  fwrite(empty, out_file)
  cat("WARN: <2 datasets retain both groups; wrote empty result, exit 0\n")
  quit(save = "no", status = 0)
}
if (length(ds_ok) < 2) write_empty_and_exit()

# --- 5. Restrict to cohorts with both groups --------------------------------
keep_ok  <- meta_s$dataset %in% ds_ok
counts_s <- counts_s[, keep_ok, drop = FALSE]
meta_s   <- meta_s[keep_ok]
stopifnot(all(meta_s$sample_id == colnames(counts_s)))
cat("Fitting on n =", ncol(counts_s), "across", length(ds_ok), "datasets:",
    paste(ds_ok, collapse = ", "), "\n")

# --- 6. Run each active method on the SAME resampled split ------------------
run_one_method <- function(m) {
  tryCatch({
    if (m == "dream") {
      run_dream(counts_s, meta_s, bp)
    } else if (m == "deseq2") {
      run_deseq2(counts_s, meta_s, bp)
    } else if (m == "metafor") {
      # RECOMPUTE per-study on the RESAMPLED samples — mandatory. Unlike the
      # LOO harness (where a held-out cohort's per-study fit is fixed and can be
      # cached), bootstrap CHANGES the sample membership of every cohort each
      # iteration, so the per-study logFC/SE inputs to metafor must be refit on
      # exactly this subsample. Reusing precomputed per_study/*.csv would
      # silently mix the full-data effect sizes into a resampled meta-analysis.
      per_study_list <- list()
      for (ds in ds_ok) {
        sel <- meta_s$dataset == ds
        ps <- tryCatch(
          run_per_study_voom(ds, counts_s[, sel, drop = FALSE],
                             meta_s[sel], contrast_mode = "binary"),
          error = function(e) {
            cat("  metafor: per-study", ds, "failed:", conditionMessage(e), "\n")
            NULL
          })
        if (!is.null(ps) && nrow(ps) > 0) per_study_list[[ds]] <- ps
      }
      K <- length(per_study_list)
      if (K < 2) {
        cat("  metafor: <2 per-study fits succeeded; skipping metafor\n")
        return(NULL)
      }
      run_metafor(per_study_list, K = K, bp)
    } else {
      cat("  Unknown method '", m, "' — skipping\n", sep = "")
      NULL
    }
  }, error = function(e) {
    cat("  METHOD", m, "ERROR:", conditionMessage(e), "\n")
    NULL
  })
}

results <- list()
for (m in methods) {
  cat("-- running", m, "--\n")
  r <- run_one_method(m)
  if (!is.null(r) && nrow(r) > 0) {
    # Each method returns gene, logFC, stat, padj, method — keep that subset.
    results[[m]] <- r[, .(method, gene, logFC, stat, padj)]
    cat("   ", m, ":", nrow(r), "genes\n")
  } else {
    cat("   ", m, ": no rows\n")
  }
}

# --- 7. Stack long, annotate, write -----------------------------------------
if (length(results) == 0) {
  # All methods failed on a valid split: still emit a valid header so the
  # aggregator does not choke. Treat as empty (0 rows per method handled above).
  write_empty_and_exit()
}

stacked <- rbindlist(results, fill = TRUE)
stacked[, `:=`(iter       = ITER,
               n_sub      = ncol(counts_s),
               n_datasets = length(ds_ok),
               ds_ok      = paste(ds_ok, collapse = ";"))]
setcolorder(stacked, OUT_COLS)
fwrite(stacked, out_file)
cat(sprintf("Elapsed %.1f min -- iter %d done (%d rows, %d methods)\n",
            (proc.time() - t0)["elapsed"] / 60, ITER, nrow(stacked),
            length(results)))
