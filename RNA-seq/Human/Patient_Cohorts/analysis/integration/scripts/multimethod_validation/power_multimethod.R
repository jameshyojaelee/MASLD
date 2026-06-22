#!/usr/bin/env Rscript
# power_multimethod.R
# ===========================================================================
# Multi-method DE validation harness — Pillar C (POWER): STEP 2 of 2.
#
# Known-truth NB simulation of K = 5 cohorts. For each grid cell (sample size x
# effect size x between-cohort heterogeneity) we simulate n_reps independent
# count matrices, run EVERY active method (dream, deseq2; metafor gated) on the
# SAME simulated matrix within a rep, and measure:
#   power (TPR @ FDR < 0.05), observed FDR, AUROC.
#
# The heterogeneity axis (tau2) directly tests whether metafor's random-effects
# modeling costs power vs the pooled dream / DESeq2 estimands at K = 5: under
# tau2 = 0 the cohorts are homogeneous (RE should pay an unnecessary variance
# penalty); under large tau2 the RE model should be better calibrated.
#
# GRID (48 cells): n_per_group {10,25,50,100} x true_lfc {0.25,0.5,1.0,2.0}
#                  x tau2 {0, 0.0375, 0.15}.  n_reps = 20 per cell.
#   tau2 is on the log2-fold-change variance scale (between-cohort SD of the
#   per-cohort true effect): 0 = homogeneous, 0.0375 ~ SD 0.19 log2FC,
#   0.15 ~ SD 0.39 log2FC.
#
# REPRODUCIBILITY: set.seed(1000 + GRID_TASK*100 + rep) is set BEFORE any RNG
# for that rep, so the simulated matrix is byte-identical across all methods
# within a rep (and deterministic across reruns / across the method set).
#
# Env:  GRID_TASK (1-48; one array task per grid cell)
#       VALIDATION_METHODS (comma list; default "dream,deseq2" — metafor gated)
#       SLURM_CPUS_PER_TASK (thread count for BiocParallel)
# Out:  results/integration/multimethod_validation/power/power_grid_{GRID_TASK}.csv
#         long: one row per (method, rep). Resumable: skip if file exists.
# ===========================================================================

t0 <- proc.time()

# --- 1. Source the shared (TESTED) helpers; do not edit them ----------------
PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HELPERS <- file.path(PROJECT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts",
  "multimethod_validation/de_validation_helpers.R")
source(HELPERS)

suppressPackageStartupMessages(library(pROC))

# --- 2. Read cached NB params (FILE 1 must have run first) ------------------
out_dir   <- file.path(RDIR, "multimethod_validation/power")
param_file <- file.path(out_dir, "nb_params.rds")
if (!file.exists(param_file)) {
  stop("NB parameter cache not found:\n  ", param_file,
       "\nRun power_estimate_nb_params.R first (FILE 1 of the power harness).")
}
nb <- readRDS(param_file)
MU0  <- nb$mu          # baseline expected counts at ref lib size L
PHI0 <- nb$phi         # NB tagwise dispersion (size = 1/phi)
ALL_GENES <- nb$genes
cat("Loaded NB params:", length(ALL_GENES), "genes; ref L =", round(nb$L), "\n")

# --- 3. Grid + constants ----------------------------------------------------
GRID <- expand.grid(
  n_per_group = c(10, 25, 50, 100),
  true_lfc    = c(0.25, 0.5, 1.0, 2.0),
  tau2        = c(0, 0.0375, 0.15),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(nrow(GRID) == 48L)

GRID_TASK <- as.integer(Sys.getenv("GRID_TASK", "1"))
if (is.na(GRID_TASK) || GRID_TASK < 1 || GRID_TASK > nrow(GRID))
  stop("GRID_TASK must be an integer in 1..", nrow(GRID))
row <- GRID[GRID_TASK, ]
cat(sprintf("=== power GRID_TASK %d/48: n_per_group=%d true_lfc=%.2f tau2=%g ===\n",
            GRID_TASK, row$n_per_group, row$true_lfc, row$tau2))

n_reps  <- 20L
K       <- 5L          # cohorts
PI_DE   <- 0.10        # fraction of genes truly DE
PI_SEX  <- 0.05        # fraction of genes with a fixed sex effect
SEX_LFC <- 0.5         # log2FC of the (fixed) sex effect
COHORT_OFFSET_SD <- 0.2 # log2-scale per-cohort baseline batch shift SD
FDR_THRESH <- 0.05

methods <- ACTIVE_METHODS()
bp      <- bp_param()
cat("ACTIVE_METHODS =", paste(methods, collapse = ", "), "\n")

# POWER_OUT_SUFFIX lets a separate single-method run (e.g. VALIDATION_METHODS=metafor,
# POWER_OUT_SUFFIX=_metafor) write power_grid_metafor_{N}.csv on the SAME seeds without
# colliding with a co-running dream+deseq2 array (power_grid_{N}.csv). Deterministic
# seeds => identical sims => own-universe metrics merge cleanly. Default "" = unchanged.
OUT_SUFFIX <- Sys.getenv("POWER_OUT_SUFFIX", "")
out_file <- file.path(out_dir, sprintf("power_grid%s_%d.csv", OUT_SUFFIX, GRID_TASK))
if (file.exists(out_file)) {
  cat("[power GRID_TASK", GRID_TASK, "] Already done — skipping\n")
  quit(save = "no", status = 0)
}

# ===========================================================================
# Simulator: one rep -> list(counts integer genes x samples, meta dt, truth)
#   Counts model for sample j in cohort k with disease status x_j in {0,1}:
#     mean = mu_g * 2^offset_k * 2^(delta_{g,k} * x_j) * 2^(sexeffect_g * sex_j)
#     Y ~ rnbinom(mu = mean, size = 1/phi_g)
#   delta_{g,k} = true_lfc_g + eps_{g,k},  eps ~ N(0, tau2)  (DE genes only)
#   true_lfc_g  = +/- true_lfc (random sign) for DE genes, 0 otherwise.
# ===========================================================================
simulate_rep <- function(n_per_group, true_lfc, tau2) {
  G  <- length(ALL_GENES)
  mu <- MU0; phi <- PHI0

  # --- truth + per-gene signed effect ---
  truth_g  <- rbinom(G, 1, PI_DE)                          # 1 = truly DE
  sign_g   <- sample(c(-1, 1), G, replace = TRUE)
  lfc_g    <- truth_g * sign_g * true_lfc                  # base log2FC per gene

  # --- fixed sex effect on a small fraction of genes ---
  sex_gene <- rbinom(G, 1, PI_SEX)
  sexeff_g <- sex_gene * sample(c(-1, 1), G, replace = TRUE) * SEX_LFC

  # --- per-cohort baseline batch offset (log2 scale) ---
  offset_k <- rnorm(K, 0, COHORT_OFFSET_SD)

  # --- per-cohort per-DE-gene heterogeneous effect delta_{g,k} ---
  # eps drawn for ALL genes x cohorts; only DE genes carry true_lfc, so non-DE
  # genes get delta = 0 (no spurious effect injected).
  eps  <- matrix(if (tau2 > 0) rnorm(G * K, 0, sqrt(tau2)) else 0,
                 nrow = G, ncol = K)
  delta <- (lfc_g + eps) * truth_g                         # genes x cohorts

  # --- assemble samples: K cohorts x 2 arms x n_per_group ---
  n_tot <- K * 2L * n_per_group
  counts <- matrix(0L, nrow = G, ncol = n_tot,
                   dimnames = list(ALL_GENES, NULL))
  sample_id    <- character(n_tot)
  dataset      <- character(n_tot)
  group_binary <- character(n_tot)
  inferred_sex <- character(n_tot)

  col <- 0L
  for (k in seq_len(K)) {
    for (grp in c(0L, 1L)) {                 # 0 = Control, 1 = Disease
      for (r in seq_len(n_per_group)) {
        col <- col + 1L
        sex01 <- rbinom(1, 1, 0.5)           # 50/50 sex assignment
        # log2-scale linear predictor -> expected counts
        log2mean <- log2(pmax(mu, 1e-8)) + offset_k[k] +
                    delta[, k] * grp + sexeff_g * sex01
        mean_g   <- 2^log2mean
        counts[, col] <- rnbinom(G, mu = mean_g, size = 1 / phi)
        sample_id[col]    <- sprintf("sim_C%d_%s_%d",
                                     k, if (grp == 1L) "D" else "C", r)
        dataset[col]      <- paste0("simC", k)
        group_binary[col] <- if (grp == 1L) "Disease" else "Control"
        inferred_sex[col] <- if (sex01 == 1L) "male" else "female"
      }
    }
  }
  storage.mode(counts) <- "integer"
  colnames(counts) <- sample_id

  meta <- data.table(
    sample_id    = sample_id,
    dataset      = dataset,
    group_binary = group_binary,
    inferred_sex = inferred_sex,
    condition    = group_binary           # sim has only Control/Disease
  )
  list(counts = counts, meta = meta, truth = setNames(truth_g, ALL_GENES))
}

# ===========================================================================
# Run one method on a simulated (counts, meta).
#   metafor: refit per-study on each cohort's columns (contrast_mode="binary",
#   matching dream/DESeq2 estimand), then run_metafor over K cohorts.
# ===========================================================================
run_one_method <- function(m, counts, meta) {
  tryCatch({
    if (m == "dream") {
      run_dream(counts, meta, bp)
    } else if (m == "deseq2") {
      run_deseq2(counts, meta, bp)
    } else if (m == "metafor") {
      per_study_list <- list()
      for (ds in sort(unique(meta$dataset))) {
        sel <- meta$dataset == ds
        ps <- tryCatch(
          run_per_study_voom(ds, counts[, sel, drop = FALSE], meta[sel],
                             contrast_mode = "binary"),
          error = function(e) {
            cat("  metafor: per-study", ds, "failed:", conditionMessage(e), "\n")
            NULL
          })
        if (!is.null(ps) && nrow(ps) > 0) per_study_list[[ds]] <- ps
      }
      Kfit <- length(per_study_list)
      if (Kfit < 2) {
        cat("  metafor: <2 per-study fits succeeded; skipping\n")
        return(NULL)
      }
      run_metafor(per_study_list, K = Kfit, bp)
    } else if (m == "limma_voom") {
      run_limma_voom(counts, meta, bp)
    } else if (m == "limma_voom_qw") {
      run_limma_voom_qw(counts, meta, bp)
    } else if (m == "limma_trend") {
      run_limma_trend(counts, meta, bp)
    } else if (m == "edger_qlf") {
      run_edger_qlf(counts, meta, bp)
    } else if (m == "edger_qlf_robust") {
      run_edger_qlf_robust(counts, meta, bp)
    } else if (m == "edger_lrt") {
      run_edger_lrt(counts, meta, bp)
    } else {
      cat("  Unknown method '", m, "' — skipping\n", sep = "")
      NULL
    }
  }, error = function(e) {
    cat("  METHOD", m, "ERROR:", conditionMessage(e), "\n")
    NULL
  })
}

# ===========================================================================
# Metric helpers (on a given tested-gene universe).
#   res : method dt with gene, stat, padj (NA-padj rows = tested-but-not-called
#         OR not-in-universe; we restrict to non-NA padj as the tested set).
#   truth: named 0/1 vector over ALL_GENES.
#   universe: character vector of genes to score on (NULL = method's own tested).
# ===========================================================================
score_on <- function(res, truth, universe = NULL) {
  tested <- res[!is.na(padj), gene]
  if (!is.null(universe)) tested <- intersect(tested, universe)
  if (length(tested) < 3)
    return(list(power = NA_real_, fdr = NA_real_, auroc = NA_real_,
                n_tested = length(tested)))
  sub  <- res[gene %in% tested]
  setkey(sub, gene)
  tr   <- truth[sub$gene]
  call <- as.integer(sub$padj < FDR_THRESH)

  n_called <- sum(call == 1L)
  # power = TPR among truly-DE & tested
  is_de  <- tr == 1L
  power  <- if (sum(is_de) > 0) mean(call[is_de] == 1L) else NA_real_
  # observed FDR = false discoveries / total discoveries.
  # If the method called NOTHING, FDR is undefined -> NA (not 0): a no-call cell
  # at low power must not masquerade as "perfectly FDR-controlled" (R3 P2-2).
  fdr    <- if (n_called == 0L) NA_real_ else sum(tr == 0L & call == 1L) / n_called
  # AUROC: |stat| ranks truth; need both classes present
  auroc <- NA_real_
  if (length(unique(tr)) == 2L) {
    pred <- abs(sub$stat)
    ok   <- is.finite(pred)
    if (sum(ok) >= 3 && length(unique(tr[ok])) == 2L) {
      auroc <- tryCatch(
        as.numeric(pROC::auc(response = tr[ok], predictor = pred[ok],
                             direction = "<", quiet = TRUE)),
        error = function(e) NA_real_)
    }
  }
  list(power = power, fdr = fdr, auroc = auroc, n_tested = length(tested))
}

# ===========================================================================
# Main rep loop
# ===========================================================================
all_rows <- list()
for (rep in seq_len(n_reps)) {
  set.seed(1000 + GRID_TASK * 100 + rep)   # BEFORE any RNG for this rep
  sim    <- simulate_rep(row$n_per_group, row$true_lfc, row$tau2)
  truth  <- sim$truth
  n_truth_de <- sum(truth == 1L)

  # 1. Run every active method on the SAME simulated matrix.
  res_list <- list()
  runtimes <- setNames(rep(NA_real_, length(methods)), methods)
  for (m in methods) {
    tm <- proc.time()
    r  <- run_one_method(m, sim$counts, sim$meta)
    runtimes[m] <- (proc.time() - tm)["elapsed"] / 60
    if (!is.null(r) && nrow(r) > 0) res_list[[m]] <- r
  }
  if (length(res_list) == 0) {
    cat("  rep", rep, ": all methods failed — skipping\n")
    next
  }

  # 2. Common universe = intersection of the tested (non-NA padj) sets of all
  #    ACTIVE methods that produced output this rep.
  tested_sets <- lapply(res_list, function(r) r[!is.na(padj), gene])
  common_universe <- Reduce(intersect, tested_sets)

  # 3. Score each method on its OWN universe and the COMMON universe.
  for (m in names(res_list)) {
    r       <- res_list[[m]]
    own     <- score_on(r, truth, universe = NULL)
    comm    <- score_on(r, truth, universe = common_universe)
    all_rows[[length(all_rows) + 1L]] <- data.table(
      grid_task       = GRID_TASK,
      rep             = rep,
      n_per_group     = row$n_per_group,
      true_lfc        = row$true_lfc,
      tau2            = row$tau2,
      method          = m,
      n_truth_de      = n_truth_de,
      n_tested        = own$n_tested,
      n_tested_common = length(common_universe),
      power_own       = own$power,
      fdr_own         = own$fdr,
      auroc_own       = own$auroc,
      power_common    = comm$power,
      fdr_common      = comm$fdr,
      auroc_common    = comm$auroc,
      runtime_min     = round(runtimes[m], 3)
    )
  }
  cat(sprintf("  rep %2d/%d done | DE=%d | common=%d | %s\n",
              rep, n_reps, n_truth_de, length(common_universe),
              paste(names(res_list), collapse = ",")))
}

# --- Write long output ------------------------------------------------------
OUT_COLS <- c("grid_task", "rep", "n_per_group", "true_lfc", "tau2", "method",
              "n_truth_de", "n_tested", "n_tested_common",
              "power_own", "fdr_own", "auroc_own",
              "power_common", "fdr_common", "auroc_common", "runtime_min")

if (length(all_rows) == 0) {
  out <- data.table(grid_task = integer(0), rep = integer(0),
                    n_per_group = integer(0), true_lfc = numeric(0),
                    tau2 = numeric(0), method = character(0),
                    n_truth_de = integer(0), n_tested = integer(0),
                    n_tested_common = integer(0),
                    power_own = numeric(0), fdr_own = numeric(0),
                    auroc_own = numeric(0), power_common = numeric(0),
                    fdr_common = numeric(0), auroc_common = numeric(0),
                    runtime_min = numeric(0))
} else {
  out <- rbindlist(all_rows, fill = TRUE)
}
setcolorder(out, OUT_COLS)
fwrite(out, out_file)
cat(sprintf("Elapsed %.1f min -- GRID_TASK %d done (%d rows)\n",
            (proc.time() - t0)["elapsed"] / 60, GRID_TASK, nrow(out)))
