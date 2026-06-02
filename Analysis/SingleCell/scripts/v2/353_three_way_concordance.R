#!/usr/bin/env Rscript
# ============================================================================
# 353_three_way_concordance.R
#
# Addresses critique #9: v1↔v2 LMM Jaccard top-50 = 0.27 across two
# scVI-based atlas builds. To distinguish "this LR ranking is method-
# dependent" from "this LR ranking is real", refit the LMM on F-stage
# predictions from THREE independent integration latents (scVI, Harmony,
# Scanorama from 352a + 352b), then compute three-way Jaccard.
#
# A robust LR pair = appears in top-50 ranking for at least 2 of 3 methods.
#
# Inputs (from 352b):
#   - results_gpu_v2_phase05/ccc/stage_trajectory_v3/donor_fstage_{scvi,harmony,scanorama}.tsv
#   - results_gpu_v2_phase05/ccc/stage_trajectory_v2/all_donor_lr_scores_v2.tsv.gz
#   - results_gpu_v2_phase05/mcp/inputs/donor_metadata_v2.tsv (covariates)
#
# Output: stage_trajectory_v3/three_way_concordance_v3.tsv
#   Long format: lr_key, ct_pair, lr_pair, scvi_rank, harmony_rank,
#                scanorama_rank, in_top50_scvi, in_top50_harmony,
#                in_top50_scanorama, n_methods_in_top50
# Plus: stage_trajectory_v3/three_way_jaccard_v3.tsv (summary metrics)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(parallel)
})

N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
options(mc.cores = N_CORES)
TOP_N <- 50L
SMALL_COHORTS <- c("GSE174748", "GSE189600")
LATENTS <- c("scvi", "harmony", "scanorama")

cat(sprintf("[config] TOP_N=%d  cores=%d\n", TOP_N, N_CORES))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
V2_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2")
V3_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3")
META_V2 <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/mcp/inputs/donor_metadata_v2.tsv")
MERGED_TSV <- file.path(V2_DIR, "all_donor_lr_scores_v2.tsv.gz")

# Load LR scores once
lr_long <- fread(MERGED_TSV)
lr_long[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr_long[, ct_pair := paste(source, target, sep = "->")]
lr_long[, score := -log10(pmax(magnitude_rank, 1e-4))]

# Load covariates (age, sex)
meta <- fread(META_V2)
meta[, dataset_pooled := ifelse(dataset %in% SMALL_COHORTS, "Other_small",
                                as.character(dataset))]
meta_covar <- meta[, .(sample, dataset, dataset_pooled, age, sex_numeric)]

fit_axis_lmm <- function(lr_dt, axis_col) {
  splits <- split(lr_dt, by = c("ct_pair", "lr_pair"), drop = TRUE)
  cat(sprintf("[lmm] %d groups\n", length(splits)))
  wrapper <- function(d) {
    if (nrow(d) < 30) return(NULL)
    use_random <- length(unique(d$dataset_pooled)) >= 3L
    rhs_terms <- c(axis_col, "log10(pmax(n_source_cells,1))",
                   "log10(pmax(n_target_cells,1))")
    # Master review M-P0-3: dataset_pooled must NOT be both fixed and random.
    if (!use_random && length(unique(d$dataset_pooled)) > 1)
      rhs_terms <- c(rhs_terms, "dataset_pooled")
    if ("age" %in% names(d) && sum(!is.na(d$age)) >= 10 &&
        length(unique(stats::na.omit(d$age))) > 1) rhs_terms <- c(rhs_terms, "age")
    if ("sex_numeric" %in% names(d) && sum(!is.na(d$sex_numeric)) >= 10 &&
        length(unique(stats::na.omit(d$sex_numeric))) > 1)
      rhs_terms <- c(rhs_terms, "sex_numeric")
    rhs <- paste(rhs_terms, collapse = " + ")
    if (use_random) rhs <- paste(rhs, "+ (1|dataset_pooled)")
    formula_str <- paste("score ~", rhs)
    fit <- tryCatch(
      if (use_random) lmerTest::lmer(stats::as.formula(formula_str),
                                     data = d, REML = TRUE)
      else stats::lm(stats::as.formula(formula_str), data = d),
      error = function(e) NULL, warning = function(w) NULL)
    if (is.null(fit)) return(NULL)
    co <- try(summary(fit)$coefficients, silent = TRUE)
    if (inherits(co, "try-error")) return(NULL)
    rn <- rownames(co)
    if (!(axis_col %in% rn)) return(NULL)
    data.table(
      lr_key = paste(d$ct_pair[1], d$lr_pair[1], sep = "||"),
      ct_pair = d$ct_pair[1], lr_pair = d$lr_pair[1],
      ligand_complex = d$ligand_complex[1], receptor_complex = d$receptor_complex[1],
      Estimate = co[axis_col, 1], pval = co[axis_col, ncol(co)],
      n_donors = nrow(d)
    )
  }
  res <- parallel::mclapply(splits, wrapper, mc.cores = N_CORES)
  rbindlist(res[!sapply(res, is.null)])
}

results_per_latent <- list()
for (latent in LATENTS) {
  fpath <- file.path(V3_DIR, sprintf("donor_fstage_%s.tsv", latent))
  if (!file.exists(fpath)) {
    cat(sprintf("[warn] %s missing -> skipping\n", fpath))
    next
  }
  cat(sprintf("\n[%s] loading %s\n", latent, fpath))
  pred <- fread(fpath)
  if (!"F_stage_predicted_argmax" %in% names(pred)) {
    cat(sprintf("[warn] %s missing F_stage_predicted_argmax\n", fpath))
    next
  }
  pred[, F_stage_pred := as.numeric(F_stage_predicted_argmax)]
  pred_use <- pred[, .(sample, F_stage_pred)]
  d <- merge(lr_long, pred_use, by = "sample", all.x = FALSE)
  d <- merge(d, meta_covar, by = "sample", all.x = TRUE)
  d <- d[!is.na(F_stage_pred)]
  cat(sprintf("[%s] %d donor-LR rows after join\n", latent, nrow(d)))
  res <- fit_axis_lmm(d, "F_stage_pred")
  if (nrow(res) == 0) {
    cat(sprintf("[%s] no LMM fits\n", latent))
    next
  }
  res[, latent := latent]
  # Rank by abs(Estimate) * -log10(pval)
  res[, rank_score := abs(Estimate) * (-log10(pmax(pval, 1e-30)))]
  setorder(res, -rank_score)
  res[, rank := seq_len(.N)]
  results_per_latent[[latent]] <- res
  cat(sprintf("[%s] %d LMM fits; top hit: %s | est=%.3f p=%.2e\n",
              latent, nrow(res), res$lr_pair[1], res$Estimate[1], res$pval[1]))
}

if (length(results_per_latent) < 2) {
  stop("Fewer than 2 latents produced LMM results; cannot compute concordance")
}

# All unique LR keys across latents
all_keys <- unique(unlist(lapply(results_per_latent, function(r) r$lr_key)))
out <- data.table(lr_key = all_keys)
for (latent in names(results_per_latent)) {
  r <- results_per_latent[[latent]]
  m <- r[, .(lr_key, rank)]
  setnames(m, "rank", paste0(latent, "_rank"))
  out <- merge(out, m, by = "lr_key", all.x = TRUE)
  topk <- r[rank <= TOP_N, lr_key]
  out[, paste0("in_top", TOP_N, "_", latent) := lr_key %in% topk]
}
in_top_cols <- grep(sprintf("^in_top%d_", TOP_N), names(out), value = TRUE)
out[, n_methods_in_topN := rowSums(.SD), .SDcols = in_top_cols]
out[, ct_pair := sub("\\|\\|.*", "", lr_key)]
out[, lr_pair := sub("^[^|]+\\|\\|", "", lr_key)]
setcolorder(out, c("ct_pair", "lr_pair", "lr_key"))

OUT_LONG <- file.path(V3_DIR, "three_way_concordance_v3.tsv")
fwrite(out, OUT_LONG, sep = "\t")
cat(sprintf("\n[output] %d LR-key rows -> %s\n", nrow(out), OUT_LONG))

# Pairwise + three-way Jaccard
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
top_sets <- list()
for (latent in names(results_per_latent)) {
  top_sets[[latent]] <- results_per_latent[[latent]][rank <= TOP_N, lr_key]
}
pair_rows <- list()
nm <- names(top_sets)
for (i in seq_along(nm)) {
  for (j in seq_along(nm)) {
    if (j <= i) next
    pair_rows[[length(pair_rows) + 1]] <- data.table(
      method_A = nm[i], method_B = nm[j],
      n_A = length(top_sets[[nm[i]]]), n_B = length(top_sets[[nm[j]]]),
      n_AB = length(intersect(top_sets[[nm[i]]], top_sets[[nm[j]]])),
      jaccard = jaccard(top_sets[[nm[i]]], top_sets[[nm[j]]])
    )
  }
}
pair_dt <- rbindlist(pair_rows)
three_way <- if (length(top_sets) == 3) {
  Reduce(intersect, top_sets)
} else character(0)
union_all <- Reduce(union, top_sets)
summary_dt <- rbind(
  pair_dt,
  data.table(method_A = "ALL3", method_B = "ALL3",
             n_A = NA, n_B = NA,
             n_AB = length(three_way),
             jaccard = ifelse(length(union_all) > 0,
                              length(three_way) / length(union_all),
                              NA_real_)),
  fill = TRUE)

# Master review M-P0-10: permutation null for Jaccard. The observed three-way
# Jaccard of 0.052 was reported as "low concordance confirms method-dependence",
# but this framing requires comparison to the random-pairing null. Shuffle each
# method's top-N labels within the per-method universe; recompute pairwise +
# three-way Jaccard; report empirical p-value and enrichment fold.
N_PERMS <- as.integer(Sys.getenv("JACCARD_N_PERMS", "1000"))
cat(sprintf("\n[perm] computing %d-permutation null for Jaccard\n", N_PERMS))
set.seed(20260516L)

method_universes <- lapply(results_per_latent, function(r) r$lr_key)
method_topN <- lapply(method_universes, function(u) min(TOP_N, length(u)))

perm_pair <- matrix(0.0, nrow = N_PERMS, ncol = nrow(pair_dt))
perm_all3 <- numeric(N_PERMS)
for (p in seq_len(N_PERMS)) {
  perm_top <- mapply(function(u, k) sample(u, k), method_universes, method_topN,
                     SIMPLIFY = FALSE)
  for (k in seq_len(nrow(pair_dt))) {
    a <- perm_top[[pair_dt$method_A[k]]]
    b <- perm_top[[pair_dt$method_B[k]]]
    perm_pair[p, k] <- jaccard(a, b)
  }
  if (length(perm_top) == 3) {
    inter <- Reduce(intersect, perm_top)
    uni <- Reduce(union, perm_top)
    perm_all3[p] <- if (length(uni) > 0) length(inter) / length(uni) else 0
  }
}

# Empirical p-value: fraction of permutations with Jaccard >= observed.
pair_dt[, perm_mean_jaccard := colMeans(perm_pair)]
pair_dt[, perm_p_value := vapply(seq_len(.N),
                                 function(k) mean(perm_pair[, k] >= jaccard[k]),
                                 numeric(1))]
pair_dt[, enrichment_fold := jaccard / pmax(perm_mean_jaccard, .Machine$double.eps)]

all3_observed <- summary_dt[method_A == "ALL3", jaccard]
all3_mean <- mean(perm_all3)
all3_p <- mean(perm_all3 >= all3_observed)
all3_fold <- all3_observed / max(all3_mean, .Machine$double.eps)
cat(sprintf("[perm] ALL3 observed=%.4f  null_mean=%.4f  fold=%.2fx  empirical p=%.4f\n",
            all3_observed, all3_mean, all3_fold, all3_p))

summary_dt[, perm_mean_jaccard := NA_real_]
summary_dt[, perm_p_value := NA_real_]
summary_dt[, enrichment_fold := NA_real_]
for (k in seq_len(nrow(pair_dt))) {
  summary_dt[method_A == pair_dt$method_A[k] & method_B == pair_dt$method_B[k],
             c("perm_mean_jaccard", "perm_p_value", "enrichment_fold") :=
               .(pair_dt$perm_mean_jaccard[k], pair_dt$perm_p_value[k],
                 pair_dt$enrichment_fold[k])]
}
summary_dt[method_A == "ALL3",
           c("perm_mean_jaccard", "perm_p_value", "enrichment_fold") :=
             .(all3_mean, all3_p, all3_fold)]

OUT_JC <- file.path(V3_DIR, "three_way_jaccard_v3.tsv")
fwrite(summary_dt, OUT_JC, sep = "\t")
cat(sprintf("[output] -> %s\n", OUT_JC))

OUT_NULL <- file.path(V3_DIR, "three_way_jaccard_perm_null_v3.tsv")
fwrite(data.table(perm_idx = seq_len(N_PERMS),
                  all3_jaccard = perm_all3),
       OUT_NULL, sep = "\t")
cat(sprintf("[output] permutation null -> %s\n", OUT_NULL))

cat("\n[summary] pairwise + three-way Jaccard with permutation null:\n")
print(summary_dt)
cat("[done] three-way concordance complete\n")
