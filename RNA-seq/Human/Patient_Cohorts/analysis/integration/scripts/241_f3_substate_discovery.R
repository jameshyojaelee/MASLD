# 241_f3_substate_discovery.R
# Phase 1.1 — Per-cohort + pooled F3 sub-state discovery.
#
# Inputs:
#   results/integration/merged_dge.rds  (DGEList, 34453 × 1444)
#   results/integration/meta_matched.rds  (sample metadata)
#
# Outputs (to results/granular_staging/):
#   f3_substate_per_cohort_labels.csv  — per-sample sub-state per cohort
#   f3_substate_pooled_labels.csv      — per-sample sub-state from pooled clustering
#   f3_substate_stability.csv          — per-cohort + pooled cophenetic / silhouette / PAC
#   f3_substate_cohort_mi.csv          — cluster ↔ cohort mutual information
#   per_cohort/{COHORT}_consensus_k*.pdf  — consensus heatmap PDFs
#   pooled_consensus_k*.pdf

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ConsensusClusterPlus)
  library(matrixStats)
  library(cluster)        # silhouette
})

set.seed(42)

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")
FIG_DIR <- file.path(PROJECT_ROOT, "figures/supplementary/figS_granular_staging")
PER_COHORT_FIG_DIR <- file.path(FIG_DIR, "per_cohort_consensus")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PER_COHORT_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("== Phase 1.1 F3 Sub-state Discovery ==\n")
cat("Output:", OUT_DIR, "\n\n")

# ---------------------------------------------------------------------------
# Load and subset
# ---------------------------------------------------------------------------
dge <- readRDS(file.path(INT_DIR, "merged_dge.rds"))
meta <- readRDS(file.path(INT_DIR, "meta_matched.rds"))

stopifnot(identical(colnames(dge), meta$sample_id) || all(colnames(dge) %in% meta$sample_id))
meta <- meta[match(colnames(dge), meta$sample_id), ]
stopifnot(identical(colnames(dge), meta$sample_id))

f3_idx <- which(!is.na(meta$fibrosis_stage) & meta$fibrosis_stage == 3)
cat(sprintf("F3 samples: %d (across %d cohorts)\n",
            length(f3_idx), length(unique(meta$dataset[f3_idx]))))
f3_table <- table(meta$dataset[f3_idx])
print(f3_table)

DISCOVERY_COHORTS <- c("GSE213621", "GSE135251", "GSE193066")
VALIDATION_COHORTS <- setdiff(names(f3_table)[f3_table > 0], DISCOVERY_COHORTS)
cat("\nDiscovery (n>=30):", paste(DISCOVERY_COHORTS, collapse = ", "), "\n")
cat("Validation (n<30):", paste(VALIDATION_COHORTS, collapse = ", "), "\n\n")

# ---------------------------------------------------------------------------
# Helper: per-cohort consensus clustering + stability
# ---------------------------------------------------------------------------
voom_top_var <- function(dge_sub, n_top = 2000) {
  keep <- filterByExpr(dge_sub, group = factor(rep("F3", ncol(dge_sub))), min.count = 5)
  dge_sub <- dge_sub[keep, , keep.lib.sizes = FALSE]
  dge_sub <- calcNormFactors(dge_sub)
  v <- voom(dge_sub, design = NULL)
  expr <- v$E
  vars <- rowVars(expr)
  top <- order(vars, decreasing = TRUE)[seq_len(min(n_top, nrow(expr)))]
  expr[top, , drop = FALSE]
}

run_ccp <- function(expr, max_k = 4, title_prefix = "ccp", out_dir, plot_only_best = TRUE) {
  # ConsensusClusterPlus uses pearson by default; use spearman+ward for robustness.
  res <- ConsensusClusterPlus(
    d = expr,
    maxK = max_k,
    reps = 1000,
    pItem = 0.8,
    pFeature = 1,
    clusterAlg = "hc",
    distance = "spearman",
    seed = 42,
    title = title_prefix,
    plot = "pdf",
    writeTable = FALSE
  )
  # Move PDF (ConsensusClusterPlus writes to current dir under title prefix).
  # Filename includes the cohort/scope tag from `title_prefix` to avoid collisions
  # when multiple cohorts share `consensus.pdf` as default name.
  if (file.exists(title_prefix)) {
    pdfs <- list.files(title_prefix, pattern = "\\.pdf$", full.names = TRUE)
    tag <- basename(title_prefix)
    for (p in pdfs) {
      out_name <- paste0(tag, "_", basename(p))
      file.copy(p, file.path(out_dir, out_name), overwrite = TRUE)
    }
    unlink(title_prefix, recursive = TRUE)
  }
  res
}

cluster_quality <- function(ccp_res, expr, k_grid = 2:4) {
  # Compute cophenetic, PAC, silhouette per k.
  rows <- list()
  for (k in k_grid) {
    cls <- ccp_res[[k]]$consensusClass
    cm <- ccp_res[[k]]$consensusMatrix
    # PAC: proportion of consensus values in [0.1, 0.9]
    cm_off <- cm[upper.tri(cm)]
    pac <- mean(cm_off >= 0.1 & cm_off <= 0.9)
    # Cophenetic correlation: cor between consensus matrix (as distance) and the
    # cluster-derived ultrametric.
    coph <- tryCatch({
      d_cm <- as.dist(1 - cm)
      hc <- hclust(d_cm, method = "average")
      cor(cophenetic(hc), d_cm)
    }, error = function(e) NA_real_)
    # Silhouette over Spearman distance
    d_expr <- as.dist(1 - cor(expr, method = "spearman"))
    sil_mat <- tryCatch(silhouette(cls, d_expr), error = function(e) NULL)
    sil_avg <- if (!is.null(sil_mat)) mean(sil_mat[, "sil_width"]) else NA_real_
    rows[[as.character(k)]] <- data.frame(
      k = k, n_samples = length(cls),
      cophenetic = coph, pac = pac, silhouette_avg = sil_avg,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Per-cohort loop (discovery + validation tiers)
# ---------------------------------------------------------------------------
all_cohorts <- c(DISCOVERY_COHORTS, VALIDATION_COHORTS)
per_cohort_labels <- list()
per_cohort_stability <- list()

for (cohort in all_cohorts) {
  cohort_idx <- f3_idx[meta$dataset[f3_idx] == cohort]
  n <- length(cohort_idx)
  if (n < 6) {
    cat(sprintf("Skipping %s: n=%d too small\n", cohort, n))
    next
  }
  cat(sprintf("\n--- %s (n=%d) ---\n", cohort, n))
  dge_sub <- dge[, cohort_idx]
  expr <- voom_top_var(dge_sub, n_top = 2000)
  cat(sprintf("voom expression: %d genes × %d samples\n", nrow(expr), ncol(expr)))

  prefix <- file.path(tempdir(), paste0("ccp_", cohort))
  ccp_res <- tryCatch(
    run_ccp(expr, max_k = 4, title_prefix = prefix, out_dir = PER_COHORT_FIG_DIR),
    error = function(e) { cat("CCP error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(ccp_res)) next

  qual <- cluster_quality(ccp_res, expr)
  qual$cohort <- cohort
  qual$score <- qual$silhouette_avg + qual$cophenetic - qual$pac
  per_cohort_stability[[cohort]] <- qual

  best_k <- qual$k[which.max(qual$score)]
  cat(sprintf("Stability per k:\n")); print(qual)
  cat(sprintf("Best k for %s: %d\n", cohort, best_k))

  cls <- ccp_res[[best_k]]$consensusClass
  per_cohort_labels[[cohort]] <- data.frame(
    sample_id = colnames(expr),
    dataset = cohort,
    n_f3_in_cohort = n,
    k_chosen = best_k,
    cluster = paste0("F3.", as.character(cls)),
    stringsAsFactors = FALSE
  )
}

# Write per-cohort results
all_labels <- do.call(rbind, per_cohort_labels)
all_stability <- do.call(rbind, per_cohort_stability)
write.csv(all_labels, file.path(OUT_DIR, "f3_substate_per_cohort_labels.csv"), row.names = FALSE)
write.csv(all_stability, file.path(OUT_DIR, "f3_substate_per_cohort_stability.csv"), row.names = FALSE)
cat(sprintf("\nWrote per-cohort labels (%d rows) and stability (%d rows)\n",
            nrow(all_labels), nrow(all_stability)))

# ---------------------------------------------------------------------------
# Pooled analysis (all 278 F3 samples, batch-corrected)
# ---------------------------------------------------------------------------
cat("\n--- Pooled F3 (n=", length(f3_idx), ") ---\n", sep = "")
dge_pool <- dge[, f3_idx]
keep <- filterByExpr(dge_pool, group = factor(rep("F3", ncol(dge_pool))), min.count = 5)
dge_pool <- dge_pool[keep, , keep.lib.sizes = FALSE]
dge_pool <- calcNormFactors(dge_pool)
v_pool <- voom(dge_pool, design = NULL)
expr_pool_raw <- v_pool$E

# Batch correction: cohort as batch, no biological covariate within F3
expr_pool <- removeBatchEffect(expr_pool_raw, batch = factor(meta$dataset[f3_idx]))
vars <- rowVars(expr_pool)
top <- order(vars, decreasing = TRUE)[seq_len(min(2000, nrow(expr_pool)))]
expr_pool_top <- expr_pool[top, ]
cat(sprintf("Pooled expression: %d genes × %d samples (after batch correction)\n",
            nrow(expr_pool_top), ncol(expr_pool_top)))

prefix <- file.path(tempdir(), "ccp_pooled")
ccp_pool <- run_ccp(expr_pool_top, max_k = 4, title_prefix = prefix, out_dir = FIG_DIR)
# Canonical filename for the pooled diagnostic
auto_pdf <- file.path(FIG_DIR, paste0(basename(prefix), "_consensus.pdf"))
canonical_pdf <- file.path(FIG_DIR, "figS_granular_consensus_pooled.pdf")
if (file.exists(auto_pdf)) {
  file.rename(auto_pdf, canonical_pdf)
}
qual_pool <- cluster_quality(ccp_pool, expr_pool_top)
qual_pool$cohort <- "POOLED"
qual_pool$score <- qual_pool$silhouette_avg + qual_pool$cophenetic - qual_pool$pac
best_k_pool <- qual_pool$k[which.max(qual_pool$score)]
cat("Pooled stability:\n"); print(qual_pool)
cat(sprintf("Best k pooled: %d\n", best_k_pool))

cls_pool <- ccp_pool[[best_k_pool]]$consensusClass
pooled_labels <- data.frame(
  sample_id = colnames(expr_pool_top),
  dataset = meta$dataset[f3_idx],
  k_chosen = best_k_pool,
  cluster_pooled = paste0("F3.pool.", as.character(cls_pool)),
  stringsAsFactors = FALSE
)
write.csv(pooled_labels, file.path(OUT_DIR, "f3_substate_pooled_labels.csv"), row.names = FALSE)
write.csv(rbind(all_stability, qual_pool), file.path(OUT_DIR, "f3_substate_stability.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Cluster ↔ cohort mutual information diagnostic
# ---------------------------------------------------------------------------
mi_metric <- function(cluster, cohort) {
  tab <- table(cluster, cohort)
  N <- sum(tab)
  pxy <- tab / N
  px <- rowSums(pxy)
  py <- colSums(pxy)
  mi <- 0
  for (i in seq_along(px)) {
    for (j in seq_along(py)) {
      if (pxy[i, j] > 0) {
        mi <- mi + pxy[i, j] * log2(pxy[i, j] / (px[i] * py[j]))
      }
    }
  }
  Hx <- -sum(px[px > 0] * log2(px[px > 0]))
  Hy <- -sum(py[py > 0] * log2(py[py > 0]))
  list(MI = mi, NMI_max = mi / max(Hx, Hy), NMI_geom = mi / sqrt(Hx * Hy),
       saturated_MI = min(Hx, Hy))
}

mi <- mi_metric(cls_pool, meta$dataset[f3_idx])
mi_df <- data.frame(
  test = "cluster_cohort_MI_pooled",
  k = best_k_pool,
  MI = mi$MI,
  NMI_max = mi$NMI_max,
  NMI_geom = mi$NMI_geom,
  saturated_MI = mi$saturated_MI,
  cohort_confound_ratio = mi$MI / mi$saturated_MI,
  passes_gate = mi$MI / mi$saturated_MI < 0.5,
  stringsAsFactors = FALSE
)
write.csv(mi_df, file.path(OUT_DIR, "f3_substate_cohort_mi.csv"), row.names = FALSE)
cat("\nCluster-cohort MI diagnostic:\n"); print(mi_df)

cat("\n== Phase 1.1 discovery complete. Inspect outputs in", OUT_DIR, "==\n")
