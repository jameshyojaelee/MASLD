#!/usr/bin/env Rscript
# pillar_B_transferability.R
# ---------------------------------------------------------------------------
# Pillar B — Cross-cohort biological transferability of dream DEGs.
# For each held-out cohort in the 5 canonical mega cohorts:
#   * Train signature = dream LOO results from `dream_loo_{cohort}.csv`
#     (top 500 genes by |t-stat|; logFC sign as weight).
#   * Score test cohort with two methods:
#       (1) dot-product activation (z-score expression × signed weight) — same
#           as cross_study_auroc.R; classic batch-robust score.
#       (2) singscore signed-signature score (Foroutan 2018, BMC Bioinformatics)
#           — rank-based, fully batch-robust. Uses up + down sub-sets.
#   * Random-label baseline: shuffle test-cohort disease labels B_NULL=100
#     times, recompute AUROC under each scoring method.
# Out: results/audit_sensitivity/pillar_B_loco_prediction.csv
# ---------------------------------------------------------------------------

t0 <- proc.time()
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(pROC); library(yaml)
})
# Rank-based signed-signature score uses GSVA::gsva(method='ssgsea') applied
# to up + down DEG sub-sets separately, then signed_score = up - down.
# (Equivalent in spirit to singscore but stays within the rnaseq env.)
have_gsva <- requireNamespace("GSVA", quietly = TRUE)
if (have_gsva) suppressPackageStartupMessages(library(GSVA))

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
LOO_DIR <- file.path(RDIR, "loo_cv")
OUT_DIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TOP_N    <- as.integer(Sys.getenv("PILLAR_B_TOPN", "500"))
B_NULL   <- as.integer(Sys.getenv("PILLAR_B_BNULL", "100"))
SET_SEED <- 42L

cat(sprintf("=== Pillar B: LOCO disease prediction (top-%d, B_null=%d) ===\n",
            TOP_N, B_NULL))

# --- Identify the 5 canonical mega cohorts (yaml include_in_mega) ---
ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("Mega cohorts (k =", length(mega_cohorts), "):", paste(mega_cohorts, collapse = ", "), "\n")

# --- Expression data (log-CPM) ---
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
keep <- dge$samples$dataset %in% mega_cohorts
dge <- dge[, keep]
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
sm <- data.table(
  sample_id = colnames(dge),
  dataset   = dge$samples$dataset,
  group_binary = dge$samples$group_binary)

# ---------------------------------------------------------------------------
# Scoring helpers
# ---------------------------------------------------------------------------

# Method 1: dot-product activation (matches cross_study_auroc.R)
score_dotprod <- function(sig_dt, expr, samples) {
  sig <- sig_dt[order(-abs(t))][1:min(TOP_N, nrow(sig_dt))]
  avail <- intersect(sig$gene, rownames(expr))
  if (length(avail) < 10) return(rep(NA_real_, length(samples)))
  sig <- sig[gene %in% avail]
  w <- setNames(sig$logFC, sig$gene)
  m <- expr[avail, samples, drop = FALSE]
  z <- t(scale(t(m)));  z[is.na(z)] <- 0
  as.numeric(t(z) %*% w[avail])
}

# Method 2: signed-signature ssGSEA (GSVA) — rank-based, batch-robust.
# Score = ssGSEA(up_DEGs) - ssGSEA(down_DEGs).
score_signed_ssgsea <- function(sig_dt, expr, samples) {
  if (!have_gsva) return(rep(NA_real_, length(samples)))
  sig <- sig_dt[order(-abs(t))][1:min(TOP_N, nrow(sig_dt))]
  up_genes   <- intersect(sig[logFC > 0, gene], rownames(expr))
  down_genes <- intersect(sig[logFC < 0, gene], rownames(expr))
  if (length(up_genes) < 5 || length(down_genes) < 5) return(rep(NA_real_, length(samples)))
  m <- expr[, samples, drop = FALSE]
  # GSVA 1.50+ uses parameter objects; fall back to legacy API otherwise.
  ssg <- tryCatch({
    if (exists("ssgseaParam", where = asNamespace("GSVA"))) {
      param <- GSVA::ssgseaParam(exprData = m, geneSets = list(up = up_genes, down = down_genes))
      GSVA::gsva(param, verbose = FALSE)
    } else {
      GSVA::gsva(m, list(up = up_genes, down = down_genes), method = "ssgsea", verbose = FALSE)
    }
  }, error = function(e) { cat("    GSVA error:", conditionMessage(e), "\n"); NULL })
  if (is.null(ssg)) return(rep(NA_real_, length(samples)))
  as.numeric(ssg["up", ] - ssg["down", ])
}

# AUROC wrapper
compute_auc <- function(scores, labels) {
  if (length(unique(labels)) < 2 || all(is.na(scores))) return(NA_real_)
  tryCatch({
    ro <- roc(response = as.numeric(labels == "Disease"), predictor = scores, quiet = TRUE)
    as.numeric(auc(ro))
  }, error = function(e) NA_real_)
}

# Random-label null
random_baseline <- function(scores, labels, B = B_NULL) {
  if (length(unique(labels)) < 2 || all(is.na(scores))) return(c(NA_real_, NA_real_))
  set.seed(SET_SEED)
  null <- vapply(seq_len(B), function(i) {
    compute_auc(scores, sample(labels))
  }, numeric(1))
  c(mean(null, na.rm = TRUE), sd(null, na.rm = TRUE))
}

# ---------------------------------------------------------------------------
# Loop over held-out cohorts
# ---------------------------------------------------------------------------
out <- data.table(
  held_out = character(0), method = character(0),
  auroc = numeric(0), null_mean = numeric(0), null_sd = numeric(0),
  z_vs_null = numeric(0), n_test = integer(0),
  n_disease = integer(0), n_control = integer(0))

for (cohort in mega_cohorts) {
  loo_f <- file.path(LOO_DIR, paste0("dream_loo_", cohort, ".csv"))
  if (!file.exists(loo_f)) {
    cat(sprintf("  %s: SKIP (missing %s)\n", cohort, basename(loo_f))); next
  }
  sig_dt <- fread(loo_f)
  test_samples <- sm[dataset == cohort, sample_id]
  test_labels  <- sm[dataset == cohort, group_binary]
  if (length(unique(test_labels)) < 2) {
    cat(sprintf("  %s: SKIP (single-class test cohort)\n", cohort)); next
  }

  # Method 1: dot-product
  s1 <- score_dotprod(sig_dt, logcpm, test_samples)
  a1 <- compute_auc(s1, test_labels)
  b1 <- random_baseline(s1, test_labels)

  # Method 2: signed ssGSEA
  s2 <- score_signed_ssgsea(sig_dt, logcpm, test_samples)
  a2 <- compute_auc(s2, test_labels)
  b2 <- random_baseline(s2, test_labels)

  cat(sprintf("  %s: dot-product AUROC=%.3f null=%.3f±%.3f  |  signed_ssGSEA AUROC=%.3f null=%.3f±%.3f  (n=%d, D=%d, C=%d)\n",
              cohort, a1, b1[1], b1[2], a2, b2[1], b2[2],
              length(test_samples), sum(test_labels == "Disease"), sum(test_labels == "Control")))

  out <- rbind(out,
    data.table(held_out = cohort, method = "dotprod",
               auroc = round(a1, 4), null_mean = round(b1[1], 4), null_sd = round(b1[2], 4),
               z_vs_null = round(if (!is.na(b1[2]) && b1[2] > 0) (a1 - b1[1]) / b1[2] else NA_real_, 3),
               n_test = length(test_samples),
               n_disease = sum(test_labels == "Disease"),
               n_control = sum(test_labels == "Control")),
    data.table(held_out = cohort, method = "signed_ssgsea",
               auroc = round(a2, 4), null_mean = round(b2[1], 4), null_sd = round(b2[2], 4),
               z_vs_null = round(if (!is.na(b2[2]) && b2[2] > 0) (a2 - b2[1]) / b2[2] else NA_real_, 3),
               n_test = length(test_samples),
               n_disease = sum(test_labels == "Disease"),
               n_control = sum(test_labels == "Control")))
}

fwrite(out, file.path(OUT_DIR, "pillar_B_loco_prediction.csv"))
cat("\n=== Pillar B summary ===\n"); print(out)

cat(sprintf("\nMean AUROC by method:\n"))
print(out[, .(mean_auroc = round(mean(auroc, na.rm = TRUE), 4),
              min_auroc  = round(min(auroc, na.rm = TRUE), 4),
              max_auroc  = round(max(auroc, na.rm = TRUE), 4),
              mean_null  = round(mean(null_mean, na.rm = TRUE), 4)), by = method])

cat(sprintf("\nElapsed %.1f min\n", (proc.time() - t0)["elapsed"] / 60))
cat("Done Pillar B.\n")
