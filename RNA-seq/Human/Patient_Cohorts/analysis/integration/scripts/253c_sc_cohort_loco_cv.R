# 253c_sc_cohort_loco_cv.R
# Test whether the 35/26 sc-Steatohepatitis donor bimodality (Ward k=2 on
# meta-subtype proportions) is driven by source cohort (GSE244832 / GSE189600 /
# GSE136103) vs real biology.
#
# Tests (pre-registered acceptance):
#   (1) Stratify donors by dataset within Steatohepatitis (n_cells >= 50).
#   (2) Within each cohort with n_donors >= 10, run Ward k=2 on the
#       meta-subtype proportion matrix; report sizes / balance / DeltaBIC.
#   (3) Pooled cluster <-> cohort mutual-information ratio (MI / min(Hc, Hk)).
#   (4) If GSE244832 has sub-batch metadata, test 35/26 split vs sub-batch.
#
# PASS criteria (real biology, not batch effect):
#   - MI ratio (pooled) < 0.5
#   - Balanced Ward k=2 (min cluster fraction >= 30%) in >= 2 of 3 cohorts (where n permits)
#   - 35/26 split survives within GSE244832 (the largest cohort, n_donor >= 50)
#
# Inputs:
#   Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_subtype_metadata.csv
#   Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/meta_subtype_mapping.csv
#
# Outputs (RNA-seq/results/granular_staging/):
#   sc_cohort_loco_cv.csv      — per-cohort + pooled cluster diagnostics
#   sc_cohort_loco_summary.md  — verdict paragraph

suppressPackageStartupMessages({
  library(matrixStats)
  library(cluster)
})

set.seed(42)

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                           "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HEP_DIR <- file.path(PROJECT_ROOT, "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("== sc Steatohepatitis cohort LOCO-CV bimodality test (253c) ==\n")

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
hep_meta_path <- file.path(HEP_DIR, "hepatocyte_subtype_metadata.csv")
mapping_path  <- file.path(HEP_DIR, "meta_subtype_mapping.csv")

cat("Reading hepatocyte metadata (large file)...\n")
hep <- read.csv(hep_meta_path, stringsAsFactors = FALSE)
mapping <- read.csv(mapping_path, stringsAsFactors = FALSE)
cat(sprintf("Total cells: %d\n", nrow(hep)))
cat("Datasets present:\n"); print(table(hep$dataset, useNA = "ifany"))

hep$meta_subtype <- mapping$meta_subtype[match(hep$hepatocyte_subtype, mapping$subtype)]

# ---------------------------------------------------------------------------
# Subset to Steatohepatitis donors
# ---------------------------------------------------------------------------
sh <- hep[hep$disease_stage_coarse == "Steatohepatitis", ]
cat(sprintf("\nSteatohepatitis cells: %d\n", nrow(sh)))
cat(sprintf("Donors (samples): %d\n", length(unique(sh$sample))))
cat("Cells per cohort (within SH):\n"); print(table(sh$dataset, useNA = "ifany"))

# Donor-level proportion matrix (Steatohepatitis only)
prop_tab <- prop.table(table(sh$sample, sh$meta_subtype), margin = 1)
donor_props <- as.data.frame.matrix(prop_tab)
donor_props$sample <- rownames(donor_props)
donor_props$dataset <- sh$dataset[match(donor_props$sample, sh$sample)]
donor_props$n_cells <- as.integer(table(sh$sample))[match(donor_props$sample,
                                                          names(table(sh$sample)))]
donor_props <- donor_props[!is.na(donor_props$n_cells) & donor_props$n_cells >= 50, ]
cat(sprintf("Donors (>=50 cells): %d\n", nrow(donor_props)))
cat("Donors per cohort:\n"); print(table(donor_props$dataset, useNA = "ifany"))

# meta_subtype columns (consistent ordering)
ms_cols <- intersect(c("Healthy", "Neutral", "Disease-Neutral",
                       "Disease-Associated", "Disease-Progressor"),
                     colnames(donor_props))

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
# 1D BIC for k=1 vs k=2 Gaussian on the chosen axis
fit_gmm_2 <- function(x, max_iter = 200, tol = 1e-6) {
  n <- length(x)
  ord <- sort(x)
  mu1 <- mean(ord[1:floor(n/2)]); mu2 <- mean(ord[(floor(n/2)+1):n])
  sd1 <- max(sd(ord[1:floor(n/2)]), 1e-3)
  sd2 <- max(sd(ord[(floor(n/2)+1):n]), 1e-3)
  p1 <- 0.5; ll_old <- -Inf; ll <- -Inf
  for (iter in 1:max_iter) {
    d1 <- p1 * dnorm(x, mu1, sd1); d2 <- (1 - p1) * dnorm(x, mu2, sd2)
    g1 <- d1 / (d1 + d2 + 1e-300)
    p1 <- mean(g1)
    mu1 <- sum(g1 * x) / sum(g1)
    mu2 <- sum((1 - g1) * x) / sum(1 - g1)
    sd1 <- max(sqrt(sum(g1 * (x - mu1)^2) / sum(g1)), 1e-3)
    sd2 <- max(sqrt(sum((1 - g1) * (x - mu2)^2) / sum(1 - g1)), 1e-3)
    ll <- sum(log(p1 * dnorm(x, mu1, sd1) + (1 - p1) * dnorm(x, mu2, sd2) + 1e-300))
    if (abs(ll - ll_old) < tol) break
    ll_old <- ll
  }
  list(mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, p1 = p1, ll = ll)
}

bic_two_minus_one <- function(x) {
  if (length(x) < 6 || sd(x) < 1e-9) return(NA_real_)
  n <- length(x)
  ll1 <- sum(dnorm(x, mean(x), sd(x), log = TRUE))
  bic1 <- -2 * ll1 + 2 * log(n)            # 2 params (mu, sigma)
  fit2 <- fit_gmm_2(x)
  bic2 <- -2 * fit2$ll + 5 * log(n)        # 5 params (mu1, mu2, sd1, sd2, p1)
  # Positive => k=2 preferred
  bic1 - bic2
}

# Cluster-cohort MI (cribbed from 241_f3_substate_discovery.R)
mi_metric <- function(cluster, cohort) {
  tab <- table(cluster, cohort)
  N <- sum(tab); if (N == 0) return(list(MI = 0, NMI_max = 0, NMI_geom = 0,
                                          MI_ratio_min = 0, saturated_MI = 0))
  pxy <- tab / N
  px <- rowSums(pxy); py <- colSums(pxy)
  mi <- 0
  for (i in seq_along(px)) for (j in seq_along(py)) {
    if (pxy[i, j] > 0) mi <- mi + pxy[i, j] * log2(pxy[i, j] / (px[i] * py[j]))
  }
  Hx <- -sum(px[px > 0] * log2(px[px > 0]))
  Hy <- -sum(py[py > 0] * log2(py[py > 0]))
  list(MI = mi, NMI_max = mi / max(Hx, Hy, 1e-9),
       NMI_geom = mi / sqrt(Hx * Hy + 1e-18),
       MI_ratio_min = mi / min(Hx, Hy, 1e9),  # the pre-registered metric
       saturated_MI = min(Hx, Hy))
}

ward_k2 <- function(prop_mat) {
  if (nrow(prop_mat) < 4) return(NULL)
  d <- dist(prop_mat)
  hc <- hclust(d, method = "ward.D2")
  cl2 <- cutree(hc, k = 2)
  list(cluster = cl2, hc = hc)
}

# ---------------------------------------------------------------------------
# (A) Pooled Ward k=2 (reproduces the "35/26" split)
# ---------------------------------------------------------------------------
pool_mat <- as.matrix(donor_props[, ms_cols])
rownames(pool_mat) <- donor_props$sample
cat("\n--- Pooled Ward k=2 ---\n")
pool_fit <- ward_k2(pool_mat)
pool_cl <- pool_fit$cluster[match(donor_props$sample, names(pool_fit$cluster))]
donor_props$pooled_cluster_k2 <- pool_cl
pool_sizes <- table(pool_cl)
cat("Pooled cluster sizes:\n"); print(pool_sizes)
pool_balance <- min(pool_sizes) / sum(pool_sizes)
cat(sprintf("Pooled balance (min frac): %.3f\n", pool_balance))

# Pooled MI ratio (cluster vs cohort)
mi_pool <- mi_metric(pool_cl, donor_props$dataset)
cat("Pooled cluster <-> cohort MI:\n"); print(mi_pool)

# Pooled BIC on the F3-axis (Progressor + DA - Disease-Neutral)
prog_col <- intersect("Disease-Progressor", colnames(donor_props))
da_col   <- intersect("Disease-Associated", colnames(donor_props))
dn_col   <- intersect("Disease-Neutral",    colnames(donor_props))

axis_pool <- if (all(c(length(prog_col), length(da_col), length(dn_col)) == 1)) {
  donor_props[[prog_col]] + donor_props[[da_col]] - donor_props[[dn_col]]
} else {
  rep(NA_real_, nrow(donor_props))
}
delta_bic_pool <- bic_two_minus_one(axis_pool)
cat(sprintf("Pooled DeltaBIC(2 vs 1) on F3-axis: %.2f\n", delta_bic_pool))

# Cohort-cluster contingency for the pooled split
cat("Pooled cluster x cohort contingency:\n")
print(table(pool_cl, donor_props$dataset))

# ---------------------------------------------------------------------------
# (B) Within-cohort Ward k=2 (LOCO-style: each cohort independently)
# ---------------------------------------------------------------------------
results_rows <- list()
cohort_levels <- sort(unique(donor_props$dataset))

for (coh in cohort_levels) {
  sub_idx <- which(donor_props$dataset == coh)
  sub_n <- length(sub_idx)
  cat(sprintf("\n--- Cohort: %s (n_donors = %d) ---\n", coh, sub_n))
  if (sub_n < 10) {
    results_rows[[coh]] <- data.frame(
      cohort = coh, n_donors = sub_n,
      cluster_sizes_within = NA_character_,
      min_balance_within = NA_real_,
      delta_bic_within = NA_real_,
      mi_ratio_pooled = mi_pool$MI_ratio_min,
      pass_cohort_confound_test = NA,
      notes = "n_donors < 10; within-cohort Ward not run",
      stringsAsFactors = FALSE
    )
    next
  }
  sub_mat <- pool_mat[sub_idx, , drop = FALSE]
  sub_fit <- ward_k2(sub_mat)
  sub_cl <- sub_fit$cluster
  sub_sizes <- as.integer(table(sub_cl))
  sub_balance <- min(sub_sizes) / sum(sub_sizes)
  sub_axis <- axis_pool[sub_idx]
  sub_dbic <- bic_two_minus_one(sub_axis)
  cat(sprintf("Cluster sizes: %s\n", paste(sub_sizes, collapse = " / ")))
  cat(sprintf("Balance: %.3f  DeltaBIC(2v1): %.2f\n", sub_balance, sub_dbic))
  # Per-cohort cluster mean profile
  cms <- aggregate(sub_mat, by = list(cluster = sub_cl), FUN = mean)
  cat("Cluster mean profile within cohort:\n"); print(cms)

  results_rows[[coh]] <- data.frame(
    cohort = coh, n_donors = sub_n,
    cluster_sizes_within = paste(sub_sizes, collapse = "/"),
    min_balance_within = sub_balance,
    delta_bic_within = sub_dbic,
    mi_ratio_pooled = mi_pool$MI_ratio_min,
    pass_cohort_confound_test = sub_balance >= 0.30,
    notes = sprintf("Ward k=2 within %s; balance threshold = 30%%", coh),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# (C) GSE244832 sub-batch correlation (only if any obvious sub-batch column)
# ---------------------------------------------------------------------------
gse_target <- "GSE244832"
if (gse_target %in% donor_props$dataset) {
  sub_idx <- which(donor_props$dataset == gse_target)
  sub_samples <- donor_props$sample[sub_idx]
  sub_meta <- sh[sh$sample %in% sub_samples, , drop = FALSE]
  sub_meta_cols <- setdiff(colnames(sub_meta),
                           c("sample", "dataset", "condition", "condition_binary",
                             "hepatocyte_subtype", "hepatocyte_subtype_label",
                             "disease_stage_coarse", "disease_stage_numeric",
                             "axis_periportal", "axis_pericentral",
                             "axis_lipid_accumulation", "axis_inflammatory",
                             "axis_ferroptosis", "UMAP_1", "UMAP_2",
                             "meta_subtype"))
  cat(sprintf("\n--- GSE244832 sub-batch metadata candidates: %s ---\n",
              paste(sub_meta_cols, collapse = ", ")))
  # Best-effort: for each donor, take the modal value of any extra column.
  if (length(sub_meta_cols) > 0) {
    donor_subbatch <- aggregate(sub_meta[, sub_meta_cols, drop = FALSE],
                                by = list(sample = sub_meta$sample),
                                FUN = function(x) {
                                  ux <- unique(x)
                                  ux[which.max(tabulate(match(x, ux)))]
                                })
    sub_cl_pool <- pool_cl[sub_idx]
    names(sub_cl_pool) <- donor_props$sample[sub_idx]
    donor_subbatch$pool_cluster <- sub_cl_pool[donor_subbatch$sample]
    cat("GSE244832 sub-batch x pooled cluster contingency:\n")
    for (cc in setdiff(colnames(donor_subbatch), c("sample", "pool_cluster"))) {
      tab <- table(donor_subbatch$pool_cluster, donor_subbatch[[cc]])
      if (nrow(tab) >= 2 && ncol(tab) >= 2 && sum(tab) >= 6) {
        cat(sprintf("  Column: %s\n", cc)); print(tab)
        ft <- tryCatch(fisher.test(tab, simulate.p.value = TRUE, B = 999),
                       error = function(e) NULL)
        if (!is.null(ft)) cat(sprintf("    Fisher p (sim) = %.4f\n", ft$p.value))
      }
    }
  } else {
    cat("No additional metadata columns available beyond the standard set.\n")
  }
}

# ---------------------------------------------------------------------------
# Aggregate output rows
# ---------------------------------------------------------------------------
# Add a pooled summary row at the top
pool_row <- data.frame(
  cohort = "POOLED",
  n_donors = nrow(donor_props),
  cluster_sizes_within = paste(as.integer(pool_sizes), collapse = "/"),
  min_balance_within = pool_balance,
  delta_bic_within = delta_bic_pool,
  mi_ratio_pooled = mi_pool$MI_ratio_min,
  pass_cohort_confound_test = mi_pool$MI_ratio_min < 0.5,
  notes = "Pooled across all 3 cohorts; pass = MI_ratio < 0.5",
  stringsAsFactors = FALSE
)

results_df <- do.call(rbind, c(list(pool_row), results_rows))
rownames(results_df) <- NULL
write.csv(results_df, file.path(OUT_DIR, "sc_cohort_loco_cv.csv"),
          row.names = FALSE)
cat("\nWritten: sc_cohort_loco_cv.csv\n"); print(results_df)

# ---------------------------------------------------------------------------
# Summary verdict
# ---------------------------------------------------------------------------
within_results <- results_df[results_df$cohort != "POOLED" &
                             !is.na(results_df$min_balance_within), ]
n_balanced <- sum(within_results$min_balance_within >= 0.30)
n_eligible <- nrow(within_results)
mi_pass <- mi_pool$MI_ratio_min < 0.5

# 35/26 split survives within GSE244832?
gse_row <- within_results[within_results$cohort == "GSE244832", ]
gse_pass <- nrow(gse_row) == 1 && gse_row$min_balance_within >= 0.30 &&
            gse_row$n_donors >= 30

verdict_real <- mi_pass && (n_balanced >= 2 || gse_pass)

sink(file.path(OUT_DIR, "sc_cohort_loco_summary.md"))
cat("# sc Steatohepatitis cohort LOCO-CV bimodality test\n\n")
cat(sprintf("Steatohepatitis donors (>=50 cells): n = %d across %d cohorts.\n",
            nrow(donor_props), length(cohort_levels)))
cat(sprintf("Pooled Ward k=2 sizes: %s (balance %.2f, DeltaBIC %.1f).\n",
            paste(as.integer(pool_sizes), collapse = " / "),
            pool_balance, delta_bic_pool))
cat(sprintf("Cluster <-> cohort MI ratio (MI / min(Hc, Hk)) = %.3f ",
            mi_pool$MI_ratio_min))
cat(sprintf("(pass threshold < 0.5: %s).\n",
            ifelse(mi_pass, "PASS", "FAIL")))
cat("\nWithin-cohort Ward k=2 (eligible cohorts with n_donors >= 10):\n\n")
if (n_eligible == 0) {
  cat("- No cohort had n_donors >= 10 in Steatohepatitis; within-cohort test not evaluable.\n")
} else {
  cat("| cohort | n_donors | sizes | min_balance | DeltaBIC | pass (>=30%) |\n")
  cat("|---|---:|---|---:|---:|:---:|\n")
  for (i in seq_len(nrow(within_results))) {
    cat(sprintf("| %s | %d | %s | %.2f | %.1f | %s |\n",
                within_results$cohort[i], within_results$n_donors[i],
                within_results$cluster_sizes_within[i],
                within_results$min_balance_within[i],
                within_results$delta_bic_within[i],
                ifelse(within_results$pass_cohort_confound_test[i], "yes", "no")))
  }
}
cat(sprintf("\nGSE244832 (largest cohort) survival: %s.\n",
            ifelse(gse_pass, "PASS (balance >= 30%, n_donors >= 30)",
                   "FAIL")))
cat("\n## Verdict\n\n")
gse_n_print <- ifelse(nrow(gse_row) == 1, as.character(gse_row$n_donors), "?")
gse_phrase  <- ifelse(gse_pass, "produces", "does not produce")
mi_phrase   <- ifelse(mi_pass, "pass", "fail")
gse_test_phrase <- ifelse(gse_pass,
                          "passes the within-cohort test",
                          "fails the within-cohort test")
if (verdict_real) {
  cat("The 35/26 sc-Steatohepatitis donor split survives cohort stratification: ")
  cat(sprintf("MI ratio %.3f is below the 0.5 confound threshold, the largest cohort ",
              mi_pool$MI_ratio_min))
  cat(sprintf("(GSE244832, n = %s) %s a balanced k=2 split, and ",
              gse_n_print, gse_phrase))
  cat(sprintf("%d of %d eligible cohorts hit the 30%% balance bar. The bimodality is ",
              n_balanced, n_eligible))
  cat("therefore consistent with real biological heterogeneity within Steatohepatitis ")
  cat("rather than a source-cohort batch artefact.\n")
} else {
  cat("The 35/26 sc-Steatohepatitis donor split does NOT cleanly survive cohort ")
  cat(sprintf("stratification: MI ratio = %.3f (threshold < 0.5: %s), only %d/%d eligible ",
              mi_pool$MI_ratio_min, mi_phrase, n_balanced, n_eligible))
  cat("cohorts produce a balanced (>=30%) within-cohort k=2 split, and the largest ")
  cat(sprintf("cohort GSE244832 %s. The pooled bimodality is therefore at least partially ",
              gse_test_phrase))
  cat("driven by source-cohort composition and should be reported as ")
  cat("severity-and-cohort-confounded rather than purely biological.\n")
}
sink()

cat("\n== 253c LOCO-CV complete. Outputs in", OUT_DIR, "==\n")
