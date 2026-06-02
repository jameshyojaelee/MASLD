# 245_phase1_gates_check.R
# Phase 1.7 — Aggregate all gates into a single GO/NO-GO decision CSV.
#
# Gates (per plan):
#   G1. Per-cohort recovery: ≥3 of 4 discovery cohorts show stable k≥2 partition
#       (cophenetic >0.85, silhouette >0.25)
#   G2. Bimodality: BIC favors k=2 mixture by ≥10 in pooled FPI density,
#       OR dip-MC p<0.05 at sep ≥1σ
#   G3. LOCO-CV: held-out cluster prediction AUROC > 0.65
#       (this gate is computed in 245 itself via held-out cohort prediction)
#   G4. Cohort-confounding: cluster ↔ cohort NMI < 0.5 of saturated
#   G5. Effect-size floor: F3 sub-state effect ≥ 30% of F1→F2 effect
#
# Inputs (must exist):
#   results/granular_staging/f3_substate_per_cohort_stability.csv (from 241)
#   results/granular_staging/f3_substate_pooled_labels.csv (from 241)
#   results/granular_staging/f3_substate_cohort_mi.csv (from 241)
#   results/granular_staging/continuous_fibrosis_index_diagnostics.csv (from 242)
#   results/granular_staging/f3_substate_signature_concordance.csv (from 243)
#   results/granular_staging/transition_effect_sizes.csv (from 244)
#
# Output:
#   results/granular_staging/phase1_gates.csv  — per-gate pass/fail + decision

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(matrixStats)
  library(pROC)
})

set.seed(42)

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")

cat("== Phase 1.7 Gate evaluation ==\n")

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

read_csv_safe <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Missing: %s — run prerequisite script first", path))
  }
  read.csv(path, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# G1 — Per-cohort stability AND cluster balance (min cluster fraction ≥ 15%)
# ---------------------------------------------------------------------------
cat("\nG1: Per-cohort stability + balance\n")
stab <- read_csv_safe(file.path(OUT_DIR, "f3_substate_per_cohort_stability.csv"))
per_cohort_lbl <- read_csv_safe(file.path(OUT_DIR, "f3_substate_per_cohort_labels.csv"))
discovery <- c("GSE213621", "GSE135251", "GSE193066")

g1_per_cohort <- list()
for (coh in discovery) {
  s <- stab[stab$cohort == coh, ]
  if (nrow(s) == 0) {
    g1_per_cohort[[coh]] <- list(pass = FALSE, reason = "no stability data")
    next
  }
  best <- s[which.max(s$silhouette_avg + s$cophenetic - s$pac), ]
  # Compute min cluster fraction at chosen k from labels
  cohort_lbls <- per_cohort_lbl[per_cohort_lbl$dataset == coh, ]
  if (nrow(cohort_lbls) > 0) {
    cluster_table <- table(cohort_lbls$cluster)
    min_frac <- min(cluster_table) / sum(cluster_table)
    n_clusters <- length(cluster_table)
  } else {
    min_frac <- NA_real_; n_clusters <- NA_integer_
  }
  stability_pass <- !is.na(best$cophenetic) && best$cophenetic > 0.85 &&
                    !is.na(best$silhouette_avg) && best$silhouette_avg > 0.25
  balance_pass <- !is.na(min_frac) && min_frac >= 0.15
  pass <- stability_pass && balance_pass
  g1_per_cohort[[coh]] <- list(pass = pass,
                                cophenetic = best$cophenetic,
                                silhouette = best$silhouette_avg,
                                k = best$k,
                                min_cluster_fraction = min_frac,
                                n_clusters = n_clusters,
                                stability_pass = stability_pass,
                                balance_pass = balance_pass)
}
n_pass_g1 <- sum(sapply(g1_per_cohort, function(x) x$pass))
g1_overall <- n_pass_g1 >= 3
cat(sprintf("  Cohorts passing per-cohort stability: %d / 4 (need ≥ 3)  [%s]\n",
            n_pass_g1, ifelse(g1_overall, "PASS", "FAIL")))
for (coh in discovery) {
  x <- g1_per_cohort[[coh]]
  cat(sprintf("  %s: coph=%s sil=%s k=%s minFrac=%s [%s, stab=%s bal=%s]\n",
              coh,
              ifelse(is.null(x$cophenetic), "NA", sprintf("%.2f", x$cophenetic %||% NA_real_)),
              ifelse(is.null(x$silhouette), "NA", sprintf("%.2f", x$silhouette %||% NA_real_)),
              ifelse(is.null(x$k), "NA", as.character(x$k)),
              ifelse(is.null(x$min_cluster_fraction), "NA",
                     sprintf("%.2f", x$min_cluster_fraction %||% NA_real_)),
              ifelse(x$pass, "PASS", "FAIL"),
              ifelse(x$stability_pass %||% FALSE, "y", "n"),
              ifelse(x$balance_pass %||% FALSE, "y", "n")))
}

# ---------------------------------------------------------------------------
# G2 — Bimodality
# ---------------------------------------------------------------------------
cat("\nG2: Bimodality\n")
diag <- read_csv_safe(file.path(OUT_DIR, "continuous_fibrosis_index_diagnostics.csv"))
bic_diff <- diag$value[diag$metric == "BIC_diff_2Gauss_pref"]
sep_sigma <- diag$value[diag$metric == "GMM_separation_sigma"]
dip_p <- diag$value[diag$metric == "DipMC_pvalue"]

g2_bic_pass <- !is.na(bic_diff) && bic_diff >= 10
g2_dip_pass <- !is.na(dip_p) && dip_p < 0.05 && sep_sigma >= 1.0
g2_overall <- g2_bic_pass || g2_dip_pass
cat(sprintf("  BIC diff: %.1f (gate: ≥10)  [%s]\n", bic_diff,
            ifelse(g2_bic_pass, "pass", "fail")))
cat(sprintf("  Dip MC p: %.3f, sep: %.2fσ (gate: p<0.05 AND sep≥1)  [%s]\n",
            dip_p, sep_sigma, ifelse(g2_dip_pass, "pass", "fail")))
cat(sprintf("  Overall G2 (OR): [%s]\n", ifelse(g2_overall, "PASS", "FAIL")))

# ---------------------------------------------------------------------------
# G3 — LOCO-CV cluster prediction AUROC
# ---------------------------------------------------------------------------
cat("\nG3: LOCO-CV cluster prediction AUROC\n")

dge <- readRDS(file.path(INT_DIR, "merged_dge.rds"))
meta <- readRDS(file.path(INT_DIR, "meta_matched.rds"))
meta <- meta[match(colnames(dge), meta$sample_id), ]
pooled_labels <- read_csv_safe(file.path(OUT_DIR, "f3_substate_pooled_labels.csv"))

f3_idx <- which(meta$fibrosis_stage == 3)
samples_f3 <- meta$sample_id[f3_idx]
cohorts_f3 <- meta$dataset[f3_idx]
labels_f3 <- pooled_labels$cluster_pooled[match(samples_f3, pooled_labels$sample_id)]
stopifnot(all(!is.na(labels_f3)))

# Use top-2000 variable batch-corrected genes
dge_f3 <- dge[, f3_idx]
keep <- filterByExpr(dge_f3, group = factor(rep("F3", ncol(dge_f3))), min.count = 5)
dge_f3 <- dge_f3[keep, , keep.lib.sizes = FALSE]
dge_f3 <- calcNormFactors(dge_f3)
v_f3 <- voom(dge_f3, design = NULL)
expr_f3 <- removeBatchEffect(v_f3$E, batch = factor(cohorts_f3))
vars <- rowVars(expr_f3)
top <- order(vars, decreasing = TRUE)[seq_len(min(2000, nrow(expr_f3)))]
expr_top <- expr_f3[top, ]

# Two-class label (largest 2 clusters)
cl_counts <- sort(table(labels_f3), decreasing = TRUE)
if (length(cl_counts) < 2) {
  cat("  Only one F3 cluster found — G3 N/A\n")
  g3_overall <- FALSE
  g3_auroc <- NA
} else {
  cl1 <- names(cl_counts)[1]; cl2 <- names(cl_counts)[2]
  binary_mask <- labels_f3 %in% c(cl1, cl2)
  expr_use <- expr_top[, binary_mask]
  y <- as.integer(labels_f3[binary_mask] == cl1)
  cohort_use <- cohorts_f3[binary_mask]

  # LOCO-CV: hold out each cohort, train PCA + LDA on rest, predict held-out
  unique_cohorts <- unique(cohort_use)
  loco_aucs <- c()
  loco_per <- list()
  for (held in unique_cohorts) {
    test_idx <- which(cohort_use == held)
    train_idx <- which(cohort_use != held)
    if (length(test_idx) < 4 || length(unique(y[test_idx])) < 2) {
      loco_per[[held]] <- list(auc = NA, n_test = length(test_idx))
      next
    }
    # Center on training
    X_train <- t(expr_use[, train_idx])
    X_test  <- t(expr_use[, test_idx])
    means <- colMeans(X_train)
    X_train <- sweep(X_train, 2, means)
    X_test  <- sweep(X_test, 2, means)
    # PCA on train
    pc <- prcomp(X_train, scale. = FALSE)
    n_comp <- min(20, ncol(pc$x))
    train_scores <- pc$x[, 1:n_comp]
    test_scores  <- X_test %*% pc$rotation[, 1:n_comp]
    # Logistic regression on PCA scores
    lr <- suppressWarnings(glm(y[train_idx] ~ ., data = data.frame(train_scores),
                              family = binomial()))
    pred <- predict(lr, newdata = data.frame(test_scores), type = "response")
    if (length(unique(y[test_idx])) >= 2) {
      auc <- suppressMessages(pROC::auc(y[test_idx], pred))
      loco_aucs <- c(loco_aucs, as.numeric(auc))
      loco_per[[held]] <- list(auc = as.numeric(auc), n_test = length(test_idx))
    }
  }
  g3_auroc <- mean(loco_aucs, na.rm = TRUE)
  g3_overall <- !is.na(g3_auroc) && g3_auroc > 0.65
  cat(sprintf("  LOCO-CV mean AUROC: %.3f across %d cohorts  [%s]\n",
              g3_auroc, length(loco_aucs), ifelse(g3_overall, "PASS", "FAIL")))
  for (h in names(loco_per)) {
    cat(sprintf("    %-15s n_test=%d AUC=%.3f\n", h,
                loco_per[[h]]$n_test, loco_per[[h]]$auc %||% NA))
  }
}

# ---------------------------------------------------------------------------
# G4 — Cohort-confounding NMI
# ---------------------------------------------------------------------------
cat("\nG4: Cohort-confounding NMI\n")
mi <- read_csv_safe(file.path(OUT_DIR, "f3_substate_cohort_mi.csv"))
ratio <- mi$cohort_confound_ratio[1]
g4_overall <- !is.na(ratio) && ratio < 0.5
cat(sprintf("  cluster-cohort MI / saturated: %.3f (gate: <0.5)  [%s]\n",
            ratio, ifelse(g4_overall, "PASS", "FAIL")))

# ---------------------------------------------------------------------------
# G5 — Effect-size floor
# ---------------------------------------------------------------------------
cat("\nG5: Effect-size floor\n")
eff <- read_csv_safe(file.path(OUT_DIR, "transition_effect_sizes.csv"))
f1f2 <- eff$median_abs_logFC_top1000[eff$contrast == "F2_vs_F1"]
sub_max <- max(eff$median_abs_logFC_top1000[grepl("F3\\.pool\\.", eff$contrast)],
               na.rm = TRUE)
ratio_eff <- sub_max / f1f2
g5_overall <- !is.na(ratio_eff) && ratio_eff >= 0.30
cat(sprintf("  F1→F2 effect: %.3f, F3 sub-state max: %.3f, ratio: %.2f (gate: ≥0.30)  [%s]\n",
            f1f2, sub_max, ratio_eff, ifelse(g5_overall, "PASS", "FAIL")))

# ---------------------------------------------------------------------------
# Aggregate decision
# ---------------------------------------------------------------------------
gates <- list(
  G1_per_cohort_recovery = list(pass = g1_overall, n_pass_of_4 = n_pass_g1),
  G2_bimodality = list(pass = g2_overall, bic_diff = bic_diff,
                        dip_p = dip_p, sep_sigma = sep_sigma),
  G3_loco_cv_auroc = list(pass = g3_overall, mean_auroc = g3_auroc),
  G4_cohort_confound = list(pass = g4_overall, mi_ratio = ratio),
  G5_effect_size_floor = list(pass = g5_overall,
                              f1f2_effect = f1f2, sub_state_max = sub_max,
                              ratio = ratio_eff)
)

gates_df <- data.frame(
  gate = names(gates),
  passes = sapply(gates, function(x) x$pass),
  diagnostic = sapply(gates, function(x) {
    paste(setdiff(names(x), "pass"), "=",
          round(unlist(x[setdiff(names(x), "pass")]), 3),
          collapse = ", ")
  }),
  stringsAsFactors = FALSE
)

all_pass <- all(gates_df$passes, na.rm = TRUE)
gates_df <- rbind(gates_df,
                  data.frame(gate = "OVERALL", passes = all_pass,
                             diagnostic = sprintf("Phase 2: %s",
                                                  ifelse(all_pass, "GO", "NO-GO")),
                             stringsAsFactors = FALSE))

write.csv(gates_df, file.path(OUT_DIR, "phase1_gates.csv"), row.names = FALSE)

cat("\n== Phase 1.7 gate decision ==\n")
print(gates_df)
cat(sprintf("\n>>> Phase 2 decision: %s <<<\n", ifelse(all_pass, "GO", "NO-GO")))
if (!all_pass) {
  failing <- gates_df$gate[!gates_df$passes & gates_df$gate != "OVERALL"]
  cat("Failing gate(s):", paste(failing, collapse = ", "), "\n")
  cat("Action: fall back to supporting-evidence supplementary in Approach A scope.\n")
  cat("Document failure mode in docs/manuscript/NUMBERS.md.\n")
}
cat("==\n")
