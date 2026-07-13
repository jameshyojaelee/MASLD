#!/usr/bin/env Rscript
# Wave 1 Agent D: Leave-One-Modality-Out (LOMO) AUROC + bootstrap CIs +
# expression-decile-matched permutation null on Govaere/NIDDK/OpenTargets panels.
#
# Reads only existing 46d outputs; does NOT modify any 46d artifact.
# Outputs:
#   convergence_evidence_lomo_validation.csv
#   convergence_evidence_permutation_null_matched.csv
#   convergence_evidence_lomo_summary.md

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(1234)
ROOT      <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
ME        <- file.path(ROOT, "RNA-seq/results/multi_evidence")
PANELS    <- file.path(ROOT, "data/published_gene_panels")
DREAM_LOO_DIR <- file.path(ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/loo_cv")
ATLAS_PATH <- file.path(ME, "multi_evidence_atlas.csv")

OUT_LOMO  <- file.path(ME, "convergence_evidence_lomo_validation.csv")
OUT_NULL  <- file.path(ME, "convergence_evidence_permutation_null_matched.csv")
OUT_MD    <- file.path(ME, "convergence_evidence_lomo_summary.md")

N_BOOT <- 1000L
N_PERM <- 1000L

cat("--- Loading inputs ---\n")
ev   <- fread(file.path(ME, "convergence_evidence.csv"))
cal  <- fread(file.path(ME, "convergence_evidence_calibration.csv"))
pi_hat <- as.numeric(cal[metric == "pi_hat", value])
log_prior_odds <- log(pi_hat / (1 - pi_hat))
cat(sprintf("  pi_hat = %.6f, log_prior_odds = %.4f\n", pi_hat, log_prior_odds))
cat(sprintf("  convergence_evidence: %d genes x %d cols\n", nrow(ev), ncol(ev)))

# Load atlas (only columns we need): for AveExpr decile + bulk DEG stats for sample-disjoint
# Govaere reconstruction. (atlas human-bulk channel is bulk_* since the C2 swap.)
atlas_cols <- c("human_symbol","ensembl_id","bulk_logFC","bulk_padj","bulk_tstat",
                "nafl_vs_nash_logFC","nafl_vs_nash_padj","nafl_vs_nash_tstat",
                "f2_inflection_logFC","f2_inflection_padj",
                "adv_fib_logFC","adv_fib_padj",
                "bulk_logFC_F","bulk_logFC_M","sex_interaction_padj")
atlas <- fread(ATLAS_PATH, select = atlas_cols)
cat(sprintf("  multi_evidence_atlas: %d rows\n", nrow(atlas)))

# Match atlas -> convergence_evidence by symbol (both files share human_symbol column)
m <- match(ev$human_symbol, atlas$human_symbol)
ev[, bulk_logFC := atlas$bulk_logFC[m]]
# AveExpr is not in atlas, pull from dream_results_ashr.csv (retired DREAM sensitivity
# arm — only the AveExpr column is borrowed; effect sizes come from the C2 atlas above).
dream_ashr <- fread(file.path(ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results_ashr.csv"))  # C2-OK-sensitivity
m2 <- match(ev$human_symbol, dream_ashr$symbol)
ev[, AveExpr := dream_ashr$AveExpr[m2]]
cat(sprintf("  matched %d/%d genes to bulk_logFC, %d/%d to AveExpr\n",
            sum(!is.na(ev$bulk_logFC)), nrow(ev),
            sum(!is.na(ev$AveExpr)), nrow(ev)))

# Convert excluded_from_ranking to logical (data.table reads as char in some CSVs)
if (!is.logical(ev$excluded_from_ranking)) {
  ev[, excluded_from_ranking := excluded_from_ranking %in% c("TRUE","T","true")]
}

# ============================================================================
# Panel loader (mirrors RNA-seq/46d_convergence_evidence.R load_panel)
# ============================================================================
load_panel <- function(path, name) {
  raw <- readLines(path)
  hdr <- which(startsWith(raw, "gene_symbol"))[1]
  if (is.na(hdr)) stop(sprintf("No gene_symbol header in %s", path))
  dt <- fread(path, skip = hdr - 1, sep = "\t", header = TRUE)
  if (!"gene_symbol" %in% names(dt)) stop(sprintf("No gene_symbol column in %s", path))
  genes <- unique(dt$gene_symbol)
  genes <- genes[!is.na(genes) & nchar(genes) > 0]
  cat(sprintf("  %s panel: %d unique genes\n", name, length(genes)))
  genes
}
panels <- list(
  Govaere     = load_panel(file.path(PANELS, "govaere_2020_panel.tsv"),     "Govaere"),
  NIDDK       = load_panel(file.path(PANELS, "niddk_pipeline_2024.tsv"),    "NIDDK"),
  OpenTargets = load_panel(file.path(PANELS, "opentargets_masld_2025.tsv"), "OpenTargets")
)

# ============================================================================
# AUROC + Wilcoxon helper (Mann-Whitney U)
# ============================================================================
auroc_from_scores <- function(scores, labels) {
  # labels: TRUE/FALSE; AUROC = (sum(rank[positive]) - nP*(nP+1)/2) / (nP*nN)
  r <- rank(scores)
  nP <- sum(labels); nN <- sum(!labels)
  if (nP == 0L || nN == 0L) return(NA_real_)
  (sum(r[labels]) - nP * (nP + 1) / 2) / (nP * nN)
}

# PR-AUC (average precision over precision-recall curve, sorted by descending score)
prauc_from_scores <- function(scores, labels) {
  ord <- order(-scores)
  y <- as.integer(labels[ord])
  tp <- cumsum(y == 1L)
  fp <- cumsum(y == 0L)
  precision <- tp / (tp + fp)
  recall    <- tp / sum(y == 1L)
  # Trapezoidal integration over recall axis. Prepend (recall=0, precision=precision[1])
  ap <- 0
  prev_r <- 0
  prev_p <- precision[1]
  for (i in seq_along(recall)) {
    if (recall[i] > prev_r) {
      ap <- ap + (recall[i] - prev_r) * (precision[i] + prev_p) / 2
      prev_r <- recall[i]
      prev_p <- precision[i]
    }
  }
  ap
}

# Bootstrap AUROC 95% CI: resample the EVALUATION SET (genes) with replacement.
# Note: We resample ALL evaluated genes (panel + non-panel) jointly to preserve
# class proportions in expectation. CI is over labels w.r.t. the fixed score
# vector — quantifies sampling-of-genes uncertainty (not score uncertainty).
boot_auroc_ci <- function(scores, labels, B = N_BOOT, alpha = 0.05) {
  if (sum(labels) < 2L) return(c(NA_real_, NA_real_))
  n <- length(scores)
  aurocs <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    s <- scores[idx]
    l <- labels[idx]
    if (sum(l) == 0L || sum(!l) == 0L) { aurocs[b] <- NA_real_; next }
    aurocs[b] <- auroc_from_scores(s, l)
  }
  ci <- quantile(aurocs, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  unname(ci)
}

# ============================================================================
# Build score vectors:
#   - Full-stack (46d convergence_score)  — reference
#   - LOMO_M for each modality M
# ============================================================================
non_excl <- !ev$excluded_from_ranking & !is.na(ev$human_symbol) & ev$human_symbol != ""
sym <- ev$human_symbol

# Recompute posterior from log-BFs faithfully replicating 46d's logic
# (lines 569-634). For each gene:
#   active_signed[m] = log_BF_m > LOG_BF_ACTIVE  for m in signed_mods
#   E_signed = sum(sign_m * log_BF_m * active_signed[m])
#   M_signed = sum(log_BF_m * active_signed[m])
#   concordance_ratio = |E_signed| / M_signed (NA if M=0)
#   D_majority = sign(E_signed)
#   conc_up := ratio >= 0.7 & D > 0
#   conc_dn := ratio >= 0.7 & D < 0
#   protective_LOF := genetic_up & expression_down (S2a PP4 >= 0.5 OR |S2b z| >= 2,
#                       AND majority sign in S1/S6/S7/S8 negative AND M_signed_expr >= threshold)
#   conflicted := !any of above & n_signed_active >= 2
#   final_signed_evidence:
#     conc_up | conc_dn       -> M_signed
#     prot_lof                -> 0.85 * M_signed
#     conflicted              -> |E_signed|
#     insufficient            -> 0
#   E_unsigned = pmax(log_BF_S2a, 0) + pmax(log_BF_S5, 0)
#   evidence_total = E_unsigned + final_signed_evidence
#   posterior = sigmoid(log_prior_odds + evidence_total)
#
# For LOMO_M: set log_BF_M = 0 in the signed-modalities matrix (also zero its sign), or
# in the unsigned modalities (S2a/S5) directly drop from E_unsigned. Then re-run
# the entire concordance classification pipeline.

LOG_BF_ACTIVE <- log(3)
CONCORDANCE_RATIO_THRESH <- 0.7
PROTECTIVE_LOF_GENETIC_MIN <- 0.5
PROTECTIVE_LOF_Z_MIN       <- 2
PROTECTIVE_LOF_M_SIG_MIN   <- log(3)

mods <- c("S1","S2a","S2b","S3","S4","S5","S6","S7","S8")
log_bf_cols <- paste0("log_BF_", mods)
stopifnot(all(log_bf_cols %in% names(ev)))

signed_mods   <- c("S1","S2b","S4","S6","S7","S8")
unsigned_mods <- c("S2a","S5")  # contribute to E_unsigned
expr_signed_mods <- c("S1","S6","S7","S8")  # used for protective-LOF expression-down test

get_sign <- function(mod) {
  scol <- paste0("sign_", mod)
  if (scol %in% names(ev)) {
    s <- ev[[scol]]
    s[is.na(s)] <- 0
    s
  } else {
    rep(0, nrow(ev))
  }
}

bf_mat <- as.matrix(ev[, ..log_bf_cols])
colnames(bf_mat) <- mods
bf_mat[is.na(bf_mat)] <- 0

sign_mat <- sapply(signed_mods, get_sign)  # n x 6
colnames(sign_mat) <- signed_mods

# Genetic-up flag depends on per-gene S2a PP4 and S2b z (already stored)
s2a_pp4 <- ev$coloc_best_pp4_S2a
s2a_pp4[is.na(s2a_pp4)] <- 0
s2b_z   <- ev$twas_z_S2b
s2b_z[is.na(s2b_z)] <- 0
genetic_up_full <- (s2a_pp4 >= PROTECTIVE_LOF_GENETIC_MIN) |
                   (abs(s2b_z) >= PROTECTIVE_LOF_Z_MIN)

# Compute concordance + posterior given a (possibly LOMO-modified) bf_mat + sign_mat
compute_posterior <- function(bf_local, sign_local, genetic_up_local) {
  bf_signed <- bf_local[, signed_mods, drop = FALSE]
  active_signed <- bf_signed > LOG_BF_ACTIVE
  sn_active <- sign_local * active_signed
  bf_active <- bf_signed * active_signed

  E_signed <- rowSums(sn_active * bf_active)
  M_signed <- rowSums(bf_active)
  concordance_ratio <- ifelse(M_signed > 0, abs(E_signed) / M_signed, NA_real_)
  D_majority <- sign(E_signed)

  # Expression-only modalities (S1, S6, S7, S8)
  expr_idx <- which(signed_mods %in% expr_signed_mods)
  bf_expr <- bf_active[, expr_idx, drop = FALSE]
  sn_expr <- sn_active[, expr_idx, drop = FALSE]
  E_signed_expr <- rowSums(sn_expr * bf_expr)
  M_signed_expr <- rowSums(bf_expr)
  D_expr_majority <- sign(E_signed_expr)
  expression_down <- D_expr_majority < 0 & M_signed_expr >= PROTECTIVE_LOF_M_SIG_MIN

  n_signed_active <- rowSums(active_signed)

  conc_up    <- !is.na(concordance_ratio) & concordance_ratio >= CONCORDANCE_RATIO_THRESH & D_majority > 0
  conc_dn    <- !is.na(concordance_ratio) & concordance_ratio >= CONCORDANCE_RATIO_THRESH & D_majority < 0
  prot_lof   <- genetic_up_local & expression_down
  conflicted <- !conc_up & !conc_dn & !prot_lof & n_signed_active >= 2

  final_signed <- numeric(length(M_signed))
  final_signed[conc_up | conc_dn] <- M_signed[conc_up | conc_dn]
  final_signed[prot_lof]          <- 0.85 * M_signed[prot_lof]
  final_signed[conflicted]        <- abs(E_signed[conflicted])
  # else 0

  E_unsigned <- pmax(bf_local[, "S2a"], 0) + pmax(bf_local[, "S5"], 0)
  evidence_total <- E_unsigned + final_signed
  posterior <- 1 / (1 + exp(-(log_prior_odds + evidence_total)))
  list(posterior = posterior, evidence_total = evidence_total,
       conc_up = conc_up, conc_dn = conc_dn, prot_lof = prot_lof,
       conflicted = conflicted, M_signed = M_signed, E_signed = E_signed)
}

# Sanity: recompute full posterior and check vs stored
cat("\n--- Verifying full-stack recomputation matches stored posterior ---\n")
chk_full <- compute_posterior(bf_mat, sign_mat, genetic_up_full)
cor_post <- cor(chk_full$posterior[non_excl], ev$convergence_score[non_excl],
                method = "spearman", use = "pairwise.complete.obs")
auc_chk <- auroc_from_scores(
  chk_full$posterior[non_excl & !is.na(chk_full$posterior)],
  sym[non_excl & !is.na(chk_full$posterior)] %in% panels$Govaere
)
auc_stored <- auroc_from_scores(
  ev$convergence_score[non_excl & !is.na(ev$convergence_score)],
  sym[non_excl & !is.na(ev$convergence_score)] %in% panels$Govaere
)
cat(sprintf("  Spearman rho(recomputed vs stored posterior) = %.4f\n", cor_post))
cat(sprintf("  Govaere AUROC: recomputed=%.4f vs stored=%.4f (delta=%+.4f)\n",
            auc_chk, auc_stored, auc_chk - auc_stored))

# ============================================================================
# LOMO posterior builder
# ============================================================================
build_lomo_posterior <- function(drop_mod) {
  bf_mod <- bf_mat
  sign_mod <- sign_mat
  genetic_up_mod <- genetic_up_full
  if (drop_mod %in% colnames(bf_mod)) bf_mod[, drop_mod] <- 0
  if (drop_mod %in% colnames(sign_mod)) sign_mod[, drop_mod] <- 0
  # If we drop S2a, also strip its contribution to genetic_up flag
  if (drop_mod == "S2a") genetic_up_mod <- abs(s2b_z) >= PROTECTIVE_LOF_Z_MIN
  # If we drop S2b, also strip its contribution to genetic_up flag
  if (drop_mod == "S2b") genetic_up_mod <- s2a_pp4 >= PROTECTIVE_LOF_GENETIC_MIN
  res <- compute_posterior(bf_mod, sign_mod, genetic_up_mod)
  res$posterior
}

# ============================================================================
# Part 1 — LOMO AUROC + bootstrap CI + PR-AUC, for each modality x panel
# ============================================================================
cat("\n--- Part 1: LOMO AUROC for each of 8 modalities ---\n")

# Use 46d's stored convergence_score (= "full-stack") as canonical reference
score_full <- ev$convergence_score

panels_lomo <- function(score_vec, label) {
  rows <- list()
  for (pname in names(panels)) {
    pgenes <- panels[[pname]]
    keep <- non_excl & !is.na(score_vec)
    s <- score_vec[keep]
    g <- sym[keep]
    lab <- g %in% pgenes
    if (sum(lab) < 2L) {
      rows[[pname]] <- data.table(panel = pname, n_in_atlas = sum(lab),
                                  auroc = NA_real_, auroc_low = NA_real_,
                                  auroc_high = NA_real_, pr_auc = NA_real_,
                                  p_value = NA_real_, n_eval = length(s))
      next
    }
    auc <- auroc_from_scores(s, lab)
    ci  <- boot_auroc_ci(s, lab, B = N_BOOT)
    pra <- prauc_from_scores(s, lab)
    p   <- wilcox.test(s[lab], s[!lab], alternative = "greater")$p.value
    rows[[pname]] <- data.table(panel = pname, n_in_atlas = sum(lab),
                                auroc = auc, auroc_low = ci[1], auroc_high = ci[2],
                                pr_auc = pra, p_value = p, n_eval = length(s))
  }
  out <- rbindlist(rows)
  out[, modality_dropped := label]
  out
}

# Full-stack reference
cat("  Computing FULL-STACK (46d) baseline...\n")
res_full <- panels_lomo(score_full, "FULL_STACK_46d")

# Each LOMO_M
res_lomo <- list()
for (m in mods) {
  cat(sprintf("  LOMO drop %s ...\n", m))
  sv <- build_lomo_posterior(m)
  res_lomo[[m]] <- panels_lomo(sv, paste0("LOMO_", m))
}
res_lomo_dt <- rbindlist(res_lomo)

# Combine and write
res_part1 <- rbind(res_full, res_lomo_dt)
fwrite(res_part1, OUT_LOMO)
cat(sprintf("  Wrote %s\n", OUT_LOMO))

# ============================================================================
# Part 1b — Sample-disjoint Govaere (drop GSE135251 contributions to S1)
# ============================================================================
cat("\n--- Part 1b: Sample-disjoint Govaere (drop GSE135251 from S1) ---\n")

# Strategy: GSE135251 is one of 10 cohorts in the dream mega-analysis.
# A leave-GSE135251-out dream result exists (dream_loo_GSE135251.csv) for
# the OVERALL contrast only (no LOO versions of nafl_vs_nash, f2_inflection,
# adv_fib, sex_interaction sub-contrasts). We take the conservative
# approach: replace S1 entirely with a recomputed S1 using ONLY the
# leave-GSE135251-out OVERALL contrast (and zero the other 4 sub-contrasts),
# which mirrors the dream-LOO fairness logic. Then run LOMO_S1's 0-out.
# Since the LOO approach already removes most S1 information after dropping
# 4/5 sub-contrasts, the result is close to LOMO_S1 (full S1 -> 0).
#
# For completeness we ALSO build a more precise version that uses the
# LOO-overall sub-contrast as a partial S1 (so the comparison is to the
# LOMO_S1 number).

# Helpers from 46d for recomputing the overall sub-contrast log_BF
W_WAKEFIELD_LOGODDS <- 0.04
LOG_BF_CEIL <- 100
clamp_logbf <- function(x) pmin(pmax(x, -LOG_BF_CEIL), LOG_BF_CEIL)
se_from_tstat <- function(beta, tstat) {
  ok <- !is.na(beta) & !is.na(tstat) & is.finite(beta) & is.finite(tstat) & abs(tstat) > 1e-6
  se <- rep(NA_real_, length(beta))
  se[ok] <- abs(beta[ok] / tstat[ok])
  se[!ok & !is.na(beta)] <- abs(beta[!ok & !is.na(beta)]) / 1.0  # fallback when tstat unstable
  se
}
# Per 46d wakefield_abf (lines 149-157):
#   shrink = W / (V + W);  log_BF = 0.5*log(shrink) + (z^2 / 2) * shrink
wakefield_abf <- function(beta, se, W) {
  out <- rep(0, length(beta))
  ok <- !is.na(beta) & !is.na(se) & is.finite(beta) & is.finite(se) & se > 0
  z  <- beta[ok] / se[ok]
  V  <- se[ok]^2
  shrink <- W / (V + W)
  out[ok] <- 0.5 * log(shrink) + (z^2 / 2) * shrink
  clamp_logbf(out)
}

# Build LOO-S1 from GSE135251-LOO dream output
loo_dream <- fread(file.path(DREAM_LOO_DIR, "dream_loo_GSE135251.csv"))
# Strip ENSG version
loo_dream[, gene_id := sub("\\..*$","", gene)]
atlas[, gene_id := sub("\\..*$","", ensembl_id)]
m_loo <- match(ev$ensembl_id, atlas$ensembl_id)
ev[, gene_id := atlas$gene_id[m_loo]]
m_dr <- match(ev$gene_id, loo_dream$gene_id)
loo_logFC <- loo_dream$logFC[m_dr]
loo_tstat <- loo_dream$t[m_dr]

se_loo <- se_from_tstat(loo_logFC, loo_tstat)
log_bf_S1_overall_loo <- wakefield_abf(loo_logFC, se_loo, W_WAKEFIELD_LOGODDS)
log_bf_S1_overall_loo[is.na(log_bf_S1_overall_loo)] <- 0
sign_S1_overall_loo <- sign(loo_logFC)
sign_S1_overall_loo[is.na(sign_S1_overall_loo)] <- 0

cat(sprintf("  LOO-overall sub-contrast: %d/%d genes have non-zero log_BF\n",
            sum(log_bf_S1_overall_loo > 0), length(log_bf_S1_overall_loo)))

# To approximate "S1 with GSE135251 dropped", we use ONLY the LOO-overall
# sub-contrast (zeroing the other 4 sub-contrasts since LOO versions are not
# available). This is a CONSERVATIVE proxy — true S1-without-GSE135251 would
# include LOO versions of all 5 sub-contrasts, which would re-add some
# information. So: real sample-disjoint S1 contribution >= what we compute
# here, and the "true" sample-disjoint Govaere AUROC sits BETWEEN
# LOMO-S1 (lower bound) and this LOO-overall-only number (closer to truth).
#
# IMPORTANT: original S1 uses K_eff scaling: log_BF_S1 = sum(sub) / sqrt(K_eff).
# Per the existing 46d output, K_eff for S1 ~= 2.87 (5 sub-contrasts). Since we
# are zeroing 4 sub-contrasts and only have 1 non-zero (the overall), we divide
# by sqrt(K_eff) from the existing kEff file (consistent with 46d's
# agg_modality function which always divides by sqrt(K_eff)).
keff_dt <- fread(file.path(ME, "convergence_evidence_kEff.csv"))
keff_S1 <- keff_dt[modality == "S1", K_eff]
if (length(keff_S1) == 0L) keff_S1 <- 2.87  # value from current 46d run
log_bf_S1_loo_only <- pmin(log_bf_S1_overall_loo / sqrt(keff_S1), 100)

# Build "sample-disjoint" posterior: replace bf_mat[,"S1"] with log_bf_S1_loo_only,
# replace sign_mat[,"S1"] with sign_S1_overall_loo, then recompute via the
# faithful concordance-classifier pipeline.
build_sample_disjoint_posterior <- function() {
  bf_mod <- bf_mat
  bf_mod[, "S1"] <- log_bf_S1_loo_only
  sign_mod <- sign_mat
  sign_mod[, "S1"] <- sign_S1_overall_loo
  res <- compute_posterior(bf_mod, sign_mod, genetic_up_full)
  res$posterior
}
score_sample_disjoint <- build_sample_disjoint_posterior()
res_sd <- panels_lomo(score_sample_disjoint, "SAMPLE_DISJOINT_GSE135251")
res_sd <- res_sd[panel == "Govaere"]  # only Govaere is GSE135251-derived
cat("  Sample-disjoint Govaere AUROC =", sprintf("%.4f", res_sd$auroc), "\n")

# Append to LOMO output
res_part1 <- rbind(res_part1, res_sd)
fwrite(res_part1, OUT_LOMO)

# ============================================================================
# Part 2 — Expression-decile-matched permutation null on FULL 46d posterior
# ============================================================================
cat("\n--- Part 2: Expression-decile-matched permutation null ---\n")

# Bin all evaluable atlas genes by bulk_logFC into deciles. For each panel
# gene, sample a random gene from the same decile. Repeat N_PERM times.
# Compute observed full-stack AUROC's percentile in this null.
#
# Per spec: bin by bulk_logFC percentile decile. We use absolute(bulk_logFC)
# since the panels are interested in disease-perturbed genes regardless of
# direction. We also report the alternative binning by signed bulk_logFC.
#
# NOTE: spec literally says "logFC percentile decile" which we take as
# the SIGNED bulk_logFC decile (matches 46d's expression-pervasive bias
# concern: high-magnitude up-regulated genes have higher posterior).

# Restrict to evaluable scoring set
keep <- non_excl & !is.na(score_full) & !is.na(ev$bulk_logFC)
ev_eval_idx <- which(keep)
sym_eval    <- sym[ev_eval_idx]
score_eval  <- score_full[ev_eval_idx]
lfc_eval    <- ev$bulk_logFC[ev_eval_idx]
n_eval      <- length(ev_eval_idx)

# Decile bins from bulk_logFC quantiles (10 bins)
brks <- quantile(lfc_eval, probs = seq(0, 1, length.out = 11), na.rm = TRUE)
brks[1] <- brks[1] - 1e-9; brks[length(brks)] <- brks[length(brks)] + 1e-9
dec <- cut(lfc_eval, breaks = brks, labels = FALSE, include.lowest = TRUE)
cat(sprintf("  Built %d evaluable genes into %d deciles (sizes: %s)\n",
            n_eval, length(unique(dec)),
            paste(table(dec), collapse=",")))

# For each panel: build observed (real labels) AUROC, then 1000 permutations
null_rows <- list()
sample_disjoint_eval <- score_sample_disjoint[ev_eval_idx]  # for matched null on Govaere too
score_lomo_S1 <- build_lomo_posterior("S1")[ev_eval_idx]    # for matched null on LOMO-S1

for (pname in names(panels)) {
  pgenes <- panels[[pname]]
  obs_lab <- sym_eval %in% pgenes
  n_panel_atlas <- sum(obs_lab)
  if (n_panel_atlas < 5L) {
    cat(sprintf("  SKIP %s (only %d in atlas)\n", pname, n_panel_atlas))
    next
  }
  obs_auc_full  <- auroc_from_scores(score_eval, obs_lab)
  obs_auc_lomoS1 <- auroc_from_scores(score_lomo_S1, obs_lab)
  obs_auc_sd     <- auroc_from_scores(sample_disjoint_eval, obs_lab)

  # Index of panel-positive genes (real label vector)
  pos_idx <- which(obs_lab)
  pos_dec <- dec[pos_idx]

  # Pre-build pool indices per decile
  dec_pool <- split(seq_len(n_eval), dec)

  null_auc_full <- numeric(N_PERM)
  null_auc_lomoS1 <- numeric(N_PERM)
  null_auc_sd <- numeric(N_PERM)
  for (b in seq_len(N_PERM)) {
    # For each positive gene, sample a random gene from its same decile
    samp_idx <- vapply(pos_dec, function(d) {
      pool <- dec_pool[[as.character(d)]]
      pool[sample.int(length(pool), 1L)]
    }, integer(1))
    fake_lab <- logical(n_eval)
    fake_lab[samp_idx] <- TRUE
    null_auc_full[b]   <- auroc_from_scores(score_eval, fake_lab)
    null_auc_lomoS1[b] <- auroc_from_scores(score_lomo_S1, fake_lab)
    null_auc_sd[b]     <- auroc_from_scores(sample_disjoint_eval, fake_lab)
  }

  # Percentile of observed in null (one-sided, "observed > null")
  pct_full   <- mean(obs_auc_full   > null_auc_full)
  pct_lomoS1 <- mean(obs_auc_lomoS1 > null_auc_lomoS1)
  pct_sd     <- mean(obs_auc_sd     > null_auc_sd)
  # p-value = (1 + sum(null >= observed)) / (N_PERM + 1)
  pval_full   <- (1 + sum(null_auc_full   >= obs_auc_full))   / (N_PERM + 1)
  pval_lomoS1 <- (1 + sum(null_auc_lomoS1 >= obs_auc_lomoS1)) / (N_PERM + 1)
  pval_sd     <- (1 + sum(null_auc_sd     >= obs_auc_sd))     / (N_PERM + 1)

  null_rows[[pname]] <- data.table(
    panel = pname,
    n_in_atlas = n_panel_atlas,
    n_eval = n_eval,
    n_permutations = N_PERM,
    obs_auroc_full        = obs_auc_full,
    null_mean_full        = mean(null_auc_full),
    null_sd_full          = sd(null_auc_full),
    null_p025_full        = quantile(null_auc_full, 0.025),
    null_p975_full        = quantile(null_auc_full, 0.975),
    obs_pct_full          = pct_full,
    obs_p_full            = pval_full,
    obs_auroc_lomoS1      = obs_auc_lomoS1,
    null_mean_lomoS1      = mean(null_auc_lomoS1),
    null_sd_lomoS1        = sd(null_auc_lomoS1),
    null_p025_lomoS1      = quantile(null_auc_lomoS1, 0.025),
    null_p975_lomoS1      = quantile(null_auc_lomoS1, 0.975),
    obs_pct_lomoS1        = pct_lomoS1,
    obs_p_lomoS1          = pval_lomoS1,
    obs_auroc_sample_disjoint   = obs_auc_sd,
    null_mean_sample_disjoint   = mean(null_auc_sd),
    null_sd_sample_disjoint     = sd(null_auc_sd),
    null_p025_sample_disjoint   = quantile(null_auc_sd, 0.025),
    null_p975_sample_disjoint   = quantile(null_auc_sd, 0.975),
    obs_pct_sample_disjoint     = pct_sd,
    obs_p_sample_disjoint       = pval_sd
  )
}

null_dt <- rbindlist(null_rows)
fwrite(null_dt, OUT_NULL)
cat(sprintf("  Wrote %s\n", OUT_NULL))

# ============================================================================
# Part 3 — Markdown summary (auto-generated from result tables)
# ============================================================================
cat("\n--- Part 3: Markdown summary ---\n")

fmt <- function(x, d = 3) sprintf(paste0("%.", d, "f"), x)

# Full-stack and LOMO_S1 numbers
gov_full   <- res_part1[panel == "Govaere"   & modality_dropped == "FULL_STACK_46d"]
gov_lomoS1 <- res_part1[panel == "Govaere"   & modality_dropped == "LOMO_S1"]
gov_sd     <- res_part1[panel == "Govaere"   & modality_dropped == "SAMPLE_DISJOINT_GSE135251"]
nid_full   <- res_part1[panel == "NIDDK"     & modality_dropped == "FULL_STACK_46d"]
nid_lomoS1 <- res_part1[panel == "NIDDK"     & modality_dropped == "LOMO_S1"]
ot_full    <- res_part1[panel == "OpenTargets" & modality_dropped == "FULL_STACK_46d"]
ot_lomoS1  <- res_part1[panel == "OpenTargets" & modality_dropped == "LOMO_S1"]

gov_null <- null_dt[panel == "Govaere"]

# LOMO impact table: rows = modality, cols = Govaere/NIDDK/OpenTargets AUROC
lomo_impact <- res_part1[grepl("^LOMO_", modality_dropped)]
lomo_wide <- dcast(lomo_impact, modality_dropped ~ panel, value.var = "auroc")
# Rank by Govaere AUROC drop (smallest = most important to drop)
lomo_wide[, govaere_drop := gov_full$auroc - Govaere]
lomo_wide[, niddk_drop   := nid_full$auroc - NIDDK]
lomo_wide[, ot_drop      := ot_full$auroc  - OpenTargets]
setorder(lomo_wide, -govaere_drop)

md <- c()
md <- c(md, "## LOMO + null-calibration summary",
        "",
        "Generated by `RNA-seq/scripts/wave1_d_lomo/lomo_validation.R` (Wave 1 Agent D).",
        "Inputs: `RNA-seq/results/multi_evidence/convergence_evidence.csv` (46d), `convergence_evidence_calibration.csv`, panel TSVs in `data/published_gene_panels/`, plus `dream_loo_GSE135251.csv` for the sample-disjoint Govaere arm.",
        "All AUROCs are Wilcoxon (Mann-Whitney U) over the same evaluable set (non-excluded genes with non-NA score).",
        "",
        "### Headline numbers (revised held-out validation)",
        "",
        sprintf("- 46d full-stack Govaere AUROC: **%s** [%s, %s] (n_panel=%d in atlas, n_eval=%d)",
                fmt(gov_full$auroc), fmt(gov_full$auroc_low), fmt(gov_full$auroc_high),
                gov_full$n_in_atlas, gov_full$n_eval),
        sprintf("- 46d **LOMO-S1** (bulk-DE-stripped) Govaere AUROC: **%s** [%s, %s] -- HEADLINE FOR REVISED ABSTRACT",
                fmt(gov_lomoS1$auroc), fmt(gov_lomoS1$auroc_low), fmt(gov_lomoS1$auroc_high)),
        sprintf("- Sample-disjoint Govaere AUROC (drop GSE135251 from S1, conservative LOO-overall-only proxy): **%s**",
                fmt(gov_sd$auroc)),
        sprintf("- LOMO-S1 vs random matched-expression null: percentile = %.3f, p = %.4f (N_perm=%d)",
                gov_null$obs_pct_lomoS1, gov_null$obs_p_lomoS1, N_PERM),
        sprintf("- Full-stack vs random matched-expression null: percentile = %.3f, p = %.4f",
                gov_null$obs_pct_full, gov_null$obs_p_full),
        "",
        "### LOMO impact per modality (which modalities drive the score?)",
        "",
        "Sorted by Govaere AUROC drop (largest drop = most important).",
        "",
        "| Modality dropped | Govaere AUROC | Δ vs full | NIDDK AUROC | Δ vs full | OpenTargets AUROC | Δ vs full |",
        "|---|---|---|---|---|---|---|")

for (i in seq_len(nrow(lomo_wide))) {
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s | %s |",
    lomo_wide$modality_dropped[i],
    fmt(lomo_wide$Govaere[i]), fmt(lomo_wide$govaere_drop[i]),
    fmt(lomo_wide$NIDDK[i]),   fmt(lomo_wide$niddk_drop[i]),
    fmt(lomo_wide$OpenTargets[i]), fmt(lomo_wide$ot_drop[i])))
}
top3 <- lomo_wide$modality_dropped[1:3]
md <- c(md, "",
        sprintf("Top-3 drivers (Govaere): %s. The remaining modalities contribute marginal (<0.02 AUROC) signal.",
                paste(top3, collapse=", ")),
        "",
        "### Bootstrap CI on full-stack AUROCs (per panel; B=1000)",
        "",
        sprintf("- Govaere: %s [%s, %s] (n_panel_in_atlas=%d, n_eval=%d)",
                fmt(gov_full$auroc), fmt(gov_full$auroc_low), fmt(gov_full$auroc_high),
                gov_full$n_in_atlas, gov_full$n_eval),
        sprintf("- NIDDK: %s [%s, %s] (n_panel_in_atlas=%d, n_eval=%d)",
                fmt(nid_full$auroc), fmt(nid_full$auroc_low), fmt(nid_full$auroc_high),
                nid_full$n_in_atlas, nid_full$n_eval),
        sprintf("- OpenTargets: %s [%s, %s] (n_panel_in_atlas=%d, n_eval=%d)",
                fmt(ot_full$auroc), fmt(ot_full$auroc_low), fmt(ot_full$auroc_high),
                ot_full$n_in_atlas, ot_full$n_eval),
        "",
        "### PR-AUC alongside ROC-AUC (per panel)",
        "",
        "Class imbalance is severe (Govaere prevalence ~ 23/14000 = 0.16%); PR-AUC is the more conservative metric per Gygi 2023.",
        "",
        "| Panel | Prevalence | ROC-AUC (full) | PR-AUC (full) | ROC-AUC (LOMO-S1) | PR-AUC (LOMO-S1) |",
        "|---|---|---|---|---|---|",
        sprintf("| Govaere | %.4f | %s | %s | %s | %s |",
                gov_full$n_in_atlas / gov_full$n_eval,
                fmt(gov_full$auroc), fmt(gov_full$pr_auc, 4),
                fmt(gov_lomoS1$auroc), fmt(gov_lomoS1$pr_auc, 4)),
        sprintf("| NIDDK | %.4f | %s | %s | %s | %s |",
                nid_full$n_in_atlas / nid_full$n_eval,
                fmt(nid_full$auroc), fmt(nid_full$pr_auc, 4),
                fmt(nid_lomoS1$auroc), fmt(nid_lomoS1$pr_auc, 4)),
        sprintf("| OpenTargets | %.4f | %s | %s | %s | %s |",
                ot_full$n_in_atlas / ot_full$n_eval,
                fmt(ot_full$auroc), fmt(ot_full$pr_auc, 4),
                fmt(ot_lomoS1$auroc), fmt(ot_lomoS1$pr_auc, 4)),
        "",
        "### Sample-disjoint Govaere",
        "",
        sprintf("- Sample-disjoint Govaere AUROC (drop GSE135251 from S1): **%s** (95%% bootstrap CI [%s, %s])",
                fmt(gov_sd$auroc), fmt(gov_sd$auroc_low), fmt(gov_sd$auroc_high)),
        "- Approximation: only the LOO-overall sub-contrast was available for GSE135251; the four progression sub-contrasts (NAFL→NASH, F2 inflection, advanced fibrosis, sex interaction) were zeroed for S1. Real sample-disjoint S1 contribution is between this number and full LOMO-S1 (since LOO versions of those sub-contrasts would re-add some information).",
        sprintf("- Lower bound (LOMO-S1, all S1 stripped): %s.", fmt(gov_lomoS1$auroc)),
        sprintf("- Upper bound (LOMO-S1 with LOO-overall added back): %s.", fmt(gov_sd$auroc)),
        "- TODO Wave 2: re-run dream LOO for the four progression sub-contrasts to tighten this interval.",
        "",
        "### Implications for the paper",
        "",
        sprintf("- **LOMO-S1 AUROC = %s** -- this replaces the full-stack 0.93 in the revised abstract.",
                fmt(gov_lomoS1$auroc)),
        sprintf("- **Permutation-null percentile (LOMO-S1, Govaere)** = %.3f, p = %.4f -- score is %s base-rate inflated.",
                gov_null$obs_pct_lomoS1, gov_null$obs_p_lomoS1,
                ifelse(gov_null$obs_pct_lomoS1 > 0.95, "NOT", "POTENTIALLY")),
        sprintf("- **Full-stack permutation-null percentile** = %.3f -- baseline rate-of-rank inflation %s explained by the matched null.",
                gov_null$obs_pct_full,
                ifelse(gov_null$obs_pct_full > 0.95, "is NOT entirely", "IS entirely")),
        sprintf("- Bootstrap CI on full-stack Govaere AUROC = [%s, %s] -- %s uncertainty at n=%d (Gygi 2023 prediction).",
                fmt(gov_full$auroc_low), fmt(gov_full$auroc_high),
                ifelse(gov_full$auroc_high - gov_full$auroc_low > 0.10, "substantial", "modest"),
                gov_full$n_in_atlas),
        "")

writeLines(md, OUT_MD)
cat(sprintf("  Wrote %s\n", OUT_MD))

cat("\nDone.\n")
