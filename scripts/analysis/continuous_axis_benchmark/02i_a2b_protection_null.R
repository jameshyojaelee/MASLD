#!/usr/bin/env Rscript
# 02i: the two controls without which arm A2b means nothing.
#
# CONTROL 1, permuted-protection null. Protecting fibrosis stage in ComBat's mod
# preserves fibrosis-associated variance by construction, so an A2b axis that
# correlates with fibrosis is partly measuring what was protected. The null
# protects a WITHIN-COHORT SHUFFLE of the same labels. That leaves exactly the
# same amount of variance unremoved, with identical cohort-by-level margins and
# identical design rank, but attaches it to the wrong donors. If the shuffled
# arm recovers fibrosis just as well, A2b's recovery is an artifact of leaving
# variance in, not of protecting the right variance.
#
# CONTROL 2, leakage-free masking. A2b's ComBat saw every donor's stage,
# including donors any held-out evaluation later scores. A leave-one-cohort-out
# C-index computed on that matrix is contaminated. Here the held-out cohort is
# masked to `unstaged` before ComBat runs, so its labels never enter the
# transform, and it is then scored.
#
# The axis fit below is a reduced copy of fit_axis() from 02a: primary setting
# only (k=8, nPC=10, seed=master), no identifiability grid. A reduced copy can
# drift from its original, so the first thing this script does is refit the real
# A2b axis with it and assert Spearman >= 0.999 against arms/A2b_axis.tsv. If
# that assertion fails, nothing else here is trustworthy and the script stops.

suppressPackageStartupMessages({
  library(data.table)
  library(sva)
  library(slingshot)
  library(SingleCellExperiment)
})

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")
ARMS <- file.path(OUT, "arms")
N_PERM <- as.integer(Sys.getenv("CAB_N_PERM", "20"))
K <- 8L; NPC <- 10L

pre <- read_prespec()
meta <- load_manifest()
E_qn <- readRDS(file.path(ARMS, "E_qn.rds"))
pca_genes <- fread(file.path(ARMS, "pca_input_genes.tsv"))$gene_id
assert_true(identical(colnames(E_qn), meta$sample_id), "E_qn columns are not in manifest order")

STAGE_LEVELS <- c("F0", "F1", "F2", "F3", "F4", "unstaged")
as_protected <- function(stage_int) factor(
  ifelse(is.na(stage_int), "unstaged", paste0("F", stage_int)), levels = STAGE_LEVELS)

# ------------------------------------------------------------------ helpers --

# Masking a held-out cohort to `unstaged` can make the protected design rank
# deficient, and ComBat then aborts with a bare "covariate is confounded with
# batch". That happens when `unstaged` ends up covering whole cohorts and nothing
# else: the indicator becomes an exact sum of batch dummies. It is not a bug to
# route around. It means that after masking, the held-out cohort carries no stage
# information, so its batch effect and the unstaged effect are not separable and
# the fold is not estimable. Check it up front so the record says that, in those
# words, instead of dying.
design_full_rank <- function(label_int) {
  info <- droplevels(copy(meta)[, stage_protected := as_protected(label_int)])
  d <- cbind(model.matrix(~ -1 + info$dataset),
             model.matrix(~ inferred_sex + stage_protected, data = info)[, -1, drop = FALSE])
  qr(d)$rank == ncol(d)
}

combat_protected <- function(label_int) {
  info <- copy(meta)[, stage_protected := as_protected(label_int)]
  m <- ComBat(dat = E_qn, batch = as.character(info$dataset),
              mod = model.matrix(~ inferred_sex + stage_protected, data = info),
              par.prior = TRUE)
  dimnames(m) <- dimnames(E_qn)
  m
}

# Reduced copy of 02a::fit_axis, primary setting, root rule "histology" driven by
# whichever stage vector was protected. Kept line-for-line faithful on purpose.
fit_primary <- function(E, root_stage, seed) {
  set.seed(seed)
  p <- prcomp(t(E[pca_genes, ]), center = TRUE, scale. = FALSE)
  rd <- p$x[, seq_len(NPC), drop = FALSE]
  cl <- kmeans(rd, centers = K, nstart = 25, iter.max = 100)$cluster
  score <- vapply(split(root_stage, cl), function(v) mean(v, na.rm = TRUE), numeric(1))
  score[!is.finite(score)] <- Inf
  start <- names(score)[which.min(score)]
  sce <- SingleCellExperiment(assays = list(logcounts = t(rd)),
                              reducedDims = list(PCA = rd))
  sce <- tryCatch(
    slingshot(sce, clusterLabels = cl, reducedDim = "PCA", start.clus = start),
    error = function(e1) tryCatch(
      slingshot(sce, clusterLabels = cl, reducedDim = "PCA", start.clus = start,
                dist.method = "simple"),
      error = function(e2) NULL))
  if (is.null(sce)) return(list(pt = rep(NA_real_, nrow(rd)), fitted = FALSE))
  ptm <- slingPseudotime(sce); w <- slingCurveWeights(sce)
  ptm[!is.finite(ptm)] <- NA_real_; w[!is.finite(w)] <- 0; w[is.na(ptm)] <- 0
  denom <- rowSums(w)
  pt <- ifelse(denom > 0, rowSums(ptm * w, na.rm = TRUE) / denom, NA_real_)
  if (cor(pt, p$x[, 1], use = "complete.obs") < 0) pt <- -pt
  if (mean(pt[cl == start], na.rm = TRUE) > mean(pt, na.rm = TRUE)) pt <- -pt
  list(pt = pt, fitted = TRUE, n_lineages = ncol(ptm))
}

# Concordance of axis with stage among pairs of unequal stage. Ties in the axis
# count half, as in Harrell's C.
c_index <- function(axis, stage) {
  ok <- is.finite(axis) & !is.na(stage)
  a <- axis[ok]; s <- stage[ok]
  if (length(a) < 3L || length(unique(s)) < 2L) return(NA_real_)
  n <- length(a); conc <- 0; tied <- 0; tot <- 0
  for (i in seq_len(n - 1L)) {
    j <- (i + 1L):n
    d <- s[j] - s[i]
    use <- d != 0
    if (!any(use)) next
    sgn <- sign(d[use]) * sign(a[j][use] - a[i])
    conc <- conc + sum(sgn > 0); tied <- tied + sum(sgn == 0); tot <- tot + sum(use)
  }
  (conc + 0.5 * tied) / tot
}

score_axis <- function(pt, label, extra = list()) {
  ok <- !is.na(meta$fibrosis_stage)
  as.data.table(c(list(
    label = label,
    n_scored = sum(is.finite(pt)),
    spearman_vs_true_fibrosis = cor(pt[ok], meta$fibrosis_stage[ok],
                                    method = "spearman", use = "complete.obs"),
    c_index_all_staged = c_index(pt[ok], meta$fibrosis_stage[ok]),
    kruskal_eta2_on_cohort = {
      kw <- kruskal.test(pt ~ meta$dataset)
      (unname(kw$statistic) - uniqueN(meta$dataset) + 1) / (length(pt) - uniqueN(meta$dataset))
    }), extra))
}

# --------------------------------------------------------- drift assertion ---
log_step("refitting the real A2b axis with the reduced implementation")
E_cbd <- readRDS(file.path(ARMS, "E_cbd.rds"))
real <- fit_primary(E_cbd, meta$fibrosis_stage, pre$seeds$master)
assert_true(isTRUE(real$fitted), "The reduced implementation could not fit the real A2b axis")

a2b <- fread(file.path(ARMS, "A2b_axis.tsv"))
assert_true(identical(a2b$sample_id, meta$sample_id), "A2b_axis.tsv is not in manifest order")
drift <- cor(real$pt, a2b$axis_raw, method = "spearman", use = "complete.obs")
log_step("reduced-vs-arm agreement: Spearman ", format(drift, digits = 6))
assert_true(is.finite(drift) && drift >= 0.999, paste0(
  "The reduced fit does not reproduce arm A2b (Spearman ", format(drift, digits = 6),
  "). 02a and this script have drifted; nothing below is interpretable."))

results <- list(score_axis(real$pt, "A2b_real_protection",
                           list(arm = "real", perm_index = NA_integer_,
                                held_out_cohort = NA_character_)))

# ------------------------------------------------- control 1: permuted null ---
log_step("permuted-protection null: ", N_PERM, " within-cohort shuffles")
for (i in seq_len(N_PERM)) {
  set.seed(pre$seeds$perm_base + i)
  perm <- meta$fibrosis_stage
  for (g in unique(meta$dataset)) {
    idx <- which(meta$dataset == g)
    perm[idx] <- perm[sample(idx)]
  }
  assert_true(identical(sort(table(addNA(perm), meta$dataset)),
                        sort(table(addNA(meta$fibrosis_stage), meta$dataset))),
              "Within-cohort permutation changed the cohort-by-stage margins")
  Ep <- combat_protected(perm)
  fp <- fit_primary(Ep, perm, pre$seeds$master)
  results[[length(results) + 1L]] <- score_axis(
    fp$pt, sprintf("perm_%02d", i),
    list(arm = "permuted_protection", perm_index = i, held_out_cohort = NA_character_))
  log_step(sprintf("  perm %2d/%d  rho_vs_true_fibrosis = %+.4f", i, N_PERM,
                   results[[length(results)]]$spearman_vs_true_fibrosis))
  rm(Ep); gc(verbose = FALSE)
}

# ------------------------------------------- control 2: leakage-free masking ---
log_step("leakage-free masking: one POP-C cohort at a time")
loco <- list()
for (co in POP_C_COHORTS) {
  masked <- meta$fibrosis_stage
  masked[meta$dataset == co] <- NA_integer_
  idx0 <- which(meta$dataset == co & !is.na(meta$fibrosis_stage))
  if (!design_full_rank(masked)) {
    loco[[co]] <- data.table(
      held_out_cohort = co, n_held_out = length(idx0),
      spearman_within_heldout = NA_real_, c_index_within_heldout = NA_real_,
      estimable = FALSE,
      reason = paste0("Masking ", co, " leaves `unstaged` covering whole cohorts only, so the ",
                      "indicator is an exact sum of batch dummies and the protected design is ",
                      "rank deficient. After masking, ", co, " carries no stage information, so ",
                      "its batch effect and the unstaged effect are not separable. The fold is ",
                      "not estimable; it is recorded, not routed around."))
    log_step("  masked ", co, " NOT ESTIMABLE (rank-deficient protected design)")
    next
  }
  Em <- combat_protected(masked)
  fm <- fit_primary(Em, masked, pre$seeds$master)
  idx <- idx0
  # Rank within the held-out cohort before scoring, as every other evaluation in
  # this workstream does, so a cohort mean-shift cannot masquerade as signal.
  r <- rank_within(fm$pt[idx], rep(co, length(idx)))
  loco[[co]] <- data.table(
    held_out_cohort = co, n_held_out = length(idx),
    spearman_within_heldout = cor(r, meta$fibrosis_stage[idx],
                                  method = "spearman", use = "complete.obs"),
    c_index_within_heldout = c_index(r, meta$fibrosis_stage[idx]),
    estimable = TRUE, reason = NA_character_)
  results[[length(results) + 1L]] <- score_axis(
    fm$pt, paste0("masked_", co),
    list(arm = "leakage_free_masked", perm_index = NA_integer_, held_out_cohort = co))
  log_step(sprintf("  masked %-10s n=%3d  rho=%+.4f  C=%.4f", co, length(idx),
                   loco[[co]]$spearman_within_heldout, loco[[co]]$c_index_within_heldout))
  rm(Em); gc(verbose = FALSE)
}

# ------------------------------------------------------------------ verdict ---
all_dt <- rbindlist(results, fill = TRUE)
write_tsv_once(all_dt, file.path(ARMS, "A2b_protection_null_axes.tsv"))
write_tsv_once(rbindlist(loco), file.path(ARMS, "A2b_leakage_free_loco.tsv"))

nullrho <- all_dt[arm == "permuted_protection", spearman_vs_true_fibrosis]
realrho <- all_dt[arm == "real", spearman_vs_true_fibrosis]
# One-sided empirical p with the +1 correction, so p is never reported as 0 from
# a finite permutation set.
emp_p <- (1 + sum(nullrho >= realrho)) / (1 + length(nullrho))

verdict <- data.table(
  arm = "A2b",
  status = "post-hoc diagnostic; sealed Q1 gate not reopened",
  n_perm = length(nullrho),
  real_spearman_vs_fibrosis = realrho,
  null_median_spearman = median(nullrho),
  null_q95_spearman = quantile(nullrho, 0.95, names = FALSE),
  null_max_spearman = max(nullrho),
  empirical_p_one_sided = emp_p,
  exceeds_permuted_protection = realrho > quantile(nullrho, 0.95, names = FALSE),
  a1_spearman_vs_fibrosis_reference = 0.0519155215026156,
  a3_spearman_vs_fibrosis_reference = 0.568054879376621,
  n_loco_folds_estimable = sum(rbindlist(loco)$estimable),
  n_loco_folds_attempted = length(loco),
  mean_leakage_free_c_index = mean(rbindlist(loco)$c_index_within_heldout, na.rm = TRUE),
  interpretation = paste(
    "Protecting stage in ComBat preserves stage variance by construction, so the",
    "real arm is only interpretable relative to the permuted-protection null.",
    "A real rho above the null q95 means the ORIGINAL ComBat was removing real",
    "between-cohort disease signal. A real rho inside the null means it was not,",
    "and the Q1 null is not an over-correction artifact."))
write_tsv_once(verdict, file.path(ARMS, "A2b_verdict.tsv"))
print(t(verdict))
log_step("A2B_NULL_COMPLETE")
