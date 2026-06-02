#!/usr/bin/env Rscript
# sex_v3/14_calibration_50seed.R
# ---------------------------------------------------------------------------
# Pillar 5 — 50-seed calibration. Forks 07a2_v5_synthetic.R but bumps SEEDS
# from 10 -> 50 (4 effect_mag values × 50 seeds = 200 fits).
#
# Logic is IDENTICAL to 07a2_v5_synthetic.R: count-scale NB injection of
# F_only / M_only / concordant / anti_correlated synthetic genes, refit
# dream M2 (current v5 formula F5 = `(1 + group_binary | dataset)`),
# apply the interaction-test v5 decision tree, score TPR/FDR per pattern.
#
# Output: intermediates/calibration_v6/synth_seed<S>_eff<E>.rds   (NEW DIR)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(reformulas); library(lme4); library(data.table); library(edgeR)
})
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({ unlockBinding(fn, ns_lme4)
          assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
          lockBinding(fn, ns_lme4) }, silent = TRUE)
  }
}
suppressPackageStartupMessages({
  library(variancePartition); library(BiocParallel)
})

SEED       <- as.integer(Sys.getenv("SEED", "1"))
EFFECT_MAG <- as.numeric(Sys.getenv("EFFECT_MAG", "0.8"))
N_PER_PATT <- as.integer(Sys.getenv("N_PER_PATT", "250"))
N_CARRIER  <- as.integer(Sys.getenv("N_CARRIER", "1000"))
set.seed(SEED)

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SEXV3 <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                   "results/integration/sex_v3")
IDIR  <- file.path(SEXV3, "intermediates")

# Contrast routing (sex_v3_utils.R::contrast_paths) — overrides SEXV3 / IDIR
# when CONTRAST_NAME != "disease_vs_ctrl"; default preserves legacy layout.
if (!exists("contrast_paths")) source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT",
             "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3/sex_v3_utils.R"))
.cpaths <- contrast_paths()
SEXV3 <- .cpaths$sexv3
IDIR  <- .cpaths$idir
dir.create(IDIR, recursive = TRUE, showWarnings = FALSE)
OUTDIR <- file.path(IDIR, "calibration_v6")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
out_rds <- file.path(OUTDIR, sprintf("synth_seed%d_eff%.2f.rds", SEED, EFFECT_MAG))

cat("=== 14_calibration_50seed (Pillar 5) ===\n")
cat("SEED:", SEED, "  EFFECT_MAG:", EFFECT_MAG, "  K:", N_PER_PATT, "\n")

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
param <- if (ncpus > 1) MulticoreParam(workers = ncpus, RNGseed = SEED) else SerialParam()

inp    <- readRDS(file.path(IDIR, "sex_v3_input.rds"))
sva    <- readRDS(file.path(IDIR, "sva_factors.rds"))
age_mi <- readRDS(file.path(IDIR, "age_mi.rds"))

info0 <- as.data.table(inp$meta, keep.rownames = "sample_id")
sv_mat <- sva$sv
sv_df  <- as.data.frame(sv_mat[match(info0$sample_id, rownames(sv_mat)), ,
                               drop = FALSE])
colnames(sv_df) <- paste0("SV", seq_len(ncol(sv_df)))
info0[, age_imputed := age_mi$age_list[[1]][match(sample_id, age_mi$sample_ids)]]
info0 <- cbind(info0, sv_df)
info <- as.data.frame(info0); rownames(info) <- info$sample_id
dge <- inp$dge[, info$sample_id]
dge <- calcNormFactors(dge)
n_samp <- ncol(dge)
cat("Samples:", n_samp, "  Genes:", nrow(dge), "\n")

cat("Estimating per-gene NB dispersion...\n")
group_for_disp <- paste(info$inferred_sex, info$group_binary, sep = ":")
design_disp <- model.matrix(~ group_for_disp)
dge <- estimateDisp(dge, design_disp, robust = TRUE)
phi <- dge$tagwise.dispersion

real_genes <- rownames(dge$counts)
set.seed(SEED * 7919L)
carrier_idx <- sample.int(length(real_genes),
                          size = min(N_CARRIER, length(real_genes)))
carrier_genes <- real_genes[carrier_idx]

patterns <- c("F_only", "M_only", "concordant", "anti_correlated")
set.seed(SEED * 31L)
donor_idx <- sample.int(length(real_genes), size = N_PER_PATT * 4)
donor_genes <- real_genes[donor_idx]

truth <- data.table(synth_id = paste0("synth_", seq_len(N_PER_PATT * 4)),
                    donor_gene = donor_genes,
                    pattern = rep(patterns, each = N_PER_PATT))
draw_beta <- function(pat) {
  bF <- runif(1, 0.5 * EFFECT_MAG, EFFECT_MAG) * sample(c(-1, 1), 1)
  bM <- runif(1, 0.5 * EFFECT_MAG, EFFECT_MAG) * sample(c(-1, 1), 1)
  # Concordant pattern: F and M have the same sign and similar magnitude.
  # Earlier `sign(bM) *` multiplier collapsed half of "concordant" truth rows
  # to anti-correlated when the independently-drawn bM happened to be negative.
  switch(pat,
         F_only          = c(bF, 0),
         M_only          = c(0, bM),
         concordant      = c(bF, bF * runif(1, 0.7, 1.3)),
         anti_correlated = c(bF, -bF * runif(1, 0.7, 1.3)))
}
betas <- t(vapply(truth$pattern, draw_beta, numeric(2)))
truth[, beta_F_true := betas[, 1]]
truth[, beta_M_true := betas[, 2]]

cat("Injecting", nrow(truth), "synthetic genes at count scale...\n")
real_cnt  <- dge$counts
synth_cnt <- matrix(0L, nrow = nrow(truth), ncol = n_samp,
                    dimnames = list(truth$synth_id, colnames(real_cnt)))

is_F <- info$inferred_sex == "F"; is_M <- info$inferred_sex == "M"
is_dis <- info$group_binary == "Disease"

for (i in seq_len(nrow(truth))) {
  gd <- truth$donor_gene[i]
  mu_real <- pmax(real_cnt[gd, ], 0.5)
  log2fc_per_sample <- rep(0, n_samp)
  log2fc_per_sample[is_F & is_dis] <- truth$beta_F_true[i]
  log2fc_per_sample[is_M & is_dis] <- truth$beta_M_true[i]
  mu_inj <- mu_real * (2 ^ log2fc_per_sample)
  phi_g <- phi[match(gd, real_genes)]
  if (!is.finite(phi_g) || phi_g <= 0) phi_g <- 0.1
  synth_cnt[i, ] <- as.integer(rnbinom(n_samp, mu = mu_inj, size = 1 / phi_g))
}

carrier_cnt <- real_cnt[carrier_genes, , drop = FALSE]
all_cnt <- rbind(synth_cnt, carrier_cnt)
dge_synth <- DGEList(all_cnt, samples = dge$samples)
dge_synth$samples$norm.factors <- dge$samples$norm.factors
dge_synth <- calcNormFactors(dge_synth)

sv_terms <- paste(colnames(sv_df), collapse = " + ")
form_str_slope <- paste0(
  "~ group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed + ", sv_terms,
  " + (1 + group_binary | dataset)"
)
form_str_int <- paste0(
  "~ group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed + ", sv_terms,
  " + (1 | dataset)"
)
form_slope     <- as.formula(form_str_slope)
form_intercept <- as.formula(form_str_int)

cat("voomWithDreamWeights + dream on synthetic DGEList...\n")
t0 <- Sys.time()
v <- tryCatch(
  suppressWarnings(voomWithDreamWeights(dge_synth, form_slope, info,
                                        BPPARAM = param, useWeights = TRUE)),
  error = function(e) {
    cat("  slope voom fail; fallback intercept:", conditionMessage(e), "\n")
    suppressWarnings(voomWithDreamWeights(dge_synth, form_intercept, info,
                                          BPPARAM = param, useWeights = TRUE))
  })
fit <- tryCatch(
  suppressWarnings(dream(v, form_slope, info, BPPARAM = param,
                         useWeights = TRUE)),
  error = function(e) {
    cat("  slope dream fail; fallback intercept:", conditionMessage(e), "\n")
    suppressWarnings(dream(v, form_intercept, info, BPPARAM = param,
                           useWeights = TRUE))
  })
cat("  dream elapsed:", format(Sys.time() - t0), "\n")

all_coefs <- colnames(fit$coefficients)
gcoef <- setdiff(grep("^group_binary", all_coefs, value = TRUE),
                 grep(":", all_coefs, value = TRUE))[1]
icoef <- grep(":", grep("^group_binary", all_coefs, value = TRUE), value = TRUE)[1]
bF <- fit$coefficients[, gcoef]
bI <- fit$coefficients[, icoef]
bM <- bF + bI
raw_se <- fit$stdev.unscaled * fit$sigma
sF <- raw_se[, gcoef]; sI <- raw_se[, icoef]
# dream MArrayLM2 stores per-gene vcov in cov.coefficients.list (single
# cov.coefficients is the limma scalar API and returns NULL for dream LMMs,
# silently zeroing cov_FI and biasing Var(β_M) downward).
if (!"cov.coefficients.list" %in% names(fit) ||
    length(fit$cov.coefficients.list) != length(bF)) {
  stop("fit$cov.coefficients.list missing or wrong length: ",
       length(fit$cov.coefficients.list), " vs ", length(bF))
}
ccl <- fit$cov.coefficients.list
cov_FI <- vapply(seq_along(ccl), function(g) {
  M <- ccl[[g]]
  if (gcoef %in% rownames(M) && icoef %in% colnames(M))
    as.numeric(M[gcoef, icoef])
  else NA_real_
}, numeric(1))
if (any(is.na(cov_FI))) stop(sum(is.na(cov_FI)), " genes have missing Cov(β_F, β_int)")
sM <- sqrt(pmax(sI^2 + sF^2 + 2 * cov_FI, .Machine$double.eps))

t_int <- bI / sI
p_int <- 2 * pnorm(-abs(t_int))

all_genes <- rownames(fit$coefficients)
g_is_synth <- all_genes %in% truth$synth_id
q_int <- rep(NA_real_, length(all_genes))
q_int[g_is_synth]  <- p.adjust(p_int[g_is_synth],  "BH")
q_int[!g_is_synth] <- p.adjust(p_int[!g_is_synth], "BH")

# posterior-P-like statistic (mirrors P(call ∈ class | data)) — used by the
# panel-C calibration deciles. Compute 1 - q_int for in-class genes and
# q_int for out-of-class as a rough monotone confidence transform.
post_p <- 1 - q_int
post_p[!is.finite(post_p)] <- NA_real_

sig_F <- abs(bF / sF) > 1.96
sig_M <- abs(bM / sM) > 1.96

n_g <- length(all_genes)
v5_class <- rep("Not_DEG", n_g)
for (i in seq_len(n_g)) {
  qi <- q_int[i]; bFi <- bF[i]; bMi <- bM[i]
  if (!is.finite(qi)) next
  if (qi >= 0.20) {
    if (sig_F[i] && sig_M[i] && sign(bFi) == sign(bMi)) {
      v5_class[i] <- "Concordant"
    } else if (sig_F[i] || sig_M[i]) {
      v5_class[i] <- "Concordant"
    } else {
      v5_class[i] <- "Not_DEG"
    }
  } else {
    if (sign(bFi) != sign(bMi) && sig_F[i] && sig_M[i]) {
      v5_class[i] <- "Divergent"
    } else if (sign(bFi) != sign(bMi) && (sig_F[i] || sig_M[i])) {
      v5_class[i] <- "Divergent_one_sided"
    } else if (sign(bFi) == sign(bMi)) {
      if (abs(bFi) >= 2 * abs(bMi) && sig_F[i] && !sig_M[i]) {
        v5_class[i] <- "Female_biased"
      } else if (abs(bMi) >= 2 * abs(bFi) && sig_M[i] && !sig_F[i]) {
        v5_class[i] <- "Male_biased"
      } else {
        v5_class[i] <- "Sex_modifier"
      }
    } else {
      v5_class[i] <- "Not_DEG"
    }
  }
}

pat_call <- rep("uncertain", n_g)
pat_call[v5_class == "Female_biased"] <- "F_only"
pat_call[v5_class == "Male_biased"]   <- "M_only"
pat_call[v5_class %in% c("Divergent", "Divergent_one_sided")] <- "anti_correlated"
pat_call[v5_class %in% c("Concordant", "Sex_modifier")] <- "concordant"
pat_call[v5_class == "Not_DEG"] <- "uncertain"

out <- data.table(gene = all_genes,
                  beta_F_est = bF, beta_M_est = bM, se_F = sF, se_M = sM,
                  beta_int = bI, se_int = sI,
                  t_int = t_int, p_int = p_int, q_int = q_int, post_p = post_p,
                  sig_F = sig_F, sig_M = sig_M,
                  v5_class = v5_class, gated = pat_call)
out <- merge(out, truth[, .(gene = synth_id, true_pattern = pattern,
                            beta_F_true, beta_M_true)],
             by = "gene", all.x = TRUE)
out[is.na(true_pattern), true_pattern := "carrier_null"]

# Per-pattern TPR / FDR
metrics <- list()
for (pat in patterns) {
  truth_pat <- out[true_pattern == pat]
  if (nrow(truth_pat) == 0) next
  tp <- sum(truth_pat$gated == pat)
  fn <- sum(truth_pat$gated != pat & truth_pat$gated != "uncertain")
  un <- sum(truth_pat$gated == "uncertain")
  call_pat <- out[gated == pat]
  fp <- sum(call_pat$true_pattern != pat & call_pat$true_pattern != "carrier_null")
  metrics[[pat]] <- list(
    pattern = pat,
    n_true = nrow(truth_pat),
    n_called = nrow(call_pat),
    tpr = tp / nrow(truth_pat),
    fdr = if (nrow(call_pat) > 0) fp / nrow(call_pat) else 0,
    uncertain_rate = un / nrow(truth_pat))
}
cat("\n== Per-pattern metrics (v6 50-seed calibration) ==\n")
print(rbindlist(metrics))

# Confusion matrix true x called (4x5 with uncertain bucket)
conf <- with(out, table(true_pattern, gated))
cat("\n== Confusion matrix (true x called) ==\n")
print(conf)

saveRDS(list(seed = SEED, effect_mag = EFFECT_MAG,
             metrics = rbindlist(metrics), out = out, truth = truth,
             confusion = conf,
             classifier = "v6_50seed_calibration"),
        out_rds)
cat("\nWrote:", out_rds, "\n")
