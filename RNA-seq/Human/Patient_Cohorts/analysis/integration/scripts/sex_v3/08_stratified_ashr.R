#!/usr/bin/env Rscript
# sex_v3/08_stratified_ashr.R
# ---------------------------------------------------------------------------
# Layer 2 — Orthogonal sex-stratified classifier (mashr-independent).
# Fits dream separately on F-only and M-only subsets, applies ashr per arm,
# computes retrospective power, joins with cached interaction-term test
# (sex_interaction_dream_v3.csv) → class_strat per gene.
#
# Reads:  intermediates/sex_v3_input.rds (Module 01)
#         intermediates/sva_factors.rds  (Module 02)
#         intermediates/age_mi.rds       (Module 03; uses first imputation)
#         sex_interaction_dream_v3.csv   (cached full-sample interaction)
# Writes: stratified_ashr_v4.csv         (per-gene class_strat + diagnostics)
# ---------------------------------------------------------------------------
set.seed(42)
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
  library(variancePartition); library(BiocParallel); library(ashr)
})

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SEXV3 <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                   "results/integration/sex_v3")
IDIR  <- file.path(SEXV3, "intermediates")
OUT_CSV <- file.path(SEXV3, "stratified_ashr_v4.csv")

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
param <- if (ncpus > 1) MulticoreParam(workers = ncpus, RNGseed = 42L) else SerialParam()
cat("Using", ncpus, "CPU cores\n")

inp    <- readRDS(file.path(IDIR, "sex_v3_input.rds"))
sva    <- readRDS(file.path(IDIR, "sva_factors.rds"))
age_mi <- readRDS(file.path(IDIR, "age_mi.rds"))

dge   <- inp$dge
info0 <- as.data.table(inp$meta, keep.rownames = "sample_id")
sv_mat <- sva$sv
sv_df  <- as.data.frame(sv_mat[match(info0$sample_id, rownames(sv_mat)), , drop = FALSE])
colnames(sv_df) <- paste0("SV", seq_len(ncol(sv_df)))
info0[, age_imputed := age_mi$age_list[[1]][match(sample_id, age_mi$sample_ids)]]
info0 <- cbind(info0, sv_df)
stopifnot(identical(info0$sample_id, colnames(dge)))

cat("\nSample counts per sex × group × cohort:\n")
print(ftable(info0$dataset, info0$inferred_sex, info0$group_binary))

# Per-arm formula: drop inferred_sex (constant within arm).
sv_terms <- paste(colnames(sv_df), collapse = " + ")
form_str <- paste0(
  "~ group_binary",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed + ", sv_terms,
  " + (1 + group_binary | dataset)"
)
form_str_int <- paste0(
  "~ group_binary",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed + ", sv_terms,
  " + (1 | dataset)"
)
form_slope <- as.formula(form_str)
form_int   <- as.formula(form_str_int)
cat("\nPer-arm formula (slope, primary):\n  ", form_str, "\n")

fit_one_arm <- function(sex_label) {
  cat("\n========== Arm:", sex_label, "==========\n")
  keep <- info0$inferred_sex == sex_label
  info_a <- as.data.frame(info0[keep])
  rownames(info_a) <- info_a$sample_id
  dge_a <- dge[, info_a$sample_id]
  dge_a <- calcNormFactors(dge_a)
  cat("  N samples:", ncol(dge_a),
      "  N Disease:", sum(info_a$group_binary == "Disease"),
      "  N Control:", sum(info_a$group_binary == "Control"), "\n")

  cat("  voomWithDreamWeights (slope)...\n")
  t0 <- Sys.time()
  v <- tryCatch(
    suppressWarnings(voomWithDreamWeights(dge_a, form_slope, info_a,
                                          BPPARAM = param, useWeights = TRUE)),
    error = function(e) {
      cat("    slope voom failed; intercept fallback:", conditionMessage(e), "\n")
      suppressWarnings(voomWithDreamWeights(dge_a, form_int, info_a,
                                            BPPARAM = param, useWeights = TRUE))
    })
  cat("    voom elapsed:", format(Sys.time() - t0), "\n")
  cat("  dream (slope)...\n")
  t0 <- Sys.time()
  fit <- tryCatch(
    suppressWarnings(dream(v, form_slope, info_a, BPPARAM = param, useWeights = TRUE)),
    error = function(e) {
      cat("    slope dream failed; intercept fallback:", conditionMessage(e), "\n")
      suppressWarnings(dream(v, form_int, info_a, BPPARAM = param, useWeights = TRUE))
    })
  cat("    dream elapsed:", format(Sys.time() - t0), "\n")

  coef_name <- "group_binaryDisease"
  beta <- fit$coefficients[, coef_name]
  se   <- (fit$stdev.unscaled * fit$sigma)[, coef_name]
  list(gene = rownames(fit$coefficients), beta = beta, se = se,
       n_dis = sum(info_a$group_binary == "Disease"),
       n_ctl = sum(info_a$group_binary == "Control"),
       pooled_sd = mean(fit$sigma, na.rm = TRUE))
}

res_F <- fit_one_arm("F")
res_M <- fit_one_arm("M")

# ashr per arm
cat("\nashr per arm...\n")
ash_F <- ash(res_F$beta, res_F$se, mixcompdist = "normal", mode = "estimate")
ash_M <- ash(res_M$beta, res_M$se, mixcompdist = "normal", mode = "estimate")
lfsr_F <- ash_F$result$lfsr; pm_F <- ash_F$result$PosteriorMean
lfsr_M <- ash_M$result$lfsr; pm_M <- ash_M$result$PosteriorMean

# Retrospective power per gene (Cohen's d ≈ |beta_F| / pooled_sd; t-test power)
cat("Retrospective power per gene...\n")
cohens_d <- abs(res_F$beta) / max(res_F$pooled_sd, 1e-6)
# Closed-form two-sample t power (Welch-like, equal-variance approx) per gene.
# ncp = d * sqrt(n1*n2/(n1+n2)); reject if |t| > qt(1-α/2, df). Avoid pwr dep.
n1 <- res_M$n_ctl; n2 <- res_M$n_dis
df <- n1 + n2 - 2
crit <- qt(1 - 0.05 / 2, df = df)
power_M <- vapply(cohens_d, function(d) {
  if (!is.finite(d) || d <= 0) return(NA_real_)
  ncp <- d * sqrt(n1 * n2 / (n1 + n2))
  1 - pt(crit, df = df, ncp = ncp) + pt(-crit, df = df, ncp = ncp)
}, numeric(1))

# Join cached interaction term
int_csv <- file.path(SEXV3, "sex_interaction_dream_v3.csv")
if (file.exists(int_csv)) {
  int_dt <- fread(int_csv)
  setnames(int_dt, old = intersect(c("gene", "gene_name", "gene_id"), names(int_dt)),
           new = rep("gene", length(intersect(c("gene", "gene_name", "gene_id"), names(int_dt)))),
           skip_absent = TRUE)
  padj_int_col <- intersect(c("padj_interaction_ashr", "padj_interaction", "padj_int"),
                            names(int_dt))[1]
  if (is.na(padj_int_col)) padj_int_col <- intersect(c("padj_BH", "padj", "FDR"),
                                                      names(int_dt))[1]
  cat("  joined interaction term from", int_csv, "via column", padj_int_col, "\n")
  int_dt <- int_dt[, .(gene, padj_interaction = get(padj_int_col))]
} else {
  cat("  WARNING: cached interaction csv missing; padj_interaction set NA\n")
  int_dt <- data.table(gene = res_F$gene, padj_interaction = NA_real_)
}

# Build per-gene table
stopifnot(identical(res_F$gene, res_M$gene))
dt <- data.table(
  gene = res_F$gene,
  beta_F_strat = res_F$beta, se_F_strat = res_F$se,
  beta_M_strat = res_M$beta, se_M_strat = res_M$se,
  lfsr_F_strat = lfsr_F, lfsr_M_strat = lfsr_M,
  PM_F_strat = pm_F, PM_M_strat = pm_M,
  cohens_d_F = cohens_d, power_M_at_F = power_M
)
dt <- merge(dt, int_dt, by = "gene", all.x = TRUE)

# Decision rule
same_sign <- sign(dt$beta_F_strat) == sign(dt$beta_M_strat) &
             dt$beta_F_strat != 0 & dt$beta_M_strat != 0
opp_sign <- sign(dt$beta_F_strat) != sign(dt$beta_M_strat) &
            dt$beta_F_strat != 0 & dt$beta_M_strat != 0
padj_int <- ifelse(is.na(dt$padj_interaction), 1, dt$padj_interaction)

class_strat <- rep("uncertain", nrow(dt))
class_strat[dt$lfsr_F_strat < 0.05 & dt$lfsr_M_strat < 0.05 & same_sign] <- "concordant"
class_strat[dt$lfsr_F_strat < 0.05 & dt$lfsr_M_strat < 0.05 & opp_sign]  <- "divergent"
class_strat[dt$lfsr_F_strat < 0.05 & dt$lfsr_M_strat > 0.20 & padj_int < 0.10] <- "F_only"
class_strat[dt$lfsr_M_strat < 0.05 & dt$lfsr_F_strat > 0.20 & padj_int < 0.10] <- "M_only"
class_strat[class_strat == "uncertain" &
            dt$lfsr_M_strat > 0.05 & dt$lfsr_M_strat < 0.20 &
            !is.na(dt$power_M_at_F) & dt$power_M_at_F < 0.5] <- "M_underpowered"
dt[, class_strat := class_strat]

cat("\n== Layer 2 stratified-ashr classification ==\n")
print(table(dt$class_strat))

fwrite(dt, OUT_CSV)
cat("\nWrote:", OUT_CSV, "  rows:", nrow(dt), "\n")
