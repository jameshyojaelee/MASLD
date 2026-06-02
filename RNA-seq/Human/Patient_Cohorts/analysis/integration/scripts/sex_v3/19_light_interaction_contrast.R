#!/usr/bin/env Rscript
# sex_v3/19_light_interaction_contrast.R
# ---------------------------------------------------------------------------
# Tier B — light fixed-effect sex × disease interaction for the two control-
# bearing contrasts (MASL-vs-Ctrl and MASH-vs-Ctrl).
#
# Why this exists (per `/docs/superpowers/specs/2026-05-18-...`):
# These contrasts collapse to 3 usable cohorts (GSE126848 has 0 F-Ctrl) with
# severe sex imbalance (up to 17.7× F:M in the Ctrl arm). N=152 (MASL-vs-Ctrl)
# and N=294 (MASH-vs-Ctrl) place TPR(M) at δ=0.8 around 0.07 and 0.10
# respectively — running the full 9-pillar cross-pillar pipeline would
# produce a calibrated NULL with no candidate-list survivors and burn ~8 hr
# of SLURM for known-null output.
#
# Instead this script runs one dream fit with a SIMPLE-intercept random
# effect on `dataset` (the F2 random-slope-on-sex×disease formula used in
# Pillar 1 is not estimable when controls of one sex are absent in any
# cohort) plus per-cohort fixed-effect Cochran's Q heterogeneity. Outputs a
# single CSV (gene-level β_int, padj, cohort-Q) that the stage-attribution
# step in Module 20 reads alongside the canonical Disease-vs-Ctrl pipeline
# and the MASH-vs-MASL Tier A pipeline.
#
# Reads (from contrast-specific intermediates dir):
#   - sex_v3_input.rds, sva_factors.rds, age_mi.rds (from Modules 01–03)
# Writes (to contrast-specific SEXV3 dir):
#   - interaction_fixed.csv  (gene, beta_int, se_int, t_int, p_int,
#                              padj_int_5k, padj_int_full, cohort_Q_pval,
#                              I2_pct, n_cohorts_identifiable)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma); library(BiocParallel)
  library(variancePartition)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE,
                 "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                 "scripts/sex_v3/sex_v3_utils.R"))

CSPEC  <- contrast_spec()
CPATHS <- contrast_paths()
SEXV3  <- CPATHS$sexv3
IDIR   <- CPATHS$idir

if (CSPEC$name %in% c("disease_vs_ctrl", "mash_vs_masl")) {
  stop("19_light_interaction_contrast.R is for the underpowered contrasts ",
       "only (masl_vs_ctrl, mash_vs_ctrl). For '", CSPEC$name,
       "' run the full Tier A pipeline via run_v6_pipeline.sh.")
}

cat("============================================================\n")
cat("Tier B light interaction — contrast = ", CSPEC$pretty_name, "\n", sep = "")
cat("  IDIR   = ", IDIR, "\n", sep = "")
cat("  SEXV3  = ", SEXV3, "\n", sep = "")
cat("============================================================\n")

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------
inp    <- readRDS(file.path(IDIR, "sex_v3_input.rds"))
sva    <- readRDS(file.path(IDIR, "sva_factors.rds"))
age_mi <- readRDS(file.path(IDIR, "age_mi.rds"))

dge   <- inp$dge
info0 <- inp$meta
sv_mat <- sva$sv
n_sv   <- ncol(sv_mat)
sv_df  <- as.data.frame(sv_mat); colnames(sv_df) <- paste0("SV", seq_len(n_sv))

info <- info0
info$age_imputed <- age_mi$age_list[[1]]
info <- cbind(info, sv_df)
rownames(info) <- info$sample_id
stopifnot(identical(rownames(info), colnames(dge)))

cat("Samples:", ncol(dge), "  Genes:", nrow(dge), "  SVs:", n_sv, "\n\n")
cat("Sex × group_binary contingency:\n")
print(table(info$inferred_sex, info$group_binary))
cat("\nCohort × group_binary contingency:\n")
print(table(info$dataset, info$group_binary))

# ---------------------------------------------------------------------------
# Identifiable cohorts: at least one F + one M in each arm
# ---------------------------------------------------------------------------
xt_per_cohort <- with(info,
                      tapply(seq_along(sample_id), dataset,
                             function(idx) {
                               x <- table(info$inferred_sex[idx],
                                          info$group_binary[idx])
                               all(x > 0) && all(dim(x) == c(2, 2))
                             }))
cohort_identifiable <- names(xt_per_cohort)[xt_per_cohort == TRUE]
n_id <- length(cohort_identifiable)
cat(sprintf("\nIdentifiable cohorts (≥1 F + ≥1 M in both arms): %d / %d\n",
            n_id, length(xt_per_cohort)))
cat("  ", paste(cohort_identifiable, collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# Pooled fixed-effect dream fit (~ group_binary * inferred_sex + ... + (1 | dataset))
# ---------------------------------------------------------------------------
cat("\n[1] voom + dream (simple intercept random effect on dataset)...\n")
sv_terms <- paste(colnames(sv_df), collapse = " + ")
form_str <- paste0("~ group_binary * inferred_sex",
                   " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
                   " + age_imputed + ", sv_terms,
                   " + (1 | dataset)")
form <- as.formula(form_str)
cat("  Formula: ", form_str, "\n\n")

param <- SerialParam()
ncpu <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
if (!is.na(ncpu) && ncpu > 1L) param <- SnowParam(ncpu)

# Pre-filter low-expressed genes (mean-variance trend fails when too many
# near-zero counts populate the lowess input; lowess delta=0 when all sx
# values are tied at the minimum count level).
keep_genes <- filterByExpr(dge, design = NULL, group = info$group_binary,
                            min.count = 10, min.total.count = 15)
cat("filterByExpr retained:", sum(keep_genes), "/", length(keep_genes), "genes\n")
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

t0 <- Sys.time()
v_voom <- tryCatch(
  suppressWarnings(
    voomWithDreamWeights(dge, form, info, BPPARAM = param, useWeights = TRUE)),
  error = function(e) {
    cat("  voomWithDreamWeights FAILED (", conditionMessage(e),
        ") — falling back to plain voom\n", sep = "")
    NULL
  })
if (is.null(v_voom)) {
  # Fallback path: plain voom() with the FIXED-effect design matrix
  # (random-effect term skipped — single-cohort or near-singular case).
  form_fixed_str <- sub("\\s*\\+\\s*\\(1\\s*\\|\\s*dataset\\)", "", form_str)
  if (form_fixed_str == form_str) {
    form_fixed_str <- paste0(form_str, " + dataset")  # absorb cohort fixed
  }
  cat("  Fallback formula (fixed-effects): ", form_fixed_str, "\n")
  design_mm <- model.matrix(as.formula(form_fixed_str), info)
  v_voom <- tryCatch(voom(dge, design = design_mm, span = 0.5),
                     error = function(e) {
                       cat("  voom() also FAILED (", conditionMessage(e), ")\n",
                           sep = "")
                       NULL
                     })
  if (is.null(v_voom)) {
    # Last resort: log-CPM transformation, no precision weights.
    cat("  Falling back to log2-CPM transform (no precision weights)\n")
    v_voom <- list(E = edgeR::cpm(dge, log = TRUE, prior.count = 2),
                   weights = matrix(1, nrow(dge), ncol(dge)),
                   design = design_mm,
                   targets = dge$samples)
    class(v_voom) <- "EList"
  }
  cat("  voom-fallback elapsed:", format(Sys.time() - t0), "\n")

  # When fallback was used, dream's random-effect path may also fail; refit
  # with plain lmFit on the fixed-effect design.
  cat("\n[1b] Fitting via limma::lmFit (fallback path, no random effect)...\n")
  fit <- eBayes(lmFit(v_voom, design = design_mm))
} else {
  cat("  voom elapsed:", format(Sys.time() - t0), "\n")
  t0 <- Sys.time()
  fit <- suppressWarnings(
    dream(v_voom, form, info, BPPARAM = param, useWeights = TRUE))
  cat("  dream elapsed:", format(Sys.time() - t0), "\n")
}

# ---------------------------------------------------------------------------
# Extract interaction coefficient (Wald)
# ---------------------------------------------------------------------------
all_coefs <- colnames(fit$coefficients)
gcoef <- setdiff(grep("^group_binary", all_coefs, value = TRUE),
                 grep(":", all_coefs, value = TRUE))[1]
icoef <- grep(":", grep("group_binary", all_coefs, value = TRUE), value = TRUE)[1]
cat("\ngroup-coef:", gcoef, "  interaction-coef:", icoef, "\n")

beta_int <- fit$coefficients[, icoef]
raw_se   <- fit$stdev.unscaled * fit$sigma
se_int   <- raw_se[, icoef]
t_int    <- beta_int / pmax(se_int, .Machine$double.eps)
p_int    <- 2 * pnorm(-abs(t_int))

# Tier-1 universe filter for the conditional BH (matches Disease-vs-Ctrl).
# Reads the per-contrast disease-effect CSV per CSPEC$tier1_csv.
tier1_path <- file.path(BASE, CSPEC$tier1_csv)
if (!file.exists(tier1_path)) {
  stop("Tier-1 CSV not found at ", tier1_path)
}
tier1 <- fread(tier1_path)
# Standard column names: logFC + adj.P.Val
tier1_gene_col <- intersect(c("gene", "Gene", "gene_id"), names(tier1))[1]
if (is.na(tier1_gene_col)) {
  # Last resort: use rownames-like first column
  tier1_gene_col <- names(tier1)[1]
}
tier1_genes <- tier1[!is.na(get(CSPEC$tier1_padj_col)) &
                       get(CSPEC$tier1_padj_col) < 0.05 &
                       abs(get(CSPEC$tier1_logFC_col)) > 0.3,
                     get(tier1_gene_col)]
cat("Tier-1 universe (padj<0.05 & |logFC|>0.3) from",
    basename(tier1_path), ":", length(tier1_genes), "genes\n")

padj_int_5k <- rep(NA_real_, length(beta_int))
in_tier1 <- names(beta_int) %in% tier1_genes
padj_int_5k[in_tier1] <- p.adjust(p_int[in_tier1], method = "BH")
padj_int_full <- p.adjust(p_int, method = "BH")

cat("padj_int_5k < 0.05:", sum(padj_int_5k < 0.05, na.rm = TRUE),
    "  < 0.10:", sum(padj_int_5k < 0.10, na.rm = TRUE),
    "  < 0.20:", sum(padj_int_5k < 0.20, na.rm = TRUE), "\n")

# ---------------------------------------------------------------------------
# Per-cohort fixed-effects Cochran's Q (heterogeneity)
# ---------------------------------------------------------------------------
cat("\n[2] Per-cohort fixed-effects β_int for Cochran's Q...\n")
percoh <- vector("list", length(cohort_identifiable))
names(percoh) <- cohort_identifiable

for (ch in cohort_identifiable) {
  keep_c <- info$dataset == ch
  info_c <- info[keep_c, , drop = FALSE]
  dge_c <- dge[, keep_c]
  dge_c$samples$norm.factors <- calcNormFactors(dge_c)$samples$norm.factors

  # Drop SVs with zero within-cohort variance
  sv_keep <- colnames(sv_df)[vapply(colnames(sv_df),
                                     function(s) var(info_c[[s]], na.rm = TRUE) > 1e-10,
                                     logical(1))]
  sv_terms_c <- if (length(sv_keep) > 0) paste("+", paste(sv_keep, collapse = " + ")) else ""
  form_c <- as.formula(paste0(
    "~ group_binary * inferred_sex",
    " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
    " + age_imputed", sv_terms_c))
  mm <- tryCatch(model.matrix(form_c, info_c), error = function(e) NULL)
  if (is.null(mm)) next
  v_c <- tryCatch(voom(dge_c, mm), error = function(e) NULL)
  if (is.null(v_c)) next
  fit_c <- tryCatch(eBayes(lmFit(v_c, mm)), error = function(e) NULL)
  if (is.null(fit_c)) next

  ic <- grep(":", grep("group_binary", colnames(fit_c$coefficients), value = TRUE), value = TRUE)[1]
  if (is.na(ic)) next
  raw_se_c <- fit_c$stdev.unscaled * fit_c$sigma
  percoh[[ch]] <- data.table(
    gene = rownames(fit_c$coefficients),
    cohort = ch,
    beta_int_c = fit_c$coefficients[, ic],
    se_int_c   = raw_se_c[, ic]
  )
  cat("  ", ch, ": ", sum(!is.na(percoh[[ch]]$beta_int_c)), " estimable genes\n", sep = "")
}

percoh_dt <- rbindlist(percoh, fill = TRUE, use.names = TRUE)
percoh_dt <- percoh_dt[!is.na(beta_int_c) & !is.na(se_int_c) & se_int_c > 0]
percoh_dt[, w := 1 / (se_int_c^2)]

pooled <- percoh_dt[, {
  if (.N < 2) {
    .(beta_pool = NA_real_, Q = NA_real_, df = NA_integer_,
      p_Q = NA_real_, I2_pct = NA_real_, k = .N)
  } else {
    wsum <- sum(w)
    bp   <- sum(w * beta_int_c) / wsum
    Q    <- sum(w * (beta_int_c - bp)^2)
    df   <- .N - 1L
    pQ   <- pchisq(Q, df, lower.tail = FALSE)
    I2   <- max(0, (Q - df) / Q) * 100
    .(beta_pool = bp, Q = Q, df = df, p_Q = pQ, I2_pct = I2, k = .N)
  }
}, by = gene]
pooled[, Q_padj_BH := p.adjust(p_Q, "BH")]

# ---------------------------------------------------------------------------
# Final per-gene table
# ---------------------------------------------------------------------------
out <- data.table(
  gene          = names(beta_int),
  beta_int      = beta_int,
  se_int        = se_int,
  t_int         = t_int,
  p_int         = p_int,
  padj_int_5k   = padj_int_5k,
  padj_int_full = padj_int_full,
  in_tier1      = in_tier1
)
out <- merge(out, pooled[, .(gene, cohort_Q_pval = p_Q, I2_pct,
                              n_cohorts_identifiable = k)],
             by = "gene", all.x = TRUE)
# Genes that weren't in the per-cohort merge have n_cohorts_identifiable = NA → set to 0
out[is.na(n_cohorts_identifiable), n_cohorts_identifiable := 0L]

OUT_CSV <- file.path(SEXV3, "interaction_fixed.csv")
write_atomic_csv(out, OUT_CSV)
cat("\nWrote: ", OUT_CSV, "  rows:", nrow(out), "  cols:", ncol(out), "\n")
cat("Identifiable cohorts in pooled set:", length(unique(percoh_dt$cohort)), "\n")
cat("Genes with |β_int|>0.5 & padj_int_5k<0.20:",
    sum(abs(out$beta_int) > 0.5 & out$padj_int_5k < 0.20, na.rm = TRUE), "\n")

dump_session_info(IDIR, "19")
cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
