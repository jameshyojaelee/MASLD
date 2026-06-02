#!/usr/bin/env Rscript
# sex_v3/02_sva_factors.R
# ---------------------------------------------------------------------------
# Module 02 — Estimate hidden Surrogate Variables (SVA) from null-model
# residuals of voomWithDreamWeights, conditioning on sex + dataset +
# composition.
#
# Strategy:
#   - Fit voomWithDreamWeights with NULL model
#         ~ inferred_sex + dataset + Hep + Mac + Endo + Chol
#     (no group_binary — SVs should capture residual structure independent
#      of biology of interest).
#   - Convert voom log2-CPM to a matrix and run sva::sva() with full vs null
#     models so SVs are independent of biological signal (group_binary
#     enters via the full model).
#   - num.sv(method = "leek") for k; cap at 5.
#
# Sanity gate:
#   - No SV should correlate with sex/disease at |r| > 0.3
#   - Save factor_loadings_qc.pdf with correlation heatmap
#
# Reads:  intermediates/sex_v3_input.rds
# Writes: intermediates/sva_factors.rds
#         factor_loadings_qc.pdf
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
  library(sva)
})
# Force injection into lme4 namespace BEFORE loading variancePartition
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({
      unlockBinding(fn, ns_lme4)
      assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
      lockBinding(fn, ns_lme4)
    }, silent = TRUE)
  }
}
suppressPackageStartupMessages({
  library(variancePartition)
  library(BiocParallel)
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RDIR   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                    "results/integration")
SEXV3  <- file.path(RDIR, "sex_v3")
IDIR   <- file.path(SEXV3, "intermediates")

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
dir.create(IDIR, recursive = TRUE, showWarnings = FALSE)

# Shared utilities (atomic writes + sessionInfo dump) -- R5 Issue 2 fix
source(file.path(BASE,
                 "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                 "scripts/sex_v3/sex_v3_utils.R"))

cat("============================================================\n")
cat("Module 02 — SVA factor estimation\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")

# ---------------------------------------------------------------------------
# Parallel setup
# ---------------------------------------------------------------------------
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus, RNGseed = 42L) else SerialParam()

# ---------------------------------------------------------------------------
# Load input bundle
# ---------------------------------------------------------------------------
in_path <- file.path(IDIR, "sex_v3_input.rds")
cat("Loading:", in_path, "\n")
inp <- readRDS(in_path)
dge  <- inp$dge
info <- inp$meta
cat("  Samples:", ncol(dge), "  Genes:", nrow(dge), "\n")

# ---------------------------------------------------------------------------
# voom log2-CPM for SVA input
# ---------------------------------------------------------------------------
# R1 Issue 3 / 7 fix: BOTH null and full models must include the group_binary *
# inferred_sex interaction so SVA cannot absorb interaction-driven variance
# (= the biology Module 04 is going to model). Old design had only main effects
# in the null model — an SV that captured the disease-by-sex interaction would
# pass the |r|>0.3 marginal-correlation gate and silently regress out sex-
# specific disease signal.
#
# Concretely:
#   mod_full = full biological model (interaction + composition + dataset)
#   mod_null = same MINUS group_binary main + interaction (so SVA estimates SVs
#              that are orthogonal to the disease/interaction direction).
biology_form <- ~ group_binary * inferred_sex + dataset + Hepatocytes +
                  Macrophages + Endothelial + Cholangiocytes
null_form    <- ~                 inferred_sex + dataset + Hepatocytes +
                  Macrophages + Endothelial + Cholangiocytes

cat("\n[1] voomWithDreamWeights with null formula (sex + dataset + composition)...\n")
cat("   null_form: ", deparse(null_form), "\n")
cat("   full_form: ", deparse(biology_form), "\n")
v_null <- suppressWarnings(
  voomWithDreamWeights(dge, null_form, info, BPPARAM = param)
)
expr <- v_null$E   # log2-CPM
stopifnot(ncol(expr) == nrow(info))
cat("  voom expression matrix:", nrow(expr), "x", ncol(expr), "\n")

# ---------------------------------------------------------------------------
# Build full and null model matrices for sva()
# ---------------------------------------------------------------------------
mod_full <- model.matrix(biology_form, data = info)   # includes interaction
mod_null <- model.matrix(null_form,    data = info)   # excludes group_binary + interaction
# Alias for backward compatibility with the saveRDS list at the bottom
full_form <- biology_form
cat("  mod_full:", paste(dim(mod_full), collapse = " x "),
    "  mod_null:", paste(dim(mod_null), collapse = " x "), "\n")

# ---------------------------------------------------------------------------
# Estimate n.sv via leek method, cap at 5
# ---------------------------------------------------------------------------
cat("\n[2] Estimating n.sv via num.sv(method = 'leek')...\n")
n_sv_leek <- tryCatch(
  num.sv(expr, mod_full, method = "leek"),
  error = function(e) {
    cat("  num.sv(leek) failed:", conditionMessage(e), "\n")
    NA_integer_
  }
)
n_sv_be <- tryCatch(
  num.sv(expr, mod_full, method = "be"),
  error = function(e) NA_integer_
)
cat("  num.sv(leek) =", n_sv_leek, "  num.sv(be) =", n_sv_be, "\n")

# Cap at 5 (keep model identifiable); use leek by default but fall back if NA
chosen_n <- if (!is.na(n_sv_leek) && n_sv_leek > 0) n_sv_leek else
            if (!is.na(n_sv_be)   && n_sv_be   > 0) n_sv_be   else 1L
chosen_n <- min(as.integer(chosen_n), 5L)
chosen_n <- max(chosen_n, 1L)
cat("  Capped n_sv =", chosen_n, "\n")

# ---------------------------------------------------------------------------
# Fit SVA
# ---------------------------------------------------------------------------
cat("\n[3] Fitting sva()...\n")
sv_obj <- tryCatch(
  sva(expr, mod_full, mod_null, n.sv = chosen_n),
  error = function(e) {
    cat("  sva() failed at n.sv =", chosen_n, ":", conditionMessage(e), "\n")
    fallback_n <- max(1L, chosen_n %/% 2L)
    cat("  Falling back to n.sv =", fallback_n, "\n")
    sva(expr, mod_full, mod_null, n.sv = fallback_n)
  }
)
final_n <- ncol(sv_obj$sv)
cat("  Final n.sv =", final_n, "\n")

sv_mat <- sv_obj$sv
colnames(sv_mat) <- paste0("SV", seq_len(final_n))
rownames(sv_mat) <- rownames(info)

# ---------------------------------------------------------------------------
# Sanity: correlate each SV with sex, disease, dataset, composition, age
# ---------------------------------------------------------------------------
cat("\n[4] SV correlation with covariates (Pearson |r|):\n")

# Build a numeric covariate matrix for correlation testing.
# R1 Issue 3 / 7 fix: also test against the 4-level group_binary * inferred_sex
# interaction (Disease & Male combination) — an SV absorbing the interaction
# signal would slip past the old marginal-only sex/disease gate.
sex_is_male  <- as.integer(info$inferred_sex == "M" | info$inferred_sex == "Male")
disease_int  <- as.integer(info$group_binary == "Disease")
covar_df <- data.frame(
  group_disease  = disease_int,
  sex_female     = as.integer(info$inferred_sex == "F" | info$inferred_sex == "Female"),
  group_x_sex    = disease_int * sex_is_male,   # interaction indicator (Disease-Male cell)
  Hepatocytes    = info$Hepatocytes,
  Macrophages    = info$Macrophages,
  Endothelial    = info$Endothelial,
  Cholangiocytes = info$Cholangiocytes
)
# Encode dataset by one-hot (max abs r per SV across dataset dummies)
ds_oh <- model.matrix(~ 0 + dataset, data = info)

cor_mat <- matrix(NA_real_, nrow = final_n,
                  ncol = ncol(covar_df),
                  dimnames = list(colnames(sv_mat), colnames(covar_df)))
for (j in seq_len(ncol(covar_df))) {
  for (i in seq_len(final_n)) {
    cor_mat[i, j] <- suppressWarnings(cor(sv_mat[, i], covar_df[[j]],
                                          use = "pairwise.complete.obs"))
  }
}
# Per-SV max |r| against dataset one-hot
ds_maxcor <- sapply(seq_len(final_n), function(i) {
  max(abs(suppressWarnings(cor(sv_mat[, i], ds_oh,
                               use = "pairwise.complete.obs"))), na.rm = TRUE)
})

cor_full <- cbind(cor_mat, dataset_maxabs = ds_maxcor)
cat("\nCorrelation matrix (rows = SVs):\n")
print(round(cor_full, 3))

# Sanity flags (R1 Issue 3 / 7 fix: gate also on the interaction)
flag_sex     <- abs(cor_mat[, "sex_female"])    > 0.3
flag_disease <- abs(cor_mat[, "group_disease"]) > 0.3
flag_interx  <- abs(cor_mat[, "group_x_sex"])   > 0.3
cat("\nSanity gate (|r| > 0.3 vs sex, disease, OR group_x_sex interaction):\n")
if (any(flag_sex)) {
  cat("  ! SV(s) correlated with sex (problematic):",
      paste(rownames(cor_mat)[flag_sex], collapse = ","), "\n")
}
if (any(flag_disease)) {
  cat("  ! SV(s) correlated with disease (problematic):",
      paste(rownames(cor_mat)[flag_disease], collapse = ","), "\n")
}
if (any(flag_interx)) {
  cat("  ! SV(s) correlated with group_x_sex INTERACTION (problematic — would mask sex-specific disease effects):",
      paste(rownames(cor_mat)[flag_interx], collapse = ","), "\n")
}
if (!any(flag_sex) && !any(flag_disease) && !any(flag_interx)) {
  cat("  OK — no SV correlates with sex, disease, or group_x_sex at |r| > 0.3\n")
}

# ---------------------------------------------------------------------------
# QC PDF — correlation heatmap
# ---------------------------------------------------------------------------
cat("\n[5] Writing QC PDF...\n")
cor_long <- reshape2::melt(cor_full)
colnames(cor_long) <- c("SV", "covariate", "r")
p <- ggplot(cor_long, aes(x = covariate, y = SV, fill = r)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
  scale_fill_gradient2(low = "#1F77B4", mid = "white", high = "#D62728",
                       midpoint = 0, limits = c(-1, 1)) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "SV factor loadings (Pearson r with covariates)",
       subtitle = "Sanity gate: no SV should hit |r| > 0.3 with sex, disease, or group_x_sex",
       x = NULL, y = NULL)

ggsave(file.path(SEXV3, "factor_loadings_qc.pdf"), p,
       width = 8, height = max(4, 1.0 + 0.4 * final_n))
cat("  Saved:", file.path(SEXV3, "factor_loadings_qc.pdf"), "\n")

# ---------------------------------------------------------------------------
# Save SV matrix
# ---------------------------------------------------------------------------
out_path <- file.path(IDIR, "sva_factors.rds")
write_atomic_rds(list(
  sv          = sv_mat,
  n_sv        = final_n,
  n_sv_leek   = n_sv_leek,
  n_sv_be     = n_sv_be,
  cor_matrix  = cor_full,
  flagged     = list(sex      = which(flag_sex),
                     disease  = which(flag_disease),
                     interx   = which(flag_interx)),
  full_form   = full_form,
  null_form   = null_form,
  pprob.gam   = sv_obj$pprob.gam,
  pprob.b     = sv_obj$pprob.b
), out_path)
cat("Saved (atomic):", out_path, "\n")

# R5 Issue 2 fix: sessionInfo dump for per-module reproducibility audit
dump_session_info(IDIR, "02")

cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
