#!/usr/bin/env Rscript
# power_estimate_nb_params.R
# ===========================================================================
# Multi-method DE validation harness — Pillar C (POWER): STEP 1 of 2.
#
# Estimate realistic per-gene negative-binomial simulation parameters from the
# REAL mega-analysis CONTROL samples, then cache them. This is run ONCE so the
# 48-task power array never races on (or redundantly recomputes) the parameter
# estimation. The array (power_multimethod.R) only READS the cache.
#
# WHY CONTROLS ONLY: the simulation injects a *known* disease effect on top of a
# disease-free baseline. Estimating mu/phi from controls gives the null baseline
# (mean expression + biological dispersion) without contaminating it with the
# real disease signal we are about to simulate.
#
# NB PARAMETERIZATION (must match power_multimethod.R rnbinom calls):
#   edgeR models Var = mu + phi * mu^2 (phi = tagwise.dispersion).
#   R's rnbinom(mu=, size=) has   Var = mu + mu^2 / size.
#   => size = 1/phi. The simulator uses rnbinom(mu = mean_gk, size = 1/phi_g).
#   mu_g (expected counts at reference lib size L) is recovered from edgeR's
#   aveLogCPM:  mu_g = (2^aveLogCPM_g / 1e6) * L,  L = median library size.
#
# Env:  MASLD_PROJECT_ROOT, SLURM_CPUS_PER_TASK (only for edgeR threading hints)
# Out:  results/integration/multimethod_validation/power/nb_params.rds
#         list(mu = named numeric, phi = named numeric, L = numeric,
#              genes = character, n_control_samples, n_cohorts, est_design)
# ===========================================================================

t0 <- proc.time()

# --- 1. Source the shared (TESTED) helpers; do not edit them ----------------
PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HELPERS <- file.path(PROJECT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts",
  "multimethod_validation/de_validation_helpers.R")
source(HELPERS)

cat("=== power_estimate_nb_params: NB parameter estimation from controls ===\n")

out_dir  <- file.path(RDIR, "multimethod_validation/power")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "nb_params.rds")

# --- 2. Load mega data; restrict to Control samples across the 5 cohorts -----
d      <- load_mega_data()
counts <- d$counts
meta   <- d$meta

is_ctrl <- meta$group_binary == "Control"
if (sum(is_ctrl) < 10) stop("Too few Control samples (", sum(is_ctrl), ") to estimate NB params")
counts <- counts[, is_ctrl, drop = FALSE]
meta   <- meta[is_ctrl]
stopifnot(all(meta$sample_id == colnames(counts)))
cat("Control samples:", ncol(counts), "across",
    length(unique(meta$dataset)), "cohorts\n")
print(table(meta$dataset))

# --- 3. edgeR: filterByExpr + TMM + estimateDisp (design ~ dataset) ----------
# Use ~dataset so cohort baseline differences in the controls do not inflate the
# tagwise dispersion (they would be soaked up as residual overdispersion under
# an intercept-only model).
dge <- DGEList(counts = counts)
dge$samples$dataset <- droplevels(factor(meta$dataset))

design <- model.matrix(~ dataset, data = dge$samples)
keep   <- filterByExpr(dge, design = design)
dge    <- dge[keep, , keep.lib.sizes = FALSE]
cat("Genes after filterByExpr(design = ~dataset):", nrow(dge), "\n")

dge <- calcNormFactors(dge, method = "TMM")
dge <- estimateDisp(dge, design = design)

# --- 4. Extract per-gene NB parameters --------------------------------------
# phi_g: tagwise dispersion (fall back to trended/common if any are NA).
phi <- dge$tagwise.dispersion
if (is.null(phi)) phi <- dge$trended.dispersion
if (is.null(phi)) phi <- rep(dge$common.dispersion, nrow(dge))
phi[!is.finite(phi) | phi <= 0] <- dge$common.dispersion
names(phi) <- rownames(dge)

# mu_g: expected counts at reference library size L (= median TMM-effective lib
# size of the control samples). aveLogCPM is on log2 CPM scale.
L  <- as.numeric(stats::median(dge$samples$lib.size * dge$samples$norm.factors))
av <- edgeR::aveLogCPM(dge)                 # log2 CPM per gene
mu <- (2^av / 1e6) * L                      # expected counts at lib size L
mu[!is.finite(mu) | mu < 0] <- 0
names(mu) <- rownames(dge)

genes <- rownames(dge)
stopifnot(identical(names(mu), genes), identical(names(phi), genes))

# --- 5. Save + summary ------------------------------------------------------
params <- list(
  mu               = mu,
  phi              = phi,
  L                = L,
  genes            = genes,
  n_control_samples = ncol(counts),
  n_cohorts        = nlevels(dge$samples$dataset),
  est_design       = "~ dataset",
  common_dispersion = dge$common.dispersion
)
saveRDS(params, out_file)

cat("\n--- NB parameter summary ---\n")
cat("n genes (filtered) :", length(genes), "\n")
cat("ref lib size L     :", format(round(L), big.mark = ","), "\n")
cat("mu  quantiles      :",
    paste(round(stats::quantile(mu,  c(0, .25, .5, .75, 1)), 2), collapse = "  "), "\n")
cat("phi quantiles      :",
    paste(round(stats::quantile(phi, c(0, .25, .5, .75, 1)), 4), collapse = "  "), "\n")
cat("common dispersion  :", round(dge$common.dispersion, 4), "\n")
cat("Saved ->", out_file, "\n")
cat(sprintf("Elapsed %.1f min\n", (proc.time() - t0)["elapsed"] / 60))
