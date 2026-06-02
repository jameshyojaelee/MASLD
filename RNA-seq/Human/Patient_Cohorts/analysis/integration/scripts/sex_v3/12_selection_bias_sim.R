#!/usr/bin/env Rscript
# sex_v3/12_selection_bias_sim.R
# ---------------------------------------------------------------------------
# Pillar 3 — Conditional-universe selection-bias simulation.
#
# Forks 07a2_v5_synthetic.R. Per cell of the (delta_main x delta_int x seed)
# grid:
#   * inject 500 sex-modulated genes whose disease-main effect is drawn
#     INDEPENDENTLY of the interaction effect (delta_main, delta_int);
#   * inject 1000 carrier-null genes (no signal);
#   * refit dream M2 (current v5 production formula F5
#     = `(1 + group_binary | dataset)`);
#   * apply the canonical Tier-1 filter `padj_main < 0.05 AND |logFC|>0.3`
#     using the dream M2 main-effect padj computed JUST on the synthetic
#     genes (so Tier-1 BH stays internally consistent for the simulation);
#   * compute retention_rate(|β_main| decile, δ_int) and write per-cell rds.
#
# Key design distinction vs the calibration sim (07a2_v5_synthetic.R):
#   * β_main and β_int drawn INDEPENDENTLY -> the only thing that can make
#     retention rate depend on |β_main| is the Tier-1 main-effect padj cut.
#   * If retention curves are flat across |β_main| deciles, no selection
#     bias. If they slope downward at low |β_main|, the conditional-universe
#     approach depletes interaction signal at low main-effect magnitudes.
#
# CLI knobs (env vars): TID -> (delta_main_idx, delta_int_idx, seed). Decoder:
#   TID in [1..120]
#   delta_main_idx = ((TID-1) %/% 20) + 1    # in [1..6]
#   delta_int_idx  = (((TID-1) %% 20) %/% 5) + 1   # in [1..4]
#   seed           = ((TID-1) %% 5) + 1            # in [1..5]
#   delta_main = c(0.0, 0.3, 0.5, 0.8, 1.2, 2.0)[delta_main_idx]
#   delta_int  = c(0.0, 0.3, 0.5, 0.8)[delta_int_idx]
#
# Output: intermediates/selection_bias_v6/cell_<TID>.rds
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

TID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID",
                              Sys.getenv("TID", "1")))
if (is.na(TID) || TID < 1L || TID > 120L) {
  stop("TID env var (or SLURM_ARRAY_TASK_ID) must be in [1, 120]; got ", TID)
}

DM_GRID <- c(0.0, 0.3, 0.5, 0.8, 1.2, 2.0)
DI_GRID <- c(0.0, 0.3, 0.5, 0.8)
dm_idx  <- ((TID - 1L) %/% 20L) + 1L
di_idx  <- (((TID - 1L) %% 20L) %/% 5L) + 1L
SEED    <- ((TID - 1L) %% 5L) + 1L
DELTA_MAIN <- DM_GRID[dm_idx]
DELTA_INT  <- DI_GRID[di_idx]

N_SEX_MOD <- as.integer(Sys.getenv("N_SEX_MOD", "500"))
N_CARRIER <- as.integer(Sys.getenv("N_CARRIER", "1000"))
set.seed(SEED * 1009L + TID * 13L)

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
OUTDIR <- file.path(IDIR, "selection_bias_v6")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
out_rds <- file.path(OUTDIR, sprintf("cell_%d.rds", TID))

cat("=== 12_selection_bias_sim ===\n")
cat("TID:", TID, "  SEED:", SEED,
    "  DELTA_MAIN:", DELTA_MAIN, "  DELTA_INT:", DELTA_INT, "\n")
cat("N_SEX_MOD:", N_SEX_MOD, "  N_CARRIER:", N_CARRIER, "\n")

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
param <- if (ncpus > 1) MulticoreParam(workers = ncpus, RNGseed = SEED) else SerialParam()

# Load real data (same scaffolding as 07a2_v5_synthetic.R)
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

# Estimate per-gene NB dispersion
cat("Estimating per-gene NB dispersion...\n")
group_for_disp <- paste(info$inferred_sex, info$group_binary, sep = ":")
design_disp <- model.matrix(~ group_for_disp)
dge <- estimateDisp(dge, design_disp, robust = TRUE)
phi <- dge$tagwise.dispersion
real_genes <- rownames(dge$counts)

# Pick donor genes for sex-modulated + carrier-null sets
set.seed(SEED * 7919L + TID * 17L)
all_idx <- sample.int(length(real_genes),
                      size = N_SEX_MOD + N_CARRIER, replace = FALSE)
sex_mod_genes <- real_genes[all_idx[seq_len(N_SEX_MOD)]]
carrier_genes  <- real_genes[all_idx[(N_SEX_MOD + 1L):(N_SEX_MOD + N_CARRIER)]]

# Truth table: per sex-modulated gene, draw β_main and β_int INDEPENDENTLY.
# This is the critical knob. β_main lives in N(0, DELTA_MAIN^2 / 2) so that
# |β_main| spans 0..~2*DELTA_MAIN; |β_int| similarly with DELTA_INT.
# (carrier-null genes get β_main = β_int = 0.)
draw_betas <- function(n, dm, di) {
  # symmetric normal with sd scaled so RMS = dm/sqrt(2) for non-zero grids
  # When dm == 0 -> exact zeros.
  bm <- if (dm > 0) rnorm(n, 0, dm) else rep(0, n)
  bi <- if (di > 0) rnorm(n, 0, di) else rep(0, n)
  cbind(beta_main_true = bm, beta_int_true = bi)
}
truth_sex <- data.table(
  synth_id = paste0("synthSM_", seq_len(N_SEX_MOD)),
  donor_gene = sex_mod_genes,
  beta_main_true = draw_betas(N_SEX_MOD, DELTA_MAIN, DELTA_INT)[, 1],
  beta_int_true  = draw_betas(N_SEX_MOD, DELTA_MAIN, DELTA_INT)[, 2]
)
# Replay to make sure both columns are independent draws.
# (The above call drew twice and took col 1 once, col 2 the second time.)
# Make this explicit:
truth_sex[, beta_main_true := rnorm(.N, 0, max(DELTA_MAIN, 1e-9))]
truth_sex[, beta_int_true  := rnorm(.N, 0, max(DELTA_INT,  1e-9))]
if (DELTA_MAIN == 0) truth_sex[, beta_main_true := 0]
if (DELTA_INT  == 0) truth_sex[, beta_int_true  := 0]
truth_sex[, pattern := "sex_modulated"]

truth_car <- data.table(
  synth_id = paste0("synthCAR_", seq_len(N_CARRIER)),
  donor_gene = carrier_genes,
  beta_main_true = 0, beta_int_true = 0, pattern = "carrier_null"
)
truth <- rbind(truth_sex, truth_car)

cat("Injecting", nrow(truth), "synthetic genes at count scale...\n")
real_cnt  <- dge$counts
synth_cnt <- matrix(0L, nrow = nrow(truth), ncol = n_samp,
                    dimnames = list(truth$synth_id, colnames(real_cnt)))

is_F   <- info$inferred_sex == "F"
is_M   <- info$inferred_sex == "M"
is_dis <- info$group_binary == "Disease"

for (i in seq_len(nrow(truth))) {
  gd <- truth$donor_gene[i]
  mu_real <- pmax(real_cnt[gd, ], 0.5)
  # β_F = beta_main_true; β_M = beta_main_true + beta_int_true
  # Disease effect lives in is_dis; sex-specific lift via inferred_sex.
  log2fc <- rep(0, n_samp)
  bF <- truth$beta_main_true[i]
  bM <- truth$beta_main_true[i] + truth$beta_int_true[i]
  log2fc[is_F & is_dis] <- bF
  log2fc[is_M & is_dis] <- bM
  mu_inj <- mu_real * (2 ^ log2fc)
  phi_g <- phi[match(gd, real_genes)]
  if (!is.finite(phi_g) || phi_g <= 0) phi_g <- 0.1
  synth_cnt[i, ] <- as.integer(rnbinom(n_samp, mu = mu_inj, size = 1 / phi_g))
}

# Build fresh DGEList with synth-only matrix (simulation universe).
dge_synth <- DGEList(synth_cnt, samples = dge$samples)
dge_synth$samples$norm.factors <- dge$samples$norm.factors
dge_synth <- calcNormFactors(dge_synth)

# v6 primary formula F2 (diagonal random slope on sex × disease). Earlier
# P3 used F5 (single random slope on group_binary), which is the v5 fallback,
# not v6's primary. Meta-review A10 noted Spearman(F5, F2)=0.832 is NOT a
# license — at the boundary cells F2 produces 0/50 Branch-B fires where F5
# produces thousands. Re-run under F2 so the Bourgon-condition claim binds
# the right model.
sv_terms <- paste(colnames(sv_df), collapse = " + ")
form_str_F2 <- paste0(
  "~ group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed + ", sv_terms,
  " + (1 + group_binary + inferred_sex + group_binary:inferred_sex || dataset)"
)
form_str_F5 <- paste0(
  "~ group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed + ", sv_terms,
  " + (1 + group_binary | dataset)"
)
form_str_F6 <- paste0(
  "~ group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed + ", sv_terms,
  " + (1 | dataset)"
)
form_slope     <- as.formula(form_str_F2)
form_F5        <- as.formula(form_str_F5)
form_intercept <- as.formula(form_str_F6)

cat("voomWithDreamWeights + dream (F2 → F5 → F6 fallback)...\n")
t0 <- Sys.time()
v <- tryCatch(
  suppressWarnings(voomWithDreamWeights(dge_synth, form_slope, info,
                                        BPPARAM = param, useWeights = TRUE)),
  error = function(e) {
    cat("  voom F2 fail; trying F5:", conditionMessage(e), "\n")
    tryCatch(
      suppressWarnings(voomWithDreamWeights(dge_synth, form_F5, info,
                                            BPPARAM = param, useWeights = TRUE)),
      error = function(e2) {
        cat("  voom F5 fail; falling back to intercept:", conditionMessage(e2), "\n")
        suppressWarnings(voomWithDreamWeights(dge_synth, form_intercept, info,
                                              BPPARAM = param, useWeights = TRUE))
      })
  })
fit <- tryCatch(
  suppressWarnings(dream(v, form_slope, info, BPPARAM = param, useWeights = TRUE)),
  error = function(e) {
    cat("  dream F2 fail; trying F5:", conditionMessage(e), "\n")
    tryCatch(
      suppressWarnings(dream(v, form_F5, info, BPPARAM = param, useWeights = TRUE)),
      error = function(e2) {
        cat("  dream F5 fail; falling back to intercept:", conditionMessage(e2), "\n")
        suppressWarnings(dream(v, form_intercept, info,
                               BPPARAM = param, useWeights = TRUE))
      })
  })
cat("  dream elapsed:", format(Sys.time() - t0), "\n")

# Extract coefficients
all_coefs <- colnames(fit$coefficients)
gcoef <- setdiff(grep("^group_binary", all_coefs, value = TRUE),
                 grep(":", all_coefs, value = TRUE))[1]
icoef <- grep(":", grep("^group_binary", all_coefs, value = TRUE), value = TRUE)[1]

bF <- fit$coefficients[, gcoef]
bI <- fit$coefficients[, icoef]
bM <- bF + bI
raw_se <- fit$stdev.unscaled * fit$sigma
sF <- raw_se[, gcoef]; sI <- raw_se[, icoef]
# Per-gene Cov(β_F, β_int) — variancePartition MArrayLM2 stores per-gene
# vcov in cov.coefficients.list; scalar cov.coefficients is the limma API
# and returns NULL for dream LMMs, silently zeroing 2*Cov.
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

# Main-effect Wald p-value (treat β_F as main effect proxy for Tier-1 filter,
# matching the production "group_binaryDisease" coefficient identity)
t_main <- bF / sF
p_main <- 2 * pnorm(-abs(t_main))
padj_main <- p.adjust(p_main, method = "BH")
t_int <- bI / sI
p_int <- 2 * pnorm(-abs(t_int))

# Tier-1 filter (production rule, exactly as 07_interaction_test_v5.R lines 42-44)
# applied to the synthetic universe.
all_genes <- rownames(fit$coefficients)
tier1 <- (!is.na(padj_main) & padj_main < 0.05 & abs(bF) > 0.3)
names(tier1) <- all_genes

# Per-gene retention table
gene_tbl <- data.table(
  gene = all_genes,
  beta_F_est = bF, se_F = sF, p_main = p_main, padj_main = padj_main,
  beta_int_est = bI, se_int = sI, t_int = t_int, p_int = p_int,
  tier1_retained = tier1
)
gene_tbl <- merge(gene_tbl,
                  truth[, .(gene = synth_id, pattern,
                            beta_main_true, beta_int_true)],
                  by = "gene", all.x = TRUE)
gene_tbl[is.na(pattern), pattern := "unknown"]

# Retention decile table (sex-modulated only)
sm <- gene_tbl[pattern == "sex_modulated"]
sm[, abs_bmain := abs(beta_main_true)]
# breaks: 10 quantile deciles of |β_main_true|. If DELTA_MAIN==0, all
# true main-effects are 0 -> single bin.
if (DELTA_MAIN == 0 || length(unique(sm$abs_bmain)) < 2) {
  sm[, decile := 1L]
  cuts <- c(0, 1e-9)
} else {
  cuts <- unique(quantile(sm$abs_bmain, probs = seq(0, 1, by = 0.1),
                          na.rm = TRUE))
  sm[, decile := as.integer(cut(abs_bmain, breaks = cuts,
                                include.lowest = TRUE, labels = FALSE))]
}

ret <- sm[, .(
  n_genes        = .N,
  n_retained     = sum(tier1_retained, na.rm = TRUE),
  retention_rate = sum(tier1_retained, na.rm = TRUE) / .N,
  mean_abs_bmain = mean(abs_bmain, na.rm = TRUE),
  mean_abs_bint  = mean(abs(beta_int_true), na.rm = TRUE)
), by = decile][order(decile)]
ret[, TID := TID][, delta_main := DELTA_MAIN][, delta_int := DELTA_INT][, seed := SEED]

cat("\n== retention by |β_main| decile ==\n")
print(ret)

# Tier-1-vs-genome-wide power side-by-side
# Genome-wide BH (on synthetic universe = sex_modulated + carrier_null)
padj_int_gw <- p.adjust(p_int, method = "BH")
# Tier-1-restricted BH
padj_int_t1 <- rep(NA_real_, length(all_genes))
if (any(tier1)) {
  padj_int_t1[tier1] <- p.adjust(p_int[tier1], method = "BH")
}
gene_tbl[, padj_int_gw := padj_int_gw]
gene_tbl[, padj_int_t1 := padj_int_t1]

# Power: only on sex-modulated truth.
power_tab <- data.table(
  TID = TID, delta_main = DELTA_MAIN, delta_int = DELTA_INT, seed = SEED,
  n_sm        = sum(gene_tbl$pattern == "sex_modulated"),
  n_sm_tier1  = sum(gene_tbl$pattern == "sex_modulated" & gene_tbl$tier1_retained),
  n_sig_gw    = sum(gene_tbl$pattern == "sex_modulated" &
                      gene_tbl$padj_int_gw < 0.05, na.rm = TRUE),
  n_sig_t1    = sum(gene_tbl$pattern == "sex_modulated" &
                      gene_tbl$padj_int_t1 < 0.05, na.rm = TRUE)
)
cat("\n== Tier-1 vs genome-wide BH power ==\n")
print(power_tab)

saveRDS(list(
  TID = TID, seed = SEED,
  delta_main = DELTA_MAIN, delta_int = DELTA_INT,
  retention = ret, power = power_tab, gene_tbl = gene_tbl,
  cuts = cuts
), out_rds)
cat("\nWrote:", out_rds, "\n")
