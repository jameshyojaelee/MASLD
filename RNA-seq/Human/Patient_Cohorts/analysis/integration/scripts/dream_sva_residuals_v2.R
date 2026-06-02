#!/usr/bin/env Rscript
# dream_sva_residuals_v2.R
# ---------------------------------------------------------------------------
# Pillar D — SVA-augmented dream sensitivity, v2: SVA estimated on dream
# residuals rather than raw counts. The original supervised svaseq() on raw
# counts double-corrected for cohort (since dream's RE already captures it),
# producing a dramatic over-shrinkage (Jaccard 0.35).  Per Leek 2014 NAR
# svaseq vignette, the proper mixed-model setup runs SVA on the residual
# matrix after the known design has been removed.
#
# Workflow:
#   1. Fast limma::lmFit with `dataset` as a FIXED effect (= matches dream's
#      structure on cohort but in a closed-form fit so we can extract clean
#      residuals quickly without per-gene lme4).
#   2. Run continuous sva() on those residuals — NOT svaseq() — because the
#      residuals are already on the log-CPM scale (continuous, mean-zero).
#   3. Cap n.sv at 5; the cohort fixed effect already removed the dominant
#      structure so few SVs should remain.
#   4. Refit dream with the SVs as additional fixed effects.
#   5. Compare DEG list to primary.
#
# Out: pillar_D_sva_v2_concordance.csv         (per-gene)
#      pillar_D_sva_v2_summary.csv             (1 row, headline)
# ---------------------------------------------------------------------------

t0 <- proc.time()
suppressPackageStartupMessages({
  library(reformulas); library(lme4); library(data.table); library(edgeR); library(yaml); library(sva); library(limma)
})
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({ unlockBinding(fn, ns_lme4)
          assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
          lockBinding(fn, ns_lme4) }, silent = TRUE)
  }
}
suppressPackageStartupMessages({ library(variancePartition); library(BiocParallel) })

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
OUT_DIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Load + filter to mega cohorts ---
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
dge <- dge[, dge$samples$dataset %in% mega_cohorts]
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge), meta_new$sample_id)]
info <- data.frame(
  group_binary = factor(dge$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = droplevels(factor(dge$samples$dataset)),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE)
rownames(info) <- colnames(dge)
na_sex <- is.na(info$inferred_sex)
if (any(na_sex)) { dge <- dge[, !na_sex]; info <- info[!na_sex, , drop = FALSE] }

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# Step 1: voom + limma::lmFit with dataset as fixed effect
cat("[step 1] voom + lmFit with dataset as fixed effect...\n")
design_full <- model.matrix(~ group_binary + inferred_sex + dataset, data = info)
v <- voom(dge, design_full)
fit_lm <- lmFit(v, design_full)

# Step 2: extract fixed-effect residuals = v$E - X*β
cat("[step 2] extract residuals...\n")
fitted_mat <- t(design_full %*% t(fit_lm$coefficients))
resid_mat  <- v$E - fitted_mat
cat(sprintf("    residual matrix: %d genes x %d samples\n", nrow(resid_mat), ncol(resid_mat)))

# Step 3: continuous sva on residuals
# mod has just group_binary (the test variable); mod0 = intercept only
# (sex and dataset have already been absorbed into the residuals).
cat("[step 3] continuous sva on residuals...\n")
mod  <- model.matrix(~ group_binary, data = info)
mod0 <- model.matrix(~ 1,           data = info)
n_sv_be <- tryCatch(num.sv(resid_mat, mod, method = "be"),
                    error = function(e) { cat("num.sv failed:", conditionMessage(e), "\n"); NA_integer_ })
n_sv <- min(if (!is.na(n_sv_be)) n_sv_be else 3L, 5L)
n_sv <- max(n_sv, 1L)
cat(sprintf("    num.sv(be) = %s, capped at %d\n", as.character(n_sv_be), n_sv))
sv_obj <- sva(resid_mat, mod, mod0, n.sv = n_sv)
n_sv_final <- ncol(sv_obj$sv)
cat(sprintf("    Final n.sv used = %d\n", n_sv_final))

# Step 4: refit dream with SVs as additional fixed effects
cat("[step 4] refit dream with SVs as fixed effects...\n")
sv_df <- as.data.frame(sv_obj$sv); colnames(sv_df) <- paste0("SV", seq_len(ncol(sv_df)))
info_sva <- cbind(info, sv_df)
sv_terms <- paste(colnames(sv_df), collapse = " + ")
form_sva <- as.formula(paste("~ group_binary + inferred_sex +", sv_terms, "+ (1|dataset)"))
cat(sprintf("    Refit formula: %s\n", deparse(form_sva)))

v_sva   <- suppressWarnings(voomWithDreamWeights(dge, form_sva, info_sva, BPPARAM = param))
fit_sva <- suppressWarnings(dream(v_sva, form_sva, info_sva, BPPARAM = param))
res_sva <- topTable(fit_sva, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res_sva$gene <- rownames(res_sva); res_sva_dt <- as.data.table(res_sva)
setnames(res_sva_dt, "adj.P.Val", "padj")
fwrite(res_sva_dt, file.path(OUT_DIR, "dream_results_sva_v2.csv"))

# Step 5: compare DEG lists vs primary
PRIMARY_PADJ <- 0.05; PRIMARY_LFC <- 0.5
primary <- fread(file.path(RDIR, "dream_results.csv"))
deg_primary <- primary[padj < PRIMARY_PADJ & abs(logFC) > PRIMARY_LFC, gene]
deg_sva     <- res_sva_dt[padj < PRIMARY_PADJ & abs(logFC) > PRIMARY_LFC, gene]

inter <- length(intersect(deg_primary, deg_sva))
union_n <- length(union(deg_primary, deg_sva))
jaccard <- if (union_n > 0) inter / union_n else NA_real_

m <- merge(primary[, .(gene, logFC_primary = logFC, padj_primary = padj)],
           res_sva_dt[, .(gene, logFC_sva = logFC, padj_sva = padj)], by = "gene")
m[, sig_primary := padj_primary < PRIMARY_PADJ & abs(logFC_primary) > PRIMARY_LFC]
m[, sig_sva     := padj_sva     < PRIMARY_PADJ & abs(logFC_sva)     > PRIMARY_LFC]
m[, sva_concordant := sig_primary & sig_sva & sign(logFC_primary) == sign(logFC_sva)]
fwrite(m, file.path(OUT_DIR, "pillar_D_sva_v2_concordance.csv"))

sum_dt <- data.table(
  n_sv = n_sv_final,
  n_deg_primary = length(deg_primary), n_deg_sva = length(deg_sva),
  intersect = inter, union = union_n, jaccard = round(jaccard, 4),
  rho_logFC = round(cor(m$logFC_primary, m$logFC_sva, use = "pairwise.complete.obs",
                        method = "spearman"), 4))
fwrite(sum_dt, file.path(OUT_DIR, "pillar_D_sva_v2_summary.csv"))
cat("\n=== Pillar D SVA v2 (residuals-based) summary ===\n"); print(sum_dt)
cat(sprintf("\nPass (Jaccard >= 0.85): %s\n",
            ifelse(!is.na(jaccard) && jaccard >= 0.85, "PASS", "FAIL")))
cat(sprintf("Elapsed %.1f min\n", (proc.time() - t0)["elapsed"] / 60))
