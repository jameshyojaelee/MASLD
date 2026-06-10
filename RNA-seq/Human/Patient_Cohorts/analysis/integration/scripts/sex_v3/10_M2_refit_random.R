#!/usr/bin/env Rscript
# sex_v3/10_M2_refit_random.R
# ---------------------------------------------------------------------------
# Pillar 1 — M2 refit (LVQW fixed-effects engine).
# LVQW re-engineering 2026-06-08: dataset random->fixed.
#   Forks 04_dream_M2_extended.R. The dream-era random-slope-on-interaction
#   `(1 + group_binary + inferred_sex + group_binary:inferred_sex || dataset)`
#   (F2 diagonal) and its F5/F6 fallback ladder are replaced by a single LVQW
#   fixed-effects design with `dataset` as a FIXED effect. Interaction,
#   composition, MI age, and SVA surrogate variables are preserved exactly.
#
# Single imputation (k=1) is used (per Pillar 2 design — across-MI variance
# is dominated by cohort effects, not age; this matches calibration design).
#
# Inputs (cached):
#   intermediates/sex_v3_input.rds (Module 01)
#   intermediates/sva_factors.rds  (Module 02)
#   intermediates/age_mi.rds       (Module 03)
#   dream_M2_v3.rds                (current; for vcov benchmark)
#
# Outputs:
#   dream_M2_v3_random.rds                  -- raw fit metadata (β_int, se_int)
#   interaction_classifier_v5_random.csv    -- per-gene v5 schema + *_random cols + class_v5_random
# ---------------------------------------------------------------------------
set.seed(42)
# LVQW re-engineering 2026-06-08: dataset random->fixed.
suppressPackageStartupMessages({
  library(limma); library(data.table); library(edgeR); library(ashr)
})

BASE  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
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
# Contrast-aware Tier-1 universe (the disease-effect DEG list the conditional
# interaction-padj filter sits on top of). For disease_vs_ctrl this is the
# canonical dream_results_ashr.csv at the integration root; for the new
# contrasts (mash_vs_masl etc.) it's the per-contrast disease-side dream
# output documented in `sex_v3_utils.R::contrast_spec()`.
.cspec <- contrast_spec()
TIER1 <- file.path(BASE, .cspec$tier1_csv)
TIER1_PADJ_COL  <- .cspec$tier1_padj_col   # "adj.P.Val" or "padj"
TIER1_LOGFC_COL <- .cspec$tier1_logFC_col  # "logFC"

OUT_RDS <- file.path(SEXV3, "dream_M2_v3_random.rds")
OUT_CSV <- file.path(SEXV3, "interaction_classifier_v5_random.csv")

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
# LVQW re-engineering 2026-06-08: dataset random->fixed.
# limma LVQW is single-threaded; cpu count retained for log parity only.
cat("Using", ncpus, "CPU core(s) reported by SLURM (limma LVQW is single-threaded)\n")

inp    <- readRDS(file.path(IDIR, "sex_v3_input.rds"))
sva    <- readRDS(file.path(IDIR, "sva_factors.rds"))
age_mi <- readRDS(file.path(IDIR, "age_mi.rds"))
dge   <- inp$dge
info0 <- inp$meta
sv_mat <- sva$sv
n_sv  <- ncol(sv_mat)
sv_df <- as.data.frame(sv_mat); colnames(sv_df) <- paste0("SV", seq_len(n_sv))
info_template <- info0
info_template$age_imputed <- age_mi$age_list[[1]]
info_template <- cbind(info_template, sv_df)
rownames(info_template) <- rownames(info0)

cat("Samples:", ncol(dge), "  Genes:", nrow(dge), "  SVs:", n_sv, "\n")
stopifnot(identical(rownames(info_template), colnames(dge)))

# LVQW re-engineering 2026-06-08: dataset random->fixed.
# Single LVQW fixed-effects design (dataset fixed); no random-slope ladder.
sv_terms <- paste(colnames(sv_df), collapse = " + ")
form_str <- paste0(
  "~ dataset",
  " + group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed + ", sv_terms
)
cat("\nFormula (LVQW fixed-effects design, dataset fixed):\n  ", form_str, "\n")
form_fixed <- as.formula(form_str)

design <- model.matrix(form_fixed, data = info_template)
stopifnot(nrow(design) == ncol(dge))
if (qr(design)$rank < ncol(design)) {
  stop("LVQW design is rank-deficient; inspect dataset/SV/composition collinearity.")
}

re_variant_used <- "lvqw_fixed"
cat("\n[1] voomWithQualityWeights (fixed-effects design)...\n")
t0 <- Sys.time()
v_fixed <- suppressWarnings(voomWithQualityWeights(dge, design))
cat("  voom elapsed:", format(Sys.time() - t0), "\n")

cat("\n[2] lmFit + eBayes...\n")
t0 <- Sys.time()
fit <- eBayes(lmFit(v_fixed, design))
cat("  lmFit+eBayes elapsed:", format(Sys.time() - t0), "\n")
stopifnot(!is.null(fit))
cat("Final variant used:", re_variant_used, "\n")

# Extract β_F, β_M, β_int, SEs (mirror 04 logic)
all_coefs <- colnames(fit$coefficients)
group_coef <- setdiff(grep("^group_binary", all_coefs, value = TRUE),
                     grep(":", all_coefs, value = TRUE))[1]
int_coef <- grep(":", grep("^group_binary", all_coefs, value = TRUE), value = TRUE)[1]
cat("group_coef (=β_F):", group_coef, "  int_coef:", int_coef, "\n")

beta_F   <- fit$coefficients[, group_coef]
beta_int <- fit$coefficients[, int_coef]
beta_M   <- beta_F + beta_int
# RAW (unmoderated) Wald SE — match dream-era choice. sigma is unmoderated.
raw_se <- fit$stdev.unscaled * fit$sigma   # genes x K, sigma recycled per row
se_F   <- raw_se[, group_coef]
se_int <- raw_se[, int_coef]
# LVQW re-engineering 2026-06-08: dataset random->fixed.
# Per-gene Cov(β_F, β_int) = cov.coefficients[F,int] (shared UNSCALED) * sigma[g]^2.
ccm <- fit$cov.coefficients
stopifnot(!is.null(ccm), group_coef %in% rownames(ccm), int_coef %in% colnames(ccm))
cov_F_int <- as.numeric(ccm[group_coef, int_coef]) * (fit$sigma)^2
names(cov_F_int) <- rownames(fit$coefficients)
cov_F_int <- cov_F_int[names(beta_F)]
var_M <- se_int^2 + se_F^2 + 2 * cov_F_int
se_M  <- sqrt(pmax(var_M, .Machine$double.eps))
t_int <- beta_int / se_int
p_int <- 2 * pnorm(-abs(t_int))

cat("\nPer-gene estimates ready:", length(beta_F), "genes\n")
cat("  median |t_int|:", round(median(abs(t_int), na.rm=TRUE), 3), "\n")
cat("  median se_int:", round(median(se_int, na.rm=TRUE), 4), "\n")

# Tier-1 universe → padj_int_5k
tier1 <- fread(TIER1)
# Use contrast-aware column names (padj for Disease-vs-Ctrl, adj.P.Val for new contrasts)
tier1_genes <- tier1[!is.na(get(TIER1_PADJ_COL)) & get(TIER1_PADJ_COL) < 0.05 &
                      abs(get(TIER1_LOGFC_COL)) > 0.3, gene]
cat("Tier-1 DEG universe:", length(tier1_genes), "\n")
padj_int_5k <- rep(NA_real_, length(beta_F))
in_tier1 <- names(beta_F) %in% tier1_genes
padj_int_5k[in_tier1] <- p.adjust(p_int[in_tier1], method = "BH")
padj_int_full <- p.adjust(p_int, method = "BH")

cat("padj_int_5k < 0.05:", sum(padj_int_5k < 0.05, na.rm=TRUE),
    "  < 0.10:", sum(padj_int_5k < 0.10, na.rm=TRUE),
    "  < 0.20:", sum(padj_int_5k < 0.20, na.rm=TRUE), "\n")
cat("padj_int_full < 0.05:", sum(padj_int_full < 0.05, na.rm=TRUE), "\n")

# Compare with current v5 (optional — only present for disease_vs_ctrl canonical
# contrast; for other contrasts the legacy v5 outputs don't exist)
v5_int_csv <- file.path(SEXV3, "interaction_classifier_v5.csv")
v5_int <- NULL
spearman <- NA_real_; se_ratio <- NA_real_; acceptance_pass <- NA
if (file.exists(v5_int_csv)) {
  v5_int <- fread(v5_int_csv)
  common <- intersect(names(beta_F), v5_int$gene)
  cmp_idx <- match(common, v5_int$gene)
  v5_t <- v5_int$t_int[cmp_idx]
  new_t <- t_int[match(common, names(beta_F))]
  spearman <- cor(v5_t, new_t, method = "spearman", use = "complete.obs")
  v5_se <- v5_int$se_int[cmp_idx]
  new_se <- se_int[match(common, names(beta_F))]
  se_ratio <- median(new_se / v5_se, na.rm = TRUE)
  cat("\n=== Comparison vs v5 ===\n")
  cat("  Spearman(t_int_v5, t_int_random):", round(spearman, 3), "\n")
  cat("  Median se_int_random / se_int_v5:", round(se_ratio, 3), "\n")
  acceptance_pass <- spearman >= 0.85 && se_ratio >= 1.0 && se_ratio <= 1.5
  cat("  Acceptance gate (Spearman >= 0.85 AND ratio in [1.0, 1.5]):",
      if (isTRUE(acceptance_pass)) "PASS" else "FAIL", "\n")
} else {
  cat("\n=== Skipping v5 comparison (interaction_classifier_v5.csv not present) ===\n")
}

# Apply v5 decision tree on new beta_int + padj_int_5k -> class_v5_random
# (optional — only run when stratified_ashr_v4.csv + v3 classification exist;
# these are legacy artifacts produced by the disease_vs_ctrl pipeline only)
strat_csv  <- file.path(SEXV3, "stratified_ashr_v4.csv")
v3meta_csv <- file.path(SEXV3, "sex_deg_classification_v3.csv")
have_legacy <- file.exists(strat_csv) && file.exists(v3meta_csv)
if (have_legacy) {
  strat <- fread(strat_csv)
  v3meta <- fread(v3meta_csv,
                  select = c("gene", "gene_symbol", "ensembl_base", "chr",
                              "chr_category", "gene_biotype"))
} else {
  cat("\n=== Skipping v5 decision-tree classification (legacy CSVs missing) ===\n")
  strat <- data.table(gene = names(beta_F),
                      beta_F_strat = NA_real_, se_F_strat = NA_real_, lfsr_F_strat = NA_real_,
                      beta_M_strat = NA_real_, se_M_strat = NA_real_, lfsr_M_strat = NA_real_,
                      power_M_at_F = NA_real_)
  v3meta <- data.table(gene = names(beta_F),
                       gene_symbol  = NA_character_,
                       ensembl_base = NA_character_,
                       chr          = NA_character_,
                       chr_category = NA_character_,
                       gene_biotype = NA_character_)
}
out <- data.table(
  gene = names(beta_F),
  beta_int_random = beta_int,
  se_int_random   = se_int,
  t_int_random    = t_int,
  p_int_random    = p_int,
  padj_int_5k_random   = padj_int_5k,
  padj_int_full_random = padj_int_full,
  re_variant_used = re_variant_used
)
out <- merge(out, strat[, .(gene, beta_F_strat, se_F_strat, lfsr_F_strat,
                             beta_M_strat, se_M_strat, lfsr_M_strat,
                             power_M_at_F)], by = "gene", all.x = TRUE)
out <- merge(out, v3meta, by = "gene", all.x = TRUE)

# Decision tree (mirrors 07_interaction_test_v5.R)
classify_v5 <- function(bF, bM, pI, lF, lM, pow) {
  if (is.na(pI) || is.na(bF) || is.na(bM)) {
    if (!is.na(lF) && !is.na(lM)) {
      if (lF < 0.05 && lM < 0.05 && sign(bF) == sign(bM)) return("Concordant")
      if (lF < 0.05 && lM < 0.05 && sign(bF) != sign(bM)) return("Divergent")
    }
    return("Not_DEG_or_NotInUniverse")
  }
  same_sign <- (sign(bF) == sign(bM)) && bF != 0 && bM != 0
  opp_sign  <- (sign(bF) != sign(bM)) && bF != 0 && bM != 0
  if (pI >= 0.20) {
    if (!is.na(lF) && !is.na(lM)) {
      if (lF < 0.05 && lM < 0.05 && same_sign) return("Concordant")
      if (lF < 0.05 || lM < 0.05)              return("Concordant_single_arm")
    }
    return("Not_DEG")
  }
  if (opp_sign) {
    if (!is.na(lF) && !is.na(lM) && lF < 0.05 && lM < 0.05) return("Divergent")
    if ((!is.na(lF) && lF < 0.05) || (!is.na(lM) && lM < 0.05))
      return("Divergent_one_sided")
  }
  if (same_sign) {
    if (abs(bF) >= 2*abs(bM) && !is.na(lF) && lF < 0.05 && !is.na(lM)) {
      if (lM > 0.20) return("Female_biased")
      if (lM >= 0.05 && lM <= 0.20 && !is.na(pow) && pow < 0.5)
        return("Female_biased_M_underpowered")
    }
    if (abs(bM) >= 2*abs(bF) && !is.na(lM) && lM < 0.05 && !is.na(lF)) {
      if (lF > 0.20) return("Male_biased")
      if (lF >= 0.05 && lF <= 0.20) return("Male_biased_F_underpowered")
    }
    return("Sex_modifier")
  }
  return("Uncertain")
}
# Use Layer-2 ashr lfsr from stratified (UNCHANGED — they don't depend on random slope)
out[, beta_F := beta_F_strat]
out[, beta_M := beta_M_strat]
out[, class_v5_random := mapply(classify_v5,
                                 beta_F_strat, beta_M_strat,
                                 padj_int_5k_random, lfsr_F_strat, lfsr_M_strat,
                                 power_M_at_F)]

# chrY F-only → Male_biased reclass
out[chr == "chrY" & class_v5_random == "Female_biased",
    class_v5_random := "Male_biased"]

cat("\n== class_v5_random distribution ==\n")
print(table(out$class_v5_random, useNA = "ifany"))

# Sanity check: compare counts (only when legacy v5 csv loaded)
if (!is.null(v5_int)) {
  v5_cls <- v5_int$class_v5_interaction
  concord <- sum(out$class_v5_random == v5_cls[match(out$gene, v5_int$gene)],
                  na.rm = TRUE)
  total <- sum(!is.na(out$class_v5_random) & !is.na(v5_cls[match(out$gene, v5_int$gene)]))
  cat("Class concordance v5 vs v5_random:",
      round(concord / total * 100, 1), "%\n")
}

# LVQW re-engineering 2026-06-08: dataset random->fixed.
# Save the LVQW voom object + the fixed-effects design so P2 (perm-FDR) can
# reuse them. P2 reads $v_voom, $form_used, $info_template, $re_variant_used;
# we additionally stash $design (the model.matrix) for the LVQW lmFit refit.
v_voom_used <- v_fixed
form_used   <- form_fixed
stopifnot(!is.null(v_voom_used))

OUT_RDS_TMP <- paste0(OUT_RDS, ".tmp")
saveRDS(list(fit = fit, re_variant_used = re_variant_used,
             engine = "lvqw_fixed",
             beta_F = beta_F, beta_M = beta_M, beta_int = beta_int,
             se_F = se_F, se_M = se_M, se_int = se_int,
             cov_F_int = cov_F_int, t_int = t_int,
             v_voom = v_voom_used,
             form_used = form_used,
             design = design,
             info_template = info_template,
             tier1_genes = tier1_genes,
             spearman_vs_v5 = spearman,
             se_ratio_vs_v5 = se_ratio,
             acceptance_pass = acceptance_pass),
        OUT_RDS_TMP)
file.rename(OUT_RDS_TMP, OUT_RDS)
cat("\nWrote:", OUT_RDS, "  voom variant:", re_variant_used, "\n")

OUT_CSV_TMP <- paste0(OUT_CSV, ".tmp")
fwrite(out, OUT_CSV_TMP)
file.rename(OUT_CSV_TMP, OUT_CSV)
cat("Wrote:", OUT_CSV, "  rows:", nrow(out), "  cols:", ncol(out), "\n")
