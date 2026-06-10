#!/usr/bin/env Rscript
# ===========================================================================
# 14.3c2_age_steatosis_sensitivity_lvqw.R
# C2-canonical (limma_voom_qw__C2) re-run of the age + steatosis/NAS confounding
# sensitivity analyses. Replaces the dream comparator with the C2 LVQW engine
# (de_engine_lvqw.R, fit_lvqw) and compares vs canonical_deg_results.csv.
#
# WHAT CHANGED vs the dream-era scripts (14.3_age_sensitivity.R + the steatosis
# sensitivity producer):
#   - Engine: dream() -> fit_lvqw() (voomWithQualityWeights -> lmFit -> eBayes).
#   - Design: C2 fixed-effects ~ dataset + inferred_sex + group_binary
#             (dataset is a FIXED effect, matching 05h, NOT a (1|dataset) RE).
#   - Sensitivity covariate is appended to that SAME C2 design:
#             base   ~ dataset + inferred_sex + group_binary
#             +age   ~ dataset + inferred_sex + group_binary + age_scaled
#             +steat ~ dataset + inferred_sex + group_binary + steatosis_grade  (within cohort)
#             +nas   ~ dataset + inferred_sex + group_binary + nas_scaled
#   - Comparator for the age-associated-gene cross-reference is the C2 canonical
#     DEG set (canonical_deg_results.csv) NOT dream_results.csv.
#   - DEG significance threshold for set-overlap kept at padj<0.1 (the dream-era
#     scripts used 0.1) so OLD->NEW is apples-to-apples; Tier-1 (padj<0.05 &
#     |logFC|>0.5) overlap is ALSO reported for the C2-canonical cross-ref.
#
# ENGINE INPUT NOTE: 05h fits voomWithQualityWeights directly on dge_mega
#   (merged_dge.rds already carries the 9-cohort filterByExpr universe + RLE
#   norm from Script 03 — 05h does NOT re-filter/re-normalize the full mega set).
#   For the MATCHED SUBSET runs here we DO re-run filterByExpr + calcNormFactors
#   (RLE) on the subset, exactly as the dream-era subset scripts did, because the
#   9-cohort universe is not appropriate for a 1-2 cohort subset.
#
# Output: RNA-seq/results/audit_sensitivity/age_confounding_c2/
#         RNA-seq/results/audit_sensitivity/steatosis_sensitivity_c2/
#
# SLURM: io+interactive, <=128G, 8 CPUs, rnaseq env.
# ===========================================================================

suppressPackageStartupMessages({
  library(edgeR); library(limma); library(ashr)
  library(data.table); library(yaml)
})
set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
SCR  <- file.path(INT, "scripts")
source(file.path(SCR, "de_engine_lvqw.R"))   # fit_lvqw(), library-only

OUT_AGE  <- file.path(BASE, "RNA-seq/results/audit_sensitivity/age_confounding_c2")
OUT_STEA <- file.path(BASE, "RNA-seq/results/audit_sensitivity/steatosis_sensitivity_c2")
dir.create(OUT_AGE,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_STEA, showWarnings = FALSE, recursive = TRUE)

cat("=== 14.3c2: Age + Steatosis/NAS sensitivity on the C2 LVQW engine ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ---------------------------------------------------------------------------
# Load data (identical loading to 05h / 14.3)
# ---------------------------------------------------------------------------
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
meta_matched <- readRDS(file.path(RDIR, "meta_matched.rds"))
meta_unified <- fread(file.path(INT, "metadata/unified_metadata.csv"))
qc <- fread(file.path(INT, "qc/sample_qc_report.csv"))
pass_samples <- qc[pass_technical == TRUE, sample_id]

ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("Mega cohorts:", paste(mega_cohorts, collapse = ", "), "\n")

# C2 canonical comparator
canon <- fread(file.path(RDIR, "canonical_deg_results.csv"))
if (!"padj" %in% names(canon) && "adj.P.Val" %in% names(canon))
  setnames(canon, "adj.P.Val", "padj")
canon_deg_p01  <- canon[padj < 0.1, gene]
canon_tier1    <- canon[padj < 0.05 & abs(logFC) > 0.5, gene]   # C2 raw Tier-1 = 1,853
cat("C2 canonical: genes =", nrow(canon),
    "| padj<0.1 DEGs =", length(canon_deg_p01),
    "| Tier-1 (padj<.05 & |lfc|>.5) =", length(canon_tier1), "\n\n")

# ---------------------------------------------------------------------------
# Helper: build a C2 design on an arbitrary info subset + extra covariate(s)
# and fit the LVQW engine for the disease coefficient. Returns the result dt.
# rhs_extra = character vector of extra covariate column names to append.
# ---------------------------------------------------------------------------
fit_c2 <- function(dge_sub, info, rhs_extra = NULL) {
  terms <- c("dataset", "inferred_sex", "group_binary")
  # drop any single-level factor term so model.matrix stays full-rank
  drop_single <- function(t) {
    x <- info[[t]]
    if ((is.factor(x) || is.character(x)) && nlevels(droplevels(as.factor(x))) < 2L) FALSE else TRUE
  }
  terms <- terms[vapply(terms, drop_single, logical(1))]
  rhs   <- c(terms, rhs_extra)
  form  <- as.formula(paste("~", paste(rhs, collapse = " + ")))
  design <- model.matrix(form, data = info)
  coef_name <- "group_binaryDisease"
  stopifnot(coef_name %in% colnames(design))
  list(res = fit_lvqw(dge_sub, design, coef_name, do_ashr = FALSE),
       formula = deparse(form))
}

# ---------------------------------------------------------------------------
# Concordance metric block shared by both arms (base vs +covariate)
# ---------------------------------------------------------------------------
concordance <- function(res_base, res_cov, padj_thr = 0.1) {
  g  <- intersect(res_base$gene, res_cov$gene)
  b  <- res_base[match(g, gene)]
  c2 <- res_cov[match(g, gene)]
  rho   <- cor(b$logFC, c2$logFC, method = "spearman", use = "complete.obs")
  pear  <- cor(b$logFC, c2$logFC, method = "pearson",  use = "complete.obs")
  rho_t <- cor(b$t,     c2$t,     method = "spearman", use = "complete.obs")
  dirc  <- mean(sign(b$logFC) == sign(c2$logFC), na.rm = TRUE) * 100
  d_b <- b[padj < padj_thr, gene]; d_c <- c2[padj < padj_thr, gene]
  jac <- length(intersect(d_b, d_c)) / max(1L, length(union(d_b, d_c)))
  shift <- c2$logFC - b$logFC
  list(n_common = length(g),
       n_deg_base = length(d_b), n_deg_cov = length(d_c),
       rho = rho, pear = pear, rho_t = rho_t, dir = dirc, jaccard = jac,
       gained = length(setdiff(d_c, d_b)), lost = length(setdiff(d_b, d_c)),
       mean_abs_shift = mean(abs(shift), na.rm = TRUE),
       median_abs_shift = median(abs(shift), na.rm = TRUE),
       max_abs_shift = max(abs(shift), na.rm = TRUE),
       base = b, cov = c2, shift = shift)
}

# ===========================================================================
# ARM 1 — AGE CONFOUNDING (matched subset, C2 +/- age_scaled)
# ===========================================================================
cat("=== ARM 1: AGE CONFOUNDING (C2 +/- age) ===\n")
AGE_DATASETS <- intersect(c("GSE130970", "GSE162694", "GSE174478", "GSE193066"), mega_cohorts)
cat("Age-eligible mega cohorts:", paste(AGE_DATASETS, collapse = ", "), "\n")

meta_age <- meta_unified[dataset %in% AGE_DATASETS & !is.na(age) & age != "" &
                           sample_id %in% pass_samples]
meta_age[, age := as.numeric(age)]
meta_age <- meta_age[!is.na(age)]
age_samples <- intersect(meta_age$sample_id, colnames(dge))
cat("Age-subset samples:", length(age_samples), "\n")

dge_age <- dge[, age_samples]
mm <- meta_age[match(colnames(dge_age), meta_age$sample_id)]
sex_age <- meta_matched$inferred_sex[match(colnames(dge_age), meta_matched$sample_id)]
info_age <- data.frame(
  group_binary = factor(mm$group_binary, levels = c("Control", "Disease")),
  dataset      = factor(mm$dataset),
  inferred_sex = factor(sex_age),
  age_scaled   = scale(as.numeric(mm$age))[, 1],
  row.names    = colnames(dge_age))
# subset-appropriate filter + RLE norm (NOT the 9-cohort universe)
keep <- filterByExpr(dge_age, group = info_age$group_binary)
dge_age <- dge_age[keep, , keep.lib.sizes = FALSE]
dge_age <- calcNormFactors(dge_age, method = "RLE")
cat("Age subset: genes after filter =", nrow(dge_age),
    "| datasets =", nlevels(droplevels(info_age$dataset)), "\n")

f_base <- fit_c2(dge_age, info_age, rhs_extra = NULL)
f_age  <- fit_c2(dge_age, info_age, rhs_extra = "age_scaled")
cat("Base formula :", f_base$formula, "\n")
cat("+age formula :", f_age$formula, "\n")
cm_age <- concordance(f_base$res, f_age$res, padj_thr = 0.1)
cat(sprintf("  rho=%.4f pear=%.4f rho_t=%.4f dir=%.1f%% jaccard=%.4f | DEG base=%d +age=%d (gain %d / lost %d)\n",
            cm_age$rho, cm_age$pear, cm_age$rho_t, cm_age$dir, cm_age$jaccard,
            cm_age$n_deg_base, cm_age$n_deg_cov, cm_age$gained, cm_age$lost))

fwrite(f_base$res, file.path(OUT_AGE, "c2_subset_noage.csv"))
fwrite(f_age$res,  file.path(OUT_AGE, "c2_subset_withage.csv"))
age_metrics <- data.table(
  metric = c("n_subset_samples","n_datasets","n_genes_tested","n_deg_noage","n_deg_age",
             "spearman_rho_logFC","pearson_r_logFC","spearman_rho_tstat",
             "direction_concordance_pct","jaccard_padj01",
             "n_gained_with_age","n_lost_with_age",
             "mean_abs_lfc_shift","median_abs_lfc_shift","max_abs_lfc_shift"),
  value = c(length(age_samples), nlevels(droplevels(info_age$dataset)), cm_age$n_common,
            cm_age$n_deg_base, cm_age$n_deg_cov,
            round(cm_age$rho,4), round(cm_age$pear,4), round(cm_age$rho_t,4),
            round(cm_age$dir,1), round(cm_age$jaccard,4),
            cm_age$gained, cm_age$lost,
            round(cm_age$mean_abs_shift,4), round(cm_age$median_abs_shift,4),
            round(cm_age$max_abs_shift,4)))
fwrite(age_metrics, file.path(OUT_AGE, "age_sensitivity_metrics_c2.csv"))

# --- Age-associated genes vs C2 canonical DEGs (Analysis C) ---
cat("\n--- Age-associated genes vs C2 canonical DEGs ---\n")
lcpm <- cpm(dge_age, log = TRUE)
nds  <- nlevels(droplevels(info_age$dataset))
age_v  <- info_age$age_scaled; ds_v <- info_age$dataset; sex_v <- info_age$inferred_sex
age_res <- rbindlist(lapply(seq_len(nrow(lcpm)), function(i) {
  df <- data.frame(expr = lcpm[i, ], age = age_v, dataset = ds_v, sex = sex_v)
  fit <- tryCatch(
    if (nds > 1) lm(expr ~ age + dataset + sex, data = df) else lm(expr ~ age + sex, data = df),
    error = function(e) NULL)
  if (is.null(fit)) return(data.table(gene = rownames(lcpm)[i], age_coef = NA_real_,
                                       age_se = NA_real_, age_t = NA_real_, age_p = NA_real_))
  s <- summary(fit)$coefficients
  if ("age" %in% rownames(s))
    data.table(gene = rownames(lcpm)[i], age_coef = s["age","Estimate"],
               age_se = s["age","Std. Error"], age_t = s["age","t value"], age_p = s["age","Pr(>|t|)"])
  else data.table(gene = rownames(lcpm)[i], age_coef = NA_real_, age_se = NA_real_,
                  age_t = NA_real_, age_p = NA_real_)
}))
age_res[, age_padj := p.adjust(age_p, method = "BH")]
age_genes_005 <- age_res[age_padj < 0.05, gene]

# Cross-ref vs C2 canonical DEGs (padj<0.1 set, to match dream-era apples-to-apples)
overlap_p01 <- length(intersect(canon_deg_p01, age_genes_005))
n_total <- nrow(age_res)
n_canon <- sum(age_res$gene %in% canon_deg_p01)
n_age   <- length(age_genes_005)
n_both  <- length(intersect(canon_deg_p01, age_genes_005))
fmat <- matrix(c(n_both, n_canon - n_both, n_age - n_both,
                 n_total - n_canon - n_age + n_both), nrow = 2)
fres <- fisher.test(fmat)
# ALSO cross-ref vs C2 Tier-1 (padj<.05 & |lfc|>.5) for the canonical claim
overlap_tier1 <- length(intersect(canon_tier1, age_genes_005))
n_canonT <- sum(age_res$gene %in% canon_tier1)
n_bothT  <- overlap_tier1
fmatT <- matrix(c(n_bothT, n_canonT - n_bothT, n_age - n_bothT,
                  n_total - n_canonT - n_age + n_bothT), nrow = 2)
fresT <- fisher.test(fmatT)

age_res[, is_c2_deg_p01 := gene %in% canon_deg_p01]
age_res[, is_c2_tier1   := gene %in% canon_tier1]
fwrite(age_res, file.path(OUT_AGE, "age_associated_genes_c2.csv"))

xref <- data.table(
  metric = c("n_c2_degs_padj01","n_c2_tier1","n_age_genes_padj005",
             "overlap_c2_p01_x_age","pct_c2_p01_degs_age","fisher_OR_p01","fisher_p_p01",
             "overlap_c2_tier1_x_age","pct_c2_tier1_degs_age","fisher_OR_tier1","fisher_p_tier1"),
  value = c(length(canon_deg_p01), length(canon_tier1), n_age,
            overlap_p01, round(overlap_p01/max(1L,length(canon_deg_p01))*100,1),
            round(fres$estimate,3), signif(fres$p.value,4),
            overlap_tier1, round(overlap_tier1/max(1L,length(canon_tier1))*100,1),
            round(fresT$estimate,3), signif(fresT$p.value,4)))
fwrite(xref, file.path(OUT_AGE, "age_c2_crossref_summary.csv"))
cat(sprintf("  age genes(padj<.05)=%d | C2 DEG(p01)=%d overlap=%d OR=%.2f p=%.2g | C2 Tier-1=%d overlap=%d OR=%.2f\n",
            n_age, length(canon_deg_p01), overlap_p01, fres$estimate, fres$p.value,
            length(canon_tier1), overlap_tier1, fresT$estimate))

# ===========================================================================
# ARM 2 — STEATOSIS / NAS CONFOUNDING (C2 +/- steatosis_grade or nas_scaled)
# ===========================================================================
cat("\n=== ARM 2: STEATOSIS / NAS CONFOUNDING (C2 +/- covariate) ===\n")

# covar_map: optional data.table(sample_id, covar_val) sourced EXTERNALLY (e.g.
# GSE130970 steatosis subscore from the SraRunTable, which is NOT in
# unified_metadata). When NULL, the covariate is read from meta_unified[[covar_col]].
run_grade_arm <- function(label, datasets, covar_col, scale_covar, covar_map = NULL) {
  cat(sprintf("\n--- %s (datasets: %s) ---\n", label, paste(datasets, collapse=", ")))
  if (is.null(covar_map)) {
    ms <- meta_unified[dataset %in% datasets & sample_id %in% pass_samples &
                         !is.na(get(covar_col)) & get(covar_col) != ""]
    ms[[covar_col]] <- as.numeric(ms[[covar_col]])
    ms <- ms[!is.na(get(covar_col))]
  } else {
    ms <- meta_unified[dataset %in% datasets & sample_id %in% pass_samples]
    ms[[covar_col]] <- as.numeric(covar_map$covar_val[match(ms$sample_id, covar_map$sample_id)])
    ms <- ms[!is.na(get(covar_col))]
  }
  samp <- intersect(ms$sample_id, colnames(dge))
  if (length(samp) < 20) { cat("  too few samples; skipping\n"); return(NULL) }
  dge_s <- dge[, samp]
  mmx <- ms[match(colnames(dge_s), ms$sample_id)]
  sexx <- meta_matched$inferred_sex[match(colnames(dge_s), meta_matched$sample_id)]
  cov_raw <- as.numeric(mmx[[covar_col]])
  info <- data.frame(
    group_binary = factor(mmx$group_binary, levels = c("Control", "Disease")),
    dataset      = factor(mmx$dataset),
    inferred_sex = factor(sexx),
    covar        = if (scale_covar) scale(cov_raw)[,1] else cov_raw,
    row.names    = colnames(dge_s))
  # need both control + disease present to estimate the disease coef
  if (length(unique(info$group_binary)) < 2) { cat("  single group; skipping\n"); return(NULL) }
  keep <- filterByExpr(dge_s, group = info$group_binary)
  dge_s <- dge_s[keep, , keep.lib.sizes = FALSE]
  dge_s <- calcNormFactors(dge_s, method = "RLE")
  cat("  samples =", length(samp), "| datasets =", nlevels(droplevels(info$dataset)),
      "| genes =", nrow(dge_s), "\n")
  fb <- fit_c2(dge_s, info, rhs_extra = NULL)
  fc <- fit_c2(dge_s, info, rhs_extra = "covar")
  cm <- concordance(fb$res, fc$res, padj_thr = 0.1)
  cat(sprintf("  base: %s\n  +cov: %s\n", fb$formula, fc$formula))
  cat(sprintf("  rho=%.4f pear=%.4f rho_t=%.4f dir=%.1f%% jaccard=%.4f | DEG base=%d +cov=%d (gain %d/lost %d)\n",
              cm$rho, cm$pear, cm$rho_t, cm$dir, cm$jaccard,
              cm$n_deg_base, cm$n_deg_cov, cm$gained, cm$lost))
  fwrite(fb$res, file.path(OUT_STEA, sprintf("c2_%s_unadjusted.csv", label)))
  fwrite(fc$res, file.path(OUT_STEA, sprintf("c2_%s_adjusted.csv", label)))
  data.table(
    analysis = label, n_samples = length(samp),
    n_datasets = nlevels(droplevels(info$dataset)), n_genes_tested = cm$n_common,
    n_deg_unadjusted = cm$n_deg_base, n_deg_adjusted = cm$n_deg_cov,
    spearman_rho_logFC = round(cm$rho,4), pearson_r_logFC = round(cm$pear,4),
    spearman_rho_tstat = round(cm$rho_t,4), direction_concordance_pct = round(cm$dir,1),
    jaccard_padj01 = round(cm$jaccard,4), n_gained = cm$gained, n_lost = cm$lost,
    mean_abs_lfc_shift = round(cm$mean_abs_shift,4),
    median_abs_lfc_shift = round(cm$median_abs_shift,4),
    max_abs_lfc_shift = round(cm$max_abs_shift,4))
}

# steatosis_grade is NOT in unified_metadata; harvest the GSE130970 subscore
# from its SraRunTable (column `steatosis_grade`, keyed by Run = sample_id),
# matching Script 89_extract_nas_components.R.
gse130970_srr <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/archive/old_data_metadata/metadata/GSE130970_SraRunTable.csv")
steat_map <- NULL
if (file.exists(gse130970_srr)) {
  srr <- fread(gse130970_srr)
  steat_map <- data.table(sample_id = srr$Run,
                          covar_val = suppressWarnings(as.numeric(srr$steatosis_grade)))
  steat_map <- steat_map[!is.na(covar_val)]
  cat(sprintf("\nGSE130970 steatosis_grade harvested for %d samples\n", nrow(steat_map)))
} else {
  cat("\nWARNING: GSE130970 SraRunTable not found; steatosis arm skipped\n")
}

steat_dt <- if (!is.null(steat_map))
  run_grade_arm("steatosis_grade", "GSE130970", "steatosis_grade",
                scale_covar = FALSE, covar_map = steat_map) else NULL
nas_dt   <- run_grade_arm("nas_score",
                          intersect(c("GSE135251","GSE130970","GSE162694"), mega_cohorts),
                          "nas_score", scale_covar = TRUE)
steat_metrics <- rbindlist(list(steat_dt, nas_dt), use.names = TRUE, fill = TRUE)
fwrite(steat_metrics, file.path(OUT_STEA, "steatosis_sensitivity_metrics_c2.csv"))

# ===========================================================================
# SUMMARY
# ===========================================================================
cat("\n============================================================\n")
cat("SUMMARY (C2 LVQW engine)\n")
cat("============================================================\n")
cat(sprintf("AGE      : rho=%.4f jaccard(p01)=%.4f dir=%.1f%% | age-genes(.05) x C2-DEG(p01) overlap=%d OR=%.2f\n",
            cm_age$rho, cm_age$jaccard, cm_age$dir, overlap_p01, fres$estimate))
if (!is.null(steat_dt))
  cat(sprintf("STEATOSIS: rho=%.4f jaccard(p01)=%.4f dir=%.1f%% DEG %d->%d\n",
              steat_dt$spearman_rho_logFC, steat_dt$jaccard_padj01,
              steat_dt$direction_concordance_pct, steat_dt$n_deg_unadjusted, steat_dt$n_deg_adjusted))
if (!is.null(nas_dt))
  cat(sprintf("NAS      : rho=%.4f jaccard(p01)=%.4f dir=%.1f%% DEG %d->%d\n",
              nas_dt$spearman_rho_logFC, nas_dt$jaccard_padj01,
              nas_dt$direction_concordance_pct, nas_dt$n_deg_unadjusted, nas_dt$n_deg_adjusted))
cat("\nOutputs:\n  ", file.path(OUT_AGE, "age_sensitivity_metrics_c2.csv"),
    "\n  ", file.path(OUT_AGE, "age_c2_crossref_summary.csv"),
    "\n  ", file.path(OUT_STEA, "steatosis_sensitivity_metrics_c2.csv"), "\n")
cat("\n=== 14.3c2 completed:", as.character(Sys.time()), "===\n")
