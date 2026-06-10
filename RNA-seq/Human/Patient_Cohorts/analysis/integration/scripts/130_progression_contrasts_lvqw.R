#!/usr/bin/env Rscript
# 130_progression_contrasts_lvqw.R
# ---------------------------------------------------------------------------
# LVQW migration of 130_progression_contrasts_dream.R.
#
# Binary progression contrasts re-fit on the canonical limma-voom-quality-
# weighted (LVQW) engine (de_engine_lvqw.R), replacing the dream + (1|dataset)
# random-intercept engine. Per-contrast cohort/sample selection is preserved
# VERBATIM from the dream script — these are disease-internal / progression
# contrasts, so include_in_mega is NOT imposed (EXCEPT for C4/C11, which involve
# Control samples and where the dream script already enforced the canonical
# mega cohort set; that selection is preserved as-is).
#
# Engine swap, design-by-design:
#   dream:  ~ <group> + sex_covar + (1 | dataset)   [dataset = random intercept]
#   LVQW :  ~ dataset + sex_covar + <group>         [dataset = FIXED effect]
#   sex_covar = NUISANCE; dataset = FIXED nuisance; <group> = tested coef.
#   C13 keeps its extra fibrosis nuisance: ~ dataset + sex_covar + fib_numeric + nafl_nash
#
# Every design is built via build_design_guarded(), which reproduces the dream
# script's protective behavior: it droplevels() factors and DROPS the `dataset`
# term when <2 datasets survive the per-contrast filter (the dream script did
# the same by stripping `(1|dataset)` for single-dataset subsets). The caller
# then verifies the target coef survived; if not, the contrast is SKIPPED.
#
# Input template — follows 05h_limma_voom_qw_canonical.R EXACTLY:
#   Load pre-normalized + pre-filtered merged_dge.rds (27,638-gene 9-cohort
#   universe; TMM norm.factors already set). SUBSET columns per contrast.
#   NO re-TMM, NO re-filterByExpr (engine parity with the canonical 05h call).
#
# Tier-2 thresholding: padj < 0.05, NO LFC filter. ALL genes written.
#
# Output schema per contrast (matches dream_results_ashr.csv consumers):
#   gene, logFC, SE, t, P.Value, padj, shrunk_logFC, lfsr, AveExpr, method
#   method = "limma_voom_qw_C2"
#
# STAGING ONLY — writes results/progression/{contrast_id}_lvqw.csv (9 files).
# Backs up each existing {contrast_id}_dream.csv -> *_premigration_backup.csv.
# DOES NOT overwrite any canonical file. Also writes a sanity gate vs the dream
# reference: progression_lvqw_sanity.csv.
#
# PROVENANCE: migrated from 130_progression_contrasts_dream.R (dream engine).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ashr)
  library(data.table)
  library(yaml)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/integration")
ODIR <- file.path(INT, "results/progression")
SCRIPTS <- file.path(INT, "scripts")
dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)

# --- Shared vetted LVQW engine (library-only; no side effects) ---
source(file.path(SCRIPTS, "de_engine_lvqw.R"))

# ===========================================================================
# Input loading — 05h template: pre-normalized + pre-filtered merged_dge.rds.
# The dream script loaded merged_counts_raw.rds and re-ran TMM + filterByExpr
# per contrast; the canonical LVQW engine instead subsets the already-vetted
# merged_dge (fixed 27,638-gene 9-cohort universe, TMM norm.factors set).
# ===========================================================================
cat("Loading merged_dge.rds (pre-normalized, pre-filtered)...\n")
dge_all <- readRDS(file.path(RDIR, "merged_dge.rds"))
cat(sprintf("[in] %d genes x %d samples; %d cohorts\n",
            nrow(dge_all), ncol(dge_all), nlevels(factor(dge_all$samples$dataset))))

# --- Metadata: meta_matched.rds + modeling_metadata NA-fill (verbatim from 130) ---
meta <- as.data.table(readRDS(file.path(RDIR, "meta_matched.rds")))

# Restrict meta to samples present in merged_dge (the QC-passing 9-cohort set).
meta <- meta[sample_id %in% colnames(dge_all)]

# Gene annotations from unified disease signatures (for symbol/gene_type if present;
# NOTE: not part of the fixed LVQW output schema, but joined for downstream parity).
annot_file <- file.path(INT, "results/disease_signatures/unified_disease_signatures.csv")
if (file.exists(annot_file)) {
  gene_annot <- fread(annot_file,
    select = c("gene", "symbol", "gene_type", "mouse_gene_id", "mouse_symbol"))
  gene_annot <- unique(gene_annot, by = "gene")
} else {
  cat("WARNING: Gene annotation file not found — outputs will lack symbol/gene_type columns\n")
  gene_annot <- data.table(gene = character(), symbol = character(),
    gene_type = character(), mouse_gene_id = character(), mouse_symbol = character())
}

# Load modeling metadata for NAS/fibrosis NA-fill (verbatim from 130).
modeling_meta <- fread(file.path(INT, "results/staging_classifier/modeling_metadata.csv"))
meta <- merge(meta, modeling_meta[, .(sample_id, fibrosis_stage, nas_score,
  diagnosis_harmonized, nas_group, fib_ge3, nas_ge5)],
  by = "sample_id", all.x = TRUE, suffixes = c("", ".mm"))
for (col in c("fibrosis_stage", "nas_score", "diagnosis_harmonized")) {
  mm_col <- paste0(col, ".mm")
  if (mm_col %in% names(meta)) {
    na_mask <- is.na(meta[[col]]) | meta[[col]] == ""
    if (any(na_mask)) meta[[col]][na_mask] <- meta[[mm_col]][na_mask]
    meta[, (mm_col) := NULL]
  }
}

# Prepare sex covariate (verbatim from 130).
meta[, sex_covar := inferred_sex]
na_sex <- is.na(meta$sex_covar) | meta$sex_covar == ""
if (any(na_sex) && "sex" %in% names(meta)) {
  meta$sex_covar[na_sex] <- meta$sex[na_sex]
}
meta[, sex_covar := factor(sex_covar)]

# group_binary lives on dge_all$samples; attach to meta for the Control contrasts.
gb <- data.table(sample_id = colnames(dge_all),
                 group_binary = as.character(dge_all$samples$group_binary))
meta <- merge(meta, gb, by = "sample_id", all.x = TRUE, suffixes = c("", ".dge"))
if ("group_binary.dge" %in% names(meta)) {
  na_gb <- is.na(meta$group_binary) | meta$group_binary == ""
  if (any(na_gb)) meta$group_binary[na_gb] <- meta$group_binary.dge[na_gb]
  meta[, group_binary.dge := NULL]
}

cat(sprintf("Total samples (merged_dge ∩ meta): %d\n", nrow(meta)))
cat(sprintf("  with fibrosis staging: %d\n", sum(!is.na(meta$fibrosis_stage))))
cat(sprintf("  with NAS score: %d\n", sum(!is.na(meta$nas_score))))
cat(sprintf("  with diagnosis_harmonized: %d\n",
  sum(!is.na(meta$diagnosis_harmonized) & meta$diagnosis_harmonized != "")))

METHOD_TAG <- "limma_voom_qw_C2"

# ===========================================================================
# Helper: run the LVQW engine for one binary contrast.
#   contrast_id  : output stem
#   meta_sub     : per-contrast metadata subset (selection already applied)
#   group_col    : the binary grouping column in meta_sub
#   group_levels : c(reference, tested)
#   coef_name    : the design column to test (= paste0(group_col, group_levels[2]))
#   extra_terms  : additional nuisance terms to put on the RHS BEFORE group_col
#                  (e.g. "fib_numeric" for C13). Default none.
#
# Design built via build_design_guarded(~ dataset + sex_covar + extra + group),
# which drops `dataset` (and any rank-deficient column) when <2 datasets survive
# — reproducing the dream script's single-dataset (1|dataset)-stripping fallback.
# Returns the per-contrast result data.table, or NULL if SKIPPED.
# ===========================================================================
run_lvqw_contrast <- function(contrast_id, meta_sub, group_col, group_levels,
                              coef_name, description, extra_terms = character(0)) {
  cat(sprintf("\n%s\n", strrep("=", 60)))
  cat(sprintf("  CONTRAST %s: %s\n", contrast_id, description))
  cat(sprintf("%s\n", strrep("=", 60)))

  meta_sub <- copy(meta_sub)
  meta_sub[[group_col]] <- factor(meta_sub[[group_col]], levels = group_levels)

  # Per-dataset distribution (pre-filter).
  cat("  Per-dataset distribution:\n")
  for (ds in sort(unique(meta_sub$dataset))) {
    sub <- meta_sub[dataset == ds]
    n_per_group <- table(sub[[group_col]])
    cat(sprintf("    %s: %s\n", ds,
      paste(paste0(names(n_per_group), "=", n_per_group), collapse = ", ")))
  }

  # Per-contrast dataset filter: >= 2 samples in EACH group (verbatim from 130).
  ds_counts <- meta_sub[, .(
    n1 = sum(.SD[[group_col]] == group_levels[1]),
    n2 = sum(.SD[[group_col]] == group_levels[2])
  ), by = dataset, .SDcols = group_col]
  valid_ds <- ds_counts[n1 >= 2 & n2 >= 2, dataset]
  if (length(valid_ds) == 0) {
    cat("  SKIPPED: no valid datasets (no dataset has >=2 per group).\n")
    return(NULL)
  }
  if (length(valid_ds) < 2) {
    cat(sprintf("  NOTE: only %d valid dataset — `dataset` term will be dropped by the guard.\n",
                length(valid_ds)))
  }
  meta_sub <- meta_sub[dataset %in% valid_ds]
  cat(sprintf("  Valid datasets: %d (%s)\n", length(valid_ds), paste(valid_ds, collapse = ", ")))

  # Subset merged_dge — NO re-TMM, NO re-filter (05h template).
  keep_samples <- intersect(meta_sub$sample_id, colnames(dge_all))
  meta_sub <- meta_sub[sample_id %in% keep_samples]
  dge <- dge_all[, keep_samples]

  # Order metadata to match dge columns; attach modeling covariates.
  m_ord <- meta_sub[match(colnames(dge), meta_sub$sample_id)]
  info <- data.frame(
    dataset   = factor(as.character(m_ord$dataset)),
    sex_covar = factor(as.character(m_ord$sex_covar)),
    stringsAsFactors = FALSE
  )
  info[[group_col]] <- factor(as.character(m_ord[[group_col]]), levels = group_levels)
  for (et in extra_terms) info[[et]] <- m_ord[[et]]
  rownames(info) <- colnames(dge)

  n1 <- sum(info[[group_col]] == group_levels[1])
  n2 <- sum(info[[group_col]] == group_levels[2])
  cat(sprintf("  Samples: %d (%s=%d, %s=%d)\n",
    ncol(dge), group_levels[1], n1, group_levels[2], n2))
  cat(sprintf("  Genes (fixed 9-cohort universe, no re-filter): %d\n", nrow(dge)))

  # Build guarded design: dataset (fixed) + sex_covar + [extra nuisance] + group.
  rhs_terms <- c("dataset", "sex_covar", extra_terms, group_col)
  gd <- build_design_guarded(info, rhs_terms)
  design <- gd$design
  cat(sprintf("  Design: %s\n", gd$formula_used))
  if (length(gd$dropped)) {
    cat(sprintf("  Guard dropped: %s\n", paste(gd$dropped, collapse = ", ")))
  }

  # CALLER CONTRACT: verify the tested coef survived design construction.
  if (!coef_name %in% colnames(design)) {
    cat(sprintf("  SKIPPED: target coef '%s' did not survive design construction (cols: %s).\n",
                coef_name, paste(colnames(design), collapse = ", ")))
    return(NULL)
  }

  # Fit the vetted LVQW engine (voomWithQualityWeights -> eBayes -> ashr).
  res <- fit_lvqw(dge, design, coef = coef_name, do_ashr = TRUE)
  res[, method := METHOD_TAG]

  # Tier-2 summary: padj<0.05, no LFC filter.
  n_sig_05 <- sum(res$padj < 0.05, na.rm = TRUE)
  n_sig_01 <- sum(res$padj < 0.1, na.rm = TRUE)
  n_up   <- sum(res$padj < 0.05 & res$logFC > 0, na.rm = TRUE)
  n_down <- sum(res$padj < 0.05 & res$logFC < 0, na.rm = TRUE)
  cat(sprintf("  DEGs (padj<0.05): %d (Up: %d, Down: %d)\n", n_sig_05, n_up, n_down))
  cat(sprintf("  DEGs (padj<0.1): %d\n", n_sig_01))

  # Back up existing dream canonical, then write LVQW staging output.
  dream_file <- file.path(ODIR, sprintf("%s_dream.csv", contrast_id))
  if (file.exists(dream_file)) {
    bk <- file.path(ODIR, sprintf("%s_premigration_backup.csv", contrast_id))
    if (!file.exists(bk)) file.copy(dream_file, bk, overwrite = FALSE)
  }
  out_file <- file.path(ODIR, sprintf("%s_lvqw.csv", contrast_id))
  fwrite(res, out_file)
  cat(sprintf("  Saved: %s\n", basename(out_file)))

  gc()
  res[, contrast_id := contrast_id]
  res[]
}

# ===========================================================================
# Per-contrast sample selection — VERBATIM from 130 (only the engine differs).
# ===========================================================================

# --- C4/C11 mega cohort set (Control-bearing contrasts) — preserved as-is ---
ycfg <- yaml::read_yaml(file.path(ROOT, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
excl_dvc <- setdiff(unique(meta$dataset), mega_cohorts)
cat(sprintf("\nControl-bearing contrasts (C4/C11) restrict to mega cohorts: %s\n",
            paste(mega_cohorts, collapse = ", ")))

# ---- C3: F3-F4 vs F0-F2 (Advanced vs Early Fibrosis, F2 in Early) ----
meta_c3 <- meta[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta_c3[, fib_group := fifelse(fibrosis_stage >= 3, "Advanced", "Early")]
c3_res <- run_lvqw_contrast(
  contrast_id = "c3_adv_vs_early_fib", meta_sub = meta_c3,
  group_col = "fib_group", group_levels = c("Early", "Advanced"),
  coef_name = "fib_groupAdvanced",
  description = "F3-F4 vs F0-F2 (Advanced vs Early Fibrosis, F2 included in Early)")

# ---- C4: NAFL vs Control (Disease Initiation) ----
stopifnot("group_binary" %in% names(meta))
meta_c4 <- meta[!is.na(diagnosis_harmonized) & diagnosis_harmonized != ""]
meta_c4 <- meta_c4[diagnosis_harmonized == "NAFL" | group_binary == "Control"]
meta_c4 <- meta_c4[!dataset %in% excl_dvc]
meta_c4[, nafl_ctrl := fifelse(group_binary == "Control", "Control", "NAFL")]
c4_res <- run_lvqw_contrast(
  contrast_id = "c4_nafl_vs_ctrl", meta_sub = meta_c4,
  group_col = "nafl_ctrl", group_levels = c("Control", "NAFL"),
  coef_name = "nafl_ctrlNAFL",
  description = "NAFL vs Control (Disease Initiation)")

# ---- C5: NAS >= 5 vs NAS < 5 (Clinical NASH Threshold) ----
meta_c5 <- meta[!is.na(nas_score)]
meta_c5[, nas_binary := fifelse(nas_score >= 5, "NAS_high", "NAS_low")]
c5_res <- run_lvqw_contrast(
  contrast_id = "c5_nas_ge5_vs_lt5", meta_sub = meta_c5,
  group_col = "nas_binary", group_levels = c("NAS_low", "NAS_high"),
  coef_name = "nas_binaryNAS_high",
  description = "NAS >= 5 vs NAS < 5 (Clinical NASH Threshold)")

# ---- C6: NASH + F3-F4 vs NAFL + F0-F1 (Extreme Endpoints) ----
meta_c6 <- meta[!is.na(diagnosis_harmonized) & !is.na(fibrosis_stage)]
meta_c6 <- meta_c6[
  (diagnosis_harmonized %in% c("NASH", "Borderline") & fibrosis_stage >= 3) |
  (diagnosis_harmonized == "NAFL" & fibrosis_stage <= 1)
]
meta_c6[, extreme_group := fifelse(
  diagnosis_harmonized == "NAFL" & fibrosis_stage <= 1, "Early_NAFL", "Advanced_NASH")]
c6_res <- run_lvqw_contrast(
  contrast_id = "c6_extreme_endpoints", meta_sub = meta_c6,
  group_col = "extreme_group", group_levels = c("Early_NAFL", "Advanced_NASH"),
  coef_name = "extreme_groupAdvanced_NASH",
  description = "NASH + F3-F4 vs NAFL + F0-F1 (Extreme Endpoints)")

# ---- C8: F4 vs F0-F3 (Cirrhosis Binary) ----
meta_c8 <- meta[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta_c8[, cirrhosis := fifelse(fibrosis_stage == 4, "Cirrhotic", "Non_cirrhotic")]
c8_res <- run_lvqw_contrast(
  contrast_id = "c8_cirrhosis", meta_sub = meta_c8,
  group_col = "cirrhosis", group_levels = c("Non_cirrhotic", "Cirrhotic"),
  coef_name = "cirrhosisCirrhotic",
  description = "F4 vs F0-F3 (Cirrhosis Binary)")

# ---- C9: F2-F4 vs F0-F1 (F2 Inflection Point) ----
meta_c9 <- meta[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta_c9[, f2_group := fifelse(fibrosis_stage >= 2, "Significant_fib", "Minimal_fib")]
c9_res <- run_lvqw_contrast(
  contrast_id = "c9_f2_inflection", meta_sub = meta_c9,
  group_col = "f2_group", group_levels = c("Minimal_fib", "Significant_fib"),
  coef_name = "f2_groupSignificant_fib",
  description = "F2-F4 vs F0-F1 (F2 Inflection Point)")

# ---- C11: NASH (incl. Borderline) vs Control ----
stopifnot("group_binary" %in% names(meta))
meta_c11 <- meta[!is.na(diagnosis_harmonized) & diagnosis_harmonized != ""]
meta_c11 <- meta_c11[diagnosis_harmonized %in% c("NASH", "Borderline") | group_binary == "Control"]
meta_c11 <- meta_c11[!dataset %in% excl_dvc]
meta_c11[, nash_ctrl := fifelse(group_binary == "Control", "Control", "NASH")]
c11_res <- run_lvqw_contrast(
  contrast_id = "c11_nash_vs_ctrl", meta_sub = meta_c11,
  group_col = "nash_ctrl", group_levels = c("Control", "NASH"),
  coef_name = "nash_ctrlNASH",
  description = "NASH (including Borderline) vs Control")

# ---- C12: Early NASH (F0-F2) vs Late NASH (F3-F4) ----
meta_c12 <- meta[diagnosis_harmonized %in% c("NASH", "Borderline") &
                 !is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta_c12[, nash_fib_group := fifelse(fibrosis_stage <= 2, "Early_NASH", "Late_NASH")]
c12_res <- run_lvqw_contrast(
  contrast_id = "c12_early_vs_late_nash", meta_sub = meta_c12,
  group_col = "nash_fib_group", group_levels = c("Early_NASH", "Late_NASH"),
  coef_name = "nash_fib_groupLate_NASH",
  description = "Early NASH (F0-F2) vs Late NASH (F3-F4) — within-NASH fibrosis progression")

# ---- C13: NASH vs NAFL (fibrosis-adjusted) — extra fib_numeric nuisance ----
meta_c13 <- meta[diagnosis_harmonized %in% c("NAFL", "NASH", "Borderline") &
                 !is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta_c13[, nafl_nash := fifelse(diagnosis_harmonized == "NAFL", "NAFL", "NASH")]
meta_c13[, fib_numeric := as.numeric(fibrosis_stage)]
c13_res <- run_lvqw_contrast(
  contrast_id = "c13_nash_vs_nafl_fib_adj", meta_sub = meta_c13,
  group_col = "nafl_nash", group_levels = c("NAFL", "NASH"),
  coef_name = "nafl_nashNASH",
  description = "NASH vs NAFL (fibrosis-adjusted; ~ dataset + sex_covar + fib_numeric + nafl_nash)",
  extra_terms = "fib_numeric")

# ===========================================================================
# SANITY GATE vs the dream reference: per contrast n DEG, direction concordance,
# logFC Spearman rho.
#
# JOIN KEY: the dream {contrast}_dream.csv files key on UNVERSIONED ENSG
# (ENSG00000241860) while the LVQW outputs key on VERSIONED ENSG
# (ENSG00000241860.8) inherited from merged_dge. Joining on the raw `gene`
# column gives n_shared=0. We strip the version suffix on BOTH sides and join
# on the stripped key.
#
# DEGENERATE-REFERENCE FLAG: the dream script re-ran filterByExpr() on raw
# counts per contrast, producing gene universes from ~26.6k to 40.6k — several
# LARGER than the curated 27,638-gene merged_dge universe the LVQW fits use.
# A dream reference with dream_n_genes far above the universe tested many genes
# LVQW never includes, so low overlap is structural, not a regression. Flagged
# when dream_n_genes > 1.10 * universe. C13 is additionally flagged as a
# covariate-adjusted (fibrosis-adjusted) contrast where fewer DEGs / weaker rho
# are expected by design.
# ===========================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("  SANITY GATE (LVQW vs dream reference)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

UNIVERSE_N <- nrow(dge_all)            # merged_dge fixed universe used by all fits
strip_ver  <- function(x) sub("[.][0-9]+$", "", x)

all_res <- list(
  c3_adv_vs_early_fib      = c3_res,
  c4_nafl_vs_ctrl          = c4_res,
  c5_nas_ge5_vs_lt5        = c5_res,
  c6_extreme_endpoints     = c6_res,
  c8_cirrhosis             = c8_res,
  c9_f2_inflection         = c9_res,
  c11_nash_vs_ctrl         = c11_res,
  c12_early_vs_late_nash   = c12_res,
  c13_nash_vs_nafl_fib_adj = c13_res
)

sanity <- rbindlist(lapply(names(all_res), function(cid) {
  lv <- all_res[[cid]]
  dream_file <- file.path(ODIR, sprintf("%s_dream.csv", cid))
  if (is.null(lv)) {
    return(data.table(contrast_id = cid, status = "LVQW_SKIPPED"))
  }
  lv <- copy(lv); lv[, eb := strip_ver(gene)]
  setorder(lv, padj); lv <- unique(lv, by = "eb")
  lvqw_n_genes <- nrow(lv); lvqw_n_deg05 <- sum(lv$padj < 0.05, na.rm = TRUE)

  if (!file.exists(dream_file)) {
    return(data.table(contrast_id = cid, status = "NO_DREAM_REF",
      lvqw_n_genes = lvqw_n_genes, lvqw_n_deg05 = lvqw_n_deg05))
  }
  dr <- fread(dream_file); dr[, eb := strip_ver(gene)]
  setorder(dr, padj); dr <- unique(dr, by = "eb")
  dream_n_genes <- nrow(dr); dream_n_deg05 <- sum(dr$padj < 0.05, na.rm = TRUE)

  m <- merge(lv[, .(eb, lv_lfc = logFC, lv_padj = padj)],
             dr[, .(eb, dr_lfc = logFC, dr_padj = padj)], by = "eb")
  n_shared <- nrow(m)
  if (n_shared > 2) {
    rho_s <- cor(m$lv_lfc, m$dr_lfc, method = "spearman", use = "complete.obs")
    r_p   <- cor(m$lv_lfc, m$dr_lfc, method = "pearson",  use = "complete.obs")
    dir_shared <- mean(sign(m$lv_lfc) == sign(m$dr_lfc), na.rm = TRUE)
    bs <- m[lv_padj < 0.05 & dr_padj < 0.05]
    n_both_sig <- nrow(bs)
    dir_bothsig <- if (n_both_sig > 0) mean(sign(bs$lv_lfc) == sign(bs$dr_lfc), na.rm = TRUE) else NA_real_
  } else {
    rho_s <- r_p <- dir_shared <- dir_bothsig <- NA_real_; n_both_sig <- NA_integer_
  }

  flags <- character(0)
  if (dream_n_genes > 1.10 * UNIVERSE_N) {
    flags <- c(flags, sprintf("dream_overfiltered(n=%d>univ*1.1)", dream_n_genes))
  }
  if (cid == "c13_nash_vs_nafl_fib_adj") flags <- c(flags, "fib_adjusted_low_DEG_expected")

  data.table(
    contrast_id   = cid, status = "OK",
    lvqw_n_genes  = lvqw_n_genes, lvqw_n_deg05 = lvqw_n_deg05,
    dream_n_genes = dream_n_genes, dream_n_deg05 = dream_n_deg05,
    n_shared = n_shared, shared_frac_of_lvqw = round(n_shared / lvqw_n_genes, 4),
    logFC_spearman = round(rho_s, 4), logFC_pearson = round(r_p, 4),
    direction_concordance_shared = round(dir_shared, 4),
    n_both_sig = n_both_sig,
    direction_concordance_bothsig = round(dir_bothsig, 4),
    degenerate_flag = if (length(flags)) paste(flags, collapse = "; ") else ""
  )
}), fill = TRUE)

fwrite(sanity, file.path(ODIR, "progression_lvqw_sanity.csv"))
cat("\nSANITY TABLE:\n")
print(sanity, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", file.path(ODIR, "progression_lvqw_sanity.csv")))

cat("\n=== Script 130 LVQW migration complete (STAGING — no canonical overwritten) ===\n")
