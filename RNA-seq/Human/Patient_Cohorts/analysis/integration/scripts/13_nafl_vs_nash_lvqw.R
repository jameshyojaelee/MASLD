#!/usr/bin/env Rscript
# ===========================================================================
# 13_nafl_vs_nash_lvqw.R — Section C (pooled NAFL-vs-NASH contrast) harmonized
# to the canonical limma-voom-quality-weighted (LVQW) engine.
#
# This REPLACES ONLY Section C of 13_nafl_vs_nash_de.R (the dream pooled
# contrast). Sections A (per-study limma) and B (metafor meta) of the original
# script remain the sensitivity arms and are untouched.
#
# Engine: the shared de_engine_lvqw.R helper (fit_lvqw + build_design_guarded),
# validated here on real NAFL-vs-NASH data.
#
# Design (sex = NUISANCE covariate, matching 05h canonical pattern):
#   ~ dataset + inferred_sex + diagnosis_harmonized
# where diagnosis_harmonized is collapsed to the NAFL-vs-NASH grouping
# (Borderline grouped with NASH per the existing Section C logic). coef tested
# = the NASH level of that factor.
#
# COHORT SET: NAFL-vs-NASH is a DISEASE-INTERNAL contrast (both arms are
# disease), so it does NOT need healthy controls and is NOT restricted to
# include_in_mega (that filter is only for the Disease-vs-Control canonical).
# This replicates the ORIGINAL Section C sample selection EXACTLY
# (13_nafl_vs_nash_de.R lines 30-44): meta_matched.rds filtered to
# pass_technical, PRJNA512027 dropped, all NAFL/NASH/Borderline-labeled
# cohorts kept (~743 samples, 7 cohorts) -> a true apples-to-apples engine swap
# (dream -> fit_lvqw) vs nafl_vs_nash_dream.csv.
#
# Input mechanics follow 05h: subset the canonical merged_dge.rds (already
# TMM-normalized + filterByExpr'd; 27,638-gene universe) -> NO re-TMM, NO
# re-filter. Sex + diagnosis pulled from meta_matched.rds.
#
# STAGING ONLY. Writes:
#   results/disease_signatures/nafl_vs_nash_lvqw.csv         (all genes)
#   results/disease_signatures/nafl_vs_nash_lvqw_sanity.csv  (gate vs dream)
# Backs up the canonical dream file to *_premigration_backup.csv (copy only).
# Does NOT overwrite any canonical file.
# ===========================================================================
suppressMessages({
  library(edgeR)
  library(limma)
  library(data.table)
})

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
BASE <- file.path(ROOT, "RNA-seq/Human/Patient_Cohorts")
INT  <- file.path(BASE, "analysis/integration")
RDIR_INT <- file.path(INT, "results/integration")
RDIR_SIG <- file.path(INT, "results/disease_signatures")
SCRIPTS  <- file.path(INT, "scripts")

source(file.path(SCRIPTS, "de_engine_lvqw.R"))

# ---------------------------------------------------------------------------
# 1) Input loading. Sample selection REPLICATES original Section C exactly
#    (13_nafl_vs_nash_de.R lines 30-44); mechanics follow 05h (subset the
#    canonical merged_dge.rds; NO re-TMM / re-filter).
# ---------------------------------------------------------------------------
dge  <- readRDS(file.path(RDIR_INT, "merged_dge.rds"))
meta <- readRDS(file.path(RDIR_INT, "meta_matched.rds"))
qc   <- fread(file.path(INT, "qc/sample_qc_report.csv"))

# --- ORIGINAL Section C selection ---
meta <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]
meta <- meta[dataset != "PRJNA512027"]   # permanently removed 2026-05-15
meta_nn <- meta[diagnosis_harmonized %in% c("NAFL", "NASH", "Borderline")]
# Borderline (NAS 3-4) grouped with NASH, exactly as original Section C
meta_nn[, nafl_nash := fifelse(diagnosis_harmonized == "NAFL", "NAFL", "NASH")]

# Keep only samples present in the canonical merged_dge.rds (all should be)
sel_ids <- intersect(meta_nn$sample_id, colnames(dge))
stopifnot(length(sel_ids) == nrow(meta_nn))   # no original-C sample dropped
keep_cols <- colnames(dge) %in% sel_ids
dge_nn <- dge[, keep_cols]                     # subset only — no re-TMM / re-filter

# Align sample-level covariates to dge column order
mm <- meta_nn[match(colnames(dge_nn), meta_nn$sample_id)]
stopifnot(all(mm$sample_id == colnames(dge_nn)))

info <- data.frame(
  diagnosis_harmonized = factor(mm$nafl_nash, levels = c("NAFL", "NASH")),
  dataset              = factor(as.character(mm$dataset)),
  inferred_sex         = factor(as.character(mm$inferred_sex))
)
rownames(info) <- colnames(dge_nn)

n_nafl <- sum(info$diagnosis_harmonized == "NAFL")
n_nash <- sum(info$diagnosis_harmonized == "NASH")
cat(sprintf("[in] %d genes x %d samples | NAFL=%d NASH=%d (Borderline->NASH) | %d cohorts: %s\n",
            nrow(dge_nn), ncol(dge_nn), n_nafl, n_nash,
            nlevels(droplevels(info$dataset)),
            paste(levels(droplevels(info$dataset)), collapse = ",")))
cat("[in] sex levels:", paste(levels(droplevels(info$inferred_sex)), collapse = ","), "\n")
cat("[in] per-cohort NAFL/NASH breakdown:\n")
print(table(info$dataset, info$diagnosis_harmonized))

# ---------------------------------------------------------------------------
# 2) Design via build_design_guarded; sex is a NUISANCE covariate.
#    ~ dataset + inferred_sex + diagnosis_harmonized
#    diagnosis_harmonized is LAST so it is NOT dropped by the guard; verify coef.
# ---------------------------------------------------------------------------
dg <- build_design_guarded(
  info,
  rhs_terms = c("dataset", "inferred_sex", "diagnosis_harmonized")
)
design <- dg$design
cat("[design] formula:", dg$formula_used, "\n")
if (length(dg$dropped)) cat("[design] dropped:", paste(dg$dropped, collapse = ", "), "\n")
cat("[design] columns:", paste(colnames(design), collapse = ", "), "\n")

coef_name <- "diagnosis_harmonizedNASH"
stopifnot(coef_name %in% colnames(design))
cat("[design] coef tested:", coef_name, "\n")

# ---------------------------------------------------------------------------
# 3) Fit the LVQW engine (validates the shared helper on real data).
# ---------------------------------------------------------------------------
cat("[fit] running fit_lvqw() ...\n")
res <- fit_lvqw(dge_nn, design, coef = coef_name, do_ashr = TRUE)
res[, method := "limma_voom_qw_C2"]

# Schema: gene, logFC, SE, t, P.Value, padj, shrunk_logFC, lfsr, AveExpr, method
setcolorder(res, c("gene", "logFC", "SE", "t", "P.Value", "padj",
                   "shrunk_logFC", "lfsr", "AveExpr", "method"))

n_deg     <- sum(res$padj < 0.05, na.rm = TRUE)
n_deg_up  <- sum(res$padj < 0.05 & res$logFC > 0, na.rm = TRUE)
n_deg_dn  <- sum(res$padj < 0.05 & res$logFC < 0, na.rm = TRUE)
cat(sprintf("[lvqw] genes=%d | DEGs(padj<0.05)=%d (up=%d down=%d) | ashr non-NA=%d\n",
            nrow(res), n_deg, n_deg_up, n_deg_dn, sum(!is.na(res$lfsr))))

# ---------------------------------------------------------------------------
# 4) STAGING write (NEVER overwrite canonical) + backup of dream file.
# ---------------------------------------------------------------------------
out_lvqw <- file.path(RDIR_SIG, "nafl_vs_nash_lvqw.csv")
fwrite(res, out_lvqw)
cat("[write] staged:", out_lvqw, "\n")

dream_path  <- file.path(RDIR_SIG, "nafl_vs_nash_dream.csv")
backup_path <- file.path(RDIR_SIG, "nafl_vs_nash_dream_premigration_backup.csv")
if (file.exists(dream_path) && !file.exists(backup_path)) {
  ok <- file.copy(dream_path, backup_path, overwrite = FALSE)
  cat(sprintf("[backup] dream -> premigration_backup (copy ok=%s)\n", ok))
} else if (file.exists(backup_path)) {
  cat("[backup] premigration_backup already exists; left as-is\n")
} else {
  cat("[backup] WARNING: dream file not found; no backup made\n")
}

# ---------------------------------------------------------------------------
# 5) SANITY GATE vs canonical dream (nafl_vs_nash_dream.csv) — SAME samples now.
#    Metrics on shared genes: n DEG each, direction concordance among
#    both-significant genes, logFC Spearman rho. Expect dir>=0.95, rho>=0.85.
# ---------------------------------------------------------------------------
dream <- fread(dream_path)
setnames(dream, "adj.P.Val", "dream_padj", skip_absent = TRUE)  # C2-OK-sensitivity
setnames(dream, "logFC", "dream_logFC", skip_absent = TRUE)  # C2-OK-sensitivity

mg <- merge(
  res[, .(gene, lvqw_logFC = logFC, lvqw_padj = padj)],
  dream[, .(gene, dream_logFC, dream_padj)],  # C2-OK-sensitivity
  by = "gene"
)
cat(sprintf("[sanity] shared genes: %d (lvqw=%d, dream=%d)\n",
            nrow(mg), nrow(res), nrow(dream)))

rho_all <- suppressWarnings(cor(mg$lvqw_logFC, mg$dream_logFC,  # C2-OK-sensitivity
                                method = "spearman", use = "complete.obs"))

both_sig <- mg[lvqw_padj < 0.05 & dream_padj < 0.05]  # C2-OK-sensitivity
dir_conc_both <- if (nrow(both_sig) > 0) {
  mean(sign(both_sig$lvqw_logFC) == sign(both_sig$dream_logFC))  # C2-OK-sensitivity
} else NA_real_

lvqw_sig <- mg[lvqw_padj < 0.05]
dir_conc_lvqw_sig <- if (nrow(lvqw_sig) > 0) {
  mean(sign(lvqw_sig$lvqw_logFC) == sign(lvqw_sig$dream_logFC))  # C2-OK-sensitivity
} else NA_real_

dir_conc_all <- mean(sign(mg$lvqw_logFC) == sign(mg$dream_logFC), na.rm = TRUE)  # C2-OK-sensitivity

sanity <- data.table(
  metric = c("n_genes_lvqw", "n_genes_dream", "n_shared_genes",
             "n_deg_lvqw_padj05", "n_deg_dream_padj05", "n_deg_both_sig",
             "direction_concordance_both_sig",
             "direction_concordance_lvqw_sig",
             "direction_concordance_all_shared",
             "logFC_spearman_rho_all_shared",
             "n_cohorts", "n_nafl", "n_nash", "n_samples"),
  value = c(nrow(res), nrow(dream), nrow(mg),
            n_deg, sum(dream$dream_padj < 0.05, na.rm = TRUE), nrow(both_sig),  # C2-OK-sensitivity
            round(dir_conc_both, 4),
            round(dir_conc_lvqw_sig, 4),
            round(dir_conc_all, 4),
            round(rho_all, 4),
            nlevels(droplevels(info$dataset)),
            n_nafl, n_nash, ncol(dge_nn))
)
sanity[, note := ""]
sanity[metric == "n_genes_lvqw", note := coef_name]

out_sanity <- file.path(RDIR_SIG, "nafl_vs_nash_lvqw_sanity.csv")
fwrite(sanity, out_sanity)
cat("[write] staged sanity:", out_sanity, "\n")

cat("\n========== SANITY GATE ==========\n")
print(sanity)
cat(sprintf("\n[gate] direction(both-sig)=%.4f (expect >=0.95) | logFC rho=%.4f (expect >=0.85)\n",
            dir_conc_both, rho_all))
pass_dir <- !is.na(dir_conc_both) && dir_conc_both >= 0.95
pass_rho <- !is.na(rho_all) && rho_all >= 0.85
cat(sprintf("[gate] PASS direction=%s | PASS rho=%s\n", pass_dir, pass_rho))

cat("\n=== 13_nafl_vs_nash_lvqw.R complete (STAGING ONLY; canonical untouched) ===\n")
