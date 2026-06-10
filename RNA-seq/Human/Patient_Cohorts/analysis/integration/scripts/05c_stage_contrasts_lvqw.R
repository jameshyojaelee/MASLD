#!/usr/bin/env Rscript
# ===========================================================================
# 05c_stage_contrasts_lvqw.R
# ---------------------------------------------------------------------------
# HARMONIZED coarse disease-stage contrasts vs Healthy, on the canonical
# limma-voom-quality-weighted (LVQW) engine (de_engine_lvqw.R). This replaces
# the dream() engine in 05c_dream_stage_contrasts.R while PRESERVING 05c's
# disease_stage_coarse derivation for the bulk samples.
#
# These three coarse contrasts (Steatosis / Steatohepatitis / Cirrhosis vs
# Healthy) are the WS3 cell-type-routing coarse anchors, so correctness is
# load-bearing.
#
# What changed vs 05c_dream_stage_contrasts.R
# -------------------------------------------
#   * ENGINE: dream() (mixed model, per-contrast binary fit) -> LVQW
#     (voomWithQualityWeights -> lmFit -> eBayes -> ashr), via fit_lvqw().
#   * DESIGN: ONE multi-level factor design with Healthy as reference, fit
#     once; one coef extracted per non-Healthy stage. Method-specific design:
#         ~ dataset + inferred_sex + disease_stage_coarse
#     (sex is a NUISANCE covariate here, not the term of interest). dataset is
#     a fixed effect (replaces dream's (1|dataset_subbatch) random intercept).
#   * INPUT: merged_dge.rds subset to yaml include_in_mega cohorts (the 05h
#     canonical pattern) -- NO re-TMM, NO re-filterByExpr. The stage-defining
#     fields (diagnosis_harmonized / fibrosis_stage / nas_score / condition)
#     are joined from meta_matched.rds by sample_id.
#
# What is PRESERVED from 05c
# --------------------------
#   * Stage definitions verbatim (steatosis_target / sh_target /
#     cirrhosis_target / is_healthy), including the cirrhosis-coded condition
#     regex. Only difference: the single-factor design needs each sample in
#     exactly ONE level, so the 3 NAFL-low-NAS-but-F4 samples that 05c's
#     separate binary fits placed in BOTH steatosis and cirrhosis are resolved
#     to Cirrhosis here (fibrosis_stage==4 is 05c's own "highest-fidelity"
#     cirrhosis signal). Samples matching NO coarse stage (e.g. GSE213621
#     fibrosis-only disease with no NAS/diagnosis) are dropped from the
#     contrast (NA stage), exactly as a coarse anchor requires.
#
# STAGING ONLY -- does NOT overwrite the canonical dream stage files. Writes:
#   results/integration/stage_steatosis_lvqw.csv
#   results/integration/stage_sh_lvqw.csv
#   results/integration/stage_cirrhosis_lvqw.csv
#   results/integration/stage_lvqw_sanity.csv
# and backs up dream_results_stage_*.csv -> *_premigration_backup.csv (copy).
#
# Usage (SLURM, never login node):
#   env -u SLURM_JOB_ID sbatch --partition=io --qos=interactive --mem=24G \
#     --cpus-per-task=2 --time=48:00:00 --job-name=limma \
#     --output=<scripts>/logs/05c_lvqw_%j.log \
#     --wrap="<rnaseq>/bin/Rscript 05c_stage_contrasts_lvqw.R"
# ===========================================================================

set.seed(42)
suppressMessages({
  library(edgeR); library(limma); library(ashr); library(data.table); library(yaml)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR  <- file.path(INT, "results/integration")
SCRDIR<- file.path(INT, "scripts")
stopifnot(dir.exists(RDIR))

# Source the canonical shared engine (library-only: no side effects on source).
source(file.path(SCRDIR, "de_engine_lvqw.R"))

METHOD_TAG <- "limma_voom_qw_C2"

cat("=== 05c LVQW: coarse disease-stage contrasts vs Healthy ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ---------------------------------------------------------------------------
# 1. Input loading -- 05h canonical pattern (subset merged_dge.rds; no re-TMM,
#    no re-filterByExpr). Stage fields joined from meta_matched.rds.
# ---------------------------------------------------------------------------
dge  <- readRDS(file.path(RDIR, "merged_dge.rds"))
ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
keep <- dge$samples$dataset %in% mega
dge_mega <- dge[, keep]
cat(sprintf("[in] merged_dge subset to include_in_mega: %d genes x %d samples; cohorts: %s\n",
            nrow(dge_mega), ncol(dge_mega),
            paste(sort(unique(dge_mega$samples$dataset)), collapse = ", ")))

meta <- readRDS(file.path(RDIR, "meta_matched.rds"))
m <- meta[match(colnames(dge_mega), meta$sample_id)]
stopifnot(!anyNA(m$sample_id))           # every dge column must match a meta row
stopifnot(identical(m$sample_id, colnames(dge_mega)))

# ---------------------------------------------------------------------------
# 2. Stage derivation -- PRESERVED VERBATIM from 05c_dream_stage_contrasts.R.
# ---------------------------------------------------------------------------
m[, steatosis_target := (
  diagnosis_harmonized == "NAFL" |
  (is.na(diagnosis_harmonized) & !is.na(nas_score) & nas_score < 3 &
   group_binary == "Disease"))]

m[, sh_target := (
  (diagnosis_harmonized %in% c("NASH", "Borderline") |
   (is.na(diagnosis_harmonized) & !is.na(nas_score) & nas_score >= 3 &
    group_binary == "Disease")) &
  (is.na(fibrosis_stage) | fibrosis_stage < 4))]

cirrhosis_labels <- unique(m$condition[
  grepl("cirrho|Cirrho|CIRRHO|F4", m$condition, ignore.case = FALSE)])
cat("Condition labels detected as cirrhosis-coded:",
    paste(cirrhosis_labels, collapse = ", "), "\n")
m[, cirrhosis_target := (
  (!is.na(fibrosis_stage) & fibrosis_stage == 4) |
  (condition %in% cirrhosis_labels & !is.na(condition)))]

m[, is_healthy := (diagnosis_harmonized == "Control")]

for (cc in c("steatosis_target", "sh_target", "cirrhosis_target", "is_healthy")) {
  m[which(is.na(get(cc))), (cc) := FALSE]
}

# Collapse the four 05c boolean flags into ONE coarse factor for the
# single-design fit. Assignment is by precedence so each sample lands in
# exactly one level: Cirrhosis > Steatohepatitis > Steatosis > Healthy. The
# only flags that ever co-fire are NAFL-low-NAS samples with fibrosis_stage==4
# (3 samples); fibrosis_stage==4 is 05c's stated highest-fidelity cirrhosis
# signal, so Cirrhosis wins. Disease samples that match no coarse stage stay NA.
m[, disease_stage_coarse := NA_character_]
m[which(m$is_healthy),       disease_stage_coarse := "Healthy"]
m[which(m$steatosis_target), disease_stage_coarse := "Steatosis"]
m[which(m$sh_target),        disease_stage_coarse := "Steatohepatitis"]
m[which(m$cirrhosis_target), disease_stage_coarse := "Cirrhosis"]

cat("\nCoarse stage assignment (cirrhosis precedence over the 3 NAFL-F4 overlaps):\n")
print(table(m$disease_stage_coarse, useNA = "always"))
cat("\nCohort x stage:\n")
print(table(m$dataset, m$disease_stage_coarse, useNA = "ifany"))

# ---------------------------------------------------------------------------
# 3. Restrict to staged + Healthy samples and fit ONE LVQW design.
# ---------------------------------------------------------------------------
STAGE_LEVELS <- c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis")
in_contrast <- !is.na(m$disease_stage_coarse)
m_fit  <- m[in_contrast]
dge_fit <- dge_mega[, in_contrast]
stopifnot(identical(m_fit$sample_id, colnames(dge_fit)))
cat(sprintf("\n[fit] %d samples enter the staged contrast (%d dropped as NA-stage)\n",
            ncol(dge_fit), sum(!in_contrast)))

# Healthy is the reference level (first).
info <- data.frame(
  dataset             = factor(dge_fit$samples$dataset),
  inferred_sex        = factor(m_fit$inferred_sex),
  disease_stage_coarse = factor(m_fit$disease_stage_coarse, levels = STAGE_LEVELS),
  row.names           = colnames(dge_fit),
  check.names         = FALSE)
cat("\ndisease_stage_coarse reference level:", levels(info$disease_stage_coarse)[1], "\n")

# Build the design via the guarded helper (drops rank-deficient / single-level
# terms). Stage is the term of interest; dataset + inferred_sex are nuisance.
gd <- build_design_guarded(info,
        rhs_terms = c("dataset", "inferred_sex", "disease_stage_coarse"))
design <- gd$design
cat("\n[build_design_guarded] formula_used:", gd$formula_used, "\n")
cat("[build_design_guarded] dropped:",
    if (length(gd$dropped)) paste(gd$dropped, collapse = ", ") else "<none>", "\n")
cat("[build_design_guarded] design columns:\n")
print(colnames(design))

# Verify all three stage coefs survived design construction (caller contract).
COEFS <- c(
  Steatosis       = "disease_stage_coarseSteatosis",
  Steatohepatitis = "disease_stage_coarseSteatohepatitis",
  Cirrhosis       = "disease_stage_coarseCirrhosis")
for (cf in COEFS) {
  stopifnot(cf %in% colnames(design))
}
cat("\n[ok] all 3 stage coefficients present in design:\n  ",
    paste(COEFS, collapse = "\n   "), "\n")

# ---------------------------------------------------------------------------
# 4. Fit LVQW once per stage coef. (Design + voom weights are identical across
#    coefs; fit_lvqw() re-runs voom/eBayes per call, which is cheap here and
#    keeps each call self-contained and identical to the canonical engine.)
# ---------------------------------------------------------------------------
sym_gm <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))
sym_gm[, eb := sub("[.][0-9]+$", "", gene_id)]

OUT_FILES <- c(
  Steatosis       = "stage_steatosis_lvqw.csv",
  Steatohepatitis = "stage_sh_lvqw.csv",
  Cirrhosis       = "stage_cirrhosis_lvqw.csv")

results <- list()
for (stage in names(COEFS)) {
  cf <- COEFS[[stage]]
  cat("\n==============================================================\n")
  cat(sprintf("LVQW contrast: %s vs Healthy   (coef = %s)\n", stage, cf))
  cat("==============================================================\n")
  res <- fit_lvqw(dge_fit, design, coef = cf, do_ashr = TRUE)

  # add method tag + symbol (drop-in parity with dream stage files)
  res[, method := METHOD_TAG]
  res[, symbol := sym_gm[match(sub("[.][0-9]+$", "", gene), eb), gene_name]]
  setcolorder(res, c("gene", "symbol", "logFC", "SE", "t", "P.Value", "padj",
                     "shrunk_logFC", "lfsr", "AveExpr", "method"))

  sig    <- sum(res$padj < 0.05, na.rm = TRUE)
  sig_up <- sum(res$padj < 0.05 & res$logFC > 0, na.rm = TRUE)
  sig_dn <- sum(res$padj < 0.05 & res$logFC < 0, na.rm = TRUE)
  cat(sprintf("genes tested: %d | DEG padj<.05: %d (up %d / dn %d)\n",
              nrow(res), sig, sig_up, sig_dn))

  fwrite(res, file.path(RDIR, OUT_FILES[[stage]]))
  cat("Saved:", file.path(RDIR, OUT_FILES[[stage]]), "\n")
  results[[stage]] <- res
}

# ---------------------------------------------------------------------------
# 5. Back up the existing dream stage files (copy, leave originals in place).
# ---------------------------------------------------------------------------
DREAM_FILES <- c(
  Steatosis       = "dream_results_stage_steatosis.csv",
  Steatohepatitis = "dream_results_stage_sh.csv",
  Cirrhosis       = "dream_results_stage_cirrhosis.csv")
for (stage in names(DREAM_FILES)) {
  src <- file.path(RDIR, DREAM_FILES[[stage]])
  if (file.exists(src)) {
    dst <- file.path(RDIR, sub("\\.csv$", "_premigration_backup.csv",
                               DREAM_FILES[[stage]]))
    ok <- file.copy(src, dst, overwrite = TRUE)
    cat(sprintf("[backup] %s -> %s  (%s)\n", basename(src), basename(dst),
                if (ok) "ok" else "FAILED"))
  } else {
    cat(sprintf("[backup] WARN: dream file missing, no backup: %s\n", src))
  }
}

# ---------------------------------------------------------------------------
# 6. SANITY GATE vs the existing dream stage files.
#    Per contrast: n DEG (padj<.05) for both engines, direction concordance,
#    and logFC Spearman rho on shared genes.
# ---------------------------------------------------------------------------
cat("\n=== SANITY GATE vs dream stage files ===\n")
san_rows <- list()
for (stage in names(COEFS)) {
  lv <- results[[stage]]
  dr <- fread(file.path(RDIR, DREAM_FILES[[stage]]))
  mm <- merge(lv[, .(gene, lv_lfc = logFC, lv_padj = padj)],
              dr[, .(gene, dr_lfc = logFC, dr_padj = padj)], by = "gene")
  rho  <- cor(mm$lv_lfc, mm$dr_lfc, method = "spearman", use = "complete.obs")
  dirA <- mean(sign(mm$lv_lfc) == sign(mm$dr_lfc), na.rm = TRUE)
  n_lv <- sum(lv$padj < 0.05, na.rm = TRUE)
  n_dr <- sum(dr$padj < 0.05, na.rm = TRUE)
  san_rows[[stage]] <- data.table(
    contrast        = stage,
    coef            = COEFS[[stage]],
    n_genes_lvqw    = nrow(lv),
    n_genes_dream   = nrow(dr),
    n_shared        = nrow(mm),
    nDEG_lvqw       = n_lv,
    nDEG_dream      = n_dr,
    direction_concordance = round(dirA, 4),
    logFC_spearman  = round(rho, 4),
    pass_direction  = dirA >= 0.95,
    pass_rho        = rho  >= 0.85)
  cat(sprintf("[%s] nDEG lvqw=%d dream=%d | shared=%d | dir=%.4f rho=%.4f | %s\n",
              stage, n_lv, n_dr, nrow(mm), dirA, rho,
              if (dirA >= 0.95 && rho >= 0.85) "PASS" else "REVIEW"))
}
san <- rbindlist(san_rows)
fwrite(san, file.path(RDIR, "stage_lvqw_sanity.csv"))
cat("\nSaved sanity table:", file.path(RDIR, "stage_lvqw_sanity.csv"), "\n")
print(san, row.names = FALSE)

cat("\n=== 05c LVQW complete:", as.character(Sys.time()),
    "(NO canonical file overwritten) ===\n")
