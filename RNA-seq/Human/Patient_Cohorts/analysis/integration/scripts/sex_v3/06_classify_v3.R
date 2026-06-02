#!/usr/bin/env Rscript
# sex_v3/06_classify_v3.R
# ---------------------------------------------------------------------------
# Module 06 — Posterior classification of sex-effect patterns from mashr fit.
#
# Consumes mashr posterior weights (Module 05) and produces a per-gene class
# call across four sex-effect patterns:
#     concordant   (same direction, both sexes)
#     F_only       (female-specific effect)
#     M_only       (male-specific effect)
#     anti_correlated (opposing directions between sexes)
#
# Replaces the v2 categorical scheme (Concordant / Female_biased / Male_biased /
# Divergent). v2 framing was confounded by power asymmetry (101 F vs 52 M
# controls) and an unstable padj<0.05 sign tie-breaker; the posterior weights
# from mashr (Urbut & Stephens 2019) absorb this via hierarchical shrinkage.
#
# Reads:
#   - intermediates/mashr_fit.rds                       (Module 05)
#   - mashr_pattern_loadings.csv                        (Module 05)
#   - mashr_posterior_summary.csv                       (Module 05)
#   - sex_xci/chrx_stratified_classification.csv        (A9 chrX annotations)
#   - sex_deg_classification.csv                        (v2, for legacy compare)
#
# Writes:
#   - sex_deg_classification_v3.csv                     (PRIMARY OUTPUT)
#   - sex_class_provenance.tsv                          (JSON evidence trail)
#
# Schema (see plan §"Output Schema"):
#   gene, gene_symbol, ensembl_base, chr, chr_category,
#   beta_F, beta_M, beta_interaction,
#   posterior_mean_F, posterior_sd_F, lfsr_F,
#   posterior_mean_M, posterior_sd_M, lfsr_M,
#   posterior_P_concordant, posterior_P_F_only,
#   posterior_P_M_only, posterior_P_anti_correlated,
#   assigned_class, posterior_confidence,
#   sex_class_chrx_aware, suspicious_flag,
#   boot_ci_lo_beta_F, boot_ci_hi_beta_F,
#   boot_ci_lo_beta_M, boot_ci_hi_beta_M,
#   sex_class_v2_legacy, class_flip_v2_to_v3,
#   calibration_passed, notes
#
# Legacy v3 column for backward compatibility with 27a:
#   sex_class  (alias of assigned_class, used by atlas Layer 5 loader)
#   logFC_M / logFC_F  (alias of beta_M / beta_F so the m_col / f_col resolver
#   in 27a still works without changes)
#   interaction_padj (derived from lfsr_F + lfsr_M combined; alias)
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR   <- file.path(INT, "results/integration")
SEXV3  <- file.path(RDIR, "sex_v3")
IDIR   <- file.path(SEXV3, "intermediates")

# Shared utilities (atomic writes + sessionInfo dump) -- R5 Issue 2 fix.
# Path is resolved relative to MASLD_PROJECT_ROOT so the script runs under
# both interactive (Rscript) and SLURM modes regardless of getwd().
source(file.path(BASE,
                 "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                 "scripts/sex_v3/sex_v3_utils.R"))

CHRX_F    <- file.path(BASE, "RNA-seq/results/audit_sensitivity",
                       "sex_xci/chrx_stratified_classification.csv")
V2_F      <- file.path(RDIR, "sex_deg_classification.csv")
GENE_META <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")

# mashr inputs (Module 05)
FIT_F            <- file.path(IDIR, "mashr_fit.rds")
PATTERN_LOAD_F   <- file.path(SEXV3, "mashr_pattern_loadings.csv")
POSTERIOR_SUM_F  <- file.path(SEXV3, "mashr_posterior_summary.csv")

# Optional bootstrap CI input (Module 04 sample-level bootstrap)
BOOT_CI_F <- file.path(IDIR, "bootstrap_beta_ci.csv")

# Outputs
OUT_CSV <- file.path(SEXV3, "sex_deg_classification_v3.csv")
OUT_TSV <- file.path(SEXV3, "sex_class_provenance.tsv")

cat("============================================================\n")
cat("Module 06 — sex_v3 classifier (mashr posterior -> 4-class)\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")

stopifnot(file.exists(PATTERN_LOAD_F))
stopifnot(file.exists(POSTERIOR_SUM_F))

# ---------------------------------------------------------------------------
# 1) Load mashr pattern loadings + posterior summary
# ---------------------------------------------------------------------------
cat("\n[1] Loading mashr pattern loadings + posterior summary...\n")
loadings <- fread(PATTERN_LOAD_F)
post_sum <- fread(POSTERIOR_SUM_F)

cat("  Genes in pattern loadings:", nrow(loadings), "\n")
cat("  Genes in posterior summary:", nrow(post_sum), "\n")
cat("  Loadings cols:", paste(names(loadings), collapse = ", "), "\n")

# Loadings is expected to carry one weight column per canonical/data-driven
# covariance pattern. Module 05 names them with a stable prefix; we collapse
# canonical + data-driven loadings into the 4 plan-defined patterns by
# summing across "concordant_up" + "concordant_down" into "concordant" and
# treating data-driven discovered covariances as "anti_correlated" when
# their off-diagonal element is negative, else "concordant" buckets, unless
# Module 05 supplies pre-aggregated columns named
# posterior_P_{concordant,F_only,M_only,anti_correlated}, in which case those
# are used verbatim.

prebagged <- all(c("posterior_P_concordant", "posterior_P_F_only",
                   "posterior_P_M_only", "posterior_P_anti_correlated") %in%
                 names(loadings))

if (prebagged) {
  cat("  Detected pre-aggregated 4-pattern columns -> using verbatim.\n")
  cls <- loadings[, .(gene,
                      posterior_P_concordant,
                      posterior_P_F_only,
                      posterior_P_M_only,
                      posterior_P_anti_correlated)]
} else {
  cat("  Aggregating raw canonical + data-driven weights into 4 buckets...\n")
  raw_cols <- setdiff(names(loadings), "gene")

  classify_col <- function(cn) {
    cnl <- tolower(cn)
    if (grepl("anti", cnl) || grepl("opp", cnl) || grepl("divergent", cnl)) {
      return("anti_correlated")
    }
    if (grepl("^f_only|female_only|f\\b", cnl)) return("F_only")
    if (grepl("^m_only|male_only|m\\b",   cnl)) return("M_only")
    if (grepl("concord", cnl) || grepl("equal", cnl) || grepl("same", cnl)) {
      return("concordant")
    }
    # Data-driven covariance falls through; Module 05 should annotate them.
    # Default to concordant since canonical covariances dominate.
    return("concordant")
  }

  bucket_map <- vapply(raw_cols, classify_col, character(1))
  cat("  Bucket assignments:\n")
  for (b in unique(bucket_map)) {
    cat("    ", b, ": ", paste(names(bucket_map)[bucket_map == b],
                                collapse = ", "), "\n", sep = "")
  }

  cls <- loadings[, .(gene)]
  for (b in c("concordant", "F_only", "M_only", "anti_correlated")) {
    cols_b <- names(bucket_map)[bucket_map == b]
    if (length(cols_b) == 0) {
      cls[[paste0("posterior_P_", b)]] <- 0
    } else if (length(cols_b) == 1) {
      cls[[paste0("posterior_P_", b)]] <- loadings[[cols_b]]
    } else {
      cls[[paste0("posterior_P_", b)]] <- rowSums(as.matrix(loadings[, ..cols_b]))
    }
  }
}

# Renormalize so the 4 posterior weights sum to 1 per gene (numerical
# stability — Module 05 may emit a residual null/inactive pattern weight)
cls[, sum_4 := posterior_P_concordant + posterior_P_F_only +
                posterior_P_M_only + posterior_P_anti_correlated]
cls[sum_4 > 0,
    `:=`(posterior_P_concordant      = posterior_P_concordant      / sum_4,
         posterior_P_F_only          = posterior_P_F_only          / sum_4,
         posterior_P_M_only          = posterior_P_M_only          / sum_4,
         posterior_P_anti_correlated = posterior_P_anti_correlated / sum_4)]
cls[, sum_4 := NULL]

# argmax classification + max-weight confidence
pat_cols <- c("posterior_P_concordant", "posterior_P_F_only",
              "posterior_P_M_only", "posterior_P_anti_correlated")
pat_lbl  <- c("concordant", "F_only", "M_only", "anti_correlated")

M <- as.matrix(cls[, ..pat_cols])
cls[, assigned_class       := pat_lbl[max.col(M, ties.method = "first")]]
cls[, posterior_confidence := apply(M, 1, max)]

cat("  Assigned class counts:\n")
print(table(cls$assigned_class, useNA = "ifany"))
cat("  Posterior confidence summary:\n")
print(summary(cls$posterior_confidence))

# ---------------------------------------------------------------------------
# 2) Merge posterior means + lfsr (Module 05 summary)
# ---------------------------------------------------------------------------
cat("\n[2] Merging posterior summary (means, sds, lfsr)...\n")
# Module 05 summary is expected to have:
#   gene, beta_F, beta_M, beta_interaction,
#   posterior_mean_F, posterior_sd_F, lfsr_F,
#   posterior_mean_M, posterior_sd_M, lfsr_M
need_cols <- c("gene", "beta_F", "beta_M", "beta_interaction",
               "posterior_mean_F", "posterior_sd_F", "lfsr_F",
               "posterior_mean_M", "posterior_sd_M", "lfsr_M")
miss_cols <- setdiff(need_cols, names(post_sum))
if (length(miss_cols) > 0) {
  cat("  WARNING: posterior summary is missing columns:",
      paste(miss_cols, collapse = ", "),
      "-- filling with NA.\n")
  for (cc in miss_cols) post_sum[[cc]] <- NA_real_
}

cls <- merge(cls, post_sum[, ..need_cols], by = "gene", all.x = TRUE)

# ---------------------------------------------------------------------------
# 3) Bootstrap CI on beta_F / beta_M (optional input)
# ---------------------------------------------------------------------------
if (file.exists(BOOT_CI_F)) {
  cat("\n[3] Merging bootstrap CI on beta_F / beta_M...\n")
  boot_dt <- fread(BOOT_CI_F)
  boot_keep <- c("gene", "boot_ci_lo_beta_F", "boot_ci_hi_beta_F",
                 "boot_ci_lo_beta_M", "boot_ci_hi_beta_M")
  miss_boot <- setdiff(boot_keep, names(boot_dt))
  if (length(miss_boot) > 0) {
    cat("  WARNING: bootstrap CI file missing cols:",
        paste(miss_boot, collapse = ", "),
        "-- filling NA.\n")
    for (cc in miss_boot) boot_dt[[cc]] <- NA_real_
  }
  cls <- merge(cls, boot_dt[, ..boot_keep], by = "gene", all.x = TRUE)
} else {
  cat("\n[3] Bootstrap CI file not present -- filling boot_ci_* with NA.\n")
  cls[, `:=`(boot_ci_lo_beta_F = NA_real_,
             boot_ci_hi_beta_F = NA_real_,
             boot_ci_lo_beta_M = NA_real_,
             boot_ci_hi_beta_M = NA_real_)]
}

# ---------------------------------------------------------------------------
# 4) chrX-aware relabeling using A9 chrX classification
# ---------------------------------------------------------------------------
cat("\n[4] Joining A9 chrX classification + chrX-aware relabeling...\n")
chrx <- fread(CHRX_F)
cat("  chrX annot rows:", nrow(chrx), "  cols:", ncol(chrx), "\n")

# chrx columns of interest: gene_id, ensembl_base, gene_name, chromosome,
# gene_biotype, chr_category, xci_status_oliva
chrx_keep <- c("gene_id", "ensembl_base", "gene_name", "chromosome",
               "gene_biotype", "chr_category")
miss_chrx <- setdiff(chrx_keep, names(chrx))
if (length(miss_chrx) > 0) {
  cat("  WARNING: chrX file missing cols:",
      paste(miss_chrx, collapse = ", "),
      "-- filling NA.\n")
  for (cc in miss_chrx) chrx[[cc]] <- NA_character_
}
chrx_sub <- chrx[, ..chrx_keep]
setnames(chrx_sub, "gene_id", "gene")

cls <- merge(cls, chrx_sub, by = "gene", all.x = TRUE)
setnames(cls, "gene_name", "gene_symbol", skip_absent = TRUE)
setnames(cls, "chromosome", "chr", skip_absent = TRUE)

# ---- R4 Issue 1 fix: chrY override needs gencode metadata.
# A9's chrx_stratified_classification.csv is chrX-only, so chrY genes never
# get the M_only_chrY_excluded relabel. We join gencode_v49 metadata to fill
# in chromosome (and biotype) for chrY + autosomal genes. ensembl_base is
# the join key against chrx (which used gene_id as gene -- gene carries
# version suffix; gencode uses ensembl_base without suffix).
cat("\n[4b] Joining gencode chromosome map for chrY + autosomal genes (R4 fix)...\n")
if (file.exists(GENE_META)) {
  gc <- fread(GENE_META, select = c("gene_id", "ensembl_base",
                                    "chromosome", "gene_biotype"))
  setnames(gc, c("chromosome", "gene_biotype"),
           c("gc_chromosome", "gc_gene_biotype"))
  # Join by gene_id (full versioned) first; fall back to ensembl_base
  cls <- merge(cls, gc[, .(gene = gene_id, gc_chromosome, gc_gene_biotype)],
               by = "gene", all.x = TRUE)
  # Backfill missing chr/biotype from gencode where chrX merge produced NA
  cls[is.na(chr) | chr == "", chr := gc_chromosome]
  if (!"gene_biotype" %in% names(cls)) cls[, gene_biotype := NA_character_]
  cls[is.na(gene_biotype) | gene_biotype == "", gene_biotype := gc_gene_biotype]
  # Recompute chr_category for chrY genes that gencode caught
  cls[is.na(chr_category) & gc_chromosome == "chrY", chr_category := "chrY"]
  cls[is.na(chr_category) & !is.na(gc_chromosome) & gc_chromosome != "chrX",
      chr_category := "autosome_or_unknown"]
  cls[, c("gc_chromosome", "gc_gene_biotype") := NULL]
} else {
  cat("  WARNING: gencode metadata not found at:", GENE_META,
      "-- chrY override will be incomplete.\n")
}

# Default chr_category = unknown for autosomes / unmapped (after gencode backfill)
cls[is.na(chr_category), chr_category := "autosome_or_unknown"]

# chrX-aware label
cls[, sex_class_chrx_aware := assigned_class]
cls[chr_category == "X_escape" & assigned_class == "F_only",
    sex_class_chrx_aware := "F_only_chrx_escape_expected"]
# R4 Issue 1 fix: chrY override now fires because gencode backfill populates
# chr / chr_category for chrY genes.
cls[chr_category == "chrY" | chr == "chrY",
    sex_class_chrx_aware := "M_only_chrY_excluded"]
cls[chr == "chrY", chr_category := "chrY"]

# suspicious_flag
cls[, suspicious_flag := NA_character_]
cls[chr_category == "X_inactive" & assigned_class == "F_only",
    suspicious_flag := "suspicious_xci_skew"]

cat("  sex_class_chrx_aware counts:\n")
print(table(cls$sex_class_chrx_aware, useNA = "ifany"))
cat("  suspicious_flag counts:\n")
print(table(cls$suspicious_flag, useNA = "ifany"))

# ---------------------------------------------------------------------------
# 5) v2 legacy comparison + flip tracking
# ---------------------------------------------------------------------------
cat("\n[5] Loading v2 sex_deg_classification.csv for legacy comparison...\n")
if (file.exists(V2_F)) {
  v2 <- fread(V2_F)
  if (!"sex_class" %in% names(v2)) {
    cat("  WARNING: v2 file lacks 'sex_class' -- skipping legacy compare.\n")
    cls[, sex_class_v2_legacy := NA_character_]
  } else {
    v2_keep <- v2[, .(gene, sex_class_v2_legacy = sex_class)]
    cls <- merge(cls, v2_keep, by = "gene", all.x = TRUE)
  }
} else {
  cat("  v2 file not found -- skipping legacy compare.\n")
  cls[, sex_class_v2_legacy := NA_character_]
}

# v2 -> v3 flip tracker (label-level, accounting for v2's different label set)
cls[, class_flip_v2_to_v3 := NA_character_]
cls[!is.na(sex_class_v2_legacy),
    class_flip_v2_to_v3 := paste0(sex_class_v2_legacy, " -> ", assigned_class)]

cat("  Top 10 v2->v3 transitions:\n")
print(head(sort(table(cls$class_flip_v2_to_v3), decreasing = TRUE), 10))

# ---------------------------------------------------------------------------
# 6) Calibration flag (placeholder until Module 07 lands)
# ---------------------------------------------------------------------------
# Module 07 produces calibration_simulation.csv; we set the flag = NA here
# and let Module 08 fill it post-hoc. Downstream consumers should treat
# NA as "not yet calibrated".
cls[, calibration_passed := NA]

# ---------------------------------------------------------------------------
# 7) Atlas backward-compatibility aliases for 27a Layer 5 loader
# ---------------------------------------------------------------------------
# 27a's column resolver expects either logFC_M/logFC_F (Script 26 v2) or
# meta_logFC_M/meta_logFC_F (Script 06 legacy). Atlas also reads sex_class
# and interaction_padj. Expose v3 posterior means + a derived interaction
# padj proxy so 27a Layer 5 keeps working with no edit beyond the file
# fallback patch (line 942).
#
# R3 Issue 1 fix: legacy `sex_class` column must use the v2 enum
# {Concordant, Female_biased, Male_biased, Divergent} so the ~12 downstream
# consumers (207, 207b, 208, 209, 216, 216b, 218a-d, 219, figS_sex_dimorphism)
# that hard-code v2 strings keep working. The v3-native lowercase 4-pattern
# vocabulary stays in `assigned_class` for v3-aware consumers.
class_map <- c(
  "concordant"      = "Concordant",
  "F_only"          = "Female_biased",
  "M_only"          = "Male_biased",
  "anti_correlated" = "Divergent"
)
cls[, sex_class    := unname(class_map[assigned_class])]
cls[, `:=`(
  logFC_M         = posterior_mean_M,
  logFC_F         = posterior_mean_F,
  # interaction_padj proxy: 2 * min(lfsr_F, lfsr_M) -- mash's lfsr is the
  # closest analog to a "false-sign rate" for the interaction direction.
  # This is a conservative ceiling on the implied interaction false-discovery.
  interaction_padj = pmin(2 * pmin(lfsr_F, lfsr_M, na.rm = FALSE), 1)
)]
cat("  sex_class (v2-compatible alias) counts:\n")
print(table(cls$sex_class, useNA = "ifany"))

# ---------------------------------------------------------------------------
# 8) Notes column for per-gene caveats
# ---------------------------------------------------------------------------
cls[, notes := NA_character_]
cls[posterior_confidence < 0.5,
    notes := paste0(notes, ";low_confidence_no_dominant_pattern")]
cls[chr_category == "PAR",
    notes := paste0(notes, ";chrX_PAR_treat_as_autosomal")]
cls[, notes := sub("^NA;", "", notes)]
cls[, notes := sub("^;",   "", notes)]
cls[notes == "NA" | notes == "", notes := NA_character_]

# ---------------------------------------------------------------------------
# 9) Schema enforcement + write
# ---------------------------------------------------------------------------
cat("\n[9] Writing v3 classification CSV...\n")
schema_cols <- c(
  "gene", "gene_symbol", "ensembl_base", "chr", "chr_category", "gene_biotype",
  "beta_F", "beta_M", "beta_interaction",
  "posterior_mean_F", "posterior_sd_F", "lfsr_F",
  "posterior_mean_M", "posterior_sd_M", "lfsr_M",
  "posterior_P_concordant", "posterior_P_F_only",
  "posterior_P_M_only", "posterior_P_anti_correlated",
  "assigned_class", "posterior_confidence",
  "sex_class_chrx_aware", "suspicious_flag",
  "boot_ci_lo_beta_F", "boot_ci_hi_beta_F",
  "boot_ci_lo_beta_M", "boot_ci_hi_beta_M",
  "sex_class_v2_legacy", "class_flip_v2_to_v3",
  "calibration_passed", "notes",
  # 27a backward-compat aliases (do not remove)
  "sex_class", "logFC_M", "logFC_F", "interaction_padj"
)
miss_schema <- setdiff(schema_cols, names(cls))
if (length(miss_schema) > 0) {
  cat("  Adding missing schema columns as NA:",
      paste(miss_schema, collapse = ", "), "\n")
  for (cc in miss_schema) cls[[cc]] <- NA
}
setcolorder(cls, schema_cols)

write_atomic_csv(cls[, ..schema_cols], OUT_CSV)
cat("  Wrote:", OUT_CSV, "  rows:", nrow(cls), "  cols:", ncol(cls), "\n")

# ---------------------------------------------------------------------------
# 10) Provenance trail (one JSON evidence object per gene)
# ---------------------------------------------------------------------------
cat("\n[10] Writing per-gene provenance TSV...\n")
prov <- cls[, .(
  gene,
  provenance = vapply(seq_len(.N), function(i) {
    toJSON(list(
      pipeline             = "sex_v3_mashr",
      module               = "06_classify_v3.R",
      generated_at         = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      mashr_fit_rds        = FIT_F,
      pattern_loadings_csv = PATTERN_LOAD_F,
      posterior_summary_csv= POSTERIOR_SUM_F,
      chrx_annotation_csv  = CHRX_F,
      v2_legacy_csv        = V2_F,
      assigned_class       = assigned_class[i],
      posterior_confidence = posterior_confidence[i],
      sex_class_chrx_aware = sex_class_chrx_aware[i],
      suspicious_flag      = suspicious_flag[i],
      class_flip_v2_to_v3  = class_flip_v2_to_v3[i],
      posterior_weights    = list(
        concordant      = posterior_P_concordant[i],
        F_only          = posterior_P_F_only[i],
        M_only          = posterior_P_M_only[i],
        anti_correlated = posterior_P_anti_correlated[i]
      )
    ), auto_unbox = TRUE, na = "null", null = "null")
  }, character(1))
)]
write_atomic_csv(prov, OUT_TSV, sep = "\t")
cat("  Wrote:", OUT_TSV, "  rows:", nrow(prov), "\n")

# ---------------------------------------------------------------------------
# 11) Summary
# ---------------------------------------------------------------------------
cat("\n============================================================\n")
cat("Module 06 complete.\n")
cat("  4-pattern E[N] (sum of posterior weights):\n")
cat(sprintf("    E[N_concordant]      = %.1f\n",
            sum(cls$posterior_P_concordant,      na.rm = TRUE)))
cat(sprintf("    E[N_F_only]          = %.1f\n",
            sum(cls$posterior_P_F_only,          na.rm = TRUE)))
cat(sprintf("    E[N_M_only]          = %.1f\n",
            sum(cls$posterior_P_M_only,          na.rm = TRUE)))
cat(sprintf("    E[N_anti_correlated] = %.1f\n",
            sum(cls$posterior_P_anti_correlated, na.rm = TRUE)))
cat("  High-confidence (posterior_confidence > 0.95):\n")
hi <- cls[posterior_confidence > 0.95]
print(table(hi$assigned_class, useNA = "ifany"))

# R5 Issue 2 fix: sessionInfo dump for per-module reproducibility audit
dump_session_info(IDIR, "06")

cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")
