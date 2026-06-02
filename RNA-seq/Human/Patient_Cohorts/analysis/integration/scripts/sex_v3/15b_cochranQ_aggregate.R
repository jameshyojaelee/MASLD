#!/usr/bin/env Rscript
# sex_v3/15b_cochranQ_aggregate.R
# ---------------------------------------------------------------------------
# Pillar 6 aggregator — fold the 5 per-cohort fits into:
#   * per_cohort_forest_data_v6.csv   (gene × cohort: β_int, SE, CI, p_int)
#   * cochranQ_v6.csv                 (per-gene Q, df, p_Q, I²)
#
# Cochran's Q (fixed-effects pooled effect):
#   wᵢ = 1 / SE_i^2 ;  β̂_pool = Σ wᵢ β_i / Σ wᵢ
#   Q  = Σ wᵢ (β_i − β̂_pool)^2 ;  df = k − 1 (k=5 cohorts -> df=4)
#   p_Q = pchisq(Q, df, lower.tail = FALSE)
#   I²  = max(0, (Q − df) / Q) · 100
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table) })

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SEXV3 <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                   "results/integration/sex_v3")

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
UTILS <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                   "scripts/sex_v3/sex_v3_utils.R")
if (file.exists(UTILS)) source(UTILS)
fwrite_x <- if (exists("write_atomic_csv")) write_atomic_csv else fwrite
QDIR <- file.path(SEXV3, "intermediates/cochranQ_v6")
OUT_FOREST <- file.path(SEXV3, "per_cohort_forest_data_v6.csv")
OUT_Q      <- file.path(SEXV3, "cochranQ_v6.csv")

rds_files <- list.files(QDIR, pattern = "^cohort_.*\\.rds$",
                        full.names = TRUE)
cat("Found", length(rds_files), "cochranQ_v6 per-cohort fits\n")
if (length(rds_files) == 0) stop("No per-cohort fits in ", QDIR)

# Long-form: gene × cohort × {β_int, se_int, ci_lo, ci_hi, p_int, n_samp,
#                              status}
long <- rbindlist(lapply(rds_files, function(f) {
  o <- readRDS(f)
  if (is.null(o$loco_dt) || nrow(o$loco_dt) == 0) return(NULL)
  dt <- copy(o$loco_dt)
  dt[, cohort   := o$cohort]
  dt[, idx      := o$idx]
  dt[, n_samp   := o$n_samp]
  dt[, status   := o$status]
  dt
}), fill = TRUE)
cat("Long-form rows:", nrow(long), "\n")
fwrite_x(long, OUT_FOREST)
cat("Wrote:", OUT_FOREST, "  rows:", nrow(long),
    "  cohorts:", length(unique(long$cohort)), "\n")

# Per-gene Cochran's Q (fixed-effects)
ok <- long[!is.na(beta_int) & !is.na(se_int) & se_int > 0]
ok[, w := 1 / (se_int^2)]
pooled <- ok[, {
  if (.N < 2) {
    .(beta_pool = NA_real_, se_pool = NA_real_, k = .N,
      Q = NA_real_, df = NA_integer_,
      p_Q = NA_real_, I2_pct = NA_real_)
  } else {
    wsum <- sum(w, na.rm = TRUE)
    bp   <- sum(w * beta_int, na.rm = TRUE) / wsum
    sep  <- sqrt(1 / wsum)
    Q    <- sum(w * (beta_int - bp)^2, na.rm = TRUE)
    df   <- .N - 1L
    pQ   <- pchisq(Q, df, lower.tail = FALSE)
    I2   <- max(0, (Q - df) / Q) * 100
    .(beta_pool = bp, se_pool = sep, k = .N,
      Q = Q, df = df, p_Q = pQ, I2_pct = I2)
  }
}, by = gene]
pooled[, Q_padj_BH := p.adjust(p_Q, "BH")]
pooled[, cohort_heterogeneous := (!is.na(Q_padj_BH) & Q_padj_BH < 0.05) |
                                  (!is.na(I2_pct) & I2_pct > 75)]
fwrite_x(pooled, OUT_Q)
cat("Wrote:", OUT_Q, "  rows:", nrow(pooled), "\n")

cat("\n== Q summary ==\n")
cat("  median I² across genes:", round(median(pooled$I2_pct, na.rm = TRUE), 1),
    "%\n")
cat("  p_Q < 0.05:", sum(pooled$p_Q < 0.05, na.rm = TRUE), "\n")
cat("  Q_padj_BH < 0.05:", sum(pooled$Q_padj_BH < 0.05, na.rm = TRUE), "\n")
cat("  cohort_heterogeneous:", sum(pooled$cohort_heterogeneous, na.rm = TRUE),
    "/", nrow(pooled), "\n")

# Cross-reference Strong / Moderate v5 tier calls if available
v5_class_file <- file.path(SEXV3, "interaction_classifier_v5.csv")
if (file.exists(v5_class_file)) {
  v5 <- fread(v5_class_file,
              select = c("gene", "gene_symbol",
                         "class_v5_interaction", "evidence_tier_A_B_C"))
  m <- merge(v5, pooled, by = "gene", all.x = TRUE)
  strong <- m[evidence_tier_A_B_C %in% c("A", "B")]
  cat("\n== Strong/Moderate (tier A/B) v5 calls: heterogeneity ==\n")
  cat("  n_strong_moderate:", nrow(strong), "\n")
  cat("  median I² (Strong/Moderate):",
      round(median(strong$I2_pct, na.rm = TRUE), 1), "%\n")
  cat("  cohort_heterogeneous in Strong/Moderate:",
      sum(strong$cohort_heterogeneous, na.rm = TRUE), "/",
      nrow(strong), "\n")
  # Per-call sanity: any Strong call where ALL 5 cohort CIs cross 0?
  strong_genes <- strong$gene
  flag_all_cross_zero <- character(0)
  for (g in strong_genes) {
    sub <- long[gene == g & !is.na(beta_int)]
    if (nrow(sub) >= 1 && all(sub$ci_lo < 0 & sub$ci_hi > 0)) {
      flag_all_cross_zero <- c(flag_all_cross_zero, g)
    }
  }
  cat("  Strong/Moderate calls with ALL cohort CIs crossing 0:",
      length(flag_all_cross_zero), "\n")
  if (length(flag_all_cross_zero) > 0) {
    head_sym <- v5[gene %in% flag_all_cross_zero, gene_symbol][1:min(10,
                  length(flag_all_cross_zero))]
    cat("    (first 10 symbols):", paste(head_sym, collapse = ", "), "\n")
  }
}
