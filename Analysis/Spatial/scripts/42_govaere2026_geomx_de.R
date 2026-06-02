#!/usr/bin/env Rscript
# ============================================================================
# Govaere 2026 Nat Genet GeoMx WTA -- Region-level differential expression
# ----------------------------------------------------------------------------
# Per paper Methods (page 14): "We fit linear mixed-effect models with the
# GeoMxTools Bioconductor R package (v3.0.1)" using lmerTest::lmer with
# patient as the random effect. We replicate that strategy here, fitting
# lmer(expr ~ group + (1|Patient)) per gene on the Q3-normalized log
# expression matrix.
#
# CONTRASTS:
#   1. SH vs PT  (CD68 only)         -- paper Sup Table 10, expect ~354 sig
#   2. SH vs LS  (CD68 only)         -- paper Sup Table 9,  expect ~215 sig
#   3. panCK vs CD68 (all regions)   -- segment-marker comparison
#   4. CD45  vs CD68 (all regions)   -- segment-marker comparison
#
# OUTPUTS:
#   geomx_de_sh_vs_pt.csv
#   geomx_de_sh_vs_ls.csv
#   geomx_de_panck_vs_cd68.csv
#   geomx_de_cd45_vs_cd68.csv
#   geomx_de_summary.tsv
#
# Each contrast CSV columns: gene_symbol, logFC, AveExpr, t_or_z_stat, pval,
# padj_bh, n_segs_group1, n_segs_group2
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lmerTest)
  library(BiocParallel)
})

# ---------------- PATHS ----------------------------------------------------
project_root <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
out_dir <- file.path(project_root, "Analysis/Spatial/results/govaere2026")
expr_path <- file.path(out_dir, "geomx_expression.tsv")
meta_path <- file.path(out_dir, "geomx_metadata.tsv")
stopifnot(file.exists(expr_path), file.exists(meta_path))

n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
cat("[", format(Sys.time(), "%H:%M:%S"), "] BiocParallel workers:", n_workers, "\n")

# ---------------- LOAD INPUTS ----------------------------------------------
cat("[", format(Sys.time(), "%H:%M:%S"), "] Reading expression + metadata\n")
expr_dt <- fread(expr_path, header = TRUE, sep = "\t", check.names = FALSE)
meta    <- fread(meta_path, header = TRUE, sep = "\t", check.names = FALSE)

# expr_dt: first column "segment" = segment ID rownames, remaining cols = genes
stopifnot("segment" %in% colnames(expr_dt))
seg_ids <- expr_dt$segment
expr_mat_raw <- as.matrix(expr_dt[, -1, with = FALSE])
rownames(expr_mat_raw) <- seg_ids
# expr_mat_raw: segments x genes (linear-scale Q3-normalized target counts).
# The Govaere loader (`40_govaere2026_geomx_load.R`) outputs Q3-normalized
# values on the linear scale (GeomxTools `normalize(norm_method="quant",
# desiredQuantile=0.75)` returns multiplicative scaling, NOT log). We apply
# log2(x + 1) here so the per-gene lmer coefficient is interpretable as a
# log2 fold change, matching paper Methods (page 14 / Sup Table 9, 10).
expr_mat <- log2(expr_mat_raw + 1)
gene_ids <- colnames(expr_mat)

cat("Expression matrix dim (segments x genes):",
    paste(dim(expr_mat), collapse = " x "), "\n")
cat("Raw Q3 value range:",
    paste(round(range(expr_mat_raw, na.rm = TRUE), 3), collapse = " - "), "\n")
cat("log2(Q3+1) value range:",
    paste(round(range(expr_mat, na.rm = TRUE), 3), collapse = " - "), "\n")
cat("Metadata rows:", nrow(meta), "\n")

# Align metadata to expression segments
meta <- meta[match(seg_ids, meta$segment), ]
stopifnot(identical(meta$segment, seg_ids))

# ---------------- DE ENGINE ------------------------------------------------
# Fit lmer(expr ~ group + (1|Patient)) per gene.
# Group must be a factor with reference level = group2 (the "denominator"),
# so the coefficient `groupXXX` represents group1 - group2 (i.e. group1 up
# when positive).
fit_lmer_gene <- function(expr_vec, group_fac, patient_vec) {
  # NA-safe: drop NA expression entries (Q3 normalization rarely produces them
  # but be defensive).
  ok <- is.finite(expr_vec)
  if (sum(ok) < length(unique(group_fac)) + 1L) {
    return(c(estimate = NA_real_, t = NA_real_, p = NA_real_,
             AveExpr = NA_real_))
  }
  df_use <- data.frame(
    y      = expr_vec[ok],
    group  = group_fac[ok],
    Patient = patient_vec[ok]
  )
  if (length(unique(df_use$Patient)) < 2L) {
    # lmer needs >=2 patients for a random intercept; fall back to lm.
    fit <- try(lm(y ~ group, data = df_use), silent = TRUE)
    if (inherits(fit, "try-error")) {
      return(c(estimate = NA_real_, t = NA_real_, p = NA_real_,
               AveExpr = mean(df_use$y)))
    }
    co <- summary(fit)$coefficients
    if (nrow(co) < 2L) {
      return(c(estimate = NA_real_, t = NA_real_, p = NA_real_,
               AveExpr = mean(df_use$y)))
    }
    return(c(estimate = co[2, 1], t = co[2, 3], p = co[2, 4],
             AveExpr = mean(df_use$y)))
  }
  fit <- try(suppressMessages(suppressWarnings(
    lmerTest::lmer(y ~ group + (1 | Patient), data = df_use,
                   REML = TRUE,
                   control = lmerControl(check.conv.singular = .makeCC(action = "ignore",
                                                                        tol = 1e-4)))
  )), silent = TRUE)
  if (inherits(fit, "try-error")) {
    # Singular RE fit, fall back to lm
    fit2 <- try(lm(y ~ group, data = df_use), silent = TRUE)
    if (inherits(fit2, "try-error")) {
      return(c(estimate = NA_real_, t = NA_real_, p = NA_real_,
               AveExpr = mean(df_use$y)))
    }
    co <- summary(fit2)$coefficients
    if (nrow(co) < 2L) {
      return(c(estimate = NA_real_, t = NA_real_, p = NA_real_,
               AveExpr = mean(df_use$y)))
    }
    return(c(estimate = co[2, 1], t = co[2, 3], p = co[2, 4],
             AveExpr = mean(df_use$y)))
  }
  co <- summary(fit)$coefficients
  if (nrow(co) < 2L) {
    return(c(estimate = NA_real_, t = NA_real_, p = NA_real_,
             AveExpr = mean(df_use$y)))
  }
  # lmerTest summary columns: Estimate, Std.Error, df, t value, Pr(>|t|)
  return(c(estimate = co[2, 1], t = co[2, "t value"],
           p = co[2, "Pr(>|t|)"], AveExpr = mean(df_use$y)))
}

run_contrast <- function(name, segs_g1, segs_g2, label_g1, label_g2,
                         workers = n_workers) {
  cat("\n========================================================\n")
  cat("[", format(Sys.time(), "%H:%M:%S"), "] CONTRAST:", name, "\n")
  cat("  group1 (", label_g1, "):", length(segs_g1), "segments\n")
  cat("  group2 (", label_g2, "):", length(segs_g2), "segments\n")
  segs_use <- c(segs_g1, segs_g2)
  group_vec <- factor(c(rep(label_g1, length(segs_g1)),
                        rep(label_g2, length(segs_g2))),
                      levels = c(label_g2, label_g1))  # ref = g2 -> coef = g1 - g2
  patient_vec <- meta$Patient[match(segs_use, meta$segment)]
  cat("  Patients (group1):",
      paste(sort(unique(meta$Patient[match(segs_g1, meta$segment)])),
            collapse = ","), "\n")
  cat("  Patients (group2):",
      paste(sort(unique(meta$Patient[match(segs_g2, meta$segment)])),
            collapse = ","), "\n")
  cat("  Total patients in contrast:",
      length(unique(patient_vec)), "\n")

  E <- expr_mat[segs_use, , drop = FALSE]  # segments x genes
  cat("  Submatrix dim:", paste(dim(E), collapse = " x "), "\n")

  # Use BiocParallel for the per-gene loop. MulticoreParam shares memory so we
  # don't pay a serialization cost for E.
  bp <- if (.Platform$OS.type == "unix") {
    MulticoreParam(workers = workers, RNGseed = 42L,
                   progressbar = FALSE, log = FALSE)
  } else {
    SnowParam(workers = workers, RNGseed = 42L, progressbar = FALSE)
  }

  t0 <- Sys.time()
  res_list <- bplapply(seq_len(ncol(E)),
                       function(j) fit_lmer_gene(E[, j], group_vec, patient_vec),
                       BPPARAM = bp)
  t1 <- Sys.time()
  cat("  Fit", ncol(E), "genes in",
      round(as.numeric(difftime(t1, t0, units = "secs")), 1), "sec\n")

  mat <- do.call(rbind, res_list)
  out <- data.frame(
    gene_symbol      = colnames(E),
    logFC            = mat[, "estimate"],
    AveExpr          = mat[, "AveExpr"],
    t_or_z_stat      = mat[, "t"],
    pval             = mat[, "p"],
    padj_bh          = p.adjust(mat[, "p"], method = "BH"),
    n_segs_group1    = length(segs_g1),
    n_segs_group2    = length(segs_g2),
    stringsAsFactors = FALSE
  )
  # Sort by p ascending
  out <- out[order(out$pval, na.last = TRUE), ]
  return(out)
}

# ---------------- DEFINE CONTRASTS -----------------------------------------
seg_sh_cd68 <- meta$segment[meta$segment_marker == "CD68" &
                              meta$region_type == "steatohepatitis"]
seg_pt_cd68 <- meta$segment[meta$segment_marker == "CD68" &
                              meta$region_type == "portal_tract"]
seg_ls_cd68 <- meta$segment[meta$segment_marker == "CD68" &
                              meta$region_type == "low_steatosis"]
seg_panck   <- meta$segment[meta$segment_marker == "PANCK"]
seg_cd68    <- meta$segment[meta$segment_marker == "CD68"]
seg_cd45    <- meta$segment[meta$segment_marker == "CD45"]

cat("\nSegment availability:\n")
cat("  SH x CD68:", length(seg_sh_cd68),
    " PT x CD68:", length(seg_pt_cd68),
    " LS x CD68:", length(seg_ls_cd68), "\n")
cat("  panCK total:", length(seg_panck),
    " CD68 total:", length(seg_cd68),
    " CD45 total:", length(seg_cd45), "\n")

# ---------------- RUN CONTRASTS --------------------------------------------
res_sh_pt <- run_contrast(
  name = "SH vs PT (CD68 only)",
  segs_g1 = seg_sh_cd68, segs_g2 = seg_pt_cd68,
  label_g1 = "SH", label_g2 = "PT"
)
fwrite(res_sh_pt, file.path(out_dir, "geomx_de_sh_vs_pt.csv"))

res_sh_ls <- run_contrast(
  name = "SH vs LS (CD68 only)",
  segs_g1 = seg_sh_cd68, segs_g2 = seg_ls_cd68,
  label_g1 = "SH", label_g2 = "LS"
)
fwrite(res_sh_ls, file.path(out_dir, "geomx_de_sh_vs_ls.csv"))

res_panck_cd68 <- run_contrast(
  name = "panCK vs CD68 (all regions)",
  segs_g1 = seg_panck, segs_g2 = seg_cd68,
  label_g1 = "panCK", label_g2 = "CD68"
)
fwrite(res_panck_cd68, file.path(out_dir, "geomx_de_panck_vs_cd68.csv"))

res_cd45_cd68 <- run_contrast(
  name = "CD45 vs CD68 (all regions)",
  segs_g1 = seg_cd45, segs_g2 = seg_cd68,
  label_g1 = "CD45", label_g2 = "CD68"
)
fwrite(res_cd45_cd68, file.path(out_dir, "geomx_de_cd45_vs_cd68.csv"))

# ---------------- SUMMARY ---------------------------------------------------
summarize_one <- function(name, df) {
  data.frame(
    contrast        = name,
    n_segs_group1   = unique(df$n_segs_group1),
    n_segs_group2   = unique(df$n_segs_group2),
    n_genes_tested  = sum(is.finite(df$pval)),
    n_sig_padj05    = sum(df$padj_bh < 0.05, na.rm = TRUE),
    n_sig_padj01    = sum(df$padj_bh < 0.01, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
summary_df <- rbind(
  summarize_one("SH_vs_PT_CD68",   res_sh_pt),
  summarize_one("SH_vs_LS_CD68",   res_sh_ls),
  summarize_one("panCK_vs_CD68",   res_panck_cd68),
  summarize_one("CD45_vs_CD68",    res_cd45_cd68)
)
fwrite(summary_df, file.path(out_dir, "geomx_de_summary.tsv"), sep = "\t")
cat("\n========================================================\n")
cat("Summary:\n")
print(summary_df)

# ---------------- VALIDATION ------------------------------------------------
cat("\n========================================================\n")
cat("VALIDATION CHECKS\n")
cat("========================================================\n")

# Top 10 by logFC (SH-up = positive logFC since group1 = SH, group2 = PT)
sh_pt_up <- res_sh_pt[res_sh_pt$logFC > 0 & is.finite(res_sh_pt$logFC), ]
sh_pt_up <- sh_pt_up[order(-sh_pt_up$logFC), ]
cat("\n[SH vs PT] Top 10 up-regulated by logFC (SH-up):\n")
print(head(sh_pt_up[, c("gene_symbol", "logFC", "pval", "padj_bh")], 10))

# Top 10 by significance + positive logFC for SH vs PT
sh_pt_sig_up <- res_sh_pt[res_sh_pt$logFC > 0 & is.finite(res_sh_pt$pval), ]
sh_pt_sig_up <- sh_pt_sig_up[order(sh_pt_sig_up$pval), ]
cat("\n[SH vs PT] Top 10 SH-up by p-value:\n")
print(head(sh_pt_sig_up[, c("gene_symbol", "logFC", "pval", "padj_bh")], 10))

# Marker recovery: at least 6 of these 7 in top 50 by logFC of SH-up
expected_markers <- c("GPNMB", "LPL", "FABP5", "HS3ST2", "MSR1", "CD36", "HLA-DRA")
top50_up <- head(sh_pt_up$gene_symbol, 50)
in_top50 <- expected_markers %in% top50_up
n_found <- sum(in_top50)
missing  <- expected_markers[!in_top50]
cat("\nMarker recovery (SH-vs-PT top 50 by logFC):\n")
cat("  Found:", n_found, "/ 7  (", paste(expected_markers[in_top50], collapse = ","), ")\n")
if (length(missing) > 0) {
  cat("  Missing:", paste(missing, collapse = ","), "\n")
  # Show their actual rank for diagnostics
  for (g in missing) {
    rank_up <- which(sh_pt_up$gene_symbol == g)
    if (length(rank_up) == 1L) {
      cat("    ", g, " is at SH-up rank by logFC:", rank_up, "\n")
    } else {
      cat("    ", g, " not present in SH-up table (length =",
          length(rank_up), ")\n")
    }
  }
}
pass_fail <- if (n_found >= 6L) "PASS" else "FAIL"
cat("Marker recovery check:", pass_fail, "(threshold: >=6/7)\n")

# Top 10 up-regulated for SH-vs-LS
sh_ls_up <- res_sh_ls[res_sh_ls$logFC > 0 & is.finite(res_sh_ls$logFC), ]
sh_ls_up <- sh_ls_up[order(-sh_ls_up$logFC), ]
cat("\n[SH vs LS] Top 10 up-regulated by logFC (SH-up):\n")
print(head(sh_ls_up[, c("gene_symbol", "logFC", "pval", "padj_bh")], 10))

# Top 5 by logFC for the other two contrasts
panck_up <- res_panck_cd68[res_panck_cd68$logFC > 0 & is.finite(res_panck_cd68$logFC), ]
panck_up <- panck_up[order(-panck_up$logFC), ]
cat("\n[panCK vs CD68] Top 5 panCK-up by logFC:\n")
print(head(panck_up[, c("gene_symbol", "logFC", "pval", "padj_bh")], 5))

cd45_up <- res_cd45_cd68[res_cd45_cd68$logFC > 0 & is.finite(res_cd45_cd68$logFC), ]
cd45_up <- cd45_up[order(-cd45_up$logFC), ]
cat("\n[CD45 vs CD68] Top 5 CD45-up by logFC:\n")
print(head(cd45_up[, c("gene_symbol", "logFC", "pval", "padj_bh")], 5))

# Compare significant counts to paper expectations
cat("\nComparison to paper-reported sig gene counts (padj<0.05):\n")
cat("  SH vs PT:  observed =", summary_df$n_sig_padj05[summary_df$contrast == "SH_vs_PT_CD68"],
    " | paper expected ~354\n")
cat("  SH vs LS:  observed =", summary_df$n_sig_padj05[summary_df$contrast == "SH_vs_LS_CD68"],
    " | paper expected ~215\n")

cat("\n[", format(Sys.time(), "%H:%M:%S"), "] DONE\n")
