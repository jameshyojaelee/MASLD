#!/usr/bin/env Rscript
# ============================================================================
# Govaere 2026 Nat Genet GeoMx WTA -- DE v2: more powerful statistical methods
# ----------------------------------------------------------------------------
# Re-runs SH-vs-PT and SH-vs-LS (CD68 segments) with three additional methods
# to attempt recovery of the paper-reported sig gene counts (~354 SH-vs-PT,
# ~215 SH-vs-LS). The baseline lmer arm (Script 42_govaere2026_geomx_de.R)
# returned n_sig_padj05 = 0 for both contrasts (BH burden on 13,969 genes x
# 7 patients).
#
# APPROACHES:
#   A) limma + duplicateCorrelation: lmFit(block=Patient, correlation=
#      duplicateCorrelation$consensus); eBayes(trend=TRUE) -- canonical
#      GeoMxTools fallback when n is small. Shares variance across genes
#      via empirical Bayes.
#   B) Permutation-based FDR: permute group labels WITHIN patient blocks 1000x;
#      empirical p = (# perms with |t| >= |t_obs|) / 1000; BH adjust.
#   C) 1-sided directional tests for the canonical SH macrophage panel
#      {GPNMB, LPL, FABP5, HS3ST2, MSR1, CD36, HLA-DRA, CD163, SPP1, TREM2}
#      via lmerTest::lmer + 1-sided p (alternative = "greater").
#
# DOES NOT overwrite existing geomx_de_*.csv (lmer baseline already in atlas).
#
# OUTPUTS (under Analysis/Spatial/results/govaere2026/):
#   geomx_de_sh_vs_pt_limma.csv
#   geomx_de_sh_vs_ls_limma.csv
#   geomx_de_sh_vs_pt_permFDR.csv
#   geomx_de_sh_vs_ls_permFDR.csv
#   geomx_de_canonical_markers_1sided.tsv
#   geomx_de_v2_summary.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(lmerTest)
  library(BiocParallel)
})

# ---------------- PATHS ----------------------------------------------------
project_root <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
out_dir   <- file.path(project_root, "Analysis/Spatial/results/govaere2026")
expr_path <- file.path(out_dir, "geomx_expression.tsv")
meta_path <- file.path(out_dir, "geomx_metadata.tsv")
stopifnot(file.exists(expr_path), file.exists(meta_path))

n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
n_perm    <- as.integer(Sys.getenv("GEOMX_N_PERM", "1000"))
seed0     <- 42L
set.seed(seed0)
cat("[", format(Sys.time(), "%H:%M:%S"), "] BiocParallel workers:", n_workers, "\n")
cat("[", format(Sys.time(), "%H:%M:%S"), "] Permutations per contrast:", n_perm, "\n")

# ---------------- LOAD INPUTS ----------------------------------------------
cat("[", format(Sys.time(), "%H:%M:%S"), "] Reading expression + metadata\n")
expr_dt <- fread(expr_path, header = TRUE, sep = "\t", check.names = FALSE)
meta    <- fread(meta_path, header = TRUE, sep = "\t", check.names = FALSE)

# File layout: first col = "segment" (segment IDs), remaining cols = genes.
# So expr_mat_raw is segments x genes; we transpose to genes x segments for
# limma + per-gene loops.
stopifnot("segment" %in% colnames(expr_dt))
seg_ids      <- expr_dt$segment
expr_mat_raw <- as.matrix(expr_dt[, -1, with = FALSE])  # segments x genes (linear Q3)
rownames(expr_mat_raw) <- seg_ids

# log2(Q3 + 1) -- matches paper Methods (page 14) and Script 42.
expr_mat_seg_gene <- log2(expr_mat_raw + 1)
expr_mat <- t(expr_mat_seg_gene)  # genes x segments (limma layout)

gene_ids <- rownames(expr_mat)
cat("Expression matrix (segments x genes):",
    paste(dim(expr_mat_seg_gene), collapse = " x "), "\n")
cat("limma-layout (genes x segments):",
    paste(dim(expr_mat), collapse = " x "), "\n")
cat("log2(Q3+1) value range:",
    paste(round(range(expr_mat, na.rm = TRUE), 3), collapse = " - "), "\n")

# Align metadata to expression segments
meta <- meta[match(seg_ids, meta$segment), ]
stopifnot(identical(meta$segment, seg_ids))

# ---------------- SEGMENT GROUPS -------------------------------------------
seg_sh_cd68 <- meta$segment[meta$segment_marker == "CD68" &
                              meta$region_type == "steatohepatitis"]
seg_pt_cd68 <- meta$segment[meta$segment_marker == "CD68" &
                              meta$region_type == "portal_tract"]
seg_ls_cd68 <- meta$segment[meta$segment_marker == "CD68" &
                              meta$region_type == "low_steatosis"]

cat("\nSegment availability (CD68 only):\n")
cat("  SH:", length(seg_sh_cd68),
    "  PT:", length(seg_pt_cd68),
    "  LS:", length(seg_ls_cd68), "\n")

# ============================================================================
# APPROACH A: limma + duplicateCorrelation
# ============================================================================
# Q3-normalized log expression. lmFit treats Patient as a within-block
# correlation structure (duplicateCorrelation), then eBayes(trend=TRUE) shares
# variance across genes. This is the canonical GeoMxTools small-n fallback.
run_limma_dupcor <- function(name, segs_g1, segs_g2, label_g1, label_g2) {
  cat("\n----- limma + dupCor:", name, "-----\n")
  segs_use   <- c(segs_g1, segs_g2)
  group_vec  <- factor(c(rep(label_g1, length(segs_g1)),
                         rep(label_g2, length(segs_g2))),
                       levels = c(label_g2, label_g1))   # ref = g2
  patient_vec <- meta$Patient[match(segs_use, meta$segment)]
  E <- expr_mat[, segs_use, drop = FALSE]   # genes x segs (limma layout)
  cat("  Genes x segs:", paste(dim(E), collapse = " x "), "\n")
  cat("  Patients:", length(unique(patient_vec)), "  Groups:", nlevels(group_vec), "\n")
  cat("  Patients(g1):", paste(sort(unique(meta$Patient[match(segs_g1, meta$segment)])), collapse=","), "\n")
  cat("  Patients(g2):", paste(sort(unique(meta$Patient[match(segs_g2, meta$segment)])), collapse=","), "\n")

  design <- model.matrix(~group_vec)
  colnames(design) <- c("Intercept", paste0("group", label_g1))

  # duplicateCorrelation expects block = patient. trend=TRUE on log Q3 -- the
  # mean-variance trend on continuous-valued normalized data is the right
  # choice (voom not needed because input is not counts).
  dc <- limma::duplicateCorrelation(E, design, block = patient_vec)
  cat("  duplicateCorrelation consensus:", round(dc$consensus, 4), "\n")

  fit  <- limma::lmFit(E, design, block = patient_vec, correlation = dc$consensus)
  fit2 <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit2, coef = 2, number = Inf, sort.by = "P")

  out <- data.frame(
    gene_symbol    = rownames(tt),
    logFC          = tt$logFC,
    AveExpr        = tt$AveExpr,
    t_or_z_stat    = tt$t,
    pval           = tt$P.Value,
    padj_bh        = tt$adj.P.Val,
    B              = tt$B,
    consensus_corr = dc$consensus,
    n_segs_group1  = length(segs_g1),
    n_segs_group2  = length(segs_g2),
    stringsAsFactors = FALSE
  )
  out
}

res_sh_pt_limma <- run_limma_dupcor(
  "SH vs PT (CD68 only)",
  seg_sh_cd68, seg_pt_cd68, "SH", "PT"
)
fwrite(res_sh_pt_limma, file.path(out_dir, "geomx_de_sh_vs_pt_limma.csv"))

res_sh_ls_limma <- run_limma_dupcor(
  "SH vs LS (CD68 only)",
  seg_sh_cd68, seg_ls_cd68, "SH", "LS"
)
fwrite(res_sh_ls_limma, file.path(out_dir, "geomx_de_sh_vs_ls_limma.csv"))

# ============================================================================
# APPROACH B: Permutation-based FDR
# ============================================================================
# Observed test stat = limma + dupCor t value (already computed). For each
# permutation, permute group labels WITHIN each patient block (so the
# patient-level random effect is preserved). Recompute lmFit + dupCor +
# eBayes -> per-gene |t| under null. Empirical p = (# perms with
# |t_perm| >= |t_obs|) / n_perm; BH adjust the empirical p.
#
# Note: we re-fit duplicateCorrelation each permutation because the consensus
# depends on the group labels through `design`. This is slower but correct.
# To bound runtime we pass design w/o intercept change between perms.
perm_within_block <- function(group_vec, patient_vec, rng_seed) {
  set.seed(rng_seed)
  gnew <- group_vec
  for (p in unique(patient_vec)) {
    idx <- which(patient_vec == p)
    gnew[idx] <- sample(group_vec[idx])
  }
  gnew
}

compute_t_one_design <- function(E, group_vec, patient_vec) {
  design <- model.matrix(~group_vec)
  dc <- tryCatch(
    limma::duplicateCorrelation(E, design, block = patient_vec),
    error = function(e) list(consensus = 0)
  )
  cons <- dc$consensus
  if (!is.finite(cons)) cons <- 0
  fit  <- limma::lmFit(E, design, block = patient_vec, correlation = cons)
  fit2 <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  fit2$t[, 2]
}

run_permFDR <- function(name, segs_g1, segs_g2, label_g1, label_g2,
                        n_perm_use = n_perm, workers = n_workers) {
  cat("\n----- permFDR:", name, "-----\n")
  segs_use    <- c(segs_g1, segs_g2)
  group_vec   <- factor(c(rep(label_g1, length(segs_g1)),
                          rep(label_g2, length(segs_g2))),
                        levels = c(label_g2, label_g1))
  patient_vec <- meta$Patient[match(segs_use, meta$segment)]
  E <- expr_mat[, segs_use, drop = FALSE]
  cat("  Genes x segs:", paste(dim(E), collapse = " x "), "\n")
  cat("  Patients:", length(unique(patient_vec)), "  n_perm:", n_perm_use, "\n")

  # How many of the within-patient permutations are non-trivial? If a patient
  # has only g1 or only g2 segments, permuting within that patient is a no-op.
  pat_grp <- table(patient_vec, group_vec)
  cat("  Patient x group table:\n")
  print(pat_grp)
  swap_pats <- rownames(pat_grp)[apply(pat_grp, 1, function(r) all(r > 0))]
  cat("  Patients with BOTH groups (effective swap units):",
      length(swap_pats), "  IDs:", paste(swap_pats, collapse = ","), "\n")
  if (length(swap_pats) == 0L) {
    cat("  WARN: zero patients carry both groups. permFDR null = empty;\n")
    cat("  falling back to UNRESTRICTED permutation (ignores patient block).\n")
  }

  # Observed t
  t_obs <- compute_t_one_design(E, group_vec, patient_vec)
  abs_t_obs <- abs(t_obs)

  # Permutation null
  bp <- if (.Platform$OS.type == "unix") {
    MulticoreParam(workers = workers, RNGseed = NULL, progressbar = FALSE, log = FALSE)
  } else {
    SnowParam(workers = workers, RNGseed = NULL, progressbar = FALSE)
  }
  t0 <- Sys.time()
  perm_seeds <- seed0 + seq_len(n_perm_use)
  perm_t_list <- bplapply(perm_seeds, function(s) {
    if (length(swap_pats) > 0L) {
      gnew <- perm_within_block(group_vec, patient_vec, s)
    } else {
      set.seed(s); gnew <- sample(group_vec)
    }
    # If permutation reproduced original labels (small swap space), keep --
    # this only inflates empirical p slightly toward the null and is safer
    # than skipping permutations.
    compute_t_one_design(E, gnew, patient_vec)
  }, BPPARAM = bp)
  t1 <- Sys.time()
  cat("  Permutations done in",
      round(as.numeric(difftime(t1, t0, units = "secs")), 1), "sec\n")

  # Build perm matrix genes x n_perm
  perm_mat <- do.call(cbind, perm_t_list)
  abs_perm <- abs(perm_mat)

  # Empirical p (two-sided) per gene
  # Use (count + 1) / (n + 1) to avoid p = 0.
  count_geq <- rowSums(abs_perm >= abs_t_obs)
  emp_p <- (count_geq + 1) / (n_perm_use + 1)
  emp_p[!is.finite(abs_t_obs)] <- NA_real_

  # Pull logFC + AveExpr from a single observed limma fit for the output
  design <- model.matrix(~group_vec)
  dc <- limma::duplicateCorrelation(E, design, block = patient_vec)
  fit  <- limma::lmFit(E, design, block = patient_vec, correlation = dc$consensus)
  fit2 <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  tt   <- limma::topTable(fit2, coef = 2, number = Inf, sort.by = "none")

  out <- data.frame(
    gene_symbol    = rownames(tt),
    logFC          = tt$logFC,
    AveExpr        = tt$AveExpr,
    t_obs          = t_obs,
    parametric_p   = tt$P.Value,
    emp_p          = emp_p,
    padj_emp_bh    = p.adjust(emp_p, method = "BH"),
    n_perm         = n_perm_use,
    n_swap_patients = length(swap_pats),
    n_segs_group1  = length(segs_g1),
    n_segs_group2  = length(segs_g2),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$emp_p, na.last = TRUE), ]
  out
}

res_sh_pt_perm <- run_permFDR(
  "SH vs PT (CD68 only)",
  seg_sh_cd68, seg_pt_cd68, "SH", "PT"
)
fwrite(res_sh_pt_perm, file.path(out_dir, "geomx_de_sh_vs_pt_permFDR.csv"))

res_sh_ls_perm <- run_permFDR(
  "SH vs LS (CD68 only)",
  seg_sh_cd68, seg_ls_cd68, "SH", "LS"
)
fwrite(res_sh_ls_perm, file.path(out_dir, "geomx_de_sh_vs_ls_permFDR.csv"))

# ============================================================================
# APPROACH C: 1-sided directional lmer for canonical SH macrophage markers
# ============================================================================
markers <- c("GPNMB", "LPL", "FABP5", "HS3ST2", "MSR1", "CD36",
             "HLA-DRA", "CD163", "SPP1", "TREM2")

fit_lmer_1sided <- function(g, segs_g1, segs_g2, label_g1, label_g2) {
  if (!(g %in% rownames(expr_mat))) {
    return(list(logFC = NA_real_, t_stat = NA_real_,
                p_2sided = NA_real_, p_1sided_greater = NA_real_,
                df = NA_real_, status = "gene_not_in_matrix"))
  }
  segs_use   <- c(segs_g1, segs_g2)
  group_vec  <- factor(c(rep(label_g1, length(segs_g1)),
                         rep(label_g2, length(segs_g2))),
                       levels = c(label_g2, label_g1))   # ref = g2
  patient_vec <- meta$Patient[match(segs_use, meta$segment)]
  y <- expr_mat[g, segs_use]
  df_use <- data.frame(y = y, group = group_vec, Patient = patient_vec)
  fit <- try(suppressMessages(suppressWarnings(
    lmerTest::lmer(y ~ group + (1 | Patient), data = df_use,
                   REML = TRUE,
                   control = lmerControl(check.conv.singular = .makeCC(action = "ignore",
                                                                        tol = 1e-4)))
  )), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(list(logFC = NA_real_, t_stat = NA_real_,
                p_2sided = NA_real_, p_1sided_greater = NA_real_,
                df = NA_real_, status = "lmer_fail"))
  }
  co <- summary(fit)$coefficients
  if (nrow(co) < 2L) {
    return(list(logFC = NA_real_, t_stat = NA_real_,
                p_2sided = NA_real_, p_1sided_greater = NA_real_,
                df = NA_real_, status = "no_group_coef"))
  }
  est  <- co[2, "Estimate"]
  tval <- co[2, "t value"]
  df_v <- co[2, "df"]
  p2   <- co[2, "Pr(>|t|)"]
  # 1-sided "greater" (group1 > group2)
  p1   <- if (is.finite(tval) && is.finite(df_v)) {
    pt(tval, df = df_v, lower.tail = FALSE)
  } else {
    NA_real_
  }
  list(logFC = est, t_stat = tval, p_2sided = p2,
       p_1sided_greater = p1, df = df_v, status = "ok")
}

run_canonical_1sided <- function(contrast_name, segs_g1, segs_g2, label_g1, label_g2) {
  rows <- lapply(markers, function(g) {
    r <- fit_lmer_1sided(g, segs_g1, segs_g2, label_g1, label_g2)
    data.frame(contrast = contrast_name, gene = g,
               logFC = r$logFC, t_stat = r$t_stat, df = r$df,
               p_2sided = r$p_2sided, p_1sided_greater = r$p_1sided_greater,
               status = r$status, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

cat("\n----- 1-sided canonical markers -----\n")
canon_sh_pt <- run_canonical_1sided("SH_vs_PT_CD68",
  seg_sh_cd68, seg_pt_cd68, "SH", "PT")
canon_sh_ls <- run_canonical_1sided("SH_vs_LS_CD68",
  seg_sh_cd68, seg_ls_cd68, "SH", "LS")
canon_out <- rbind(canon_sh_pt, canon_sh_ls)
fwrite(canon_out, file.path(out_dir, "geomx_de_canonical_markers_1sided.tsv"), sep = "\t")
cat("Canonical marker results:\n"); print(canon_out)

# ============================================================================
# SUMMARY
# ============================================================================
n_sig <- function(p, thr) sum(p < thr, na.rm = TRUE)

summary_df <- data.frame(
  contrast = c("SH_vs_PT_CD68", "SH_vs_LS_CD68"),
  method   = NA_character_,
  n_segs_group1 = NA_integer_,
  n_segs_group2 = NA_integer_,
  n_genes_tested = NA_integer_,
  n_sig_padj05 = NA_integer_,
  n_sig_padj01 = NA_integer_,
  stringsAsFactors = FALSE
)
summary_rows <- list()
add_row <- function(contrast, method, df, padj_col) {
  data.frame(
    contrast = contrast,
    method   = method,
    n_segs_group1 = unique(df$n_segs_group1)[1],
    n_segs_group2 = unique(df$n_segs_group2)[1],
    n_genes_tested = sum(is.finite(df[[padj_col]])),
    n_sig_padj05 = n_sig(df[[padj_col]], 0.05),
    n_sig_padj01 = n_sig(df[[padj_col]], 0.01),
    stringsAsFactors = FALSE
  )
}
summary_rows[["sh_pt_limma"]] <- add_row("SH_vs_PT_CD68", "limma_dupCor", res_sh_pt_limma, "padj_bh")
summary_rows[["sh_ls_limma"]] <- add_row("SH_vs_LS_CD68", "limma_dupCor", res_sh_ls_limma, "padj_bh")
summary_rows[["sh_pt_perm"]]  <- add_row("SH_vs_PT_CD68", "permFDR_BH",   res_sh_pt_perm,  "padj_emp_bh")
summary_rows[["sh_ls_perm"]]  <- add_row("SH_vs_LS_CD68", "permFDR_BH",   res_sh_ls_perm,  "padj_emp_bh")

# Also bring forward original lmer counts so the report is self-contained.
lmer_summary_path <- file.path(out_dir, "geomx_de_summary.tsv")
if (file.exists(lmer_summary_path)) {
  lmer_sum <- fread(lmer_summary_path)
  for (cn in c("SH_vs_PT_CD68", "SH_vs_LS_CD68")) {
    row <- lmer_sum[contrast == cn]
    if (nrow(row) == 1) {
      summary_rows[[paste0(cn, "_lmer")]] <- data.frame(
        contrast = cn,
        method = "lmer_BH (baseline)",
        n_segs_group1 = row$n_segs_group1,
        n_segs_group2 = row$n_segs_group2,
        n_genes_tested = row$n_genes_tested,
        n_sig_padj05 = row$n_sig_padj05,
        n_sig_padj01 = row$n_sig_padj01,
        stringsAsFactors = FALSE
      )
    }
  }
}
summary_df <- do.call(rbind, summary_rows)
fwrite(summary_df, file.path(out_dir, "geomx_de_v2_summary.tsv"), sep = "\t")

cat("\n=================== SUMMARY ===================\n")
print(summary_df)

# ============================================================================
# VALIDATION
# ============================================================================
cat("\n=================== VALIDATION ===================\n")
cat("Paper expected: SH-vs-PT ~354 sig | SH-vs-LS ~215 sig (FDR<0.05).\n")

cat("\n[SH vs PT] limma+dupCor top 10 by p:\n")
print(head(res_sh_pt_limma[order(res_sh_pt_limma$pval),
            c("gene_symbol", "logFC", "AveExpr", "t_or_z_stat", "pval", "padj_bh")], 10))

cat("\n[SH vs PT] permFDR top 10 by emp_p:\n")
print(head(res_sh_pt_perm[order(res_sh_pt_perm$emp_p),
            c("gene_symbol", "logFC", "t_obs", "emp_p", "padj_emp_bh")], 10))

cat("\n[SH vs LS] limma+dupCor top 10 by p:\n")
print(head(res_sh_ls_limma[order(res_sh_ls_limma$pval),
            c("gene_symbol", "logFC", "AveExpr", "t_or_z_stat", "pval", "padj_bh")], 10))

# Canonical markers in limma top-N?
cat("\nCanonical SH markers in limma SH-vs-PT (sorted by p):\n")
m_idx <- match(markers, res_sh_pt_limma$gene_symbol)
canon_in_limma <- data.frame(
  gene = markers,
  in_matrix = !is.na(m_idx),
  logFC = res_sh_pt_limma$logFC[m_idx],
  pval  = res_sh_pt_limma$pval[m_idx],
  padj_bh = res_sh_pt_limma$padj_bh[m_idx],
  rank_by_p = match(markers,
                    res_sh_pt_limma$gene_symbol[order(res_sh_pt_limma$pval)])
)
print(canon_in_limma)

cat("\n[", format(Sys.time(), "%H:%M:%S"), "] DONE\n")
