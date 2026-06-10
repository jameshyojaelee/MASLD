#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Surgical S1 K_eff recompute against the C2 canonical bulk axis.
#
# Faithfully reproduces the S1 (human bulk DE, 5 sub-contrasts) Empirical K_eff
# computation in RNA-seq/46d_convergence_evidence.R (lines 134-396, 824-867)
# WITHOUT re-running the heavy permutation-FDR / full convergence ranking.
#
# All inputs are the same atlas columns 46d reads:
#   S1.1 overall  : bulk_logFC / bulk_tstat  (== C2 canonical_deg_results.csv)
#   S1.2 NAFL->NASH: nafl_vs_nash_logFC / nafl_vs_nash_tstat  (LVQW arm)
#   S1.3 F2 inflection: f2_inflection_logFC / f2_inflection_padj
#   S1.4 advanced fibrosis: adv_fib_logFC / adv_fib_padj
#   S1.5 sex interaction: (bulk_logFC_F - bulk_logFC_M) / sex_interaction_padj
#
# The bulk_* columns in the atlas were verified identical (corr 1.0, max abs
# diff 0.0 on ENSG merge) to canonical_deg_results.csv, so this IS the C2 axis.
#
# Output: RNA-seq/results/multi_evidence/convergence_evidence_kEff_C2.csv
#         (does NOT overwrite the canonical kEff produced by 46d).
# ---------------------------------------------------------------------------
suppressMessages({ library(data.table) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ME   <- file.path(BASE, "RNA-seq/results/multi_evidence")

# --- constants (copied verbatim from 46d) ----------------------------------
LOG_BF_FLOOR        <- log(0.1)   # -2.303
LOG_BF_CEIL         <- log(1e6)   # 13.816
LOG_BF_ACTIVE       <- log(3)
W_WAKEFIELD_LOGODDS <- 0.04

# --- BF primitives (copied verbatim from 46d) ------------------------------
clamp_logbf <- function(x, floor = LOG_BF_FLOOR, ceil = LOG_BF_CEIL) {
  x[is.na(x)] <- 0
  pmin(pmax(x, floor), ceil)
}
wakefield_abf <- function(effect, se, W) {
  out <- rep(0, length(effect))
  ok <- !is.na(effect) & !is.na(se) & se > 0 & is.finite(se)
  V <- se[ok]^2
  z <- effect[ok] / se[ok]
  shrink <- W / (V + W)
  out[ok] <- 0.5 * log(1 - shrink) + (z^2 / 2) * shrink
  clamp_logbf(out)
}
se_from_tstat <- function(effect, tstat) {
  ok <- !is.na(effect) & !is.na(tstat) & abs(tstat) > 0
  se <- rep(NA_real_, length(effect))
  se[ok] <- abs(effect[ok]) / abs(tstat[ok])
  se
}
se_from_padj <- function(effect, padj) {
  se <- rep(NA_real_, length(effect))
  ok <- !is.na(effect) & !is.na(padj) & padj > 0 & padj < 1
  idx <- which(ok)
  if (length(idx) > 0) {
    z <- abs(qnorm(padj[idx] / 2))
    usable <- z > 0.01 & is.finite(z)
    se[idx[usable]] <- abs(effect[idx[usable]]) / z[usable]
  }
  se
}

# --- load atlas ------------------------------------------------------------
atlas <- fread(file.path(ME, "multi_evidence_atlas.csv"))
gm    <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))
atlas[, ensembl_base := sub("\\..*", "", ensembl_id)]
gm_u <- unique(gm, by = "ensembl_base")
atlas[, gene_biotype_gm := gm_u$gene_biotype[match(ensembl_base, gm_u$ensembl_base)]]
atlas[, chromosome_gm  := gm_u$chromosome[match(ensembl_base, gm_u$ensembl_base)]]

# --- confounder filter (copied verbatim from 46d lines 178-197) ------------
sym <- atlas$human_symbol; sym[is.na(sym)] <- ""
drop_chrY  <- !is.na(atlas$chromosome_gm) & atlas$chromosome_gm == "chrY"
drop_chrM  <- !is.na(atlas$chromosome_gm) & atlas$chromosome_gm %in% c("chrM","chrMT")
drop_rRNA  <- !is.na(atlas$gene_biotype_gm) & atlas$gene_biotype_gm %in%
              c("rRNA","Mt_rRNA","rRNA_pseudogene")
drop_rpsl  <- grepl("^RPS[0-9]+[A-Z]*$|^RPL[0-9]+[A-Z]*$", sym) |
              sym %in% c("RPSA","RPLP0","RPLP1","RPLP2") |
              grepl("^RPL[0-9]+P[0-9]+$|^RPS[0-9]+P[0-9]+$", sym)
drop_mrpsl <- grepl("^MRPS[0-9]+[A-Z]*$|^MRPL[0-9]+[A-Z]*$", sym)
drop_ig    <- grepl("^IG[HKL][VJCD][0-9]", sym)
drop_hb    <- sym %in% c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2",
                         "HBM","HBQ1","HBZ")
drop_xesc  <- sym %in% c("XIST","TSIX","KDM6A","DDX3X","EIF1AX","UTX","UBA1",
                          "RLIM","ZFX","RPS4X","EIF2S3","TXLNG")
drop_hla   <- grepl("^HLA-", sym)
excluded <- drop_chrY | drop_chrM | drop_rRNA | drop_rpsl | drop_mrpsl |
            drop_ig | drop_hb | drop_xesc | drop_hla
cat(sprintf("Confounder filter: %d/%d genes excluded\n", sum(excluded), nrow(atlas)))

# --- S1 sub-contrast log-BFs (copied verbatim from 46d lines 356-394) ------
se_dream <- se_from_tstat(atlas$bulk_logFC, atlas$bulk_tstat)
overall_bf <- wakefield_abf(atlas$bulk_logFC, se_dream, W_WAKEFIELD_LOGODDS)

se_nn <- se_from_tstat(atlas$nafl_vs_nash_logFC, atlas$nafl_vs_nash_tstat)
nn_bf <- wakefield_abf(atlas$nafl_vs_nash_logFC, se_nn, W_WAKEFIELD_LOGODDS)

se_f2 <- se_from_padj(atlas$f2_inflection_logFC, atlas$f2_inflection_padj)
f2_bf <- wakefield_abf(atlas$f2_inflection_logFC, se_f2, W_WAKEFIELD_LOGODDS)

se_af <- se_from_padj(atlas$adv_fib_logFC, atlas$adv_fib_padj)
af_bf <- wakefield_abf(atlas$adv_fib_logFC, se_af, W_WAKEFIELD_LOGODDS)

sex_effect <- atlas$bulk_logFC_F - atlas$bulk_logFC_M
se_sex <- se_from_padj(sex_effect, atlas$sex_interaction_padj)
sex_bf <- wakefield_abf(sex_effect, se_sex, W_WAKEFIELD_LOGODDS)

# --- est_keff (copied verbatim from 46d lines 831-840) ---------------------
est_keff <- function(sub_contrast_mat, label) {
  mat <- sub_contrast_mat[!excluded, , drop = FALSE]
  R <- suppressWarnings(cor(mat, method = "spearman", use = "pairwise.complete.obs"))
  R[is.na(R)] <- 0
  lam <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  lam <- pmax(lam, 0)
  keff <- if (sum(lam^2) > 0) sum(lam)^2 / sum(lam^2) else 1
  list(K = ncol(mat), K_eff = keff, eigenvalues = lam, cor_mat = R, label = label)
}

s1_mat <- cbind(overall_bf, nn_bf, f2_bf, af_bf, sex_bf)
colnames(s1_mat) <- c("overall","nafl_vs_nash","f2_inflection","adv_fib","sex_interaction")
keff_s1 <- est_keff(s1_mat, "S1")

cat(sprintf("\nS1: K=%d, K_eff=%.4f\n", keff_s1$K, keff_s1$K_eff))
cat("S1 eigenvalues:", paste(round(keff_s1$eigenvalues, 3), collapse=";"), "\n")
cat("\nS1 sub-contrast log-BF correlation matrix (spearman):\n")
print(round(keff_s1$cor_mat, 3))
cat(sprintf("\nMean non-zero log-BF per sub-contrast:\n"))
for (cn in colnames(s1_mat)) {
  v <- s1_mat[!excluded, cn]; v <- v[v > LOG_BF_ACTIVE]
  cat(sprintf("  %-15s active=%6d  mean=%.3f\n", cn, length(v),
              if (length(v)) mean(v) else NA_real_))
}

keff_out <- data.table(
  modality = "S1",
  K = keff_s1$K,
  K_eff = keff_s1$K_eff,
  eigenvalues = paste(round(keff_s1$eigenvalues, 3), collapse=";")
)
out_f <- file.path(ME, "convergence_evidence_kEff_C2.csv")
fwrite(keff_out, out_f)
cat(sprintf("\nWrote %s\n", out_f))
