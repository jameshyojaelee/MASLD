#!/usr/bin/env Rscript
# ===========================================================================
# recompute_integration_only_cohort_support_C2.R
#
# Recomputes the manuscript's Fig 3 "cohort-support breakdown" (A1) and
# "integration-only DEGs" (A2) numbers against the C2 canonical
# (limma_voom_qw, canonical_deg_results.csv), replacing the stale dream/STAR-s0
# values (1,877 supported; 593/506/268/83; 4,653 integration-only; 49.5% vs
# 13.4% directional concordance, Fisher OR=6.4).
#
# Method (C2-consistent):
#   * Integrated set  = canonical_deg_results.csv (pooled LVQW, 5 mega cohorts).
#   * Per-cohort DE   = re-fit with the SAME engine (de_engine_lvqw::fit_lvqw,
#                       voomWithQualityWeights -> eBayes), one cohort at a time,
#                       Disease-vs-Control (group_binary), design
#                       ~ inferred_sex + group_binary (sex dropped if degenerate),
#                       filterByExpr + calcNormFactors(RLE) WITHIN each cohort
#                       (mirrors lvqw_loo_cv_C2.R / 05h producer filter).
#
# Definitions (locked to the manuscript prose):
#   per-cohort "called" bar = adj.P.Val<0.05 & |logFC|>0.5 & same direction as
#                             the integrated effect.
#   A1 cohort-support : over integrated Tier-1 (padj<0.05 & |logFC|>0.5),
#                       n_support = # of 5 cohorts meeting the per-cohort bar.
#   A2 integration-only : over integrated padj<0.05 (no LFC filter), genes with
#                       n_support==0 (fail the per-cohort bar in EVERY cohort).
#       directional concordance: among integration-only genes, fraction whose
#       per-cohort logFC sign matches the integrated sign in >=4 of 5 cohorts,
#       versus an AveExpr-size-matched set of integrated-nonsignificant genes
#       (padj>=0.05). Fisher exact -> OR + p.
#
# Output: results/integration/integration_only_C2/
# ===========================================================================
set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(yaml)
})

BASE    <- Sys.getenv("MASLD_PROJECT_ROOT",
             "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT     <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR    <- file.path(INT, "results/integration")
SCRIPTS <- file.path(INT, "scripts")
OUT     <- file.path(RDIR, "integration_only_C2")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source(file.path(SCRIPTS, "de_engine_lvqw.R"))

strip_ver <- function(x) sub("\\.[0-9]+$", "", x)

# ---------------------------------------------------------------------------
# 1. Integrated set (C2 canonical)
# ---------------------------------------------------------------------------
cat("== Loading C2 canonical integrated DEGs ==\n")
intg <- fread(file.path(RDIR, "canonical_deg_results.csv"))
stopifnot(all(c("gene","logFC","padj","AveExpr") %in% names(intg)))
intg[, gid := strip_ver(gene)]
intg_sig05  <- intg[!is.na(padj) & padj < 0.05]                       # 13,043
intg_tier1  <- intg_sig05[abs(logFC) > 0.5]                           # 1,853
cat(sprintf("  integrated genes tested: %d\n", nrow(intg)))
cat(sprintf("  integrated padj<0.05 (no LFC): %d\n", nrow(intg_sig05)))
cat(sprintf("  integrated Tier-1 padj<0.05 & |logFC|>0.5: %d\n", nrow(intg_tier1)))

# ---------------------------------------------------------------------------
# 2. Per-cohort LVQW DE (Disease vs Control), 5 mega cohorts
# ---------------------------------------------------------------------------
ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("\n== Mega cohorts:", paste(mega_cohorts, collapse = ", "), "==\n")

cat("Loading raw counts + metadata...\n")
merged_raw <- readRDS(file.path(RDIR, "merged_counts_raw.rds"))
meta_all   <- readRDS(file.path(RDIR, "meta_matched.rds")); setDT(meta_all)
qc <- fread(file.path(INT, "qc/sample_qc_report.csv"))
pass_ids <- qc[pass_technical == TRUE, sample_id]
merged_raw <- merged_raw[, colnames(merged_raw) %in% pass_ids]
meta_all   <- meta_all[sample_id %in% colnames(merged_raw)]
meta_all   <- meta_all[match(colnames(merged_raw), sample_id)]
stopifnot(all(meta_all$sample_id == colnames(merged_raw)))

per_cohort <- list()
percohort_counts <- list()
for (co in mega_cohorts) {
  idx <- meta_all$dataset == co
  cat(sprintf("\n-- %s: %d samples --\n", co, sum(idx)))
  m   <- meta_all[idx]
  dge <- DGEList(counts = merged_raw[, idx])
  gb  <- factor(m$group_binary, levels = c("Control", "Disease"))
  cat("   group_binary:", paste(names(table(gb)), table(gb), collapse = " "), "\n")
  if (nlevels(droplevels(gb)) < 2) { cat("   SKIP (single group)\n"); next }

  # filterByExpr within cohort (same producer filter style as 05h/LOO)
  keep <- filterByExpr(dge, design = model.matrix(~ 0 + gb))
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge, method = "RLE")
  cat("   genes post-filterByExpr:", sum(keep), "\n")

  info <- data.frame(
    inferred_sex = factor(m$inferred_sex),
    group_binary = gb,
    stringsAsFactors = FALSE
  )
  gd  <- build_design_guarded(info, c("inferred_sex", "group_binary"))
  design <- gd$design
  stopifnot("group_binaryDisease" %in% colnames(design))
  cat("   design:", gd$formula_used, "\n")

  res <- fit_lvqw(dge, design, coef = "group_binaryDisease",
                  do_ashr = FALSE, weights = TRUE)
  res[, gid := strip_ver(gene)]
  per_cohort[[co]] <- res[, .(gid, logFC, padj, AveExpr)]
  fwrite(res, file.path(OUT, paste0("percohort_lvqw_", co, ".csv")))

  n_deg <- res[!is.na(padj) & padj < 0.05 & abs(logFC) > 0.5, .N]
  percohort_counts[[co]] <- data.table(cohort = co, n_genes = nrow(res),
                                        n_deg_tier1 = n_deg)
  cat(sprintf("   per-cohort DEGs (padj<0.05 & |logFC|>0.5): %d\n", n_deg))
}
pc_counts <- rbindlist(percohort_counts)
fwrite(pc_counts, file.path(OUT, "per_cohort_deg_counts_C2.csv"))
cohorts <- names(per_cohort); K <- length(cohorts)
cat(sprintf("\n== %d cohorts fit; per-cohort Tier-1 DEG range: %d - %d ==\n",
            K, min(pc_counts$n_deg_tier1), max(pc_counts$n_deg_tier1)))

# ---------------------------------------------------------------------------
# 3. Build gene x cohort logFC / padj matrices (aligned on stripped gene id)
# ---------------------------------------------------------------------------
all_gids <- unique(intg$gid)
lfc_mat  <- matrix(NA_real_, nrow = length(all_gids), ncol = K,
                   dimnames = list(all_gids, cohorts))
padj_mat <- lfc_mat
for (co in cohorts) {
  d <- per_cohort[[co]]
  mi <- match(d$gid, all_gids); ok <- !is.na(mi)
  lfc_mat[mi[ok], co]  <- d$logFC[ok]
  padj_mat[mi[ok], co] <- d$padj[ok]
}

# per-cohort "called" bar: padj<0.05 & |logFC|>0.5 & same direction as integrated
intg_sign <- sign(intg$logFC[match(all_gids, intg$gid)])
called_mat <- (padj_mat < 0.05) & (abs(lfc_mat) > 0.5) &
              (sign(lfc_mat) == matrix(intg_sign, length(all_gids), K))
called_mat[is.na(called_mat)] <- FALSE
n_support <- rowSums(called_mat)
# per-cohort sign concordance (regardless of significance)
concord_mat <- sign(lfc_mat) == matrix(intg_sign, length(all_gids), K)
concord_mat[is.na(concord_mat)] <- FALSE
n_concord <- rowSums(concord_mat)
names(n_support) <- names(n_concord) <- all_gids

# ---------------------------------------------------------------------------
# 4. A1 — cohort-support breakdown over integrated Tier-1 (1,853)
# ---------------------------------------------------------------------------
t1_ns <- n_support[intg_tier1$gid]
supp_tab <- table(factor(t1_ns, levels = 0:K))
a1 <- data.table(
  n_cohorts_supporting = 0:K,
  n_genes = as.integer(supp_tab)
)
a1_summary <- data.table(
  metric = c("integrated_Tier1_total","supported_by_ge1","supported_by_ge2",
             paste0("exactly_", 0:K, "_cohorts")),
  value = c(nrow(intg_tier1), sum(t1_ns >= 1), sum(t1_ns >= 2),
            as.integer(supp_tab))
)
fwrite(a1, file.path(OUT, "A1_cohort_support_breakdown_C2.csv"))
fwrite(a1_summary, file.path(OUT, "A1_cohort_support_summary_C2.csv"))

cat("\n", strrep("=", 64), "\n  A1 COHORT-SUPPORT (over", nrow(intg_tier1),
    "Tier-1 integrated DEGs)\n", strrep("=", 64), "\n")
cat(sprintf("  supported by >=1 cohort: %d (%.1f%%)\n",
            sum(t1_ns >= 1), 100*mean(t1_ns >= 1)))
for (k in 1:K) cat(sprintf("  exactly %d cohort(s): %d\n", k, sum(t1_ns == k)))
cat(sprintf("  in NO cohort (Tier-1 but unsupported): %d\n", sum(t1_ns == 0)))

# ---------------------------------------------------------------------------
# 5. A2 — integration-only over integrated padj<0.05 (13,043)
# ---------------------------------------------------------------------------
sig_ns   <- n_support[intg_sig05$gid]
io_mask  <- sig_ns == 0
io_genes <- intg_sig05$gid[io_mask]
n_io     <- length(io_genes)
io_lfc05 <- sum(abs(intg_sig05$logFC[io_mask]) > 0.5)

cat("\n", strrep("=", 64), "\n  A2 INTEGRATION-ONLY (over", nrow(intg_sig05),
    "integrated padj<0.05 genes)\n", strrep("=", 64), "\n")
cat(sprintf("  integration-only (n_support==0): %d (%.1f%%)\n",
            n_io, 100*n_io/nrow(intg_sig05)))
cat(sprintf("  ... remaining after |logFC|>0.5: %d\n", io_lfc05))

# directional concordance: integration-only vs AveExpr-size-matched non-DEGs
io_ge4 <- n_concord[io_genes] >= 4
nonsig <- intg[!is.na(padj) & padj >= 0.05]           # integrated non-significant
# size-match on AveExpr deciles of the integration-only set
io_ave  <- intg_sig05$AveExpr[io_mask]
brk     <- unique(quantile(io_ave, probs = seq(0,1,0.1), na.rm = TRUE))
io_bin  <- cut(io_ave, breaks = brk, include.lowest = TRUE)
ns_bin  <- cut(nonsig$AveExpr, breaks = brk, include.lowest = TRUE)
matched <- character(0)
for (b in levels(io_bin)) {
  need <- sum(io_bin == b, na.rm = TRUE)
  pool <- nonsig$gid[which(ns_bin == b)]
  if (length(pool) == 0 || need == 0) next
  matched <- c(matched, sample(pool, need, replace = length(pool) < need))
}
m_ge4 <- n_concord[matched] >= 4

p_io <- mean(io_ge4); p_m <- mean(m_ge4)
ct <- matrix(c(sum(io_ge4), sum(!io_ge4), sum(m_ge4), sum(!m_ge4)),
             nrow = 2, byrow = TRUE,
             dimnames = list(c("integration_only","matched_nonDEG"),
                             c("ge4of5_concordant","lt4")))
ft <- fisher.test(ct)
cat(sprintf("\n  directional concordance (>=4 of %d cohorts, same sign as integrated):\n", K))
cat(sprintf("    integration-only:      %.1f%% (%d/%d)\n",
            100*p_io, sum(io_ge4), length(io_ge4)))
cat(sprintf("    size-matched non-DEGs: %.1f%% (%d/%d)\n",
            100*p_m, sum(m_ge4), length(m_ge4)))
cat(sprintf("    Fisher OR = %.2f   p = %.3g\n", ft$estimate, ft$p.value))

a2_summary <- data.table(
  metric = c("integrated_padj05_total","integration_only_n","integration_only_pct",
             "integration_only_after_lfc05",
             "dirconcord_io_pct","dirconcord_matched_pct","fisher_OR","fisher_p",
             "n_cohorts"),
  value = c(nrow(intg_sig05), n_io, round(100*n_io/nrow(intg_sig05),1), io_lfc05,
            round(100*p_io,1), round(100*p_m,1),
            round(unname(ft$estimate),2), signif(ft$p.value,3), K)
)
fwrite(a2_summary, file.path(OUT, "A2_integration_only_summary_C2.csv"))
fwrite(data.table(gid = io_genes,
                  logFC = intg_sig05$logFC[io_mask],
                  padj  = intg_sig05$padj[io_mask],
                  n_concord = n_concord[io_genes]),
       file.path(OUT, "A2_integration_only_genes_C2.csv"))

cat("\n== DONE. Outputs in", OUT, "==\n")
