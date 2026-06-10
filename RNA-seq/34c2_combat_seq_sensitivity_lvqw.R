#!/usr/bin/env Rscript
# ===========================================================================
# 34c2_combat_seq_sensitivity_lvqw.R
# C2-canonical (limma_voom_qw__C2) re-run of the ComBat-seq batch-correction
# sensitivity analysis. Replaces dream with the C2 LVQW engine and compares the
# ComBat-corrected DE vs canonical_deg_results.csv.
#
# WHAT CHANGED vs 34_combat_seq_sensitivity.R:
#   - Engine on the ComBat-corrected counts: dream() -> fit_lvqw()
#     (voomWithQualityWeights -> lmFit -> eBayes).
#   - ComBat-seq is applied to the SAME mega DGE the C2 canonical was fit on
#     (merged_dge.rds subset to include_in_mega cohorts), NOT the full 9-cohort
#     merged set, so the comparison vs canonical_deg_results.csv is on the same
#     sample universe (n=846, 5 cohorts).
#   - Design after ComBat: ~ inferred_sex + group_binary. The `dataset` term is
#     DROPPED because ComBat-seq has already removed the batch (dataset) effect;
#     keeping it would double-correct. The C2 canonical instead adjusts `dataset`
#     as a fixed covariate (~ dataset + inferred_sex + group_binary) — this is
#     exactly the batch-handling contrast the sensitivity analysis exists to test.
#   - Comparator: canonical_deg_results.csv (C2), NOT dream_results.csv.
#   - Significance threshold for the set-overlap / r reported at padj<0.1 to keep
#     OLD->NEW apples-to-apples with the dream-era Script 34 (which used 0.1).
#
# Output: RNA-seq/results/audit_sensitivity/combat_seq_c2/
# SLURM: bigmem+interactive, --mem=500G (ComBat-seq is RAM-heavy), 8 CPUs, rnaseq.
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table); library(sva); library(edgeR); library(limma)
  library(ashr); library(ggplot2); library(yaml)
})
set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
source(file.path(INT, "scripts/de_engine_lvqw.R"))   # fit_lvqw(), library-only
OUTDIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity/combat_seq_c2")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 34c2: ComBat-seq sensitivity on the C2 LVQW engine ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# ---------------------------------------------------------------------------
# 1. Load mega DGE (subset to include_in_mega cohorts — same universe as 05h)
# ---------------------------------------------------------------------------
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
dge  <- dge[, dge$samples$dataset %in% mega]
cat("  Mega cohorts:", paste(mega, collapse = ", "), "\n")
cat("  Samples:", ncol(dge), " Genes:", nrow(dge), "\n")

meta <- as.data.table(dge$samples)
condition_var <- "group_binary"
stopifnot(condition_var %in% names(meta))

# C2-parity sex (from meta_matched, as in 05h)
meta_matched <- readRDS(file.path(RDIR, "meta_matched.rds"))
sex <- meta_matched$inferred_sex[match(colnames(dge), meta_matched$sample_id)]

# ---------------------------------------------------------------------------
# 2. ComBat-seq (batch = dataset, biological group preserved)
# ---------------------------------------------------------------------------
cat("\n--- ComBat-seq (batch=dataset, group=group_binary preserved) ---\n")
batch <- as.factor(meta$dataset)
group <- as.factor(meta[[condition_var]])
counts_adj <- tryCatch(
  ComBat_seq(counts = as.matrix(dge$counts), batch = batch, group = group),
  error = function(e) {
    cat("  WARNING ComBat_seq w/ group failed:", conditionMessage(e), "-> retry no covar\n")
    ComBat_seq(counts = as.matrix(dge$counts), batch = batch)
  })
cat("  ComBat-seq adjusted counts:", paste(dim(counts_adj), collapse = " x "), "\n")

# ---------------------------------------------------------------------------
# 3. C2 LVQW engine on ComBat-corrected counts (dataset term dropped — already
#    removed by ComBat). Design: ~ inferred_sex + group_binary.
# ---------------------------------------------------------------------------
cat("\n--- C2 LVQW engine on ComBat-corrected counts ---\n")
dge_adj <- DGEList(counts = counts_adj, samples = meta)
dge_adj <- calcNormFactors(dge_adj, method = "RLE")
keep <- filterByExpr(dge_adj, group = meta[[condition_var]])
dge_adj <- dge_adj[keep, , keep.lib.sizes = FALSE]
cat("  Genes after filterByExpr:", nrow(dge_adj), "\n")

info <- data.frame(
  group_binary = factor(meta$group_binary, levels = c("Control", "Disease")),
  inferred_sex = factor(sex),
  row.names    = colnames(dge_adj))
design <- model.matrix(~ inferred_sex + group_binary, data = info)
coef_name <- "group_binaryDisease"
stopifnot(coef_name %in% colnames(design))
cat("  Design: ~ inferred_sex + group_binary (dataset removed by ComBat)\n")

res_combat <- fit_lvqw(dge_adj, design, coef_name, do_ashr = FALSE)
n_sig_combat <- sum(res_combat$padj < 0.1, na.rm = TRUE)
cat("  ComBat+LVQW DEGs (padj<0.1):", n_sig_combat, "\n")
fwrite(res_combat, file.path(OUTDIR, "c2_combat_seq_results.csv"))

# ---------------------------------------------------------------------------
# 4. Concordance vs C2 canonical (canonical_deg_results.csv)
# ---------------------------------------------------------------------------
cat("\n--- Concordance vs C2 canonical ---\n")
canon <- fread(file.path(RDIR, "canonical_deg_results.csv"))
if (!"padj" %in% names(canon) && "adj.P.Val" %in% names(canon))
  setnames(canon, "adj.P.Val", "padj")
canon[,      gene_clean := sub("\\..*", "", gene)]
res_combat[, gene_clean := sub("\\..*", "", gene)]

merged <- merge(
  canon[, .(gene_clean, logFC_canon = logFC, padj_canon = padj)],
  res_combat[, .(gene_clean, logFC_combat = logFC, padj_combat = padj)],
  by = "gene_clean")
cat("  Genes in both:", nrow(merged), "\n")

lfc_r   <- cor(merged$logFC_canon, merged$logFC_combat, use = "complete.obs")  # pearson
lfc_rho <- cor(merged$logFC_canon, merged$logFC_combat, method = "spearman", use = "complete.obs")
sig_canon  <- merged[padj_canon  < 0.1, gene_clean]
sig_combat <- merged[padj_combat < 0.1, gene_clean]
ov  <- length(intersect(sig_canon, sig_combat))
jac <- ov / max(1L, length(union(sig_canon, sig_combat)))
ov_dt <- merged[gene_clean %in% intersect(sig_canon, sig_combat)]
same_dir <- if (nrow(ov_dt)) sum(sign(ov_dt$logFC_canon) == sign(ov_dt$logFC_combat)) else NA_integer_
dir_pct  <- if (nrow(ov_dt)) 100 * same_dir / nrow(ov_dt) else NA_real_

# Tier-1 cross-ref (C2 canonical Tier-1 = padj<.05 & |lfc|>.5)
tier1_canon  <- merged[padj_canon  < 0.05 & abs(logFC_canon)  > 0.5, gene_clean]
tier1_combat <- merged[padj_combat < 0.05 & abs(logFC_combat) > 0.5, gene_clean]
jac_tier1 <- length(intersect(tier1_canon, tier1_combat)) /
             max(1L, length(union(tier1_canon, tier1_combat)))

cat(sprintf("  logFC pearson r = %.4f | spearman rho = %.4f\n", lfc_r, lfc_rho))
cat(sprintf("  DEGs padj<0.1 — canonical: %d, ComBat: %d, overlap: %d, jaccard: %.4f\n",
            length(sig_canon), length(sig_combat), ov, jac))
cat(sprintf("  Tier-1 — canonical: %d, ComBat: %d, jaccard: %.4f\n",
            length(tier1_canon), length(tier1_combat), jac_tier1))
cat(sprintf("  Direction concordance in overlap: %s\n",
            if (!is.na(dir_pct)) sprintf("%d/%d (%.1f%%)", same_dir, nrow(ov_dt), dir_pct) else "NA"))

conc <- data.table(
  metric = c("n_genes_both","lfc_pearson_r","lfc_spearman_rho",
             "sig_canon_p01","sig_combat_p01","sig_overlap_p01","jaccard_p01",
             "tier1_canon","tier1_combat","jaccard_tier1","direction_concordance"),
  value = c(nrow(merged), round(lfc_r,4), round(lfc_rho,4),
            length(sig_canon), length(sig_combat), ov, round(jac,4),
            length(tier1_canon), length(tier1_combat), round(jac_tier1,4),
            if (!is.na(dir_pct)) round(dir_pct/100,4) else NA_real_))
fwrite(conc, file.path(OUTDIR, "combat_seq_concordance_c2.csv"))

pdf(file.path(OUTDIR, "lfc_scatter_canon_vs_combat_c2.pdf"), width = 6, height = 6)
print(ggplot(merged, aes(logFC_canon, logFC_combat)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = sprintf("r = %.3f", lfc_r), color = "red", size = 4) +
  labs(title = "logFC: C2 canonical vs ComBat-seq + LVQW",
       x = "logFC (C2 canonical)", y = "logFC (ComBat-seq + LVQW)") +
  theme_bw(base_size = 11))
dev.off()

cat("\n=== 34c2 complete:", format(Sys.time()), "===\n")
cat("Outputs:\n  ", file.path(OUTDIR, "c2_combat_seq_results.csv"),
    "\n  ", file.path(OUTDIR, "combat_seq_concordance_c2.csv"),
    "\n  ", file.path(OUTDIR, "lfc_scatter_canon_vs_combat_c2.pdf"), "\n")
