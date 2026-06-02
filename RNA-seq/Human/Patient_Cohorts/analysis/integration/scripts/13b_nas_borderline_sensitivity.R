#!/usr/bin/env Rscript
# 13b_nas_borderline_sensitivity.R
# ---------------------------------------------------------------------------
# NAS Borderline Sensitivity Analysis [R3 P1-1]
#
# Reviewer concern: NAS 3-4 patients are grouped with NASH, but clinically
# NAS 3-4 is indeterminate (Kleiner et al. 2005), not definite NASH.
#
# Current classification (Script 00 line 298):
#   NAS < 3  -> NAFL
#   NAS 3-4  -> Borderline (grouped with NASH in Script 13)
#   NAS >= 5 -> NASH
#
# Two sensitivity arms:
#   Arm 1 — Strict: NAS >= 5 only as NASH; NAS 3-4 EXCLUDED entirely
#   Arm 2 — Borderline separate: NAS 3-4 as own category, excluded from
#           both NAFL and NASH groups (NAFL-vs-NASH DE on NAS<3 vs NAS>=5)
#
# Both arms should yield identical samples (NAS<3 vs NAS>=5) — they differ
# only semantically and are collapsed into a single re-analysis to avoid
# redundant compute. The label "Arm 1/2" is retained in the report for
# reviewer traceability.
#
# Compares against canonical Script 13 results (NAS 3-4 grouped with NASH).
#
# Usage: Rscript 13b_nas_borderline_sensitivity.R
# SLURM: bigmem partition, 16 CPUs, 500G RAM, 48h
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
})

# Force injection into lme4 namespace BEFORE loading variancePartition
# (same pattern as 14.3_age_sensitivity.R)
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({
      unlockBinding(fn, ns_lme4)
      assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
      lockBinding(fn, ns_lme4)
    }, silent = TRUE)
  }
}

suppressPackageStartupMessages({
  library(variancePartition)
  library(BiocParallel)
  library(ggplot2)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
DSIG <- file.path(INT, "results/disease_signatures")
OUTDIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity/nas_borderline_sensitivity")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 13b: NAS Borderline Sensitivity Analysis [R3 P1-1] ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# --- Setup parallel backend ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus, progressbar = TRUE) else SerialParam()

# ============================================================
#  Load data (same pattern as Script 13)
# ============================================================
cat("Loading merged counts...\n")
counts <- readRDS(file.path(RDIR, "merged_counts_raw.rds"))
cat("  Count matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")

cat("Loading matched metadata...\n")
meta <- readRDS(file.path(RDIR, "meta_matched.rds"))
setDT(meta)

cat("Loading QC report...\n")
qc <- fread(file.path(INT, "qc/sample_qc_report.csv"))
meta <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]
cat("  QC-passing samples:", nrow(meta), "\n")

# PRJNA512027 permanently removed (Script 00 defensive filter)
meta <- meta[dataset != "PRJNA512027"]
cat("  After PRJNA512027 removal:", nrow(meta), "\n")

# ============================================================
#  Load canonical NAFL-vs-NASH dream results
# ============================================================
cat("\nLoading canonical NAFL-vs-NASH dream results...\n")
canon_file <- file.path(DSIG, "nafl_vs_nash_dream.csv")
if (!file.exists(canon_file)) stop("nafl_vs_nash_dream.csv not found — run Script 13 first")
canon <- fread(canon_file)
cat("  Canonical genes:", nrow(canon), "\n")
cat("  Canonical DEGs (padj<0.05):", sum(canon$adj.P.Val < 0.05), "\n")

# ============================================================
#  Sample inventory: canonical vs strict
# ============================================================
cat("\n========== SAMPLE INVENTORY ==========\n")

# Canonical: NAFL + NASH + Borderline (Borderline grouped with NASH)
meta_canon <- meta[diagnosis_harmonized %in% c("NAFL", "NASH", "Borderline")]
meta_canon[, nafl_nash := fifelse(diagnosis_harmonized == "NAFL", "NAFL", "NASH")]

# Strict / Borderline-separate: NAS < 3 (NAFL) vs NAS >= 5 (NASH)
# Both arms exclude NAS 3-4, so the sample set is identical
meta_strict <- meta[diagnosis_harmonized %in% c("NAFL", "NASH")]
meta_strict[, nafl_nash := diagnosis_harmonized]

# Also include NAFL/NASH from label-only datasets (GSE126848, GSE167523)
# These have no NAS score -> diagnosis_harmonized comes from condition labels
# They are already "NAFL" or "NASH" (not "Borderline"), so they appear in both

n_borderline <- nrow(meta[diagnosis_harmonized == "Borderline"])
n_canon_nafl <- sum(meta_canon$nafl_nash == "NAFL")
n_canon_nash <- sum(meta_canon$nafl_nash == "NASH")
n_strict_nafl <- sum(meta_strict$nafl_nash == "NAFL")
n_strict_nash <- sum(meta_strict$nafl_nash == "NASH")

cat(sprintf("Borderline (NAS 3-4) samples: %d\n", n_borderline))
cat(sprintf("\nCanonical (Script 13): NAFL=%d, NASH=%d, total=%d\n",
  n_canon_nafl, n_canon_nash, nrow(meta_canon)))
cat(sprintf("Strict (NAS 3-4 excluded): NAFL=%d, NASH=%d, total=%d\n",
  n_strict_nafl, n_strict_nash, nrow(meta_strict)))
cat(sprintf("Samples dropped: %d (%.1f%%)\n",
  n_borderline, 100 * n_borderline / nrow(meta_canon)))

cat("\nBorderline samples by dataset:\n")
print(meta[diagnosis_harmonized == "Borderline", .N, by = dataset][order(dataset)])

cat("\nCanonical grouping by dataset:\n")
print(meta_canon[, .N, by = .(dataset, nafl_nash)][order(dataset, nafl_nash)])

cat("\nStrict grouping by dataset:\n")
print(meta_strict[, .N, by = .(dataset, nafl_nash)][order(dataset, nafl_nash)])

# Check minimum group sizes per dataset
strict_sizes <- meta_strict[, .(
  n_nafl = sum(nafl_nash == "NAFL"),
  n_nash = sum(nafl_nash == "NASH")
), by = dataset]
cat("\nPer-dataset group sizes (strict):\n")
print(strict_sizes)

# Flag datasets that lose a group entirely
lost_datasets <- strict_sizes[n_nafl < 3 | n_nash < 3, dataset]
if (length(lost_datasets) > 0) {
  cat("WARNING: Datasets with <3 in a group after exclusion:", paste(lost_datasets, collapse = ", "), "\n")
}

# ============================================================
#  Run dream NAFL-vs-NASH on strict subset
# ============================================================
cat("\n========== DREAM: STRICT (NAS 3-4 excluded) ==========\n")

# Factor levels
meta_strict[, nafl_nash := factor(nafl_nash, levels = c("NAFL", "NASH"))]

# Align counts
idx_strict <- colnames(counts) %in% meta_strict$sample_id
dge_strict <- DGEList(counts = counts[, idx_strict])

# Attach metadata
dge_strict$samples <- cbind(dge_strict$samples,
  meta_strict[match(colnames(dge_strict), meta_strict$sample_id),
    .(nafl_nash, dataset, sex, inferred_sex, age)])

# Dataset subbatch (PRJNA512027 already removed)
dge_strict$samples$dataset_subbatch <- as.character(dge_strict$samples$dataset)

# Sex handling (same as Script 13)
dge_strict$samples$sex_for_model <- dge_strict$samples$sex
na_sex <- is.na(dge_strict$samples$sex_for_model) | dge_strict$samples$sex_for_model == ""
dge_strict$samples$sex_for_model[na_sex] <- as.character(dge_strict$samples$inferred_sex[na_sex])
dge_strict$samples$sex_for_model <- factor(dge_strict$samples$sex_for_model)

# TMM normalization
dge_strict <- calcNormFactors(dge_strict, method = "TMM")

# Filter genes
keep_strict <- filterByExpr(dge_strict, group = dge_strict$samples$nafl_nash)
dge_strict <- dge_strict[keep_strict, , keep.lib.sizes = FALSE]

cat("Strict dream samples:", ncol(dge_strict), "\n")
cat("Strict dream genes:", nrow(dge_strict), "\n")
cat("Batch levels:", paste(sort(unique(dge_strict$samples$dataset_subbatch)), collapse = ", "), "\n")
cat("Sex levels:", paste(levels(dge_strict$samples$sex_for_model), collapse = ", "), "\n")

# Dream formula (same as Script 13)
form <- ~ nafl_nash + sex_for_model + (1 | dataset_subbatch)
cat("Formula:", deparse(form), "\n")

cat("Running voomWithDreamWeights...\n")
vobjStrict <- voomWithDreamWeights(dge_strict, form, dge_strict$samples, BPPARAM = param)

cat("Running dream...\n")
fitStrict <- dream(vobjStrict, form, dge_strict$samples, BPPARAM = param)
# Do NOT call eBayes() after dream() — dream already applies moderated t-stats

strict_tt <- topTable(fitStrict, coef = "nafl_nashNASH", number = Inf, sort.by = "none")
strict_tt$gene <- rownames(strict_tt)
strict_dt <- as.data.table(strict_tt)

n_strict_degs <- sum(strict_dt$adj.P.Val < 0.05)
n_strict_up <- sum(strict_dt$adj.P.Val < 0.05 & strict_dt$logFC > 0)
n_strict_down <- sum(strict_dt$adj.P.Val < 0.05 & strict_dt$logFC < 0)
cat(sprintf("Strict DEGs (padj<0.05): %d (Up: %d, Down: %d)\n",
  n_strict_degs, n_strict_up, n_strict_down))
cat(sprintf("Mean LFC: %.4f\n", mean(strict_dt$logFC)))

fwrite(strict_dt, file.path(OUTDIR, "dream_strict_nas5_only.csv"))
cat("Saved: dream_strict_nas5_only.csv\n")

# ============================================================
#  Concordance: canonical vs strict
# ============================================================
cat("\n========== CONCORDANCE ==========\n")

# Merge on shared genes
merged <- merge(
  canon[, .(gene, canon_lfc = logFC, canon_padj = adj.P.Val)],
  strict_dt[, .(gene, strict_lfc = logFC, strict_padj = adj.P.Val)],
  by = "gene"
)
cat("Shared genes:", nrow(merged), "\n")

# Spearman rho (all shared genes)
rho_all <- cor(merged$canon_lfc, merged$strict_lfc, method = "spearman")
rho_pearson <- cor(merged$canon_lfc, merged$strict_lfc, method = "pearson")
cat(sprintf("Spearman rho (all genes): %.4f\n", rho_all))
cat(sprintf("Pearson r (all genes): %.4f\n", rho_pearson))

# Direction concordance (all genes)
dir_concord_all <- mean(sign(merged$canon_lfc) == sign(merged$strict_lfc))
cat(sprintf("Direction concordance (all genes): %.4f\n", dir_concord_all))

# DEG sets
canon_degs_03 <- canon[adj.P.Val < 0.05 & abs(logFC) > 0.3, gene]
strict_degs_03 <- strict_dt[adj.P.Val < 0.05 & abs(logFC) > 0.3, gene]
canon_degs_05 <- canon[adj.P.Val < 0.05 & abs(logFC) > 0.5, gene]
strict_degs_05 <- strict_dt[adj.P.Val < 0.05 & abs(logFC) > 0.5, gene]
canon_degs_nolfc <- canon[adj.P.Val < 0.05, gene]
strict_degs_nolfc <- strict_dt[adj.P.Val < 0.05, gene]

jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

j_03 <- jaccard(canon_degs_03, strict_degs_03)
j_05 <- jaccard(canon_degs_05, strict_degs_05)
j_nolfc <- jaccard(canon_degs_nolfc, strict_degs_nolfc)

cat(sprintf("\nJaccard (padj<0.05, no LFC filter): %.4f (intersect %d / union %d)\n",
  j_nolfc, length(intersect(canon_degs_nolfc, strict_degs_nolfc)),
  length(union(canon_degs_nolfc, strict_degs_nolfc))))
cat(sprintf("Jaccard (padj<0.05, |LFC|>0.3): %.4f (intersect %d / union %d)\n",
  j_03, length(intersect(canon_degs_03, strict_degs_03)),
  length(union(canon_degs_03, strict_degs_03))))
cat(sprintf("Jaccard (padj<0.05, |LFC|>0.5): %.4f (intersect %d / union %d)\n",
  j_05, length(intersect(canon_degs_05, strict_degs_05)),
  length(union(canon_degs_05, strict_degs_05))))

# Spearman rho on shared DEGs
shared_degs <- intersect(canon_degs_nolfc, strict_degs_nolfc)
if (length(shared_degs) > 10) {
  rho_degs <- cor(merged[gene %in% shared_degs, canon_lfc],
                  merged[gene %in% shared_degs, strict_lfc],
                  method = "spearman")
  cat(sprintf("Spearman rho (shared DEGs): %.4f\n", rho_degs))
} else {
  rho_degs <- NA_real_
  cat("Too few shared DEGs for correlation\n")
}

# Direction concordance on shared DEGs
if (length(shared_degs) > 0) {
  dir_shared <- merged[gene %in% shared_degs,
    mean(sign(canon_lfc) == sign(strict_lfc))]
  cat(sprintf("Direction concordance (shared DEGs): %.4f\n", dir_shared))
} else {
  dir_shared <- NA_real_
}

# Genes gained/lost
lost_03 <- setdiff(canon_degs_03, strict_degs_03)
gained_03 <- setdiff(strict_degs_03, canon_degs_03)
cat(sprintf("\nDEGs lost when excluding NAS 3-4 (|LFC|>0.3): %d\n", length(lost_03)))
cat(sprintf("DEGs gained when excluding NAS 3-4 (|LFC|>0.3): %d\n", length(gained_03)))

lost_05 <- setdiff(canon_degs_05, strict_degs_05)
gained_05 <- setdiff(strict_degs_05, canon_degs_05)
cat(sprintf("DEGs lost when excluding NAS 3-4 (|LFC|>0.5): %d\n", length(lost_05)))
cat(sprintf("DEGs gained when excluding NAS 3-4 (|LFC|>0.5): %d\n", length(gained_05)))

# ============================================================
#  Per-gene comparison table
# ============================================================
per_gene <- merge(
  canon[, .(gene, canon_lfc = logFC, canon_padj = adj.P.Val, canon_t = t)],
  strict_dt[, .(gene, strict_lfc = logFC, strict_padj = adj.P.Val, strict_t = t)],
  by = "gene", all = TRUE
)
per_gene[, `:=`(
  lfc_diff = strict_lfc - canon_lfc,
  canon_sig = canon_padj < 0.05,
  strict_sig = strict_padj < 0.05,
  direction_match = sign(canon_lfc) == sign(strict_lfc)
)]
per_gene[, status := fcase(
  canon_sig & strict_sig & direction_match, "shared",
  canon_sig & !strict_sig, "lost",
  !canon_sig & strict_sig, "gained",
  canon_sig & strict_sig & !direction_match, "flipped",
  default = "NS_both"
)]

fwrite(per_gene, file.path(OUTDIR, "per_gene_canon_vs_strict.csv"))
cat("\nSaved: per_gene_canon_vs_strict.csv\n")

# Status breakdown
cat("\nPer-gene status:\n")
print(per_gene[, .N, by = status][order(-N)])

# ============================================================
#  NAS 3-4 characterization: what do these patients look like?
# ============================================================
cat("\n========== NAS 3-4 CHARACTERIZATION ==========\n")

borderline_meta <- meta[diagnosis_harmonized == "Borderline"]
cat("Borderline (NAS 3-4) patients:", nrow(borderline_meta), "\n")

if (nrow(borderline_meta) > 0) {
  # NAS score breakdown
  cat("\nNAS score distribution:\n")
  print(borderline_meta[, .N, by = nas_score][order(nas_score)])

  # Fibrosis stage distribution
  cat("\nFibrosis stage distribution (Borderline):\n")
  print(borderline_meta[, .N, by = fibrosis_stage][order(fibrosis_stage)])

  # Compare fibrosis distributions
  cat("\nFibrosis stage by diagnosis group:\n")
  fib_tab <- meta[diagnosis_harmonized %in% c("NAFL", "Borderline", "NASH"),
    .N, by = .(diagnosis_harmonized, fibrosis_stage)][order(diagnosis_harmonized, fibrosis_stage)]
  print(fib_tab)

  fwrite(fib_tab, file.path(OUTDIR, "fibrosis_by_diagnosis_group.csv"))
}

# ============================================================
#  Volcano comparison plot
# ============================================================
cat("\nGenerating comparison plots...\n")

pdf(file.path(OUTDIR, "nas_borderline_sensitivity.pdf"), width = 14, height = 10)

# Panel 1: LFC scatter (canonical vs strict)
p1 <- ggplot(merged, aes(x = canon_lfc, y = strict_lfc)) +
  geom_point(size = 0.3, alpha = 0.2, color = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = min(merged$canon_lfc) + 0.5, y = max(merged$strict_lfc) - 0.2,
    label = sprintf("rho = %.3f\nr = %.3f", rho_all, rho_pearson),
    hjust = 0, size = 4) +
  labs(
    title = "LFC: Canonical (NAS 3-4 with NASH) vs Strict (NAS >=5 only)",
    x = "Canonical logFC", y = "Strict logFC"
  ) +
  theme_minimal(base_size = 11)
print(p1)

# Panel 2: -log10(padj) scatter
p2 <- ggplot(merged, aes(x = -log10(canon_padj), y = -log10(strict_padj))) +
  geom_point(size = 0.3, alpha = 0.2, color = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "-log10(padj): Canonical vs Strict",
    x = "Canonical -log10(padj)", y = "Strict -log10(padj)"
  ) +
  theme_minimal(base_size = 11)
print(p2)

# Panel 3: Side-by-side volcanos
vdat <- rbind(
  merged[, .(gene, lfc = canon_lfc, padj = canon_padj, arm = "Canonical\n(NAS 3-4 grouped with NASH)")],
  merged[, .(gene, lfc = strict_lfc, padj = strict_padj, arm = "Strict\n(NAS >= 5 only)")]
)
vdat[, sig := padj < 0.05 & abs(lfc) > 0.5]

p3 <- ggplot(vdat, aes(x = lfc, y = -log10(padj))) +
  geom_point(aes(color = sig), size = 0.3, alpha = 0.3) +
  scale_color_manual(values = c("grey70", "#E91E63"), guide = "none") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  facet_wrap(~ arm) +
  labs(
    title = "NAFL vs NASH Volcanos: NAS Borderline Sensitivity",
    x = "log2 Fold Change (NASH vs NAFL)", y = "-log10(padj)"
  ) +
  theme_minimal(base_size = 11)
print(p3)

# Panel 4: Highlight gained/lost DEGs
if (nrow(per_gene[status != "NS_both"]) > 0) {
  pdat <- per_gene[status != "NS_both"]
  p4 <- ggplot(pdat, aes(x = canon_lfc, y = strict_lfc, color = status)) +
    geom_point(size = 0.8, alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c(
      "shared" = "#4CAF50", "lost" = "#F44336",
      "gained" = "#2196F3", "flipped" = "#FF9800"
    )) +
    labs(
      title = "DEG status: canonical vs strict",
      subtitle = sprintf("shared=%d, lost=%d, gained=%d, flipped=%d",
        sum(pdat$status == "shared"), sum(pdat$status == "lost"),
        sum(pdat$status == "gained"), sum(pdat$status == "flipped")),
      x = "Canonical logFC", y = "Strict logFC", color = "Status"
    ) +
    theme_minimal(base_size = 11)
  print(p4)
}

dev.off()
cat("Saved: nas_borderline_sensitivity.pdf\n")

# ============================================================
#  Summary metrics table
# ============================================================
n_canon_degs_05_count <- sum(canon$adj.P.Val < 0.05 & abs(canon$logFC) > 0.5)
n_strict_degs_05_count <- sum(strict_dt$adj.P.Val < 0.05 & abs(strict_dt$logFC) > 0.5)
n_canon_degs_03_count <- sum(canon$adj.P.Val < 0.05 & abs(canon$logFC) > 0.3)
n_strict_degs_03_count <- sum(strict_dt$adj.P.Val < 0.05 & abs(strict_dt$logFC) > 0.3)

summary_metrics <- data.table(
  metric = c(
    "canonical_n_samples", "strict_n_samples", "n_borderline_excluded",
    "pct_samples_excluded",
    "canonical_n_degs_padj05", "strict_n_degs_padj05",
    "canonical_n_degs_lfc03", "strict_n_degs_lfc03",
    "canonical_n_degs_lfc05", "strict_n_degs_lfc05",
    "spearman_rho_all", "pearson_r_all",
    "spearman_rho_shared_degs",
    "direction_concordance_all", "direction_concordance_shared_degs",
    "jaccard_padj05_nolfc", "jaccard_padj05_lfc03", "jaccard_padj05_lfc05",
    "n_degs_lost_lfc03", "n_degs_gained_lfc03",
    "n_degs_lost_lfc05", "n_degs_gained_lfc05"
  ),
  value = c(
    nrow(meta_canon), nrow(meta_strict), n_borderline,
    round(100 * n_borderline / nrow(meta_canon), 1),
    sum(canon$adj.P.Val < 0.05), sum(strict_dt$adj.P.Val < 0.05),
    n_canon_degs_03_count, n_strict_degs_03_count,
    n_canon_degs_05_count, n_strict_degs_05_count,
    round(rho_all, 4), round(rho_pearson, 4),
    round(rho_degs, 4),
    round(dir_concord_all, 4), round(dir_shared, 4),
    round(j_nolfc, 4), round(j_03, 4), round(j_05, 4),
    length(lost_03), length(gained_03),
    length(lost_05), length(gained_05)
  )
)

fwrite(summary_metrics, file.path(OUTDIR, "summary_metrics.csv"))
cat("\nSaved: summary_metrics.csv\n")

# ============================================================
#  Generate REPORT.md
# ============================================================
cat("\nGenerating REPORT.md...\n")

# Verdict logic: Jaccard > 0.70 and rho > 0.90 = PASS
verdict_pass <- j_03 >= 0.70 && rho_all >= 0.90
verdict_str <- if (verdict_pass) {
  sprintf("**PASS**: Jaccard = %.3f (>= 0.70), Spearman rho = %.3f (>= 0.90). Excluding NAS 3-4 does not materially change the NAFL-vs-NASH DE results. The borderline grouping decision is robust.",
    j_03, rho_all)
} else {
  sprintf("**FAIL**: Jaccard = %.3f and/or Spearman rho = %.3f below threshold. NAS 3-4 classification materially affects DE results; consider reporting both analyses.",
    j_03, rho_all)
}

report <- sprintf(
'# NAS Borderline Sensitivity Analysis [R3 P1-1]

**Date:** %s

## Design

- **Reviewer concern**: NAS 3-4 patients are grouped with NASH, but clinically
  NAS 3-4 is indeterminate (Kleiner et al. 2005), not definite NASH.
- **Current classification** (Script 00, line 298):
  NAS < 3 = NAFL; NAS 3-4 = Borderline (grouped with NASH); NAS >= 5 = NASH
- **Sensitivity arms**:
  - **Arm 1 (Strict)**: NAS >= 5 only as NASH; NAS 3-4 excluded
  - **Arm 2 (Borderline separate)**: NAS 3-4 as own category, excluded from both groups
  - Both arms yield identical NAFL-vs-NASH samples (NAS < 3 vs NAS >= 5)
- **Model**: `~ nafl_nash + sex_for_model + (1|dataset_subbatch)` (dream)
- **Reference**: Canonical Script 13 (NAS 3-4 grouped with NASH, N=%d)
- **Test**: Strict (NAS 3-4 excluded, N=%d)

## Sample changes

| Group | Canonical | Strict | Change |
|-------|-----------|--------|--------|
| NAFL | %d | %d | 0 |
| NASH | %d | %d | -%d |
| **Total** | **%d** | **%d** | **-%d** |

- Borderline (NAS 3-4) samples excluded: **%d** (%.1f%% of canonical)
- Datasets with Borderline samples: %s

## DEG counts (padj < 0.05)

| Threshold | Canonical | Strict | Delta |
|-----------|-----------|--------|-------|
| No LFC filter | %d | %d | %+d |
| \\|logFC\\| > 0.3 | %d | %d | %+d |
| \\|logFC\\| > 0.5 | %d | %d | %+d |

## Concordance metrics

- **Spearman rho (all genes)**: %.4f
- **Pearson r (all genes)**: %.4f
- **Spearman rho (shared DEGs)**: %s
- **Direction concordance (all genes)**: %.4f
- **Direction concordance (shared DEGs)**: %s
- **Jaccard (padj<0.05, no LFC)**: %.4f (intersect %d / union %d)
- **Jaccard (padj<0.05, |LFC|>0.3)**: %.4f (intersect %d / union %d)
- **Jaccard (padj<0.05, |LFC|>0.5)**: %.4f (intersect %d / union %d)

## Genes affected (padj < 0.05, |LFC| > 0.3)

- DEGs lost when excluding NAS 3-4: %d
- DEGs gained when excluding NAS 3-4: %d
- DEGs shared: %d
- DEGs with flipped direction: %d

## Verdict

%s

## Files

- `dream_strict_nas5_only.csv` -- full per-gene strict dream results
- `per_gene_canon_vs_strict.csv` -- matched logFC/padj/status for canonical vs strict
- `summary_metrics.csv` -- all concordance metrics
- `fibrosis_by_diagnosis_group.csv` -- fibrosis distribution across NAFL/Borderline/NASH
- `nas_borderline_sensitivity.pdf` -- comparison plots (4 panels)
',
  as.character(Sys.time()),
  nrow(meta_canon), nrow(meta_strict),
  n_canon_nafl, n_strict_nafl,
  n_canon_nash, n_strict_nash, n_borderline,
  nrow(meta_canon), nrow(meta_strict), n_borderline,
  n_borderline, 100 * n_borderline / nrow(meta_canon),
  paste(meta[diagnosis_harmonized == "Borderline", unique(dataset)], collapse = ", "),
  # DEG counts
  sum(canon$adj.P.Val < 0.05), sum(strict_dt$adj.P.Val < 0.05),
  as.integer(sum(strict_dt$adj.P.Val < 0.05) - sum(canon$adj.P.Val < 0.05)),
  n_canon_degs_03_count, n_strict_degs_03_count,
  as.integer(n_strict_degs_03_count - n_canon_degs_03_count),
  n_canon_degs_05_count, n_strict_degs_05_count,
  as.integer(n_strict_degs_05_count - n_canon_degs_05_count),
  # Concordance
  rho_all, rho_pearson,
  if (is.na(rho_degs)) "N/A (too few)" else sprintf("%.4f", rho_degs),
  dir_concord_all,
  if (is.na(dir_shared)) "N/A" else sprintf("%.4f", dir_shared),
  j_nolfc, length(intersect(canon_degs_nolfc, strict_degs_nolfc)),
  length(union(canon_degs_nolfc, strict_degs_nolfc)),
  j_03, length(intersect(canon_degs_03, strict_degs_03)),
  length(union(canon_degs_03, strict_degs_03)),
  j_05, length(intersect(canon_degs_05, strict_degs_05)),
  length(union(canon_degs_05, strict_degs_05)),
  # Genes affected
  length(lost_03), length(gained_03),
  length(intersect(canon_degs_03, strict_degs_03)),
  sum(per_gene$status == "flipped", na.rm = TRUE),
  verdict_str
)

writeLines(report, file.path(OUTDIR, "REPORT.md"))
cat("Saved: REPORT.md\n")

cat("\n=== Script 13b complete ===\n")
cat("Finished:", as.character(Sys.time()), "\n")
