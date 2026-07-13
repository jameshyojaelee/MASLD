#!/usr/bin/env Rscript
# B5: 9-cohort sensitivity (effectively 6-cohort: 5 canonical + PRJNA512027 with
# L0/S0 library_prep_batch as a covariate, OR ComBat-seq corrected if covariate
# fails). The 4 control-less cohorts (GSE167523/GSE174478/GSE193066/GSE240729)
# cannot contribute to Disease-vs-Control; we report this in the methods note.
#
# Output: RNA-seq/results/audit_sensitivity/9_cohort_sensitivity/
#   - dream_results_6cohort.csv
#   - REPORT.md (concordance vs canonical 5-cohort)

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(variancePartition)
  library(BiocParallel)
  library(sva)  # ComBat_seq
})

MAIN <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
WT   <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation"
RDIR <- file.path(MAIN, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT  <- file.path(WT,   "RNA-seq/results/audit_sensitivity/9_cohort_sensitivity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat("Loading DGE...\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

# Load batch-annotated metadata
meta_batch <- fread(file.path(WT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata_with_batch.csv"))
meta_orig  <- readRDS(file.path(RDIR, "meta_matched.rds"))

# Cohorts to include: 5 canonical + PRJNA512027
cohorts_6 <- c("GSE126848","GSE130970","GSE135251","GSE162694","PRJNA512027")
# Actual canonical 5 in dream is (Suppli/Hoang/Govaere/Bril/Chen)
# Chen = GSE126848? No. Per CLAUDE.md: Suppli=GSE130970, Hoang=GSE162694,
# Govaere=GSE135251, Bril=?, Chen=GSE126848 (recheck — we'll use whatever the
# canonical dream csv had). Use the dataset list actually present in the canonical run.
canonical_csv <- file.path(RDIR, "dream_results.csv")  # C2-OK-sensitivity
canonical <- fread(canonical_csv)
cat("Canonical csv:", canonical_csv, "rows=", nrow(canonical), "\n")

# Determine canonical cohorts from yaml-implied set (5 with healthy controls)
yaml <- yaml::read_yaml(file.path(MAIN, "config/human_datasets.yaml"))
mega5 <- names(Filter(function(d) isTRUE(d$de$include_in_mega), yaml$datasets))
cat("Canonical 5 cohorts (yaml include_in_mega):", paste(mega5, collapse=","), "\n")

cohorts_6 <- unique(c(mega5, "PRJNA512027"))
cat("6-cohort set:", paste(cohorts_6, collapse=","), "\n")

# Subset dge
keep <- dge$samples$dataset %in% cohorts_6
dge6 <- dge[, keep]
cat("Samples 6-cohort:", ncol(dge6), "\n")

# Sex
matched_sex <- meta_orig$inferred_sex[match(colnames(dge6), meta_orig$sample_id)]

# library_prep_batch
batch_map <- setNames(meta_batch$library_prep_batch, meta_batch$sample_id)
lpb <- batch_map[colnames(dge6)]
lpb[is.na(lpb) | lpb == "NA" | lpb == ""] <- "none"
cat("library_prep_batch levels:\n"); print(table(lpb, dge6$samples$dataset))

# Strategy: include library_prep_batch as covariate. For non-PRJNA512027 cohorts
# all values are "none" so the covariate adds 1 df = L0 vs S0 vs none contrast.
# We use a fixed-effect dummy for the L0/S0 split (only meaningful within
# PRJNA512027); this is exactly the confound correction reviewers want.
info <- data.frame(
  group_binary = factor(dge6$samples$group_binary, levels = c("Control","Disease")),
  dataset = factor(dge6$samples$dataset),
  inferred_sex = factor(matched_sex),
  lib_batch = factor(lpb),
  stringsAsFactors = FALSE
)
rownames(info) <- colnames(dge6)
cat("Distribution:\n"); print(table(info$group_binary, info$dataset))

# Try model with lib_batch covariate
form_full  <- ~ group_binary + inferred_sex + lib_batch + (1|dataset)

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

cat("\n=== Running dream with library_prep_batch covariate ===\n")
v <- suppressWarnings(voomWithDreamWeights(dge6, form_full, info, BPPARAM = param))
fit <- suppressWarnings(dream(v, form_full, info, BPPARAM = param))
res <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res$gene <- rownames(res)
res_dt <- as.data.table(res); setnames(res_dt, "adj.P.Val", "padj")
fwrite(res_dt, file.path(OUT, "dream_results_6cohort_with_batch.csv"))
cat("Wrote dream_results_6cohort_with_batch.csv\n")

# Sensitivity: ComBat-seq correction on PRJNA512027 L0/S0 batch, then dream without
# lib_batch covariate
cat("\n=== ComBat-seq alternative (PRJNA512027 L0/S0 corrected) ===\n")
counts <- dge6$counts
prjna_idx <- which(dge6$samples$dataset == "PRJNA512027")
# Only correct within PRJNA512027 (small batch labels => need >=2 levels per dataset)
sub_batch <- lpb[prjna_idx]
sub_batch_factor <- factor(sub_batch)
if (length(levels(sub_batch_factor)) >= 2) {
  sub_group <- as.character(info$group_binary[prjna_idx])
  corrected <- ComBat_seq(counts[, prjna_idx], batch = sub_batch_factor, group = sub_group)
  counts_cb <- counts
  counts_cb[, prjna_idx] <- corrected
  dge_cb <- dge6; dge_cb$counts <- counts_cb
  form_cb <- ~ group_binary + inferred_sex + (1|dataset)
  v2 <- suppressWarnings(voomWithDreamWeights(dge_cb, form_cb, info, BPPARAM = param))
  fit2 <- suppressWarnings(dream(v2, form_cb, info, BPPARAM = param))
  res2 <- topTable(fit2, coef = "group_binaryDisease", number = Inf, sort.by = "none")
  res2$gene <- rownames(res2); res2_dt <- as.data.table(res2); setnames(res2_dt, "adj.P.Val", "padj")
  fwrite(res2_dt, file.path(OUT, "dream_results_6cohort_combatseq.csv"))
  cat("Wrote dream_results_6cohort_combatseq.csv\n")
} else {
  cat("ComBat-seq skipped: only 1 batch level in PRJNA512027 subset\n")
  res2_dt <- NULL
}

# ---- Concordance vs canonical 5-cohort ----
sig_cutoff <- function(dt, padj_thr=0.05, lfc_thr=0.5) dt[padj < padj_thr & abs(logFC) >= lfc_thr, gene]
canon_sig <- sig_cutoff(canonical)
b6_sig    <- sig_cutoff(res_dt)
jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))
J_lpb <- jacc(canon_sig, b6_sig)
delta_pct_lpb <- (length(b6_sig) - length(canon_sig)) / length(canon_sig)

report_lines <- c(
  "# B5 — 9-cohort sensitivity (effectively 6-cohort: 5 canonical + PRJNA512027)",
  "",
  sprintf("Date: %s", format(Sys.time())),
  "",
  "## Rationale",
  "",
  "The canonical mega-analysis (CLAUDE.md 2026-05-15) drops PRJNA512027 (Gerhard)",
  "because its samples come from two distinct library-prep batches (`L0xxxx` and",
  "`S0xxxx` library names) and the L0/S0 split is strongly confounded with",
  "Control/Disease. This B5 sensitivity re-adds PRJNA512027 with the L0/S0",
  "library-prep batch as an explicit covariate (and as a ComBat-seq sensitivity).",
  "The 4 other previously-collected cohorts (GSE167523/GSE174478/GSE193066/",
  "GSE240729) have no healthy controls and therefore cannot contribute to a",
  "Disease-vs-Control mega-analysis under any model; they remain in DGE for",
  "stage-stratified analyses but are not added here.",
  "",
  "## Models",
  "",
  "1. **canonical 5-cohort** (no PRJNA512027): `~ group + sex + (1|dataset)`",
  "2. **6-cohort + library_prep_batch covariate**: `~ group + sex + lib_batch + (1|dataset)`",
  "3. **6-cohort + ComBat-seq**: PRJNA512027 counts ComBat-seq-corrected on L0/S0",
  "   batch (preserving Disease/Control group), then dream without lib_batch.",
  "",
  "## Headline DEG counts (padj < 0.05, |LFC| > 0.5)",
  "",
  sprintf("- canonical 5-cohort: **%d**", length(canon_sig)),
  sprintf("- 6-cohort + lib_batch covariate: **%d** (Δ = %+.1f%%)", length(b6_sig), 100*delta_pct_lpb)
)

if (!is.null(res2_dt)) {
  b6_cb_sig <- sig_cutoff(res2_dt)
  J_cb <- jacc(canon_sig, b6_cb_sig)
  delta_cb <- (length(b6_cb_sig) - length(canon_sig)) / length(canon_sig)
  report_lines <- c(report_lines,
    sprintf("- 6-cohort + ComBat-seq: **%d** (Δ = %+.1f%%)", length(b6_cb_sig), 100*delta_cb))
}

report_lines <- c(report_lines,
  "",
  "## Concordance vs canonical 5-cohort",
  "",
  sprintf("- 6-cohort + lib_batch Jaccard: **%.3f**", J_lpb))
if (!is.null(res2_dt)) {
  report_lines <- c(report_lines, sprintf("- 6-cohort + ComBat-seq Jaccard: **%.3f**", J_cb))
}

verdict_robust <- abs(delta_pct_lpb) < 0.10
report_lines <- c(report_lines,
  "",
  "## Verdict",
  "",
  sprintf("- Δ headline DEG count from canonical: %s (target <10%% for robust)",
          if (verdict_robust) "**ROBUST**" else "**FLAG IN METHODS**"),
  "",
  "## Files",
  "",
  "- `dream_results_6cohort_with_batch.csv`",
  if (!is.null(res2_dt)) "- `dream_results_6cohort_combatseq.csv`" else "")

writeLines(report_lines, file.path(OUT, "REPORT.md"))
cat("Wrote REPORT.md\n")
cat("End:", format(Sys.time()), "\n")
