#!/usr/bin/env Rscript
# 06a_bariatric_sensitivity.R
# ---------------------------------------------------------------------------
# R3 P0-2: Bariatric sensitivity arm
# GSE162694 (Bril) consists entirely of bariatric surgery wedge biopsies
# (BMI > 35, different tissue procurement). Run dream on 4 cohorts EXCLUDING
# GSE162694 and compare with the 5-cohort canonical kallisto dream results.
#
# Model: ~ group_binary + inferred_sex + (1|dataset)
# Input: RNA-seq/results/kallisto/all_cohorts_gene_counts.tsv.gz
# Reference: RNA-seq/results/kallisto/dream_results_kallisto.csv
# Output: RNA-seq/results/audit_sensitivity/bariatric_sensitivity/
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(reformulas); library(lme4); library(data.table); library(edgeR); library(yaml)
})
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) try({
    unlockBinding(fn, ns_lme4)
    assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
    lockBinding(fn, ns_lme4)
  }, silent = TRUE)
}
suppressPackageStartupMessages({ library(variancePartition); library(BiocParallel) })

PROJECT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR    <- file.path(PROJECT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
KALL    <- file.path(PROJECT, "RNA-seq/results/kallisto")
OUTDIR  <- file.path(PROJECT, "RNA-seq/results/audit_sensitivity/bariatric_sensitivity")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

EXCL_DATASET <- "GSE162694"  # Bril bariatric cohort

# ---- Load kallisto counts ---------------------------------------------------
cat("[", as.character(Sys.time()), "] Loading kallisto counts\n")
kc <- fread(file.path(KALL, "all_cohorts_gene_counts.tsv.gz"))
gene_ids <- kc$gene_id
kc[, gene_id := NULL]
cnts <- as.matrix(kc)
rownames(cnts) <- gene_ids
storage.mode(cnts) <- "double"
cnts[is.na(cnts)] <- 0
cnts <- round(cnts)
rownames(cnts) <- sub("\\..*$", "", rownames(cnts))
cat("  kallisto matrix:", nrow(cnts), "genes x", ncol(cnts), "samples\n")

# ---- Canonical metadata from merged_dge.rds --------------------------------
cat("[", as.character(Sys.time()), "] Loading canonical metadata\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))

yaml_path <- file.path(PROJECT, "config/human_datasets.yaml")
ycfg <- yaml::read_yaml(yaml_path)$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))

# Exclude GSE162694 (bariatric)
mega_4 <- setdiff(mega_cohorts, EXCL_DATASET)
cat("  5-cohort canonical:", paste(mega_cohorts, collapse = ", "), "\n")
cat("  4-cohort (excl bariatric):", paste(mega_4, collapse = ", "), "\n")

# ---- Subset samples ---------------------------------------------------------
ds_samples <- dge$samples
ds_samples$sample_id <- rownames(ds_samples)
keep_4 <- ds_samples$dataset %in% mega_4
canonical_4 <- ds_samples$sample_id[keep_4]
common_4 <- intersect(canonical_4, colnames(cnts))
cat("  4-cohort samples:", length(common_4), "\n")

cnts_4 <- cnts[, common_4]
ds_sub <- ds_samples[match(common_4, ds_samples$sample_id), ]

# ---- Build DGEList + filter --------------------------------------------------
group <- factor(ds_sub$group_binary, levels = c("Control", "Disease"))
y <- DGEList(counts = cnts_4, samples = data.frame(
  sample_id    = common_4,
  dataset      = ds_sub$dataset,
  group_binary = ds_sub$group_binary,
  stringsAsFactors = FALSE))
keep_g <- filterByExpr(y, group = group)
y <- y[keep_g, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
cat("  After filterByExpr+TMM:", nrow(y), "genes x", ncol(y), "samples\n")

matched_sex <- meta_new$inferred_sex[match(common_4, meta_new$sample_id)]
info <- data.frame(
  group_binary = factor(ds_sub$group_binary, levels = c("Control", "Disease")),
  dataset      = droplevels(factor(ds_sub$dataset)),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE)
rownames(info) <- common_4

cat("\n4-cohort Group x Dataset:\n"); print(table(info$group_binary, info$dataset))
cat("Sex distribution:\n"); print(table(info$inferred_sex, useNA = "always"))

# ---- Dream (4-cohort, RE) ---------------------------------------------------
form <- ~ group_binary + inferred_sex + (1|dataset)
cat("\nFormula:", deparse(form), "\n")

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

cat("[", as.character(Sys.time()), "] voomWithDreamWeights...\n")
v <- suppressWarnings(voomWithDreamWeights(y, form, info, BPPARAM = param))

cat("[", as.character(Sys.time()), "] dream()...\n")
fit <- suppressWarnings(dream(v, form, info, BPPARAM = param))

res <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res$gene <- rownames(res)
res_dt <- as.data.table(res)
setnames(res_dt, "adj.P.Val", "padj")

cat("\n===== 4-COHORT (EXCL BARIATRIC) RESULTS =====\n")
cat("Total genes tested:", nrow(res_dt), "\n")
cat("DEGs padj<0.05 & |logFC|>=0.3:", nrow(res_dt[padj < 0.05 & abs(logFC) >= 0.3]), "\n")
cat("DEGs padj<0.05 & |logFC|>=0.5:", nrow(res_dt[padj < 0.05 & abs(logFC) >= 0.5]), "\n")

fwrite(res_dt, file.path(OUTDIR, "dream_4cohort_nobariatric_kallisto.csv"))
cat("Saved 4-cohort results\n")

# ---- Load 5-cohort canonical reference (kallisto) ----------------------------
cat("\n[", as.character(Sys.time()), "] Loading 5-cohort canonical kallisto reference\n")
ref <- fread(file.path(KALL, "dream_results_kallisto.csv"))
setnames(ref, "adj.P.Val", "padj", skip_absent = TRUE)

# ---- Comparison metrics ------------------------------------------------------
ref_sub <- ref[, .(gene, logFC_5c = logFC, padj_5c = padj)]
new_sub <- res_dt[, .(gene, logFC_4c = logFC, padj_4c = padj)]
comp <- merge(ref_sub, new_sub, by = "gene")
cat("Genes in common (5c vs 4c):", nrow(comp), "\n")

# DEG sets at padj<0.05, |logFC|>0.3
deg_5c <- comp[padj_5c < 0.05 & abs(logFC_5c) > 0.3, gene]
deg_4c <- comp[padj_4c < 0.05 & abs(logFC_4c) > 0.3, gene]

# Also at padj<0.05, |logFC|>0.5
deg_5c_t1 <- comp[padj_5c < 0.05 & abs(logFC_5c) > 0.5, gene]
deg_4c_t1 <- comp[padj_4c < 0.05 & abs(logFC_4c) > 0.5, gene]

jacc <- function(a, b) length(intersect(a, b)) / max(length(union(a, b)), 1)

# Spearman rho
rho_all <- cor(comp$logFC_5c, comp$logFC_4c, method = "spearman")
rho_sig <- tryCatch({
  shared_sig <- intersect(deg_5c, deg_4c)
  sub <- comp[gene %in% shared_sig]
  cor(sub$logFC_5c, sub$logFC_4c, method = "spearman")
}, error = function(e) NA_real_)

# Direction concordance among shared DEGs
shared <- intersect(deg_5c, deg_4c)
dir_conc <- if (length(shared) > 0) {
  sub <- comp[gene %in% shared]
  mean(sign(sub$logFC_5c) == sign(sub$logFC_4c))
} else NA_real_

# Genes lost/gained
lost_genes  <- setdiff(deg_5c, deg_4c)
gained_genes <- setdiff(deg_4c, deg_5c)

# Save per-gene comparison
fwrite(comp, file.path(OUTDIR, "per_gene_5c_vs_4c.csv"))

# ---- Write REPORT.md --------------------------------------------------------
report_path <- file.path(OUTDIR, "REPORT.md")
sink(report_path)
cat("# Bariatric Sensitivity Arm [R3 P0-2]\n\n")
cat("**Date:**", as.character(Sys.time()), "\n\n")

cat("## Design\n\n")
cat("- **Question**: Does excluding GSE162694 (Bril, N=142 bariatric surgery wedge biopsies,\n")
cat("  BMI>35, different tissue procurement) materially change the mega-analysis?\n")
cat("- **Quantifier**: kallisto (tximport gene-level counts)\n")
cat("- **Model**: `~ group_binary + inferred_sex + (1|dataset)` (dream)\n")
cat("- **Reference**: 5-cohort canonical kallisto dream (k=5, N=", length(canonical_4) + sum(ds_samples$dataset == EXCL_DATASET & ds_samples$sample_id %in% colnames(cnts)), ")\n", sep = "")
cat("- **Test**: 4-cohort (excl GSE162694) kallisto dream (k=4, N=", ncol(y), ")\n\n", sep = "")

cat("## Cohort composition\n\n")
cat("| Cohort | Role | N |\n")
cat("|--------|------|---|\n")
for (coh in mega_cohorts) {
  role <- if (coh == EXCL_DATASET) "EXCLUDED (bariatric)" else "Included"
  n <- sum(ds_samples$dataset == coh & ds_samples$sample_id %in% colnames(cnts))
  cat(sprintf("| %s | %s | %d |\n", coh, role, n))
}

cat("\n## Results\n\n")
cat("### DEG counts (padj < 0.05)\n\n")
cat("| Threshold | 5-cohort | 4-cohort (excl bariatric) |\n")
cat("|-----------|----------|---------------------------|\n")
cat(sprintf("| \\|logFC\\| > 0.3 | %d | %d |\n", length(deg_5c), length(deg_4c)))
cat(sprintf("| \\|logFC\\| > 0.5 | %d | %d |\n", length(deg_5c_t1), length(deg_4c_t1)))

cat("\n### Concordance metrics\n\n")
cat(sprintf("- **Spearman rho (all genes)**: %.4f\n", rho_all))
cat(sprintf("- **Spearman rho (shared DEGs)**: %.4f\n", rho_sig))
cat(sprintf("- **Direction concordance (shared DEGs)**: %.4f\n", dir_conc))
cat(sprintf("- **Jaccard (padj<0.05, |LFC|>0.3)**: %.4f (intersect %d / union %d)\n",
            jacc(deg_5c, deg_4c), length(intersect(deg_5c, deg_4c)), length(union(deg_5c, deg_4c))))
cat(sprintf("- **Jaccard (padj<0.05, |LFC|>0.5)**: %.4f (intersect %d / union %d)\n",
            jacc(deg_5c_t1, deg_4c_t1), length(intersect(deg_5c_t1, deg_4c_t1)), length(union(deg_5c_t1, deg_4c_t1))))

cat("\n### Genes affected\n\n")
cat(sprintf("- DEGs lost when excluding bariatric: %d\n", length(lost_genes)))
cat(sprintf("- DEGs gained when excluding bariatric: %d\n", length(gained_genes)))

cat("\n### Verdict\n\n")
j03 <- jacc(deg_5c, deg_4c)
if (j03 >= 0.85) {
  cat(sprintf("**PASS**: Jaccard = %.3f (>= 0.85 threshold). Inclusion of GSE162694 is defensible;\n", j03))
  cat("bariatric tissue procurement does not materially alter the mega-analysis.\n")
} else if (j03 >= 0.70) {
  cat(sprintf("**CONDITIONAL**: Jaccard = %.3f (0.70-0.85 range). Moderate perturbation;\n", j03))
  cat("recommend Methods note documenting bariatric tissue procurement and sensitivity.\n")
} else {
  cat(sprintf("**FAIL**: Jaccard = %.3f (< 0.70). GSE162694 substantially drives the mega-analysis;\n", j03))
  cat("consider excluding or running parallel analyses.\n")
}

cat("\n## Files\n\n")
cat("- `dream_4cohort_nobariatric_kallisto.csv` — full per-gene results\n")
cat("- `per_gene_5c_vs_4c.csv` — matched logFC/padj for 5c vs 4c\n")

sink()
cat("\nWrote:", report_path, "\n")
cat("[", as.character(Sys.time()), "] Done.\n")
