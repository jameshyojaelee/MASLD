#!/usr/bin/env Rscript
# 04_variance_partition.R
# ---------------------------------------------------------------------------
# Variance decomposition with variancePartition to quantify sources of
# variation (condition, sex, dataset/batch).
# Output: results/integration/variance_partition.{csv,pdf}
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
})

# Force injection into lme4 namespace BEFORE loading variancePartition
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
  library(ggplot2)
  library(BiocParallel)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/integration")

# --- Load merged DGE ---
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

# --- Update metadata with inferred sex ---
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge), meta_new$sample_id)]

info <- data.frame(
  sample_id  = colnames(dge),
  dataset    = factor(dge$samples$dataset),
  condition  = factor(dge$samples$condition),
  group_binary = factor(dge$samples$group_binary),
  sex        = factor(matched_sex),
  stringsAsFactors = FALSE
)
rownames(info) <- info$sample_id

cat("Samples:", nrow(info), "\n")
cat("Datasets:", nlevels(info$dataset), "\n")
cat("Sex levels:", levels(info$sex), "\n")

# --- voom weights ---
# Use group_binary as the design for voom so mean-variance weights are estimated
# per condition group rather than with a single intercept. This produces better
# precision weights for the downstream variance partition model.
design_voom <- model.matrix(~ group_binary, data = info)
v <- voom(dge, design_voom)

# --- Variance partition ---
# R1-C13 fix: Use group_binary (Disease/Control, 2 levels) instead of condition
# (6+ levels including NAFL/NASH/Borderline) — the variable of interest should
# match the primary contrast used in dream mega-analysis.
# Also use inferred_sex (from k-means XIST/DDX3Y) for consistency with Script 26.
form <- ~ (1|group_binary) + (1|sex) + (1|dataset)
cat("\nFitting variance partition model...\n")
cat("Formula:", deparse(form), "\n")

# Use serial processing to avoid memory issues
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores for variance partition\n")
register(if (ncpus > 1) MulticoreParam(ncpus) else SerialParam())

# Fix 6: Use ALL genes (or large subset)
# Using all genes
v_sub <- v

vp <- tryCatch({
  suppressWarnings(fitExtractVarPartModel(v_sub, form, info))
}, error = function(e) {
  cat("Full model failed:", conditionMessage(e), "\n")
  cat("Trying simpler model without sex...\n")
  form2 <- ~ (1|group_binary) + (1|dataset)
  suppressWarnings(fitExtractVarPartModel(v_sub, form2, info))
})

cat("\nMedian variance fractions:\n")
print(apply(vp, 2, median))

# --- Plot ---
pdf(file.path(RDIR, "variance_partition.pdf"), width = 8, height = 6)
p <- plotVarPart(vp) +
  ggtitle("Variance Partition: All Genes") +
  theme(plot.title = element_text(size = 14, face = "bold"))
print(p)
dev.off()

# --- Save ---
vp_dt <- as.data.table(vp, keep.rownames = "gene")
fwrite(vp_dt, file.path(RDIR, "variance_partition.csv"))
cat("\nSaved: variance_partition.csv and variance_partition.pdf\n")
