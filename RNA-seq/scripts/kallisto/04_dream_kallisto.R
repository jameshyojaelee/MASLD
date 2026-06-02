#!/usr/bin/env Rscript
# B1 step 4 - dream mega-analysis on kallisto gene counts
# Mirrors 05_dream_mega_analysis.R (formula, model, voomWithDreamWeights, no
# eBayes, BiocParallel) but reads from RNA-seq/results/kallisto/all_cohorts_gene_counts.tsv.gz
# Output: RNA-seq/results/kallisto/dream_results_kallisto.csv

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

WT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation"
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
RDIR <- file.path(BASE, "analysis/integration/results/integration")
KALL <- file.path(WT, "RNA-seq/results/kallisto")
OUT  <- file.path(KALL, "dream_results_kallisto.csv")

# Load kallisto counts (pooled)
cat("[", as.character(Sys.time()), "] loading kallisto counts\n")
kc <- fread(file.path(KALL, "all_cohorts_gene_counts.tsv.gz"))
gene_ids <- kc$gene_id
kc[, gene_id := NULL]
cnts <- as.matrix(kc)
rownames(cnts) <- gene_ids
storage.mode(cnts) <- "double"
cnts[is.na(cnts)] <- 0
cnts <- round(cnts)
cat("  kallisto matrix:", nrow(cnts), "genes x", ncol(cnts), "samples\n")

# Strip ENSG version suffix to align with canonical DGE gene_ids
rownames(cnts) <- sub("\\..*$", "", rownames(cnts))

# Pull canonical sample metadata from merged_dge.rds
cat("[", as.character(Sys.time()), "] loading canonical metadata\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))

# Mega cohort filter via yaml
yaml_path <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/config/human_datasets.yaml"
ycfg <- yaml::read_yaml(yaml_path)$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("  mega cohorts:", paste(mega_cohorts, collapse=", "), "\n")

# Restrict to samples in canonical mega-set AND present in kallisto matrix
ds_samples <- dge$samples
ds_samples$sample_id <- rownames(ds_samples)
keep_canonical <- ds_samples$dataset %in% mega_cohorts
canonical_ids <- ds_samples$sample_id[keep_canonical]
common_ids <- intersect(canonical_ids, colnames(cnts))
cat("  canonical mega samples:", length(canonical_ids),
    " | kallisto samples:", ncol(cnts),
    " | common (intersect):", length(common_ids), "\n")

cnts <- cnts[, common_ids]
ds_samples <- ds_samples[match(common_ids, ds_samples$sample_id), ]

# Filter genes by edgeR filterByExpr against group_binary
group <- factor(ds_samples$group_binary, levels = c("Control", "Disease"))
y <- DGEList(counts = cnts, samples = data.frame(
  sample_id    = common_ids,
  dataset      = ds_samples$dataset,
  group_binary = ds_samples$group_binary,
  stringsAsFactors = FALSE))
keep_g <- filterByExpr(y, group = group)
y <- y[keep_g, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
cat("  after filterByExpr+TMM:", nrow(y), "genes x", ncol(y), "samples\n")

matched_sex <- meta_new$inferred_sex[match(common_ids, meta_new$sample_id)]
info <- data.frame(
  group_binary = factor(ds_samples$group_binary, levels = c("Control", "Disease")),
  dataset      = factor(ds_samples$dataset),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE)
rownames(info) <- common_ids
cat("Group distribution by dataset:\n"); print(table(info$group_binary, info$dataset))
cat("Sex distribution:\n");              print(table(info$inferred_sex, useNA="always"))

form <- ~ group_binary + inferred_sex + (1|dataset)
cat("\nFormula:", deparse(form), "\n")

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

cat("voomWithDreamWeights...\n")
v <- suppressWarnings(voomWithDreamWeights(y, form, info, BPPARAM = param))

cat("dream()...\n")
fit <- suppressWarnings(dream(v, form, info, BPPARAM = param))
# No eBayes() (matches canonical)

res <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res$gene <- rownames(res)
res_dt <- as.data.table(res)
setnames(res_dt, "adj.P.Val", "padj")

sig_lfc <- res_dt[padj < 0.05 & abs(logFC) >= 0.5]
cat("\n===== KALLISTO DREAM RESULTS =====\n")
cat("Total genes tested:", nrow(res_dt), "\n")
cat("DEGs padj<0.05 & |logFC|>=0.5:", nrow(sig_lfc),
    " (up:", sum(sig_lfc$logFC > 0), "down:", sum(sig_lfc$logFC < 0), ")\n")
cat("DEGs padj<0.1:", sum(res_dt$padj < 0.1, na.rm = TRUE), "\n")

fwrite(res_dt, OUT)
cat("Saved:", OUT, "\n")
