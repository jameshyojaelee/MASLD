#!/usr/bin/env Rscript
# Enhanced dream mega-analysis on kallisto gene counts
# Extends 04_dream_kallisto.R with SE, df.total, isSingular columns (Phase 1.4)
# Reads from main repo RNA-seq/results/kallisto/ (post-worktree-merge)
# Output: RNA-seq/Human/.../dream_results.csv (CANONICAL — replaces STAR)

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

PROJECT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
BASE <- file.path(PROJECT, "RNA-seq/Human/Patient_Cohorts")
RDIR <- file.path(BASE, "analysis/integration/results/integration")
KALL <- file.path(PROJECT, "RNA-seq/results/kallisto")
OUT  <- file.path(RDIR, "dream_results.csv")

cat("[", as.character(Sys.time()), "] Enhanced kallisto dream (Phase 1.4)\n")

# Load kallisto counts
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

rownames(cnts) <- sub("\\..*$", "", rownames(cnts))

# Pull canonical sample metadata
cat("[", as.character(Sys.time()), "] loading canonical metadata\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))

yaml_path <- file.path(PROJECT, "config/human_datasets.yaml")
ycfg <- yaml::read_yaml(yaml_path)$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("  mega cohorts:", paste(mega_cohorts, collapse=", "), "\n")

ds_samples <- dge$samples
ds_samples$sample_id <- rownames(ds_samples)
keep_canonical <- ds_samples$dataset %in% mega_cohorts
canonical_ids <- ds_samples$sample_id[keep_canonical]
common_ids <- intersect(canonical_ids, colnames(cnts))
cat("  canonical mega samples:", length(canonical_ids),
    " | kallisto samples:", ncol(cnts),
    " | common:", length(common_ids), "\n")

cnts <- cnts[, common_ids]
ds_samples <- ds_samples[match(common_ids, ds_samples$sample_id), ]

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

# Phase 1.4: Extract SE, df.total, isSingular
cat("[", as.character(Sys.time()), "] Extracting SE, df.total, isSingular...\n")
coef_idx <- which(colnames(fit$coefficients) == "group_binaryDisease")
se_vec <- fit$stdev.unscaled[, coef_idx] * fit$sigma.post
res_dt[, SE := se_vec[match(gene, rownames(fit$coefficients))]]
if (!is.null(fit$df.total)) {
  res_dt[, df.total := fit$df.total[match(gene, rownames(fit$coefficients))]]
} else {
  res_dt[, df.total := NA_real_]
}
vc <- try(variancePartition::extractVarPart(fit), silent = TRUE)
if (!inherits(vc, "try-error") && "dataset" %in% colnames(vc)) {
  dataset_var <- vc[match(res_dt$gene, rownames(vc)), "dataset"]
  res_dt[, isSingular := dataset_var < 1e-8]
} else {
  res_dt[, isSingular := NA]
}

sig_lfc <- res_dt[padj < 0.05 & abs(logFC) >= 0.5]
cat("\n===== KALLISTO DREAM RESULTS (ENHANCED) =====\n")
cat("Total genes tested:", nrow(res_dt), "\n")
cat("DEGs padj<0.05 & |logFC|>=0.5:", nrow(sig_lfc),
    " (up:", sum(sig_lfc$logFC > 0), "down:", sum(sig_lfc$logFC < 0), ")\n")
cat("DEGs padj<0.1:", sum(res_dt$padj < 0.1, na.rm = TRUE), "\n")
cat("Columns:", paste(names(res_dt), collapse=", "), "\n")
cat("isSingular count:", sum(res_dt$isSingular, na.rm = TRUE), "/", nrow(res_dt), "\n")

fwrite(res_dt, OUT)
cat("[", as.character(Sys.time()), "] Saved:", OUT, "\n")
