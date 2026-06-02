#!/usr/bin/env Rscript
# 131_nas_component_ordinal_dream.R
# ---------------------------------------------------------------------------
# Ordinal dream models for individual NAS components.
#
# For each NAS component (treated as continuous ordinal covariate):
#   C7a: Steatosis grade (0-3)
#   C7b: Lobular inflammation grade (0-3)
#   C7c: Ballooning grade (0-2)
#
# The coefficient of interest is the slope per unit increase in the component
# score, capturing dose-response transcriptomic changes.
#
# Formula: ~ component_score + sex_covar + (1 | dataset)
# Falls back to fixed effects (~ component_score + sex_covar + dataset) when
# only 1 dataset has valid component annotations.
#
# NOTE: Currently only GSE130970 (N=76) provides individual NAS component
# grades. All other datasets have NAS total scores but lack the component
# breakdown. If future datasets add component grades, the random-intercept
# model will engage automatically (>= 2 datasets).
#
# NOTE: This script correctly avoids calling eBayes() after dream().
# dream() already computes moderated t-statistics via the Satterthwaite
# approximation.
#
# Output: results/progression/c7{a,b,c}_{component}_ordinal_dream.csv
#         results/progression/c7_nas_components_ordinal_combined.csv
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
  library(BiocParallel)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
ODIR <- file.path(INT, "results/progression")
dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)

# --- Parallel setup ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
cat("Using", ncpus, "CPU cores\n")
BPPARAM <- if (ncpus > 1) MulticoreParam(ncpus, progressbar = TRUE) else SerialParam()

# --- Load data ---
cat("Loading data...\n")
counts   <- readRDS(file.path(INT, "results/integration/merged_counts_raw.rds"))
meta     <- readRDS(file.path(INT, "results/integration/meta_matched.rds"))
qc       <- fread(file.path(INT, "qc/sample_qc_report.csv"))
meta     <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

# Gene annotations (dream_results.csv lacks annotations; use unified_disease_signatures)
annot_file <- file.path(INT, "results/disease_signatures/unified_disease_signatures.csv")
if (file.exists(annot_file)) {
  gene_annot <- fread(annot_file,
    select = c("gene", "symbol", "gene_type", "mouse_gene_id", "mouse_symbol"))
  gene_annot <- unique(gene_annot, by = "gene")
} else {
  gene_annot <- data.table(gene = character(), symbol = character(),
    gene_type = character(), mouse_gene_id = character(), mouse_symbol = character())
}

# Load modeling metadata for NAS component columns
modeling_meta <- fread(file.path(INT, "results/staging_classifier/modeling_metadata.csv"))

# Merge NAS component columns into meta
merge_cols <- c("sample_id", "steatosis_grade", "lobular_inflammation_grade", "ballooning_grade")
merge_cols <- intersect(merge_cols, names(modeling_meta))
meta <- merge(meta, modeling_meta[, ..merge_cols],
  by = "sample_id", all.x = TRUE)

# Prepare sex covariate
meta[, sex_covar := inferred_sex]
na_sex <- is.na(meta$sex_covar) | meta$sex_covar == ""
if (any(na_sex) && "sex" %in% names(meta)) {
  meta$sex_covar[na_sex] <- meta$sex[na_sex]
}
meta[, sex_covar := factor(sex_covar)]

cat(sprintf("Total QC-passing samples: %d\n", nrow(meta)))

# ============================================================
# Define NAS components
# ============================================================
components <- list(
  list(
    id       = "c7a",
    name     = "steatosis",
    col      = "steatosis_grade",
    label    = "Steatosis Grade (0-3)"
  ),
  list(
    id       = "c7b",
    name     = "inflammation",
    col      = "lobular_inflammation_grade",
    label    = "Lobular Inflammation Grade (0-3)"
  ),
  list(
    id       = "c7c",
    name     = "ballooning",
    col      = "ballooning_grade",
    label    = "Ballooning Grade (0-2)"
  )
)

# ============================================================
# Helper: run dream for ordinal continuous covariate
# ============================================================
run_dream_ordinal <- function(comp, meta_all, counts_all, gene_annot, BPPARAM) {
  id    <- comp$id
  name  <- comp$name
  col   <- comp$col
  label <- comp$label

  cat(sprintf("\n%s\n", strrep("=", 60)))
  cat(sprintf("  CONTRAST %s: %s (ordinal)\n", toupper(id), label))
  cat(sprintf("%s\n", strrep("=", 60)))

  # --- Subset to samples with valid component score ---
  # Sentinel value -1 = missing; also exclude NA
  meta_sub <- copy(meta_all[!is.na(get(col)) & get(col) >= 0])

  if (nrow(meta_sub) < 10) {
    cat(sprintf("  SKIPPED: Only %d samples with valid %s values (need >= 10).\n",
      nrow(meta_sub), col))
    return(NULL)
  }

  # Treat component score as numeric (continuous ordinal)
  meta_sub[, component_score := as.numeric(get(col))]

  cat(sprintf("  Samples with valid %s: %d\n", col, nrow(meta_sub)))
  cat("  Score distribution:\n")
  print(table(meta_sub$component_score))

  # Per-dataset sample counts
  cat("  Per-dataset distribution:\n")
  for (ds in sort(unique(meta_sub$dataset))) {
    sub <- meta_sub[dataset == ds]
    n_per_score <- table(sub$component_score)
    cat(sprintf("    %s (N=%d): %s\n", ds, nrow(sub),
      paste(paste0(names(n_per_score), "=", n_per_score), collapse = ", ")))
  }

  # --- Filter datasets with < 2 samples total ---
  ds_counts <- meta_sub[, .N, by = dataset]
  valid_ds <- ds_counts[N >= 2, dataset]

  # Also require at least 2 distinct score values per dataset
  ds_diversity <- meta_sub[dataset %in% valid_ds,
    .(n_levels = uniqueN(component_score)), by = dataset]
  valid_ds <- ds_diversity[n_levels >= 2, dataset]

  if (length(valid_ds) == 0) {
    cat("  SKIPPED: No datasets with >= 2 samples and >= 2 distinct score levels.\n")
    return(NULL)
  }
  meta_sub <- meta_sub[dataset %in% valid_ds]
  cat(sprintf("  Valid datasets: %d (%s)\n", length(valid_ds),
    paste(valid_ds, collapse = ", ")))

  # --- Build DGEList from raw counts ---
  keep_samples <- intersect(meta_sub$sample_id, colnames(counts_all))
  meta_sub <- meta_sub[sample_id %in% keep_samples]
  dge <- DGEList(counts = counts_all[, keep_samples])

  # Attach metadata
  m_ordered <- meta_sub[match(colnames(dge), meta_sub$sample_id)]
  dge$samples <- cbind(dge$samples, m_ordered[, .(dataset, sex_covar, component_score)])
  dge$samples$dataset   <- factor(dge$samples$dataset)
  dge$samples$sex_covar <- factor(dge$samples$sex_covar)

  # TMM normalization
  dge <- calcNormFactors(dge, method = "TMM")

  # Filter low-expression genes
  keep_genes <- filterByExpr(dge)
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

  cat(sprintf("  Samples: %d\n", ncol(dge)))
  cat(sprintf("  Genes after filterByExpr: %d\n", nrow(dge)))

  # --- Determine formula (random vs fixed) ---
  has_random <- length(valid_ds) >= 2
  if (has_random) {
    formula_str <- "~ component_score + sex_covar + (1 | dataset)"
  } else {
    # Single dataset — drop the dataset term entirely
    formula_str <- "~ component_score + sex_covar"
  }
  form <- as.formula(formula_str)
  cat(sprintf("  Formula: %s\n", formula_str))

  # --- Run dream or limma ---
  if (has_random) {
    vobj <- suppressWarnings(voomWithDreamWeights(dge, form, dge$samples, BPPARAM = BPPARAM))
    fit  <- suppressWarnings(dream(vobj, form, dge$samples, BPPARAM = BPPARAM))
    # NOTE: do NOT call eBayes() after dream() — dream() already computes
    # moderated t-statistics via Satterthwaite approximation.
  } else {
    design <- model.matrix(form, data = dge$samples)
    vobj <- voom(dge, design, plot = FALSE)
    fit  <- lmFit(vobj, design)
    fit  <- eBayes(fit)
  }

  # --- Extract results for component_score coefficient ---
  res <- topTable(fit, coef = "component_score", number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  res_dt <- as.data.table(res)

  # Add gene annotations
  res_dt <- merge(res_dt, gene_annot, by = "gene", all.x = TRUE)

  # Rename for consistency with Script 05 output convention
  setnames(res_dt, "adj.P.Val", "padj", skip_absent = TRUE)

  # Add component metadata
  res_dt[, component := name]
  res_dt[, contrast_id := id]
  res_dt[, n_samples := ncol(dge)]
  res_dt[, n_datasets := length(valid_ds)]

  # --- Summary ---
  n_sig_01 <- sum(res_dt$padj < 0.1, na.rm = TRUE)
  n_sig_05 <- sum(res_dt$padj < 0.05, na.rm = TRUE)
  n_up   <- sum(res_dt$padj < 0.1 & res_dt$logFC > 0, na.rm = TRUE)
  n_down <- sum(res_dt$padj < 0.1 & res_dt$logFC < 0, na.rm = TRUE)
  cat(sprintf("  DEGs (padj<0.1): %d (Up: %d, Down: %d)\n", n_sig_01, n_up, n_down))
  cat(sprintf("  DEGs (padj<0.05): %d\n", n_sig_05))

  # logFC here = slope per unit increase in component score
  if (n_sig_01 > 0) {
    top5 <- head(res_dt[order(padj)], 5)
    cat("  Top 5 by padj:\n")
    for (i in seq_len(nrow(top5))) {
      cat(sprintf("    %s (slope=%.3f, padj=%.2e)\n",
        top5$symbol[i], top5$logFC[i], top5$padj[i]))
    }
  }

  # --- Save ---
  out_file <- file.path(ODIR, sprintf("%s_%s_ordinal_dream.csv", id, name))
  fwrite(res_dt, out_file)
  cat(sprintf("  Saved: %s\n", basename(out_file)))

  gc()
  return(res_dt)
}

# ============================================================
# Run each NAS component
# ============================================================
all_results <- list()
for (comp in components) {
  res <- run_dream_ordinal(comp, meta, counts, gene_annot, BPPARAM)
  all_results[[comp$id]] <- res
}

# ============================================================
# Combined output
# ============================================================
valid_results <- Filter(Negate(is.null), all_results)
if (length(valid_results) > 0) {
  combined <- rbindlist(valid_results, fill = TRUE)
  fwrite(combined, file.path(ODIR, "c7_nas_components_ordinal_combined.csv"))
  cat(sprintf("\nSaved combined: c7_nas_components_ordinal_combined.csv (%d rows)\n",
    nrow(combined)))
}

# ============================================================
# Summary across all components
# ============================================================
cat(sprintf("\n%s\n", strrep("=", 60)))
cat("  NAS COMPONENT ORDINAL SUMMARY\n")
cat(sprintf("%s\n", strrep("=", 60)))

summary_table <- rbindlist(lapply(names(all_results), function(nm) {
  r <- all_results[[nm]]
  comp <- components[[which(sapply(components, `[[`, "id") == nm)]]
  if (is.null(r)) {
    cat(sprintf("  %s (%s): SKIPPED\n", nm, comp$label))
    return(data.table(contrast_id = nm, component = comp$name,
      label = comp$label, n_samples = 0, n_datasets = 0,
      n_genes = 0, n_deg_01 = 0, n_deg_05 = 0, n_up = 0, n_down = 0))
  }
  data.table(
    contrast_id = nm,
    component   = comp$name,
    label       = comp$label,
    n_samples   = r$n_samples[1],
    n_datasets  = r$n_datasets[1],
    n_genes     = nrow(r),
    n_deg_01    = sum(r$padj < 0.1, na.rm = TRUE),
    n_deg_05    = sum(r$padj < 0.05, na.rm = TRUE),
    n_up        = sum(r$padj < 0.1 & r$logFC > 0, na.rm = TRUE),
    n_down      = sum(r$padj < 0.1 & r$logFC < 0, na.rm = TRUE)
  )
}))

print(summary_table)
fwrite(summary_table, file.path(ODIR, "c7_nas_component_ordinal_summary.csv"))
cat("Saved: c7_nas_component_ordinal_summary.csv\n")

# ============================================================
# C17: Fibrosis Ordinal (continuous, all datasets)
# Dose-response model: fibrosis_stage 0-4 as numeric covariate.
# Uses all 8 datasets with fibrosis staging (~1,128 samples).
# Fixes the eBayes-after-dream bug in Script 14's version and
# includes more datasets (8 vs 3 in the original).
# ============================================================
cat(sprintf("\n%s\n", strrep("=", 60)))
cat("  C17: FIBROSIS ORDINAL (CONTINUOUS, ALL DATASETS)\n")
cat(sprintf("%s\n", strrep("=", 60)))

# Load unified metadata for fibrosis staging (broader than modeling_metadata)
meta_fib <- readRDS(file.path(INT, "results/integration/meta_matched.rds"))
qc_fib <- fread(file.path(INT, "qc/sample_qc_report.csv"))
meta_fib <- meta_fib[sample_id %in% qc_fib[pass_technical == TRUE, sample_id]]

# Merge fibrosis from modeling_metadata
mm <- fread(file.path(INT, "results/staging_classifier/modeling_metadata.csv"),
  select = c("sample_id", "fibrosis_stage"))
meta_fib <- merge(meta_fib, mm, by = "sample_id", all.x = TRUE, suffixes = c("", ".mm"))
if ("fibrosis_stage.mm" %in% names(meta_fib)) {
  na_mask <- is.na(meta_fib$fibrosis_stage) | meta_fib$fibrosis_stage == ""
  if (any(na_mask)) meta_fib$fibrosis_stage[na_mask] <- meta_fib$fibrosis_stage.mm[na_mask]
  meta_fib[, fibrosis_stage.mm := NULL]
}

meta_c17 <- meta_fib[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta_c17[, fib_numeric := as.numeric(fibrosis_stage)]

# Prepare sex covariate
meta_c17[, sex_covar := inferred_sex]
na_sex <- is.na(meta_c17$sex_covar) | meta_c17$sex_covar == ""
if (any(na_sex) && "sex" %in% names(meta_c17)) {
  meta_c17$sex_covar[na_sex] <- meta_c17$sex[na_sex]
}
meta_c17[, sex_covar := factor(sex_covar)]

cat(sprintf("  Samples with F0-F4: %d\n", nrow(meta_c17)))
cat("  Per stage:\n")
print(table(meta_c17$fibrosis_stage))
cat("  Per dataset:\n")
print(meta_c17[, .N, by = dataset][order(dataset)])

# Filter datasets with >= 2 distinct fibrosis levels
ds_levels <- meta_c17[, .(n_levels = length(unique(fibrosis_stage))), by = dataset]
valid_ds <- ds_levels[n_levels >= 2, dataset]
meta_c17 <- meta_c17[dataset %in% valid_ds]
cat(sprintf("  Valid datasets (>= 2 F levels): %d\n", length(valid_ds)))

# Build DGEList
keep_s <- intersect(meta_c17$sample_id, colnames(counts))
meta_c17 <- meta_c17[sample_id %in% keep_s]
dge_c17 <- DGEList(counts = counts[, keep_s])
m_c17 <- meta_c17[match(colnames(dge_c17), meta_c17$sample_id)]
dge_c17$samples <- cbind(dge_c17$samples, m_c17[, .(dataset, sex_covar, fib_numeric)])
dge_c17$samples$dataset <- factor(dge_c17$samples$dataset)
dge_c17$samples$sex_covar <- factor(dge_c17$samples$sex_covar)

dge_c17 <- calcNormFactors(dge_c17, method = "TMM")
keep_genes <- filterByExpr(dge_c17)
dge_c17 <- dge_c17[keep_genes, , keep.lib.sizes = FALSE]

cat(sprintf("  Samples: %d, Genes: %d\n", ncol(dge_c17), nrow(dge_c17)))

# Dream with fibrosis as continuous ordinal
form_c17 <- ~ fib_numeric + sex_covar + (1 | dataset)
cat(sprintf("  Formula: %s\n", deparse(form_c17)))

vobj_c17 <- suppressWarnings(voomWithDreamWeights(dge_c17, form_c17, dge_c17$samples, BPPARAM = BPPARAM))
fit_c17 <- suppressWarnings(dream(vobj_c17, form_c17, dge_c17$samples, BPPARAM = BPPARAM))

res_c17 <- topTable(fit_c17, coef = "fib_numeric", number = Inf, sort.by = "none")
res_c17$gene <- rownames(res_c17)
c17_dt <- as.data.table(res_c17)
setnames(c17_dt, "adj.P.Val", "padj", skip_absent = TRUE)
c17_dt <- merge(c17_dt, gene_annot, by = "gene", all.x = TRUE)
c17_dt[, contrast_id := "c17"]
c17_dt[, component := "fibrosis"]
c17_dt[, n_samples := ncol(dge_c17)]
c17_dt[, n_datasets := length(valid_ds)]

n17_01 <- sum(c17_dt$padj < 0.1, na.rm = TRUE)
n17_05 <- sum(c17_dt$padj < 0.05, na.rm = TRUE)
n17_up <- sum(c17_dt$padj < 0.1 & c17_dt$logFC > 0, na.rm = TRUE)
n17_dn <- sum(c17_dt$padj < 0.1 & c17_dt$logFC < 0, na.rm = TRUE)
cat(sprintf("  DEGs (padj<0.1): %d (Up: %d, Down: %d)\n", n17_01, n17_up, n17_dn))
cat(sprintf("  DEGs (padj<0.05): %d\n", n17_05))

fwrite(c17_dt, file.path(ODIR, "c17_fibrosis_ordinal_dream.csv"))
cat("  Saved: c17_fibrosis_ordinal_dream.csv\n")

cat("\n=== Script 131 complete ===\n")
