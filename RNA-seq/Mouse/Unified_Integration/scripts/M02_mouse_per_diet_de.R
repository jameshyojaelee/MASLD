#!/usr/bin/env Rscript
# M02_mouse_per_diet_de.R
# ---------------------------------------------------------------------------
# Per-diet-type differential expression using limma-voom.
# Comparisons built dynamically from config/mouse_datasets.yaml — each
# dataset's diet_model field determines which comparison group it belongs to.
# Adding a new diet model = editing the YAML, not this script.
# Input:  merged_counts_raw.rds, meta_matched.rds, sample_qc_report.csv,
#         config/mouse_datasets.yaml
# Output: results/per_diet/<diet_model>_de_results.csv
# ---------------------------------------------------------------------------

# ---- Seed pinning (T2.4, 2026-04-22) -----
set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(yaml)
})

PROJECT_ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
MOUSE  <- file.path(PROJECT_ROOT, "RNA-seq/Mouse")
INT    <- file.path(MOUSE, "Unified_Integration")
RDIR   <- file.path(INT, "results")
DEDIR  <- file.path(RDIR, "per_diet")
dir.create(DEDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== M02: Per-Diet DE (limma-voom) ===\n\n")

# --- Load data ---
merged <- readRDS(file.path(RDIR, "merged_counts_raw.rds"))
meta   <- readRDS(file.path(RDIR, "meta_matched.rds"))
qc     <- fread(file.path(INT, "qc/sample_qc_report.csv"))

# --- Load tximport gene-length offsets if available (Kallisto pipeline) ---
# When counts come from Kallisto via tximport(countsFromAbundance="no"),
# gene-level effective lengths vary per sample (isoform-usage dependent).
# Passing log(length) as voom offsets corrects for condition-dependent isoform
# switching that would otherwise bias fold-change estimates.
lengths_file <- file.path(RDIR, "merged_gene_lengths.rds")
HAS_TX_OFFSETS <- file.exists(lengths_file)
if (HAS_TX_OFFSETS) {
  gene_lengths_all <- readRDS(lengths_file)
  cat("Loaded tximport gene-length matrix:", nrow(gene_lengths_all), "genes x",
      ncol(gene_lengths_all), "samples\n")
  cat("  -> Will apply log(length) offsets to voom\n\n")
} else {
  cat("No tximport gene-length matrix found (featureCounts pipeline) — no offsets\n\n")
}

# Filter to QC-passing samples
pass <- qc[pass_qc == TRUE, sample_id]
merged <- merged[, colnames(merged) %in% pass]
meta   <- meta[sample_id %in% pass]
cat("QC-passing samples:", ncol(merged), "/", nrow(qc), "\n\n")

# --- Build comparisons list from YAML config ---
# Each diet model type gets its own DE analysis.
# Controls are identified from the same datasets as the disease samples.
cfg_path <- file.path(PROJECT_ROOT, "config/mouse_datasets.yaml")
if (!file.exists(cfg_path)) stop("Mouse config not found: ", cfg_path)
mouse_cfg <- yaml::read_yaml(cfg_path)
cat("Loaded config:", cfg_path, "\n")

# Group datasets by diet_model from the YAML.
# CDAHFD_FPC entries register in both CDAHFD and FPC groups.
comparisons <- list()
for (ds_id in names(mouse_cfg$datasets)) {
  ds <- mouse_cfg$datasets[[ds_id]]
  models <- ds$diet_model
  # Expand dual-model entries (e.g. "CDAHFD_FPC" -> c("CDAHFD", "FPC"))
  if (grepl("_", models) && !models %in% c("AMLN_ob")) {
    models <- unlist(strsplit(models, "_"))
  }
  for (dm in models) {
    if (is.null(comparisons[[dm]])) {
      comparisons[[dm]] <- list(datasets = character(0))
    }
    comparisons[[dm]]$datasets <- c(comparisons[[dm]]$datasets, ds_id)
  }
}

# Derive has_batch (multiple datasets) and has_sex (determined at runtime from
# metadata, so default TRUE — the DE loop below checks actual sex availability).
for (dm in names(comparisons)) {
  comparisons[[dm]]$has_batch <- length(comparisons[[dm]]$datasets) > 1
  comparisons[[dm]]$has_sex   <- TRUE
}

# Log the comparisons built from YAML
cat("\nComparisons built from YAML:\n")
for (dm in names(comparisons)) {
  cat("  ", dm, ": datasets =", paste(comparisons[[dm]]$datasets, collapse = ", "),
      " has_batch =", comparisons[[dm]]$has_batch, "\n")
}
cat("\n")

results_summary <- data.table()

for (diet_name in names(comparisons)) {
  comp <- comparisons[[diet_name]]
  
  cat("============================================================\n")
  cat("  ", diet_name, "\n")
  cat("============================================================\n")
  
  # Select samples for this comparison
  idx <- meta$dataset %in% comp$datasets &
         (meta$diet_model == diet_name | meta$group_binary == "Control")
  sub_meta <- meta[idx]
  sub_counts <- merged[, sub_meta$sample_id]
  
  cat("  Samples:", nrow(sub_meta), "\n")
  cat("    Disease:", sum(sub_meta$group_binary == "Disease"), "\n")
  cat("    Control:", sum(sub_meta$group_binary == "Control"), "\n")
  
  if (sum(sub_meta$group_binary == "Disease") < 2 |
      sum(sub_meta$group_binary == "Control") < 2) {
    cat("  SKIPPING: insufficient samples\n\n")
    next
  }
  
  # Create DGEList and filter
  # Relaxed filterByExpr thresholds (2026-05-21, perdiet_de_code_review P1):
  # defaults (min.count=10, min.total.count=15) silently drop canonical liver
  # fibrosis markers (Acta2, Timp1, Mmp2, Mmp9, Saa3) in CDAHFD/FPC. Relaxing
  # to min.count=5, min.total.count=10 recovers ~500-1000 biologically real
  # low-expression-in-control genes critical for Cas13 screen positive-control
  # representation.
  dge <- DGEList(counts = sub_counts)
  keep <- filterByExpr(dge, group = sub_meta$group_binary,
                       min.count = 5, min.total.count = 10)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  cat("  Genes after filterByExpr:", nrow(dge), "\n")
  
  # Build design matrix
  group <- factor(sub_meta$group_binary, levels = c("Control", "Disease"))
  
  if (comp$has_batch && comp$has_sex && 
      length(unique(sub_meta[!is.na(sex) & sex != "unspecified", sex])) > 1) {
    # Full model: group + sex + batch
    batch <- factor(sub_meta$dataset)
    sex <- factor(sub_meta$sex)
    # Check if sex has enough variation
    sex_valid <- !is.na(sex) & sex != "unspecified"
    if (sum(sex_valid) > nrow(sub_meta) * 0.5) {
      # Impute missing sex as the majority
      majority_sex <- names(sort(table(sex[sex_valid]), decreasing = TRUE))[1]
      sex[!sex_valid] <- majority_sex
      design <- model.matrix(~ 0 + group + sex + batch)
    } else {
      design <- model.matrix(~ 0 + group + batch)
    }
  } else if (comp$has_batch) {
    batch <- factor(sub_meta$dataset)
    design <- model.matrix(~ 0 + group + batch)
  } else if (comp$has_sex &&
             length(unique(sub_meta[!is.na(sex) & sex != "unspecified", sex])) > 1) {
    sex <- factor(sub_meta$sex)
    design <- model.matrix(~ 0 + group + sex)
  } else {
    design <- model.matrix(~ 0 + group)
  }
  
  colnames(design) <- gsub("^group", "", colnames(design))
  
  cat("  Design cols:", paste(colnames(design), collapse = ", "), "\n")
  
  # Fit limma-voom
  v <- voom(dge, design, plot = FALSE)

  # Apply tximport transcript-length offsets when available.
  # v$offset is on log scale; limma::lmFit() uses it alongside lib-size
  # normalization to correct for per-gene, per-sample effective length
  # differences caused by differential isoform usage across conditions.
  if (HAS_TX_OFFSETS) {
    # Subset length matrix to current genes and samples
    common_g <- intersect(rownames(v), rownames(gene_lengths_all))
    common_s <- intersect(colnames(v), colnames(gene_lengths_all))
    if (length(common_g) < nrow(v) || length(common_s) < ncol(v)) {
      cat("    [offset] Subsetting lengths:", length(common_g), "/", nrow(v),
          "genes,", length(common_s), "/", ncol(v), "samples\n")
    }
    len_sub <- gene_lengths_all[common_g, common_s]
    # Guard against zero or NA lengths (would produce -Inf in log)
    len_sub[is.na(len_sub) | len_sub <= 0] <- 1
    # Align voom object and design to genes/samples present in length matrix
    v <- v[common_g, common_s]
    v$offset <- log(len_sub)
    # If any samples were dropped, realign the design matrix
    if (length(common_s) < nrow(design)) {
      design <- design[match(common_s, rownames(design)), , drop = FALSE]
    }
    cat("    [offset] Applied tximport length offsets\n")
  }

  fit <- lmFit(v, design)
  
  # Contrast: Disease - Control
  contrasts <- makeContrasts(Disease - Control, levels = design)
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)
  
  # Extract results
  res <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  # Add unmoderated SE (stdev.unscaled * sigma) for downstream meta-analysis.
  # fit2$stdev.unscaled[,1] = unscaled SD per gene for contrast coef 1;
  # fit2$sigma = residual SD per gene BEFORE eBayes shrinkage.
  # Their product gives the ordinary-least-squares SE, which is the correct
  # input to metafor::rma() because rma() estimates its own between-study τ².
  # Using logFC/t (moderated SE) would bias τ² downward.
  res$SE_unmoderated <- as.numeric(fit2$stdev.unscaled[, 1]) * fit2$sigma
  res$df_total <- fit2$df.total
  res <- as.data.table(res)
  setnames(res, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "gene",
                   "SE_unmoderated", "df_total"))
  setcolorder(res, "gene")
  
  sig <- res[adj.P.Val < 0.05]
  sig_up <- sig[logFC > 0]
  sig_down <- sig[logFC < 0]
  
  cat("  DEGs (padj < 0.05):", nrow(sig),
      " (Up:", nrow(sig_up), "Down:", nrow(sig_down), ")\n")
  
  # Save
  out_file <- file.path(DEDIR, paste0(diet_name, "_de_results.csv"))
  fwrite(res[order(adj.P.Val)], out_file)
  cat("  Saved:", out_file, "\n\n")
  
  results_summary <- rbindlist(list(results_summary, data.table(
    diet_model = diet_name,
    n_samples = nrow(sub_meta),
    n_disease = sum(sub_meta$group_binary == "Disease"),
    n_control = sum(sub_meta$group_binary == "Control"),
    n_genes = nrow(dge),
    degs_005 = nrow(sig),
    degs_up = nrow(sig_up),
    degs_down = nrow(sig_down)
  )))
}

# --- Summary table ---
cat("============================================================\n")
cat("  CROSS-DIET DE SUMMARY\n")
cat("============================================================\n\n")
print(results_summary)

fwrite(results_summary, file.path(DEDIR, "de_summary.csv"))
cat("\nPer-diet DE complete.\n")
