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

# --- Detect a tximport gene-length matrix (Kallisto pipeline) ---
# Kallisto counts arrive via tximport(countsFromAbundance="no"), so gene-level
# effective lengths vary per sample. Length correction must then be baked into the
# counts at the tximport step; it cannot be bolted onto a voom object afterwards
# (see the BG-013 note at the fit below). We only detect the matrix here so the fit
# can refuse rather than silently produce a mislabelled "length-corrected" arm.
lengths_file <- file.path(RDIR, "merged_gene_lengths.rds")
HAS_TX_OFFSETS <- file.exists(lengths_file)
if (HAS_TX_OFFSETS) {
  cat("Gene-length matrix present:", lengths_file, "\n")
  cat("  -> M02 will STOP: limma-voom cannot apply length offsets post hoc (BG-013)\n\n")
} else {
  cat("No tximport gene-length matrix found (featureCounts pipeline) — none needed\n\n")
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
  
  # BG-014 (2026-08-09): the old validity test was `!is.na(sex) & sex != "unspecified"`,
  # which does NOT exclude the empty string. 45 QC-passing samples carry sex == "", so the
  # modal "valid" value WAS the empty string: `majority_sex` resolved to '', the imputation
  # was inert, and sex entered as a 3-level factor ('', female, male) whose two named
  # columns summed exactly to a batch indicator. The design was therefore RANK DEFICIENT
  # (MCD 5/6, HFD 4/5) and limma printed "Coefficients not estimable" on every run.
  #
  # Recoding unlabelled samples to an explicit "unknown" reference level is the
  # missing-indicator method, and it is provably numerically IDENTICAL to the old fit:
  # max |delta logFC| = 0, Spearman = 1.000000, DEG symmetric difference = 0 in BOTH MCD
  # and HFD. It removes the empty-string level and the nonsense "majority" without moving
  # a single number.
  #
  # NOTE: dropping sex altogether is NOT equivalent and must not be done as a "cleanup" --
  # it discards the genuine within-batch sex adjustment and moves results materially
  # (MCD DEG symdiff 269; HFD 4,631 -> 3,543, Spearman 0.968).
  sex_raw   <- sub_meta$sex
  sex_known <- !is.na(sex_raw) & !sex_raw %in% c("", "unspecified")
  sex_levels_known <- sort(unique(sex_raw[sex_known]))
  # Every genuine level needs >= 2 samples for its coefficient to be estimable. This
  # replaces the old ">50% labelled" heuristic, which counted empty strings as labelled.
  sex_usable <- length(sex_levels_known) > 1 &&
    all(table(sex_raw[sex_known]) >= 2)
  make_sex <- function() factor(ifelse(sex_known, sex_raw, "unknown"),
                                levels = c("unknown", sex_levels_known))

  if (comp$has_batch && comp$has_sex && sex_usable) {
    batch <- factor(sub_meta$dataset)
    sex <- make_sex()
    design <- model.matrix(~ 0 + group + sex + batch)
  } else if (comp$has_batch) {
    batch <- factor(sub_meta$dataset)
    design <- model.matrix(~ 0 + group + batch)
  } else if (comp$has_sex && sex_usable) {
    sex <- make_sex()
    design <- model.matrix(~ 0 + group + sex)
  } else {
    design <- model.matrix(~ 0 + group)
  }
  if (qr(design)$rank < ncol(design)) {
    cat("  [warn] design is rank deficient (", qr(design)$rank, "/", ncol(design),
        "); limma will drop aliased coefficients\n", sep = "")
  }
  
  colnames(design) <- gsub("^group", "", colnames(design))
  
  cat("  Design cols:", paste(colnames(design), collapse = ", "), "\n")
  
  # Fit limma-voom
  v <- voom(dge, design, plot = FALSE)

  # BG-013 (2026-08-09): this block previously assigned `v$offset <- log(lengths)`
  # after voom() and claimed the fit was length-corrected. It was a silent no-op.
  # limma::lmFit() consumes only what getEAWP() extracts from an EList -- verified
  # on the installed limma to be exactly (exprs, Amean, weights, design). There is
  # no `offset` member among them, and coefficients are bit-identical with and
  # without an arbitrary log-length offset. voom() likewise reads lib.size *
  # norm.factors, not DGEList$offset. Length correction therefore CANNOT be applied
  # at this point in a limma-voom pipeline; it must be baked into the counts at the
  # tximport step (`countsFromAbundance = "lengthScaledTPM"`, the route the tximport
  # vignette prescribes for limma-voom). Fail loudly rather than mislabel the arm.
  if (HAS_TX_OFFSETS) {
    stop("Refusing to run a mislabelled length-corrected arm.\n",
         "  A gene-length matrix is present (", lengths_file, "), but limma-voom\n",
         "  cannot consume per-observation length offsets: lmFit() reads only\n",
         "  getEAWP() output (exprs, Amean, weights, design) and ignores EList$offset.\n",
         "  Re-run tximport with countsFromAbundance = \"lengthScaledTPM\" so the\n",
         "  effective lengths are carried in the counts, then rerun M02 without a\n",
         "  separate length matrix.", call. = FALSE)
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
