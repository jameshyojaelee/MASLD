#!/usr/bin/env Rscript
# ============================================================================
# 348b_celltype_fstage_dream.R
#
# Per-cell-type dream model across donor F-stage. Replaces the bulk-derived
# signature-scored heatmap (celltype_stage_mean_consensus.csv) with a
# ground-truth per-donor measurement.
#
# For each cell type:
#   pseudobulk_<ct>.tsv.gz (genes x donors, raw counts)
#   metadata: F_stage_augmented (0..4), dataset, sex, n_cells
#
# Model: voomWithDreamWeights + dream(~ F_stage_factor + dataset + sex + (1|sample))
# Contrasts: F1-F0, F2-F1, F3-F2, F4-F3 (per-transition logFC vectors)
#
# Outputs (under disease_signatures/celltype_fstage_dream/):
#   per_celltype_per_transition_logFC_<COHORT_TAG>.csv   (long: ct, transition, gene, logFC, t, p)
#   celltype_stage_perdonor_summary_<COHORT_TAG>.csv     (ct x transition mean |logFC|)
#   celltype_stage_perdonor_z_<COHORT_TAG>.csv           (ct x transition z-scaled mean |logFC|)
#
# COHORT_TAG: "augmented" (default) uses F_stage_augmented for all 260 donors;
#             "documented" restricts to Andrews documented (58 donors, sensitivity).
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(BiocParallel)
})

BASE       <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
COHORT_TAG <- Sys.getenv("COHORT_TAG", "augmented")  # "augmented" or "documented"
N_THREADS  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
MIN_DONORS_PER_STAGE <- as.integer(Sys.getenv("MIN_DONORS_PER_STAGE", "3"))
MIN_CELLS_PER_DONOR  <- as.integer(Sys.getenv("MIN_CELLS_PER_DONOR", "30"))

PB_DIR  <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk")
OUT_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_dream")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

stage_col <- switch(COHORT_TAG,
  augmented  = "F_stage_augmented_clean",
  documented = "F_stage_documented",
  stop("COHORT_TAG must be 'augmented' or 'documented'")
)

cat(sprintf("[cfg] BASE=%s\n[cfg] COHORT_TAG=%s -> stage column %s\n",
            BASE, COHORT_TAG, stage_col))
cat(sprintf("[cfg] N_THREADS=%d  MIN_DONORS_PER_STAGE=%d  MIN_CELLS_PER_DONOR=%d\n",
            N_THREADS, MIN_DONORS_PER_STAGE, MIN_CELLS_PER_DONOR))

register(SerialParam())  # dream sets its own BPPARAM via param arg

# ----------------------------------------------------------------------------
# Load per-(donor, celltype) metadata
# ----------------------------------------------------------------------------
meta_path <- file.path(PB_DIR, "sample_celltype_metadata.tsv")
stopifnot(file.exists(meta_path))
meta <- fread(meta_path)
cat(sprintf("[meta] loaded %d (donor x celltype) rows\n", nrow(meta)))

# Backwards-compat: protocol-contamination flag and clean F_stage column may
# be missing if Phase 1 atlas refresh hasn't completed yet.
if (!"exclude_stage_analysis" %in% names(meta))
  meta[, exclude_stage_analysis := FALSE]
if (!"F_stage_augmented_clean" %in% names(meta))
  meta[, F_stage_augmented_clean := if ("F_stage_augmented" %in% names(meta))
       F_stage_augmented else NA_real_]

# Drop protocol-contaminated donors (GSE136103, Liver_Atlas). Belt-and-
# suspenders: the upstream pseudobulk script already excludes them.
n_before_excl <- nrow(meta)
n_donors_before_excl <- uniqueN(meta$sample)
meta <- meta[is.na(exclude_stage_analysis) | exclude_stage_analysis == FALSE]
cat(sprintf("[filter] excluded %d rows (%d donors) flagged exclude_stage_analysis; %d rows (%d donors) remain\n",
            n_before_excl - nrow(meta),
            n_donors_before_excl - uniqueN(meta$sample),
            nrow(meta), uniqueN(meta$sample)))

# Restrict to donors with the chosen stage column
meta <- meta[!is.na(get(stage_col))]
cat(sprintf("[meta] after %s filter: %d rows, %d donors, %d celltypes\n",
            stage_col, nrow(meta),
            uniqueN(meta$sample), uniqueN(meta$cell_type)))

if (nrow(meta) == 0) stop("no donors carry ", stage_col)

celltypes <- sort(unique(meta$cell_type))

# ----------------------------------------------------------------------------
# Helper: per-cell-type dream
# ----------------------------------------------------------------------------
fit_one <- function(ct) {
  ct_safe <- gsub("[^A-Za-z0-9]+", "_", ct)
  pb_file <- file.path(PB_DIR, sprintf("pseudobulk_%s.tsv.gz", ct_safe))
  if (!file.exists(pb_file)) {
    cat(sprintf("[skip %s] no pseudobulk file\n", ct))
    return(NULL)
  }

  pb <- fread(pb_file, header = TRUE)
  setnames(pb, 1, "gene")
  genes <- pb$gene
  M <- as.matrix(pb[, -1, with = FALSE])
  rownames(M) <- genes
  storage.mode(M) <- "integer"

  ct_meta <- meta[cell_type == ct & sample %in% colnames(M)]
  ct_meta <- ct_meta[n_cells >= MIN_CELLS_PER_DONOR]
  if (nrow(ct_meta) < 6) {
    cat(sprintf("[skip %s] only %d donors\n", ct, nrow(ct_meta)))
    return(NULL)
  }

  ct_meta[, F_factor := factor(get(stage_col),
                               levels = sort(unique(get(stage_col))))]
  if (nlevels(ct_meta$F_factor) < 2) {
    cat(sprintf("[skip %s] only one stage level\n", ct))
    return(NULL)
  }

  # Drop stages with too few donors
  stage_n <- ct_meta[, .N, by = F_factor]
  keep_stages <- stage_n[N >= MIN_DONORS_PER_STAGE, F_factor]
  ct_meta <- ct_meta[F_factor %in% keep_stages]
  ct_meta[, F_factor := droplevels(F_factor)]
  if (nlevels(ct_meta$F_factor) < 2) {
    cat(sprintf("[skip %s] <2 stages after donor-N filter\n", ct))
    return(NULL)
  }

  M_ct <- M[, ct_meta$sample, drop = FALSE]

  cat(sprintf("[%s] donors: %d, stages: %s\n",
              ct, ncol(M_ct), paste(levels(ct_meta$F_factor), collapse = ",")))

  # Filter low-expressed genes (CPM > 1 in >= 3 samples)
  dge <- DGEList(counts = M_ct)
  keep_g <- rowSums(cpm(dge) > 1) >= 3
  dge <- dge[keep_g, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  cat(sprintf("[%s] genes after CPM filter: %d\n", ct, nrow(dge)))

  # Design — fixed effect F_factor + dataset; donor=sample is the unit so no
  # repeated-measures random effect (each pseudobulk is one donor x one CT).
  has_dataset <- uniqueN(ct_meta$dataset) > 1
  has_sex     <- uniqueN(ct_meta$sex)     > 1 && all(!is.na(ct_meta$sex))
  fmla_terms  <- "~ F_factor"
  if (has_dataset) fmla_terms <- paste(fmla_terms, "+ dataset")
  if (has_sex)     fmla_terms <- paste(fmla_terms, "+ sex")
  fmla <- as.formula(fmla_terms)
  cat(sprintf("[%s] formula: %s\n", ct, fmla_terms))

  param <- BiocParallel::SnowParam(N_THREADS, "SOCK", progressbar = FALSE)

  vobj <- tryCatch(
    voomWithDreamWeights(dge, fmla, ct_meta, BPPARAM = param),
    error = function(e) {
      cat(sprintf("[%s] voom failed: %s -- falling back to plain voom\n",
                  ct, conditionMessage(e)))
      voom(dge, model.matrix(fmla, ct_meta))
    }
  )

  fit <- tryCatch(
    dream(vobj, fmla, ct_meta, BPPARAM = param),
    error = function(e) {
      cat(sprintf("[%s] dream failed: %s -- falling back to lmFit\n",
                  ct, conditionMessage(e)))
      lmFit(vobj, model.matrix(fmla, ct_meta))
    }
  )
  fit <- tryCatch(eBayes(fit), error = function(e) fit)

  # ----- Contrasts: adjacent stage transitions present in this CT -----
  stage_levels <- levels(ct_meta$F_factor)
  pairs <- list()
  for (i in seq_len(length(stage_levels) - 1)) {
    a <- stage_levels[i]; b <- stage_levels[i + 1]
    pairs[[paste0("F", b, "_vs_F", a)]] <- list(a = a, b = b)
  }

  coef_names <- colnames(fit$coefficients)
  intercept_name <- coef_names[1]      # (Intercept) reflects F_factor baseline

  results <- list()
  for (nm in names(pairs)) {
    a <- pairs[[nm]]$a; b <- pairs[[nm]]$b
    coef_a <- if (a == stage_levels[1]) NULL else paste0("F_factor", a)
    coef_b <- if (b == stage_levels[1]) NULL else paste0("F_factor", b)

    # Contrast vector: zeros, +1 for coef_b, -1 for coef_a
    cont <- setNames(rep(0, length(coef_names)), coef_names)
    if (!is.null(coef_b) && coef_b %in% coef_names) cont[coef_b] <- cont[coef_b] + 1
    if (!is.null(coef_a) && coef_a %in% coef_names) cont[coef_a] <- cont[coef_a] - 1
    if (all(cont == 0)) next

    # Apply via contrasts.fit (dream returns MArrayLM-like)
    cfit <- tryCatch(
      contrasts.fit(fit, contrasts = matrix(cont, ncol = 1,
                                            dimnames = list(coef_names, nm))),
      error = function(e) NULL
    )
    if (is.null(cfit)) next
    cfit <- tryCatch(eBayes(cfit), error = function(e) cfit)

    tab <- topTable(cfit, coef = 1, number = Inf, sort.by = "none")
    tab$gene       <- rownames(tab)
    tab$cell_type  <- ct
    tab$transition <- nm
    results[[nm]] <- as.data.table(tab)
  }
  if (!length(results)) return(NULL)
  rbindlist(results, fill = TRUE)
}

# ----------------------------------------------------------------------------
# Run all cell types
# ----------------------------------------------------------------------------
all_results <- list()
for (ct in celltypes) {
  cat("\n=========================================================\n")
  res <- fit_one(ct)
  if (!is.null(res)) all_results[[ct]] <- res
}

if (!length(all_results)) stop("no celltype results to save")

big <- rbindlist(all_results, fill = TRUE)

out_logfc <- file.path(OUT_DIR,
  sprintf("per_celltype_per_transition_logFC_%s.csv", COHORT_TAG))
fwrite(big, out_logfc)
cat(sprintf("\n[out] wrote %s (%d rows)\n", out_logfc, nrow(big)))

# Summary: per (cell_type, transition) mean LFC (signed) and mean |LFC|
summ <- big[, .(
    n_genes      = .N,
    mean_lfc     = mean(logFC, na.rm = TRUE),
    mean_abs_lfc = mean(abs(logFC), na.rm = TRUE),
    median_abs_lfc = median(abs(logFC), na.rm = TRUE),
    n_padj_lt_05 = sum(adj.P.Val < 0.05, na.rm = TRUE)
  ), by = .(cell_type, transition)]

# Wide tables for the figure
wide_signed <- dcast(summ, cell_type ~ transition, value.var = "mean_lfc")
wide_abs    <- dcast(summ, cell_type ~ transition, value.var = "mean_abs_lfc")

out_summ <- file.path(OUT_DIR,
  sprintf("celltype_stage_perdonor_summary_%s.csv", COHORT_TAG))
fwrite(summ, out_summ)
cat(sprintf("[out] wrote %s\n", out_summ))

out_signed <- file.path(OUT_DIR,
  sprintf("celltype_stage_perdonor_meanLFC_%s.csv", COHORT_TAG))
fwrite(wide_signed, out_signed)
cat(sprintf("[out] wrote %s\n", out_signed))

out_abs <- file.path(OUT_DIR,
  sprintf("celltype_stage_perdonor_meanabs_%s.csv", COHORT_TAG))
fwrite(wide_abs, out_abs)
cat(sprintf("[out] wrote %s\n", out_abs))

cat("\nDone.\n")
