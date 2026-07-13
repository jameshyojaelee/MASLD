#!/usr/bin/env Rscript
# de_validation_helpers.R
# ===========================================================================
# Shared module for the multi-method DE validation harness (LOO-CV / bootstrap
# / power) comparing dream, DESeq2, and metafor on IDENTICAL data splits.
#
# Sourced (not executed standalone) by:
#   loocv_multimethod.R, bootstrap_multimethod.R, power_multimethod.R
#
# Exposes:
#   - mega_cohorts()                  cohort list from yaml include_in_mega
#   - de_cfg(ds)                      per-cohort de: block from yaml
#   - load_mega_data()                QC-filtered, mega-only counts+meta
#   - run_dream(counts, meta, bp)     -> dt(gene, logFC, stat, padj, method)
#   - run_deseq2(counts, meta, bp)    -> dt(gene, logFC, stat, padj, method)
#   - run_per_study_voom(ds, c, m, contrast_mode) -> dt(gene, logFC, SE_unmoderated, df.total, dataset)
#   - run_metafor(per_study_list, K, bp)          -> dt(gene, logFC, stat, padj, method)  [PROVISIONAL]
#   - storey_pi1(pvals)               Storey pi1 via qvalue (with fallback)
#   - compute_loocv_metrics(method_dt, ho_dt, method_name) -> 1-row metrics dt
#
# COMMON METHOD SCHEMA: every run_*() that represents a "method result" returns
#   data.table(gene, logFC, stat, padj, method)
# with NA-padj rows RETAINED (so each method's tested universe is explicit).
#
# ACTIVE METHODS: controlled by env VALIDATION_METHODS (comma list).
#   Default "dream,deseq2" — metafor is GATED OFF pending the metafor-recipe
#   investigation (the run_metafor recipe below is PROVISIONAL, mirroring the
#   canonical 06_meta_analysis.R: REML + SE_unmoderated, NO HKSJ, min_K=3).
# ===========================================================================

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
  library(limma)
  library(yaml)
})

# --- Force reformulas into lme4 namespace BEFORE variancePartition loads -----
# (verbatim from dream_loo_cv_v2.R lines 36-45)
local({
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
})

suppressPackageStartupMessages({
  library(variancePartition)
  library(BiocParallel)
})

# ===========================================================================
# Paths / constants
# ===========================================================================
PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
PERSTUDY_DIR <- file.path(INT, "results/per_study")
YAML_PATH    <- file.path(PROJECT_ROOT, "config/human_datasets.yaml")

`%||%` <- function(x, y) if (!is.null(x)) x else y

# Active methods (metafor gated off by default — see header)
ACTIVE_METHODS <- function() {
  m <- Sys.getenv("VALIDATION_METHODS", "dream,deseq2")
  trimws(strsplit(m, ",")[[1]])
}

bp_param <- function() {
  ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
  if (is.na(ncpus) || ncpus < 1) ncpus <- 1L
  if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()
}

# ===========================================================================
# Config helpers
# ===========================================================================
mega_cohorts <- function() {
  ycfg <- yaml::read_yaml(YAML_PATH)$datasets
  names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
}

de_cfg <- function(ds) {
  ycfg <- yaml::read_yaml(YAML_PATH)$datasets
  cfg  <- ycfg[[ds]]
  # Unknown dataset (e.g. simulated cohorts simC1..simC5 in the power harness):
  # return defaults. Binary-mode run_per_study_voom only needs norm/qw/sex, so
  # the NULL formula/contrast/levels are never dereferenced for these.
  # (R1 P1-1: previously stop()'d → crashed the metafor power arm.)
  if (is.null(cfg)) {
    return(list(formula = NULL, contrast = NULL, condition_levels = NULL,
                norm_method = "TMM", quality_weights = FALSE, sex_annotated = TRUE))
  }
  list(
    formula          = cfg$de$formula,
    contrast         = cfg$de$contrast,
    condition_levels = cfg$de$condition_levels,
    norm_method      = cfg$de$norm_method %||% "TMM",
    quality_weights  = isTRUE(cfg$de$quality_weights),
    sex_annotated    = isTRUE(cfg$sex_annotated %||% TRUE)
  )
}

# ===========================================================================
# Data loading — raw counts + metadata, QC-filtered, mega cohorts only.
# Filtering/normalization happen INSIDE each method runner (leakage-safe).
# NEVER load merged_dge.rds (pre-filtered) — use merged_counts_raw.rds.
# ===========================================================================
load_mega_data <- function() {
  counts <- readRDS(file.path(RDIR, "merged_counts_raw.rds"))
  meta   <- readRDS(file.path(RDIR, "meta_matched.rds"))
  setDT(meta)
  for (col in names(meta)) if (is.character(meta[[col]])) meta[get(col) == "", (col) := NA]
  qc     <- fread(file.path(INT, "qc/sample_qc_report.csv"))
  pass   <- qc[pass_technical == TRUE, sample_id]
  mega   <- mega_cohorts()

  ds_by_sample <- meta$dataset[match(colnames(counts), meta$sample_id)]
  keep <- (colnames(counts) %in% pass) & (ds_by_sample %in% mega)
  counts <- counts[, keep, drop = FALSE]
  meta   <- meta[match(colnames(counts), sample_id)]
  stopifnot(all(meta$sample_id == colnames(counts)))
  list(counts = counts, meta = meta)
}

# ===========================================================================
# METHOD RUNNER: dream
#   Mirrors dream_loo_cv_v2.R: independent filterByExpr + calcNormFactors(RLE)
#   inside the split, form = ~ group_binary + inferred_sex + (1|dataset),
#   coef "group_binaryDisease", NO eBayes after dream().
#   `counts` = raw integer matrix (genes x samples) for THIS split.
#   `meta`   = data.table aligned to colnames(counts) with group_binary,
#              dataset, inferred_sex.
# ===========================================================================
run_dream <- function(counts, meta, bp = bp_param()) {
  stopifnot(all(meta$sample_id == colnames(counts)))
  dge <- DGEList(counts = counts)
  dge$samples$dataset      <- meta$dataset
  dge$samples$group_binary <- factor(meta$group_binary, levels = c("Control", "Disease"))

  design_filter <- model.matrix(~ 0 + group_binary, data = dge$samples)
  keep <- filterByExpr(dge, design = design_filter)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge, method = "RLE")

  info <- data.frame(
    group_binary = factor(meta$group_binary, levels = c("Control", "Disease")),
    dataset      = droplevels(factor(meta$dataset)),
    inferred_sex = factor(meta$inferred_sex),
    row.names    = colnames(dge),
    stringsAsFactors = FALSE
  )
  if (nlevels(info$dataset) < 2)
    stop("run_dream: need >= 2 datasets for (1|dataset); got ", nlevels(info$dataset))

  # Drop inferred_sex from the model only if degenerate (single level) — keeps
  # bootstrap subsamples from hard-failing; LOO never triggers this.
  form <- if (nlevels(droplevels(info$inferred_sex)) >= 2)
            ~ group_binary + inferred_sex + (1 | dataset)
          else ~ group_binary + (1 | dataset)

  v   <- suppressWarnings(voomWithDreamWeights(dge, form, info, BPPARAM = bp))
  fit <- suppressWarnings(dream(v, form, info, BPPARAM = bp))   # NO eBayes
  tt  <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), logFC = tt$logFC, stat = tt$t,
             padj = tt$adj.P.Val, method = "dream")
}

# ===========================================================================
# METHOD RUNNER: DESeq2
#   Mirrors mega_validation/06_DESeq2_cohortFE.R: prefilter rowSums(>=10)>=10,
#   design ~ dataset + inferred_sex + group_binary, Wald,
#   results(name="group_binary_Disease_vs_Control"). NA-padj rows retained.
# ===========================================================================
run_deseq2 <- function(counts, meta, bp = bp_param()) {
  suppressPackageStartupMessages(library(DESeq2))
  stopifnot(all(meta$sample_id == colnames(counts)))
  cm <- as.matrix(counts); storage.mode(cm) <- "integer"
  keep <- rowSums(cm >= 10) >= 10
  cm   <- cm[keep, , drop = FALSE]

  inferred_sex <- droplevels(factor(meta$inferred_sex))
  dataset      <- droplevels(factor(meta$dataset))
  coldata <- DataFrame(
    dataset      = dataset,
    inferred_sex = inferred_sex,
    group_binary = factor(meta$group_binary, levels = c("Control", "Disease")),
    row.names    = colnames(cm)
  )
  design <- if (nlevels(inferred_sex) >= 2)
              ~ dataset + inferred_sex + group_binary
            else ~ dataset + group_binary
  dds <- DESeqDataSetFromMatrix(countData = cm, colData = coldata, design = design)
  ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
  if (!is.na(ncpus) && ncpus > 1) {
    dds <- DESeq(dds, test = "Wald", parallel = TRUE, BPPARAM = bp)
  } else {
    dds <- DESeq(dds, test = "Wald")
  }
  r <- results(dds, name = "group_binary_Disease_vs_Control")
  data.table(gene = rownames(r), logFC = r$log2FoldChange, stat = r$stat,
             padj = r$padj, method = "deseq2")   # keep NA padj rows
}

# ===========================================================================
# Per-study limma-voom DE (mirror of 02_per_study_de.R).
#   Returns gene, logFC, SE_unmoderated (= stdev.unscaled * sigma), df.total.
#   contrast_mode:
#     "config" — yaml multi-level averaged-substage contrast (canonical input
#                to the production metafor / meta_analysis_results.csv).
#     "binary" — clean Disease-vs-Control (group_binary), matching dream/DESeq2
#                estimand. [Estimand choice PENDING metafor investigation.]
#   `counts_ds` raw counts for this cohort's samples; `meta_ds` aligned rows.
# ===========================================================================
run_per_study_voom <- function(ds, counts_ds, meta_ds, contrast_mode = "config") {
  cfg <- de_cfg(ds)
  m   <- as.data.frame(meta_ds)
  stopifnot(all(m$sample_id == colnames(counts_ds)))

  if (contrast_mode == "config") {
    m$condition <- factor(m$condition, levels = cfg$condition_levels)
    if ("age" %in% names(m) && !all(is.na(m$age))) m$age <- as.numeric(m$age)
    if ("sex" %in% names(m) && !all(is.na(m$sex))) m$sex <- factor(m$sex)
    if (!cfg$sex_annotated && "inferred_sex" %in% names(m))
      m$inferred_sex <- factor(m$inferred_sex)
    form          <- as.formula(cfg$formula)
    contrast_name <- cfg$contrast
  } else if (contrast_mode == "binary") {
    m$group_binary <- factor(m$group_binary, levels = c("Control", "Disease"))
    # cohort-appropriate covariate: annotated sex if present else inferred_sex
    has_sex <- "sex" %in% names(m) && cfg$sex_annotated && !all(is.na(m$sex))
    if (has_sex) { m$sex <- factor(m$sex); covar <- "+ sex" }
    else if ("inferred_sex" %in% names(m) && !all(is.na(m$inferred_sex))) {
      m$inferred_sex <- factor(m$inferred_sex); covar <- "+ inferred_sex"
    } else covar <- ""
    form          <- as.formula(paste("~ 0 + group_binary", covar))
    contrast_name <- "group_binaryDisease - group_binaryControl"
  } else stop("contrast_mode must be 'config' or 'binary'")

  # Drop NA-covariate rows (matches 02_per_study_de.R)
  fv <- setdiff(all.vars(form), "0")
  for (v in fv) if (v %in% names(m)) {
    ok <- !is.na(m[[v]]); m <- m[ok, , drop = FALSE]; counts_ds <- counts_ds[, ok, drop = FALSE]
  }
  # Drop degenerate single-level factors from the formula (bootstrap safety)
  drop_terms <- fv[vapply(fv, function(v)
    is.factor(m[[v]]) && nlevels(droplevels(m[[v]])) < 2, logical(1))]
  if (length(drop_terms)) {
    for (v in drop_terms) form <- update(form, as.formula(paste(". ~ . -", v)))
  }

  dge    <- DGEList(counts = counts_ds)
  design <- model.matrix(form, data = m)
  keep   <- filterByExpr(dge, design = design)
  dge    <- dge[keep, , keep.lib.sizes = FALSE]
  dge    <- calcNormFactors(dge, method = cfg$norm_method)
  v      <- if (cfg$quality_weights) voomWithQualityWeights(dge, design, plot = FALSE)
            else voom(dge, design, plot = FALSE)
  fit    <- lmFit(v, design)
  contr  <- makeContrasts(contrasts = contrast_name, levels = design)
  fit2   <- eBayes(contrasts.fit(fit, contr))
  su     <- if (is.matrix(fit2$stdev.unscaled)) fit2$stdev.unscaled[, 1] else fit2$stdev.unscaled
  tt     <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), logFC = tt$logFC,
             SE_unmoderated = su * fit2$sigma,
             df.total = fit2$df.total, dataset = ds)
}

# ===========================================================================
# METHOD RUNNER: metafor  [FINALIZED 2026-06-03 per the metafor investigation]
#   Per-gene random-effects meta-analysis of per-study logFC + SE_unmoderated.
#   Recipe (investigation verdict):
#     - REML tau^2 estimator (Veroniki 2016 / Langan 2019: preferred over PM
#       under the cohort size-imbalance here).
#     - SE = SE_unmoderated (M1: moderated-vs-unmoderated empirically a no-op
#       but unmoderated is correct by principle; matches 06_meta_analysis.R).
#     - BINARY per-study estimand: callers feed run_per_study_voom(...,"binary")
#       so metafor targets the SAME Disease-vs-Control coefficient as
#       dream/DESeq2 (M3: removes the averaged-substage estimand mismatch).
#     - TIERED inference (M2): the PRIMARY result `padj` is REML WITHOUT HKSJ
#       (matches the 1,848-DEG comparator in meta_analysis_results.csv, gives
#       non-empty sets for rank/overlap metrics, point estimates concord with
#       dream at rho~0.89). `padj_hksj` (REML + Hartung-Knapp test="knha") and
#       `padj_fe` (fixed/equal-effects floor) are emitted as EXTRA columns for
#       the power FDR-control panel — HKSJ is type-I-honest but ~0 power at K=5,
#       which is itself the argument for dream/DESeq2 being canonical.
#     - Coverage: gene must appear in >= min(3, K) cohorts (matches canonical
#       MIN_DATASETS=3; degrades gracefully for small bootstrap subsamples).
#   per_study_list: data.tables (gene, logFC, SE_unmoderated, dataset).
#   The metric/aggregator functions consume gene/logFC/stat/padj; the extra
#   padj_hksj/padj_fe/I2/tau2 columns are ignored unless explicitly used.
# ===========================================================================
run_metafor <- function(per_study_list, K = length(per_study_list), bp = bp_param()) {
  suppressPackageStartupMessages(library(metafor))
  empty <- data.table(gene = character(), logFC = numeric(), stat = numeric(),
                       padj = numeric(), padj_hksj = numeric(), padj_fe = numeric(),
                       I2 = numeric(), tau2 = numeric(), method = character())
  per_study_list <- Filter(function(x) !is.null(x) && nrow(x) > 0, per_study_list)
  if (length(per_study_list) < min(3L, K)) return(empty)

  ps <- rbindlist(per_study_list, fill = TRUE)
  cohorts <- vapply(per_study_list, function(x) as.character(x$dataset[1]), "")
  genes   <- sort(unique(ps$gene))
  B <- S <- matrix(NA_real_, length(genes), length(cohorts),
                   dimnames = list(genes, cohorts))
  for (i in seq_along(cohorts)) {
    sub <- ps[dataset == cohorts[i]]
    idx <- match(sub$gene, genes)
    B[idx, i] <- sub$logFC
    S[idx, i] <- sub$SE_unmoderated
  }
  min_K    <- min(3L, K)
  idx_keep <- which(rowSums(!is.na(B)) >= min_K)

  run1 <- function(i) {
    yi <- B[i, ]; sei <- S[i, ]
    ok <- !is.na(yi) & !is.na(sei) & sei > 0
    if (sum(ok) < min_K) return(NULL)
    tryCatch({
      ctrl <- list(stepadj = 0.5, maxiter = 1000)
      r  <- rma(yi = yi[ok], sei = sei[ok], method = "REML", control = ctrl)
      rh <- tryCatch(rma(yi = yi[ok], sei = sei[ok], method = "REML",
                         test = "knha", control = ctrl), error = function(e) NULL)
      rf <- tryCatch(rma(yi = yi[ok], sei = sei[ok], method = "EE"),
                     error = function(e) NULL)
      data.table(
        gene      = genes[i],
        logFC     = r$b[1, 1],
        stat      = r$zval,
        pval      = r$pval,                                   # REML, no HKSJ (primary)
        pval_hksj = if (!is.null(rh)) rh$pval else NA_real_,  # REML + Hartung-Knapp
        pval_fe   = if (!is.null(rf)) rf$pval else NA_real_,  # fixed/equal-effects floor
        I2        = r$I2,
        tau2      = r$tau2
      )
    }, error = function(e) NULL)
  }
  res <- rbindlist(Filter(Negate(is.null), bplapply(idx_keep, run1, BPPARAM = bp)))
  if (nrow(res) == 0) return(empty)
  res[, padj      := p.adjust(pval,      method = "BH")]
  res[, padj_hksj := p.adjust(pval_hksj, method = "BH")]
  res[, padj_fe   := p.adjust(pval_fe,   method = "BH")]
  res[, .(gene, logFC, stat, padj, padj_hksj, padj_fe, I2, tau2, method = "metafor")]
}

# ===========================================================================
# METHOD RUNNERS: limma family (voom / voom+QW / trend) and edgeR family
# (QLF / QLF-robust / LRT). All share the canonical C2 FIXED-effect design
#   ~ dataset + inferred_sex + group_binary,  coef "group_binaryDisease"
# fit INSIDE the training split (leakage-safe filterByExpr + TMM), matching the
# DESeq2 runner's cohort-fixed-effect estimand (so they are head-to-head with
# dream's (1|dataset) on the identical split). Each returns the COMMON SCHEMA
#   data.table(gene, logFC, stat, padj, method)  with NA-padj rows retained.
#
# Shared scaffold: build a TMM-normalized, filterByExpr-ed DGEList + C2 design
# on the training split. `inferred_sex` is dropped only if degenerate (single
# level) — LOO never triggers this; bootstrap subsamples might.
# ===========================================================================
.c2_dge_design <- function(counts, meta) {
  stopifnot(all(meta$sample_id == colnames(counts)))
  dge <- DGEList(counts = counts)
  dge$samples$dataset      <- droplevels(factor(meta$dataset))
  dge$samples$inferred_sex <- droplevels(factor(meta$inferred_sex))
  dge$samples$group_binary <- factor(meta$group_binary, levels = c("Control", "Disease"))
  if (nlevels(dge$samples$dataset) < 2)
    stop(".c2_dge_design: need >= 2 datasets; got ", nlevels(dge$samples$dataset))

  design <- if (nlevels(dge$samples$inferred_sex) >= 2)
              model.matrix(~ dataset + inferred_sex + group_binary, data = dge$samples)
            else
              model.matrix(~ dataset + group_binary, data = dge$samples)
  stopifnot("group_binaryDisease" %in% colnames(design))

  keep <- filterByExpr(dge, design = design, group = dge$samples$group_binary)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge, method = "TMM")
  list(dge = dge, design = design, coef = "group_binaryDisease")
}

# --- limma-voom (plain voom -> lmFit -> eBayes) -----------------------------
run_limma_voom <- function(counts, meta, bp = bp_param()) {
  cd  <- .c2_dge_design(counts, meta)
  v   <- voom(cd$dge, cd$design, plot = FALSE)
  fit <- eBayes(lmFit(v, cd$design))
  tt  <- topTable(fit, coef = cd$coef, number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), logFC = tt$logFC, stat = tt$t,
             padj = tt$adj.P.Val, method = "limma_voom")
}

# --- limma-voom quality-weighted (CANONICAL C2 engine) ----------------------
run_limma_voom_qw <- function(counts, meta, bp = bp_param()) {
  cd  <- .c2_dge_design(counts, meta)
  v   <- voomWithQualityWeights(cd$dge, cd$design, plot = FALSE)
  fit <- eBayes(lmFit(v, cd$design))
  tt  <- topTable(fit, coef = cd$coef, number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), logFC = tt$logFC, stat = tt$t,
             padj = tt$adj.P.Val, method = "limma_voom_qw")
}

# --- limma-trend (logCPM -> lmFit -> eBayes(trend=TRUE)) ---------------------
run_limma_trend <- function(counts, meta, bp = bp_param()) {
  cd     <- .c2_dge_design(counts, meta)
  logcpm <- edgeR::cpm(cd$dge, log = TRUE, prior.count = 3)
  fit    <- eBayes(lmFit(logcpm, cd$design), trend = TRUE)
  tt     <- topTable(fit, coef = cd$coef, number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), logFC = tt$logFC, stat = tt$t,
             padj = tt$adj.P.Val, method = "limma_trend")
}

# --- edgeR QLF (glmQLFit + glmQLFTest; estimateDisp robust) -----------------
run_edger_qlf <- function(counts, meta, bp = bp_param(), robust = FALSE) {
  cd  <- .c2_dge_design(counts, meta)
  dge <- estimateDisp(cd$dge, cd$design, robust = robust)
  fit <- glmQLFit(dge, cd$design, robust = robust)
  qlf <- glmQLFTest(fit, coef = cd$coef)
  tt  <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  data.table(gene = rownames(tt), logFC = tt$logFC, stat = tt$F,
             padj = tt$FDR, method = if (robust) "edger_qlf_robust" else "edger_qlf")
}

run_edger_qlf_robust <- function(counts, meta, bp = bp_param())
  run_edger_qlf(counts, meta, bp, robust = TRUE)

# --- edgeR LRT (glmFit + glmLRT) --------------------------------------------
run_edger_lrt <- function(counts, meta, bp = bp_param()) {
  cd  <- .c2_dge_design(counts, meta)
  dge <- estimateDisp(cd$dge, cd$design, robust = TRUE)
  fit <- glmFit(dge, cd$design)
  lrt <- glmLRT(fit, coef = cd$coef)
  tt  <- edgeR::topTags(lrt, n = Inf, sort.by = "none")$table
  data.table(gene = rownames(tt), logFC = tt$logFC, stat = tt$LR,
             padj = tt$FDR, method = "edger_lrt")
}

# ===========================================================================
# Storey pi1 (= 1 - pi0). Uses qvalue smoother; falls back to lambda=0.5 closed
# form when qvalue errors or n < 20 (matches dream_loo_cv_v2.R fallback).
# ===========================================================================
storey_pi1 <- function(pvals) {
  pvals <- pvals[is.finite(pvals)]
  if (length(pvals) < 20) return(NA_real_)
  pi0 <- tryCatch({
    suppressWarnings(qvalue::qvalue(pvals)$pi0)
  }, error = function(e) {
    lambda <- 0.5
    min(1, sum(pvals > lambda) / (length(pvals) * (1 - lambda)))
  })
  1 - pi0
}

# ===========================================================================
# LOO-CV metrics: one method's training result vs the held-out cohort's
# per-study DE. Mirrors the 6 metrics in dream_loo_cv_v2.R (Fix B).
#   method_dt: dt(gene, logFC, padj[, stat])  — this method on the training set
#   ho_dt:     dt(gene, ho_logFC, ho_padj, ho_pval) — held-out cohort per-study
#   Merge on each method's tested universe ∩ held-out (own-universe metrics).
# ===========================================================================
compute_loocv_metrics <- function(method_dt, ho_dt, method_name) {
  md <- method_dt[!is.na(gene), .(gene, train_logFC = logFC, train_padj = padj)]
  common <- merge(md, ho_dt, by = "gene")
  n_common <- nrow(common)

  if (n_common < 10) {
    return(data.table(method = method_name, n_genes_common = n_common,
                      n_degs_train_005 = NA_integer_, n_degs_train_01 = NA_integer_,
                      n_degs_ho_005 = NA_integer_, jaccard_01 = NA_real_,
                      lfc_spearman = NA_real_, lfc_spearman_sig = NA_real_,
                      direction_concordance = NA_real_, auc_replication = NA_real_,
                      fisher_or = NA_real_, fisher_pval = NA_real_, pi1_heldout = NA_real_))
  }

  train_sig01 <- common[!is.na(train_padj) & train_padj < 0.1, gene]
  ho_sig01    <- common[!is.na(ho_padj)    & ho_padj    < 0.1, gene]

  # 1. Jaccard (padj<0.1)
  un <- length(union(train_sig01, ho_sig01))
  jaccard_01 <- if (un > 0) length(intersect(train_sig01, ho_sig01)) / un else 0

  # 2. LFC Spearman (all common + training-DEG-restricted)
  lfc_spearman <- cor(common$train_logFC, common$ho_logFC,
                      method = "spearman", use = "pairwise.complete.obs")
  tsig <- common[!is.na(train_padj) & train_padj < 0.1]
  lfc_spearman_sig <- if (nrow(tsig) > 5)
    cor(tsig$train_logFC, tsig$ho_logFC, method = "spearman", use = "pairwise.complete.obs")
    else NA_real_

  # 3. Direction concordance among training DEGs
  dir_conc <- if (nrow(tsig) > 0)
    mean(sign(tsig$train_logFC) == sign(tsig$ho_logFC)) * 100 else NA_real_

  # 4. AUC: training |logFC| predicts held-out padj<0.05 (Wilcoxon U)
  # Guard NA in ho_padj / train_logFC so the denominator counts only valid pairs
  # (R1 P1-4: NA in the predictor silently deflates AUC).
  ok_auc  <- !is.na(common$ho_padj) & is.finite(common$train_logFC)
  ho_bin  <- as.integer(common$ho_padj[ok_auc] < 0.05)
  lfc_auc <- abs(common$train_logFC[ok_auc])
  pos <- lfc_auc[ho_bin == 1]
  neg <- lfc_auc[ho_bin == 0]
  auc_rep <- if (length(pos) > 0 && length(neg) > 0) {
    wt <- suppressWarnings(wilcox.test(pos, neg))
    as.numeric(wt$statistic) / (length(pos) * length(neg))
  } else NA_real_

  # 5. Fisher OR: training DEGs (padj<0.05) enriched in held-out DEGs (padj<0.05)
  ts5 <- !is.na(common$train_padj) & common$train_padj < 0.05
  hs5 <- !is.na(common$ho_padj)    & common$ho_padj    < 0.05
  ct  <- matrix(c(sum(ts5 & hs5), sum(ts5 & !hs5),
                  sum(!ts5 & hs5), sum(!ts5 & !hs5)), nrow = 2)
  ft  <- tryCatch(fisher.test(ct), error = function(e) list(estimate = NA_real_, p.value = NA_real_))

  # 6. pi1 among training DEGs' held-out p-values
  pi1 <- storey_pi1(common[!is.na(train_padj) & train_padj < 0.1, ho_pval])

  data.table(
    method                = method_name,
    n_genes_common        = n_common,
    n_degs_train_005      = length(common[!is.na(train_padj) & train_padj < 0.05, gene]),
    n_degs_train_01       = length(train_sig01),
    n_degs_ho_005         = length(common[!is.na(ho_padj) & ho_padj < 0.05, gene]),
    jaccard_01            = round(jaccard_01, 4),
    lfc_spearman          = round(lfc_spearman, 4),
    lfc_spearman_sig      = round(lfc_spearman_sig, 4),
    direction_concordance = round(dir_conc, 1),
    auc_replication       = round(auc_rep, 4),
    fisher_or             = round(as.numeric(ft$estimate), 3),
    fisher_pval           = ft$p.value,
    pi1_heldout           = round(pi1, 4)
  )
}

# ===========================================================================
# Cohort-stratified 50% subsample WITHOUT replacement.
#   Identical seed scheme to dream_bootstrap_iter.R: caller sets
#   set.seed(42 + ITER*7919) BEFORE calling this (and before any other RNG),
#   so the subsample indices are byte-identical across all three methods.
#   Returns integer column indices into the (already NA-sex-dropped) sample set.
# ===========================================================================
stratified_half_indices <- function(meta) {
  idx <- integer(0)
  for (ds in unique(meta$dataset)) {
    ds_idx <- which(meta$dataset == ds)
    idx <- c(idx, sample(ds_idx, floor(length(ds_idx) / 2)))
  }
  sort(idx)
}

# Datasets in a sample set that retain BOTH Control and Disease (needed for a
# valid contrast). Used by the bootstrap failure guard.
datasets_with_both_groups <- function(meta) {
  tab <- table(meta$dataset, meta$group_binary)
  rn  <- rownames(tab)
  rn[tab[, "Control"] > 0 & tab[, "Disease"] > 0]
}

cat("[de_validation_helpers] loaded. ACTIVE_METHODS =",
    paste(ACTIVE_METHODS(), collapse = ", "), "\n")
