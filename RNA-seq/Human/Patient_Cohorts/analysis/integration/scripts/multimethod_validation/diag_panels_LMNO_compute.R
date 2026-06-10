#!/usr/bin/env Rscript
# ===========================================================================
# diag_panels_LMNO_compute.R
# ---------------------------------------------------------------------------
# Heavy compute for the cross-method DE statistical-diagnostic panels
# (L = variance partition, M = residual PVCA, N = genomic inflation lambda + QQ,
#  O = p-value distribution + Storey pi1).
#
# Methods (5): limma-voom QW (canonical), dream, DESeq2, metafor (RE), edgeR-QLF.
#   - canonical / dream / metafor p-value vectors are READ from disk.
#   - DESeq2 + edgeR-QLF are FIT FRESH on the 846-sample mega set (not on disk).
#
# Writes CSV intermediates into
#   results/integration/multimethod_validation/diagnostics/
# consumed by scripts/figures/figS_multimethod_diagnostics.R.
#
# Reuses de_validation_helpers.R (load_mega_data, run_deseq2 recipe, storey_pi1,
# bp_param) and mirrors the batch-correction code in
# scripts/figures/figS_multimethod_batch_pca.R for the residual matrices.
# Does NOT edit the shared helper / 05h / the stale variance_partition.csv.
# ===========================================================================

PROJECT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(PROJECT, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
OUTD <- file.path(RDIR, "multimethod_validation/diagnostics")
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)

# de_validation_helpers.R sets up reformulas->lme4, variancePartition, edgeR,
# limma, data.table, BiocParallel, yaml and exposes load_mega_data/run_deseq2/
# storey_pi1/bp_param. Source it for the engines + the env setup.
source(file.path(INT, "scripts/multimethod_validation/de_validation_helpers.R"))
suppressPackageStartupMessages({ library(ggplot2) })

NCPU <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(NCPU) || NCPU < 1) NCPU <- 1L
register(if (NCPU > 1) MulticoreParam(NCPU) else SerialParam())
cat(sprintf("[diag] %d CPU(s); OUT=%s\n", NCPU, OUTD))

base_id <- function(x) sub("[.][0-9]+$", "", x)   # strip Ensembl version

# ===========================================================================
# 1. edgeR quasi-likelihood runner (NOT in the shared helper; canonical-style
#    design ~ dataset + inferred_sex + group_binary, glmQLFit/glmQLFTest).
# ===========================================================================
run_edger_qlf <- function(counts, meta, bp = bp_param()) {
  stopifnot(all(meta$sample_id == colnames(counts)))
  dge <- DGEList(counts = counts)
  info <- data.frame(
    group_binary = factor(meta$group_binary, levels = c("Control", "Disease")),
    dataset      = droplevels(factor(meta$dataset)),
    inferred_sex = droplevels(factor(meta$inferred_sex)),
    row.names    = colnames(counts))
  design <- if (nlevels(info$inferred_sex) >= 2)
              model.matrix(~ dataset + inferred_sex + group_binary, data = info)
            else model.matrix(~ dataset + group_binary, data = info)
  stopifnot("group_binaryDisease" %in% colnames(design))
  keep <- filterByExpr(dge, design = design)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge, method = "RLE")
  dge  <- estimateDisp(dge, design)
  fit  <- glmQLFit(dge, design)
  qlf  <- glmQLFTest(fit, coef = "group_binaryDisease")
  tt   <- topTags(qlf, n = Inf, sort.by = "none")$table
  data.table(gene = rownames(tt), logFC = tt$logFC, stat = tt$F,
             P.Value = tt$PValue, padj = tt$FDR, method = "edger_qlf")
}

# ===========================================================================
# 2. DESeq2 runner with raw P.Value retained.
#    Mirrors de_validation_helpers::run_deseq2 EXACTLY but adds r$pvalue
#    (the shared helper returns stat/padj only). Most-faithful option (b) in
#    the plan; avoids any independent-filtering / Cook's-outlier NA mismatch.
# ===========================================================================
run_deseq2_pv <- function(counts, meta, bp = bp_param()) {
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
    row.names    = colnames(cm))
  design <- if (nlevels(inferred_sex) >= 2)
              ~ dataset + inferred_sex + group_binary
            else ~ dataset + group_binary
  dds <- DESeqDataSetFromMatrix(countData = cm, colData = coldata, design = design)
  dds <- if (NCPU > 1) DESeq(dds, test = "Wald", parallel = TRUE, BPPARAM = bp)
         else          DESeq(dds, test = "Wald")
  r <- results(dds, name = "group_binary_Disease_vs_Control")
  data.table(gene = rownames(r), logFC = r$log2FoldChange, stat = r$stat,
             P.Value = r$pvalue, padj = r$padj, method = "deseq2")
}

# ===========================================================================
# 3. Fresh fits (DESeq2 + edgeR-QLF) on the 846 mega set.
# ===========================================================================
cat("[diag] loading mega data (846)...\n")
md <- load_mega_data()
cat(sprintf("[diag] mega: %d genes x %d samples; cohorts %s\n",
            nrow(md$counts), ncol(md$counts),
            paste(sort(unique(md$meta$dataset)), collapse = ", ")))

cat("[diag] fitting DESeq2 (Wald, ~dataset+inferred_sex+group_binary)...\n")
deseq2 <- run_deseq2_pv(md$counts, md$meta)
fwrite(deseq2, file.path(OUTD, "deseq2_mega_results.csv"))
# sanity: option-(a) derivation vs DESeq2's own pvalue (verification log)
.dz <- deseq2[is.finite(stat) & is.finite(P.Value)]
.corr_pa <- suppressWarnings(cor(2 * pnorm(-abs(.dz$stat)), .dz$P.Value, use = "complete.obs"))
cat(sprintf("[diag] DESeq2 done: %d genes; corr(2*pnorm(-|z|), pvalue) = %.5f\n",
            nrow(deseq2), .corr_pa))

cat("[diag] fitting edgeR quasi-likelihood...\n")
edger <- run_edger_qlf(md$counts, md$meta)
fwrite(edger, file.path(OUTD, "edger_qlf_mega_results.csv"))
cat(sprintf("[diag] edgeR-QLF done: %d genes\n", nrow(edger)))

# ===========================================================================
# 4. Assemble per-method p-value vectors (keyed by Ensembl BASE id, since
#    metafor uses unversioned IDs while the rest are versioned).
# ===========================================================================
canon <- fread(file.path(RDIR, "canonical_deg_results.csv"))[, .(gene, P.Value)]
dream <- fread(file.path(RDIR, "dream_results.csv"))[, .(gene, P.Value)]
mfor  <- fread(file.path(RDIR, "meta_results_ashr.csv"))[, .(gene, P.Value)]

pv_list <- list(
  `limma-voom QW (canonical)` = canon,
  dream                       = dream,
  DESeq2                      = deseq2[, .(gene, P.Value)],
  `metafor (RE)`              = mfor,
  `edgeR-QLF`                 = edger[, .(gene, P.Value)])
for (m in names(pv_list)) {
  d <- pv_list[[m]]
  d <- d[is.finite(P.Value)]
  d[, eb := base_id(gene)]
  d <- d[!duplicated(eb)]                         # one row per base id
  pv_list[[m]] <- d
}
cat("[diag] per-method tested genes:\n")
for (m in names(pv_list)) cat(sprintf("    %-26s %d\n", m, nrow(pv_list[[m]])))

# common universe across all 5 methods
common_eb <- Reduce(intersect, lapply(pv_list, function(d) d$eb))
cat(sprintf("[diag] common universe (intersection of 5): %d genes\n", length(common_eb)))

# ===========================================================================
# 5. lambda_GC + Storey pi1 (both on common set and each method's full set).
#    lambda computed via the p-value -> 1df chi-square route so the 5 methods'
#    different null reference distributions (moderated-t / Wald-z / RE-z / QL-F)
#    are directly comparable.
# ===========================================================================
lambda_gc <- function(p) {
  n_floor <- sum(p == 0, na.rm = TRUE)
  p[p == 0] <- .Machine$double.xmin
  p <- p[is.finite(p) & p > 0 & p <= 1]
  list(lambda = median(qchisq(1 - p, df = 1)) / qchisq(0.5, df = 1),
       n = length(p), n_floor = n_floor)
}

lambda_rows <- list(); pi1_rows <- list()
for (m in names(pv_list)) {
  d <- pv_list[[m]]
  full <- lambda_gc(d$P.Value)
  comm <- lambda_gc(d[eb %in% common_eb, P.Value])
  lambda_rows[[m]] <- data.table(
    method = m, lambda_common = comm$lambda, lambda_full = full$lambda,
    n_common = comm$n, n_full = full$n, n_p_floored = full$n_floor)
  pi1_rows[[m]] <- data.table(
    method = m, pi1_full = storey_pi1(d$P.Value),
    pi1_common = storey_pi1(d[eb %in% common_eb, P.Value]), n_full = full$n)
}
lambda_table <- rbindlist(lambda_rows)
pi1_table    <- rbindlist(pi1_rows)
fwrite(lambda_table, file.path(OUTD, "lambda_table.csv"))
fwrite(pi1_table,    file.path(OUTD, "pi1_table.csv"))
cat("[diag] lambda_table:\n"); print(lambda_table)
cat("[diag] pi1_table:\n");    print(pi1_table)

# ===========================================================================
# 6. QQ data (common universe), per method. Beta order-statistic CI band is a
#    function of (rank i, m) only, so it is shared across methods. Thin the
#    null bulk to keep the PDF small (keep top-2000 by significance + every 20th).
# ===========================================================================
m_common <- length(common_eb)
qq_list <- lapply(names(pv_list), function(m) {
  p <- sort(pv_list[[m]][eb %in% common_eb, P.Value])
  p[p == 0] <- .Machine$double.xmin
  i <- seq_along(p)
  keepr <- i <= 2000 | (i %% 20 == 0) | i >= (m_common - 50)
  data.table(
    method = m, idx = i[keepr],
    obs = -log10(p[keepr]),
    exp = -log10((i[keepr] - 0.5) / m_common),
    lo  = -log10(qbeta(0.975, i[keepr], m_common - i[keepr] + 1)),
    hi  = -log10(qbeta(0.025, i[keepr], m_common - i[keepr] + 1)))
})
qq_data <- rbindlist(qq_list)
fwrite(qq_data, file.path(OUTD, "qq_data.csv"))
cat(sprintf("[diag] qq_data: %d rows (m_common=%d)\n", nrow(qq_data), m_common))

# Full per-method p-value vectors for the Panel O histograms (own universe).
fwrite(rbindlist(lapply(names(pv_list), function(m)
         data.table(method = m, P.Value = pv_list[[m]]$P.Value))),
       file.path(OUTD, "pvalue_long.csv"))

# ===========================================================================
# 7. Raw variance partition on the 846 mega subset (Panel L).
#    Mirror 04_variance_partition.R but subset merged_dge.rds to mega (like 05h).
#    Writes variance_partition_mega846.csv — does NOT touch variance_partition.csv.
# ===========================================================================
cat("[diag] raw variance partition on 846 mega subset...\n")
dge_all <- readRDS(file.path(RDIR, "merged_dge.rds"))
ycfg <- yaml::read_yaml(file.path(PROJECT, "config/human_datasets.yaml"))$datasets
mega <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
dge_mega <- dge_all[, dge_all$samples$dataset %in% mega]
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
sex_v <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]
if (any(is.na(sex_v))) { cat(sprintf("  WARN %d NA inferred_sex -> 'Unknown'\n", sum(is.na(sex_v)))); sex_v[is.na(sex_v)] <- "Unknown" }
info <- data.frame(
  group_binary = factor(dge_mega$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = factor(dge_mega$samples$dataset),
  inferred_sex = factor(sex_v),
  row.names    = colnames(dge_mega))
cat(sprintf("  %d genes x %d samples; %d cohorts; sex %s\n",
            nrow(dge_mega), ncol(dge_mega), nlevels(info$dataset),
            paste(levels(info$inferred_sex), collapse = "/")))
v_vp <- voom(dge_mega, model.matrix(~ group_binary, data = info))
form_vp <- ~ (1 | group_binary) + (1 | inferred_sex) + (1 | dataset)
vp <- fitExtractVarPartModel(v_vp, form_vp, info)
vp_dt <- as.data.table(as.data.frame(vp), keep.rownames = "gene")
fwrite(vp_dt, file.path(OUTD, "variance_partition_mega846.csv"))
cat("[diag] varpart median fractions:\n"); print(round(apply(vp, 2, median), 5))

# ===========================================================================
# 8. Residual PVCA, per method's batch model (Panel M).
#    Reuse the figS_multimethod_batch_pca.R correction recipe: top-2000
#    within-cohort HVG log-CPM -> X (Raw); X_fixed = removeBatchEffect (~grp,
#    fixed dataset = limma-voom-QW / DESeq2 / edgeR-QLF); X_rand = per-gene
#    lmer (1|dataset) shrunken-BLUP subtraction (dream). metafor = n/a (no
#    pooled matrix) -> handled in the figure.
#    PVCA = eigenvalue-weighted average of per-PC variance components from a
#    random-effects model lmer(PC ~ (1|dataset)+(1|group_binary)+(1|inferred_sex)),
#    which is order-independent (Bushel-style).
# ===========================================================================
cat("[diag] residual PVCA (per-method batch correction)...\n")
logcpm <- edgeR::cpm(dge_mega, log = TRUE, prior.count = 1)
logcpm_wc <- logcpm
for (d in unique(info$dataset)) {
  idx <- which(info$dataset == d)
  logcpm_wc[, idx] <- logcpm[, idx] - rowMeans(logcpm[, idx, drop = FALSE])
}
rv  <- matrixStats::rowVars(logcpm_wc)
top <- order(rv, decreasing = TRUE)[seq_len(min(2000, length(rv)))]
X   <- logcpm[top, ]
grp <- info$group_binary
ds  <- info$dataset

X_fixed <- limma::removeBatchEffect(X, batch = ds, design = model.matrix(~ grp))

ctrl_lmer <- lme4::lmerControl(optimizer = "nloptwrap", check.conv.singular = "ignore")
ds_chr <- as.character(ds); X_rand <- X
for (i in seq_len(nrow(X))) {
  y <- X[i, ]
  b <- tryCatch({
    fit <- suppressWarnings(suppressMessages(lme4::lmer(y ~ grp + (1 | ds), control = ctrl_lmer)))
    re <- lme4::ranef(fit)$ds[, 1]; names(re) <- rownames(lme4::ranef(fit)$ds); re[ds_chr]
  }, error = function(e) { mns <- tapply(y, ds_chr, mean); (mns - mean(mns))[ds_chr] })
  X_rand[i, ] <- y - as.numeric(b)
}

# PVCA on a corrected matrix Xc (genes x samples)
covars <- data.frame(dataset = ds, group_binary = grp, inferred_sex = info$inferred_sex)
pvca_one <- function(Xc, cum_cut = 0.90) {
  pr  <- prcomp(t(Xc), center = TRUE, scale. = FALSE)
  eig <- pr$sdev^2
  k   <- max(1L, which(cumsum(eig) / sum(eig) >= cum_cut)[1])
  w   <- eig[seq_len(k)] / sum(eig[seq_len(k)])          # normalise over kept PCs
  comps <- c("dataset", "group_binary", "inferred_sex", "Residuals")
  acc <- setNames(numeric(length(comps)), comps)
  for (j in seq_len(k)) {
    df <- cbind(score = pr$x[, j], covars)
    fit <- suppressWarnings(lme4::lmer(
      score ~ (1 | dataset) + (1 | group_binary) + (1 | inferred_sex),
      data = df, control = ctrl_lmer))
    vc  <- as.data.frame(lme4::VarCorr(fit))
    vv  <- setNames(vc$vcov, ifelse(vc$grp == "Residual", "Residuals", vc$grp))
    fr  <- vv / sum(vv)
    for (cc in comps) acc[cc] <- acc[cc] + w[j] * (if (cc %in% names(fr)) fr[cc] else 0)
  }
  list(frac = acc, k = k)
}
states <- list(Raw = X, `Fixed (limma-QW/DESeq2/edgeR)` = X_fixed, `Random (dream)` = X_rand)
pvca_dt <- rbindlist(lapply(names(states), function(s) {
  r <- pvca_one(states[[s]])
  cat(sprintf("  %-30s PCs=%d  dataset=%.3f group=%.3f sex=%.3f resid=%.3f\n",
              s, r$k, r$frac["dataset"], r$frac["group_binary"],
              r$frac["inferred_sex"], r$frac["Residuals"]))
  data.table(correction = s, covariate = names(r$frac), pct_variance = 100 * as.numeric(r$frac))
}))
fwrite(pvca_dt, file.path(OUTD, "pvca_before_after.csv"))

cat("\n[diag] DONE. Wrote intermediates to", OUTD, "\n")
cat("  deseq2_mega_results.csv, edger_qlf_mega_results.csv, lambda_table.csv,\n")
cat("  pi1_table.csv, qq_data.csv, pvalue_long.csv, variance_partition_mega846.csv,\n")
cat("  pvca_before_after.csv\n")
