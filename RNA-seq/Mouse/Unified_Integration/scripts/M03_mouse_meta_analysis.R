#!/usr/bin/env Rscript
# M03_mouse_meta_analysis.R
# ---------------------------------------------------------------------------
# Two analyses:
#  (A) Per-diet-type meta-analysis — random-effects rma() across diet types
#  (B) Pooled mega-analysis — limma-voom quality-weighted across ALL disease vs
#      control samples, dataset as FIXED effect (mirrors the human
#      limma_voom_qw__C2 canonical; dream retired for mouse 2026-06-16), ashr-shrunk
# Input:  per_diet/*_de_results.csv, merged_counts_raw.rds, meta_matched.rds
# Output: meta_per_diet.csv, lvqw_pooled_results.csv
# ---------------------------------------------------------------------------

# ---- Seed pinning (T2.4, 2026-04-22) -----
set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(metafor)
  library(ashr)
  library(BiocParallel)
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
METADIR <- file.path(RDIR, "meta_analysis")
dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)

cat("=== M03: Meta-Analysis + Pooled limma-voom-qw ===\n\n")

# ============================================================
# (A) Per-Diet-Type Meta-Analysis (random-effects rma)
# ============================================================
cat("===== (A) PER-DIET META-ANALYSIS =====\n\n")

# BG-012 fix (2026-08-07): consume the explicit diet_models manifest from
# config/mouse_datasets.yaml instead of globbing per_diet/. Globbing silently
# ingested retired models (e.g. LIDPAD, AMLN_ob) or a stale/overwritten file,
# changing the canonical mouse meta with no schema or provenance failure.
# Fail closed on BOTH missing and extra models.
cfg_path <- file.path(PROJECT_ROOT, "config/mouse_datasets.yaml")
if (!file.exists(cfg_path)) stop("Mouse config not found: ", cfg_path)
mouse_cfg <- yaml::read_yaml(cfg_path)
diet_models <- as.character(mouse_cfg$diet_models)
if (!length(diet_models) || anyNA(diet_models) || anyDuplicated(diet_models)) {
  stop("config/mouse_datasets.yaml diet_models must be a nonempty unique list")
}
cat("Expected", length(diet_models), "diet models from config:",
    paste(diet_models, collapse = ", "), "\n")

expected_files <- file.path(DEDIR, paste0(diet_models, "_de_results.csv"))
missing_models <- diet_models[!file.exists(expected_files)]
if (length(missing_models)) {
  stop("Missing per-diet DE results for expected diet model(s): ",
       paste(missing_models, collapse = ", "),
       "\n  Expected in ", DEDIR,
       "\n  Run M02_mouse_per_diet_de.R first.")
}
observed_models <- sub("_de_results\\.csv$", "",
                       basename(list.files(DEDIR, pattern = "_de_results\\.csv$")))
extra_models <- setdiff(observed_models, diet_models)
if (length(extra_models)) {
  stop("Refusing unexpected per-diet DE results in ", DEDIR, ": ",
       paste(sort(extra_models), collapse = ", "),
       "\n  These are not in config/mouse_datasets.yaml diet_models.",
       "\n  Archive retired models before rerunning M03.")
}

de_list <- list()
for (dm in diet_models) {
  f <- file.path(DEDIR, paste0(dm, "_de_results.csv"))
  de_list[[dm]] <- fread(f)
  cat("  Loaded:", dm, "-", nrow(de_list[[dm]]), "genes\n")
}

# Find genes present in >= 3 diet models.
# Rationale: Balances gene coverage with measurement reliability
# Aligns with human pipeline (>= 2 datasets) and Tier 3 logic (>= 3 diets)
# NOTE (2026-08-07): the comment here previously read "out of 5 total". The
# registry has held 4 models since LIDPAD was removed (2026-05-28), so the
# retained threshold of 3 is now >=3/4, not >=3/5 — a stringency change that
# happened silently. The threshold is left at 3 deliberately; revisit it as a
# scientific decision, not as part of this robustness fix.
all_unique_genes <- unique(unlist(lapply(de_list, function(x) x$gene)))
gene_diet_counts <- sapply(all_unique_genes, function(g) {
  sum(sapply(de_list, function(x) g %in% x$gene))
})

min_diets_required <- 3
all_genes <- all_unique_genes[gene_diet_counts >= min_diets_required]

cat("\nGene coverage across diet models:\n")
cat("  Total unique genes:", length(all_unique_genes), "\n")
cat("  Genes in >=", min_diets_required, "models:", length(all_genes),
    sprintf("(%.1f%%)\n", 100 * length(all_genes) / length(all_unique_genes)))
cat("  Distribution:\n")
print(table(gene_diet_counts))
cat("\n")

# ---- BG-011 (2026-08-09): shared-control sampling covariance --------------------
# The registry declares GSE162876 as diet_model "CDAHFD_FPC" with a shared control group,
# and M02's selection rule (`diet_model == dm | group_binary == "Control"`) therefore puts
# the SAME control samples in both comparisons. Passing those two effects to rma(yi, sei)
# asserts independent sampling errors, which they are not: they share a control mean, so
# their errors are positively correlated and the pooled SE / I-squared come out
# overconfident.
#
# Gleser & Olkin (2009), Handbook of Research Synthesis & Meta-Analysis, ch.19, give the
# shared-control covariance for mean-difference contrasts:
#     Cov(d_i, d_j) = sigma^2 * n_shared / (n_Ci * n_Cj)
#     Var(d_i)      = sigma^2 * (1/n_i + 1/n_Ci)
# =>  r_ij = [n_s/(n_Ci n_Cj)] / sqrt[(1/n_i + 1/n_Ci)(1/n_j + 1/n_Cj)]
# This assumes a common sigma^2 across the two fits. CDAHFD and FPC are fit on the same
# cohort through the same pipeline, so the approximation is close; it is stated, not hidden.
#
# The structure is DERIVED from the registry + metadata, never hardcoded, so adding or
# removing a shared-control model updates it automatically. When no diet pair shares
# controls the matrix is the identity and we fall through to the original rma() path, so
# independent designs are numerically untouched.
meta_matched_path <- file.path(RDIR, "meta_matched.rds")
qc_path <- file.path(INT, "qc/sample_qc_report.csv")
shared_R <- diag(length(diet_models))
dimnames(shared_R) <- list(diet_models, diet_models)
if (file.exists(meta_matched_path) && file.exists(qc_path)) {
  .m  <- as.data.table(readRDS(meta_matched_path))
  .qc <- fread(qc_path)
  .m  <- .m[sample_id %in% .qc[pass_qc == TRUE, sample_id]]
  .comps <- list()
  for (.id in names(mouse_cfg$datasets)) {
    .dm <- mouse_cfg$datasets[[.id]]$diet_model
    if (grepl("_", .dm) && !.dm %in% c("AMLN_ob")) .dm <- unlist(strsplit(.dm, "_"))
    for (.d in .dm) .comps[[.d]] <- c(.comps[[.d]], .id)
  }
  .sel <- lapply(diet_models, function(dm) {
    if (is.null(.comps[[dm]])) return(character(0))
    .m[dataset %in% .comps[[dm]] & (diet_model == dm | group_binary == "Control"), sample_id]
  })
  names(.sel) <- diet_models
  .ctl <- lapply(diet_models, function(dm) .m[sample_id %in% .sel[[dm]] & group_binary == "Control", sample_id])
  .ndis <- sapply(diet_models, function(dm) .m[sample_id %in% .sel[[dm]] & group_binary == "Disease", .N])
  names(.ctl) <- diet_models
  for (i in seq_along(diet_models)) for (j in seq_along(diet_models)) if (i < j) {
    a <- diet_models[i]; b <- diet_models[j]
    ns <- length(intersect(.ctl[[a]], .ctl[[b]]))
    if (ns > 0 && .ndis[a] > 0 && .ndis[b] > 0) {
      r <- (ns / (length(.ctl[[a]]) * length(.ctl[[b]]))) /
        sqrt((1 / .ndis[a] + 1 / length(.ctl[[a]])) * (1 / .ndis[b] + 1 / length(.ctl[[b]])))
      r <- max(-0.99, min(0.99, r))
      shared_R[a, b] <- shared_R[b, a] <- r
      cat(sprintf("  [BG-011] %s/%s share %d control samples -> sampling correlation r = %.4f\n",
                  a, b, ns, r))
    }
  }
} else {
  warning("BG-011: meta_matched.rds or sample_qc_report.csv absent; cannot derive the ",
          "shared-control structure. Falling back to the independence assumption.")
}
USE_SHARED_COV <- any(shared_R[upper.tri(shared_R)] != 0)
if (USE_SHARED_COV) {
  if (min(eigen(shared_R, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("BG-011: derived sampling-correlation matrix is not positive definite")
  }
  cat("  [BG-011] using rma.mv with V = D R D (correlated sampling errors)\n")
} else {
  cat("  [BG-011] no diet pair shares controls; using independent rma()\n")
}

# Run rma() for each gene
cat("Running random-effects meta-analysis across diet types...\n")

meta_results <- rbindlist(lapply(all_genes, function(g) {
  # Use scalar extraction; return NA if gene absent in a diet model
  effects <- sapply(de_list, function(x) {
    v <- x[gene == g, logFC]
    if (length(v) == 0L) NA_real_ else v[1L]
  })
  ses <- sapply(de_list, function(x) {
    v <- x[gene == g, SE_unmoderated]  # unmoderated SE from M02 (stdev.unscaled * sigma)
    if (length(v) == 0L) NA_real_ else v[1L]
  })

  # Remove NAs and non-finite SEs (e.g. t=0 gives Inf)
  valid <- !is.na(effects) & !is.na(ses) & is.finite(ses) & ses > 0
  if (sum(valid, na.rm = TRUE) < 2) return(NULL)

  tryCatch({
    y  <- effects[valid]
    s  <- ses[valid]
    dm <- names(effects)[valid]
    # BG-011: when this gene's surviving diets include a shared-control pair, meta-analyse
    # with the correlated sampling-error matrix V = D R D; otherwise the original rma().
    Rg <- shared_R[dm, dm, drop = FALSE]
    if (USE_SHARED_COV && any(Rg[upper.tri(Rg)] != 0)) {
      V   <- diag(s, nrow = length(s)) %*% Rg %*% diag(s, nrow = length(s))
      fit <- rma.mv(yi = y, V = V, random = ~ 1 | factor(dm), method = "REML")
      # rma.mv reports no I^2, so keep the column populated with the Higgins & Thompson
      # (2002) generalization: I^2 = sigma^2 / (sigma^2 + s^2_typical), where the typical
      # within-study variance uses w = 1/v_i. Returning NA here would silently drop a
      # column downstream consumers already read.
      w  <- 1 / diag(V)
      s2_typ <- (length(y) - 1) * sum(w) / (sum(w)^2 - sum(w^2))
      sig2   <- sum(fit$sigma2)
      i2 <- if (is.finite(s2_typ) && (sig2 + s2_typ) > 0) 100 * sig2 / (sig2 + s2_typ) else NA_real_
    } else {
      fit <- rma(yi = y, sei = s, method = "REML")
      i2  <- as.numeric(fit$I2)
    }
    data.table(
      gene = g,
      meta_logFC = as.numeric(fit$b),
      meta_se    = as.numeric(fit$se),
      meta_pval  = as.numeric(fit$pval),
      meta_I2    = i2,
      n_diets    = sum(valid),
      diet_types = paste(dm, collapse = ";")
    )
  }, error = function(e) NULL)
}))

meta_results[, meta_padj := p.adjust(meta_pval, method = "BH")]
meta_results <- meta_results[order(meta_padj)]

sig_meta <- meta_results[meta_padj < 0.05]
cat("\nMeta-analysis DEGs (padj<0.05):", nrow(sig_meta),
    " (Up:", sum(sig_meta$meta_logFC > 0),
    "Down:", sum(sig_meta$meta_logFC < 0), ")\n")

fwrite(meta_results, file.path(METADIR, "meta_per_diet.csv"))
cat("Saved: meta_per_diet.csv\n\n")

# ============================================================
# (B) Pooled limma-voom-qw Mega-Analysis  (replaces dream, 2026-06-16)
# ============================================================
cat("===== (B) POOLED limma-voom-qw MEGA-ANALYSIS =====\n\n")

# Load raw counts and metadata
merged <- readRDS(file.path(RDIR, "merged_counts_raw.rds"))
meta   <- readRDS(file.path(RDIR, "meta_matched.rds"))
qc     <- fread(file.path(INT, "qc/sample_qc_report.csv"))

# --- Detect a tximport gene-length matrix (Kallisto pipeline) ---
# See the BG-013 note at the fit below: length correction cannot be applied to a
# voom object post hoc, so we only detect the matrix here in order to refuse.
lengths_file <- file.path(RDIR, "merged_gene_lengths.rds")
HAS_TX_OFFSETS <- file.exists(lengths_file)
if (HAS_TX_OFFSETS) {
  cat("Gene-length matrix present:", lengths_file, "\n")
  cat("  -> M03 will STOP: limma-voom cannot apply length offsets post hoc (BG-013)\n\n")
} else {
  cat("No tximport gene-length matrix found — none needed (featureCounts)\n\n")
}

# Filter to QC-passing
pass <- qc[pass_qc == TRUE, sample_id]
merged <- merged[, colnames(merged) %in% pass]
meta   <- meta[sample_id %in% pass]

cat("Samples for pooled mega-analysis:", ncol(merged), "\n")

# Disease vs Control distribution
cat("\nGroup distribution:\n")
print(table(meta$group_binary, meta$dataset))

# Set up for dream
meta$group_binary <- factor(meta$group_binary, levels = c("Control", "Disease"))
meta$dataset <- factor(meta$dataset)

# Create DGE and filter
dge <- DGEList(counts = merged)
keep <- filterByExpr(dge, group = meta$group_binary)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

cat("\nGenes after filterByExpr:", nrow(dge), "\n")

# C2-style fixed-effect design (2026-06-16): mouse now mirrors the human
# limma_voom_qw__C2 canonical. `dataset` is a FIXED effect (batch correction);
# the retired dream random intercept (1|dataset) is no longer used. Mouse sex is
# too sparsely annotated to include as a covariate, so it is omitted.
# Align metadata rows to the count-matrix column order before building the design.
meta <- meta[match(colnames(dge), meta$sample_id)]
design <- model.matrix(~ dataset + group_binary, data = meta)
cat("Design columns:", paste(colnames(design), collapse = ", "), "\n")

# voom with quality weights (limma); replaces voomWithDreamWeights
cat("Running voomWithQualityWeights...\n")
vobj <- voomWithQualityWeights(dge, design)

# BG-013 (2026-08-09): this block previously assigned `vobj$offset <- log(lengths)`
# after voomWithQualityWeights() and reported the arm as length-corrected. It was a
# silent no-op -- lmFit() consumes only getEAWP() output, verified on the installed
# limma to be exactly (exprs, Amean, weights, design), with no `offset` member;
# coefficients are bit-identical with and without an arbitrary log-length offset.
# Length correction belongs at the tximport step, not here.
if (HAS_TX_OFFSETS) {
  stop("Refusing to run a mislabelled length-corrected arm.\n",
       "  A gene-length matrix is present (", lengths_file, "), but limma-voom\n",
       "  cannot consume per-observation length offsets: lmFit() reads only\n",
       "  getEAWP() output (exprs, Amean, weights, design) and ignores EList$offset.\n",
       "  Re-run tximport with countsFromAbundance = \"lengthScaledTPM\" so the\n",
       "  effective lengths are carried in the counts, then rerun M03 without a\n",
       "  separate length matrix.", call. = FALSE)
}

cat("Running lmFit + eBayes (limma-voom quality-weighted)...\n")
fit <- lmFit(vobj, design)
fit <- eBayes(fit)

# Extract the Disease-vs-Control coefficient
lvqw_res <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
lvqw_res$gene <- rownames(lvqw_res)
lvqw_res <- as.data.table(lvqw_res)
setcolorder(lvqw_res, "gene")

# ashr adaptive shrinkage (mirrors the human canonical: adds shrunk_logFC + lfsr)
se_hat  <- lvqw_res$logFC / lvqw_res$t
ash_fit <- ashr::ash(as.numeric(lvqw_res$logFC), as.numeric(se_hat),
                     mixcompdist = "normal")
lvqw_res[, shrunk_logFC := ashr::get_pm(ash_fit)]
lvqw_res[, lfsr := ashr::get_lfsr(ash_fit)]

sig_lvqw <- lvqw_res[adj.P.Val < 0.05]
cat("\nPooled limma-voom-qw DEGs (padj<0.05):", nrow(sig_lvqw),
    " (Up:", sum(sig_lvqw$logFC > 0),
    "Down:", sum(sig_lvqw$logFC < 0), ")\n")

fwrite(lvqw_res[order(adj.P.Val)], file.path(METADIR, "lvqw_pooled_results.csv"))
cat("Saved: lvqw_pooled_results.csv\n")

# Save DGE object for downstream
saveRDS(dge, file.path(RDIR, "merged_dge.rds"))
cat("Saved: merged_dge.rds\n")

cat("\nM03 complete.\n")
