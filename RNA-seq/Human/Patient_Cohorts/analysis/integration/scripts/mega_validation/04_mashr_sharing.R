#!/usr/bin/env Rscript
# 04_mashr_sharing.R
# ---------------------------------------------------------------------------
# Mega-validation Arm 4: multivariate adaptive shrinkage (mashr; Urbut 2019)
# on the per-cohort limma-voom posterior estimates. Treats K cohorts as K
# "conditions" (a la GTEx 44 tissues), learns canonical + data-driven
# covariance structure, and classifies each gene into a sharing pattern.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mashr); library(data.table); library(yaml)
})

set.seed(42)   # ExtremeDeconvolution in cov_ed uses random init

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
PERSTUDY <- file.path(INT, "results/per_study")
RDIR_OUT <- file.path(INT, "results/mega_validation/mashr")
dir.create(RDIR_OUT, recursive = TRUE, showWarnings = FALSE)

ycfg <- read_yaml(file.path(PROJECT_ROOT, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
K <- length(mega_cohorts)
cat("Mega cohorts (K =", K, "):", paste(mega_cohorts, collapse=", "), "\n")

per_study <- rbindlist(lapply(mega_cohorts, function(ds) {
  p <- file.path(PERSTUDY, paste0(ds, "_de_results.csv"))
  if (!file.exists(p)) stop("Missing per-study DE for ", ds)
  dt <- fread(p)
  stopifnot("SE" %in% names(dt) || "SE_unmoderated" %in% names(dt))
  dt
}), fill = TRUE)

genes <- sort(unique(per_study$gene))
B <- matrix(NA_real_, nrow = length(genes), ncol = K,
            dimnames = list(genes, mega_cohorts))
S <- B
# P1 fix 2026-05-28: use unmoderated SE for mashr second-level model (consistent
# with Script 06 metafor). The legacy `SE` column is the eBayes-MODERATED SE
# (Script 02 sets SE := SE_moderated); for a second-level hierarchical model
# (mashr) the UNMODERATED SE is the consistent choice. Fall back to `SE` only if
# the per-study file predates the SE_unmoderated column.
se_col <- if ("SE_unmoderated" %in% names(per_study)) "SE_unmoderated" else "SE"
cat("mashr Shat source column:", se_col, "\n")
for (i in seq_len(K)) {
  sub <- per_study[dataset == mega_cohorts[i]]
  idx <- match(sub$gene, genes)
  B[idx, i] <- sub$logFC
  S[idx, i] <- sub[[se_col]]
}

# mashr requires complete cases — drop genes with any NA / non-positive SE.
complete <- rowSums(is.na(B) | is.na(S) | S <= 0) == 0
B <- B[complete, ]; S <- S[complete, ]
cat("Genes with complete per-cohort data:", nrow(B), "/", length(genes), "\n")

# --- mash setup ---
data_all <- mash_set_data(Bhat = B, Shat = S)

# Strong subset: top 1000 genes by max |Z|
Z <- B / S
strong_idx <- order(apply(abs(Z), 1, max), decreasing = TRUE)[seq_len(min(1000, nrow(B)))]
data_strong <- mash_set_data(Bhat = B[strong_idx, ], Shat = S[strong_idx, ])

# --- Covariance estimation ---
U_c   <- cov_canonical(data_all)
U_pca <- cov_pca(data_strong, npc = min(3, K - 1))
U_ed  <- tryCatch(
  cov_ed(data_strong, U_pca),
  error = function(e) {
    cat("cov_ed failed:", conditionMessage(e),
        "\nFalling back to canonical-only U.\n")
    NULL
  })
U <- if (is.null(U_ed)) U_c else c(U_c, U_ed)
cat("Using", length(U), "covariance matrices\n")

# --- Fit mash ---
cat("Fitting mash on", nrow(B), "genes ×", K, "cohorts...\n")
m <- mash(data_all, U)

# --- Posterior summaries ---
post_mean <- get_pm(m)
post_lfsr <- get_lfsr(m)

# Lock matrix orientation: mash preserves the column order from Bhat (which we
# built with dimnames=...,mega_cohorts). This assertion guards the long-format
# flatten below — if mashr ever reorders columns, the pivot would silently break.
stopifnot(identical(colnames(post_mean), mega_cohorts))
stopifnot(identical(colnames(post_lfsr), mega_cohorts))

# --- Sharing classes ---
n_sig <- rowSums(post_lfsr < 0.05)
sig_signs <- sign(post_mean) * (post_lfsr < 0.05)
same_sign <- apply(sig_signs, 1, function(x) {
  s <- x[x != 0]
  if (length(s) == 0) TRUE else length(unique(s)) == 1
})

majority_thresh <- floor(K / 2) + 1
sharing_class <- ifelse(n_sig == 0,                              "null",
                ifelse(n_sig == K & same_sign,                   "pan_cohort",
                ifelse(n_sig >= majority_thresh & same_sign,     "majority_shared",
                ifelse(n_sig == 1,                               "cohort_specific",
                                                                 "divergent"))))

# --- Pack outputs ---
gene_ids <- rownames(post_mean)
posterior_dt <- data.table(
  gene       = rep(gene_ids, each = K),
  cohort     = rep(mega_cohorts, times = length(gene_ids)),
  post_logFC = as.vector(t(post_mean)),
  lfsr       = as.vector(t(post_lfsr))
)

sharing_dt <- data.table(
  gene               = gene_ids,
  n_sig_cohorts      = n_sig,
  sharing_class      = factor(sharing_class,
                              levels = c("pan_cohort","majority_shared",
                                         "cohort_specific","divergent","null")),
  pan_cohort         = sharing_class == "pan_cohort",
  max_lfsr_among_sig = apply(post_lfsr, 1, function(x) {
                          s <- x[x < 0.05]; if (length(s)) max(s) else NA_real_
                       }),
  max_abs_logFC      = apply(abs(post_mean), 1, max),
  min_abs_logFC      = apply(abs(post_mean), 1, min)
)

cat("\n===== mashr sharing-class breakdown =====\n")
print(table(sharing_class, useNA = "ifany"))
cat("Total genes:", length(gene_ids), "\n")

fwrite(posterior_dt, file.path(RDIR_OUT, "posterior_per_cohort.csv"))
fwrite(sharing_dt,   file.path(RDIR_OUT, "sharing_classes.csv"))
fwrite(as.data.table(post_lfsr, keep.rownames = "gene"),
       file.path(RDIR_OUT, "lfsr_matrix.csv"))
saveRDS(m, file.path(RDIR_OUT, "mash_fit.rds"))

cat("\nSaved:", file.path(RDIR_OUT, "sharing_classes.csv"), "\n")
