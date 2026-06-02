#!/usr/bin/env Rscript
# 03_metafor_REML_HKSJ.R
# ---------------------------------------------------------------------------
# Mega-validation Arm 3: per-gene REML random-effects meta-analysis with
# Hartung-Knapp-Sidik-Jonkman small-K correction (K = number of mega cohorts).
#
# Inputs:  per-cohort limma-voom CSVs from 02_per_study_de.R (must have
#          logFC, SE, df.total, gene, dataset)
# Outputs: per-gene tau^2, I^2, Q, HKSJ-corrected meta logFC + p + padj
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(metafor); library(data.table); library(yaml); library(BiocParallel)
})

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
PERSTUDY <- file.path(INT, "results/per_study")
RDIR_OUT <- file.path(INT, "results/mega_validation/metafor")
dir.create(RDIR_OUT, recursive = TRUE, showWarnings = FALSE)

ycfg <- read_yaml(file.path(PROJECT_ROOT, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
K <- length(mega_cohorts)
cat("Mega cohorts (K =", K, "):", paste(mega_cohorts, collapse=", "), "\n")

# --- Load per-cohort DE results ---
per_study <- rbindlist(lapply(mega_cohorts, function(ds) {
  p <- file.path(PERSTUDY, paste0(ds, "_de_results.csv"))
  if (!file.exists(p)) stop("Missing per-study DE for ", ds,
                            " — re-run 02_per_study_de.R")
  dt <- fread(p)
  stopifnot("SE" %in% names(dt), "df.total" %in% names(dt))
  dt
}), fill = TRUE)
cat("Per-study rows:", nrow(per_study), "\n")

# --- Pivot to gene × cohort matrices ---
genes <- sort(unique(per_study$gene))
B <- matrix(NA_real_, nrow = length(genes), ncol = K,
            dimnames = list(genes, mega_cohorts))
S <- B; DF <- B
for (i in seq_len(K)) {
  ds <- mega_cohorts[i]
  sub <- per_study[dataset == ds]
  idx <- match(sub$gene, genes)
  B[idx, i]  <- sub$logFC
  S[idx, i]  <- sub$SE
  DF[idx, i] <- sub$df.total
}

# Require >= floor(K/2)+1 cohorts contributing per gene (3+ of 5 at K=5)
n_obs <- rowSums(!is.na(B))
min_K <- floor(K/2) + 1
keep <- n_obs >= min_K
cat("Genes meeting >=", min_K, "cohort coverage:", sum(keep), "/",
    length(keep), "\n")

# --- Per-gene REML+HKSJ via metafor ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
bp <- if (ncpus > 1) MulticoreParam(ncpus, progressbar = FALSE) else SerialParam()

run_rma <- function(i) {
  yi <- B[i, ]; sei <- S[i, ]
  ok <- !is.na(yi) & !is.na(sei) & sei > 0
  if (sum(ok) < min_K) return(NULL)
  tryCatch({
    r <- rma(yi = yi[ok], sei = sei[ok], method = "REML",
             test = "knha",
             control = list(stepadj = 0.5, maxiter = 1000))
    data.table(
      gene       = genes[i],
      meta_logFC = r$b[1, 1],
      meta_SE    = r$se,
      meta_z     = r$zval,           # under knha this is tval
      HKSJ_p     = r$pval,
      tau2       = r$tau2,
      I2         = r$I2,
      Q          = r$QE,
      Q_p        = r$QEp,
      n_studies  = sum(ok)
    )
  }, error = function(e) NULL)
}

idx_keep <- which(keep)
cat("Running REML+HKSJ on", length(idx_keep), "genes...\n")
results_list <- bplapply(idx_keep, run_rma, BPPARAM = bp)
res_dt <- rbindlist(Filter(Negate(is.null), results_list))
conv_pct <- 100 * nrow(res_dt) / length(idx_keep)
cat("Convergence:", round(conv_pct, 2), "%\n")

res_dt[, HKSJ_padj := p.adjust(HKSJ_p, method = "BH")]

cat("\n===== metafor REML+HKSJ RESULTS =====\n")
cat("Genes meta-analyzed:", nrow(res_dt), "\n")
cat("DEGs HKSJ_padj < 0.05:", res_dt[HKSJ_padj < 0.05, .N], "\n")
cat("Median tau^2:", round(median(res_dt$tau2), 4), "\n")
cat("Median I^2:",   round(median(res_dt$I2),   2), "%\n")
cat("Genes with I^2 > 75%:", res_dt[I2 > 75, .N], "\n")

fwrite(res_dt, file.path(RDIR_OUT, "per_gene_meta.csv"))

het <- data.table(
  metric = c("n_genes_meta", "n_DEG_HKSJ_padj005", "median_tau2",
             "median_I2_pct", "n_genes_I2_gt75", "convergence_pct",
             "K", "min_K_required"),
  value  = c(nrow(res_dt), res_dt[HKSJ_padj < 0.05, .N],
             round(median(res_dt$tau2), 4),
             round(median(res_dt$I2), 2),
             res_dt[I2 > 75, .N],
             round(conv_pct, 2), K, min_K)
)
fwrite(het, file.path(RDIR_OUT, "heterogeneity_summary.csv"))
