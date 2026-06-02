#!/usr/bin/env Rscript
##############################################################################
# fig3_cross_ancestry_classifier.R  (new 2026-04-15)
# Fig 3 panel: Cross-ancestry COLOC genetic elastic-net classifier.
#
# Train an elastic-net logistic classifier on the cross-ancestry COLOC panel
# (EUR + EAS validated genes) to predict disease status (case/control).
# Benchmark against:
#   - Random gene panels of equal size (permutation null)
#   - EUR-only COLOC panel
#
# Features: gene-level expression on the validated panel; standardized.
# CV: leave-one-cohort-out (same split as Script 192).
#
# Output: figures/main/fig3_regulatory_architecture/panels/
#         fig3_cross_ancestry_classifier.pdf (ROC + AUROC comparison)
#
# Compute note: single-fold EN on 10 cohorts; fits on CPU in <10 min.
# Submit via run_fig3_regulatory.sbatch (NOT a login-node job).
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(ggplot2)
  library(pROC)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

dir.create(file.path(FIG3_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)
OUT   <- file.path(FIG3_DIR, "panels", "fig3_cross_ancestry_classifier.pdf")
OUTCSV <- file.path(FIG3_DIR, "panels", "fig3_cross_ancestry_classifier_metrics.csv")

# -----------------------------------------------------------------------------
# Inputs: cross-ancestry validated gene panel + expression matrix + metadata
# -----------------------------------------------------------------------------
CA <- fread(file.path(BASE,
  "RNA-seq/results/causal_inference/cross_ancestry/cross_ancestry_validated_targets.csv"))
panel_xa  <- unique(CA[n_validated_ancestry >= 2, gene])     # multi-ancestry
panel_eur <- unique(CA[best_pp4_eur >= 0.9, gene])           # EUR-only
message("Cross-ancestry panel size: ", length(panel_xa))
message("EUR-only panel size:       ", length(panel_eur))

# Normalized expression + metadata via atlas helper
EXPR_RDS <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration",
  "corrected_logcpm.rds")
META_CSV <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata",
  "unified_metadata.csv")
stopifnot(file.exists(EXPR_RDS), file.exists(META_CSV))

expr <- readRDS(EXPR_RDS)      # genes × samples, ENSG IDs as rownames
meta <- fread(META_CSV)
meta <- meta[sample_id %in% colnames(expr)]
meta <- meta[!is.na(group_binary)]
expr <- expr[, meta$sample_id]
y    <- as.integer(meta$group_binary == "Disease")

# Map ENSG rownames to HGNC symbols so the COLOC panels (HGNC) can be intersected.
# Use load_gene_map() from load_figure_data.R; strip version suffix to match atlas.
gene_map <- load_gene_map()                                   # cols: ensembl_clean, symbol
rownames(expr) <- sub("\\..*$", "", rownames(expr))
lookup <- setNames(gene_map$symbol, gene_map$ensembl_clean)
sym    <- lookup[rownames(expr)]
keep   <- !is.na(sym) & !duplicated(sym)
expr   <- expr[keep, , drop = FALSE]
rownames(expr) <- sym[keep]
message("Expression matrix after symbol mapping: ",
        nrow(expr), " genes × ", ncol(expr), " samples")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
loco_auc <- function(features) {
  present <- intersect(features, rownames(expr))
  if (length(present) < 3) return(NA_real_)
  X  <- t(expr[present, , drop = FALSE])
  folds <- split(seq_along(y), meta$dataset)
  preds <- rep(NA_real_, length(y))
  for (ds in names(folds)) {
    te <- folds[[ds]]; tr <- setdiff(seq_along(y), te)
    if (length(unique(y[tr])) < 2) next
    fit <- cv.glmnet(X[tr, , drop = FALSE], y[tr],
                     family = "binomial", alpha = 0.5, nfolds = 5)
    preds[te] <- as.numeric(predict(fit, X[te, , drop = FALSE],
                                    s = "lambda.min", type = "response"))
  }
  ok <- !is.na(preds)
  if (sum(ok) < 20 || length(unique(y[ok])) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(y[ok], preds[ok], quiet = TRUE)))
}

# -----------------------------------------------------------------------------
# Evaluate panels
# -----------------------------------------------------------------------------
auc_xa  <- loco_auc(panel_xa)
auc_eur <- loco_auc(panel_eur)

n_rand  <- 100
auc_rand <- replicate(n_rand, {
  pool <- rownames(expr)
  s <- sample(pool, size = length(panel_xa))
  loco_auc(s)
})
auc_rand <- auc_rand[!is.na(auc_rand)]
p_emp <- mean(auc_rand >= auc_xa)

metrics <- data.table(
  panel  = c("Cross-ancestry (EUR+EAS)", "EUR-only",
             sprintf("Random (n=%d)", length(auc_rand))),
  n_genes = c(length(panel_xa), length(panel_eur), length(panel_xa)),
  auroc  = c(auc_xa, auc_eur, mean(auc_rand)),
  auroc_sd = c(NA_real_, NA_real_, sd(auc_rand)),
  p_vs_random = c(p_emp, NA_real_, NA_real_))
fwrite(metrics, OUTCSV)

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
plot_dt <- data.table(
  auc = c(auc_xa, auc_eur, auc_rand),
  panel = c("Cross-ancestry", "EUR-only", rep("Random", length(auc_rand))))
plot_dt[, panel := factor(panel, levels = c("Cross-ancestry", "EUR-only", "Random"))]

p <- ggplot(plot_dt, aes(x = panel, y = auc, color = panel)) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.35, linewidth = 0.5,
               color = "black") +
  scale_color_manual(values = c("Cross-ancestry" = "#1b7837",
                                "EUR-only"       = "#4575b4",
                                "Random"         = "grey60"),
                     guide = "none") +
  labs(x = NULL, y = "AUROC (leave-one-cohort-out)",
       title = "Cross-ancestry COLOC panel elastic-net classifier",
       subtitle = sprintf("Cross-ancestry AUROC = %.3f (p_emp vs random = %.3f, n = %d genes)",
                          auc_xa, p_emp, length(panel_xa))) +
  theme_masld()

ggsave(OUT, p, width = 6.0, height = 4.0, device = cairo_pdf)
message("Wrote: ", OUT)
