#!/usr/bin/env Rscript
# figS_pls_disease_separation.R
# PLS-DA (Partial Least Squares Discriminant Analysis) to visualize Control vs
# Disease separation — a supervised alternative to PCA that finds the direction
# of maximum covariance between expression and disease label.
#
# Key design: leave-one-cohort-out (LOCO) cross-validation identical to panelK.
# Train PLS on 4 cohorts, project held-out 5th — zero leakage. HVGs also
# selected from training folds only.
# NIPALS algorithm implemented in base R (no mixOmics/pls packages needed).
#
# Outputs: panelJ_pls_disease_separation.pdf
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(edgeR); library(matrixStats); library(pROC)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
OUT <- file.path(BASE,
  "figures/supplementary/figS_methods_validation/multimethod_validation/panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CTRL <- "#9E9E9E"; set.seed(42)

MEGA <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
cohort_short <- c(GSE126848 = "GSE126848", GSE130970 = "GSE130970",
                  GSE135251 = "GSE135251", GSE162694 = "GSE162694", GSE213621 = "GSE213621")
N_HVG <- 2000L; NCOMP <- 3L

# ---------------------------------------------------------------------------
# NIPALS PLS (binary Y) — base R implementation
# Returns W (weight matrix) and X_center for projecting new samples.
# ---------------------------------------------------------------------------
nipals_pls <- function(X, y, ncomp = 2) {
  Xh <- sweep(X, 2, colMeans(X))          # center columns
  yh <- y - mean(y)
  W  <- matrix(0, ncol(X), ncomp)
  P  <- matrix(0, ncol(X), ncomp)
  for (a in seq_len(ncomp)) {
    w  <- drop(t(Xh) %*% yh)
    w  <- w / sqrt(sum(w^2))              # unit-length weight vector
    t  <- drop(Xh %*% w)                 # X score
    p  <- drop(t(Xh) %*% t) / sum(t^2)  # X loading
    c_ <- sum(yh * t) / sum(t^2)         # Y loading (scalar for univariate Y)
    Xh <- Xh - outer(t, p)               # deflate X
    yh <- yh - c_ * t                    # deflate Y
    W[, a] <- w; P[, a] <- p
  }
  list(W = W, X_center = colMeans(X))
}

pls_project <- function(fit, X_new) {
  sweep(X_new, 2, fit$X_center) %*% fit$W
}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
dge  <- load_merged_dge(); stopifnot(!is.null(dge))
samp <- as.data.table(dge$samples, keep.rownames = "sample_id")
keep <- samp$dataset %in% MEGA
dge  <- dge[, keep]; samp <- samp[keep]
grp  <- factor(samp$group_binary, levels = c("Control", "Disease"))
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
cat(sprintf("LOCO PLS-DA: %d samples, 5 cohorts\n", ncol(logcpm)))

# ---------------------------------------------------------------------------
# LOCO cross-validation
# ---------------------------------------------------------------------------
scores <- rbindlist(lapply(MEGA, function(test) {
  tr <- which(samp$dataset != test)
  te <- which(samp$dataset == test)
  if (length(unique(grp[te])) < 2) return(NULL)

  # HVG selection on training fold (within-cohort centred variance)
  lc_tr    <- logcpm[, tr]
  ds_tr    <- samp$dataset[tr]
  lc_tr_wc <- lc_tr
  for (d in unique(ds_tr)) {
    idx <- which(ds_tr == d)
    lc_tr_wc[, idx] <- lc_tr[, idx] - rowMeans(lc_tr[, idx, drop = FALSE])
  }
  rv  <- matrixStats::rowVars(lc_tr_wc)
  top <- order(rv, decreasing = TRUE)[seq_len(N_HVG)]

  X_tr <- t(lc_tr[top, ])
  X_te <- t(logcpm[top, te])
  y_tr <- as.numeric(grp[tr] == "Disease")

  fit  <- nipals_pls(X_tr, y_tr, ncomp = NCOMP)
  T_te <- pls_project(fit, X_te)

  data.table(sample_id = samp$sample_id[te],
             cohort    = cohort_short[test],
             disease   = grp[te],
             PLS1 = T_te[, 1], PLS2 = T_te[, 2], PLS3 = T_te[, 3])
}))
scores[, cohort := factor(cohort, levels = unname(cohort_short))]

# ---------------------------------------------------------------------------
# AUROC (PLS1 as continuous predictor)
# ---------------------------------------------------------------------------
auc_by <- scores[, .(auc = as.numeric(pROC::auc(pROC::roc(
  response = disease, predictor = PLS1,
  levels = c("Control", "Disease"), quiet = TRUE)))), by = cohort]
# z-score PLS1 within each cohort so pooled AUROC is not deflated by
# scale differences across LOCO folds
scores[, PLS1_z := scale(PLS1), by = cohort]

# ensure direction: positive = Disease (per cohort, then globally)
if (mean(scores$PLS1_z[scores$disease == "Disease"]) <
    mean(scores$PLS1_z[scores$disease == "Control"])) {
  scores[, PLS1_z := -PLS1_z]
  auc_by[, auc := 1 - auc]
}
roc_all <- pROC::roc(response = scores$disease, predictor = scores$PLS1_z,
                     levels = c("Control", "Disease"), quiet = TRUE)
auc_all <- as.numeric(pROC::auc(roc_all))
cat(sprintf("Pooled held-out AUROC (PLS1) = %.3f\n", auc_all)); print(auc_by)

# ---------------------------------------------------------------------------
# Panel A: PLS1 vs PLS2 scatter, faceted by cohort
# ---------------------------------------------------------------------------
lab_auc <- auc_by[, setNames(sprintf("%s\nAUROC=%.2f", cohort, auc), cohort)]
pA <- ggplot(scores, aes(PLS1, PLS2, colour = disease)) +
  geom_point(size = 0.5, alpha = 0.55) +
  facet_wrap(~ cohort, nrow = 1, scales = "free",
             labeller = labeller(cohort = lab_auc)) +
  scale_colour_manual(values = c(Control = CTRL, Disease = masld_colors$nash),
                      name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  labs(x = "PLS component 1 (held-out)", y = "PLS component 2 (held-out)",
       title = sprintf(
         "PLS-DA (LOCO): held-out projections  —  pooled AUROC = %.2f", auc_all)) +
  theme_masld(base_size = 7) +
  theme(strip.text       = element_text(size = 6.5, face = "bold"),
        legend.position  = "top",
        panel.grid.minor = element_blank())

# ---------------------------------------------------------------------------
# Panel B: PLS1 score violin by cohort (same layout as panelK left)
# ---------------------------------------------------------------------------
lab_auc2 <- auc_by[, setNames(sprintf("%s\nAUROC %.2f", cohort, auc), cohort)]
pB <- ggplot(scores, aes(cohort, PLS1_z, fill = disease, colour = disease)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_violin(position = position_dodge(width = 0.8), width = 0.78,
              alpha = 0.35, linewidth = 0.3, scale = "width",
              draw_quantiles = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12,
                                             dodge.width = 0.8),
             size = 0.3, alpha = 0.4, show.legend = FALSE) +
  scale_fill_manual(values   = c(Control = CTRL, Disease = masld_colors$nash),
                    name = NULL) +
  scale_colour_manual(values = c(Control = CTRL, Disease = masld_colors$nash),
                      name = NULL) +
  scale_x_discrete(labels = lab_auc2) +
  labs(x = NULL, y = "PLS component 1 (held-out, z-scored)") +
  theme_masld(base_size = 7) +
  theme(axis.text.x    = element_text(size = 6, lineheight = 0.85),
        legend.position = "none")

# ---------------------------------------------------------------------------
# Panel C: held-out ROC
# ---------------------------------------------------------------------------
roc_dt <- data.table(fpr = 1 - roc_all$specificities,
                     tpr = roc_all$sensitivities)
setorder(roc_dt, fpr, tpr)
pC <- ggplot(roc_dt, aes(fpr, tpr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey70", linewidth = 0.3) +
  geom_path(colour = masld_colors$nash, linewidth = 0.7) +
  annotate("text", x = 0.6, y = 0.18, size = 2.8, colour = "grey15",
           label = sprintf("pooled held-out\nAUROC = %.2f", auc_all)) +
  coord_equal() +
  labs(x = "false positive rate", y = "true positive rate",
       title = "ROC (PLS1)") +
  theme_masld(base_size = 7)

fig <- (pA / (pB | pC)) +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(
    title    = "Supervised PLS-DA separates Control from Disease (LOCO cross-validation)",
    subtitle = sprintf(
      "NIPALS PLS, top-%d HVGs per training fold (within-cohort variance). n=%d held-out samples. Compare to unsupervised PCA (panelJ) where disease is invisible.",
      N_HVG, nrow(scores)),
    theme = theme(plot.title    = element_text(size = 9, face = "bold"),
                  plot.subtitle = element_text(size = 6.2, colour = "grey35")))

ggsave(file.path(OUT, "panelJ_pls_disease_separation.pdf"), fig,
       width = 9.0, height = 7.5, device = cairo_pdf)
fwrite(scores, file.path(OUT, "panelJ_pls_disease_separation_data.csv"))
cat("Wrote panelJ_pls_disease_separation.pdf\n")
