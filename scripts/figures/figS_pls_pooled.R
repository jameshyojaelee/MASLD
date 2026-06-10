#!/usr/bin/env Rscript
# figS_pls_pooled.R
# PLS-DA on all 5 cohorts pooled (batch-corrected), mirroring the layout of
# panelJ_batch_model_pca.pdf. Rows: cohort coloured | disease coloured.
# Columns: PCA (unsupervised) vs PLS-DA (supervised) — direct comparison.
# Note: PLS here is in-sample (no held-out), so it is a visualization of the
# supervised axis, not a cross-validated classifier. See panelK/panelJ_pls_*
# for held-out performance numbers.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(edgeR); library(limma); library(matrixStats)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
OUT <- file.path(BASE,
  "figures/supplementary/figS_methods_validation/multimethod_validation/panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CTRL <- "#9E9E9E"

MEGA <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
cohort_short <- c(GSE126848 = "Suppli", GSE130970 = "Hoang",
                  GSE135251 = "Govaere", GSE162694 = "Bril", GSE213621 = "Chen")
cohort_pal   <- c(Suppli = "#1F77B4", Hoang = "#FF7F0E", Govaere = "#2CA02C",
                  Bril = "#D62728", Chen = "#9467BD")
N_HVG <- 2000L; NCOMP <- 2L

# ---------------------------------------------------------------------------
# NIPALS PLS (base R, no extra packages)
# ---------------------------------------------------------------------------
nipals_pls <- function(X, y, ncomp = 2) {
  Xh <- sweep(X, 2, colMeans(X))
  yh <- y - mean(y)
  W  <- matrix(0, ncol(X), ncomp)
  P  <- matrix(0, ncol(X), ncomp)
  for (a in seq_len(ncomp)) {
    w  <- drop(t(Xh) %*% yh); w <- w / sqrt(sum(w^2))
    t  <- drop(Xh %*% w)
    p  <- drop(t(Xh) %*% t) / sum(t^2)
    c_ <- sum(yh * t) / sum(t^2)
    Xh <- Xh - outer(t, p); yh <- yh - c_ * t
    W[, a] <- w; P[, a] <- p
  }
  list(W = W, X_center = colMeans(X))
}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
dge  <- load_merged_dge(); stopifnot(!is.null(dge))
samp <- as.data.table(dge$samples, keep.rownames = "sample_id")
keep <- samp$dataset %in% MEGA
dge  <- dge[, keep]; samp <- samp[keep]
grp  <- factor(samp$group_binary, levels = c("Control", "Disease"))
ds   <- factor(samp$dataset)
cat(sprintf("Pooled: %d samples\n", ncol(dge)))

logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

# HVG selection — within-cohort centred (same as panelJ)
logcpm_wc <- logcpm
for (d in unique(samp$dataset)) {
  idx <- which(samp$dataset == d)
  logcpm_wc[, idx] <- logcpm[, idx] - rowMeans(logcpm[, idx, drop = FALSE])
}
rv  <- matrixStats::rowVars(logcpm_wc)
top <- order(rv, decreasing = TRUE)[seq_len(N_HVG)]
X   <- logcpm[top, ]

# batch-corrected matrix (same as panelJ)
X_corr <- limma::removeBatchEffect(X, batch = ds,
                                   design = model.matrix(~ grp))

# ---------------------------------------------------------------------------
# PCA (batch-corrected)
# ---------------------------------------------------------------------------
pr_pca <- prcomp(t(X_corr), center = TRUE, scale. = FALSE)
pve    <- round(100 * pr_pca$sdev^2 / sum(pr_pca$sdev^2), 1)
pca_dt <- data.table(samp[, .(sample_id, dataset, group_binary)],
                     Dim1 = pr_pca$x[, 1], Dim2 = pr_pca$x[, 2],
                     method = "PCA (unsupervised)",
                     xlab = sprintf("PC1 (%.1f%% var)", pve[1]),
                     ylab = sprintf("PC2 (%.1f%% var)", pve[2]))

# ---------------------------------------------------------------------------
# PLS-DA (batch-corrected, in-sample)
# ---------------------------------------------------------------------------
y    <- as.numeric(grp == "Disease")
fit  <- nipals_pls(t(X_corr), y, ncomp = NCOMP)
T_pls <- sweep(t(X_corr), 2, fit$X_center) %*% fit$W

# flip so Disease is positive on PLS1
if (mean(T_pls[grp == "Disease", 1]) < mean(T_pls[grp == "Control", 1]))
  T_pls[, 1] <- -T_pls[, 1]

pls_dt <- data.table(samp[, .(sample_id, dataset, group_binary)],
                     Dim1 = T_pls[, 1], Dim2 = T_pls[, 2],
                     method = "PLS-DA (supervised)",
                     xlab = "PLS component 1",
                     ylab = "PLS component 2")

both <- rbind(pca_dt, pls_dt)
both[, cohort  := factor(cohort_short[dataset], levels = unname(cohort_short))]
both[, disease := factor(group_binary, levels = c("Control", "Disease"))]
both[, method  := factor(method, levels = c("PCA (unsupervised)",
                                             "PLS-DA (supervised)"))]

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
base_theme <- function() {
  theme_masld(base_size = 7) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          strip.text = element_text(size = 7, face = "bold"),
          strip.background = element_blank(),
          legend.position = "right",
          legend.title = element_text(size = 6.5, face = "bold"),
          legend.text  = element_text(size = 6),
          legend.key.size = unit(0.25, "cm"))
}

# row 1 — coloured by cohort
row_cohort <- ggplot(both, aes(Dim1, Dim2, colour = cohort)) +
  geom_point(size = 0.4, alpha = 0.7) +
  facet_wrap(~ method, nrow = 1, scales = "free") +
  scale_colour_manual(values = cohort_pal, name = "Cohort") +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  labs(x = NULL, y = NULL,
       title = "Coloured by cohort — PCA mixes cohorts after correction; PLS focuses on disease") +
  base_theme()

# row 2 — coloured by disease
row_dis <- ggplot(both, aes(Dim1, Dim2, colour = disease)) +
  geom_point(size = 0.4, alpha = 0.7) +
  facet_wrap(~ method, nrow = 1, scales = "free") +
  scale_colour_manual(values = c(Control = CTRL, Disease = masld_colors$nash),
                      name = "Disease") +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  labs(x = NULL, y = NULL,
       title = "Coloured by disease — PLS-DA explicitly separates Control from Disease on axis 1") +
  base_theme()

fig <- (row_cohort / row_dis) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title    = "PCA (unsupervised) vs PLS-DA (supervised) — batch-corrected, 5 cohorts pooled",
    subtitle = sprintf(
      "Top-%d HVGs, within-cohort variance selection. Batch removed via limma::removeBatchEffect. PLS-DA is in-sample (not cross-validated) — see panelK for held-out AUROC.",
      N_HVG),
    theme = theme(plot.title    = element_text(size = 9, face = "bold"),
                  plot.subtitle = element_text(size = 6.2, colour = "grey35")))

ggsave(file.path(OUT, "panelJ_pls_pooled.pdf"), fig,
       width = 8.0, height = 5.8, device = cairo_pdf)
fwrite(both[, .(sample_id, dataset, cohort, disease, method, Dim1, Dim2)],
       file.path(OUT, "panelJ_pls_pooled_data.csv"))
cat("Wrote panelJ_pls_pooled.pdf\n")
