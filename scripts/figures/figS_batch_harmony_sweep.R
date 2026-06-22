#!/usr/bin/env Rscript
# figS_batch_harmony_sweep.R — Harmony parameter sweep over theta to test
# whether more aggressive batch correction can improve cohort iLISI without
# destroying disease biology.
#
# Sweep:
#   theta in {0 (no correction), 0.5, 1, 2 (current default), 4, 8, 16, 32}
#   max.iter.harmony = 30 (more than the default 10/20 to ensure convergence)
#
# Per-theta metrics on the corrected 30-D embedding (k=30 nearest neighbors):
#   - Cohort iLISI (max = 10)              -- higher = better cross-cohort mixing
#   - Disease iLISI (max = 2)              -- LOWER = cleaner Control/Disease separation
#   - Disease AUROC (logistic, 5-fold CV)  -- HIGHER = biology preserved
#   - Avg disease silhouette               -- HIGHER = biology preserved
#   - Cohort silhouette                    -- LOWER = better mixing
#
# Outputs:
#   figS_batch_harmony_sweep_metrics.csv      per-theta summary
#   figS_batch_harmony_sweep_tradeoff.pdf     bivariate trade-off plot
#   figS_batch_harmony_sweep_umaps.pdf        UMAP grid (one row per theta x cohort/disease)
#   figS_batch_harmony_sweep_perTheta_*.csv   per-sample iLISI export per theta

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(harmony)
  library(uwot)
  library(RANN)
  library(cluster)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- file.path(FIGS_BATCH_DIR, "panels")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

COHORT_LABEL <- c(
  GSE126848   = "GSE126848",   GSE130970   = "GSE130970",
  GSE135251   = "GSE135251",   GSE162694   = "GSE162694",
  GSE167523   = "GSE167523",   GSE174478   = "GSE174478",
  GSE193066   = "GSE193066",   GSE213621   = "GSE213621",
  GSE240729   = "GSE240729",   PRJNA512027 = "PRJNA512027"
)
COHORT_COLORS <- c(
  "GSE126848"="#1F77B4","GSE130970"="#FF7F0E","GSE135251"="#2CA02C","GSE162694"="#D62728",
  "GSE167523"="#9467BD","GSE174478"="#8C564B","GSE193066"="#E377C2","GSE213621"="#7F7F7F",
  "GSE240729"="#BCBD22","PRJNA512027"="#17BECF"
)

KNN_K   <- 30
N_PCS   <- 30
N_HVG   <- 2000
THETAS  <- c(0, 0.5, 1, 2, 4, 8, 16, 32)
MAX_IT  <- 30

# ----------------------------------------------------------------------------
# Load + pre-process
# ----------------------------------------------------------------------------
message("Loading merged DGE...")
dge <- load_merged_dge()
n_samples <- ncol(dge)
n_genes   <- nrow(dge)
message(sprintf("DGE: %d samples x %d genes", n_samples, n_genes))

logcpm <- cpm(dge, log = TRUE, prior.count = 1)
vars   <- matrixStats::rowVars(logcpm)
top_g  <- order(vars, decreasing = TRUE)[seq_len(min(N_HVG, length(vars)))]

message(sprintf("Computing PCA on top %d variable genes...", N_HVG))
pca   <- prcomp(t(logcpm[top_g, ]), scale. = TRUE, center = TRUE)
pcs0  <- pca$x[, seq_len(N_PCS)]

batch    <- as.character(dge$samples$dataset)
disease  <- factor(dge$samples$group_binary, levels = c("Control", "Disease"))
cohort_f <- factor(COHORT_LABEL[batch], levels = unname(COHORT_LABEL))

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
ilisi <- function(coords, labels, k = KNN_K) {
  nn  <- RANN::nn2(coords, query = coords, k = k + 1)$nn.idx[, -1]
  lab <- as.integer(labels)
  apply(nn, 1, function(idx) {
    p <- tabulate(lab[idx]) / k
    1 / sum(p^2)
  })
}

cv_auroc <- function(X, y, k_fold = 5, seed = 42) {
  # Simple CV logistic regression -> AUROC, on standardized columns
  set.seed(seed)
  n <- nrow(X); folds <- sample(rep_len(seq_len(k_fold), n))
  preds <- numeric(n)
  for (f in seq_len(k_fold)) {
    tr <- folds != f; te <- !tr
    fit <- suppressWarnings(glm.fit(cbind(1, X[tr, , drop = FALSE]),
                                    as.integer(y[tr]) - 1L,
                                    family = binomial()))
    coefs <- fit$coefficients; coefs[is.na(coefs)] <- 0
    preds[te] <- as.vector(cbind(1, X[te, , drop = FALSE]) %*% coefs)
  }
  # AUROC
  ord <- order(preds, decreasing = TRUE)
  yo  <- as.integer(y)[ord] - 1L
  pos <- sum(yo); neg <- length(yo) - pos
  if (pos == 0 || neg == 0) return(NA_real_)
  cum_pos <- cumsum(yo); cum_neg <- cumsum(1 - yo)
  tpr <- cum_pos / pos; fpr <- cum_neg / neg
  sum(diff(c(0, fpr)) * (c(0, tpr[-length(tpr)]) + tpr) / 2)
}

avg_silhouette <- function(coords, labels, max_n = 1500, seed = 42) {
  # subsample for speed (silhouette is O(n^2))
  set.seed(seed)
  if (nrow(coords) > max_n) {
    idx <- sample.int(nrow(coords), max_n)
    coords <- coords[idx, , drop = FALSE]; labels <- labels[idx]
  }
  d  <- dist(coords)
  cl <- as.integer(labels)
  if (length(unique(cl)) < 2) return(NA_real_)
  s  <- cluster::silhouette(cl, d)
  mean(s[, "sil_width"])
}

# ----------------------------------------------------------------------------
# Sweep
# ----------------------------------------------------------------------------
metrics <- list()
embeds  <- list()  # store UMAP coords for plotting

for (th in THETAS) {
  tag <- if (th == 0) "raw" else sprintf("theta%g", th)
  message(sprintf("\n=== theta = %s ===", tag))

  if (th == 0) {
    H <- pcs0
  } else {
    H <- HarmonyMatrix(
      data_mat = pcs0,
      meta_data = data.frame(batch = batch),
      vars_use = "batch",
      do_pca = FALSE,
      theta = th,
      max.iter.harmony = MAX_IT,
      verbose = FALSE
    )
  }

  message("  iLISI cohort/disease...")
  li_c <- ilisi(H, cohort_f,  KNN_K)
  li_d <- ilisi(H, disease,   KNN_K)

  message("  Disease classification AUROC (5-fold CV logistic on PCs)...")
  auc <- cv_auroc(H, disease, k_fold = 5)

  message("  Silhouettes (subsampled)...")
  sil_d <- avg_silhouette(H, disease,  max_n = 1500)
  sil_c <- avg_silhouette(H, cohort_f, max_n = 1500)

  message("  UMAP for visualization...")
  set.seed(42)
  um <- uwot::umap(H, n_neighbors = 30, min_dist = 0.3, metric = "cosine",
                   verbose = FALSE)

  metrics[[tag]] <- data.table(
    theta             = th,
    cohort_ilisi_med  = median(li_c),
    cohort_ilisi_iqr  = IQR(li_c),
    disease_ilisi_med = median(li_d),
    disease_ilisi_iqr = IQR(li_d),
    disease_auroc     = auc,
    silhouette_disease= sil_d,
    silhouette_cohort = sil_c
  )
  embeds[[tag]] <- data.table(
    UMAP1 = um[, 1], UMAP2 = um[, 2],
    cohort = cohort_f, disease = disease,
    theta = th, theta_lab = tag
  )

  # per-sample export
  fwrite(data.table(sample_id = colnames(dge),
                    theta = th, cohort = cohort_f, disease = disease,
                    cohort_ilisi = li_c, disease_ilisi = li_d),
         file.path(OUT_DIR, sprintf("figS_batch_harmony_sweep_perTheta_%s.csv", tag)))
}

mt <- rbindlist(metrics)
fwrite(mt, file.path(OUT_DIR, "figS_batch_harmony_sweep_metrics.csv"))
print(mt)

# ----------------------------------------------------------------------------
# Plot 1: trade-off curve
# ----------------------------------------------------------------------------
mt[, theta_lab := ifelse(theta == 0, "raw (no correction)", sprintf("theta = %g", theta))]
mt[, theta_lab := factor(theta_lab, levels = c("raw (no correction)",
                                              sprintf("theta = %g", THETAS[THETAS > 0])))]

p_a <- ggplot(mt, aes(x = cohort_ilisi_med, y = disease_ilisi_med)) +
  geom_path(color = "grey60", linewidth = 0.4) +
  geom_point(aes(color = theta_lab), size = 3) +
  ggrepel::geom_text_repel(aes(label = theta_lab), size = 2.4,
                          family = "Helvetica", max.overlaps = Inf,
                          segment.size = 0.2) +
  scale_color_viridis_d(option = "plasma", end = 0.9, name = NULL, guide = "none") +
  scale_x_continuous(limits = c(1, 10),
                     breaks = c(1, 2, 4, 6, 8, 10)) +
  scale_y_continuous(limits = c(1, 2),
                     breaks = c(1, 1.25, 1.5, 1.75, 2)) +
  labs(x = "Cohort iLISI (max 10; higher = better cohort mixing)",
       y = "Disease iLISI (max 2; lower = cleaner biology)",
       title = "Harmony theta sweep: cohort mixing vs biology preservation",
       subtitle = "Top-left = good (cohorts mixed, disease/control preserved). Top-right = biology destroyed.") +
  theme_masld()

p_b <- ggplot(mt, aes(x = cohort_ilisi_med, y = disease_auroc)) +
  geom_path(color = "grey60", linewidth = 0.4) +
  geom_point(aes(color = theta_lab), size = 3) +
  ggrepel::geom_text_repel(aes(label = theta_lab), size = 2.4,
                          family = "Helvetica", max.overlaps = Inf,
                          segment.size = 0.2) +
  scale_color_viridis_d(option = "plasma", end = 0.9, name = NULL, guide = "none") +
  scale_x_continuous(limits = c(1, 10),
                     breaks = c(1, 2, 4, 6, 8, 10)) +
  scale_y_continuous(limits = c(0.5, 1),
                     breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  labs(x = "Cohort iLISI (higher = better cohort mixing)",
       y = "Disease vs control AUROC (logistic on 30 PCs, 5-fold CV)",
       title = "Disease classifiability after batch correction",
       subtitle = "Drop in AUROC at high theta = batch correction over-corrects biology.") +
  theme_masld()

p_c <- ggplot(mt, aes(x = cohort_ilisi_med, y = silhouette_disease)) +
  geom_path(color = "grey60", linewidth = 0.4) +
  geom_point(aes(color = theta_lab), size = 3) +
  ggrepel::geom_text_repel(aes(label = theta_lab), size = 2.4,
                          family = "Helvetica", max.overlaps = Inf,
                          segment.size = 0.2) +
  scale_color_viridis_d(option = "plasma", end = 0.9, name = NULL, guide = "none") +
  geom_hline(yintercept = 0, color = "grey80", linetype = "dashed", linewidth = 0.3) +
  labs(x = "Cohort iLISI (higher = better cohort mixing)",
       y = "Mean silhouette width (disease label)",
       title = "Disease silhouette after batch correction",
       subtitle = "Negative silhouette = disease & control inseparable in the embedding.") +
  theme_masld()

tradeoff <- p_a / p_b / p_c +
  plot_annotation(tag_levels = list(c("a", "b", "c"))) &
  theme(plot.tag = element_text(face = "bold", size = 10))

save_fig(tradeoff,
         file.path(OUT_DIR, "figS_batch_harmony_sweep_tradeoff.pdf"),
         width = 6.5, height = 12)

# ----------------------------------------------------------------------------
# Plot 2: UMAP grid (cohort/disease columns x theta rows)
# ----------------------------------------------------------------------------
all_em <- rbindlist(embeds)
all_em[, theta_lab := factor(
  ifelse(theta == 0, "raw (no correction)", sprintf("theta = %g", theta)),
  levels = c("raw (no correction)", sprintf("theta = %g", THETAS[THETAS > 0]))
)]
mt_lookup <- mt[, .(theta_lab, cohort_ilisi_med, disease_ilisi_med, disease_auroc)]
all_em <- merge(all_em, mt_lookup, by = "theta_lab")

# Add iLISI annotation as facet strip text
all_em[, cohort_strip := sprintf("%s | cohort iLISI = %.2f", theta_lab, cohort_ilisi_med)]
all_em[, disease_strip := sprintf("%s | disease iLISI = %.2f", theta_lab, disease_ilisi_med)]

p_um_cohort <- ggplot(all_em, aes(UMAP1, UMAP2, color = cohort)) +
  rasterize_layer(geom_point(size = 0.18, alpha = 0.7, shape = 16)) +
  scale_color_manual(values = COHORT_COLORS, name = NULL) +
  facet_wrap(~ cohort_strip, ncol = 2, scales = "free") +
  guides(color = guide_legend(override.aes = list(size = 1.6, alpha = 1),
                              ncol = 2)) +
  labs(title = "UMAP by cohort across theta sweep") +
  theme_masld() +
  theme(strip.text = element_text(size = 7),
        axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 7))

p_um_disease <- ggplot(all_em, aes(UMAP1, UMAP2, color = disease)) +
  rasterize_layer(geom_point(size = 0.18, alpha = 0.7, shape = 16)) +
  scale_color_manual(values = c(Control = masld_colors$control,
                                Disease = masld_colors$nash), name = NULL) +
  facet_wrap(~ disease_strip, ncol = 2, scales = "free") +
  guides(color = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
  labs(title = "UMAP by disease state across theta sweep") +
  theme_masld() +
  theme(strip.text = element_text(size = 7),
        axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 7))

save_fig(p_um_cohort,
         file.path(OUT_DIR, "figS_batch_harmony_sweep_umaps_cohort.pdf"),
         width = 8, height = 13)
save_fig(p_um_disease,
         file.path(OUT_DIR, "figS_batch_harmony_sweep_umaps_disease.pdf"),
         width = 8, height = 13)

message("\nWrote:")
message("  ", file.path(OUT_DIR, "figS_batch_harmony_sweep_metrics.csv"))
message("  ", file.path(OUT_DIR, "figS_batch_harmony_sweep_tradeoff.pdf"))
message("  ", file.path(OUT_DIR, "figS_batch_harmony_sweep_umaps_cohort.pdf"))
message("  ", file.path(OUT_DIR, "figS_batch_harmony_sweep_umaps_disease.pdf"))
message("\nDone.")
